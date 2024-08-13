// ---------------------------------------------------------------
//
// History:
//
// Created by Evgueni Goudzovski (eg@hep.ph.bham.ac.uk) 2016-03-25
// Extended functionality by Lubos Bician (lubos.bician@cern.ch) 2023-07-12
//
// ---------------------------------------------------------------

/// \class MyKmu2Selection
/// \ingroup ph_analyzers gr_physics
/// \Brief
/// Kmu2 decay selection
/// \EndBrief
/// \Detailed
/// A simple Kmu2 decay selection producing basic monitoring plots:
/// missing mass, momentum, vertex position, etc.
/// By default, the main selection does not use the GTK information.
/// However, an additional selection (which is stricter than it could be),
/// relying on kaon momentum measurement with the GTK,
/// is employed on top of the standard selection for GTK performance monitoring.
/// Both selections output their acceptances at step 2. <br>
/// Update on 2023-07-12: <br>
/// Using parameter UseGTK (false by default), one can force the
/// selection to provide outputs obtained with using the GTK kaon.
/// Parameter TightPID (true by default) can be used to enable/disable
/// the E/p and MUV3 association cuts used for positive muon PID.
/// Both control and physics L0 data types are considered.
/// For the physics data, the L0 trigger mask to be used can be specified
/// from the command line (0xFF by default).
/// For example, to select events with L0 trigger bit 5 up, one can call
/// \code
/// ./MyApplication ... -p "MyKmu2Selection:TriggerMask=0x20"
/// \endcode
/// or equivalently, using decimal notation,
/// \code
/// ./MyApplication ... -p "MyKmu2Selection:TriggerMask=32"
/// \endcode
/// \EndDetailed
/// \author Evgueni Goudzovski (eg@hep.ph.bham.ac.uk)

#include "MyKmu2Selection.hh"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "KaonDecayConstants.hh"
#include "LAVMatching.hh"
#include "Persistency.hh"
#include "SAVMatching.hh"
#include "TriggerConditions.hh"

#include <TLatex.h>

#include <iostream>
#include <stdlib.h>

using namespace NA62Analysis;

MyKmu2Selection::MyKmu2Selection(Core::BaseAnalysis *ba):
    Analyzer(ba, "MyKmu2Selection"), fHZtrue(nullptr), fHPhysicsEventsPerBurst(nullptr),
    fHPhysicsEventsGoodQualityPerBurst(nullptr), fHPhysicsEventsGoodQualityGTKPerBurst(nullptr),
    fHPhysicsEventsGoodQualitywGTKPerBurst(nullptr), fHKmu2EventsPerBurst(nullptr), fHMass(nullptr),
    fHEOP(nullptr), fHZvtx(nullptr), fHZvtxGTK(nullptr), fHxyMUV3(nullptr), fHMassVsP(nullptr) 
{
  RequestTree("Cedar", new TRecoCedarEvent, "Reco");
  RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
  RequestTree("CHOD", new TRecoCHODEvent, "Reco");
  RequestTree("RICH", new TRecoRICHEvent, "Reco");
  RequestTree("LKr", new TRecoLKrEvent, "Reco");
  RequestTree("LAV", new TRecoLAVEvent, "Reco");
  RequestTree("IRC", new TRecoIRCEvent, "Reco");
  RequestTree("SAC", new TRecoSACEvent, "Reco");
  RequestL0Data();

  // Acceptances of the two selections, Ztrue range (105-180)m, MC v2.1.6/v2.2.2 non-overlaid
  fAcceptance = 0.4413;
  fAcceptanceGTK = 0.3388;

  fCDAcomp = new TwoLinesCDA();
  fMatchingRG = new MatchingRG(ba, this, "MatchingRG");

  AddParam("TriggerMask", &fTriggerMask, 0xFF);
  AddParam("TightPID", &fTightPID, true);
  AddParam("UseGTK", &fUseGTK, false);
  AddParam("MaxNBursts", &fMaxNBursts, 5000);  // max number of bins in histograms
  AddParam("SkipWrongType", &fSkipWrongType, false);
}

void MyKmu2Selection::StartOfBurstUser() {
}

void MyKmu2Selection::EndOfBurstUser() {
}

void MyKmu2Selection::InitHist() {
  fReadingData = GetIsTree();

  if(fReadingData) {
    // Histograms of MC true quantities
    if(GetWithMC()) {
      BookHisto("mctrue/hZvertex", new TH1F("Zvertex_true", "True Zvertex; z [m]", 300, 0, 300));
    }

    BookHisto("hNTracks", new TH1F("hNTracks", "Number of tracks", 11, -0.5, 10.5));
    BookHisto("hNGoodTracks", new TH1F("hNGoodTracks", "Number of good tracks", 11, -0.5, 10.5));
    BookHisto("TrackTime", new TH1F("TrackTime", "TrackTime", 100, -25, 25));
    BookHisto("TrackNchambers", new TH1F("TrackNchambers", "TrackNchambers", 5, -0.5, 4.5));
    BookHisto("QChi2Track", new TH1F("QChi2Track", "QChi2Track", 100, -100, 100));

    BookHisto("hEoP", new TH1F("EoP", "Track E/p; E/p", 150, 0.0, 1.5));
    BookHisto("hdTTrackCedar", new TH1F("dTTrackCedar", "dTTrackCedar;#Deltat [ns]", 100, -25, 25));
    BookHisto("hZvtx", new TH1F("Zvtx", "Z of track-beam axis vertex;Vertex z [m]", 200, 50, 250));
    BookHisto("hZvtxSelected", new TH1F("ZvtxSelected", "Z of track-beam axis vertex;Vertex z [m]", 200, 50, 250));
    BookHisto("hZvtxSelected_WithGTK", new TH1F("ZvtxSelected_WithGTK", "Z of track-beam axis vertex;Vertex z [m]", 200, 50, 250));
    BookHisto("hCDA", new TH1F("CDA", "CDA of the track-beam axis vertex;CDA [mm]", 200, 0, 200));
    BookHisto("hZvtxCDA", new TH2F("ZvtxCDA", "CDA vs Z of the track-beam axis vertex;Vertex z [m];CDA [mm]", 75, 50, 200, 100, 0, 200));
    BookHisto("hMMiss2Mu", new TH1F("MMiss2Mu",
                                    "Squared missing mass in muon hypothesis (beam average "
                                    "momentum);M_{miss}^{2}(#mu) [GeV^{2}/c^{4}]",
                                    300, -0.15, 0.15));
    BookHisto("hPMMiss2Mu", new TH2F("PMMiss2Mu",
                                     "Squared missing mass in muon hypothesis vs momentum; Track "
                                     "momentum [GeV/c];M_{miss}^{2}(#mu) [GeV^{2}/c^{4}]",
                                     160, 0, 80, 100, -0.15, 0.15));
    BookHisto("hMMiss2Mu_WithGTK", new TH1F( "MMiss2Mu_WithGTK",
        "Squared missing mass in muon hypothesis (GTK track);M_{miss}^{2}(#mu) [GeV^{2}/c^{4}]",
        300, -0.15, 0.15));
    BookHisto("hTrackXYMUV3", new TH2F("TrackXYMUV3", "Track (x,y) at MUV3 plane;x [mm];y [mm]", 65,
                                       -1300, 1300, 65, -1300, 1300));

    BookHisto("hPhysicsEventsPerBurst",
              new TH1F("PhysicsEventsPerBurst", "Physics events per burst;Burst ID", fMaxNBursts,
                       -0.5, fMaxNBursts - 0.5));

    BookHisto("hPhysicsEventsGoodQualityPerBurst",
              new TH1F("PhysicsEventsGoodQualityPerBurst",
                       "Physics events with Good Event Quality Mask per burst ; Burst ID",
                       fMaxNBursts, -0.5, fMaxNBursts - 0.5));
    BookHisto("hPhysicsEventsGoodQualityGTKPerBurst",
              new TH1F("PhysicsEventsGoodQualityGTKPerBurst",
                       "Physics events with no GTK Error Mask per burst ; Burst ID", fMaxNBursts,
                       -0.5, fMaxNBursts - 0.5));
    BookHisto(
      "hPhysicsEventsGoodQualitywGTKPerBurst",
      new TH1F(
        "PhysicsEventsGoodQualitywGTKPerBurst",
        "Physics events with Good Event Quality mask and No GTK Error Mask per burst ; Burst ID",
        fMaxNBursts, -0.5, fMaxNBursts - 0.5));

    BookHisto("hKmu2EventsPerBurst",
              new TH1F("Kmu2EventsPerBurst",
                       "Kmu2 candidates per burst (beam average momentum);Burst ID", fMaxNBursts,
                       -0.5, fMaxNBursts - 0.5));
    BookHisto("hKmu2Events_WithGTK_PerBurst",
              new TH1F("Kmu2Events_WithGTK_PerBurst",
                       "Kmu2 candidates per burst (GTK track);Burst ID", fMaxNBursts, -0.5,
                       fMaxNBursts - 0.5));

    BookHisto("hNCHODHits",
              new TH1F("NCHODHits", "Number of CHOD hits;Number of hits", 50, -0.5, 49.5));
    BookHisto("hNRICHHits",
              new TH1F("NRICHHits", "Number of RICH hits;Number of hits", 100, -0.5, 99.5));
    BookHisto("hNLKrCells",
              new TH1F("NLKrCells", "Number of LKr cells with any signal;Number of cells", 100,
                       -0.5, 99.5));
    BookHisto("hNLKrGoodCells",
              new TH1F("NLKrGoodCells", "Number of LKr cells with E>40MeV;Number of cells", 100,
                       -0.5, 99.5));
    BookHisto("hNLKrGoodCellsIn",
              new TH1F("NLKrGoodCellsIn", "Number of LKr cells with E>40MeV;Number of cells", 100,
                       -0.5, 99.5));
    BookHisto("hNLKrGoodCellsOut",
              new TH1F("NLKrGoodCellsOut", "Number of LKr cells with E>40MeV;Number of cells", 100,
                       -0.5, 99.5));
    BookHisto(
      "hCellTrackDistance",
      new TH1F("CellTrackDistance", "Good cell - track distance;Distance [mm]", 200, 0, 2000));
    BookHisto("hNLKrClusters",
              new TH1F("NLKrClusters", "Number of LKr clusters;Number of clusters", 10, -0.5, 9.5));
    BookHisto("hLKrClusterEnergy",
              new TH1F("LKrClusterEnergy", "LKr cluster energy;Energy [GeV]", 100, 0, 50));
    BookHisto("hLKrCellTotalEnergy",
              new TH1F("LKrCellTotalEnergy",
                       "LKr total cell (E>40MeV) energy;Total cell energy [GeV]", 70, 0, 70));
    BookHisto("hLKrCellTotalEnergyIn",
              new TH1F("LKrCellTotalEnergyIn",
                       "LKr total cell (E>40MeV) energy near track;Total cell energy [GeV]", 70, 0,
                       70));
    BookHisto("hLKrCellTotalEnergyOut",
              new TH1F("LKrCellTotalEnergyOut",
                       "LKr total cell (E>40MeV) energy far from track;Total cell energy [GeV]", 70,
                       0, 70));
    BookHisto(
      "hLKrCellTotalEnergyInOut",
      new TH2F("LKrCellTotalEnergyInOut",
               "LKr total cell energies;Energy near track [GeV];Energy far from track [GeV]", 100,
               0, 50, 100, 0, 50));

    BookHisto("hLKrCellClusterTotalEnergy",
              new TH2F("LKrCellClusterTotalEnergy",
                       "LKr total cluster energy vs cell energy;Total cell (>40MeV) energy "
                       "[GeV];Total cluster energy [GeV]",
                       70, 0, 70, 70, 0, 70));
  }
  else {  // step 2
    std::cout << user_normal() << "Reading my own output" << std::endl;
    fHZtrue = static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "mctrue/Zvertex_true", true));
    fHPhysicsEventsPerBurst =
      static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "PhysicsEventsPerBurst", true));
    fHPhysicsEventsGoodQualityPerBurst = static_cast<TH1F *>(
      RequestHistogram(fAnalyzerName, "PhysicsEventsGoodQualityPerBurst", true));
    fHPhysicsEventsGoodQualityGTKPerBurst = static_cast<TH1F *>(
      RequestHistogram(fAnalyzerName, "PhysicsEventsGoodQualityGTKPerBurst", true));
    fHPhysicsEventsGoodQualitywGTKPerBurst = static_cast<TH1F *>(
      RequestHistogram(fAnalyzerName, "PhysicsEventsGoodQualitywGTKPerBurst", true));
    fHKmu2EventsPerBurst =
      static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "Kmu2EventsPerBurst", true));
    fHMass = static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "MMiss2Mu", true));
    fHEOP = static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "EoP", true));
    fHZvtx = static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "ZvtxSelected", true));
    fHZvtxGTK = static_cast<TH1F *>(RequestHistogram(fAnalyzerName, "ZvtxSelected_WithGTK", true));
    fHxyMUV3 = static_cast<TH2F *>(RequestHistogram(fAnalyzerName, "TrackXYMUV3", true));
    fHMassVsP = static_cast<TH2F *>(RequestHistogram(fAnalyzerName, "PMMiss2Mu", true));
  }
}

void MyKmu2Selection::InitOutput() {
  RegisterOutput("EventSelected", &fEventSelected);
  RegisterOutput("MuonMomentum", &fMuonMomentum);
  RegisterOutput("MissingMass2", &fMissingMass2);
  RegisterOutput("VertexZ", &fVertexZ);
  RegisterOutput("GTKTrackID", &fGTKTrackID);
  RegisterOutput("Kmu2TrackID", &fKmu2TrackID);
}

void MyKmu2Selection::Process(Int_t) {

  // Initialize the outputs
  SetOutputState("EventSelected", kOValid);
  SetOutputState("MuonMomentum", kOValid);
  SetOutputState("Kmu2TrackID", kOValid);
  fEventSelected = false;
  fMuonMomentum = 0.0;
  fMissingMass2 = 9999.0;
  fVertexZ = 0.0;
  fKmu2TrackID = 0;
  fGTKTrackID = 0;

  if(!fReadingData)
    return;  // no action if reading its own output in --histo mode
  if(GetWithMC() && fSkipWrongType && GetMCEvent()->GetEventBoundary(0)->GetStreamID() % 1000 != 31)
    return;  // not a Kmu2 event

  Bool_t PhysicsData = TriggerConditions::GetInstance()->IsPhysicsTrigger(GetL0Data());
  Int_t L0TriggerWord = GetL0Data()->GetTriggerFlags();
  Bool_t TriggerOK = TriggerConditions::GetInstance()->IsControlTrigger(GetL0Data());
  Int_t TrigCtrl1 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-Q1");
  if(GetRunID() > 10000)  // Run 2
    TriggerOK = TriggerConditions::GetInstance()->L0TriggerOn(GetRunID(), GetL0Data(), TrigCtrl1);
  TriggerOK |= (PhysicsData && (L0TriggerWord & fTriggerMask));
  if(!TriggerOK)
    return;

  FillHisto("hPhysicsEventsPerBurst", GetBurstID());

  if(GetWithMC()) {
    Event *evt = GetMCEvent();
    EventBoundary *evt0 = static_cast<EventBoundary *>(evt->GetEventBoundary(0));
    if(evt0->GetNKineParts()) {
      FillHisto("mctrue/hZvertex", 0.001 * evt->GetKinePart(0)->GetEndPos().Z());  // [m]
    }
  }

  TRecoCedarEvent *CEDARevent = GetEvent<TRecoCedarEvent>();
  TRecoGigaTrackerEvent *GTKevent = GetEvent<TRecoGigaTrackerEvent>();
  TRecoCHODEvent *CHODevent = GetEvent<TRecoCHODEvent>();
  TRecoRICHEvent *RICHevent = GetEvent<TRecoRICHEvent>();
  TRecoLKrEvent *LKRevent = GetEvent<TRecoLKrEvent>();
  TRecoLAVEvent *LAVevent = GetEvent<TRecoLAVEvent>();
  TRecoIRCEvent *IRCevent = GetEvent<TRecoIRCEvent>();
  TRecoSACEvent *SACevent = GetEvent<TRecoSACEvent>();

  // Check impact of event quality checks
  if(GetEventHeader()->GetEventQualityMask() == 0) {
    FillHisto("hPhysicsEventsGoodQualityPerBurst", GetBurstID());
    if(!GTKevent->GetErrorMask())
      FillHisto("hPhysicsEventsGoodQualitywGTKPerBurst", GetBurstID());
  }
  if(!GTKevent->GetErrorMask())
    FillHisto("hPhysicsEventsGoodQualityGTKPerBurst", GetBurstID());

  // This is the trigger time
  Double_t RefTime = GetEventHeader()->GetFineTime() * TriggerCalib;

  //////////////////////////
  // Find a good STRAW track

  std::vector<DownstreamTrack> Tracks = *GetOutput<std::vector<DownstreamTrack>>("DownstreamTrackBuilder.Output");
  FillHisto("hNTracks", Tracks.size());
  if(Tracks.size() > 10)
    return;

  Bool_t TrackIsGood[10];
  for(UInt_t i = 0; i < Tracks.size(); i++) {
    TrackIsGood[i] = false;

    if(!Tracks[i].CHODTimeExists() && !Tracks[i].NewCHODAssociationExists())
      continue;

    Double_t TrackTime = (Tracks[i].CHODTimeExists()) ? Tracks[i].GetCHODTime() : Tracks[i].GetNewCHODTime();
    TrackTime -= RefTime;  // with respect to the trigger time
    FillHisto("TrackTime", TrackTime);
    if(fabs(TrackTime) > 10.0)
      continue;  // out of time with the trigger

    FillHisto("TrackNchambers", Tracks[i].GetNChambers());
    if(Tracks[i].GetNChambers() != 4)
      continue;
    FillHisto("QChi2Track", Tracks[i].GetCharge() * Tracks[i].GetChi2());
    if(Tracks[i].GetChi2() > 20.0)
      continue;

    TrackIsGood[i] = true;

    // Forbid tracks forming two-track vertices with other non-fake tracks
    fCDAcomp->SetLine1PointDir(TVector3(Tracks[i].xAt(0), Tracks[i].yAt(0), 0.0),
                               Tracks[i].GetMomentumBeforeMagnet());
    for(UInt_t j = 0; j < Tracks.size(); j++) {
      if(i == j)
        continue;
      if(Tracks[j].GetIsFake())
        continue;
      fCDAcomp->SetLine2PointDir(TVector3(Tracks[j].xAt(0), Tracks[j].yAt(0), 0.0),
                                 Tracks[j].GetMomentumBeforeMagnet());
      fCDAcomp->ComputeVertexCDA();
      double cv = fCDAcomp->GetCDA();
      double zv = 1e-3 * fCDAcomp->GetVertex().Z();
      if(cv < 50.0 && zv > 60.0 && zv < 200.0)
        TrackIsGood[i] = false;
    }
  }

  // Require exactly one good track
  Int_t itr = -1, NGoodTracks = 0;
  for(UInt_t i = 0; i < Tracks.size(); i++) {
    if(TrackIsGood[i]) {
      NGoodTracks++;
      itr = i;
    }
  }
  FillHisto("hNGoodTracks", NGoodTracks);
  if(NGoodTracks != 1)
    return;

  // Selection conditions based on the track properties
  Int_t Q = Tracks[itr].GetCharge();
  Double_t Ptrack = Tracks[itr].GetMomentum();  // spectrometer calibration included
  Double_t Ttrack = Tracks[itr].GetCHODTime();
  Bool_t UseCHODTime = true;
  if(!Tracks[itr].CHODTimeExists()) {
    Ttrack = Tracks[itr].GetNewCHODTime();
    UseCHODTime = false;
  }
  if(Q != 1 || Ptrack < 5000 || Ptrack > 70000)
    return;
  TRecoSpectrometerCandidate *Scand = Tracks[itr].GetSpectrometerCandidate();
  if(!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 0))
    return;
  if(!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 1))
    return;
  if(!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 2))
    return;
  if(!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 3))
    return;
  if(!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kMUV3))
    return;

  Double_t eop = Tracks[itr].GetLKrEoP();
  FillHisto("hEoP", eop);

  ////////////////////////////////////////////////////////////////////
  // Track-Cedar timing: look for the closest Cedar candidate to track

  Double_t dT_Track_Cedar = 999.999;
  Double_t CedarTime = -999.999;
  for(Int_t i = 0; i < CEDARevent->GetNCandidates(); i++) {
    TRecoCedarCandidate *Ccand = static_cast<TRecoCedarCandidate *>(CEDARevent->GetCandidate(i));
    if(Ccand->GetNSectors() < 5)
      continue;
    Double_t dT = Ttrack - Ccand->GetTime();
    if(fabs(dT) < fabs(dT_Track_Cedar)) {
      CedarTime = Ccand->GetTime();
      dT_Track_Cedar = dT;
    }
  }
  FillHisto("hdTTrackCedar", dT_Track_Cedar);
  if((UseCHODTime && fabs(dT_Track_Cedar) > 2.0) || (!UseCHODTime && fabs(dT_Track_Cedar) > 5.0))
    return;

  /////////////////////////////////////////
  // Zvertex & CDA: track wrt the beam axis

  Double_t cda_rundep = Tracks[itr].GetBeamAxisCDA();
  Double_t Zvtx_rundep = Tracks[itr].GetBeamAxisVertex().Z();

  TVector3 KaonThreeMomentum = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  TLorentzVector Kaon;
  Kaon.SetVectM(KaonThreeMomentum, MKCH);
  TLorentzVector Muon;
  Muon.SetVectM(Scand->GetThreeMomentumBeforeMagnet(), MMU);
  Double_t Mmiss2Mu_rundep = (Kaon - Muon).M2();

  ////////////////////////////////////
  // Zvertex & CDA: track wrt GTK kaon

  fMatchingRG->Process(GTKevent, Tracks[itr].GetSpectrometerCandidate(), CedarTime, CedarTime,
                       Ttrack, 0, "");                    // ref time, Cedar time, RICH time
  fMatchingRG->FinalSelection(CedarTime, Ttrack, 0, "");  // Cedar time, RICH time

  std::vector<Int_t> matchedGTKIDs = fMatchingRG->GetMatchedGTKIDs();
  Bool_t MatchedRG = (matchedGTKIDs.size() > 0 && matchedGTKIDs.at(0) != -1);
  Double_t Mmiss2Mu_GTK{9999.}, Zvtx_GTK{0.}, cda_GTK{9999.};
  if(fUseGTK && !MatchedRG)
    return;
  if(MatchedRG) {
    std::vector<TVector3> matchedGTKMomenta = fMatchingRG->GetGTKMomentaAtVertices();
    std::vector<TVector3> matchedSTRAWMomenta = fMatchingRG->GetTrackMomentaAtVertices();
    std::vector<TVector3> matchedVertices = fMatchingRG->GetVertices();
    std::vector<double> matchedCDA = fMatchingRG->GetCDA();
    Zvtx_GTK = matchedVertices.at(0).Z();
    cda_GTK = matchedCDA.at(0);
    TLorentzVector KaonRG;  // beam kaon
    KaonRG.SetVectM(matchedGTKMomenta.at(0), MKCH);
    TLorentzVector MuonRG;  // downstream muon
    MuonRG.SetVectM(matchedSTRAWMomenta.at(0), MMU);
    Mmiss2Mu_GTK = (KaonRG - MuonRG).M2();
  }

  //////////////////////////////////////////////////
  // Set observables based on the selection settings

  Double_t Zvtx = (fUseGTK) ? Zvtx_GTK : Zvtx_rundep;
  Double_t cda = (fUseGTK) ? cda_GTK : cda_rundep;
  Double_t Mmiss2Mu = (fUseGTK) ? Mmiss2Mu_GTK : Mmiss2Mu_rundep;
  Double_t Mmiss2Mu_cut = (fUseGTK) ? 0.01 : 0.02;
  Double_t PassedMmiss2MuCut = fabs(Mmiss2Mu * 1e-6) < Mmiss2Mu_cut;

  FillHisto("hZvtx", 0.001 * Zvtx);
  FillHisto("hCDA", cda);
  FillHisto("hZvtxCDA", 0.001 * Zvtx, cda);

  // cut on Zvtx and CDA
  if(Zvtx < 120000.0 || Zvtx > 180000.0 || cda > 40.0)
    return;

  /////////////////////////
  // LAV veto (with timing)

  LAVMatching *pLAVMatching = *GetOutput<LAVMatching *>("PhotonVetoHandler.LAVMatching");
  pLAVMatching->SetReferenceTime(CedarTime);
  if(pLAVMatching->LAVHasTimeMatching(LAVevent))
    return;

  /////////////////////////////////
  // IRC and SAC veto (with timing)

  SAVMatching *pSAVMatching = *GetOutput<SAVMatching *>("PhotonVetoHandler.SAVMatching");
  pSAVMatching->SetReferenceTime(CedarTime);
  pSAVMatching->SetIRCTimeCuts(10.0, 10.0);  // half time window; default = 5ns
  pSAVMatching->SetSACTimeCuts(10.0, 10.0);  // half time window; default = 5ns
  Bool_t SAVmatched = pSAVMatching->SAVHasTimeMatching(IRCevent, SACevent);
  if(SAVmatched)
    return;

  /////////////////////////////////////////////////////////////////////////
  // MUV3 trigger performance check: all cuts except Track-MUV3 association

  if((!fTightPID || eop < 0.2) && PassedMmiss2MuCut) {
    Double_t Zmuv3 = GeometricAcceptance::GetInstance()->GetZMUV3();
    FillHisto("hTrackXYMUV3", Tracks[itr].xAtAfterMagnet(Zmuv3), Tracks[itr].yAtAfterMagnet(Zmuv3));
  }

  /////////////////////////////////////////////////////////////////////
  // Muon-Cedar timing: looking for the closest muon in time with Cedar

  if(fTightPID) {
    if(!Tracks[itr].MUV3AssociationExists())
      return;

    Double_t dTmin = 9999;
    for(Int_t iMu = 0; iMu < Tracks[itr].GetNMUV3AssociationRecords(); iMu++) {
      Double_t dT = Tracks[itr].GetMUV3Time(iMu) - CedarTime;
      if(fabs(dT) < fabs(dTmin))
        dTmin = dT;
    }
    if(fabs(dTmin) > 1.5)
      return;

    // Track-MUV3 timing:
    Double_t dTtrkmin = 9999;
    for(Int_t iMu = 0; iMu < Tracks[itr].GetNMUV3AssociationRecords(); iMu++) {
      Double_t dT = Tracks[itr].GetMUV3Time(iMu) - Ttrack;
      if(fabs(dT) < fabs(dTtrkmin))
        dTtrkmin = dT;
    }
    if((UseCHODTime && fabs(dTtrkmin) > 2.) || (!UseCHODTime && fabs(dTtrkmin) > 5.))
      return;
  }

  ///////////////////
  // Muon ID with E/p

  if(fTightPID && eop > 0.2)
    return;

  FillHisto("hMMiss2Mu", Mmiss2Mu * 1e-6);  // [GeV^2]
  FillHisto("hPMMiss2Mu", 0.001 * Ptrack, Mmiss2Mu * 1e-6);
  if(PassedMmiss2MuCut)
    FillHisto("hKmu2EventsPerBurst", GetBurstID());

  ////////////////////////////////////////////////////////////////////
  // This is mainly to check that the GTK information looks plausible.

  if(MatchedRG) {
    FillHisto("hMMiss2Mu_WithGTK", Mmiss2Mu_GTK * 1e-6);  // [GeV^2]
    FillHisto("hKmu2Events_WithGTK_PerBurst", GetBurstID());
    if(fabs(Mmiss2Mu_GTK * 1e-6) < 0.01)
      FillHisto("hZvtxSelected_WithGTK", 0.001 * Zvtx_GTK);
  }

  if(!PassedMmiss2MuCut)
    return;

  FillHisto("hZvtxSelected", 0.001 * Zvtx);

  // CHOD response studies
  Int_t NCHODHits = CHODevent->GetNHits();
  FillHisto("hNCHODHits", NCHODHits);

  // RICH response studies
  Int_t NRICHHitsAll = RICHevent->GetNHits();  // including super-cells
  Int_t NRICHHits = 0;
  for(Int_t i = 0; i < NRICHHitsAll; i++) {
    TRecoRICHHit *hit = static_cast<TRecoRICHHit *>(RICHevent->GetHit(i));
    if(hit->GetOrSuperCellID() == 0)
      NRICHHits++;  // no super-cells
  }
  FillHisto("hNRICHHits", NRICHHits);

  // LKr response studies
  if(GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[itr], NA62::kLKr)) {
    Int_t NCells = LKRevent->GetNHits();
    Int_t NClusters = LKRevent->GetNCandidates();
    FillHisto("hNLKrCells", NCells);
    FillHisto("hNLKrClusters", NClusters);

    Double_t TotalCellEnergy = 0.0;
    Double_t TotalCellEnergy_in = 0.0;
    Double_t TotalCellEnergy_out = 0.0;
    Double_t TotalClusterEnergy = 0.0;
    Int_t NGoodCells = 0, NGoodCells_in = 0, NGoodCells_out = 0;
    for(Int_t i = 0; i < NCells; i++) {
      TRecoLKrHit *hit = static_cast<TRecoLKrHit *>(LKRevent->GetHit(i));
      Double_t energy = hit->GetEnergy();
      if(energy < 40.0)
        continue;
      NGoodCells++;
      TotalCellEnergy += energy;
      Double_t Zlkr = GeometricAcceptance::GetInstance()->GetZLKr();
      Double_t dx = hit->GetPosition().x() - Tracks[itr].xAt(Zlkr);
      Double_t dy = hit->GetPosition().y() - Tracks[itr].yAt(Zlkr);
      Double_t dist = sqrt(dx * dx + dy * dy);
      FillHisto("hCellTrackDistance", dist);
      if(dist < 80.0) {
        NGoodCells_in++;
        TotalCellEnergy_in += energy;
      }
      else {
        NGoodCells_out++;
        TotalCellEnergy_out += energy;
      }
    }
    FillHisto("hNLKrGoodCells", NGoodCells);
    FillHisto("hNLKrGoodCellsIn", NGoodCells_in);
    FillHisto("hNLKrGoodCellsOut", NGoodCells_out);
    for(Int_t i = 0; i < NClusters; i++) {
      TRecoLKrCandidate *Lcand = static_cast<TRecoLKrCandidate *>(LKRevent->GetCandidate(i));
      Double_t energy = Lcand->GetEnergy();
      FillHisto("hLKrClusterEnergy", 0.001 * energy);
      TotalClusterEnergy += energy;
    }
    FillHisto("hLKrCellTotalEnergy", 0.001 * TotalCellEnergy);
    FillHisto("hLKrCellTotalEnergyIn", 0.001 * TotalCellEnergy_in);
    FillHisto("hLKrCellTotalEnergyOut", 0.001 * TotalCellEnergy_out);
    FillHisto("hLKrCellTotalEnergyInOut", 0.001 * TotalCellEnergy_in, 0.001 * TotalCellEnergy_out);
    FillHisto("hLKrCellClusterTotalEnergy", 0.001 * TotalCellEnergy, 0.001 * TotalClusterEnergy);
  }

  // Save the outputs
  fEventSelected = true;
  fMuonMomentum = Ptrack;
  fMissingMass2 = Mmiss2Mu;
  fVertexZ = Zvtx;
  fGTKTrackID = (MatchedRG) ? matchedGTKIDs.at(0) : -1;
  fKmu2TrackID = itr;
}

void MyKmu2Selection::EndOfJobUser() {
  if(fReadingData) {  // Data mode: save output
    SaveAllPlots();
    return;
  }
  if(!fHPhysicsEventsPerBurst) {  // Histo mode required but no histograms found
    std::cout << user_normal() << "Asked to read my own output but cannot find it" << std::endl;
    return;
  }

  /////////////////////////////////////
  // Histo mode: analyze the histograms

  // Print out acceptance (this is used for revision metrics)
  if(GetWithMC() && fHZtrue) {
    Double_t n1 = fHZvtx->Integral();
    Double_t n2 = fHZvtxGTK->Integral();
    Double_t N = fHZtrue->Integral(106, 180);  // 105 m < Ztrue < 180 m
    fAcceptanceNew = n1 / N;
    fAcceptanceNewGTK = n2 / N;
    Double_t dAcc = sqrt(fAcceptanceNew * (1.0 - fAcceptanceNew) / N);
    Double_t dAccGTK = sqrt(fAcceptanceNewGTK * (1.0 - fAcceptanceNewGTK) / N);
    std::cout << user_normal() << Form("MC events read: %d\n", (Int_t)fHZtrue->Integral());
    std::cout << user_normal()
              << Form("MC acceptance = %d/%d = %7.5f +- %7.5f\n", (Int_t)n1, (Int_t)N,
                      fAcceptanceNew, dAcc);
    std::cout << user_normal()
              << Form("MC acceptance (with GTK) = %d/%d = %7.5f +- %7.5f\n", (Int_t)n2, (Int_t)N,
                      fAcceptanceNewGTK, dAccGTK);
    std::cout << user_normal()
              << Form("##CI_DASH::Acceptance.Kmu2=%7.5f+-%7.5f", fAcceptanceNew, dAcc) << std::endl;
    std::cout << user_normal()
              << Form("##CI_DASH::Acceptance.Kmu2GTK=%7.5f+-%7.5f", fAcceptanceNewGTK, dAccGTK)
              << std::endl;
  }

  BuildPDFReport();
}

void MyKmu2Selection::BuildPDFReport() {

  TString OutputPDFFileName = fAnalyzerName + ".pdf";
  gErrorIgnoreLevel = 5000;  // suppress messages generated for each page printed
  gStyle->SetOptStat(11);

  TCanvas *Canvas = new TCanvas("Kmu2Canvas");
  Canvas->Print(Form(OutputPDFFileName + "["), "pdf");  // open file

  Canvas->Divide(2, 2);
  for(Int_t i = 1; i <= 4; i++) {
    Canvas->GetPad(i)->SetLeftMargin(0.1);
    Canvas->GetPad(i)->SetRightMargin(0.03);
    Canvas->GetPad(i)->SetTopMargin(0.06);
    Canvas->GetPad(i)->SetBottomMargin(0.10);
  }
  fHMass->SetLineColor(kBlue);
  fHMass->SetFillColor(kYellow);
  fHEOP->SetLineColor(kBlue);
  fHEOP->SetFillColor(kYellow);
  fHMassVsP->SetMarkerColor(kBlue);

  Canvas->cd(1);
  gPad->SetLogy();
  fHEOP->Draw();
  Canvas->cd(2);
  fHMass->Fit("gaus", "Q", "", -0.01, 0.01);
  fHMass->Draw();
  Canvas->cd(3);
  fHMassVsP->Draw();
  Canvas->cd(4);
  Int_t MaxNonEmptyBurstID = 0;
  Int_t MaxY = -99;
  for(Int_t i = 0; i < fMaxNBursts; i++) {
    if(fHKmu2EventsPerBurst->GetBinContent(i) > 0) {
      MaxNonEmptyBurstID = i;
      if(fHKmu2EventsPerBurst->GetBinContent(i) > MaxY)
        MaxY = fHKmu2EventsPerBurst->GetBinContent(i);
    }
  }
  TH1F *Kmu2EventsPerBurst = new TH1F("Kmu2EventsPerBurst", "Kmu2 Events per Burst",
                                      MaxNonEmptyBurstID, -0.5, MaxNonEmptyBurstID - 0.5);
  for(Int_t iBin = 1; iBin <= fHKmu2EventsPerBurst->GetNbinsX(); iBin++) {
    Kmu2EventsPerBurst->SetBinContent(iBin, fHKmu2EventsPerBurst->GetBinContent(iBin));
  }
  Kmu2EventsPerBurst->SetLineColor(kBlue);
  Kmu2EventsPerBurst->SetFillColor(kYellow);
  Kmu2EventsPerBurst->SetAxisRange(0, MaxY + 30, "Y");
  Kmu2EventsPerBurst->GetXaxis()->SetTitle("Burst ID");
  Kmu2EventsPerBurst->Draw();

  Canvas->Print(OutputPDFFileName, "pdf");

  Canvas->Clear();
  Canvas->Divide(1, 2);
  for(Int_t i = 1; i <= 2; i++) {
    Canvas->GetPad(i)->SetLeftMargin(0.04);
    Canvas->GetPad(i)->SetRightMargin(0.01);
    Canvas->GetPad(i)->SetTopMargin(0.07);
    Canvas->GetPad(i)->SetBottomMargin(0.07);
  }
  fHZvtx->SetLineColor(kBlue);
  fHZvtx->SetFillColor(kYellow);

  Canvas->cd(1);
  fHZvtx->Draw();

  if(GetWithMC()) {  // Acceptance vs Zvertex plot
    Canvas->cd(2);
    TH1F *fHZtrue2 =
      new TH1F("Zvertex position for MC true",
               "Zvertex position for MC true with different Bin ranges", 200, 50, 250);
    for(Int_t iBin = 1; iBin <= fHZtrue->GetNbinsX(); iBin++) {
      fHZtrue2->SetBinContent(fHZtrue2->FindBin(iBin), fHZtrue->GetBinContent(iBin));
    }
    TH1F *fHAcceptanceZvtx = new TH1F("AcceptanceatZvtx", "Acceptance vs Z position", 200, 50, 250);
    fHAcceptanceZvtx->Divide(fHZvtx, fHZtrue2, 1., 1., "B");
    fHAcceptanceZvtx->GetXaxis()->SetTitle("Vertex z [m]");
    fHAcceptanceZvtx->SetAxisRange(0, 1.0, "Y");
    fHAcceptanceZvtx->Draw();
    TLatex text;
    text.SetTextSize(0.04);
    text.DrawLatex(105., .8, Form("Current Acceptance: %5.3f", fAcceptanceNew));
    text.SetNDC(kTRUE);
  }

  Canvas->Print(OutputPDFFileName, "pdf");

  Canvas->Clear();
  Canvas->GetPad(0)->SetLeftMargin(0.12);
  Canvas->GetPad(0)->SetRightMargin(0.04);
  Canvas->GetPad(0)->SetTopMargin(0.06);
  Canvas->GetPad(0)->SetBottomMargin(0.10);
  fHxyMUV3->SetMarkerColor(kBlue);
  fHxyMUV3->Draw("col");
  Canvas->Print(OutputPDFFileName, "pdf");

  // monitoring for event quality mask & GTK error mask
  fHPhysicsEventsPerBurst->SetLineColor(kBlack);              // all events
  fHPhysicsEventsGoodQualityPerBurst->SetLineColor(kBlue);    // no EventQualityMask (good events)
  fHPhysicsEventsGoodQualityGTKPerBurst->SetLineColor(kRed);  // no GTKErrorMask
  fHPhysicsEventsGoodQualitywGTKPerBurst->SetLineColor(
    kGreen + 2);  // no EventQualityMask and no GTKErrorMask

  fHPhysicsEventsPerBurst->SetTitle(" ; Burst ID ; Number of events");
  fHPhysicsEventsPerBurst->GetXaxis()->SetRangeUser(0., MaxNonEmptyBurstID);
  fHPhysicsEventsPerBurst->GetYaxis()->SetRangeUser(0.,
                                                    fHPhysicsEventsPerBurst->GetMaximum() * 1.2);
  fHPhysicsEventsPerBurst->SetStats(0);
  fHPhysicsEventsPerBurst->Draw();
  fHPhysicsEventsGoodQualityPerBurst->Draw("SAME");
  fHPhysicsEventsGoodQualityGTKPerBurst->Draw("SAME");
  fHPhysicsEventsGoodQualitywGTKPerBurst->Draw("SAME");

  TLegend *Legend = new TLegend(0.15, 0.90, 0.99, 0.99);
  Legend->SetNColumns(4);
  Legend->AddEntry(fHPhysicsEventsPerBurst, "All physics events", "l");
  Legend->AddEntry(fHPhysicsEventsGoodQualityPerBurst, "Good quality events", "l");
  Legend->AddEntry(fHPhysicsEventsGoodQualityGTKPerBurst, "GTK good quality events", "l");
  Legend->AddEntry(fHPhysicsEventsGoodQualitywGTKPerBurst, "Good quality (check all)", "l");
  Legend->Draw();

  Canvas->Print(OutputPDFFileName, "pdf");

  Canvas->Print(Form(OutputPDFFileName + "]"), "pdf");  // close file
  gErrorIgnoreLevel = -1;                               // restore the default

  delete Canvas;
  delete Legend;
  // PrintStatisticsPerBurst();
}

void MyKmu2Selection::PrintStatisticsPerBurst() {
  for(Int_t i = 1; i <= fHPhysicsEventsPerBurst->GetNbinsX(); i++) {
    Double_t N = fHPhysicsEventsPerBurst->GetBinContent(i);
    if(!N)
      continue;
    Double_t n = fHKmu2EventsPerBurst->GetBinContent(i);
    Double_t e = n / N;
    Double_t de = sqrt(e * (1.0 - e) / N);
    std::cout << user_standard() << "@@Kmu2 " << i - 1 << " " << n << " " << N << " " << e << " "
              << de << std::endl;
  }
}
