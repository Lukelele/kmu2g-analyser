#include "MyK3piSelection.hh"

#include "Event.hh"
#include "Persistency.hh"
#include "functions.hh"

#include <iostream>
#include <stdlib.h>



using namespace NA62Analysis;


MyK3piSelection::MyK3piSelection(Core::BaseAnalysis *ba) :
  Analyzer(ba, "MyK3piSelection"),
  fRunID(0),
  fBurstID(0)
{
  RequestTree("Cedar", new TRecoCedarEvent, "Reco");
  RequestTree("CHOD",  new TRecoCHODEvent,  "Reco");
  RequestTree("RICH",  new TRecoRICHEvent,  "Reco");
  RequestTree("LKr",   new TRecoLKrEvent,   "Reco");
  RequestTree("LAV",   new TRecoLAVEvent,   "Reco");
  RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
 
  RequestL0Data();
  RequestL1Data();
 
  fSGMatching = new SpectrometerGigaTrackerMatching();
  fSGMatching->SetMatchingTimingCuts(-0.5, 0.5);
}

void MyK3piSelection::InitOutput() {

}

void MyK3piSelection::InitHist(){
	if (GetWithMC()) { // Book those histograms only if the input file contains MC
    BookHisto("mctrue/hZvertex", new TH1F("Zvertex_true", "True Zvertex; z [m]", 300, 0, 300));
    BookHisto("mctrue/hMomentum", new TH1F("Momentum_true", "True kaon momentum; Momentum [GeV/c]", 60, 72, 78));
    BookHisto("mctrue/hTrackMomentum", new TH1F("TrackMomentum_true", "True track momentum;Momentum [GeV/c]", 50, 0, 50));
    BookHisto("mctrue/hTrackX1", new TH1F("TrackX1_true", "True track x at Straw1;x [mm]", 50, -500, 500));
    BookHisto("mctrue/hTrackY1", new TH1F("TrackY1_true", "True track y at Straw1;y [mm]", 50, -500, 500));
    BookHisto("mctrue/hdxdz", new TH1F("dxdz_true", "True kaon dx/dz", 100, -0.0025, 0.0025));
    BookHisto("mctrue/hdydz", new TH1F("dydz_true", "True kaon dy/dz", 100, -0.0025, 0.0025));
    BookHisto("mctrue/hXYstart", new TH2F("XYstart_true", "True kaon (x,y) at GTK3 plane (z=102.4m);x [mm]; y[mm]", 50, -50, 50, 50, -50, 50));
    BookHisto("mctrue/hXstart", new TH1F("Xstart_true", "True kaon x at GTK3 plane (z=102.4m);x [mm]", 100, -50, 50));
    BookHisto("mctrue/hYstart", new TH1F("Ystart_true", "True kaon y at GTK3 plane (z=102.4m);y [mm]", 100, -50, 50));
    BookHisto("mctrue/hdxdzVsXstart", new TH2F("dxdzVsX_true", "True dx/dz vs x; x [mm]; dx/dz", 20, -50, 50, 50, -0.0025, 0.0025));
    BookHisto("mctrue/hdydzVsXstart", new TH2F("dydzVsX_true", "True dy/dz vs x; x [mm]; dy/dz", 20, -50, 50, 50, -0.0025, 0.0025));
    BookHisto("mctrue/hdxdzVsYstart", new TH2F("mctrue/dxdzVsY_true", "True dx/dz vs y; y [mm]; dx/dz", 20, -25, 25, 50, -0.0025, 0.0025));
    BookHisto("mctrue/hdydzVsYstart", new TH2F("dydzVsY_true", "True dy/dz vs y; y [mm]; dy/dz", 20, -25, 25, 50, -0.0025, 0.0025));
    BookHisto("mctrue/pdxdzVsXstart", new TProfile("dxdzVsXP_true", "True dx/dz vs x; x [mm]; dx/dz", 20, -50, 50));
    BookHisto("mctrue/pdydzVsXstart", new TProfile("dydzVsXP_true", "True dy/dz vs x; x [mm]; dy/dz", 20, -50, 50));
    BookHisto("mctrue/pdxdzVsYstart", new TProfile("dxdzVsYP_true", "True dx/dz vs y; y [mm]; dx/dz", 20, -25, 25));
    BookHisto("mctrue/pdydzVsYstart", new TProfile("dydzVsYP_true", "True dy/dz vs y; y [mm]; dy/dz", 20, -25, 25));
 
    // M3pi vs the phi angle: true MC events, this is a pure blue field effect
    BookHisto("mctrue/hM3piVsPhi", new TH2F("M3piVsPhi_true", "True M(3#pi) vs #varphi: pure blue field effect; #varphi/#pi; M(3#pi) [MeV]", 40, -1, 1, 30, 492, 495));
 
    BookHisto("mctrue/CedarTimeMC", new TH1F("CedarTimeMC", "Cedar Time #minus Fine Time;Time [ns]", 200, -2, 2));
  }
 
  // Histograms of reconstructed quantities: general monitoring
  	BookHisto("general/hNTracks", new TH1F("NTracks", "Number of tracks", 11, -0.5, 10.5));
  	BookHisto("general/hNVertices", new TH1F("NVertices", "Number of 3-track vertices", 11, -0.5, 10.5));
  	BookHisto("general/hVertexTimeControl", new TH1F("VertexTimeControl", "Vertex CHOD time wrt trigger (control triggers);time [ns]", 100, -50, 50));
  	BookHisto("general/hVertexTimePhysics", new TH1F("VertexTimePhysics", "Vertex CHOD time wrt trigger (physics triggers);time [ns]", 100, -50, 50));
	
  	BookHisto("general/hTrackXY_straw1", new TH2F("TrackXY_straw1", "Track (x,y) at Straw chamber 1;x [mm];y [mm]", 60, -1200, 1200, 60, -1200, 1200));
  	BookHisto("general/hTrackXY_straw2", new TH2F("TrackXY_straw2", "Track (x,y) at Straw chamber 2;x [mm];y [mm]", 60, -1200, 1200, 60, -1200, 1200));
  	BookHisto("general/hTrackXY_straw3", new TH2F("TrackXY_straw3", "Track (x,y) at Straw chamber 3;x [mm];y [mm]", 60, -1200, 1200, 60, -1200, 1200));
  	BookHisto("general/hTrackXY_straw4", new TH2F("TrackXY_straw4", "Track (x,y) at Straw chamber 4;x [mm];y [mm]", 60, -1200, 1200, 60, -1200, 1200));
  	BookHisto("general/DistanceStraw1", new TH1D("DistanceStraw1", "Distance between track pairs in Straw1;[mm]", 500, 0, 1000));
  	BookHisto("general/DistanceLKr", new TH1D("DistanceLKr", "Distance between track pairs in LKr;[mm]", 500, 0, 2000));
  	BookHisto("general/hVertexChi2", new TH1F("VertexChi2", "#chi^{2} of the vertices;#chi^{2}", 100, 0.0, 100.0));
  	BookHisto("general/hZvertex0", new TH1F("Zvertex0", "Zvertex; z [m]", 100, 90, 190));
  	BookHisto("general/hMomentum0", new TH1F("Momentum0", "Kaon momentum;Momentum [GeV/c]", 200, 50, 100));
  	BookHisto("general/hPt0", new TH1F("Pt0", "Kaon transverse momentum wrt beam axis from the DB;Transverse momentum [MeV/c]", 200, 0, 200));
  	BookHisto("general/hM3pi",          new TH1F("M3pi",     "M(3#pi); M(3#pi) [MeV]", 120, 480, 510));
	
  	BookHisto("selected/Cuts",           new TH1F("Cuts",     "Cuts", 10, -0.5, 9.5));
	
	BookHisto("selected/hM3pi",          new TH1F("M3pi",     "M(3#pi); M(3#pi) [MeV]", 120, 480, 510));
	BookHisto("selected/hZvertex",       new TH1F("Zvertex",  "Zvertex; z [m]", 40, 100, 180));
	BookHisto("selected/hMomentum",      new TH1F("Momentum", "Kaon momentum; Momentum [GeV/c]", 60, 72, 78));
	BookHisto("selected/hPt", new TH1F("Pt", "Kaon transverse momentum wrt beam axis from the DB;Momentum [MeV/c]", 200, 0, 200));
	BookHisto("selected/hM3piControlTrigger", new TH1F("M3piControlTrigger", "M(3#pi); M(3#pi) [MeV]", 120, 480, 510));
	BookHisto("selected/hZvertexControlTrigger", new TH1F("ZvertexControlTrigger",  "Zvertex; z [m]", 40, 100, 180));
	BookHisto("selected/hMomentumControlTrigger", new TH1F("MomentumControlTrigger", "Kaon momentum; Momentum [GeV/c]", 60, 72, 78));
	BookHisto("selected/hTrackMomentum", new TH1F("TrackMomentum", "Track momentum;Momentum [GeV/c]", 50, 0, 50));
	
	BookHisto("selected/hK3piEventsPerBurst", new TH1F("K3piEventsPerBurst", "K3pi candidates per burst;Burst ID",cMaxNBursts, -0.5, cMaxNBursts-0.5));
	BookHisto("selected/hK3piEventsPerBurstControlTrigger", new TH1F("K3piEventsPerBurstControlTrigger", "K3pi candidates per burst (control trigger)*DS;Burst ID",cMaxNBursts, -0.5, cMaxNBursts-0.5));
	BookHisto("selected/hK3piEventsPerBurstControlTriggerQM0", new TH1F("K3piEventsPerBurstControlTriggerQM0", "K3pi candidates per burst (control trigger, quality mask = 0)*DS;Burst ID",cMaxNBursts, -0.5, cMaxNBursts-0.5));
	BookHisto("selected/hK3piEventsPerRunControlTrigger", new TH1F("K3piEventsPerRunControlTrigger", "K3pi candidates per run (control trigger)*DS;Run ID",cMaxRunID-cMinRunID+1, cMinRunID-0.5, cMaxRunID+0.5));
	
	// CHOD RICH Cedar studies
	BookHisto("selected/hNCHODHits", new TH1F("NCHODHits", "Number of CHOD hits;Number of hits", 50, -0.5, 49.5));
	BookHisto("selected/hNRICHHits", new TH1F("NRICHHits", "Number of RICH hits;Number of hits", 100, -0.5, 99.5));
	BookHisto("selected/CedarTime", new TH1F("CedarTime", "Cedar Time #minus Mean vertex CHOD track time;Time [ns]", 200, -50, 50));
	BookHisto("selected/NCedarOctants", new TH1F("NCedarOctants", "Number of Cedar Octants in a candidate;N(octants)", 9, -0.5, 8.5));
	
	// LKr-related quantities
	BookHisto("LKr/hNLKrCells", new TH1F("NLKrCells", "Number of LKr cells with signal;Number of cells", 125, -0.5, 249.5));
	BookHisto("LKr/hNLKrClusters", new TH1F("NLKrClusters", "Number of LKr clusters;Number of clusters", 10, -0.5, 9.5));
	BookHisto("LKr/hLKrClusterEnergy", new TH1F("LKrClusterEnergy", "LKr cluster energy;Energy [GeV]", 100, 0, 50));
	BookHisto("LKr/hLKrClusterTime", new TH1F("LKrClusterTime", "LKr cluster time wrt trigger;Time [ns]", 200, -50, 50));
	BookHisto("LKr/hLKrCellTotalEnergy", new TH1F("LKrCellTotalEnergy","LKr total cell (E>40MeV) energy;Total cell energy [GeV]", 70, 0, 70));
	BookHisto("LKr/hLKrCellClusterTotalEnergy", new TH2F("LKrCellClusterTotalEnergy","LKr total cluster energy vs cell energy;Total cell (>40MeV) energy [GeV];Total cluster energy [GeV]",70, 0, 70, 70, 0, 70));
	
	// Straw-GTK missing mass resolution studies
	BookHisto("selected/hMMiss", new TH2F("hMMiss", "Mmiss;From spectrometer [MeV];From recoil [MeV]",100, 270, 370, 100, 270, 370));
	BookHisto("selected/hMMiss_res", new TH2F("hMMiss_res", "Mmiss_res;From spectrometer [MeV];From recoil - from spectrometer [MeV]",50, 270, 370, 150, -15, 15));
	BookHisto("selected/hMMissGTK", new TH2F("hMMissGTK", "MmissGTK;From spectrometer [MeV];From recoil [MeV]",100, 270, 370, 100, 270, 370));
	BookHisto("selected/hMMissGTK_res", new TH2F("hMMissGTK_res", "MmissGTK_res;From spectrometer [MeV];From recoil - from spectrometer [MeV]",50, 270, 370, 150, -15, 15));
	
	// Straw-GTK squared missing mass resolution studies
	BookHisto("selected/hMMiss2", new TH2F("hMMiss2", "Mmiss2;Mmiss2 From Beam and Spectrometer [MeV^{2}];From recoil [MeV^{2}]",100, 72900, 136900, 100, 72900, 136900));
	BookHisto("selected/hMMiss2_res", new TH2F("hMMiss2_res", "Mmiss2_res;Mmiss2 From Beam and Spectrometer [MeV^{2}];From recoil - from spectrometer [MeV^{2}]",100, 72900, 136900, 160, -8000, 8000));
	BookHisto("selected/hMMiss2GTK", new TH2F("hMMiss2GTK", "Mmiss2GTK;From GTK and Spectrometer [MeV^{2}];From recoil [MeV^{2}]",100, 72900, 136900, 100, 72900, 136900));
	BookHisto("selected/hMMiss2GTK_res", new TH2F("hMMiss2GTK_res", "Mmiss2GTK_res;From GTK and Spectrometer [MeV^{2}];From recoil - from spectrometer [MeV^{2}]",100, 72900, 136900, 160, -8000, 8000));
}


void MyK3piSelection::StartOfRunUser() {
	fRunID = GetRunID();
}

void MyK3piSelection::StartOfBurstUser() {
  fBurstID = GetBurstID();
}

void MyK3piSelection::ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType) {
  /// \MemberDescr
  /// variable to be filled, RequestL0SpecialTrigger must be called in the constructor.\n
  /// Process method for special triggers. Called on each special trigger event after each start of burst.
  /// \EndMemberDescr
  /// \param iEvent : Special event number
  /// \param triggerType : Special trigger type (-1 if not known). For this
}

void MyK3piSelection::Process(int /*iEvent*/){
  // Trigger reference time
  const double TdcCalib = 1.0;
  
  Double_t refTime = GetEventHeader()->GetFineTime() * TdcCalib;
  FillHisto("selected/Cuts", 0.0);
 
  if (GetWithMC()) {
    // When running on MC files, we are able to extract the MC truth information from
    // the MC event.
    Event *evt = GetMCEvent();
 
    // The EventBoundary concept comes from the overlay. In such case, multiple MC events
    // are overlayed on each other, forming one big event. The EventBoundary class keeps a
    // record of where each event starts and ends. The boundary of the original (before
    // overlay) event is retrieved here:
    EventBoundary* evt0 = static_cast<EventBoundary*>(evt->GetEventBoundary(0));
 
    if (evt0->GetNKineParts()) { // We check that we do have kineparts in the boundary
      FillHisto("mctrue/hZvertex", 0.001*evt->GetKinePart(0)->GetEndPos().Z()); // [m]
 
      // Checkpoints are particle positions recorded regularly during the tracking
      // in the simulation. This allows to easily retrieve the position of the particle
      // at specific places as it travelled along the detector (not trivial because
      // subjected to magnetic fields).
      // Here we retrieve the checkpoint at the end of GTK and verify that it exists
      // for the beam kaon (meaning it did not decay or interact before)
      Int_t iGTKExit = GetMCInfo()->FindCheckPoint("GTKExit");
      if (evt->GetKinePart(0)->GetNCheckPoints()>0 && iGTKExit>=0 &&
          evt->GetKinePart(0)->GetPosAtCheckPoint(iGTKExit).size()) {
        // If it does exist, we extract the momentum information and fill the relavant
        // histograms.
        TVector3       pos = evt->GetKinePart(0)->GetPosAtCheckPoint(iGTKExit).at(0);
        TLorentzVector mom = evt->GetKinePart(0)->GetMomAtCheckPoint(iGTKExit).at(0);
        Double_t dxdz   = mom.X()/mom.Z();
        Double_t dydz   = mom.Y()/mom.Z();
        Double_t xStart = pos.X();
        Double_t yStart = pos.Y();
        FillHisto("mctrue/hMomentum", 0.001*mom.Vect().Mag()); // [GeV/c]
        // Fill momentum information for daughters of the beam kaon
        for (Int_t i=1; i<evt0->GetNKineParts(); i++) {
          if (evt->GetKinePart(i)->GetParentIndex()==0) { // K+ daughter
            FillHisto("mctrue/hTrackMomentum", 0.001*evt->GetKinePart(i)->GetInitialMomentum().Mag()); // [GeV/c]
            FillHisto("mctrue/hTrackX1",
                      evt->GetKinePart(i)->xAt(GeometricAcceptance::GetInstance()->GetZStraw(0)));
            FillHisto("mctrue/hTrackY1",
                      evt->GetKinePart(i)->yAt(GeometricAcceptance::GetInstance()->GetZStraw(0)));
          }
        }
        FillHisto("mctrue/hdxdz", dxdz);
        FillHisto("mctrue/hdydz", dydz);
        FillHisto("mctrue/hXYstart", xStart, yStart);
        FillHisto("mctrue/hXstart", xStart);
        FillHisto("mctrue/hYstart", yStart);
        FillHisto("mctrue/hdxdzVsXstart", xStart, dxdz);
        FillHisto("mctrue/hdydzVsXstart", xStart, dydz);
        FillHisto("mctrue/hdxdzVsYstart", yStart, dxdz);
        FillHisto("mctrue/hdydzVsYstart", yStart, dydz);
        FillHisto("mctrue/pdxdzVsXstart", xStart, dxdz);
        FillHisto("mctrue/pdydzVsXstart", xStart, dydz);
        FillHisto("mctrue/pdxdzVsYstart", yStart, dxdz);
        FillHisto("mctrue/pdydzVsYstart", yStart, dydz);
      }
    }
 
    ////////////////////////////////////////////////////////////////////////////////////
    // Angular dependence of the true 3pi mass: this is purely a blue field effect.
    // Acceptance (mainly chamber hole) is not taken into account, therefore
    // the observed phi-variation of the mass is smaller than in the reconstructed data.
 
    // Here we are doing some K3pi specific operation, so we make sure that this is what
    // we have in the input (decay type 10 is K3pi).
    if (GetMCInfo()->GetDecayType()==10 && evt0->GetNKineParts()==4) {
      // Get the checkpoint and make sure it exists for all three pions
      Int_t iSpectEntry = GetMCInfo()->FindCheckPoint("STRAWCH0Entry");
      if (iSpectEntry>=0) {
        if (evt->GetKinePart(1)->GetMomAtCheckPoint(iSpectEntry).size() &&
            evt->GetKinePart(2)->GetMomAtCheckPoint(iSpectEntry).size() &&
            evt->GetKinePart(3)->GetMomAtCheckPoint(iSpectEntry).size()) {
            // Compute the invariant mass of the 3 pions
            TLorentzVector p3pi = TLorentzVector(0.0, 0.0, 0.0, 0.0);
            for (Int_t i=1; i<=3; i++) {
              p3pi += evt->GetKinePart(i)->GetMomAtCheckPoint(iSpectEntry).at(0);
            }
            Double_t phi = evt->GetKinePart(1)->GetInitialMomentum().Phi() / TMath::Pi(); // odd pion
            FillHisto("mctrue/hM3piVsPhi", phi, p3pi.M());
        }
      }
    }
  }
 
  // The TriggerConditions class is the interface with the trigger conditions database.
  // We want to know whether this event comes from a Physics trigger, or from a control
  // trigger (implications on timing, downscaling, bias)
  Bool_t physicsTrigger = TriggerConditions::GetInstance()->IsPhysicsTrigger(GetL0Data());
  Bool_t controlTrigger = TriggerConditions::GetInstance()->IsControlTrigger(GetL0Data());
 
  // This is kind of the equivalent of what we have done in the previous step of the tutorial.
  // Instead of printing the information, we plot it.
  TRecoCedarEvent* cedarEvent = GetEvent<TRecoCedarEvent>();
  if (GetWithMC()) {
    for (Int_t i=0; i<cedarEvent->GetNCandidates(); i++) {
      TRecoCedarCandidate* cand = static_cast<TRecoCedarCandidate*>(cedarEvent->GetCandidate(i));
      // We align the plot such that the trigger is at 0.
      FillHisto("mctrue/CedarTimeMC", cand->GetTime() - refTime);
    }
  }
 

  /////////////////////////
  // Run the vertexing tool
  // The vertexing tool will take the downstream tracks (the 3 pions) and reconstruct their
  // vertex. This is not trivial as we might have spurious tracks, magnetic fields are involved
  // and it also takes into account timing.
  // This tool is actually another analyzer that will run previous to this one on every event.
  // It is providing some output variables that we retrieve with the GetOutput method (the return
  // type must be checked in the documentation of the tool itself).
  std::vector<SpectrometerTrackVertex> vertices = *(std::vector<SpectrometerTrackVertex>*)GetOutput("SpectrometerVertexBuilder.Output3");
  // This other tool on the other hand, builds coherent downstream track by assembling the information
  // provided by the spectrometer tracks and those provided by the other downstream detectors, and
  // supposedly belonging to the same track.
  std::vector<DownstreamTrack> tracks = *(std::vector<DownstreamTrack>*)GetOutput("DownstreamTrackBuilder.Output");
 
  FillHisto("general/hNTracks", tracks.size());
  FillHisto("general/hNVertices", vertices.size());
  // Fill the track positions at the 4 Straw stations (extrapolated)
  for (UInt_t i=0; i<tracks.size(); i++) {
    FillHisto("general/hTrackXY_straw1", tracks[i].xAt(183508.), tracks[i].yAt(183508.));
    FillHisto("general/hTrackXY_straw2", tracks[i].xAt(194066.), tracks[i].yAt(194066.));
    FillHisto("general/hTrackXY_straw3", tracks[i].xAt(204459.), tracks[i].yAt(204459.));
    FillHisto("general/hTrackXY_straw4", tracks[i].xAt(218885.), tracks[i].yAt(218885.));
  }
 
  // Compute the number of good vertices
  // For now a good vertex is simply a vertex composed of exactly 3 tracks.
  Int_t nGoodVertex = 0;
  Int_t vtx_index = -1;
  for (UInt_t i=0; i<vertices.size(); i++) {
    if (vertices[i].GetNTracks()!=3) continue;

    Double_t vtxTimeDelta = vertices[i].GetCHODTime()-refTime;
    if (controlTrigger) FillHisto("general/hVertexTimeControl", vtxTimeDelta);
    if (physicsTrigger) FillHisto("general/hVertexTimePhysics", vtxTimeDelta);
	if (fabs(vtxTimeDelta) < 5.0) {
    	vtx_index = i;
    	nGoodVertex++;
	}

	// // We want the total charge to be +1 as the beam kaons are K+
	// Int_t totalCharge = vertices[vtx_index].GetCharge();
	// if (totalCharge!=1) return;
	// FillHisto("selected/Cuts", 3.0);
	
	// // We are also interest in the distribution of the vertex quality and z position
	// Double_t chi2 = vertices[vtx_index].GetChi2();
	// FillHisto("general/hVertexChi2", chi2);
	// if (chi2>25.0) return;
	// FillHisto("selected/Cuts", 4.0);
	// Double_t zVertex = vertices[vtx_index].GetPosition().z();
	// FillHisto("general/hZvertex0", 0.001*zVertex); // [m]
	// // This is our definition of the fiducial decay region
	// if (zVertex<104000 || zVertex>180000) return;
	// FillHisto("selected/Cuts", 5.0);
  }


 
  // If the event does not contain exactly one vertex with 3 tracks, we skip it.
  if (nGoodVertex!=1) return;
  // If this is a good event, we record that it passed cut 1.0
  FillHisto("selected/Cuts", 1.0);
 
  // We verify that all the tracks are inside the geometric acceptance of the spectrometer and NewCHOD
  for (Int_t i=0; i<3; i++) {
    TRecoSpectrometerCandidate* Scand = vertices[vtx_index].GetSpectrometerCandidate(i);
    if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 0)) return;
    if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 1)) return;
    if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 2)) return;
    if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kSpectrometer, 3)) return;
    if (!GeometricAcceptance::GetInstance()->InAcceptance(Scand, NA62::kNewCHOD))         return;
  }
  // All tracks succeeded the check, we record that it passed our cut 2.0
  FillHisto("selected/Cuts", 2.0);
 
  // We are now interested in the track separations in STRAW1 and LKr planes
  // We might want to cut on this in the future.
  for (Int_t i=0; i<3; i++) {
    TRecoSpectrometerCandidate* cand = vertices[vtx_index].GetSpectrometerCandidate(i);
    double x1s = cand->xAt(GeometricAcceptance::GetInstance()->GetZStraw(0)); // Straw1
    double y1s = cand->yAt(GeometricAcceptance::GetInstance()->GetZStraw(0));
    double x1c = cand->xAt(241093.0); // LKr
    double y1c = cand->yAt(241093.0);
    for (Int_t j=i+1; j<3; j++) {
      TRecoSpectrometerCandidate* cand2 = vertices[vtx_index].GetSpectrometerCandidate(j);
      double x2s = cand2->xAt(GeometricAcceptance::GetInstance()->GetZStraw(0)); // Straw1
      double y2s = cand2->yAt(GeometricAcceptance::GetInstance()->GetZStraw(0));
      double x2c = cand2->xAt(241093.0); // LKr
      double y2c = cand2->yAt(241093.0);
      double rs  = sqrt((x1s-x2s)*(x1s-x2s)+(y1s-y2s)*(y1s-y2s)); // Straw1
      double rc  = sqrt((x1c-x2c)*(x1c-x2c)+(y1c-y2c)*(y1c-y2c)); // LKr
      // We fill the histogram with those distances
      FillHisto("general/DistanceStraw1", rs);
      FillHisto("general/DistanceLKr", rc);
    }
  }
 
  // We are also interest in the distribution of the vertex quality and z position
  Double_t chi2 = vertices[vtx_index].GetChi2();
  FillHisto("general/hVertexChi2", chi2);
  Double_t zVertex = vertices[vtx_index].GetPosition().z();
  FillHisto("general/hZvertex0", 0.001*zVertex); // [m]
 
  // Finally we will here check the reconstructed total momentum, assuming 3 pions
  // We might want to cut on the total and transverse momentum as well as the
  // reconstructed invariant mass.
  TLorentzVector v[3];
  for (Int_t i=0; i<3; i++) {
    v[i].SetVectM(vertices[vtx_index].GetTrackThreeMomentum(i), MPI);
  }
 
  TLorentzVector k3pi_FourMomentum  = v[0] + v[1] + v[2];
  TVector3       k3pi_Momentum      = k3pi_FourMomentum.Vect();
  Double_t       m3pi               = k3pi_FourMomentum.M();
   
  // Longitudinal and transverse momenta with respect to the beam axis
  TVector3 beamAxis = BeamParameters::GetInstance()->GetBeamThreeMomentum();
  Double_t pl       = (beamAxis*k3pi_Momentum) / beamAxis.Mag();
  Double_t pt       = sqrt(k3pi_Momentum*k3pi_Momentum-pl*pl);
   
  FillHisto("general/hMomentum0", 0.001*k3pi_Momentum.Mag()); // [GeV/c]
  if (fabs(k3pi_Momentum.Mag()-75000.0)>3000.0) return;
  FillHisto("selected/Cuts", 6.0);
   
  FillHisto("general/hPt0", pt); // [MeV/c]
  if (pt>30.0) return;
  FillHisto("selected/Cuts", 7.0);

  FillHisto("general/hM3pi",     m3pi);
  if (m3pi<490.0 || m3pi>497.0) return; // last cut
  FillHisto("selected/Cuts", 8.0);


  // K3pi yields per trigger: control and multi-track triggers.
// Runs with variable or unknown control trigger downscaling cannot be processed: request DS>0.
FillHisto("selected/hK3piEventsPerBurst", fBurstID);
Int_t controlDownscaling = TriggerConditions::GetInstance()->GetControlTriggerDownscaling(fRunID);
if (controlTrigger && controlDownscaling>0) {
  FillHisto("selected/hK3piEventsPerBurstControlTrigger", fBurstID, controlDownscaling);
  FillHisto("selected/hK3piEventsPerRunControlTrigger", fRunID, controlDownscaling);
  if (!GetEventHeader()->GetEventQualityMask()) {
    FillHisto("selected/hK3piEventsPerBurstControlTriggerQM0", fBurstID, controlDownscaling);
  }
}
 
// Check that the event passes L0 and L1 triggers (We use L0 multi-track trigger which is defined by
// RICH-QX
Int_t multiTrackTriggerID = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-QX");
Int_t multiTrackDownscaling = TriggerConditions::GetInstance()->GetControlTriggerDownscaling(fRunID);
if (TriggerConditions::GetInstance()->L0TriggerOn(fRunID, GetL0Data(), multiTrackTriggerID) && multiTrackDownscaling > 0
  && TriggerConditions::GetInstance()->L1TriggerOnIgnoringFlagging(fRunID, GetL1Data(), multiTrackTriggerID)) {
  FillHisto("selected/hK3piEventsPerBurstMultiTrackTrigger", fBurstID, 1.0 * multiTrackDownscaling);
  FillHisto("selected/hK3piEventsPerRunMultiTrackTrigger", fRunID, 1.0 * multiTrackDownscaling);
}
 
// Same as before, but we have only selected events
FillHisto("selected/hZvertex", 0.001*zVertex); // [m]
FillHisto("selected/hMomentum", 0.001*k3pi_Momentum.Mag()); // [GeV/c]
FillHisto("selected/hPt", pt); // [MeV/c]
// And here we separate further events from control trigger
if (controlTrigger) {
  FillHisto("selected/hZvertexControlTrigger", 0.001*zVertex); // [m]
  FillHisto("selected/hMomentumControlTrigger", 0.001*k3pi_Momentum.Mag()); // [GeV/c]
}
// We look at the individual pions momenta
for (Int_t i=0; i<3; i++) {
  FillHisto("selected/hTrackMomentum", 0.001*vertices[vtx_index].GetTrackThreeMomentum(i).Mag());
}
 
//Study detector responses: number of hits per selected event
// CHOD response studies
TRecoCHODEvent* chodEvent = GetEvent<TRecoCHODEvent>();
FillHisto("selected/hNCHODHits", chodEvent->GetNHits());
 
// RICH response studies
TRecoRICHEvent* richEvent = GetEvent<TRecoRICHEvent>();
Int_t nRICHHitsAll = richEvent->GetNHits(); // including super-cells
Int_t nRICHHits = 0;
for (Int_t i=0; i<nRICHHitsAll; i++) {
  TRecoRICHHit* hit = static_cast<TRecoRICHHit*>(richEvent->GetHit(i));
  if (hit->GetOrSuperCellID()==0) nRICHHits++; // no super-cells
}
FillHisto("selected/hNRICHHits", nRICHHits);
 
// This is the weighted mean CHOD+NewCHOD time of the vertex tracks,
// with up to one inconsistent time measurement excluded.
Double_t vtxTime = vertices[vtx_index].GetCHODTime();
 
// Cedar/KTAG response studies: timing and number of sectors in coincidence
Double_t cedarTime = -999.0, dT_CHOD_Cedar = 999.0;
for (Int_t i=0; i<cedarEvent->GetNCandidates(); i++) {
  TRecoCedarCandidate* cand = static_cast<TRecoCedarCandidate*>(cedarEvent->GetCandidate(i));
  Double_t dT = cand->GetTime() - vtxTime;
  FillHisto("selected/CedarTime", dT);
  // Check number of octants only for candidates close enough to the vertex time
  if (fabs(dT)<2.0) FillHisto("selected/NCedarOctants", cand->GetNSectors());
  // Select the candidate with more than 5 sectors closest in time with the vertex
  if (cand->GetNSectors()>=5 && fabs(dT)<fabs(dT_CHOD_Cedar)) {
    cedarTime = cand->GetTime();
    dT_CHOD_Cedar = dT;
  }
}
 
// LKr response studies
TRecoLKrEvent* lkrEvent = GetEvent<TRecoLKrEvent>();
Double_t totalCellEnergy = 0.0; // [MeV]
Double_t totalInTimeClusterEnergy = 0.0; // [MeV]
 
Int_t nCells = lkrEvent->GetNHits();
Int_t nClusters = lkrEvent->GetNCandidates();
FillHisto("LKr/hNLKrCells", nCells);
FillHisto("LKr/hNLKrClusters", nClusters);
 
// Compute the total LKr deposited energy (in cells with at least 40MeV - zero suppression)
for (Int_t i=0; i<nCells; i++) {
  TRecoLKrHit *hit = static_cast<TRecoLKrHit*>(lkrEvent->GetHit(i));
  Double_t energy = hit->GetEnergy();
  if (energy>40.0) totalCellEnergy += energy;
}
FillHisto("LKr/hLKrCellTotalEnergy", 0.001*totalCellEnergy);
FillHisto("LKr/hLKrCellClusterTotalEnergy", 0.001*totalCellEnergy, 0.001*lkrEvent->GetEnergyTotal());
 
// Now check the energy and time of each individual reconstructed cluster
for (Int_t i=0; i<nClusters; i++) {
  TRecoLKrCandidate* cand = static_cast<TRecoLKrCandidate*>(lkrEvent->GetCandidate(i));
  FillHisto("LKr/hLKrClusterEnergy", 0.001*cand->GetEnergy());
  FillHisto("LKr/hLKrClusterTime", cand->GetTime()-refTime);
  // Count the total energy within 6ns of the reference time
  if (fabs(cand->GetTime()-refTime)<6.0)
    totalInTimeClusterEnergy += cand->GetEnergy();
}
 
// Squared missing mass and missing mass resolution studies
Double_t mpipi = (v[0] + v[1]).M(); // two random pions
TLorentzVector kaon;
kaon.SetVectM(BeamParameters::GetInstance()->GetBeamThreeMomentum(), MKCH);
Double_t mMiss = (kaon - v[2]).M(); // the third pion
FillHisto("selected/hMMiss", mpipi, mMiss);
FillHisto("selected/hMMiss_res", mpipi, mMiss - mpipi);
 
Double_t mMiss2 = (kaon - v[2]).M2(); // Mmiss2 = (PK - Ppi3)^2
Double_t m2pipi = (v[0] + v[1]).M2(); // M12^2 = (Ppi1 + Ppi2)^2
FillHisto("selected/hMMiss2", mMiss2, m2pipi);
FillHisto("selected/hMMiss2_res", mMiss2, mMiss2 - m2pipi);
 
// Match the third STRAW pion track to a GTK track
if (fabs(dT_CHOD_Cedar) < 5.0) { // matching Cedar candidate found
  TRecoGigaTrackerEvent *gtkEvent = GetEvent<TRecoGigaTrackerEvent>();
  fSGMatching->Match(gtkEvent, tracks[vertices[vtx_index].GetTrackIndex(2)].GetSpectrometerCandidate(), cedarTime, NA62::kCedar);
  if (fSGMatching->BestGTKTrackFound() && fSGMatching->GetBestDiscriminant() > 0.01 && fSGMatching->GetBestCDA() < 7.0) {
    TLorentzVector kaonGTK; // beam kaon
    kaonGTK.SetVectM(fSGMatching->GetBestCorrectedBeamMomentum(), MKCH);
    TLorentzVector pionGTK; // downstream pion
    pionGTK.SetVectM(fSGMatching->GetBestCorrectedTrackMomentum(), MPI);
    Double_t mMissGTK = (kaonGTK - pionGTK).M();
    FillHisto("selected/hMMissGTK", mpipi, mMissGTK);
    FillHisto("selected/hMMissGTK_res", mpipi, mMissGTK - mpipi);
    Double_t mMiss2GTK = (kaonGTK - pionGTK).M2();
    FillHisto("selected/hMMiss2GTK", mMiss2GTK, m2pipi);
    FillHisto("selected/hMMiss2GTK_res", mMiss2GTK, mMiss2GTK - m2pipi);
  }
}
}

void MyK3piSelection::PostProcess() {
  /// \MemberDescr
  /// This function is called after an event has been processed by all analyzers. It could be used to free some memory allocated
  /// during the Process.
  /// \EndMemberDescr
}

void MyK3piSelection::EndOfBurstUser() {
  /// \MemberDescr
  /// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
  /// Do here your start/end of burst processing if any.
  /// Be careful: this is called after the event/file has changed.
  /// \EndMemberDescr
}

void MyK3piSelection::EndOfRunUser() {
  /// \MemberDescr
  /// This method is called at the end of the processing (corresponding to a end of run in the normal NA62 data taking)\n
  /// Do here your end of run processing if any
  /// \EndMemberDescr
}

void MyK3piSelection::EndOfJobUser(){
  SaveAllPlots();
}

void MyK3piSelection::DrawPlot() {
  /// \MemberDescr
  /// This method is called at the end of processing to draw plots when the -g option is used.\n
  /// If you want to draw all the plots, just call
  /// \code
  ///   DrawAllPlots();
  /// \endcode
  /// Or get the pointer to the histogram with
  /// \code
  ///   fHisto.GetTH1("histoName");// for TH1
  ///   fHisto.GetTH2("histoName");// for TH2
  ///   fHisto.GetGraph("graphName");// for TGraph and TGraphAsymmErrors
  ///   fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
  /// \endcode
  /// and manipulate it as usual (TCanvas, Draw, ...)
  /// \EndMemberDescr
}

MyK3piSelection::~MyK3piSelection() {
  	if (fSGMatching) delete fSGMatching;
}
