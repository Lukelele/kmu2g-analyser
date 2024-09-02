#define KAON_REST_MASS 0.49368
#define MUON_REST_MASS 0.10566

#include "Kmu2gSelection.hh"

#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixDfwd.h"
#include "TMatrixDUtils.h"

#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "Event.hh"
#include "KaonDecayConstants.hh"
#include "LAVMatching.hh"
#include "Persistency.hh"
#include "SAVMatching.hh"
#include "TriggerConditions.hh"

#include "Event.hh"
#include "Persistency.hh"
#include "functions.hh"

#include <iostream>
#include <stdlib.h>


using namespace NA62Analysis;
using namespace TMath;


Kmu2gSelection::Kmu2gSelection(Core::BaseAnalysis *ba): Analyzer(ba, "Kmu2gSelection") {
	RequestTree("Cedar", new TRecoCedarEvent, "Reco");
	RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
	RequestTree("CHOD", new TRecoCHODEvent, "Reco");
	RequestTree("RICH", new TRecoRICHEvent, "Reco");
	RequestTree("LKr", new TRecoLKrEvent, "Reco");
	RequestTree("LAV", new TRecoLAVEvent, "Reco");
	RequestTree("IRC", new TRecoIRCEvent, "Reco");
	RequestTree("SAC", new TRecoSACEvent, "Reco");
	RequestL0Data();

	LKrZ = GeometricAcceptance::GetInstance()->GetZLKr();

	modeCanvas = new TCanvas("c1", "Decay Mode Scatter", 200, 10, 700, 500);
	modeScatter = new TGraph();

	NTracksBeforeCut = 0;
	NTracksAfterCut = 0;
	
	m_eventCount = 0;

    filenumber = 1;
    
	initCSV("/eos/user/y/yuanye/Analysis/data/run3");
}

void Kmu2gSelection::InitOutput() {

}

void Kmu2gSelection::InitHist() {
	BookHisto("GTKMomentum", new TH1F("GTKMomentum", "GTK Momentum;Momentum [GeV/c]", 300, 70, 80));
	BookHisto("GTKEnergy", new TH1F("GTKEnergy", "GTK Energy;Energy [GeV]", 100, 70, 80));
	BookHisto("GTKBeta", new TH1F("GTKBeta", "GTK Beta;Beta", 1000, 0.99996, 1));
	BookHisto("GTKGamma", new TH1F("GTKGamma", "GTK Gamma;Gamma", 300, 130, 180));
	BookHisto("NLKrCells", new TH1F("hNLKrCells", "Number of LKr cells", 100, -10, 100));
	BookHisto("NLKrClusters", new TH1F("hNLKrClusters", "Number of LKr clusters", 100, -10, 15));
    BookHisto("LKrClusterEnergy", new TH1F("LKrClusterEnergy", "LKr Cluster Energy;Energy [GeV]", 100, -5, 60));
	BookHisto("NVertices", new TH1F("hNVertices", "Number of vertices", 100, -1, 3));
	BookHisto("NTracks", new TH1F("hNTracks", "Number of tracks", 11, -0.5, 10.5));
	BookHisto("CDA", new TH1F("CDA", "CDA", 100, -0.5, 50));
	BookHisto("TrackTime", new TH1F("TrackTimes", "Time", 100, -3, 3));
	BookHisto("NGoodTracks", new TH1F("hNGoodTracks", "Number of good tracks", 11, -0.5, 10.5));
	BookHisto("QChi2Track", new TH1F("QChi2Track", "QChi2Track", 100, -2, 8));
	BookHisto("EoP", new TH1F("EoP", "Track E/p; E/p", 150, -1.0, 1.5));
    BookHisto("trackLKrTime", new TH1F("trackLKrTime", "Track LKr Time;Time [ns]", 100, -10, 10));
	BookHisto("MuonMomentum", new TH1F("MuonMomentum", "Muon Momentum;Momentum [GeV/c]", 300, 0, 100));
    BookHisto("MuonTrueEnergy", new TH1F("MuonTrueEnergy", "Muon Energy;Energy [GeV]", 100, -5, 100));
	BookHisto("MuonLKrClusterEnergy", new TH1F("MuonLKrClusterEnergy", "Muon LKr Cluster Energy;Energy [GeV]", 100, -5, 20));
	BookHisto("MuonMomentumPrime", new TH1F("MuonMomentumPrime", "Muon Momentum Prime;Momentum [GeV/c]", 300, 0, 0.5));
	BookHisto("MuonEnergyPrime", new TH1F("MuonEnergyPrime", "Muon Energy Prime;Energy [GeV]", 100, 0, 0.5));
	BookHisto("PhotonEnergy", new TH1F("PhotonEnergy", "Photon Energy / Momentum;Energy [GeV]", 100, -10, 75));
    BookHisto("PhotonEnergyPrime", new TH1F("PhotonEnergyPrime", "Photon Energy / Momentum Prime;Energy [GeV]", 100, 0, 0.5));
	BookHisto("GammaMuAngle", new TH1F("GammaMuAngle", "Gamma Mu Angle;Angle [rad]", 100, -0.05, 3.5));
	BookHisto("CosGammaMuAngle", new TH1F("CosGammaMuAngle", "Cos(Gamma Mu Angle);Cos(Gamma Mu Angle)", 100, -1, 1));
	BookHisto("MissingEnergy", new TH1F("MissingEnergy", "Missing Energy;Energy [GeV]", 100, -0.05, 0.5));
	BookHisto("MissingMomentum", new TH1F("MissingMomentum", "Missing Momentum;Momentum [GeV/c]", 100, -0.05, 0.5));
	BookHisto("MissingMassSquared", new TH1F("MissingMassSquared", "Missing Mass Squared;Mass [GeV/c^{2}]", 300, -0.05, 0.15));
	BookHisto("KaonMass", new TH1F("KaonMass", "Kaon Mass;Mass [GeV/c^{2}]", 300, 0.2, 0.7));
	BookHisto("CDAAfterCut", new TH1F("CDAAfterCut", "CDA After Cut;CDA", 100, -0.5, 50));
	BookHisto("XGamma", new TH1F("XGamma", "XGamma", 100, -0.5, 1.5));
	BookHisto("XMu", new TH1F("XMu", "XMu", 100, -0.5, 1.5));
}

void Kmu2gSelection::StartOfRunUser() {

}

void Kmu2gSelection::StartOfBurstUser() {

}

void Kmu2gSelection::ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType) {

}


void Kmu2gSelection::Process(int iEvent) {
	m_eventCount++;

	// Get the GTK event
	TRecoGigaTrackerEvent *GTKevent = GetEvent<TRecoGigaTrackerEvent>();
	Int_t NGTKCandidates = GTKevent->GetNCandidates();

	std::vector<TRecoGigaTrackerCandidate*> GTKCandidates;

	for (int i = 0; i < NGTKCandidates; i++) {
		TRecoGigaTrackerCandidate* candidate = static_cast<TRecoGigaTrackerCandidate*>(GTKevent->GetCandidate(i));
		GTKCandidates.push_back(candidate);
	}

	// filter out any events which did not come from a kaon decay
	if (NGTKCandidates == 0) return;

	// Get kaon momentum
	TVector3 GTKMomentum = GTKCandidates[0]->GetMomentum();
	GTKMomentum = TVector3(GTKMomentum.x() / 1000, GTKMomentum.y() / 1000, GTKMomentum.z() / 1000);
	FillHisto("GTKMomentum", GTKMomentum.Mag());

	// Get kaon energy
	Double_t GTKEnergy = sqrt(GTKMomentum.Mag2() + KAON_REST_MASS * KAON_REST_MASS);
	FillHisto("GTKEnergy", GTKEnergy);

	Double_t GTKBeta = GTKMomentum.Mag() / GTKEnergy;
	Double_t GTKGamma = 1 / sqrt(1 - GTKBeta * GTKBeta);

	FillHisto("GTKBeta", GTKBeta);
	FillHisto("GTKGamma", GTKGamma);

	// rotate the axis such that the kaon momentum is along the z-axis
	// Rodrigues' formula
	TVector3 r = GTKMomentum.Cross(TVector3(0, 0, 1)).Unit();
	Double_t theta = GTKMomentum.Angle(TVector3(0, 0, 1));

	// skew-symmetric matrix
	TMatrixD K(3, 3);
	K[0][0] = 0;
	K[1][1] = 0;
	K[2][2] = 0;
	K[0][1] = -r.z();
	K[0][2] = r.y();
	K[1][0] = r.z();
	K[1][2] = -r.x();
	K[2][0] = -r.y();
	K[2][1] = r.x();

	TMatrixD I(3, 3);
	I.UnitMatrix();

	// Transform matrix
	TMatrixD R(3, 3);
	R = I + (Sin(theta) * K + (1 - Cos(theta)) * K * K);

	// apply the transformation
	TVector3 standardGTKMomentum = R * GTKMomentum;


	// Get the LKr events for photon detection
	TRecoLKrEvent *LKRevent = GetEvent<TRecoLKrEvent>();
	Int_t NCells = LKRevent->GetNHits();
    Int_t NClusters = LKRevent->GetNCandidates();
    FillHisto("NLKrCells", NCells);
    FillHisto("NLKrClusters", NClusters);

    // Make sure there are 2 clusters in the LKr corresponding to a muon and a photon
	if (NClusters != 2) return;

	std::vector<TRecoLKrCandidate*> LKrClusters;

	for (int i = 0; i < 2; i++) {
		TRecoLKrCandidate *cluster = static_cast<TRecoLKrCandidate*>(LKRevent->GetCandidate(i));
		LKrClusters.push_back(cluster);
		Double_t clusterEnergy = cluster->GetEnergy() / 1000;    // convert to GeV
		FillHisto("LKrClusterEnergy", clusterEnergy);
	}

	// Retrieve the spectrometer vertices
	std::vector<SpectrometerTrackVertex> vertices = *(std::vector<SpectrometerTrackVertex>*)GetOutput("SpectrometerVertexBuilder.Output3");
	FillHisto("NVertices", vertices.size());

	// Retrieve the downstream tracks
	std::vector<DownstreamTrack> tracks = *GetOutput<std::vector<DownstreamTrack>>("DownstreamTrackBuilder.Output");
	FillHisto("NTracks", tracks.size());
	NTracksBeforeCut++;

	// Get the quality tracks (1 track vertex)
	if (tracks.size() != 1) return;
	DownstreamTrack track = tracks[0];

	Double_t cda = track.GetBeamAxisCDA();
	FillHisto("CDA", cda);
	if (cda > 30) return;

	// Get track stats
	Int_t trackCharge = track.GetCharge();
	Double_t trackChi2 = track.GetChi2();
	Int_t trackNChambers = track.GetNChambers();
	Double_t eop = track.GetLKrEoP();

	// Get the track time
	Double_t trackTime = (track.CHODTimeExists()) ? track.GetCHODTime() : track.GetNewCHODTime();
	Double_t refTime = GetEventHeader()->GetFineTime() * TriggerCalib;
	trackTime -= refTime;
	FillHisto("TrackTime", trackTime);

	// make sure the muon charge is positive 1
	if(trackCharge != 1) return;

	// select for tracks within 10 ns of trigger time
	if(fabs(trackTime) > 10.0) return;

	// check if track passes through all the chambers
    if(trackNChambers != 4) return;

	// get the chi squared
	FillHisto("QChi2Track", trackCharge * trackChi2);
    if(trackChi2 > 4.0) return;

	// plot the number of good quality tracks
	FillHisto("NGoodTracks", 1);

	// get the e/p for the muon
	FillHisto("EoP", track.GetLKrEoP());
    // since the muons should have a very low E/p, we can cut on this
    if(eop > 0.15) return;

	// get the muon momentum after selection
	Double_t muonMomentumMag = track.GetMomentum() / 1000;     // convert to GeV
	TVector3 muonMomentum = TVector3(track.GetMomentumBeforeMagnet().x() / 1000, track.GetMomentumBeforeMagnet().y() / 1000, track.GetMomentumBeforeMagnet().z() / 1000);
	Double_t muonEnergy = sqrt(muonMomentumMag * muonMomentumMag + MUON_REST_MASS * MUON_REST_MASS);

	// if there are 0 or more than 1 clusters, remove the event
	if (track.GetLKrNAssociatedClusters() != 1) return;

	TRecoLKrCandidate* trackCandidate = track.GetLKrCandidate(0);
	Double_t muonClusterEnergy = trackCandidate->GetEnergy() / 1000;    // convert to GeV
	FillHisto("MuonLKrClusterEnergy", muonClusterEnergy);


    FillHisto("trackLKrTime", trackCandidate->GetTime() - refTime);
    // if time if greater than 5 ns, remove the event
    if (fabs(trackCandidate->GetTime() - refTime) > 5.0) return;

	// Find the photon by finding the LKr cluster that is not associated with the track
	for (int i = 0; i < LKrClusters.size(); i++) {
		if (LKrClusters[i] != trackCandidate) {
			Double_t photonEnergy = LKrClusters[i]->GetEnergy() / 1000;
			FillHisto("PhotonEnergy", photonEnergy);

			// change the muon momentum and energy from lab frame to the kaon rest frame
			TVector3 standardMuonMomentum = R * muonMomentum;
			TVector3 standardMuonMomentumPrime = TVector3(standardMuonMomentum.x(), standardMuonMomentum.y(), GTKGamma * (standardMuonMomentum.z() - GTKBeta * muonEnergy));
			Double_t muonEnergyPrime = GTKGamma * (muonEnergy - GTKBeta * standardMuonMomentum.z());

            FillHisto("MuonTrueEnergy", muonEnergy);
			FillHisto("MuonMomentum", standardMuonMomentum.Mag());
			FillHisto("MuonMomentumPrime", standardMuonMomentumPrime.Mag());
			FillHisto("MuonEnergyPrime", muonEnergyPrime);


			// change the photon momentum and energy from lab frame to the kaon rest frame
			TVector3 vertexPosition = track.GetBeamAxisVertex();
            TVector3 photonMomentum = (TVector3(LKrClusters[i]->GetX(), LKrClusters[i]->GetY(), LKrZ) - vertexPosition).Unit() * photonEnergy;
            TVector3 standardPhotonMomentum = R * photonMomentum;
            TVector3 standardPhotonMomentumPrime = TVector3(standardPhotonMomentum.x(), standardPhotonMomentum.y(), GTKGamma * (standardPhotonMomentum.z() - GTKBeta * photonEnergy));
            Double_t photonEnergyPrime = GTKGamma * (photonEnergy - GTKBeta * standardPhotonMomentum.z());

			FillHisto("PhotonEnergyPrime", photonEnergyPrime);

			// Get the angle between the photon and the kaon momentum
			Double_t angle = standardPhotonMomentumPrime.Angle(standardMuonMomentumPrime);
			FillHisto("GammaMuAngle", angle);

			Double_t cosAngle = standardPhotonMomentumPrime.Dot(standardMuonMomentumPrime) / (standardPhotonMomentumPrime.Mag() * standardMuonMomentumPrime.Mag());
			FillHisto("CosGammaMuAngle", cosAngle);

            // calculate the missing energy and momentum
			Double_t missingEnergy = KAON_REST_MASS - muonEnergyPrime - photonEnergyPrime;
			TVector3 missingMomentum = TVector3(0, 0, 0) - standardMuonMomentumPrime - standardPhotonMomentumPrime;
			Double_t missingMassSquared = missingEnergy * missingEnergy - missingMomentum.Mag2();
			FillHisto("MissingEnergy", missingEnergy);
			FillHisto("MissingMomentum", missingMomentum.Mag());
			FillHisto("MissingMassSquared", missingMassSquared);

			// Calculate the Kaon mass assuming the neutrino mass is zero
			TVector3 neutrinoMomentum = -standardMuonMomentumPrime - standardPhotonMomentumPrime;
			Double_t kaonMass = neutrinoMomentum.Mag() + photonEnergyPrime + muonEnergyPrime;
			FillHisto("KaonMass", kaonMass);

            // Cut on the missing mass (kaon mass)
            if (kaonMass < 0.48 || kaonMass > 0.51) return;

			Double_t cdaAfterCut = track.GetBeamAxisCDA();
			FillHisto("CDAAfterCut", cdaAfterCut);


			// calculate the x_gamma and x_mu, plot on histograms and scatter plot
			Double_t x_gamma = 2 * photonEnergyPrime / KAON_REST_MASS;
			Double_t x_mu = 2 * muonEnergyPrime / KAON_REST_MASS;

			FillHisto("XGamma", x_gamma);
			FillHisto("XMu", x_mu);


            // Final Processing
			modeScatter->AddPoint(x_mu, x_gamma);

			saveCSV(m_eventCount, GTKMomentum.x(), GTKMomentum.y(), GTKMomentum.z(), GTKEnergy, GTKGamma, GTKBeta, LKrClusters.size(), LKrClusters.size(), muonClusterEnergy, vertices.size(), tracks.size(), cdaAfterCut, trackTime, trackChi2, eop, muonMomentum.x(), muonMomentum.y(), muonMomentum.z(), muonEnergy, muonClusterEnergy, standardMuonMomentumPrime.x(), standardMuonMomentumPrime.y(), standardMuonMomentumPrime.z(), muonEnergyPrime, photonEnergy, photonEnergyPrime, angle, cosAngle, missingEnergy, missingMomentum.x(), missingMomentum.y(), missingMomentum.z(), missingMassSquared, kaonMass, cdaAfterCut, x_gamma, x_mu);

			NTracksAfterCut++;
			m_eventCount = 0;
        }
	}
}

void Kmu2gSelection::PostProcess() {

}

void Kmu2gSelection::EndOfBurstUser() {

}

void Kmu2gSelection::EndOfRunUser() {

}

void Kmu2gSelection::EndOfJobUser() {
	SaveAllPlots();


	modeCanvas->cd();
	modeCanvas->SetGrid();
	modeScatter->SetTitle("Decay Mode Scatter;x_{#gamma};x_{#mu}");
	modeScatter->SetMarkerStyle(5);
	modeScatter->SetMarkerSize(0.5);

	modeScatter->GetXaxis()->SetLimits(-0.2, 1.2); // Set the limits for x-axis
    modeScatter->GetXaxis()->SetRangeUser(0, 1.2); // Set the viewing range for x-axis

    modeScatter->GetYaxis()->SetLimits(-0.2, 1.2); // Set the limits for y-axis
    modeScatter->GetYaxis()->SetRangeUser(0, 1.2); // Set the viewing range for y-axis

    modeCanvas->Update();

	modeScatter->Draw("AP");
	modeCanvas->SaveAs("plots/modeScatter.jpg");

}

void Kmu2gSelection::DrawPlot() {
}

Kmu2gSelection::~Kmu2gSelection() {

}


void Kmu2gSelection::initCSV(const char* directory) {
    // std::ifstream inFile(filename, std::ios::in);
    // std::string firstline;

    // // Check if the file was opened successfully
    // if (inFile.is_open()) {
    //     std::getline(inFile, firstline);
    // }
    // inFile.close();

    // std::string header = "EventCoumt,GTKMomentumX,GTKMomentumY,GTKMomentumZ,GTKEnergy,GTKGamma,GTKBeta,NLKrCells,NLKrClusters,LKrClusterEnergy,NVertices,NTracks,CDA,TrackTime,QChi2Track,EoP,MuonMomentumX,MuonMomentumY,MuonMomentumZ,MuonTrueEnergy,MuonLKrClusterEnergy,MuonMomentumPrimeX,MuonMomentumPrimeY,MuonMomentumPrimeZ,MuonEnergyPrime,PhotonEnergy,PhotonEnergyPrime,GammaMuAngle,CosGammaMuAngle,MissingEnergy,MissingMomentumX,MissingMomentumY,MissingMomentumZ,MissingMassSquared,KaonMass,CDAAfterCut,XGamma,XMu";

    // if (firstline != header) {
    //     file.open(filename, std::ofstream::trunc);
    //     file << header << std::endl;
    // }
    // else {
    //     file.open(filename, std::ofstream::app);
    // }

    while (true) {
        std::string filename = directory + std::to_string(filenumber) + ".csv";
        if (!std::filesystem::exists(filename)) {
            file.open(filename, std::ofstream::trunc);
            file << "EventCount,GTKMomentumX,GTKMomentumY,GTKMomentumZ,GTKEnergy,GTKGamma,GTKBeta,NLKrCells,NLKrClusters,LKrClusterEnergy,NVertices,NTracks,CDA,TrackTime,QChi2Track,EoP,MuonMomentumX,MuonMomentumY,MuonMomentumZ,MuonTrueEnergy,MuonLKrClusterEnergy,MuonMomentumPrimeX,MuonMomentumPrimeY,MuonMomentumPrimeZ,MuonEnergyPrime,PhotonEnergy,PhotonEnergyPrime,GammaMuAngle,CosGammaMuAngle,MissingEnergy,MissingMomentumX,MissingMomentumY,MissingMomentumZ,MissingMassSquared,KaonMass,CDAAfterCut,XGamma,XMu" << std::endl;
            break;
        }
        else {
            filenumber++;
        }
    }
}

void Kmu2gSelection::saveCSV(Int_t eventCount, Double_t GTKMomentumX, Double_t GTKMomentumY, Double_t GTKMomentumZ, Double_t GTKEnergy, Double_t GTKGamma, Double_t GTKBeta, Int_t NLKrCells, Int_t NLKrClusters, Double_t LKrClusterEnergy, Int_t NVertices, Int_t NTracks, Double_t CDA, Double_t TrackTime, Double_t QChi2Track, Double_t EoP, Double_t MuonMomentumX, Double_t MuonMomentumY, Double_t MuonMomentumZ, Double_t MuonTrueEnergy, Double_t MuonLKrClusterEnergy, Double_t MuonMomentumPrimeX, Double_t MuonMomentumPrimeY, Double_t MuonMomentumPrimeZ, Double_t MuonEnergyPrime, Double_t PhotonEnergy, Double_t PhotonEnergyPrime, Double_t GammaMuAngle, Double_t CosGammaMuAngle, Double_t MissingEnergy, Double_t MissingMomentumX, Double_t MissingMomentumY, Double_t MissingMomentumZ, Double_t MissingMassSquared, Double_t KaonMass, Double_t CDAAfterCut, Double_t XGamma, Double_t XMu) {
	if (!file.is_open()) {
		std::cout << "File is not open" << std::endl;
		return;
	}
	file << eventCount << "," << GTKMomentumX << "," << GTKMomentumY << "," << GTKMomentumZ << "," << GTKEnergy << "," << GTKGamma << "," << GTKBeta << "," << NLKrCells << "," << NLKrClusters << "," << LKrClusterEnergy << "," << NVertices << "," << NTracks << "," << CDA << "," << TrackTime  << "," << QChi2Track << "," << EoP << "," << MuonMomentumX << "," << MuonMomentumY << "," << MuonMomentumZ << "," << MuonTrueEnergy << "," << MuonLKrClusterEnergy << "," << MuonMomentumPrimeX << "," << MuonMomentumPrimeY << "," << MuonMomentumPrimeZ << "," << MuonEnergyPrime << "," << PhotonEnergy << "," << PhotonEnergyPrime << "," << GammaMuAngle << "," << CosGammaMuAngle << "," << MissingEnergy << "," << MissingMomentumX << "," << MissingMomentumY << "," << MissingMomentumZ << "," << MissingMassSquared << "," << KaonMass << "," << CDAAfterCut << "," << XGamma << "," << XMu << std::endl;
}
