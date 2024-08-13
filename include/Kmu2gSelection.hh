#ifndef KMU2GSELECTION_HH
#define KMU2GSELECTION_HH


#include <fstream>
#include <filesystem>

#include "Analyzer.hh"
#include "GeometricAcceptance.hh"
#include "MatchingRG.hh"
#include "SpectrometerTrackVertex.hh"
#include "TwoLinesCDA.hh"

#include <stdlib.h>
#include <vector>
#include <TCanvas.h>

class TH1I;
class TH2F;
class TGraph;
class TTree;


class Kmu2gSelection : public NA62Analysis::Analyzer
{
public:
    explicit Kmu2gSelection(NA62Analysis::Core::BaseAnalysis *ba);
    ~Kmu2gSelection();
    void InitHist();
    void InitOutput();
    void ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType);
    void Process(int iEvent);
    void PostProcess();
    void StartOfBurstUser();
    void EndOfBurstUser();
    void StartOfRunUser();
    void EndOfRunUser();
    void EndOfJobUser();
    void DrawPlot();

private:
    Double_t LKrZ;

	TCanvas* modeCanvas;
	TGraph* modeScatter;

	Int_t NTracksBeforeCut;
	Int_t NTracksAfterCut;
	Double_t acceptance;

	std::ofstream file;
    int filenumber;

    Int_t m_eventCount;
	
	void initCSV(const char* directory);
	void saveCSV(Int_t eventCount, Double_t GTKMomentumX, Double_t GTKMomentumY, Double_t GTKMomentumZ, Double_t GTKEnergy, Double_t GTKGamma, Double_t GTKBeta, Int_t NLKrCells, Int_t NLKrClusters, Double_t LKrClusterEnergy, Int_t NVertices, Int_t NTracks, Double_t CDA, Double_t TrackTime, Double_t QChi2Track, Double_t EoP, Double_t MuonMomentumX, Double_t MuonMomentumY, Double_t MuonMomentumZ, Double_t MuonTrueEnergy, Double_t MuonLKrClusterEnergy, Double_t MuonMomentumPrimeX, Double_t MuonMomentumPrimeY, Double_t MuonMomentumPrimeZ, Double_t MuonEnergyPrime, Double_t PhotonEnergy, Double_t PhotonEnergyPrime, Double_t GammaMuAngle, Double_t CosGammaMuAngle, Double_t MissingEnergy, Double_t MissingMomentumX, Double_t MissingMomentumY, Double_t MissingMomentumZ, Double_t MissingMassSquared, Double_t KaonMass, Double_t CDAAfterCut, Double_t XGamma, Double_t XMu);
};
#endif
