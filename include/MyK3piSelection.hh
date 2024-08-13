#ifndef MYK3PISELECTION_HH
#define MYK3PISELECTION_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>

#include <TProfile.h>
#include "DownstreamTrack.hh"
#include "SpectrometerTrackVertex.hh"
#include "TriggerConditions.hh"
#include "GeometricAcceptance.hh"
#include "NA62Global.hh"
#include "BeamParameters.hh"
#include "KaonDecayConstants.hh"
#include "SpectrometerGigaTrackerMatching.hh"


class TH1I;
class TH2F;
class TGraph;
class TTree;


class MyK3piSelection : public NA62Analysis::Analyzer
{
public:
  explicit MyK3piSelection(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyK3piSelection();
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
protected:

private:
  Int_t       fRunID;
  Int_t       fBurstID;
  const Int_t cMaxNBursts = 5000;
  const Int_t cMinRunID   = 5435;
  const Int_t cMaxRunID   = 9462;
  SpectrometerGigaTrackerMatching *fSGMatching;
};
#endif
