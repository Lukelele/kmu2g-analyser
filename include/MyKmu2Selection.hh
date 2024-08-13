// ---------------------------------------------------------------
//
// History:
//
// Created by Evgueni Goudzovski (eg@hep.ph.bham.ac.uk) 2016-03-25
//
// ---------------------------------------------------------------

#ifndef MYKMU2SELECTION_HH
#define MYKMU2SELECTION_HH

#include "Analyzer.hh"
#include "GeometricAcceptance.hh"
#include "MatchingRG.hh"
#include "SpectrometerTrackVertex.hh"
#include "TwoLinesCDA.hh"

#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <stdlib.h>

class MyKmu2Selection: public NA62Analysis::Analyzer {

public:
  explicit MyKmu2Selection(NA62Analysis::Core::BaseAnalysis *ba);
  ~MyKmu2Selection() {}
  void InitHist();
  void InitOutput();
  void Process(Int_t);
  void StartOfBurstUser();
  void EndOfBurstUser();
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser(); 
  void BuildPDFReport();
  void PostProcess() {}
  void DrawPlot() {}
  void PrintStatisticsPerBurst();

private:
  Int_t fTriggerMask;   ///< Definition of the data sample by L0 trigger mask
  Bool_t fTightPID;     ///< Tight PID: E/p and MUV3 association
  Bool_t fUseGTK;       ///< Whether to use GTK kaon momentum
  Bool_t fReadingData;  ///< Reading data or my own output?
  Bool_t
    fSkipWrongType;  ///< If true, do not process MC events not of Kmu2 type (used by automatic revision metrics)
  Double_t fMaxNBursts;  ///< Number of bins in the histograms of counts vs burst ID, default = 5000
  Double_t fAcceptance;  ///< Standard acceptance for on-overlaid MC
  Double_t fAcceptanceGTK;  ///< Standard acceptance for on-overlaid MC (for selection with GTK)
  Double_t fAcceptanceNew;  ///< Acceptance for the input MC
  Double_t
    fAcceptanceNewGTK;  ///< Acceptance for the input MC (with a selection involving a GTK candidate)

  TH1F *fHZtrue;
  TH1F *fHPhysicsEventsPerBurst;
  TH1F *fHPhysicsEventsGoodQualityPerBurst;
  TH1F *fHPhysicsEventsGoodQualityGTKPerBurst;
  TH1F *fHPhysicsEventsGoodQualitywGTKPerBurst;
  TH1F *fHKmu2EventsPerBurst;
  TH1F *fHMass;
  TH1F *fHEOP;
  TH1F *fHZvtx;
  TH1F *fHZvtxGTK;
  TH2F *fHxyMUV3;
  TH2F *fHMassVsP;

  MatchingRG *fMatchingRG;  ///< STRAW-GTK track matching algorithm
  TwoLinesCDA *fCDAcomp;    ///< A tool for computation of distance between track pairs

  // Outputs
  Bool_t fEventSelected;
  Double_t fMuonMomentum;
  Double_t fMissingMass2;
  Double_t fVertexZ;
  Int_t fGTKTrackID;
  Int_t fKmu2TrackID;
};
#endif
