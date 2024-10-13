#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TChain.h>
#include <TClass.h>
#include <TCollection.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "FlowAnalysis_Fitting.h"
#include "Framework/Logger.h"

using namespace std;

void LoadDataRun2(double *&x, double *&y, double *&ex, double *&ey,
                  double *&ey_sys, int flag);

//______________________________________________________________________________
void FlowAnalysis_CompareRun2(std::string FileName = "AnalysisResults.root") {

  // Load data from AnalysisResults.root
  TFile *Input_File = TFile::Open(FileName.c_str());
  THashList *list_hist =
      (THashList *)Input_File->Get("analysis-same-event-pairing/output");
  TList *sublist = (TList *)list_hist->FindObject(
      Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));

  // Get TProfiles of correlations
  TProfile3D *tp_Corr2Ref =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr2REF");
  TProfile3D *tp_Corr4Ref =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr4REF");
  TProfile3D *tp_Corr2Poi =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr2POI");
  TProfile3D *tp_Corr4Poi =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr4POI");

  double *x_run2_1, *y_run2_1, *ex_run2_1, *ey_run2_1, *eysys_run2_1;
  double *x_run2_2, *y_run2_2, *ex_run2_2, *ey_run2_2, *eysys_run2_2;
}