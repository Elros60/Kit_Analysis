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
#include <THStack.h>
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

#include "FlowAnalysis_Helper.h"
#include "Framework/Logger.h"

using namespace std;

//______________________________________________________________________________
void FlowAnalysis_SingleMuV2(std::string FileName = "AnalysisResults.root") {

  TProfile2D *tp_corr2ref = new TProfile2D();
  TProfile2D *tp_corr4ref = new TProfile2D();
  TProfile2D *tp_corr2poi = new TProfile2D();
  TProfile2D *tp_corr4poi = new TProfile2D();

  TFile *Input_File = TFile::Open(FileName.c_str());
  THashList *list_hist =
      (THashList *)Input_File->Get("analysis-muon-selection/output");
  TList *sublist =
      (TList *)list_hist->FindObject("TrackMuon_muonLowPt210SigmaPDCA");
  tp_corr2ref =
      (TProfile2D *)sublist->FindObject("Pt_centrFT0C_Corr2REFsingle");
  tp_corr4ref =
      (TProfile2D *)sublist->FindObject("Pt_centrFT0C_Corr4REFsingle");
  tp_corr2poi =
      (TProfile2D *)sublist->FindObject("Pt_centrFT0C_Corr2POIsingle");
  tp_corr4poi =
      (TProfile2D *)sublist->FindObject("Pt_centrFT0C_Corr4POIsingle");

  FlowAnalysis_Helper *helper = new FlowAnalysis_Helper();
  TAxis *ptAxis = new TAxis();
  TAxis *centAxis = new TAxis();
  ptAxis = tp_corr2ref->GetXaxis();
  centAxis = tp_corr2ref->GetYaxis();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_pt = helper->CreateBinsFromAxis(ptAxis);
  double *Bin_cent = helper->CreateBinsFromAxis(centAxis);

  TH2D *hist_v22single = new TH2D("hist_v22single", "hist_v22single", Nbins_pt,
                                  Bin_pt, Nbins_cent, Bin_cent);
  TH2D *hist_v24single = new TH2D("hist_v24single", "hist_v24single", Nbins_pt,
                                  Bin_pt, Nbins_cent, Bin_cent);
  TH2D *hist_mult = new TH2D("Multiplicity_singlemu", "Multiplicity_singlemu",
                             Nbins_pt, Bin_pt, Nbins_cent, Bin_cent);
  for (int i = 0; i < Nbins_pt; i++) {
    for (int j = 0; j < Nbins_cent; j++) {
      int global_id = tp_corr2ref->GetBin(i + 1, j + 1);
      double val_mult = tp_corr2ref->GetBinEntries(global_id);
      double val_corr2ref = tp_corr2ref->GetBinContent(i + 1, j + 1);
      double val_corr4ref = tp_corr4ref->GetBinContent(i + 1, j + 1);
      double val_corr2poi = tp_corr2poi->GetBinContent(i + 1, j + 1);
      double val_corr4poi = tp_corr4poi->GetBinContent(i + 1, j + 1);

      double d22 = val_corr2poi;
      double v22 = d22 / pow(val_corr2ref, 0.5);

      double c24 = val_corr4ref - 2.0 * pow(val_corr2ref, 2.0);
      double d24 = val_corr4poi - 2.0 * val_corr2poi * val_corr2ref;
      double v24 = d24 / pow(-1. * c24, 3. / 4.);

      v22 = isnan(v22) || isinf(v22) ? 1000. : v22;
      v24 = isnan(v24) || isinf(v24) ? 1000. : v24;
      hist_v22single->SetBinContent(i + 1, j + 1, v22);
      //   hist_v22single->SetBinError(i + 1, j + 1, v22 <= 0. ? -1. * v22 :
      //   0.);
      hist_v24single->SetBinContent(i + 1, j + 1,
                                    isnan(v24) || isinf(v24) ? 1000. : v24);
      hist_mult->SetBinContent(i + 1, j + 1, val_mult);
    }
  }

  TFile *f = new TFile("FlowAnalysisResults_SingleMuV2.root", "RECREATE");
  f->cd();
  hist_mult->Write();
  hist_v22single->Write();
  hist_v24single->Write();
  f->Close();

  Input_File->Close();
}