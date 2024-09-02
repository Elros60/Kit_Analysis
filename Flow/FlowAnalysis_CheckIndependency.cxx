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
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMatrixD.h>
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

#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include "Framework/Logger.h"

using namespace std;

double *CreateBinsFromAxis(TAxis *axis);

void FlowAnalysis_CheckIndependency(
    std::string muonCut = "muonLowPt210SigmaPDCA",
    std::string dimuonCut = "pairRapidityForward",
    std::string FileName = "AnalysisResults.root") {

  // Load data from AnalysisResults.root
  TFile *Input_File = TFile::Open(FileName.c_str());
  THashList *list_hist =
      (THashList *)Input_File->Get("analysis-same-event-pairing/output");
  TList *sublist = (TList *)list_hist->FindObject(
      Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));

  // Get histograms of mult-dimuon vs. correlators
  TH2D *histCor2_xy = (TH2D *)sublist->FindObject("MultiDimuons_Corr2POI");
  TH2D *histCor4_xy = (TH2D *)sublist->FindObject("MultiDimuons_Corr4POI");

  /// Step 1.
  // Normalize 2D histograms
  histCor2_xy->Scale(1. / histCor2_xy->Integral());
  histCor4_xy->Scale(1. / histCor4_xy->Integral());

  // Step 2.
  // Projection onto 2 axes to get marginal normalized pdfs
  TH1D *histCor2_x = histCor2_xy->ProjectionX("cor2_x"); // axis of multi-dimuon
  TH1D *histCor2_y = histCor2_xy->ProjectionY("cor2_y"); // axis of correlator

  TH1D *histCor4_x = histCor4_xy->ProjectionX("cor4_x"); // axis of multi-dimuon
  TH1D *histCor4_y = histCor4_xy->ProjectionY("cor4_y"); // axis of correlator

  // Step 3.
  // Get P_xy/(P_x*P_y) histogram
  double *BinCor2_x = CreateBinsFromAxis(histCor2_xy->GetXaxis());
  double *BinCor2_y = CreateBinsFromAxis(histCor2_xy->GetYaxis());
  int NbinCor2_x = histCor2_xy->GetXaxis()->GetNbins();
  int NbinCor2_y = histCor2_xy->GetYaxis()->GetNbins();

  double *BinCor4_x = CreateBinsFromAxis(histCor4_xy->GetXaxis());
  double *BinCor4_y = CreateBinsFromAxis(histCor4_xy->GetYaxis());
  int NbinCor4_x = histCor4_xy->GetXaxis()->GetNbins();
  int NbinCor4_y = histCor4_xy->GetYaxis()->GetNbins();

  TH2D *histCor2_ratio = new TH2D("histCor2_ratio", "histCor2_ratio",
                                  NbinCor2_x, BinCor2_x, NbinCor2_y, BinCor2_y);
  TH2D *histCor4_ratio = new TH2D("histCor4_ratio", "histCor4_ratio",
                                  NbinCor4_x, BinCor4_x, NbinCor4_y, BinCor4_y);
  for (int i = 0; i < NbinCor2_x; i++) {
    for (int j = 0; j < NbinCor2_y; j++) {

      double pCor2_x = histCor2_x->GetBinContent(i + 1);
      double pCor2_y = histCor2_y->GetBinContent(j + 1);
      double pCor4_x = histCor4_x->GetBinContent(i + 1);
      double pCor4_y = histCor4_y->GetBinContent(j + 1);
      double pCor2_xy = histCor2_xy->GetBinContent(i + 1, j + 1);
      double pCor4_xy = histCor4_xy->GetBinContent(i + 1, j + 1);

      double ratioCor2 = pCor2_x * pCor2_y * pCor2_xy == 0
                             ? 1.
                             : pCor2_xy / (pCor2_x * pCor2_y);
      double ratioCor4 = pCor4_x * pCor4_y * pCor4_xy == 0
                             ? 1.
                             : pCor4_xy / (pCor4_x * pCor4_y);

      histCor2_ratio->SetBinContent(i + 1, j + 1, ratioCor2);
      histCor4_ratio->SetBinContent(i + 1, j + 1, ratioCor4);
    }
  }

  TList *lsCor2 = new TList();
  TList *lsCor4 = new TList();
  TFile f("FlowAnalysisResults_CheckIndependency.root", "RECREATE");

  lsCor2->Add(histCor2_ratio);
  lsCor2->Add(histCor2_xy);
  lsCor2->Add(histCor2_x);
  lsCor2->Add(histCor2_y);
  lsCor2->Write("Results_Corr2POI", TObject::kSingleKey);

  lsCor4->Add(histCor4_ratio);
  lsCor4->Add(histCor4_xy);
  lsCor4->Add(histCor4_x);
  lsCor4->Add(histCor4_y);
  lsCor4->Write("Results_Corr4POI", TObject::kSingleKey);

  f.Close();
}

double *CreateBinsFromAxis(TAxis *axis) {
  int Nbins = axis->GetNbins();
  double *Bins = new double[Nbins + 1];
  axis->GetLowEdge(Bins);
  Bins[Nbins] = axis->GetBinUpEdge(Nbins);
  return Bins;
}