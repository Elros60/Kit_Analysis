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

#include "FlowAnalysis_Fitting.h"
#include "FlowAnalysis_Helper.h"
#include "Framework/Logger.h"

using namespace std;

//______________________________________________________________________________
void FlowAnalysis_MassOnly(int flag_sig, int flag_bkg, double mass_min = 2.3,
                           double mass_max = 4.3, double chi2max_mass = 2.,
                           std::string label = "template",
                           std::string FileName = "AnalysisResults.root",
                           std::string muonCut = "muonLowPt510SigmaPDCA") {
  // Init Helper class
  FlowAnalysis_Helper *helper = new FlowAnalysis_Helper();

  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist = new THashList();
  TList *sublist = new TList();
  TFile *Input_File = TFile::Open(FileName.c_str());
  list_hist =
      (THashList *)Input_File->Get("analysis-same-event-pairing/output");
  sublist =
      (TList *)list_hist->FindObject(Form("PairsMuonSEPM_%s", muonCut.c_str()));
  TH2F *hist_mass = new TH2F();
  hist_mass = (TH2F *)sublist->FindObject("Mass_Pt");
  Input_File->Close();
  // delete list_hist;

  // Initialize fitter
  FlowAnalysis_Fitting *fitter = new FlowAnalysis_Fitting();
  fitter->init();
  fitter->setChi2MaxMass(chi2max_mass);
  fitter->setCentRange(0, 100);

  // Create output file
  TFile *f = new TFile(
      Form("FlowAnalysisResults_MassOnly_%s.root", label.c_str()), "RECREATE");

  // Pt bins
  vector<double> Bin_pt = {0, 2, 4, 6, 8, 10, 12, 20};
  vector<double> BinCenter_pt = {1, 3, 5, 7, 9, 11, 16};
  vector<double> BinError_pt = {1, 1, 1, 1, 1, 1, 4};
  vector<double> Mean_pt;
  vector<double> MeanError_pt;
  vector<double> Width_pt;
  vector<double> WidthError_pt;

  // Run over all Pt bins for fitting
  for (int i = 0; i < int(Bin_pt.size()) - 1; i++) {
    // Get data
    TH1D *hs_mass_proj = new TH1D();
    TH2F *hs_mass_copy = new TH2F();

    hs_mass_copy = dynamic_cast<TH2F *>(hist_mass->Clone("Mass_Copy"));
    hs_mass_copy->GetXaxis()->SetRangeUser(mass_min + 1E-5, mass_max);
    hs_mass_copy->GetYaxis()->SetRangeUser(Bin_pt[i] + 1E-5, Bin_pt[i + 1]);
    hs_mass_proj = hs_mass_copy->ProjectionX("ProjX");
    // Setup fitter
    fitter->setModel(flag_sig, flag_bkg);
    fitter->setMassRange(mass_min, mass_max);
    fitter->setPtRange(Bin_pt[i], Bin_pt[i + 1]);
    fitter->setMode(0);

    // Run fitting
    TList *l_diff_fit = new TList();
    vector<double> results =
        fitter->runFittingMassOnly(hs_mass_proj, l_diff_fit);

    Mean_pt.emplace_back(results[0]);
    MeanError_pt.emplace_back(results[1]);
    Width_pt.emplace_back(results[2]);
    WidthError_pt.emplace_back(results[3]);

    f->cd();
    l_diff_fit->SetOwner();
    l_diff_fit->Write(Form("Fit_%g_%g", Bin_pt[i], Bin_pt[i + 1]),
                      TObject::kSingleKey);

    delete l_diff_fit;
  }

  TList *l_final = new TList();
  TH1D *hist_peak = new TH1D("PeakVsPt", "PeakVsPt", Bin_pt.size() - 1,
                             reinterpret_cast<double *>(Bin_pt.data()));
  TH1D *hist_width = new TH1D("WidthVsPt", "WidthVsPt", Bin_pt.size() - 1,
                              reinterpret_cast<double *>(Bin_pt.data()));
  for (int i = 0; i < Bin_pt.size() - 1; i++) {
    hist_peak->SetBinContent(i + 1, Mean_pt[i]);
    hist_peak->SetBinError(i + 1, MeanError_pt[i]);
    hist_width->SetBinContent(i + 1, Width_pt[i]);
    hist_width->SetBinError(i + 1, WidthError_pt[i]);
  }
  l_final->Add(hist_peak);
  l_final->Add(hist_width);
  f->cd();
  l_final->SetOwner();
  l_final->Write("FinalResults", TObject::kSingleKey);

  f->Close();
  delete f;
}