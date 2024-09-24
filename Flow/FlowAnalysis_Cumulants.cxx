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
double *CreateBinsFromAxis(TAxis *axis);
void CreateBins(double *axis, double min, double max, int Nbins = 10);
vector<TH1D *> GetVD2(double pt_min, double pt_max, TProfile *Corr2Ref_mass,
                      TProfile *Corr4Ref_mass, TProfile *Corr2Poi_mass,
                      TProfile *Corr4Poi_mass);
vector<TProfile *> GetProfiles(double ptmin, double ptmax, double massmin,
                               double massmax, double centmin, double centmax,
                               TProfile3D *tp_Corr2Ref, TProfile3D *tp_Corr4Ref,
                               TProfile3D *tp_Corr2Poi,
                               TProfile3D *tp_Corr4Poi);
TH1D *GetMass(double ptmin, double ptmax, double massmin, double massmax,
              double centmin, double centmax, THnSparse *hist_mass);

//______________________________________________________________________________
void FlowAnalysis_Cumulants(int flag_sig, int flag_bkg,
                            std::string FileName = "AnalysisResults.root",
                            double cent_min = 10., double cent_max = 50.,
                            double chi2max = 2., bool sys = false,
                            std::string muonCut = "muonLowPt210SigmaPDCA",
                            std::string dimuonCut = "pairRapidityForward") {
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

  // Get histograms for J/Psi invariant mass
  THnSparse *hist_mass =
      (THnSparse *)sublist->FindObject("Mass_Pt_Rapidity_CentFT0C");

  // Get binning information
  TAxis *massAxis = tp_Corr2Ref->GetXaxis();
  TAxis *ptAxis = tp_Corr2Ref->GetYaxis();
  TAxis *centAxis = tp_Corr2Ref->GetZaxis();
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_mass = CreateBinsFromAxis(massAxis);
  double *Bin_pt = CreateBinsFromAxis(ptAxis);
  double *Bin_cent = CreateBinsFromAxis(centAxis);

  // Initialization for fitting
  FlowAnalysis_Fitting fitter;
  fitter.init();
  fitter.setChi2Max(chi2max);
  fitter.setCentRange(cent_min, cent_max);

  // Define variables' range for analysis
  double Bin_pt_mass[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
  double mass_min = 2.3;
  double mass_max = 4.3;

  // Define the pool for systematics: 36
  // combinationss
  double mass_min_sys[3] = {2.2, 2.3, 2.4};
  double mass_max_sys[3] = {4.2, 4.3, 4.4};
  string sig_enum[5] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev", "VWG"};
  string bkg_v2_enum[2] = {"Pol2", "Chebychev"};
  int sig_mass[1] = {0};    // CB2(MC,data) NA60
  int bkg_mass[2] = {3, 4}; // Chebychev VWG
  int bkg_v2[1] = {0};      // Pol2 and Chebychev

  // Create output file
  TFile f(sys ? Form("FlowAnalysisResults_"
                     "Cumulants_%g_%"
                     "g_withSys.root",
                     cent_min, cent_max)
              : Form("FlowAnalysisResults_"
                     "Cumulants_%g_%"
                     "g.root",
                     cent_min, cent_max),
          "RECREATE");

  ///////////////////////////////////////
  ///                                 ///
  ///   Analysis for Reference Flow   ///
  ///                                 ///
  ///////////////////////////////////////

  // Define projected profiles w.r.t above
  // defined variables' ranges
  TList *l_refflow = new TList();
  TProfile *tp_Corr2Ref_cent = tp_Corr2Ref->Project3DProfile("xz")->ProfileX();
  TProfile *tp_Corr4Ref_cent = tp_Corr4Ref->Project3DProfile("xz")->ProfileX();

  TH1D *hist_c22ref =
      new TH1D("c22ref", "c^{REF}_{2}{2}", Nbins_cent, Bin_cent);
  hist_c22ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_c22ref->GetYaxis()->SetTitle("c^{REF}_{2}{2}");
  TH1D *hist_c24ref =
      new TH1D("c24ref", "c^{REF}_{2}{4}", Nbins_cent, Bin_cent);
  hist_c24ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_c24ref->GetYaxis()->SetTitle("c^{REF}_{2}{4}");
  TH1D *hist_v22ref =
      new TH1D("v22ref", "v^{REF}_{2}{2}", Nbins_cent, Bin_cent);
  hist_v22ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_v22ref->GetYaxis()->SetTitle("v^{REF}_{2}{2}");
  TH1D *hist_v24ref =
      new TH1D("v24ref", "v^{REF}_{2}{4}", Nbins_cent, Bin_cent);
  hist_v24ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_v24ref->GetYaxis()->SetTitle("v^{REF}_{2}{4}");

  // Evaluation for Cn{n} and Vn{n}
  for (int i = 0; i < Nbins_cent; i++) {
    double corr2 = tp_Corr2Ref_cent->GetBinContent(i + 1);
    double corr2e = tp_Corr2Ref_cent->GetBinError(i + 1);
    double corr4 = tp_Corr4Ref_cent->GetBinContent(i + 1);
    double corr4e = tp_Corr4Ref_cent->GetBinError(i + 1);
    double c22 = corr2;
    double c22e = corr2e;
    double c24 = corr4 - 2 * TMath::Power(c22, 2);
    double c24e = corr4e * corr4e + 16 * corr2 * corr2 * corr2e * corr2e;
    c24e = c24e > 0 ? TMath::Sqrt(c24e) : 0;
    hist_c22ref->SetBinContent(i + 1, isnan(c22) || isinf(c22) ? 0 : c22);
    hist_c22ref->SetBinError(i + 1, isnan(c22e) || isinf(c22e) ? 0 : c22e);
    hist_c24ref->SetBinContent(i + 1, isnan(c24) || isinf(c24) ? 0 : c24);
    hist_c24ref->SetBinError(i + 1, isnan(c24e) || isinf(c24e) ? 0 : c24e);

    double v22 = c22 > 0 ? TMath::Sqrt(c22) : 0;
    double v22e = c22 > 0 ? 0.5 * c22e / TMath::Sqrt(c22) : 0;
    double v24 = c24 < 0 ? TMath::Power(-c24, 1. / 4.) : 0;
    double v24e = c24 < 0 ? TMath::Power(-c24, (-3. / 4)) * c24e / 4 : 0;
    hist_v22ref->SetBinContent(i + 1, isnan(v22) || isinf(v22) ? 0 : v22);
    hist_v22ref->SetBinError(i + 1, isnan(v22e) || isinf(v22e) ? 0 : v22e);
    hist_v24ref->SetBinContent(i + 1, isnan(v24) || isinf(v24) ? 0 : v24);
    hist_v24ref->SetBinError(i + 1, isnan(v24e) || isinf(v24e) ? 0 : v24e);
  }

  l_refflow->Add(hist_c22ref);
  l_refflow->Add(hist_c24ref);
  l_refflow->Add(hist_v22ref);
  l_refflow->Add(hist_v24ref);
  f.cd();
  l_refflow->Write("ReferenceFlow", TObject::kSingleKey);

  /////////////////////////////////////////////////
  ///                                           ///
  ///   Analysis for Differential Flow of POI ///
  ///                                           ///
  /////////////////////////////////////////////////

  // Create histogram for pt-differential v2
  TList *l_results = new TList();
  const int nbCombo = int(size(sig_mass)) * int(size(bkg_mass)) *
                      int(size(bkg_v2)) * int(size(mass_max_sys));
  double x_yield[10], y_yield[10], ex_yield[10], ey_yield[10];
  double x_v22pt[10], y_v22pt[10], ex_v22pt[10], ey_v22pt[10];
  double x_v24pt[10], y_v24pt[10], ex_v24pt[10], ey_v24pt[10];
  double *x_run2_1, *y_run2_1, *ex_run2_1, *ey_run2_1, *eysys_run2_1;
  double *x_run2_2, *y_run2_2, *ex_run2_2, *ey_run2_2, *eysys_run2_2;

  vector<double *> x_sys_yield, y_sys_yield, ey_sys_yield, x_sys_v22, y_sys_v22,
      ey_sys_v22, x_sys_v24, y_sys_v24, ey_sys_v24;

  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    x_sys_yield.emplace_back(new double[nbCombo]);
    y_sys_yield.emplace_back(new double[nbCombo]);
    ey_sys_yield.emplace_back(new double[nbCombo]);

    x_sys_v22.emplace_back(new double[nbCombo]);
    y_sys_v22.emplace_back(new double[nbCombo]);
    ey_sys_v22.emplace_back(new double[nbCombo]);

    x_sys_v24.emplace_back(new double[nbCombo]);
    y_sys_v24.emplace_back(new double[nbCombo]);
    ey_sys_v24.emplace_back(new double[nbCombo]);
  }

  LoadDataRun2(x_run2_1, y_run2_1, ex_run2_1, ey_run2_1, eysys_run2_1, 1);
  LoadDataRun2(x_run2_2, y_run2_2, ex_run2_2, ey_run2_2, eysys_run2_2, 2);

  TGraphMultiErrors *graph_v2pt_run2_1 =
      new TGraphMultiErrors("v2_pt_run2_1", "", 10, x_run2_1, y_run2_1,
                            ex_run2_1, ex_run2_1, ey_run2_1, ey_run2_1);
  graph_v2pt_run2_1->AddYError(10, eysys_run2_1, eysys_run2_1);
  graph_v2pt_run2_1->SetTitle("#sqrt{#it{s}_{NN}} = 5.02 TeV, 10-30%");
  graph_v2pt_run2_1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt_run2_1->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

  TGraphMultiErrors *graph_v2pt_run2_2 =
      new TGraphMultiErrors("v2_pt_run2_2", "", 10, x_run2_2, y_run2_2,
                            ex_run2_2, ex_run2_2, ey_run2_2, ey_run2_2);
  graph_v2pt_run2_2->AddYError(10, eysys_run2_2, eysys_run2_2);
  graph_v2pt_run2_2->SetTitle("#sqrt{#it{s}_{NN}} = 5.02 TeV, 30-50%");
  graph_v2pt_run2_2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt_run2_2->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    // Get mass and v2 to fit
    TH1D *hist_mass_proj = dynamic_cast<TH1D *>(
        GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                cent_min, cent_max, hist_mass));

    vector<TProfile *> results_profile = GetProfiles(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, tp_Corr2Ref, tp_Corr4Ref, tp_Corr2Poi, tp_Corr4Poi);
    TProfile *tp_Corr2Ref_mass = dynamic_cast<TProfile *>(results_profile[0]);
    TProfile *tp_Corr4Ref_mass = dynamic_cast<TProfile *>(results_profile[1]);
    TProfile *tp_Corr2Poi_mass = dynamic_cast<TProfile *>(results_profile[2]);
    TProfile *tp_Corr4Poi_mass = dynamic_cast<TProfile *>(results_profile[3]);

    vector<TH1D *> results =
        GetVD2(Bin_pt_mass[i], Bin_pt_mass[i + 1], tp_Corr2Ref_mass,
               tp_Corr4Ref_mass, tp_Corr2Poi_mass, tp_Corr4Poi_mass);
    TH1D *hist_d22 = dynamic_cast<TH1D *>(results[0]);
    TH1D *hist_d24 = dynamic_cast<TH1D *>(results[1]);
    TH1D *hist_vd22 = dynamic_cast<TH1D *>(results[2]);
    TH1D *hist_vd24 = dynamic_cast<TH1D *>(results[3]);
    /// Do fitting
    // Configuration for fitting
    fitter.setModel(flag_sig, flag_bkg);
    fitter.setMassRange(mass_min, mass_max);

    // Fit invariant mass + v2
    fitter.setOrder(2);
    fitter.setMode(0);
    TList *l_diff_fit2 = new TList();
    vector<double> results_v22 =
        fitter.runFitting(hist_mass_proj, hist_vd22, l_diff_fit2,
                          Bin_pt_mass[i], Bin_pt_mass[i + 1]);

    fitter.setOrder(4);
    fitter.setMode(0);
    TList *l_diff_fit4 = new TList();
    vector<double> results_v24 =
        fitter.runFitting(hist_mass_proj, hist_vd24, l_diff_fit4,
                          Bin_pt_mass[i], Bin_pt_mass[i + 1]);

    // Run fittings for systematics
    if (sys) {
      int index_sys = 0;
      for (int i1 = 0; i1 < int(size(mass_max_sys)); i1++) {
        for (int i2 = 0; i2 < int(size(sig_mass)); i2++) {
          for (int i3 = 0; i3 < int(size(bkg_mass)); i3++) {
            for (int i4 = 0; i4 < int(size(bkg_v2)); i4++) {
              // Get mass and v2 to fit
              TH1D *hist_mass_proj_sys = dynamic_cast<TH1D *>(
                  GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                          mass_max_sys[i1], cent_min, cent_max, hist_mass));

              vector<TProfile *> results_profile_sys = GetProfiles(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                  mass_max_sys[i1], cent_min, cent_max, tp_Corr2Ref,
                  tp_Corr4Ref, tp_Corr2Poi, tp_Corr4Poi);
              TProfile *tp_Corr2Ref_mass_sys =
                  dynamic_cast<TProfile *>(results_profile_sys[0]);
              TProfile *tp_Corr4Ref_mass_sys =
                  dynamic_cast<TProfile *>(results_profile_sys[1]);
              TProfile *tp_Corr2Poi_mass_sys =
                  dynamic_cast<TProfile *>(results_profile_sys[2]);
              TProfile *tp_Corr4Poi_mass_sys =
                  dynamic_cast<TProfile *>(results_profile_sys[3]);

              vector<TH1D *> results_sys =
                  GetVD2(Bin_pt_mass[i], Bin_pt_mass[i + 1],
                         tp_Corr2Ref_mass_sys, tp_Corr4Ref_mass_sys,
                         tp_Corr2Poi_mass_sys, tp_Corr4Poi_mass_sys);
              TH1D *hist_d22 = dynamic_cast<TH1D *>(results[0]);
              TH1D *hist_d24 = dynamic_cast<TH1D *>(results[1]);
              TH1D *hist_vd22 = dynamic_cast<TH1D *>(results[2]);
              TH1D *hist_vd24 = dynamic_cast<TH1D *>(results[3]);

              LOG(info) << Form(
                  "{%s,%s}+{%s}+[%g-%g]", sig_enum[sig_mass[i2]].c_str(),
                  sig_enum[bkg_mass[i3]].c_str(), bkg_v2_enum[i4].c_str(),
                  mass_min_sys[i1], mass_max_sys[i1]);

              // Do evaluation for v22
              fitter.setOrder(2);
              fitter.setMassRange(mass_min_sys[i1], mass_max_sys[i1]);
              fitter.setModel(sig_mass[i2], bkg_mass[i3]);
              fitter.setModelV2(bkg_v2[i4]);
              fitter.setMode(1);
              TList *l_diff_sys2 = new TList();
              vector<double> results_sys_v22 =
                  fitter.runFitting(hist_mass_proj_sys, hist_vd22, l_diff_sys2,
                                    Bin_pt_mass[i], Bin_pt_mass[i + 1]);

              // Do evaluation for v24
              fitter.setOrder(4);
              fitter.setMassRange(mass_min_sys[i1], mass_max_sys[i1]);
              fitter.setModel(sig_mass[i2], bkg_mass[i3]);
              fitter.setModelV2(bkg_v2[i4]);
              fitter.setMode(1);
              TList *l_diff_sys4 = new TList();
              vector<double> results_sys_v24 =
                  fitter.runFitting(hist_mass_proj_sys, hist_vd24, l_diff_sys4,
                                    Bin_pt_mass[i], Bin_pt_mass[i + 1]);

              // Fill pT-differential v2 and jpsi
              // yields
              x_sys_yield[i][index_sys] = index_sys;
              x_sys_v22[i][index_sys] = index_sys;
              x_sys_v24[i][index_sys] = index_sys;

              y_sys_yield[i][index_sys] =
                  results_sys_v22[2] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
              y_sys_v22[i][index_sys] = results_sys_v22[0];
              y_sys_v24[i][index_sys] = results_sys_v24[0];

              ey_sys_yield[i][index_sys] =
                  results_sys_v22[3] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
              ey_sys_v22[i][index_sys] = results_sys_v22[1];
              ey_sys_v24[i][index_sys] = results_sys_v24[1];

              f.cd();
              l_diff_sys2->Write(Form("FitSys_Combo%d_Fit2_%g_%"
                                      "g",
                                      index_sys, Bin_pt_mass[i],
                                      Bin_pt_mass[i + 1]),
                                 TObject::kSingleKey);
              l_diff_sys4->Write(Form("FitSys_Combo%d_Fit4_%g_%"
                                      "g",
                                      index_sys, Bin_pt_mass[i],
                                      Bin_pt_mass[i + 1]),
                                 TObject::kSingleKey);

              index_sys++;
              delete l_diff_sys2;
              delete l_diff_sys4;
            }
          }
        }
      }
    }

    // Fill pT-differential v2 and jpsi yields
    x_v22pt[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    y_v22pt[i] = results_v22[0];
    ex_v22pt[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;
    ey_v22pt[i] = results_v22[1];

    x_v24pt[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    y_v24pt[i] = results_v24[0];
    ex_v24pt[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;
    ey_v24pt[i] = results_v24[1];

    x_yield[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    y_yield[i] = results_v22[2] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
    ex_yield[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;
    ey_yield[i] = results_v22[3] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);

    // Save results
    TList *l_diff_hist = new TList();
    l_diff_hist->Add(hist_d22);
    l_diff_hist->Add(hist_d24);
    l_diff_hist->Add(hist_vd22);
    l_diff_hist->Add(hist_vd24);
    f.cd();
    l_diff_hist->Write(
        Form("DifferentialFlow_Hist_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    l_diff_fit2->Write(
        Form("DifferentialFlow_Fit2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    l_diff_fit4->Write(
        Form("DifferentialFlow_Fit4_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    delete l_diff_hist;
    delete l_diff_fit2;
    delete l_diff_fit4;
  }

  // Compare with Run2 data: 5.02 TeV, 10-30% and
  // 30-50%
  TGraphMultiErrors *graph_v22pt =
      new TGraphMultiErrors("v22_pt", "", 10, x_v22pt, y_v22pt, ex_v22pt,
                            ex_v22pt, ey_v22pt, ey_v22pt);
  graph_v22pt->SetTitle("#it{v}^{J/#psi}_{2}{2} "
                        "#sqrt{#it{s}_{NN}} = "
                        "5.36 TeV, 10-50%");
  graph_v22pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v22pt->GetYaxis()->SetTitle("v^{J/#psi}_{2}{2}}");

  TGraphMultiErrors *graph_v24pt =
      new TGraphMultiErrors("v24_pt", "", 10, x_v24pt, y_v24pt, ex_v24pt,
                            ex_v24pt, ey_v24pt, ey_v24pt);
  graph_v24pt->SetTitle("#it{v}^{J/#psi}_{2}{4} "
                        "#sqrt{#it{s}_{NN}} = "
                        "5.36 TeV, 10-50%");
  graph_v24pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v24pt->GetYaxis()->SetTitle("v^{J/#psi}_{2}{4}");

  TGraphMultiErrors *graph_yield =
      new TGraphMultiErrors("yields_pt", "", 10, x_yield, y_yield, ex_yield,
                            ex_yield, ey_yield, ey_yield);
  graph_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_yield->GetYaxis()->SetTitle("Raw counts / d#it{p}_{T} (GeV/c)^{-1}");

  // Save final results into graphs
  graph_yield->SetMarkerStyle(20);
  graph_yield->SetMarkerSize(1.);
  graph_yield->SetMarkerColor(kBlue);
  graph_yield->SetLineColor(kBlue);
  graph_yield->SetLineWidth(2);
  graph_yield->SetFillStyle(0);
  // graph_yield->Draw("A P Z");
  l_results->Add(graph_yield);
  /*
  TLatex *text_yield = new TLatex();
  text_yield->SetTextSize(0.04);
  text_yield->SetTextFont(42);
  text_yield->DrawLatexNDC(
      .3, .82, "ALICE Performance, Pb-Pb
  #sqrt{#it{s}_{NN}} = 5.36 TeV");
  gPad->ModifiedUpdate();
  text_yield->DrawLatexNDC(.3, .77,
                           "J/#psi#rightarrow#mu^{+}#mu^{-},
  2.5 < y < 4"); gPad->ModifiedUpdate();
  */

  auto mg_2vs4 = new TMultiGraph("v2_pt_2vs4", "");
  graph_v22pt->SetMarkerStyle(20);
  graph_v22pt->SetMarkerSize(1.);
  graph_v22pt->SetMarkerColor(kBlue);
  graph_v22pt->SetLineColor(kBlue);
  graph_v22pt->SetLineWidth(2);
  graph_v22pt->SetFillStyle(0);
  graph_v24pt->SetMarkerStyle(20);
  graph_v24pt->SetMarkerSize(1.);
  graph_v24pt->SetMarkerColor(kRed);
  graph_v24pt->SetLineColor(kRed);
  graph_v24pt->SetLineWidth(2);
  graph_v24pt->SetFillStyle(0);
  mg_2vs4->Add(graph_v22pt);
  mg_2vs4->Add(graph_v24pt);
  // mg_2vs4->Draw("A P Z ; Z ; 5 s=0.5");
  l_results->Add(graph_v22pt);
  l_results->Add(graph_v24pt);
  l_results->Add(mg_2vs4);
  /*
  TLatex *text_pt_2vs4 = new TLatex();
  text_pt_2vs4->SetTextSize(0.04);
  text_pt_2vs4->SetTextFont(42);
  text_pt_2vs4->DrawLatexNDC(
      .18, .82, "ALICE Performance, Pb-Pb
  #sqrt{#it{s}_{NN}} = 5.36 TeV");
  gPad->ModifiedUpdate();
  text_pt_2vs4->DrawLatexNDC(.18, .77,
                             "J/#psi#rightarrow#mu^{+}#mu^{-},
  2.5 < y < 4"); gPad->ModifiedUpdate();
  */

  auto mg_run2 = new TMultiGraph("v2_pt_run2", "");
  graph_v2pt_run2_1->SetMarkerStyle(20);
  graph_v2pt_run2_1->SetMarkerSize(1.);
  graph_v2pt_run2_1->SetMarkerColor(kOrange);
  graph_v2pt_run2_1->SetLineColor(kOrange);
  graph_v2pt_run2_1->SetLineWidth(2);
  graph_v2pt_run2_1->SetFillStyle(0);
  graph_v2pt_run2_2->SetMarkerStyle(20);
  graph_v2pt_run2_2->SetMarkerSize(1.);
  graph_v2pt_run2_2->SetMarkerColor(kBlack);
  graph_v2pt_run2_2->SetLineColor(kBlack);
  graph_v2pt_run2_2->SetLineWidth(2);
  graph_v2pt_run2_2->SetFillStyle(0);
  mg_run2->Add(graph_v22pt);
  mg_run2->Add(graph_v24pt);
  mg_run2->Add(graph_v2pt_run2_1);
  mg_run2->Add(graph_v2pt_run2_2);
  // mg_run2->Draw("A P Z ; Z ; 5 s=0.5");
  l_results->Add(mg_run2);

  /*
  TLatex *text_pt_run2 = new TLatex();
  text_pt_run2->SetTextSize(0.04);
  text_pt_run2->SetTextFont(42);
  text_pt_run2->DrawLatexNDC(
      .18, .82, "ALICE Performance, Pb-Pb
  #sqrt{#it{s}_{NN}} = 5.36 TeV");
  gPad->ModifiedUpdate();
  text_pt_run2->DrawLatexNDC(.18, .77,
                             "J/#psi#rightarrow#mu^{+}#mu^{-},
  2.5 < y < 4"); gPad->ModifiedUpdate();
  */
  if (sys) {
    vector<TGraphMultiErrors *> graph_sys_yield;
    vector<TGraphMultiErrors *> graph_sys_v22;
    vector<TGraphMultiErrors *> graph_sys_v24;
    for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
      graph_sys_yield.emplace_back(new TGraphMultiErrors(
          Form("sys_yield_%d", i), "", nbCombo, x_sys_yield[i], y_sys_yield[i],
          nullptr, nullptr, ey_sys_yield[i], ey_sys_yield[i]));
      graph_sys_v22.emplace_back(new TGraphMultiErrors(
          Form("sys_v22_%d", i), "", nbCombo, x_sys_v22[i], y_sys_v22[i],
          nullptr, nullptr, ey_sys_v22[i], ey_sys_v22[i]));
      graph_sys_v24.emplace_back(new TGraphMultiErrors(
          Form("sys_v24_%d", i), "", nbCombo, x_sys_v24[i], y_sys_v24[i],
          nullptr, nullptr, ey_sys_v24[i], ey_sys_v24[i]));
    }

    for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
      graph_sys_v22[i]->SetMarkerStyle(20);
      graph_sys_v22[i]->SetMarkerSize(1.);
      graph_sys_v22[i]->SetMarkerColor(kBlue);
      graph_sys_v22[i]->SetLineColor(kBlue);
      graph_sys_v22[i]->SetLineWidth(2);
      graph_sys_v22[i]->SetFillStyle(0);

      graph_sys_v24[i]->SetMarkerStyle(20);
      graph_sys_v24[i]->SetMarkerSize(1.);
      graph_sys_v24[i]->SetMarkerColor(kBlue);
      graph_sys_v24[i]->SetLineColor(kBlue);
      graph_sys_v24[i]->SetLineWidth(2);
      graph_sys_v24[i]->SetFillStyle(0);

      graph_sys_yield[i]->SetMarkerStyle(20);
      graph_sys_yield[i]->SetMarkerSize(1.);
      graph_sys_yield[i]->SetMarkerColor(kBlue);
      graph_sys_yield[i]->SetLineColor(kBlue);
      graph_sys_yield[i]->SetLineWidth(2);
      graph_sys_yield[i]->SetFillStyle(0);

      l_results->Add(graph_sys_v22[i]);
      l_results->Add(graph_sys_v24[i]);
      l_results->Add(graph_sys_yield[i]);
    }
  }

  f.cd();
  l_results->Write("FinalResults", TObject::kSingleKey);

  f.Close();

  f.Close();
  Input_File->Close();
}

//______________________________________________________________________________
double *CreateBinsFromAxis(TAxis *axis) {
  int Nbins = axis->GetNbins();
  double *Bins = new double[Nbins + 1];
  axis->GetLowEdge(Bins);
  Bins[Nbins] = axis->GetBinUpEdge(Nbins);
  return Bins;
}

//______________________________________________________________________________
void CreateBins(double *axis, double min, double max, int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}

//______________________________________________________________________________
TH1D *GetMass(double ptmin, double ptmax, double massmin, double massmax,
              double centmin, double centmax, THnSparse *hist_mass) {
  THnSparse *hist_mass_cp = dynamic_cast<THnSparse *>(hist_mass->Clone(
      Form("Mass_Pt_Rapidity_centrFT0C_Copy_%g_%g", ptmin, ptmax)));
  hist_mass_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hist_mass_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hist_mass_cp->GetAxis(3)->SetRangeUser(centmin, centmax);
  TH1D *hist_mass_proj = hist_mass_cp->Projection(0);
  return hist_mass_proj;
}

//______________________________________________________________________________
vector<TProfile *> GetProfiles(double ptmin, double ptmax, double massmin,
                               double massmax, double centmin, double centmax,
                               TProfile3D *tp_Corr2Ref, TProfile3D *tp_Corr4Ref,
                               TProfile3D *tp_Corr2Poi,
                               TProfile3D *tp_Corr4Poi) {
  vector<TProfile *> results;
  // Copy original profiles for projections
  TProfile3D *tp_Corr2Ref_cp = dynamic_cast<TProfile3D *>(tp_Corr2Ref->Clone(
      Form("Mass_Pt_centrFT0C_Corr2REF_Copy_%g_%g", ptmin, ptmax)));
  TProfile3D *tp_Corr4Ref_cp = dynamic_cast<TProfile3D *>(tp_Corr4Ref->Clone(
      Form("Mass_Pt_centrFT0C_Corr4REF_Copy_%g_%g", ptmin, ptmax)));
  TProfile3D *tp_Corr2Poi_cp = dynamic_cast<TProfile3D *>(tp_Corr2Poi->Clone(
      Form("Mass_Pt_centrFT0C_Corr2POI_Copy_%g_%g", ptmin, ptmax)));
  TProfile3D *tp_Corr4Poi_cp = dynamic_cast<TProfile3D *>(tp_Corr4Poi->Clone(
      Form("Mass_Pt_centrFT0C_Corr4POI_Copy_%g_%g", ptmin, ptmax)));

  // Set axes' ranges for mass-differential study
  tp_Corr2Ref_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_Corr2Ref_cp->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_Corr2Ref_cp->GetZaxis()->SetRangeUser(centmin, centmax);

  tp_Corr4Ref_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_Corr4Ref_cp->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_Corr4Ref_cp->GetZaxis()->SetRangeUser(centmin, centmax);

  tp_Corr2Poi_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_Corr2Poi_cp->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_Corr2Poi_cp->GetZaxis()->SetRangeUser(centmin, centmax);

  tp_Corr4Poi_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_Corr4Poi_cp->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_Corr4Poi_cp->GetZaxis()->SetRangeUser(centmin, centmax);

  // Define projected profiles w.r.t above
  // defined Variables' ranges
  TProfile *tp_Corr2Ref_mass =
      tp_Corr2Ref_cp->Project3DProfile("yx")->ProfileX();
  TProfile *tp_Corr4Ref_mass =
      tp_Corr4Ref_cp->Project3DProfile("yx")->ProfileX();
  TProfile *tp_Corr2Poi_mass =
      tp_Corr2Poi_cp->Project3DProfile("yx")->ProfileX();
  TProfile *tp_Corr4Poi_mass =
      tp_Corr4Poi_cp->Project3DProfile("yx")->ProfileX();

  results.emplace_back(tp_Corr2Ref_mass);
  results.emplace_back(tp_Corr4Ref_mass);
  results.emplace_back(tp_Corr2Poi_mass);
  results.emplace_back(tp_Corr4Poi_mass);
  return results;
}

//______________________________________________________________________________
vector<TH1D *> GetVD2(double pt_min, double pt_max, TProfile *Corr2Ref_mass,
                      TProfile *Corr4Ref_mass, TProfile *Corr2Poi_mass,
                      TProfile *Corr4Poi_mass) {
  vector<TH1D *> results;
  double *Bin_mass_new = CreateBinsFromAxis(Corr2Poi_mass->GetXaxis());
  int NBins_mass_new = Corr2Poi_mass->GetXaxis()->GetNbins();

  // Define histograms
  TH1D *hist_d22 = new TH1D(Form("d22_%g_%g", pt_min, pt_max),
                            Form("d^{#mu#mu}_{2}{2}_%g_%g", pt_min, pt_max),
                            NBins_mass_new, Bin_mass_new);
  hist_d22->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_d22->GetYaxis()->SetTitle("d^{#mu#mu}_{2}{2}");

  TH1D *hist_d24 = new TH1D(Form("d24_%g_%g", pt_min, pt_max),
                            Form("d^{#mu#mu}_{2}{4}_%g_%g", pt_min, pt_max),
                            NBins_mass_new, Bin_mass_new);
  hist_d24->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_d24->GetYaxis()->SetTitle("d^{#mu#mu}_{2}{4}");

  TH1D *hist_vd22 = new TH1D(Form("vd22_%g_%g", pt_min, pt_max),
                             Form("v'^{#mu#mu}_{2}{2}_%g_%g", pt_min, pt_max),
                             NBins_mass_new, Bin_mass_new);
  hist_vd22->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_vd22->GetYaxis()->SetTitle("v'^{#mu#mu}_{2}{2}");

  TH1D *hist_vd24 = new TH1D(Form("vd24_%g_%g", pt_min, pt_max),
                             Form("v'^{#mu#mu}_{2}{4}_%g_%g", pt_min, pt_max),
                             NBins_mass_new, Bin_mass_new);
  hist_vd24->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_vd24->GetYaxis()->SetTitle("v'^{#mu#mu}_{2}{4}");

  for (int j = 0; j < NBins_mass_new; j++) {
    double corr2_ref = Corr2Ref_mass->GetBinContent(j + 1);
    double corr2e_ref = Corr2Ref_mass->GetBinError(j + 1);
    double corr4_ref = Corr4Ref_mass->GetBinContent(j + 1);
    double corr4e_ref = Corr4Ref_mass->GetBinError(j + 1);
    double corr2_poi = Corr2Poi_mass->GetBinContent(j + 1);
    double corr2e_poi = Corr2Poi_mass->GetBinError(j + 1);
    double corr4_poi = Corr4Poi_mass->GetBinContent(j + 1);
    double corr4e_poi = Corr4Poi_mass->GetBinError(j + 1);

    double c22 = corr2_ref;
    double c22e = corr2e_ref;
    double c24 = corr4_ref - 2 * TMath::Power(c22, 2);
    double c24e = corr4e_ref * corr4e_ref +
                  16 * corr2_ref * corr2_ref * corr2e_ref * corr2e_ref;
    c24e = c24e > 0 ? TMath::Sqrt(c24e) : 0;

    double d22 = corr2_poi;
    double d22e = corr2e_poi;
    double d24 = corr4_poi - 2 * corr2_ref * corr2_poi;
    double d24e = corr4e_poi * corr4e_ref +
                  4 * corr2_ref * corr2_ref * corr2e_poi * corr2e_poi +
                  4 * corr2_poi * corr2_poi * corr2e_ref * corr2e_ref;
    d24e = d24e > 0 ? TMath::Sqrt(d24e) : 0;
    hist_d22->SetBinContent(j + 1, isnan(d22) || isinf(d22) ? 0 : d22);
    hist_d22->SetBinError(j + 1, isnan(d22e) || isinf(d22e) ? 0 : d22e);
    hist_d24->SetBinContent(j + 1, isnan(d24) || isinf(d24) ? 0 : d24);
    hist_d24->SetBinError(j + 1, isnan(d24e) || isinf(d24e) ? 0 : d24e);

    double vd22 = c22 > 0 ? d22 / TMath::Sqrt(c22) : 0;
    double vd22e = corr2e_poi * corr2e_poi / corr2_ref +
                   0.25 * corr2_poi * corr2_poi * corr2e_ref * corr2e_ref /
                       (corr2_ref * corr2_ref * corr2_ref);
    vd22e = vd22e > 0 ? TMath::Sqrt(vd22e) : 0;
    double vd24 = c24 < 0 ? -d24 / TMath::Power(-c24, (3. / 4)) : 0;
    double vd24e = c24 < 0
                       ? TMath::Sqrt(TMath::Power(-c24, -6. / 4) * d24e * d24e +
                                     TMath::Power(-c24, -14. / 4) * d24 * d24 *
                                         c24e * c24e * 9. / 16)
                       : 0;
    hist_vd22->SetBinContent(j + 1, isnan(vd22) || isinf(vd22) ? 0 : vd22);
    hist_vd22->SetBinError(j + 1, isnan(vd22e) || isinf(vd22e) ? 0 : vd22e);
    hist_vd24->SetBinContent(j + 1, isnan(vd24) || isinf(vd24) ? 0 : vd24);
    hist_vd24->SetBinError(j + 1, isnan(vd24e) || isinf(vd24e) ? 0 : vd24e);
  }
  results.emplace_back(hist_d22);
  results.emplace_back(hist_d24);
  results.emplace_back(hist_vd22);
  results.emplace_back(hist_vd24);
  return results;
}

//______________________________________________________________________________
void LoadDataRun2(double *&x, double *&y, double *&ex, double *&ey,
                  double *&ey_sys, int flag) {
  x = new double[10];
  y = new double[10];
  ex = new double[10];
  ey = new double[10];
  ey_sys = new double[10];
  if (flag == 1) {
    // 10-30%
    x[0] = 0.64;
    x[1] = 1.49;
    x[2] = 2.47;
    x[3] = 3.46;
    x[4] = 4.45;
    x[5] = 5.45;
    x[6] = 6.819;
    x[7] = 8.835;
    x[8] = 10.84;
    x[9] = 14.25;

    y[0] = 0.011;
    y[1] = 0.043;
    y[2] = 0.074;
    y[3] = 0.088;
    y[4] = 0.085;
    y[5] = 0.103;
    y[6] = 0.083;
    y[7] = 0.1;
    y[8] = 0.049;
    y[9] = 0.022;

    ex[0] = 0.5;
    ex[1] = 0.5;
    ex[2] = 0.5;
    ex[3] = 0.5;
    ex[4] = 0.5;
    ex[5] = 0.5;
    ex[6] = 1.;
    ex[7] = 1.;
    ex[8] = 1.;
    ex[9] = 1.5;

    ey[0] = 0.0085;
    ey[1] = 0.0069;
    ey[2] = 0.0069;
    ey[3] = 0.0077;
    ey[4] = 0.009;
    ey[5] = 0.011;
    ey[6] = 0.011;
    ey[7] = 0.018;
    ey[8] = 0.028;
    ey[9] = 0.032;

    ey_sys[0] = 0.0038391;
    ey_sys[1] = 0.0036633;
    ey_sys[2] = 0.004898;
    ey_sys[3] = 0.0035068;
    ey_sys[4] = 0.0037855;
    ey_sys[5] = 0.0029726;
    ey_sys[6] = 0.0036802;
    ey_sys[7] = 0.0075789;
    ey_sys[8] = 0.0093488;
    ey_sys[9] = 0.0091828;
  } else if (flag == 2) {
    // 30-50%
    x[0] = 0.64;
    x[1] = 1.49;
    x[2] = 2.47;
    x[3] = 3.46;
    x[4] = 4.45;
    x[5] = 5.45;
    x[6] = 6.819;
    x[7] = 8.835;
    x[8] = 10.84;
    x[9] = 14.25;

    y[0] = 0.0008;
    y[1] = 0.029;
    y[2] = 0.067;
    y[3] = 0.099;
    y[4] = 0.098;
    y[5] = 0.101;
    y[6] = 0.098;
    y[7] = 0.092;
    y[8] = 0.055;
    y[9] = 0.026;

    ex[0] = 0.5;
    ex[1] = 0.5;
    ex[2] = 0.5;
    ex[3] = 0.5;
    ex[4] = 0.5;
    ex[5] = 0.5;
    ex[6] = 1.;
    ex[7] = 1.;
    ex[8] = 1.;
    ex[9] = 1.5;

    ey[0] = 0.011;
    ey[1] = 0.0091;
    ey[2] = 0.0089;
    ey[3] = 0.0095;
    ey[4] = 0.011;
    ey[5] = 0.013;
    ey[6] = 0.023;
    ey[7] = 0.022;
    ey[8] = 0.037;
    ey[9] = 0.039;

    ey_sys[0] = 0.0032389;
    ey_sys[1] = 0.0032904;
    ey_sys[2] = 0.0033594;
    ey_sys[3] = 0.0043267;
    ey_sys[4] = 0.0065719;
    ey_sys[5] = 0.0066256;
    ey_sys[6] = 0.0065651;
    ey_sys[7] = 0.0067724;
    ey_sys[8] = 0.0075293;
    ey_sys[9] = 0.0093145;
  } else {
    // 20-40%
    x[0] = 0.64;
    x[1] = 1.49;
    x[2] = 2.47;
    x[3] = 3.46;
    x[4] = 4.45;
    x[5] = 5.45;
    x[6] = 6.819;
    x[7] = 8.835;
    x[8] = 10.84;
    x[9] = 14.25;

    y[0] = 0.017;
    y[1] = 0.044;
    y[2] = 0.081;
    y[3] = 0.1;
    y[4] = 0.107;
    y[5] = 0.115;
    y[6] = 0.096;
    y[7] = 0.099;
    y[8] = 0.042;
    y[9] = 0.047;

    ex[0] = 0.5;
    ex[1] = 0.5;
    ex[2] = 0.5;
    ex[3] = 0.5;
    ex[4] = 0.5;
    ex[5] = 0.5;
    ex[6] = 1.;
    ex[7] = 1.;
    ex[8] = 1.;
    ex[9] = 1.5;

    ey[0] = 0.0094;
    ey[1] = 0.0076;
    ey[2] = 0.0075;
    ey[3] = 0.0095;
    ey[4] = 0.0083;
    ey[5] = 0.0097;
    ey[6] = 0.012;
    ey[7] = 0.011;
    ey[8] = 0.031;
    ey[9] = 0.035;

    ey[0] = 0.0031529;
    ey[1] = 0.0038394;
    ey[2] = 0.0040278;
    ey[3] = 0.0039357;
    ey[4] = 0.0033191;
    ey[5] = 0.003598;
    ey[6] = 0.0042368;
    ey[7] = 0.006164;
    ey[8] = 0.0072915;
    ey[9] = 0.010211;
  }
}