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

#include "Framework/Logger.h"

using namespace std;
using namespace RooFit;

void LoadDataRun2(double *&x, double *&y, double *&ex, double *&ey,
                  double *&ey_sys, int flag);
double *CreateBinsFromAxis(TAxis *axis);
void CreateBins(double *axis, double min, double max, int Nbins = 10);
void LoadData(std::string FileName, THnSparse *&hs_V2, TH2F *&hs_R2SPAB,
              TH2F *&hs_R2SPAC, TH2F *&hs_R2SPBC, std::string muonCut,
              std::string dimuonCut);
TH1D *GetMass(double ptmin, double ptmax, double massmin, double massmax,
              double centmin, double centmax, THnSparse *hist_V2);
TH1D *GetV2(double ptmin, double ptmax, double massmin, double massmax,
            double centmin, double centmax, THnSparse *hist_V2, double R2SP);
vector<double> GetStats(int size, double *sample, double *sample_error);
void PlotSNRvsRun2(int size_ptbin, double *pt_bins, int size_run3,
                   double *x_run3, double *snr_run3, int size_run2,
                   double *x_run2, double *snr_run2, TList *ls);
void PlotSystematics(int index, TCanvas *c_sys_yield, TCanvas *c_sys_v2,
                     TH1D *hist_sys_yield, TH1D *hist_sys_v2,
                     double *bins_sys_yield, double *bins_sys_v2,
                     double *chi2_yield, double *chi2_v2, int nbCombo_yield,
                     int nbCombo_v2, vector<double> stats_yield,
                     vector<double> stats_v2, int size_ptbin, double *pt_bins,
                     TList *ls_sys_yield, TList *ls_sys_v2);
void PlotFinalResults(int size_ptbin, double *pt_bins, double *x_v2pt,
                      double *y_v2pt, double *ex_v2pt, double *ey_v2pt,
                      double *eysys_v2pt, double *x_run2, double *y_run2,
                      double *ex_run2, double *ey_run2, double *eysys_run2,
                      double *x_yield, double *y_yield, double *ex_yield,
                      double *ey_yield, double *eysys_yield,
                      double *x_yield_run2, double *y_yield_run2,
                      double *ex_yield_run2, double *ey_yield_run2,
                      double *eysys_yield_run2, TList *ls);

//______________________________________________________________________________
void FlowAnalysis_ScalarProduct(
    int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    std::string FileName = "AnalysisResults.root", double mass_min = 2.3,
    double mass_max = 4.3, double cent_min = 10., double cent_max = 50.,
    double chi2max_mass = 2., double chi2max_v2 = 2., bool sys = false,
    std::string muonCut = "muonLowPt210SigmaPDCA", std::string dimuonCut = "") {

  // Load input data for analysis
  THnSparse *hs_V2;
  TH2F *hs_R2SPAB, *hs_R2SPAC, *hs_R2SPBC;
  LoadData(FileName, hs_V2, hs_R2SPAB, hs_R2SPAC, hs_R2SPBC, muonCut,
           dimuonCut);

  // Get binning information
  TAxis *massAxis = hs_V2->GetAxis(0);
  TAxis *ptAxis = hs_V2->GetAxis(1);
  TAxis *centAxis = hs_V2->GetAxis(3);
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_mass = CreateBinsFromAxis(massAxis);
  double *Bin_pt = CreateBinsFromAxis(ptAxis);
  double *Bin_cent = CreateBinsFromAxis(centAxis);

  // Initialize fitter
  FlowAnalysis_Fitting fitter;
  fitter.init();
  fitter.setChi2MaxMass(chi2max_mass);
  fitter.setChi2MaxV2(chi2max_v2);
  fitter.setCentRange(cent_min, cent_max);

  // Define variables' range for analysis
  // double Bin_pt_mass[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
  double Bin_pt_mass[16] = {0., 0.3, 1., 2.,  3.,  4.,  5.,  6.,
                            7., 8.,  9., 10., 11., 12., 15., 20.};

  // Define the pool for systematics: 36
  // combinationss
  double mass_min_sys[3] = {2.2, 2.3, 2.4};
  double mass_max_sys[3] = {4.2, 4.3, 4.4};
  string sig_enum[5] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev", "VWG"};
  string bkg_v2_enum[2] = {"Pol2", "Chebychev"};
  int sig_mass[3] = {0, 1, 2}; // CB2(MC,data) NA60
  int bkg_mass[2] = {3, 4};    // Chebychev VWG
  int bkg_v2[2] = {0, 1};      // Pol2 and Chebychev

  // Create output file
  TFile f(sys ? Form("FlowAnalysisResults_"
                     "ScalarProduct_%s_%g_%"
                     "g_%dBinPt_withSys.root",
                     muonCut.c_str(), cent_min, cent_max,
                     int(size(Bin_pt_mass)) - 1)
              : Form("FlowAnalysisResults_"
                     "ScalarProduct_%s_%g_%"
                     "g_%dBinPt.root",
                     muonCut.c_str(), cent_min, cent_max,
                     int(size(Bin_pt_mass)) - 1),
          "RECREATE");

  ///////////////////////////////////////////////////
  ///                                             ///
  ///   Analysis for Differential Flow of J/Psi   ///
  ///                                             ///
  ///////////////////////////////////////////////////

  LOG(info) << "Processing analysis for differential flow ...";
  // Resolution factors
  LOG(info) << "Processing evaluation for resolution factors ...";
  TH1D *hist_r2sp =
      new TH1D("R2SP_Cent", "R_{2}{SP} 3-Sub", Nbins_cent, Bin_cent);
  hist_r2sp->GetXaxis()->SetTitle("Centrality FT0 %");
  hist_r2sp->GetYaxis()->SetTitle("R_{2}{SP}");
  hist_r2sp->SetMarkerStyle(20);
  hist_r2sp->SetMarkerSize(1.3);
  hist_r2sp->SetMarkerColor(kRed);
  hist_r2sp->SetLineColor(kRed);
  TH1D *hist_r2sp_AB =
      new TH1D("R2SPAB_Cent", "R_{2}{SP} TPC-FT0A", Nbins_cent, Bin_cent);
  hist_r2sp_AB->GetXaxis()->SetTitle("Centrality FT0 %");
  hist_r2sp_AB->GetYaxis()->SetTitle("R^{AB}_{2}{SP}");
  hist_r2sp_AB->SetMarkerStyle(21);
  hist_r2sp_AB->SetMarkerSize(1.3);
  hist_r2sp_AB->SetMarkerColor(kBlue);
  hist_r2sp_AB->SetLineColor(kBlue);
  TH1D *hist_r2sp_AC =
      new TH1D("R2SPAC_Cent", "R_{2}{SP} TPC-FT0C", Nbins_cent, Bin_cent);
  hist_r2sp_AC->GetXaxis()->SetTitle("Centrality FT0 %");
  hist_r2sp_AC->GetYaxis()->SetTitle("R^{AC}_{2}{SP}");
  hist_r2sp_AC->SetMarkerStyle(22);
  hist_r2sp_AC->SetMarkerSize(1.3);
  hist_r2sp_AC->SetMarkerColor(kOrange);
  hist_r2sp_AC->SetLineColor(kOrange);
  TH1D *hist_r2sp_BC =
      new TH1D("R2SPBC_Cent", "R_{2}{SP} FT0A-FT0C", Nbins_cent, Bin_cent);
  hist_r2sp_BC->GetXaxis()->SetTitle("Centrality FT0 %");
  hist_r2sp_BC->GetYaxis()->SetTitle("R^{BC}_{2}{SP}");
  hist_r2sp_BC->SetMarkerStyle(23);
  hist_r2sp_BC->SetMarkerSize(1.3);
  hist_r2sp_BC->SetMarkerColor(kGreen);
  hist_r2sp_BC->SetLineColor(kGreen);

  for (int i = 0; i < Nbins_cent; i++) {
    // Resolution factor for SP
    TH2F *hs_R2SPAB_cp = dynamic_cast<TH2F *>(hs_R2SPAB->Clone(
        Form("R2SPAB_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
    TH2F *hs_R2SPAC_cp = dynamic_cast<TH2F *>(hs_R2SPAC->Clone(
        Form("R2SPAC_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
    TH2F *hs_R2SPBC_cp = dynamic_cast<TH2F *>(hs_R2SPBC->Clone(
        Form("R2SPBC_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
    hs_R2SPAB_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);
    hs_R2SPAC_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);
    hs_R2SPBC_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);

    double R2SPAB = hs_R2SPAB_cp->GetMean(2);
    double R2SPAC = hs_R2SPAC_cp->GetMean(2);
    double R2SPBC = hs_R2SPBC_cp->GetMean(2);
    double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
    double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
    hist_r2sp->SetBinContent(i + 1, R2SP);
    hist_r2sp_AB->SetBinContent(i + 1, R2SPAB > 0 ? TMath::Sqrt(R2SPAB) : 0.0);
    hist_r2sp_AC->SetBinContent(i + 1, R2SPAC > 0 ? TMath::Sqrt(R2SPAC) : 0.0);
    hist_r2sp_BC->SetBinContent(i + 1, R2SPBC > 0 ? TMath::Sqrt(R2SPBC) : 0.0);

    delete hs_R2SPAB_cp;
    delete hs_R2SPAC_cp;
    delete hs_R2SPBC_cp;
  }
  // Save resolution factor histograms
  auto hs_r2stack = new THStack("R2SPALL", "R_{2}{SP}(3-Sub,AB,AC,BC)");
  hs_r2stack->Add(hist_r2sp);
  hs_r2stack->Add(hist_r2sp_AB);
  hs_r2stack->Add(hist_r2sp_AC);
  hs_r2stack->Add(hist_r2sp_BC);
  TCanvas *c_r2 = new TCanvas("R2All", "R2All");
  c_r2->cd();
  hs_r2stack->Draw("nostack HIST PL");
  TList *l_r2 = new TList();
  l_r2->Add(hist_r2sp);
  l_r2->Add(hist_r2sp_AB);
  l_r2->Add(hist_r2sp_AC);
  l_r2->Add(hist_r2sp_BC);
  l_r2->Add(c_r2);
  f.cd();
  l_r2->Write("ResolutionFactor", TObject::kSingleKey);

  // Evaluation of resolution factors in sub-bins of given centrality range
  double Bin_cent_r2[5] = {10., 20., 30., 40., 50.};
  double R2SP_sub[4];
  for (int i = 0; i < int(size(Bin_cent_r2)) - 1; i++) {
    // Resolution factor for SP
    TH2F *hs_R2SPAB_cp = dynamic_cast<TH2F *>(hs_R2SPAB->Clone(
        Form("R2SPAB_Cent_Copy_%g_%g", Bin_cent_r2[i], Bin_cent_r2[i + 1])));
    TH2F *hs_R2SPAC_cp = dynamic_cast<TH2F *>(hs_R2SPAC->Clone(
        Form("R2SPAC_Cent_Copy_%g_%g", Bin_cent_r2[i], Bin_cent_r2[i + 1])));
    TH2F *hs_R2SPBC_cp = dynamic_cast<TH2F *>(hs_R2SPBC->Clone(
        Form("R2SPBC_Cent_Copy_%g_%g", Bin_cent_r2[i], Bin_cent_r2[i + 1])));

    hs_R2SPAB_cp->GetXaxis()->SetRangeUser(Bin_cent_r2[i], Bin_cent_r2[i + 1]);
    hs_R2SPAC_cp->GetXaxis()->SetRangeUser(Bin_cent_r2[i], Bin_cent_r2[i + 1]);
    hs_R2SPBC_cp->GetXaxis()->SetRangeUser(Bin_cent_r2[i], Bin_cent_r2[i + 1]);

    double R2SPAB = hs_R2SPAB_cp->GetMean(2);
    double R2SPAC = hs_R2SPAC_cp->GetMean(2);
    double R2SPBC = hs_R2SPBC_cp->GetMean(2);
    double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
    double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;

    R2SP_sub[i] = R2SP;

    delete hs_R2SPAB_cp;
    delete hs_R2SPAC_cp;
    delete hs_R2SPBC_cp;
  }

  // Calculation of resolution factor in the given centrality bin
  hs_R2SPAB->GetXaxis()->SetRangeUser(cent_min, cent_max);
  hs_R2SPAC->GetXaxis()->SetRangeUser(cent_min, cent_max);
  hs_R2SPBC->GetXaxis()->SetRangeUser(cent_min, cent_max);
  double R2SPAB = hs_R2SPAB->GetMean(2);
  double R2SPAC = hs_R2SPAC->GetMean(2);
  double R2SPBC = hs_R2SPBC->GetMean(2);
  double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
  double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
  if (R2SP == 0.0) {
    LOG(fatal) << "Inconsistent value of resolution factor!";
  }

  // Create histogram for pt-differential v2
  TList *l_results = new TList();
  double *x_yield = new double[int(size(Bin_pt_mass)) - 1];
  double *y_yield = new double[int(size(Bin_pt_mass)) - 1];
  double *ex_yield = new double[int(size(Bin_pt_mass)) - 1];
  double *ey_yield = new double[int(size(Bin_pt_mass)) - 1];
  double *eysys_yield = new double[int(size(Bin_pt_mass)) - 1];
  double *SNR = new double[int(size(Bin_pt_mass)) - 1];
  double *x_v2pt = new double[int(size(Bin_pt_mass)) - 1];
  double *y_v2pt = new double[int(size(Bin_pt_mass)) - 1];
  double *ex_v2pt = new double[int(size(Bin_pt_mass)) - 1];
  double *ey_v2pt = new double[int(size(Bin_pt_mass)) - 1];
  double *eysys_v2pt = new double[int(size(Bin_pt_mass)) - 1];

  double x_yield_run2[15] = {0.15, 0.65, 1.5, 2.5,  3.5,  4.5,  5.5, 6.5,
                             7.5,  8.5,  9.5, 10.5, 11.5, 13.5, 17.5};
  double ex_yield_run2[15] = {0.15, 0.35, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                              0.5,  0.5,  0.5, 0.5, 0.5, 1.5, 2.5};
  double y_yield_run2[15] = {
      6584 / 0.3, 38709 / 0.7, 66303, 49791, 30467, 17566,    9805,    5789,
      3203,       1781,        1023,  638,   457,   509 / 3., 249 / 5.};
  double ey_yield_run2[15] = {401 / 0.3, 984 / 0.7, 1241, 972,     491,
                              306,       206,       140,  97,      70,
                              54,        41,        36,   37 / 3., 24 / 5.};
  double eysys_yield_run2[15] = {206 / 0.3, 1052 / 0.7, 1748, 1272,    892,
                                 447,       189,        112,  67,      49,
                                 21,        13,         14,   23 / 3., 6 / 5.};
  double SNR_run2[15] = {0.18, 0.14, 0.18, 0.3,  0.44, 0.67, 0.87, 1.19,
                         1.53, 1.76, 2.01, 2.10, 2.59, 1.92, 2.81};

  // Load Run2 data for comparaison
  double *x_run2, *y_run2, *ex_run2, *ey_run2, *eysys_run2;
  LoadDataRun2(x_run2, y_run2, ex_run2, ey_run2, eysys_run2, flag_run2);

  // Initialize arrays for each trial of systematic study
  const int nbCombo_v2 = int(size(sig_mass)) * int(size(bkg_mass)) *
                         int(size(bkg_v2)) * int(size(mass_max_sys));
  const int nbCombo_yield =
      int(size(sig_mass)) * int(size(bkg_mass)) * int(size(mass_max_sys));

  vector<double *> y_sys_yield, ey_sys_yield, y_sys_v2, ey_sys_v2, chi2_yield,
      chi2_v2;
  double *bins_sys_v2 = new double[nbCombo_v2 + 1];
  double *bins_sys_yield = new double[nbCombo_yield + 1];
  for (int i = 0; i < nbCombo_v2 + 1; i++) {
    bins_sys_v2[i] = i;
  }
  for (int i = 0; i < nbCombo_yield + 1; i++) {
    bins_sys_yield[i] = i;
  }

  vector<TH1D *> hist_sys_yield, hist_sys_v2;
  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    y_sys_yield.emplace_back(new double[nbCombo_yield]);
    ey_sys_yield.emplace_back(new double[nbCombo_yield]);

    y_sys_v2.emplace_back(new double[nbCombo_v2]);
    ey_sys_v2.emplace_back(new double[nbCombo_v2]);

    chi2_yield.emplace_back(new double[nbCombo_yield]);
    chi2_v2.emplace_back(new double[nbCombo_v2]);

    hist_sys_yield.emplace_back(new TH1D(
        Form("hist_sys_yield_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        Form("hist_sys_yield_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        nbCombo_yield, bins_sys_yield));
    hist_sys_v2.emplace_back(
        new TH1D(Form("hist_sys_v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                 Form("hist_sys_v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                 nbCombo_v2, bins_sys_v2));
  }

  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    // Get mass and v2 to fit
    TH1D *hs_mass_proj = GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                 mass_max, cent_min, cent_max, hs_V2);
    TH1D *hs_v2sp = GetV2(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                          mass_max, cent_min, cent_max, hs_V2, R2SP);

    /// Do fitting
    // Configuration for fitting
    fitter.setModel(flag_sig, flag_bkg);
    fitter.setModelV2(flag_v2);
    fitter.setMassRange(mass_min, mass_max);
    fitter.setPtRange(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
    fitter.setOrder(2);
    fitter.setMode(0);

    // Fit invariant mass + v2
    TList *l_diff_fit = new TList();
    vector<double> results_v2 =
        fitter.runFitting(hs_mass_proj, hs_v2sp, l_diff_fit);

    SNR[i] = results_v2[4];
    x_yield[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    ex_yield[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;
    x_v2pt[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    ex_v2pt[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;

    // Save results
    TList *l_diff_hist = new TList();
    l_diff_hist->Add(hs_v2sp);
    f.cd();
    l_diff_hist->Write(
        Form("DifferentialFlow_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    l_diff_fit->Write(
        Form("DifferentialFlow_Fit_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);

    delete hs_mass_proj;
    delete hs_v2sp;
    delete l_diff_hist;
    delete l_diff_fit;

    // Run fittings for systematics
    if (sys) {
      LOG(info) << "Processing systematics from fitting ...";
      int index_sys_v2 = 0;
      int index_sys_yield = 0;
      for (int i1 = 0; i1 < int(size(mass_max_sys)); i1++) {
        for (int i2 = 0; i2 < int(size(sig_mass)); i2++) {
          for (int i3 = 0; i3 < int(size(bkg_mass)); i3++) {
            TString combo_yield =
                Form("{%s;%s}[%g-%g]", sig_enum[sig_mass[i2]].c_str(),
                     sig_enum[bkg_mass[i3]].c_str(), mass_min_sys[i1],
                     mass_max_sys[i1]);
            for (int i4 = 0; i4 < int(size(bkg_v2)); i4++) {
              // Get mass and v2 to fit
              TH1D *hist_mass_proj_sys =
                  GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                          mass_max_sys[i1], cent_min, cent_max, hs_V2);

              TH1D *hs_v2sp_sys =
                  GetV2(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                        mass_max_sys[i1], cent_min, cent_max, hs_V2, R2SP);

              TString combo_v2 =
                  Form("{%s;%s}{%s}[%g-%g]", sig_enum[sig_mass[i2]].c_str(),
                       sig_enum[bkg_mass[i3]].c_str(), bkg_v2_enum[i4].c_str(),
                       mass_min_sys[i1], mass_max_sys[i1]);
              LOG(info) << Form("Processing combination: {%s,%s}+{%s}+[%g-%g]",
                                sig_enum[sig_mass[i2]].c_str(),
                                sig_enum[bkg_mass[i3]].c_str(),
                                bkg_v2_enum[i4].c_str(), mass_min_sys[i1],
                                mass_max_sys[i1]);

              // Do evaluation for v2
              fitter.setOrder(2);
              fitter.setMassRange(mass_min_sys[i1], mass_max_sys[i1]);
              fitter.setPtRange(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
              fitter.setModel(sig_mass[i2], bkg_mass[i3]);
              fitter.setModelV2(bkg_v2[i4]);
              fitter.setMode(1);
              TList *l_diff_sys = new TList();
              LOG(info) << "Processing fitting(systematic) for v2{SP} ...";
              vector<double> results_sys_v2 = fitter.runFitting(
                  hist_mass_proj_sys, hs_v2sp_sys, l_diff_sys);

              // Fill pT-differential v2 and jpsi
              // yields
              y_sys_yield[i][index_sys_yield] = results_sys_v2[2];
              y_sys_v2[i][index_sys_v2] = results_sys_v2[0];

              ey_sys_yield[i][index_sys_yield] = results_sys_v2[3];
              ey_sys_v2[i][index_sys_v2] = results_sys_v2[1];

              chi2_yield[i][index_sys_yield] = results_sys_v2[5];
              chi2_v2[i][index_sys_v2] = results_sys_v2[6];

              hist_sys_yield[i]->GetXaxis()->SetBinLabel(index_sys_yield + 1,
                                                         combo_yield);
              hist_sys_v2[i]->GetXaxis()->SetBinLabel(index_sys_v2 + 1,
                                                      combo_v2);
              hist_sys_yield[i]->SetBinContent(index_sys_yield + 1,
                                               y_sys_yield[i][index_sys_yield]);
              hist_sys_v2[i]->SetBinContent(index_sys_v2 + 1,
                                            y_sys_v2[i][index_sys_v2]);
              hist_sys_yield[i]->SetBinError(index_sys_yield + 1,
                                             ey_sys_yield[i][index_sys_yield]);
              hist_sys_v2[i]->SetBinError(index_sys_v2 + 1,
                                          ey_sys_v2[i][index_sys_v2]);

              f.cd();
              l_diff_sys->Write(Form("FitSys_%g_%g_%s", Bin_pt_mass[i],
                                     Bin_pt_mass[i + 1], combo_v2.Data()),
                                TObject::kSingleKey);

              delete l_diff_sys;

              delete hist_mass_proj_sys;
              delete hs_v2sp_sys;

              index_sys_v2++;
            }
            index_sys_yield++;
          }
        }
      }
    }
  }

  // Print resolution factors
  LOG(info) << "Resolution factors: ";
  for (int i = 0; i < int(size(Bin_cent_r2)) - 1; i++) {
    LOG(info) << Form("R2SP [%g-%g%%] = %g", Bin_cent_r2[i], Bin_cent_r2[i + 1],
                      R2SP_sub[i]);
  }
  LOG(info) << Form("R2SP [10-50%%] = %g", R2SP);
  for (int i = 0; i < int(size(Bin_cent_r2)) - 1; i++) {
    LOG(info) << Form("1/R2SP [%g-%g%%] = %g", Bin_cent_r2[i],
                      Bin_cent_r2[i + 1], 1. / R2SP_sub[i]);
  }
  LOG(info) << Form("1/R2SP [10-50%%] = %g", 1. / R2SP);

  // Saving plot for J/psi yields SNR as function of pT
  // compared with Run2
  PlotSNRvsRun2(int(size(Bin_pt_mass)) - 1, Bin_pt_mass,
                int(size(Bin_pt_mass)) - 1, x_yield, SNR,
                int(size(x_yield_run2)), x_yield_run2, SNR_run2, l_results);

  // Save plots for systematics
  if (sys) {
    TList *l_results_sys_yield = new TList();
    TList *l_results_sys_v2 = new TList();
    vector<TCanvas *> c_sys_yield;
    vector<TCanvas *> c_sys_v2;
    for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
      c_sys_yield.emplace_back(new TCanvas(
          Form("Sys_yield_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
          Form("Sys_yield_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1])));
      c_sys_v2.emplace_back(new TCanvas(
          Form("Sys_v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
          Form("Sys_v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1])));
    }
    for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
      vector<double> stats_yield =
          GetStats(nbCombo_yield, y_sys_yield[i], ey_sys_yield[i]);
      vector<double> stats_v2 = GetStats(nbCombo_v2, y_sys_v2[i], ey_sys_v2[i]);

      // Fill pT-differential v2 and jpsi yields
      y_v2pt[i] = stats_v2[0];
      ey_v2pt[i] = stats_v2[1];
      eysys_v2pt[i] = stats_v2[2];

      y_yield[i] = stats_yield[0] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
      ey_yield[i] = stats_yield[1] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
      eysys_yield[i] = stats_yield[2] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);

      // Saving results for systematics
      PlotSystematics(i, c_sys_yield[i], c_sys_v2[i], hist_sys_yield[i],
                      hist_sys_v2[i], bins_sys_yield, bins_sys_v2,
                      chi2_yield[i], chi2_v2[i], nbCombo_yield, nbCombo_v2,
                      stats_yield, stats_v2, int(size(Bin_pt_mass)) - 1,
                      Bin_pt_mass, l_results_sys_yield, l_results_sys_v2);
    }
    f.cd();
    l_results_sys_yield->Write("FitYieldSystematics", TObject::kSingleKey);
    l_results_sys_v2->Write("FitV2Systematics", TObject::kSingleKey);

    // Saving final results
    PlotFinalResults(int(size(Bin_pt_mass)) - 1, Bin_pt_mass, x_v2pt, y_v2pt,
                     ex_v2pt, ey_v2pt, eysys_v2pt, x_run2, y_run2, ex_run2,
                     ey_run2, eysys_run2, x_yield, y_yield, ex_yield, ey_yield,
                     eysys_yield, x_yield_run2, y_yield_run2, ex_yield_run2,
                     ey_yield_run2, eysys_yield_run2, l_results);
  }
  f.cd();
  l_results->Write("FinalResults", TObject::kSingleKey);
  f.Close();
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
              double centmin, double centmax, THnSparse *hist_V2) {

  // Copy original profiles for projections
  THnSparse *hist_V2_cp = dynamic_cast<THnSparse *>(
      hist_V2->Clone(Form("Mass_Pt_centrFT0C_V2_Copy_%g_%g", ptmin, ptmax)));

  // Set axes' ranges for mass-differential study
  hist_V2_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hist_V2_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hist_V2_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  // Get mass
  TH1D *hist_proj = hist_V2_cp->Projection(0);
  TH1D *hist_mass_proj = dynamic_cast<TH1D *>(
      hist_proj->Clone(Form("Proj_%s", hist_proj->GetName())));

  delete hist_V2_cp;
  delete hist_proj;

  return hist_mass_proj;
}

//______________________________________________________________________________
TH1D *GetV2(double ptmin, double ptmax, double massmin, double massmax,
            double centmin, double centmax, THnSparse *hist_V2, double R2SP) {

  // Copy original profiles for projections
  THnSparse *hist_V2_cp = dynamic_cast<THnSparse *>(
      hist_V2->Clone(Form("Mass_Pt_centrFT0C_V2_Copy_%g_%g", ptmin, ptmax)));

  // Set axes' ranges for mass-differential study
  hist_V2_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hist_V2_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hist_V2_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  // Get v2
  TH2D *hs_v2_sp_proj = hist_V2_cp->Projection(4, 0);

  // Define histograms
  double *Bin_mass_new = CreateBinsFromAxis(hs_v2_sp_proj->GetXaxis());
  int NBins_mass_new = hs_v2_sp_proj->GetXaxis()->GetNbins();
  TH1D *hist_v2sp = new TH1D(Form("v2sp_%g_%g", ptmin, ptmax),
                             Form("v^{#mu#mu}_{2}{SP}_%g_%g", ptmin, ptmax),
                             NBins_mass_new, Bin_mass_new);
  hist_v2sp->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_v2sp->GetYaxis()->SetTitle("v^{#mu#mu}_{2}{SP}");

  // Evaluation of differential flow as function of invariant mass
  for (int i = 0; i < NBins_mass_new; i++) {
    TH2D *hs_v2_sp_proj_cp =
        dynamic_cast<TH2D *>(hs_v2_sp_proj->Clone("Mass_V2SP_Copy"));
    hs_v2_sp_proj_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                               Bin_mass_new[i + 1]);
    double v2sp = hs_v2_sp_proj_cp->GetMean(2);
    double v2spe = hs_v2_sp_proj_cp->GetMean(12);
    hist_v2sp->SetBinContent(i + 1, v2sp / R2SP);
    hist_v2sp->SetBinError(i + 1, v2spe / R2SP);

    delete hs_v2_sp_proj_cp;
  }

  delete hist_V2_cp;
  delete hs_v2_sp_proj;

  return hist_v2sp;
}

//______________________________________________________________________________
vector<double> GetStats(int size, double *sample, double *sample_error) {
  vector<double> results;
  double mean = 0.;
  double mean_error = 0.;
  for (int i = 0; i < size; i++) {
    mean += sample[i] / size;
    mean_error += sample_error[i] / size;
  }
  results.emplace_back(mean);
  results.emplace_back(mean_error);

  double sum2 = 0;
  for (int i = 0; i < size; i++) {
    sum2 += pow(sample[i] - mean, 2.) / size;
  }
  double rms = pow(sum2, 0.5);
  results.emplace_back(rms);
  return results;
}

//______________________________________________________________________________
void PlotSystematics(int index, TCanvas *c_sys_yield, TCanvas *c_sys_v2,
                     TH1D *hist_sys_yield, TH1D *hist_sys_v2,
                     double *bins_sys_yield, double *bins_sys_v2,
                     double *chi2_yield, double *chi2_v2, int nbCombo_yield,
                     int nbCombo_v2, vector<double> stats_yield,
                     vector<double> stats_v2, int size_ptbin, double *pt_bins,
                     TList *ls_sys_yield, TList *ls_sys_v2) {

  // Plotting for yields systematics
  c_sys_yield->cd();
  c_sys_yield->SetBottomMargin(0);
  c_sys_yield->SetCanvasSize(1000, 400);
  TPad *pad_sys_yield =
      new TPad(Form("pad_sys_yield_%d", index), Form("pad_sys_yield_%d", index),
               0, 0.3, 1, 1.0);
  pad_sys_yield->SetBottomMargin(0);
  pad_sys_yield->Draw();
  pad_sys_yield->cd();
  hist_sys_yield->SetMarkerStyle(20);
  hist_sys_yield->SetMarkerSize(1.);
  hist_sys_yield->SetMarkerColor(kBlack);
  hist_sys_yield->SetLineColor(kBlack);
  hist_sys_yield->SetLineWidth(2);
  hist_sys_yield->SetFillStyle(0);
  hist_sys_yield->SetStats(0);
  hist_sys_yield->SetTitle("");
  hist_sys_yield->GetYaxis()->SetTitle("N_{J/#psi}");
  hist_sys_yield->GetYaxis()->SetTitleSize(0.05);
  hist_sys_yield->GetYaxis()->SetTitleOffset(0.5);
  hist_sys_yield->GetXaxis()->SetLabelOffset(999);
  hist_sys_yield->GetXaxis()->SetLabelSize(0);
  hist_sys_yield->Draw("HIST EP");
  TF1 *lyield_mean = new TF1("meanyield", "[0]", bins_sys_yield[0],
                             bins_sys_yield[nbCombo_yield]);
  lyield_mean->SetParameter(0, stats_yield[0]);
  lyield_mean->SetLineColor(kBlue);
  lyield_mean->SetLineWidth(3);
  lyield_mean->SetLineStyle(1);
  lyield_mean->Draw("same");
  TF1 *lyield_meanerrorp =
      new TF1("meanerrorpyield", "[0]+[1]", bins_sys_yield[0],
              bins_sys_yield[nbCombo_yield]);
  lyield_meanerrorp->SetParameter(0, stats_yield[0]);
  lyield_meanerrorp->SetParameter(1, stats_yield[1]);
  lyield_meanerrorp->SetLineColor(kBlue);
  lyield_meanerrorp->SetLineWidth(3);
  lyield_meanerrorp->SetLineStyle(7);
  lyield_meanerrorp->Draw("same");
  TF1 *lyield_meanerrorm =
      new TF1("meanerrormyield", "[0]-[1]", bins_sys_yield[0],
              bins_sys_yield[nbCombo_yield]);
  lyield_meanerrorm->SetParameter(0, stats_yield[0]);
  lyield_meanerrorm->SetParameter(1, stats_yield[1]);
  lyield_meanerrorm->SetLineColor(kBlue);
  lyield_meanerrorm->SetLineWidth(3);
  lyield_meanerrorm->SetLineStyle(7);
  lyield_meanerrorm->Draw("same");
  TF1 *lyield_rmsp = new TF1("rmspyield", "[0]+[1]", bins_sys_yield[0],
                             bins_sys_yield[nbCombo_yield]);
  lyield_rmsp->SetParameter(0, stats_yield[0]);
  lyield_rmsp->SetParameter(1, stats_yield[2]);
  lyield_rmsp->SetLineColor(kBlue);
  lyield_rmsp->SetLineWidth(3);
  lyield_rmsp->SetLineStyle(9);
  lyield_rmsp->Draw("same");
  TF1 *lyield_rmsm = new TF1("rmsmyield", "[0]-[1]", bins_sys_yield[0],
                             bins_sys_yield[nbCombo_yield]);
  lyield_rmsm->SetParameter(0, stats_yield[0]);
  lyield_rmsm->SetParameter(1, stats_yield[2]);
  lyield_rmsm->SetLineColor(kBlue);
  lyield_rmsm->SetLineWidth(3);
  lyield_rmsm->SetLineStyle(9);
  lyield_rmsm->Draw("same");
  TLatex *text_sys_yield = new TLatex();
  text_sys_yield->SetTextSize(0.05);
  text_sys_yield->SetTextFont(42);
  text_sys_yield->SetTextColor(kBlue);
  text_sys_yield->DrawLatexNDC(
      .12, .85,
      Form("N_{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
           "(sys)",
           pt_bins[index], pt_bins[index + 1], stats_yield[0], stats_yield[1],
           100. * stats_yield[1] / stats_yield[0], stats_yield[2],
           100. * stats_yield[2] / stats_yield[0]));
  pad_sys_yield->ModifiedUpdate();
  c_sys_yield->cd();
  TPad *pad_sys_yield_chi =
      new TPad(Form("pad_sys_yield_chi2_%d", index),
               Form("pad_sys_yield_chi2_%d", index), 0, 0., 1, 0.3);
  pad_sys_yield_chi->SetTopMargin(0);
  pad_sys_yield_chi->SetBottomMargin(0.6);
  pad_sys_yield_chi->Draw();
  pad_sys_yield_chi->cd();
  TH1D *hist_chi2_yield =
      (TH1D *)hist_sys_yield->Clone(Form("hist_yield_chi2_%d", index));
  for (int j = 0; j < hist_chi2_yield->GetNbinsX(); j++) {
    hist_chi2_yield->SetBinContent(j + 1, chi2_yield[j]);
    hist_chi2_yield->SetBinError(j + 1, 0.);
  }
  hist_chi2_yield->SetTitle("");
  hist_chi2_yield->GetYaxis()->SetLabelSize(0.05);
  hist_chi2_yield->GetYaxis()->SetRangeUser(0.5, 3.);
  hist_chi2_yield->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hist_chi2_yield->GetYaxis()->SetTitleSize(0.08);
  hist_chi2_yield->GetYaxis()->SetTitleOffset(0.25);
  hist_chi2_yield->GetXaxis()->SetLabelSize(0.14);
  hist_chi2_yield->GetXaxis()->SetLabelOffset(0.02);
  hist_chi2_yield->Draw("HIST P");
  TF1 *lchi2_yield1 = new TF1("lchi2_yield1", "[0]", bins_sys_yield[0],
                              bins_sys_yield[nbCombo_yield]);
  lchi2_yield1->SetParameter(0, 1.0);
  lchi2_yield1->SetLineColor(kBlue);
  lchi2_yield1->SetLineWidth(2);
  lchi2_yield1->SetLineStyle(9);
  lchi2_yield1->Draw("same");
  TF1 *lchi2_yield2 = new TF1("lchi2_yield2", "[0]", bins_sys_yield[0],
                              bins_sys_yield[nbCombo_yield]);
  lchi2_yield2->SetParameter(0, 2.0);
  lchi2_yield2->SetLineColor(kRed);
  lchi2_yield2->SetLineWidth(2);
  lchi2_yield2->SetLineStyle(1);
  lchi2_yield2->Draw("same");

  // Plotting for v2 systematics
  c_sys_v2->cd();
  c_sys_v2->SetBottomMargin(0);
  c_sys_v2->SetCanvasSize(1200, 400);
  TPad *pad_sys_v2 = new TPad(Form("pad_sys_v2_%d", index),
                              Form("pad_sys_v2_%d", index), 0, 0.3, 1, 1.0);
  pad_sys_v2->SetBottomMargin(0);
  pad_sys_v2->Draw();
  pad_sys_v2->cd();
  hist_sys_v2->SetMarkerStyle(20);
  hist_sys_v2->SetMarkerSize(1.);
  hist_sys_v2->SetMarkerColor(kBlack);
  hist_sys_v2->SetLineColor(kBlack);
  hist_sys_v2->SetLineWidth(2);
  hist_sys_v2->SetFillStyle(0);
  hist_sys_v2->SetStats(0);
  hist_sys_v2->SetTitle("");
  hist_sys_v2->GetYaxis()->SetTitle("#it{v}_{2}{SP}");
  hist_sys_v2->GetYaxis()->SetTitleSize(0.05);
  hist_sys_v2->GetYaxis()->SetTitleOffset(0.5);
  hist_sys_v2->GetXaxis()->SetLabelOffset(999);
  hist_sys_v2->GetXaxis()->SetLabelSize(0);
  hist_sys_v2->Draw("HIST EP");
  TF1 *lv2_mean =
      new TF1("meanv2", "[0]", bins_sys_v2[0], bins_sys_v2[nbCombo_v2]);
  lv2_mean->SetParameter(0, stats_v2[0]);
  lv2_mean->SetLineColor(kBlue);
  lv2_mean->SetLineWidth(3);
  lv2_mean->SetLineStyle(1);
  lv2_mean->Draw("same");
  TF1 *lv2_meanerrorp = new TF1("meanerrorpv2", "[0]+[1]", bins_sys_v2[0],
                                bins_sys_v2[nbCombo_v2]);
  lv2_meanerrorp->SetParameter(0, stats_v2[0]);
  lv2_meanerrorp->SetParameter(1, stats_v2[1]);
  lv2_meanerrorp->SetLineColor(kBlue);
  lv2_meanerrorp->SetLineWidth(3);
  lv2_meanerrorp->SetLineStyle(7);
  lv2_meanerrorp->Draw("same");
  TF1 *lv2_meanerrorm = new TF1("meanerrormv2", "[0]-[1]", bins_sys_v2[0],
                                bins_sys_v2[nbCombo_v2]);
  lv2_meanerrorm->SetParameter(0, stats_v2[0]);
  lv2_meanerrorm->SetParameter(1, stats_v2[1]);
  lv2_meanerrorm->SetLineColor(kBlue);
  lv2_meanerrorm->SetLineWidth(3);
  lv2_meanerrorm->SetLineStyle(7);
  lv2_meanerrorm->Draw("same");
  TF1 *lv2_rmsp =
      new TF1("rmspv2", "[0]+[1]", bins_sys_v2[0], bins_sys_v2[nbCombo_v2]);
  lv2_rmsp->SetParameter(0, stats_v2[0]);
  lv2_rmsp->SetParameter(1, stats_v2[2]);
  lv2_rmsp->SetLineColor(kBlue);
  lv2_rmsp->SetLineWidth(3);
  lv2_rmsp->SetLineStyle(9);
  lv2_rmsp->Draw("same");
  TF1 *lv2_rmsm =
      new TF1("rmsmv2", "[0]-[1]", bins_sys_v2[0], bins_sys_v2[nbCombo_v2]);
  lv2_rmsm->SetParameter(0, stats_v2[0]);
  lv2_rmsm->SetParameter(1, stats_v2[2]);
  lv2_rmsm->SetLineColor(kBlue);
  lv2_rmsm->SetLineWidth(3);
  lv2_rmsm->SetLineStyle(9);
  lv2_rmsm->Draw("same");
  TLatex *text_sys_v2 = new TLatex();
  text_sys_v2->SetTextSize(0.05);
  text_sys_v2->SetTextFont(42);
  text_sys_v2->SetTextColor(kBlue);
  text_sys_v2->DrawLatexNDC(
      .12, .85,
      Form("N_{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
           "(sys)",
           pt_bins[index], pt_bins[index + 1], stats_v2[0], stats_v2[1],
           100. * stats_v2[1] / stats_v2[0], stats_v2[2],
           100. * stats_v2[2] / stats_v2[0]));
  pad_sys_v2->ModifiedUpdate();
  c_sys_v2->cd();
  TPad *pad_sys_v2_chi =
      new TPad(Form("pad_sys_v2_chi2_%d", index),
               Form("pad_sys_v2_chi2_%d", index), 0, 0., 1, 0.3);
  pad_sys_v2_chi->SetTopMargin(0);
  pad_sys_v2_chi->SetBottomMargin(0.6);
  pad_sys_v2_chi->Draw();
  pad_sys_v2_chi->cd();
  TH1D *hist_chi2_v2 =
      (TH1D *)hist_sys_v2->Clone(Form("hist_v2_chi2_%d", index));
  for (int j = 0; j < hist_chi2_v2->GetNbinsX(); j++) {
    hist_chi2_v2->SetBinContent(j + 1, chi2_v2[j]);
    hist_chi2_v2->SetBinError(j + 1, 0.);
  }
  hist_chi2_v2->SetTitle("");
  hist_chi2_v2->GetYaxis()->SetLabelSize(0.05);
  hist_chi2_v2->GetYaxis()->SetRangeUser(0.5, 3.);
  hist_chi2_v2->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hist_chi2_v2->GetYaxis()->SetTitleSize(0.08);
  hist_chi2_v2->GetYaxis()->SetTitleOffset(0.25);
  hist_chi2_v2->GetXaxis()->SetLabelSize(0.13);
  hist_chi2_v2->GetXaxis()->SetLabelOffset(0.02);
  hist_chi2_v2->Draw("HIST P");
  TF1 *lchi2_v21 =
      new TF1("lchi2_v21", "[0]", bins_sys_v2[0], bins_sys_v2[nbCombo_v2]);
  lchi2_v21->SetParameter(0, 1.0);
  lchi2_v21->SetLineColor(kBlue);
  lchi2_v21->SetLineWidth(2);
  lchi2_v21->SetLineStyle(9);
  lchi2_v21->Draw("same");
  TF1 *lchi2_v22 =
      new TF1("lchi2_v22", "[0]", bins_sys_v2[0], bins_sys_v2[nbCombo_v2]);
  lchi2_v22->SetParameter(0, 2.0);
  lchi2_v22->SetLineColor(kRed);
  lchi2_v22->SetLineWidth(2);
  lchi2_v22->SetLineStyle(1);
  lchi2_v22->Draw("same");

  ls_sys_v2->Add(c_sys_v2);
  ls_sys_yield->Add(c_sys_yield);
}

//______________________________________________________________________________
void PlotSNRvsRun2(int size_ptbin, double *pt_bins, int size_run3,
                   double *x_run3, double *snr_run3, int size_run2,
                   double *x_run2, double *snr_run2, TList *ls) {

  // Saving plot for J/psi yields SNR as function of pT
  // compared with Run2
  TGraph *graph_SNR_Run3 = new TGraph(size_run3, x_run3, snr_run3);
  graph_SNR_Run3->SetNameTitle("graph_snr_run3", "Run3");
  TGraph *graph_SNR_Run2 = new TGraph(size_run2, x_run2, snr_run2);
  graph_SNR_Run2->SetNameTitle("graph_snr_run2", "Run2");
  TCanvas *c_SNR = new TCanvas("jpsi_SNR_pT", "jpsi_SNR_pT");
  c_SNR->cd();
  if (!(size_run2 == size_run3)) {
    auto mg_SNR = new TMultiGraph();
    graph_SNR_Run3->SetMarkerStyle(20);
    graph_SNR_Run3->SetMarkerSize(1.);
    graph_SNR_Run3->SetMarkerColor(kBlue);
    graph_SNR_Run3->SetLineColor(kBlue);
    graph_SNR_Run3->SetLineWidth(2);
    graph_SNR_Run3->SetFillStyle(0);
    graph_SNR_Run2->SetMarkerStyle(20);
    graph_SNR_Run2->SetMarkerSize(1.);
    graph_SNR_Run2->SetMarkerColor(kRed);
    graph_SNR_Run2->SetLineColor(kRed);
    graph_SNR_Run2->SetLineWidth(2);
    graph_SNR_Run2->SetFillStyle(0);
    mg_SNR->Add(graph_SNR_Run3);
    mg_SNR->Add(graph_SNR_Run2);
    mg_SNR->Draw("A P Z ");
  } else {
    c_SNR->SetBottomMargin(0);
    TPad *pad_snr = new TPad("pad_snr", "pad_snr", 0, 0.3, 1, 1.0);
    pad_snr->SetBottomMargin(0);
    pad_snr->Draw();
    pad_snr->cd();
    auto mg_SNR = new TMultiGraph();
    graph_SNR_Run3->SetMarkerStyle(20);
    graph_SNR_Run3->SetMarkerSize(1.);
    graph_SNR_Run3->SetMarkerColor(kBlue);
    graph_SNR_Run3->SetLineColor(kBlue);
    graph_SNR_Run3->SetLineWidth(2);
    graph_SNR_Run3->SetFillStyle(0);
    graph_SNR_Run2->SetMarkerStyle(20);
    graph_SNR_Run2->SetMarkerSize(1.);
    graph_SNR_Run2->SetMarkerColor(kRed);
    graph_SNR_Run2->SetLineColor(kRed);
    graph_SNR_Run2->SetLineWidth(2);
    graph_SNR_Run2->SetFillStyle(0);
    mg_SNR->Add(graph_SNR_Run3);
    mg_SNR->Add(graph_SNR_Run2);
    mg_SNR->GetXaxis()->SetLimits(0., 20.);
    mg_SNR->GetYaxis()->SetTitle("(S/B)_{3#sigma}");
    mg_SNR->Draw("A P Z ");
    pad_snr->BuildLegend();
    pad_snr->ModifiedUpdate();

    c_SNR->cd();
    TPad *pad_snr_ratio =
        new TPad("pad_snr_ratio", "pad_snr_ratio", 0, 0., 1, 0.3);
    pad_snr_ratio->SetTopMargin(0);
    pad_snr_ratio->SetBottomMargin(0.22);
    pad_snr_ratio->Draw();
    pad_snr_ratio->cd();
    TH1D *hs_snr_ratio = new TH1D("hist_snr_ratio", "", size_ptbin, pt_bins);
    for (int i = 0; i < size_ptbin; i++) {
      hs_snr_ratio->SetBinContent(i + 1, snr_run2[i] / snr_run3[i]);
    }
    hs_snr_ratio->SetStats(0);
    hs_snr_ratio->SetTitle("");
    hs_snr_ratio->SetMarkerStyle(20);
    hs_snr_ratio->SetMarkerSize(0.8);
    hs_snr_ratio->GetYaxis()->SetTitleSize(0.1);
    hs_snr_ratio->GetYaxis()->SetTitle("ratio");
    hs_snr_ratio->GetYaxis()->SetTitleOffset(0.25);
    hs_snr_ratio->GetYaxis()->SetLabelSize(0.1);
    hs_snr_ratio->GetYaxis()->SetRangeUser(0., 10.);
    hs_snr_ratio->GetXaxis()->SetLabelSize(0.1);
    hs_snr_ratio->GetXaxis()->SetLabelOffset();
    hs_snr_ratio->GetXaxis()->SetTitleSize(0.1);
    hs_snr_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hs_snr_ratio->Draw("HIST P");
    TF1 *lratio = new TF1("lratio", "[0]", pt_bins[0], pt_bins[size_ptbin]);
    lratio->SetParameter(0, 1.);
    lratio->SetLineColor(kBlue);
    lratio->SetLineWidth(3);
    lratio->SetLineStyle(1);
    lratio->Draw("same");
    pad_snr_ratio->ModifiedUpdate();
  }
  ls->Add(c_SNR);
}

//______________________________________________________________________________
void PlotFinalResults(int size_ptbin, double *pt_bins, double *x_v2pt,
                      double *y_v2pt, double *ex_v2pt, double *ey_v2pt,
                      double *eysys_v2pt, double *x_run2, double *y_run2,
                      double *ex_run2, double *ey_run2, double *eysys_run2,
                      double *x_yield, double *y_yield, double *ex_yield,
                      double *ey_yield, double *eysys_yield,
                      double *x_yield_run2, double *y_yield_run2,
                      double *ex_yield_run2, double *ey_yield_run2,
                      double *eysys_yield_run2, TList *ls) {

  // Compare with Run2 data: 5.02 TeV
  TGraphMultiErrors *graph_v2pt =
      new TGraphMultiErrors("graph_v2_pt", "", size_ptbin, x_v2pt, y_v2pt,
                            ex_v2pt, ex_v2pt, ey_v2pt, ey_v2pt);
  graph_v2pt->SetTitle("#sqrt{#it{s}_{NN}} = 5.36 TeV, 10-50%");
  graph_v2pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");
  graph_v2pt->AddYError(size_ptbin, eysys_v2pt, eysys_v2pt);

  TGraphMultiErrors *graph_v2pt_run2 =
      new TGraphMultiErrors("graph_v2_pt_run2", "", 10, x_run2, y_run2, ex_run2,
                            ex_run2, ey_run2, ey_run2);
  graph_v2pt_run2->AddYError(10, eysys_run2, eysys_run2);
  graph_v2pt_run2->SetTitle("#sqrt{#it{s}_{NN}} = 5.02 TeV, 10-30%");
  graph_v2pt_run2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt_run2->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

  TGraphMultiErrors *graph_yield =
      new TGraphMultiErrors("graph_yields_pt", "Run3", size_ptbin, x_yield,
                            y_yield, ex_yield, ex_yield, ey_yield, ey_yield);
  graph_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_yield->GetYaxis()->SetTitle("dN_{J/#psi}/d#it{p}_{T} (GeV/c)^{-1}");
  graph_yield->AddYError(size_ptbin, eysys_yield, eysys_yield);

  TGraphMultiErrors *graph_yield_run2 = new TGraphMultiErrors(
      "graph_yields_run2_pt", "Run2", 15, x_yield_run2, y_yield_run2,
      ex_yield_run2, ex_yield_run2, ey_yield_run2, ey_yield_run2);
  graph_yield_run2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_yield_run2->GetYaxis()->SetTitle(
      "dN_{J/#psi}/d#it{p}_{T} (GeV/c)^{-1}");
  graph_yield_run2->AddYError(15, eysys_yield_run2, eysys_yield_run2);

  // Save final results
  // Saving plot for J/psi yields as function of pT
  TCanvas *c_yield = new TCanvas("jpsi_yield_pT", "jpsi_yield_pT");
  c_yield->cd();
  graph_yield->SetMarkerStyle(20);
  graph_yield->SetMarkerSize(1.);
  graph_yield->SetMarkerColor(kBlue);
  graph_yield->SetLineColor(kBlue);
  graph_yield->SetLineWidth(2);
  graph_yield->SetFillStyle(0);
  graph_yield->Draw("A P Z ; Z ; 5 s=0.5");
  // gPad->ModifiedUpdate();
  TLatex *text_yield = new TLatex();
  text_yield->SetTextSize(0.04);
  text_yield->SetTextFont(42);
  text_yield->DrawLatexNDC(.3, .82,
                           "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                           "= 5.36 TeV");
  gPad->ModifiedUpdate();
  text_yield->DrawLatexNDC(.3, .77,
                           "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                           "4");
  gPad->ModifiedUpdate();
  ls->Add(c_yield);

  // Saving plot for J/psi yields as function of pT
  // compared with Run2
  TCanvas *c_yield_vsRun2 =
      new TCanvas("jpsi_yield_raw_pT", "jpsi_yield_raw_pT");
  if (size_ptbin == 15) {
    c_yield_vsRun2->cd();
    c_yield_vsRun2->SetBottomMargin(0);
    TPad *pad_yield = new TPad("pad_yield", "pad_yield", 0, 0.3, 1, 1.0);
    pad_yield->SetBottomMargin(0);
    pad_yield->Draw();
    pad_yield->cd();
    auto mg_raw = new TMultiGraph();
    graph_yield->SetMarkerStyle(20);
    graph_yield->SetMarkerSize(1.);
    graph_yield->SetMarkerColor(kBlue);
    graph_yield->SetLineColor(kBlue);
    graph_yield->SetLineWidth(2);
    graph_yield->SetFillStyle(0);
    graph_yield_run2->SetMarkerStyle(20);
    graph_yield_run2->SetMarkerSize(1.);
    graph_yield_run2->SetMarkerColor(kRed);
    graph_yield_run2->SetLineColor(kRed);
    graph_yield_run2->SetLineWidth(2);
    graph_yield_run2->SetFillStyle(0);
    mg_raw->Add(graph_yield);
    mg_raw->Add(graph_yield_run2);
    mg_raw->Draw("A P Z ; Z ; 5 s=0.5");
    pad_yield->BuildLegend();

    c_yield_vsRun2->cd();
    TPad *pad_yield_ratio =
        new TPad("pad_yield_ratio", "pad_yield_ratio", 0, 0., 1, 0.3);
    pad_yield_ratio->SetTopMargin(0);
    pad_yield_ratio->SetBottomMargin(0.22);
    pad_yield_ratio->Draw();
    pad_yield_ratio->cd();
    TH1D *hs_yield_ratio =
        new TH1D("hist_yield_ratio", "", size_ptbin, pt_bins);
    for (int i = 0; i < size_ptbin; i++) {
      hs_yield_ratio->SetBinContent(i + 1, y_yield[i] / y_yield_run2[i]);
    }
    hs_yield_ratio->SetStats(0);
    hs_yield_ratio->SetTitle("");
    hs_yield_ratio->SetMarkerStyle(20);
    hs_yield_ratio->SetMarkerSize(0.8);
    hs_yield_ratio->GetYaxis()->SetTitleSize(0.1);
    hs_yield_ratio->GetYaxis()->SetTitle("ratio");
    hs_yield_ratio->GetYaxis()->SetTitleOffset(0.25);
    hs_yield_ratio->GetYaxis()->SetLabelSize(0.1);
    hs_yield_ratio->GetYaxis()->SetRangeUser(0., 10.);
    hs_yield_ratio->GetXaxis()->SetLabelSize(0.1);
    hs_yield_ratio->GetXaxis()->SetLabelOffset();
    hs_yield_ratio->GetXaxis()->SetTitleSize(0.1);
    hs_yield_ratio->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c2)");
    hs_yield_ratio->Draw("HIST P");
    TF1 *lratio_yield =
        new TF1("lratio_yield", "[0]", pt_bins[0], pt_bins[size_ptbin]);
    lratio_yield->SetParameter(0, 1.);
    lratio_yield->SetLineColor(kBlue);
    lratio_yield->SetLineWidth(3);
    lratio_yield->SetLineStyle(1);
    lratio_yield->Draw("same");
    pad_yield_ratio->ModifiedUpdate();
  } else {
    c_yield_vsRun2->cd();
    auto mg_raw = new TMultiGraph();
    graph_yield->SetMarkerStyle(20);
    graph_yield->SetMarkerSize(1.);
    graph_yield->SetMarkerColor(kBlue);
    graph_yield->SetLineColor(kBlue);
    graph_yield->SetLineWidth(2);
    graph_yield->SetFillStyle(0);
    graph_yield_run2->SetMarkerStyle(20);
    graph_yield_run2->SetMarkerSize(1.);
    graph_yield_run2->SetMarkerColor(kRed);
    graph_yield_run2->SetLineColor(kRed);
    graph_yield_run2->SetLineWidth(2);
    graph_yield_run2->SetFillStyle(0);
    mg_raw->Add(graph_yield);
    mg_raw->Add(graph_yield_run2);
    mg_raw->Draw("A P Z ; Z ; 5 s=0.5");
    gPad->BuildLegend();
  }
  ls->Add(c_yield_vsRun2);

  // Saving plot for J/psi v2 as function of pT and
  // compared with Run2
  TCanvas *c_pt = new TCanvas("v2_pT", "v2_pT");
  c_pt->cd();
  auto mg = new TMultiGraph();
  graph_v2pt->SetMarkerStyle(20);
  graph_v2pt->SetMarkerSize(1.);
  graph_v2pt->SetMarkerColor(kBlue);
  graph_v2pt->SetLineColor(kBlue);
  graph_v2pt->SetLineWidth(2);
  graph_v2pt->SetFillStyle(0);
  graph_v2pt_run2->SetMarkerStyle(20);
  graph_v2pt_run2->SetMarkerSize(1.);
  graph_v2pt_run2->SetMarkerColor(kRed);
  graph_v2pt_run2->SetLineColor(kRed);
  graph_v2pt_run2->SetLineWidth(2);
  graph_v2pt_run2->SetFillStyle(0);
  mg->Add(graph_v2pt);
  mg->Add(graph_v2pt_run2);
  mg->Draw("A P Z ; Z ; 5 s=0.5");
  gPad->BuildLegend();
  // gPad->ModifiedUpdate();
  TLatex *text_pt = new TLatex();
  text_pt->SetTextSize(0.04);
  text_pt->SetTextFont(42);
  text_pt->DrawLatexNDC(.18, .82,
                        "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                        "= 5.36 TeV");
  gPad->ModifiedUpdate();
  text_pt->DrawLatexNDC(.18, .77,
                        "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                        "4");
  gPad->ModifiedUpdate();
  ls->Add(c_pt);
}

//______________________________________________________________________________
void LoadData(std::string FileName, THnSparse *&hs_V2, TH2F *&hs_R2SPAB,
              TH2F *&hs_R2SPAC, TH2F *&hs_R2SPBC, std::string muonCut,
              std::string dimuonCut) {
  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist_v2, *list_hist_r2;
  TList *sublist_v2, *sublist_r2;
  if (filePath.extension() == ".root") {
    // Load data from AnalysisResults.root
    TFile *Input_File = TFile::Open(FileName.c_str());
    list_hist_v2 =
        (THashList *)Input_File->Get("analysis-same-event-pairing/output");
    sublist_v2 = (TList *)list_hist_v2->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonSEPM_%s", muonCut.c_str())
            : Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    list_hist_r2 =
        (THashList *)Input_File->Get("analysis-event-selection/output");
    sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

    // Get histograms of correlations and resolution factors
    hs_V2 = (THnSparse *)sublist_v2->FindObject("Mass_Pt_centrFT0C_V2");
    hs_R2SPAB = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C");
    hs_R2SPAC = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C");
    hs_R2SPBC = (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C");
    Input_File->Close();
  } else {
    // Load data from a list of AnalysisResults.root
    fstream InputFiles;
    InputFiles.open(FileName, ios::in);
    if (InputFiles.is_open()) {
      string File;
      cout << "Start Loading input AnalysisResults in list..." << endl;
      bool first = true;
      while (getline(InputFiles, File)) {
        cout << "Loading input from: " << File << endl;
        TFile *inFile = TFile::Open(File.c_str());
        list_hist_v2 =
            (THashList *)inFile->Get("analysis-same-event-pairing/output");
        sublist_v2 = (TList *)list_hist_v2->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPM_%s", muonCut.c_str())
                            : Form("PairsMuonSEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        list_hist_r2 =
            (THashList *)inFile->Get("analysis-event-selection/output");
        sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

        if (first) {
          hs_V2 = (THnSparse *)sublist_v2->FindObject("Mass_Pt_centrFT0C_V2");
          hs_R2SPAB = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C");
          hs_R2SPAC = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C");
          hs_R2SPBC = (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C");
        } else {
          hs_V2->Add(
              (THnSparse *)sublist_v2->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_R2SPAB->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C"));
          hs_R2SPAC->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C"));
          hs_R2SPBC->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C"));
        }

        inFile->Close();
        first = false;
      }
      InputFiles.close();
    }
  }
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