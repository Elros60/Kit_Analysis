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
void FlowAnalysis_ScalarProduct(
    int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    std::string FileName = "AnalysisResults.root", double mass_min = 2.3,
    double mass_max = 4.3, double cent_min = 10., double cent_max = 50.,
    double chi2max_mass = 2., double chi2max_v2 = 2., bool sys = false,
    std::string muonCut = "muonLowPt210SigmaPDCA", std::string dimuonCut = "") {

  // Init Helper class
  FlowAnalysis_Helper helper;

  // Load input data for analysis
  THnSparse *hs_V2;
  TH2F *hs_R2SPAB, *hs_R2SPAC, *hs_R2SPBC;
  helper.LoadData(FileName, hs_V2, hs_R2SPAB, hs_R2SPAC, hs_R2SPBC, muonCut,
                  dimuonCut);

  // Get binning information
  TAxis *massAxis = hs_V2->GetAxis(0);
  TAxis *ptAxis = hs_V2->GetAxis(1);
  TAxis *centAxis = hs_V2->GetAxis(3);
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_mass = helper.CreateBinsFromAxis(massAxis);
  double *Bin_pt = helper.CreateBinsFromAxis(ptAxis);
  double *Bin_cent = helper.CreateBinsFromAxis(centAxis);

  // Initialize fitter
  FlowAnalysis_Fitting fitter;
  fitter.init();
  fitter.setChi2MaxMass(chi2max_mass);
  fitter.setChi2MaxV2(chi2max_v2);
  fitter.setCentRange(cent_min, cent_max);

  // Define variables' range for analysis
  double Bin_pt_mass[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
  // double Bin_pt_mass[16] = {0., 0.3, 1., 2.,  3.,  4.,  5.,  6.,
  //                           7., 8.,  9., 10., 11., 12., 15., 20.};

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
  helper.LoadDataRun2(x_run2, y_run2, ex_run2, ey_run2, eysys_run2, flag_run2);

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
    TH1D *hs_mass_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2);
    TH1D *hs_v2sp = helper.GetV2(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
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
              TH1D *hist_mass_proj_sys = helper.GetMass(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                  mass_max_sys[i1], cent_min, cent_max, hs_V2);

              TH1D *hs_v2sp_sys = helper.GetV2(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
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
  helper.PlotSNRvsRun2(
      int(size(Bin_pt_mass)) - 1, Bin_pt_mass, int(size(Bin_pt_mass)) - 1,
      x_yield, SNR, int(size(x_yield_run2)), x_yield_run2, SNR_run2, l_results);

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
          helper.GetStats(nbCombo_yield, y_sys_yield[i], ey_sys_yield[i]);
      vector<double> stats_v2 =
          helper.GetStats(nbCombo_v2, y_sys_v2[i], ey_sys_v2[i]);

      // Fill pT-differential v2 and jpsi yields
      y_v2pt[i] = stats_v2[0];
      ey_v2pt[i] = stats_v2[1];
      eysys_v2pt[i] = stats_v2[2];

      y_yield[i] = stats_yield[0] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
      ey_yield[i] = stats_yield[1] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
      eysys_yield[i] = stats_yield[2] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);

      // Saving results for systematics
      helper.PlotSystematics(i, c_sys_yield[i], c_sys_v2[i], hist_sys_yield[i],
                             hist_sys_v2[i], bins_sys_yield, bins_sys_v2,
                             chi2_yield[i], chi2_v2[i], nbCombo_yield,
                             nbCombo_v2, stats_yield, stats_v2, Bin_pt_mass,
                             l_results_sys_yield, l_results_sys_v2);
    }
    f.cd();
    l_results_sys_yield->Write("FitYieldSystematics", TObject::kSingleKey);
    l_results_sys_v2->Write("FitV2Systematics", TObject::kSingleKey);

    // Saving final results
    helper.PlotFinalResults(int(size(Bin_pt_mass)) - 1, Bin_pt_mass, x_v2pt,
                            y_v2pt, ex_v2pt, ey_v2pt, eysys_v2pt, x_run2,
                            y_run2, ex_run2, ey_run2, eysys_run2, x_yield,
                            y_yield, ex_yield, ey_yield, eysys_yield,
                            x_yield_run2, y_yield_run2, ex_yield_run2,
                            ey_yield_run2, eysys_yield_run2, l_results);
  }
  f.cd();
  l_results->Write("FinalResults", TObject::kSingleKey);
  f.Close();
}