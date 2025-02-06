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
void FlowAnalysis_EventMixing_CentDiff(
    int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    std::string FileName = "AnalysisResults.root", double mass_min = 2.3,
    double mass_max = 4.3, double pt_min = 0., double pt_max = 5.,
    double chi2max_mass = 2., double chi2max_v2 = 2., bool sys = false,
    std::string muonCut = "muonLowPt210SigmaPDCA", std::string dimuonCut = "") {
  // Init Helper class
  FlowAnalysis_Helper helper;

  // Load input data for analysis
  TProfile3D *tp_V2SEPM, *tp_V2SEPP, *tp_V2SEMM;
  TProfile3D *tp_V2MEPM, *tp_V2MEPP, *tp_V2MEMM;

  helper.LoadDataMEProfile(FileName, tp_V2SEPM, tp_V2SEPP, tp_V2SEMM, tp_V2MEPM,
                           tp_V2MEPP, tp_V2MEMM, muonCut, dimuonCut);

  // Get general binning information
  TAxis *massAxis = tp_V2SEPM->GetXaxis();
  TAxis *ptAxis = tp_V2SEPM->GetYaxis();
  TAxis *centAxis = tp_V2SEPM->GetZaxis();
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
  fitter.setPtRange(pt_min, pt_max);

  // Define variables' range for analysis
  double Bin_cent_mass[16] = {0,  5,  10, 15, 20, 25, 30, 35,
                              40, 45, 50, 55, 60, 70, 80, 90};

  // Define the pool for systematics: 36
  // combinations
  double mass_min_sys[3] = {2.2, 2.3, 2.4};
  double mass_max_sys[3] = {4.2, 4.3, 4.4};
  string sig_enum[5] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev",
                        "EventMixing"};
  string bkg_v2_enum[1] = {"EventMixing"};
  int sig_mass[3] = {0, 1, 2}; // CB2(MC,data) NA60
  int bkg_mass[2] = {3, 4};    // Chebychev Event-Mixing
  int bkg_v2[1] = {2};         // Event-Mixing

  // Create output file
  TFile f(
      sys ? Form("FlowAnalysisResults_"
                 "EventMixing_CentDiff_%s_%g_%"
                 "g_%dBinPt_withSys.root",
                 muonCut.c_str(), pt_min, pt_max, int(size(Bin_cent_mass)) - 1)
          : Form("FlowAnalysisResults_"
                 "EventMixing_CentDiff_%s_%g_%"
                 "g_%dBinPt.root",
                 muonCut.c_str(), pt_min, pt_max, int(size(Bin_cent_mass)) - 1),
      "RECREATE");

  ///////////////////////////////////////////////////
  ///                                             ///
  ///   Analysis for Differential Flow of J/Psi   ///
  ///                                             ///
  ///////////////////////////////////////////////////

  LOG(info) << "Processing analysis for differential flow ...";

  // Calculate R factors and F factors
  LOG(info) << "Processing analysis for R and F factors ...";
  vector<double> ffactor;
  for (int i = 0; i < int(size(Bin_cent_mass)) - 1; i++) {
    TH1D *hist_rfactor = helper.GetRfactorProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
    double F_value = helper.GetFfactorProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
    LOG(info) << Form("F factor for [%g - %g] (%%): %g", Bin_cent_mass[i],
                      Bin_cent_mass[i + 1], F_value);
    ffactor.emplace_back(F_value);
    TList *l_rfactor = new TList();
    l_rfactor->Add(hist_rfactor);
    f.cd();
    l_rfactor->Write(Form("RNormalizationFactor_%g_%g", Bin_cent_mass[i],
                          Bin_cent_mass[i + 1]),
                     TObject::kSingleKey);

    delete hist_rfactor;
    delete l_rfactor;
  }

  // Create histogram for pt-differential v2
  TList *l_results = new TList();
  double *x_yield = new double[int(size(Bin_cent_mass)) - 1];
  double *y_yield = new double[int(size(Bin_cent_mass)) - 1];
  double *ex_yield = new double[int(size(Bin_cent_mass)) - 1];
  double *ey_yield = new double[int(size(Bin_cent_mass)) - 1];
  double *eysys_yield = new double[int(size(Bin_cent_mass)) - 1];
  double *SNR = new double[int(size(Bin_cent_mass)) - 1];
  double *x_v2cent = new double[int(size(Bin_cent_mass)) - 1];
  double *y_v2cent = new double[int(size(Bin_cent_mass)) - 1];
  double *ex_v2cent = new double[int(size(Bin_cent_mass)) - 1];
  double *ey_v2cent = new double[int(size(Bin_cent_mass)) - 1];
  double *eysys_v2cent = new double[int(size(Bin_cent_mass)) - 1];

  // Load Run2 data for comparaison
  double *x_run2, *y_run2, *ex_run2, *ey_run2, *eysys_run2;
  helper.LoadDataRun2Cent(x_run2, y_run2, ex_run2, ey_run2, eysys_run2,
                          flag_run2);

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
  for (int i = 0; i < int(size(Bin_cent_mass)) - 1; i++) {
    y_sys_yield.emplace_back(new double[nbCombo_yield]);
    ey_sys_yield.emplace_back(new double[nbCombo_yield]);

    y_sys_v2.emplace_back(new double[nbCombo_v2]);
    ey_sys_v2.emplace_back(new double[nbCombo_v2]);

    chi2_yield.emplace_back(new double[nbCombo_yield]);
    chi2_v2.emplace_back(new double[nbCombo_v2]);

    hist_sys_yield.emplace_back(new TH1D(
        Form("hist_sys_yield_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        Form("hist_sys_yield_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        nbCombo_yield, bins_sys_yield));
    hist_sys_v2.emplace_back(new TH1D(
        Form("hist_sys_v2_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        Form("hist_sys_v2_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        nbCombo_v2, bins_sys_v2));
  }

  for (int i = 0; i < int(size(Bin_cent_mass)) - 1; i++) {
    // Same-event profiles: mass
    TH1D *hs_mass_sepm_proj = helper.GetMassProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEPM, "SEPM");
    TH1D *hs_mass_sepp_proj = helper.GetMassProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEPP, "SEPP");
    TH1D *hs_mass_semm_proj = helper.GetMassProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEMM, "SEMM");

    // Same-event profiles: v2
    TH1D *hs_v2_sepm_proj = helper.GetV2Profile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEPM, "SEPM");
    TH1D *hs_v2_sepp_proj = helper.GetV2Profile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEPP, "SEPP");
    TH1D *hs_v2_semm_proj = helper.GetV2Profile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2SEMM, "SEMM");

    // Mixed-event profiles: mass
    TH1D *hs_mass_mepm_proj = helper.GetMassProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEPM, "MEPM");
    TH1D *hs_mass_mepp_proj = helper.GetMassProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEPP, "MEPP");
    TH1D *hs_mass_memm_proj = helper.GetMassProfile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEMM, "MEMM");

    // Mixed-event profiles: v2
    TH1D *hs_v2_mepm_proj = helper.GetV2Profile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEPM, "MEPM");
    TH1D *hs_v2_mepp_proj = helper.GetV2Profile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEPP, "MEPP");
    TH1D *hs_v2_memm_proj = helper.GetV2Profile(
        pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
        Bin_cent_mass[i + 1], tp_V2MEMM, "MEMM");

    // Scale mixed-event spectra with F factor
    hs_mass_mepm_proj->Scale(ffactor[i]);
    hs_mass_mepp_proj->Scale(ffactor[i]);
    hs_mass_memm_proj->Scale(ffactor[i]);

    // Save plots for invariant mass
    TList *l_SE_ME = new TList();
    helper.PlotSEME("PM", pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
                    Bin_cent_mass[i + 1], hs_mass_sepm_proj, hs_mass_mepm_proj,
                    l_SE_ME);
    helper.PlotSEME("PP", pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
                    Bin_cent_mass[i + 1], hs_mass_sepp_proj, hs_mass_mepp_proj,
                    l_SE_ME);
    helper.PlotSEME("MM", pt_min, pt_max, mass_min, mass_max, Bin_cent_mass[i],
                    Bin_cent_mass[i + 1], hs_mass_semm_proj, hs_mass_memm_proj,
                    l_SE_ME);
    f.cd();
    l_SE_ME->Write(
        Form("Mass_SEME_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        TObject::kSingleKey);
    delete l_SE_ME;

    TList *l_SE_ME_V2 = new TList();
    l_SE_ME_V2->Add(hs_v2_sepm_proj);
    l_SE_ME_V2->Add(hs_v2_sepp_proj);
    l_SE_ME_V2->Add(hs_v2_semm_proj);
    l_SE_ME_V2->Add(hs_v2_mepm_proj);
    l_SE_ME_V2->Add(hs_v2_mepp_proj);
    l_SE_ME_V2->Add(hs_v2_memm_proj);
    f.cd();
    l_SE_ME_V2->Write(
        Form("V2_SEME_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        TObject::kSingleKey);
    delete l_SE_ME_V2;

    /// Do fitting
    // Configuration for fitting
    fitter.setModel(flag_sig, flag_bkg);
    fitter.setModelV2(flag_v2);
    fitter.setMassRange(mass_min, mass_max);
    fitter.setCentRange(Bin_cent_mass[i], Bin_cent_mass[i + 1]);
    fitter.setOrder(2);
    fitter.setMode(0); // standard mode, no systematics

    // Fit invariant mass + v2
    TList *l_diff_fit = new TList();
    vector<double> results_v2 = fitter.runFittingEMNoMeanPt(
        hs_mass_sepm_proj, hs_mass_mepm_proj, hs_v2_sepm_proj, hs_v2_mepm_proj,
        l_diff_fit);

    SNR[i] = results_v2[4];
    x_yield[i] = (Bin_cent_mass[i] + Bin_cent_mass[i + 1]) / 2;
    ex_yield[i] = (Bin_cent_mass[i + 1] - Bin_cent_mass[i]) / 2;
    x_v2cent[i] = (Bin_cent_mass[i] + Bin_cent_mass[i + 1]) / 2;
    ex_v2cent[i] = (Bin_cent_mass[i + 1] - Bin_cent_mass[i]) / 2;

    f.cd();
    l_diff_fit->SetOwner();
    l_diff_fit->Write(Form("DifferentialFlow_Fit_%g_%g", Bin_cent_mass[i],
                           Bin_cent_mass[i + 1]),
                      TObject::kSingleKey);
    delete l_diff_fit;

    delete hs_mass_sepm_proj;
    delete hs_mass_sepp_proj;
    delete hs_mass_semm_proj;
    delete hs_mass_mepm_proj;
    delete hs_mass_mepp_proj;
    delete hs_mass_memm_proj;
    delete hs_v2_sepm_proj;
    delete hs_v2_sepp_proj;
    delete hs_v2_semm_proj;
    delete hs_v2_mepm_proj;
    delete hs_v2_mepp_proj;
    delete hs_v2_memm_proj;

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
              // Get F factors
              // Calculate R factors and F factors
              LOG(info) << "Processing analysis for R and F factors ...";
              vector<double> ffactor_sys;
              for (int k = 0; k < int(size(Bin_cent_mass)) - 1; k++) {
                TH1D *hist_rfactor_sys = helper.GetRfactorProfile(
                    pt_min, pt_max, mass_min_sys[i1], mass_max_sys[i1],
                    Bin_cent_mass[k], Bin_cent_mass[k + 1], tp_V2MEPM,
                    tp_V2MEPP, tp_V2MEMM);
                double F_value = helper.GetFfactorProfile(
                    pt_min, pt_max, mass_min_sys[i1], mass_max_sys[i1],
                    Bin_cent_mass[k], Bin_cent_mass[k + 1], tp_V2SEPP,
                    tp_V2SEMM, tp_V2MEPM, hist_rfactor_sys);

                ffactor_sys.emplace_back(F_value);
                delete hist_rfactor_sys;
              }

              // Get mass and v2 to fit
              // Same-event profiles: mass
              TH1D *hs_mass_sepm_proj_sys = helper.GetMassProfile(
                  pt_min, pt_max, mass_min_sys[i1], mass_max_sys[i1],
                  Bin_cent_mass[i], Bin_cent_mass[i + 1], tp_V2SEPM, "SEPM");

              // Same-event profiles: v2
              TH1D *hs_v2_sepm_proj_sys = helper.GetV2Profile(
                  pt_min, pt_max, mass_min_sys[i1], mass_max_sys[i1],
                  Bin_cent_mass[i], Bin_cent_mass[i + 1], tp_V2SEPM, "SEPM");

              // Mixed-event profiles: mass
              TH1D *hs_mass_mepm_proj_sys = helper.GetMassProfile(
                  pt_min, pt_max, mass_min_sys[i1], mass_max_sys[i1],
                  Bin_cent_mass[i], Bin_cent_mass[i + 1], tp_V2MEPM, "MEPM");

              // Mixed-event profiles: v2
              TH1D *hs_v2_mepm_proj_sys = helper.GetV2Profile(
                  pt_min, pt_max, mass_min_sys[i1], mass_max_sys[i1],
                  Bin_cent_mass[i], Bin_cent_mass[i + 1], tp_V2MEPM, "MEPM");

              // Scale mixed-event spectra with F factor
              hs_mass_mepm_proj_sys->Scale(ffactor_sys[i]);

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
              fitter.setCentRange(Bin_cent_mass[i], Bin_cent_mass[i + 1]);
              fitter.setModel(sig_mass[i2], bkg_mass[i3]);
              fitter.setModelV2(bkg_v2[i4]);
              fitter.setMode(1);
              TList *l_diff_sys = new TList();
              LOG(info) << "Processing fitting(systematic) for v2{SP} ...";
              vector<double> results_sys_v2 = fitter.runFittingEMNoMeanPt(
                  hs_mass_sepm_proj_sys, hs_mass_mepm_proj_sys,
                  hs_v2_sepm_proj_sys, hs_v2_mepm_proj_sys, l_diff_sys);

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
              l_diff_sys->SetOwner();
              l_diff_sys->Write(Form("FitSys_%g_%g_%s", Bin_cent_mass[i],
                                     Bin_cent_mass[i + 1], combo_v2.Data()),
                                TObject::kSingleKey);

              delete l_diff_sys;
              delete hs_mass_sepm_proj_sys;
              delete hs_mass_mepm_proj_sys;
              delete hs_v2_sepm_proj_sys;
              delete hs_v2_mepm_proj_sys;

              index_sys_v2++;
            }
            index_sys_yield++;
          }
        }
      }
    }
  }

// Save plots for systematics
if (sys) {
  TList *l_results_sys_yield = new TList();
  TList *l_results_sys_v2 = new TList();
  vector<TCanvas *> c_sys_yield;
  vector<TCanvas *> c_sys_v2;
  for (int i = 0; i < int(size(Bin_cent_mass)) - 1; i++) {
    c_sys_yield.emplace_back(new TCanvas(
        Form("Sys_yield_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        Form("Sys_yield_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1])));
    c_sys_v2.emplace_back(new TCanvas(
        Form("Sys_v2_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1]),
        Form("Sys_v2_%g_%g", Bin_cent_mass[i], Bin_cent_mass[i + 1])));
  }
  for (int i = 0; i < int(size(Bin_cent_mass)) - 1; i++) {
    vector<double> stats_yield =
        helper.GetStats(nbCombo_yield, y_sys_yield[i], ey_sys_yield[i]);
    vector<double> stats_v2 =
        helper.GetStats(nbCombo_v2, y_sys_v2[i], ey_sys_v2[i]);

    // Fill pT-differential v2 and jpsi yields
    y_v2cent[i] = stats_v2[0];
    ey_v2cent[i] = stats_v2[1];
    eysys_v2cent[i] = stats_v2[2];

    y_yield[i] = stats_yield[0] / (Bin_cent_mass[i + 1] - Bin_cent_mass[i]);
    ey_yield[i] = stats_yield[1] / (Bin_cent_mass[i + 1] - Bin_cent_mass[i]);
    eysys_yield[i] = stats_yield[2] / (Bin_cent_mass[i + 1] - Bin_cent_mass[i]);

    // Saving results for systematics
    helper.PlotSystematicsNoMeanPt(
        i, c_sys_yield[i], c_sys_v2[i], hist_sys_yield[i], hist_sys_v2[i],
        bins_sys_yield, bins_sys_v2, chi2_yield[i], chi2_v2[i], nbCombo_yield,
        nbCombo_v2, stats_yield, stats_v2, Bin_cent_mass, l_results_sys_yield,
        l_results_sys_v2, "cent");
  }
  f.cd();
  l_results_sys_yield->Write("FitYieldSystematics", TObject::kSingleKey);
  l_results_sys_v2->Write("FitV2Systematics", TObject::kSingleKey);

  // Saving final results
  helper.PlotFinalResultsCent(int(size(Bin_cent_mass)) - 1, pt_min, pt_max,
                              Bin_cent_mass, x_v2cent, y_v2cent, ex_v2cent,
                              ey_v2cent, eysys_v2cent, x_run2, y_run2, ex_run2,
                              ey_run2, eysys_run2, x_yield, y_yield, ex_yield,
                              ey_yield, eysys_yield, l_results);
}
f.cd();
l_results->SetOwner();
l_results->Write("FinalResults", TObject::kSingleKey);
delete l_results;
f.Close();
}