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

// Predefined binnings
vector<double> Bin_pt_mass_10 = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 20};
vector<double> Bin_pt_mass_11 = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20};
vector<double> Bin_pt_mass_15 = {0., 0.3, 1., 2.,  3.,  4.,  5.,  6.,
                                 7., 8.,  9., 10., 11., 12., 15., 20.};
vector<double> Bin_pt_mass_18 = {0,   1, 1.5, 2, 2.5, 3,  3.5, 4,  4.5, 5,
                                 5.5, 6, 7,   8, 9,   10, 12,  15, 20};
vector<double> Bin_PoolEM = {0.,  5.,  10., 20., 30., 40.,
                             50., 60., 70., 80., 90.};

//______________________________________________________________________________
void FlowAnalysis_EventMixing(
    int flag_binning, int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    int flag_run2yield, std::string FileName = "AnalysisResults.root",
    double mass_min = 2.3, double mass_max = 4.3, double cent_min = 10.,
    double cent_max = 50., double chi2max_mass = 2., double chi2max_v2 = 2.,
    bool sys = false, bool meanPt = false, bool SaveSys = false,
    std::string inputFlag = "goodmedium",
    std::string muonCut = "muonLowPt210SigmaPDCA", std::string dimuonCut = "") {

  LOG(info) << "Start flow analysis...";
  // Init Helper class
  FlowAnalysis_Helper *helper = new FlowAnalysis_Helper();

  // Load input data for analysis
  TProfile3D *tp_V2SEPM, *tp_V2SEPP, *tp_V2SEMM;
  TProfile3D *tp_V2MEPM, *tp_V2MEPP, *tp_V2MEMM;

  helper->LoadDataMEProfile(FileName, tp_V2SEPM, tp_V2SEPP, tp_V2SEMM,
                            tp_V2MEPM, tp_V2MEPP, tp_V2MEMM, muonCut,
                            dimuonCut);

  // Get general binning information
  TAxis *massAxis = tp_V2SEPM->GetXaxis();
  TAxis *ptAxis = tp_V2SEPM->GetYaxis();
  TAxis *centAxis = tp_V2SEPM->GetZaxis();
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_mass = helper->CreateBinsFromAxis(massAxis);
  double *Bin_pt = helper->CreateBinsFromAxis(ptAxis);
  double *Bin_cent = helper->CreateBinsFromAxis(centAxis);

  // Initialize fitter
  FlowAnalysis_Fitting *fitter = new FlowAnalysis_Fitting();
  fitter->init();
  fitter->setChi2MaxMass(chi2max_mass);
  fitter->setChi2MaxV2(chi2max_v2);
  fitter->setCentRange(cent_min, cent_max);

  // Define variables' range for analysis
  vector<double> Bin_pt_mass;
  switch (flag_binning) {
  case 10:
    Bin_pt_mass = Bin_pt_mass_10;
    break;
  case 11:
    Bin_pt_mass = Bin_pt_mass_11;
    break;
  case 15:
    Bin_pt_mass = Bin_pt_mass_15;
    break;
  case 18:
    Bin_pt_mass = Bin_pt_mass_11;
    break;
  default:
    Bin_pt_mass = Bin_pt_mass_10;
    break;
  }

  // Define the pool for systematics:
  // combinations
  double mass_min_sys[3] = {2.2, 2.3, 2.4};
  double mass_max_sys[3] = {4.2, 4.3, 4.4};
  string sig_enum[5] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev",
                        "EventMixing"};
  string bkg_v2_enum[2] = {"EventMixing(beta free)", "EventMixing(beta fix)"};
  int sig_mass[3] = {0, 1, 2}; // CB2(MC,data) NA60
  int bkg_mass[2] = {3, 4};    // Chebychev Event-Mixing
  int bkg_v2[2] = {0, 1};      // Event-Mixing, beta fix or free
  int nb_trials = int(size(sig_mass)) * int(size(bkg_mass)) *
                  int(size(bkg_v2)) * int(size(mass_min_sys));

  // Create output file
  TFile *f = new TFile(
      sys ? Form("FlowAnalysisResults_%s_"
                 "EventMixing%d_%s_%g_%"
                 "g_%dBinPt_%s_withSys.root",
                 inputFlag.c_str(), nb_trials, muonCut.c_str(), cent_min,
                 cent_max, int(Bin_pt_mass.size()) - 1,
                 meanPt ? "MeanPt" : "NoMeanPt")
          : Form("FlowAnalysisResults_%s_"
                 "EventMixing_%s_%g_%"
                 "g_%dBinPt_%s.root",
                 inputFlag.c_str(), muonCut.c_str(), cent_min, cent_max,
                 int(Bin_pt_mass.size()) - 1, meanPt ? "MeanPt" : "NoMeanPt"),
      "RECREATE");

  ///////////////////////////////////////////////////
  ///                                             ///
  ///   Analysis for Differential Flow of J/Psi   ///
  ///                                             ///
  ///////////////////////////////////////////////////

  LOG(info) << "Processing analysis for differential flow ...";

  // Get index of interested centrality bin
  int itmin = std::find(Bin_PoolEM.begin(), Bin_PoolEM.end(), cent_min) -
              Bin_PoolEM.begin();
  int itmax = std::find(Bin_PoolEM.begin(), Bin_PoolEM.end(), cent_max) -
              Bin_PoolEM.begin();

  // Calculate pT-integrated R factors and F factors
  LOG(info) << "Processing analysis for pT-integrated R and F factors ...";
  vector<double> ffactor;
  for (int i = 0; i < int(Bin_PoolEM.size()) - 1; i++) {
    TH1D *hist_rfactor = helper->GetRfactorProfile(
        0., 20., 0., 5., Bin_PoolEM[i], Bin_PoolEM[i + 1], tp_V2MEPM, tp_V2MEPP,
        tp_V2MEMM);
    double F_value = helper->GetFfactorProfile(
        0., 20., mass_min, mass_max, Bin_PoolEM[i], Bin_PoolEM[i + 1],
        tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
    LOG(info) << Form("F factor for [%g - %g] (%%): %g", Bin_PoolEM[i],
                      Bin_PoolEM[i + 1], F_value);
    ffactor.emplace_back(F_value);
    TList *l_rfactor = new TList();
    l_rfactor->Add(hist_rfactor);
    f->cd();
    l_rfactor->Write(
        Form("RNormalizationFactor_%g_%g", Bin_PoolEM[i], Bin_PoolEM[i + 1]),
        TObject::kSingleKey);

    delete hist_rfactor;
    delete l_rfactor;
  }

  LOG(info) << "Processing analysis for pT-differential R and F factors ...";
  vector<double *> ffactor_ptDiff;
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
    ffactor_ptDiff.emplace_back(new double[itmax - itmin]);
  }
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
    for (int j = itmin; j < itmax; j++) {
      TH1D *hist_rfactor = helper->GetRfactorProfile(
          Bin_pt_mass[i], Bin_pt_mass[i + 1], 0., 5., Bin_PoolEM[j],
          Bin_PoolEM[j + 1], tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
      double F_value = helper->GetFfactorProfile(
          Bin_pt_mass[i], Bin_pt_mass[i + 1], 2., 5., Bin_PoolEM[j],
          Bin_PoolEM[j + 1], tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
      LOG(info) << Form("F factor for [%g - %g] (%%) [%g - %g] (GeV/c): %g",
                        Bin_PoolEM[j], Bin_PoolEM[j + 1], Bin_pt_mass[i],
                        Bin_pt_mass[i + 1], F_value);
      ffactor_ptDiff[i][j - itmin] = F_value;
      delete hist_rfactor;
    }
  }
  cout << "Flag1" << endl;

  // Create histogram for pt-differential v2
  double *x_yield = new double[int(Bin_pt_mass.size()) - 1];
  double *y_yield = new double[int(Bin_pt_mass.size()) - 1];
  double *ex_yield = new double[int(Bin_pt_mass.size()) - 1];
  double *ey_yield = new double[int(Bin_pt_mass.size()) - 1];
  double *eysys_yield = new double[int(Bin_pt_mass.size()) - 1];
  double *SNR = new double[int(Bin_pt_mass.size()) - 1];
  double *x_v2pt = new double[int(Bin_pt_mass.size()) - 1];
  double *y_v2pt = new double[int(Bin_pt_mass.size()) - 1];
  double *ex_v2pt = new double[int(Bin_pt_mass.size()) - 1];
  double *ey_v2pt = new double[int(Bin_pt_mass.size()) - 1];
  double *eysys_v2pt = new double[int(Bin_pt_mass.size()) - 1];

  // Load Run2 data for comparaison
  double *x_run2, *y_run2, *ex_run2, *ey_run2, *eysys_run2;
  helper->LoadDataRun2(x_run2, y_run2, ex_run2, ey_run2, eysys_run2, flag_run2);

  double *x_yield_run2, *ex_yield_run2, *y_yield_run2, *ey_yield_run2,
      *eysys_yield_run2, *SNR_run2;
  int nbins_run2yield = 15;
  helper->LoadDataYieldRun2(x_yield_run2, y_yield_run2, ex_yield_run2,
                            ey_yield_run2, eysys_yield_run2, SNR_run2,
                            flag_run2yield);

  // Initialize arrays for each trial of systematic study
  const int nbCombo_v2 = int(size(sig_mass)) * int(size(bkg_mass)) *
                         int(size(bkg_v2)) * int(size(mass_max_sys));
  const int nbCombo_yield =
      int(size(sig_mass)) * int(size(bkg_mass)) * int(size(mass_max_sys));

  vector<double *> x_sys_pt, ex_sys_pt, y_sys_yield, ey_sys_yield, y_sys_v2,
      ey_sys_v2, chi2_yield, chi2_v2, chi2_meanPt;
  double *bins_sys_v2 = new double[nbCombo_v2 + 1];
  double *bins_sys_yield = new double[nbCombo_yield + 1];
  for (int i = 0; i < nbCombo_v2 + 1; i++) {
    bins_sys_v2[i] = i;
  }
  for (int i = 0; i < nbCombo_yield + 1; i++) {
    bins_sys_yield[i] = i;
  }

  vector<TH1D *> hist_sys_yield, hist_sys_v2, hist_sys_meanPt;
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
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

    if (meanPt) {
      x_sys_pt.emplace_back(new double[nbCombo_yield]);
      ex_sys_pt.emplace_back(new double[nbCombo_yield]);
      chi2_meanPt.emplace_back(new double[nbCombo_yield]);
      hist_sys_meanPt.emplace_back(new TH1D(
          Form("hist_sys_meanPt_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
          Form("hist_sys_meanPt_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
          nbCombo_yield, bins_sys_yield));
    }
  }
  cout << "Flag2" << endl;
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {

    // Normalisation for mixed-event spectra
    TH1D *hs_mass_mepm_proj = new TH1D();
    TH1D *hs_mass_mepp_proj = new TH1D();
    TH1D *hs_mass_memm_proj = new TH1D();
    cout << "Flag3" << endl;
    for (int j = itmin; j < itmax; j++) {
      if (j == itmin) {
        hs_mass_mepm_proj = helper->GetMassProfile(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
            Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPM, "MEPM");
        hs_mass_mepm_proj->Scale(ffactor_ptDiff[i][j - itmin]);

        hs_mass_mepp_proj = helper->GetMassProfile(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
            Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPP, "MEPP");
        hs_mass_mepp_proj->Scale(ffactor_ptDiff[i][j - itmin]);

        hs_mass_memm_proj = helper->GetMassProfile(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
            Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEMM, "MEMM");
        hs_mass_memm_proj->Scale(ffactor_ptDiff[i][j - itmin]);

      } else {
        hs_mass_mepm_proj->Add(
            helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                   mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1],
                                   tp_V2MEPM, "MEPM"),
            ffactor_ptDiff[i][j - itmin]);

        hs_mass_mepp_proj->Add(
            helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                   mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1],
                                   tp_V2MEPP, "MEPP"),
            ffactor_ptDiff[i][j - itmin]);

        hs_mass_memm_proj->Add(
            helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                   mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1],
                                   tp_V2MEMM, "MEMM"),
            ffactor_ptDiff[i][j - itmin]);
      }
      cout << "Flag3" << endl;
    }

    // Same-event profiles: mass
    TH1D *hs_mass_sepm_proj =
        helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                               mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    TH1D *hs_mass_sepp_proj =
        helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                               mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP");
    TH1D *hs_mass_semm_proj =
        helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                               mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM");

    // Same-event profiles: v2
    TH1D *hs_v2_sepm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    TH1D *hs_v2_sepp_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP");
    TH1D *hs_v2_semm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM");

    // Mixed-event profiles: v2
    TH1D *hs_v2_mepm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2MEPM, "MEPM");
    TH1D *hs_v2_mepp_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2MEPP, "MEPP");
    TH1D *hs_v2_memm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2MEMM, "MEMM");

    // Mean pT profile
    TH1D *hs_mass_sepm_proj_meanPt = new TH1D();
    if (meanPt) {
      hs_mass_sepm_proj_meanPt =
          helper->GetMeanPt(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    }
    cout << "Flag4" << endl;
    // Save plots for invariant mass
    TList *l_SE_ME = new TList();
    helper->PlotSEME("PM", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_mass_sepm_proj,
                     hs_mass_mepm_proj, l_SE_ME);
    helper->PlotSEME("PP", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_mass_sepp_proj,
                     hs_mass_mepp_proj, l_SE_ME);
    helper->PlotSEME("MM", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_mass_semm_proj,
                     hs_mass_memm_proj, l_SE_ME);
    cout << "Flag4.5" << endl;
    TList *l_SE_ME_V2 = new TList();
    helper->PlotSEME("PM", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_v2_sepm_proj,
                     hs_v2_mepm_proj, l_SE_ME_V2);
    helper->PlotSEME("PP", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_v2_sepp_proj,
                     hs_v2_mepp_proj, l_SE_ME_V2);
    helper->PlotSEME("MM", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_v2_semm_proj,
                     hs_v2_memm_proj, l_SE_ME_V2);
    cout << "Flag5" << endl;
    /// Do fitting
    // Configuration for fitting
    fitter->setModel(flag_sig, flag_bkg);
    fitter->setModelV2(flag_v2);
    fitter->setMassRange(mass_min, mass_max);
    fitter->setPtRange(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
    fitter->setOrder(2);
    fitter->setMode(0); // standard mode, no systematics

    // Fit invariant mass + v2
    TList *l_diff_fit = new TList();
    vector<double> results_v2 =
        meanPt ? fitter->runFittingEM(hs_mass_sepm_proj, hs_mass_mepm_proj,
                                      hs_v2_sepm_proj, hs_v2_mepm_proj,
                                      hs_mass_sepm_proj_meanPt, l_diff_fit)
               : fitter->runFittingEM(hs_mass_sepm_proj, hs_mass_mepm_proj,
                                      hs_v2_sepm_proj, hs_v2_mepm_proj, nullptr,
                                      l_diff_fit);

    SNR[i] = results_v2[4];
    cout << "Flag6" << endl;
    f->cd();
    l_SE_ME->SetOwner();
    l_SE_ME->Write(Form("Mass_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                   TObject::kSingleKey);
    delete l_SE_ME;

    f->cd();
    l_SE_ME_V2->SetOwner();
    l_SE_ME_V2->Write(Form("V2_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                      TObject::kSingleKey);
    delete l_SE_ME_V2;

    f->cd();
    l_diff_fit->SetOwner();
    l_diff_fit->Write(
        Form("DifferentialFlow_Fit_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    delete l_diff_fit;
    cout << "Flag7" << endl;
    delete hs_mass_sepm_proj;
    delete hs_mass_sepp_proj;
    delete hs_mass_semm_proj;
    delete hs_mass_mepm_proj;
    delete hs_mass_mepp_proj;
    delete hs_mass_memm_proj;
    delete hs_v2_sepm_proj;
    delete hs_v2_sepp_proj;
    delete hs_v2_semm_proj;
    if (hs_v2_mepm_proj) {
      delete hs_v2_mepm_proj;
    }
    if (hs_v2_mepp_proj) {
      delete hs_v2_mepp_proj;
    }
    if (hs_v2_memm_proj) {
      delete hs_v2_memm_proj;
    }
    if (hs_mass_sepm_proj_meanPt) {
      delete hs_mass_sepm_proj_meanPt;
    }
    cout << "Flag8" << endl;
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
            TString combo_meanPt = Form("{%s;%s}{Chebychev}[%g-%g]",
                                        sig_enum[sig_mass[i2]].c_str(),
                                        sig_enum[bkg_mass[i3]].c_str(),
                                        mass_min_sys[i1], mass_max_sys[i1]);
            for (int i4 = 0; i4 < int(size(bkg_v2)); i4++) {
              // Get F factors
              // Calculate R factors and F factors
              LOG(info) << "Processing analysis for R and F factors ...";
              vector<double> ffactor_sys;
              for (int k = 0; k < int(Bin_pt_mass.size()) - 1; k++) {
                TH1D *hist_rfactor_sys = helper->GetRfactorProfile(
                    Bin_pt_mass[k], Bin_pt_mass[k + 1], mass_min_sys[i1],
                    mass_max_sys[i1], cent_min, cent_max, tp_V2MEPM, tp_V2MEPP,
                    tp_V2MEMM);
                double F_value = helper->GetFfactorProfile(
                    Bin_pt_mass[k], Bin_pt_mass[k + 1], mass_min_sys[i1],
                    mass_max_sys[i1], cent_min, cent_max, tp_V2SEPP, tp_V2SEMM,
                    tp_V2MEPM, hist_rfactor_sys);

                ffactor_sys.emplace_back(F_value);
                delete hist_rfactor_sys;
              }

              TH1D *hist_rfactor_sys = helper->GetRfactorProfile(
                  Bin_pt_mass[0], Bin_pt_mass[int(Bin_pt_mass.size()) - 1],
                  mass_min_sys[i1], mass_max_sys[i1], cent_min, cent_max,
                  tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
              double F_value_sys = helper->GetFfactorProfile(
                  Bin_pt_mass[0], Bin_pt_mass[int(Bin_pt_mass.size()) - 1],
                  mass_min_sys[i1], mass_max_sys[i1], cent_min, cent_max,
                  tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor_sys);
              delete hist_rfactor_sys;

              // Get mass and v2 to fit
              // Same-event profiles: mass
              TH1D *hs_mass_sepm_proj_sys = helper->GetMassProfile(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                  mass_max_sys[i1], cent_min, cent_max, tp_V2SEPM, "SEPM");

              // Same-event profiles: v2
              TH1D *hs_v2_sepm_proj_sys = helper->GetV2Profile(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                  mass_max_sys[i1], cent_min, cent_max, tp_V2SEPM, "SEPM");

              // Mixed-event profiles: mass
              TH1D *hs_mass_mepm_proj_sys = helper->GetMassProfile(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                  mass_max_sys[i1], cent_min, cent_max, tp_V2MEPM, "MEPM");

              // Mixed-event profiles: v2
              TH1D *hs_v2_mepm_proj_sys = helper->GetV2Profile(
                  Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min_sys[i1],
                  mass_max_sys[i1], cent_min, cent_max, tp_V2MEPM, "MEPM");

              // Scale mixed-event spectra with F factor
              // hs_mass_mepm_proj_sys->Scale(ffactor_sys[i]);
              hs_mass_mepm_proj_sys->Scale(F_value_sys);

              // Mean pT profile
              TH1D *hs_mass_sepm_proj_meanPt_sys = new TH1D();
              if (meanPt) {
                hs_mass_sepm_proj_meanPt_sys = helper->GetMeanPt(
                    Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                    cent_min, cent_max, tp_V2SEPM, "MeanPt");
              }

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
              fitter->setOrder(2);
              fitter->setMassRange(mass_min_sys[i1], mass_max_sys[i1]);
              fitter->setPtRange(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
              fitter->setModel(sig_mass[i2], bkg_mass[i3]);
              fitter->setModelV2(bkg_v2[i4]);
              fitter->setMode(1);
              TList *l_diff_sys = new TList();
              LOG(info) << "Processing fitting(systematic) for v2{SP} ...";
              vector<double> results_sys_v2 =
                  meanPt ? fitter->runFittingEM(
                               hs_mass_sepm_proj_sys, hs_mass_mepm_proj_sys,
                               hs_v2_sepm_proj_sys, hs_v2_mepm_proj_sys,
                               hs_mass_sepm_proj_meanPt_sys, l_diff_sys)
                         : fitter->runFittingEM(
                               hs_mass_sepm_proj_sys, hs_mass_mepm_proj_sys,
                               hs_v2_sepm_proj_sys, hs_v2_mepm_proj_sys,
                               nullptr, l_diff_sys);

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
              if (meanPt) {
                x_sys_pt[i][index_sys_yield] = results_sys_v2[7];
                ex_sys_pt[i][index_sys_yield] = results_sys_v2[8];
                chi2_meanPt[i][index_sys_v2] = results_sys_v2[9];
                hist_sys_meanPt[i]->GetXaxis()->SetBinLabel(index_sys_yield + 1,
                                                            combo_meanPt);
                hist_sys_meanPt[i]->SetBinContent(index_sys_yield + 1,
                                                  x_sys_pt[i][index_sys_yield]);
                hist_sys_meanPt[i]->SetBinError(index_sys_yield + 1,
                                                ex_sys_pt[i][index_sys_yield]);
              }
              f->cd();
              l_diff_sys->SetOwner();
              l_diff_sys->Write(Form("FitSys_%g_%g_%s", Bin_pt_mass[i],
                                     Bin_pt_mass[i + 1], combo_v2.Data()),
                                TObject::kSingleKey);
              delete l_diff_sys;
              delete hs_mass_sepm_proj_sys;
              delete hs_mass_mepm_proj_sys;
              delete hs_v2_sepm_proj_sys;
              delete hs_v2_mepm_proj_sys;
              if (hs_mass_sepm_proj_meanPt_sys) {
                delete hs_mass_sepm_proj_meanPt_sys;
              }
              index_sys_v2++;
            }
            index_sys_yield++;
          }
        }
      }
    }
  }

  // Save plots for systematics
  TList *l_results = new TList();
  TList *l_results_sys_yield = new TList();
  TList *l_results_sys_v2 = new TList();
  TList *l_results_sys_meanPt = new TList();
  if (sys) {
    for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
      vector<double> stats_yield =
          helper->GetStats(nbCombo_yield, y_sys_yield[i], ey_sys_yield[i]);
      vector<double> stats_v2 =
          helper->GetStats(nbCombo_v2, y_sys_v2[i], ey_sys_v2[i]);
      vector<double> stats_meanPt;
      if (meanPt) {
        stats_meanPt =
            helper->GetStats(nbCombo_yield, x_sys_pt[i], ex_sys_pt[i]);
      }

      // Fill pT-differential v2 and jpsi yields
      x_yield[i] =
          meanPt ? stats_meanPt[0] : (Bin_pt_mass[i + 1] + Bin_pt_mass[i]) / 2.;
      ex_yield[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2.;
      x_v2pt[i] =
          meanPt ? stats_meanPt[0] : (Bin_pt_mass[i + 1] + Bin_pt_mass[i]) / 2.;
      // ex_v2pt[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2.;
      ex_v2pt[i] = 0.25;

      y_v2pt[i] = stats_v2[0];
      ey_v2pt[i] = stats_v2[1];
      eysys_v2pt[i] = stats_v2[2];

      y_yield[i] = stats_yield[0] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
      ey_yield[i] = stats_yield[1] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
      eysys_yield[i] = stats_yield[2] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);

      // Saving results for systematics
      if (meanPt) {
        helper->PlotSystematics(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], cent_min, cent_max,
            int(Bin_pt_mass.size()) - 1, i, hist_sys_yield[i], hist_sys_v2[i],
            hist_sys_meanPt[i], bins_sys_yield, bins_sys_v2, chi2_yield[i],
            chi2_v2[i], chi2_meanPt[i], nbCombo_yield, nbCombo_v2, stats_yield,
            stats_v2, stats_meanPt,
            reinterpret_cast<double *>(Bin_pt_mass.data()), l_results_sys_yield,
            l_results_sys_v2, l_results_sys_meanPt, SaveSys);
      } else {
        helper->PlotSystematics(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], cent_min, cent_max,
            int(Bin_pt_mass.size()) - 1, i, hist_sys_yield[i], hist_sys_v2[i],
            nullptr, bins_sys_yield, bins_sys_v2, chi2_yield[i], chi2_v2[i],
            nullptr, nbCombo_yield, nbCombo_v2, stats_yield, stats_v2,
            stats_meanPt, reinterpret_cast<double *>(Bin_pt_mass.data()),
            l_results_sys_yield, l_results_sys_v2, nullptr, SaveSys);
      }
    }
    f->cd();
    l_results_sys_yield->SetOwner();
    l_results_sys_yield->Write("FitYieldSystematics", TObject::kSingleKey);

    if (meanPt) {
      f->cd();
      l_results_sys_meanPt->SetOwner();
      l_results_sys_meanPt->Write("FitMeanPtSystematics", TObject::kSingleKey);
    }

    f->cd();
    l_results_sys_v2->SetOwner();
    l_results_sys_v2->Write("FitV2Systematics", TObject::kSingleKey);

    // Saving final results
    helper->PlotFinalResults(int(Bin_pt_mass.size()) - 1, cent_min, cent_max,
                             reinterpret_cast<double *>(Bin_pt_mass.data()),
                             x_v2pt, y_v2pt, ex_v2pt, ey_v2pt, eysys_v2pt,
                             x_run2, y_run2, ex_run2, ey_run2, eysys_run2,
                             x_yield, y_yield, ex_yield, ey_yield, eysys_yield,
                             x_yield_run2, y_yield_run2, ex_yield_run2,
                             ey_yield_run2, eysys_yield_run2, l_results);
    f->cd();
    l_results->SetOwner();
    l_results->Write("FinalResults", TObject::kSingleKey);
  }
  f->Close();
  LOG(info) << "Analysis done."; // this is a temporary solution to get rid of
                                 // destructor issue with in-list tcanvas
  delete f;
}