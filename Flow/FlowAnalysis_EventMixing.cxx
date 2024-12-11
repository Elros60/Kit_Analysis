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
void FlowAnalysis_EventMixing(
    int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    std::string FileName = "AnalysisResults.root", double mass_min = 2.3,
    double mass_max = 4.3, double cent_min = 10., double cent_max = 50.,
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
  fitter.setCentRange(cent_min, cent_max);

  // Define variables' range for analysis
  double Bin_pt_mass[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};

  // Create output file
  TFile f(sys ? Form("FlowAnalysisResults_"
                     "EventMixing%s_%g_%"
                     "g_%dBinPt_withSys.root",
                     muonCut.c_str(), cent_min, cent_max,
                     int(size(Bin_pt_mass)) - 1)
              : Form("FlowAnalysisResults_"
                     "EventMixing%s_%g_%"
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

  // Calculate R factors and F factors
  LOG(info) << "Processing analysis for R and F factors ...";
  vector<double> ffactor;
  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    TH1D *hist_rfactor = helper.GetRfactorProfile(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
    double F_value = helper.GetFfactorProfile(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
    LOG(info) << Form("F factor for [%g - %g] (Gev/c): %g", Bin_pt_mass[i],
                      Bin_pt_mass[i + 1], F_value);
    ffactor.emplace_back(F_value);
    TList *l_rfactor = new TList();
    l_rfactor->Add(hist_rfactor);
    f.cd();
    l_rfactor->Write(
        Form("RNormalizationFactor_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);

    delete hist_rfactor;
    delete l_rfactor;
  }

  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    // Same-event profiles: mass
    TH1D *hs_mass_sepm_proj =
        helper.GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                              mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    TH1D *hs_mass_sepp_proj =
        helper.GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                              mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP");
    TH1D *hs_mass_semm_proj =
        helper.GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                              mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM");

    // Same-event profiles: v2
    TH1D *hs_v2_sepm_proj =
        helper.GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    TH1D *hs_v2_sepp_proj =
        helper.GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP");
    TH1D *hs_v2_semm_proj =
        helper.GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM");

    // Mixed-event profiles: mass
    TH1D *hs_mass_mepm_proj =
        helper.GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                              mass_max, cent_min, cent_max, tp_V2MEPM, "MEPM");
    TH1D *hs_mass_mepp_proj =
        helper.GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                              mass_max, cent_min, cent_max, tp_V2MEPP, "MEPP");
    TH1D *hs_mass_memm_proj =
        helper.GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                              mass_max, cent_min, cent_max, tp_V2MEMM, "MEMM");

    // Mixed-event profiles: v2
    TH1D *hs_v2_mepm_proj =
        helper.GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2MEPM, "MEPM");
    TH1D *hs_v2_mepp_proj =
        helper.GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2MEPP, "MEPP");
    TH1D *hs_v2_memm_proj =
        helper.GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2MEMM, "MEMM");

    // Scale mixed-event spectra with F factor
    hs_mass_mepm_proj->Scale(ffactor[i]);
    hs_mass_mepp_proj->Scale(ffactor[i]);
    hs_mass_memm_proj->Scale(ffactor[i]);
    // hs_v2_mepm_proj->Scale(ffactor[i]);
    // hs_v2_mepp_proj->Scale(ffactor[i]);
    // hs_v2_memm_proj->Scale(ffactor[i]);

    // Save plots for invariant mass
    TList *l_SE_ME = new TList();
    helper.PlotSEME("PM", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                    mass_max, cent_min, cent_max, hs_mass_sepm_proj,
                    hs_mass_mepm_proj, l_SE_ME);
    helper.PlotSEME("PP", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                    mass_max, cent_min, cent_max, hs_mass_sepp_proj,
                    hs_mass_mepp_proj, l_SE_ME);
    helper.PlotSEME("MM", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                    mass_max, cent_min, cent_max, hs_mass_semm_proj,
                    hs_mass_memm_proj, l_SE_ME);
    f.cd();
    l_SE_ME->Write(Form("Mass_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
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
    l_SE_ME_V2->Write(Form("V2_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                      TObject::kSingleKey);
    delete l_SE_ME_V2;

    /// Do fitting
    // Configuration for fitting
    fitter.setModel(flag_sig, flag_bkg);
    fitter.setModelV2(flag_v2);
    fitter.setMassRange(mass_min, mass_max);
    fitter.setPtRange(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
    fitter.setOrder(2);
    fitter.setMode(0); // standard mode, no systematics

    // Fit invariant mass + v2
    TList *l_diff_fit = new TList();
    vector<double> results_v2 = fitter.runFittingEM(
        hs_mass_mepm_proj->Integral(), hs_mass_sepm_proj, hs_mass_mepm_proj,
        hs_v2_sepm_proj, hs_v2_mepm_proj, l_diff_fit);

    f.cd();
    l_diff_fit->Write(
        Form("DifferentialFlow_Fit_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);

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
    delete l_diff_fit;
  }

  f.Close();
}