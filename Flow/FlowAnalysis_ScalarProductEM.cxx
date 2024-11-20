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
void FlowAnalysis_ScalarProductEM(
    int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    std::string FileName = "AnalysisResults.root", double mass_min = 2.3,
    double mass_max = 4.3, double cent_min = 10., double cent_max = 50.,
    double chi2max_mass = 2., double chi2max_v2 = 2., bool sys = false,
    std::string muonCut = "muonLowPt210SigmaPDCA", std::string dimuonCut = "") {
  // Init Helper class
  FlowAnalysis_Helper helper;

  // Load input data for analysis
  THnSparse *hs_V2SEPM, *hs_V2SEPP, *hs_V2SEMM;
  THnSparse *hs_V2MEPM, *hs_V2MEPP, *hs_V2MEMM;
  TH2F *hs_R2SPAB, *hs_R2SPAC, *hs_R2SPBC;
  THnSparse *hs_u2q2_cosDeltaPhi_MEPM1, *hs_u2q2_cosDeltaPhi_MEPP1,
      *hs_u2q2_cosDeltaPhi_MEMM1;
  THnSparse *hs_u2q2_cosDeltaPhi_MEPM2, *hs_u2q2_cosDeltaPhi_MEPP2,
      *hs_u2q2_cosDeltaPhi_MEMM2;
  TH3F *hs_r2spABMEPM1, *hs_r2spABMEPP1, *hs_r2spABMEMM1;
  TH3F *hs_r2spACMEPM1, *hs_r2spACMEPP1, *hs_r2spACMEMM1;
  TH3F *hs_r2spBCMEPM1, *hs_r2spBCMEPP1, *hs_r2spBCMEMM1;
  TH3F *hs_r2spABMEPM2, *hs_r2spABMEPP2, *hs_r2spABMEMM2;
  TH3F *hs_r2spACMEPM2, *hs_r2spACMEPP2, *hs_r2spACMEMM2;
  TH3F *hs_r2spBCMEPM2, *hs_r2spBCMEPP2, *hs_r2spBCMEMM2;

  helper.LoadDataME(
      FileName, hs_V2SEPM, hs_V2SEPP, hs_V2SEMM, hs_V2MEPM, hs_V2MEPP,
      hs_V2MEMM, hs_R2SPAB, hs_R2SPAC, hs_R2SPBC, hs_u2q2_cosDeltaPhi_MEPM1,
      hs_u2q2_cosDeltaPhi_MEPP1, hs_u2q2_cosDeltaPhi_MEMM1,
      hs_u2q2_cosDeltaPhi_MEPM2, hs_u2q2_cosDeltaPhi_MEPP2,
      hs_u2q2_cosDeltaPhi_MEMM2, hs_r2spABMEPM1, hs_r2spABMEPP1, hs_r2spABMEMM1,
      hs_r2spACMEPM1, hs_r2spACMEPP1, hs_r2spACMEMM1, hs_r2spBCMEPM1,
      hs_r2spBCMEPP1, hs_r2spBCMEMM1, hs_r2spABMEPM2, hs_r2spABMEPP2,
      hs_r2spABMEMM2, hs_r2spACMEPM2, hs_r2spACMEPP2, hs_r2spACMEMM2,
      hs_r2spBCMEPM2, hs_r2spBCMEPP2, hs_r2spBCMEMM2, muonCut, dimuonCut);

  // Get general binning information
  TAxis *massAxis = hs_V2SEPM->GetAxis(0);
  TAxis *ptAxis = hs_V2SEPM->GetAxis(1);
  TAxis *centAxis = hs_V2SEPM->GetAxis(3);
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
                     "ScalarProductEM_%s_%g_%"
                     "g_%dBinPt_withSys.root",
                     muonCut.c_str(), cent_min, cent_max,
                     int(size(Bin_pt_mass)) - 1)
              : Form("FlowAnalysisResults_"
                     "ScalarProductEM_%s_%g_%"
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

  // Calculate R factors and F factors
  LOG(info) << "Processing analysis for R and F factors ...";
  vector<double> ffactor;
  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {
    TH1D *hist_rfactor = helper.GetRfactor(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, hs_V2MEPM, hs_V2MEPP, hs_V2MEMM);
    double F_value = helper.GetFfactor(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, hs_V2SEPP, hs_V2SEMM, hs_V2MEPM, hist_rfactor);
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
    TList *l_SE_ME = new TList();
    // Same-event profiles: mass
    TH1D *hs_mass_sepm_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2SEPM, "SEPM");
    TH1D *hs_mass_sepp_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2SEPP, "SEPP");
    TH1D *hs_mass_semm_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2SEMM, "SEMM");

    // Same-event profiles: v2
    TH1D *hs_v2_sepm_proj =
        helper.GetV2(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                     cent_min, cent_max, hs_V2SEPM, R2SP, "SEPM");
    TH1D *hs_v2_sepp_proj =
        helper.GetV2(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                     cent_min, cent_max, hs_V2SEPP, R2SP, "SEPP");
    TH1D *hs_v2_semm_proj =
        helper.GetV2(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                     cent_min, cent_max, hs_V2SEMM, R2SP, "SEMM");

    // Mixed-event profiles: mass
    TH1D *hs_mass_mepm_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2MEPM, "MEPM");
    TH1D *hs_mass_mepp_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2MEPP, "MEPP");
    TH1D *hs_mass_memm_proj =
        helper.GetMass(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
                       cent_min, cent_max, hs_V2MEMM, "MEMM");

    // Mixed-event profiles: v2
    TH1D *hs_v2_mepm_proj = helper.GetV2EM(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, hs_u2q2_cosDeltaPhi_MEPM1, hs_u2q2_cosDeltaPhi_MEPM2,
        hs_r2spABMEPM1, hs_r2spACMEPM1, hs_r2spBCMEPM1, hs_r2spABMEPM2,
        hs_r2spACMEPM2, hs_r2spBCMEPM2, "MEPM");
    TH1D *hs_v2_mepp_proj = helper.GetV2EM(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, hs_u2q2_cosDeltaPhi_MEPP1, hs_u2q2_cosDeltaPhi_MEPP2,
        hs_r2spABMEPP1, hs_r2spACMEPP1, hs_r2spBCMEPP1, hs_r2spABMEPP2,
        hs_r2spACMEPP2, hs_r2spBCMEPP2, "MEPP");
    TH1D *hs_v2_memm_proj = helper.GetV2EM(
        Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min,
        cent_max, hs_u2q2_cosDeltaPhi_MEMM1, hs_u2q2_cosDeltaPhi_MEMM2,
        hs_r2spABMEMM1, hs_r2spACMEMM1, hs_r2spBCMEMM1, hs_r2spABMEMM2,
        hs_r2spACMEMM2, hs_r2spBCMEMM2, "MEMM");

    // Scale mixed-event spectra with F factor
    hs_mass_mepm_proj->Scale(ffactor[i]);
    hs_mass_mepp_proj->Scale(ffactor[i]);
    hs_mass_memm_proj->Scale(ffactor[i]);
    hs_v2_mepm_proj->Scale(ffactor[i]);
    hs_v2_mepp_proj->Scale(ffactor[i]);
    hs_v2_memm_proj->Scale(ffactor[i]);

    // Save plots for invariant mass
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