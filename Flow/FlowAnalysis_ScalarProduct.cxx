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

//______________________________________________________________________________
void FlowAnalysis_ScalarProduct(int flag_sig, int flag_bkg,
                                std::string FileName = "AnalysisResults.root",
                                double cent_min = 10., double cent_max = 50.,
                                double chi2max = 2., bool sys = false,
                                std::string muonCut = "muonLowPt210SigmaPDCA",
                                std::string dimuonCut = "pairRapidityForward") {

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
  fitter.setChi2Max(chi2max);

  // Define variables' range for analysis
  double Bin_pt_mass[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
  double mass_min = 2.3;
  double mass_max = 4.3;

  // Create output file
  TFile f(
      Form("FlowAnalysisResults_ScalarProduct_%g_%g.root", cent_min, cent_max),
      "RECREATE");

  ///////////////////////////////////////////////////
  ///                                             ///
  ///   Analysis for Differential Flow of J/Psi   ///
  ///                                             ///
  ///////////////////////////////////////////////////

  // Resolution factors
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
  double x_yield[10], y_yield[10], ex_yield[10], ey_yield[10];
  double x_v2pt[10], y_v2pt[10], ex_v2pt[10], ey_v2pt[10];
  double *x_run2_1, *y_run2_1, *ex_run2_1, *ey_run2_1, *eysys_run2_1;
  double *x_run2_2, *y_run2_2, *ex_run2_2, *ey_run2_2, *eysys_run2_2;

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

    TList *l_diff = new TList();

    // Copy original profiles for projections
    THnSparse *hs_V2_cp = dynamic_cast<THnSparse *>(
        hs_V2->Clone(Form("Mass_Pt_centrFT0C_V2_Copy_%g_%g", Bin_pt_mass[i],
                          Bin_pt_mass[i + 1])));

    // Set axes' ranges for mass-differential study
    hs_V2_cp->GetAxis(0)->SetRangeUser(mass_min, mass_max);
    hs_V2_cp->GetAxis(1)->SetRangeUser(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
    hs_V2_cp->GetAxis(3)->SetRangeUser(cent_min, cent_max);

    // Get mass and v2
    TH1D *hs_mass_proj = hs_V2_cp->Projection(0);
    TH2D *hs_v2_sp_proj = hs_V2_cp->Projection(4, 0);

    // Define histograms
    double *Bin_mass_new = CreateBinsFromAxis(hs_mass_proj->GetXaxis());
    int NBins_mass_new = hs_mass_proj->GetXaxis()->GetNbins();
    TH1D *hist_v2sp = new TH1D(
        Form("v2sp_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        Form("v^{#mu#mu}_{2}{SP}_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        NBins_mass_new, Bin_mass_new);
    hist_v2sp->GetXaxis()->SetTitle("mass (GeV/c2)");
    hist_v2sp->GetYaxis()->SetTitle("v^{#mu#mu}_{2}{SP}");

    // Evaluation of differential flow as function of invariant mass
    for (int j = 0; j < NBins_mass_new; j++) {
      TH2D *hs_v2_sp_proj_cp =
          dynamic_cast<TH2D *>(hs_v2_sp_proj->Clone("Mass_V2SP_Copy"));
      hs_v2_sp_proj_cp->GetXaxis()->SetRangeUser(Bin_mass_new[j],
                                                 Bin_mass_new[j + 1]);
      double v2sp = hs_v2_sp_proj_cp->GetMean(2);
      double v2spe = hs_v2_sp_proj_cp->GetMean(12);

      hist_v2sp->SetBinContent(j + 1, v2sp / R2SP);
      hist_v2sp->SetBinError(j + 1, v2spe / R2SP);

      delete hs_v2_sp_proj_cp;
    }

    /// Do fitting

    // Configuration for fitting
    fitter.setModel(flag_sig, flag_bkg);
    fitter.setMassRange(mass_min, mass_max);
    fitter.setCentRange(cent_min, cent_max);

    // Fit invariant mass + v2
    TList *l_diff_fit = new TList();
    vector<double> results_v2 =
        fitter.runFitting(hs_mass_proj, hist_v2sp, l_diff_fit, Bin_pt_mass[i],
                          Bin_pt_mass[i + 1]);

    // Fill pT-differential v2 and jpsi yields
    x_v2pt[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    y_v2pt[i] = results_v2[0];
    ex_v2pt[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;
    ey_v2pt[i] = results_v2[1];

    x_yield[i] = (Bin_pt_mass[i] + Bin_pt_mass[i + 1]) / 2;
    y_yield[i] = results_v2[2] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);
    ex_yield[i] = (Bin_pt_mass[i + 1] - Bin_pt_mass[i]) / 2;
    ey_yield[i] = results_v2[3] / (Bin_pt_mass[i + 1] - Bin_pt_mass[i]);

    // Save results
    TList *l_diff_hist = new TList();
    l_diff_hist->Add(hist_v2sp);
    f.cd();
    l_diff_hist->Write(
        Form("DifferentialFlow_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    l_diff_fit->Write(
        Form("DifferentialFlow_Fit_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);

    delete hs_mass_proj;
    delete hs_v2_sp_proj;
    delete hist_v2sp;
    delete hs_V2_cp;
    delete l_diff_hist;
    delete l_diff_fit;
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

  // Compare with Run2 data: 5.02 TeV, 10-30% and 30-50%
  TGraphMultiErrors *graph_v2pt = new TGraphMultiErrors(
      "v2_pt", "", 10, x_v2pt, y_v2pt, ex_v2pt, ex_v2pt, ey_v2pt, ey_v2pt);
  graph_v2pt->SetTitle("#sqrt{#it{s}_{NN}} = 5.36 TeV, 10-50%");
  graph_v2pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

  TGraphMultiErrors *graph_yield =
      new TGraphMultiErrors("yields_pt", "", 10, x_yield, y_yield, ex_yield,
                            ex_yield, ey_yield, ey_yield);
  graph_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_yield->GetYaxis()->SetTitle("Raw counts / d#it{p}_{T} (GeV/c)^{-1}");

  // Save final results
  TCanvas *c_yield = new TCanvas("jpsi_yield_pT", "jpsi_yield_pT");
  c_yield->cd();
  graph_yield->SetMarkerStyle(20);
  graph_yield->SetMarkerSize(1.);
  graph_yield->SetMarkerColor(kBlue);
  graph_yield->SetLineColor(kBlue);
  graph_yield->SetLineWidth(2);
  graph_yield->SetFillStyle(0);
  graph_yield->Draw("A P Z");
  // gPad->ModifiedUpdate();
  TLatex *text_yield = new TLatex();
  text_yield->SetTextSize(0.04);
  text_yield->SetTextFont(42);
  text_yield->DrawLatexNDC(
      .3, .82, "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  gPad->ModifiedUpdate();
  text_yield->DrawLatexNDC(.3, .77,
                           "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");
  gPad->ModifiedUpdate();
  l_results->Add(c_yield);

  TCanvas *c_pt = new TCanvas("v2_pT", "v2_pT");
  c_pt->cd();
  auto mg = new TMultiGraph();
  graph_v2pt->SetMarkerStyle(20);
  graph_v2pt->SetMarkerSize(1.);
  graph_v2pt->SetMarkerColor(kBlue);
  graph_v2pt->SetLineColor(kBlue);
  graph_v2pt->SetLineWidth(2);
  graph_v2pt->SetFillStyle(0);
  graph_v2pt_run2_1->SetMarkerStyle(20);
  graph_v2pt_run2_1->SetMarkerSize(1.);
  graph_v2pt_run2_1->SetMarkerColor(kRed);
  graph_v2pt_run2_1->SetLineColor(kRed);
  graph_v2pt_run2_1->SetLineWidth(2);
  graph_v2pt_run2_1->SetFillStyle(0);
  graph_v2pt_run2_2->SetMarkerStyle(20);
  graph_v2pt_run2_2->SetMarkerSize(1.);
  graph_v2pt_run2_2->SetMarkerColor(kOrange);
  graph_v2pt_run2_2->SetLineColor(kOrange);
  graph_v2pt_run2_2->SetLineWidth(2);
  graph_v2pt_run2_2->SetFillStyle(0);
  mg->Add(graph_v2pt);
  mg->Add(graph_v2pt_run2_1);
  mg->Add(graph_v2pt_run2_2);
  mg->Draw("A P Z ; Z ; 5 s=0.5");
  gPad->BuildLegend();
  // gPad->ModifiedUpdate();
  TLatex *text_pt = new TLatex();
  text_pt->SetTextSize(0.04);
  text_pt->SetTextFont(42);
  text_pt->DrawLatexNDC(
      .18, .82, "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  gPad->ModifiedUpdate();
  text_pt->DrawLatexNDC(.18, .77,
                        "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");
  gPad->ModifiedUpdate();
  l_results->Add(c_pt);

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
        Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
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
            Form("PairsMuonSEPM_%s", muonCut.c_str()));
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