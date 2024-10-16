#include "FlowAnalysis_Helper.h"

//______________________________________________________________________________
double *FlowAnalysis_Helper::CreateBinsFromAxis(TAxis *axis) {
  int Nbins = axis->GetNbins();
  double *Bins = new double[Nbins + 1];
  axis->GetLowEdge(Bins);
  Bins[Nbins] = axis->GetBinUpEdge(Nbins);
  return Bins;
}

//______________________________________________________________________________
void FlowAnalysis_Helper::CreateBins(double *axis, double min, double max,
                                     int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Helper::GetMass(double ptmin, double ptmax, double massmin,
                                   double massmax, double centmin,
                                   double centmax, THnSparse *hist_V2) {

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
TH1D *FlowAnalysis_Helper::GetRfactor(double ptmin, double ptmax,
                                      double massmin, double massmax,
                                      double centmin, double centmax,
                                      THnSparse *hs_V2MEPM,
                                      THnSparse *hs_V2MEPP,
                                      THnSparse *hs_V2MEMM) {
  // Copy original profiles for projections
  THnSparse *hs_V2MEPM_cp = dynamic_cast<THnSparse *>(
      hs_V2MEPM->Clone(Form("MEPM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  THnSparse *hs_V2MEPP_cp = dynamic_cast<THnSparse *>(
      hs_V2MEPP->Clone(Form("MEPP_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  THnSparse *hs_V2MEMM_cp = dynamic_cast<THnSparse *>(
      hs_V2MEMM->Clone(Form("MEMM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));

  // Set axes' ranges for mass-differential study
  hs_V2MEPM_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_V2MEPM_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_V2MEPM_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  hs_V2MEPP_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_V2MEPP_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_V2MEPP_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  hs_V2MEMM_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_V2MEMM_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_V2MEMM_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  TH1D *hs_V2MEPM_cp_proj = hs_V2MEPM_cp->Projection(0);
  TH1D *hs_V2MEPP_cp_proj = hs_V2MEPP_cp->Projection(0);
  TH1D *hs_V2MEMM_cp_proj = hs_V2MEMM_cp->Projection(0);

  // Define resulting histogram
  double *Bin_mass_new = CreateBinsFromAxis(hs_V2MEPM_cp_proj->GetXaxis());
  int NBins_mass_new = hs_V2MEPM_cp_proj->GetXaxis()->GetNbins();
  TH1D *hist_rfactor = new TH1D(Form("Rfactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                Form("Rfactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                NBins_mass_new, Bin_mass_new);
  for (int i = 0; i < NBins_mass_new; i++) {
    double Npm = hs_V2MEPM_cp_proj->GetBinContent(i + 1);
    double Npp = hs_V2MEPP_cp_proj->GetBinContent(i + 1);
    double Nmm = hs_V2MEMM_cp_proj->GetBinContent(i + 1);

    double val = Npp * Nmm <= 0. ? 0. : 0.5 * Npm / pow(Npp * Nmm, 0.5);
    hist_rfactor->SetBinContent(i + 1, val);
  }

  delete hs_V2MEPM_cp;
  delete hs_V2MEPP_cp;
  delete hs_V2MEMM_cp;
  delete hs_V2MEPM_cp_proj;
  delete hs_V2MEPP_cp_proj;
  delete hs_V2MEMM_cp_proj;

  return hist_rfactor;
}

//______________________________________________________________________________
double
FlowAnalysis_Helper::GetFfactor(double ptmin, double ptmax, double massmin,
                                double massmax, double centmin, double centmax,
                                THnSparse *hs_V2SEPP, THnSparse *hs_V2SEMM,
                                THnSparse *hs_V2MEPM, TH1D *hist_rfactor) {
  // Copy original profiles for projections
  THnSparse *hs_V2SEPP_cp = dynamic_cast<THnSparse *>(
      hs_V2SEPP->Clone(Form("SEPP_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  THnSparse *hs_V2SEMM_cp = dynamic_cast<THnSparse *>(
      hs_V2SEMM->Clone(Form("SEMM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  THnSparse *hs_V2MEPM_cp = dynamic_cast<THnSparse *>(
      hs_V2MEPM->Clone(Form("MEPM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));

  // Set axes' ranges for mass-differential study
  hs_V2SEPP_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_V2SEPP_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_V2SEPP_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  hs_V2SEMM_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_V2SEMM_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_V2SEMM_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  hs_V2MEPM_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_V2MEPM_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_V2MEPM_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  TH1D *hs_V2SEPP_cp_proj = hs_V2SEPP_cp->Projection(0);
  TH1D *hs_V2SEMM_cp_proj = hs_V2SEMM_cp->Projection(0);
  TH1D *hs_V2MEPM_cp_proj = hs_V2MEPM_cp->Projection(0);

  double *Bin_mass_new = CreateBinsFromAxis(hs_V2MEPM_cp_proj->GetXaxis());
  int NBins_mass_new = hs_V2MEPM_cp_proj->GetXaxis()->GetNbins();
  TH1D *hist_ffactor = new TH1D(Form("Ffactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                Form("Ffactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                NBins_mass_new, Bin_mass_new);
  for (int i = 0; i < NBins_mass_new; i++) {
    double R_val = hist_rfactor->GetBinContent(i + 1);
    double N_SEPP = hs_V2SEPP_cp_proj->GetBinContent(i + 1);
    double N_SEMM = hs_V2SEMM_cp_proj->GetBinContent(i + 1);

    double F_val =
        N_SEPP * N_SEMM < 0. ? 0. : 2. * R_val * pow(N_SEPP * N_SEMM, 0.5);
    hist_ffactor->SetBinContent(i + 1, F_val);
  }

  return hist_ffactor->Integral("width") / hs_V2MEPM_cp_proj->Integral("width");
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Helper::GetV2(double ptmin, double ptmax, double massmin,
                                 double massmax, double centmin, double centmax,
                                 THnSparse *hist_V2, double R2SP) {

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
vector<double> FlowAnalysis_Helper::GetStats(int size, double *sample,
                                             double *sample_error) {
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
void FlowAnalysis_Helper::PlotSystematics(
    int index, TCanvas *c_sys_yield, TCanvas *c_sys_v2, TH1D *hist_sys_yield,
    TH1D *hist_sys_v2, double *bins_sys_yield, double *bins_sys_v2,
    double *chi2_yield, double *chi2_v2, int nbCombo_yield, int nbCombo_v2,
    vector<double> stats_yield, vector<double> stats_v2, double *pt_bins,
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
void FlowAnalysis_Helper::PlotSNRvsRun2(int size_ptbin, double *pt_bins,
                                        int size_run3, double *x_run3,
                                        double *snr_run3, int size_run2,
                                        double *x_run2, double *snr_run2,
                                        TList *ls) {

  // Saving plot for J/psi yields SNR as function of pT
  // compared with Run2
  TGraph *graph_SNR_Run3 = new TGraph(size_run3, x_run3, snr_run3);
  graph_SNR_Run3->SetNameTitle("graph_snr_run3", "Run3");
  TGraph *graph_SNR_Run2 = new TGraph(size_run2, x_run2, snr_run2);
  graph_SNR_Run2->SetNameTitle("graph_snr_run2", "Run2");
  TCanvas *c_SNR = new TCanvas("jpsi_SNR_pT", "jpsi_SNR_pT");
  c_SNR->cd();
  if (!(size_run2 == size_run3)) {
    TPad *pad_snr = new TPad("pad_snr", "pad_snr", 0, 0., 1, 1.0);
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
void FlowAnalysis_Helper::PlotFinalResults(
    int size_ptbin, double *pt_bins, double *x_v2pt, double *y_v2pt,
    double *ex_v2pt, double *ey_v2pt, double *eysys_v2pt, double *x_run2,
    double *y_run2, double *ex_run2, double *ey_run2, double *eysys_run2,
    double *x_yield, double *y_yield, double *ex_yield, double *ey_yield,
    double *eysys_yield, double *x_yield_run2, double *y_yield_run2,
    double *ex_yield_run2, double *ey_yield_run2, double *eysys_yield_run2,
    TList *ls) {
  // Compare with Run2 data: 5.02 TeV
  TCanvas *c_yield = new TCanvas("jpsi_yield_pT", "jpsi_yield_pT");
  TCanvas *c_yield_vsRun2 =
      new TCanvas("jpsi_yield_raw_pT", "jpsi_yield_raw_pT");
  TCanvas *c_pt = new TCanvas("v2_pT", "v2_pT");

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
  c_yield->cd();
  TPad *pad_yield_final =
      new TPad("pad_yield_final", "pad_yield_final", 0, 0, 1, 1);
  pad_yield_final->Draw();
  pad_yield_final->cd();
  graph_yield->SetMarkerStyle(20);
  graph_yield->SetMarkerSize(1.);
  graph_yield->SetMarkerColor(kBlue);
  graph_yield->SetLineColor(kBlue);
  graph_yield->SetLineWidth(2);
  graph_yield->SetFillStyle(0);
  graph_yield->GetXaxis()->SetLimits(pt_bins[0], pt_bins[size_ptbin]);
  graph_yield->Draw("A P Z ; Z ; 5 s=0.5");
  TLatex *text_yield = new TLatex();
  text_yield->SetTextSize(0.04);
  text_yield->SetTextFont(42);
  text_yield->DrawLatexNDC(.3, .82,
                           "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                           "= 5.36 TeV");
  pad_yield_final->ModifiedUpdate();
  text_yield->DrawLatexNDC(.3, .77,
                           "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                           "4");
  pad_yield_final->ModifiedUpdate();
  ls->Add(c_yield);

  // Saving plot for J/psi yields as function of pT
  // compared with Run2
  if (size_ptbin == 15) {
    c_yield_vsRun2->cd();
    c_yield_vsRun2->SetBottomMargin(0);
    TPad *pad_yield_final_run2 = new TPad(
        "pad_yield_final_run2", "pad_yield_final_run2", 0, 0.3, 1, 1.0);
    pad_yield_final_run2->SetBottomMargin(0);
    pad_yield_final_run2->Draw();
    pad_yield_final_run2->cd();
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
    mg_raw->GetXaxis()->SetLimits(pt_bins[0], pt_bins[size_ptbin]);
    mg_raw->Draw("A P Z ; Z ; 5 s=0.5");
    pad_yield_final_run2->BuildLegend();

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
    mg_raw->GetXaxis()->SetLimits(0., 20.);
    mg_raw->Draw("A P Z ; Z ; 5 s=0.5");
  }
  ls->Add(c_yield_vsRun2);

  // Saving plot for J/psi v2 as function of pT and
  // compared with Run2
  c_pt->cd();
  TPad *pad_pt_final = new TPad("pad_pt_final", "pad_pt_final", 0, 0, 1, 1);
  pad_pt_final->Draw();
  pad_pt_final->cd();
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
  mg->GetXaxis()->SetLimits(pt_bins[0], pt_bins[size_ptbin]);
  mg->Draw("A P Z ; Z ; 5 s=0.5");
  pad_pt_final->BuildLegend();
  TLatex *text_pt = new TLatex();
  text_pt->SetTextSize(0.04);
  text_pt->SetTextFont(42);
  text_pt->DrawLatexNDC(.18, .82,
                        "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                        "= 5.36 TeV");
  pad_pt_final->ModifiedUpdate();
  text_pt->DrawLatexNDC(.18, .77,
                        "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                        "4");
  pad_pt_final->ModifiedUpdate();
  TF1 *lv2_pt = new TF1("lv2_pt", "[0]", pt_bins[0], pt_bins[size_ptbin]);
  lv2_pt->SetParameter(0, 0.);
  lv2_pt->SetLineColor(18);
  lv2_pt->SetLineWidth(3);
  lv2_pt->SetLineStyle(9);
  lv2_pt->Draw("same");
  pad_pt_final->ModifiedUpdate();
  ls->Add(c_pt);
}

//______________________________________________________________________________
void FlowAnalysis_Helper::LoadData(std::string FileName, THnSparse *&hs_V2,
                                   TH2F *&hs_R2SPAB, TH2F *&hs_R2SPAC,
                                   TH2F *&hs_R2SPBC, std::string muonCut,
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
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C"));
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
void FlowAnalysis_Helper::LoadDataME(
    std::string FileName, THnSparse *&hs_V2SEPM, THnSparse *&hs_V2SEPP,
    THnSparse *&hs_V2SEMM, THnSparse *&hs_V2MEPM, THnSparse *&hs_V2MEPP,
    THnSparse *&hs_V2MEMM, TH2F *&hs_R2SPAB, TH2F *&hs_R2SPAC, TH2F *&hs_R2SPBC,
    TH3F *&hs_u2q2MEPM1, TH3F *&hs_u2q2MEPP1, TH3F *&hs_u2q2MEMM1,
    TH3F *&hs_u2q2MEPM2, TH3F *&hs_u2q2MEPP2, TH3F *&hs_u2q2MEMM2,
    TH3F *&hs_r2spMEPM1, TH3F *&hs_r2spMEPP1, TH3F *&hs_r2spMEMM1,
    TH3F *&hs_r2spMEPM2, TH3F *&hs_r2spMEPP2, TH3F *&hs_r2spMEMM2,
    TH2F *&hs_cosDeltaPhiMEPM1, TH2F *&hs_cosDeltaPhiMEPP1,
    TH2F *&hs_cosDeltaPhiMEMM1, TH2F *&hs_cosDeltaPhiMEPM2,
    TH2F *&hs_cosDeltaPhiMEPP2, TH2F *&hs_cosDeltaPhiMEMM2, std::string muonCut,
    std::string dimuonCut) {
  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist_v2se, *list_hist_v2me, *list_hist_r2;
  TList *sublist_v2sepm, *sublist_v2sepp, *sublist_v2semm, *sublist_r2;
  TList *sublist_v2mepm, *sublist_v2mepp, *sublist_v2memm;
  if (filePath.extension() == ".root") {
    // Load data from AnalysisResults.root
    TFile *Input_File = TFile::Open(FileName.c_str());
    list_hist_v2se =
        (THashList *)Input_File->Get("analysis-same-event-pairing/output");
    list_hist_v2me =
        (THashList *)Input_File->Get("analysis-event-mixing/output");
    sublist_v2sepm = (TList *)list_hist_v2se->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonSEPM_%s", muonCut.c_str())
            : Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    sublist_v2sepp = (TList *)list_hist_v2se->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonSEPP_%s", muonCut.c_str())
            : Form("PairsMuonSEPP_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    sublist_v2semm = (TList *)list_hist_v2se->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonSEMM_%s", muonCut.c_str())
            : Form("PairsMuonSEMM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    sublist_v2mepm = (TList *)list_hist_v2me->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonMEPM_%s", muonCut.c_str())
            : Form("PairsMuonMEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    sublist_v2mepp = (TList *)list_hist_v2me->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonMEPP_%s", muonCut.c_str())
            : Form("PairsMuonMEPP_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    sublist_v2memm = (TList *)list_hist_v2me->FindObject(
        dimuonCut == ""
            ? Form("PairsMuonMEMM_%s", muonCut.c_str())
            : Form("PairsMuonMEMM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
    list_hist_r2 =
        (THashList *)Input_File->Get("analysis-event-selection/output");
    sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

    // Get histograms of correlations and resolution factors
    hs_V2SEPM = (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2");
    hs_V2SEPP = (THnSparse *)sublist_v2sepp->FindObject("Mass_Pt_centrFT0C_V2");
    hs_V2SEMM = (THnSparse *)sublist_v2semm->FindObject("Mass_Pt_centrFT0C_V2");
    hs_V2MEPM = (THnSparse *)sublist_v2mepm->FindObject("Mass_Pt_centrFT0C_V2");
    hs_V2MEPP = (THnSparse *)sublist_v2mepp->FindObject("Mass_Pt_centrFT0C_V2");
    hs_V2MEMM = (THnSparse *)sublist_v2memm->FindObject("Mass_Pt_centrFT0C_V2");

    hs_R2SPAB = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C");
    hs_R2SPAC = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C");
    hs_R2SPBC = (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C");

    hs_u2q2MEPM1 = (TH3F *)sublist_v2mepm->FindObject("U2Q2_CentFT0C_ev1");
    hs_u2q2MEPM2 = (TH3F *)sublist_v2mepm->FindObject("U2Q2_CentFT0C_ev2");
    hs_u2q2MEPP1 = (TH3F *)sublist_v2mepp->FindObject("U2Q2_CentFT0C_ev1");
    hs_u2q2MEPP2 = (TH3F *)sublist_v2mepp->FindObject("U2Q2_CentFT0C_ev2");
    hs_u2q2MEMM1 = (TH3F *)sublist_v2memm->FindObject("U2Q2_CentFT0C_ev1");
    hs_u2q2MEMM2 = (TH3F *)sublist_v2memm->FindObject("U2Q2_CentFT0C_ev2");
    hs_r2spMEPM1 = (TH3F *)sublist_v2mepm->FindObject("R2SP1_CentFT0C");
    hs_r2spMEPM2 = (TH3F *)sublist_v2mepm->FindObject("R2SP2_CentFT0C");
    hs_r2spMEPP1 = (TH3F *)sublist_v2mepp->FindObject("R2SP1_CentFT0C");
    hs_r2spMEPP2 = (TH3F *)sublist_v2mepp->FindObject("R2SP2_CentFT0C");
    hs_r2spMEMM1 = (TH3F *)sublist_v2memm->FindObject("R2SP1_CentFT0C");
    hs_r2spMEMM2 = (TH3F *)sublist_v2memm->FindObject("R2SP2_CentFT0C");

    hs_cosDeltaPhiMEPM1 =
        (TH2F *)sublist_v2mepm->FindObject("Mass_cos2DeltaPhiMu1");
    hs_cosDeltaPhiMEPM2 =
        (TH2F *)sublist_v2mepm->FindObject("Mass_cos2DeltaPhiMu2");
    hs_cosDeltaPhiMEPP1 =
        (TH2F *)sublist_v2mepp->FindObject("Mass_cos2DeltaPhiMu1");
    hs_cosDeltaPhiMEPP2 =
        (TH2F *)sublist_v2mepp->FindObject("Mass_cos2DeltaPhiMu2");
    hs_cosDeltaPhiMEMM1 =
        (TH2F *)sublist_v2memm->FindObject("Mass_cos2DeltaPhiMu1");
    hs_cosDeltaPhiMEMM2 =
        (TH2F *)sublist_v2memm->FindObject("Mass_cos2DeltaPhiMu2");
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
        cout << "Flag1" << endl;
        list_hist_v2se =
            (THashList *)inFile->Get("analysis-same-event-pairing/output");
        cout << "Flag1" << endl;
        list_hist_v2me =
            (THashList *)inFile->Get("analysis-event-mixing/output");
        cout << "Flag1" << endl;
        sublist_v2sepm = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPM_%s", muonCut.c_str())
                            : Form("PairsMuonSEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        cout << "Flag1" << endl;
        sublist_v2sepp = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPP_%s", muonCut.c_str())
                            : Form("PairsMuonSEPP_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        cout << "Flag1" << endl;
        sublist_v2semm = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEMM_%s", muonCut.c_str())
                            : Form("PairsMuonSEMM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        cout << "Flag1" << endl;
        sublist_v2mepm = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEPM_%s", muonCut.c_str())
                            : Form("PairsMuonMEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        cout << "Flag1" << endl;
        sublist_v2mepp = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEPP_%s", muonCut.c_str())
                            : Form("PairsMuonMEPP_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        cout << "Flag1" << endl;
        sublist_v2memm = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEMM_%s", muonCut.c_str())
                            : Form("PairsMuonMEMM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        cout << "Flag1" << endl;
        list_hist_r2 =
            (THashList *)inFile->Get("analysis-event-selection/output");
        sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");
        cout << "Flag1" << endl;

        if (first) {
          hs_V2SEPM =
              (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2");
          cout << "Flag1" << endl;
          hs_V2SEPP =
              (THnSparse *)sublist_v2sepp->FindObject("Mass_Pt_centrFT0C_V2");
          cout << "Flag1" << endl;
          hs_V2SEMM =
              (THnSparse *)sublist_v2semm->FindObject("Mass_Pt_centrFT0C_V2");
          cout << "Flag1" << endl;
          hs_V2MEPM =
              (THnSparse *)sublist_v2mepm->FindObject("Mass_Pt_centrFT0C_V2");
          cout << "Flag1" << endl;
          hs_V2MEPP =
              (THnSparse *)sublist_v2mepp->FindObject("Mass_Pt_centrFT0C_V2");
          cout << "Flag1" << endl;
          hs_V2MEMM =
              (THnSparse *)sublist_v2memm->FindObject("Mass_Pt_centrFT0C_V2");
          cout << "Flag1" << endl;

          hs_R2SPAB = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C");
          hs_R2SPAC = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C");
          hs_R2SPBC = (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C");
          cout << "Flag1" << endl;

          hs_u2q2MEPM1 =
              (TH3F *)sublist_v2mepm->FindObject("U2Q2_CentFT0C_ev1");
          hs_u2q2MEPM2 =
              (TH3F *)sublist_v2mepm->FindObject("U2Q2_CentFT0C_ev2");
          hs_u2q2MEPP1 =
              (TH3F *)sublist_v2mepp->FindObject("U2Q2_CentFT0C_ev1");
          hs_u2q2MEPP2 =
              (TH3F *)sublist_v2mepp->FindObject("U2Q2_CentFT0C_ev2");
          hs_u2q2MEMM1 =
              (TH3F *)sublist_v2memm->FindObject("U2Q2_CentFT0C_ev1");
          hs_u2q2MEMM2 =
              (TH3F *)sublist_v2memm->FindObject("U2Q2_CentFT0C_ev2");
          cout << "Flag1" << endl;

          hs_r2spMEPM1 = (TH3F *)sublist_v2mepm->FindObject("R2SP1_CentFT0C");
          hs_r2spMEPM2 = (TH3F *)sublist_v2mepm->FindObject("R2SP2_CentFT0C");
          hs_r2spMEPP1 = (TH3F *)sublist_v2mepp->FindObject("R2SP1_CentFT0C");
          hs_r2spMEPP2 = (TH3F *)sublist_v2mepp->FindObject("R2SP2_CentFT0C");
          hs_r2spMEMM1 = (TH3F *)sublist_v2memm->FindObject("R2SP1_CentFT0C");
          hs_r2spMEMM2 = (TH3F *)sublist_v2memm->FindObject("R2SP2_CentFT0C");
          cout << "Flag1" << endl;

          hs_cosDeltaPhiMEPM1 =
              (TH2F *)sublist_v2mepm->FindObject("Mass_cos2DeltaPhiMu1");
          hs_cosDeltaPhiMEPM2 =
              (TH2F *)sublist_v2mepm->FindObject("Mass_cos2DeltaPhiMu2");
          hs_cosDeltaPhiMEPP1 =
              (TH2F *)sublist_v2mepp->FindObject("Mass_cos2DeltaPhiMu1");
          hs_cosDeltaPhiMEPP2 =
              (TH2F *)sublist_v2mepp->FindObject("Mass_cos2DeltaPhiMu2");
          hs_cosDeltaPhiMEMM1 =
              (TH2F *)sublist_v2memm->FindObject("Mass_cos2DeltaPhiMu1");
          hs_cosDeltaPhiMEMM2 =
              (TH2F *)sublist_v2memm->FindObject("Mass_cos2DeltaPhiMu2");
          cout << "Flag1" << endl;
        } else {
          hs_V2SEPM->Add(
              (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2SEPP->Add(
              (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2SEMM->Add(
              (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2MEPM->Add(
              (THnSparse *)sublist_v2mepm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2MEPP->Add(
              (THnSparse *)sublist_v2mepp->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2MEMM->Add(
              (THnSparse *)sublist_v2memm->FindObject("Mass_Pt_centrFT0C_V2"));
          cout << "Flag1" << endl;

          hs_R2SPAB->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C"));
          hs_R2SPAC->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C"));
          hs_R2SPBC->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C"));
          cout << "Flag1" << endl;

          hs_u2q2MEPM1->Add(
              (TH3F *)sublist_v2mepm->FindObject("U2Q2_CentFT0C_ev1"));
          hs_u2q2MEPM2->Add(
              (TH3F *)sublist_v2mepm->FindObject("U2Q2_CentFT0C_ev2"));
          hs_u2q2MEPP1->Add(
              (TH3F *)sublist_v2mepp->FindObject("U2Q2_CentFT0C_ev1"));
          hs_u2q2MEPP2->Add(
              (TH3F *)sublist_v2mepp->FindObject("U2Q2_CentFT0C_ev2"));
          hs_u2q2MEMM1->Add(
              (TH3F *)sublist_v2memm->FindObject("U2Q2_CentFT0C_ev1"));
          hs_u2q2MEMM2->Add(
              (TH3F *)sublist_v2memm->FindObject("U2Q2_CentFT0C_ev2"));
          cout << "Flag1" << endl;

          hs_r2spMEPM1->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SP1_CentFT0C"));
          hs_r2spMEPM2->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SP2_CentFT0C"));
          hs_r2spMEPP1->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SP1_CentFT0C"));
          hs_r2spMEPP2->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SP2_CentFT0C"));
          hs_r2spMEMM1->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SP1_CentFT0C"));
          hs_r2spMEMM2->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SP2_CentFT0C"));
          cout << "Flag1" << endl;

          hs_cosDeltaPhiMEPM1->Add(
              (TH2F *)sublist_v2mepm->FindObject("Mass_cos2DeltaPhiMu1"));
          hs_cosDeltaPhiMEPM2->Add(
              (TH2F *)sublist_v2mepm->FindObject("Mass_cos2DeltaPhiMu2"));
          hs_cosDeltaPhiMEPP1->Add(
              (TH2F *)sublist_v2mepp->FindObject("Mass_cos2DeltaPhiMu1"));
          hs_cosDeltaPhiMEPP2->Add(
              (TH2F *)sublist_v2mepp->FindObject("Mass_cos2DeltaPhiMu2"));
          hs_cosDeltaPhiMEMM1->Add(
              (TH2F *)sublist_v2memm->FindObject("Mass_cos2DeltaPhiMu1"));
          hs_cosDeltaPhiMEMM2->Add(
              (TH2F *)sublist_v2memm->FindObject("Mass_cos2DeltaPhiMu2"));
          cout << "Flag1" << endl;
        }

        inFile->Close();
        first = false;
      }
      InputFiles.close();
    }
  }
}

//______________________________________________________________________________
void FlowAnalysis_Helper::LoadDataRun2(double *&x, double *&y, double *&ex,
                                       double *&ey, double *&ey_sys, int flag) {
  x = new double[10];
  y = new double[10];
  ex = new double[10];
  ey = new double[10];
  ey_sys = new double[10];
  if (flag == 0) {
    // 0-10%
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

    y[0] = 0.03;
    y[1] = 0.034;
    y[2] = 0.022;
    y[3] = 0.053;
    y[4] = 0.043;
    y[5] = 0.05;
    y[6] = 0.045;
    y[7] = 0.0006;
    y[8] = 0.068;
    y[9] = 0.002;

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

    ey[0] = 0.013;
    ey[1] = 0.01;
    ey[2] = 0.011;
    ey[3] = 0.013;
    ey[4] = 0.016;
    ey[5] = 0.019;
    ey[6] = 0.019;
    ey[7] = 0.028;
    ey[8] = 0.046;
    ey[9] = 0.059;

    ey_sys[0] = 0.0041243;
    ey_sys[1] = 0.0031563;
    ey_sys[2] = 0.0040117;
    ey_sys[3] = 0.0030156;
    ey_sys[4] = 0.0017944;
    ey_sys[5] = 0.0025338;
    ey_sys[6] = 0.0033808;
    ey_sys[7] = 0.0051196;
    ey_sys[8] = 0.0059211;
    ey_sys[9] = 0.0051197;
  } else if (flag == 1) {
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

//______________________________________________________________________________
void FlowAnalysis_Helper::LoadDataYieldRun2(double *&x, double *&y, double *&ex,
                                            double *&ey, double *&ey_sys,
                                            double *&SNR, int flag) {
  x = new double[15];
  y = new double[15];
  ex = new double[15];
  ey = new double[15];
  ey_sys = new double[15];
  SNR = new double[15];

  if (flag == 0) {
    // 0-20%
    x[0] = 0.15;
    x[1] = 0.65;
    x[2] = 1.5;
    x[3] = 2.5;
    x[4] = 3.5;
    x[5] = 4.5;
    x[6] = 5.5;
    x[7] = 6.5;
    x[8] = 7.5;
    x[9] = 8.5;
    x[10] = 9.5;
    x[11] = 10.5;
    x[12] = 11.5;
    x[13] = 13.5;
    x[14] = 17.5;

    y[0] = 13890 / 0.3;
    y[1] = 106158 / 0.7;
    y[2] = 183306;
    y[3] = 131450;
    y[4] = 73631;
    y[5] = 37874;
    y[6] = 20346;
    y[7] = 10440;
    y[8] = 5796;
    y[9] = 3125;
    y[10] = 1826;
    y[11] = 1148;
    y[12] = 701;
    y[13] = 1050 / 3.;
    y[14] = 393 / 5.;

    ex[0] = 0.15;
    ex[1] = 0.35;
    ex[2] = 0.5;
    ex[3] = 0.5;
    ex[4] = 0.5;
    ex[5] = 0.5;
    ex[6] = 0.5;
    ex[7] = 0.5;
    ex[8] = 0.5;
    ex[9] = 0.5;
    ex[10] = 0.5;
    ex[11] = 0.5;
    ex[12] = 0.5;
    ex[13] = 1.5;
    ex[14] = 2.5;

    ey[0] = 824 / 0.3;
    ey[1] = 2240 / 0.7;
    ey[2] = 2389;
    ey[3] = 1602;
    ey[4] = 1340;
    ey[5] = 764;
    ey[6] = 432;
    ey[7] = 260;
    ey[8] = 179;
    ey[9] = 121;
    ey[10] = 92;
    ey[11] = 68;
    ey[12] = 51;
    ey[13] = 63 / 3.;
    ey[14] = 37 / 5.;

    ey_sys[0] = 552 / 0.3;
    ey_sys[1] = 2985 / 0.7;
    ey_sys[2] = 3277;
    ey_sys[3] = 3950;
    ey_sys[4] = 2069;
    ey_sys[5] = 1064;
    ey_sys[6] = 329;
    ey_sys[7] = 188;
    ey_sys[8] = 99;
    ey_sys[9] = 46;
    ey_sys[10] = 70;
    ey_sys[11] = 34;
    ey_sys[12] = 23;
    ey_sys[13] = 30 / 3.;
    ey_sys[14] = 23 / 5.;

    SNR[0] = 0.08;
    SNR[1] = 0.07;
    SNR[2] = 0.09;
    SNR[3] = 0.14;
    SNR[4] = 0.21;
    SNR[5] = 0.27;
    SNR[6] = 0.37;
    SNR[7] = 0.49;
    SNR[8] = 0.59;
    SNR[9] = 0.72;
    SNR[10] = 0.74;
    SNR[11] = 0.95;
    SNR[12] = 0.98;
    SNR[13] = 0.95;
    SNR[14] = 1.07;
  } else if (flag == 1) {
    // 20-40%
    x[0] = 0.15;
    x[1] = 0.65;
    x[2] = 1.5;
    x[3] = 2.5;
    x[4] = 3.5;
    x[5] = 4.5;
    x[6] = 5.5;
    x[7] = 6.5;
    x[8] = 7.5;
    x[9] = 8.5;
    x[10] = 9.5;
    x[11] = 10.5;
    x[12] = 11.5;
    x[13] = 13.5;
    x[14] = 17.5;

    y[0] = 6584 / 0.3;
    y[1] = 38709 / 0.7;
    y[2] = 66303;
    y[3] = 49791;
    y[4] = 30467;
    y[5] = 17566;
    y[6] = 9805;
    y[7] = 5789;
    y[8] = 3203;
    y[9] = 1781;
    y[10] = 1023;
    y[11] = 638;
    y[12] = 457;
    y[13] = 509 / 3.;
    y[14] = 249 / 5.;

    ex[0] = 0.15;
    ex[1] = 0.35;
    ex[2] = 0.5;
    ex[3] = 0.5;
    ex[4] = 0.5;
    ex[5] = 0.5;
    ex[6] = 0.5;
    ex[7] = 0.5;
    ex[8] = 0.5;
    ex[9] = 0.5;
    ex[10] = 0.5;
    ex[11] = 0.5;
    ex[12] = 0.5;
    ex[13] = 1.5;
    ex[14] = 2.5;

    ey[0] = 401 / 0.3;
    ey[1] = 984 / 0.7;
    ey[2] = 1241;
    ey[3] = 972;
    ey[4] = 491;
    ey[5] = 306;
    ey[6] = 206;
    ey[7] = 140;
    ey[8] = 97;
    ey[9] = 70;
    ey[10] = 54;
    ey[11] = 41;
    ey[12] = 36;
    ey[13] = 37 / 3.;
    ey[14] = 24 / 5.;

    ey_sys[0] = 206 / 0.3;
    ey_sys[1] = 1052 / 0.7;
    ey_sys[2] = 1748;
    ey_sys[3] = 1272;
    ey_sys[4] = 892;
    ey_sys[5] = 447;
    ey_sys[6] = 189;
    ey_sys[7] = 112;
    ey_sys[8] = 67;
    ey_sys[9] = 49;
    ey_sys[10] = 21;
    ey_sys[11] = 13;
    ey_sys[12] = 14;
    ey_sys[13] = 23 / 3.;
    ey_sys[14] = 6 / 5.;

    SNR[0] = 0.18;
    SNR[1] = 0.14;
    SNR[2] = 0.18;
    SNR[3] = 0.3;
    SNR[4] = 0.44;
    SNR[5] = 0.67;
    SNR[6] = 0.87;
    SNR[7] = 1.19;
    SNR[8] = 1.53;
    SNR[9] = 1.76;
    SNR[10] = 2.01;
    SNR[11] = 2.10;
    SNR[12] = 2.59;
    SNR[13] = 1.92;
    SNR[14] = 2.81;
  } else {
    // 40-90%
    x[0] = 0.15;
    x[1] = 0.65;
    x[2] = 1.5;
    x[3] = 2.5;
    x[4] = 3.5;
    x[5] = 4.5;
    x[6] = 5.5;
    x[7] = 6.5;
    x[8] = 7.5;
    x[9] = 8.5;
    x[10] = 9.5;
    x[11] = 10.5;
    x[12] = 11.5;
    x[13] = 13.5;
    x[14] = 17.5;

    y[0] = 6379 / 0.3;
    y[1] = 14349 / 0.7;
    y[2] = 25545;
    y[3] = 19354;
    y[4] = 13759;
    y[5] = 8491;
    y[6] = 5252;
    y[7] = 3151;
    y[8] = 1743;
    y[9] = 1062;
    y[10] = 584;
    y[11] = 347;
    y[12] = 220;
    y[13] = 381 / 3.;
    y[14] = 123 / 5.;

    ex[0] = 0.15;
    ex[1] = 0.35;
    ex[2] = 0.5;
    ex[3] = 0.5;
    ex[4] = 0.5;
    ex[5] = 0.5;
    ex[6] = 0.5;
    ex[7] = 0.5;
    ex[8] = 0.5;
    ex[9] = 0.5;
    ex[10] = 0.5;
    ex[11] = 0.5;
    ex[12] = 0.5;
    ex[13] = 1.5;
    ex[14] = 2.5;

    ey[0] = 156 / 0.3;
    ey[1] = 389 / 0.7;
    ey[2] = 436;
    ey[3] = 321;
    ey[4] = 215;
    ey[5] = 148;
    ey[6] = 111;
    ey[7] = 78;
    ey[8] = 55;
    ey[9] = 42;
    ey[10] = 31;
    ey[11] = 23;
    ey[12] = 18;
    ey[13] = 24 / 3.;
    ey[14] = 15 / 5.;

    ey_sys[0] = 139 / 0.3;
    ey_sys[1] = 416 / 0.7;
    ey_sys[2] = 609;
    ey_sys[3] = 458;
    ey_sys[4] = 438;
    ey_sys[5] = 221;
    ey_sys[6] = 95;
    ey_sys[7] = 56;
    ey_sys[8] = 41;
    ey_sys[9] = 35;
    ey_sys[10] = 9;
    ey_sys[11] = 21;
    ey_sys[12] = 6;
    ey_sys[13] = 31 / 3.;
    ey_sys[14] = 13 / 5.;

    SNR[0] = 1.23;
    SNR[1] = 0.39;
    SNR[2] = 0.5;
    SNR[3] = 0.75;
    SNR[4] = 1.33;
    SNR[5] = 1.92;
    SNR[6] = 2.73;
    SNR[7] = 3.71;
    SNR[8] = 4.53;
    SNR[9] = 4.67;
    SNR[10] = 4.48;
    SNR[11] = 3.86;
    SNR[12] = 5.87;
    SNR[13] = 5.72;
    SNR[14] = 4.51;
  }
}