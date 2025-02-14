#include "FlowAnalysis_Helper.h"

//______________________________________________________________________________
vector<string> FlowAnalysis_Helper::tokenize(string input_string) {
  vector<string> output_string;
  string temp;
  istringstream ss(input_string);
  while (getline(ss, temp, '/')) {
    output_string.push_back(temp);
  }
  return output_string;
}

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
                                   double centmax, THnSparse *hist_V2,
                                   std::string flag) {

  // Copy original profiles for projections
  THnSparse *hist_V2_cp = dynamic_cast<THnSparse *>(hist_V2->Clone(
      Form("Mass_Pt_centrFT0C_V2_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));

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
TH1D *FlowAnalysis_Helper::GetMassProfile(double ptmin, double ptmax,
                                          double massmin, double massmax,
                                          double centmin, double centmax,
                                          TProfile3D *tp_V2, std::string flag) {

  // Copy original profiles for projections
  TProfile3D *tp_V2_cp = dynamic_cast<TProfile3D *>(tp_V2->Clone(
      Form("Mass_Pt_centrFT0C_V2_%s_Copy_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax)));

  TH3D *tp_V2_cp_projxyz = tp_V2_cp->ProjectionXYZ(
      Form("Mass_Pt_centrFT0C_V2_%s_Pxyz_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2_cp_proj = dynamic_cast<TH1D *>(tp_V2_cp_projxyz->Project3D("x"));
  TH1D *hist_mass_proj = dynamic_cast<TH1D *>(
      tp_V2_cp_proj->Clone(Form("Proj_%s", tp_V2_cp_proj->GetName())));

  delete tp_V2_cp;
  delete tp_V2_cp_projxyz;
  delete tp_V2_cp_proj;

  return hist_mass_proj;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Helper::GetMeanPt(double ptmin, double ptmax, double massmin,
                                     double massmax, double centmin,
                                     double centmax, TProfile3D *tp_V2,
                                     std::string flag) {
  // Copy original profiles for projections
  TProfile3D *tp_V2_cp = dynamic_cast<TProfile3D *>(tp_V2->Clone(
      Form("Mass_Pt_centrFT0C_V2_%s_Copy_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax)));

  TH3D *tp_V2_cp_projxyz = tp_V2_cp->ProjectionXYZ(
      Form("Mass_Pt_centrFT0C_V2_%s_Pxyz_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);

  TH2D *tp_V2_cp_projxy =
      dynamic_cast<TH2D *>(tp_V2_cp_projxyz->Project3D("yx"));
  TProfile *tp_V2_cp_projx = tp_V2_cp_projxy->ProfileX(
      Form("Mass_Pt_centrFT0C_V2_%s_Projx_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax));

  double *Bin_mass_new = CreateBinsFromAxis(tp_V2_cp_projx->GetXaxis());
  int NBins_mass_new = tp_V2_cp_projx->GetXaxis()->GetNbins();
  TH1D *hist_meanPt = new TH1D(Form("ProjPt_%s", tp_V2_cp_projx->GetName()),
                               Form("ProjPt_%s", tp_V2_cp_projx->GetName()),
                               NBins_mass_new, Bin_mass_new);
  hist_meanPt->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_meanPt->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");

  for (int i = 0; i < NBins_mass_new; i++) {
    hist_meanPt->SetBinContent(i + 1, tp_V2_cp_projx->GetBinContent(i + 1));
    hist_meanPt->SetBinError(i + 1, tp_V2_cp_projx->GetBinError(i + 1));
  }

  delete tp_V2_cp;
  delete tp_V2_cp_projxyz;
  delete tp_V2_cp_projxy;
  delete tp_V2_cp_projx;

  return hist_meanPt;
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
TH1D *FlowAnalysis_Helper::GetRfactorProfile(double ptmin, double ptmax,
                                             double massmin, double massmax,
                                             double centmin, double centmax,
                                             TProfile3D *tp_V2MEPM,
                                             TProfile3D *tp_V2MEPP,
                                             TProfile3D *tp_V2MEMM) {
  // Copy original profiles for projections
  TProfile3D *tp_V2MEPM_cp = dynamic_cast<TProfile3D *>(
      tp_V2MEPM->Clone(Form("MEPM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  TProfile3D *tp_V2MEPP_cp = dynamic_cast<TProfile3D *>(
      tp_V2MEPP->Clone(Form("MEPP_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  TProfile3D *tp_V2MEMM_cp = dynamic_cast<TProfile3D *>(
      tp_V2MEMM->Clone(Form("MEMM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));

  TH3D *tp_V2MEPM_cp_projxyz = tp_V2MEPM_cp->ProjectionXYZ(
      Form("MEPM_Mass_Pt_centrFT0C_V2_Pxyz_%g_%g_%g_%g_%g_%g", massmin, massmax,
           ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2MEPM_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2MEPM_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2MEPM_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2MEPM_cp_proj =
      dynamic_cast<TH1D *>(tp_V2MEPM_cp_projxyz->Project3D("x"));

  TH3D *tp_V2MEPP_cp_projxyz = tp_V2MEPP_cp->ProjectionXYZ(
      Form("MEPP_Mass_Pt_centrFT0C_V2_Pxyz_%g_%g_%g_%g_%g_%g", massmin, massmax,
           ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2MEPP_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2MEPP_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2MEPP_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2MEPP_cp_proj =
      dynamic_cast<TH1D *>(tp_V2MEPP_cp_projxyz->Project3D("x"));

  TH3D *tp_V2MEMM_cp_projxyz = tp_V2MEMM_cp->ProjectionXYZ(
      Form("MEMM_Mass_Pt_centrFT0C_V2_Pxyz_%g_%g_%g_%g_%g_%g", massmin, massmax,
           ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2MEMM_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2MEMM_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2MEMM_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2MEMM_cp_proj =
      dynamic_cast<TH1D *>(tp_V2MEMM_cp_projxyz->Project3D("x"));

  // Define resulting histogram
  double *Bin_mass_new = CreateBinsFromAxis(tp_V2MEPM_cp_proj->GetXaxis());
  int NBins_mass_new = tp_V2MEPM_cp_proj->GetXaxis()->GetNbins();
  TH1D *hist_rfactor = new TH1D(Form("Rfactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                Form("Rfactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                NBins_mass_new, Bin_mass_new);
  for (int i = 0; i < NBins_mass_new; i++) {
    double Npm = tp_V2MEPM_cp_proj->GetBinContent(i + 1);
    double Npp = tp_V2MEPP_cp_proj->GetBinContent(i + 1);
    double Nmm = tp_V2MEMM_cp_proj->GetBinContent(i + 1);

    double val = Npp * Nmm <= 0. ? 0. : 0.5 * Npm / pow(Npp * Nmm, 0.5);
    hist_rfactor->SetBinContent(i + 1, val);
  }

  delete tp_V2MEPM_cp;
  delete tp_V2MEPP_cp;
  delete tp_V2MEMM_cp;
  delete tp_V2MEPM_cp_projxyz;
  delete tp_V2MEPP_cp_projxyz;
  delete tp_V2MEMM_cp_projxyz;
  delete tp_V2MEPM_cp_proj;
  delete tp_V2MEPP_cp_proj;
  delete tp_V2MEMM_cp_proj;

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

  double int_V2MEPM = hs_V2MEPM_cp_proj->Integral("width");

  delete hs_V2SEPP_cp;
  delete hs_V2SEMM_cp;
  delete hs_V2MEPM_cp;
  delete hs_V2SEPP_cp_proj;
  delete hs_V2SEMM_cp_proj;
  delete hs_V2MEPM_cp_proj;

  return hist_ffactor->Integral("width") / int_V2MEPM;
}

//______________________________________________________________________________
double FlowAnalysis_Helper::GetFfactorProfile(
    double ptmin, double ptmax, double massmin, double massmax, double centmin,
    double centmax, TProfile3D *tp_V2SEPP, TProfile3D *tp_V2SEMM,
    TProfile3D *tp_V2MEPM, TH1D *hist_rfactor) {
  // Copy original profiles for projections
  TProfile3D *tp_V2SEPP_cp = dynamic_cast<TProfile3D *>(
      tp_V2SEPP->Clone(Form("SEPP_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  TProfile3D *tp_V2SEMM_cp = dynamic_cast<TProfile3D *>(
      tp_V2SEMM->Clone(Form("SEMM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));
  TProfile3D *tp_V2MEPM_cp = dynamic_cast<TProfile3D *>(
      tp_V2MEPM->Clone(Form("MEPM_Mass_Pt_centrFT0C_V2_Copy_%g_%g_%g_%g_%g_%g",
                            massmin, massmax, ptmin, ptmax, centmin, centmax)));

  TH3D *tp_V2SEPP_cp_projxyz = tp_V2SEPP_cp->ProjectionXYZ(
      Form("SEPP_Mass_Pt_centrFT0C_V2_Pxyz_%g_%g_%g_%g_%g_%g", massmin, massmax,
           ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2SEPP_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2SEPP_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2SEPP_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2SEPP_cp_proj =
      dynamic_cast<TH1D *>(tp_V2SEPP_cp_projxyz->Project3D("x"));

  TH3D *tp_V2SEMM_cp_projxyz = tp_V2SEMM_cp->ProjectionXYZ(
      Form("SEMM_Mass_Pt_centrFT0C_V2_Pxyz_%g_%g_%g_%g_%g_%g", massmin, massmax,
           ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2SEMM_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2SEMM_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2SEMM_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2SEMM_cp_proj =
      dynamic_cast<TH1D *>(tp_V2SEMM_cp_projxyz->Project3D("x"));

  TH3D *tp_V2MEPM_cp_projxyz = tp_V2MEPM_cp->ProjectionXYZ(
      Form("MEPM_Mass_Pt_centrFT0C_V2_Pxyz_%g_%g_%g_%g_%g_%g", massmin, massmax,
           ptmin, ptmax, centmin, centmax),
      "B");
  tp_V2MEPM_cp_projxyz->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2MEPM_cp_projxyz->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2MEPM_cp_projxyz->GetZaxis()->SetRangeUser(centmin, centmax);
  TH1D *tp_V2MEPM_cp_proj =
      dynamic_cast<TH1D *>(tp_V2MEPM_cp_projxyz->Project3D("x"));

  double *Bin_mass_new = CreateBinsFromAxis(tp_V2SEPP_cp_proj->GetXaxis());
  int NBins_mass_new = tp_V2SEPP_cp_proj->GetXaxis()->GetNbins();
  TH1D *hist_ffactor = new TH1D(Form("Ffactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                Form("Ffactor_%g_%g_%g_%g_%g_%g", massmin,
                                     massmax, ptmin, ptmax, centmin, centmax),
                                NBins_mass_new, Bin_mass_new);
  for (int i = 0; i < NBins_mass_new; i++) {
    double R_val = hist_rfactor->GetBinContent(i + 1);
    double N_SEPP = tp_V2SEPP_cp_proj->GetBinContent(i + 1);
    double N_SEMM = tp_V2SEMM_cp_proj->GetBinContent(i + 1);

    double F_val =
        N_SEPP * N_SEMM < 0. ? 0. : 2. * R_val * pow(N_SEPP * N_SEMM, 0.5);
    hist_ffactor->SetBinContent(i + 1, F_val);
  }

  double int_V2MEPM = tp_V2MEPM_cp_proj->Integral("width");
  double val_ffactor = hist_ffactor->Integral("width") / int_V2MEPM;

  delete tp_V2SEPP_cp;
  delete tp_V2SEMM_cp;
  delete tp_V2MEPM_cp;
  delete tp_V2SEPP_cp_projxyz;
  delete tp_V2SEMM_cp_projxyz;
  delete tp_V2MEPM_cp_projxyz;
  delete tp_V2SEPP_cp_proj;
  delete tp_V2SEMM_cp_proj;
  delete tp_V2MEPM_cp_proj;
  delete hist_ffactor;

  return val_ffactor;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Helper::GetV2(double ptmin, double ptmax, double massmin,
                                 double massmax, double centmin, double centmax,
                                 THnSparse *hist_V2, double R2SP,
                                 std::string flag) {

  // Copy original profiles for projections
  THnSparse *hist_V2_cp = dynamic_cast<THnSparse *>(hist_V2->Clone(
      Form("Mass_Pt_centrFT0C_V2_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));

  // Set axes' ranges for mass-differential study
  hist_V2_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hist_V2_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hist_V2_cp->GetAxis(3)->SetRangeUser(centmin, centmax);

  // Get v2
  TH2D *hs_v2_sp_proj = hist_V2_cp->Projection(4, 0);

  // Define histograms
  double *Bin_mass_new = CreateBinsFromAxis(hs_v2_sp_proj->GetXaxis());
  int NBins_mass_new = hs_v2_sp_proj->GetXaxis()->GetNbins();
  TH1D *hist_v2sp =
      new TH1D(Form("v2sp_%s_%g_%g", flag.c_str(), ptmin, ptmax),
               Form("v^{#mu#mu}_{2}{SP}_%s_%g_%g", flag.c_str(), ptmin, ptmax),
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
TH1D *FlowAnalysis_Helper::GetV2Profile(double ptmin, double ptmax,
                                        double massmin, double massmax,
                                        double centmin, double centmax,
                                        TProfile3D *tp_V2, std::string flag) {

  // Copy original profiles for projections
  TProfile3D *tp_V2_cp = dynamic_cast<TProfile3D *>(tp_V2->Clone(
      Form("Mass_Pt_centrFT0C_V2_%s_Copy_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax)));
  tp_V2_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  tp_V2_cp->GetYaxis()->SetRangeUser(ptmin, ptmax);
  tp_V2_cp->GetZaxis()->SetRangeUser(centmin, centmax);

  TProfile2D *tp_V2_cp_projxy = tp_V2_cp->Project3DProfile("yx");
  TProfile *tp_V2_cp_projx = tp_V2_cp_projxy->ProfileX(
      Form("Mass_Pt_centrFT0C_V2_%s_Projx_%g_%g_%g_%g_%g_%g", flag.c_str(),
           massmin, massmax, ptmin, ptmax, centmin, centmax));

  double *Bin_mass_new = CreateBinsFromAxis(tp_V2_cp_projx->GetXaxis());
  int NBins_mass_new = tp_V2_cp_projx->GetXaxis()->GetNbins();
  TH1D *hist_v2sp = new TH1D(Form("Proj_%s", tp_V2_cp_projx->GetName()),
                             Form("Proj_%s", tp_V2_cp_projx->GetName()),
                             NBins_mass_new, Bin_mass_new);
  hist_v2sp->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_v2sp->GetYaxis()->SetTitle("v^{#mu#mu}_{2}{SP}");

  for (int i = 0; i < NBins_mass_new; i++) {
    hist_v2sp->SetBinContent(i + 1, tp_V2_cp_projx->GetBinContent(i + 1));
    hist_v2sp->SetBinError(i + 1, tp_V2_cp_projx->GetBinError(i + 1));
  }

  delete tp_V2_cp;
  delete tp_V2_cp_projxy;
  delete tp_V2_cp_projx;

  return hist_v2sp;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Helper::GetV2EM(
    double ptmin, double ptmax, double massmin, double massmax, double centmin,
    double centmax, THnSparse *hs_u2q2_cosDeltaPhi_ME1,
    THnSparse *hs_u2q2_cosDeltaPhi_ME2, TH3F *hs_r2spAB1, TH3F *hs_r2spAC1,
    TH3F *hs_r2spBC1, TH3F *hs_r2spAB2, TH3F *hs_r2spAC2, TH3F *hs_r2spBC2,
    std::string flag) {

  // Copy original profiles for projections
  THnSparse *hs_u2q2_cosDeltaPhi_ME1_cp =
      dynamic_cast<THnSparse *>(hs_u2q2_cosDeltaPhi_ME1->Clone(
          Form("Mass_Pt_centrFT0C_u2q2_cosDeltaPhi_1_Copy_%s_%g_%g",
               flag.c_str(), ptmin, ptmax)));
  THnSparse *hs_u2q2_cosDeltaPhi_ME2_cp =
      dynamic_cast<THnSparse *>(hs_u2q2_cosDeltaPhi_ME2->Clone(
          Form("Mass_Pt_centrFT0C_u2q2_cosDeltaPhi_2_Copy_%s_%g_%g",
               flag.c_str(), ptmin, ptmax)));
  TH3F *hs_r2spAB1_cp = dynamic_cast<TH3F *>(hs_r2spAB1->Clone(Form(
      "Mass_Pt_centrFT0C_R2SPAB_1_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));
  TH3F *hs_r2spAC1_cp = dynamic_cast<TH3F *>(hs_r2spAC1->Clone(Form(
      "Mass_Pt_centrFT0C_R2SPAC_1_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));
  TH3F *hs_r2spBC1_cp = dynamic_cast<TH3F *>(hs_r2spBC1->Clone(Form(
      "Mass_Pt_centrFT0C_R2SPBC_1_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));
  TH3F *hs_r2spAB2_cp = dynamic_cast<TH3F *>(hs_r2spAB2->Clone(Form(
      "Mass_Pt_centrFT0C_R2SPAB_2_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));
  TH3F *hs_r2spAC2_cp = dynamic_cast<TH3F *>(hs_r2spAC2->Clone(Form(
      "Mass_Pt_centrFT0C_R2SPAC_2_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));
  TH3F *hs_r2spBC2_cp = dynamic_cast<TH3F *>(hs_r2spBC2->Clone(Form(
      "Mass_Pt_centrFT0C_R2SPBC_2_Copy_%s_%g_%g", flag.c_str(), ptmin, ptmax)));

  // Set axes' ranges for mass-differential study
  hs_u2q2_cosDeltaPhi_ME1_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_u2q2_cosDeltaPhi_ME1_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_u2q2_cosDeltaPhi_ME1_cp->GetAxis(2)->SetRangeUser(centmin, centmax);

  hs_u2q2_cosDeltaPhi_ME2_cp->GetAxis(0)->SetRangeUser(massmin, massmax);
  hs_u2q2_cosDeltaPhi_ME2_cp->GetAxis(1)->SetRangeUser(ptmin, ptmax);
  hs_u2q2_cosDeltaPhi_ME2_cp->GetAxis(2)->SetRangeUser(centmin, centmax);

  hs_r2spAB1_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  hs_r2spAB1_cp->GetYaxis()->SetRangeUser(centmin, centmax);
  hs_r2spAC1_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  hs_r2spAC1_cp->GetYaxis()->SetRangeUser(centmin, centmax);
  hs_r2spBC1_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  hs_r2spBC1_cp->GetYaxis()->SetRangeUser(centmin, centmax);

  hs_r2spAB2_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  hs_r2spAB2_cp->GetYaxis()->SetRangeUser(centmin, centmax);
  hs_r2spAC2_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  hs_r2spAC2_cp->GetYaxis()->SetRangeUser(centmin, centmax);
  hs_r2spBC2_cp->GetXaxis()->SetRangeUser(massmin, massmax);
  hs_r2spBC2_cp->GetYaxis()->SetRangeUser(centmin, centmax);

  // Get u2q2 and cosDeltaPhi
  TH2D *hs_u2q2_proj1 = hs_u2q2_cosDeltaPhi_ME1_cp->Projection(3, 0);
  TH2D *hs_u2q2_proj2 = hs_u2q2_cosDeltaPhi_ME2_cp->Projection(3, 0);
  TH2D *hs_cosDeltaPhi_proj1 = hs_u2q2_cosDeltaPhi_ME1_cp->Projection(4, 0);
  TH2D *hs_cosDeltaPhi_proj2 = hs_u2q2_cosDeltaPhi_ME2_cp->Projection(4, 0);
  TH2D *hs_r2spAB_proj1 = dynamic_cast<TH2D *>(hs_r2spAB1_cp->Project3D("zx"));
  TH2D *hs_r2spAC_proj1 = dynamic_cast<TH2D *>(hs_r2spAC1_cp->Project3D("zx"));
  TH2D *hs_r2spBC_proj1 = dynamic_cast<TH2D *>(hs_r2spBC1_cp->Project3D("zx"));
  TH2D *hs_r2spAB_proj2 = dynamic_cast<TH2D *>(hs_r2spAB2_cp->Project3D("zx"));
  TH2D *hs_r2spAC_proj2 = dynamic_cast<TH2D *>(hs_r2spAC2_cp->Project3D("zx"));
  TH2D *hs_r2spBC_proj2 = dynamic_cast<TH2D *>(hs_r2spBC2_cp->Project3D("zx"));

  // Define histograms
  double *Bin_mass_new = CreateBinsFromAxis(hs_u2q2_proj1->GetXaxis());
  int NBins_mass_new = hs_u2q2_proj1->GetXaxis()->GetNbins();
  TH1D *hist_v2spme =
      new TH1D(Form("v2sp_%s_%g_%g", flag.c_str(), ptmin, ptmax),
               Form("v^{#mu#mu}_{2}{SP}_%s_%g_%g", flag.c_str(), ptmin, ptmax),
               NBins_mass_new, Bin_mass_new);
  hist_v2spme->GetXaxis()->SetTitle("mass (GeV/c2)");
  hist_v2spme->GetYaxis()->SetTitle("v^{#mu#mu}_{2}{SP}");

  // Evaluation of differential flow as function of invariant mass
  for (int i = 0; i < NBins_mass_new; i++) {
    TH2D *hs_u2q2_proj1_cp =
        dynamic_cast<TH2D *>(hs_u2q2_proj1->Clone("u2q2_1_Copy"));
    hs_u2q2_proj1_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                               Bin_mass_new[i + 1]);
    TH2D *hs_u2q2_proj2_cp =
        dynamic_cast<TH2D *>(hs_u2q2_proj2->Clone("u2q2_2_Copy"));
    hs_u2q2_proj2_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                               Bin_mass_new[i + 1]);

    TH2D *hs_cosDeltaPhi_proj1_cp =
        dynamic_cast<TH2D *>(hs_cosDeltaPhi_proj1->Clone("cosDeltaPhi_1_Copy"));
    hs_cosDeltaPhi_proj1_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                      Bin_mass_new[i + 1]);
    TH2D *hs_cosDeltaPhi_proj2_cp =
        dynamic_cast<TH2D *>(hs_cosDeltaPhi_proj2->Clone("cosDeltaPhi_2_Copy"));
    hs_cosDeltaPhi_proj2_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                      Bin_mass_new[i + 1]);

    TH2D *hs_r2spAB_proj1_cp =
        dynamic_cast<TH2D *>(hs_r2spAB_proj1->Clone("r2spAB_1_Copy"));
    hs_r2spAB_proj1_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                 Bin_mass_new[i + 1]);
    TH2D *hs_r2spAC_proj1_cp =
        dynamic_cast<TH2D *>(hs_r2spAC_proj1->Clone("r2spAC_1_Copy"));
    hs_r2spAC_proj1_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                 Bin_mass_new[i + 1]);
    TH2D *hs_r2spBC_proj1_cp =
        dynamic_cast<TH2D *>(hs_r2spBC_proj1->Clone("r2spBC_1_Copy"));
    hs_r2spBC_proj1_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                 Bin_mass_new[i + 1]);

    TH2D *hs_r2spAB_proj2_cp =
        dynamic_cast<TH2D *>(hs_r2spAB_proj2->Clone("r2spAB_2_Copy"));
    hs_r2spAB_proj2_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                 Bin_mass_new[i + 1]);
    TH2D *hs_r2spAC_proj2_cp =
        dynamic_cast<TH2D *>(hs_r2spAC_proj2->Clone("r2spAC_2_Copy"));
    hs_r2spAC_proj2_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                 Bin_mass_new[i + 1]);
    TH2D *hs_r2spBC_proj2_cp =
        dynamic_cast<TH2D *>(hs_r2spBC_proj2->Clone("r2spBC_2_Copy"));
    hs_r2spBC_proj2_cp->GetXaxis()->SetRangeUser(Bin_mass_new[i],
                                                 Bin_mass_new[i + 1]);

    double u2q2_1 = hs_u2q2_proj1_cp->GetMean(2);
    double u2q2e_1 = hs_u2q2_proj1_cp->GetMean(12);
    double u2q2_2 = hs_u2q2_proj2_cp->GetMean(2);
    double u2q2e_2 = hs_u2q2_proj2_cp->GetMean(12);

    double cosDeltaPhi_1 = hs_cosDeltaPhi_proj1_cp->GetMean(2);
    double cosDeltaPhie_1 = hs_cosDeltaPhi_proj1_cp->GetMean(12);
    double cosDeltaPhi_2 = hs_cosDeltaPhi_proj2_cp->GetMean(2);
    double cosDeltaPhie_2 = hs_cosDeltaPhi_proj2_cp->GetMean(12);

    double r2spAB_1 = hs_r2spAB_proj1_cp->GetMean(2);
    double r2spAC_1 = hs_r2spAC_proj1_cp->GetMean(2);
    double r2spBC_1 = hs_r2spBC_proj1_cp->GetMean(2);

    double r2spAB_2 = hs_r2spAB_proj2_cp->GetMean(2);
    double r2spAC_2 = hs_r2spAC_proj2_cp->GetMean(2);
    double r2spBC_2 = hs_r2spBC_proj2_cp->GetMean(2);

    double r2sp_1 = r2spBC_1 != 0 ? r2spAB_1 * r2spAC_1 / r2spBC_1 : 0.0;
    r2sp_1 = r2sp_1 > 0 ? TMath::Sqrt(r2sp_1) : 0.0;
    double r2sp_2 = r2spBC_2 != 0 ? r2spAB_2 * r2spAC_2 / r2spBC_2 : 0.0;
    r2sp_2 = r2sp_2 > 0 ? TMath::Sqrt(r2sp_2) : 0.0;
    double v2 =
        r2sp_1 == 0 || r2sp_2 == 0
            ? 0.0
            : u2q2_1 / r2sp_1 * cosDeltaPhi_1 + u2q2_2 / r2sp_2 * cosDeltaPhi_2;
    double v2e =
        r2sp_1 == 0 || r2sp_2 == 0
            ? 0.0
            : TMath::Sqrt(
                  pow(cosDeltaPhi_1, 2.) * pow(u2q2e_1, 2.) / pow(r2sp_1, 2.) +
                  pow(cosDeltaPhie_1, 2.) * pow(u2q2_1, 2.) / pow(r2sp_1, 2.) +
                  pow(cosDeltaPhi_2, 2.) * pow(u2q2e_2, 2.) / pow(r2sp_2, 2.) +
                  pow(cosDeltaPhie_2, 2.) * pow(u2q2_2, 2.) / pow(r2sp_2, 2.));

    hist_v2spme->SetBinContent(i + 1, v2);
    hist_v2spme->SetBinError(i + 1, v2e);

    delete hs_u2q2_proj1_cp;
    delete hs_u2q2_proj2_cp;
    delete hs_cosDeltaPhi_proj1_cp;
    delete hs_cosDeltaPhi_proj2_cp;
    delete hs_r2spAB_proj1_cp;
    delete hs_r2spAB_proj2_cp;
    delete hs_r2spAC_proj1_cp;
    delete hs_r2spAC_proj2_cp;
    delete hs_r2spBC_proj1_cp;
    delete hs_r2spBC_proj2_cp;
  }

  delete hs_u2q2_cosDeltaPhi_ME1_cp;
  delete hs_u2q2_cosDeltaPhi_ME2_cp;
  delete hs_r2spAB1_cp;
  delete hs_r2spAB2_cp;
  delete hs_r2spAC1_cp;
  delete hs_r2spAC2_cp;
  delete hs_r2spBC1_cp;
  delete hs_r2spBC2_cp;
  delete hs_u2q2_proj1;
  delete hs_u2q2_proj2;
  delete hs_cosDeltaPhi_proj1;
  delete hs_cosDeltaPhi_proj2;
  delete hs_r2spAB_proj1;
  delete hs_r2spAB_proj2;
  delete hs_r2spAC_proj1;
  delete hs_r2spAC_proj2;
  delete hs_r2spBC_proj1;
  delete hs_r2spBC_proj2;

  return hist_v2spme;
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
    double ptmin, double ptmax, double centmin, double centmax, int nBins,
    int index, TCanvas *c_sys_yield, TCanvas *c_sys_v2, TCanvas *c_sys_meanPt,
    TH1D *hist_sys_yield, TH1D *hist_sys_v2, TH1D *hist_sys_meanPt,
    double *bins_sys_yield, double *bins_sys_v2, double *chi2_yield,
    double *chi2_v2, double *chi2_meanPt, int nbCombo_yield, int nbCombo_v2,
    vector<double> stats_yield, vector<double> stats_v2,
    vector<double> stats_meanPt, double *pt_bins, TList *ls_sys_yield,
    TList *ls_sys_v2, TList *ls_sys_meanPt, std::string flag) {

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
  if (flag == "pt") {
    text_sys_yield->DrawLatexNDC(
        .12, .85,
        Form("N_{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
             "(sys)",
             pt_bins[index], pt_bins[index + 1], stats_yield[0], stats_yield[1],
             100. * stats_yield[1] / stats_yield[0], stats_yield[2],
             100. * stats_yield[2] / stats_yield[0]));
  } else {
    text_sys_yield->DrawLatexNDC(
        .12, .85,
        Form("N_{J/#psi} [%g-%g] (%%) = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
             "(sys)",
             pt_bins[index], pt_bins[index + 1], stats_yield[0], stats_yield[1],
             100. * stats_yield[1] / stats_yield[0], stats_yield[2],
             100. * stats_yield[2] / stats_yield[0]));
  }
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

  // Plotting for mean Pt systematics
  c_sys_meanPt->cd();
  c_sys_meanPt->SetBottomMargin(0);
  c_sys_meanPt->SetCanvasSize(1000, 400);
  TPad *pad_sys_meanPt =
      new TPad(Form("pad_sys_meanPt_%d", index),
               Form("pad_sys_meanPt_%d", index), 0, 0.3, 1, 1.0);
  pad_sys_meanPt->SetBottomMargin(0);
  pad_sys_meanPt->Draw();
  pad_sys_meanPt->cd();
  hist_sys_meanPt->SetMarkerStyle(20);
  hist_sys_meanPt->SetMarkerSize(1.);
  hist_sys_meanPt->SetMarkerColor(kBlack);
  hist_sys_meanPt->SetLineColor(kBlack);
  hist_sys_meanPt->SetLineWidth(2);
  hist_sys_meanPt->SetFillStyle(0);
  hist_sys_meanPt->SetStats(0);
  hist_sys_meanPt->SetTitle("");
  hist_sys_meanPt->GetYaxis()->SetTitle("<#it{p}_{T}> (GeV/c)");
  hist_sys_meanPt->GetYaxis()->SetTitleSize(0.05);
  hist_sys_meanPt->GetYaxis()->SetTitleOffset(0.5);
  hist_sys_meanPt->GetXaxis()->SetLabelOffset(999);
  hist_sys_meanPt->GetXaxis()->SetLabelSize(0);
  hist_sys_meanPt->Draw("HIST EP");
  TF1 *lmeanPt_mean = new TF1("meanmeanPt", "[0]", bins_sys_yield[0],
                              bins_sys_yield[nbCombo_yield]);
  lmeanPt_mean->SetParameter(0, stats_meanPt[0]);
  lmeanPt_mean->SetLineColor(kBlue);
  lmeanPt_mean->SetLineWidth(3);
  lmeanPt_mean->SetLineStyle(1);
  lmeanPt_mean->Draw("same");
  TF1 *lmeanPt_meanerrorp =
      new TF1("meanerrorpmeanPt", "[0]+[1]", bins_sys_yield[0],
              bins_sys_yield[nbCombo_yield]);
  lmeanPt_meanerrorp->SetParameter(0, stats_meanPt[0]);
  lmeanPt_meanerrorp->SetParameter(1, stats_meanPt[1]);
  lmeanPt_meanerrorp->SetLineColor(kBlue);
  lmeanPt_meanerrorp->SetLineWidth(3);
  lmeanPt_meanerrorp->SetLineStyle(7);
  lmeanPt_meanerrorp->Draw("same");
  TF1 *lmeanPt_meanerrorm =
      new TF1("meanerrormmeanPt", "[0]-[1]", bins_sys_yield[0],
              bins_sys_yield[nbCombo_yield]);
  lmeanPt_meanerrorm->SetParameter(0, stats_meanPt[0]);
  lmeanPt_meanerrorm->SetParameter(1, stats_meanPt[1]);
  lmeanPt_meanerrorm->SetLineColor(kBlue);
  lmeanPt_meanerrorm->SetLineWidth(3);
  lmeanPt_meanerrorm->SetLineStyle(7);
  lmeanPt_meanerrorm->Draw("same");
  TF1 *lmeanPt_rmsp = new TF1("rmspmeanPt", "[0]+[1]", bins_sys_yield[0],
                              bins_sys_yield[nbCombo_yield]);
  lmeanPt_rmsp->SetParameter(0, stats_meanPt[0]);
  lmeanPt_rmsp->SetParameter(1, stats_meanPt[2]);
  lmeanPt_rmsp->SetLineColor(kBlue);
  lmeanPt_rmsp->SetLineWidth(3);
  lmeanPt_rmsp->SetLineStyle(9);
  lmeanPt_rmsp->Draw("same");
  TF1 *lmeanPt_rmsm = new TF1("rmsmmeanPt", "[0]-[1]", bins_sys_yield[0],
                              bins_sys_yield[nbCombo_yield]);
  lmeanPt_rmsm->SetParameter(0, stats_meanPt[0]);
  lmeanPt_rmsm->SetParameter(1, stats_meanPt[2]);
  lmeanPt_rmsm->SetLineColor(kBlue);
  lmeanPt_rmsm->SetLineWidth(3);
  lmeanPt_rmsm->SetLineStyle(9);
  lmeanPt_rmsm->Draw("same");
  TLatex *text_sys_meanPt = new TLatex();
  text_sys_meanPt->SetTextSize(0.05);
  text_sys_meanPt->SetTextFont(42);
  text_sys_meanPt->SetTextColor(kBlue);
  text_sys_meanPt->DrawLatexNDC(
      .12, .85,
      Form("<#it{p}_{T}>^{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm "
           "%g[%g%%] "
           "(sys)",
           pt_bins[index], pt_bins[index + 1], stats_meanPt[0], stats_meanPt[1],
           100. * stats_meanPt[1] / stats_meanPt[0], stats_meanPt[2],
           100. * stats_meanPt[2] / stats_meanPt[0]));
  pad_sys_meanPt->ModifiedUpdate();
  c_sys_meanPt->cd();
  TPad *pad_sys_meanPt_chi =
      new TPad(Form("pad_sys_meanPt_chi2_%d", index),
               Form("pad_sys_meanPt_chi2_%d", index), 0, 0., 1, 0.3);
  pad_sys_meanPt_chi->SetTopMargin(0);
  pad_sys_meanPt_chi->SetBottomMargin(0.6);
  pad_sys_meanPt_chi->Draw();
  pad_sys_meanPt_chi->cd();
  TH1D *hist_chi2_meanPt =
      (TH1D *)hist_sys_meanPt->Clone(Form("hist_meanPt_chi2_%d", index));
  for (int j = 0; j < hist_chi2_meanPt->GetNbinsX(); j++) {
    hist_chi2_meanPt->SetBinContent(j + 1, chi2_meanPt[j]);
    hist_chi2_meanPt->SetBinError(j + 1, 0.);
  }
  hist_chi2_meanPt->SetTitle("");
  hist_chi2_meanPt->GetYaxis()->SetLabelSize(0.05);
  hist_chi2_meanPt->GetYaxis()->SetRangeUser(0.5, 3.);
  hist_chi2_meanPt->GetYaxis()->SetTitle("#chi^{2}/ndf");
  hist_chi2_meanPt->GetYaxis()->SetTitleSize(0.08);
  hist_chi2_meanPt->GetYaxis()->SetTitleOffset(0.25);
  hist_chi2_meanPt->GetXaxis()->SetLabelSize(0.14);
  hist_chi2_meanPt->GetXaxis()->SetLabelOffset(0.02);
  hist_chi2_meanPt->Draw("HIST P");
  TF1 *lchi2_meanPt1 = new TF1("lchi2_meanPt1", "[0]", bins_sys_yield[0],
                               bins_sys_yield[nbCombo_yield]);
  lchi2_meanPt1->SetParameter(0, 1.0);
  lchi2_meanPt1->SetLineColor(kBlue);
  lchi2_meanPt1->SetLineWidth(2);
  lchi2_meanPt1->SetLineStyle(9);
  lchi2_meanPt1->Draw("same");
  TF1 *lchi2_meanPt2 = new TF1("lchi2_meanPt2", "[0]", bins_sys_yield[0],
                               bins_sys_yield[nbCombo_yield]);
  lchi2_meanPt2->SetParameter(0, 2.0);
  lchi2_meanPt2->SetLineColor(kRed);
  lchi2_meanPt2->SetLineWidth(2);
  lchi2_meanPt2->SetLineStyle(1);
  lchi2_meanPt2->Draw("same");

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
  if (flag == "pt") {
    text_sys_v2->DrawLatexNDC(
        .12, .85,
        Form("v_{2}^{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm "
             "%g[%g%%] "
             "(sys)",
             pt_bins[index], pt_bins[index + 1], stats_v2[0], stats_v2[1],
             100. * stats_v2[1] / stats_v2[0], stats_v2[2],
             100. * stats_v2[2] / stats_v2[0]));
  } else {
    text_sys_v2->DrawLatexNDC(
        .12, .85,
        Form(
            "v_{2}^{J/#psi} [%g-%g] (%%) = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
            "(sys)",
            pt_bins[index], pt_bins[index + 1], stats_v2[0], stats_v2[1],
            100. * stats_v2[1] / stats_v2[0], stats_v2[2],
            100. * stats_v2[2] / stats_v2[0]));
  }
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

  c_sys_v2->SaveAs(Form("FitSys_v2EM%d_%g_%g_%g_%g_%d.pdf", nbCombo_v2, ptmin,
                        ptmax, centmin, centmax, nBins));
  ls_sys_v2->Add(c_sys_v2);

  c_sys_yield->SaveAs(Form("FitSys_massEM%d_%g_%g_%g_%g_%d.pdf", nbCombo_yield,
                           ptmin, ptmax, centmin, centmax, nBins));
  ls_sys_yield->Add(c_sys_yield);

  ls_sys_meanPt->SaveAs(Form("FitSys_meanPtEM%d_%g_%g_%g_%g_%d.pdf",
                             nbCombo_yield, ptmin, ptmax, centmin, centmax,
                             nBins));
  ls_sys_meanPt->Add(c_sys_meanPt);
}

//______________________________________________________________________________
void FlowAnalysis_Helper::PlotSystematicsNoMeanPt(
    double ptmin, double ptmax, double centmin, double centmax, int nBins,
    int index, TCanvas *c_sys_yield, TCanvas *c_sys_v2, TH1D *hist_sys_yield,
    TH1D *hist_sys_v2, double *bins_sys_yield, double *bins_sys_v2,
    double *chi2_yield, double *chi2_v2, int nbCombo_yield, int nbCombo_v2,
    vector<double> stats_yield, vector<double> stats_v2, double *pt_bins,
    TList *ls_sys_yield, TList *ls_sys_v2, std::string flag) {

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
  if (flag == "pt") {
    text_sys_yield->DrawLatexNDC(
        .12, .85,
        Form("N_{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
             "(sys)",
             pt_bins[index], pt_bins[index + 1], stats_yield[0], stats_yield[1],
             100. * stats_yield[1] / stats_yield[0], stats_yield[2],
             100. * stats_yield[2] / stats_yield[0]));
  } else {
    text_sys_yield->DrawLatexNDC(
        .12, .85,
        Form("N_{J/#psi} [%g-%g] (%%) = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
             "(sys)",
             pt_bins[index], pt_bins[index + 1], stats_yield[0], stats_yield[1],
             100. * stats_yield[1] / stats_yield[0], stats_yield[2],
             100. * stats_yield[2] / stats_yield[0]));
  }
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
  if (flag == "pt") {
    text_sys_v2->DrawLatexNDC(
        .12, .85,
        Form("v_{2}^{J/#psi} [%g-%g] GeV/c = %g #pm %g[%g%%] (stat) #pm "
             "%g[%g%%] "
             "(sys)",
             pt_bins[index], pt_bins[index + 1], stats_v2[0], stats_v2[1],
             100. * stats_v2[1] / stats_v2[0], stats_v2[2],
             100. * stats_v2[2] / stats_v2[0]));
  } else {
    text_sys_v2->DrawLatexNDC(
        .12, .85,
        Form(
            "v_{2}^{J/#psi} [%g-%g] (%%) = %g #pm %g[%g%%] (stat) #pm %g[%g%%] "
            "(sys)",
            pt_bins[index], pt_bins[index + 1], stats_v2[0], stats_v2[1],
            100. * stats_v2[1] / stats_v2[0], stats_v2[2],
            100. * stats_v2[2] / stats_v2[0]));
  }
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

  c_sys_v2->SaveAs(Form("FitSys_v2EM%d_%g_%g_%g_%g_%d.pdf", nbCombo_v2, ptmin,
                        ptmax, centmin, centmax, nBins));
  ls_sys_v2->Add(c_sys_v2);

  c_sys_yield->SaveAs(Form("FitSys_massEM%d_%g_%g_%g_%g_%d.pdf", nbCombo_yield,
                           ptmin, ptmax, centmin, centmax, nBins));
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
    int size_ptbin, double cent_min, double cent_max, double *pt_bins,
    double *x_v2pt, double *y_v2pt, double *ex_v2pt, double *ey_v2pt,
    double *eysys_v2pt, double *x_run2, double *y_run2, double *ex_run2,
    double *ey_run2, double *eysys_run2, double *x_yield, double *y_yield,
    double *ex_yield, double *ey_yield, double *eysys_yield,
    double *x_yield_run2, double *y_yield_run2, double *ex_yield_run2,
    double *ey_yield_run2, double *eysys_yield_run2, TList *ls) {
  // Compare with Run2 data: 5.02 TeV
  TCanvas *c_yield = new TCanvas("jpsi_yield_pT", "jpsi_yield_pT");
  TCanvas *c_yield_vsRun2 =
      new TCanvas("jpsi_yield_raw_pT", "jpsi_yield_raw_pT");
  TCanvas *c_pt = new TCanvas("v2_pT", "v2_pT");

  TGraphMultiErrors *graph_v2pt =
      new TGraphMultiErrors("graph_v2_pt", "", size_ptbin, x_v2pt, y_v2pt,
                            ex_v2pt, ex_v2pt, ey_v2pt, ey_v2pt);
  graph_v2pt->SetTitle(Form("Run3 #sqrt{#it{s}_{NN}} = 5.36 TeV, %d-%d%%",
                            int(cent_min), int(cent_max)));
  graph_v2pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");
  graph_v2pt->AddYError(size_ptbin, eysys_v2pt, eysys_v2pt);

  TGraphMultiErrors *graph_v2pt_run2 =
      new TGraphMultiErrors("graph_v2_pt_run2", "", 10, x_run2, y_run2, ex_run2,
                            ex_run2, ey_run2, ey_run2);
  graph_v2pt_run2->AddYError(10, eysys_run2, eysys_run2);
  graph_v2pt_run2->SetTitle(Form("Run2 #sqrt{#it{s}_{NN}} = 5.02 TeV, %d-%d%%",
                                 int(cent_min), int(cent_max)));
  graph_v2pt_run2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_v2pt_run2->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

  TGraphMultiErrors *graph_yield =
      new TGraphMultiErrors("graph_yields_pt", "Run3", size_ptbin, x_yield,
                            y_yield, ex_yield, ex_yield, ey_yield, ey_yield);
  graph_yield->SetTitle(Form("Run3 #sqrt{#it{s}_{NN}} = 5.36 TeV, %d-%d%%",
                             int(cent_min), int(cent_max)));
  graph_yield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  graph_yield->GetYaxis()->SetTitle("dN_{J/#psi}/d#it{p}_{T} (GeV/c)^{-1}");
  graph_yield->AddYError(size_ptbin, eysys_yield, eysys_yield);

  TGraphMultiErrors *graph_yield_run2 = new TGraphMultiErrors(
      "graph_yields_run2_pt", "Run2", 15, x_yield_run2, y_yield_run2,
      ex_yield_run2, ex_yield_run2, ey_yield_run2, ey_yield_run2);
  graph_yield_run2->SetTitle(Form("Run2 #sqrt{#it{s}_{NN}} = 5.02 TeV, %d-%d%%",
                                  int(cent_min), int(cent_max)));
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
    mg_raw->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    mg_raw->GetYaxis()->SetTitle("dN_{J/#psi}/d#it{p}_{T} (GeV/c)^{-1}");
    mg_raw->SetTitle("");
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
    mg_raw->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    mg_raw->GetYaxis()->SetTitle("dN_{J/#psi}/d#it{p}_{T} (GeV/c)^{-1}");
    mg_raw->SetTitle("");
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
  mg->Add(graph_v2pt_run2);
  mg->Add(graph_v2pt);
  mg->GetXaxis()->SetRangeUser(0, mg->GetXaxis()->GetXmax());
  mg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  mg->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
  mg->SetTitle("");
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
void FlowAnalysis_Helper::PlotFinalResultsCent(
    int size_centbin, double pt_min, double pt_max, double *cent_bins,
    double *x_v2cent, double *y_v2cent, double *ex_v2cent, double *ey_v2cent,
    double *eysys_v2cent, double *x_run2, double *y_run2, double *ex_run2,
    double *ey_run2, double *eysys_run2, double *x_yield, double *y_yield,
    double *ex_yield, double *ey_yield, double *eysys_yield, TList *ls) {
  // Compare with Run2 data: 5.02 TeV
  TCanvas *c_yield = new TCanvas("jpsi_yield_cent", "jpsi_yield_cent");
  TCanvas *c_cent = new TCanvas("v2_cent", "v2_cent");

  TGraphMultiErrors *graph_v2cent =
      new TGraphMultiErrors("graph_v2_pt", "", size_centbin, x_v2cent, y_v2cent,
                            ex_v2cent, ex_v2cent, ey_v2cent, ey_v2cent);
  graph_v2cent->SetTitle(Form("Run3 #sqrt{#it{s}_{NN}} = 5.36 TeV, %d-%d GeV/c",
                              int(pt_min), int(pt_max)));
  graph_v2cent->GetXaxis()->SetTitle("Centrality (%)");
  graph_v2cent->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");
  graph_v2cent->AddYError(size_centbin, eysys_v2cent, eysys_v2cent);

  TGraphMultiErrors *graph_v2cent_run2 =
      new TGraphMultiErrors("graph_v2_cent_run2", "", 7, x_run2, y_run2,
                            ex_run2, ex_run2, ey_run2, ey_run2);
  graph_v2cent_run2->AddYError(7, eysys_run2, eysys_run2);
  graph_v2cent_run2->SetTitle(
      Form("Run2 #sqrt{#it{s}_{NN}} = 5.02 TeV, %d-%d GeV/c", int(pt_min),
           int(pt_max)));
  graph_v2cent_run2->GetXaxis()->SetTitle("Centrality (%)");
  graph_v2cent_run2->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

  TGraphMultiErrors *graph_yield =
      new TGraphMultiErrors("graph_yields_cent", "Run3", size_centbin, x_yield,
                            y_yield, ex_yield, ex_yield, ey_yield, ey_yield);
  graph_yield->SetTitle(Form("Run3 #sqrt{#it{s}_{NN}} = 5.36 TeV, %d-%d GeV/c",
                             int(pt_min), int(pt_max)));
  graph_yield->GetXaxis()->SetTitle("Centrality (%)");
  graph_yield->GetYaxis()->SetTitle("N_{J/#psi}");
  graph_yield->AddYError(size_centbin, eysys_yield, eysys_yield);

  // Save final results
  // Saving plot for J/psi yields as function of centrality
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
  graph_yield->GetXaxis()->SetLimits(cent_bins[0], cent_bins[size_centbin]);
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

  // Saving plot for J/psi v2 as function of centrality and
  // compared with Run2
  c_cent->cd();
  TPad *pad_cent_final =
      new TPad("pad_cent_final", "pad_cent_final", 0, 0, 1, 1);
  pad_cent_final->Draw();
  pad_cent_final->cd();
  auto mg = new TMultiGraph();
  graph_v2cent->SetMarkerStyle(20);
  graph_v2cent->SetMarkerSize(1.);
  graph_v2cent->SetMarkerColor(kBlue);
  graph_v2cent->SetLineColor(kBlue);
  graph_v2cent->SetLineWidth(2);
  graph_v2cent->SetFillStyle(0);
  graph_v2cent_run2->SetMarkerStyle(20);
  graph_v2cent_run2->SetMarkerSize(1.);
  graph_v2cent_run2->SetMarkerColor(kRed);
  graph_v2cent_run2->SetLineColor(kRed);
  graph_v2cent_run2->SetLineWidth(2);
  graph_v2cent_run2->SetFillStyle(0);
  mg->Add(graph_v2cent);
  mg->Add(graph_v2cent_run2);
  mg->GetXaxis()->SetLimits(cent_bins[0], cent_bins[size_centbin]);
  mg->GetXaxis()->SetTitle("Centrality (%)");
  mg->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
  mg->SetTitle("");
  mg->Draw("A P Z ; Z ; 5 s=0.5");
  pad_cent_final->BuildLegend();
  TLatex *text_cent = new TLatex();
  text_cent->SetTextSize(0.04);
  text_cent->SetTextFont(42);
  text_cent->DrawLatexNDC(.18, .82,
                          "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                          "= 5.36 TeV");
  pad_cent_final->ModifiedUpdate();
  text_cent->DrawLatexNDC(.18, .77,
                          "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                          "4");
  text_cent->DrawLatexNDC(.18, .72,
                          Form("%g < #it{p}_{T} < %g GeV/c", pt_min, pt_max));
  pad_cent_final->ModifiedUpdate();
  TF1 *lv2_cent =
      new TF1("lv2_cent", "[0]", cent_bins[0], cent_bins[size_centbin]);
  lv2_cent->SetParameter(0, 0.);
  lv2_cent->SetLineColor(18);
  lv2_cent->SetLineWidth(3);
  lv2_cent->SetLineStyle(9);
  lv2_cent->Draw("same");
  pad_cent_final->ModifiedUpdate();
  ls->Add(c_cent);
}

//______________________________________________________________________________
void FlowAnalysis_Helper::PlotSEME(std::string flag, double ptmin, double ptmax,
                                   double massmin, double massmax,
                                   double centmin, double centmax,
                                   TH1D *hist_SE, TH1D *hist_ME, TList *ls) {
  TCanvas *c_SEME =
      new TCanvas(Form("hist_Mass_%s_%g_%g_%g_%g_%g_%g", flag.c_str(), ptmin,
                       ptmax, massmin, massmax, centmin, centmax),
                  Form("hist_Mass_%s_%g_%g_%g_%g_%g_%g", flag.c_str(), ptmin,
                       ptmax, massmin, massmax, centmin, centmax));
  c_SEME->cd();

  hist_SE->SetStats(0);
  hist_SE->SetTitle(Form("Mass spectra with event-mixing: %s", flag.c_str()));
  hist_SE->SetMarkerStyle(20);
  hist_SE->SetMarkerSize(0.8);
  hist_SE->SetMarkerColor(kBlue);
  hist_SE->GetYaxis()->SetTitle("Counts");
  hist_SE->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c2)");
  hist_SE->Draw("HIST EP");

  hist_ME->SetStats(0);
  hist_ME->SetTitle("");
  hist_ME->SetMarkerStyle(20);
  hist_ME->SetMarkerSize(0.8);
  hist_ME->SetMarkerColor(kRed);
  hist_ME->GetYaxis()->SetTitle("Counts");
  hist_ME->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c2)");
  hist_ME->Draw("HIST EP same");
  gPad->BuildLegend();
  ls->Add(c_SEME);
}

//______________________________________________________________________________
void FlowAnalysis_Helper::LoadReso(std::string FileName, TH2F *&hs_R2SPAB,
                                   TH2F *&hs_R2SPAC, TH2F *&hs_R2SPBC) {
  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist_r2;
  TList *sublist_r2;
  if (filePath.extension() == ".root") {
    // Load data from AnalysisResults.root
    TFile *Input_File = TFile::Open(FileName.c_str());
    list_hist_r2 =
        (THashList *)Input_File->Get("analysis-event-selection/output");
    sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

    // Get histograms of correlations and resolution factors
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
        list_hist_r2 =
            (THashList *)inFile->Get("analysis-event-selection/output");
        sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

        if (first) {
          hs_R2SPAB = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C");
          hs_R2SPAC = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C");
          hs_R2SPBC = (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C");
        } else {
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
void FlowAnalysis_Helper::LoadResoProfile(
    std::string FileName, TProfile *&tp_R2SPAB, TProfile *&tp_R2SPAC,
    TProfile *&tp_R2SPBC, TProfile *&tp_R2EPAB, TProfile *&tp_R2EPAC,
    TProfile *&tp_R2EPBC, TProfile *&tp_R2SPAB_Im, TProfile *&tp_R2SPAC_Im,
    TProfile *&tp_R2SPBC_Im, TProfile *&tp_R2EPAB_Im, TProfile *&tp_R2EPAC_Im,
    TProfile *&tp_R2EPBC_Im) {
  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist_r2;
  TList *sublist_r2;
  if (filePath.extension() == ".root") {
    // Load data from AnalysisResults.root
    TFile *Input_File = TFile::Open(FileName.c_str());
    list_hist_r2 =
        (THashList *)Input_File->Get("analysis-event-selection/output");
    sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

    // Get histograms of correlations and resolution factors
    tp_R2SPAB =
        (TProfile *)sublist_r2->FindObject("Profile_R2SP_TPCFT0A_CentFT0C");
    tp_R2SPAC =
        (TProfile *)sublist_r2->FindObject("Profile_R2SP_TPCFT0C_CentFT0C");
    tp_R2SPBC =
        (TProfile *)sublist_r2->FindObject("Profile_R2SP_FT0AFT0C_CentFT0C");
    tp_R2EPAB =
        (TProfile *)sublist_r2->FindObject("Profile_R2EP_TPCFT0A_CentFT0C");
    tp_R2EPAC =
        (TProfile *)sublist_r2->FindObject("Profile_R2EP_TPCFT0C_CentFT0C");
    tp_R2EPBC =
        (TProfile *)sublist_r2->FindObject("Profile_R2EP_FT0AFT0C_CentFT0C");

    tp_R2SPAB_Im =
        (TProfile *)sublist_r2->FindObject("Profile_R2SP_Im_TPCFT0A_CentFT0C");
    tp_R2SPAC_Im =
        (TProfile *)sublist_r2->FindObject("Profile_R2SP_Im_TPCFT0C_CentFT0C");
    tp_R2SPBC_Im =
        (TProfile *)sublist_r2->FindObject("Profile_R2SP_Im_FT0AFT0C_CentFT0C");
    tp_R2EPAB_Im =
        (TProfile *)sublist_r2->FindObject("Profile_R2EP_Im_TPCFT0A_CentFT0C");
    tp_R2EPAC_Im =
        (TProfile *)sublist_r2->FindObject("Profile_R2EP_Im_TPCFT0C_CentFT0C");
    tp_R2EPBC_Im =
        (TProfile *)sublist_r2->FindObject("Profile_R2EP_Im_FT0AFT0C_CentFT0C");
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
        list_hist_r2 =
            (THashList *)inFile->Get("analysis-event-selection/output");
        sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

        if (first) {
          tp_R2SPAB = (TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_TPCFT0A_CentFT0C");
          tp_R2SPAC = (TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_TPCFT0C_CentFT0C");
          tp_R2SPBC = (TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_FT0AFT0C_CentFT0C");
          tp_R2EPAB = (TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_TPCFT0A_CentFT0C");
          tp_R2EPAC = (TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_TPCFT0C_CentFT0C");
          tp_R2EPBC = (TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_FT0AFT0C_CentFT0C");

          tp_R2SPAB_Im = (TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_Im_TPCFT0A_CentFT0C");
          tp_R2SPAC_Im = (TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_Im_TPCFT0C_CentFT0C");
          tp_R2SPBC_Im = (TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_Im_FT0AFT0C_CentFT0C");
          tp_R2EPAB_Im = (TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_Im_TPCFT0A_CentFT0C");
          tp_R2EPAC_Im = (TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_Im_TPCFT0C_CentFT0C");
          tp_R2EPBC_Im = (TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_Im_FT0AFT0C_CentFT0C");
        } else {
          tp_R2SPAB->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_TPCFT0A_CentFT0C"));
          tp_R2SPAC->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_TPCFT0C_CentFT0C"));
          tp_R2SPBC->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_FT0AFT0C_CentFT0C"));
          tp_R2EPAB->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_TPCFT0A_CentFT0C"));
          tp_R2EPAC->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_TPCFT0C_CentFT0C"));
          tp_R2EPBC->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_FT0AFT0C_CentFT0C"));

          tp_R2SPAB_Im->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_Im_TPCFT0A_CentFT0C"));
          tp_R2SPAC_Im->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_Im_TPCFT0C_CentFT0C"));
          tp_R2SPBC_Im->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2SP_Im_FT0AFT0C_CentFT0C"));
          tp_R2EPAB_Im->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_Im_TPCFT0A_CentFT0C"));
          tp_R2EPAC_Im->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_Im_TPCFT0C_CentFT0C"));
          tp_R2EPBC_Im->Add((TProfile *)sublist_r2->FindObject(
              "Profile_R2EP_Im_FT0AFT0C_CentFT0C"));
        }

        inFile->Close();
        first = false;
      }
      InputFiles.close();
    }
  }
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
    THnSparse *&hs_u2q2_cosDeltaPhi_MEPM1,
    THnSparse *&hs_u2q2_cosDeltaPhi_MEPP1,
    THnSparse *&hs_u2q2_cosDeltaPhi_MEMM1,
    THnSparse *&hs_u2q2_cosDeltaPhi_MEPM2,
    THnSparse *&hs_u2q2_cosDeltaPhi_MEPP2,
    THnSparse *&hs_u2q2_cosDeltaPhi_MEMM2, TH3F *&hs_r2spABMEPM1,
    TH3F *&hs_r2spABMEPP1, TH3F *&hs_r2spABMEMM1, TH3F *&hs_r2spACMEPM1,
    TH3F *&hs_r2spACMEPP1, TH3F *&hs_r2spACMEMM1, TH3F *&hs_r2spBCMEPM1,
    TH3F *&hs_r2spBCMEPP1, TH3F *&hs_r2spBCMEMM1, TH3F *&hs_r2spABMEPM2,
    TH3F *&hs_r2spABMEPP2, TH3F *&hs_r2spABMEMM2, TH3F *&hs_r2spACMEPM2,
    TH3F *&hs_r2spACMEPP2, TH3F *&hs_r2spACMEMM2, TH3F *&hs_r2spBCMEPM2,
    TH3F *&hs_r2spBCMEPP2, TH3F *&hs_r2spBCMEMM2, std::string muonCut,
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

    hs_u2q2_cosDeltaPhi_MEPM1 = (THnSparse *)sublist_v2mepm->FindObject(
        "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1");
    hs_u2q2_cosDeltaPhi_MEPM2 = (THnSparse *)sublist_v2mepm->FindObject(
        "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2");
    hs_u2q2_cosDeltaPhi_MEPP1 = (THnSparse *)sublist_v2mepp->FindObject(
        "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1");
    hs_u2q2_cosDeltaPhi_MEPP2 = (THnSparse *)sublist_v2mepp->FindObject(
        "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2");
    hs_u2q2_cosDeltaPhi_MEMM1 = (THnSparse *)sublist_v2memm->FindObject(
        "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1");
    hs_u2q2_cosDeltaPhi_MEMM2 = (THnSparse *)sublist_v2memm->FindObject(
        "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2");

    hs_r2spABMEPM1 = (TH3F *)sublist_v2mepm->FindObject("R2SPAB1_CentFT0C");
    hs_r2spABMEPM2 = (TH3F *)sublist_v2mepm->FindObject("R2SPAB2_CentFT0C");
    hs_r2spABMEPP1 = (TH3F *)sublist_v2mepp->FindObject("R2SPAB1_CentFT0C");
    hs_r2spABMEPP2 = (TH3F *)sublist_v2mepp->FindObject("R2SPAB2_CentFT0C");
    hs_r2spABMEMM1 = (TH3F *)sublist_v2memm->FindObject("R2SPAB1_CentFT0C");
    hs_r2spABMEMM2 = (TH3F *)sublist_v2memm->FindObject("R2SPAB2_CentFT0C");

    hs_r2spACMEPM1 = (TH3F *)sublist_v2mepm->FindObject("R2SPAC1_CentFT0C");
    hs_r2spACMEPM2 = (TH3F *)sublist_v2mepm->FindObject("R2SPAC2_CentFT0C");
    hs_r2spACMEPP1 = (TH3F *)sublist_v2mepp->FindObject("R2SPAC1_CentFT0C");
    hs_r2spACMEPP2 = (TH3F *)sublist_v2mepp->FindObject("R2SPAC2_CentFT0C");
    hs_r2spACMEMM1 = (TH3F *)sublist_v2memm->FindObject("R2SPAC1_CentFT0C");
    hs_r2spACMEMM2 = (TH3F *)sublist_v2memm->FindObject("R2SPAC2_CentFT0C");

    hs_r2spBCMEPM1 = (TH3F *)sublist_v2mepm->FindObject("R2SPBC1_CentFT0C");
    hs_r2spBCMEPM2 = (TH3F *)sublist_v2mepm->FindObject("R2SPBC2_CentFT0C");
    hs_r2spBCMEPP1 = (TH3F *)sublist_v2mepp->FindObject("R2SPBC1_CentFT0C");
    hs_r2spBCMEPP2 = (TH3F *)sublist_v2mepp->FindObject("R2SPBC2_CentFT0C");
    hs_r2spBCMEMM1 = (TH3F *)sublist_v2memm->FindObject("R2SPBC1_CentFT0C");
    hs_r2spBCMEMM2 = (TH3F *)sublist_v2memm->FindObject("R2SPBC2_CentFT0C");

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
        list_hist_v2se =
            (THashList *)inFile->Get("analysis-same-event-pairing/output");
        list_hist_v2me =
            (THashList *)inFile->Get("analysis-event-mixing/output");

        sublist_v2sepm = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPM_%s", muonCut.c_str())
                            : Form("PairsMuonSEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2sepp = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPP_%s", muonCut.c_str())
                            : Form("PairsMuonSEPP_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2semm = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEMM_%s", muonCut.c_str())
                            : Form("PairsMuonSEMM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2mepm = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEPM_%s", muonCut.c_str())
                            : Form("PairsMuonMEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2mepp = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEPP_%s", muonCut.c_str())
                            : Form("PairsMuonMEPP_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2memm = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEMM_%s", muonCut.c_str())
                            : Form("PairsMuonMEMM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        list_hist_r2 =
            (THashList *)inFile->Get("analysis-event-selection/output");
        sublist_r2 = (TList *)list_hist_r2->FindObject("Event_AfterCuts");

        if (first) {
          hs_V2SEPM =
              (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2");
          hs_V2SEPP =
              (THnSparse *)sublist_v2sepp->FindObject("Mass_Pt_centrFT0C_V2");
          hs_V2SEMM =
              (THnSparse *)sublist_v2semm->FindObject("Mass_Pt_centrFT0C_V2");
          hs_V2MEPM =
              (THnSparse *)sublist_v2mepm->FindObject("Mass_Pt_centrFT0C_V2");
          hs_V2MEPP =
              (THnSparse *)sublist_v2mepp->FindObject("Mass_Pt_centrFT0C_V2");
          hs_V2MEMM =
              (THnSparse *)sublist_v2memm->FindObject("Mass_Pt_centrFT0C_V2");

          hs_R2SPAB = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C");
          hs_R2SPAC = (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C");
          hs_R2SPBC = (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C");

          hs_u2q2_cosDeltaPhi_MEPM1 = (THnSparse *)sublist_v2mepm->FindObject(
              "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1");
          hs_u2q2_cosDeltaPhi_MEPM2 = (THnSparse *)sublist_v2mepm->FindObject(
              "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2");
          hs_u2q2_cosDeltaPhi_MEPP1 = (THnSparse *)sublist_v2mepp->FindObject(
              "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1");
          hs_u2q2_cosDeltaPhi_MEPP2 = (THnSparse *)sublist_v2mepp->FindObject(
              "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2");
          hs_u2q2_cosDeltaPhi_MEMM1 = (THnSparse *)sublist_v2memm->FindObject(
              "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1");
          hs_u2q2_cosDeltaPhi_MEMM2 = (THnSparse *)sublist_v2memm->FindObject(
              "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2");

          hs_r2spABMEPM1 =
              (TH3F *)sublist_v2mepm->FindObject("R2SPAB1_CentFT0C");
          hs_r2spABMEPM2 =
              (TH3F *)sublist_v2mepm->FindObject("R2SPAB2_CentFT0C");
          hs_r2spABMEPP1 =
              (TH3F *)sublist_v2mepp->FindObject("R2SPAB1_CentFT0C");
          hs_r2spABMEPP2 =
              (TH3F *)sublist_v2mepp->FindObject("R2SPAB2_CentFT0C");
          hs_r2spABMEMM1 =
              (TH3F *)sublist_v2memm->FindObject("R2SPAB1_CentFT0C");
          hs_r2spABMEMM2 =
              (TH3F *)sublist_v2memm->FindObject("R2SPAB2_CentFT0C");

          hs_r2spACMEPM1 =
              (TH3F *)sublist_v2mepm->FindObject("R2SPAC1_CentFT0C");
          hs_r2spACMEPM2 =
              (TH3F *)sublist_v2mepm->FindObject("R2SPAC2_CentFT0C");
          hs_r2spACMEPP1 =
              (TH3F *)sublist_v2mepp->FindObject("R2SPAC1_CentFT0C");
          hs_r2spACMEPP2 =
              (TH3F *)sublist_v2mepp->FindObject("R2SPAC2_CentFT0C");
          hs_r2spACMEMM1 =
              (TH3F *)sublist_v2memm->FindObject("R2SPAC1_CentFT0C");
          hs_r2spACMEMM2 =
              (TH3F *)sublist_v2memm->FindObject("R2SPAC2_CentFT0C");

          hs_r2spBCMEPM1 =
              (TH3F *)sublist_v2mepm->FindObject("R2SPBC1_CentFT0C");
          hs_r2spBCMEPM2 =
              (TH3F *)sublist_v2mepm->FindObject("R2SPBC2_CentFT0C");
          hs_r2spBCMEPP1 =
              (TH3F *)sublist_v2mepp->FindObject("R2SPBC1_CentFT0C");
          hs_r2spBCMEPP2 =
              (TH3F *)sublist_v2mepp->FindObject("R2SPBC2_CentFT0C");
          hs_r2spBCMEMM1 =
              (TH3F *)sublist_v2memm->FindObject("R2SPBC1_CentFT0C");
          hs_r2spBCMEMM2 =
              (TH3F *)sublist_v2memm->FindObject("R2SPBC2_CentFT0C");
        } else {
          hs_V2SEPM->Add(
              (THnSparse *)sublist_v2sepm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2SEPP->Add(
              (THnSparse *)sublist_v2sepp->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2SEMM->Add(
              (THnSparse *)sublist_v2semm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2MEPM->Add(
              (THnSparse *)sublist_v2mepm->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2MEPP->Add(
              (THnSparse *)sublist_v2mepp->FindObject("Mass_Pt_centrFT0C_V2"));
          hs_V2MEMM->Add(
              (THnSparse *)sublist_v2memm->FindObject("Mass_Pt_centrFT0C_V2"));

          hs_R2SPAB->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0A_CentFT0C"));
          hs_R2SPAC->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_TPCFT0C_CentFT0C"));
          hs_R2SPBC->Add(
              (TH2F *)sublist_r2->FindObject("R2SP_FT0AFT0C_CentFT0C"));

          hs_u2q2_cosDeltaPhi_MEPM1->Add(
              (THnSparse *)sublist_v2mepm->FindObject(
                  "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1"));
          hs_u2q2_cosDeltaPhi_MEPM2->Add(
              (THnSparse *)sublist_v2mepm->FindObject(
                  "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2"));
          hs_u2q2_cosDeltaPhi_MEPP1->Add(
              (THnSparse *)sublist_v2mepp->FindObject(
                  "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1"));
          hs_u2q2_cosDeltaPhi_MEPP2->Add(
              (THnSparse *)sublist_v2mepp->FindObject(
                  "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2"));
          hs_u2q2_cosDeltaPhi_MEMM1->Add(
              (THnSparse *)sublist_v2memm->FindObject(
                  "Mass_Pt_CentFT0C_U2Q2ev1_cos2DeltaPhiMu1"));
          hs_u2q2_cosDeltaPhi_MEMM2->Add(
              (THnSparse *)sublist_v2memm->FindObject(
                  "Mass_Pt_CentFT0C_U2Q2ev2_cos2DeltaPhiMu2"));

          hs_r2spABMEPM1->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SPAB1_CentFT0C"));
          hs_r2spABMEPM2->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SPAB2_CentFT0C"));
          hs_r2spABMEPP1->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SPAB1_CentFT0C"));
          hs_r2spABMEPP2->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SPAB2_CentFT0C"));
          hs_r2spABMEMM1->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SPAB1_CentFT0C"));
          hs_r2spABMEMM2->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SPAB2_CentFT0C"));

          hs_r2spACMEPM1->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SPAC1_CentFT0C"));
          hs_r2spACMEPM2->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SPAC2_CentFT0C"));
          hs_r2spACMEPP1->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SPAC1_CentFT0C"));
          hs_r2spACMEPP2->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SPAC2_CentFT0C"));
          hs_r2spACMEMM1->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SPAC1_CentFT0C"));
          hs_r2spACMEMM2->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SPAC2_CentFT0C"));

          hs_r2spBCMEPM1->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SPBC1_CentFT0C"));
          hs_r2spBCMEPM2->Add(
              (TH3F *)sublist_v2mepm->FindObject("R2SPBC2_CentFT0C"));
          hs_r2spBCMEPP1->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SPBC1_CentFT0C"));
          hs_r2spBCMEPP2->Add(
              (TH3F *)sublist_v2mepp->FindObject("R2SPBC2_CentFT0C"));
          hs_r2spBCMEMM1->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SPBC1_CentFT0C"));
          hs_r2spBCMEMM2->Add(
              (TH3F *)sublist_v2memm->FindObject("R2SPBC2_CentFT0C"));
        }

        inFile->Close();
        first = false;
      }
      InputFiles.close();
    }
  }
}

//______________________________________________________________________________
void FlowAnalysis_Helper::LoadDataMEProfile(
    std::string FileName, TProfile3D *&tp_V2SEPM, TProfile3D *&tp_V2SEPP,
    TProfile3D *&tp_V2SEMM, TProfile3D *&tp_V2MEPM, TProfile3D *&tp_V2MEPP,
    TProfile3D *&tp_V2MEMM, std::string muonCut, std::string dimuonCut) {
  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist_v2se, *list_hist_v2me;
  TList *sublist_v2sepm, *sublist_v2sepp, *sublist_v2semm;
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

    // Get histograms of v2
    tp_V2SEPM =
        (TProfile3D *)sublist_v2sepm->FindObject("Mass_Pt_CentFT0C_V2SPwR");
    tp_V2SEPP =
        (TProfile3D *)sublist_v2sepp->FindObject("Mass_Pt_CentFT0C_V2SPwR");
    tp_V2SEMM =
        (TProfile3D *)sublist_v2semm->FindObject("Mass_Pt_CentFT0C_V2SPwR");

    tp_V2MEPM =
        (TProfile3D *)sublist_v2mepm->FindObject("Mass_Pt_CentFT0C_V2ME_SP");
    tp_V2MEPP =
        (TProfile3D *)sublist_v2mepp->FindObject("Mass_Pt_CentFT0C_V2ME_SP");
    tp_V2MEMM =
        (TProfile3D *)sublist_v2memm->FindObject("Mass_Pt_CentFT0C_V2ME_SP");

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
        list_hist_v2se =
            (THashList *)inFile->Get("analysis-same-event-pairing/output");
        list_hist_v2me =
            (THashList *)inFile->Get("analysis-event-mixing/output");

        sublist_v2sepm = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPM_%s", muonCut.c_str())
                            : Form("PairsMuonSEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2sepp = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEPP_%s", muonCut.c_str())
                            : Form("PairsMuonSEPP_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2semm = (TList *)list_hist_v2se->FindObject(
            dimuonCut == "" ? Form("PairsMuonSEMM_%s", muonCut.c_str())
                            : Form("PairsMuonSEMM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));

        sublist_v2mepm = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEPM_%s", muonCut.c_str())
                            : Form("PairsMuonMEPM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2mepp = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEPP_%s", muonCut.c_str())
                            : Form("PairsMuonMEPP_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));
        sublist_v2memm = (TList *)list_hist_v2me->FindObject(
            dimuonCut == "" ? Form("PairsMuonMEMM_%s", muonCut.c_str())
                            : Form("PairsMuonMEMM_%s_%s", muonCut.c_str(),
                                   dimuonCut.c_str()));

        if (first) {
          // Get histograms of v2
          tp_V2SEPM = (TProfile3D *)sublist_v2sepm->FindObject(
              "Mass_Pt_CentFT0C_V2SPwR");
          tp_V2SEPP = (TProfile3D *)sublist_v2sepp->FindObject(
              "Mass_Pt_CentFT0C_V2SPwR");
          tp_V2SEMM = (TProfile3D *)sublist_v2semm->FindObject(
              "Mass_Pt_CentFT0C_V2SPwR");

          tp_V2MEPM = (TProfile3D *)sublist_v2mepm->FindObject(
              "Mass_Pt_CentFT0C_V2ME_SP");
          tp_V2MEPP = (TProfile3D *)sublist_v2mepp->FindObject(
              "Mass_Pt_CentFT0C_V2ME_SP");
          tp_V2MEMM = (TProfile3D *)sublist_v2memm->FindObject(
              "Mass_Pt_CentFT0C_V2ME_SP");
        } else {
          // Get histograms of v2
          tp_V2SEPM->Add((TProfile3D *)sublist_v2sepm->FindObject(
              "Mass_Pt_CentFT0C_V2SPwR"));
          tp_V2SEPP->Add((TProfile3D *)sublist_v2sepp->FindObject(
              "Mass_Pt_CentFT0C_V2SPwR"));
          tp_V2SEMM->Add((TProfile3D *)sublist_v2semm->FindObject(
              "Mass_Pt_CentFT0C_V2SPwR"));

          tp_V2MEPM->Add((TProfile3D *)sublist_v2mepm->FindObject(
              "Mass_Pt_CentFT0C_V2ME_SP"));
          tp_V2MEPP->Add((TProfile3D *)sublist_v2mepp->FindObject(
              "Mass_Pt_CentFT0C_V2ME_SP"));
          tp_V2MEMM->Add((TProfile3D *)sublist_v2memm->FindObject(
              "Mass_Pt_CentFT0C_V2ME_SP"));
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

    ey_sys[0] = 0.0031529;
    ey_sys[1] = 0.0038394;
    ey_sys[2] = 0.0040278;
    ey_sys[3] = 0.0039357;
    ey_sys[4] = 0.0033191;
    ey_sys[5] = 0.003598;
    ey_sys[6] = 0.0042368;
    ey_sys[7] = 0.006164;
    ey_sys[8] = 0.0072915;
    ey_sys[9] = 0.010211;
  }
}

//______________________________________________________________________________
void FlowAnalysis_Helper::LoadDataRun2Cent(double *&x, double *&y, double *&ex,
                                           double *&ey, double *&ey_sys,
                                           int flag) {
  x = new double[7];
  y = new double[7];
  ex = new double[7];
  ey = new double[7];
  ey_sys = new double[7];

  if (flag == 0) {
    // 0-5 GeV/c
    x[0] = 5.;
    x[1] = 15.;
    x[2] = 25.;
    x[3] = 35.;
    x[4] = 45.;
    x[5] = 55.;
    x[6] = 75.;

    y[0] = 0.036;
    y[1] = 0.05;
    y[2] = 0.056;
    y[3] = 0.055;
    y[4] = 0.031;
    y[5] = 0.065;
    y[6] = 0.045;

    ex[0] = 5.;
    ex[1] = 5.;
    ex[2] = 5.;
    ex[3] = 5.;
    ex[4] = 5.;
    ex[5] = 5.;
    ex[6] = 15.;

    ey[0] = 0.0093;
    ey[1] = 0.0081;
    ey[2] = 0.0084;
    ey[3] = 0.0099;
    ey[4] = 0.011;
    ey[5] = 0.014;
    ey[6] = 0.019;

    ey_sys[0] = 0.0031087;
    ey_sys[1] = 0.0044159;
    ey_sys[2] = 0.0057615;
    ey_sys[3] = 0.0053348;
    ey_sys[4] = 0.0065608;
    ey_sys[5] = 0.0065192;
    ey_sys[6] = 0.0082189;
  } else {
    // 5-20 GeV/c
    x[0] = 5.;
    x[1] = 15.;
    x[2] = 25.;
    x[3] = 35.;
    x[4] = 45.;
    x[5] = 55.;
    x[6] = 75.;

    y[0] = 0.04;
    y[1] = 0.064;
    y[2] = 0.105;
    y[3] = 0.096;
    y[4] = 0.092;
    y[5] = 0.106;
    y[6] = 0.1;

    ex[0] = 5.;
    ex[1] = 5.;
    ex[2] = 5.;
    ex[3] = 5.;
    ex[4] = 5.;
    ex[5] = 5.;
    ex[6] = 15.;

    ey[0] = 0.02;
    ey[1] = 0.015;
    ey[2] = 0.015;
    ey[3] = 0.016;
    ey[4] = 0.019;
    ey[5] = 0.025;
    ey[6] = 0.036;

    ey_sys[0] = 0.0044385;
    ey_sys[1] = 0.0044843;
    ey_sys[2] = 0.0065419;
    ey_sys[3] = 0.0064786;
    ey_sys[4] = 0.0083954;
    ey_sys[5] = 0.0085108;
    ey_sys[6] = 0.011057;
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