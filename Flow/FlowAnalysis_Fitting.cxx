#include "FlowAnalysis_Fitting.h"

double FlowAnalysis_Fitting::massmin = 0.;
double FlowAnalysis_Fitting::massmax = 10.;
double FlowAnalysis_Fitting::centmin = 0.;
double FlowAnalysis_Fitting::centmax = 0.;
int FlowAnalysis_Fitting::mflag_sig = 0.;
int FlowAnalysis_Fitting::mflag_bkg = 0.;
int FlowAnalysis_Fitting::mflag_bkg_v2 = 0.;
int FlowAnalysis_Fitting::norder = 2;
int FlowAnalysis_Fitting::nhar = 2;
int FlowAnalysis_Fitting::mode = 0;
string FlowAnalysis_Fitting::mode_string[2] = {"Std", "Sys"};

//______________________________________________________________________________
void FlowAnalysis_Fitting::init() {
  FlowAnalysis_Fitting::mflag_sig = 0;
  FlowAnalysis_Fitting::mflag_bkg = 1;
  mchi2max = 2.0;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setModel(int flag_sig, int flag_bkg) {
  FlowAnalysis_Fitting::mflag_sig = flag_sig;
  FlowAnalysis_Fitting::mflag_bkg = flag_bkg;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setModelV2(int flag_bkg_v2) {
  FlowAnalysis_Fitting::mflag_bkg_v2 = flag_bkg_v2;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setMassRange(double mass_min, double mass_max) {
  FlowAnalysis_Fitting::massmax = mass_max;
  FlowAnalysis_Fitting::massmin = mass_min;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setCentRange(double cent_min, double cent_max) {
  FlowAnalysis_Fitting::centmax = cent_max;
  FlowAnalysis_Fitting::centmin = cent_min;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setOrder(int order) {
  FlowAnalysis_Fitting::norder = order;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setHarmonic(int har) {
  FlowAnalysis_Fitting::nhar = har;
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::DoubleSidedCB2(double x, double mu, double width,
                                            double a1, double p1, double a2,
                                            double p2) {
  double u = (x - mu) / width;
  double A1 = pow(p1 / abs(a1), p1) * exp(-a1 * a1 / 2);
  double A2 = pow(p2 / abs(a2), p2) * exp(-a2 * a2 / 2);
  double B1 = p1 / abs(a1) - abs(a1);
  double B2 = p2 / abs(a2) - abs(a2);

  double result = 0.;
  if (u < -a1) {
    result = A1 * pow(B1 - u, -p1);
  } else if (u >= -a1 && u <= a2) {
    result = exp(-u * u / 2);
  } else {
    result = A2 * pow(B2 + u, -p2);
  }
  return result;
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::DoubleSidedCB(double *x, double *par) {
  return FlowAnalysis_Fitting::DoubleSidedCB2(x[0], par[0], par[1], par[2],
                                              par[3], par[4], par[5]);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::Cheby7(double *x, double *par) {
  double t =
      -1. + 2. * (x[0] - FlowAnalysis_Fitting::massmin) /
                (FlowAnalysis_Fitting::massmax - FlowAnalysis_Fitting::massmin);
  double T0 = 1;
  double T1 = t;
  double T2 = 2. * pow(t, 2.) - 1.;
  double T3 = 4. * pow(t, 3.) - 3. * t;
  double T4 = 8. * pow(t, 4.) - 8. * pow(t, 2.) + 1;
  double T5 = 16. * pow(t, 5.) - 20. * pow(t, 3.) + 5. * t;
  double T6 = 32. * pow(t, 6.) - 48. * pow(t, 4.) + 18. * pow(t, 2) - 1;
  double T7 = 2. * t * T6 - T5;
  return (T0 + par[0] * T1 + par[1] * T2 + par[2] * T3 + par[3] * T4 +
          par[4] * T5 + par[5] * T6 + par[6] * T7);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::Cheby3(double *x, double *par) {
  double t =
      -1. + 2. * (x[0] - FlowAnalysis_Fitting::massmin) /
                (FlowAnalysis_Fitting::massmax - FlowAnalysis_Fitting::massmin);
  double T0 = 1;
  double T1 = t;
  double T2 = 2. * pow(t, 2.) - 1.;
  double T3 = 4. * pow(t, 3.) - 3. * t;
  return (T0 + par[0] * T1 + par[1] * T2 + par[2] * T3);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::VariableWidthGauss(double *x, double *par) {
  double sigma = par[1] + par[2] * (x[0] - par[0]) / par[0];
  double mu = (x[0] - par[0]) / sigma;
  return exp(-1. / 2 * pow(mu, 2.));
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::DoubleExp(double *x, double *par) {
  return par[0] * exp(par[1] * (x[0] - par[2])) +
         par[3] * exp(par[4] * (x[0] - par[5]));
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::PolyExp(double *x, double *par) {
  return (par[0] + par[1] * x[0] + par[3] * pow(x[0], 2.) +
          par[4] * pow(x[0], 3.) + par[5] * pow(x[0], 4.)) *
         exp(par[2] * x[0]);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::FittedSignal(double *x, double *par) {
  double value = 0.;
  if (FlowAnalysis_Fitting::mflag_sig == CB2) {
    value = FlowAnalysis_Fitting::DoubleSidedCB(x, &par[2]);
  }
  return par[0] * value / par[1];
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::FittedBkg(double *x, double *par) {
  double value = 0.;
  if (FlowAnalysis_Fitting::mflag_bkg == Chebychev) {
    value = FlowAnalysis_Fitting::Cheby7(x, &par[2]);
  }
  if (FlowAnalysis_Fitting::mflag_bkg == PolExp) {
    value = FlowAnalysis_Fitting::DoubleExp(x, &par[2]);
  }
  if (FlowAnalysis_Fitting::mflag_bkg == VWG) {
    value = FlowAnalysis_Fitting::VariableWidthGauss(x, &par[2]);
  }
  if (FlowAnalysis_Fitting::mflag_bkg == PolExp) {
    value = FlowAnalysis_Fitting::PolyExp(x, &par[2]);
  }
  return par[0] * value / par[1];
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::CreateModel(TF1 *&model, int flag) {
  if (flag == CB2) {
    model = new TF1("CB2", FlowAnalysis_Fitting::DoubleSidedCB,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 6);
    model->SetParameter(0, 3.097);      // mean
    model->SetParLimits(0, 3.05, 3.13); // mean
    model->SetParName(0, "#mu");        // mean
    model->SetParameter(1, 0.08);       // sigma
    model->SetParLimits(1, 0.05, 0.12); // sigma
    model->SetParName(1, "#sigma");     // sigma
    model->SetParameter(2, 0.883);      // alphaL
    model->SetParName(2, "#alpha_{L}"); // alphaL
    model->SetParameter(3, 9.940);      // nL
    model->SetParName(3, "n_{L}");      // nL
    model->SetParameter(4, 1.832);      // alphaR
    model->SetParName(4, "#alpha_{R}"); // alphaR
    model->SetParameter(5, 15.323);     // nR
    model->SetParName(5, "n_{R}");      // nR
  }
  if (flag == VWG) {
    model = new TF1("VWG2", FlowAnalysis_Fitting::VariableWidthGauss,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 3);
    model->SetParameter(0, 0.1);     // a0
    model->SetParLimits(0, -1., 5.); // a0
    model->SetParName(0, "a0");      // a0
    model->SetParameter(1, 0.2);     // a1
    model->SetParLimits(1, 0., 3.);  // a1
    model->SetParName(1, "a1");      // a1
    model->SetParameter(2, 0.1);     // a2
    model->SetParLimits(2, -1., 1.); // a2
    model->SetParName(2, "a2");      // a2
  }
  if (flag == PolExp) {
    model = new TF1("Exp2", FlowAnalysis_Fitting::PolyExp,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 6);
    model->SetParameter(0, 1.);         // a0
    model->SetParLimits(0, -20., 500.); // a0
    model->SetParName(0, "a0");         // a0
    model->SetParameter(1, 1.);         // a1
    model->SetParLimits(1, -10., 10.);  // a1
    model->SetParName(1, "a1");         // a1
    model->SetParameter(2, -1.0);       // a2
    model->SetParLimits(2, -10., 10.);  // a2
    model->SetParName(2, "a2");         // a2
    model->SetParameter(3, 0.);         // a3
    model->SetParLimits(3, -10., 10.);  // a3
    model->SetParName(3, "a3");         // a3
    model->SetParameter(4, 0.);         // a4
    model->SetParLimits(4, -10., 10.);  // a4
    model->SetParName(4, "a4");         // a4
    model->SetParameter(5, 0.);         // a5
    model->SetParLimits(5, -10., 10.);  // a5
    model->SetParName(5, "a5");         // a5
  }
  if (flag == Exp2) {
    model = new TF1("Exp2", FlowAnalysis_Fitting::DoubleExp,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 6);
    model->SetParameter(0, 0.0);       // a0
    model->SetParLimits(0, -10., 10.); // a0
    model->SetParName(0, "a0");        // a0
    model->SetParameter(1, -0.1);      // a1
    model->SetParLimits(1, -5., 5.);   // a1
    model->SetParName(1, "a1");        // a1
    model->SetParameter(2, 0.1);       // a2
    model->SetParLimits(2, -5., 5.);   // a2
    model->SetParName(2, "a2");        // a2
    model->SetParameter(3, 0.0);       // a3
    model->SetParLimits(3, -10., 10.); // a3
    model->SetParName(3, "a3");        // a3
    model->SetParameter(4, -0.2);      // a4
    model->SetParLimits(4, -5., 5.);   // a4
    model->SetParName(4, "a4");        // a4
    model->SetParameter(5, 0.2);       // a5
    model->SetParLimits(5, -5., 5.);   // a5
    model->SetParName(5, "a5");        // a5
  }
  if (flag == Chebychev) {
    model = new TF1("Chebychev", FlowAnalysis_Fitting::Cheby7,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 7);
    model->SetParameter(0, 0.0);     // a0
    model->SetParLimits(0, -8., 8.); // a0
    model->SetParName(0, "a0");      // a0
    model->SetParameter(1, 0.0);     // a1
    model->SetParLimits(1, -2., 2.); // a1
    model->SetParName(1, "a1");      // a1
    model->SetParameter(2, 0.0);     // a2
    model->SetParLimits(2, -2., 2.); // a2
    model->SetParName(2, "a2");      // a2
    model->SetParameter(3, 0.0);     // a3
    model->SetParLimits(3, -2., 2.); // a3
    model->SetParName(3, "a3");      // a3
    model->SetParameter(4, 0.0);     // a4
    model->SetParLimits(4, -2., 2.); // a4
    model->SetParName(4, "a4");      // a4
    model->SetParameter(5, 0.0);     // a5
    model->SetParLimits(5, -2., 2.); // a5
    model->SetParName(5, "a5");      // a5
    model->SetParameter(6, 0.0);     // a6
    model->SetParLimits(6, -2., 2.); // a6
    model->SetParName(6, "a6");      // a6
  }
}

//______________________________________________________________________________
vector<double> FlowAnalysis_Fitting::runFitting(TH1D *hs_input,
                                                TH1D *hs_v2_input, TList *ls,
                                                double ptmin, double ptmax) {
  vector<double> results;
  // Fitting for dimuon invariant mass + v2 signal
  cout << ">>>>>>>>>>>>>> Start processing Pt range: " << ptmin << " " << ptmax
       << endl;
  // Make copies for input histograms for safety
  TH1D *hs = dynamic_cast<TH1D *>(hs_input->Clone());
  TH1D *hs_v2 = dynamic_cast<TH1D *>(hs_v2_input->Clone());
  //////////////////////////////////////////////////////////////////////////////
  ///      INVARIANT MASS FIT
  //////////////////////////////////////////////////////////////////////////////
  // Setting up fit model for invariant mass fit
  cout << ">>>>>>>>>>>>>> Start processing dimuon invariant mass fit..."
       << endl;
  // Construct a combined signal+background model
  TF1 *sig, *bkg;
  CreateModel(sig, ModelType(FlowAnalysis_Fitting::mflag_sig));
  CreateModel(bkg, ModelType(FlowAnalysis_Fitting::mflag_bkg));
  int nsig = 5.E2;
  int nbkg = 5.E5;
  TF1NormSum *sum_model = new TF1NormSum(sig, bkg, nsig, nbkg);
  TF1 *model = new TF1("model", *sum_model, FlowAnalysis_Fitting::massmin,
                       FlowAnalysis_Fitting::massmax, sum_model->GetNpar());
  model->SetParameters(sum_model->GetParameters().data());
  model->SetParName(0, "N_{J/#psi}");
  model->SetParName(1, "N_{bkg}");
  model->SetParLimits(0, 1., 1.E12);
  model->SetParLimits(1, 1., 1.E12);
  int nPar_model = model->GetNpar();
  int nPar_sig = sig->GetNpar();
  int nPar_bkg = bkg->GetNpar();
  for (int i = 2; i < nPar_model; i++) {
    if (i < 2 + nPar_sig) {
      model->SetParName(i, sig->GetParName(i - 2));
      double min, max;
      sig->GetParLimits(i - 2, min, max);
      model->SetParLimits(i, min, max);
      if (i - 2 >= 2) {
        model->FixParameter(i, model->GetParameter(i));
      }
    } else {
      model->SetParName(i, bkg->GetParName(i - 2 - nPar_sig));
      double min, max;
      bkg->GetParLimits(i - 2 - nPar_sig, min, max);
      model->SetParLimits(i, min, max);
      if (i > 2 + nPar_sig) {
        model->FixParameter(i, model->GetParameter(i));
      }
    }
  }

  // Do fitting for invariant mass
  // Configuration for fitting
  MinimizerOptions::SetDefaultMinimizer("Minuit2");

  // Iterative fitting
  hs->Scale(1., "width");
  int nfree_bkg = model->GetNumberFreeParameters() - 4;
  for (int i = 0; i <= nPar_bkg - nfree_bkg; i++) {
    auto result = hs->Fit("model", "SQ0");
    result->Print();
    double chi2ndf = model->GetChisquare() / model->GetNDF();
    cout << "chi2/ndf: " << chi2ndf << endl;
    if (chi2ndf > mchi2max && i != nPar_bkg) {
      model->ReleaseParameter(i + 9);
    } else {
      break;
    }
  }

  // Getting components of fitted model
  int nsig_fitted = model->GetParameter(0);
  int nbkg_fitted = model->GetParameter(1);

  // Fitted signal function
  TF1 *sig_fitted = new TF1("sig_fitted", FlowAnalysis_Fitting::FittedSignal,
                            FlowAnalysis_Fitting::massmin,
                            FlowAnalysis_Fitting::massmax, nPar_sig + 2);
  for (int i = 0; i < nPar_sig + 2; i++) {
    if (i == 0) {
      sig_fitted->SetParameter(i, model->GetParameter(0));
    } else if (i == 1) {
      sig_fitted->SetParameter(i, 1.0);
    } else {
      sig_fitted->SetParameter(i, model->GetParameter(i));
    }
  }
  double int_sig = sig_fitted->Integral(FlowAnalysis_Fitting::massmin,
                                        FlowAnalysis_Fitting::massmax) /
                   model->GetParameter(0);
  sig_fitted->SetParameter(1, int_sig);

  // Fitted background function
  TF1 *bkg_fitted = new TF1("bkg_fitted", FlowAnalysis_Fitting::FittedBkg,
                            FlowAnalysis_Fitting::massmin,
                            FlowAnalysis_Fitting::massmax, nPar_bkg + 2);
  for (int i = 0; i < nPar_bkg + 2; i++) {
    if (i == 0) {
      bkg_fitted->SetParameter(i, model->GetParameter(1));
    } else if (i == 1) {
      bkg_fitted->SetParameter(i, 1.0);
    } else {
      bkg_fitted->SetParameter(i, model->GetParameter(i + nPar_sig));
    }
  }
  double int_bkg = bkg_fitted->Integral(FlowAnalysis_Fitting::massmin,
                                        FlowAnalysis_Fitting::massmax) /
                   model->GetParameter(1);
  bkg_fitted->SetParameter(1, int_bkg);

  // Plotting
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  TCanvas *c = new TCanvas(
      Form(
          "Fitted_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%d_%d_%d",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder, ptmin,
          ptmax, FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::mflag_sig, FlowAnalysis_Fitting::mflag_bkg,
          FlowAnalysis_Fitting::mflag_bkg_v2),
      Form(
          "Fitted_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%d_%d_%d",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder, ptmin,
          ptmax, FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::mflag_sig, FlowAnalysis_Fitting::mflag_bkg,
          FlowAnalysis_Fitting::mflag_bkg_v2));
  c->cd();
  hs->SetMarkerStyle(20);
  hs->SetMarkerSize(0.8);
  hs->SetTitle(Form("J/#psi invariant mass: %g - %g GeV/c, %g - %g %%", ptmin,
                    ptmax, FlowAnalysis_Fitting::centmin,
                    FlowAnalysis_Fitting::centmax));
  hs->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c2)");
  hs->GetYaxis()->SetTitle(Form("Events / ( %g GeV/c )", hs->GetBinWidth(1)));
  hs->Draw("HIST EP");
  gPad->ModifiedUpdate();
  model->SetLineWidth(5.0);
  model->SetLineColor(kBlue);
  model->Draw("same");
  gPad->ModifiedUpdate();
  bkg_fitted->SetLineWidth(5.0);
  bkg_fitted->SetLineColor(kBlue);
  bkg_fitted->SetLineStyle(kDashed);
  bkg_fitted->Draw("same");
  gPad->ModifiedUpdate();
  sig_fitted->SetLineWidth(5.0);
  sig_fitted->SetLineColor(kRed);
  sig_fitted->Draw("same");
  gPad->ModifiedUpdate();

  double mean_fitted = model->GetParameter(2);
  double sigma_fitted = model->GetParameter(3);
  double S_3sigma = sig_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                         mean_fitted + 3. * sigma_fitted);
  double B_3sigma = bkg_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                         mean_fitted + 3. * sigma_fitted);

  TPaveStats *sb = (TPaveStats *)gPad->GetPrimitive("stats");
  sb->SetName("J/#psi invariant mass fit");
  sb->SetX1NDC(0.7);
  sb->SetX2NDC(0.95);
  sb->SetY1NDC(0.88);
  sb->SetY2NDC(0.38);
  sb->AddText(
      TString::Format("%s = %f", "(S/B)_{3#sigma}", S_3sigma / B_3sigma));
  hs->SetStats(0);
  sb->Draw();
  gPad->ModifiedUpdate();

  TLatex *text_info = new TLatex();
  text_info->SetTextSize(0.04);
  text_info->SetTextFont(42);
  text_info->DrawLatexNDC(
      .18, .82, "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  gPad->ModifiedUpdate();
  text_info->DrawLatexNDC(.18, .77,
                          "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");
  gPad->ModifiedUpdate();
  text_info->DrawLatexNDC(.18, .72,
                          Form("%g < #it{p}_{T} < %g GeV/c", ptmin, ptmax));
  gPad->ModifiedUpdate();
  ls->Add(c);

  //////////////////////////////////////////////////////////////////////////////
  ///      V2 FIT
  //////////////////////////////////////////////////////////////////////////////
  // Setting up fit model for v2 signal fit
  cout << ">>>>>>>>>>>>>> Start processing v2 fit..." << endl;

  // Constructing alpha function
  auto fct_alpha = [sig_fitted, bkg_fitted](double *x, double *) {
    double S = sig_fitted->Eval(x[0]);
    double B = bkg_fitted->Eval(x[0]);
    return S / (S + B);
  };
  TF1 *alpha = new TF1("alpha", fct_alpha, FlowAnalysis_Fitting::massmin,
                       FlowAnalysis_Fitting::massmax, 0);

  // Constructing combined v2 fit function
  auto fct_v2 = [alpha](double *x, double *par) {
    double val_alpha = alpha->Eval(x[0]);
    double value = 0.;
    if (FlowAnalysis_Fitting::mflag_bkg_v2 == 0) {
      // Using Pol2
      value = (par[0] * val_alpha + (1 - val_alpha) * (par[1] + par[2] * x[0] +
                                                       par[3] * x[0] * x[0]));
    }
    if (FlowAnalysis_Fitting::mflag_bkg_v2 == 1) {
      // Using Chebychev
      value = (par[0] * val_alpha +
               (1 - val_alpha) * FlowAnalysis_Fitting::Cheby3(x, &par[1]));
    }
    return value;
  };
  TF1 *model_v2 = new TF1("model_v2", fct_v2, FlowAnalysis_Fitting::massmin,
                          FlowAnalysis_Fitting::massmax, 4);
  model_v2->SetParameter(0, 0.);
  model_v2->SetParLimits(0, -1., 1.);
  model_v2->SetParName(0, "v^{J/#psi}_{2}");
  model_v2->SetParameter(1, 0.);
  model_v2->SetParLimits(1, -5., 5.);
  model_v2->SetParName(1, "a0");
  model_v2->SetParameter(2, 0.);
  model_v2->SetParLimits(2, -2., 2.);
  model_v2->SetParName(2, "a2");
  model_v2->SetParameter(3, 0.);
  model_v2->SetParLimits(3, -2., 2.);
  model_v2->SetParName(3, "a3");

  // Do fitting
  hs_v2->Rebin();
  hs_v2->Scale(0.5);
  auto result_v2 = hs_v2->Fit("model_v2", "SQ0");
  result_v2->Print();
  double chi2ndf_v2 = model_v2->GetChisquare() / model_v2->GetNDF();
  cout << "chi2/ndf: " << chi2ndf_v2 << endl;

  // Getting fitted background
  auto fitted_v2bkg = [](double *x, double *par) {
    return (par[0] + par[1] * x[0] + par[2] * x[0] * x[0]);
  };
  TF1 *v2bkg_fitted =
      new TF1("v2bkg", fitted_v2bkg, FlowAnalysis_Fitting::massmin,
              FlowAnalysis_Fitting::massmax, 3);
  v2bkg_fitted->SetParameter(0, model_v2->GetParameter(1));
  v2bkg_fitted->SetParameter(1, model_v2->GetParameter(2));
  v2bkg_fitted->SetParameter(2, model_v2->GetParameter(3));
  // Do plotting
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  TCanvas *c_v2 = new TCanvas(
      Form(
          "Fitted_%s_v%d%d_%g_%g_%g_%g_%g_%g_%d_%d_%d",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder, ptmin,
          ptmax, FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::mflag_sig, FlowAnalysis_Fitting::mflag_bkg,
          FlowAnalysis_Fitting::mflag_bkg_v2),
      Form(
          "Fitted_%s_v%d%d_%g_%g_%g_%g_%g_%g_%d_%d_%d",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder, ptmin,
          ptmax, FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::mflag_sig, FlowAnalysis_Fitting::mflag_bkg,
          FlowAnalysis_Fitting::mflag_bkg_v2));
  c_v2->cd();
  hs_v2->SetMarkerStyle(20);
  hs_v2->SetMarkerSize(0.8);
  hs_v2->SetTitle(Form("J/#psi #it{v}_{2}{%d}: %g - %g GeV/c, %g - %g %%",
                       FlowAnalysis_Fitting::norder, ptmin, ptmax,
                       FlowAnalysis_Fitting::centmin,
                       FlowAnalysis_Fitting::centmax));
  hs_v2->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c2)");
  hs_v2->GetYaxis()->SetTitle("#it{v}^{#mu#mu}_{2}");
  hs_v2->Draw("HIST EP");
  gPad->ModifiedUpdate();
  model_v2->SetLineWidth(5.0);
  model_v2->SetLineColor(kBlue);
  model_v2->Draw("same");
  gPad->ModifiedUpdate();
  v2bkg_fitted->SetLineWidth(5.0);
  v2bkg_fitted->SetLineColor(kBlue);
  v2bkg_fitted->SetLineStyle(kDashed);
  v2bkg_fitted->Draw("same");
  gPad->ModifiedUpdate();
  TPaveStats *sb_v2 = (TPaveStats *)gPad->GetPrimitive("stats");
  sb_v2->SetName("J/#psi v2 fit");
  sb_v2->SetX1NDC(0.7);
  sb_v2->SetX2NDC(0.95);
  sb_v2->SetY1NDC(0.88);
  sb_v2->SetY2NDC(0.38);
  hs_v2->SetStats(0);
  sb_v2->Draw();
  gPad->ModifiedUpdate();
  ls->Add(c_v2);
  cout << endl;
  cout << endl;
  results.emplace_back(model_v2->GetParameter(0));
  results.emplace_back(model_v2->GetParError(0));
  results.emplace_back(model->GetParameter(0));
  results.emplace_back(model->GetParError(0));
  return results;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::Print() {
  cout << "Info for fitter: " << endl;
  cout << "Signal model flag: " << FlowAnalysis_Fitting::mflag_sig << endl;
  cout << "Background model flag: " << FlowAnalysis_Fitting::mflag_bkg << endl;
  cout << "Chi2 threshold: " << mchi2max << endl;
}