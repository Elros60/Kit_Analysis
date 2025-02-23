#include "FlowAnalysis_Fitting.h"

double FlowAnalysis_Fitting::massmin = 0.;
double FlowAnalysis_Fitting::massmax = 10.;
double FlowAnalysis_Fitting::centmin = 0.;
double FlowAnalysis_Fitting::centmax = 0.;
double FlowAnalysis_Fitting::ptmin = 0.;
double FlowAnalysis_Fitting::ptmax = 10.;
int FlowAnalysis_Fitting::mflag_sig = 0.;
int FlowAnalysis_Fitting::mflag_bkg = 0.;
int FlowAnalysis_Fitting::mflag_bkg_v2 = 0.;
int FlowAnalysis_Fitting::norder = 2;
int FlowAnalysis_Fitting::nhar = 2;
int FlowAnalysis_Fitting::mode = 0;
string FlowAnalysis_Fitting::mode_string[2] = {"Std", "Sys"};
string FlowAnalysis_Fitting::model_string[8] = {
    "CB2(data)",   "CB2(MC)", "NA60", "Cheby",
    "EventMixing", "VWG",     "Exp2", "PolExp"};
string FlowAnalysis_Fitting::v2bkg_string[2] = {"EventMixing(beta free)",
                                                "EventMixing(beta fix)"};

//______________________________________________________________________________
void FlowAnalysis_Fitting::init() {
  FlowAnalysis_Fitting::mflag_sig = 0;
  FlowAnalysis_Fitting::mflag_bkg = 1;
  mchi2max_mass = 2.0;
  mchi2max_v2 = 2.0;
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
void FlowAnalysis_Fitting::setPtRange(double pt_min, double pt_max) {
  FlowAnalysis_Fitting::ptmax = pt_max;
  FlowAnalysis_Fitting::ptmin = pt_min;
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
double FlowAnalysis_Fitting::NA60Function(double *x, double *par) {
  double t = (x[0] - par[0]) / par[1];
  double t0 = 0;
  if (t <= par[2]) {
    t0 =
        1. + pow(par[4] * (par[2] - t), par[5] - par[6] * pow(par[2] - t, 0.5));
  } else if (t >= par[3]) {
    t0 =
        1. + pow(par[7] * (t - par[3]), par[8] - par[9] * pow(t - par[3], 0.5));
  } else {
    t0 = 1.;
  }
  return exp(-0.5 * pow(t / t0, 2));
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
  return (par[0] * T0 + par[1] * T1 + par[2] * T2 + par[3] * T3);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::Cheby4(double *x, double *par) {
  double t =
      -1. + 2. * (x[0] - FlowAnalysis_Fitting::massmin) /
                (FlowAnalysis_Fitting::massmax - FlowAnalysis_Fitting::massmin);
  double T0 = 1;
  double T1 = t;
  double T2 = 2. * pow(t, 2.) - 1.;
  double T3 = 4. * pow(t, 3.) - 3. * t;
  double T4 = 8. * pow(t, 4.) - 8. * pow(t, 2.) + 1;
  return (par[0] * T0 + par[1] * T1 + par[2] * T2 + par[3] * T3 + par[4] * T4);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::Cheby5(double *x, double *par) {
  double t =
      -1. + 2. * (x[0] - FlowAnalysis_Fitting::massmin) /
                (FlowAnalysis_Fitting::massmax - FlowAnalysis_Fitting::massmin);
  double T0 = 1;
  double T1 = t;
  double T2 = 2. * pow(t, 2.) - 1.;
  double T3 = 4. * pow(t, 3.) - 3. * t;
  double T4 = 8. * pow(t, 4.) - 8. * pow(t, 2.) + 1;
  double T5 = 16. * pow(t, 5.) - 20. * pow(t, 3.) + 5. * t;
  return (par[0] * T0 + par[1] * T1 + par[2] * T2 + par[3] * T3 + par[4] * T4 +
          par[5] * T5);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::Cheby6(double *x, double *par) {
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
  return (par[0] * T0 + par[1] * T1 + par[2] * T2 + par[3] * T3 + par[4] * T4 +
          par[5] * T5 + par[6] * T6);
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::VariableWidthGauss(double *x, double *par) {
  double sigma = par[1] + par[2] * (x[0] - par[0]) / par[0] +
                 par[3] * pow((x[0] - par[0]) / par[0], 2.);
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
  if (FlowAnalysis_Fitting::mflag_sig == CB2Data ||
      FlowAnalysis_Fitting::mflag_sig == CB2MC) {
    value = FlowAnalysis_Fitting::DoubleSidedCB(x, &par[2]);
  }
  if (FlowAnalysis_Fitting::mflag_sig == NA60) {
    value = FlowAnalysis_Fitting::NA60Function(x, &par[2]);
  }
  return par[0] * value / par[1];
}

//______________________________________________________________________________
double FlowAnalysis_Fitting::FittedBkg(double *x, double *par) {
  double value = 0.;
  if (FlowAnalysis_Fitting::mflag_bkg == Chebychev) {
    value = FlowAnalysis_Fitting::Cheby7(x, &par[2]);
  }
  if (FlowAnalysis_Fitting::mflag_bkg == Exp2) {
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
  if (flag == CB2Data) {
    model = new TF1("CB2Data", FlowAnalysis_Fitting::DoubleSidedCB,
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
  if (flag == CB2MC) {
    model = new TF1("CB2MC", FlowAnalysis_Fitting::DoubleSidedCB,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 6);
    model->SetParameter(0, 3.097);      // mean
    model->SetParLimits(0, 3.05, 3.13); // mean
    model->SetParName(0, "#mu");        // mean
    model->SetParameter(1, 0.08);       // sigma
    model->SetParLimits(1, 0.05, 0.12); // sigma
    model->SetParName(1, "#sigma");     // sigma
    model->SetParameter(2, 0.993);      // alphaL
    // model->SetParameter(2, 1.078);      // alphaL
    model->SetParName(2, "#alpha_{L}"); // alphaL
    model->SetParameter(3, 2.9075);     // nL
    // model->SetParameter(3, 2.402); // nL
    model->SetParName(3, "n_{L}"); // nL
    model->SetParameter(4, 2.182); // alphaR
    // model->SetParameter(4, 2.283);      // alphaR
    model->SetParName(4, "#alpha_{R}"); // alphaR
    model->SetParameter(5, 3.122);      // nR
    // model->SetParameter(5, 5.268); // nR
    model->SetParName(5, "n_{R}"); // nR
  }
  if (flag == NA60) {
    model = new TF1("NA60", FlowAnalysis_Fitting::NA60Function,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 10);
    model->SetParameter(0, 3.097);      // mean
    model->SetParLimits(0, 3.05, 3.13); // mean
    model->SetParName(0, "#mu");        // mean
    model->SetParameter(1, 0.08);       // sigma
    model->SetParLimits(1, 0.05, 0.12); // sigma
    model->SetParName(1, "#sigma");     // sigma
    model->SetParameter(2, -0.6);       // t1
    model->SetParName(2, "t1");         // t1
    model->SetParameter(3, 1.8);        // t2
    model->SetParName(3, "t2");         // t2
    model->SetParameter(4, 0.2);        // p1
    model->SetParName(4, "p1");         // p1
    model->SetParameter(5, 1.5);        // p2
    model->SetParName(5, "p2");         // p2
    model->SetParameter(6, 0.1);        // p3
    model->SetParName(6, "p3");         // p3
    model->SetParameter(7, 0.18);       // p4
    model->SetParName(7, "p4");         // p4
    model->SetParameter(8, 1.54);       // p5
    model->SetParName(8, "p5");         // p5
    model->SetParameter(9, 0.17);       // p6
    model->SetParName(9, "p6");         // p6
  }
  if (flag == VWG) {
    model = new TF1("VWG", FlowAnalysis_Fitting::VariableWidthGauss,
                    FlowAnalysis_Fitting::massmin,
                    FlowAnalysis_Fitting::massmax, 4);
    double bin_center =
        (FlowAnalysis_Fitting::ptmin + FlowAnalysis_Fitting::ptmax) / 2;
    if (bin_center >= 0 && bin_center < 2) {
      model->SetParameter(0, 2.5);      // a0
      model->SetParLimits(0, -2.5, 5.); // a0
      model->SetParName(0, "a0");       // a0
      model->SetParameter(1, 1.);       // a1
      model->SetParLimits(1, 0., 5.);   // a1
      model->SetParName(1, "a1");       // a1
      model->SetParameter(2, 0.);       // a2
      model->SetParLimits(2, -2., 2.);  // a2
      model->SetParName(2, "a2");       // a2
      model->SetParameter(3, 0.);       // a2
      model->SetParLimits(3, -1., 1.);  // a3
      model->SetParName(3, "a3");       // a3
    } else if (bin_center >= 2 && bin_center < 3) {
      model->SetParameter(0, 1.5);      // a0
      model->SetParLimits(0, -2.0, 5.); // a0
      model->SetParName(0, "a0");       // a0
      model->SetParameter(1, 2.);       // a1
      model->SetParLimits(1, 0.5, 5.);  // a1
      model->SetParName(1, "a1");       // a1
      model->SetParameter(2, 0.);       // a2
      model->SetParLimits(2, -2., 2.);  // a2
      model->SetParName(2, "a2");       // a2
      model->SetParameter(3, 0.);       // a2
      model->SetParLimits(3, -1., 1.);  // a3
      model->SetParName(3, "a3");       // a3
    } else if (bin_center >= 3 && bin_center < 5) {
      model->SetParameter(0, 2.);       // a0
      model->SetParLimits(0, 0.1, 3.2); // a0
      model->SetParName(0, "a0");       // a0
      model->SetParameter(1, 0.65);     // a1
      model->SetParLimits(1, 0.5, 2.);  // a1
      model->SetParName(1, "a1");       // a1
      model->SetParameter(2, 0.);       // a2
      model->SetParLimits(2, -2., 2.);  // a2
      model->SetParName(2, "a2");       // a2
      model->SetParameter(3, 0.);       // a2
      model->SetParLimits(3, -1., 1.);  // a3
      model->SetParName(3, "a3");       // a3
    } else if (bin_center >= 5 && bin_center < 8) {
      model->SetParameter(0, 2.3);      // a0
      model->SetParLimits(0, 0.1, 3.2); // a0
      model->SetParName(0, "a0");       // a0
      model->SetParameter(1, 1.);       // a1
      model->SetParLimits(1, 0.5, 2.);  // a1
      model->SetParName(1, "a1");       // a1
      model->SetParameter(2, 0.);       // a2
      model->SetParLimits(2, -2., 2.);  // a2
      model->SetParName(2, "a2");       // a2
      model->SetParameter(3, 0.);       // a2
      model->SetParLimits(3, -1., 1.);  // a3
      model->SetParName(3, "a3");       // a3
    } else if (bin_center >= 8 && bin_center < 12) {
      model->SetParameter(0, 3.);      // a0
      model->SetParLimits(0, 2., 5.);  // a0
      model->SetParName(0, "a0");      // a0
      model->SetParameter(1, 2.5);     // a1
      model->SetParLimits(1, 1.5, 4.); // a1
      model->SetParName(1, "a1");      // a1
      model->SetParameter(2, 0.);      // a2
      model->SetParLimits(2, -2., 2.); // a2
      model->SetParName(2, "a2");      // a2
      model->SetParameter(3, 0.);      // a2
      model->SetParLimits(3, -1., 1.); // a3
      model->SetParName(3, "a3");      // a3
    } else {
      model->SetParameter(0, 3.8);     // a0
      model->SetParLimits(0, 3., 4.5); // a0
      model->SetParName(0, "a0");      // a0
      model->SetParameter(1, 8.);      // a1
      model->SetParLimits(1, 5., 15.); // a1
      model->SetParName(1, "a1");      // a1
      model->SetParameter(2, 0.);      // a2
      model->SetParLimits(2, -2., 2.); // a2
      model->SetParName(2, "a2");      // a2
      model->SetParameter(3, 0.);      // a2
      model->SetParLimits(3, -1., 1.); // a3
      model->SetParName(3, "a3");      // a3
    }
  }
  if (flag == PolExp) {
    model = new TF1("PolExp", FlowAnalysis_Fitting::PolyExp,
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
                                                TH1D *hs_v2_input, TList *ls) {
  vector<double> results;
  // Fitting for dimuon invariant mass + v2 signal
  cout << ">>>>>>>>>>>>>> Start processing Pt range: "
       << FlowAnalysis_Fitting::ptmin << " " << FlowAnalysis_Fitting::ptmax
       << endl;
  cout << ">>>>>>>>>>>>>> Start processing centrality range: "
       << FlowAnalysis_Fitting::centmin << " " << FlowAnalysis_Fitting::centmax
       << endl;
  // Make copies for input histograms for safety
  TH1D *hs = dynamic_cast<TH1D *>(hs_input->Clone(Form(
      "Proj_Mass_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));
  TH1D *hs_v2 = dynamic_cast<TH1D *>(hs_v2_input->Clone(Form(
      "Proj_v2_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));

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
  model->SetParLimits(0, 1., 1.E10);
  model->SetParLimits(1, 1., 1.E10);
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
  MinimizerOptions::SetDefaultMaxIterations(10000);
  // MinimizerOptions::SetDefaultTolerance(1.);
  // MinimizerOptions::SetDefaultErrorDef(0.5);
  // MinimizerOptions::SetDefaultPrecision(1.E-9);
  // MinimizerOptions::SetDefaultPrintLevel(1);
  // IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
  // IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);

  // Iterative fitting
  // hs->Scale(1., "width");
  double chi2ndf_mass;
  int nfree_bkg = model->GetNumberFreeParameters() - 4;
  for (int i = nfree_bkg - 1; i < nPar_bkg; i++) {
    auto result = hs->Fit("model", "S Q B 0");
    int fitStatus = result;
    int bit_improve = int(fitStatus / 1000);
    int bit_minos = int((fitStatus - 1000 * bit_improve) / 100);
    int bit_hesse =
        int((fitStatus - 1000 * bit_improve - 100 * bit_minos) / 10);
    int bit_migrad = int(
        (fitStatus - 1000 * bit_improve - 100 * bit_minos - 10 * bit_hesse) /
        1);
    result->Print();
    cout << "Fit result = " << fitStatus << endl;
    chi2ndf_mass = model->GetChisquare() / model->GetNDF();
    cout << "chi2/ndf: " << chi2ndf_mass << endl;
    if ((chi2ndf_mass > mchi2max_mass || (bit_migrad != 0 && bit_migrad != 1) ||
         bit_minos != 0 || bit_hesse != 0) &&
        i < nPar_bkg - 1) {
      double min, max;
      bkg->GetParLimits(i + 1, min, max);
      model->ReleaseParameter(i + 3 + nPar_sig);
      model->SetParLimits(i + 3 + nPar_sig, min, max);
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
          "Fitted_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      Form(
          "Fitted_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()));
  c->cd();
  c->SetBottomMargin(0);
  TPad *pad1_yield = new TPad(
      Form(
          "Pad1_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0.3, 1, 1.0);
  pad1_yield->SetBottomMargin(0);
  pad1_yield->Draw();
  pad1_yield->cd();
  hs->SetMarkerStyle(20);
  hs->SetMarkerSize(0.8);
  hs->SetTitle(Form("J/#psi invariant mass: %g - %g GeV/c, %g - %g %%",
                    FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
                    FlowAnalysis_Fitting::centmin,
                    FlowAnalysis_Fitting::centmax));
  hs->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs->GetYaxis()->SetTitle(Form("Counts per %g GeV/c^{2}", hs->GetBinWidth(1)));
  hs->Draw("HIST EP");
  pad1_yield->ModifiedUpdate();
  model->SetLineWidth(3.0);
  model->SetLineColor(kBlue);
  model->Draw("same");
  pad1_yield->ModifiedUpdate();
  bkg_fitted->SetLineWidth(3.0);
  bkg_fitted->SetLineColor(kBlue);
  bkg_fitted->SetLineStyle(kDashed);
  bkg_fitted->Draw("same");
  pad1_yield->ModifiedUpdate();
  sig_fitted->SetLineWidth(3.0);
  sig_fitted->SetLineColor(kRed);
  sig_fitted->Draw("same");
  pad1_yield->ModifiedUpdate();

  double mean_fitted = model->GetParameter(2);
  double sigma_fitted = model->GetParameter(3);
  double S_3sigma = sig_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                         mean_fitted + 3. * sigma_fitted);
  double B_3sigma = bkg_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                         mean_fitted + 3. * sigma_fitted);

  TPaveStats *sb = (TPaveStats *)pad1_yield->GetPrimitive("stats");
  sb->SetName(Form(
      "Stats_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str()));
  sb->SetX1NDC(0.7);
  sb->SetX2NDC(0.95);
  sb->SetY1NDC(0.88);
  sb->SetY2NDC(0.38);
  TList *sb_list = sb->GetListOfLines();
  TText *tconst1 = sb->GetLineWith("N_{J/#psi}");
  TText *tconst2 = sb->GetLineWith("N_{bkg}");
  sb_list->Remove(tconst1);
  sb_list->Remove(tconst2);
  sb->AddText(TString::Format("%s = %d #pm %d", "N_{J#psi}",
                              int(model->GetParameter(0) / hs->GetBinWidth(1)),
                              int(model->GetParError(0) / hs->GetBinWidth(1))));
  sb->AddText(TString::Format("%s = %d #pm %d", "N_{bkg}",
                              int(model->GetParameter(1) / hs->GetBinWidth(1)),
                              int(model->GetParError(1) / hs->GetBinWidth(1))));
  sb->AddText(
      TString::Format("%s = %f", "(S/B)_{3#sigma}", S_3sigma / B_3sigma));
  sb->AddText(TString::Format("%s = %f", "(S / #sqrt{S+B})_{3#sigma}",
                              S_3sigma / pow(S_3sigma + B_3sigma, 0.5)));
  hs->SetStats(0);
  sb->Draw();
  pad1_yield->ModifiedUpdate();

  TLatex *text_info = new TLatex();
  text_info->SetTextSize(0.05);
  text_info->SetTextFont(42);
  text_info->DrawLatexNDC(.18, .82,
                          "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                          "= 5.36 TeV");
  pad1_yield->ModifiedUpdate();
  text_info->DrawLatexNDC(.18, .77,
                          "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                          "4");
  pad1_yield->ModifiedUpdate();
  text_info->DrawLatexNDC(.18, .72,
                          Form("%g < #it{p}_{T} < %g GeV/c",
                               FlowAnalysis_Fitting::ptmin,
                               FlowAnalysis_Fitting::ptmax));
  pad1_yield->ModifiedUpdate();
  c->cd();
  TPad *pad2_yield = new TPad(
      Form(
          "Pad2_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%"
          "s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0., 1, 0.3);
  pad2_yield->SetTopMargin(0);
  pad2_yield->SetBottomMargin(0.22);
  pad2_yield->Draw();
  pad2_yield->cd();

  TH1D *hs_pull_yield =
      dynamic_cast<TH1D *>(FlowAnalysis_Fitting::GetPull(hs, model, "yield"));
  hs_pull_yield->SetStats(0);
  hs_pull_yield->SetTitle("");
  hs_pull_yield->SetMarkerStyle(20);
  hs_pull_yield->SetMarkerSize(0.8);
  hs_pull_yield->GetYaxis()->SetTitleSize(0.1);
  hs_pull_yield->GetYaxis()->SetTitle("pull");
  hs_pull_yield->GetYaxis()->SetTitleOffset(0.25);
  hs_pull_yield->GetYaxis()->SetLabelSize(0.1);
  hs_pull_yield->GetYaxis()->SetRangeUser(-5., 5.);
  hs_pull_yield->GetXaxis()->SetLabelSize(0.1);
  hs_pull_yield->GetXaxis()->SetLabelOffset();
  hs_pull_yield->GetXaxis()->SetTitleSize(0.1);
  hs_pull_yield->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_pull_yield->Draw("HIST P");
  TF1 *lyield1 = new TF1("lyield1", "[0]", FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax);
  lyield1->SetParameter(0, 0.);
  lyield1->SetLineColor(kBlue);
  lyield1->SetLineWidth(3);
  lyield1->SetLineStyle(1);
  lyield1->Draw("same");
  TF1 *lyield2 = new TF1("lyield2", "[0]", FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax);
  lyield2->SetParameter(0, 1.);
  lyield2->SetLineColor(kBlue);
  lyield2->SetLineWidth(3);
  lyield2->SetLineStyle(9);
  lyield2->Draw("same");
  TF1 *lyield3 = new TF1("lyield3", "[0]", FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax);
  lyield3->SetParameter(0, -1.);
  lyield3->SetLineColor(kBlue);
  lyield3->SetLineWidth(3);
  lyield3->SetLineStyle(9);
  lyield3->Draw("same");

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
                          FlowAnalysis_Fitting::massmax,
                          FlowAnalysis_Fitting::mflag_bkg_v2 == 0 ? 4 : 5);
  model_v2->SetParameter(0, 0.01);
  model_v2->SetParLimits(0, -0.1, 0.5);
  model_v2->SetParName(0, "v^{J/#psi}_{2}");
  hs_v2->Rebin();
  hs_v2->Scale(0.5);
  double chi2ndf_v2;
  if (FlowAnalysis_Fitting::mflag_bkg_v2 == 0) {
    model_v2->SetParameter(1, 0.);
    model_v2->SetParLimits(1, -5., 5.);
    model_v2->SetParName(1, "a0");
    model_v2->SetParameter(2, 0.);
    model_v2->FixParameter(2, 0.0);
    model_v2->SetParName(2, "a2");
    model_v2->SetParameter(3, 0.);
    model_v2->FixParameter(3, 0.0);
    model_v2->SetParName(3, "a3");

    for (int i = 1; i <= 3; i++) {
      auto result_v2 = hs_v2->Fit("model_v2", "S Q B 0");
      int fitStatus_v2 = result_v2;
      int bit_improve_v2 = int(fitStatus_v2 / 1000);
      int bit_minos_v2 = int((fitStatus_v2 - 1000 * bit_improve_v2) / 100);
      int bit_hesse_v2 =
          int((fitStatus_v2 - 1000 * bit_improve_v2 - 100 * bit_minos_v2) / 10);
      int bit_migrad_v2 = int((fitStatus_v2 - 1000 * bit_improve_v2 -
                               100 * bit_minos_v2 - 10 * bit_hesse_v2) /
                              1);
      result_v2->Print();
      cout << "Fit result = " << fitStatus_v2 << endl;
      chi2ndf_v2 = model_v2->GetChisquare() / model_v2->GetNDF();
      cout << "chi2/ndf: " << chi2ndf_v2 << endl;
      if ((chi2ndf_v2 > mchi2max_v2 ||
           (bit_migrad_v2 != 0 && bit_migrad_v2 != 1) || bit_minos_v2 != 0 ||
           bit_hesse_v2 != 0) &&
          i < 3) {
        model_v2->ReleaseParameter(i + 1);
        model_v2->SetParLimits(i + 1, -2., 2.);
      } else {
        break;
      }
    }
  } else {
    model_v2->SetParameter(1, 0.0);
    model_v2->SetParLimits(1, -8., 8.);
    model_v2->SetParName(1, "a0");
    model_v2->SetParameter(2, 0.0);
    model_v2->SetParLimits(2, -2., 2.);
    model_v2->SetParName(2, "a1");
    model_v2->FixParameter(2, 0.0);
    model_v2->SetParameter(3, 0.0);
    model_v2->SetParLimits(3, -2., 2.);
    model_v2->SetParName(3, "a2");
    model_v2->FixParameter(3, 0.0);
    model_v2->SetParameter(4, 0.0);
    model_v2->SetParLimits(4, -2., 2.);
    model_v2->SetParName(4, "a3");
    model_v2->FixParameter(4, 0.0);

    for (int i = 1; i <= 4; i++) {
      auto result_v2 = hs_v2->Fit("model_v2", "S Q B 0");
      int fitStatus_v2 = result_v2;
      int bit_improve_v2 = int(fitStatus_v2 / 1000);
      int bit_minos_v2 = int((fitStatus_v2 - 1000 * bit_improve_v2) / 100);
      int bit_hesse_v2 =
          int((fitStatus_v2 - 1000 * bit_improve_v2 - 100 * bit_minos_v2) / 10);
      int bit_migrad_v2 = int((fitStatus_v2 - 1000 * bit_improve_v2 -
                               100 * bit_minos_v2 - 10 * bit_hesse_v2) /
                              1);
      result_v2->Print();
      cout << "Fit result = " << fitStatus_v2 << endl;
      chi2ndf_v2 = model_v2->GetChisquare() / model_v2->GetNDF();
      cout << "chi2/ndf: " << chi2ndf_v2 << endl;
      if ((chi2ndf_v2 > mchi2max_v2 ||
           (bit_migrad_v2 != 0 && bit_migrad_v2 != 1) || bit_minos_v2 != 0 ||
           bit_hesse_v2 != 0) &&
          i < 4) {
        model_v2->ReleaseParameter(i + 1);
        model_v2->SetParLimits(i + 1, -2., 2.);
      } else {
        break;
      }
    }
  }

  // Getting fitted background
  auto fitted_v2bkg = [](double *x, double *par) {
    double val = 0.;
    if (FlowAnalysis_Fitting::mflag_bkg_v2 == 0) {
      val = (par[0] + par[1] * x[0] + par[2] * x[0] * x[0]);
    } else {
      val = FlowAnalysis_Fitting::Cheby3(x, par);
    }
    return val;
  };
  TF1 *v2bkg_fitted =
      new TF1("v2bkg", fitted_v2bkg, FlowAnalysis_Fitting::massmin,
              FlowAnalysis_Fitting::massmax,
              FlowAnalysis_Fitting::mflag_bkg_v2 == 0 ? 3 : 4);
  if (FlowAnalysis_Fitting::mflag_bkg_v2 == 0) {
    v2bkg_fitted->SetParameter(0, model_v2->GetParameter(1));
    v2bkg_fitted->SetParameter(1, model_v2->GetParameter(2));
    v2bkg_fitted->SetParameter(2, model_v2->GetParameter(3));

  } else {
    v2bkg_fitted->SetParameter(0, model_v2->GetParameter(1));
    v2bkg_fitted->SetParameter(1, model_v2->GetParameter(2));
    v2bkg_fitted->SetParameter(2, model_v2->GetParameter(3));
    v2bkg_fitted->SetParameter(3, model_v2->GetParameter(4));
  }

  // Do plotting
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  TCanvas *c_v2 = new TCanvas(
      Form(
          "Fitted_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      Form(
          "Fitted_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()));
  c_v2->cd();
  c_v2->SetBottomMargin(0);
  TPad *pad1_v2 = new TPad(
      Form(
          "Pad1_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0.3, 1, 1.0);
  pad1_v2->SetBottomMargin(0);
  pad1_v2->Draw();
  pad1_v2->cd();
  hs_v2->SetMarkerStyle(20);
  hs_v2->SetMarkerSize(0.8);
  hs_v2->SetTitle(Form("J/#psi #it{v}_{2}{%d}: %g - %g GeV/c, %g - %g %%",
                       FlowAnalysis_Fitting::norder,
                       FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
                       FlowAnalysis_Fitting::centmin,
                       FlowAnalysis_Fitting::centmax));
  hs_v2->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_v2->GetYaxis()->SetTitle("#it{v}^{#mu#mu}_{2}");
  hs_v2->Draw("HIST EP");
  pad1_v2->ModifiedUpdate();
  model_v2->SetLineWidth(3.0);
  model_v2->SetLineColor(kBlue);
  model_v2->Draw("same");
  pad1_v2->ModifiedUpdate();
  v2bkg_fitted->SetLineWidth(3.0);
  v2bkg_fitted->SetLineColor(kBlue);
  v2bkg_fitted->SetLineStyle(kDashed);
  v2bkg_fitted->Draw("same");
  pad1_v2->ModifiedUpdate();
  TPaveStats *sb_v2 = (TPaveStats *)pad1_v2->GetPrimitive("stats");
  sb_v2->SetName("J/#psi v2 fit");
  sb_v2->SetX1NDC(0.7);
  sb_v2->SetX2NDC(0.95);
  sb_v2->SetY1NDC(0.88);
  sb_v2->SetY2NDC(0.38);
  hs_v2->SetStats(0);
  sb_v2->Draw();
  pad1_v2->ModifiedUpdate();
  c_v2->cd();
  TPad *pad2_v2 = new TPad(
      Form(
          "Pad2_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0., 1, 0.3);
  pad2_v2->SetTopMargin(0);
  pad2_v2->SetBottomMargin(0.22);
  pad2_v2->Draw();
  pad2_v2->cd();
  TH1D *hs_pull_v2 = dynamic_cast<TH1D *>(
      FlowAnalysis_Fitting::GetPull(hs_v2, model_v2, "v2"));
  hs_pull_v2->SetStats(0);
  hs_pull_v2->SetTitle("");
  hs_pull_v2->SetMarkerStyle(20);
  hs_pull_v2->SetMarkerSize(0.8);
  hs_pull_v2->GetYaxis()->SetTitleSize(0.1);
  hs_pull_v2->GetYaxis()->SetTitle("pull");
  hs_pull_v2->GetYaxis()->SetTitleOffset(0.25);
  hs_pull_v2->GetYaxis()->SetLabelSize(0.1);
  hs_pull_v2->GetYaxis()->SetRangeUser(-5., 5.);
  hs_pull_v2->GetXaxis()->SetLabelSize(0.1);
  hs_pull_v2->GetXaxis()->SetLabelOffset();
  hs_pull_v2->GetXaxis()->SetTitleSize(0.1);
  hs_pull_v2->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_pull_v2->Draw("HIST P");
  TF1 *lv21 = new TF1("lv21", "[0]", FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax);
  lv21->SetParameter(0, 0.);
  lv21->SetLineColor(kBlue);
  lv21->SetLineWidth(3);
  lv21->SetLineStyle(1);
  lv21->Draw("same");
  TF1 *lv22 = new TF1("lv22", "[0]", FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax);
  lv22->SetParameter(0, 1.);
  lv22->SetLineColor(kBlue);
  lv22->SetLineWidth(3);
  lv22->SetLineStyle(9);
  lv22->Draw("same");
  TF1 *lv23 = new TF1("lv23", "[0]", FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax);
  lv23->SetParameter(0, -1.);
  lv23->SetLineColor(kBlue);
  lv23->SetLineWidth(3);
  lv23->SetLineStyle(9);
  lv23->Draw("same");

  ls->Add(c_v2);

  cout << endl;
  cout << endl;
  results.emplace_back(model_v2->GetParameter(0));
  results.emplace_back(model_v2->GetParError(0));
  results.emplace_back(model->GetParameter(0) / hs->GetBinWidth(1));
  results.emplace_back(model->GetParError(0) / hs->GetBinWidth(1));
  results.emplace_back(S_3sigma / B_3sigma);
  results.emplace_back(chi2ndf_mass);
  results.emplace_back(chi2ndf_v2);

  return results;
}

//______________________________________________________________________________
vector<double>
FlowAnalysis_Fitting::runFittingEM(TH1D *hs_mse_input, TH1D *hs_mme_input,
                                   TH1D *hs_v2se_input, TH1D *hs_v2me_input,
                                   TH1D *hs_meanPt_input, TList *ls) {
  vector<double> results;
  // Fitting for dimuon invariant mass + v2 signal
  cout << ">>>>>>>>>>>>>> Start processing Pt range: "
       << FlowAnalysis_Fitting::ptmin << " " << FlowAnalysis_Fitting::ptmax
       << endl;
  cout << ">>>>>>>>>>>>>> Start processing centrality range: "
       << FlowAnalysis_Fitting::centmin << " " << FlowAnalysis_Fitting::centmax
       << endl;
  // Make copies for input histograms for safety
  TH1D *hs_se = dynamic_cast<TH1D *>(hs_mse_input->Clone(Form(
      "Proj_MassSE_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));
  TH1D *hs_me = dynamic_cast<TH1D *>(hs_mme_input->Clone(Form(
      "Proj_MassME_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));
  TH1D *hs_v2se = dynamic_cast<TH1D *>(hs_v2se_input->Clone(Form(
      "Proj_v2SE_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));
  TH1D *hs_v2me = dynamic_cast<TH1D *>(hs_v2me_input->Clone(Form(
      "Proj_v2ME_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));
  TH1D *hs_meanPt = new TH1D();
  if (hs_meanPt_input) {
    hs_meanPt = dynamic_cast<TH1D *>(hs_meanPt_input->Clone(Form(
        "Proj_meanPt_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
        FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
        FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
        FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
        FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
        FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
        FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
            .c_str(),
        FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
            .c_str(),
        FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
            .c_str())));
  }
  TH1D *hs_se_res = dynamic_cast<TH1D *>(hs_mse_input->Clone(Form(
      "Proj_MassSEResidual_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str())));
  hs_se_res->Add(hs_me, -1.0);

  //////////////////////////////////////////////////////////////////////////////
  ///      INVARIANT MASS FIT
  //////////////////////////////////////////////////////////////////////////////
  // Setting up fit model for invariant mass fit
  cout << ">>>>>>>>>>>>>> Start processing dimuon invariant mass fit..."
       << endl;
  // Construct a combined signal+background model
  TF1 *sig = new TF1();
  TF1 *bkg = new TF1();
  CreateModel(sig, ModelType(FlowAnalysis_Fitting::mflag_sig));
  auto fct_bkgResidual = [](double *x, double *par) {
    double value = 0.0;
    value = exp(par[0] * x[0]);
    return value;
  };
  if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
      "EventMixing") {
    bkg = new TF1("bkgResidual", fct_bkgResidual, FlowAnalysis_Fitting::massmin,
                  FlowAnalysis_Fitting::massmax, 1);
    int idx_massmin = hs_me->FindBin(FlowAnalysis_Fitting::massmin);
    double initMin = hs_me->GetBinContent(idx_massmin);
    bkg->SetParameter(0, 0.);
    bkg->SetParLimits(0, -2., 0.);
    bkg->SetParName(0, "a0");
  } else {
    CreateModel(bkg, ModelType(FlowAnalysis_Fitting::mflag_bkg));
  }
  int nsig = 5.E2;
  int nbkg = 5.E5;
  TF1NormSum *sum_model = new TF1NormSum(sig, bkg, nsig, nbkg);
  TF1 *model = new TF1("model", *sum_model, FlowAnalysis_Fitting::massmin,
                       FlowAnalysis_Fitting::massmax, sum_model->GetNpar());
  model->SetParameters(sum_model->GetParameters().data());
  model->SetParName(0, "N_{J/#psi}");
  model->SetParName(1, "N_{bkg}");
  model->SetParLimits(0, 1., 1.E10);
  model->SetParLimits(1, 1., 1.E10);
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
  MinimizerOptions::SetDefaultMaxIterations(10000);
  // MinimizerOptions::SetDefaultTolerance(1.);
  // MinimizerOptions::SetDefaultErrorDef(0.5);
  // MinimizerOptions::SetDefaultPrecision(1.E-9);
  // MinimizerOptions::SetDefaultPrintLevel(1);
  // IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
  // IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);

  // Iterative fitting
  double chi2ndf_mass;
  int nfree_bkg = model->GetNumberFreeParameters() - 4;
  for (int i = nfree_bkg - 1; i < nPar_bkg; i++) {
    auto result =
        FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
                "EventMixing"
            ? hs_se_res->Fit("model", "S Q B 0")
            : hs_se->Fit("model", "S Q B 0");
    int fitStatus = result;
    int bit_improve = int(fitStatus / 1000);
    int bit_minos = int((fitStatus - 1000 * bit_improve) / 100);
    int bit_hesse =
        int((fitStatus - 1000 * bit_improve - 100 * bit_minos) / 10);
    int bit_migrad = int(
        (fitStatus - 1000 * bit_improve - 100 * bit_minos - 10 * bit_hesse) /
        1);
    result->Print();
    cout << "Fit result = " << fitStatus << endl;
    chi2ndf_mass = model->GetChisquare() / model->GetNDF();
    cout << "chi2/ndf: " << chi2ndf_mass << endl;
    if ((chi2ndf_mass > mchi2max_mass || (bit_migrad != 0 && bit_migrad != 1) ||
         bit_minos != 0 || bit_hesse != 0) &&
        i < nPar_bkg - 1) {
      double min, max;
      bkg->GetParLimits(i + 1, min, max);
      model->ReleaseParameter(i + 3 + nPar_sig);
      model->SetParLimits(i + 3 + nPar_sig, min, max);
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
  auto fct_bkgResidual_norm = [](double *x, double *par) {
    double value = 0.0;
    value = par[0] * exp(par[2] * x[0]) / par[1];
    return value;
  };
  TF1 *bkg_fitted = new TF1();
  if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
      "EventMixing") {
    bkg_fitted = new TF1("bkgResidual_norm", fct_bkgResidual_norm,
                         FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax, 3);
  } else {
    bkg_fitted = new TF1("bkg_fitted", FlowAnalysis_Fitting::FittedBkg,
                         FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax, nPar_bkg + 2);
  }
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
          "Fitted_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      Form(
          "Fitted_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()));
  c->cd();
  c->SetBottomMargin(0);
  TPad *pad1_yield = new TPad(
      Form(
          "Pad1_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0.3, 1, 1.0);
  pad1_yield->SetBottomMargin(0);
  pad1_yield->Draw();
  pad1_yield->cd();
  if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
      "EventMixing") {
    hs_se_res->SetStats(1);
    hs_se_res->SetTitle("");
    hs_se_res->SetMarkerStyle(4);
    hs_se_res->SetMarkerSize(0.6);
    hs_se_res->SetMarkerColor(kBlue);
    hs_se_res->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    hs_se_res->GetYaxis()->SetTitle(
        Form("Counts per %g GeV/c^{2}", hs_se->GetBinWidth(1)));
    hs_se_res->Draw("HIST EP");
    hs_se_res->GetYaxis()->SetRangeUser(
        0.1, 1.2 * hs_se->GetBinContent(hs_se->GetMaximumBin()));
    pad1_yield->ModifiedUpdate();
    hs_se->SetStats(0);
    hs_se->SetMarkerStyle(20);
    hs_se->SetMarkerSize(0.6);
    hs_se->SetMarkerColor(kBlack);
    hs_se->SetTitle("");
    hs_se->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    hs_se->GetYaxis()->SetTitle(
        Form("Counts per %g GeV/c^{2}", hs_se->GetBinWidth(1)));
    hs_se->Draw("HIST EP same");
    pad1_yield->ModifiedUpdate();
  } else {
    hs_se->SetStats(1);
    hs_se->SetMarkerStyle(20);
    hs_se->SetMarkerSize(0.6);
    hs_se->SetMarkerColor(kBlack);
    hs_se->SetTitle("");
    hs_se->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    hs_se->GetYaxis()->SetTitle(
        Form("Counts per %g GeV/c^{2}", hs_se->GetBinWidth(1)));
    hs_se->GetYaxis()->SetRangeUser(
        0.1, 1.2 * hs_se->GetBinContent(hs_se->GetMaximumBin()));
    hs_se->Draw("HIST EP");
    pad1_yield->ModifiedUpdate();
  }
  hs_me->SetStats(0);
  hs_me->SetTitle("");
  hs_me->SetMarkerStyle(4);
  hs_me->SetMarkerSize(0.6);
  hs_me->SetMarkerColor(kBlack);
  hs_me->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_me->GetYaxis()->SetTitle(
      Form("Counts per %g GeV/c^{2}", hs_se->GetBinWidth(1)));
  hs_me->Draw("HIST EP same");
  pad1_yield->ModifiedUpdate();
  model->SetLineWidth(2.0);
  model->SetLineColor(kBlue);
  model->SetTitle("");
  model->Draw("same");
  pad1_yield->ModifiedUpdate();
  bkg_fitted->SetLineWidth(2.0);
  bkg_fitted->SetLineColor(kBlue);
  bkg_fitted->SetLineStyle(kDashed);
  bkg_fitted->SetTitle("");
  bkg_fitted->Draw("same");
  pad1_yield->ModifiedUpdate();
  sig_fitted->SetLineWidth(2.0);
  sig_fitted->SetLineColor(kRed);
  sig_fitted->SetTitle("");
  sig_fitted->Draw("same");
  pad1_yield->ModifiedUpdate();

  double mean_fitted = model->GetParameter(2);
  double sigma_fitted = model->GetParameter(3);
  double S_3sigma = sig_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                         mean_fitted + 3. * sigma_fitted);
  double B_3sigma =
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
              "EventMixing"
          ? bkg_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                 mean_fitted + 3. * sigma_fitted) +
                hs_me->GetBinWidth(1) *
                    hs_me->Integral(
                        hs_me->FindBin(mean_fitted - 3. * sigma_fitted),
                        hs_me->FindBin(mean_fitted + 3. * sigma_fitted))
          : bkg_fitted->Integral(mean_fitted - 3. * sigma_fitted,
                                 mean_fitted + 3. * sigma_fitted);
  TPaveStats *sb = (TPaveStats *)pad1_yield->GetPrimitive("stats");
  sb->SetName(Form(
      "Stats_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
      FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
      FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
      FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
      FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
      FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
          .c_str(),
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
          .c_str(),
      FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
          .c_str()));
  sb->SetX1NDC(0.64);
  sb->SetX2NDC(0.89);
  sb->SetY1NDC(0.88);
  sb->SetY2NDC(0.38);
  TList *sb_list = sb->GetListOfLines();
  TText *tconst1 = sb->GetLineWith("N_{J/#psi}");
  TText *tconst2 = sb->GetLineWith("N_{bkg}");
  sb_list->Remove(tconst1);
  sb_list->Remove(tconst2);
  sb->AddText(
      TString::Format("%s = %d #pm %d", "N_{J/#psi}",
                      int(model->GetParameter(0) / hs_se->GetBinWidth(1)),
                      int(model->GetParError(0) / hs_se->GetBinWidth(1))));
  sb->AddText(TString::Format(
      "%s = %d #pm %d", "N_{bkg}",
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
              "EventMixing"
          ? int(model->GetParameter(1) / hs_se->GetBinWidth(1) +
                hs_me->Integral())
          : int(model->GetParameter(1) / hs_se->GetBinWidth(1)),
      int(model->GetParError(1) / hs_se->GetBinWidth(1))));
  sb->AddText(
      TString::Format("%s = %f", "(S/B)_{3#sigma}", S_3sigma / B_3sigma));
  sb->AddText(TString::Format("%s = %f", "(S / #sqrt{S+B})_{3#sigma}",
                              S_3sigma / pow(S_3sigma + B_3sigma, 0.5)));
  if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
      "EventMixing") {
    hs_se_res->SetStats(0);
  } else {
    hs_se->SetStats(0);
  }
  sb->Draw();
  pad1_yield->ModifiedUpdate();

  TLatex *text_info_yield = new TLatex();
  /*
  text_info->SetTextSize(0.04);
  text_info->SetTextFont(42);
  text_info->DrawLatexNDC(.12, .85,
                          "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
                          "= 5.36 TeV");
  pad1_yield->ModifiedUpdate();
  text_info->DrawLatexNDC(.12, .80,
                          "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < "
                          "4");
  pad1_yield->ModifiedUpdate();
  */
  text_info_yield->SetTextSize(0.05);
  text_info_yield->SetTextFont(62);
  text_info_yield->DrawLatexNDC(
      .12, .82,
      Form("%g < #it{p}_{T} < %g GeV/c    %g-%g%%", FlowAnalysis_Fitting::ptmin,
           FlowAnalysis_Fitting::ptmax, FlowAnalysis_Fitting::centmin,
           FlowAnalysis_Fitting::centmax));
  pad1_yield->ModifiedUpdate();
  TLegend *legend_yield = new TLegend(0.13, 0.05, 0.35, 0.45);
  legend_yield->SetBorderSize(0);
  legend_yield->SetFillStyle(0);
  legend_yield->AddEntry(hs_se, "Data", "P");
  legend_yield->AddEntry(hs_me, "Event mixing", "P");
  if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
      "EventMixing") {
    legend_yield->AddEntry(hs_se_res, "Residual", "P");
  }
  legend_yield->AddEntry(model, "Total");
  legend_yield->AddEntry(sig_fitted, "Signal");
  legend_yield->AddEntry(bkg_fitted, "Background");
  legend_yield->Draw("same");
  pad1_yield->ModifiedUpdate();

  c->cd();
  TPad *pad2_yield = new TPad(
      Form(
          "Pad2_%s_Mass_v%d%d_%g_%g_%g_%g_%g_%g_%"
          "s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0., 1, 0.3);
  pad2_yield->SetTopMargin(0);
  pad2_yield->SetBottomMargin(0.22);
  pad2_yield->Draw();
  pad2_yield->cd();
  TH1D *hs_pull_yield =
      FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
              "EventMixing"
          ? dynamic_cast<TH1D *>(
                FlowAnalysis_Fitting::GetPull(hs_se_res, model, "yield"))
          : dynamic_cast<TH1D *>(
                FlowAnalysis_Fitting::GetPull(hs_se, model, "yield"));
  hs_pull_yield->SetStats(0);
  hs_pull_yield->SetTitle("");
  hs_pull_yield->SetMarkerStyle(20);
  hs_pull_yield->SetMarkerSize(0.6);
  hs_pull_yield->GetYaxis()->SetTitleSize(0.1);
  hs_pull_yield->GetYaxis()->SetTitle("pull");
  hs_pull_yield->GetYaxis()->CenterTitle();
  hs_pull_yield->GetYaxis()->SetTitleOffset(0.25);
  hs_pull_yield->GetYaxis()->SetLabelSize(0.1);
  hs_pull_yield->GetYaxis()->SetRangeUser(-4.9, 4.9);
  hs_pull_yield->GetXaxis()->SetLabelSize(0.1);
  hs_pull_yield->GetXaxis()->SetLabelOffset();
  hs_pull_yield->GetXaxis()->SetTitleSize(0.1);
  hs_pull_yield->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_pull_yield->Draw("HIST P");
  TF1 *lyield1 = new TF1("lyield1", "[0]", FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax);
  lyield1->SetParameter(0, 0.);
  lyield1->SetLineColor(kBlue);
  lyield1->SetLineWidth(2);
  lyield1->SetLineStyle(1);
  lyield1->Draw("same");
  TF1 *lyield2 = new TF1("lyield2", "[0]", FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax);
  lyield2->SetParameter(0, 1.);
  lyield2->SetLineColor(kBlue);
  lyield2->SetLineWidth(2);
  lyield2->SetLineStyle(9);
  lyield2->Draw("same");
  TF1 *lyield3 = new TF1("lyield3", "[0]", FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax);
  lyield3->SetParameter(0, -1.);
  lyield3->SetLineColor(kBlue);
  lyield3->SetLineWidth(2);
  lyield3->SetLineStyle(9);
  lyield3->Draw("same");
  pad2_yield->ModifiedUpdate();

  ls->Add(c);
  // Constructing alpha function
  auto fct_alpha = [sig_fitted, bkg_fitted, hs_me](double *x, double *) {
    double S = sig_fitted->Eval(x[0]);
    double B = 0.;
    if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
        "EventMixing") {
      double B_res = bkg_fitted->Eval(x[0]);
      double B_em = hs_me->GetBinContent(hs_me->FindBin(x[0]));
      B = B_res + B_em;
    } else {
      B = bkg_fitted->Eval(x[0]);
    }
    return S / (S + B);
  };
  TF1 *alpha = new TF1("alpha", fct_alpha, FlowAnalysis_Fitting::massmin,
                       FlowAnalysis_Fitting::massmax, 0);

  //////////////////////////////////////////////////////////////////////////////
  ///      Mean pT FIT
  //////////////////////////////////////////////////////////////////////////////
  double meanPt_fitted = 0.;
  double meanPterror_fitted = 0.;
  double chi2ndf_meanPt = 0.;
  TF1 *model_meanPt = new TF1();
  TGraph *hs_model_meanPt = new TGraph();
  TLegend *legend_meanPt = new TLegend(0.13, 0.45, 0.28, 0.65);
  if (hs_meanPt_input) {
    // Setting up fit model for mean pT fit
    cout << ">>>>>>>>>>>>>> Start processing mean pT fit..." << endl;
    auto fct_meanPt = [alpha](double *x, double *par) {
      double value = 0.0;
      double val_alpha = alpha->Eval(x[0]);
      value = (par[0] * val_alpha +
               (1 - val_alpha) * FlowAnalysis_Fitting::Cheby6(x, &par[1]));

      return value;
    };

    model_meanPt =
        FlowAnalysis_Fitting::ptmax <= 3
            ? new TF1("model_meanPt", fct_meanPt, 2.5, 3.8, 8)
            : new TF1("model_meanPt", fct_meanPt, FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax, 8);

    model_meanPt->SetParameter(
        0, (FlowAnalysis_Fitting::ptmin + FlowAnalysis_Fitting::ptmax) / 2.);
    model_meanPt->SetParLimits(0, FlowAnalysis_Fitting::ptmin,
                               FlowAnalysis_Fitting::ptmax);
    model_meanPt->SetParName(0, "<#it{p}_{T}>^{J/#psi}");
    model_meanPt->SetParameter(
        1, (FlowAnalysis_Fitting::ptmin + FlowAnalysis_Fitting::ptmax) / 2.);
    model_meanPt->SetParLimits(1, FlowAnalysis_Fitting::ptmin,
                               FlowAnalysis_Fitting::ptmax);
    model_meanPt->SetParName(1, "a0");
    model_meanPt->SetParameter(2, 0.0);
    model_meanPt->SetParLimits(2, -2., 2.);
    model_meanPt->SetParName(2, "a1");
    model_meanPt->FixParameter(2, 0.0);
    model_meanPt->SetParameter(3, 0.0);
    model_meanPt->SetParLimits(3, -2., 2.);
    model_meanPt->SetParName(3, "a2");
    model_meanPt->FixParameter(3, 0.0);
    model_meanPt->SetParameter(4, 0.0);
    model_meanPt->SetParLimits(4, -2., 2.);
    model_meanPt->SetParName(4, "a3");
    model_meanPt->FixParameter(4, 0.0);
    model_meanPt->SetParameter(5, 0.0);
    model_meanPt->SetParLimits(5, -1., 1.);
    model_meanPt->SetParName(5, "a4");
    model_meanPt->FixParameter(5, 0.0);
    model_meanPt->SetParameter(6, 0.0);
    model_meanPt->SetParLimits(6, -1., 1.);
    model_meanPt->SetParName(6, "a5");
    model_meanPt->FixParameter(6, 0.0);
    model_meanPt->SetParameter(7, 0.0);
    model_meanPt->SetParLimits(7, -1., 1.);
    model_meanPt->SetParName(7, "a6");
    model_meanPt->FixParameter(7, 0.0);

    hs_meanPt->Rebin();
    hs_meanPt->Scale(0.5);

    for (int i = 1; i <= 7; i++) {
      auto result_meanPt = FlowAnalysis_Fitting::ptmax <= 3
                               ? hs_meanPt->Fit("model_meanPt", "S Q B R 0")
                               : hs_meanPt->Fit("model_meanPt", "S Q B 0");
      int fitStatus_meanPt = result_meanPt;
      int bit_improve_meanPt = int(fitStatus_meanPt / 1000);
      int bit_minos_meanPt =
          int((fitStatus_meanPt - 1000 * bit_improve_meanPt) / 100);
      int bit_hesse_meanPt = int((fitStatus_meanPt - 1000 * bit_improve_meanPt -
                                  100 * bit_minos_meanPt) /
                                 10);
      int bit_migrad_meanPt =
          int((fitStatus_meanPt - 1000 * bit_improve_meanPt -
               100 * bit_minos_meanPt - 10 * bit_hesse_meanPt) /
              1);
      result_meanPt->Print();
      cout << "Fit result = " << fitStatus_meanPt << endl;
      chi2ndf_meanPt = model_meanPt->GetChisquare() / model_meanPt->GetNDF();
      cout << "chi2/ndf: " << chi2ndf_meanPt << endl;
      if ((chi2ndf_meanPt > mchi2max_mass ||
           (bit_migrad_meanPt != 0 && bit_migrad_meanPt != 1) ||
           bit_minos_meanPt != 0 || bit_hesse_meanPt != 0) &&
          i < 7) {
        model_meanPt->ReleaseParameter(i + 1);
        model_meanPt->SetParLimits(i + 1, -1., 1.);
      } else {
        break;
      }
    }

    meanPt_fitted = model_meanPt->GetParameter(0);
    meanPterror_fitted = model_meanPt->GetParError(0);

    // Do plotting
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();
    TCanvas *c_meanPt = new TCanvas(
        Form("Fitted_%s_MeanPt_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
             FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode]
                 .c_str(),
             FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
                 .c_str(),
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
                 .c_str(),
             FlowAnalysis_Fitting::v2bkg_string
                 [FlowAnalysis_Fitting::mflag_bkg_v2]
                     .c_str()),
        Form("Fitted_%s_MeanPt_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
             FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode]
                 .c_str(),
             FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
                 .c_str(),
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
                 .c_str(),
             FlowAnalysis_Fitting::v2bkg_string
                 [FlowAnalysis_Fitting::mflag_bkg_v2]
                     .c_str()));
    c_meanPt->cd();
    c_meanPt->SetBottomMargin(0);
    TPad *pad1_meanPt = new TPad(
        Form("Pad1_%s_MeanPt_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
             FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode]
                 .c_str(),
             FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
                 .c_str(),
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
                 .c_str(),
             FlowAnalysis_Fitting::v2bkg_string
                 [FlowAnalysis_Fitting::mflag_bkg_v2]
                     .c_str()),
        "", 0, 0.3, 1, 1.0);
    pad1_meanPt->SetBottomMargin(0);
    pad1_meanPt->Draw();
    pad1_meanPt->cd();
    hs_meanPt->SetStats(1);
    hs_meanPt->SetTitle("");
    hs_meanPt->SetMarkerStyle(20);
    hs_meanPt->SetMarkerSize(0.8);
    hs_meanPt->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    hs_meanPt->GetYaxis()->SetTitle("<#it{p}^{#mu#mu}_{T}> (GeV/c)");
    hs_meanPt->Draw("HIST EP");
    pad1_meanPt->ModifiedUpdate();
    /*
    hs_model_meanPt = dynamic_cast<TH1D *>(FlowAnalysis_Fitting::GetHistFromTF(
        hs_meanPt, model_meanPt, "ModelMeanPt"));
    */
    hs_model_meanPt = new TGraph(model_meanPt);
    hs_model_meanPt->SetLineWidth(3.0);
    hs_model_meanPt->SetLineColor(kBlue);
    hs_model_meanPt->SetTitle("");
    hs_model_meanPt->Draw("C same");
    pad1_meanPt->ModifiedUpdate();
    TPaveStats *sb_meanPt = (TPaveStats *)pad1_meanPt->GetPrimitive("stats");
    sb_meanPt->SetName("J/#psi <#it{p}_{T}> fit");
    sb_meanPt->SetX1NDC(0.64);
    sb_meanPt->SetX2NDC(0.89);
    sb_meanPt->SetY1NDC(0.88);
    sb_meanPt->SetY2NDC(0.38);
    hs_meanPt->SetStats(0);
    sb_meanPt->Draw();
    pad1_meanPt->ModifiedUpdate();
    TLatex *text_info_meanPt = new TLatex();
    text_info_meanPt->SetTextSize(0.05);
    text_info_meanPt->SetTextFont(62);
    text_info_meanPt->DrawLatexNDC(
        .12, .82,
        Form("%g < #it{p}_{T} < %g GeV/c    %g-%g%%",
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax));
    pad1_meanPt->ModifiedUpdate();
    legend_meanPt->SetBorderSize(0);
    legend_meanPt->SetFillStyle(0);
    legend_meanPt->AddEntry(hs_meanPt, "Data", "P");
    legend_meanPt->AddEntry(hs_model_meanPt, "Fit");
    legend_meanPt->Draw("same");
    pad1_meanPt->ModifiedUpdate();

    c_meanPt->cd();
    TPad *pad2_meanPt = new TPad(
        Form("Pad2_%s_MeanPt_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
             FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode]
                 .c_str(),
             FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
                 .c_str(),
             FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
                 .c_str(),
             FlowAnalysis_Fitting::v2bkg_string
                 [FlowAnalysis_Fitting::mflag_bkg_v2]
                     .c_str()),
        "", 0, 0., 1, 0.3);
    pad2_meanPt->SetTopMargin(0);
    pad2_meanPt->SetBottomMargin(0.22);
    pad2_meanPt->Draw();
    pad2_meanPt->cd();
    TH1D *hs_pull_meanPt = dynamic_cast<TH1D *>(
        FlowAnalysis_Fitting::GetPull(hs_meanPt, model_meanPt, "MeanPt"));
    hs_pull_meanPt->SetStats(0);
    hs_pull_meanPt->SetTitle("");
    hs_pull_meanPt->SetMarkerStyle(20);
    hs_pull_meanPt->SetMarkerSize(0.6);
    hs_pull_meanPt->GetYaxis()->SetTitleSize(0.1);
    hs_pull_meanPt->GetYaxis()->SetTitle("pull");
    hs_pull_meanPt->GetYaxis()->CenterTitle();
    hs_pull_meanPt->GetYaxis()->SetTitleOffset(0.25);
    hs_pull_meanPt->GetYaxis()->SetLabelSize(0.1);
    hs_pull_meanPt->GetYaxis()->SetRangeUser(-4.9, 4.9);
    hs_pull_meanPt->GetXaxis()->SetLabelSize(0.1);
    hs_pull_meanPt->GetXaxis()->SetLabelOffset();
    hs_pull_meanPt->GetXaxis()->SetTitleSize(0.1);
    hs_pull_meanPt->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
    hs_pull_meanPt->Draw("HIST P");
    TF1 *lmeanPt1 = new TF1("lmeanPt1", "[0]", FlowAnalysis_Fitting::massmin,
                            FlowAnalysis_Fitting::massmax);
    lmeanPt1->SetParameter(0, 0.);
    lmeanPt1->SetLineColor(kBlue);
    lmeanPt1->SetLineWidth(2);
    lmeanPt1->SetLineStyle(1);
    lmeanPt1->Draw("same");
    TF1 *lmeanPt2 = new TF1("lmeanPt2", "[0]", FlowAnalysis_Fitting::massmin,
                            FlowAnalysis_Fitting::massmax);
    lmeanPt2->SetParameter(0, 1.);
    lmeanPt2->SetLineColor(kBlue);
    lmeanPt2->SetLineWidth(2);
    lmeanPt2->SetLineStyle(9);
    lmeanPt2->Draw("same");
    TF1 *lmeanPt3 = new TF1("lmeanPt3", "[0]", FlowAnalysis_Fitting::massmin,
                            FlowAnalysis_Fitting::massmax);
    lmeanPt3->SetParameter(0, -1.);
    lmeanPt3->SetLineColor(kBlue);
    lmeanPt3->SetLineWidth(2);
    lmeanPt3->SetLineStyle(9);
    lmeanPt3->Draw("same");
    pad2_meanPt->ModifiedUpdate();
    ls->Add(c_meanPt);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///      V2 FIT
  //////////////////////////////////////////////////////////////////////////////
  // Setting up fit model for v2 signal fit
  cout << ">>>>>>>>>>>>>> Start processing v2 fit..." << endl;
  TH1D *hs_v2bkg = new TH1D();
  if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
      "EventMixing") {
    auto fct_bkgSum = [bkg_fitted, hs_me](double *x, double *) {
      int idx = hs_me->FindBin(x[0]);
      double val_me = hs_me->GetBinContent(idx);
      double val_fitted = bkg_fitted->Eval(x[0]);
      return val_me + val_fitted;
    };
    TF1 *bkg_sum =
        new TF1("fct_bkgSum", fct_bkgSum, FlowAnalysis_Fitting::massmin,
                FlowAnalysis_Fitting::massmax, 0);
    hs_v2bkg = FlowAnalysis_Fitting::GetV2BkgCorrected(hs_v2me, hs_me, bkg_sum);
  } else {
    hs_v2bkg =
        FlowAnalysis_Fitting::GetV2BkgCorrected(hs_v2me, hs_me, bkg_fitted);
  }

  hs_v2bkg->Rebin();
  hs_v2bkg->Scale(0.5);
  hs_v2me->Rebin();
  hs_v2me->Scale(0.5);

  auto fct_v2 = [alpha, hs_v2bkg, hs_v2me](double *x, double *par) {
    double val_fct = 0.0;
    double val_alpha = alpha->Eval(x[0]);
    int idx = hs_v2bkg->FindBin(x[0]);
    double val_bkg = hs_v2bkg->GetBinContent(idx);
    double val_v2me = hs_v2me->GetBinContent(idx);

    if (FlowAnalysis_Fitting::mode == 0) {
      val_fct = (par[0] * val_alpha + (1. - val_alpha) * val_bkg);
    } else {
      if (FlowAnalysis_Fitting::v2bkg_string
              [FlowAnalysis_Fitting::mflag_bkg_v2] ==
          "EventMixing(beta free)") {
        val_fct =
            (par[0] * val_alpha +
             (1. - val_alpha) * (val_bkg + par[1] * (val_v2me - val_bkg)));
      } else {
        val_fct = (par[0] * val_alpha + (1. - val_alpha) * val_bkg);
      }
    }

    return val_fct;
  };

  // Setup for model
  TF1 *model_v2 = new TF1();
  if (FlowAnalysis_Fitting::mode == 0) {
    model_v2 = new TF1("model_v2", fct_v2, FlowAnalysis_Fitting::massmin,
                       FlowAnalysis_Fitting::massmax, 1);
    model_v2->SetParameter(0, 0.01);
    model_v2->SetParLimits(0, -0.1, 0.5);
    model_v2->SetParName(0, "v^{J/#psi}_{2}");
  } else {
    if (FlowAnalysis_Fitting::v2bkg_string
            [FlowAnalysis_Fitting::mflag_bkg_v2] == "EventMixing(beta free)") {
      model_v2 = new TF1("model_v2", fct_v2, FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax, 2);
      model_v2->SetParameter(0, 0.01);
      model_v2->SetParLimits(0, -0.1, 0.5);
      model_v2->SetParName(0, "v^{J/#psi}_{2}");
      model_v2->SetParameter(1, 0.01);
      model_v2->SetParLimits(1, 0., 2.);
      model_v2->SetParName(1, "#beta");
    } else {
      model_v2 = new TF1("model_v2", fct_v2, FlowAnalysis_Fitting::massmin,
                         FlowAnalysis_Fitting::massmax, 1);
      model_v2->SetParameter(0, 0.01);
      model_v2->SetParLimits(0, -0.1, 0.5);
      model_v2->SetParName(0, "v^{J/#psi}_{2}");
    }
  }
  hs_v2se->Rebin();
  hs_v2se->Scale(0.5);

  // Do fitting
  double chi2ndf_v2;
  auto result_v2 = hs_v2se->Fit("model_v2", "S Q B 0");
  int fitStatus_v2 = result_v2;
  int bit_improve_v2 = int(fitStatus_v2 / 1000);
  int bit_minos_v2 = int((fitStatus_v2 - 1000 * bit_improve_v2) / 100);
  int bit_hesse_v2 =
      int((fitStatus_v2 - 1000 * bit_improve_v2 - 100 * bit_minos_v2) / 10);
  int bit_migrad_v2 = int((fitStatus_v2 - 1000 * bit_improve_v2 -
                           100 * bit_minos_v2 - 10 * bit_hesse_v2) /
                          1);
  result_v2->Print();
  cout << "Fit result = " << fitStatus_v2 << endl;
  chi2ndf_v2 = model_v2->GetChisquare() / model_v2->GetNDF();
  cout << "chi2/ndf: " << chi2ndf_v2 << endl;

  // Do plotting
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  TCanvas *c_v2 = new TCanvas(
      Form(
          "Fitted_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      Form(
          "Fitted_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()));
  c_v2->cd();
  c_v2->SetBottomMargin(0);
  TPad *pad1_v2 = new TPad(
      Form(
          "Pad1_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0.3, 1, 1.0);
  pad1_v2->SetBottomMargin(0);
  pad1_v2->Draw();
  pad1_v2->cd();
  hs_v2se->SetStats(1);
  hs_v2se->SetMarkerStyle(20);
  hs_v2se->SetMarkerSize(0.8);
  hs_v2se->SetTitle("");
  hs_v2se->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_v2se->GetYaxis()->SetTitle("#it{v}^{#mu#mu}_{2}{SP}");
  hs_v2se->Draw("HIST EP");
  pad1_v2->ModifiedUpdate();
  TF1 *lv2 = new TF1("lv2", "[0]", FlowAnalysis_Fitting::massmin,
                     FlowAnalysis_Fitting::massmax);
  lv2->SetParameter(0, 0.);
  lv2->SetLineColor(18);
  lv2->SetLineWidth(3);
  lv2->SetLineStyle(9);
  lv2->Draw("same");
  pad1_v2->ModifiedUpdate();
  hs_v2bkg->SetStats(0);
  hs_v2bkg->SetTitle("");
  hs_v2bkg->SetMarkerStyle(4);
  hs_v2bkg->SetMarkerSize(0.8);
  hs_v2bkg->SetMarkerColor(kBlue);
  hs_v2bkg->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_v2bkg->GetYaxis()->SetTitle("#it{v}^{#mu#mu}_{2}");
  hs_v2bkg->Draw("HIST EP same");
  pad1_v2->ModifiedUpdate();
  TH1D *hs_model_v2 = dynamic_cast<TH1D *>(
      FlowAnalysis_Fitting::GetHistFromTF(hs_v2se, model_v2, "ModelV2"));
  hs_model_v2->SetLineWidth(3.0);
  hs_model_v2->SetLineColor(kBlue);
  hs_model_v2->SetTitle("");
  hs_model_v2->Draw("HIST same");
  pad1_v2->ModifiedUpdate();
  TPaveStats *sb_v2 = (TPaveStats *)pad1_v2->GetPrimitive("stats");
  sb_v2->SetName("J/#psi v2 fit");
  sb_v2->SetX1NDC(0.64);
  sb_v2->SetX2NDC(0.89);
  sb_v2->SetY1NDC(0.88);
  sb_v2->SetY2NDC(0.38);
  hs_v2se->SetStats(0);
  sb_v2->Draw();
  pad1_v2->ModifiedUpdate();
  TLatex *text_info_v2 = new TLatex();
  text_info_v2->SetTextSize(0.05);
  text_info_v2->SetTextFont(62);
  text_info_v2->DrawLatexNDC(
      .12, .82,
      Form("%g < #it{p}_{T} < %g GeV/c    %g-%g%%", FlowAnalysis_Fitting::ptmin,
           FlowAnalysis_Fitting::ptmax, FlowAnalysis_Fitting::centmin,
           FlowAnalysis_Fitting::centmax));
  pad1_v2->ModifiedUpdate();
  TLegend *legend_v2 = new TLegend(0.13, 0.45, 0.35, 0.7);
  legend_v2->SetBorderSize(0);
  legend_v2->SetFillStyle(0);
  legend_v2->AddEntry(hs_v2se, "Data", "P");
  legend_v2->AddEntry(hs_v2bkg, "Event mixing", "P");
  legend_v2->AddEntry(hs_model_v2, "Fit");
  legend_v2->Draw("same");
  pad1_v2->ModifiedUpdate();

  c_v2->cd();
  TPad *pad2_v2 = new TPad(
      Form(
          "Pad2_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", 0, 0., 1, 0.3);
  pad2_v2->SetTopMargin(0);
  pad2_v2->SetBottomMargin(0.22);
  pad2_v2->Draw();
  pad2_v2->cd();
  TH1D *hs_pull_v2 = dynamic_cast<TH1D *>(
      FlowAnalysis_Fitting::GetPull(hs_v2se, model_v2, "v2"));
  hs_pull_v2->SetStats(0);
  hs_pull_v2->SetTitle("");
  hs_pull_v2->SetMarkerStyle(20);
  hs_pull_v2->SetMarkerSize(0.6);
  hs_pull_v2->GetYaxis()->SetTitleSize(0.1);
  hs_pull_v2->GetYaxis()->SetTitle("pull");
  hs_pull_v2->GetYaxis()->CenterTitle();
  hs_pull_v2->GetYaxis()->SetTitleOffset(0.25);
  hs_pull_v2->GetYaxis()->SetLabelSize(0.1);
  hs_pull_v2->GetYaxis()->SetRangeUser(-4.9, 4.9);
  hs_pull_v2->GetXaxis()->SetLabelSize(0.1);
  hs_pull_v2->GetXaxis()->SetLabelOffset();
  hs_pull_v2->GetXaxis()->SetTitleSize(0.1);
  hs_pull_v2->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  hs_pull_v2->Draw("HIST P");
  TF1 *lv21 = new TF1("lv21", "[0]", FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax);
  lv21->SetParameter(0, 0.);
  lv21->SetLineColor(kBlue);
  lv21->SetLineWidth(2);
  lv21->SetLineStyle(1);
  lv21->Draw("same");
  TF1 *lv22 = new TF1("lv22", "[0]", FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax);
  lv22->SetParameter(0, 1.);
  lv22->SetLineColor(kBlue);
  lv22->SetLineWidth(2);
  lv22->SetLineStyle(9);
  lv22->Draw("same");
  TF1 *lv23 = new TF1("lv23", "[0]", FlowAnalysis_Fitting::massmin,
                      FlowAnalysis_Fitting::massmax);
  lv23->SetParameter(0, -1.);
  lv23->SetLineColor(kBlue);
  lv23->SetLineWidth(2);
  lv23->SetLineStyle(9);
  lv23->Draw("same");
  pad2_v2->ModifiedUpdate();
  ls->Add(c_v2);

  //////////////////////////////////////////////////////////////////////////////
  ///      Combined plot
  //////////////////////////////////////////////////////////////////////////////
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle();
  TCanvas *c_all = new TCanvas(
      Form(
          "Fitted_%s_Total_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      Form(
          "Fitted_%s_Total_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()));
  if (hs_meanPt_input) {
    c_all->SetCanvasSize(600, 800);
    c_all->cd();
    c_all->Divide(1, 3, 0, 0);
    c_all->cd(1);
    TH1D *hs_se_res_cp = new TH1D();
    TH1D *hs_se_cp = new TH1D();
    if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
        "EventMixing") {
      hs_se_res->Copy(*hs_se_res_cp);
      hs_se_cp = dynamic_cast<TH1D *>(hs_se->Clone());
      hs_se_res_cp->SetStats(0);
      hs_se_res_cp->GetYaxis()->SetTitleSize(0.06);
      hs_se_res_cp->GetYaxis()->SetLabelSize(0.05);
      hs_se_res_cp->GetYaxis()->SetTitleOffset(0.8);
      hs_se_res_cp->GetYaxis()->CenterTitle();
      hs_se_res_cp->Draw("HIST EP");
      hs_se_cp->Draw("HIST EP same");
    } else {
      hs_se->Copy(*hs_se_cp);
      hs_se_cp->SetStats(0);
      hs_se_cp->GetYaxis()->SetTitleSize(0.06);
      hs_se_cp->GetYaxis()->SetLabelSize(0.05);
      hs_se_cp->GetYaxis()->SetTitleOffset(0.8);
      hs_se_cp->GetYaxis()->CenterTitle();
      hs_se_cp->Draw("HIST EP");
    }
    TH1D *hs_me_cp = dynamic_cast<TH1D *>(hs_me->Clone());
    TGraph *gr_model = new TGraph(model);
    gr_model->SetLineWidth(2.0);
    gr_model->SetLineColor(kBlue);
    gr_model->SetTitle("");
    TGraph *gr_bkg_fitted = new TGraph(bkg_fitted);
    gr_bkg_fitted->SetLineWidth(2.0);
    gr_bkg_fitted->SetLineColor(kBlue);
    gr_bkg_fitted->SetLineStyle(kDashed);
    gr_bkg_fitted->SetTitle("");
    TGraph *gr_sig_fitted = new TGraph(bkg_fitted);
    gr_sig_fitted->SetLineWidth(2.0);
    gr_sig_fitted->SetLineColor(kRed);
    gr_sig_fitted->SetTitle("");
    hs_me_cp->Draw("HIST EP same");
    gr_model->Draw("C same");
    gr_bkg_fitted->Draw("C same");
    gr_sig_fitted->Draw("C same");
    TLatex *text_info_yield_all = new TLatex();
    text_info_yield_all->SetTextSize(0.06);
    text_info_yield_all->SetTextFont(42);
    text_info_yield_all->DrawLatexNDC(
        .55, .90,
        "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
        "= 5.36 TeV");
    text_info_yield_all->DrawLatexNDC(
        .55, .83,
        "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < "
        "4");
    text_info_yield_all->SetTextSize(0.07);
    text_info_yield_all->SetTextFont(62);
    text_info_yield_all->DrawLatexNDC(
        .15, .9,
        Form("%g < #it{p}_{T} < %g GeV/c    %g-%g%%",
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax));
    text_info_yield_all->DrawLatexNDC(
        .65, .70,
        Form("%s = %d #pm %d", "N_{J/#psi}",
             int(model->GetParameter(0) / hs_se->GetBinWidth(1)),
             int(model->GetParError(0) / hs_se->GetBinWidth(1))));
    text_info_yield_all->DrawLatexNDC(
        .65, .60, Form("%s = %.2f", "(S/B)_{3#sigma}", S_3sigma / B_3sigma));

    TLegend *legend_yield_all = new TLegend(0.13, 0.05, 0.35, 0.45);
    legend_yield_all->SetBorderSize(0);
    legend_yield_all->SetFillStyle(0);
    legend_yield_all->AddEntry(hs_se_cp, "Data", "P");
    legend_yield_all->AddEntry(hs_me_cp, "Event mixing", "P");
    if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
        "EventMixing") {
      legend_yield_all->AddEntry(hs_se_res_cp, "Residual", "P");
    }
    legend_yield_all->AddEntry(gr_model, "Total");
    legend_yield_all->AddEntry(gr_sig_fitted, "Signal");
    legend_yield_all->AddEntry(gr_bkg_fitted, "Background");
    legend_yield_all->Draw("same");

    c_all->cd(2);
    TH1D *hs_v2se_cp = new TH1D();
    hs_v2se->Copy(*hs_v2se_cp);
    TH1D *hs_v2bkg_cp = dynamic_cast<TH1D *>(hs_v2bkg->Clone());
    TH1D *hs_model_v2_cp = dynamic_cast<TH1D *>(hs_model_v2->Clone());
    hs_v2se_cp->SetStats(0);
    hs_v2se_cp->GetYaxis()->SetTitleSize(0.06);
    hs_v2se_cp->GetYaxis()->SetLabelSize(0.05);
    hs_v2se_cp->GetYaxis()->SetTitleOffset(0.8);
    hs_v2se_cp->GetYaxis()->CenterTitle();
    hs_v2se_cp->Draw("HIST EP");
    lv2->DrawClone("same");
    hs_v2bkg_cp->Draw("HIST EP same");
    hs_model_v2_cp->Draw("HIST same");
    TLegend *legend_v2_all = new TLegend(0.13, 0.45, 0.35, 0.7);
    legend_v2_all->SetBorderSize(0);
    legend_v2_all->SetFillStyle(0);
    legend_v2_all->AddEntry(hs_v2se_cp, "Data", "P");
    legend_v2_all->AddEntry(hs_v2bkg_cp, "Event mixing", "P");
    legend_v2_all->AddEntry(hs_model_v2_cp, "Fit");
    legend_v2_all->Draw("same");

    c_all->cd(3);
    TH1D *hs_meanPt_cp = new TH1D();
    hs_meanPt->Copy(*hs_meanPt_cp);
    TGraph *hs_model_meanPt_cp =
        dynamic_cast<TGraph *>(hs_model_meanPt->Clone());
    hs_meanPt_cp->SetStats(0);
    hs_meanPt_cp->GetYaxis()->SetTitleSize(0.06);
    hs_meanPt_cp->GetYaxis()->SetLabelSize(0.05);
    hs_meanPt_cp->GetYaxis()->SetTitleOffset(0.8);
    hs_meanPt_cp->GetYaxis()->CenterTitle();
    hs_meanPt_cp->GetXaxis()->SetTitleSize(0.05);
    hs_meanPt_cp->GetXaxis()->SetLabelSize(0.05);
    hs_meanPt_cp->Draw("HIST EP");
    hs_model_meanPt_cp->Draw("C same");
    TLegend *legend_meanPt_all = new TLegend(0.13, 0.45, 0.28, 0.65);
    legend_meanPt_all->SetBorderSize(0);
    legend_meanPt_all->SetFillStyle(0);
    legend_meanPt_all->AddEntry(hs_meanPt_cp, "Data", "P");
    legend_meanPt_all->AddEntry(hs_model_meanPt_cp, "Fit");
    legend_meanPt_all->Draw("same");
  } else {
    c_all->SetCanvasSize(600, 600);
    c_all->cd();
    c_all->Divide(1, 2, 0, 0);
    c_all->cd(1);
    TH1D *hs_se_res_cp = new TH1D();
    TH1D *hs_se_cp = new TH1D();
    if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
        "EventMixing") {
      hs_se_res->Copy(*hs_se_res_cp);
      hs_se_cp = dynamic_cast<TH1D *>(hs_se->Clone());
      hs_se_res_cp->SetStats(0);
      hs_se_res_cp->GetYaxis()->SetTitleSize(0.06);
      hs_se_res_cp->GetYaxis()->SetLabelSize(0.05);
      hs_se_res_cp->GetYaxis()->SetTitleOffset(0.8);
      hs_se_res_cp->Draw("HIST EP");
      hs_se_cp->Draw("HIST EP same");
    } else {
      hs_se->Copy(*hs_se_cp);
      hs_se_cp->SetStats(0);
      hs_se_cp->GetYaxis()->SetTitleSize(0.06);
      hs_se_cp->GetYaxis()->SetLabelSize(0.05);
      hs_se_cp->GetYaxis()->SetTitleOffset(0.8);
      hs_se_cp->Draw("HIST EP");
    }
    TH1D *hs_me_cp = dynamic_cast<TH1D *>(hs_me->Clone());
    TGraph *gr_model = new TGraph(model);
    gr_model->SetLineWidth(2.0);
    gr_model->SetLineColor(kBlue);
    gr_model->SetTitle("");
    TGraph *gr_bkg_fitted = new TGraph(bkg_fitted);
    gr_bkg_fitted->SetLineWidth(2.0);
    gr_bkg_fitted->SetLineColor(kBlue);
    gr_bkg_fitted->SetLineStyle(kDashed);
    gr_bkg_fitted->SetTitle("");
    TGraph *gr_sig_fitted = new TGraph(bkg_fitted);
    gr_sig_fitted->SetLineWidth(2.0);
    gr_sig_fitted->SetLineColor(kRed);
    gr_sig_fitted->SetTitle("");
    hs_me_cp->Draw("HIST EP same");
    gr_model->Draw("C same");
    gr_bkg_fitted->Draw("C same");
    gr_sig_fitted->Draw("C same");
    TLatex *text_info_yield_all = new TLatex();
    text_info_yield_all->SetTextSize(0.06);
    text_info_yield_all->SetTextFont(42);
    text_info_yield_all->DrawLatexNDC(
        .55, .9,
        "ALICE Performance, Pb-Pb #sqrt{#it{s}_{NN}} "
        "= 5.36 TeV");
    text_info_yield_all->DrawLatexNDC(
        .55, .83,
        "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < "
        "4");
    text_info_yield_all->SetTextSize(0.07);
    text_info_yield_all->SetTextFont(62);
    text_info_yield_all->DrawLatexNDC(
        .15, .9,
        Form("%g < #it{p}_{T} < %g GeV/c    %g-%g%%",
             FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
             FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax));
    text_info_yield_all->DrawLatexNDC(
        .65, .70,
        Form("%s = %d #pm %d", "N_{J/#psi}",
             int(model->GetParameter(0) / hs_se->GetBinWidth(1)),
             int(model->GetParError(0) / hs_se->GetBinWidth(1))));
    text_info_yield_all->DrawLatexNDC(
        .65, .60, Form("%s = %.2f", "(S/B)_{3#sigma}", S_3sigma / B_3sigma));
    TLegend *legend_yield_all = new TLegend(0.13, 0.05, 0.35, 0.45);
    legend_yield_all->SetBorderSize(0);
    legend_yield_all->SetFillStyle(0);
    legend_yield_all->AddEntry(hs_se_cp, "Data", "P");
    legend_yield_all->AddEntry(hs_me_cp, "Event mixing", "P");
    if (FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg] ==
        "EventMixing") {
      legend_yield_all->AddEntry(hs_se_res_cp, "Residual", "P");
    }
    legend_yield_all->AddEntry(gr_model, "Total");
    legend_yield_all->AddEntry(gr_sig_fitted, "Signal");
    legend_yield_all->AddEntry(gr_bkg_fitted, "Background");
    legend_yield_all->Draw("same");

    c_all->cd(2);
    TH1D *hs_v2se_cp = new TH1D();
    hs_v2se->Copy(*hs_v2se_cp);
    TH1D *hs_v2bkg_cp = dynamic_cast<TH1D *>(hs_v2bkg->Clone());
    TH1D *hs_model_v2_cp = dynamic_cast<TH1D *>(hs_model_v2->Clone());
    hs_v2se_cp->SetStats(0);
    hs_v2se_cp->GetYaxis()->SetTitleSize(0.06);
    hs_v2se_cp->GetYaxis()->SetLabelSize(0.05);
    hs_v2se_cp->GetYaxis()->SetTitleOffset(0.8);
    hs_v2se_cp->GetYaxis()->CenterTitle();
    hs_v2se_cp->Draw("HIST EP");
    lv2->DrawClone("same");
    hs_v2bkg_cp->Draw("HIST EP same");
    hs_model_v2_cp->Draw("HIST same");
    TLegend *legend_v2_all = new TLegend(0.13, 0.45, 0.35, 0.7);
    legend_v2_all->SetBorderSize(0);
    legend_v2_all->SetFillStyle(0);
    legend_v2_all->AddEntry(hs_v2se_cp, "Data", "P");
    legend_v2_all->AddEntry(hs_v2bkg_cp, "Event mixing", "P");
    legend_v2_all->AddEntry(hs_model_v2_cp, "Fit");
    legend_v2_all->Draw("same");
  }
  ls->Add(c_all);

  cout << endl;
  results.emplace_back(model_v2->GetParameter(0));
  results.emplace_back(model_v2->GetParError(0));
  results.emplace_back(model->GetParameter(0) / hs_se->GetBinWidth(1));
  results.emplace_back(model->GetParError(0) / hs_se->GetBinWidth(1));
  results.emplace_back(S_3sigma / B_3sigma);
  results.emplace_back(chi2ndf_mass);
  results.emplace_back(chi2ndf_v2);
  if (hs_meanPt_input) {
    results.emplace_back(meanPt_fitted);
    results.emplace_back(meanPterror_fitted);
    results.emplace_back(chi2ndf_meanPt);
  }
  return results;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Fitting::GetPull(TH1D *hs, TF1 *model, string fit_case) {
  double xmin, xmax;
  model->GetRange(xmin, xmax);
  TH1D *hs_cp = new TH1D();
  hs->Copy(*hs_cp);
  int idx_massmin = hs_cp->FindBin(xmin);
  int idx_massmax = hs_cp->FindBin(xmax);
  int Nbins = hs_cp->GetXaxis()->GetNbins();
  double *Bins = new double[Nbins + 1];
  hs_cp->GetXaxis()->GetLowEdge(Bins);
  Bins[Nbins] = hs_cp->GetXaxis()->GetBinUpEdge(Nbins);
  TH1D *hs_pull = new TH1D(
      Form(
          "Pull_%s_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s", fit_case.c_str(),
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", Nbins, Bins);
  for (int i = 0; i < Nbins; i++) {
    if ((i + 1) >= idx_massmin && (i + 1) <= idx_massmax) {
      double val = (hs_cp->GetBinContent(i + 1) -
                    model->Eval(hs_cp->GetBinCenter(i + 1))) /
                   hs_cp->GetBinError(i + 1);
      hs_pull->SetBinContent(i + 1, val);
    } else {
      hs_pull->SetBinContent(i + 1, 1000.);
    }
  }
  delete hs_cp;
  return hs_pull;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Fitting::GetHistFromTF(TH1D *hs, TF1 *model, string label) {
  int Nbins = hs->GetXaxis()->GetNbins();
  double *Bins = new double[Nbins + 1];
  hs->GetXaxis()->GetLowEdge(Bins);
  Bins[Nbins] = hs->GetXaxis()->GetBinUpEdge(Nbins);
  TH1D *hs_TF = new TH1D(
      Form(
          "Hist_%s_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s", label.c_str(),
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", Nbins, Bins);

  for (int i = 0; i < Nbins; i++) {
    double val = model->Eval(hs->GetBinCenter(i + 1));
    hs_TF->SetBinContent(i + 1, val);
  }
  return hs_TF;
}

//______________________________________________________________________________
TH1D *FlowAnalysis_Fitting::GetV2BkgCorrected(TH1D *hs, TH1D *hs_mepm,
                                              TF1 *bkg) {
  int Nbins = hs->GetXaxis()->GetNbins();
  double *Bins = new double[Nbins + 1];
  hs->GetXaxis()->GetLowEdge(Bins);
  Bins[Nbins] = hs->GetXaxis()->GetBinUpEdge(Nbins);
  TH1D *hs_v2bkg = new TH1D(
      Form(
          "Hist_V2BkgCorrected_%s_v%d%d_%g_%g_%g_%g_%g_%g_%s_%s_%s",
          FlowAnalysis_Fitting::mode_string[FlowAnalysis_Fitting::mode].c_str(),
          FlowAnalysis_Fitting::nhar, FlowAnalysis_Fitting::norder,
          FlowAnalysis_Fitting::ptmin, FlowAnalysis_Fitting::ptmax,
          FlowAnalysis_Fitting::massmin, FlowAnalysis_Fitting::massmax,
          FlowAnalysis_Fitting::centmin, FlowAnalysis_Fitting::centmax,
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_sig]
              .c_str(),
          FlowAnalysis_Fitting::model_string[FlowAnalysis_Fitting::mflag_bkg]
              .c_str(),
          FlowAnalysis_Fitting::v2bkg_string[FlowAnalysis_Fitting::mflag_bkg_v2]
              .c_str()),
      "", Nbins, Bins);

  for (int i = 0; i < Nbins; i++) {
    double N_mepm = hs_mepm->GetBinContent(i + 1);
    double N_bkg = bkg->Eval(hs->GetBinCenter(i + 1));
    double val = hs->GetBinContent(i + 1) * N_mepm / N_bkg;
    double error = hs->GetBinError(i + 1) * N_mepm / N_bkg;
    hs_v2bkg->SetBinContent(i + 1, val);
    hs_v2bkg->SetBinError(i + 1, error);
  }
  return hs_v2bkg;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::Print() {
  cout << "Info for fitter: " << endl;
  cout << "Signal model flag: " << FlowAnalysis_Fitting::mflag_sig << endl;
  cout << "Background model flag: " << FlowAnalysis_Fitting::mflag_bkg << endl;
  cout << "Chi2 threshold for mass fit: " << mchi2max_mass << endl;
  cout << "Chi2 threshold for v2 fit: " << mchi2max_v2 << endl;
}