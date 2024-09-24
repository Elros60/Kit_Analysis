#ifndef FLOWANALYSIS_FITTING_H
#define FLOWANALYSIS_FITTING_H

#include <Math/IntegratorOptions.h>
#include <Math/MinimizerOptions.h>
#include <iostream>

#include <TCanvas.h>
#include <TF1.h>
#include <TF1NormSum.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

using namespace std;
using namespace ROOT::Math;

// CB2 functions will be parametrised with 2 sets of tails: MC and data
enum ModelType { CB2 = 0, CB2Bis, NA60, Chebychev, VWG, Exp2, PolExp };

class FlowAnalysis_Fitting {
public:
  void init();
  void setModel(int flag_sig, int flag_bkg);
  void setModelV2(int flag_bkg_v2);
  void setChi2Max(double chi2) { mchi2max = chi2; };
  void setMassRange(double mass_min, double mass_max);
  void setCentRange(double cent_min, double cent_max);
  void setHarmonic(int har);
  void setMode(int mode_flag) { mode = mode_flag; };
  void setOrder(int order);
  vector<double> runFitting(TH1D *hs_input, TH1D *hs_v2_input, TList *ls,
                            double ptmin, double ptmax);
  void Print();

private:
  void CreateModel(TF1 *&model, int flag);
  static double DoubleSidedCB2(double x, double mu, double width, double a1,
                               double p1, double a2, double p2);
  static double DoubleSidedCB(double *x, double *par);
  static double Cheby7(double *x, double *par);
  static double Cheby3(double *x, double *par);
  static double VariableWidthGauss(double *x, double *par);
  static double DoubleExp(double *x, double *par);
  static double PolyExp(double *x, double *par);
  static double FittedSignal(double *x, double *par);
  static double FittedBkg(double *x, double *par);

  static double massmin;
  static double massmax;
  static double centmin;
  static double centmax;
  static int mflag_sig;
  static int mflag_bkg;
  static int mflag_bkg_v2;
  static int norder;
  static int nhar;
  static int mode;
  static string mode_string[2];

  double mchi2max{1.};
};
#endif