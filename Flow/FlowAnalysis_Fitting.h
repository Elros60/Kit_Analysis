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

enum ModelType { CB2 = 0, Chebychev, VWG, Exp2, PolExp };

class FlowAnalysis_Fitting {
public:
  void init();
  void setModel(int flag_sig, int flag_bkg);
  void setChi2Max(double chi2) { mchi2max = chi2; };
  void setMassRange(double mass_min, double mass_max);
  void setCentRange(double cent_min, double cent_max);
  vector<double> runFitting(TH1D *hs, TH1D *hs_v2, TList *ls, double ptmin,
                            double ptmax);
  void Print();

private:
  void CreateModel(TF1 *&model, int flag);
  static double DoubleSidedCB2(double x, double mu, double width, double a1,
                               double p1, double a2, double p2);
  static double DoubleSidedCB(double *x, double *par);
  static double Cheby(double *x, double *par);
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

  double mchi2max{1.};
};
#endif