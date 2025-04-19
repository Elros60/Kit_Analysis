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
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

using namespace std;
using namespace ROOT::Math;

// CB2 functions will be parametrised with 2 sets of tails: MC and data
enum ModelType {
  CB2Data = 0,
  CB2MC,
  NA60,
  Chebychev,
  EventMixing_Landau,
  EventMixing_GausExp,
  VWG,
  Exp2,
  PolExp
};

class FlowAnalysis_Fitting {
public:
  void init();
  void setModel(int flag_sig, int flag_bkg);
  void setModelV2(int flag_bkg_v2);
  void setChi2MaxMass(double chi2) { mchi2max_mass = chi2; };
  void setChi2MaxV2(double chi2) { mchi2max_v2 = chi2; };
  void setMassRange(double mass_min, double mass_max);
  void setCentRange(double cent_min, double cent_max);
  void setPtRange(double pt_min, double pt_max);
  void setHarmonic(int har);
  void setMode(int mode_flag) { mode = mode_flag; };
  void setOrder(int order);
  TH1D *GetPull(TH1D *hs, TF1 *model, string fit_case);
  TH1D *GetHistFromTF(TH1D *hs, TF1 *model, string label);
  TH1D *GetV2BkgCorrected(TH1D *hs, TH1D *hs_mepm, TF1 *bkg);
  vector<double> runFitting(TH1D *hs_input, TH1D *hs_v2_input, TList *ls);
  vector<double> runFittingEM(TH1D *hs_mse_input, TH1D *hs_mme_input,
                              TH1D *hs_v2se_input, TH1D *hs_v2me_input,
                              TH1D *hs_meanPt_input, TList *ls);
  vector<double> runFittingMassOnly(TH1D *hs_input, TList *ls);
  void Print();

private:
  void CreateModel(TF1 *&model, int flag);
  void GetNparFullModel(int &nParSig, int &nParBkg);
  void InitFullModel(TF1 *&model);
  static double DoubleSidedCB2(double x, double mu, double width, double a1,
                               double p1, double a2, double p2);
  static double DoubleSidedCB(double *x, double *par);
  static double NA60Function(double *x, double *par);
  static double Cheby7(double *x, double *par);
  static double Cheby3(double *x, double *par);
  static double Cheby4(double *x, double *par);
  static double Cheby5(double *x, double *par);
  static double Cheby6(double *x, double *par);
  static double VariableWidthGauss(double *x, double *par);
  static double DoubleExp(double *x, double *par);
  static double PolyExp(double *x, double *par);
  static double FittedSignal(double *x, double *par);
  static double FittedBkg(double *x, double *par);
  static double FullModelWithPsi2s(double *x, double *par);

  static double massmin;
  static double massmax;
  static double centmin;
  static double centmax;
  static double ptmin;
  static double ptmax;
  static int mflag_sig;
  static int mflag_bkg;
  static int mflag_bkg_v2;
  static int norder;
  static int nhar;
  static int mode;
  static string mode_string[2];
  static string model_string[8];
  static string v2bkg_string[2];

  double mchi2max_mass{1.};
  double mchi2max_v2{1.};
};
#endif