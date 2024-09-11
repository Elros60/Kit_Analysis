#ifndef FLOWANALYSIS_FITTING_H
#define FLOWANALYSIS_FITTING_H

#include <TCanvas.h>
#include <TH1F.h>
#include <TList.h>

#include <RooAbsData.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooChebychev.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooPowerSum.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

using namespace std;
using namespace RooFit;

enum ModelType { CB2 = 0, Chebychev, VWG, POL, Exp2, PolExp };

class FlowAnalysis_Fitting {
public:
  FlowAnalysis_Fitting() : mWS("w", "workspace") {};
  void init();
  void setModel(int flag_sig, int flag_bkg);
  void setChi2Max(double chi2) { mchi2max = chi2; };
  void runFitting(TH1D *hs, TList *ls, double ptmin, double ptmax,
                  double massmin, double massmax);

private:
  void CreateModel(RooRealVar x, int flag);
  int mflag_sig{};
  int mflag_bkg{};
  double mchi2max{};
  RooWorkspace mWS{};
};
#endif