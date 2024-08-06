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
#include <TH1F.h>
#include <TH2F.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMatrixD.h>
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

#include <RooAbsData.h>
#include <RooAddPdf.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include "Framework/Logger.h"

using namespace std;
using namespace RooFit;

double *CreateBinsFromAxis(TAxis *axis);
void CreateBins(double *axis, double min, double max, int Nbins = 10);

void FlowAnalysis_FitTemplate(
    std::string FileName = "AnalysisResults_Flow.root") {

  // Load data from AnalysisResults.root
  TFile *Input_File = TFile::Open(FileName.c_str());
  TList *list_hist = (TList *)Input_File->Get("DifferentialFlow_2_3");
  TH1D *hist_mass = (TH1D *)list_hist->FindObject(
      "Mass_Pt_Rapidity_centrFT0C_Copy_2_3_proj_0");

  // Setting up RooFit
  RooRealVar m("m", "mass", 2.5, 4.5);
  RooDataHist dh("dh", "dh", m, Import(*hist_mass));
  RooPlot *frame =
      m.frame(Name("mFrame"), Title("Imported Invariant Mass histogram"));
  dh.plotOn(frame, DataError(RooAbsData::Poisson), Name("data"));

  // Setting up fit model
  // Signal model: double-sided crytalball
  RooRealVar m0("m0", "m0", 3.097, 3.05, 3.12);
  RooRealVar sigma("sigma", "sigma", 0.1, 0., 3.0);
  RooRealVar alphaL("alphaL", "alphaL", 10.0, 5.0, 30.0);
  RooRealVar alphaR("alphaR", "alphaR", 15.0, 5.0, 30.0);
  RooRealVar nL("nL", "nL", 1.0, 0.1, 10.0);
  RooRealVar nR("nR", "nR", 1.0, 0.1, 10.0);
  RooCrystalBall sig("Signal", "CB2", m, m0, sigma, alphaL, nL, alphaR, nR);
  // Background model: VWG
  RooRealVar A("A", "A", -0.1, -2.5, 2.5);
  RooRealVar B("B", "B", -0.1, -5.0, 5.0);
  RooRealVar C("C", "C", -0.1, -5.0, 5.0);
  RooFormulaVar sigma_VWG("sigma_VWG", "sigma_VWG", "B+C*((m-A)/A)",
                          RooArgList(m, A, B, C));

  RooGenericPdf bkg("Bkg", "VWG", "exp(-(m-A)*(m-A)/(2.*sigma_VWG*sigma_VWG))",
                    RooArgSet(m, A, sigma_VWG));
  // Construct a combined signal+background model
  RooRealVar nsig("nsig", "#signal events", 200, 0., 1000000);
  RooRealVar nbkg("nbkg", "#background events", 800, 0., 1000000);
  RooAddPdf model("model", "CB2+VWG", {sig, bkg}, {nsig, nbkg});

  // Do fitting
  RooFitResult *fit_result = model.fitTo(dh, Extended(kTRUE), Save(kTRUE));
  model.plotOn(frame, LineColor(kBlue), Name("model"));
  model.plotOn(frame, Components(sig), LineColor(kRed), Name("Signal"));
  model.plotOn(frame, Components(bkg), LineStyle(ELineStyle::kDashed),
               LineColor(kBlue), Name("Bkg"));
  double chi2 = frame->chiSquare("model", "data");
  TString chi2text = TString::Format("#chi^{2}/ndf = %f ", chi2);
  model.paramOn(frame, Label(chi2text), Layout(0.63, 0.9, 0.9), Format("NE"));

  frame->Draw();
}

double *CreateBinsFromAxis(TAxis *axis) {
  int Nbins = axis->GetNbins();
  double *Bins = new double[Nbins + 1];
  axis->GetLowEdge(Bins);
  Bins[Nbins] = axis->GetBinUpEdge(Nbins);
  return Bins;
}

void CreateBins(double *axis, double min, double max, int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}