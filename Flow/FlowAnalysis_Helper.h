#ifndef FLOWANALYSIS_HELPER_H
#define FLOWANALYSIS_HELPER_H

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
#include <TLatex.h>
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

using namespace std;
using namespace ROOT::Math;

class FlowAnalysis_Helper {
public:
  FlowAnalysis_Helper() = default;

  void LoadDataRun2(double *&x, double *&y, double *&ex, double *&ey,
                    double *&ey_sys, int flag);
  double *CreateBinsFromAxis(TAxis *axis);
  void CreateBins(double *axis, double min, double max, int Nbins = 10);
  void LoadData(std::string FileName, THnSparse *&hs_V2, TH2F *&hs_R2SPAB,
                TH2F *&hs_R2SPAC, TH2F *&hs_R2SPBC, std::string muonCut,
                std::string dimuonCut);
  void LoadDataME(std::string FileName, THnSparse *&hs_V2SEPM,
                  THnSparse *&hs_V2SEPP, THnSparse *&hs_V2SEMM,
                  THnSparse *&hs_V2MEPM, THnSparse *&hs_V2MEPP,
                  THnSparse *&hs_V2MEMM, TH2F *&hs_R2SPAB, TH2F *&hs_R2SPAC,
                  TH2F *&hs_R2SPBC, TH3F *&hs_u2q2MEPM1, TH3F *&hs_u2q2MEPP1,
                  TH3F *&hs_u2q2MEMM1, TH3F *&hs_u2q2MEPM2, TH3F *&hs_u2q2MEPP2,
                  TH3F *&hs_u2q2MEMM2, TH3F *&hs_r2spMEPM1, TH3F *&hs_r2spMEPP1,
                  TH3F *&hs_r2spMEMM1, TH3F *&hs_r2spMEPM2, TH3F *&hs_r2spMEPP2,
                  TH3F *&hs_r2spMEMM2, TH2F *&hs_cosDeltaPhiMEPM1,
                  TH2F *&hs_cosDeltaPhiMEPP1, TH2F *&hs_cosDeltaPhiMEMM1,
                  TH2F *&hs_cosDeltaPhiMEPM2, TH2F *&hs_cosDeltaPhiMEPP2,
                  TH2F *&hs_cosDeltaPhiMEMM2, std::string muonCut,
                  std::string dimuonCut);
  TH1D *GetMass(double ptmin, double ptmax, double massmin, double massmax,
                double centmin, double centmax, THnSparse *hist_V2);
  TH1D *GetV2(double ptmin, double ptmax, double massmin, double massmax,
              double centmin, double centmax, THnSparse *hist_V2, double R2SP);
  TH1D *GetRfactor(double ptmin, double ptmax, double massmin, double massmax,
                   double centmin, double centmax, THnSparse *hs_V2MEPM,
                   THnSparse *hs_V2MEPP, THnSparse *hs_V2MEMM);
  double GetFfactor(double ptmin, double ptmax, double massmin, double massmax,
                    double centmin, double centmax, THnSparse *hs_V2SEPP,
                    THnSparse *hs_V2SEMM, THnSparse *hs_V2MEPM,
                    TH1D *hist_rfactor);
  vector<double> GetStats(int size, double *sample, double *sample_error);
  void PlotSNRvsRun2(int size_ptbin, double *pt_bins, int size_run3,
                     double *x_run3, double *snr_run3, int size_run2,
                     double *x_run2, double *snr_run2, TList *ls);
  void PlotSystematics(int index, TCanvas *c_sys_yield, TCanvas *c_sys_v2,
                       TH1D *hist_sys_yield, TH1D *hist_sys_v2,
                       double *bins_sys_yield, double *bins_sys_v2,
                       double *chi2_yield, double *chi2_v2, int nbCombo_yield,
                       int nbCombo_v2, vector<double> stats_yield,
                       vector<double> stats_v2, double *pt_bins,
                       TList *ls_sys_yield, TList *ls_sys_v2);
  void PlotFinalResults(int size_ptbin, double *pt_bins, double *x_v2pt,
                        double *y_v2pt, double *ex_v2pt, double *ey_v2pt,
                        double *eysys_v2pt, double *x_run2, double *y_run2,
                        double *ex_run2, double *ey_run2, double *eysys_run2,
                        double *x_yield, double *y_yield, double *ex_yield,
                        double *ey_yield, double *eysys_yield,
                        double *x_yield_run2, double *y_yield_run2,
                        double *ex_yield_run2, double *ey_yield_run2,
                        double *eysys_yield_run2, TList *ls);
};
#endif