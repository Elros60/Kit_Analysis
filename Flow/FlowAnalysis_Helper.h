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

  vector<string> tokenize(string input_string);
  double *CreateBinsFromAxis(TAxis *axis);
  void CreateBins(double *axis, double min, double max, int Nbins = 10);
  void LoadReso(std::string FileName, TH2F *&hs_R2SPAB, TH2F *&hs_R2SPAC,
                TH2F *&hs_R2SPBC);
  void LoadResoProfile(std::string FileName, TProfile *&tp_R2SPAB,
                       TProfile *&tp_R2SPAC, TProfile *&tp_R2SPBC,
                       TProfile *&tp_R2EPAB, TProfile *&tp_R2EPAC,
                       TProfile *&tp_R2EPBC, TProfile *&tp_R2SPAB_Im,
                       TProfile *&tp_R2SPAC_Im, TProfile *&tp_R2SPBC_Im,
                       TProfile *&tp_R2EPAB_Im, TProfile *&tp_R2EPAC_Im,
                       TProfile *&tp_R2EPBC_Im);
  void LoadData(std::string FileName, THnSparse *&hs_V2, TH2F *&hs_R2SPAB,
                TH2F *&hs_R2SPAC, TH2F *&hs_R2SPBC, std::string muonCut,
                std::string dimuonCut);
  void LoadDataME(
      std::string FileName, THnSparse *&hs_V2SEPM, THnSparse *&hs_V2SEPP,
      THnSparse *&hs_V2SEMM, THnSparse *&hs_V2MEPM, THnSparse *&hs_V2MEPP,
      THnSparse *&hs_V2MEMM, TH2F *&hs_R2SPAB, TH2F *&hs_R2SPAC,
      TH2F *&hs_R2SPBC, THnSparse *&hs_u2q2_cosDeltaPhi_MEPM1,
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
      std::string dimuonCut);
  void LoadDataMEProfile(std::string FileName, TProfile3D *&tp_V2SEPM,
                         TProfile3D *&tp_V2SEPP, TProfile3D *&tp_V2SEMM,
                         TProfile3D *&tp_V2MEPM, TProfile3D *&tp_V2MEPP,
                         TProfile3D *&tp_V2MEMM, std::string muonCut,
                         std::string dimuonCut);
  void LoadDataRun2(double *&x, double *&y, double *&ex, double *&ey,
                    double *&ey_sys, int flag);
  void LoadDataRun2Cent(double *&x, double *&y, double *&ex, double *&ey,
                        double *&ey_sys, int flag);
  void LoadDataYieldRun2(double *&x, double *&y, double *&ex, double *&ey,
                         double *&ey_sys, double *&SNR, int flag);
  TH1D *GetMass(double ptmin, double ptmax, double massmin, double massmax,
                double centmin, double centmax, THnSparse *hist_V2,
                std::string flag = "SEPM");
  TH1D *GetMassProfile(double ptmin, double ptmax, double massmin,
                       double massmax, double centmin, double centmax,
                       TProfile3D *tp_V2, std::string flag = "SEPM");
  TH1D *GetMeanPt(double ptmin, double ptmax, double massmin, double massmax,
                  double centmin, double centmax, TProfile3D *tp_V2,
                  std::string flag = "SEPM");
  TH1D *GetV2(double ptmin, double ptmax, double massmin, double massmax,
              double centmin, double centmax, THnSparse *hist_V2, double R2SP,
              std::string flag = "SEPM");
  TH1D *GetV2Profile(double ptmin, double ptmax, double massmin, double massmax,
                     double centmin, double centmax, TProfile3D *tp_V2,
                     std::string flag = "SEPM");
  TH1D *GetV2EM(double ptmin, double ptmax, double massmin, double massmax,
                double centmin, double centmax,
                THnSparse *hs_u2q2_cosDeltaPhi_MEPM1,
                THnSparse *hs_u2q2_cosDeltaPhi_MEPM2, TH3F *hs_r2spAB1,
                TH3F *hs_r2spAC1, TH3F *hs_r2spBC1, TH3F *hs_r2spAB2,
                TH3F *hs_r2spAC2, TH3F *hs_r2spBC2, std::string flag);
  TH1D *GetRfactor(double ptmin, double ptmax, double massmin, double massmax,
                   double centmin, double centmax, THnSparse *hs_V2MEPM,
                   THnSparse *hs_V2MEPP, THnSparse *hs_V2MEMM);
  TH1D *GetRfactorProfile(double ptmin, double ptmax, double massmin,
                          double massmax, double centmin, double centmax,
                          TProfile3D *tp_V2MEPM, TProfile3D *tp_V2MEPP,
                          TProfile3D *tp_V2MEMM);
  double GetFfactor(double ptmin, double ptmax, double massmin, double massmax,
                    double centmin, double centmax, THnSparse *hs_V2SEPP,
                    THnSparse *hs_V2SEMM, THnSparse *hs_V2MEPM,
                    TH1D *hist_rfactor);
  double GetFfactorProfile(double ptmin, double ptmax, double massmin,
                           double massmax, double centmin, double centmax,
                           TProfile3D *tp_V2SEPP, TProfile3D *tp_V2SEMM,
                           TProfile3D *tp_V2MEPM, TH1D *hist_rfactor);
  vector<double> GetStats(int size, double *sample, double *sample_error);
  void PlotSNRvsRun2(int size_ptbin, double *pt_bins, int size_run3,
                     double *x_run3, double *snr_run3, int size_run2,
                     double *x_run2, double *snr_run2, TList *ls);
  void PlotSystematics(double ptmin, double ptmax, double centmin,
                       double centmax, int nBins, int index,
                       TH1D *hist_sys_yield, TH1D *hist_sys_v2,
                       TH1D *hist_sys_meanPt, double *bins_sys_yield,
                       double *bins_sys_v2, double *chi2_yield, double *chi2_v2,
                       double *chi2_meanPt, int nbCombo_yield, int nbCombo_v2,
                       vector<double> stats_yield, vector<double> stats_v2,
                       vector<double> stats_meanPt, double *pt_bins,
                       TList *ls_sys_yield, TList *ls_sys_v2,
                       TList *ls_sys_meanPt, bool Save,
                       std::string flag = "pt");
  void PlotFinalResults(int size_ptbin, double cent_min, double cent_max,
                        double *pt_bins, double *x_v2pt, double *y_v2pt,
                        double *ex_v2pt, double *ey_v2pt, double *eysys_v2pt,
                        double *x_run2, double *y_run2, double *ex_run2,
                        double *ey_run2, double *eysys_run2, double *x_yield,
                        double *y_yield, double *ex_yield, double *ey_yield,
                        double *eysys_yield, double *x_yield_run2,
                        double *y_yield_run2, double *ex_yield_run2,
                        double *ey_yield_run2, double *eysys_yield_run2,
                        TList *ls);
  void PlotFinalResultsCent(int size_centbin, double pt_min, double pt_max,
                            double *cent_bins, double *x_v2cent,
                            double *y_v2cent, double *ex_v2cent,
                            double *ey_v2cent, double *eysys_v2cent,
                            double *x_run2, double *y_run2, double *ex_run2,
                            double *ey_run2, double *eysys_run2,
                            double *x_yield, double *y_yield, double *ex_yield,
                            double *ey_yield, double *eysys_yield, TList *ls);
  void PlotSEME(std::string flag, std::string type, double ptmin, double ptmax,
                double massmin, double massmax, double centmin, double centmax,
                TH1D *hist_SE, TH1D *hist_ME, TList *ls);
};
#endif