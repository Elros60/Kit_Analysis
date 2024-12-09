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

#include "FlowAnalysis_Fitting.h"
#include "FlowAnalysis_Helper.h"
#include "Framework/Logger.h"

using namespace std;
vector<double> GetMeanError(int size, double *sample);

void FlowAnalysis_RunByRunResolution(
    string Label, string FileName = "input_AnalysisResults.txt") {

  FlowAnalysis_Helper helper;
  fstream InputFiles;
  vector<vector<double>> stats;
  THStack *hs_stack = new THStack("Run_by_Run_ResolutionFactors", "");
  InputFiles.open(FileName, ios::in);
  TList *ls = new TList();
  TList *ls_run_sp = new TList();
  TList *ls_run_ep = new TList();
  int Nrow = 0;
  int Ncol = 0;
  int NBins_cent = 0;
  int NBins_cent_profile = 0;
  double *Bin_cent = nullptr;
  double *Bin_cent_profile = nullptr;
  if (InputFiles.is_open()) {
    string File;
    cout << "Start reading input AnalysisResults.root ..." << endl;
    while (getline(InputFiles, File)) {
      cout << "Reading input from: " << File << endl;
      TH2F *hs_R2SPAB, *hs_R2SPAC, *hs_R2SPBC;
      TProfile *tp_R2SPAB, *tp_R2SPAC, *tp_R2SPBC;
      TProfile *tp_R2SPAB_Im, *tp_R2SPAC_Im, *tp_R2SPBC_Im;
      TProfile *tp_R2EPAB, *tp_R2EPAC, *tp_R2EPBC;
      TProfile *tp_R2EPAB_Im, *tp_R2EPAC_Im, *tp_R2EPBC_Im;
      vector<string> File_string = helper.tokenize(File);
      string run_number = File_string[File_string.size() - 2];
      helper.LoadReso(File.c_str(), hs_R2SPAB, hs_R2SPAC, hs_R2SPBC);
      helper.LoadResoProfile(File.c_str(), tp_R2SPAB, tp_R2SPAC, tp_R2SPBC,
                             tp_R2EPAB, tp_R2EPAC, tp_R2EPBC, tp_R2SPAB_Im,
                             tp_R2SPAC_Im, tp_R2SPBC_Im, tp_R2EPAB_Im,
                             tp_R2EPAC_Im, tp_R2EPBC_Im);

      // From histograms
      Bin_cent = helper.CreateBinsFromAxis(hs_R2SPAB->GetXaxis());
      NBins_cent = hs_R2SPAB->GetXaxis()->GetNbins();
      Ncol = NBins_cent;
      TH1D *hist_r2sp =
          new TH1D(Form("R2SP_Cent_%s", run_number.c_str()),
                   Form("R_{2}{SP} for run: %s", run_number.c_str()),
                   NBins_cent, Bin_cent);
      hist_r2sp->GetXaxis()->SetTitle("Centrality FT0C(%)");
      hist_r2sp->GetYaxis()->SetTitle("R_{2}{SP}");

      vector<double> row;
      for (int i = 0; i < NBins_cent; i++) {
        // Resolution factor for SP
        TH2F *hs_R2SPAB_cp = dynamic_cast<TH2F *>(hs_R2SPAB->Clone(
            Form("R2SPAB_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
        TH2F *hs_R2SPAC_cp = dynamic_cast<TH2F *>(hs_R2SPAC->Clone(
            Form("R2SPAC_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
        TH2F *hs_R2SPBC_cp = dynamic_cast<TH2F *>(hs_R2SPBC->Clone(
            Form("R2SPBC_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
        hs_R2SPAB_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);
        hs_R2SPAC_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);
        hs_R2SPBC_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);

        double R2SPAB = hs_R2SPAB_cp->GetMean(2);
        double R2SPAC = hs_R2SPAC_cp->GetMean(2);
        double R2SPBC = hs_R2SPBC_cp->GetMean(2);
        double R2SPABe = hs_R2SPAB_cp->GetMean(12);
        double R2SPACe = hs_R2SPAC_cp->GetMean(12);
        double R2SPBCe = hs_R2SPBC_cp->GetMean(12);

        double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
        double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
        double R2SPe =
            R2SPAB * R2SPAC * R2SPBC == 0
                ? 0.0
                : TMath::Sqrt(
                      1. / 4 * (R2SPAC / (R2SPAB * R2SPBC)) * pow(R2SPABe, 2.) +
                      1. / 4 * (R2SPAB / (R2SPAC * R2SPBC)) * pow(R2SPACe, 2.) +
                      1. / 4 * R2SPAC * R2SPAB / pow(R2SPBC, 3.) *
                          pow(R2SPBCe, 2.));
        R2SPe = isnan(R2SPe) || isinf(R2SPe) ? 0. : R2SPe;

        hist_r2sp->SetBinContent(i + 1, R2SP);
        hist_r2sp->SetBinError(i + 1, R2SPe);
        row.emplace_back(R2SP);
        delete hs_R2SPAB_cp;
        delete hs_R2SPAC_cp;
        delete hs_R2SPBC_cp;
      }

      // From profile histograms
      Bin_cent_profile = helper.CreateBinsFromAxis(tp_R2SPAB->GetXaxis());
      NBins_cent_profile = tp_R2SPAB->GetXaxis()->GetNbins();
      TH1D *hist_r2sp_profile =
          new TH1D(Form("%s_R2SP_Cent_FromProfile", run_number.c_str()),
                   Form("R_{2}{SP} for run: %s of %s", run_number.c_str(),
                        Label.c_str()),
                   NBins_cent_profile, Bin_cent_profile);
      hist_r2sp_profile->GetXaxis()->SetTitle("Centrality FT0C(%)");
      hist_r2sp_profile->GetYaxis()->SetTitle("R_{2}{SP}");
      TH1D *hist_r2ep_profile =
          new TH1D(Form("%s_R2EP_Cent_FromProfile", run_number.c_str()),
                   Form("R_{2}{EP} for run: %s of %s", run_number.c_str(),
                        Label.c_str()),
                   NBins_cent_profile, Bin_cent_profile);
      hist_r2ep_profile->GetXaxis()->SetTitle("Centrality FT0C(%)");
      hist_r2ep_profile->GetYaxis()->SetTitle("R_{2}{EP}");

      for (int i = 0; i < NBins_cent_profile; i++) {
        // Resolution factor for SP/EP
        TProfile *tp_R2SPAB_cp = dynamic_cast<TProfile *>(tp_R2SPAB->Clone(
            Form("Profile_R2SPAB_Cent_Copy_%g_%g", Bin_cent_profile[i],
                 Bin_cent_profile[i + 1])));
        TProfile *tp_R2SPAC_cp = dynamic_cast<TProfile *>(tp_R2SPAC->Clone(
            Form("Profile_R2SPAC_Cent_Copy_%g_%g", Bin_cent_profile[i],
                 Bin_cent_profile[i + 1])));
        TProfile *tp_R2SPBC_cp = dynamic_cast<TProfile *>(tp_R2SPBC->Clone(
            Form("Profile_R2SPBC_Cent_Copy_%g_%g", Bin_cent_profile[i],
                 Bin_cent_profile[i + 1])));
        TProfile *tp_R2EPAB_cp = dynamic_cast<TProfile *>(tp_R2EPAB->Clone(
            Form("Profile_R2EPAB_Cent_Copy_%g_%g", Bin_cent_profile[i],
                 Bin_cent_profile[i + 1])));
        TProfile *tp_R2EPAC_cp = dynamic_cast<TProfile *>(tp_R2EPAC->Clone(
            Form("Profile_R2EPAC_Cent_Copy_%g_%g", Bin_cent_profile[i],
                 Bin_cent_profile[i + 1])));
        TProfile *tp_R2EPBC_cp = dynamic_cast<TProfile *>(tp_R2EPBC->Clone(
            Form("Profile_R2EPBC_Cent_Copy_%g_%g", Bin_cent_profile[i],
                 Bin_cent_profile[i + 1])));

        double R2SPAB = tp_R2SPAB_cp->GetBinContent(i + 1);
        double R2SPAC = tp_R2SPAC_cp->GetBinContent(i + 1);
        double R2SPBC = tp_R2SPBC_cp->GetBinContent(i + 1);
        double R2SPABe = tp_R2SPAB_cp->GetBinError(i + 1);
        double R2SPACe = tp_R2SPAC_cp->GetBinError(i + 1);
        double R2SPBCe = tp_R2SPBC_cp->GetBinError(i + 1);

        double R2EPAB = tp_R2EPAB_cp->GetBinContent(i + 1);
        double R2EPAC = tp_R2EPAC_cp->GetBinContent(i + 1);
        double R2EPBC = tp_R2EPBC_cp->GetBinContent(i + 1);
        double R2EPABe = tp_R2EPAB_cp->GetBinError(i + 1);
        double R2EPACe = tp_R2EPAC_cp->GetBinError(i + 1);
        double R2EPBCe = tp_R2EPBC_cp->GetBinError(i + 1);

        double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
        double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
        double R2SPe =
            R2SPAB * R2SPAC * R2SPBC == 0
                ? 0.0
                : TMath::Sqrt(
                      1. / 4 * (R2SPAC / (R2SPAB * R2SPBC)) * pow(R2SPABe, 2.) +
                      1. / 4 * (R2SPAB / (R2SPAC * R2SPBC)) * pow(R2SPACe, 2.) +
                      1. / 4 * R2SPAC * R2SPAB / pow(R2SPBC, 3.) *
                          pow(R2SPBCe, 2.));
        R2SPe = isnan(R2SPe) || isinf(R2SPe) ? 0. : R2SPe;

        double R22EP = R2EPBC != 0 ? R2EPAB * R2EPAC / R2EPBC : 0.0;
        double R2EP = R22EP > 0 ? TMath::Sqrt(R22EP) : 0.0;
        double R2EPe =
            R2EPAB * R2EPAC * R2EPBC == 0
                ? 0.0
                : TMath::Sqrt(
                      1. / 4 * (R2EPAC / (R2EPAB * R2EPBC)) * pow(R2EPABe, 2.) +
                      1. / 4 * (R2EPAB / (R2EPAC * R2EPBC)) * pow(R2EPACe, 2.) +
                      1. / 4 * R2EPAC * R2EPAB / pow(R2EPBC, 3.) *
                          pow(R2EPBCe, 2.));
        R2EPe = isnan(R2EPe) || isinf(R2EPe) ? 0. : R2EPe;

        hist_r2sp_profile->SetBinContent(i + 1, R2SP);
        hist_r2sp_profile->SetBinError(i + 1, R2SPe);
        hist_r2ep_profile->SetBinContent(i + 1, R2EP);
        hist_r2ep_profile->SetBinError(i + 1, R2EPe);

        delete tp_R2SPAB_cp;
        delete tp_R2SPAC_cp;
        delete tp_R2SPBC_cp;
        delete tp_R2EPAB_cp;
        delete tp_R2EPAC_cp;
        delete tp_R2EPBC_cp;
      }

      ls_run_sp->Add(hist_r2sp_profile);
      ls_run_ep->Add(hist_r2ep_profile);
      ls->Add(hist_r2sp);
      hs_stack->Add(hist_r2sp);
      stats.emplace_back(row);
      Nrow++;

      delete hs_R2SPAB;
      delete hs_R2SPAC;
      delete hs_R2SPBC;
      delete tp_R2SPAB;
      delete tp_R2SPAC;
      delete tp_R2SPBC;
      delete tp_R2SPAB_Im;
      delete tp_R2SPAC_Im;
      delete tp_R2SPBC_Im;
      delete tp_R2EPAB;
      delete tp_R2EPAC;
      delete tp_R2EPBC;
      delete tp_R2EPAB_Im;
      delete tp_R2EPAC_Im;
      delete tp_R2EPBC_Im;
    }
  }

  vector<double *> stats_trans;
  for (int i = 0; i < Ncol; i++) {
    stats_trans.emplace_back(new double[Nrow]);
  }
  for (int i = 0; i < Nrow; i++) {
    for (int j = 0; j < Ncol; j++) {
      stats_trans[j][i] = stats[i][j];
    }
  }

  TH1D *hist_final =
      new TH1D(Form("R2SP_All_%s", Label.c_str()),
               Form("R2SP_All_%s", Label.c_str()), NBins_cent, Bin_cent);
  for (int i = 0; i < Ncol; i++) {
    vector<double> final_stats = GetMeanError(Nrow, stats_trans[i]);
    hist_final->SetBinContent(i + 1, final_stats[0]);
    hist_final->SetBinError(i + 1, final_stats[1]);
  }
  TList *ls_final = new TList();
  ls_final->Add(hist_final);
  ls_final->Add(hs_stack);

  TFile froot(Form("ResolutionsAll_%s.root", Label.c_str()), "RECREATE");
  ls->Write("RunByRun_Histogram", TObject::kSingleKey);
  ls_run_sp->Write("RunByRun_HistogramFromProfile_SP", TObject::kSingleKey);
  ls_run_ep->Write("RunByRun_HistogramFromProfile_EP", TObject::kSingleKey);
  ls_final->Write("FinalStatistics", TObject::kSingleKey);
  froot.Close();
  InputFiles.close();
}

vector<double> GetMeanError(int size, double *sample) {
  vector<double> results;
  double mean = 0.;
  for (int i = 0; i < size; i++) {
    mean += sample[i] / size;
  }
  results.emplace_back(mean);

  double sum2 = 0;
  for (int i = 0; i < size; i++) {
    sum2 += pow(sample[i] - mean, 2.) / size;
  }
  double rms = pow(sum2, 0.5);
  results.emplace_back(rms);
  return results;
}