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

#include "FlowAnalysis_Fitting.h"

#include "Framework/Logger.h"

using namespace std;
using namespace RooFit;

double *CreateBinsFromAxis(TAxis *axis);
void CreateBins(double *axis, double min, double max, int Nbins = 10);

void FlowAnalysis_ScalarProduct(nt flag_sig, int flag_bkg, double chi2max,
                                std::string FileName = "AnalysisResults.root",
                                std::string muonCut = "muonLowPt210SigmaPDCA",
                                std::string dimuonCut = "pairRapidityForward") {
  // Load data from AnalysisResults.root
  TFile *Input_File = TFile::Open(FileName.c_str());
  THashList *list_hist_v2 =
      (THashList *)Input_File->Get("analysis-same-event-pairing/output");
  TList *sublist_v2 = (TList *)list_hist->FindObject(
      Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));
  THashList *list_hist_r2 =
      (THashList *)Input_File->Get("analysis-event-selection/output");
  TList *sublist_r2 = (TList *)list_hist->FindObject("Event_AfterCuts");

  // Initialization for fitting
  FlowAnalysis_Fitting fitter;
  fitter.init();
  fitter.setModel(flag_sig, flag_bkg);
  fitter.setChi2Max(chi2max);
  fitter.Print();

  // Get histograms of correlations and resolution factors
  THnSparse *hs_V2 = (THnSparse *)sublist->FindObject("Mass_Pt_centrFT0C_V2");
  TH2F *hs_R2AB = (TH2F *)sublist->FindObject("R2SP_TPCFT0A_CentFT0C");
  TH2F *hs_R2AC = (TH2F *)sublist->FindObject("R2SP_TPCFT0C_CentFT0C");
  TH2F *hs_R2BC = (TH2F *)sublist->FindObject("R2SP_FT0AFT0C_CentFT0C");

  // Get binning information
  TAxis *massAxis = hs_V2->GetXaxis();
  TAxis *ptAxis = hs_V2->GetYaxis();
  TAxis *centAxis = hs_V2->GetZaxis();
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_mass = CreateBinsFromAxis(massAxis);
  double *Bin_pt = CreateBinsFromAxis(ptAxis);
  double *Bin_cent = CreateBinsFromAxis(centAxis);

  ///////////////////////////////////////////////////
  ///                                             ///
  ///   Analysis for Differential Flow of J/Psi   ///
  ///                                             ///
  ///////////////////////////////////////////////////

  // Define variables' range for analysis
  double Bin_pt_mass[11] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
  double cent_min = 10.0; // 20
  double cent_max = 50.0; // 50
  double mass_min = 2.3;  // same as initial setup
  double mass_max = 4.2;  // same as initial setup
  int iBin_cent_min = centAxis->FindBin(cent_min);
  int iBin_cent_max = centAxis->FindBin(cent_max);
  int iBin_mass_min = massAxis->FindBin(mass_min);
  int iBin_mass_max = massAxis->FindBin(mass_max);

  hs_R2AB->GetXaxis()->SetRange(iBin_cent_min, iBin_cent_max);
  hs_R2AC->GetXaxis()->SetRange(iBin_cent_min, iBin_cent_max);
  hs_R2BC->GetXaxis()->SetRange(iBin_cent_min, iBin_cent_max);

  // Calculation of resolution factor in the given centrality bin
  double R2AB = hs_R2AB->GetMean(2);
  double R2AC = hs_R2AB->GetMean(2);
  double R2BC = hs_R2AB->GetMean(2);
  double R22 = R2BC != 0 ? R2AB * R2AC / R2BC : 0.0;
  double R2 = R22 > 0 ? TMath::Sqrt(R22) : 0.0;
  if (R2 == 0.0) {
    LOG(fatal) << "Inconsistent value of resolution factor!";
  }

  for (int i = 0; i < int(size(Bin_pt_mass)) - 1; i++) {

    TList *l_diff = new TList();
    int iBin_pt_min = ptAxis->FindBin(Bin_pt_mass[i]);
    int iBin_pt_max = ptAxis->FindBin(Bin_pt_mass[i + 1]);

    // Copy original profiles for projections
    THnSparse *hs_V2_cp = dynamic_cast<THnSparse *>(
        hs_V2->Clone(Form("Mass_Pt_centrFT0C_V2_Copy_%g_%g", Bin_pt_mass[i],
                          Bin_pt_mass[i + 1])));

    // Set axes' ranges for mass-differential study
    hs_V2_cp->GetXaxis()->SetRange(iBin_mass_min, iBin_mass_max);
    hs_V2_cp->GetYaxis()->SetRange(iBin_pt_min, iBin_pt_max);
    hs_V2_cp->GetZaxis()->SetRange(iBin_cent_min, iBin_cent_max);

    // Get mass and v2
    TH1D *hs_mass_proj = hs_V2_cp->Projection(0);
    TH2D *hs_v2_proj = hs_V2_cp->Projection(0, 3);

    // Define histograms
    double *Bin_mass_new = CreateBinsFromAxis(hs_V2_cp->GetXaxis());
    int NBins_mass_new = hs_V2_cp->GetXaxis()->GetNbins();
    TH1D *hist_v2 = new TH1D(
        Form("v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        Form("v^{J/#psi}_{2}{SP}_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        NBins_mass_new, Bin_mass_new);
    hist_v2->GetXaxis()->SetTitle("mass (GeV/c2)");
    hist_v2->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");

    // Evaluation of differential flow as function of invariant mass
    for (int j = 0; j < NBins_mass_new; j++) {
      TH2D *hs_v2_proj_copy =
          dynamic_cast<TH2D *>(hs_v2_proj->Clone("Mass_V2_Copy"));
      hs_v2_proj_copy->GetXaxis()->SetRange(Bin_mass_new[j],
                                            Bin_mass_new[j + 1]);
      double v2 = hs_v2_proj_copy->GetMean(2);
      double v2e = hs_v2_proj_copy->GetMean(12);

      hist_v2->SetBinContent(j + 1, v2 / R2);
      hist_v2->SetBinError(j + 1, v2e / R2);
    }

    // Do fitting
    fitter.runFitting(hs_mass_proj, l_diff, Bin_pt_mass[i], Bin_pt_mass[i + 1],
                      mass_min, mass_max);

    // Save results
    l_diff->Add(hs_mass_proj);
    l_diff->Add(hist_v2);
    l_diff->Write(
        Form("DifferentialFlow_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
  }
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