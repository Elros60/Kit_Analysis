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
void runFitting(TH1D *hs, TList *ls, double ptmin, double ptmax);

void FlowAnalysis_Cumulants(std::string muonCut = "muonLowPt210SigmaPDCA",
                            std::string dimuonCut = "pairRapidityForward",
                            std::string FileName = "AnalysisResults.root") {

  // Load data from AnalysisResults.root
  TFile *Input_File = TFile::Open(FileName.c_str());
  THashList *list_hist =
      (THashList *)Input_File->Get("analysis-same-event-pairing/output");
  TList *sublist = (TList *)list_hist->FindObject(
      Form("PairsMuonSEPM_%s_%s", muonCut.c_str(), dimuonCut.c_str()));

  // Get TProfiles of correlations
  TProfile3D *tp_Corr2Ref =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr2REF");
  TProfile3D *tp_Corr4Ref =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr4REF");
  TProfile3D *tp_Corr2Poi =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr2POI");
  TProfile3D *tp_Corr4Poi =
      (TProfile3D *)sublist->FindObject("Mass_Pt_centrFT0C_Corr4POI");

  // Get histograms for J/Psi invariant mass
  THnSparse *hist_mass =
      (THnSparse *)sublist->FindObject("Mass_Pt_Rapidity_CentFT0C");

  // Get binning information
  TAxis *massAxis = tp_Corr2Ref->GetXaxis();
  TAxis *ptAxis = tp_Corr2Ref->GetYaxis();
  TAxis *centAxis = tp_Corr2Ref->GetZaxis();
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  double *Bin_mass = CreateBinsFromAxis(massAxis);
  double *Bin_pt = CreateBinsFromAxis(ptAxis);
  double *Bin_cent = CreateBinsFromAxis(centAxis);

  ///////////////////////////////////////
  ///                                 ///
  ///   Analysis for Reference Flow   ///
  ///                                 ///
  ///////////////////////////////////////

  // Define projected profiles w.r.t above defined variables' ranges
  TList *l_refflow = new TList();
  TFile f("FlowAnalysisResults_Cumulants.root", "RECREATE");
  TProfile *tp_Corr2Ref_cent = tp_Corr2Ref->Project3DProfile("xz")->ProfileX();
  TProfile *tp_Corr4Ref_cent = tp_Corr4Ref->Project3DProfile("xz")->ProfileX();

  TH1D *hist_c22ref =
      new TH1D("c22ref", "c^{REF}_{2}{2}", Nbins_cent, Bin_cent);
  hist_c22ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_c22ref->GetYaxis()->SetTitle("c^{REF}_{2}{2}");
  TH1D *hist_c24ref =
      new TH1D("c24ref", "c^{REF}_{2}{4}", Nbins_cent, Bin_cent);
  hist_c24ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_c24ref->GetYaxis()->SetTitle("c^{REF}_{2}{4}");
  TH1D *hist_v22ref =
      new TH1D("v22ref", "v^{REF}_{2}{2}", Nbins_cent, Bin_cent);
  hist_v22ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_v22ref->GetYaxis()->SetTitle("v^{REF}_{2}{2}");
  TH1D *hist_v24ref =
      new TH1D("v24ref", "v^{REF}_{2}{4}", Nbins_cent, Bin_cent);
  hist_v24ref->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_v24ref->GetYaxis()->SetTitle("v^{REF}_{2}{4}");

  // Evaluation for Cn{n} and Vn{n}
  for (int i = 0; i < Nbins_cent; i++) {
    double corr2 = tp_Corr2Ref_cent->GetBinContent(i + 1);
    double corr2e = tp_Corr2Ref_cent->GetBinError(i + 1);
    double corr4 = tp_Corr4Ref_cent->GetBinContent(i + 1);
    double corr4e = tp_Corr4Ref_cent->GetBinError(i + 1);

    double c22 = corr2;
    double c22e = corr2e;
    double c24 = corr4 - 2 * TMath::Power(c22, 2);
    double c24e = corr4e * corr4e + 16 * corr2 * corr2 * corr2e * corr2e;
    c24e = c24e > 0 ? TMath::Sqrt(c24e) : 0;
    hist_c22ref->SetBinContent(i + 1, isnan(c22) || isinf(c22) ? 0 : c22);
    hist_c22ref->SetBinError(i + 1, isnan(c22e) || isinf(c22e) ? 0 : c22e);
    hist_c24ref->SetBinContent(i + 1, isnan(c24) || isinf(c24) ? 0 : c24);
    hist_c24ref->SetBinError(i + 1, isnan(c24e) || isinf(c24e) ? 0 : c24e);

    double v22 = c22 > 0 ? TMath::Sqrt(c22) : 0;
    double v22e = c22 > 0 ? 0.5 * c22e / TMath::Sqrt(c22) : 0;
    double v24 = c24 < 0 ? TMath::Power(-c24, 1. / 4.) : 0;
    double v24e = c24 < 0 ? TMath::Power(-c24, (-3. / 4)) * c24e / 4 : 0;
    hist_v22ref->SetBinContent(i + 1, isnan(v22) || isinf(v22) ? 0 : v22);
    hist_v22ref->SetBinError(i + 1, isnan(v22e) || isinf(v22e) ? 0 : v22e);
    hist_v24ref->SetBinContent(i + 1, isnan(v24) || isinf(v24) ? 0 : v24);
    hist_v24ref->SetBinError(i + 1, isnan(v24e) || isinf(v24e) ? 0 : v24e);
  }

  l_refflow->Add(hist_c22ref);
  l_refflow->Add(hist_c24ref);
  l_refflow->Add(hist_v22ref);
  l_refflow->Add(hist_v24ref);
  l_refflow->Write("ReferenceFlow", TObject::kSingleKey);

  /////////////////////////////////////////////////
  ///                                           ///
  ///   Analysis for Differential Flow of POI   ///
  ///                                           ///
  /////////////////////////////////////////////////

  // Define variables' range for analysis
  double Bin_pt_mass[9] = {1., 2., 3., 4., 5., 6., 9., 12., 20.};
  double cent_min = 10.0;
  double cent_max = 60.0;
  double mass_min = 2.5; // same as initial setup
  double mass_max = 4.5; // same as initial setup
  int iBin_cent_min = centAxis->FindBin(cent_min);
  int iBin_cent_max = centAxis->FindBin(cent_max);
  int iBin_mass_min = massAxis->FindBin(mass_min);
  int iBin_mass_max = massAxis->FindBin(mass_max);

  for (int i = 0; i < 8; i++) {
    TList *l_diff = new TList();

    int iBin_pt_min = ptAxis->FindBin(Bin_pt_mass[i]);
    int iBin_pt_max = ptAxis->FindBin(Bin_pt_mass[i + 1]);
    // Copy original profiles for projections
    TProfile3D *tp_Corr2Ref_cp = dynamic_cast<TProfile3D *>(
        tp_Corr2Ref->Clone(Form("Mass_Pt_centrFT0C_Corr2REF_Copy_%d_%d",
                                int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1]))));
    TProfile3D *tp_Corr4Ref_cp = dynamic_cast<TProfile3D *>(
        tp_Corr4Ref->Clone(Form("Mass_Pt_centrFT0C_Corr4REF_Copy_%d_%d",
                                int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1]))));
    TProfile3D *tp_Corr2Poi_cp = dynamic_cast<TProfile3D *>(
        tp_Corr2Poi->Clone(Form("Mass_Pt_centrFT0C_Corr2POI_Copy_%d_%d",
                                int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1]))));
    TProfile3D *tp_Corr4Poi_cp = dynamic_cast<TProfile3D *>(
        tp_Corr4Poi->Clone(Form("Mass_Pt_centrFT0C_Corr4POI_Copy_%d_%d",
                                int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1]))));
    THnSparse *hist_mass_cp = dynamic_cast<THnSparse *>(
        hist_mass->Clone(Form("Mass_Pt_Rapidity_centrFT0C_Copy_%d_%d",
                              int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1]))));

    // Set axes' ranges for mass-differential study
    tp_Corr2Ref_cp->GetYaxis()->SetRange(iBin_pt_min, iBin_pt_max);
    tp_Corr2Ref_cp->GetZaxis()->SetRange(iBin_cent_min, iBin_cent_max);

    tp_Corr4Ref_cp->GetYaxis()->SetRange(iBin_pt_min, iBin_pt_max);
    tp_Corr4Ref_cp->GetZaxis()->SetRange(iBin_cent_min, iBin_cent_max);

    tp_Corr2Poi_cp->GetYaxis()->SetRange(iBin_pt_min, iBin_pt_max);
    tp_Corr2Poi_cp->GetZaxis()->SetRange(iBin_cent_min, iBin_cent_max);

    tp_Corr4Poi_cp->GetYaxis()->SetRange(iBin_pt_min, iBin_pt_max);
    tp_Corr4Poi_cp->GetZaxis()->SetRange(iBin_cent_min, iBin_cent_max);

    hist_mass_cp->GetAxis(0)->SetRange(iBin_mass_min, iBin_mass_max);
    hist_mass_cp->GetAxis(1)->SetRange(iBin_pt_min, iBin_pt_max);
    hist_mass_cp->GetAxis(3)->SetRange(iBin_cent_min, iBin_cent_max);

    // Define projected profiles w.r.t above defined variables' ranges
    TProfile *tp_Corr2Ref_mass =
        tp_Corr2Ref_cp->Project3DProfile("yx")->ProfileX();
    TProfile *tp_Corr4Ref_mass =
        tp_Corr4Ref_cp->Project3DProfile("yx")->ProfileX();
    TProfile *tp_Corr2Poi_mass =
        tp_Corr2Poi_cp->Project3DProfile("yx")->ProfileX();
    TProfile *tp_Corr4Poi_mass =
        tp_Corr4Poi_cp->Project3DProfile("yx")->ProfileX();
    TH1D *hist_mass_proj = hist_mass_cp->Projection(0);

    // Define histograms
    double *Bin_mass_new = CreateBinsFromAxis(hist_mass_proj->GetXaxis());
    int NBins_mass_new = hist_mass_proj->GetXaxis()->GetNbins();
    TH1D *hist_d22 = new TH1D(
        Form("d22_%d_%d", int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1])),
        Form("d^{J/#psi}_{2}{2}_%d_%d", int(Bin_pt_mass[i]),
             int(Bin_pt_mass[i + 1])),
        NBins_mass_new, Bin_mass_new);
    hist_d22->GetXaxis()->SetTitle("mass (GeV/c2)");
    hist_d22->GetYaxis()->SetTitle("d^{J/#psi}_{2}{2}");

    TH1D *hist_d24 = new TH1D(
        Form("d24_%d_%d", int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1])),
        Form("d^{J/#psi}_{2}{4}_%d_%d", int(Bin_pt_mass[i]),
             int(Bin_pt_mass[i + 1])),
        NBins_mass_new, Bin_mass_new);
    hist_d24->GetXaxis()->SetTitle("mass (GeV/c2)");
    hist_d24->GetYaxis()->SetTitle("d^{J/#psi}_{2}{4}");

    TH1D *hist_vd22 = new TH1D(
        Form("vd22_%d_%d", int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1])),
        Form("v'^{J/#psi}_{2}{2}_%d_%d", int(Bin_pt_mass[i]),
             int(Bin_pt_mass[i + 1])),
        NBins_mass_new, Bin_mass_new);
    hist_vd22->GetXaxis()->SetTitle("mass (GeV/c2)");
    hist_vd22->GetYaxis()->SetTitle("v'^{J/#psi}_{2}{2}");

    TH1D *hist_vd24 = new TH1D(
        Form("vd24_%d_%d", int(Bin_pt_mass[i]), int(Bin_pt_mass[i + 1])),
        Form("v'^{J/#psi}_{2}{4}_%d_%d", int(Bin_pt_mass[i]),
             int(Bin_pt_mass[i + 1])),
        NBins_mass_new, Bin_mass_new);
    hist_vd24->GetXaxis()->SetTitle("mass (GeV/c2)");
    hist_vd24->GetYaxis()->SetTitle("v'^{J/#psi}_{2}{4}");

    // Evaluation of differential flow as function of invariant mass
    for (int j = 0; j < NBins_mass_new; j++) {
      double corr2_ref = tp_Corr2Ref_mass->GetBinContent(j + 1);
      double corr2e_ref = tp_Corr2Ref_mass->GetBinError(j + 1);
      double corr4_ref = tp_Corr4Ref_mass->GetBinContent(j + 1);
      double corr4e_ref = tp_Corr4Ref_mass->GetBinError(j + 1);
      double corr2_poi = tp_Corr2Poi_mass->GetBinContent(j + 1);
      double corr2e_poi = tp_Corr2Poi_mass->GetBinError(j + 1);
      double corr4_poi = tp_Corr4Poi_mass->GetBinContent(j + 1);
      double corr4e_poi = tp_Corr4Poi_mass->GetBinError(j + 1);

      double c22 = corr2_ref;
      double c22e = corr2e_ref;
      double c24 = corr4_ref - 2 * TMath::Power(c22, 2);
      double c24e = corr4e_ref * corr4e_ref +
                    16 * corr2_ref * corr2_ref * corr2e_ref * corr2e_ref;
      c24e = c24e > 0 ? TMath::Sqrt(c24e) : 0;

      double d22 = corr2_poi;
      double d22e = corr2e_poi;
      double d24 = corr4_poi - 2 * corr2_ref * corr2_poi;
      double d24e = corr4e_poi * corr4e_ref +
                    4 * corr2_ref * corr2_ref * corr2e_poi * corr2e_poi +
                    4 * corr2_poi * corr2_poi * corr2e_ref * corr2e_ref;
      d24e = d24e > 0 ? TMath::Sqrt(d24e) : 0;
      hist_d22->SetBinContent(j + 1, isnan(d22) || isinf(d22) ? 0 : d22);
      hist_d22->SetBinError(j + 1, isnan(d22e) || isinf(d22e) ? 0 : d22e);
      hist_d24->SetBinContent(j + 1, isnan(d24) || isinf(d24) ? 0 : d24);
      hist_d24->SetBinError(j + 1, isnan(d24e) || isinf(d24e) ? 0 : d24e);

      double vd22 = c22 > 0 ? d22 / TMath::Sqrt(c22) : 0;
      double vd22e = corr2e_poi * corr2e_poi / corr2_ref +
                     0.25 * corr2_poi * corr2_poi * corr2e_ref * corr2e_ref /
                         (corr2_ref * corr2_ref * corr2_ref);
      vd22e = vd22e > 0 ? TMath::Sqrt(vd22e) : 0;
      double vd24 = c24 < 0 ? -d24 / TMath::Power(-c24, (3. / 4)) : 0;
      double vd24e =
          c24 < 0 ? TMath::Sqrt(TMath::Power(-c24, -6. / 4) * d24e * d24e +
                                TMath::Power(-c24, -14. / 4) * d24 * d24 *
                                    c24e * c24e * 9. / 16)
                  : 0;
      hist_vd22->SetBinContent(j + 1, isnan(vd22) || isinf(vd22) ? 0 : vd22);
      hist_vd22->SetBinError(j + 1, isnan(vd22e) || isinf(vd22e) ? 0 : vd22e);
      hist_vd24->SetBinContent(j + 1, isnan(vd24) || isinf(vd24) ? 0 : vd24);
      hist_vd24->SetBinError(j + 1, isnan(vd24e) || isinf(vd24e) ? 0 : vd24e);
    }

    // Do fitting
    runFitting(hist_mass_proj, l_diff, Bin_pt_mass[i], Bin_pt_mass[i + 1]);

    // Save results
    l_diff->Add(hist_mass_proj);
    l_diff->Add(hist_d22);
    l_diff->Add(hist_d24);
    l_diff->Add(hist_vd22);
    l_diff->Add(hist_vd24);
    l_diff->Write(Form("DifferentialFlow_%d_%d", int(Bin_pt_mass[i]),
                       int(Bin_pt_mass[i + 1])),
                  TObject::kSingleKey);
  }

  f.Close();
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

void runFitting(TH1D *hs, TList *ls, double ptmin, double ptmax) {

  // Setting up RooFit
  cout << "Start processing Pt range: " << ptmin << " " << ptmax << endl;
  RooRealVar m("m", "m_{#mu#mu}", 2.6, 4.5, "GeV/c2");
  RooDataHist dh("dh", "dh", m, Import(*hs));
  RooPlot *frame =
      m.frame(Name("mFrame"), Title(Form("m_{#mu#mu} signal fit (%d < pT < %d)",
                                         int(ptmin), int(ptmax))));
  dh.plotOn(frame, DataError(RooAbsData::Poisson), Name("data"));

  // Setting up fit model
  /// Signal model: double-sided crytalball
  RooRealVar m0("m0", "m0", 3.097, 3.09, 3.1);
  RooRealVar sigma("sigma", "sigma", 0.06, 0.001, 0.12);
  RooRealVar alphaL("alphaL", "alphaL", 5.0, 3.0, 10.0);
  RooRealVar alphaR("alphaR", "alphaR", 5.0, 3.0, 10.0);
  RooRealVar nL("nL", "nL", 0.05, 0.0, 1.0);
  RooRealVar nR("nR", "nR", 0.05, 0.0, 1.0);
  RooCrystalBall sig("Signal", "CB2", m, m0, sigma, alphaL, nL, alphaR, nR);
  /// Background model: VWG
  RooRealVar A("A", "A", -0.1, -2.0, 2.0);
  RooRealVar B("B", "B", 0.1, -5.5, 5.5);
  RooRealVar C("C", "C", -0.1, -1.5, 1.5);
  RooFormulaVar sigma_VWG("sigma_VWG", "sigma_VWG", "B+C*((m-A)/A)",
                          RooArgList(m, A, B, C));

  RooGenericPdf bkg("Bkg", "VWG", "exp(-(m-A)*(m-A)/(2.*sigma_VWG*sigma_VWG))",
                    RooArgSet(m, A, sigma_VWG));
  /// Construct a combined signal+background model
  RooRealVar nsig("nsig", "#signal events", 500, 0., 10000000);
  RooRealVar nbkg("nbkg", "#background events", 10000, 0., 10000000);
  RooAddPdf model("model", "CB2+VWG", {sig, bkg}, {nsig, nbkg});

  // Do fitting
  // RooFitResult *fit_result = model.fitTo(dh, Save(kTRUE));
  model.fitTo(dh, Minimizer("Minuit2", "migrad"), Strategy(0), Hesse(kFALSE));
  model.plotOn(frame, LineColor(kBlue), Name("model"));
  model.plotOn(frame, Components(sig), LineColor(kRed), Name("Signal"));
  model.plotOn(frame, Components(bkg), LineStyle(ELineStyle::kDashed),
               LineColor(kBlue), Name("Bkg"));
  double chi2 = frame->chiSquare("model", "data");
  TString chi2text = TString::Format("#chi^{2}/ndf = %f ", chi2);
  model.paramOn(frame, Label(chi2text), Layout(0.63, 0.9, 0.9),
                Format("NE", FixedPrecision(5)));
  // Saving plot
  TCanvas *c = new TCanvas(Form("Fitted_%s", hs->GetName()));
  c->cd();
  frame->Draw();
  ls->Add(c);
}