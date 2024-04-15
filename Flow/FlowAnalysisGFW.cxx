#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <THashList.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TString.h>
#include <iostream>
#include <vector>

#include "Framework/Logger.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

using namespace std;
using namespace o2;
TH1D *Prof2Hist(TProfile *inpf);

void FlowAnalysisGFW(std::string FileName = "AnalysisResults.root") {

  // Load object FlowContainer
  TFile *file;
  try {
    file = TFile::Open(FileName.c_str());
  } catch (std::exception const &e) {
    LOG(fatal) << "Cannot open input file!";
  }

  FlowContainer *fFC;
  try {
    // file->GetObject(Form("d-q-event-qvector/%s", "FlowContainer"), fFC);
    file->GetObject(Form("d-q-event-qvector/%s", "FlowContainer"), fFC);
  } catch (std::exception const &e) {
    LOG(fatal) << "Cannot get GFWWeights from inpout file!";
  }

  TCanvas *c1 = new TCanvas("2P plots");
  c1->Divide(2);
  TCanvas *c2 = new TCanvas("4P plots");
  c2->Divide(2);

  // Get 2D profile in FlowContainer
  TProfile2D *tp_2D = fFC->GetProfile();

  // Get 2,4-correlation: c_2{2,4}
  int id_corrfull22 = tp_2D->GetYaxis()->FindBin("ChFull22");
  int id_corrfull24 = tp_2D->GetYaxis()->FindBin("ChFull24");
  TProfile *tp_corrfull22 =
      tp_2D->ProfileX("CorrFull22", id_corrfull22, id_corrfull22);
  TProfile *tp_corrfull24 =
      tp_2D->ProfileX("CorrFull24", id_corrfull24, id_corrfull24);

  TH1D *hist_corrfull22 = Prof2Hist(tp_corrfull22); // <<2>>
  TH1D *hist_v22 = dynamic_cast<TH1D *>(hist_corrfull22->Clone("v2_2"));
  for (int i = 1; i <= hist_v22->GetNbinsX(); i++) {
    double d2 = hist_corrfull22->GetBinContent(i);
    double d2e = hist_corrfull22->GetBinError(i);
    if (d2 > 0) {
      hist_v22->SetBinContent(i, d2 > 0 ? TMath::Sqrt(d2) : -2);
      hist_v22->SetBinError(i, d2 > 0 ? 0.5 * d2e / TMath::Sqrt(d2) : 0);
    }
  }

  c1->cd(1);
  tp_corrfull22->Draw();

  c1->cd(2);
  hist_v22->Draw();

  TH1D *hist_corrfull24 = Prof2Hist(tp_corrfull24); // <<4>>
  TH1D *hist_c24 = dynamic_cast<TH1D *>(hist_corrfull24->Clone("c2_4"));
  for (int i = 1; i <= hist_c24->GetNbinsX(); i++) {
    double cor4v = hist_corrfull24->GetBinContent(i);
    double cor4e = hist_corrfull24->GetBinError(i);
    double cor2v = hist_corrfull22->GetBinContent(i);
    double cor2e = hist_corrfull22->GetBinError(i);
    hist_c24->SetBinContent(i, cor4v - 2 * TMath::Power(cor2v, 2));
    hist_c24->SetBinError(
        i, TMath::Sqrt(cor4e * cor4e + 16 * cor2v * cor2v * cor2e * cor2e));
  }

  hist_c24->SetBinContent(0, 0);
  hist_c24->SetBinError(0, 0);

  c2->cd(1);
  hist_c24->Draw();

  TH1D *hist_v24 = dynamic_cast<TH1D *>(hist_c24->Clone("v2_4"));
  for (int i = 1; i <= hist_v24->GetNbinsX(); i++) {
    double d4 = hist_c24->GetBinContent(i);
    double d4e = hist_c24->GetBinError(i);
    hist_v24->SetBinContent(i, d4 >= 0 ? -2 : TMath::Power(-d4, 1. / 4));
    hist_v24->SetBinError(i, d4 >= 0 ? 0.
                                     : TMath::Power(-d4, (-3. / 4)) * d4e / 4);
  }

  c2->cd(2);
  hist_v24->Draw();
}

TH1D *Prof2Hist(TProfile *inpf) {
  int nbins = inpf->GetNbinsX();
  double *xbs = new double[nbins + 1];
  inpf->GetLowEdge(xbs);
  xbs[nbins] = xbs[nbins - 1] + inpf->GetBinWidth(nbins);
  TH1D *rethist =
      new TH1D(Form("%s_hist", inpf->GetName()), inpf->GetTitle(), nbins, xbs);
  for (int i = 1; i <= rethist->GetNbinsX(); i++) {
    if (inpf->GetBinContent(i) != 0) {
      rethist->SetBinContent(i, inpf->GetBinContent(i));
      rethist->SetBinError(i, inpf->GetBinError(i));
    }
  }
  return rethist;
}
