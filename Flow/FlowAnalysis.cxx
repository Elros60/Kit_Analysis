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
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <iostream>
#include <string>

#include "Framework/Logger.h"

using namespace std;
void LoadData(TChain *fChain, std::string TreeName = "O2rerefflow",
              std::string FileName = "AO2D.root");

void CreateBins(float *axis, float min, float max, int Nbins = 10);

void FlowAnalysis(float dimuonMassMin = 0., float dimuonMassMax = 4.,
                  float dimuonPtMin = 0., float dimuonPtMax = 10.,
                  bool fPropError = true) {

  TChain *fChain_REF = new TChain();
  TChain *fChain_POI = new TChain();
  LoadData(fChain_REF, "O2rerefflow");
  LoadData(fChain_POI, "O2rtdimuonall");

  int NBinsMult = 10;
  int NBinsMass = 20;
  int NBinsPt = 20;
  float Cent[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  float *MassBins = new float[NBinsMass + 1];
  float *PtBins = new float[NBinsPt + 1];

  CreateBins(MassBins, dimuonMassMin, dimuonMassMax, NBinsMass);
  CreateBins(PtBins, dimuonPtMin, dimuonPtMax, NBinsPt);

  //////////////////////////////////////////////////////////////////////////////////////
  // Post-processing for reference flow
  float C22MultREF[10] = {0.};
  float C22EMultREF[10] = {0.};
  float V22MultREF[10] = {0.};
  float V22EMultREF[10] = {0.};
  float C24MultREF[10] = {0.};
  float C24EMultREF[10] = {0.};
  float V24MultREF[10] = {0.};
  float V24EMultREF[10] = {0.};
  float Sum22REF[10] = {0.};
  float Sum24REF[10] = {0.};
  float C22REF = 0.; //! Event averaged 2-p cumulant for reference flow
  float C24REF = 0.; //! Event averaged 4-p cumulant for reference flow
  float C22EREF =
      0.; //! Error for event averaged 2-p cumulant for reference flow
  float C24EREF =
      0.; //! Error for event averaged 4-p cumulant for reference flow
  float Sum22REFAll = 0.;
  float Sum24REFAll = 0.;
  float M11Ref;
  float M1111Ref;
  float Corr2Ref;
  float Corr4Ref;
  float CentFT0REF;

  fChain_REF->SetBranchAddress("fM11REF", &M11Ref);
  fChain_REF->SetBranchAddress("fM1111REF", &M1111Ref);
  fChain_REF->SetBranchAddress("fCORR2REF", &Corr2Ref);
  fChain_REF->SetBranchAddress("fCORR4REF", &Corr4Ref);
  fChain_REF->SetBranchAddress("fCentFT0C", &CentFT0REF);

  TH1F *hist_c22REF = new TH1F("c22REF", "c22REF", NBinsMult, Cent);
  TH1F *hist_v22REF = new TH1F("v22REF", "v22REF", NBinsMult, Cent);
  TH1F *hist_c24REF = new TH1F("c24REF", "c24REF", NBinsMult, Cent);
  TH1F *hist_v24REF = new TH1F("v24REF", "v24REF", NBinsMult, Cent);

  // Loop over all event Q-vectors
  for (int i = 0; i < fChain_REF->GetEntries(); i++) {

    fChain_REF->GetEntry(i);
    if (isnan(Corr2Ref) || isnan(Corr4Ref) || isinf(Corr2Ref) ||
        isinf(Corr4Ref)) {
      continue;
    }

    C22REF += Corr2Ref * M11Ref;
    C24REF += Corr4Ref * M1111Ref;
    C22EREF += pow(Corr2Ref, 2) * M11Ref;
    C24EREF += pow(Corr4Ref, 2) * M1111Ref;
    Sum22REFAll += M11Ref;
    Sum24REFAll += M1111Ref;

    // Centrality-differential correlations
    for (int i = 0; i < NBinsMult; i++) {
      if (CentFT0REF > Cent[i] && CentFT0REF <= Cent[i + 1]) {
        C22MultREF[i] += Corr2Ref * M11Ref;
        C24MultREF[i] += Corr4Ref * M1111Ref;
        C22EMultREF[i] += pow(Corr2Ref, 2) * M11Ref;
        C24EMultREF[i] += pow(Corr4Ref, 2) * M1111Ref;
        Sum22REF[i] += M11Ref;
        Sum24REF[i] += M1111Ref;
      }
    }
  }

  // Events averaged reference flow
  C22REF = C22REF / Sum22REFAll;
  C24REF = C24REF / Sum24REFAll - 2. * pow(C22REF, 2);
  C22EREF = C22EREF / Sum22REFAll;
  C22EREF = pow(C22EREF - C22REF * C22REF, 1. / 2);
  C24EREF = C24EREF / Sum24REFAll;
  C24EREF = pow(C24EREF - C24REF * C24REF, 1. / 2);
  C24EREF = pow(C24EREF * C24EREF + 16. * C22EREF * C22EREF * C22REF * C22REF,
                1. / 2);

  // Centrality-dependent flow for reference
  for (int i = 0; i < NBinsMult; i++) {

    if (Sum22REF[i] * Sum24REF[i] != 0.) {

      C22MultREF[i] = C22MultREF[i] / Sum22REF[i];
      C22EMultREF[i] = C22EMultREF[i] / Sum22REF[i];
      C22EMultREF[i] =
          pow(C22EMultREF[i] - C22MultREF[i] * C22MultREF[i], 1. / 2);

      V22MultREF[i] = pow(C22MultREF[i], 1. / 2);
      V22EMultREF[i] = 0.5 * C22EMultREF[i] * pow(C22MultREF[i], -1. / 2);

      C24MultREF[i] = C24MultREF[i] / Sum24REF[i] - 2. * pow(C22MultREF[i], 2);
      C24EMultREF[i] = C24EMultREF[i] / Sum24REF[i];
      C24EMultREF[i] =
          C24EMultREF[i] - C24MultREF[i] * C24MultREF[i] > 0.
              ? pow(C24EMultREF[i] - C24MultREF[i] * C24MultREF[i], 1. / 2)
              : 0.;
      C24EMultREF[i] = pow(C24EMultREF[i] * C24EMultREF[i] +
                               16. * C22MultREF[i] * C22MultREF[i] *
                                   C22EMultREF[i] * C22EMultREF[i],
                           1. / 2);

      V24MultREF[i] = pow(-1. * C24MultREF[i], 1. / 4);
      V24EMultREF[i] = C24MultREF[i] > 0 ? 0
                                         : 1. / 4 * C24EMultREF[i] *
                                               pow(-C24MultREF[i], -3. / 4);

      cout << V24MultREF[i] << endl;
      cout << V24EMultREF[i] << endl;
      cout << C24MultREF[i] << endl;
      cout << C24EMultREF[i] << endl;
      cout << V22MultREF[i] << endl;
      cout << V22EMultREF[i] << endl;
      cout << C22MultREF[i] << endl;
      cout << C22EMultREF[i] << endl;
      cout << endl;
      cout << "=================================" << endl;

      if (fPropError) {

        hist_c22REF->SetBinContent(i + 1, C22MultREF[i]);
        hist_c22REF->SetBinError(i + 1, C22EMultREF[i]);

        hist_v22REF->SetBinContent(i + 1, V22MultREF[i]);
        hist_v22REF->SetBinError(i + 1, V22EMultREF[i]);

        hist_c24REF->SetBinContent(i + 1, C24MultREF[i]);
        hist_c24REF->SetBinError(i + 1, C24EMultREF[i]);

        hist_v24REF->SetBinContent(i + 1, V24MultREF[i]);
        hist_v24REF->SetBinError(i + 1, V24EMultREF[i]);

      } else {
        hist_c22REF->SetBinContent(i + 1, C22MultREF[i]);
        hist_v22REF->SetBinContent(i + 1, V22MultREF[i]);
        hist_c24REF->SetBinContent(i + 1, C24MultREF[i]);
        hist_v24REF->SetBinContent(i + 1, V24MultREF[i]);
      }
    }
  }

  TCanvas *c1REF = new TCanvas("c22REF");
  TCanvas *c2REF = new TCanvas("v22REF");
  TCanvas *c3REF = new TCanvas("c24REF");
  TCanvas *c4REF = new TCanvas("v24REF");

  c1REF->cd();
  hist_c22REF->SetMarkerStyle(20);
  if (fPropError) {
    hist_c22REF->Draw("EP");
  } else {
    hist_c22REF->Draw("HIST P");
  }
  c2REF->cd();
  hist_v22REF->SetMarkerStyle(20);
  if (fPropError) {
    hist_v22REF->Draw("EP");
  } else {
    hist_v22REF->Draw("HIST P");
  }
  c3REF->cd();
  hist_c24REF->SetMarkerStyle(20);
  if (fPropError) {
    hist_c24REF->Draw("EP");
  } else {
    hist_c24REF->Draw("HIST P");
  }
  c4REF->cd();
  hist_v24REF->SetMarkerStyle(20);
  if (fPropError) {
    hist_v24REF->Draw("EP");
  } else {
    hist_v24REF->Draw("HIST P");
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Post-processing for flow of POI
  float *C22POI = new float[NBinsMass];
  float *C22EPOI = new float[NBinsMass];
  float *V22POI = new float[NBinsMass];
  float *V22EPOI = new float[NBinsMass];
  float *C24POI = new float[NBinsMass];
  float *C24EPOI = new float[NBinsMass];
  float *V24POI = new float[NBinsMass];
  float *V24EPOI = new float[NBinsMass];
  float *Sum22POI = new float[NBinsMass];
  float *Sum24POI = new float[NBinsMass];

  for (int i = 0; i < NBinsMass; i++) {
    C22POI[i] = 0.;
    C22EPOI[i] = 0.;
    C24POI[i] = 0.;
    C24EPOI[i] = 0.;
    V22POI[i] = 0.;
    V22EPOI[i] = 0.;
    V24POI[i] = 0.;
    V24EPOI[i] = 0.;
    Sum22POI[i] = 0.;
    Sum24POI[i] = 0.;
  }

  float M01Poi;
  float M0111Poi;
  float Corr2Poi;
  float Corr4Poi;
  float CentFT0POI;
  float fMass;
  float fPt, fPt1, fPt2;
  float fEta, fEta1, fEta2;
  int fMultDimuons;

  fChain_POI->SetBranchAddress("fM01POI", &M01Poi);
  fChain_POI->SetBranchAddress("fM0111POI", &M0111Poi);
  fChain_POI->SetBranchAddress("fCORR2POI", &Corr2Poi);
  fChain_POI->SetBranchAddress("fCORR4POI", &Corr4Poi);
  fChain_POI->SetBranchAddress("fCentFT0C", &CentFT0POI);
  fChain_POI->SetBranchAddress("fMass", &fMass);
  fChain_POI->SetBranchAddress("fPt", &fPt);
  fChain_POI->SetBranchAddress("fPt1", &fPt1);
  fChain_POI->SetBranchAddress("fPt2", &fPt2);
  fChain_POI->SetBranchAddress("fEta", &fEta);
  fChain_POI->SetBranchAddress("fEta1", &fEta1);
  fChain_POI->SetBranchAddress("fEta2", &fEta2);
  fChain_POI->SetBranchAddress("fMultDimuons", &fMultDimuons);

  TH1F *hist_c22POIMass =
      new TH1F("c22POIMass", "c22POIMass", NBinsMass, MassBins);
  TH1F *hist_v22POIMass =
      new TH1F("v22POIMass", "v22POIMass", NBinsMass, MassBins);
  TH1F *hist_c24POIMass =
      new TH1F("c24POIMass", "c24POIMass", NBinsMass, MassBins);
  TH1F *hist_v24POIMass =
      new TH1F("v24POIMass", "v24POIMass", NBinsMass, MassBins);

  // Loop over all dimuons
  for (int i = 0; i < fChain_POI->GetEntries(); i++) {

    fChain_POI->GetEntry(i);
    if (isnan(Corr2Poi) || isnan(Corr4Poi) || isinf(Corr2Poi) ||
        isinf(Corr4Poi)) {
      continue;
    }

    // Dimuons selection
    if (!(fPt > dimuonPtMin && fPt <= dimuonPtMax && fMass > dimuonMassMin &&
          fMass <= dimuonMassMax && fEta > -4. && fEta < -2.5 &&
          CentFT0POI > 0. && CentFT0POI <= 70.)) {
      continue;
    }

    // Mass-differential
    for (int i = 0; i < NBinsMass; i++) {

      if (fMass > MassBins[i] && fMass <= MassBins[i + 1]) {
        C22POI[i] += Corr2Poi * M01Poi;
        C24POI[i] += Corr4Poi * M0111Poi;
        C22EPOI[i] += pow(Corr2Poi, 2) * M01Poi;
        C24EPOI[i] += pow(Corr4Poi, 2) * M0111Poi;
        Sum22POI[i] += M01Poi / fMultDimuons;
        Sum24POI[i] += M0111Poi / fMultDimuons;
      }
    }
  }

  // Mass-dependent flow for POI
  for (int i = 0; i < NBinsMass; i++) {

    if (Sum22POI[i] * Sum24POI[i] != 0.) {

      C22POI[i] = C22POI[i] / Sum22POI[i];
      C22EPOI[i] = C22EPOI[i] / Sum22POI[i];
      C22EPOI[i] = pow(C22EPOI[i] - C22POI[i] * C22POI[i], 1. / 2);

      V22POI[i] = C22POI[i] * pow(C22REF, -1. / 2);
      V22EPOI[i] = pow(C22EPOI[i] * C22EPOI[i] / C22REF +
                           0.25 * C22POI[i] * C22POI[i] * C22EREF * C22EREF /
                               pow(C22REF, 3),
                       1. / 2);
      V22EPOI[i] = V22EPOI[i] > 0 ? V22EPOI[i] : 0.;

      C24POI[i] = C24POI[i] / Sum24POI[i] - 2. * C22POI[i] * C22REF;
      C24EPOI[i] = C24EPOI[i] / Sum24POI[i];
      C24EPOI[i] = pow(C24EPOI[i] - C24POI[i] * C24POI[i], 1. / 2);
      C24EPOI[i] = pow(C24EPOI[i] * C24EPOI[i] +
                           4. * C22REF * C22REF * C22EPOI[i] * C22EPOI[i] +
                           4. * C22POI[i] * C22POI[i] * C22EREF * C22EREF,
                       1. / 2);

      V24POI[i] = -1. * C24POI[i] * pow(-1. * C24REF, -3. / 4);
      V24EPOI[i] = pow(pow(-1. * C24REF, -6. / 4) * C24EPOI[i] * C24EPOI[i] +
                           pow(-1. * C24REF, -14. / 4) * C24POI[i] * C24POI[i] *
                               C24EREF * C24EREF * 9. / 16,
                       1. / 2);

      cout << V24POI[i] << endl;
      cout << V24EPOI[i] << endl;
      cout << C24POI[i] << endl;
      cout << C24EPOI[i] << endl;
      cout << V22POI[i] << endl;
      cout << V22EPOI[i] << endl;
      cout << C22POI[i] << endl;
      cout << C22EPOI[i] << endl;
      cout << endl;

      if (fPropError) {
        if (isnan(C22EPOI[i])) {
          hist_c22POIMass->SetBinContent(i + 1, C22POI[i]);
          hist_c22POIMass->SetBinError(i + 1, 0.);
        } else {
          hist_c22POIMass->SetBinContent(i + 1, C22POI[i]);
          hist_c22POIMass->SetBinError(i + 1, C22EPOI[i]);
        }

        if (isnan(V22EPOI[i])) {
          hist_v22POIMass->SetBinContent(i + 1, V22POI[i]);
          hist_v22POIMass->SetBinError(i + 1, 0.);
        } else {
          hist_v22POIMass->SetBinContent(i + 1, V22POI[i]);
          hist_v22POIMass->SetBinError(i + 1, V22EPOI[i]);
        }

        if (isnan(C24EPOI[i])) {
          hist_c24POIMass->SetBinContent(i + 1, C24POI[i]);
          hist_c24POIMass->SetBinError(i + 1, 0.);
        } else {
          hist_c24POIMass->SetBinContent(i + 1, C24POI[i]);
          hist_c24POIMass->SetBinError(i + 1, C24EPOI[i]);
        }

        if (isnan(V24EPOI[i])) {
          hist_v24POIMass->SetBinContent(i + 1, V24POI[i]);
          hist_v24POIMass->SetBinError(i, 0.);
        } else {
          hist_v24POIMass->SetBinContent(i + 1, V24POI[i]);
          hist_v24POIMass->SetBinError(i + 1, V24EPOI[i]);
        }
      } else {
        hist_c22POIMass->SetBinContent(i + 1, C22POI[i]);
        hist_v22POIMass->SetBinContent(i + 1, V22POI[i]);
        hist_c24POIMass->SetBinContent(i + 1, C24POI[i]);
        hist_v24POIMass->SetBinContent(i + 1, V24POI[i]);
      }
    }
  }

  TCanvas *c1POI = new TCanvas("c22POI");
  TCanvas *c2POI = new TCanvas("v22POI");
  TCanvas *c3POI = new TCanvas("c24POI");
  TCanvas *c4POI = new TCanvas("v24POI");

  c1POI->cd();
  hist_c22POIMass->SetMarkerStyle(20);
  if (fPropError) {
    hist_c22POIMass->Draw("EP");
  } else {
    hist_c22POIMass->Draw("HIST P");
  }
  c2POI->cd();
  hist_v22POIMass->SetMarkerStyle(20);
  if (fPropError) {
    hist_v22POIMass->Draw("EP");
  } else {
    hist_v22POIMass->Draw("HIST P");
  }
  c3POI->cd();
  hist_c24POIMass->SetMarkerStyle(20);
  if (fPropError) {
    hist_c24POIMass->Draw("EP");
  } else {
    hist_c24POIMass->Draw("HIST P");
  }
  c4POI->cd();
  hist_v24POIMass->SetMarkerStyle(20);
  if (fPropError) {
    hist_v24POIMass->Draw("EP");
  } else {
    hist_v24POIMass->Draw("HIST P");
  }
}

void LoadData(TChain *fChain, std::string TreeName, std::string FileName) {

  TFile *fInput = TFile::Open(FileName.c_str());
  TIter keyList(fInput->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)keyList())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    std::string dir = key->GetName();
    if (dir == "parentFiles") {
      continue;
    }
    fChain->Add(
        Form("%s?#%s/%s", FileName.c_str(), dir.c_str(), TreeName.c_str()));
  }
}

void CreateBins(float *axis, float min, float max, int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}