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
#include <THStack.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
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
#include <TTreeIndex.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>

#include "Framework/Logger.h"

using namespace std;

void LoadData(TChain *fChain, std::string TreeName = "O2rerefflow",
              std::string FileName = "AO2D.root");
double getEventPlane(int n, double x, double y);
void CreateBins(double *axis, double min, double max, int Nbins = 10);

// Resolution factor(Reference flow) evaluation with all events included
void Resolution(int NbinCent = 20) {
  TChain *fChain_Qvec = new TChain("Qvectors");
  TChain *fChain_Event = new TChain("EventExtended");
  LoadData(fChain_Qvec, "O2reqvectorctr");
  LoadData(fChain_Event, "O2reextended");

  fChain_Qvec->AddFriend(fChain_Event);
  std::string ComboLabel[9] = {"TCPPos_TPCNeg", "FT0A_TPCPos",  "FT0A_TPCNeg",
                               "FT0C_TPCPos",   "FT0C_TPCNeg",  "FT0C_FT0A",
                               "FT0A_TPCFull",  "FT0C_TPCFull", "3-Sub"};
  double Combo[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  double *Cent = new double[NbinCent + 1];
  CreateBins(Cent, 0., 90., NbinCent);
  TProfile2D *ResoProfileSPRe = new TProfile2D(
      "ResoProfileSPRe", "Profile of resolution factors for SP (real)", 8,
      Combo, NbinCent, Cent);
  TProfile2D *ResoProfileSPIm = new TProfile2D(
      "ResoProfileSPIm", "Profile of resolution factors for SP (imaginary)", 8,
      Combo, NbinCent, Cent);
  TProfile2D *ResoProfileEP =
      new TProfile2D("ResoProfileEP", "Profile of resolution factors for EP", 8,
                     Combo, NbinCent, Cent);

  TH1D *hist_RSPRe[9];
  TH1D *hist_RSPIm[9];
  TH1D *hist_REP[9];

  // Compare with current formulae
  TProfile2D *ResoProfileSPRe_old = new TProfile2D(
      "ResoProfileSPRe_old", "Profile of resolution factors for SP (real)", 8,
      Combo, NbinCent, Cent);

  TH1D *hist_RSPRe_old[9];
  TH1D *hist_RSPRe_diff[9];

  for (int i = 0; i < 9; i++) {
    hist_RSPRe[i] = new TH1D(Form("RSPRe_%s", ComboLabel[i].c_str()),
                             ComboLabel[i].c_str(), NbinCent, Cent);
    hist_RSPIm[i] = new TH1D(Form("RSPIm_%s", ComboLabel[i].c_str()),
                             ComboLabel[i].c_str(), NbinCent, Cent);
    hist_REP[i] = new TH1D(Form("REP_%s", ComboLabel[i].c_str()),
                           ComboLabel[i].c_str(), NbinCent, Cent);

    // Compare with current formulae
    hist_RSPRe_old[i] = new TH1D(Form("RSPRe_old_%s", ComboLabel[i].c_str()),
                                 ComboLabel[i].c_str(), NbinCent, Cent);
    hist_RSPRe_diff[i] = new TH1D(Form("RSPRe_diff_%s", ComboLabel[i].c_str()),
                                  ComboLabel[i].c_str(), NbinCent, Cent);
  }

  float qvecFT0CRe, qvecFT0CIm, qvecFT0ARe, qvecFT0AIm, qvecFT0MRe, qvecFT0MIm,
      qvecFV0ARe, qvecFV0AIm, qvecBPosRe, qvecBPosIm, qvecBNegRe, qvecBNegIm,
      qvecBTotRe, qvecBTotIm, centFT0C;
  int nTrkBPos, nTrkBNeg;

  fChain_Qvec->SetBranchAddress("fQvecFT0ARe", &qvecFT0ARe);
  fChain_Qvec->SetBranchAddress("fQvecFT0AIm", &qvecFT0AIm);
  fChain_Qvec->SetBranchAddress("fQvecFT0CRe", &qvecFT0CRe);
  fChain_Qvec->SetBranchAddress("fQvecFT0CIm", &qvecFT0CIm);
  fChain_Qvec->SetBranchAddress("fQvecFT0MRe", &qvecFT0MRe);
  fChain_Qvec->SetBranchAddress("fQvecFT0MIm", &qvecFT0MIm);
  fChain_Qvec->SetBranchAddress("fQvecFV0ARe", &qvecFV0ARe);
  fChain_Qvec->SetBranchAddress("fQvecFV0AIm", &qvecFV0AIm);
  fChain_Qvec->SetBranchAddress("fQvecBPosRe", &qvecBPosRe);
  fChain_Qvec->SetBranchAddress("fQvecBPosIm", &qvecBPosIm);
  fChain_Qvec->SetBranchAddress("fQvecBNegRe", &qvecBNegRe);
  fChain_Qvec->SetBranchAddress("fQvecBNegIm", &qvecBNegIm);
  fChain_Qvec->SetBranchAddress("fNTrkBPos", &nTrkBPos);
  fChain_Qvec->SetBranchAddress("fNTrkBNeg", &nTrkBNeg);
  // fChain_Qvec->SetBranchAddress("fQvecBTotRe", &qvecBTotRe);
  // fChain_Qvec->SetBranchAddress("fQvecBTotIm", &qvecBTotIm);
  fChain_Qvec->SetBranchAddress("EventExtended.fCentFT0C", &centFT0C);

  for (int i = 0; i < fChain_Qvec->GetEntries(); i++) {
    fChain_Qvec->GetEntry(i);
    for (int j = 0; j < NbinCent; j++) {
      if (centFT0C >= Cent[j] && centFT0C < Cent[j + 1]) {
        complex<double> Q1, Q2, Q12;
        double Psi1 = 0, Psi2 = 0, cosPsi12 = 0;
        for (int k = 0; k < 8; k++) {
          switch (k) {
          case 0:
            Q1 = complex<double>(qvecBPosRe, qvecBPosIm);
            Q2 = complex<double>(qvecBNegRe, qvecBNegIm);
            Psi1 = getEventPlane(2, qvecBPosRe, qvecBPosIm);
            Psi2 = getEventPlane(2, qvecBNegRe, qvecBNegIm);
            break;
          case 1:
            Q1 = complex<double>(qvecFT0ARe, qvecFT0AIm);
            Q2 = complex<double>(qvecBPosRe, qvecBPosIm);
            Psi1 = getEventPlane(2, qvecFT0ARe, qvecFT0AIm);
            Psi2 = getEventPlane(2, qvecBPosRe, qvecBPosIm);
            break;
          case 2:
            Q1 = complex<double>(qvecFT0ARe, qvecFT0AIm);
            Q2 = complex<double>(qvecBNegRe, qvecBNegIm);
            Psi1 = getEventPlane(2, qvecFT0ARe, qvecFT0AIm);
            Psi2 = getEventPlane(2, qvecBNegRe, qvecBNegIm);
            break;
          case 3:
            Q1 = complex<double>(qvecFT0CRe, qvecFT0CIm);
            Q2 = complex<double>(qvecBPosRe, qvecBPosIm);
            Psi1 = getEventPlane(2, qvecFT0CRe, qvecFT0CIm);
            Psi2 = getEventPlane(2, qvecBPosRe, qvecBPosIm);
            break;
          case 4:
            Q1 = complex<double>(qvecFT0CRe, qvecFT0CIm);
            Q2 = complex<double>(qvecBNegRe, qvecBNegIm);
            Psi1 = getEventPlane(2, qvecFT0CRe, qvecFT0CIm);
            Psi2 = getEventPlane(2, qvecBNegRe, qvecBNegIm);
            break;
          case 5:
            Q1 = complex<double>(qvecFT0CRe, qvecFT0CIm);
            Q2 = complex<double>(qvecFT0ARe, qvecFT0AIm);
            Psi1 = getEventPlane(2, qvecFT0CRe, qvecFT0CIm);
            Psi2 = getEventPlane(2, qvecFT0ARe, qvecFT0AIm);
            break;
          case 6:
            qvecBTotRe = (nTrkBNeg + nTrkBPos) == 0
                             ? 0
                             : (nTrkBPos * qvecBPosRe + nTrkBNeg * qvecBNegRe) /
                                   (nTrkBNeg + nTrkBPos);
            qvecBTotIm = (nTrkBNeg + nTrkBPos) == 0
                             ? 0
                             : (nTrkBPos * qvecBPosIm + nTrkBNeg * qvecBNegIm) /
                                   (nTrkBNeg + nTrkBPos);
            Q1 = complex<double>(qvecFT0ARe, qvecFT0AIm);
            Q2 = complex<double>(qvecBTotRe, qvecBTotIm);
            Psi1 = getEventPlane(2, qvecFT0ARe, qvecFT0AIm);
            Psi2 = getEventPlane(2, qvecBTotRe, qvecBTotIm);
            break;
          case 7:
            qvecBTotRe = (nTrkBNeg + nTrkBPos) == 0
                             ? 0
                             : (nTrkBPos * qvecBPosRe + nTrkBNeg * qvecBNegRe) /
                                   (nTrkBNeg + nTrkBPos);
            qvecBTotIm = (nTrkBNeg + nTrkBPos) == 0
                             ? 0
                             : (nTrkBPos * qvecBPosIm + nTrkBNeg * qvecBNegIm) /
                                   (nTrkBNeg + nTrkBPos);
            Q1 = complex<double>(qvecFT0CRe, qvecFT0CIm);
            Q2 = complex<double>(qvecBTotRe, qvecBTotIm);
            Psi1 = getEventPlane(2, qvecFT0CRe, qvecFT0CIm);
            Psi2 = getEventPlane(2, qvecBTotRe, qvecBTotIm);
            break;
          }

          Q12 = Q1 * conj(Q2);
          cosPsi12 = TMath::Cos(2 * (Psi1 - Psi2));
          ResoProfileEP->Fill(k, centFT0C, cosPsi12);
          ResoProfileSPRe->Fill(k, centFT0C, Q12.real());
          ResoProfileSPIm->Fill(k, centFT0C, Q12.imag());

          // Compare with current formulae
          ResoProfileSPRe_old->Fill(
              k, centFT0C, Q1.real() * Q2.real() + Q1.imag() * Q2.imag());
        }
      }
    }
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < NbinCent; j++) {
      if (i == 8) {
        double Q12Re = ResoProfileSPRe->GetBinContent(7, j + 1);
        double Q12Im = ResoProfileSPIm->GetBinContent(7, j + 1);
        double Q13Re = ResoProfileSPRe->GetBinContent(8, j + 1);
        double Q13Im = ResoProfileSPIm->GetBinContent(8, j + 1);
        double Q23Re = ResoProfileSPRe->GetBinContent(6, j + 1);
        double Q23Im = ResoProfileSPIm->GetBinContent(6, j + 1);
        complex<double> Q12(Q12Re, Q12Im);
        complex<double> Q13(Q13Re, Q13Im);
        complex<double> Q23(Q23Re, Q23Im);

        double Psi12 = ResoProfileEP->GetBinContent(7, j + 1);
        double Psi13 = ResoProfileEP->GetBinContent(8, j + 1);
        double Psi23 = ResoProfileEP->GetBinContent(6, j + 1);

        double REP = Psi23 != 0 ? TMath::Sqrt(Psi12 * Psi13 / Psi23) : 0;
        double RSPRe = abs(Q23) != 0 ? pow(Q12 * Q13 / Q23, 0.5).real() : 0;
        double RSPIm = abs(Q23) != 0 ? pow(Q12 * Q13 / Q23, 0.5).imag() : 0;
        hist_RSPRe[i]->SetBinContent(j + 1, RSPRe);
        hist_RSPIm[i]->SetBinContent(j + 1, RSPIm);
        hist_REP[i]->SetBinContent(j + 1, REP);

        // Compare with current formulae
        double Q12Re_old = ResoProfileSPRe_old->GetBinContent(7, j + 1);
        double Q13Re_old = ResoProfileSPRe_old->GetBinContent(8, j + 1);
        double Q23Re_old = ResoProfileSPRe_old->GetBinContent(6, j + 1);
        double RSPRe_old =
            Q23Re_old != 0 ? TMath::Sqrt(Q12Re_old * Q13Re_old / Q23Re_old) : 0;
        hist_RSPRe_old[i]->SetBinContent(j + 1, RSPRe_old);
      } else {
        double RSPRe = ResoProfileSPRe->GetBinContent(i + 1, j + 1);
        double RSPIm = ResoProfileSPIm->GetBinContent(i + 1, j + 1);
        complex<double> RSP(RSPRe, RSPIm);
        double REP =
            ResoProfileEP->GetBinContent(i + 1, j + 1) >= 0
                ? TMath::Sqrt(ResoProfileEP->GetBinContent(i + 1, j + 1))
                : 0;

        hist_RSPRe[i]->SetBinContent(j + 1, pow(RSP, 0.5).real());
        hist_RSPIm[i]->SetBinContent(j + 1, pow(RSP, 0.5).imag());
        hist_REP[i]->SetBinContent(j + 1, REP);

        // Compare with current formulae
        double RSPRe_old =
            ResoProfileSPRe_old->GetBinContent(i + 1, j + 1) >= 0
                ? TMath::Sqrt(ResoProfileSPRe_old->GetBinContent(i + 1, j + 1))
                : 0;
        hist_RSPRe_old[i]->SetBinContent(j + 1, RSPRe_old);
      }
    }
  }

  // Draw plots
  TCanvas *c = new TCanvas("ResoPlotsRe", "");
  TCanvas *c1 = new TCanvas("ResoPlotsIm", "");
  auto hs1 = new THStack("hs1", "");
  auto hs2 = new THStack("hs2", "");
  auto hs3 = new THStack("hs3", "");

  // Compare with current formulae
  TCanvas *c_old = new TCanvas("ResoPlotsRe_old", "");
  TCanvas *c_diff = new TCanvas("ResoPlotsRe_diff", "");
  auto hs1_old = new THStack("hs1_old", "");
  auto hs1_diff = new THStack("hs1_diff", "");

  for (int i = 1; i < 9; i++) {
    hist_RSPRe[i]->SetMarkerColor(i + 1);
    hist_RSPRe[i]->SetMarkerSize(1.5);
    hist_RSPRe[i]->SetLineWidth(2);
    hist_RSPRe[i]->SetLineColor(i + 1);
    hist_RSPRe[i]->SetMarkerStyle(i + 86);
    hs1->Add(hist_RSPRe[i]);

    hist_REP[i]->SetMarkerColor(i + 1);
    hist_REP[i]->SetMarkerSize(1.5);
    hist_REP[i]->SetLineWidth(2);
    hist_REP[i]->SetLineColor(i + 1);
    hist_REP[i]->SetMarkerStyle(i + 86);
    hs2->Add(hist_REP[i]);

    hist_RSPIm[i]->SetMarkerColor(i + 1);
    hist_RSPIm[i]->SetMarkerSize(1.5);
    hist_RSPIm[i]->SetLineWidth(2);
    hist_RSPIm[i]->SetLineColor(i + 1);
    hist_RSPIm[i]->SetMarkerStyle(i + 86);
    hs3->Add(hist_RSPIm[i]);

    // Compare with current formulae
    hist_RSPRe_old[i]->SetMarkerColor(i + 1);
    hist_RSPRe_old[i]->SetMarkerSize(1.5);
    hist_RSPRe_old[i]->SetLineWidth(2);
    hist_RSPRe_old[i]->SetLineColor(i + 1);
    hist_RSPRe_old[i]->SetMarkerStyle(i + 86);
    hs1_old->Add(hist_RSPRe_old[i]);

    hist_RSPRe_diff[i]->Add(hist_RSPRe[i], hist_RSPRe_old[i], 1., -1.);
    hs1_diff->Add(hist_RSPRe_diff[i]);
  }

  hist_REP[0]->SetMarkerColor(1);
  hist_REP[0]->SetMarkerSize(1.5);
  hist_REP[0]->SetLineWidth(2);
  hist_REP[0]->SetLineColor(1);
  hist_REP[0]->SetMarkerStyle(86);
  hs2->Add(hist_REP[0]);

  hist_RSPIm[0]->SetMarkerColor(1);
  hist_RSPIm[0]->SetMarkerSize(1.5);
  hist_RSPIm[0]->SetLineWidth(2);
  hist_RSPIm[0]->SetLineColor(1);
  hist_RSPIm[0]->SetMarkerStyle(86);
  hs3->Add(hist_RSPIm[0]);

  /*
  // Compare with current formulae
  hist_RSPRe_old[0]->SetMarkerColor(1);
  hist_RSPRe_old[0]->SetMarkerSize(1.5);
  hist_RSPRe_old[0]->SetLineWidth(2);
  hist_RSPRe_old[0]->SetLineColor(1);
  hist_RSPRe_old[0]->SetMarkerStyle(86);
  hs1_old->Add(hist_RSPRe_old[0]);
  */

  c->Divide(2, 1);
  c->cd(1);
  c->cd(1)->SetGrid();
  hs1->Draw("cp nostack");
  hs1->GetXaxis()->SetTitle("Centrality FT0C %");
  hs1->GetYaxis()->SetTitle("Re#{}{R_{2}(SP)}");
  c->cd(2);
  c->cd(2)->SetGrid();
  hs2->Draw("cp nostack");
  hs2->GetXaxis()->SetTitle("Centrality FT0C %");
  hs2->GetYaxis()->SetTitle("R_{2}(EP)");

  c1->cd();
  c1->cd()->SetGrid();
  hs3->Draw("cp nostack");
  hs3->GetXaxis()->SetTitle("Centrality FT0C %");
  hs3->GetYaxis()->SetTitle("Im#{}{R_{2}(SP)}");

  // Compare with current formulae
  c_old->Divide(2, 1);
  c_old->cd(1);
  c_old->cd(1)->SetGrid();
  hs1->Draw("cp nostack");
  hs1->GetXaxis()->SetTitle("Centrality FT0C %");
  hs1->GetYaxis()->SetTitle("Re#{}{R_{2}(SP)}");
  c_old->cd(2);
  c_old->cd(2)->SetGrid();
  hs1_old->Draw("cp nostack");
  hs1_old->GetXaxis()->SetTitle("Centrality FT0C %");
  hs1_old->GetYaxis()->SetTitle("Re#{}{R_{2}(SP)}");

  c_diff->cd();
  c_diff->cd()->SetGrid();
  hs1_diff->Draw("cp nostack");
  hs1_diff->GetXaxis()->SetTitle("Centrality FT0C %");
  hs1_diff->GetYaxis()->SetTitle("#Delta{R_{2}(SP)}");

  TFile f("histos.root", "RECREATE");
  TList *l1 = new TList();
  TList *l2 = new TList();
  for (int i = 0; i < 9; i++) {
    l1->Add(hist_RSPRe[i]);
    l1->Add(hist_RSPIm[i]);
    l1->Add(hist_REP[i]);

    // Compare with current formulae
    l2->Add(hist_RSPRe_old[i]);
  }

  l1->Add(c);
  l1->Add(c1);
  l1->Write("ResolutionHistList_New", TObject::kSingleKey);

  // Compare with current formulae
  l2->Add(c_old);
  l2->Add(c_diff);
  l2->Write("ResolutionHistList_Old", TObject::kSingleKey);

  f.ls();
  f.Close();
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

double getEventPlane(int n, double x, double y) {
  return (1.0 / n) * TMath::ATan2(y, x);
}

void CreateBins(double *axis, double min, double max, int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}