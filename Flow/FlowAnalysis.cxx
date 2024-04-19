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
#include <iostream>
#include <string>

#include "Framework/Logger.h"

using namespace std;

void LoadData(TChain *fChain, std::string TreeName = "O2rerefflow",
              std::string FileName = "AO2D.root");
void CreateBins(float *axis, float min, float max, int Nbins = 10);

void FlowAnalysis(std::array<float, 3> dimuonMassRange,
                  std::array<float, 3> dimuonPtRange,
                  std::array<float, 2> dimuonCentRange, bool fPropError = true,
                  bool fCumulant = true) {

  TChain *fChain_REF = new TChain();
  if (fCumulant) {
    LoadData(fChain_REF, "O2rerefflow");
  }
  TChain *fChain_POI = new TChain();
  LoadData(fChain_POI, "O2rtdimuonall");

  int NBinsMass = static_cast<int>(dimuonMassRange[2]);
  int NBinsPt = static_cast<int>(dimuonPtRange[2]);
  int NBinsMult = 10;
  float Cent[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  float *MassBins = new float[NBinsMass + 1];
  float *PtBins = new float[NBinsPt + 1];

  CreateBins(MassBins, dimuonMassRange[0], dimuonMassRange[1], NBinsMass);
  CreateBins(PtBins, dimuonPtRange[0], dimuonPtRange[1], NBinsPt);

  //////////////////////////////////////////////////////////////////////////////////////
  // Post-processing of reference flow for cumulants
  //////////////////////////////////////////////////////////////////////////////////////
  if (fCumulant) {

    TProfile *Corr22Ref = new TProfile(
        "Corr22Ref", "Profile of <2> with n=2 for ref", NBinsMult, Cent);
    TProfile *Corr24Ref = new TProfile(
        "Corr24Ref", "Profile of <4> with n=2 for ref", NBinsMult, Cent);

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
    hist_c22REF->GetXaxis()->SetTitle("Centrality (%)");
    hist_c22REF->GetYaxis()->SetTitle("c^{Ref}_{2}{2}");
    TH1F *hist_v22REF = new TH1F("v22REF", "v22REF", NBinsMult, Cent);
    hist_v22REF->GetXaxis()->SetTitle("Centrality (%)");
    hist_v22REF->GetYaxis()->SetTitle("v^{Ref}_{2}{2}");
    TH1F *hist_c24REF = new TH1F("c24REF", "c24REF", NBinsMult, Cent);
    hist_c24REF->GetXaxis()->SetTitle("Centrality (%)");
    hist_c24REF->GetYaxis()->SetTitle("c^{Ref}_{2}{4}");
    TH1F *hist_v24REF = new TH1F("v24REF", "v24REF", NBinsMult, Cent);
    hist_v24REF->GetXaxis()->SetTitle("Centrality (%)");
    hist_v24REF->GetYaxis()->SetTitle("v^{Ref}_{2}{4}");

    // Loop over all event Q-vectors
    for (int i = 0; i < fChain_REF->GetEntries(); i++) {

      fChain_REF->GetEntry(i);
      if (isnan(Corr2Ref) || isnan(Corr4Ref) || isinf(Corr2Ref) ||
          isinf(Corr4Ref)) {
        continue;
      }

      Corr22Ref->Fill(CentFT0REF, Corr2Ref, M11Ref);
      Corr24Ref->Fill(CentFT0REF, Corr4Ref, M1111Ref);
    }

    // Events averaged reference flow in given centrality range
    TProfile Corr22RefAll;
    Corr22Ref->Copy(Corr22RefAll);
    TProfile Corr24RefAll;
    Corr24Ref->Copy(Corr24RefAll);
    // Rebin thep profile to cover centrality range for differential study
    Corr22RefAll.SetBins(1, dimuonCentRange[0], dimuonCentRange[1]);
    Corr24RefAll.SetBins(1, dimuonCentRange[0], dimuonCentRange[1]);
    // Evaluate v2 in centrality range
    float cor2All = Corr22RefAll.GetBinContent(1);
    float cor2eAll = Corr22RefAll.GetBinError(1);
    float v22All = cor2All < 0 ? 0. : pow(cor2All, 1. / 2);
    float v22eAll = cor2All < 0 ? 0. : 0.5 * cor2eAll * pow(cor2All, 1. / 2);
    float cor4All = Corr24RefAll.GetBinContent(1);
    float cor4eAll = Corr24RefAll.GetBinError(1);
    float c24All = cor4All - 2 * pow(cor2All, 2);
    float c24eAll =
        pow(cor4eAll * cor4eAll + 16 * cor2All * cor2All * cor2eAll * cor2eAll,
            1. / 2);
    float v24All = c24All > 0 ? 0. : pow(-c24All, 1. / 4);
    float v24eAll = c24All > 0 ? 0. : pow(-c24All, (-3. / 4)) * c24eAll / 4;

    // Centrality-dependent reference flow
    for (int i = 0; i < NBinsMult; i++) {
      // Fill v2{2}
      float cor2 = Corr22Ref->GetBinContent(i + 1);
      float cor2e = Corr22Ref->GetBinError(i + 1);
      float v22 = cor2 < 0 ? 0. : pow(cor2, 1. / 2);
      float v22e = cor2 < 0 ? 0. : 0.5 * cor2e * pow(cor2, 1. / 2);
      hist_c22REF->SetBinContent(i + 1, cor2);
      hist_c22REF->SetBinError(i + 1, cor2e);
      hist_v22REF->SetBinContent(i + 1, v22);
      hist_v22REF->SetBinError(i + 1, v22e);

      // Fill v2{4}
      float cor4 = Corr24Ref->GetBinContent(i + 1);
      float cor4e = Corr24Ref->GetBinError(i + 1);
      float c24 = cor4 - 2 * pow(cor2, 2);
      float c24e =
          pow(cor4e * cor4e + 16 * cor2 * cor2 * cor2e * cor2e, 1. / 2);
      float v24 = c24 > 0 ? 0. : pow(-c24, 1. / 4);
      float v24e = c24 > 0 ? 0. : pow(-c24, (-3. / 4)) * c24e / 4;
      hist_c24REF->SetBinContent(i + 1, c24);
      hist_c24REF->SetBinError(i + 1, c24e);
      hist_v24REF->SetBinContent(i + 1, v24);
      hist_v24REF->SetBinError(i + 1, v24e);
    }

    TCanvas *c1REF = new TCanvas("c22REF");
    TCanvas *c2REF = new TCanvas("v22REF");
    TCanvas *c3REF = new TCanvas("c24REF");
    TCanvas *c4REF = new TCanvas("v24REF");

    c1REF->cd();
    hist_c22REF->SetMarkerStyle(20);
    hist_c22REF->SetStats(0);
    if (fPropError) {
      hist_c22REF->Draw("EP");
    } else {
      hist_c22REF->Draw("HIST P");
    }
    c2REF->cd();
    hist_v22REF->SetMarkerStyle(20);
    hist_v22REF->SetStats(0);
    if (fPropError) {
      hist_v22REF->Draw("EP");
    } else {
      hist_v22REF->Draw("HIST P");
    }
    c3REF->cd();
    hist_c24REF->SetMarkerStyle(20);
    hist_c24REF->SetStats(0);
    if (fPropError) {
      hist_c24REF->Draw("EP");
    } else {
      hist_c24REF->Draw("HIST P");
    }
    c4REF->cd();
    hist_v24REF->SetMarkerStyle(20);
    hist_v24REF->SetStats(0);
    if (fPropError) {
      hist_v24REF->Draw("EP");
    } else {
      hist_v24REF->Draw("HIST P");
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Post-processing for POI flow
  //////////////////////////////////////////////////////////////////////////////////////
  TProfile3D *Corr22Poi =
      new TProfile3D("Corr22Poi", "Profile of <2> with n=2 for poi", NBinsMass,
                     MassBins, NBinsPt, PtBins, NBinsMult, Cent);
  TProfile3D *Corr24Poi =
      new TProfile3D("Corr24Poi", "Profile of <4> with n=2 for poi", NBinsMass,
                     MassBins, NBinsPt, PtBins, NBinsMult, Cent);

  TProfile3D *U2Q2 = new TProfile3D("U2Q2", "U2Q2", NBinsMass, MassBins,
                                    NBinsPt, PtBins, NBinsMult, Cent);
  TProfile3D *R2SP = new TProfile3D("R2SP", "Resolution SP", NBinsMass,
                                    MassBins, NBinsPt, PtBins, NBinsMult, Cent);
  TProfile3D *Cos2DeltaPhi =
      new TProfile3D("Cos2DeltaPhi", "Cos2DeltaPhi", NBinsMass, MassBins,
                     NBinsPt, PtBins, NBinsMult, Cent);
  TProfile3D *R2EP = new TProfile3D("R2EP", "Resolution EP", NBinsMass,
                                    MassBins, NBinsPt, PtBins, NBinsMult, Cent);

  float M01Poi = 0.;
  float M0111Poi = 0.;
  float Corr2Poi = 0.;
  float Corr4Poi = 0.;
  float CentFT0POI;
  float fMass;
  float fPt, fPt1, fPt2;
  float fEta, fEta1, fEta2;
  int fMultDimuons = 0.;
  float fU2Q2, fU3Q3, fCos2DeltaPhi, fCos3DeltaPhi;
  float fR2EP, fR2SP;

  if (fCumulant) {
    fChain_POI->SetBranchAddress("fM01POI", &M01Poi);
    fChain_POI->SetBranchAddress("fM0111POI", &M0111Poi);
    fChain_POI->SetBranchAddress("fCORR2POI", &Corr2Poi);
    fChain_POI->SetBranchAddress("fCORR4POI", &Corr4Poi);
    fChain_POI->SetBranchAddress("fMultDimuons", &fMultDimuons);
  }

  fChain_POI->SetBranchAddress("fCentFT0C", &CentFT0POI);
  fChain_POI->SetBranchAddress("fMass", &fMass);
  fChain_POI->SetBranchAddress("fPt", &fPt);
  fChain_POI->SetBranchAddress("fPt1", &fPt1);
  fChain_POI->SetBranchAddress("fPt2", &fPt2);
  fChain_POI->SetBranchAddress("fEta", &fEta);
  fChain_POI->SetBranchAddress("fEta1", &fEta1);
  fChain_POI->SetBranchAddress("fEta2", &fEta2);
  fChain_POI->SetBranchAddress("fU2Q2", &fU2Q2);
  fChain_POI->SetBranchAddress("fU3Q3", &fU3Q3);
  fChain_POI->SetBranchAddress("fCos2DeltaPhi", &fCos2DeltaPhi);
  fChain_POI->SetBranchAddress("fCos3DeltaPhi", &fCos3DeltaPhi);
  fChain_POI->SetBranchAddress("fR2EP", &fR2EP);
  fChain_POI->SetBranchAddress("fR2SP", &fR2SP);

  TH1F *hist_c22POIMass =
      new TH1F("c22POIMass", "c22POIMass", NBinsMass, MassBins);
  hist_c22POIMass->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
  hist_c22POIMass->GetYaxis()->SetTitle("c^{J/#psi}_{2}{2}");
  TH1F *hist_v22POIMass =
      new TH1F("v22POIMass", "v22POIMass", NBinsMass, MassBins);
  hist_v22POIMass->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
  hist_v22POIMass->GetYaxis()->SetTitle("v^{J/#psi}_{2}{2}");
  TH1F *hist_c24POIMass =
      new TH1F("c24POIMass", "c24POIMass", NBinsMass, MassBins);
  hist_c24POIMass->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
  hist_c24POIMass->GetYaxis()->SetTitle("c^{J/#psi}_{2}{4}");
  TH1F *hist_v24POIMass =
      new TH1F("v24POIMass", "v24POIMass", NBinsMass, MassBins);
  hist_v24POIMass->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
  hist_v24POIMass->GetYaxis()->SetTitle("v^{J/#psi}_{2}{4}");

  TH1F *hist_v2SP = new TH1F("v2SP", "v2SP", NBinsMass, MassBins);
  hist_v2SP->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
  hist_v2SP->GetYaxis()->SetTitle("v^{J/#psi}_{2}{SP}");
  TH1F *hist_v2EP = new TH1F("v2EP", "v2EP", NBinsMass, MassBins);
  hist_v2EP->GetXaxis()->SetTitle("M_{#mu #mu} (GeV/c^{2})");
  hist_v2EP->GetYaxis()->SetTitle("v^{J/#psi}_{2}{EP}");

  // Loop over all dimuons
  for (int i = 0; i < fChain_POI->GetEntries(); i++) {
    fChain_POI->GetEntry(i);
    // Dimuons general selection
    // {Pt range, msass bin, centrality range}
    if (!(fPt > dimuonPtRange[0] && fPt <= dimuonPtRange[1] &&
          fMass > dimuonMassRange[0] && fMass <= dimuonMassRange[1] &&
          fEta > -4. && fEta < -2.5 && CentFT0POI > dimuonCentRange[0] &&
          CentFT0POI <= dimuonCentRange[1])) {
      continue;
    }
    // Fill (mass, pt, centrality) bins for cumulants
    if (fCumulant) {
      if (!(isnan(Corr2Poi) || isnan(Corr4Poi) || isinf(Corr2Poi) ||
            isinf(Corr4Poi))) {
        Corr22Poi->Fill(fMass, fPt, CentFT0POI, Corr2Poi, M01Poi);
        Corr24Poi->Fill(fMass, fPt, CentFT0POI, Corr4Poi, M0111Poi);
      }
    }

    // Fill (mass, pt, centrality) bins for SP and EP
    if (!(isnan(fR2SP) || isinf(fR2SP) || isnan(fR2EP) || isinf(fR2EP))) {
      U2Q2->Fill(fMass, fPt, CentFT0POI, fU2Q2);
      R2SP->Fill(fMass, fPt, CentFT0POI, fR2SP);
      Cos2DeltaPhi->Fill(fMass, fPt, CentFT0POI, fCos2DeltaPhi);
      R2EP->Fill(fMass, fPt, CentFT0POI, fR2EP);
    }
  }

  // Mass-dependent flow for POI
  for (int i = 0; i < NBinsMass; i++) {

    // Scalar-Product & Event-Plane method
    if (V2MultBinMass[i] != 0.) {
      V2SP[i] = V2SP[i] / V2MultBinMass[i];
      R2SP[i] = R2SP[i] / V2MultBinMass[i];
      V2ESP[i] = V2ESP[i] / V2MultBinMass[i];
      R2ESP[i] = R2ESP[i] / V2MultBinMass[i];
      V2ESP[i] = pow(V2ESP[i] - V2SP[i] * V2SP[i], 1. / 2);
      R2ESP[i] = pow(R2ESP[i] - R2SP[i] * R2SP[i], 1. / 2);
      V2SP[i] = R2SP[i] > 0 ? V2SP[i] / pow(R2SP[i], 1. / 2) : 0.;
      V2ESP[i] = pow(
          (V2ESP[i] * V2ESP[i] / R2SP[i] +
           0.25 * V2SP[i] * V2SP[i] * R2ESP[i] * R2ESP[i] / pow(R2SP[i], 3)) /
              V2MultBinMass[i],
          1. / 2);

      V2EP[i] = V2EP[i] / V2MultBinMass[i];
      R2EP[i] = R2EP[i] / V2MultBinMass[i];
      V2EEP[i] = V2EEP[i] / V2MultBinMass[i];
      R2EEP[i] = R2EEP[i] / V2MultBinMass[i];
      V2EEP[i] = pow(V2EEP[i] - V2EP[i] * V2EP[i], 1. / 2);
      R2EEP[i] = pow(R2EEP[i] - R2EP[i] * R2EP[i], 1. / 2);
      V2EP[i] = R2EP[i] > 0 ? V2EP[i] / pow(R2EP[i], 1. / 2) : 0.;
      V2EEP[i] = pow(
          (V2EEP[i] * V2EEP[i] / R2EP[i] +
           0.25 * V2EP[i] * V2EP[i] * R2EEP[i] * R2EEP[i] / pow(R2EP[i], 3)) /
              V2MultBinMass[i],
          1. / 2);

      if (fPropError) {
        hist_v2SP->SetBinContent(i + 1, V2SP[i]);
        hist_v2EP->SetBinContent(i + 1, V2EP[i]);
        if (isnan(V2ESP[i]) || isnan(V2EEP[i])) {
          hist_v2SP->SetBinError(i + 1, 0.);
          hist_v2EP->SetBinError(i + 1, 0.);
        } else {
          hist_v2SP->SetBinError(i + 1, V2ESP[i]);
          hist_v2EP->SetBinError(i + 1, V2EEP[i]);
        }
      } else {
        hist_v2SP->SetBinContent(i + 1, V2SP[i]);
        hist_v2EP->SetBinContent(i + 1, V2EP[i]);
      }
    }

    // Multiparticle Cumulant method
    if (fCumulant) {
      if (Sum22POI[i] * Sum24POI[i] != 0.) {
        C22POI[i] = C22POI[i] / Sum22POI[i];
        C22EPOI[i] = C22EPOI[i] / Sum22POI[i];
        C22EPOI[i] = pow(C22EPOI[i] - C22POI[i] * C22POI[i], 1. / 2);

        V22POI[i] = C22POI[i] * pow(C22REF, -1. / 2);
        V22EPOI[i] = pow(C22EPOI[i] * C22EPOI[i] / (C22REF * Sum22POI[i]) +
                             0.25 * C22POI[i] * C22POI[i] * C22EREF * C22EREF /
                                 (pow(C22REF, 3) * Sum22REFAll),
                         1. / 2);
        V22EPOI[i] = V22EPOI[i] > 0 ? V22EPOI[i] : 0.;

        C24POI[i] = C24POI[i] / Sum24POI[i] - 2. * C22POI[i] * C22REF;
        C24EPOI[i] = C24EPOI[i] / Sum24POI[i];
        C24EPOI[i] = pow(C24EPOI[i] - C24POI[i] * C24POI[i], 1. / 2);
        C24EPOI[i] = pow(
            C24EPOI[i] * C24EPOI[i] / Sum24POI[i] +
                4. * C22REF * C22REF * C22EPOI[i] * C22EPOI[i] / Sum22POI[i] +
                4. * C22POI[i] * C22POI[i] * C22EREF * C22EREF / Sum22REFAll,
            1. / 2);

        V24POI[i] = -1. * C24POI[i] * pow(-1. * C24REF, -3. / 4);
        V24EPOI[i] = pow(
            pow(-1. * C24REF, -6. / 4) * C24EPOI[i] * C24EPOI[i] / Sum24POI[i] +
                pow(-1. * C24REF, -14. / 4) * C24POI[i] * C24POI[i] * C24EREF *
                    C24EREF / Sum24REFAll * 9. / 16,
            1. / 2);
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
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Plot for results, save to histograms for fit
  //////////////////////////////////////////////////////////////////////////////////////
  if (fCumulant) {
    TCanvas *c1POI = new TCanvas("c22POI");
    TCanvas *c2POI = new TCanvas("v22POI");
    TCanvas *c3POI = new TCanvas("c24POI");
    TCanvas *c4POI = new TCanvas("v24POI");
    c1POI->cd();
    hist_c22POIMass->SetMarkerStyle(20);
    hist_c22POIMass->SetStats(0);
    if (fPropError) {
      hist_c22POIMass->Draw("EP");
    } else {
      hist_c22POIMass->Draw("HIST P");
    }
    c2POI->cd();
    hist_v22POIMass->SetMarkerStyle(20);
    hist_v22POIMass->SetStats(0);
    if (fPropError) {
      hist_v22POIMass->Draw("EP");
    } else {
      hist_v22POIMass->Draw("HIST P");
    }
    c3POI->cd();
    hist_c24POIMass->SetMarkerStyle(20);
    hist_c24POIMass->SetStats(0);
    if (fPropError) {
      hist_c24POIMass->Draw("EP");
    } else {
      hist_c24POIMass->Draw("HIST P");
    }
    c4POI->cd();
    hist_v24POIMass->SetMarkerStyle(20);
    hist_v24POIMass->SetStats(0);
    if (fPropError) {
      hist_v24POIMass->Draw("EP");
    } else {
      hist_v24POIMass->Draw("HIST P");
    }
  }

  TCanvas *cSP = new TCanvas("v2SP");
  TCanvas *cEP = new TCanvas("v2EP");
  cSP->cd();
  hist_v2SP->SetMarkerStyle(20);
  hist_v2SP->SetStats(0);
  if (fPropError) {
    hist_v2SP->Draw("EP");
  } else {
    hist_v2SP->Draw("HIST P");
  }
  cEP->cd();
  hist_v2EP->SetMarkerStyle(20);
  hist_v2EP->SetStats(0);
  if (fPropError) {
    hist_v2EP->Draw("EP");
  } else {
    hist_v2EP->Draw("HIST P");
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