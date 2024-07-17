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
#include "TGraphAsymmErrors.h"
#include <fstream>
#include <vector>
#include "TGraphErrors.h"
#include "TAxis.h"
#include <sstream>
#include <TStyle.h>

#include "Framework/Logger.h"

using namespace std;

void LoadData(TChain *fChain, std::string TreeName = "O2rerefflow",
              std::string FileName = "AOD.root");
void CreateBins(double *axis, double min, double max, int Nbins = 10);

void writeToFile(const char *filename, TH1F *hist);

void run2_vs_run3();

void plotData(const char *filename);

void FlowAnalysis(std::array<float, 3> dimuonMassRange,
                  std::array<float, 3> dimuonPtRange,
                  std::array<float, 2> dimuonCentRange,
                  bool fCumulant = false)
{

  TChain *fChain_REF = new TChain();
  if (fCumulant)
  {
    LoadData(fChain_REF, "O2rerefflow");
  }

  TChain *fChain_POI = new TChain();
  LoadData(fChain_POI, "O2rtdimuonall");

  int NBinsMass = static_cast<int>(dimuonMassRange[2]);
  int NBinsPt = static_cast<int>(dimuonPtRange[2]);
  int NBinsMult = 10;
  double Cent[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  double *MassBins = new double[NBinsMass + 1];
  double *PtBins = new double[NBinsPt + 1];

  CreateBins(MassBins, dimuonMassRange[0], dimuonMassRange[1], NBinsMass);
  CreateBins(PtBins, dimuonPtRange[0], dimuonPtRange[1], NBinsPt);

  //////////////////////////////////////////////////////////////////////////////////////
  // Post-processing of reference flow for cumulants
  //////////////////////////////////////////////////////////////////////////////////////
  if (fCumulant)
  {
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
    for (int i = 0; i < fChain_REF->GetEntries(); i++)
    {

      fChain_REF->GetEntry(i);
      if (isnan(Corr2Ref) || isnan(Corr4Ref) || isinf(Corr2Ref) ||
          isinf(Corr4Ref))
      {
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
    for (int i = 0; i < NBinsMult; i++)
    {
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
    // Write v22REF to a file
    writeToFile("v22REF_Run3_w.txt", hist_v22REF);

    // Write v24REF to a file
    writeToFile("v24REF_Run3_w.txt", hist_v24REF);

    run2_vs_run3();

    TCanvas *c1REF = new TCanvas("c22REF");
    TCanvas *c2REF = new TCanvas("v22REF");
    TCanvas *c3REF = new TCanvas("c24REF");
    TCanvas *c4REF = new TCanvas("v24REF");

    c1REF->cd();
    hist_c22REF->SetMarkerStyle(20);
    hist_c22REF->SetStats(0);
    hist_c22REF->Draw("EP");

    c2REF->cd();
    hist_v22REF->SetMarkerStyle(20);
    hist_v22REF->SetStats(0);
    hist_v22REF->Draw("EP");

    c3REF->cd();
    hist_c24REF->SetMarkerStyle(20);
    hist_c24REF->SetStats(0);
    hist_c24REF->Draw("EP");

    c4REF->cd();
    hist_v24REF->SetMarkerStyle(20);
    hist_v24REF->SetStats(0);
    hist_v24REF->Draw("EP");
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Post-processing for POI flow
  //////////////////////////////////////////////////////////////////////////////////////
  /*
  TProfile3D *Corr22Poi =
      new TProfile3D("Corr22Poi", "Profile of <2> with n=2 for poi", NBinsMass,
                     MassBins, NBinsPt, PtBins, NBinsMult, Cent);
  TProfile3D *Corr24Poi =
      new TProfile3D("Corr24Poi", "Profile of <4> with n=2 for poi", NBinsMass,
                     MassBins, NBinsPt, PtBins, NBinsMult, Cent);
  */

  TProfile2D *U2Q2 =
      new TProfile2D("U2Q2", "U2Q2", NBinsMass, MassBins, NBinsPt, PtBins);
  TProfile2D *R2SP = new TProfile2D("R2SP", "Resolution SP", NBinsMass,
                                    MassBins, NBinsPt, PtBins);
  TProfile2D *Cos2DeltaPhi = new TProfile2D(
      "Cos2DeltaPhi", "Cos2DeltaPhi", NBinsMass, MassBins, NBinsPt, PtBins);
  TProfile2D *R2EP = new TProfile2D("R2EP", "Resolution EP", NBinsMass,
                                    MassBins, NBinsPt, PtBins);

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

  if (fCumulant)
  {
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
  for (int i = 0; i < fChain_POI->GetEntries(); i++)
  {
    fChain_POI->GetEntry(i);
    // Dimuons general selection
    // {Pt range, msass bin, centrality range}
    if (!(fPt > dimuonPtRange[0] && fPt <= dimuonPtRange[1] &&
          fMass > dimuonMassRange[0] && fMass <= dimuonMassRange[1] &&
          fEta > -4. && fEta < -2.5 && CentFT0POI > dimuonCentRange[0] &&
          CentFT0POI <= dimuonCentRange[1]))
    {
      continue;
    }
    /*
    // Fill (mass, pt, centrality) bins for cumulants
    if (fCumulant) {
      if (!(isnan(Corr2Poi) || isnan(Corr4Poi) || isinf(Corr2Poi) ||
            isinf(Corr4Poi))) {
        Corr22Poi->Fill(fMass, fPt, CentFT0POI, Corr2Poi, M01Poi);
        Corr24Poi->Fill(fMass, fPt, CentFT0POI, Corr4Poi, M0111Poi);
      }
    }
    */

    // Fill (mass, pt, centrality) bins for SP and EP
    if (!(isnan(fR2SP) || isinf(fR2SP) || isnan(fR2EP) || isinf(fR2EP)))
    {
      U2Q2->Fill(fMass, fPt, fU2Q2);
      R2SP->Fill(fMass, fPt, fR2SP);
      Cos2DeltaPhi->Fill(fMass, fPt, fCos2DeltaPhi);
      R2EP->Fill(fMass, fPt, fR2EP);
    }
  }

  // Mass-dependent flow for POI
  TProfile *U2Q2Mass = U2Q2->ProfileX("u2q2mass", 1, NBinsPt);
  TProfile *R2SPMass = R2SP->ProfileX("r2spmass", 1, NBinsPt);
  TProfile *Cos2DeltaPhiMass =
      Cos2DeltaPhi->ProfileX("cos2deltaphimass", 1, NBinsPt);
  TProfile *R2EPMass = R2EP->ProfileX("r2epmass", 1, NBinsPt);

  for (int i = 0; i < NBinsMass; i++)
  {

    // Scalar-Product & Event-Plane method
    float u2q2 = U2Q2Mass->GetBinContent(i + 1);
    float u2q2e = U2Q2Mass->GetBinError(i + 1);
    float r2sp = R2SPMass->GetBinContent(i + 1);
    float r2spe = R2SPMass->GetBinError(i + 1);
    float v2sp = u2q2 / pow(r2sp, 1. / 2);
    float v2spe = pow(u2q2e * u2q2e / r2sp +
                          0.25 * pow(r2sp, -3) * u2q2 * u2q2 * r2spe * r2spe,
                      1. / 2);

    float cos2deltaphi = Cos2DeltaPhiMass->GetBinContent(i + 1);
    float cos2deltaphie = Cos2DeltaPhiMass->GetBinError(i + 1);
    float r2ep = R2EPMass->GetBinContent(i + 1);
    float r2epe = R2EPMass->GetBinError(i + 1);
    float v2ep = cos2deltaphi / pow(r2ep, 1. / 2);
    float v2epe = pow(cos2deltaphie * cos2deltaphie / r2ep +
                          0.25 * pow(r2ep, -3) * cos2deltaphi * cos2deltaphi *
                              r2epe * r2epe,
                      1. / 2);

    hist_v2SP->SetBinContent(i + 1, v2sp);
    hist_v2EP->SetBinContent(i + 1, v2ep);
    hist_v2SP->SetBinError(i + 1, isnan(v2spe) ? 0. : v2spe);
    hist_v2EP->SetBinError(i + 1, isnan(v2epe) ? 0. : v2epe);

    /*
      // Multiparticle Cumulant method
      if (fCumulant) {
        if (Sum22POI[i] * Sum24POI[i] != 0.) {
          C22POI[i] = C22POI[i] / Sum22POI[i];
          C22EPOI[i] = C22EPOI[i] / Sum22POI[i];
          C22EPOI[i] = pow(C22EPOI[i] - C22POI[i] * C22POI[i], 1. / 2);

          V22POI[i] = C22POI[i] * pow(C22REF, -1. / 2);
          V22EPOI[i] = pow(C22EPOI[i] * C22EPOI[i] / (C22REF * Sum22POI[i]) +
                               0.25 * C22POI[i] * C22POI[i] * C22EREF * C22EREF
      / (pow(C22REF, 3) * Sum22REFAll),
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
              pow(-1. * C24REF, -6. / 4) * C24EPOI[i] * C24EPOI[i] / Sum24POI[i]
      + pow(-1. * C24REF, -14. / 4) * C24POI[i] * C24POI[i] * C24EREF * C24EREF
      / Sum24REFAll * 9. / 16,
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
      */
  }

  //////////////////////////////////////////////////////////////////////////////////////
  // Plot for results, save to histograms for fit
  //////////////////////////////////////////////////////////////////////////////////////
  if (fCumulant)
  {
    TCanvas *c1POI = new TCanvas("c22POI");
    TCanvas *c2POI = new TCanvas("v22POI");
    TCanvas *c3POI = new TCanvas("c24POI");
    TCanvas *c4POI = new TCanvas("v24POI");
    c1POI->cd();
    hist_c22POIMass->SetMarkerStyle(20);
    hist_c22POIMass->SetStats(0);
    hist_c22POIMass->Draw("EP");

    c2POI->cd();
    hist_v22POIMass->SetMarkerStyle(20);
    hist_v22POIMass->SetStats(0);
    hist_v22POIMass->Draw("EP");

    c3POI->cd();
    hist_c24POIMass->SetMarkerStyle(20);
    hist_c24POIMass->SetStats(0);
    hist_c24POIMass->Draw("EP");

    c4POI->cd();
    hist_v24POIMass->SetMarkerStyle(20);
    hist_v24POIMass->SetStats(0);
    hist_v24POIMass->Draw("EP");
  }

  TCanvas *cSP = new TCanvas("v2SP");
  TCanvas *cEP = new TCanvas("v2EP");
  cSP->cd();
  hist_v2SP->SetMarkerStyle(20);
  hist_v2SP->SetStats(0);
  hist_v2SP->Draw("EP");

  cEP->cd();
  hist_v2EP->SetMarkerStyle(20);
  hist_v2EP->SetStats(0);
  hist_v2EP->Draw("EP");
}

void LoadData(TChain *fChain, std::string TreeName, std::string FileName)
{

  TFile *fInput = TFile::Open(FileName.c_str());
  TIter keyList(fInput->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)keyList()))
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    std::string dir = key->GetName();
    if (dir == "parentFiles")
    {
      continue;
    }
    fChain->Add(
        Form("%s?#%s/%s", FileName.c_str(), dir.c_str(), TreeName.c_str()));
  }
}

void CreateBins(double *axis, double min, double max, int Nbins)
{
  for (int i = 0; i < Nbins; i++)
  {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}

void writeToFile(const char *filename, TH1F *hist)
{
  std::ofstream outFile(filename);
  if (!outFile.is_open())
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  // Write the header
  outFile << "bin_range\tbin_content\tbin_error\n";

  // Write the bin contents and errors to the file
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    double binMin = hist->GetXaxis()->GetBinLowEdge(i);
    double binMax = hist->GetXaxis()->GetBinUpEdge(i);
    double binContent = hist->GetBinContent(i);
    double binError = hist->GetBinError(i);
    outFile << binMin << " " << binMax << "\t" << binContent << "\t" << binError << std::endl;
  }

  outFile.close();
}

void plotData(const char *filename)
{
  // Open the input file
  std::ifstream inputFile(filename);
  if (!inputFile.is_open())
  {
    std::cerr << "Error: Unable to open file " << filename << std::endl;
    return;
  }

  // Create arrays to hold data
  const int nPoints = 9;
  double xmin[nPoints], xmax[nPoints], y[nPoints], error[nPoints];

  // Read data from file
  for (int i = 0; i < nPoints; ++i)
  {
    inputFile >> xmin[i] >> xmax[i] >> y[i] >> error[i];
  }

  // Close the input file
  inputFile.close();

  // Create TGraphErrors object
  TGraphErrors *graph = new TGraphErrors(nPoints, xmin, y, 0, error);
  graph->SetTitle("Data Plot");
  graph->GetXaxis()->SetTitle("Centrality");
  graph->GetYaxis()->SetTitle("v_{2}{} Run 3 ");

  // Set histogram style
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);

  // Create a canvas to draw the graph
  TCanvas *canvas = new TCanvas("canvas", "Data Plot", 800, 600);
  graph->Draw("AP"); // Draw the graph with error bars

  // Save the canvas as an image file
  canvas->SaveAs((std::string(filename) + ".pdf").c_str());
}

void run2_vs_run3()
{
  double minPt = 0.2;
  double maxPt = 5;

  TCanvas c("v2_ref_cent", "", 900, 600);
  gStyle->SetOptStat(0);

  TLegend *legend = new TLegend(0.64, 0.68, 0.86, 0.83);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.03);

  // Adjusted coordinates for the second legend
  TLegend *legend2 = new TLegend(0.3, 0.88, 0.85, 0.88);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.03);

  // RUN3 VALUES: From v22_cent.txt
  std::ifstream in22("v22REF_Run3_w.txt");
  std::string skipLine22;
  std::getline(in22, skipLine22); // Skip the first line

  Float_t min, max, V2, errorStatref;
  TGraphAsymmErrors StatGraph22;
  Int_t pointIndex22 = 0;

  while (in22 >> min >> max >> V2 >> errorStatref)
  {
    Float_t binCenter = (min + max) / 2;
    Double_t binWidth = (max - min) / 2;
    StatGraph22.SetPoint(pointIndex22, binCenter, V2);
    StatGraph22.SetPointError(pointIndex22, binWidth, binWidth, errorStatref, errorStatref);
    pointIndex22++;
  }
  in22.close();

  StatGraph22.SetTitle("; Centrality; v^{REF}_{2}");
  StatGraph22.SetMarkerStyle(20);
  StatGraph22.SetMarkerColor(kBlue);
  StatGraph22.SetLineColor(kBlue);

  // RUN3 VALUES: From v24_cent.txt (new data)
  std::ifstream in24("v24REF_Run3_w.txt");
  std::string skipLine24;
  std::getline(in24, skipLine24);

  TGraphAsymmErrors StatGraph24;
  Int_t pointIndex24 = 0;

  while (in24 >> min >> max >> V2 >> errorStatref)
  {
    Float_t binCenter = (min + max) / 2;
    Double_t binWidth = (max - min) / 2;
    StatGraph24.SetPoint(pointIndex24, binCenter, V2);
    StatGraph24.SetPointError(pointIndex24, binWidth, binWidth, errorStatref, errorStatref);
    pointIndex24++;
  }
  in24.close();

  StatGraph24.SetMarkerStyle(20);
  StatGraph24.SetMarkerColor(kRed);
  StatGraph24.SetLineColor(kRed);

  // Save the canvas to a file
  auto f = TFile::Open("run2vsrun3.root", "Recreate");

  // RUN2 VALUES:
  std::ifstream file("RUN2_v24.txt");

  // Skip the header
  std::string header;
  std::getline(file, header);

  const int nPoints = 9; // NUMBER OF POINTS
  double x[nPoints], y[nPoints], statError[nPoints], systError[nPoints];

  // Read the data
  for (int i = 0; i < nPoints && !file.eof(); ++i)
  {
    file >> x[i] >> y[i] >> statError[i] >> systError[i];
  }

  file.close();

  // Create a graph for the statistical errors
  TGraphErrors *graphStat = new TGraphErrors(nPoints, x, y, nullptr, statError);
  graphStat->SetTitle("Data with Errors;X-axis;Y-axis");
  graphStat->SetMarkerStyle(20);
  graphStat->SetMarkerColor(kGreen);
  graphStat->SetLineColor(kGreen);
  graphStat->Draw("AP");

  TGraphAsymmErrors StatGraphRUN2_1;
  StatGraphRUN2_1.SetMarkerStyle(20);
  StatGraphRUN2_1.SetMarkerColor(kGreen);
  StatGraphRUN2_1.SetLineColor(kGreen);

  TGraphAsymmErrors StatGraphRUN2_2;
  StatGraphRUN2_2.SetMarkerStyle(20);
  StatGraphRUN2_2.SetMarkerColor(kViolet);
  StatGraphRUN2_2.SetLineColor(kViolet);

  std::ifstream file2("RUN2_v22.txt");
  std::getline(file2, header);
  double x2[nPoints], y2[nPoints], statError2[nPoints], systError2[nPoints];

  // Read the data
  for (int i = 0; i < nPoints && !file2.eof(); ++i)
  {
    file2 >> x2[i] >> y2[i] >> statError2[i] >> systError2[i];
  }
  file.close();

  // Create a graph for the statistical errors
  TGraphErrors *graphStat2 = new TGraphErrors(nPoints, x2, y2, nullptr, statError2);
  graphStat2->SetTitle("Data with Errors;X-axis;Y-axis");
  graphStat2->SetMarkerStyle(20);
  graphStat2->SetMarkerColor(kViolet);
  graphStat2->SetLineColor(kViolet);

  legend->AddEntry(&StatGraphRUN2_2, "v_{2}{2} Run 2", "lp");
  legend->AddEntry(&StatGraphRUN2_1, "v_{2}{4} Run 2", "lp");
  legend->AddEntry(&StatGraph22, "v_{2}^{c}{2} Run 3", "lp");
  legend->AddEntry(&StatGraph24, "v_{2}^{c}{4} Run 3", "lp");

  // Convert double to string with 1 decimal place
  std::string minPtStr = std::to_string(minPt);
  std::string maxPtStr = std::to_string(maxPt);

  // Find the position of the decimal point
  size_t decimalPosMinPt = minPtStr.find('.');
  size_t decimalPosMaxPt = maxPtStr.find('.');

  // If the decimal point is found, truncate the string after 1 decimal place
  if (decimalPosMinPt != std::string::npos)
  {
    minPtStr = minPtStr.substr(0, decimalPosMinPt + 2);
  }
  if (decimalPosMaxPt != std::string::npos)
  {
    maxPtStr = maxPtStr.substr(0, decimalPosMaxPt + 2);
  }

  // Construct the centString using minPtStr and maxPtStr
  std::string centString = "ALICE, Pb-Pb #sqrt{s} = 5.36 TeV, -0.8 < y < 0.8, p_{T} " + minPtStr + "-" + maxPtStr + " GeV/c";
  legend2->AddEntry(&StatGraph24, centString.c_str(), "");

  // Drawing both graphs
  StatGraph22.Draw("AP");
  StatGraph24.Draw("P SAME");

  // Systematic errors
  for (int i = 0; i < nPoints; i++)
  {
    TBox *box = new TBox(x[i] - 0.1, y[i] - systError[i], x[i] + 0.1, y[i] + systError[i]);
    box->SetFillStyle(1001);              // Solid fill
    box->SetFillColorAlpha(kGreen, 0.35); // Semi-transparent red fill
    box->SetLineColor(kGreen);            // Solid red outline
    box->SetLineWidth(1);                 // Outline width
    box->Draw("l SAME");                  // Draw with line to include the outline
  }

  for (int i = 0; i < nPoints; i++)
  {
    TBox *box2 = new TBox(x2[i] - 0.1, y2[i] - systError2[i], x2[i] + 0.1, y2[i] + systError2[i]);
    box2->SetFillStyle(1001);               // Solid fill
    box2->SetFillColorAlpha(kViolet, 0.35); // Semi-transparent red fill
    box2->SetLineColor(kViolet);            // Solid red outline
    box2->SetLineWidth(1);                  // Outline width
    box2->Draw("l SAME");                   // Draw with line to include the outline
  }

  graphStat2->Draw("P SAME");
  graphStat->Draw("P SAME");

  legend->Draw("P SAME");
  legend2->Draw("P SAME");

  c.Write();
  f->Close();
}