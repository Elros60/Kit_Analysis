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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Framework/Logger.h"

using namespace std;

void LoadData(TChain *fChain, std::string TreeName = "O2rerefflow",
              std::string FileName = "AO2D.root");
void CreateBins(double *axis, double min, double max, int Nbins = 10);

void RefFlowVsCent(std::string FileName = "input_data.txt") {

  // Load data from AO2Ds
  TChain *fChain_Flow = new TChain();
  fstream InputFiles;
  InputFiles.open(FileName, ios::in);
  if (InputFiles.is_open()) {
    string File;
    cout << "Start reading input AO2Ds..." << endl;
    while (getline(InputFiles, File)) {
      cout << "Reading input from: " << File << endl;
      LoadData(fChain_Flow, "O2rtdileptonflow", File);
    }
    InputFiles.close();
  }

  // Retrieve and setup binning information
  int NBinsMult = 10;
  double Cent[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};

  // Define variables to be used in post-processing
  float Mass, Pt, Eta, Phi, CentFT0C, M11REF, M1111REF, CORR2REF, CORR4REF,
      M01POI, M0111POI, CORR2POI, CORR4POI;
  int MultDimuons;

  fChain_Flow->SetBranchAddress("fMass", &Mass);
  fChain_Flow->SetBranchAddress("fPt", &Pt);
  fChain_Flow->SetBranchAddress("fEta", &Eta);
  fChain_Flow->SetBranchAddress("fPhi", &Phi);
  fChain_Flow->SetBranchAddress("fM11REF", &M11REF);
  fChain_Flow->SetBranchAddress("fM1111REF", &M1111REF);
  fChain_Flow->SetBranchAddress("fCORR2REF", &CORR2REF);
  fChain_Flow->SetBranchAddress("fCORR4REF", &CORR4REF);
  fChain_Flow->SetBranchAddress("fM01POI", &M01POI);
  fChain_Flow->SetBranchAddress("fM0111POI", &M0111POI);
  fChain_Flow->SetBranchAddress("fCORR2POI", &CORR2POI);
  fChain_Flow->SetBranchAddress("fCORR4POI", &CORR4POI);
  fChain_Flow->SetBranchAddress("fCentFT0C", &CentFT0C);
  fChain_Flow->SetBranchAddress("fMultDimuons", &MultDimuons);

  // Define TProfiles for averages
  TProfile *Corr22Ref =
      new TProfile("Corr22Ref", "Profile of <<2>>", NBinsMult, Cent);
  TProfile *Corr24Ref =
      new TProfile("Corr24Ref", "Profile of <<4>>", NBinsMult, Cent);

  // Define histograms for correlations
  TH1F *hist_c22REF = new TH1F("c22REF", "c22REF", NBinsMult, Cent);
  hist_c22REF->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_c22REF->GetYaxis()->SetTitle("c^{Ref}_{2}{2}");

  TH1F *hist_v22REF = new TH1F("v22REF", "v22REF", NBinsMult, Cent);
  hist_v22REF->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_v22REF->GetYaxis()->SetTitle("v^{Ref}_{2}{2}");

  TH1F *hist_c24REF = new TH1F("c24REF", "c24REF", NBinsMult, Cent);
  hist_c24REF->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_c24REF->GetYaxis()->SetTitle("c^{Ref}_{2}{4}");

  TH1F *hist_v24REF = new TH1F("v24REF", "v24REF", NBinsMult, Cent);
  hist_v24REF->GetXaxis()->SetTitle("Centrality FT0C %");
  hist_v24REF->GetYaxis()->SetTitle("v^{Ref}_{2}{4}");

  // Loop over all dimuons for both REF
  for (int i = 0; i < fChain_Flow->GetEntries(); i++) {
    fChain_Flow->GetEntry(i);
    // Correlations averaged over events within the same centrality bin
    // For Reference Flow
    if (!(isnan(CORR2REF) || isnan(CORR4REF) || isinf(CORR2REF) ||
          isinf(CORR4REF))) {
      if (MultDimuons != 0) {
        Corr22Ref->Fill(CentFT0C, CORR2REF, M11REF / MultDimuons);
        Corr24Ref->Fill(CentFT0C, CORR4REF, M1111REF / MultDimuons);
      }
    }
  }

  // Extract flow coefficient from correlations
  for (int i = 0; i < NBinsMult; i++) {
    // Fill v2{2} with error propagation
    float cor2_ref = Corr22Ref->GetBinContent(i + 1);
    float cor2e_ref = Corr22Ref->GetBinError(i + 1);
    float v22_ref = cor2_ref < 0 ? 0. : pow(cor2_ref, 1. / 2);
    float v22e_ref =
        cor2_ref < 0 ? 0. : 0.5 * cor2e_ref * pow(cor2_ref, 1. / 2);
    hist_c22REF->SetBinContent(i + i, cor2_ref);
    hist_c22REF->SetBinError(i + i, cor2e_ref);
    hist_v22REF->SetBinContent(i + i, v22_ref);
    hist_v22REF->SetBinError(i + i, v22e_ref);

    // Fill v2{4} with error propagation
    float cor4_ref = Corr24Ref->GetBinContent(i + 1);
    float cor4e_ref = Corr24Ref->GetBinError(i + 1);
    float c24 = cor4_ref - 2 * pow(cor2_ref, 2);
    float c24e = pow(cor4e_ref * cor4e_ref +
                         16 * cor2_ref * cor2_ref * cor2e_ref * cor2e_ref,
                     1. / 2);
    float v24_ref = c24 > 0 ? 0. : pow(-c24, 1. / 4);
    float v24e_ref = c24 > 0 ? 0. : pow(-c24, (-3. / 4)) * c24e / 4;
    hist_c24REF->SetBinContent(i + 1, c24);
    hist_c24REF->SetBinError(i + 1, c24e);
    hist_v24REF->SetBinContent(i + 1, v24_ref);
    hist_v24REF->SetBinError(i + 1, v24e_ref);
  }

  // Save outputs
  TFile f("AnalysisResults_RefFlowVsCent.root", "RECREATE");
  TList *l = new TList();
  l->Add(hist_c22REF);
  l->Add(hist_c24REF);
  l->Add(hist_v22REF);
  l->Add(hist_v24REF);
  l->Write("ReferenceFlowVsCent", TObject::kSingleKey);
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

void CreateBins(double *axis, double min, double max, int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}