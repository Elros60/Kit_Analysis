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

double *CreateBinsFromAxis(TAxis *axis);
void CreateBins(double *axis, double min, double max, int Nbins = 10);
void LoadData(std::string FileName, TH1F *&hs_Q2XA, TH1F *&hs_Q2YA,
              TH1F *&hs_Q2XB, TH1F *&hs_Q2YB, TH1F *&hs_Q2XC, TH1F *&hs_Q2YC);

//______________________________________________________________________________
void FlowAnalysis_ScalarProductQA(
    std::string FileName = "AnalysisResults.root") {
  TH1F *hs_QAX, *hs_QBX, *hs_QCX, *hs_QAY, *hs_QBY, *hs_QCY;
  LoadData(FileName, hs_QAX, hs_QAY, hs_QBX, hs_QBY, hs_QCX, hs_QCY);
  TFile f("FlowAnalysisResults_ScalarProductQA.root", "RECREATE");
  TList *ls = new TList();
  hs_QAX->GetXaxis()->SetRangeUser(-4., 4.);
  hs_QAY->GetXaxis()->SetRangeUser(-4., 4.);
  hs_QBX->GetXaxis()->SetRangeUser(-4., 4.);
  hs_QBY->GetXaxis()->SetRangeUser(-4., 4.);
  hs_QCX->GetXaxis()->SetRangeUser(-4., 4.);
  hs_QCY->GetXaxis()->SetRangeUser(-4., 4.);
  ls->Add(hs_QAX);
  ls->Add(hs_QAY);
  ls->Add(hs_QBX);
  ls->Add(hs_QBY);
  ls->Add(hs_QCX);
  ls->Add(hs_QCY);
  f.cd();
  ls->Write("QVectors", TObject::kSingleKey);
  f.Close();
}

//______________________________________________________________________________
double *CreateBinsFromAxis(TAxis *axis) {
  int Nbins = axis->GetNbins();
  double *Bins = new double[Nbins + 1];
  axis->GetLowEdge(Bins);
  Bins[Nbins] = axis->GetBinUpEdge(Nbins);
  return Bins;
}

//______________________________________________________________________________
void CreateBins(double *axis, double min, double max, int Nbins) {
  for (int i = 0; i < Nbins; i++) {
    axis[i] = min + i * (max - min) / Nbins;
  }
  axis[Nbins] = max;
}

//______________________________________________________________________________
void LoadData(std::string FileName, TH1F *&hs_Q2XA, TH1F *&hs_Q2YA,
              TH1F *&hs_Q2XB, TH1F *&hs_Q2YB, TH1F *&hs_Q2XC, TH1F *&hs_Q2YC) {
  // Load input data for analysis
  filesystem::path filePath = FileName;
  THashList *list_hist;
  TList *sublist;
  if (filePath.extension() == ".root") {
    // Load data from AnalysisResults.root
    TFile *Input_File = TFile::Open(FileName.c_str());
    list_hist = (THashList *)Input_File->Get("analysis-event-selection/output");
    sublist = (TList *)list_hist->FindObject("Event_AfterCuts");

    // Get histograms of correlations and resolution factors
    hs_Q2XA = (TH1F *)sublist->FindObject("Q2X0A");
    hs_Q2XB = (TH1F *)sublist->FindObject("Q2X0B");
    hs_Q2XC = (TH1F *)sublist->FindObject("Q2X0C");
    hs_Q2YA = (TH1F *)sublist->FindObject("Q2Y0A");
    hs_Q2YB = (TH1F *)sublist->FindObject("Q2Y0B");
    hs_Q2YC = (TH1F *)sublist->FindObject("Q2Y0C");
    Input_File->Close();
  } else {
    // Load data from a list of AnalysisResults.root
    fstream InputFiles;
    InputFiles.open(FileName, ios::in);
    if (InputFiles.is_open()) {
      string File;
      cout << "Start Loading input AnalysisResults in list..." << endl;
      bool first = true;
      while (getline(InputFiles, File)) {
        cout << "Loading input from: " << File << endl;
        TFile *inFile = TFile::Open(File.c_str());
        list_hist = (THashList *)inFile->Get("analysis-event-selection/output");
        sublist = (TList *)list_hist->FindObject("Event_AfterCuts");

        if (first) {
          hs_Q2XA = (TH1F *)sublist->FindObject("Q2X0A");
          hs_Q2XB = (TH1F *)sublist->FindObject("Q2X0B");
          hs_Q2XC = (TH1F *)sublist->FindObject("Q2X0C");
          hs_Q2YA = (TH1F *)sublist->FindObject("Q2Y0A");
          hs_Q2YB = (TH1F *)sublist->FindObject("Q2Y0B");
          hs_Q2YC = (TH1F *)sublist->FindObject("Q2Y0C");
        } else {
          hs_Q2XA->Add((TH1F *)sublist->FindObject("Q2X0A"));
          hs_Q2XB->Add((TH1F *)sublist->FindObject("Q2X0B"));
          hs_Q2XC->Add((TH1F *)sublist->FindObject("Q2X0C"));
          hs_Q2YA->Add((TH1F *)sublist->FindObject("Q2Y0A"));
          hs_Q2YB->Add((TH1F *)sublist->FindObject("Q2Y0B"));
          hs_Q2YC->Add((TH1F *)sublist->FindObject("Q2Y0C"));
        }
        inFile->Close();
        first = false;
      }
      InputFiles.close();
    }
  }
}