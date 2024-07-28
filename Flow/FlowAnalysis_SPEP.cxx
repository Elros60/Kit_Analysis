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
              std::string FileName = "AO2Dsmall.root");
void CreateBins(double *axis, double min, double max, int Nbins = 10);

void FlowAnalysis_SPEP(std::array<float, 3> dimuonMassRange,
                       std::array<float, 3> dimuonPtRange,
                       std::array<float, 2> dimuonCentRange,
                       std::string FileName = "input_data.txt") {

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
  int NBinsMass = static_cast<int>(dimuonMassRange[2]);
  int NBinsPt = static_cast<int>(dimuonPtRange[2]);
  int NBinsMult = 10;
  double Cent[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  double *MassBins = new double[NBinsMass + 1];
  double *PtBins = new double[NBinsPt + 1];
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