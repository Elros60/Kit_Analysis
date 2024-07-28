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
#include <TTreeIndex.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

void LoadData(TChain *fChain, std::string TreeName = "O2rerefflow",
              std::string FileName = "AO2D.root");

void MergeTree(std::string FileName = "AO2D.root",
               std::string TreeName = "O2rtdimuonall") {
  TChain *file = new TChain(TreeName.c_str());
  LoadData(file, TreeName, FileName);
  file->Merge(Form("AO2D_merged_%s.root", TreeName.c_str()));
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