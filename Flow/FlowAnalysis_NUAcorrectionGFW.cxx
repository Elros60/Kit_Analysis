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
#include <TH3D.h>
#include <THStack.h>
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

#include "Framework/Logger.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

using namespace std;

vector<string> tokenize(string input_string);

void FlowAnalysis_NUAcorrections(
    string FileName = "input_AnalysisResults.txt") {

  auto hs_stack = new THStack("Run_by_Run_NUA_correction", "");
  TList *ls = new TList();
  fstream InputFiles;
  InputFiles.open(FileName, ios::in);
  if (InputFiles.is_open()) {
    string File;
    cout << "Start reading input AnalysisResults.root ..." << endl;
    int count = 0;
    while (getline(InputFiles, File)) {
      cout << "Reading input from: " << File << endl;
      TFile *f = TFile::Open(File.c_str());
      vector<string> File_string = tokenize(File);
      string run_number = File_string[File_string.size() - 2];
      GFWWeights *weights = (GFWWeights *)f->Get("d-q-event-qvector/weights");
      TObjArray *ls_weights = weights->GetDataArray();
      TH3D *hs = reinterpret_cast<TH3D *>(ls_weights->At(0));
      TH1D *hs_x = hs->ProjectionX();
      hs_x->Scale(1. / hs_x->Integral());
      hs_x->SetName(run_number.c_str());
      hs_x->SetTitle(run_number.c_str());
      hs_x->SetLineColor(count + 1);
      hs_x->SetLineWidth(3);
      ls->Add(hs_x);
      hs_stack->Add(hs_x);
      count++;
    }
  }
  InputFiles.close();

  TCanvas *c = new TCanvas("NUA_Corrections_PbPb_Pass3", "NUA_Corrections");
  c->cd();
  hs_stack->Draw("HIST nostack");
  hs_stack->GetXaxis()->SetTitle("#phi (Rad)");
  hs_stack->GetYaxis()->SetTitle("Normalized counts");

  ls->Add(c);

  TFile fout("AnalysisResults_NUAcorrections.root", "RECREATE");
  ls->Write("RunByRun_histograms", TObject::kSingleKey);
  fout.Close();
}

vector<string> tokenize(string input_string) {
  vector<string> output_string;
  string temp;
  istringstream ss(input_string);
  while (getline(ss, temp, '/')) {
    output_string.push_back(temp);
  }
  return output_string;
}