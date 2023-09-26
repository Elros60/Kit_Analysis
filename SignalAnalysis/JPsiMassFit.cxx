#include <cmath>
#include <string>
#include <iostream>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TChain.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooCrystalBall.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooExponential.h>

using namespace RooFit;

void JPsiMassFit(){

	// Loading data set to be fitted from analysis results
	/*
	// Read data from binned histogram(only standard cuts applied)
	TFile *input = TFile::Open("AnalysisResults.root");
	THashList *hlist = (THashList*)input->Get("analysis-same-event-pairing/output");
	TList *flist = (TList*)hlist->FindObject("PairsMuonSEPM_matchedQualityCuts");
	TH1F *histo = (TH1F*)flist->FindObject("Mass");
	*/

	
	// Read unbinned data from AO2D
	TChain dataChain;
	TString tree_list[] = {"DF_2302029532114912", "DF_2302029532142432", "DF_2302029540324192", "DF_2302029601531360", "DF_2302029575874464", "DF_2302029590405856", "DF_2302029646758464", "DF_2302029579578944"};
	for(int i=0; i<size(tree_list);i++){
		std::string file_name(tree_list[i]);
		dataChain.Add(Form("AO2D.root?#%s/O2rtdimuonall", file_name.c_str()));
	}
	cout << "Total entries before cuts: " << dataChain.GetEntries() <<endl;
	

	
	RooRealVar x("fMass","Mass",2.6,3.5, "GeV/c2");
	RooRealVar Eta("fEta","Eta",-3.5,-2.6);
	RooRealVar Chi2MCHMFT1("fChi2MatchMCHMFT1","fChi2MatchMCHMFT1",0.0,1000.0);
	RooRealVar Chi2MCHMFT2("fChi2MatchMCHMFT2","fChi2MatchMCHMFT2",0.0,1000.0);
	RooRealVar Sign("fSign","fSign",-3,3);
	RooRealVar pT("fPt", "fPt",0.0,100);

	//RooDataHist data("data", "binned J/Psi invariant mass", x, histo);
	std::string Filter = "fPt >=0.5 && fSign==0 && fChi2MatchMCHMFT1<=40 && fChi2MatchMCHMFT2<=40";
	RooDataSet data("data", "unbinned dataset", &dataChain, RooArgSet(x, pT, Sign, Chi2MCHMFT1, Chi2MCHMFT2), Filter.c_str());
	

	
	// Initializing double-sided crystalball function
	RooRealVar x0("x0","mass0",3.,2.6,3.5);
	RooRealVar sigmaL("sigma","sigmaL",0.07,0,0.5);
	//RooRealVar sigmaR("sigmaR","sigmaR",0.07,0,1);
	RooRealVar alphaL("alpha","alphaL",1,0,10);
	RooRealVar nL("n","nL",1,0,10);
	//RooRealVar alphaR("alphaR","alphaR",1,0,10);
	//RooRealVar nR("nR","nR",1,0,10);
	RooCrystalBall signal = RooCrystalBall("signal", "double-sided crystalball function", x, x0, sigmaL, alphaL, nL, false);
	// RooCrystalBall signal = RooCrystalBall("signal", "double-sided crystalball function", x, x0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

	// Initializing background model using power function
	RooRealVar a0("a0","a0",-1,-5,5);
	RooExponential BKG = RooExponential("BKG","Power background", x, a0);

	// Signal+BKG model composition
	RooRealVar coeff1("coeff_{sig}","signal coeff",100,-100000,100000);
	RooRealVar coeff2("coeff_{bkg}","BKG coeff",100,-100000,100000);
	//fsig.setConstant(kTRUE);
	RooAddPdf model("model", "model", RooArgList(signal, BKG), RooArgList(coeff1,coeff2));

	// Fitting
	model.fitTo(data);
	RooPlot *xframe = x.frame();
	xframe->SetTitle("M_{#mu#mu} fitted with unbinned tree using LHC23K_pass1_small");
	data.plotOn(xframe);
	model.plotOn(xframe, Components(signal),LineStyle(kDashed),LineColor(kRed));
	model.plotOn(xframe, Components(BKG),LineStyle(kDashed),LineColor(kBlue));
	model.plotOn(xframe, LineColor(kBlue));
	
	model.paramOn(xframe);

	xframe->Draw();
	cout << endl;
	cout << xframe->chiSquare(7) << endl;

}