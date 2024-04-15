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
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TChain.h>
#include <TPad.h>
#include <TLegend.h>
#include <TString.h>
#include <TCollection.h>
#include <TKey.h>
#include <TClass.h>
#include <TROOT.h>

#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooCrystalBall.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooWorkspace.h>
#include <RooProdPdf.h>
#include <RooChi2Var.h>
#include "RooStats/SPlot.h"
#include <RooChebychev.h>
#include <RooHist.h>

using namespace RooFit;
using namespace RooStats;

void AddMassModel(RooWorkspace &ws);
void AddData(RooWorkspace &ws, TChain *dChain);
void DoMassFitting(RooWorkspace &ws);
void DoSplot(RooWorkspace &ws);
void DrawPlots(RooWorkspace &ws);
void GetTree(TChain *fChain, std::string FileName="AO2D.root", std::string TreeName="O2rtdimuonall");

void JPsiQAFit(){

	// Initializing workspace
	RooWorkspace wspace{"myWS"};

	// Reading unbinned data from AO2D
	TChain *dataChain = new TChain();
	/*
	TString tree_list[] = {"DF_2302029532114912", "DF_2302029532142432", "DF_2302029540324192", "DF_2302029601531360", "DF_2302029575874464", "DF_2302029590405856", "DF_2302029646758464", "DF_2302029579578944"};
	for(int i=0; i<size(tree_list);i++){
		std::string file_name(tree_list[i]);
		dataChain->Add(Form("AO2D.root?#%s/O2rtdimuonall", file_name.c_str()));
	}
	*/
	GetTree(dataChain);

	cout << "Total entries before cuts: " << dataChain->GetEntries() <<endl;

	// Processing
	AddMassModel(wspace);
	AddData(wspace, dataChain);
	DoMassFitting(wspace);
	//DoSplot(wspace);
	DrawPlots(wspace);

}


void GetTree(TChain *fChain, std::string FileName="AO2D.root", std::string TreeName="O2rtdimuonall"){
	TFile *fInput = TFile::Open(FileName.c_str());
	TIter keyList(fInput->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)keyList())) {
	    TClass *cl = gROOT->GetClass(key->GetClassName());
	    //if (!cl->InheritsFrom("TDirectoryFile")) continue;
	    std::string dir = key->GetName();
	    if(dir == "parentFiles") continue;
	    cout << dir <<endl;
	    fChain->Add(Form("%s?#%s/%s",FileName.c_str(),dir.c_str(),TreeName.c_str()));
  	}

}



void AddMassModel(RooWorkspace &ws){

	cout << "Adding model..." <<endl;

	//====================================================
	// Invariant mass
	//====================================================
	// Initializing free variables 
	RooRealVar x("fMass","M_{#mu^{+}#mu^{-}}",2.0,4.0, "GeV/c2");

	// Initializing double-sided crystalball function for signal
	RooRealVar x0("x0","m_{0}",3.097,2.95,3.15);
	RooRealVar sigmaL("sigmaL","#sigma_{L}",0.05,0.01,0.2);
	RooRealVar sigmaR("sigmaR","#sigma_{R}",0.08,0.05,0.2);
	RooRealVar alphaL("alphaL","#alpha_{L}",0.8,0.05,2.0);
	RooRealVar alphaR("alphaR","#alpha_{R}",0.8,0.05,1.5);
	RooRealVar nL("nL","n_{L}",0.05,0.0005,2.0);
	RooRealVar nR("nR","n_{R}",0.05,0.0005,2.0);
	//RooCrystalBall signal = RooCrystalBall("signal", "double-sided crystalball function", x, x0, sigmaL, alphaL, nL, false);
	RooCrystalBall signal = RooCrystalBall("signal", "double-sided crystalball function", x, x0, sigmaL, alphaL, nL, alphaR, nR);

	// Initializing exponential function for background
	//RooRealVar a0("a0","a_{0}",-0.5,-2,2);
	//RooExponential BKG = RooExponential("BKG","Power background", x, a0);

	// Initializing Chebyshev polynomials for background
	
	// Convenient to set intial values to negative better guiding fit
	RooRealVar a0("a0","a_{0}",-0.01,-1.0,1.0);
	RooRealVar a1("a1","a_{1}",-0.01,-1.0,1.0);
	RooRealVar a2("a2","a_{2}",-0.01,-1.0,1.0);
	//RooRealVar a3("a3","a_{3}",0.1,-2.0,2.0);
	//RooRealVar a4("a4","a_{4}",0.0,-2.0,2.0);
	//RooRealVar a5("a5","a_{5}",0.0,-2.0,2.0);
	//RooRealVar a6("a6","a_{6}",0.0,-2.0,2.0);
	//RooRealVar a7("a7","a_{7}",0.0,-2.0,2.0);
	RooChebychev BKG = RooChebychev("BKG","Chebyshev background",x, {a0,a1,a2});
	




	//====================================================
	// Tauz
	//====================================================
	// Initializing free fvariables 
	RooRealVar Tauz("fTauz","#tau_{z}",-0.004,0.004, "ns");

	// Toy model for signal after sPlot
	RooRealVar mean_tau1("mean_tau1","#mu_{#tau1}",0.0,-0.001,0.001);
	RooRealVar sigma_tau1("sigma_tau1","#sigma_{#tau1}",0.0001,0.0,0.0008);
	RooGaussian signal_tau1("signal_tau1","signal part in #tau_{z}",Tauz, mean_tau1, sigma_tau1);
	mean_tau1.setConstant();

	RooRealVar mean_tau2("mean_tau2","#mu_{#tau2}",0.0,-0.001,0.001);
	RooRealVar sigma_tau2("sigma_tau2","#sigma_{#tau2}",0.001,0.0,0.05);
	RooGaussian signal_tau2("signal_tau2","signal part in #tau_{z}",Tauz, mean_tau2, sigma_tau2);
	mean_tau2.setConstant();


	// Model for prompt J/Psi
	RooRealVar tau0_prompt("tau0_prompt","tau0_prompt",0.0,-0.006,0.006);
	RooRealVar sigmaL_prompt("sigmaL_prompt","sigmaL_prompt",0.005,0,0.01);
	RooRealVar sigmaR_prompt("sigmaR_prompt","sigmaR_prompt",0.005,0,0.01);
	RooRealVar alphaL_prompt("alphaL_prompt","alphaL_prompt",1,0,10);
	RooRealVar nL_prompt("nL_prompt","nL_prompt",1,0,10);
	RooRealVar alphaR_prompt("alphaR_prompt","alphaR_prompt",1,0,10);
	RooRealVar nR_prompt("nR_prompt","nR_prompt",1,0,10);
	RooCrystalBall prompt_signal = RooCrystalBall("prompt_signal", "double-sided crystalball function", Tauz, tau0_prompt, sigmaL_prompt, sigmaR_prompt, alphaL_prompt, nL_prompt, alphaR_prompt, nR_prompt);

	// Model for background
	RooRealVar tau0_BKG("tau0_BKG","tau0_BKG",0.0,-0.01,0.01);
	RooRealVar sigmaL_BKG("sigmaL_BKG","sigmaL_BKG",0.005,0,0.01);
	RooRealVar sigmaR_BKG("sigmaR_BKG","sigmaR_BKG",0.005,0,0.01);
	RooRealVar alphaL_BKG("alphaL_BKG","alphaL_BKG",1,0,10);
	RooRealVar nL_BKG("nL_BKG","nL_BKG",1.,0.5,10);
	RooRealVar alphaR_BKG("alphaR_BKG","alphaR_BKG",1,0,10);
	RooRealVar nR_BKG("nR_BKG","nR_BKG",1.,0.5,10);
	RooCrystalBall BKG_tau = RooCrystalBall("BKG_tau", "double-sided crystalball function", Tauz, tau0_BKG, sigmaL_BKG, sigmaR_BKG, alphaL_BKG, nL_BKG, alphaR_BKG, nR_BKG);



	//====================================================
	// Combinined models
	//====================================================

	// Yields for fitting
	RooRealVar sigYield("sigYield", "N_{sig}", 50, 0., 100000000);
   	RooRealVar bkgYield("bkgYield", "N_{bkg}", 100, 0., 1000000000);
   	RooRealVar fGauss1("fGauss1","fGauss1",0.,0.,10000);
   	RooRealVar fGauss2("fGauss2","fGauss2",0.,0.,10000);
   	// Combining
	//RooAddPdf tau_model("tau_model", "combined tau model", {sigModel, bkgModel}, {sigYield, bkgYield});
	RooAddPdf model("model", "combined mass model", {signal, BKG}, {sigYield, bkgYield});
	RooAddPdf double_gauss("double_gauss","double gaussian function", {signal_tau1,signal_tau2},{fGauss1,fGauss2});

	//ws.import(tau_model);
	ws.import(model);
	ws.import(double_gauss);
	ws.Print();


}



void AddData(RooWorkspace &ws, TChain *dChain){
	cout << "Adding data to workspace..." <<endl;
	// Variables
	RooRealVar* x = ws.var("fMass");
	RooRealVar* Tauz = ws.var("fTauz");
	//RooRealVar Tauz("fTauz","#tau_{z}",-0.004,0.004,"ns");
	//RooRealVar Eta("fEta","#eta",-3.5,-2.6);

	// Single muon selection
	RooRealVar Eta1("fEta1","#eta_{1}",-4.0,-2.5);
	RooRealVar Eta2("fEta2","#eta_{2}",-4.0,-2.5);
	//RooRealVar Chi2MCHMFT1("fChi2MatchMCHMFT1","fChi2MatchMCHMFT1",0.0,1000.0);
	//RooRealVar Chi2MCHMFT2("fChi2MatchMCHMFT2","fChi2MatchMCHMFT2",0.0,1000.0);
	RooRealVar Chi2MCHMID1("fChi2MatchMCHMID1","fChi2MatchMCHMID1",0.0,4.0);
	RooRealVar Chi2MCHMID2("fChi2MatchMCHMID2","fChi2MatchMCHMID2",0.0,4.0);
	RooRealVar pT1("fPt1", "fPt1",0.0,10.0,"GeV/c");
	RooRealVar pT2("fPt2", "fPt2",0.0,10.0,"GeV/c");

	// Dimuon selection
	RooRealVar Sign("fSign","fSign",-3,3);
	//RooRealVar Chi2pca("fChi2pca","fChi2pca",0.0, 4.0);
	RooRealVar pT("fPt", "fPt",2.0,10.0,"GeV/c");
	//RooRealVar Eta("fEta","#eta",-3.5,-2.6);

	// Initializing data set with selected cuts
	//std::string Filter = "fPt1>0.5 && fPt2>0.5 && fSign==0 && fChi2MatchMCHMFT1<40 && fChi2MatchMCHMFT2<40";
	std::string Filter = "fPt1>1.0 && fPt2>1.0 && fSign==0";
	RooDataSet data("data", "unbinned dataset", dChain, RooArgSet(*x, *Tauz, pT1, pT2, Eta1, Eta2, Sign, pT, Chi2MCHMID1, Chi2MCHMID2), Filter.c_str());
	data.Print();
	ws.import(data);
}



void DoSplot(RooWorkspace &ws){

	RooAbsPdf *model = ws.pdf("model");
	//RooAbsPdf *tau_model = ws.pdf("tau_model");
	RooRealVar *sigYield = ws.var("sigYield");
   	RooRealVar *bkgYield = ws.var("bkgYield");
   	RooDataSet& data = static_cast<RooDataSet&>(*ws.data("data"));


   	// Fixing parameters

   	// Fixing parameters from invariant mass fitting
   	
   	RooRealVar *x0 = ws.var("x0");
	RooRealVar *sigmaL = ws.var("sigmaL");
	RooRealVar *alphaL = ws.var("alphaL");
	RooRealVar *nL = ws.var("nL");
	x0->setConstant();
	sigmaL->setConstant();
	alphaL->setConstant();
	nL->setConstant();

	RooRealVar *a0 = ws.var("a0");
	a0->setConstant();


	/*
	RooRealVar *tau0_prompt = ws.var("tau0_prompt");
	RooRealVar *sigmaL_prompt = ws.var("sigmaL_prompt");
	RooRealVar *sigmaR_prompt = ws.var("sigmaR_prompt");
	RooRealVar *alphaL_prompt = ws.var("alphaL_prompt");
	RooRealVar *nL_prompt = ws.var("nL_prompt");
	RooRealVar *alphaR_prompt = ws.var("alphaR_prompt");
	RooRealVar *nR_prompt = ws.var("nR_prompt");
	tau0_prompt->setConstant();
	sigmaL_prompt->setConstant();
	sigmaR_prompt->setConstant();
	alphaL_prompt->setConstant();
	nL_prompt->setConstant();
	alphaR_prompt->setConstant();
	nR_prompt->setConstant();
	*/

   	//RooMsgService::instance().setSilentMode(true);
   	cout << "\n\n------------------------------------------\nThe dataset before creating sWeights:\n";
   	data.Print();

   	// Fitting with sPlot
   	//RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

   	RooStats::SPlot sData{"sData", "An SPlot", data, model, RooArgList(*sigYield, *bkgYield)};


   	std::cout << "\n\n------------------------------------------\nThe dataset after creating sWeights:\n";
    data.Print();
    cout << "Yields check after sPlot:" <<endl;
    cout << "Signal yield from mass fitting: "<< sigYield->getVal() <<endl;
    cout << "Signal yield from sWeight: " << sData.GetYieldFromSWeight("sigYield") <<endl;
    cout << "Background yield from mass fitting: "<< bkgYield->getVal() <<endl;
    cout << "Background yield from sWeight: " << sData.GetYieldFromSWeight("bkgYield") <<endl;

    ws.import(data, Rename("dataWithSWeights"));
    //RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);

}




void DoMassFitting(RooWorkspace &ws){
	cout << "Fitting invariant mass..."<<endl;
	RooAbsPdf *model = ws.pdf("model");
	RooDataSet& data = static_cast<RooDataSet&>(*ws.data("data"));
	//model->fitTo(data, Save(kTRUE),Extended(kTRUE),Minos(kTRUE),Strategy(2));
	model->fitTo(data,Extended(kTRUE));

}



void DrawPlots(RooWorkspace &ws){
	cout << "Plotting results" <<endl;
	// Initializing plot canvas
	TCanvas *c1 = new TCanvas("FitInvMass","Fit result for invariant mass");
	//TCanvas *c2 = new TCanvas("FitTauz","Fit results for tauz");
	TCanvas *c3 = new TCanvas("Resi&Pull", "Residu and pulls plot for invarian mass fit", 1200, 600);
	c3->Divide(2);

	// Get ingredients from workspace
	//RooAbsPdf *tau_model = ws.pdf("tau_model");
	RooAbsPdf *model = ws.pdf("model");
	RooAbsPdf *gauss_tau = ws.pdf("double_gauss");
	//RooAbsPdf *sigModel = ws.pdf("sigModel");
	//RooAbsPdf *bkgModel = ws.pdf("bkgModel");
	RooAbsPdf *signal = ws.pdf("signal");
	RooAbsPdf *BKG = ws.pdf("BKG");

	RooRealVar *x = ws.var("fMass");
	RooRealVar *Tauz = ws.var("fTauz");
	//Tauz->setRange("left",-0.003,0.0000001);
	Tauz->setRange("full",-0.004,0.004);
	//RooRealVar *Tauz = ws.var("fTauz");
	RooRealVar *sigYield = ws.var("sigYield");
	RooRealVar *bkgYield = ws.var("bkgYield");

	//RooDataSet& data = static_cast<RooDataSet&>(*ws.data("dataWithSWeights"));
	RooDataSet& data = static_cast<RooDataSet&>(*ws.data("data"));

	//RooDataSet dataw_sig{data.GetName(), data.GetTitle(), &data, *data.get(), nullptr, "sigYield_sw"};
    //RooDataSet dataw_bkg{data.GetName(), data.GetTitle(), &data, *data.get(), nullptr, "bkgYield_sw"};

    // Fitting tau with sWeighted data;
    //gauss_tau->fitTo(dataw_sig,Extended(kTRUE),Range("full"));
    //double renorm_scale = dataw_sig.sumEntries("fTauz>-0.003&&fTauz<0.0000001"); // normalization w.r.t left


    // Chi-square extraction
    int nBins = 30;
    int ndf = 11;
    TH1* hist_sdata = data.createHistogram("fMass", nBins);
    RooDataHist* h_sdata = new RooDataHist("h_sdata", "h_sdata", RooArgSet(*x), hist_sdata);
    RooChi2Var* chi2Var = new RooChi2Var("chi2","chi2", *model, *h_sdata, DataError(RooAbsData::Poisson));
    double chi2M = chi2Var->getVal();



	// Initializing data frame for plots
	gStyle->SetLegendFont(22);
   	//gStyle->SetLegendTextSize(0.1);
	c1->cd();
	RooPlot *xframe = x->frame(nBins=50);
	xframe->SetTitle("M_{#mu#mu} fitted with unbinned data");
	data.plotOn(xframe);
   	model->plotOn(xframe, Name("FullModel"));
   	model->plotOn(xframe, Components(*signal), LineStyle(kDashed), LineColor(kRed), Name("signal_Model"));
   	model->plotOn(xframe, Components(*BKG), LineStyle(kDashed), LineColor(kGreen), Name("bkg_Model"));
   	model->paramOn(xframe,Label(Form("#chi^{2}/ndf = %f",chi2M/ndf)),ShowConstants(true),Layout(0.63,0.9,0.9),Format("TE",AutoPrecision(3)));
   	xframe->getAttText()->SetTextSize(0.025);
   	xframe->getAttText()->SetTextFont(22);
   	xframe->Draw();
   	TLegend* lgd1 = new TLegend(0.15, 0.7, 0.4, 0.85);
    lgd1->AddEntry(xframe->findObject("signal_Model"), "J/#Psi", "L");
    lgd1->AddEntry(xframe->findObject("bkg_Model"), "Bkg", "L");
    lgd1->AddEntry(xframe->findObject("FullModel"), "Full Model", "L");
    lgd1->Draw();

    RooPlot *ResiFrame = x->frame(Title("Residu Distribution"));
    RooPlot *PullFrame = x->frame(Title("Pull Distribution"));
    RooHist *hResidu = xframe->residHist();
    RooHist *hPull = xframe->pullHist();

    ResiFrame->addPlotable(hResidu,"P");
    PullFrame->addPlotable(hPull,"P");
    c3->cd(1);
    gPad->SetLeftMargin(0.15);
    ResiFrame->Draw();
    c3->cd(2);
    gPad->SetLeftMargin(0.15);
    PullFrame->Draw();



   	/*
   	c2->cd();
	RooPlot *tauframe = Tauz->frame(nBins=100);
	tauframe->SetAxisRange(2,300,"Y");
	tauframe->SetTitle("#tau_{z}");
	cout << "Hello" <<endl;
	dataw_sig.plotOn(tauframe, LineColor(kBlue),Name("dataw_sig"));
	tauframe->getAttMarker()->SetMarkerColor(kBlue);
	dataw_bkg.plotOn(tauframe, LineColor(kGreen),Name("dataw_bkg"));
	cout << "Hello" <<endl;
	tauframe->getAttMarker()->SetMarkerColor(kGreen);
	data.plotOn(tauframe);
	tauframe->getAttMarker()->SetMarkerColor(kBlack);
	gauss_tau->plotOn(tauframe,LineStyle(kDashed),LineColor(kRed),Normalization(renorm_scale,2),Range("full"));
	gauss_tau->paramOn(tauframe,Format("TE"),Layout(0.59,0.9,0.9));
	tauframe->getAttText()->SetTextSize(0.028);
   	tauframe->getAttText()->SetTextFont(22);
	tauframe->Draw();
	TLegend* lgd2 = new TLegend(0.15, 0.7, 0.3, 0.85);
    lgd2->AddEntry(tauframe->findObject("dataw_sig"), "J/#Psi", "L");
    lgd2->AddEntry(tauframe->findObject("dataw_bkg"), "Bkg", "L");
    lgd2->Draw();
    */


}








