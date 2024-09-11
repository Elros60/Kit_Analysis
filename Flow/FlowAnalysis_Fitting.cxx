#include "FlowAnalysis_Fitting.h"

//______________________________________________________________________________
void FlowAnalysis_Fitting::init() {
  mflag_sig = 0;
  mflag_bkg = 1;
  mchi2max = 2.0;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::setModel(int flag_sig, int flag_bkg) {
  mflag_sig = flag_sig;
  mflag_bkg = flag_bkg;
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::CreateModel(RooRealVar x, int flag) {
  if (flag == CB2) {
    RooRealVar m0("mean", "mean", 3.097, 3.09, 3.1);
    RooRealVar sigma("sigma", "sigma", 0.08, 0.05, 0.12);
    RooRealVar alphaL("alphaL", "alphaL", 0.883, 0.1, 3.0);
    RooRealVar alphaR("alphaR", "alphaR", 1.832, 0.1, 3.0);
    RooRealVar nL("nL", "nL", 9.940, 1.0, 7.0);
    RooRealVar nR("nR", "nR", 15.323, 2.0, 7.0);
    alphaL.setConstant(kTRUE);
    alphaR.setConstant(kTRUE);
    nL.setConstant(kTRUE);
    nR.setConstant(kTRUE);
    RooCrystalBall model =
        RooCrystalBall("Signal", "CB2", x, m0, sigma, alphaL, nL, alphaR, nR);
    mWS.import(model);
  }
  if (flag == VWG) {
    RooRealVar a0("a0", "a0", 0.1, -1.0, 5.0);
    RooRealVar a1("a1", "a1", 0.2, 0.0, 3.0);
    RooRealVar a2("a2", "a2", 0.1, -1.0, 1.0);
    a1.setConstant(kTRUE);
    a2.setConstant(kTRUE);

    RooFormulaVar sigma_VWG("sigma_VWG", "sigma_VWG", "a1+a2*((m-a0)/a0)",
                            RooArgList(x, a0, a1, a2));
    RooGenericPdf model = RooGenericPdf(
        "Bkg", "VWG", "exp(-(m-a0)*(m-a0)/(2.*sigma_VWG*sigma_VWG))",
        RooArgSet(x, a0, sigma_VWG));
    mWS.import(model);
  }
  if (flag == POL) {
    RooRealVar a0("a0", "a0", -0.2, -10.0, 10.0);
    RooRealVar a1("a1", "a1", 0.0, -10.0, 10.0);
    RooRealVar a2("a2", "a2", 0.0, -10.0, 10.0);
    RooRealVar a3("a3", "a3", 0.0, -2.0, 2.0);
    RooRealVar a4("a4", "a4", 0.0, -1.0, 1.0);
    RooRealVar a5("a5", "a5", 0.0, -0.5, 0.5);
    RooRealVar a6("a6", "a6", 0.0, -0.1, 0.1);
    RooRealVar a7("a7", "a7", 0.0, -0.05, 0.05);
    a1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    a3.setConstant(kTRUE);
    a4.setConstant(kTRUE);
    a5.setConstant(kTRUE);
    a6.setConstant(kTRUE);
    a7.setConstant(kTRUE);
    RooArgList param_bkg = RooArgList(a0, a1, a2, a3, a4, a5, a6, a7);
    RooPolynomial model = RooPolynomial("Bkg", "Pol", x, param_bkg);
    mWS.import(model);
  }
  if (flag == PolExp) {
    RooRealVar a0("a0", "a0", 1.0, -20.0, 500.0);
    RooRealVar a1("a1", "a1", 1.0, -10.0, 10.0);
    RooRealVar a2("a2", "a2", -1.0, -10.0, 10.0);
    RooRealVar a3("a3", "a3", 0.0, -10.0, 10.0);
    RooRealVar a4("a4", "a4", 0.0, -10.0, 10.0);
    RooRealVar a5("a5", "a5", 0.0, -10.0, 10.0);
    a2.setConstant(kTRUE);
    a3.setConstant(kTRUE);
    a4.setConstant(kTRUE);
    a5.setConstant(kTRUE);
    RooGenericPdf model = RooGenericPdf(
        "Bkg", "PolExp", "(a0+a1*m+a3*m*m+a4*m*m*m+a5*m*m*m*m)*exp(a2*m)",
        RooArgSet(x, a0, a1, a2, a3, a4, a5));
    mWS.import(model);
  }
  if (flag == Exp2) {
    RooRealVar a0("a0", "a0", 0.0, -10.0, 10.0);
    RooRealVar a1("a1", "a1", -0.1, -5.0, 5.0);
    RooRealVar a2("a2", "a2", 0.1, -5.0, 5.0);
    RooRealVar a3("a3", "a3", 0.0, -10.0, 10.0);
    RooRealVar a4("a4", "a4", -0.2, -5.0, 5.0);
    RooRealVar a5("a5", "a5", 0.2, -5.0, 5.0);
    a1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    a3.setConstant(kTRUE);
    a4.setConstant(kTRUE);
    a5.setConstant(kTRUE);
    RooGenericPdf model =
        RooGenericPdf("Bkg", "Exp2", "a0*exp(a1*(m-a2))+a3*exp(a4*(m-a5))",
                      RooArgSet(x, a0, a1, a2, a3, a4, a5));
    mWS.import(model);
  }
  if (flag == Chebychev) {
    RooRealVar a0("a0", "a0", 0.0, -5.0, 5.0);
    RooRealVar a1("a1", "a1", 0.0, -2.0, 2.0);
    RooRealVar a2("a2", "a2", 0.0, -2.0, 2.0);
    RooRealVar a3("a3", "a3", 0.0, -2.0, 2.0);
    RooRealVar a4("a4", "a4", 0.0, -2.0, 2.0);
    RooRealVar a5("a5", "a5", 0.0, -2.0, 2.0);
    RooRealVar a6("a6", "a6", 0.0, -2.0, 2.0);
    RooRealVar a7("a7", "a7", 0.0, -2.0, 2.0);
    a1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    a3.setConstant(kTRUE);
    a4.setConstant(kTRUE);
    a5.setConstant(kTRUE);
    a6.setConstant(kTRUE);
    a7.setConstant(kTRUE);
    RooArgList param_bkg = RooArgList(a0, a1, a2, a3, a4, a5, a6, a7);
    RooChebychev model = RooChebychev("Bkg", "Chebychev", x, param_bkg);
    mWS.import(model);
  }
}

//______________________________________________________________________________
void FlowAnalysis_Fitting::runFitting(TH1D *hs, TList *ls, double ptmin,
                                      double ptmax, double massmin,
                                      double massmax) {

  // Setting up RooFit
  cout << "Start processing Pt range: " << ptmin << " " << ptmax << endl;
  RooRealVar m("m", "m_{#mu#mu}", massmin, massmax, "GeV/c2");
  RooDataHist dh("dh", "dh", m, Import(*hs));
  RooPlot *frame = m.frame(
      Name(Form("mFrame_%g_%g", ptmin, ptmax)),
      Title(Form("m_{#mu#mu} signal fit (%g < pT < %g)", ptmin, ptmax)));
  dh.plotOn(frame, DataError(RooAbsData::Poisson), Name("data"));

  // Setting up fit model
  CreateModel(m, ModelType(mflag_sig)); // Create model for signal
  CreateModel(m, ModelType(mflag_bkg)); // Create model for background
  auto sig = mWS.pdf("Signal");
  auto bkg = mWS.pdf("Bkg");
  RooArgList param_bkg = RooArgList(*(bkg->getParameters(m)));

  /// Construct a combined signal+background model
  RooRealVar nsig("nsig", "#signal events", 500, 0., 100000000);
  RooRealVar nbkg("nbkg", "#background events", 500000, 0., 100000000);
  RooAddPdf model("model", "model", {*sig, *bkg}, {nsig, nbkg});

  // Do fitting
  /// Iterative fitting
  RooArgList param_bkg_const =
      RooArgList(*param_bkg.selectByAttrib("Constant", kTRUE));

  for (int i = 0; i <= param_bkg_const.getSize(); i++) {
    model.fitTo(dh, Hesse(0));
    RooAbsReal *chi2var = model.createChi2(dh, DataError(RooAbsData::Poisson));
    int nfree =
        model.getParameters(dh)->selectByAttrib("Constant", kFALSE)->getSize();
    double chi2ndf = chi2var->getVal() / (dh.numEntries() - nfree);
    if (chi2ndf > mchi2max && i != param_bkg_const.getSize()) {
      static_cast<RooRealVar *>(param_bkg_const.at(i))->setConstant(kFALSE);
    } else {
      break;
    }
  }

  // Plotting
  model.plotOn(frame, LineColor(kBlue), LineWidth(5.0), Name("model"));
  model.plotOn(frame, Components(*sig), LineColor(kRed), LineWidth(5.0),
               Name("Signal"));
  model.plotOn(frame, Components(*bkg), LineStyle(ELineStyle::kDashed),
               LineColor(kBlue), LineWidth(5.0), Name("Bkg"));

  int nfree =
      model.getParameters(dh)->selectByAttrib("Constant", kFALSE)->getSize();
  double chi2 = frame->chiSquare("model", "data", nfree);
  string chi2text = Form("#chi^{2}/ndf = %f ", chi2);

  double sigma_fitted =
      model.getParameters(dh)->selectByName("sigma")->getRealValue("sigma");
  double mean_fitted =
      model.getParameters(dh)->selectByName("mean")->getRealValue("mean");
  RooRealVar m_3sigma("m3sig", "m_{#mu#mu}", mean_fitted - 3. * sigma_fitted,
                      mean_fitted + 3. * sigma_fitted, "GeV/c2");
  m.setRange("3sigma", mean_fitted - 3. * sigma_fitted,
             mean_fitted + 3. * sigma_fitted);
  double nsig_3sigma =
      nsig.getVal() * sig->createIntegral(m, m, "3sigma")->getVal();
  double nbkg_3sigma =
      nbkg.getVal() * bkg->createIntegral(m, m, "3sigma")->getVal();
  double SNR = nsig_3sigma / nbkg_3sigma;
  string SNRtext = Form("(S/B)_{3#sigma} = %f ", SNR);
  model.paramOn(frame, Label(Form("%s\n%s", chi2text.c_str(), SNRtext.c_str())),
                Layout(0.63, 0.9, 0.9), Format("NE", FixedPrecision(5)));

  // Get residual and pull histograms
  RooHist *hresid = frame->residHist("data", "model");
  RooHist *hpull = frame->pullHist("data", "model");
  RooPlot *frame_resid = m.frame(Title("Residual Distribution"));
  RooPlot *frame_pull = m.frame(Title("Pull Distribution"));
  frame_resid->addPlotable(hresid, "P");
  frame_pull->addPlotable(hpull, "P");

  // Saving plot
  TCanvas *c = new TCanvas(Form("Fitted_%s", hs->GetName()),
                           Form("Fitted_%s", hs->GetName()));
  TCanvas *c_resid = new TCanvas(Form("Residual_%s", hs->GetName()),
                                 Form("Residual_%s", hs->GetName()));
  TCanvas *c_pull = new TCanvas(Form("Pull_%s", hs->GetName()),
                                Form("Pull_%s", hs->GetName()));
  c->cd();
  frame->Draw();
  ls->Add(c);

  c_resid->cd();
  frame_resid->Draw();
  ls->Add(c_resid);

  c_pull->cd();
  frame_pull->Draw();
  ls->Add(c_pull);
}