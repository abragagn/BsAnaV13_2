#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "Math/MinimizerOptions.h"

#include "RooPlot.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussModel.h"
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooGamma.h"
#include "RooTruthModel.h"

#include "RooJohnsonLocal.cxx"

using namespace RooFit;
using namespace std;

TString datafile = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBdData2018.root";
TString mcfile   = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBdMC2018.root";
double ctau_PDG = 0.04554; // from PDG 2019
double ctau_PDG_err = 0.00012;
double ctau_mc = 0.045718350;

void ct2Dfit()
{
    // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    // Too many warnings
    auto& msg = RooMsgService::instance();
    for(int iStream = 0; iStream < msg.numStreams(); iStream++) {
        auto& stream = msg.getStream(iStream);
        if (stream.minLevel <= RooFit::WARNING) {
            stream.removeTopic(RooFit::Eval);
        }
    }

    auto fMC = new TFile(mcfile);
    auto tMC = (TTree*)fMC->Get("treeFit");
    int neMC = tMC->GetEntries();

    auto fDATA = new TFile(datafile);
    auto tDATA = (TTree*)fDATA->Get("treeFit");
    int neDATA = tDATA->GetEntries();

    cout<<"neDATA "<<neDATA<<endl;
    cout<<"neMC "<<neMC<<endl;

    vector<double> ctBins = {0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025,
                           0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100,
                           0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.350, 0.400, 0.500};

    int  nCtBins = ctBins.size() - 1;

    auto hct_MC  = new TH1D("hct_MC","ct MC; ct [cm]", nCtBins, ctBins.data());
    auto hct_GEN = new TH1D("hct_GEN","ct GEN from expo(ctau_PDG); ct [cm]", nCtBins, ctBins.data());
    auto resobulk = new TH1D("resobulk", "ct Resolution;ct_gen - ct_reco [cm]",  100, -0.01, 0.01);
    auto pullbulk = new TH1D("pullbulk", "ct Pulls; pull",  50, -5.0, 5.0);
    auto hct_DATA = new TH1D("hct_DATA","ct DATA; ct [cm]", 100, 0.007, 0.4);

    auto rnd = new TRandom3(36);
    
    float ctMC, ctMCerr, ctMCgen, ctDATA;

    tMC->SetBranchAddress("BsCt2DMC",&ctMC);
    tMC->SetBranchAddress("BsCt2DMCErr",&ctMCerr);
    tMC->SetBranchAddress("BsCt2DMC_GEN",&ctMCgen);
    tDATA->SetBranchAddress("BsCt2DMC",&ctDATA);

    // Fill mc
    for (int i=0;i<neMC;i++){
        tMC->GetEntry(i);
        hct_MC->Fill(ctMC);
        resobulk->Fill(ctMCgen-ctMC);
        pullbulk->Fill((ctMCgen-ctMC)/ctMCerr);
    }

    float kappa_val = pullbulk->GetRMS();
    cout<<" --- kappa value : "<<kappa_val<<endl;


    // Fill data
    for(int i=0;i<neDATA;i++){
        tDATA->GetEntry(i);
        hct_DATA->Fill(ctDATA);
    }

    // Fill gen
    for(int i=0;i<1e+07;i++){
        double randct = rnd->Exp(ctau_mc);
        double randreso = resobulk->GetRandom();
        hct_GEN->Fill(randct + randreso);
    }

    // Make eff
    auto hEFF = (TH1D*)hct_MC->Clone("hEFF");
    hEFF->Sumw2();
    hEFF->Divide(hct_GEN);
    
    auto modelEFF = new TF1("modelEFF","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.4);
    modelEFF->SetParameter(1,-1);
    modelEFF->FixParameter(2,1);

    auto c1 = new TCanvas();

    TFitResultPtr fitResEFF = hEFF->Fit("modelEFF","MS");
    TMatrixDSym covEFF = fitResEFF->GetCovarianceMatrix();
    auto fittedEFF = hEFF->GetFunction("modelEFF");
    hEFF->Draw();

    // ROOFIT

    // Eff model parameters
    RooRealVar ctp0("ctp0", "ctp0", fittedEFF->GetParameter(0));
    ctp0.setError(fittedEFF->GetParError(0));
    RooRealVar ctp1("ctp1", "ctp1", fittedEFF->GetParameter(1));
    ctp1.setError(fittedEFF->GetParError(1));
    RooRealVar ctp2("ctp2", "ctp2", fittedEFF->GetParameter(2));
    ctp2.setError(fittedEFF->GetParError(2));
    RooRealVar ctp3("ctp3", "ctp3", fittedEFF->GetParameter(3));
    ctp3.setError(fittedEFF->GetParError(3));
    RooRealVar ctp4("ctp4", "ctp4", fittedEFF->GetParameter(4));
    ctp4.setError(fittedEFF->GetParError(4));
    RooRealVar ctp5("ctp5", "ctp5", fittedEFF->GetParameter(5));
    ctp5.setError(fittedEFF->GetParError(5));
    RooRealVar ctp6("ctp6", "ctp6", fittedEFF->GetParameter(6));
    ctp6.setError(fittedEFF->GetParError(6));

    //Bs 2017
    // RooRealVar ctp0("ctp0", "ctp0", -0.051);
    // ctp0.setError(0.02);
    // RooRealVar ctp1("ctp1", "ctp1", -11.473);
    // ctp1.setError(0.09);
    // RooRealVar ctp2("ctp2", "ctp2", 1.);
    // RooRealVar ctp3("ctp3", "ctp3", 15.72);
    // ctp3.setError(0.04);
    // RooRealVar ctp4("ctp4", "ctp4", 1.060;
    // ctp4.setError(0.002);
    // RooRealVar ctp5("ctp5", "ctp5", 4.96);
    // ctp5.setError(0.01);
    // RooRealVar ctp6("ctp6", "ctp6", 40.126);
    // ctp6.setError(0.0002);

    //Bs2018
    // RooRealVar ctp0("ctp0", "ctp0", -2.528);
    // ctp0.setError(0.02);
    // RooRealVar ctp1("ctp1", "ctp1", -14.018);
    // ctp1.setError(0.07);
    // RooRealVar ctp2("ctp2", "ctp2", 1.);
    // RooRealVar ctp3("ctp3", "ctp3", 1022.5);
    // ctp3.setError(0.08);
    // RooRealVar ctp4("ctp4", "ctp4", 4.);
    // ctp4.setError(0.03);
    // RooRealVar ctp5("ctp5", "ctp5", 332.8);
    // ctp5.setError(0.03);
    // RooRealVar ctp6("ctp6", "ctp6", 4.32);
    // ctp6.setError(0.03);


    // Input variables
    auto svmass= new RooRealVar("svmass", "M_{B^{0}} [GeV/c^{2}]",5.1,5.5);
    auto BsCt2DMC = new RooRealVar("BsCt2DMC","ct [cm]",0.007,0.5);
    auto BsCt2DMCErr = new  RooRealVar("BsCt2DMCErr", " #sigma_{ct}", 0.0002, 0.005,"cm");

    // Dataset
    auto data = new RooDataSet("data", "Bd data", RooArgSet(*svmass,*BsCt2DMC,*BsCt2DMCErr),Import(*tDATA));

    // Sidebands
    svmass->setRange("sbL",5.05,5.15);
    svmass->setRange("sbR",5.4,5.5);

    // --- MASS MODEL

    auto mass_mu     = new RooRealVar("mass_mu", "mass_mu", 5.27964, 5.27, 5.29);
    auto mass_lambda = new RooRealVar("mass_lambda", "mass_lambda",  0.02, 0., 0.1);
    auto mass_gamma  = new RooRealVar("mass_gamma", "mass_gamma",  0., -0.5, 0.5);
    auto mass_delta  = new RooRealVar("mass_delta", "mass_delta", 1., 0., 2.);
    auto massSigPdf  = new RooJohnsonLocal("massSigPdf","massSigPdf",*svmass, *mass_mu, *mass_lambda,*mass_gamma,*mass_delta);

    auto mass_bkgSlope = new RooRealVar("mass_bkgSlope", "bkg slope", -5, -10., 0.);
    auto massBkgPdf   =  new RooExponential("massBkgPdf", "mass bkg expo", *svmass, *mass_bkgSlope);
    auto mass_promtP1 = new RooRealVar("mass_promtP1","mass_promtP1", -0.1, -0.5, 0.5);
    auto massPromtPdf = new RooPolynomial("massPromtPdf", " prompt background mass",*svmass, RooArgList(*mass_promtP1));
    auto expo_frac = new RooRealVar("expo_frac", "Fraction of exponential background", 0.9, 0., 1.);

    // auto massBkgPdf   =  new RooAddPdf("massBkgPdf", "mass bkg", RooArgList(*massBkgExpo,*massPromtPdf), RooArgList(*expo_frac),true);

    auto nSig = new RooRealVar("nSig", "Number of Signal Events", 0.6*data->sumEntries(), 0., data->sumEntries());
    auto nBkg = new RooRealVar("nBkg", "Number of Backgound Events", 0.4*data->sumEntries(), 0., data->sumEntries());

    auto massFullPdf = new RooAddPdf("massFullPdf","Full mass Model",RooArgList(*massSigPdf,*massBkgPdf), RooArgList(*nSig,*nBkg));

    RooFitResult* massPreFit = massFullPdf->fitTo(*data,Save(),NumCPU(8));

    // --- CT ERROR MODEL
    // Sgn 
    // auto ctErrSig_gammashape1 = new RooRealVar("ctErrSig_gammashape1", "ct error pdf, first gamma shape", 8.0440, 0., 16);
    // auto ctErrSig_gammashape2 = new RooRealVar("ctErrSig_gammashape2", "ct error pdf, second gamma shape", 3.798, 0., 8);
    // auto ctErrSig_beta1  = new RooRealVar("ctErrSig_beta1", "ct error pdf, first gamma scale", 0.00014473, 0., 0.0003);
    // auto ctErrSig_beta2  = new RooRealVar("ctErrSig_beta2", "ct error pdf, second gamma scale", 0.000541, 0., 0.001);
    // auto ctErrSig_mu1 = new RooRealVar("ctErrSig_mu1", "ct error pdf, first gamma displacement", 0.000001, 0., 0.0002);
    // auto ctErrSig_mu2 = new RooRealVar("ctErrSig_mu2", "ct error pdf, second gamma displacement", 0.000001, 0., 0.0002);

    // auto ctErrSig_gamma1 = new RooGamma("ctErrSig_gamma1", "ct error pdf, first gamma", *BsCt2DMCErr, *ctErrSig_gammashape1, *ctErrSig_beta1, *ctErrSig_mu1);
    // auto ctErrSig_gamma2 = new RooGamma("ctErrSig_gamma2", "ct error pdf, second gamma", *BsCt2DMCErr, *ctErrSig_gammashape2, *ctErrSig_beta2, *ctErrSig_mu2);

    // auto ctErrSig_frac1 = new RooRealVar("ctErrSig_frac1", "frac", 0.5, 0., 1.);

    // auto ctErrSigPdf = new RooAddPdf("ctErrSigPdf", "ct error pdf", RooArgList(*ctErrSig_gamma1, *ctErrSig_gamma2), *ctErrSig_frac1);

    // auto ctErrBkg_gammashape1 = new RooRealVar("ctErrBkg_gammashape1", "ct error pdf, first gamma shape", 8.0440, 0., 16);
    // auto ctErrBkg_gammashape2 = new RooRealVar("ctErrBkg_gammashape2", "ct error pdf, second gamma shape", 3.798, 0., 8);
    // auto ctErrBkg_beta1  = new RooRealVar("ctErrBkg_beta1", "ct error pdf, first gamma scale", 0.00014473, 0., 0.0003);
    // auto ctErrBkg_beta2  = new RooRealVar("ctErrBkg_beta2", "ct error pdf, second gamma scale", 0.000541, 0., 0.001);
    // auto ctErrBkg_mu1 = new RooRealVar("ctErrBkg_mu1", "ct error pdf, first gamma displacement", 0.000001, 0., 0.0002);
    // auto ctErrBkg_mu2 = new RooRealVar("ctErrBkg_mu2", "ct error pdf, second gamma displacement", 0.000001, 0., 0.0002);

    // auto ctErrBkg_gamma1 = new RooGamma("ctErrBkg_gamma1", "ct error pdf, first gamma", *BsCt2DMCErr, *ctErrBkg_gammashape1, *ctErrBkg_beta1, *ctErrBkg_mu1);
    // auto ctErrBkg_gamma2 = new RooGamma("ctErrBkg_gamma2", "ct error pdf, second gamma", *BsCt2DMCErr, *ctErrBkg_gammashape2, *ctErrBkg_beta2, *ctErrBkg_mu2);

    // auto ctErrBkg_frac1 = new RooRealVar("ctErrBkg_frac1", "frac", 0.5, 0., 1.);

    // auto ctErrBkgPdf = new RooAddPdf("ctErrBkgPdf", "ct error pdf", RooArgList(*ctErrBkg_gamma1, *ctErrBkg_gamma2), *ctErrBkg_frac1);

    // auto MassCtErr_sigPdf = new RooProdPdf("MassCtErr_sigPdf", "",RooArgList(*massSigPdf,*ctErrSigPdf));
    // auto MassCtErr_bkgPdf  = new RooProdPdf("MassCtErr_bkgPdf", "",RooArgList(*massBkgPdf,*ctErrBkgPdf));  
    // auto MassCtErr_fullPdf = new RooAddPdf("MassCtErr_fullPdf","",RooArgList(*MassCtErr_sigPdf,*MassCtErr_bkgPdf), RooArgList(*nSig,*nBkg));


    // RooFitResult* masscterrBkgPreFit = ctErrBkgPdf->fitTo(*data,Save(),NumCPU(8),Range("sbL,sbR"));
    // RooFitResult* masscterrPreFit = MassCtErr_fullPdf->fitTo(*data,Save(),NumCPU(8));

    // ctErrSig_gammashape1->setConstant();
    // ctErrSig_gammashape2->setConstant();
    // ctErrSig_beta1->setConstant();
    // ctErrSig_beta2->setConstant();
    // ctErrSig_mu1->setConstant();
    // ctErrSig_mu2->setConstant();
    // ctErrSig_frac1->setConstant();

    // ctErrBkg_gammashape1->setConstant();
    // ctErrBkg_gammashape2->setConstant();
    // ctErrBkg_beta1->setConstant();
    // ctErrBkg_beta2->setConstant();
    // ctErrBkg_mu1->setConstant();
    // ctErrBkg_mu2->setConstant();
    // ctErrBkg_frac1->setConstant();

    nSig->setConstant();
    nBkg->setConstant();

    mass_mu->setConstant();
    mass_delta->setConstant();
    mass_gamma->setConstant();
    mass_lambda->setConstant();
    mass_bkgSlope->setConstant();

    // --- CT MODEL

    auto kappa = new RooRealVar("kappa", "kappa", kappa_val);
    auto resoModel = new RooGaussModel("resoModel","", *BsCt2DMC, RooConst(0), *BsCt2DMCErr, *kappa);

    auto ctau = new RooRealVar("ctau", "ctau",  ctau_PDG, 0.04, 0.05);

    auto ct_wReso = new RooDecay("ct_wReso","decay",*BsCt2DMC, *ctau, *resoModel, RooDecay::SingleSided);
    auto ctEffFunc = new RooFormulaVar("ctEffFunc","eff Model",
                                                "exp(@0+@1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,@2,@3,@4,@5,@6)",
                                                RooArgList(ctp0, ctp1, ctp2, ctp3, ctp4, ctp5, ctp6, *BsCt2DMC));

    auto ctSigPdf = new RooEffProd("ctSigPdf","model with efficiency", *ct_wReso, *ctEffFunc);

    //ct Bkg Model
    auto ctTruth = new RooTruthModel("ctTruth", "ctTruth", *BsCt2DMC);

    auto ct_bkgTau1 = new RooRealVar("ct_bkgTau1", "Bkg c#tau_{1}",0.004, 0, 0.5);
    auto ct_bkgTau2 = new RooRealVar("ct_bkgTau2", "Bkg c#tau_{2}",0.04, 0, 0.5);
    auto ct_bkgTau3 = new RooRealVar("ct_bkgTau3", "Bkg c#tau_{3}",0.05, 0, 0.5);
    auto ctBkg1 = new RooDecay("ctBkg1", "decay1", *BsCt2DMC, *ct_bkgTau1, *ctTruth, RooDecay::SingleSided);
    auto ctBkg2 = new RooDecay("ctBkg2", "decay2", *BsCt2DMC, *ct_bkgTau2, *ctTruth, RooDecay::SingleSided);
    auto ctBkg3 = new RooDecay("ctBkg3", "decay3", *BsCt2DMC, *ct_bkgTau3, *ctTruth, RooDecay::SingleSided);

    auto ct_bkgFrac1 = new RooRealVar("ct_bkgFrac1", "ct_bkgFrac1", 0.2, 0., 1.0);
    auto ct_bkgFrac2 = new RooRealVar("ct_bkgFrac2", "ct_bkgFrac2", 0.2, 0., 1.0);
    
    // auto ct_bkgMean = new RooRealVar("ct_bkgMean","ct_bkgMean",0.);
    // auto ct_bkgSigma = new RooRealVar("ct_bkgSigma","ct_bkgSigma",0.01, 0.,0.05);
    // auto ct_bkgPromt = new RooGaussian("ct_bkgPromt","ct_bkgPromt",*BsCt2DMC,*ct_bkgMean,*ct_bkgSigma);

    auto ctBkgPdf = new RooAddPdf("ctBkgPdf", "ctBg_pdf", RooArgList(*ctBkg1, *ctBkg2), RooArgList(*ct_bkgFrac1), true);

    RooFitResult* ctBkg_Prefit = ctBkgPdf->fitTo(*data,Save(),NumCPU(8),Range("sbL,sbR"));

    // --- FULL MODEL

    auto sigPdf = new RooProdPdf("sigPdf", "",RooArgList(*massSigPdf, *ctSigPdf));
    auto bkgPdf = new RooProdPdf("bkgPdf", "",RooArgList(*massBkgPdf, *ctBkgPdf));
  
    auto fullPdf = new RooAddPdf("fullPdf","Full 3D Model",RooArgList(*sigPdf,*bkgPdf), RooArgList(*nSig,*nBkg));

    // --- FIT
    
    RooFitResult* fitRes = fullPdf->fitTo(*data,Save(),NumCPU(8), PrintLevel(3), Verbose(true));

    fitRes->Print("v");

    // --- PLOTTING

    // Eff

    // Mass
    cout<<" START PLOTTING MASS "<<endl;
    auto plotMass = svmass->frame(Title("m [GeV/c^{2}]"),Bins(100));
    data->plotOn(plotMass,DataError(RooAbsData::SumW2));
    fullPdf->plotOn(plotMass);
    auto pullframe = svmass->frame(RooFit::Title("Mass pull"));
    auto hpull1 = plotMass->pullHist();
    pullframe->addPlotable(hpull1,"P0");
    pullframe->SetMinimum(-5);
    pullframe->SetMaximum(+5);
    pullframe->SetYTitle("pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);
    double chisquare_mass = plotMass->chiSquare();
    
    fullPdf->plotOn(plotMass, RooFit::LineColor(kOrange),RooFit::Components("massSigPdf"), RooFit::Name("signal1"), LineWidth(2), LineStyle(4));
    fullPdf->plotOn(plotMass,RooFit::LineColor(kRed),RooFit::Components("massBkgPdf"), RooFit::Name("background"), LineWidth(2), LineStyle(6));

    // ct
    cout<<" START PLOTTING CT "<<endl;
    RooPlot* plotCt = BsCt2DMC->frame(Title("ct [cm]"),Bins(100));
    data->plotOn(plotCt,DataError(RooAbsData::SumW2));
    fullPdf->plotOn(plotCt);
    auto pullframect = BsCt2DMC->frame(RooFit::Title("ct pull"));
    auto hpullct = plotCt->pullHist(0,0,true);
    pullframect->addPlotable(hpullct,"P0");
    pullframect->SetMinimum(-5);
    pullframect->SetMaximum(+5);
    pullframect->SetYTitle("pull");
    pullframect->SetMarkerStyle(20);
    pullframect->SetNdivisions(10);
    double chisquare_time = plotCt->chiSquare();

    fullPdf->plotOn(plotCt,RooFit::LineColor(kGreen),RooFit::Components("ctSigPdf"), RooFit::Name("signalct"), LineWidth(2), LineStyle(4));
    fullPdf->plotOn(plotCt,RooFit::LineColor(kRed), RooFit::Components("ctBkg1"), RooFit::Name("backgroundct"), LineWidth(2), LineStyle(2));
    fullPdf->plotOn(plotCt,RooFit::LineColor(kRed), RooFit::Components("ctBkg2"), RooFit::Name("backgroundct2"), LineWidth(2), LineStyle(2));
    fullPdf->plotOn(plotCt,RooFit::LineColor(kRed), RooFit::Components("ctBkg3"), RooFit::Name("backgroundct3"), LineWidth(2), LineStyle(2));
    // fullPdf->plotOn(plotCt,RooFit::LineColor(kOrange), RooFit::Components("ct_bkgPromt"), RooFit::Name("promptct"), LineWidth(2), LineStyle(2));

    auto c2 = new TCanvas("c2", "c2",0,0,600,600);
    auto pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    auto pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
    c2->SetFillColor(0);
    c2->SetBorderSize(2);
    c2->SetLeftMargin(0.1422222);
    c2->SetRightMargin(0.04444445);
    plotMass->SetStats(0);
    plotMass->Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c2->cd();
    c2->Print("mass_data.png");
  
    auto c3 = new TCanvas("c3", "c3",0,0,600,600);
    auto pad11 = new TPad("pad1","pad1",0,0.33,1,1);
    auto pad21 = new TPad("pad2","pad2",0,0,1,0.33);
    pad11->SetBottomMargin(0.00001);
    pad11->SetBorderMode(0);
    pad21->SetTopMargin(0.00001);
    pad21->SetBottomMargin(0.1);
    pad21->SetBorderMode(0);
    pad11->Draw();
    pad21->Draw();
    pad11->cd();
    pad11->SetLogy();
    gStyle->SetOptTitle(0);
    c3->SetFillColor(0);
    c3->SetBorderSize(2);
    c3->SetLeftMargin(0.1422222);
    c3->SetRightMargin(0.04444445);
    plotCt->SetStats(0);
    plotCt->Draw();
    pad21->cd();
    pullframect->SetStats(0);
    pullframect->Draw();
    c3->cd();
    c3->Print("ct_data.png");


    // CtErr
    // cout<<" START PLOTTING CTERR "<<endl;
    // auto plotCtErr = BsCt2DMCErr->frame(Title("#sigma_{ct} [cm]"),Bins(100));
    // data->plotOn(plotCtErr,DataError(RooAbsData::SumW2));
    // fullPdf->plotOn(plotCtErr);
    // auto pullframecterr = BsCt2DMCErr->frame(RooFit::Title("#sigma_{ct} pull"));
    // auto hpullcterr = plotCtErr->pullHist();
    // pullframecterr->addPlotable(hpullcterr,"P0");
    // pullframecterr->SetMinimum(-5);
    // pullframecterr->SetMaximum(+5);
    // pullframecterr->SetYTitle("pull");
    // pullframecterr->SetMarkerStyle(20);
    // pullframecterr->SetNdivisions(10);
    // double chisquare_cterr = plotCtErr->chiSquare();
    
    // fullPdf->plotOn(plotCtErr, RooFit::LineColor(kOrange),RooFit::Components("ctErrSigPdf"), RooFit::Name("signalcterr"), LineWidth(2), LineStyle(4));
    // fullPdf->plotOn(plotCtErr,RooFit::LineColor(kRed),RooFit::Components("ctErrBkgPdf"), RooFit::Name("backgroundcterr"), LineWidth(2), LineStyle(6));

    // auto c4 = new TCanvas("c4", "c4",0,0,600,600);
    // auto padcterr1 = new TPad("padcterr1","padcterr1",0,0.33,1,1);
    // auto padcterr2 = new TPad("padcterr2","padcterr2",0,0,1,0.33);
    // padcterr1->SetBottomMargin(0.00001);
    // padcterr1->SetBorderMode(0);
    // padcterr2->SetTopMargin(0.00001);
    // padcterr2->SetBottomMargin(0.1);
    // padcterr2->SetBorderMode(0);
    // padcterr1->Draw();
    // padcterr2->Draw();
    // padcterr1->cd();
    // gStyle->SetOptTitle(0);
    // c4->SetFillColor(0);
    // c4->SetBorderSize(2);
    // c4->SetLeftMargin(0.1422222);
    // c4->SetRightMargin(0.04444445);
    // plotCtErr->SetStats(0);
    // plotCtErr->Draw();
    // padcterr2->cd();
    // pullframecterr->SetStats(0);
    // pullframecterr->Draw();
    // c4->cd();
    // c4->Print("ctErr_data.png");

    
    cout << "final value of floating parameters" << endl;
    fitRes->floatParsFinal().Print("s");

    // cout<<"Chi square of cterr  fit is :"<< chisquare_cterr<< endl;
    cout<<"Chi square of lifetime  fit is :"<< chisquare_time<< endl;
    cout<<"Chi square of mass fit is :"<< chisquare_mass<< endl;
    cout<<endl;
    cout<<"Fitted ctau: "<<1e+04*ctau->getVal()<<" +- "<<1e+04*ctau->getError()<<endl;
    cout<<"PDG ctau: "<<1e+04*ctau_PDG<<" +- "<<1e+04*ctau_PDG_err<<endl;
    cout<<"Diff: "<<1e+04*(ctau->getVal() - ctau_PDG)<<endl;
    cout<<"Pull: "<<(ctau->getVal() - ctau_PDG)/sqrt(pow(ctau_PDG_err,2)+pow(ctau->getError(),2))<<endl;
    

}
