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

TString mcfile   = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBuMC2018.root";

TString mcfile_Bd   = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBdMC2018.root";

double ctau_mc      = 5.0095320e-02;

void ct2Dfit_BpMC()
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

    cout<<"neMC "<<neMC<<endl;

    auto fMC_Bd = new TFile(mcfile_Bd);
    auto tMC_Bd = (TTree*)fMC_Bd->Get("treeFit");
    int neMC_Bd = tMC_Bd->GetEntries();

    Double_t bins[] = {0.007,0.0073,0.0076,0.0079, 0.008, 0.009, 0.01,0.011,0.012, 0.013,0.014,0.015,0.016,0.017,0.018,0.019,  0.02165, 0.02865,0.03565,0.04265, 0.04965, 0.05665,  0.06365,0.07065,0.07765,0.08465,0.09165,0.09865,0.10565,0.11265,0.11965,0.12665,0.13365,0.14065,0.14765,0.15465,0.16165,0.16865,0.17565,0.18265,0.18965,0.19665,0.20365,0.21065,0.21765,0.22465,0.23165,0.23865,0.24565,0.25265,0.25965,0.26665,0.27365,0.28065,0.28765,0.29465,0.3, 0.32,0.34,0.36,0.38,0.4};
    Int_t  binnum = sizeof(bins)/sizeof(Double_t) - 1;

    auto hct_MC  = new TH1D("hct_MC","ct MC; ct [cm]", binnum, bins);
    auto hct_MC_even  = new TH1D("hct_MC_even","ct MC; ct [cm]", 100, 0.007,0.4);
    auto hct_GEN = new TH1D("hct_GEN","ct GEN from expo(ctau_PDG); ct [cm]", binnum, bins);
    auto resobulk = new TH1D("resobulk", "ct Resolution;ct_gen - ct_reco [cm]",  100, -0.004, 0.004);
    auto pullbulk = new TH1D("pullbulk", "ct Pulls; pull",  50, -5.0, 5.0);

    auto rnd = new TRandom3(458);
    
    float ctMC, ctMCerr;
    float ctMC_Bd, ctMCerr_Bd, ctMCgen_Bd;

    tMC->SetBranchAddress("BsCt2DMC",&ctMC);
    tMC->SetBranchAddress("BsCt2DMCErr",&ctMCerr);
    tMC_Bd->SetBranchAddress("BsCt2DMC",&ctMC_Bd);
    tMC_Bd->SetBranchAddress("BsCt2DMCErr",&ctMCerr_Bd);
    tMC_Bd->SetBranchAddress("BsCt2DMC_GEN",&ctMCgen_Bd);

    // Fill mc
    for (int i=0;i<neMC_Bd;i++){
        tMC_Bd->GetEntry(i);
        if(ctMCerr_Bd != ctMCerr_Bd) continue;
        if(ctMC_Bd != ctMC_Bd) continue;
        if(ctMCgen_Bd != ctMCgen_Bd) continue;
        resobulk->Fill(ctMCgen_Bd-ctMC_Bd);
        pullbulk->Fill((ctMCgen_Bd-ctMC_Bd)/ctMCerr_Bd);
    }

    float kappa_val = pullbulk->GetRMS();
    cout<<" --- kappa value : "<<kappa_val<<endl;

    for (int i=0;i<neMC;i++){
        tMC->GetEntry(i);
        if(ctMC != ctMC) continue;
        hct_MC->Fill(ctMC);
        hct_MC_even->Fill(ctMC);
    }

    // Fill gen
    for(int i=0;i<neMC;i++){
        double randct = rnd->Exp(ctau_mc);
        double randreso = resobulk->GetRandom();
        hct_GEN->Fill(randct + randreso);
    }

    // Make eff
    auto hEFF = (TH1D*)hct_MC->Clone("hEFF");
    hEFF->Sumw2();
    hEFF->Divide(hct_GEN);
    
    auto ctform = new TF1("ctEffFn","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.4);
    ctform->SetParameter(1,-1);
    ctform->FixParameter(2,1);      

    auto c1 = new TCanvas();

    hEFF->Fit(ctform);

    auto fittedEFF = hEFF->GetFunction("ctEffFn");
    hEFF->Draw();

//    TCanvas *cc = new TCanvas();
//    TF1 *ctModel = new TF1("ctModel"," [8]*TMath::Exp(-x/[7])*(expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6]))",0.007,0.4);

//    ctModel->FixParameter(0,fittedEFF->GetParameter(0));
//    ctModel->FixParameter(1,fittedEFF->GetParameter(1));
//    ctModel->FixParameter(2,fittedEFF->GetParameter(2));
//    ctModel->FixParameter(3,fittedEFF->GetParameter(3));
//    ctModel->FixParameter(4,fittedEFF->GetParameter(4));
//    ctModel->FixParameter(5,fittedEFF->GetParameter(5));
//    ctModel->FixParameter(6,fittedEFF->GetParameter(6));
//    ctModel->SetParameter(7,0.05);

//    hct_MC_even->Fit("ctModel","MIE");
//    hct_MC_even->Draw();

    // ROOFIT

    // Eff model parameters
    RooRealVar ctp0("ctp0", "ctp0", fittedEFF->GetParameter(0));
    ctp0.setError(fittedEFF->GetParError(0));
    RooRealVar ctp1("ctp1", "ctp1", fittedEFF->GetParameter(1), "cm^{-1}");
    ctp1.setError(fittedEFF->GetParError(1));
    RooRealVar ctp2("ctp2", "ctp2", fittedEFF->GetParameter(2));
    ctp2.setError(fittedEFF->GetParError(2));
    RooRealVar ctp3("ctp3", "ctp3", fittedEFF->GetParameter(3), "cm^{-1}");
    ctp3.setError(fittedEFF->GetParError(3));
    RooRealVar ctp4("ctp4", "ctp4", fittedEFF->GetParameter(4), "cm^{-2}");
    ctp4.setError(fittedEFF->GetParError(4));
    RooRealVar ctp5("ctp5", "ctp5", fittedEFF->GetParameter(5), "cm^{-3}");
    ctp5.setError(fittedEFF->GetParError(5));
    RooRealVar ctp6("ctp6", "ctp6", fittedEFF->GetParameter(6), "cm^{-4}");
    ctp6.setError(fittedEFF->GetParError(6));

    // Input variables
    auto svmass= new RooRealVar("svmass", "M_{B^{+}} (GeV/c^{2})",5.1,5.5);
    auto BsCt2DMC = new RooRealVar("BsCt2DMC","ct [cm]",0.007,0.4);
    auto BsCt2DMCErr = new  RooRealVar("BsCt2DMCErr", " #sigma_{ct}", 0.0007, 0.008,"cm");

    // Dataset
    auto data = new RooDataSet("data", "Bu mc", RooArgSet(*svmass, *BsCt2DMC,*BsCt2DMCErr),Import(*tMC));

    // --- CT MODEL

    auto kappa = new RooRealVar("kappa", "kappa", kappa_val);
    auto bias = new RooRealVar("bias","bias",0.);
    auto resoModel = new RooGaussModel("resoModel","", *BsCt2DMC, *bias, *kappa, *BsCt2DMCErr);

    auto ctau = new RooRealVar("ctau", "ctau",  ctau_mc, 0.04, 0.06);

    auto ctTruth = new RooTruthModel("ctTruth", "ctTruth", *BsCt2DMC);

    auto ct_wReso = new RooDecay("ct_wReso","decay",*BsCt2DMC, *ctau, *resoModel, RooDecay::SingleSided);
    auto ctEffFunc = new RooFormulaVar("ctEffFunc","eff Model", "exp(@0+@1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,@2,@3,@4,@5,@6)", RooArgList(ctp0, ctp1, ctp2, ctp3, ctp4, ctp5, ctp6, *BsCt2DMC));

    auto ctSigPdf = new RooEffProd("ctSigPdf","model with efficiency", *ct_wReso, *ctEffFunc);

    // --- FIT

    RooFitResult* fitRes = ctSigPdf->fitTo(*data,Save(),NumCPU(8));
// ConditionalObservables(RooArgSet(*BsCt2DMCErr))
    fitRes->Print("v");

    // --- PLOTTING

    // ct
    cout<<" START PLOTTING CT "<<endl;
    RooPlot* plotCt = BsCt2DMC->frame(Title("ct [cm]"),Bins(100));
    data->plotOn(plotCt,DataError(RooAbsData::SumW2));
    ctSigPdf->plotOn(plotCt);
    auto pullframect = BsCt2DMC->frame(RooFit::Title("ct pull"));
    auto hpullct = plotCt->pullHist(0,0,true);
    pullframect->addPlotable(hpullct,"P0");
    pullframect->SetMinimum(-5);
    pullframect->SetMaximum(+5);
    pullframect->SetYTitle("pull");
    pullframect->SetMarkerStyle(20);
    pullframect->SetNdivisions(10);
    double chisquare_time = plotCt->chiSquare();
  
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

    
    cout << "final value of floating parameters" << endl;
    fitRes->floatParsFinal().Print("s");

    cout<<"Chi square of lifetime  fit is :"<< chisquare_time<< endl;
    cout<<endl;
    cout<<"Fitted ctau: "<<1e+04*ctau->getVal()<<" +- "<<1e+04*ctau->getError()<<endl;
    cout<<"PDG ctau: "<<ctau_mc*1e+04<<endl;
    cout<<"Diff: "<<1e+04*(ctau->getVal() - ctau_mc)<<endl;
    cout<<"Pull: "<<(ctau->getVal() - ctau_mc)/ctau->getError()<<endl;
    

}
