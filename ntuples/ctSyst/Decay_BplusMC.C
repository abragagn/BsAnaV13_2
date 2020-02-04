/// \author  - Muhammad Alibordi
// Test of RooKeyPDF has ability to discriminate the background and

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooJohnsonLocal.cxx"

using namespace RooFit ;
using namespace std;

string filename = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBuMC2018.root";
string plotString = "";
string myInputFile = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBuData2018A.root";

TString bdfile = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBdMC2018.root";

float evtfact(int n){
    int prod =1;
    for(int k=1;k<=n;k++){prod=prod*k;}
    return prod;
    
}


void Decay_Bplus()
{
   
    Int_t nbins = 100;
    //std::cout<<" Give the value of bin number "<<"\n";
    //cin>>nbins;
    
    
    auto ctform = new TF1("ctEffFn","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.4);
    ctform->SetParameter(1,-1);
    ctform->FixParameter(2,1);

    
    
    Double_t bins[] = {0.007,0.0073,0.0076,0.0079, 0.008, 0.009, 0.01,0.011,0.012, 0.013,0.014,0.015,0.016,0.017,0.018,0.019,  0.02165, 0.02865,0.03565,0.04265, 0.04965, 0.05665,  0.06365,0.07065,0.07765,0.08465,0.09165,0.09865,0.10565,0.11265,0.11965,0.12665,0.13365,0.14065,0.14765,0.15465,0.16165,0.16865,0.17565,0.18265,0.18965,0.19665,0.20365,0.21065,0.21765,0.22465,0.23165,0.23865,0.24565,0.25265,0.25965,0.26665,0.27365,0.28065,0.28765,0.29465,0.3, 0.32,0.34,0.36,0.38,0.4};
    Int_t  binnum = sizeof(bins)/sizeof(Double_t) - 1;
    
    auto resobulk = new TH1D("resobulk", "Resobulk ; From B_{d}^{0}; Events ",  nbins, -0.004, 0.004);
    auto pullbulk = new TH1D("pullbulk", "pullbulk ; From B_{d}^{0}; Events ",  nbins, -5.0, 5.0);
    auto cterror = new TH1D("cterror", "CtError; CtErrorEvents ",  nbins, 0.0005, 0.005);
    auto hctaucut = new TH1D("hctaucut","Recolifetime  ;ct[cm];#epsilon [a.u.]",  binnum, bins);
    auto hctaunocut = new TH1D("hctaunocut","Gen lifetime ;ct[cm];#epsilon [a.u.]",  binnum, bins);
    
    Float_t ctreco , npv, cterr;
    Double_t kappa_val;
    Double_t par1, par2, par3, par4, par5, par6, par7;
    Double_t time = 5.0095320e-02;//2018 pdg Bplus 491.1micron, B0 = 455.7 micron, 500.95320
    
    TTree          *fChain, *copy;
    TFile* f = new TFile(filename.c_str());
    fChain = new TTree;
    fChain=(TTree*)f->Get("treeFit");
    if (!fChain){
        cout << "No TTree found in input file, returning" << endl;
        return;
    }
    TRandom3 *decay = new TRandom3();
    fChain->SetBranchAddress("BsCt2DMC",&ctreco);
    fChain->SetBranchAddress("BsCt2DMCErr",&cterr);
    Long64_t nentries = fChain->GetEntries();
    cout << "Start Processing " << nentries << " events" <<"\n";
    for (int jentry=0; jentry<nentries;jentry++) {

        fChain->GetEntry(jentry);
        if (cterr !=cterr)continue;
        if (ctreco !=ctreco)continue;
        hctaucut->Fill(ctreco);
        cterror->Fill(cterr);
    }
    
    //==========================The resolution and the pull is from B0 17MC/18MC sample
    
    
    TFile *fb0 = new TFile(bdfile);
    
    TTree* b0file = (TTree*)fb0->Get("treeFit");
    Int_t n_entriesb0 = b0file->GetEntries();
    Float_t ctrecob0, ctgenb0, cterrb0;
    b0file->SetBranchAddress("BsCt2DMC",&ctrecob0);
    b0file->SetBranchAddress("BsCt2DMC_GEN",&ctgenb0);
    b0file->SetBranchAddress("BsCt2DMCErr",&cterrb0);
    for( int k = 0 ; k < nentries ; k++ )
    {
    b0file->GetEntry(k);
        resobulk->Fill(ctgenb0-ctrecob0);
        pullbulk->Fill((ctgenb0-ctrecob0)/cterrb0);
    }
    
    TCanvas *c4= new TCanvas("c4", "c4",0,0,800,600);
    resobulk->Draw();//hctaucut->Draw();//resobulk->Draw();//
    TCanvas *c3= new TCanvas("c3", "c3",0,0,800,600);
    pullbulk->Draw();//hctaunocut->Draw();//pullbulk->Draw();//
    
    
    for( int i = 0 ; i < nentries ; i++ )
    {
        Double_t ctaugen1 =  decay->Exp(time);
        Double_t actualctgen1 =  ctaugen1;// + cterror->GetRandom();//resobulk->GetRandom();
        hctaunocut->Fill(actualctgen1);
    }
    
    kappa_val = pullbulk->GetRMS();
    std::cout<<" kappa value for main pull:  "<<"The Value of kappa: "<<pullbulk->GetRMS()<<"\n";
    TCanvas *ceff= new TCanvas("ceff", "ceff",0,0,800,600);
    auto hdivideI = (TH1D*)hctaucut->Clone("hdivideI");
    hdivideI->Sumw2();
    hdivideI->Divide(hctaunocut);
    hdivideI->Draw("colz");    //colz ep
    hdivideI->Fit(ctform);
    TFitResultPtr r =    hdivideI->Fit(ctform, "S");
    TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
    Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
    Double_t par0   = r->Parameter(0);            // retrieve the value for the parameter 0
    Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
    r->Print("V");                                // print full information of fit including covariance matrix
    
    
    
    RooRealVar ctp0("ctp0", "ctp0", ctform->GetParameter(0));
    ctp0.setError(ctform->GetParError(0));
    RooRealVar ctp1("ctp1", "ctp1", ctform->GetParameter(1), "cm^{-1}");
    ctp1.setError(ctform->GetParError(1));
    RooRealVar ctp2("ctp2", "ctp2", ctform->GetParameter(2));
    ctp2.setError(ctform->GetParError(2));
    RooRealVar ctp3("ctp3", "ctp3", ctform->GetParameter(3), "cm^{-1}");
    ctp3.setError(ctform->GetParError(3));
    RooRealVar ctp4("ctp4", "ctp4", ctform->GetParameter(4), "cm^{-2}");
    ctp4.setError(ctform->GetParError(4));
    RooRealVar ctp5("ctp5", "ctp5", ctform->GetParameter(5), "cm^{-3}");
    ctp5.setError(ctform->GetParError(5));
    RooRealVar ctp6("ctp6", "ctp6", ctform->GetParameter(6), "cm^{-4}");
    ctp6.setError(ctform->GetParError(6));
    
    //==================
    //==Model
    //==================
    
    TChain* chain_data = new TChain("treeFit");
    //chain_data->Add("/afs/cern.ch/work/m/mumuhamm/private/CMSSW_10_2_7/src/Bplus/fittree_ntuBpdata2018A.root");
    chain_data->Add(myInputFile.c_str());
    Int_t nevt = (int)chain_data->GetEntries();
    std::cout<<"Number of total events"<<nevt<<"\n";
    
    
    RooRealVar *svmass= new RooRealVar("svmass", "M_{B^{+}} (GeV/c^{2})",5.18,5.4);
    RooRealVar *BsCt2DMC = new RooRealVar("BsCt2DMC","proper time (ct) (cm)",0.007,0.4);
    RooRealVar *BsCt2DMCErr = new  RooRealVar("BsCt2DMCErr", " #Deltact_{B0}",0.0007, 0.008,"cm");
    
    RooDataSet* data = new RooDataSet("data", "raw data1", RooArgSet(*svmass,*BsCt2DMC,*BsCt2DMCErr),Import(*chain_data));
    
    RooRealVar *conts = new RooRealVar("conts","conts", -3.0, 0.0);
    RooExponential *expoBg = new RooExponential("expoBg","expoBG",*svmass,*conts);
    
    RooRealVar *prompt_p1 =new RooRealVar("prompt_p1","prompt background mass polynomial", 0.1, -0.2, 0.3);
    RooPolynomial *prompt_mass =new RooPolynomial("prompt_mass", " prompt background mass",*svmass, RooArgList(*prompt_p1));
    
    
    RooRealVar *mu= new RooRealVar("mu", "mu", 5.27929, 5.2, 5.35);
    RooRealVar *lambda = new RooRealVar("lambda", "lambda",  0.5, 0, 1);
    RooRealVar *gamma = new RooRealVar("gamma", "gamma",  0., -1, 1);
    RooRealVar *delta = new RooRealVar("delta", "delta", 1., 0, 10);
    RooJohnsonLocal *john = new RooJohnsonLocal("john", "john", *svmass, *mu, *lambda, *gamma, *delta);
    
    
    
    
    RooRealVar *tau = new RooRealVar("tau", "tau",  5.0095320e-02, 0.04, 0.06);
    RooRealVar *kappa = new RooRealVar("kappa", "kappa", kappa_val);
    RooFormulaVar *kapreso = new RooFormulaVar("kapreso", "kappa*BsCt2DMCErr", RooArgList(*kappa, *BsCt2DMCErr));
    RooRealVar *bias = new RooRealVar("bias","bias",0);
    RooRealVar *sdscale = new RooRealVar("sdscale","per-event error scale factor",1);
    RooGaussModel *gm_sig = new RooGaussModel("gm_sig","gauss model scaled bt per-event error", *BsCt2DMC, *bias, *sdscale, *kapreso);
    RooDecay *decay_gm = new RooDecay("decay_gm","decay",*BsCt2DMC, *tau, *gm_sig, RooDecay::SingleSided);
    RooFormulaVar *ctEffFunc = new RooFormulaVar("ctEffFunc","eff", "exp(@0+@1*BsCt2DMC)*ROOT::Math::Chebyshev4(BsCt2DMC,@2,@3,@4,@5,@6)",RooArgList(ctp0, ctp1, ctp2, ctp3, ctp4, ctp5, ctp6, *BsCt2DMC));
    
    RooEffProd *lifeEff = new RooEffProd("lifeEff","model with efficiency-Full 2D ", *decay_gm, *ctEffFunc);
    
    RooRealVar *bgtau1= new RooRealVar("bgtau1", "Bg ct_{1}", 0, 0.5);
    RooRealVar *bgtau2= new RooRealVar("bgtau2", "Bg ct_{2}", 0, 0.5);
    RooRealVar *bgtau3= new RooRealVar("bgtau3", "Bg ct_{3}",3.74228e-03, 0, 0.5);
    RooDecay *bgdecay1= new RooDecay("bgdecay1", "decay1",*BsCt2DMC, *bgtau1, *gm_sig, RooDecay::SingleSided);
    RooDecay *bgdecay2= new RooDecay("bgdecay2", "decay2",*BsCt2DMC, *bgtau2, *gm_sig, RooDecay::SingleSided);
    RooDecay *bgdecay3= new RooDecay("bgdecay3", "decay3",*BsCt2DMC, *bgtau3, *gm_sig, RooDecay::SingleSided);
    
    RooRealVar *g_mean = new RooRealVar("g_mean","g_mean",0.0,0.06) ;
    RooRealVar *g_sigma =new RooRealVar("g_sigma","g_sigma",0.0,1.0) ;
    RooGaussian *ct_gauss =new RooGaussian("ct_gauss","ct_gauss",*BsCt2DMC,*g_mean,*g_sigma) ;
    RooRealVar *frac1 = new RooRealVar("frac1", "frac1", 0.1, 0.9);
    RooAddPdf *proptm_bg = new RooAddPdf("proptm_bg", "ctBg_pdf", RooArgList(*bgdecay1, *bgdecay2), RooArgList(*frac1), true);
    
   
     RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*ct*cterr",RooArgList(*john,*decay_gm));
    //RooProdPdf *sigpdf = new RooProdPdf("sigpdf", "mass*ct*cterr",RooArgList(*john,*lifeEff));
    RooProdPdf *Combi_bg = new RooProdPdf("Combi_bg", "massbg*ctbg*cterrbg",RooArgList(*expoBg,*proptm_bg));
    RooProdPdf *Prompt_bg = new RooProdPdf("Prompt_bg", "massbg*ctbg*cterrbg",RooArgList(*prompt_mass,*ct_gauss));
    
    
    
    
    RooRealVar *nSig = new RooRealVar("nSig", "Number of Signal Events in SIGNAL MC",0.7e+05,(int)chain_data->GetEntries());//e05 for all era 3 for 17D
    RooRealVar *nBkg_Combi = new RooRealVar("nBkg_Combi", "Number of Backgound Events in produced  MC", 1000 ,(int)chain_data->GetEntries());//1000 for all other era
    RooRealVar *nBkg_Prompt = new RooRealVar("nBkg_Prompt", "Number of Backgound Events in produced  MC", 500 ,(int)chain_data->GetEntries());//500 for all other era
    
    RooAddPdf *MTpdf = new RooAddPdf("MTpdf","MTpdf_inter",RooArgList(*sigpdf,*Combi_bg,*Prompt_bg), RooArgList(*nSig,*nBkg_Combi,*nBkg_Prompt));
     //RooAddPdf *MTpdf = new RooAddPdf("MTpdf","MTpdf_inter",RooArgList(*sigpdf,*Combi_bg), RooArgList(*nSig,*nBkg_Combi));//, Conditional(*BsCt2DMC,*BsCt2DMCErr));
    //RooAddPdf *MTpdf_inter = new RooAddPdf("MTpdf_inter","MTpdf_inter",RooArgList(*sigpdf,*Combi_bg,*Prompt_bg), RooArgList(*nSig,*nBkg_Combi,*nBkg_Prompt));    
    //RooEffProd *MTpdf = new RooEffProd("MTpdf","model with efficiency-Full 2D ", *MTpdf_inter, *ctEffFunc);//
    
    
    RooAbsReal* nll = MTpdf->createNLL(*data,NumCPU(8)) ;
    RooMinuit(*nll).migrad() ;
    RooMinuit(*nll).minos();
    RooMinuit(*nll).hesse();
    
    
    RooFitResult* fitRes = MTpdf->fitTo(*data,Save(),NumCPU(8));//,ConditionalObservables(*BsCt2DMCErr));//data_SigReg
    fitRes->Print("v");
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    TH2* hcorr = fitRes->correlationHist() ;
    TCanvas* ccor = new TCanvas("Corr Matrix","M/ct/cterr correaltion matrix",800,400) ;
    gPad->SetLeftMargin(0.15) ; hcorr->GetYaxis()->SetTitleOffset(1.4) ; hcorr->Draw("colz") ;ccor->Print(("correlation_Matrix_"+plotString+"_data.png").c_str(),"png");
    
    
    
    
    RooPlot* bsmass = svmass->frame(Title("M_{B^{+}} (GeV/c^{2})"),Bins(nbins));
    data->plotOn(bsmass,DataError(RooAbsData::SumW2));
    MTpdf->plotOn(bsmass) ;
    MTpdf->paramOn(bsmass);
    RooPlot* pullframe = svmass->frame(RooFit::Title("Mass pull"));
    RooHist* hpull1 = bsmass->pullHist();
    pullframe->addPlotable(hpull1,"P0") ;
    pullframe->SetMinimum(-3) ;
    pullframe->SetMaximum(+3) ;
    pullframe->SetYTitle("pull");
    pullframe->SetMarkerStyle(20);
    pullframe->SetNdivisions(10);
    Double_t chisquare_mass = bsmass->chiSquare();
    
    MTpdf->plotOn(bsmass, RooFit::LineColor(kGreen),RooFit::Components("john"), RooFit::Name("signal"), LineWidth(2), LineStyle(4));
    MTpdf->plotOn(bsmass,RooFit::LineColor(kRed),RooFit::Components("expoBg"), RooFit::Name("combinatorial"), LineWidth(2), LineStyle(6));
    MTpdf->plotOn(bsmass,RooFit::LineColor(kMagenta),RooFit::Components("prompt_mass"), RooFit::Name("prompt"), LineWidth(2), LineStyle(6));
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(bsmass->findObject("signal"),"B^{+}#rightarrow J/#psi K^{+}","l");
    leg->AddEntry(bsmass->findObject("combinatorial"),"Combinatorial","l");
    leg->AddEntry(bsmass->findObject("prompt"),"B^{+}#rightarrow J/#psi X","l");
  
    RooPlot* decaytau = BsCt2DMC->frame(Title("Proper decay "),Bins(100));
    data->plotOn(decaytau,DataError(RooAbsData::SumW2));
    MTpdf->plotOn(decaytau);
    MTpdf->paramOn(decaytau);
    RooPlot* pullframect = BsCt2DMC->frame(RooFit::Title("ct pull"));
    RooHist* hpullct = decaytau->pullHist();
    pullframect->addPlotable(hpullct,"P0") ;
    pullframect->SetMinimum(-3) ;
    pullframect->SetMaximum(+3) ;
    pullframect->SetYTitle("pull");
    pullframect->SetMarkerStyle(20);
    pullframect->SetNdivisions(10);
    Double_t chisquare_time = decaytau->chiSquare();
    
    MTpdf->plotOn(decaytau,RooFit::LineColor(kGreen),RooFit::Components("decay_gm"), RooFit::Name("signalct"), LineWidth(2), LineStyle(4));
    MTpdf->plotOn(decaytau,RooFit::LineColor(kRed), RooFit::Components("proptm_bg"), RooFit::Name("combinatorialct"), LineWidth(2), LineStyle(2));
    MTpdf->plotOn(decaytau,RooFit::LineColor(kMagenta),RooFit::Components("ct_gauss"), RooFit::Name("promptct"), LineWidth(2), LineStyle(6));

    TLegend *legct = new TLegend(0.7,0.7,0.9,0.9);
    legct->AddEntry(decaytau->findObject("signalct"),"B^{+}#rightarrow J/#psi K^{+}","l");
    legct->AddEntry(decaytau->findObject("combinatorialct"),"Combinatorial","l");
    legct->AddEntry(decaytau->findObject("promptct"),"B^{+}#rightarrow J/#psi X","l");   

 
    TCanvas *c = new TCanvas("c", "c",0,0,600,600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gStyle->SetOptTitle(0);
    c->SetFillColor(0);
    c->SetBorderSize(2);
    c->SetLeftMargin(0.1422222);
    c->SetRightMargin(0.04444445);
    bsmass->SetStats(0);
    bsmass->Draw();
    leg->Draw("same");
    auto cms1 = new TLatex(5.15, 5300, "#bf{CMS} #it{Preliminary} 2018A, #sqrt{s} = 13 TeV");
    cms1->SetNDC(false);
    cms1->SetTextColor(12);
    cms1->SetTextFont(42);
    cms1->SetTextSize(0.055);
    cms1-> Draw();
    pad2->cd();
    pullframe->SetStats(0);
    pullframe->Draw();
    c->cd();
    c->Print(("Bplus_Mass_"+plotString+"_data.root").c_str(), "root");
  
    
    TCanvas *cc = new TCanvas("cc", "cc",0,0,600,600);
    TPad *pad11 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad21 = new TPad("pad2","pad2",0,0,1,0.33);
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
    // cc->Range(-0.3786885,-19.31166,0.1745902,168.5381);
    cc->SetFillColor(0);
    cc->SetBorderSize(2);
    cc->SetLeftMargin(0.1422222);
    cc->SetRightMargin(0.04444445);
    decaytau->SetStats(0);
    decaytau->Draw();
    legct->Draw("same");
    auto cms11 = new TLatex(0.007, 26200, "#bf{CMS} #it{Preliminary} 2018A, #sqrt{s} = 13 TeV");
    cms11->SetNDC(false);
    cms11->SetTextColor(12);
    cms11->SetTextFont(42);
    cms11->SetTextSize(0.055);
    cms11-> Draw();
    pad21->cd();
    pullframect->SetStats(0);
    pullframect->Draw();
    cc->cd();
    cc->Print(("Bplus_lifetiem_"+plotString+"_data.root").c_str(), "root");
    

    Double_t taulow, tauhigh, gsigmalow, gsigmahigh;
    taulow = tau->getVal() - tau->getAsymErrorHi();
    tauhigh =tau->getVal() + tau->getAsymErrorHi();
    gsigmalow = g_sigma->getVal()-g_sigma->getAsymErrorHi();
    gsigmahigh = g_sigma->getVal()+g_sigma->getAsymErrorHi();
   
    RooPlot* frame = new RooPlot(*tau, *g_sigma, taulow, tauhigh, gsigmalow, gsigmahigh) ;
    frame->SetTitle("Covariance between #tau and #sigma_{G}") ;
    fitRes->plotOn(frame,*tau,*g_sigma,"ME12ABHV") ;
    
    cout << "final value of floating parameters" << endl ;
    fitRes->floatParsFinal().Print("s") ;

    cout<<"Chi square of lifetime  fit is :"<< chisquare_time<< endl;
    cout<<"Chi square of mass fit is :"<< chisquare_mass<< endl; 
    
}
