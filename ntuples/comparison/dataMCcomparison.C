#include <vector>
#include <cmath>
#include "Riostream.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TEventList.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TChain.h"
#include "TEventList.h"
#include "TTree.h"
#include "TVector.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooStats/SPlot.h"
#include "RooJohnsonLocal.cxx"

using namespace std;
using namespace RooFit;

TString datainput = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBsData2018.root";
TString mcinput = "/lustre/cmswork/abragagn/BPH/BsAnaV13_2/src/PDAnalysis/Ntu/bin/ntuples/fittree_ntuBsDG0MC2018.root";
TCut selection = "1";
TString outDir = "./outputDataMC/18";


void dataMCcomparison(){

  Int_t nbins = 50;

  cout << "nbins = " << nbins << endl;
  
  gStyle->SetOptTitle(0);

  // INPUT
  TFile *dataFile = new TFile(datainput);
  TTree *dataTree = (TTree*)dataFile->Get("treeFit");

  TFile *mcFile = new TFile(mcinput);
  TTree *mcTree = (TTree*)mcFile->Get("treeFit");

  TFile *outFile = new TFile("comparison.root","RECREATE");

  TTree *selectedData = dataTree->CopyTree(selection);
  TTree *selectedMC   = mcTree->CopyTree(selection);

  cout<<"Selected number of Data events: "<<selectedData->GetEntries()<<" / "<<dataTree->GetEntries()<<endl;
  cout<<"Selected number of MC events: "<<selectedMC->GetEntries()<<" / "<<mcTree->GetEntries()<<endl;

  // VARIABLES TO PLOT

  RooRealVar *svmass        = new RooRealVar("svmass","Bs mass", 5.24,5.49,"GeV");
  RooRealVar *BsCt2DMC      = new RooRealVar("BsCt2DMC","Bs ct", 0.007,0.5,"cm");
  RooRealVar *BsCt2DMCErr   = new RooRealVar("BsCt2DMCErr","Bs ct uncertainty", 0.00001,0.005,"cm");

  RooRealVar *BscosthetaMC  = new RooRealVar("BscosthetaMC","cos(#theta_{T})", -1,1);
  RooRealVar *BscospsiMC    = new RooRealVar("BscospsiMC","cos(#psi_{T})", -1,1);
  RooRealVar *BsphiMC       = new RooRealVar("BsphiMC","#phi_{T}", -TMath::Pi(),TMath::Pi(),"rad");

  RooRealVar *muonmpt       = new RooRealVar("muonmpt","p_{T} #mu_{1}",2,25,"GeV");
  RooRealVar *muonppt       = new RooRealVar("muonppt","p_{T} #mu_{2}",2,25,"GeV");
  RooRealVar *kaonppt       = new RooRealVar("kaonppt","p_{T} K_{1}",1,15,"GeV");
  RooRealVar *kaonmpt       = new RooRealVar("kaonmpt","p_{T} K_{2}",1,15,"GeV");

  RooRealVar *muonpeta       = new RooRealVar("muonpeta","#eta #mu_{1}",-2.4,2.4);
  RooRealVar *muonmeta       = new RooRealVar("muonmeta","#eta #mu_{2}",-2.4,2.4);
  RooRealVar *kaonpeta       = new RooRealVar("kaonpeta","#eta K_{1}",-2.5,2.5);
  RooRealVar *kaonmeta       = new RooRealVar("kaonmeta","#eta K_{2}",-2.5,2.5);

  RooRealVar *jpsimass       = new RooRealVar("jpsimass","J/#psi mass",2.95,3.25,"GeV");
  RooRealVar *phimass        = new RooRealVar("phimass","#phi(1020) mass",1.010,1.030,"GeV");

  RooPlot *frame_svmass       = svmass->frame();
  RooPlot *frame_BsCt2DMC     = BsCt2DMC->frame();
  RooPlot *frame_BsCt2DMCErr  = BsCt2DMCErr->frame();

  RooPlot *frame_BscosthetaMC = BscosthetaMC->frame();
  RooPlot *frame_BscospsiMC   = BscospsiMC->frame();
  RooPlot *frame_BsphiMC      = BsphiMC->frame();

  RooPlot *frame_muonmpt      = muonmpt->frame();
  RooPlot *frame_muonppt      = muonppt->frame();
  RooPlot *frame_kaonppt      = kaonppt->frame();
  RooPlot *frame_kaonmpt      = kaonmpt->frame();

  RooPlot *frame_muonpeta     = muonpeta->frame();
  RooPlot *frame_muonmeta     = muonmeta->frame();
  RooPlot *frame_kaonpeta     = kaonpeta->frame();
  RooPlot *frame_kaonmeta     = kaonmeta->frame();

  RooPlot *frame_jpsimass     = jpsimass->frame();
  RooPlot *frame_phimass      = phimass->frame();




  TList varlist;
  varlist.Add(svmass);
  varlist.Add(BsCt2DMC);
  varlist.Add(BsCt2DMCErr);

  varlist.Add(BscosthetaMC);
  varlist.Add(BscospsiMC);
  varlist.Add(BsphiMC);

  varlist.Add(muonmpt);
  varlist.Add(muonppt);
  varlist.Add(kaonppt);
  varlist.Add(kaonmpt);

  varlist.Add(muonpeta);
  varlist.Add(muonmeta);
  varlist.Add(kaonpeta);
  varlist.Add(kaonmeta);

  varlist.Add(jpsimass);
  varlist.Add(phimass);


  TList plotlist;
  plotlist.Add(frame_svmass);
  plotlist.Add(frame_BsCt2DMC);
  plotlist.Add(frame_BsCt2DMCErr);

  plotlist.Add(frame_BscosthetaMC);
  plotlist.Add(frame_BscospsiMC);
  plotlist.Add(frame_BsphiMC);

  plotlist.Add(frame_muonmpt);
  plotlist.Add(frame_muonppt);
  plotlist.Add(frame_kaonppt);
  plotlist.Add(frame_kaonmpt);

  plotlist.Add(frame_muonpeta);
  plotlist.Add(frame_muonmeta);
  plotlist.Add(frame_kaonpeta);
  plotlist.Add(frame_kaonmeta);

  plotlist.Add(frame_jpsimass);
  plotlist.Add(frame_phimass);


  // MASS FIT

  RooRealVar mass_mu("mass_mu", "mass_mu", 5.36679, 5.35, 5.37);
  RooRealVar mass_lambda("mass_lambda", "mass_lambda", 0.5, 0, 1);
  RooRealVar mass_gamma("mass_gamma", "mass_gamma", 0., -1, 1);
  RooRealVar mass_delta("mass_delta", "mass_delta", 1., 0, 10);
  RooRealVar bkgSlope("mass_bkgSlope", "bkg slope", -100, 100);

  RooJohnsonLocal sgnPdf("mass_sgn", "mass signal", *svmass, mass_mu, mass_lambda, mass_gamma, mass_delta);
  RooExponential bkgPdf("mass_bkg", "mass bkg", *svmass, bkgSlope);

  RooRealVar nSig("nSig", "Number of Signal Events", 1e+04, 0., 1e+07);
  RooRealVar nBkg("nBkg", "Number of BG events", 1e+04, 0., 1e+07);

  RooAddPdf massPdf("mass_pdf", "Total pdf", RooArgList(sgnPdf, bkgPdf), RooArgList(nSig, nBkg));

  // PLOT

  cout<<endl<<endl<<" ---- BEGIN LOOP ---- "<<endl<<endl;

  for(int t=0; t<varlist.GetSize(); ++t){
  //        for(int t=1; t<3; ++t){
    cout<<endl<<"----- It's time for: "<<varlist.At(t)->GetName()<<"-----"<<endl<<endl;
    RooArgSet arg_set(*BsCt2DMC,*BsCt2DMCErr,*BscosthetaMC,*BscospsiMC,*BsphiMC,*svmass);

    TString lname = varlist.At(t)->GetName();

    if ( lname != "BsCt2DMC" 
      && lname != "BsCt2DMCErr" 
      && lname != "BscosthetaMC" 
      && lname != "BscospsiMC" 
      && lname != "BsphiMC" ) arg_set.add(*(RooRealVar*)varlist.At(t));

    RooDataSet *dataDataset = new RooDataSet("dataDataset", "dataDataset", arg_set, Import(*selectedData));  
    cout<<"Used number of data events: "<<dataDataset->sumEntries()<<endl;

    massPdf.fitTo(*dataDataset,Extended());
    if (lname == "svmass") {
      dataDataset->plotOn(frame_svmass);
      massPdf.plotOn(frame_svmass);
      massPdf.plotOn(frame_svmass,Components("mass_sgn"),LineStyle(kDashed),LineColor(kGreen));
      massPdf.plotOn(frame_svmass,Components("mass_bkg"),LineStyle(kDashed),LineColor(kRed));
      frame_svmass->Write();
    }

    mass_mu.setConstant();
    mass_lambda.setConstant();
    mass_gamma.setConstant();
    mass_delta.setConstant();
    bkgSlope.setConstant();

    RooMsgService::instance().setSilentMode(true);
    RooStats::SPlot *sData = new RooStats::SPlot("sData","An SPlot",*dataDataset, &massPdf, RooArgList(nSig,nBkg) );

    dataDataset->Print("v");
    RooDataSet *dataw_z = new RooDataSet(dataDataset->GetName(),dataDataset->GetTitle(),dataDataset,*dataDataset->get(),0,"nSig_sw");
    dataw_z->Print("v");

    RooDataSet *mcDataset = new RooDataSet("mcDataset","mcDataset",arg_set, Import(*selectedMC));  
    cout<<"Used number of MC events: "<<mcDataset->sumEntries()<<endl;

    mcDataset->plotOn((RooPlot*)plotlist.At(t),
      Rescale(dataDataset->sumEntries()/mcDataset->sumEntries()*nSig.getVal()/(nSig.getVal()+nBkg.getVal())),
      LineColor(2),
      MarkerColor(2),
      DataError(RooAbsData::SumW2),
      Binning(nbins));

    dataw_z->plotOn((RooPlot*)plotlist.At(t), DataError(RooAbsData::SumW2) , Binning(nbins));
    ((RooPlot*)plotlist.At(t))->SetMaximum(((RooPlot*)plotlist.At(t))->GetMaximum()*1.1*dataDataset->sumEntries()/mcDataset->sumEntries()*nSig.getVal()/(nSig.getVal()+nBkg.getVal()));
    ((RooPlot*)plotlist.At(t))->Write(varlist.At(t)->GetName());


  // Construct a histogram with the residuals of the data w.r.t. montecarlo, and a pull histo
    TH1 *hmc   = mcDataset->createHistogram(lname,nbins);
    TH1 *hdata = dataw_z->createHistogram(lname,nbins);
    TH1 *hresid = (TH1*)hdata->Clone();
    TH1 *hpull  = (TH1*)hdata->Clone();
    TH1 *hdummy = (TH1*)hpull->Clone();

    hresid->Sumw2();
    hpull->Sumw2();
    // mc normalization factor
    Float_t datatomc = dataDataset->sumEntries()/mcDataset->sumEntries()*nSig.getVal()/(nSig.getVal()+nBkg.getVal());
    Float_t pullscale; 


    hresid->Add(hmc, -datatomc);

    for (Int_t i=1;i<=nbins;i++) {
      pullscale = 1.;
      if (hresid->GetBinError(i)!=0.)       pullscale = 1./hresid->GetBinError(i);
      //   if (hdata->GetBinError(i)!=0.)       pullscale = 1./hdata->GetBinError(i)/datatomc;
      Double_t pull = hresid->GetBinContent(i)*pullscale;
      Double_t pullerror = hresid->GetBinError(i)*pullscale; 
      // cout <<"i " << i <<  " error " << hresid->GetBinError(i) << " pullscale " << pullscale << " product " << hresid->GetBinError(i)*pullscale << endl; 
      hpull->SetBinContent(i,pull);
      //      hpull->SetBinError(i,pullerror);
      hpull->SetBinError(i,0.); //  zero error
      //      cout << "i " << " pull " << hpull-> GetBinContent(i) << " pull error " << hpull->GetBinError(i) << endl; 
      hdummy->SetBinContent(i,0.);
      hdummy->SetBinError(i,0.);
      hpull->GetYaxis()->SetLabelSize(12);
    }


    //    TCanvas c;
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "c",800,600);
    TPad *pad1 = new TPad("pad1", "pad1",0,0.18,1,1);
    TPad *pad2 = new TPad("pad2", "pad1",0,0,1,0.2);

    pad1->SetBottomMargin(0.1);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderSize(0);

    pad1->Draw();
    pad2->Draw();
   
    pad1->cd();

    ((RooPlot*)plotlist.At(t))->Draw();


    pad2->cd();

    //    Int_t ylabelsize = ((RooPlot*)plotlist.At(t))->GetYaxis()->GetLabelSize();
    hpull->GetYaxis()->SetLabelFont(63);
    hpull->GetYaxis()->SetLabelSize(12);
    hpull->GetXaxis()->SetLabelSize(0);

    // set pull range
    Float_t pullmax = 5.;
      if(hpull->GetMaximum()>pullmax) pullmax=hpull->GetMaximum()*1.1;
      if(hpull->GetMinimum()<-pullmax) pullmax=-hpull->GetMinimum()*1.1;
    hpull->GetYaxis()->SetRangeUser(-pullmax,pullmax);
    hpull->GetYaxis()->SetTitleSize(0);
    hpull->GetXaxis()->SetTitleSize(0);
    hpull->SetFillColor(38);
    hpull->Draw("B");
    hdummy->SetLineColor(1); 
    hdummy->Draw("same");

    //    Float_t ymax = hpull->GetMaximum();
    // TLine *line = new TLine(0,ymax,3,ymax);
    // line->SetLineColor(kRed);
    //line->Draw();

    c->Print(outDir + "/" + lname+".pdf");
    c->Print(outDir + "/" + lname+".png");

    delete hresid; 
    delete hpull; 
    delete hdummy; 
  }

  outFile->Write();
  outFile->Close();
}
