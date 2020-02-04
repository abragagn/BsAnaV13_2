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

void datadatacomparison(TString input1 = "../fittree_ntuBsData2017.root", TString input2 = "../fittree_ntuBsData2018.root", TCut selection = "1", TString outDir = "./outputData"){

  Int_t nbins = 50;

  cout << "nbins = " << nbins << endl;
  
  gStyle->SetOptTitle(0);

  // INPUT
  TFile *dataFile1 = new TFile(input1);
  TTree *dataTree1 = (TTree*)dataFile1->Get("treeFit");

  TFile *dataFile2 = new TFile(input2);
  TTree *dataTree2 = (TTree*)dataFile2->Get("treeFit");

  TFile *outFile = new TFile("comparisonData.root","RECREATE");

  TTree *selectedData1 = dataTree1->CopyTree(selection);
  TTree *selectedData2 = dataTree2->CopyTree(selection);

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

  RooRealVar *muonpeta      = new RooRealVar("muonpeta","#eta #mu_{1}",-2.4,2.4);
  RooRealVar *muonmeta      = new RooRealVar("muonmeta","#eta #mu_{2}",-2.4,2.4);
  RooRealVar *kaonpeta      = new RooRealVar("kaonpeta","#eta K_{1}",-2.5,2.5);
  RooRealVar *kaonmeta      = new RooRealVar("kaonmeta","#eta K_{2}",-2.5,2.5);

  RooRealVar *jpsimass      = new RooRealVar("jpsimass","J/#psi mass",2.95,3.25,"GeV");
  RooRealVar *phimass       = new RooRealVar("phimass","#phi(1020) mass",1.010,1.030,"GeV");

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

  RooRealVar mass_mu1("mass_mu1", "mass_mu1", 5.36679, 5.35, 5.37);
  RooRealVar mass_lambda1("mass_lambda1", "mass_lambda1", 0.5, 0, 1);
  RooRealVar mass_gamma1("mass_gamma1", "mass_gamma1", 0., -1, 1);
  RooRealVar mass_delta1("mass_delta1", "mass_delta1", 1., 0, 10);
  RooRealVar bkgSlope1("mass_bkgSlope1", "bkg_slope1", -100, 100);
  RooJohnsonLocal sgnPdf1("mass_sgn1", "mass_signal1", *svmass, mass_mu1, mass_lambda1, mass_gamma1, mass_delta1);
  RooExponential bkgPdf1("mass_bkg1", "mass_bkg1", *svmass, bkgSlope1);
  RooRealVar nSig1("nSig1", "Number of Signal Events 1", 1e+04, 0., 1e+07);
  RooRealVar nBkg1("nBkg1", "Number of BG events 1", 1e+04, 0., 1e+07);
  RooAddPdf massPdf1("mass_pdf1", "Total pdf 1", RooArgList(sgnPdf1, bkgPdf1), RooArgList(nSig1, nBkg1));

  RooRealVar mass_mu2("mass_mu2", "mass_mu2", 5.36679, 5.35, 5.37);
  RooRealVar mass_lambda2("mass_lambda2", "mass_lambda2", 0.5, 0, 1);
  RooRealVar mass_gamma2("mass_gamma2", "mass_gamma2", 0., -1, 1);
  RooRealVar mass_delta2("mass_delta2", "mass_delta2", 1., 0, 10);
  RooRealVar bkgSlope2("mass_bkgSlope2", "bkg_slope2", -100, 100);
  RooJohnsonLocal sgnPdf2("mass_sgn2", "mass_signal2", *svmass, mass_mu2, mass_lambda2, mass_gamma2, mass_delta2);
  RooExponential bkgPdf2("mass_bkg2", "mass_bkg2", *svmass, bkgSlope2);
  RooRealVar nSig2("nSig2", "Number of Signal Events 2", 1e+04, 0., 1e+07);
  RooRealVar nBkg2("nBkg2", "Number of BG events 2", 1e+04, 0., 1e+07);
  RooAddPdf massPdf2("mass_pdf2", "Total pdf 2", RooArgList(sgnPdf2, bkgPdf2), RooArgList(nSig2, nBkg2));

  // PLOT

  cout<<endl<<endl<<" ---- BEGIN LOOP ---- "<<endl<<endl;

  for(int t=0; t<varlist.GetSize(); ++t){

    cout<<endl<<"----- It's time for: "<<varlist.At(t)->GetName()<<"-----"<<endl<<endl;
    RooArgSet arg_set(*BsCt2DMC,*BsCt2DMCErr,*BscosthetaMC,*BscospsiMC,*BsphiMC,*svmass);

    TString lname = varlist.At(t)->GetName();

    if ( lname != "BsCt2DMC" 
      && lname != "BsCt2DMCErr" 
      && lname != "BscosthetaMC" 
      && lname != "BscospsiMC" 
      && lname != "BsphiMC" ) arg_set.add(*(RooRealVar*)varlist.At(t));

    RooDataSet *dataDataset1 = new RooDataSet("dataDataset1", "dataDataset1", arg_set, Import(*selectedData1));  
    RooDataSet *dataDataset2 = new RooDataSet("dataDataset2", "dataDataset2", arg_set, Import(*selectedData2));  

    massPdf1.fitTo(*dataDataset1,Extended());
    mass_mu1.setConstant();
    mass_lambda1.setConstant();
    mass_gamma1.setConstant();
    mass_delta1.setConstant();
    bkgSlope1.setConstant();

    massPdf2.fitTo(*dataDataset2,Extended());
    mass_mu2.setConstant();
    mass_lambda2.setConstant();
    mass_gamma2.setConstant();
    mass_delta2.setConstant();
    bkgSlope2.setConstant();

    RooMsgService::instance().setSilentMode(true);
    RooStats::SPlot *sData1 = new RooStats::SPlot("sData1","SPlot 1",*dataDataset1, &massPdf1, RooArgList(nSig1,nBkg1) );
    RooStats::SPlot *sData2 = new RooStats::SPlot("sData2","SPlot 2",*dataDataset2, &massPdf2, RooArgList(nSig2,nBkg2) );

    RooDataSet *dataw_z1 = new RooDataSet(dataDataset1->GetName(),dataDataset1->GetTitle(),dataDataset1,*dataDataset1->get(),0,"nSig1_sw");
    RooDataSet *dataw_z2 = new RooDataSet(dataDataset2->GetName(),dataDataset2->GetTitle(),dataDataset2,*dataDataset2->get(),0,"nSig2_sw");

    Float_t data1to2 = dataw_z1->sumEntries()/dataw_z2->sumEntries();

    dataw_z2->plotOn(
        (RooPlot*)plotlist.At(t),
        Rescale(data1to2),
        LineColor(4),
        MarkerColor(4),
        DataError(RooAbsData::SumW2),
        Binning(nbins)
      );

    dataw_z1->plotOn((RooPlot*)plotlist.At(t), DataError(RooAbsData::SumW2) , Binning(nbins), LineColor(2), MarkerColor(2));
    ((RooPlot*)plotlist.At(t))->SetMaximum(((RooPlot*)plotlist.At(t))->GetMaximum()*1.1*data1to2);
    ((RooPlot*)plotlist.At(t))->Write(varlist.At(t)->GetName());


  // Construct a histogram with the residuals of the data w.r.t. montecarlo, and a pull histo
    TH1 *hdata1 = dataw_z1->createHistogram(lname,nbins);
    TH1 *hdata2 = dataw_z2->createHistogram(lname,nbins);
    TH1 *hresid = (TH1*)hdata1->Clone();
    TH1 *hpull  = (TH1*)hdata1->Clone();
    TH1 *hdummy = (TH1*)hpull->Clone();

    hresid->Sumw2();
    hpull->Sumw2();
    // normalization factor
    Float_t pullscale; 

    hresid->Add(hdata2, -data1to2);

    for (Int_t i=1;i<=nbins;i++) {
      pullscale = 1.;
      if (hresid->GetBinError(i)!=0.) pullscale = 1./hresid->GetBinError(i);
      Double_t pull = hresid->GetBinContent(i)*pullscale;
      Double_t pullerror = hresid->GetBinError(i)*pullscale; 
      hpull->SetBinContent(i,pull);
      //      hpull->SetBinError(i,pullerror);
      hpull->SetBinError(i,0.); //  zero error
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
    Float_t pullmax = 3.;
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
