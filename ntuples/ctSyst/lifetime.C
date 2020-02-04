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

using namespace RooFit ;
using namespace std;

TString Bfilename = "/afs/cern.ch/user/a/abragagn/public/BsJpsiPhiDatasets/fittree_ntuBdMC2018.root";
Double_t ctauBd = 0.045718350;
Double_t ctauBs = 0.04412945;
Double_t ctau   = ctauBd;

void lifetime()
{
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 100000 );

    TFile *fIn1 = new TFile(Bfilename);
    TTree *tIn1 = (TTree*)fIn1->Get("treeFit");
    Int_t n_entries1 = tIn1->GetEntries();

    TFile *fGEN = new TFile("/afs/cern.ch/user/a/abragagn/public/BsJpsiPhiDatasets/ntuBsGEN.root");
    TTree *tGEN = (TTree*)fGEN->Get("OutTree");
    Int_t n_GEN = tGEN->GetEntries();

    cout<<"nReco "<<n_entries1<<endl;
    cout<<"nGen "<<n_GEN<<endl;

    vector<double> bins = {0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025,
                           0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100,
                           0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.350, 0.400, 0.500};

    Int_t  binnum = bins.size() - 1;

    auto hctau_reco   = new TH1D("hctau_reco","ct(proper time) Efficiency Plot; ct (cm); #epsilon (ct) [a.u.] ", binnum, bins.data());
    auto hctau_reco_evenBinning       = new TH1D("hctau","ct(proper time); ct (cm);  ", 100, 0.007, 0.5);
    auto hctau_GEN_expo = new TH1D("hctau_GEN_expo","ct GEN from expo; ct[cm];", binnum, bins.data());

    auto resobulk       = new TH1D("resobulk", "Resobulk ; From B_{d}^{0}; Events ",  100, -0.01, 0.01);
    auto pullbulk       = new TH1D("pullbulk", "pullbulk ; From B_{d}^{0}; Events ",  50, -5.0, 5.0);

    TRandom3 *rnd = new TRandom3(36);

    // Fill reco
    Float_t ctreco, cterr, ctgen, ct_gensample;
    tIn1->SetBranchAddress("BsCt2DMC",&ctreco);
    tIn1->SetBranchAddress("BsCt2DMCErr",&cterr);
    tIn1->SetBranchAddress("BsCt2DMC_GEN",&ctgen);
    tGEN->SetBranchAddress("ctau_GEN",&ct_gensample);

    for (Int_t in=0; in<n_entries1; in++){
        tIn1->GetEntry(in);
        hctau_reco->Fill(ctreco);
        hctau_reco_evenBinning->Fill(ctreco);
        resobulk->Fill(ctgen-ctreco);
        pullbulk->Fill((ctgen-ctreco)/cterr);
    }

    // Fill gen
    for( int i=0; i<n_GEN; i++ )
    {
        hctau_GEN_expo->Fill(rnd->Exp(ctau) + resobulk->GetRandom());
    }


    // Make eff
    auto effgraph = (TH1D*)hctau_reco->Clone("effgraph");
    effgraph->Sumw2();
    // effgraph->Divide(hctau_GEN_sample);
    effgraph->Divide(hctau_GEN_expo);
    
    TF1 *effModel = new TF1("effModel","expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6])",0.007,0.5);
    effModel->SetParameter(1,-1);
    effModel->FixParameter(2,1);

    TCanvas *c1 = new TCanvas();
    effgraph->Fit("effModel","M");
    effgraph->Draw();

    TCanvas *c2 = new TCanvas();
    c2->Divide(2,2);
    c2->cd(1);
    resobulk->Draw();
    c2->cd(2);
    pullbulk->Draw();

    TCanvas *c3 = new TCanvas();
    TF1 *ctModel = new TF1("ctModel"," [8]*TMath::Exp(-x/[7])*(expo(0)*ROOT::Math::Chebyshev4(x,[2],[3],[4],[5],[6]))",0.007,0.5);
    TF1 *fit = effgraph->GetFunction("effModel");

    ctModel->FixParameter(0,fit->GetParameter(0));
    ctModel->FixParameter(1,fit->GetParameter(1));
    ctModel->FixParameter(2,fit->GetParameter(2));
    ctModel->FixParameter(3,fit->GetParameter(3));
    ctModel->FixParameter(4,fit->GetParameter(4));
    ctModel->FixParameter(5,fit->GetParameter(5));
    ctModel->FixParameter(6,fit->GetParameter(6));
    ctModel->SetParameter(7,0.045);

    hctau_reco_evenBinning->Fit("ctModel","MIE");
    hctau_reco_evenBinning->Draw();

    TF1 *fitCt = hctau_reco_evenBinning->GetFunction("ctModel");

    cout<<endl<<"ct pull = "<<(fitCt->GetParameter(7)-ctau)/fitCt->GetParError(7)<<endl;
}


