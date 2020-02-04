#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"

using namespace std;

int nbinsx = 70;
int nbinsy = 70;
int nbinsz = 30;

int nhist = 100;

void testEff2()
{
    TFile *f = new TFile("/afs/cern.ch/user/a/abragagn/public/BsJpsiPhiDatasets/fittree_ntuBsDG0MC2018.root");
    TFile *fGen = new TFile("/afs/cern.ch/user/a/abragagn/public/BsJpsiPhiDatasets/ntuBsGEN.root");

    TTree *t = (TTree*)f->Get("treeFit");
    TTree *tGEN = (TTree*)fGen->Get("OutTree");

    TH3D *hRECO = new TH3D("hRECO","",nbinsx,-1.,1.,nbinsy,-1.,1., nbinsz, -TMath::Pi(), TMath::Pi());
    TH3D *hGEN = new TH3D("hGEN","",nbinsx,-1.,1.,nbinsy,-1.,1., nbinsz, -TMath::Pi(), TMath::Pi());

    TH3D *aRECO[nhist];
    TH3D *aGEN[nhist];
    for(int i=0; i<nhist; ++i){
        aRECO[i] = (TH3D*)hRECO->Clone(Form("aRECO_%i", i));
        aGEN[i] = (TH3D*)hGEN->Clone(Form("aGEN_%i", i));
    }

    float costheta, cospsi, phi;
    float costheta_gen, cospsi_gen, phi_gen;
    t->SetBranchAddress("BscosthetaMC",&costheta);
    t->SetBranchAddress("BscospsiMC",&cospsi);
    t->SetBranchAddress("BsphiMC",&phi);
    tGEN->SetBranchAddress("angle_costheta_GEN",&costheta_gen);
    tGEN->SetBranchAddress("angle_cospsi_GEN",&cospsi_gen);
    tGEN->SetBranchAddress("angle_phi_GEN",&phi_gen);

    for(int i=0; i<t->GetEntries(); i++){
        t->GetEntry(i);
        hRECO->Fill(costheta,cospsi,phi);
    }

    for(int i=0; i<tGEN->GetEntries(); i++){
        tGEN->GetEntry(i);
        hGEN->Fill(costheta_gen,cospsi_gen,phi_gen);
    }

    int nzero=0;
    int nzero_gen=0;
    for(int j=1; j<=nbinsx; ++j)
        for(int k=1; k<=nbinsy; ++k)
            for(int m=1; m<=nbinsz; ++m){
                if(hRECO->GetBinContent(hRECO->GetBin(j,k,m)) == 0) nzero++;
                if(hGEN->GetBinContent(hGEN->GetBin(j,k,m)) == 0) nzero_gen++;
            }

    cout<<" --- reco bin with zero content = "<<nzero<<" out of "<<nbinsx*nbinsy*nbinsz<<endl;
    cout<<" --- gen bin with zero content = "<<nzero_gen<<" out of "<<nbinsx*nbinsy*nbinsz<<endl;

    for(int i=0; i<nhist; ++i){
        cout<<" --- generating histograms #"<<i<<endl;
        aRECO[i]->FillRandom(hRECO, gRandom->Poisson(hRECO->Integral()));
        aGEN[i]->FillRandom(hGEN, gRandom->Poisson(hGEN->Integral()));
    }

    cout<<"--- histograms generated ---"<<endl;

    auto hEFF = (TH3D*)hRECO->Clone("hEFF");
    hEFF->Sumw2();
    hEFF->Divide(hGEN);

    TH3D *aEFF[nhist];
    TH3D *aDIFF_reco[nhist];
    TH3D *aDIFF_gen[nhist];
    TH3D *aDIFF_eff[nhist];

    for(int i=0; i<nhist; ++i){
        aEFF[i] = (TH3D*)hRECO->Clone(Form("aEFF_%i",i));
        aEFF[i]->Sumw2();
        aEFF[i]->Divide(aGEN[i]);

        aDIFF_reco[i] = (TH3D*)aRECO[i]->Clone(Form("aDIFF_reco_%i",i));
        aDIFF_gen[i] = (TH3D*)aGEN[i]->Clone(Form("aDIFF_gen_%i",i));
        aDIFF_eff[i] = (TH3D*)aEFF[i]->Clone(Form("aDIFF_eff_%i",i));

        for(int j=1; j<=nbinsx; ++j){
            for(int k=1; k<=nbinsy; ++k){
                for(int m=1; m<=nbinsz; ++m){
                    int gb = hRECO->GetBin(j,k,m);
                    aDIFF_reco[i]->SetBinContent(gb, (aRECO[i]->GetBinContent(gb)-hRECO->GetBinContent(gb)));
                    aDIFF_gen[i]->SetBinContent(gb, (aGEN[i]->GetBinContent(gb)-hGEN->GetBinContent(gb)));
                    aDIFF_eff[i]->SetBinContent(gb, (aEFF[i]->GetBinContent(gb)-hEFF->GetBinContent(gb)));

                }
            }
        }
    }

    cout<<"--- efficiencies computed ---"<<endl;

    int range = 100;

    TH1D *reco_diff_distribution = new TH1D("reco_diff_distribution","RECO bin diff", 2*range, -range, +range);
    TH1D *gen_diff_distribution = new TH1D("gen_diff_distribution","GEN bin diff", 2*range, -range, +range);
    TH1D *eff_diff_distribution = new TH1D("eff_diff_distribution","EFF bin diff", 1000, -0.1, 0.1);
    TH1D *eff_diff_rel_distribution = new TH1D("eff_diff_rel_distribution","EFF bin diff REL", 1000, -1, 1);

    int nbineff0 = 0;
    for(int i=0; i<nhist; ++i){
        for(int j=1; j<=nbinsx; ++j){
            for(int k=1; k<=nbinsy; ++k){
                for(int m=1; m<=nbinsz; ++m){
                    int gb = hRECO->GetBin(j,k,m);
                    reco_diff_distribution->Fill(aDIFF_reco[i]->GetBinContent(gb));
                    gen_diff_distribution->Fill(aDIFF_gen[i]->GetBinContent(gb));
                    eff_diff_distribution->Fill(aDIFF_eff[i]->GetBinContent(gb));
                    if(hEFF->GetBinContent(gb) != 0) eff_diff_rel_distribution->Fill(aDIFF_eff[i]->GetBinContent(gb)/hEFF->GetBinContent(gb));
                    else if(i==0) nbineff0++; //count only one time
                }
            }
        }
    }

    cout<<" --- eff bin with zero content = "<<nbineff0<<" out of "<<nbinsx*nbinsy*nbinsz<<endl;


    TCanvas *c = new TCanvas();
    c->Divide(2,2);

    c->cd(1);
    reco_diff_distribution->Draw();

    c->cd(2);
    gen_diff_distribution->Draw();

    c->cd(3);
    eff_diff_distribution->Draw();


    c->cd(4);
    eff_diff_rel_distribution->Draw();


    TCanvas *c2 = new TCanvas();
    c2->Divide(3,3);

    TH1D *hRECO_X = hRECO->ProjectionX();
    TH1D *hRECO_Y = hRECO->ProjectionY();
    TH1D *hRECO_Z = hRECO->ProjectionZ();

    TH1D *hGEN_X = hGEN->ProjectionX();
    TH1D *hGEN_Y = hGEN->ProjectionY();
    TH1D *hGEN_Z = hGEN->ProjectionZ();

    TH1D *hEFF_X = hEFF->ProjectionX();
    TH1D *hEFF_Y = hEFF->ProjectionY();
    TH1D *hEFF_Z = hEFF->ProjectionZ();

    c2->cd(1);
    hRECO_X->Draw();

    c2->cd(2);
    hRECO_Y->Draw();

    c2->cd(3);
    hRECO_Z->Draw();

    c2->cd(4);
    hGEN_X->Draw();

    c2->cd(5);
    hGEN_Y->Draw();

    c2->cd(6);
    hGEN_Z->Draw();

    c2->cd(7);
    hEFF_X->Draw();

    c2->cd(8);
    hEFF_Y->Draw();

    c2->cd(9);
    hEFF_Z->Draw();
}