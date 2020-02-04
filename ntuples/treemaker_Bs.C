//Authors: Alibordi, Giacomo Fedi, Alberto Bragagnolo
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include <TMath.h>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBuffer.h>

using namespace std;

void treemaker_Bs(TString name = ""){

    TFile *fIn1 = new TFile(name);

    if(!fIn1) return;
    TTree* t = (TTree*)fIn1->Get("OutTree");
    Int_t n_entries = t->GetEntries();

    Int_t isBs, MuonsPassSoftSelection, HLT_JpsiMu, HLT_JpsiTkTk, MC_Flavour, Tag, N_PV;
    Float_t angle_cospsi, angle_costheta, angle_phi, ctau, ctauErr, B_VProb, B_Mass;
    Float_t B_Pt, B_MassFromSV, Mum_Pt, Mum_Eta, Mum_Phi, Mup_Pt, Mup_Eta, Mup_Phi;
    Float_t KmK_Pt, KmK_Eta, KmK_Phi, KpPi_Pt, KpPi_Eta, KpPi_Phi; 
    Float_t angle_cospsi_GEN, angle_costheta_GEN, angle_phi_GEN, ctau_GEN, Jpsi_Mass;
    Float_t KpPi_Hits, KmK_Hits;
    Float_t PhiKstar_Mass, MisTag, MisTagCal, MisTagCalBs;
    Float_t Mup_HltPt, Mum_HltPt;
Float_t B_Eta;
    
    t->SetBranchAddress("isBs",&isBs);
    t->SetBranchAddress("HLT_JpsiMu",&HLT_JpsiMu);
    t->SetBranchAddress("HLT_JpsiTkTk",&HLT_JpsiTkTk);
    t->SetBranchAddress("MuonsPassSoftSelection",&MuonsPassSoftSelection);
    t->SetBranchAddress("angle_cospsi",&angle_cospsi);
    t->SetBranchAddress("angle_costheta",& angle_costheta);
    t->SetBranchAddress("angle_phi",&angle_phi);
    t->SetBranchAddress("ctau",&ctau);
    t->SetBranchAddress("ctauErr",&ctauErr);
    t->SetBranchAddress("B_VProb",&B_VProb);
    t->SetBranchAddress("B_Mass",&B_Mass);
    t->SetBranchAddress("B_Pt",&B_Pt);
    t->SetBranchAddress("B_Eta",&B_Eta);
    t->SetBranchAddress("B_MassFromSV",&B_MassFromSV);

    t->SetBranchAddress("Mum_Pt",&Mum_Pt);
    t->SetBranchAddress("Mum_Eta",&Mum_Eta);
    t->SetBranchAddress("Mum_Phi",&Mum_Phi);
    t->SetBranchAddress("Mup_Pt",&Mup_Pt);
    t->SetBranchAddress("Mup_Eta",&Mup_Eta);
    t->SetBranchAddress("Mup_Phi",&Mup_Phi);

    t->SetBranchAddress("Mup_HltPt",&Mup_HltPt);
    t->SetBranchAddress("Mum_HltPt",&Mum_HltPt);

    t->SetBranchAddress("KmK_Pt",&KmK_Pt);
    t->SetBranchAddress("KmK_Eta",&KmK_Eta);
    t->SetBranchAddress("KmK_Phi",& KmK_Phi);

    t->SetBranchAddress("KpPi_Pt",&KpPi_Pt);
    t->SetBranchAddress("KpPi_Eta",&KpPi_Eta);
    t->SetBranchAddress("KpPi_Phi",& KpPi_Phi);

    t->SetBranchAddress("N_PV",&N_PV);

    t->SetBranchAddress("angle_cospsi_GEN",&angle_cospsi_GEN);
    t->SetBranchAddress("angle_costheta_GEN",&angle_costheta_GEN);
    t->SetBranchAddress("angle_phi_GEN",&angle_phi_GEN);

    t->SetBranchAddress("ctau_GEN",&ctau_GEN);
    t->SetBranchAddress("Jpsi_Mass",&Jpsi_Mass);
    t->SetBranchAddress("PhiKstar_Mass",&PhiKstar_Mass);

    t->SetBranchAddress("KpPi_Hits",&KpPi_Hits);
    t->SetBranchAddress("KmK_Hits",&KmK_Hits);

    t->SetBranchAddress("MC_Flavour",&MC_Flavour);

    t->SetBranchAddress("Tag",&Tag);
    t->SetBranchAddress("MisTag",&MisTag);
    t->SetBranchAddress("MisTagCal",&MisTagCal);
    t->SetBranchAddress("MisTagCalBs",&MisTagCalBs);

       
    Float_t svmass, BsCt2DMC, BscosthetaMC, BscospsiMC, BsphiMC, BsCt2DMCErr, BsCt2DMC_GEN, BscosthetaMC_GEN, BscospsiMC_GEN, BsphiMC_GEN;
    Float_t mistag;
    Float_t BsPt, muonppt, muonmpt, muonpHLTpt, muonmHLTpt, kaonppt, kaonmpt;
    Float_t muonpeta, muonmeta, kaonpeta, kaonmeta;
    Float_t jpsimass, phimass;
    Int_t Bs_NPV, tag;   
    Float_t BsEta;
       
    TFile *f = new TFile("fittree_" + name,"RECREATE");

    TTree *treefit = new TTree("treeFit", "");

    treefit->Branch("svmass",&svmass,"svmass/F");
    treefit->Branch("BsPt",&BsPt,"BsPt/F");
    treefit->Branch("BsCt2DMC",&BsCt2DMC,"BsCt2DMC/F");
    treefit->Branch("BsCt2DMCErr",&BsCt2DMCErr,"BsCt2DMCErr/F");
    treefit->Branch("BscosthetaMC",&BscosthetaMC,"BscosthetaMC/F");
    treefit->Branch("BscospsiMC",&BscospsiMC,"BscospsiMC/F");
    treefit->Branch("BsphiMC",&BsphiMC,"BsphiMC/F");
    treefit->Branch("BsCt2DMC_GEN",&BsCt2DMC_GEN,"BsCt2DMC_GEN/F");
    treefit->Branch("BscosthetaMC_GEN",&BscosthetaMC_GEN,"BscosthetaMC_GEN/F");
    treefit->Branch("BscospsiMC_GEN",&BscospsiMC_GEN,"BscospsiMC_GEN/F");
    treefit->Branch("BsphiMC_GEN",&BsphiMC_GEN,"BsphiMC_GEN/F");
    treefit->Branch("tag",&tag,"tag/I");
    treefit->Branch("mistag", &mistag, "mistag/F");


    //selection variables
    treefit->Branch("muonmHLTpt", &muonmHLTpt, "muonmHLTpt/F");
    treefit->Branch("muonpHLTpt", &muonpHLTpt, "muonpHLTpt/F");

    treefit->Branch("muonmpt", &muonmpt, "muonmpt/F");
    treefit->Branch("muonppt", &muonppt, "muonppt/F");
    treefit->Branch("kaonppt", &kaonppt, "kaonppt/F");
    treefit->Branch("kaonmpt", &kaonmpt, "kaonmpt/F");

    treefit->Branch("muonmeta", &muonmeta, "muonmeta/F");
    treefit->Branch("muonpeta", &muonpeta, "muonpeta/F");
    treefit->Branch("kaonpeta", &kaonpeta, "kaonpeta/F");
    treefit->Branch("kaonmeta", &kaonmeta, "kaonmeta/F");

    treefit->Branch("jpsimass", &jpsimass, "jpsimass/F");
    treefit->Branch("phimass", &phimass, "phimass/F");

    treefit->Branch("Bs_NPV",&Bs_NPV,"Bs_NPV/I");
    treefit->Branch("BsEta",&BsEta,"BsEta/F");

    Int_t npass = 0;

    for (Int_t i=0;i<n_entries;i++) {

        if(i%10000 == 0) cout<<i<<" / "<<n_entries<<endl;

        t->GetEntry(i);

        if( B_MassFromSV > 5.24 && B_MassFromSV < 5.49
            && Mum_Pt > 3.5 && Mup_Pt > 3.5 && B_Pt > 11 
            && KmK_Pt > 1.2 && KpPi_Pt > 1.2 
            && fabs(Mum_Eta) < 2.5 && fabs(Mup_Eta) < 2.5
            && fabs(KmK_Eta) < 2.5 && fabs(KpPi_Eta) < 2.5
            //&& ctau > 0.007
            && abs(Jpsi_Mass-3.0969) < 0.150 
            && abs(PhiKstar_Mass-1.01946) < 0.010 
            && B_VProb > 0.001
            && HLT_JpsiMu == 1
            && ((int(KpPi_Hits)/100)%10000)%100 >= 4 
            && ((int(KmK_Hits)/100)%10000)%100 >= 4 
            && MuonsPassSoftSelection == 1){

            if(TMath::IsNaN(ctauErr)) continue;
            Float_t w = MisTagCal;
            Float_t correctedMistag =  w;
            if(abs(Tag) != 1 && w != 0.5) correctedMistag = 0.5;
            else correctedMistag = w;

            if(correctedMistag > 1.) correctedMistag = 1.;
            if(correctedMistag < 0.) correctedMistag = 0.;

            svmass          = B_MassFromSV;                                                                                               
            BsPt            = B_Pt;                                                                                               
            BsCt2DMC        = ctau;
            BsCt2DMCErr     = ctauErr;
            BscosthetaMC    = angle_costheta;
            BscospsiMC      = angle_cospsi;
            BsphiMC         = angle_phi;
            Bs_NPV          = (int)N_PV;
            BscosthetaMC_GEN= angle_costheta_GEN;
            BscospsiMC_GEN  = angle_cospsi_GEN;
            BsphiMC_GEN     = angle_phi_GEN;
            BsCt2DMC_GEN    = ctau_GEN;
            mistag          = correctedMistag;
            tag             = (int)Tag;

            muonppt     = Mup_Pt;
            muonmpt     = Mum_Pt;
            kaonppt     = KpPi_Pt;
            kaonmpt     = KmK_Pt;
            muonpeta    = Mup_Eta;
            muonmeta    = Mum_Eta;
            kaonpeta    = KpPi_Eta;
            kaonmeta    = KmK_Eta;
            jpsimass    = Jpsi_Mass;
            phimass     = PhiKstar_Mass;

            muonpHLTpt  = Mup_HltPt;
            muonmHLTpt  = Mum_HltPt;

            BsEta = B_Eta;


            treefit->Fill();

            npass++;

        }

    }

    cout<<endl<<"nselected "<<npass<<" (eff "<<(float)npass/(float)n_entries<<")"<<endl;

    treefit->Print();
    f->Write();
}
