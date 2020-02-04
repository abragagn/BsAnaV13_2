void testEff()
{
    TFile *f = new TFile("/afs/cern.ch/user/a/abragagn/public/BsJpsiPhiDatasets/fittree_ntuBsDG0MC2018.root");
    TFile *fGen = new TFile("/afs/cern.ch/user/a/abragagn/public/BsJpsiPhiDatasets/ntuBsGEN.root");

    TTree *t = (TTree*)f->Get("treeFit");
    TTree *tGEN = (TTree*)fGen->Get("OutTree");

    TH1D *hcostheta = new TH1D("hcostheta","", 100, -1, 1);
    TH1D *hcostheta_1 = new TH1D("hcostheta_1","", 100, -1, 1);
    TH1D *hcostheta_GEN = new TH1D("hcostheta_GEN","", 100, -1, 1);
    TH1D *hcostheta_GEN_1 = new TH1D("hcostheta_GEN_1","", 100, -1, 1);

    TCanvas *c = new TCanvas();

    t->Draw("BscosthetaMC>>hcostheta");
    tGEN->Draw("angle_costheta_GEN>>hcostheta_GEN");

    auto hcostheta_eff = (TH1D*)hcostheta->Clone("hcostheta_eff");
    hcostheta_eff->Sumw2();
    hcostheta_eff->Divide(hcostheta_GEN);

    hcostheta_1->FillRandom(hcostheta, gRandom->Poisson(hcostheta->Integral()));
    hcostheta_GEN_1->FillRandom(hcostheta_GEN, gRandom->Poisson(hcostheta_GEN->Integral()));

    auto hcostheta_eff_1 = (TH1D*)hcostheta_1->Clone("hcostheta_eff_1");
    hcostheta_eff_1->Sumw2();
    hcostheta_eff_1->Divide(hcostheta_GEN_1);

    hcostheta_1->SetLineColor(kRed);
    hcostheta_GEN_1->SetLineColor(kRed);
    hcostheta_eff_1->SetLineColor(kRed);

    gPad->Clear();
    c->Divide(2,2);

    c->cd(1);
    hcostheta->Draw();
    hcostheta_1->Draw("SAME");

    c->cd(2);
    hcostheta_GEN->Draw();
    hcostheta_GEN_1->Draw("SAME");

    c->cd(3);
    hcostheta_eff->Draw();
    hcostheta_eff_1->Draw("SAME");

    c->cd(4);

}