TString oldfilename = "../fittree_ntuBsDG0MC2018.root";
TString newfilename = "BsDG02018_kaonPtSF.root";
TString SFfilename = "kaonPtSF18.root";

void addSF(){

    auto oldFile = TFile::Open(oldfilename);
    auto tree = (TTree*)oldFile->Get("treeFit");

    auto sfFile = TFile::Open(SFfilename);
    auto hSF = (TH1D*)sfFile->Get("hSF");

    Float_t kaonppt, kaonmpt;

    tree->SetBranchAddress("kaonppt",&kaonppt);
    tree->SetBranchAddress("kaonmpt",&kaonmpt);
    
    auto newFile = TFile::Open(newfilename, "RECREATE");
    auto weightTree = new TNtuple("weightTree", "", "weight");

    for (int iEntry = 0; iEntry < tree->GetEntries(); iEntry++ )
    {
	tree->GetEntry(iEntry);

        double sf1 = hSF->GetBinContent(hSF->FindBin(kaonppt));
	double sf2 = hSF->GetBinContent(hSF->FindBin(kaonmpt));

        auto weight = sf1*sf2;
        weightTree->Fill(weight);
    }
    weightTree->Write();
    delete newFile;

}
