void skim(TString name)
{
   TFile oldfile(name + ".root");
   TTree *oldtree;
   oldfile.GetObject("OutTree", oldtree);

   TFile newfile(name + "_skim.root", "recreate");
   auto newtree = oldtree->CopyTree("HLT_JpsiMu == 1");
   newtree->Print();
   newfile.Write();
}