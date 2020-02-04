void skimToDataStat(TString name, int nEvt)
{
   TFile oldfile(name + ".root");
   TTree *oldtree;
   oldfile.GetObject("treeFit", oldtree);


   TFile newfile(name + "_DataStat.root", "recreate");
   auto newtree = oldtree->CloneTree(nEvt);
   newtree->Print();
   newfile.Write();
}