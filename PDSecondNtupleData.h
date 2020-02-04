#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset() { autoReset(); }

PDSecondNtupleData() {


}
virtual ~PDSecondNtupleData() {
}

void initTree() {
  treeName = "OutTree";

  setBranch("run", &Trun, "run/F", &b_Trun);
  setBranch("evt", &Tevt, "evt/l", &b_Tevt);
  setBranch("lumi", &Tlumi, "lumi/F", &b_Tlumi);

  setBranch("Tag", &Ttag, "Tag/I", &b_Ttag);
  setBranch("MisTag", &Tmtag, "MisTag/F", &b_Tmtag);
  setBranch("MisTagCal", &TmtagCal, "MisTagCal/F", &b_TmtagCal);
  setBranch("MisTagCalBs", &TmtagCalBs, "MisTagCalBs/F", &b_TmtagCalBs);

  setBranch("isBs", &TisBs, "isBs/I", &b_TisBs);
  
  setBranch("HLT_JpsiTkTk", &THLT_Jtktk, "HLT_JpsiTkTk/I", &b_THLT_Jtktk);
  setBranch("HLT_JpsiTk", &THLT_Jtk, "HLT_JpsiTk/I", &b_THLT_Jtk);
  setBranch("HLT_JpsiMu", &THLT_Jmu, "HLT_JpsiMu/I", &b_THLT_Jmu);

  setBranch("HLT_MatchedJpsi", &TmatchHLTmu, "HLT_MatchedJpsi/I", &b_TmatchHLTmu);
  setBranch("HLT_MatchedTracks", &TmatchHLTk, "HLT_MatchedTracks/I", &b_TmatchHLTk);

  setBranch("MuonsPassSoftSelection", &TmuonsPassSoft, "MuonsPassSoftSelection/I", &b_TmuonsPassSoft);

  //polarization angles
  setBranch("angle_cospsi", &Tangle_cospsi, "angle_cospsi/F", &b_Tangle_cospsi);
  setBranch("angle_costheta", &Tangle_costheta, "angle_costheta/F", &b_Tangle_costheta);
  setBranch("angle_phi", &Tangle_phi, "angle_phi/F", &b_Tangle_phi);

  setBranch("PV_cosPoint", &TcosPoint, "PV_cosPoint/F", &b_TcosPoint);
  
  setBranch("ctau", &TLxy, "ctau/F", &b_TLxy);
  setBranch("ctauPD", &TLxyPD, "ctauPD/F", &b_TLxyPD);
  setBranch("ctauANBS", &TLxyANBS, "ctauANBS/F", &b_TLxyANBS);
  setBranch("ctauAT", &TLxyAT, "ctauAT/F", &b_TLxyAT);
  setBranch("ctauCNBS", &TLxyCNBS, "ctauCNBS/F", &b_TLxyCNBS);
  setBranch("ctauCO", &TLxyCO, "ctauCO/F", &b_TLxyCO);
  setBranch("ctau3D", &TLxyz, "ctau3D/F", &b_TLxyz);

  setBranch("ctauErr", &TLxyErr, "ctauErr/F", &b_TLxyErr);
  setBranch("ctau3DErr", &TLxyzErr, "ctau3DErr/F", &b_TLxyzErr);
  
  setBranch("B_VProb", &TVprob, "B_VProb/F", &b_TVprob);
  setBranch("Jpsi_VProb", &TVprobJpsi, "Jpsi_VProb/F", &b_TVprobJpsi);
  setBranch("PhiKstar_VProb", &TVprobPhi, "PhiKstar_VProb/F", &b_TVprobPhi);
  
  setBranch("B_Mass", &TBsMass, "B_Mass/F", &b_TBsMass);
  setBranch("B_Pt", &TBsPt, "B_Pt/F", &b_TBsPt);
  setBranch("B_Eta", &TBsEta, "B_Eta/F", &b_TBsEta);
  setBranch("B_Phi", &TBsPhi, "B_Phi/F", &b_TBsPhi);

  setBranch("B_MassNoRefit", &TBsMassNoRefit, "B_MassNoRefit/F", &b_TBsMassNoRefit);
  setBranch("Jpsi_MassNoRefit", &TJpsiMassNoRefit, "Jpsi_MassNoRefit/F", &b_TJpsiMassNoRefit);
  setBranch("PhiKstar_MassNoRefit", &TPhiMassNoRefit, "PhiKstar_MassNoRefit/F", &b_TPhiMassNoRefit);

  setBranch("B_MassFromSV", &TBsMassFromSV, "B_MassFromSV/F", &b_TBsMassFromSV);
  setBranch("Jpsi_MassFromSV", &TJpsiMassFromSV, "Jpsi_MassFromSV/F", &b_TJpsiMassFromSV);
  setBranch("PhiKstar_MassFromSV", &TPhiMassFromSV, "PhiKstar_MassFromSV/F", &b_TPhiMassFromSV);

  setBranch("Jpsi_Mass", &TJpsiMass, "Jpsi_Mass/F", &b_TJpsiMass);
  setBranch("Jpsi_Pt", &TJpsiPt, "Jpsi_Pt/F", &b_TJpsiPt);
  setBranch("Jpsi_Eta", &TJpsiEta, "Jpsi_Eta/F", &b_TJpsiEta);
  setBranch("Jpsi_Phi", &TJpsiPhi, "Jpsi_Phi/F", &b_TJpsiPhi);

  setBranch("PhiKstar_Mass", &TPhiMass, "PhiKstar_Mass/F", &b_TPhiMass);
  setBranch("PhiKstar_Pt", &TPhiPt, "PhiKstar_Pt/F", &b_TPhiPt);
  setBranch("PhiKstar_Eta", &TPhiEta, "PhiKstar_Eta/F", &b_TPhiEta);
  setBranch("PhiKstar_Phi", &TPhiPhi, "PhiKstar_Phi/F", &b_TPhiPhi);
  
  setBranch("Mup_Pt", &TMupPt, "Mup_Pt/F", &b_TMupPt);
  setBranch("Mup_Eta", &TMupEta, "Mup_Eta/F", &b_TMupEta);
  setBranch("Mup_Phi", &TMupPhi, "Mup_Phi/F", &b_TMupPhi);
  setBranch("Mup_Hits", &TMupHits, "Mup_Hits/F", &b_TMupHits);
  setBranch("Mup_HltPt", &TMupHltPt, "Mup_HltPt/F", &b_TMupHltPt);

  setBranch("Mum_Pt", &TMumPt, "Mum_Pt/F", &b_TMumPt);
  setBranch("Mum_Eta", &TMumEta, "Mum_Eta/F", &b_TMumEta);
  setBranch("Mum_Phi", &TMumPhi, "Mum_Phi/F", &b_TMumPhi);
  setBranch("Mum_Hits", &TMumHits, "Mum_Hits/F", &b_TMumHits);
  setBranch("Mum_HltPt", &TMumHltPt, "Mum_HltPt/F", &b_TMumHltPt);
  
  setBranch("KmK_Pt", &TKmPt, "KmK_Pt/F", &b_TKmPt);
  setBranch("KmK_Eta", &TKmEta, "KmK_Eta/F", &b_TKmEta);
  setBranch("KmK_Phi", &TKmPhi, "KmK_Phi/F", &b_TKmPhi);
  setBranch("KmK_Hits", &TKmHits, "KmK_Hits/F", &b_TKmHits);
  
  setBranch("KpPi_Pt", &TKpPt, "KpPi_Pt/F", &b_TKpPt);
  setBranch("KpPi_Eta", &TKpEta, "KpPi_Eta/F", &b_TKpEta);
  setBranch("KpPi_Phi", &TKpPhi, "KpPi_Phi/F", &b_TKpPhi);
  setBranch("KpPi_Hits", &TKpHits, "KpPi_Hits/F", &b_TKpHits);

  setBranch("Pi_Pt", &TPiPt, "Pi_Pt/F", &b_TPiPt);
  setBranch("Pi_Eta", &TPiEta, "Pi_Eta/F", &b_TPiEta);
  setBranch("Pi_Phi", &TPiPhi, "Pi_Phi/F", &b_TPiPhi);
  setBranch("Pi_Hits", &TPiHits, "Pi_Hits/F", &b_TPiHits);
  
  setBranch("Kaon_Deltadz", &TKaon_Deltadz, "Kaon_Deltadz/F", &b_TKaon_Deltadz);
  setBranch("Kaon_DeltaR", &TKaon_DeltaR, "Kaon_DeltaR/F", &b_TKaon_DeltaR);

  setBranch("N_PV", &TNpv, "N_PV/I", &b_TNpv);

  setBranch("Bc_Mass", &TBcMass, "Bc_Mass/F", &b_TBcMass);

  // Gen

  setBranch("MC_genFound", &TgenFound, "MC_genFound/I", &b_TgenFound);

  setBranch("MC_Flavour", &TgenFlavour, "MC_Flavour/I", &b_TgenFlavour);
  setBranch("MC_matched", &Tmatched, "MC_matched/I", &b_Tmatched);

  setBranch("angle_cospsi_GEN", &Tangle_cospsiGen, "angle_cospsi_GEN/F", &b_Tangle_cospsiGen);
  setBranch("angle_costheta_GEN", &Tangle_costhetaGen, "angle_costheta_GEN/F", &b_Tangle_costhetaGen);
  setBranch("angle_phi_GEN", &Tangle_phiGen, "angle_phi_GEN/F", &b_Tangle_phiGen);

  setBranch("ctau_GEN", &TLxyGen, "ctau_GEN/F", &b_TLxyGen);
  setBranch("ctau3D_GEN", &TLxyzGen, "ctau3D_GEN/F", &b_TLxyzGen);

  setBranch("B_Mass_GEN", &TBsMassGen, "B_Mass_GEN/F", &b_TBsMassGen);
  setBranch("B_Pt_GEN", &TBsPtGen, "B_Pt_GEN/F", &b_TBsPtGen);
  setBranch("B_Eta_GEN", &TBsEtaGen, "B_Eta_GEN/F", &b_TBsEtaGen);
  setBranch("B_Phi_GEN", &TBsPhiGen, "B_Phi_GEN/F", &b_TBsPhiGen);

  setBranch("Jpsi_Mass_GEN", &TJpsiMassGen, "Jpsi_Mass_GEN/F", &b_TJpsiMassGen);
  setBranch("Jpsi_Pt_GEN", &TJpsiPtGen, "Jpsi_Pt_GEN/F", &b_TJpsiPtGen);
  setBranch("Jpsi_Eta_GEN", &TJpsiEtaGen, "Jpsi_Eta_GEN/F", &b_TJpsiEtaGen);
  setBranch("Jpsi_Phi_GEN", &TJpsiPhiGen, "Jpsi_Phi_GEN/F", &b_TJpsiPhiGen);

  setBranch("PhiKstar_Mass_GEN", &TPhiMassGen, "PhiKstar_Mass_GEN/F", &b_TPhiMassGen);
  setBranch("PhiKstar_Pt_GEN", &TPhiPtGen, "PhiKstar_Pt_GEN/F", &b_TPhiPtGen);
  setBranch("PhiKstar_Eta_GEN", &TPhiEtaGen, "PhiKstar_Eta_GEN/F", &b_TPhiEtaGen);
  setBranch("PhiKstar_Phi_GEN", &TPhiPhiGen, "PhiKstar_Phi_GEN/F", &b_TPhiPhiGen);

  setBranch("Mup_Pt_GEN", &TMupPtGen, "Mup_Pt_GEN/F", &b_TMupPtGen);
  setBranch("Mup_Eta_GEN", &TMupEtaGen, "Mup_Eta_GEN/F", &b_TMupEtaGen);
  setBranch("Mup_Phi_GEN", &TMupPhiGen, "Mup_Phi_GEN/F", &b_TMupPhiGen);

  setBranch("Mum_Pt_GEN", &TMumPtGen, "Mum_Pt_GEN/F", &b_TMumPtGen);
  setBranch("Mum_Eta_GEN", &TMumEtaGen, "Mum_Eta_GEN/F", &b_TMumEtaGen);
  setBranch("Mum_Phi_GEN", &TMumPhiGen, "Mum_Phi_GEN/F", &b_TMumPhiGen);

  setBranch("KmK_Pt_GEN", &TKmPtGen, "KmK_Pt_GEN/F", &b_TKmPtGen);
  setBranch("KmK_Eta_GEN", &TKmEtaGen, "KmK_Eta_GEN/F", &b_TKmEtaGen);
  setBranch("KmK_Phi_GEN", &TKmPhiGen, "KmK_Phi_GEN/F", &b_TKmPhiGen);

  setBranch("KpPi_Pt_GEN", &TKpPtGen, "KpPi_Pt_GEN/F", &b_TKpPtGen);
  setBranch("KpPi_Eta_GEN", &TKpEtaGen, "KpPi_Eta_GEN/F", &b_TKpEtaGen);
  setBranch("KpPi_Phi_GEN", &TKpPhiGen, "KpPi_Phi_GEN/F", &b_TKpPhiGen);


}

float TVprob, TVprobJpsi, TVprobPhi;
float Tangle_cospsiGen,  Tangle_costhetaGen, Tangle_phiGen;
float Tangle_cospsi,  Tangle_costheta, Tangle_phi;
float TLxy, TLxyPD, TLxyANBS, TLxyAT, TLxyCNBS, TLxyCO;
float  TLxyz, TLxyErr, TLxyzErr, TLxyGen, TLxyzGen;
int TnBs, Tmatched;
//the Bs
float TBsMass, TBsPt, TBsPhi, TBsEta, TBsMassGen;
float TBsPtGen, TBsPhiGen, TBsEtaGen;
float TBsMassNoRefit, TJpsiMassNoRefit, TPhiMassNoRefit;
float TBsMassFromSV, TJpsiMassFromSV, TPhiMassFromSV;
//Jpsi and Phi
float TJpsiMass, TJpsiPt, TJpsiPhi, TJpsiEta, TJpsiMassGen;
float TJpsiPtGen, TJpsiPhiGen, TJpsiEtaGen;
float TPhiMass, TPhiPt, TPhiPhi, TPhiEta, TPhiMassGen;
float TPhiPtGen, TPhiPhiGen, TPhiEtaGen;
//Muons and Kaons
float TMupPt, TMupPhi, TMupEta, TMupHits, TMupPtGen, TMupPhiGen, TMupEtaGen;
float TMumPt, TMumPhi, TMumEta, TMumHits, TMumPtGen, TMumPhiGen, TMumEtaGen;
float TKpPt, TKpPhi, TKpEta, TKpHits, TKpPtGen, TKpPhiGen, TKpEtaGen;
float TKmPt, TKmPhi, TKmEta, TKmHits, TKmPtGen, TKmPhiGen, TKmEtaGen;
float TPiPt, TPiPhi, TPiEta, TPiHits, TBcMass;
float TcosPoint;
int   TNpv;
float TMupHltPt, TMumHltPt;
//event
float Tmtag, TmtagCal, TmtagCalBs, TmtagErr;
float Trun, Tlumi;
int   THLT_Jtktk, THLT_Jmu, THLT_Jtk, TisBs;  
int   TgenFlavour,Ttag;
int   TmatchHLTmu, TmatchHLTk;
int   TgenFound, TmuonsPassSoft;
ULong64_t Tevt;
float TKaon_Deltadz,TKaon_DeltaR;


//BRANCHES
TBranch *b_TVprob, *b_TVprobJpsi, *b_TVprobPhi, *b_Tmatched;
TBranch *b_Tangle_cospsiGen,  *b_Tangle_costhetaGen, *b_Tangle_phiGen;
TBranch *b_Tangle_cospsi,  *b_Tangle_costheta, *b_Tangle_phi;
TBranch *b_TLxy, *b_TLxyPD, *b_TLxyANBS, *b_TLxyAT, *b_TLxyCNBS, *b_TLxyCO;
TBranch *b_TLxyz, *b_TLxyErr, *b_TLxyzErr, *b_TLxyGen, *b_TLxyzGen;
TBranch *b_TnBs;
//the Bs
TBranch *b_TBsMass, *b_TBsPt, *b_TBsPhi, *b_TBsEta, *b_TBsMassGen;
TBranch *b_TBsPtGen, *b_TBsPhiGen, *b_TBsEtaGen;
TBranch *b_TBsMassNoRefit, *b_TJpsiMassNoRefit, *b_TPhiMassNoRefit;
TBranch *b_TBsMassFromSV, *b_TJpsiMassFromSV, *b_TPhiMassFromSV;
//Jpsi and Phi
TBranch *b_TJpsiMass, *b_TJpsiPt, *b_TJpsiPhi, *b_TJpsiEta, *b_TJpsiMassGen;
TBranch *b_TJpsiPtGen, *b_TJpsiPhiGen, *b_TJpsiEtaGen;
TBranch *b_TPhiMass, *b_TPhiPt, *b_TPhiPhi, *b_TPhiEta, *b_TPhiMassGen;
TBranch *b_TPhiPtGen, *b_TPhiPhiGen, *b_TPhiEtaGen;
//Muons and Kaons
TBranch *b_TMupPt, *b_TMupPhi, *b_TMupEta, *b_TMupHits, *b_TMupPtGen, *b_TMupPhiGen, *b_TMupEtaGen;
TBranch *b_TMumPt, *b_TMumPhi, *b_TMumEta, *b_TMumHits, *b_TMumPtGen, *b_TMumPhiGen, *b_TMumEtaGen;
TBranch *b_TKpPt, *b_TKpPhi, *b_TKpEta, *b_TKpHits, *b_TKpPtGen, *b_TKpPhiGen, *b_TKpEtaGen;
TBranch *b_TKmPt, *b_TKmPhi, *b_TKmEta, *b_TKmHits, *b_TKmPtGen, *b_TKmPhiGen, *b_TKmEtaGen;
TBranch *b_TPiPt, *b_TPiPhi, *b_TPiEta, *b_TPiHits, *b_TNpv, *b_TBcMass;
TBranch *b_TcosPoint;
TBranch *b_TMupHltPt, *b_TMumHltPt;
//event
TBranch *b_Tmtag, *b_TmtagCal, *b_TmtagCalBs, *b_TmtagErr;
TBranch *b_Trun, *b_Tlumi, *b_THLT_Jtktk, *b_THLT_Jmu, *b_THLT_Jtk, *b_TisBs;  
TBranch *b_TgenFlavour, *b_Ttag;
TBranch *b_TmatchHLTmu, *b_TmatchHLTk;
TBranch *b_TgenFound, *b_TmuonsPassSoft;
TBranch *b_Tevt;
TBranch *b_TKaon_Deltadz, *b_TKaon_DeltaR;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

