#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>


#include "TMatrixF.h"
#include "TVectorF.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TMatrixTBase.h"
#include "Math/SMatrix.h"

#include "PDAnalyzer.h"
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"

#include "PDSecondNtupleWriter.h"

using namespace std;

const float MUMASS  = 0.105658;
const float KMASS   = 0.493677;
const float BSMASS  = 5.36689;
const float PIMASS  = 0.139570;
const float B0MASS  = 5.27963;

// pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list hist.root -v outputFile ntu.root -v process Bs -v tagCalChannel BsJPsiPhiMC2017 -v histoMode RECREATE -v use_gen t -v genOnly f -v debug f -v use_tracks t -v use_hlts t -v use_pvts t -v use_hlto t -n 1000

PDAnalyzer::PDAnalyzer() {

  std::cout << "new PDAnalyzer" << std::endl;

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor

  setUserParameter( "verbose", "f" );
  setUserParameter( "outputFile", "ntu.root" );
  setUserParameter( "process", "Bs" );
  setUserParameter( "ptCut", "4.0" );
  setUserParameter( "maxAngle", "1.0" );
  setUserParameter( "genOnly", "f" ); //to be set to true to run on GEN Only samples
  setUserParameter( "debug", "f" ); // set to true to activate debugging printout
  setUserParameter( "tagCalChannel", "BuJPsiKData2018" );
  setUserParameter( "useBestPV", "f" );
}

PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

  PDAnalyzerUtil::beginJob();

  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly seture
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );
  getUserParameter( "ptCut", ptCut );
  getUserParameter( "maxAngle", maxAngle );
  getUserParameter( "genOnly", genOnly );
  getUserParameter( "debug", debug );
  getUserParameter( "outputFile", outputFile );
  getUserParameter( "process", process );
  getUserParameter( "tagCalChannel", tagCalChannel );
  getUserParameter( "useBestPV", useBestPV );
  
  cout<<"verbose "<<verbose<<endl;
  cout<<"ptCut "<<ptCut<<endl;
  cout<<"maxAngle "<<maxAngle<<endl;
  cout<<"genOnly "<<genOnly<<endl;
  cout<<"debug "<<debug<<endl;
  cout<<"outputFile "<<outputFile<<endl;
  cout<<"process "<<process<<endl;
  cout<<"tagCalChannel "<<tagCalChannel<<endl;
  cout<<"askHltMatch "<<askHltMatch<<endl;
  cout<<"useBestPV "<<useBestPV<<endl;


  if(genOnly && !use_gen){cout<<" !!!!------------------- genOnly && !use_gen -----------------!!!! "<<endl;}

  tWriter = new PDSecondNtupleWriter; // second ntuple
  tWriter->open( getUserParameter("outputFile"), "RECREATE" );

  if(!genOnly){
    inizializeMuonMvaReader();
    inizializeOSMuonMvaReader();
    bool osInit = inizializeOSMuonCalibration(tagCalChannel, "BuJPsiKMC2018", "BsJPsiPhiMC2018");
    if(!osInit) cout<<endl<<"!!! FAILED TO INIZIALIZED TAG CALIBRATION"<<endl<<endl;
  }

  for(int i=0;i<5;++i) counter[i] = 0;

  return;
}


void PDAnalyzer::book() {

  autoSavedObject =
    hPhi          = new TH1F( "Phi mass", "Phi mass", 240, 0.95, 1.1 );
  
  autoSavedObject =
    hJpsi         = new TH1F( "Jpsi mass", "Jpsi mass", 240, 2.8, 3.3 );

  autoSavedObject =
    hPhiHLT       = new TH1F( "Phi mass HLT", "Phi mass HLT", 240, 0.95, 1.1 );
  
  autoSavedObject =
    hJpsiHLT      = new TH1F( "Jpsi mass HLT", "Jpsi mass HLT", 240, 2.8, 3.3 );
  
  autoSavedObject =
    hBs           = new TH1F( "Bs mass", "Bs mass", 50, 5.2, 5.5 );
  
  autoSavedObject =
    hAnglePhi     = new TH1F( "Phi angle", "Phi Angle", 50, -3.2, 3.2 );
  
  autoSavedObject =
    hCosPsi       = new TH1F( "CosPsi", "CosPsi", 51, -1, 1 );
  
  autoSavedObject =
    hCosTheta     = new TH1F( "CosTheta", "CosTheta", 51, -1, 1 );

  return;
}


void PDAnalyzer::reset() {
  autoReset();
  return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

  if ( (!(event_tot%10) && event_tot<100 ) || 
  (!(event_tot %100) && event_tot<1000 ) || 
  (!(event_tot %1000)&& event_tot<10000 ) || 
  (!(event_tot %10000) && event_tot<100000 ) || 
  (!(event_tot %100000) && event_tot<1000000 ) || 
  (!(event_tot %1000000) && event_tot<10000000 ) )
      cout << " == at event " << event_file << " " << event_tot << endl;

  convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
  convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
  convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
  convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

  if(!genOnly) computeMuonVar();   //compute muon variable for soft id
  if(!genOnly) inizializeOsMuonTagVars();
  tWriter->Reset();

  if(debug) cout<<endl<<" --------------------- NEW EVENT "<<event_tot<<endl;

  TLorentzVector mu_p_Bus, mu_m_Bus, k_p_Bus, k_m_Bus; //needed in gen block so declared here

  int selectedPV, selectedPVgen, selectedSV;
  TVector3 vSV, pSV;
  // float selectedCosPoint;


  if(!genOnly){
    //----- SEARCH FOR CANDIDATES
    int iSV = -1;
    float selVprob = -1;
    int MuMatchHLT = 0, KMatchHLT = 0;
    int nCandidates = 0;

    for ( int isvt = 0; isvt < nSVertices; ++isvt ){
      if(svtType->at(isvt) != PDEnumString::svtBsJPsiPhi && process == "Bs") continue;
      if(svtType->at(isvt) != PDEnumString::svtBdJPsiKx && process == "Bd") continue;
      int mumatch = 0;
      int kmatch = 0;

      int iJPsi = (subVtxFromSV(isvt)).at(0);
      int iHad  = (subVtxFromSV(isvt)).at(1);
      const vector<int>& tkJpsi = tracksFromSV(iJPsi);
      const vector<int>& tkHad  = tracksFromSV(iHad);

      if(tkJpsi.size()!=2) continue;
      if(tkHad.size()!=2) continue;

      int m1 = -1, m2 = -1; //hlt index
      int k1 = -1, k2 = -1;
      if (use_hlto) {
        m1 = HLTMatch(tkJpsi[0],1);
        m2 = HLTMatch(tkJpsi[1],1);
        k1 = HLTMatch(tkHad[0],0);
        k2 = HLTMatch(tkHad[1],0);
        
        if (m1 != m2 && m1>=0 && m2>=0) {
          mumatch = 1;
          if(process=="Bs"){
            TLorentzVector hltm1, hltm2;
            hltm1.SetPtEtaPhiM(hltPt->at(m1),hltEta->at(m1),hltPhi->at(m1),MUMASS);
            hltm2.SetPtEtaPhiM(hltPt->at(m2),hltEta->at(m2),hltPhi->at(m2),MUMASS);
            hJpsiHLT->Fill((hltm1+hltm2).M());
          }
        }

        if (k1 != k2 && k1>=0 && k2>=0) {
          kmatch = 1;
          if(process=="Bs"){
            TLorentzVector hltk1, hltk2;
            hltk1.SetPtEtaPhiM(hltPt->at(k1),hltEta->at(k1),hltPhi->at(k1),KMASS);
            hltk2.SetPtEtaPhiM(hltPt->at(k2),hltEta->at(k2),hltPhi->at(k2),KMASS);
            hPhiHLT->Fill((hltk1+hltk2).M());
          }
        }
      }

      if(mumatch == 0 && askHltMatch == true) continue;

      float tmpVprob = TMath::Prob(svtChi2->at(isvt),svtNDOF->at(isvt));
      if(tmpVprob<selVprob) continue;

      iSV = isvt;
      selVprob = tmpVprob;
      MuMatchHLT = mumatch;
      KMatchHLT = kmatch;

      nCandidates++;

    }

    (tWriter->nCand) = nCandidates;

    if(iSV < 0) return false;

    selectedSV = iSV;

    int isbs = -1;
    if(process=="Bs") isbs = 1;
    else isbs = 0;

    (tWriter->TisBs) = isbs;
    (tWriter->TmatchHLTmu) = MuMatchHLT;
    (tWriter->TmatchHLTk) = KMatchHLT;

    if(debug) cout<<"  ++++++++++++++ Candidate found "<<isbs<<" ++++++++++++++ "<<endl;

    // ----- BEGINNING OF RECONSTRUCTION

    (tWriter->Trun) = runNumber;
    (tWriter->Tevt) = ULong64_t(eventNumber);
    (tWriter->Tlumi) = lumiSection;
    (tWriter->TNpv) = nPVertices;

    if (use_hlts){
      if(debug) cout << " ++++++++++++++ Checking HLT decision ++++++++++++++ " << endl;
      bool jpsimu = false;
      bool jpsitktk = false;
      bool jpsitk = false;

      if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
      if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;
      if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;

      (tWriter->THLT_Jmu) = jpsimu;
      (tWriter->THLT_Jtktk) = jpsitktk;
      (tWriter->THLT_Jtk) = jpsitk;

      if(!jpsimu) return false;
    }

    vector<int> mu_idx_refit, k_idx_refit;
    TLorentzVector mu_p, mu_m, k_p, k_m;
    TLorentzVector mu_p_norefit, mu_m_norefit, k_p_norefit, k_m_norefit;
    TLorentzVector Bs, Jpsi, Phi, BsNoRefit, JpsiNoRefit, PhiNoRefit, pi, Bc;
    TVector3 SVpos, PVpos;
    int mu_idx_p, mu_idx_m;

    float mupH = -999, mumH = -999, kpH = -999, kmH = -999;
      
    double angle_costheta;
    double angle_phi;
    double angle_cospsi;

    // const vector<int>& tks = tracksFromSV(iSV);
    int iJPsi = (subVtxFromSV(iSV)).at(0);
    if(svtType->at(iJPsi) != PDEnumString::svtJPsi){
      cout<<"(subVtxFromSV(iSV)).at(0) NOT A JPSI but "<<svtType->at(iJPsi)<<endl;
      return false;
    }
    int iHad  = (subVtxFromSV(iSV)).at(1);
    if(process == "Bs"){
      if(svtType->at(iHad) != PDEnumString::svtPhi){
        cout<<"(subVtxFromSV(iSV)).at(1) NOT A PHI but "<<svtType->at(iHad)<<endl;
        return false;
      }
    }
    if(process == "Bd"){
      if(svtType->at(iHad) != PDEnumString::svtKx0){
        cout<<"(subVtxFromSV(iSV)).at(1) NOT A K0* but "<<svtType->at(iHad)<<endl;
        return false;
      }
    }

    const vector<int>& mu_idx = tracksFromSV(iJPsi);
    const vector<int>& k_idx  = tracksFromSV(iHad);

    float candMass = svtMass->at(iSV);
    float jpsiMass = svtMass->at(iJPsi);
    float hadMass = svtMass->at(iHad);

    hBs->Fill(candMass);
    hJpsi->Fill(jpsiMass);
    hPhi->Fill(hadMass);
    
    (tWriter->TBsMassFromSV) = candMass;
    (tWriter->TJpsiMassFromSV) = jpsiMass;
    (tWriter->TPhiMassFromSV) = hadMass;

    float VProb = TMath::Prob(svtChi2->at(iSV),svtNDOF->at(iSV));
    float VProbJpsi = TMath::Prob(svtChi2->at(iJPsi),svtNDOF->at(iJPsi));
    float VProbPhi = TMath::Prob(svtChi2->at(iHad),svtNDOF->at(iHad));

    for(auto it:mu_idx)
      mu_idx_refit.push_back(vpIndex(it,iSV));

    for(auto it:k_idx)
      k_idx_refit.push_back(vpIndex(it,iSV));

    if (k_idx.size()!=2 || mu_idx.size()!=2 || k_idx_refit.size()!=2 || mu_idx_refit.size()!=2 ){
      cout<<" PROBLEMS WITH TRACKS.SIZE()"<<endl; 
      return false;
    }
    
    SVpos.SetXYZ(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));

    if ( (k_idx_refit[0] >= int(tvpPt->size())
      ||  k_idx_refit[1] >= int(tvpPt->size()))
      ||  (mu_idx_refit[0] >= int(tvpPt->size())
      ||  mu_idx_refit[1] >= int(tvpPt->size())) ){
      cout << "PROBLEM IN THE vpIndex FUNCTION in run " << runNumber << " event " << eventNumber << " lumi " << lumiSection << endl;
      cout << "Size of the full track collection " << nTracks << " N refitted " << tvpPt->size() << endl;
      cout << " refitted indices k1 " <<  k_idx_refit[0] << " k2 " <<  k_idx_refit[1] << " mu 1 " << mu_idx_refit[0] << " m2 " <<  mu_idx_refit[1] << endl;
      return false;
    }

    //analyze the muons 
    if (trkCharge->at(mu_idx[0]) > 0 && trkCharge->at(mu_idx[1]) < 0){
      mu_idx_p = mu_idx[0];
      mu_idx_m = mu_idx[1];

      mu_p_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[0]),trkEta->at(mu_idx[0]),trkPhi->at(mu_idx[0]),MUMASS);
      mu_m_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[1]),trkEta->at(mu_idx[1]),trkPhi->at(mu_idx[1]),MUMASS);

      mupH = trkHitPattern->at(mu_idx[0]);
      mumH = trkHitPattern->at(mu_idx[1]);

      mu_p.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[0]),tvpEta->at(mu_idx_refit[0]),tvpPhi->at(mu_idx_refit[0]),MUMASS);
      mu_m.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[1]),tvpEta->at(mu_idx_refit[1]),tvpPhi->at(mu_idx_refit[1]),MUMASS);
    } else if (trkCharge->at(mu_idx[0]) < 0 && trkCharge->at(mu_idx[1]) > 0){
      mu_idx_p = mu_idx[1];
      mu_idx_m = mu_idx[0];

      mu_p_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[1]),trkEta->at(mu_idx[1]),trkPhi->at(mu_idx[1]),MUMASS);
      mu_m_norefit.SetPtEtaPhiM(trkPt->at(mu_idx[0]),trkEta->at(mu_idx[0]),trkPhi->at(mu_idx[0]),MUMASS);

      mupH = trkHitPattern->at(mu_idx[1]);
      mumH = trkHitPattern->at(mu_idx[0]);

      mu_p.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[1]),tvpEta->at(mu_idx_refit[1]),tvpPhi->at(mu_idx_refit[1]),MUMASS);
      mu_m.SetPtEtaPhiM(tvpPt->at(mu_idx_refit[0]),tvpEta->at(mu_idx_refit[0]),tvpPhi->at(mu_idx_refit[0]),MUMASS);
    }else{
      cout<<" PROBLEMS WITH MUONS CHARGES "<< endl;
      return false;
    }

    // Muon Soft Selection Flag
    int muonSoft = 1;
    for(auto it:mu_idx){
      if(!isTrkHighPurity(it)) muonSoft = 0;
      if(abs(trkDz->at(it))>20) muonSoft = 0;
      if(abs(trkDxy->at(it))>0.3) muonSoft = 0;
      if((((int(trkHitPattern->at(it))/100)%10000)%100)<5) muonSoft = 0;
      if(MuonFromTrack(it) >= 0 ){
        bool goodMuTMO = ( muoType->at( MuonFromTrack(it) ) ) & PDEnumString::tmOneStation;
        if(!goodMuTMO) muonSoft = 0;
      }else{muonSoft = 0;}
    }

    // Hlt Pt
    float hltMuPt_p = -1, hltMuPt_m = -1;

    if(MuMatchHLT){
      int mp = HLTMatch(mu_idx_p,1);
      int mm = HLTMatch(mu_idx_m,1);
      
      hltMuPt_p = hltPt->at(mp);
      hltMuPt_m = hltPt->at(mm);
    }

    //analyze the traks
    if(process == "Bs"){
      if (trkCharge->at(k_idx[0]) > 0 && trkCharge->at(k_idx[1]) < 0){
        k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
        k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);

        kpH = trkHitPattern->at(k_idx[0]);
        kmH = trkHitPattern->at(k_idx[1]);

        k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
        k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
      } else if (trkCharge->at(k_idx[0]) < 0 && trkCharge->at(k_idx[1]) > 0){
        k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
        k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);

        kpH = trkHitPattern->at(k_idx[1]);
        kmH = trkHitPattern->at(k_idx[0]);    

        k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
        k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
      }else{
        cout<<" PROBLEMS WITH KAONS CHARGES "<< endl;
        return false;
      }

    }else if(process == "Bd"){
      TLorentzVector cand1, cand2;

      k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
      k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),PIMASS);
        
      cand1 = k_p_norefit+ k_m_norefit;
        
      k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),PIMASS);
      k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
        
      cand2 = k_p_norefit+ k_m_norefit;

      if (abs(cand1.M()-svtMass->at( iSV )) <= abs(cand2.M()-svtMass->at( iSV ))){
        k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),KMASS);
        k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),PIMASS);
        kpH = trkHitPattern->at(k_idx[0]);
        kmH = trkHitPattern->at(k_idx[1]);  
      }
      else {
        k_p_norefit.SetPtEtaPhiM(trkPt->at(k_idx[0]),trkEta->at(k_idx[0]),trkPhi->at(k_idx[0]),PIMASS);
        k_m_norefit.SetPtEtaPhiM(trkPt->at(k_idx[1]),trkEta->at(k_idx[1]),trkPhi->at(k_idx[1]),KMASS);
        kpH = trkHitPattern->at(k_idx[1]);
        kmH = trkHitPattern->at(k_idx[0]);    
      }
      //select the mass hypotesis closer to PDG value
      k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
      k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),PIMASS);
      cand1 = k_p + k_m;

      k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),PIMASS);
      k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
      cand2 = k_p + k_m;

      if (abs(cand1.M()-svtMass->at( iSV )) <= abs(cand2.M()-svtMass->at( iSV ))){
        k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),KMASS);
        k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),PIMASS);
      }
      else {
        k_p.SetPtEtaPhiM(tvpPt->at(k_idx_refit[0]),tvpEta->at(k_idx_refit[0]),tvpPhi->at(k_idx_refit[0]),PIMASS);
        k_m.SetPtEtaPhiM(tvpPt->at(k_idx_refit[1]),tvpEta->at(k_idx_refit[1]),tvpPhi->at(k_idx_refit[1]),KMASS);
      }
    }

    if (debug) {
      cout << " +++++++++++++ Building Candidates ++++++++++++++ " << endl;
      cout << "From refitted tracks: J/psi mass " << (mu_p+mu_m).Mag() << " Phi/KStar mass " << (k_p+k_m).Mag() << endl;
      cout << "From native tracks: J/psi mass " << (mu_p_norefit+mu_m_norefit).Mag() << " Phi/Kstar mass " << (k_p_norefit+k_m_norefit).Mag() << endl;
      cout << endl;
    }

    Jpsi = mu_p + mu_m;
    Phi = k_p + k_m;
    Bs = Jpsi + Phi;

    JpsiNoRefit = mu_p_norefit+mu_m_norefit;
    PhiNoRefit = k_p_norefit + k_m_norefit;
    BsNoRefit = JpsiNoRefit + PhiNoRefit;

    //selecting the PV and computing Lxy and Error
    int ipv_refitANBS = -1;
    int ipv_refitAT = -1;
    int ipv_refitCNBS = -1;
    int ipv_refitCO = -1;
    int pvIndex = -1;
    float cosPoint;

    FindPV(SVpos, PVpos, Bs, pvIndex, cosPoint);

    if(pvIndex < 0){
      cout<<"PV not found"<<endl;
      return false;
    }

    int ipv = svtPVtx->at(iSV);

    selectedPV = pvIndex;
    pSV = Bs.Vect();
    vSV = SVpos;
    // selectedCosPoint = cosPoint;

    (tWriter->TcosPoint) = cosPoint;

    for(int isub=0; isub<nCompVts; ++isub){
      if(subPart->at(isub)!=iSV) continue;
      int ivtx = subSVtx->at(isub);
      if(svtType->at(ivtx)==PDEnumString::refittedPVAllTrkNoBS) ipv_refitANBS = ivtx;
      if(svtType->at(ivtx)==PDEnumString::refittedPVAllTrk) ipv_refitAT = ivtx;
      if(svtType->at(ivtx)==PDEnumString::refittedPVRCOnlyNoBS) ipv_refitCNBS = ivtx;
      if(svtType->at(ivtx)==PDEnumString::refittedPVRCOnly) ipv_refitCO = ivtx;
    }

    TVector3 PVposPD;
    PVposPD.SetXYZ(pvtX->at(ipv),pvtY->at(ipv),pvtZ->at(ipv));
    TVector3 PVSVdistPD=SVpos-PVposPD;
    TVector3 PVposANBS;
    if (ipv_refitANBS>=0) PVposANBS.SetXYZ(svtX->at(ipv_refitANBS),svtY->at(ipv_refitANBS),svtZ->at(ipv_refitANBS));
    TVector3 PVSVdistANBS=SVpos-PVposANBS;
    TVector3 PVposAT;
    if (ipv_refitAT>=0) PVposAT.SetXYZ(svtX->at(ipv_refitAT),svtY->at(ipv_refitAT),svtZ->at(ipv_refitAT));
    TVector3 PVSVdistAT=SVpos-PVposAT;
    TVector3 PVposCNBS;
    if (ipv_refitCNBS>=0) PVposCNBS.SetXYZ(svtX->at(ipv_refitCNBS),svtY->at(ipv_refitCNBS),svtZ->at(ipv_refitCNBS));
    TVector3 PVSVdistCNBS=SVpos-PVposCNBS;
    TVector3 PVposCO;
    if (ipv_refitCO>=0) PVposCO.SetXYZ(svtX->at(ipv_refitCO),svtY->at(ipv_refitCO),svtZ->at(ipv_refitCO));
    TVector3 PVSVdistCO=SVpos-PVposCO;

    // tagging
    if(debug) cout<<"  ++++++++++++++ Tagging ++++++++++++++ "<<endl;

    setVtxOsMuonTag(iSV, pvIndex); //set PV and SVT for tagging class
    bool istagged = makeOsMuonTagging();

    (tWriter->Ttag) = getOsMuonTag(); //get tag decision -1, +1 or 0 

    if(istagged){
      (tWriter->Tmtag) = getOsMuonTagMistagProbRaw();
      (tWriter->TmtagCal) = getOsMuonTagMistagProbCalProcess();
      (tWriter->TmtagCalBs) = getOsMuonTagMistagProbCalProcessBuBs();
    }else{
      (tWriter->Tmtag) = 0.5;
      (tWriter->TmtagCal) = 0.5;
      (tWriter->TmtagCalBs) = 0.5;
    }

    // lifetime
    if(debug) cout<<"  ++++++++++++++ Lifetime ++++++++++++++ "<<endl;
    float Lxy, Lxy3D;
    float LxyPD, LxyANBS, LxyAT, LxyCNBS, LxyCO;
    float BMASS;
    if(process == "Bs") BMASS = BSMASS;
    if(process == "Bd") BMASS = B0MASS;

    TVector3 BsMomentum = Bs.Vect(); //3Momentum of B candidate
    TVector3 PVSVdist = SVpos - PVpos;
      
    Lxy3D = BMASS/Bs.P()*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag();//(PVpos-SVpos).Mag();
    
    BsMomentum.SetZ(0.);
    PVSVdist.SetZ(0.);
  
    Lxy = BMASS*PVSVdist.Dot(BsMomentum)/BsMomentum.Mag2(); //Transverse ctau after setting to 0 the Z components of the vectors above
    LxyPD = BMASS*PVSVdistPD.Dot(BsMomentum)/BsMomentum.Mag2();      
    if (ipv_refitANBS>=0) { LxyANBS=BMASS/Bs.Pt()*PVSVdistANBS.Dot(BsMomentum)/BsMomentum.Mag(); } else { LxyANBS=-99999;}
    if (ipv_refitAT>=0) {LxyAT=BMASS/Bs.Pt()*PVSVdistAT.Dot(BsMomentum)/BsMomentum.Mag(); } else {LxyAT=-99999;}
    if (ipv_refitCNBS>=0) {LxyCNBS=BMASS/Bs.Pt()*PVSVdistCNBS.Dot(BsMomentum)/BsMomentum.Mag(); } else {LxyCNBS=-99999;}      
    if (ipv_refitCO>=0) {LxyCO=BMASS/Bs.Pt()*PVSVdistCO.Dot(BsMomentum)/BsMomentum.Mag(); } else {LxyCO=-99999;}      

    if (debug) cout << "Lxy " << Lxy << " Lxy 3D " << Lxy3D << endl;
    if (debug) cout << "LxyPD " << LxyPD << " LxyANBS " << LxyANBS << endl;
    if (debug) cout << "LxyAT " << LxyAT << " LxyCNBS " << LxyCNBS << endl;
    if (debug) cout << "LxyCO " << LxyCO << endl;
      
    TMatrixF covSV(3,3);
    float covSVArray[]={
      svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
      svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV),
      svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)
    };
    covSV.SetMatrixArray(covSVArray);
      
    TMatrixF covPV(3,3);
    float covPVArray[]={
      pvtSxx->at(pvIndex),pvtSxy->at(pvIndex),pvtSxz->at(pvIndex),
      pvtSxy->at(pvIndex),pvtSyy->at(pvIndex),pvtSyz->at(pvIndex),
      pvtSxz->at(pvIndex),pvtSyz->at(pvIndex),pvtSzz->at(pvIndex)
    };
    covPV.SetMatrixArray(covPVArray);
      
    TMatrixF covTot= covSV+covPV;
      
    float distArray[]={float((SVpos-PVpos).X()),float((SVpos-PVpos).Y()),float((SVpos-PVpos).Z())};
    TVectorF diff(3,distArray);

    float distArray2D[]={float((SVpos-PVpos).X()),float((SVpos-PVpos).Y()),0.};
    TVectorF diff2D(3,distArray2D);
      
    float LxyErr=-999,Lxy3DErr=-999;
      
    if (diff2D.Norm2Sqr()==0 || diff.Norm2Sqr()==0){
      cout << "secondary vertex is exactly the same as PV" << endl;
      return false; //if the secondary vertex is exactly the same as PV continue
    }
      
    TVector3 diff2DV(diff2D[0],diff2D[1],diff2D[2]); 
    LxyErr=TMath::Abs(BMASS/BsMomentum.Mag()*sqrt(covTot.Similarity(diff2D))/sqrt(diff2D.Norm2Sqr())*cos(diff2DV.Angle(BsMomentum))); 
    Lxy3DErr= BMASS/Bs.P()*sqrt(covTot.Similarity(diff))/sqrt(diff.Norm2Sqr()); 
      
    float temparray[]={float(BsMomentum[0]),float(BsMomentum[1]),float(BsMomentum[2])}; 
    TVectorF BsMomentum2DF(3,temparray);
    
    if (debug) cout << "LxyErr " << LxyErr << " Lxy 3DErr " << Lxy3DErr << endl;
        
    //Bc->Bs pi
    if(process == "Bs"){
      int PionIndex=-999;

      float dRtmp=1;
      float pttmp=1;
      for (int i=0; i<nTracks; ++i){
        TLorentzVector pi;
        pi.SetPtEtaPhiM(trkPt->at(i),trkEta->at(i),trkPhi->at(i),PIMASS);
        if (pi.Pt() < 1) continue;
        if (pi.DeltaR(mu_p) < 1.e-3 || pi.DeltaR(mu_m) < 1.e-3 || pi.DeltaR(k_p) < 1.e-3 || pi.DeltaR(k_m) < 1.e-3) continue;
        if ((Bs+pi).M() < 5.0 || (Bs+pi).M() > 6.6) continue;
        float dR=Bs.DeltaR(pi);
        if (dR < dRtmp && pi.Pt()> pttmp){
          PionIndex=i;
          pttmp=pi.Pt();
          dRtmp=dR;
        }
      }
      if (PionIndex>=0){
        TLorentzVector pi;
        pi.SetPtEtaPhiM(trkPt->at(PionIndex),trkEta->at(PionIndex),trkPhi->at(PionIndex),PIMASS);
        (tWriter->TPiPt) = trkPt->at(PionIndex);
        (tWriter->TPiPhi) = trkPhi->at(PionIndex);
        (tWriter->TPiEta) = trkEta->at(PionIndex);
        (tWriter->TPiHits) = trkHitPattern->at(PionIndex);
        (tWriter->TBcMass) = (Bs+pi).M();
      }
    }

    //----- COMPUTING ANGLES
    if(debug) cout<<"  ++++++++++++++ Computing Angles ++++++++++++++ "<<endl;
    // the betas for the boost
    TVector3 p3_JPsi;
    p3_JPsi = Jpsi.Vect();
    p3_JPsi *= -1./Jpsi.E();

    // the boost matrix
    TLorentzRotation boost_jpsi(p3_JPsi);
    TLorentzVector p_JPsi_JPsi = boost_jpsi.VectorMultiplication(Jpsi);
      
    // the different momenta in the new frame
    TLorentzVector p_JPsi_muplus = boost_jpsi.VectorMultiplication(mu_p); 
    TLorentzVector p_JPsi_Kplus = boost_jpsi.VectorMultiplication(k_p); 
    TLorentzVector p_JPsi_phi = boost_jpsi.VectorMultiplication(Phi);
      
    // the 3-momenta
    TVector3 p3_JPsi_muplus = p_JPsi_muplus.Vect();
    TVector3 p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
    TVector3 p3_JPsi_phi = p_JPsi_phi.Vect();
  
    // coordinate system
    TVector3 x,y,z;
    x = p3_JPsi_phi.Unit();
    y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
    y = y.Unit();
    z = x.Cross(y);
      
    // Transversity Basis
    angle_costheta = p3_JPsi_muplus.Unit() * z;

    double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costheta*angle_costheta);
    double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costheta*angle_costheta);
    angle_phi = TMath::ACos(cos_phi);
    if (sin_phi < 0){
      angle_phi =  -angle_phi;
    }
  
    // boosting in phi restframe
    // the betas for the boost
    TVector3 p3_phi;
    p3_phi = Phi.Vect();
    p3_phi *= -1./Phi.E();
      
    // the boost matrix
    TLorentzRotation boost_phi(p3_phi);
    TLorentzVector p_phi_phi;
    p_phi_phi = boost_phi.VectorMultiplication(Phi);
      
    // the different momenta in the new frame
    TLorentzVector p_phi_Kplus = boost_phi.VectorMultiplication(k_p);
    TLorentzVector p_phi_JPsi = boost_phi.VectorMultiplication(Jpsi);
    TLorentzVector p_phi_Bs  = boost_phi.VectorMultiplication(Bs);
      
    // the 3-momenta
    TVector3 p3_phi_Kplus = p_phi_Kplus.Vect();
    TVector3 p3_phi_JPsi = p_phi_JPsi.Vect();
    TVector3 p3_phi_Bs = p_phi_Bs.Vect();
      
    angle_cospsi = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
    
    hCosPsi->Fill(angle_cospsi); 
    hCosTheta->Fill(angle_costheta);
    hAnglePhi->Fill(angle_phi);
      
    if (debug) cout << "CosPsi " << angle_cospsi << " CosTheta " << angle_costheta << " Angle Phi " << angle_phi << endl;

    // Tree Fill

    if(debug) cout<<"  ++++++++++++++ Filling tree ++++++++++++++ "<<endl;

    (tWriter->TmuonsPassSoft) = muonSoft;

    (tWriter->Tangle_cospsi) = angle_cospsi;
    (tWriter->Tangle_costheta) = angle_costheta;
    (tWriter->Tangle_phi) = angle_phi;
      
    (tWriter->TLxy) = Lxy;
    (tWriter->TLxyz) = Lxy3D;
    (tWriter->TLxyErr) = LxyErr;
    (tWriter->TLxyzErr) = Lxy3DErr;
      
    (tWriter->TVprob) = VProb;
    (tWriter->TVprobJpsi) = VProbJpsi;
    (tWriter->TVprobPhi) = VProbPhi;
      
    (tWriter->TBsMass) = Bs.M();
    (tWriter->TBsPt) = Bs.Pt();
    (tWriter->TBsEta) = Bs.Eta();
    (tWriter->TBsPhi) = Bs.Phi();
      
    (tWriter->TBsMassNoRefit) = BsNoRefit.M();
    (tWriter->TJpsiMassNoRefit) = JpsiNoRefit.M();
    (tWriter->TPhiMassNoRefit) = PhiNoRefit.M();
      
    (tWriter->TJpsiMass) = Jpsi.M();
    (tWriter->TJpsiPt) = Jpsi.Pt();
    (tWriter->TJpsiEta) = Jpsi.Eta();
    (tWriter->TJpsiPhi) = Jpsi.Phi();
      
    (tWriter->TPhiMass) = Phi.M();
    (tWriter->TPhiPt) = Phi.Pt();
    (tWriter->TPhiEta) = Phi.Eta();
    (tWriter->TPhiPhi) = Phi.Phi();
      
    (tWriter->TMupPt) = mu_p.Pt();
    (tWriter->TMupPhi) = mu_p.Phi();
    (tWriter->TMupEta) = mu_p.Eta();
    (tWriter->TMupHits) = mupH;
    (tWriter->TMupHltPt) = hltMuPt_p;
      
    (tWriter->TMumPt) = mu_m.Pt();
    (tWriter->TMumPhi) = mu_m.Phi();
    (tWriter->TMumEta) = mu_m.Eta();
    (tWriter->TMumHits) = mumH;
    (tWriter->TMumHltPt) = hltMuPt_m;
      
    (tWriter->TKpPt) = k_p.Pt();
    (tWriter->TKpPhi) = k_p.Phi();
    (tWriter->TKpEta) = k_p.Eta();
    (tWriter->TKpHits) = kpH;
      
    (tWriter->TKmPt) = k_m.Pt();
    (tWriter->TKmPhi) = k_m.Phi();
    (tWriter->TKmEta) = k_m.Eta();
    (tWriter->TKmHits) = kmH;

    mu_p_Bus = mu_p;
    mu_m_Bus = mu_m;
    k_p_Bus = k_p;
    k_m_Bus = k_m;
  }

  //----- MC TRUTH HERE
  if (use_gen || genOnly){
    if(genOnly && process == "Bd"){
      cout<< "GEN ONLY NOT WORKING ON Bd" << endl;
    }else{
      if(debug) cout<<"  ++++++++++++++ Begin GEN INFO "<<nGenP<<" ++++++++++++++ "<<endl;
     
      int BsIndex = -999, JpsiIndex = -999, PhiIndex = -999;
      int KpIndex = -999, KmIndex = -999, MupIndex = -999, MumIndex = -999;
      int genFlavour = 0;

      int processLund[3]; // B lund, jpsi lund, phi/kstar lung
      if(process == "Bs"){
        processLund[0] = 531;
        processLund[1] = 443;
        processLund[2] = 333;
      }
      if(process == "Bd"){
        processLund[0] = 511;
        processLund[1] = 443;
        processLund[2] = 313;
      }

      for (int p=0; p<nGenP; ++p){
        BsIndex = -999;
        JpsiIndex = -999;
        PhiIndex = -999;
        KpIndex = -999;
        KmIndex = -999;
        MupIndex = -999;
        MumIndex = -999;

        if(abs(genId->at(p))!=processLund[0]) continue;
        int tagMix = GetMixStatus( p );
        if(tagMix == 2) continue; // 2 -> Bs that are going to oscillate, 1 -> Bs have oscillated, 0 -> no oscillation

        if(debug) cout<<" found gen "<<genId->at(p)<<endl;

        const vector<int>& aD = allDaughters(p);
        for(int it:aD){
          if(abs(genId->at(it))==processLund[1]) JpsiIndex = it;
          if(abs(genId->at(it))==processLund[2]) PhiIndex = it;
        }
        if( JpsiIndex==-999 || PhiIndex==-999 ) continue;

        const vector<int>& aDjpsi = allDaughters(JpsiIndex);
        const vector<int>& aDphi = allDaughters(PhiIndex);

        if(debug) cout<<" primary Daughters ok "<<endl;

        //Jpsi muons
        for(int it:aDjpsi){
          if(genId->at(it)==13) MumIndex=it;
          if(genId->at(it)==-13) MupIndex=it;
        }
        if( MumIndex==-999 || MupIndex==-999 ) continue;

        //Phi Kaons
        if(process == "Bs"){ /// kstar reco traks not present (to check)
          for(int it:aDphi){
            if(genId->at(it)==321) KpIndex=it;
            if(genId->at(it)==-321) KmIndex=it;
          }
          if( KpIndex==-999 || KmIndex==-999 ) continue;
        }

        if(debug) cout<<" secondary Daughters ok "<<endl;

        BsIndex=p;
        if(genId->at(p)>0) genFlavour=1;
        else  genFlavour=-1;

        if(tagMix == 1) genFlavour*=-1;

        break;
      }

      (tWriter->TgenFlavour) = genFlavour;

      if(BsIndex==-999) (tWriter->TgenFound) = 0;
      else (tWriter->TgenFound) = 1;

      if(genOnly && BsIndex==-999) return false;

      if (GetMixStatus(BsIndex)==0) {
        selectedPVgen = BsIndex;
      }else {
        const vector <int>& aM = allMothers(BsIndex);
        selectedPVgen = aM[0];
      }

      float pvXgen = genVx->at(selectedPVgen);
      float pvYgen = genVy->at(selectedPVgen);
      float pvZgen = genVz->at(selectedPVgen);

      counter[0]++;

      float bestDist = 999.;
      int bestPV = -1;

      for(int iPV = 0; iPV < nPVertices; ++iPV ) {
        float tmpPx = pvtX->at(iPV);
        float tmpPy = pvtY->at(iPV);
        float tmpPz = pvtZ->at(iPV);

        float dx = pvXgen-tmpPx;
        float dy = pvYgen-tmpPy;
        float dz = pvZgen-tmpPz;

        float dist = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

        if(dist < bestDist){
          bestDist = dist;
          bestPV   = iPV;
        }
      }

      if(bestPV !=-1 && bestPV != selectedPV){
        counter[1]++;
        // float dx = genVx->at(selectedPVgen)-pvtX->at(bestPV);
        // float dy = genVy->at(selectedPVgen)-pvtY->at(bestPV);
        // float dz = genVz->at(selectedPVgen)-pvtZ->at(bestPV);
        // float distErr = sqrt( (pow(dx,2)*pow(pvtSxx->at(bestPV),2)
        //               +pow(dy,2)*pow(pvtSyy->at(bestPV),2)
        //               +pow(dz,2)*pow(pvtSzz->at(bestPV),2)
        //               + 2*dx*dy*pvtSxy->at(bestPV)
        //               + 2*dx*dz*pvtSxz->at(bestPV)
        //               + 2*dy*dz*pvtSyz->at(bestPV)
        //               )/(pow(dx,2)+pow(dy,2)+pow(dz,2)) );

        (tWriter->pvtBestDist3D) = bestDist;
        (tWriter->pvtSelectedBestDist3D) = sqrt(pow(pvtX->at(bestPV)-pvtX->at(selectedPV),2)+pow(pvtY->at(bestPV)-pvtY->at(selectedPV),2)+pow(pvtY->at(bestPV)-pvtY->at(selectedPV),2));

        if(useBestPV){
          // Recomputing ct
          TVector3 newPV;
          newPV.SetXYZ(pvtX->at(bestPV),pvtY->at(bestPV),pvtZ->at(bestPV));
          TVector3 PVSVdist = vSV - newPV;
          pSV.SetZ(0.);
          PVSVdist.SetZ(0.);
          float newLxy = BSMASS*PVSVdist.Dot(pSV)/pSV.Mag2();
          (tWriter->TLxy) = newLxy;

          // Recomputing ctErr
          TMatrixF covSV(3,3);
          float covSVArray[]={
            svtSxx->at(selectedSV),svtSxy->at(selectedSV),svtSxz->at(selectedSV),
            svtSxy->at(selectedSV),svtSyy->at(selectedSV),svtSyz->at(selectedSV),
            svtSxz->at(selectedSV),svtSyz->at(selectedSV),svtSzz->at(selectedSV)
          };
          covSV.SetMatrixArray(covSVArray);
            
          TMatrixF covPV(3,3);
          float covPVArray[]={
            pvtSxx->at(bestPV),pvtSxy->at(bestPV),pvtSxz->at(bestPV),
            pvtSxy->at(bestPV),pvtSyy->at(bestPV),pvtSyz->at(bestPV),
            pvtSxz->at(bestPV),pvtSyz->at(bestPV),pvtSzz->at(bestPV)
          };
          covPV.SetMatrixArray(covPVArray);
            
          TMatrixF covTot= covSV+covPV;

          float distArray2D[]={float((vSV-newPV).X()),float((vSV-newPV).Y()),0.};
          TVectorF diff2D(3,distArray2D);
            
          float newLxyErr=-999;
            
          if(diff2D.Norm2Sqr()==0){
            cout << "secondary vertex is exactly the same as PV" << endl;
            return false; //if the secondary vertex is exactly the same as PV continue
          }
            
          TVector3 diff2DV(diff2D[0],diff2D[1],diff2D[2]); 
          newLxyErr = TMath::Abs(BSMASS/pSV.Mag()*sqrt(covTot.Similarity(diff2D))/sqrt(diff2D.Norm2Sqr())*cos(diff2DV.Angle(pSV))); 
          
          (tWriter->TLxyErr) = newLxyErr;

          // recomputing tagging
          setVtxOsMuonTag(selectedSV, bestPV); //set PV and SVT for tagging class
          bool istagged = makeOsMuonTagging();

          (tWriter->Ttag) = getOsMuonTag(); //get tag decision -1, +1 or 0 

          if(istagged){
            (tWriter->Tmtag) = getOsMuonTagMistagProbRaw();
            (tWriter->TmtagCal) = getOsMuonTagMistagProbCalProcess();
            (tWriter->TmtagCalBs) = getOsMuonTagMistagProbCalProcessBuBs();
          }else{
            (tWriter->Tmtag) = 0.5;
            (tWriter->TmtagCal) = 0.5;
            (tWriter->TmtagCalBs) = 0.5;
          }

          // cout<<"Dist: "<<bestDist<<" v "<<sqrt(pow(genVx->at(selectedPVgen)-pvtX->at(selectedPV),2)+pow(genVy->at(selectedPVgen)-pvtY->at(selectedPV),2)+pow(genVz->at(selectedPVgen)-pvtZ->at(selectedPV),2))<<endl;
          // TVector3 vPV;
          // vPV.SetXYZ(pvtX->at(bestPV),pvtY->at(bestPV),pvtZ->at(bestPV));
          // TVector3 diff= vSV-vPV;
          // float cpoint = cos(diff.Angle(pSV));
          // cout<<"Point: "<<cpoint<<" v "<<selectedCosPoint<<endl<<endl;
        }
      }else{
        (tWriter->pvtBestDist3D) = -1;
      }

      if (debug) cout << "Found MC event, tot particles " << nGenP << " indexes:" << BsIndex << " " <<  JpsiIndex << " " <<  PhiIndex << " " << MupIndex  << " " << MumIndex << " " << KpIndex << " "<< KmIndex << endl;

      TVector3 GenSV, GenPV, GenDiff, BsMomentumGen;
      GenSV.SetXYZ(genVx->at(JpsiIndex),genVy->at(JpsiIndex),genVz->at(JpsiIndex));//checked that the Jpsi SV is exactly the same as the Phi SV
      if (GetMixStatus(BsIndex)==0) {
        GenPV.SetXYZ(genVx->at(BsIndex),genVy->at(BsIndex),genVz->at(BsIndex));
      } else {
        const vector <int>& aM = allMothers(BsIndex);
        if( aM.size()>1)  cout<<"Bs with more than one mother"<<endl;
        GenPV.SetXYZ(genVx->at(aM[0]),genVy->at(aM[0]),genVz->at(aM[0]));
      }

      TLorentzVector BsGen, JpsiGen, PhiGen, MupGen, MumGen, KpGen, KmGen;
      
      BsGen.SetPtEtaPhiM(genPt->at(BsIndex), genEta->at(BsIndex), genPhi->at(BsIndex), genMass->at(BsIndex));

      JpsiGen.SetPtEtaPhiM(genPt->at(JpsiIndex), genEta->at(JpsiIndex), genPhi->at(JpsiIndex), genMass->at(JpsiIndex));
      PhiGen.SetPtEtaPhiM(genPt->at(PhiIndex), genEta->at(PhiIndex), genPhi->at(PhiIndex), genMass->at(PhiIndex));
    
      MupGen.SetPtEtaPhiM(genPt->at(MupIndex), genEta->at(MupIndex), genPhi->at(MupIndex), genMass->at(MupIndex));
      MumGen.SetPtEtaPhiM(genPt->at(MumIndex), genEta->at(MumIndex), genPhi->at(MumIndex), genMass->at(MumIndex));
    
      if(process == "Bs"){
        KpGen.SetPtEtaPhiM(genPt->at(KpIndex), genEta->at(KpIndex), genPhi->at(KpIndex), genMass->at(KpIndex));
        KmGen.SetPtEtaPhiM(genPt->at(KmIndex), genEta->at(KmIndex), genPhi->at(KmIndex), genMass->at(KmIndex));
      }
      if(process == "Bd" && !genOnly){
        KpGen.SetPtEtaPhiM(k_p_Bus.Pt(), k_p_Bus.Eta(), k_p_Bus.Phi(), k_p_Bus.M());//GEN tracks are set to the reco ones since they are not included in the collection
        KmGen.SetPtEtaPhiM(k_m_Bus.Pt(), k_m_Bus.Eta(), k_m_Bus.Phi(), k_m_Bus.M());
      }      

      if (debug) cout << " Gen Masses " << BsGen.Mag() << " " << JpsiGen.Mag() << " " << PhiGen.Mag() << " " << MupGen.Mag() << " " << MumGen.Mag() << " " << KpGen.Mag() << " " << KmGen.Mag() << endl;
      
      //mc matching
      int isMCMatched = 0;
      if(!genOnly){
        if(process == "Bs")
           if(KpGen.DeltaR(k_p_Bus) < 0.005 && KmGen.DeltaR(k_m_Bus)< 0.005 && MumGen.DeltaR(mu_m_Bus) < 0.005 &&  MupGen.DeltaR(mu_p_Bus) < 0.005)
            isMCMatched =  1;

        if(process == "Bd")
           if(PhiGen.DeltaR(k_p_Bus+k_m_Bus) < 0.005 && MumGen.DeltaR(mu_m_Bus) < 0.005 &&  MupGen.DeltaR(mu_p_Bus) < 0.005)
            isMCMatched =  1;

        (tWriter->Tmatched) = isMCMatched;
        if (debug) cout << "Distance GEN reco: Kp " << KpGen.DeltaR(k_p_Bus) << " Km " <<  KmGen.DeltaR(k_m_Bus) << " Mum " << MumGen.DeltaR(mu_m_Bus) << " Mup " << MupGen.DeltaR(mu_p_Bus) << endl;
        if (debug) cout << "Distance GEN reco: Kstar " << PhiGen.DeltaR(k_p_Bus+k_m_Bus)  << " Mum " << MumGen.DeltaR(mu_m_Bus) << " Mup " << MupGen.DeltaR(mu_p_Bus) << endl;
      }
        
      // Gen Lxy
      BsMomentumGen = BsGen.Vect();
    
      GenDiff = GenSV-GenPV;
    
      float GenLxy3D = BsGen.M()/BsGen.P()*GenDiff.Dot(BsMomentumGen)/BsMomentumGen.Mag();
    
      GenDiff.SetZ(0.);//Lxy in 2D now
      BsMomentumGen.SetZ(0.);
    
      float GenLxy = BsGen.M()/BsGen.Pt()*GenDiff.Dot(BsMomentumGen)/BsMomentumGen.Mag();
      
      if (debug) cout << "Gen Lxy3D " << GenLxy3D << " Gen Lxy " << GenLxy << " BsIndex " << BsIndex << " JpsiIndex " << JpsiIndex << endl;

      //Computing Gen angles
      double angle_costhetaGen;
      double angle_phiGen;
      double angle_cospsiGen;
    
      TVector3 p3_JPsiGen;
      p3_JPsiGen = JpsiGen.Vect();
      p3_JPsiGen *= -1./JpsiGen.E();
    
      // the boost matrix
      TLorentzRotation boost_jpsiGen(p3_JPsiGen);
      TLorentzVector p_JPsi_JPsiGen = boost_jpsiGen.VectorMultiplication(JpsiGen);
    
      // the different momenta in the new frame                                                                                                       
      TLorentzVector p_JPsi_muplusGen = boost_jpsiGen.VectorMultiplication(MupGen); 
      TLorentzVector p_JPsi_KplusGen = boost_jpsiGen.VectorMultiplication(KpGen); 
      TLorentzVector p_JPsi_phiGen = boost_jpsiGen.VectorMultiplication(PhiGen);      
    
      // the 3-momenta
      TVector3 p3_JPsi_muplusGen = p_JPsi_muplusGen.Vect();
      TVector3 p3_JPsi_KplusGen = p_JPsi_KplusGen.Vect();
      TVector3 p3_JPsi_phiGen = p_JPsi_phiGen.Vect();
    
      // coordinate syste`m
      TVector3 xg,yg,zg;
      xg = p3_JPsi_phiGen.Unit();
      yg = p3_JPsi_KplusGen.Unit() - (xg * (xg * p3_JPsi_KplusGen.Unit()));
      yg = yg.Unit();
      zg = xg.Cross(yg);
    
      // Transversity Basis
      angle_costhetaGen = p3_JPsi_muplusGen.Unit() * zg;
    
      double cos_phiGen = p3_JPsi_muplusGen.Unit() * xg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
      double sin_phiGen = p3_JPsi_muplusGen.Unit() * yg / TMath::Sqrt(1 - angle_costhetaGen*angle_costhetaGen);
      angle_phiGen = TMath::ACos(cos_phiGen);
      if (sin_phiGen < 0){
        angle_phiGen = - angle_phiGen;
      }
    
      // boosting in phi restframe                                                                                                          
      // the betas for the boost
      TVector3 p3_phiGen;
      p3_phiGen = PhiGen.Vect();
      p3_phiGen *= -1./PhiGen.E();
    
      // the boost matrix
      TLorentzRotation boost_phiGen(p3_phiGen);
      TLorentzVector p_phi_phiGen;
      p_phi_phiGen = boost_phiGen.VectorMultiplication(PhiGen);
    
      // the different momenta in the new frame
      TLorentzVector p_phi_KplusGen = boost_phiGen.VectorMultiplication(KpGen);
      TLorentzVector p_phi_JPsiGen = boost_phiGen.VectorMultiplication(JpsiGen);
      TLorentzVector p_phi_BsGen  = boost_phiGen.VectorMultiplication(BsGen);
    
      // the 3-momenta
      TVector3 p3_phi_KplusGen = p_phi_KplusGen.Vect();
      TVector3 p3_phi_JPsiGen = p_phi_JPsiGen.Vect();
      TVector3 p3_phi_BsGen = p_phi_BsGen.Vect();
    
      angle_cospsiGen = -1 * p3_phi_KplusGen.Unit() * p3_phi_JPsiGen.Unit();
    
      if (debug) cout << "CosPsi GEN " << angle_cospsiGen << " CosTheta GEN " << angle_costhetaGen << " Angle Phi GEN " << angle_phiGen << endl;

      // GEN Fill
      (tWriter->Tangle_cospsiGen) = angle_cospsiGen;
      (tWriter->Tangle_costhetaGen) = angle_costhetaGen;
      (tWriter->Tangle_phiGen) = angle_phiGen;
    
      (tWriter->TLxyGen) = GenLxy;
      (tWriter->TLxyzGen) = GenLxy3D;
    
      (tWriter->TBsMassGen) = BsGen.M();
      (tWriter->TBsPtGen) = BsGen.Pt();
      (tWriter->TBsEtaGen) = BsGen.Eta();
      (tWriter->TBsPhiGen) = BsGen.Phi();
    
      (tWriter->TJpsiMassGen) = JpsiGen.M();
      (tWriter->TJpsiPtGen) = JpsiGen.Pt();
      (tWriter->TJpsiEtaGen) = JpsiGen.Eta();
      (tWriter->TJpsiPhiGen) = JpsiGen.Phi();
    
      (tWriter->TPhiMassGen) = PhiGen.M();
      (tWriter->TPhiPtGen) = PhiGen.Pt();
      (tWriter->TPhiEtaGen) = PhiGen.Eta();
      (tWriter->TPhiPhiGen) = PhiGen.Phi();
    
      (tWriter->TMupPtGen) = MupGen.Pt();
      (tWriter->TMupPhiGen) = MupGen.Phi();
      (tWriter->TMupEtaGen) = MupGen.Eta();
      
      (tWriter->TMumPtGen) = MumGen.Pt();
      (tWriter->TMumPhiGen) = MumGen.Phi();
      (tWriter->TMumEtaGen) = MumGen.Eta();
      
      (tWriter->TKpPtGen) = KpGen.Pt();
      (tWriter->TKpPhiGen) = KpGen.Phi();
      (tWriter->TKpEtaGen) = KpGen.Eta();
      
      (tWriter->TKmPtGen) = KmGen.Pt();
      (tWriter->TKmPhiGen) = KmGen.Phi();
      (tWriter->TKmEtaGen) = KmGen.Eta();
      float dx = genVx->at(selectedPVgen)-pvtX->at(selectedPV);
      float dy = genVy->at(selectedPVgen)-pvtY->at(selectedPV);
      float dz = genVz->at(selectedPVgen)-pvtZ->at(selectedPV);
      (tWriter->pvtGenDist) = sqrt(pow(dx,2)+pow(dy,2));
      (tWriter->pvtGenDist3D) = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

    }
  }

  tWriter->fill();

  if(debug) cout<<"  ++++++++++++++ EVENT END ++++++++++++++ "<<endl;

  return true;

}



int PDAnalyzer::HLTMatch(int trackIdx,int isMu){
  
  int hltIdx= -999;

  float mass;//to be set properly for the two objects
  if (isMu==1) mass=0.105658;
  else mass=0.493677;;

  TLorentzVector trk, hltObj;
  trk.SetPtEtaPhiM(trkPt->at(trackIdx),trkEta->at(trackIdx),trkPhi->at(trackIdx),mass);
  
  float minDist=0.1; 

  // type 5 = track, 2 = muon
  for (int k=0; k< nHLTObjects; k++){

    if (hltObjType->at(k) != 2 && isMu==1) continue; // check only matches between trigger muons and muons
    if (hltObjType->at(k) != 5 && isMu==0) continue; // check only matches between trigger tracks and kaons

    hltObj.SetPtEtaPhiM(hltPt->at(k),hltEta->at(k),hltPhi->at(k),mass);
    
    float dRtmp=trk.DeltaR(hltObj);
    
    if (dRtmp < minDist){
      minDist=dRtmp;
      hltIdx=k;
    }
  }
  return hltIdx;
}


int PDAnalyzer::GetMixStatus( unsigned int genindex )
{
  int Code = genId->at( genindex );

  const vector <int>& aD = allDaughters(genindex);
  if( aD.size()>0 && genId->at(aD[0]) == -Code ) return 2;

  const vector <int>& aM = allMothers(genindex); 
  if( aM.size()>0 && genId->at(aM[0]) == -Code ) return 1;

  return 0;
}



void PDAnalyzer::FindPV(TVector3 & sv, TVector3 & pv, TLorentzVector & BsP4, int & index, float &cosPoint) {
  
  TVector3 BsP3 = BsP4.Vect();
  TVector3 tmpPV;
  float tmpCos = -999.;

  int iPV;
  for ( iPV = 0; iPV < nPVertices; ++iPV ) {
    tmpPV.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    TVector3 diff= sv-tmpPV; //vector pointing from PV to SV
    // if (abs(sv.Z()-tmpPV.Z())>0.5) continue; //if the PV has a z distance more than 0.5 cm just skip the vertex
    //diff.SetZ(0.);
    //BsP3.SetZ(0.);
    float cosP=cos(diff.Angle(BsP3));
    if (cosP > tmpCos){
      tmpCos=cosP;
      cosPoint=cosP;
      pv.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
      index=iPV;
    }
  }  

  return;
}

void PDAnalyzer::endJob() {
  tWriter->close();                              // second ntuple
  for(int i=0;i<5;++i) cout<<counter[i]<<endl;
  return;
}


void PDAnalyzer::save() {
#  if UTIL_USE == FULL
  // explicit saving not necessary for "autoSavedObjects"
  autoSave();
#elif UTIL_USE == BARE
  // explicit save histos when not using the full utility
#endif

  return;
}