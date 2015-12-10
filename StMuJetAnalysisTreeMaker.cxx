/*
  Author: Sevil Salur 2007
  Analysis code for the tree reconstruction for Jets...
  Added Mark's request April 2008
  Updated with Matuesz's "user friendly" classes See (TStarJetPico*) April 2008
*/

//  Include header files. 
#include "TFile.h"
#include "TTree.h"

#include "StMessMgr.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <TVector2.h>
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "TVector3.h"

//StEmc
#include "StEmcClusterCollection.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"

#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
#include "StEmcRawMaker/defines.h"

#include "tables/St_emcStatus_Table.h"
#include "tables/St_smdStatus_Table.h"


#include "StMuDSTMaker/COMMON/StMuEmcCollection.h" 
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StMuDSTMaker/COMMON/StMuEmcPoint.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/filters/StEmcFilter.h" 
//
#include "StEmcRawHit.h" 
#include "StEmcModule.h"  
#include "StEmcDetector.h" 
#include "StEmcClusterCollection.h"  
#include "StDaqLib/EMC/StEmcDecoder.h"

#include "tables/St_emcStatus_Table.h"
#include "tables/St_smdStatus_Table.h"

#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDebug.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"

#include "StEvent/StBTofHeader.h"

#include "StEvent/StTriggerId.h"
#include "StTriggerUtilities/StTriggerSimuMaker.h"
#include "StTriggerUtilities/StTriggerSimuResult.h"
#include "StTriggerUtilities/Bemc/StBemcTriggerSimu.h"
#include "StTriggerUtilities/L2Emulator/StL2TriggerSimu.h"

//#include "StStrangeMuDstMaker/StV0MuDst.hh"

#include "StMuJetAnalysisTreeMaker.h"

#include "TStMuCutEventJet.h"
//#include "TStMuCutV0Jet.h"
//#include "eventStructure/TParticleTow.h"

#include "JetPicoEvent/TStarJetPicoEvent.h"
#include "JetPicoEvent/TStarJetPicoEventHeader.h"
#include "JetPicoEvent/TStarJetPicoPrimaryTrack.h"
#include "JetPicoEvent/TStarJetPicoGlobalTrack.h"
//#include "JetPicoEvent/TStarJetPicoTower.h"
//#include "JetPicoEvent/TStarJetPicoV0.h"
#include "JetPicoEvent/TStarJetPicoTriggerInfo.h"

#include "JetPicoEvent/TStarJetPicoDefinitions.h"
//#include "JetPicoEvent/TStarJetPicoQAHistograms.h"
#include "JetPicoEvent/TStarJetPicoUtils.h"

#include "SystemOfUnits.h"
#include "JetPicoEvent/TStarJetPicoPairs.h"

#define TOFILLTREE 0


#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

//const int nadc_max = 4; //number of emc cluster hits

//  Prototype 
void muEventInfo(StMuEvent&, const Int_t&);


ClassImp(StMuJetAnalysisTreeMaker)
  
  // The constructor. Initialize data members here.
  StMuJetAnalysisTreeMaker::StMuJetAnalysisTreeMaker(const Char_t *name) : StMaker(name)
{ 
  mEventCounter = 0;   
  mInputEventCounter = 0;
  mFile = NULL; 
  mFile_monitors = NULL; 
  fNHits=20;
  fEta=0.7;
  //mGeom = StEmcGeom::instance("bemc");
  mFilter   = new StEmcFilter();  
  mTables   = new StBemcTables(); 
  rMatch    = 0; //counter for matched
//  fMatchTrArr = 0;
//  fMatchTrEtaArr = 0;
//  fMatchTrPhiArr = 0;
//  fPrimIndexArray = 0;
//  fGlobIndexArray = 0;
//  fEtaArray = 0;
//  fPhiArray = 0;

/*  idMuDstPrimIndexArray = 0;
  idMuDstGlobIndexArray = 0;
  idPrimIndexArray = 0;
  idGlobIndexArray = 0;
  PrimIndex = 0;
  GlobIndex = 0;*/

  //fRefMultCorr=new StRefMultCorr();

  //mDoV0 = kFALSE;
  mDoPrimTrks = kTRUE;
  mDoGlobTrks = kTRUE;
  //fTCutV0       = new TStMuCutV0Jet();
  //fTCutV0->SetVerbose(0); 
  //fNTOTMatchedTr=0;
  //fMatchedTow=0;
  
  //moved here from Init: enable setting Verbose in macro
  fTCutEvent       = new TStMuCutEventJet();
  fTCutEvent->SetVerbose(1); 
}

StMuJetAnalysisTreeMaker::~StMuJetAnalysisTreeMaker() { 
  //delete fRefMultCorr;
  delete fTCutEvent;
  //delete fTCutV0;
  delete mFilter;  
  delete mTables;
 }
//
//  Called once at the beginning.
Int_t StMuJetAnalysisTreeMaker::Init()
{

  cout << "open a root ouput file " << endl;
  //mFile = new TFile(fFilename,"RECREATE","PicoFile");
  mFile = new TFile(file_name,"RECREATE");

  mFile->cd();
   
  cout<<"!!!!!!!!!!!!!!!!!!!!!!!Setting Event Cuts!!!!!"<<endl;
    switch (nFlagData){
  case 0 :{ // U+U, year 12 all
    cout<<" U+U, year 12 all"<<endl;
    fTCutEvent->SetStandardCuts_UU_Y12_all();
    fTCutEvent->SetMult(0,5000);
    fTCutEvent->SetVertexZ(-100.,100.); 
    fTCutEvent->SetVzVpdDiff(3);
  }
  case 1: { // U+U, year 12 minb
    cout<<" U+U, year 12 minb+NPE"<<endl;
    fTCutEvent->SetStandardCuts_UU_Y12_minB_NPE();
    fTCutEvent->SetMult(0,5000);
    fTCutEvent->SetVertexZ(-30.,30.);
    fTCutEvent->SetVzVpdDiff(3);
     break;
  }
    case 2:{ // U+U, year 12 central
     cout<<" U+U, year 12 central"<<endl;
  fTCutEvent->SetStandardCuts_UU_Y12_central();
    //fTCutEvent->SetMult(0,5000);
    //fTCutEvent->SetVertexZ(-30.,30.); 
    //fTCutEvent->SetVzVpdDiff(3);
    break;
  }
  default: assert(0);
  }
  
 cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
 


  cout<<"******** doing Init  **********"<<endl; 
  //MEvent =0;
  //mFile->cd();
  //MTree=new TTree("JetTree"," Pico Tree for Jet");
  //MTree->Branch("PicoJetTree","TStarJetPicoEvent",&MEvent);//,64000,50);
  //MTree->SetAutoSave(100000);

//  initEmc();	   ///......................uncomment
//  HistAllocation();    //..................ucomment

//******multiplicity histrogams
//
//

refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();

mult_mb = new TH1F("mult_mb","mult_mb",800,0,800);
   mult_mb->Sumw2();
   centrality = new TH1F("centrality","centrality",9,0,9);
   centrality->Sumw2();
   mult_cent = new TH1F("mult_cent","mult_cent",800,0,800);
   mult_cent->Sumw2();
   mult_cent_cut = new TH1F("mult_cent_cut","mult_cent_cut",800,0,800);
   mult_cent_cut->Sumw2();
   mult_cent_cut_5_10 = new TH1F("mult_cent_cut_5_10","mult_cent_cut_5_10",800,0,800);
   mult_cent_cut_5_10->Sumw2();
   mult_cent_cut_10_20 = new TH1F("mult_cent_cut_10_20","mult_cent_cut_10_20",800,0,800);
   mult_cent_cut_10_20->Sumw2();
   mult_runid = new TH1F("mult_runid","mult_runid",10000, 13114000, 13137010);
   mult_runid->Sumw2();
   event_runid = new TH1F("event_runid","event_runid",10000, 13114000, 13137010);
   event_runid->Sumw2();
   mult_runid_cent = new TH1F("mult_runid_cent","mult_runid_cent",10000, 13114000, 13137010);
   mult_runid_cent->Sumw2();
   event_runid_cent = new TH1F("event_runid_cent","event_runid_cent",10000, 13114000, 13137010);
   event_runid_cent->Sumw2();
   mult_runid_cent_cut = new TH1F("mult_runid_cent_cut","mult_runid_cent_cut",10000, 13114000, 13137010);
   mult_runid_cent_cut->Sumw2();
   event_runid_cent_cut = new TH1F("event_runid_cent_cut","event_runid_cent_cut",10000, 13114000, 13137010);
   event_runid_cent_cut->Sumw2();

   mult_runid_zoom = new TH1F("mult_runid_zoom","mult_runid_zoom",70, 13117000, 13117070);
   mult_runid_zoom->Sumw2();
   mult_runid_zoom_corr = new TH1F("mult_runid_zoom_corr","mult_runid_zoom_corr",70, 13117000, 13117070);
   mult_runid_zoom_corr->Sumw2();
   event_runid_zoom = new TH1F("event_runid_zoom","event_runid_zoom",70, 13117000, 13117070);
   event_runid_zoom->Sumw2();
   mult_runid_cent_zoom = new TH1F("mult_runid_cent_zoom","mult_runid_cent_zoom",70, 13117000, 13117070);
   mult_runid_cent_zoom->Sumw2();
   mult_runid_cent_zoom_corr = new TH1F("mult_runid_cent_zoom_corr","mult_runid_cent_zoom_corr",70, 13117000, 13117070);
   mult_runid_cent_zoom_corr->Sumw2();
   event_runid_cent_zoom = new TH1F("event_runid_cent_zoom","event_runid_cent_zoom",70, 13117000, 13117070);
   event_runid_cent_zoom->Sumw2();
   mult_runid_cent_cut_zoom = new TH1F("mult_runid_cent_cut_zoom","mult_runid_cent_cut_zoom",70, 13117000, 13117070);
   mult_runid_cent_cut_zoom->Sumw2();
   mult_runid_cent_cut_zoom_corr = new TH1F("mult_runid_cent_cut_zoom_corr","mult_runid_cent_cut_zoom_corr",70, 13117000, 13117070);
   mult_runid_cent_cut_zoom_corr->Sumw2();
   event_runid_cent_cut_zoom = new TH1F("event_runid_cent_cut_zoom","event_runid_cent_cut_zoom",70, 13117000, 13117070);
   event_runid_cent_cut_zoom->Sumw2();
   
   multCorr = new TH1F("multCorr","multCorr",800,0,800);
   multCorr->Sumw2();
   multCorr_cent = new TH1F("multCorr_cent","multCorr_cent",800,0,800);
   multCorr_cent->Sumw2();
   multCorr_cent_cut = new TH1F("multCorr_cent_cut","multCorr_cent_cut",800,0,800);
   multCorr_cent_cut->Sumw2();
hWeight = new TH1F("hWeight","hWeight;weight",100,0,2);
//******end of multiplicity histograms

  cout<<"******** end of Init Stage **********"<<endl; 


  return StMaker::Init();
}
//
//  Called every event after Make(). 
void StMuJetAnalysisTreeMaker::Clear(Option_t *opt)
{
 // delete [] fMatchTrArr;fMatchTrArr=0;
 // delete [] fMatchTrEtaArr;fMatchTrEtaArr=0;
 // delete [] fMatchTrPhiArr;fMatchTrPhiArr=0;
//  delete [] fPrimIndexArray; fPrimIndexArray=0;
//  delete [] fGlobIndexArray; fGlobIndexArray=0;
//  delete [] fEtaArray; fEtaArray=0;
 // delete [] fPhiArray; fPhiArray=0;

  /*delete [] idMuDstPrimIndexArray; idMuDstPrimIndexArray=0;
  delete [] idMuDstGlobIndexArray; idMuDstGlobIndexArray=0;
  delete [] idPrimIndexArray; idPrimIndexArray=0;
  delete [] idGlobIndexArray; idGlobIndexArray=0;*/
 // delete PrimIndex;
 // delete GlobIndex;

  StMaker::Clear();
}
//

void   StMuJetAnalysisTreeMaker::setRootFile(const char* outpath,const char* filename)
{fFilename = outpath;
 fFilename +="/picoDst_";
 fFilename +=filename;fFilename +=".root";
 
 file_name = outpath;
 file_name +="/allMuDst_histos_";
 file_name +=filename;
 file_name +=".root";
} 
 
//  Called once at the end.
Int_t StMuJetAnalysisTreeMaker::Finish()
{
  //  Summarize the run.
  cout << "StMuJetAnalysisTreeMaker::Finish()\n";
  cout << "\tProcessed " << mEventCounter << " events." << endl;
  cout << "\tInput events: "<<mInputEventCounter<<endl;
//  mFile->cd();


//  MTree->Write(NULL,TFile::kWriteDelete); 



// mFile->cd();

//-----------
cout<<"saving histograms"<<endl;

mFile->cd();

//*********multiplicity histograms

mult_mb->Write();
   centrality->Write();
   mult_cent->Write();
   mult_cent_cut->Write();
   mult_cent_cut_5_10->Write();
   mult_cent_cut_10_20->Write();
   mult_runid->Write();
   event_runid->Write();
   mult_runid_cent->Write();
   event_runid_cent->Write();
   mult_runid_cent_cut->Write();
   event_runid_cent_cut->Write();
   mult_runid_zoom->Write();
   mult_runid_zoom_corr->Write();
   event_runid_zoom->Write();
   mult_runid_cent_zoom->Write();
   mult_runid_cent_zoom_corr->Write();
   event_runid_cent_zoom->Write();
   mult_runid_cent_cut_zoom->Write();
   mult_runid_cent_cut_zoom_corr->Write();
   event_runid_cent_cut_zoom->Write();
   multCorr->Write();
   multCorr_cent->Write();
   multCorr_cent_cut->Write();
hWeight->Write();

//*******end of multiplicity histograms

/*
Vz_before->Write();
DiffVz_before->Write();

hVz->Write();
hDiffVz->Write();

refMultiplicity->Write();
centrality->Write();
mult_cent->Write();

hdEdx->Write();
hnsigma->Write();

firstTPC->Write();
nHits_nMax_pt->Write();
ndEdx_pt->Write();
eta_pt->Write();
phi_pt->Write();
dca_pt->Write();

electron_pt->Write();
nsigma_purity_pt->Write();

hDeltaZ->Write();
hDeltaZ_unlike->Write();
hDeltaZ_like->Write();
hDeltaPhi->Write();
hDeltaPhi_unlike->Write();
hDeltaPhi_like->Write();
hpE->Write();
hpE_unlike->Write();
hpE_like->Write();
hnSMDE->Write();
smde_unlike->Write();
smde_like->Write();
hnSMDP->Write();
smdp_unlike->Write();
smdp_like->Write();

hinvmass_pt_unlike->Write();
hinvmass_pt_like->Write();
nsigma_partner_unlike->Write();
nsigma_partner_like->Write();

photonic_pt_unlike->Write();
photonic_pt_like->Write();

nsigma_primary_unlike->Write();
nsigma_primary_like->Write();

noEMC_unlike->Write();
noEMC_like->Write();
EMC_unlike->Write();
EMC_like->Write();

eta_pt_unlike->Write();
phi_pt_unlike->Write();
dca_pt_unlike->Write();
eta_pt_like->Write();
phi_pt_like->Write();
dca_pt_like->Write();

nHits_nHitsRatio_pt_partner_unlike->Write();
nHits_nHitsRatio_pt_partner_like->Write();
pT_partner_unlike->Write();
pT_partner_like->Write();
pair_dca_unlike->Write();
pair_dca_like->Write();
*/
 mFile->Close();


//  finishEmc();  //.............uncomment

  cout<<"End of StMuJetAnalysisTreeMaker::Finish"<<endl;

  return kStOK;
}
//
//  This method is called every event.
Int_t StMuJetAnalysisTreeMaker::Make()
{
  //  cout<<"EVENT = "<<mEventCounter<<endl;
  //   mEventCounter++;  // increase counter
cout<<"starting Make"<<endl;

  ++mInputEventCounter;
  DEBUGVALUE2(mEventCounter);
  
  mCurrentMu=  (StMuDst*) GetInputDS("MuDst"); 
  DEBUGVALUE2(mCurrentMu);
  
  if (!mCurrentMu){
    gMessMgr->Warning() << "StMuJetAnalysisTreeMaker::Make : No MuDst" << endm;
    return kStOK;        // if no event, we're done
  }

  
  //  Check StMuEvent branch
  StMuEvent* muEvent =  mCurrentMu->event();
 /*if (mCurrentMu->btofHeader())
   QAHist->Vpd_crosscheck->Fill(muEvent->vpdVz(),mCurrentMu->btofHeader()->vpdVz());
 else QAHist->Vpd_crosscheck->Fill(muEvent->vpdVz(),-500);*/

  if (!muEvent) {
    cout << "ReadStMuEvent() -- no event class???" << endl;
    return kStOK;
  }
  
 
  
  //// Event cuts check triggers - for saving histos only
//  if(fTCutEvent->CheckEventTriggers(muEvent) ==kFALSE) {    //........uncomment
//    if (mVerbose)cout<<" ### event not taken - triggers ### "<<endl;   //...............uncomment
//	return kStOK;    ///KATKA added
//  }     //.............uncomment
  
 
  //  if (mVerbose)cout<<" ### Petr event passed triggers ### "<<endl;

  //Event cuts - includes also the trigger cut
//  if(fTCutEvent->CheckEvent(mCurrentMu) != 1) {   //............uncomment
//    if (mVerbose)cout<<" ### event not taken ### "<<endl;     //................uncomment
//    return kStOK;   //...........................uncomment
//  }        //.................uncomment
  
/*  
  StMuPrimaryVertex* fprimaryVertex =  mCurrentMu->primaryVertex();
  PrimVertexZ=fprimaryVertex->position().z(); //for tower Eta correction due to z shift of the primary vertex

cout<<"before writing primary vertex position"<<endl; 
  if(mVerbose) {
    cout<<"PRIMARY VERTEX = "<<fprimaryVertex->position().x()<<"  "<<fprimaryVertex->position().y()<<"  "<<fprimaryVertex->position().z()<<"  "<<endl;
    cout<<"PV errors: "<<fprimaryVertex->posError().x()<<"  "<<fprimaryVertex->posError().y()<<"  "<<fprimaryVertex->posError().z()<<"  "<<endl;
  } */ 
    //backward compatibility
  if(mVerbose)cout<<"EVENT = "<<mEventCounter<<endl;
  mEventCounter++;  // increase counter
  

 

  //fill jet event
//  MEvent=new TStarJetPicoEvent();
  

  Bfield  = 0.1*muEvent->runInfo().magneticField(); //this is in tesla,ie 0.1*kilogauss
  cout<<"     ***** B Field = "<<Bfield<<endl;
  //StEventSummary &eventsum = muEvent->eventSummary();

 

  /*   //.......................uncomment
  if(!mCurrentMu->muEmcCollection()) {
    cout<<"NO EMC collection - skipping event"<<endl;
    delete MEvent;
    return kStOk;      
  }
 
  mEmcCol = (StEmcCollection*)mCurrentMu->emcCollection();  //KAtka
  assert(mEmcCol);
  //emcCollection = (StEmcCollection*)mCurrentMu->emcCollection();
  //assert(emcCollection);
  //.....................uncomment
    */

/*  if (mDoPrimTrks){
    doPrimTrks();  
    }


  if (mDoGlobTrks) doGlobTrks();  

  doPhotoElectrons();*/

 
 
  
 // rBarrelPts = mEmcCol->barrelPoints().size();
 /* 
  
  mTables->loadTables((StMaker*)this);
    //  cout<<" reaction plane  "<<GetReactionPlane()<<endl;
  rplane = GetReactionPlane();
  cout<<" reaction plane  ="<<rplane<<endl; //PVok
  
   
  //fill header
  TStarJetPicoEventHeader *picoHeader=MEvent->GetHeader();
  assert(picoHeader);
  
  //---eventlevel tree variables
  picoHeader->SetNominalTriggerId(muEvent->triggerIdCollection().nominal());
  picoHeader->SetEventId(static_cast<Int_t>(muEvent->eventId()));
  picoHeader->SetRunId(static_cast<Int_t>(muEvent->runId()));
  picoHeader->SetReferenceMultiplicity(static_cast<Int_t>(muEvent->refMult()));
  //Petr  picoHeader->SetGReferenceMultiplicity(static_cast<Int_t>(muEvent->grefmult()));
  picoHeader->SetNMuDstPrimaryTracks(static_cast<Int_t>(mCurrentMu->numberOfPrimaryTracks()));
  picoHeader->SetNMuDstGlobalTracks(static_cast<Int_t>(mCurrentMu->numberOfGlobalTracks()));
  picoHeader->SetReactionPlaneAngle(static_cast<Float_t>(rplane));
//  picoHeader->SetNOfMatchedTowers(fMatchedTow);

  //picoHeader->SetNOfTowers(MEvent->GetTowers()->GetEntriesFast());
  //picoHeader->SetNOfPrimaryTracks(MEvent->GetPrimaryTracks()->GetEntriesFast());
  //cout<<"HP 1"<<endl;
 
//  picoHeader->SetNOfMatchedTracks(fNTOTMatchedTr);
  
  picoHeader->SetNOfEMCPoints(static_cast<Int_t>(rBarrelPts));
  picoHeader->SetPrimaryVertexX(static_cast<Float_t>((eventsum.primaryVertexPosition()).x()));  
  picoHeader->SetPrimaryVertexY(static_cast<Float_t>((eventsum.primaryVertexPosition()).y()));  
  picoHeader->SetPrimaryVertexZ(static_cast<Float_t>(PrimVertexZ));  
  picoHeader->SetCTBMultiplicity(muEvent->ctbMultiplicity());      
  picoHeader->SetTOFMultiplicity(muEvent->btofTrayMultiplicity());      
  if (fprimaryVertex) {
    picoHeader->SetPrimaryVertexMeanDipAngle(fprimaryVertex->meanDip());
    picoHeader->SetPrimaryVertexRanking(fprimaryVertex->ranking());
  } else {
    picoHeader->SetPrimaryVertexMeanDipAngle(-99);
    picoHeader->SetPrimaryVertexRanking(-99);
  }
  picoHeader->SetNumberOfVertices(mCurrentMu->primaryVertices()->GetEntries());
  picoHeader->SetDSMInput(muEvent->l0Trigger().dsmInput());
  picoHeader->SetTrigMask(muEvent->eventInfo().triggerMask());
 
  //VPD
   if (mCurrentMu->btofHeader())picoHeader->SetVpdZ(mCurrentMu->btofHeader()->vpdVz());
  else picoHeader->SetVpdZ(-500);
 
  //runinfo
  picoHeader->SetZdcWestRate(static_cast<Float_t>(muEvent->runInfo().zdcWestRate()));
  picoHeader->SetZdcEastRate(static_cast<Float_t>(muEvent->runInfo().zdcEastRate()));
  picoHeader->SetZdcCoincidenceRate(static_cast<Float_t>(muEvent->runInfo().zdcCoincidenceRate()));
  picoHeader->SetBbcWestRate(static_cast<Float_t>(muEvent->runInfo().bbcWestRate()));
  picoHeader->SetBbcEastRate(static_cast<Float_t>(muEvent->runInfo().bbcEastRate()));
  picoHeader->SetBbcCoincidenceRate(static_cast<Float_t>(muEvent->runInfo().bbcCoincidenceRate()));
  picoHeader->SetBbcBlueBackgroundRate(static_cast<Float_t>(muEvent->runInfo().bbcBlueBackgroundRate()));
  picoHeader->SetBbcYellowBackgroundRate(static_cast<Float_t>(muEvent->runInfo().bbcYellowBackgroundRate()));
  picoHeader->SetBbcAdcSumEast(static_cast<Int_t>(muEvent->bbcTriggerDetector().adcSumEast()));

 //clear vectors
 idMuDstPrimIndex.clear();
 idMuDstGlobIndex.clear();
 idJetTreePrimIndex.clear();
 idJetTreeGlobIndex.clear();
*/

 // cout<<" #prim="<<mCurrentMu->numberOfPrimaryTracks()<<endl; 
  
 //MTree->Fill();
    
  
  //delete MEvent;
  //cout<<"KG 4"<<endl;


//...........multiplicity histograms
//StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
//double corr = refmultCorrUtil->getRefMultCorr();

double Vz = muEvent->primaryVertexPosition().z();
double vpdVz = muEvent->vpdVz();
double DiffVz = vpdVz-Vz;

double trigCent[4] = {400132, 400142, 400122, 400102};//..physics central5
int ntrigCent = 4;
double trig[6] = {400004, 400005, 400014, 400015, 400025, 400035}; //..MB
int ntrig = 6;
bool isTrigger = false;
bool isTriggerCent = false;

for(int i=0; i<ntrig; i++){
        isTrigger = muEvent->triggerIdCollection().nominal().isTrigger(trig[i]);
        if(isTrigger == true) break;
}

for(int iCent=0; iCent<ntrigCent; iCent++){
        isTriggerCent = muEvent->triggerIdCollection().nominal().isTrigger(trigCent[iCent]);
        if(isTriggerCent == true) break;
}

if(fabs(Vz)<30 && fabs(DiffVz)<3){// && isTrigger == true){

	refmultCorrUtil->init(muEvent->runId());
	refmultCorrUtil->initEvent(muEvent->refMult(), Vz, muEvent->runInfo().zdcCoincidenceRate());
	if (refmultCorrUtil->isBadRun(muEvent->runId()))
      	{
                return kStOK;
        }
	const double corr = refmultCorrUtil->getRefMultCorr();
	const double weight = refmultCorrUtil->getWeight();

        if(isTrigger == true){
                mult_mb->Fill(muEvent->refMult());
                mult_runid->Fill(muEvent->runId(), corr);
                event_runid->Fill(muEvent->runId());
                mult_runid_zoom->Fill(muEvent->runId(), muEvent->refMult());
                event_runid_zoom->Fill(muEvent->runId());
                mult_runid_zoom_corr->Fill(muEvent->runId(), corr);
		multCorr->Fill(corr);
                }
        if(isTriggerCent == true){

                mult_cent->Fill(muEvent->refMult());
                mult_runid_cent->Fill(muEvent->runId(), corr);
                event_runid_cent->Fill(muEvent->runId());
                mult_runid_cent_zoom->Fill(muEvent->runId(), muEvent->refMult());
                event_runid_cent_zoom->Fill(muEvent->runId());
                mult_runid_cent_zoom_corr->Fill(muEvent->runId(), corr);
		multCorr_cent->Fill(corr);

                //refmultCorrUtil->init(muEvent->runId());
                //refmultCorrUtil->initEvent(muEvent->refMult(), Vz, muEvent->runInfo().zdcCoincidenceRate());
                const Int_t cent9 = refmultCorrUtil->getCentralityBin9();
                centrality->Fill(cent9);

                if(cent9 == 7){mult_cent_cut_5_10->Fill(muEvent->refMult());}//..5-10%
                if(cent9 == 6){mult_cent_cut_10_20->Fill(muEvent->refMult());}//..10-20%

                if(cent9 == 8){ //..0-5% centrality
			hWeight->Fill(weight);
                        mult_cent_cut->Fill(muEvent->refMult());
                        mult_runid_cent_cut->Fill(muEvent->runId(), corr);
                        event_runid_cent_cut->Fill(muEvent->runId());
                        mult_runid_cent_cut_zoom->Fill(muEvent->runId(),muEvent->refMult());
                        event_runid_cent_cut_zoom->Fill(muEvent->runId());
                        mult_runid_cent_cut_zoom_corr->Fill(muEvent->runId(), corr);
			multCorr_cent_cut->Fill(corr);
                }
        }

}

//..........end of multiplicity histograms

/*

//BEGINS PART TO ONLY SAVE HISTOGRAMS----------------------------
StRefMultCorr* refmultCorrUtil  = CentralityMaker::instance()->getRefMultCorr() ;

cout<<"begin of part for saving only histograms"<<endl;
Int_t nTracks = mCurrentMu->primaryTracks()->GetEntries();
Int_t nGlobTracks = mCurrentMu->globalTracks()->GetEntries();
const Double_t ELECMASS=0.00051099906;

float Vz = muEvent->primaryVertexPosition().z();
//Float_t vpdVz = mCurrentMu->btofHeader()->vpdVz();
float vpdVz = muEvent->vpdVz();
float DiffVz = vpdVz-Vz;

Vz_before->Fill(Vz);
DiffVz_before->Fill(DiffVz);


 if(Vz<30 && Vz>-30 && DiffVz<3 && DiffVz>-3){//if passes event cuts

   float multiplicity = muEvent->refMult();
   refMultiplicity->Fill(multiplicity);

   refmultCorrUtil->init(muEvent->runId());
   refmultCorrUtil->initEvent(muEvent->refMult(), Vz, muEvent->runInfo().zdcCoincidenceRate());

   const Int_t cent9 = refmultCorrUtil->getCentralityBin9();
   centrality->Fill(cent9);

   if(cent9 == 8){ //..0-5% centrality

      mult_cent->Fill(muEvent->refMult());

      hVz->Fill(Vz);
      hDiffVz->Fill(DiffVz);

      for(int i=0;i<nTracks;i++){//..loop over primary tracks

	StMuTrack* track = mCurrentMu->primaryTracks(i);
  	StMuTrack* gl_track = (StMuTrack*)track->globalTrack();
	
	assert(track);

	float DCA = gl_track->dca().mag();
	double pT = gl_track->pt();
	double nHits = gl_track->nHits();
	double possHits = gl_track->nHitsPoss();
	double hitsRatio = nHits/possHits;
	double dedxHits = gl_track->nHitsDedx();
	double pseudorapidity = gl_track->eta();
	double phi = gl_track->phi();
	StThreeVectorF firstpoint = gl_track->firstPoint();
	double firstPoint_x = firstpoint.x();
	double firstPoint_y = firstpoint.y();
	double firstPoint = TMath::Sqrt(firstPoint_x*firstPoint_x + firstPoint_y*firstPoint_y);

	int id; int adc0; float e[2]; float dist[2]; int nhit[2]; int ntow;
        getBEMC(gl_track, &id, &adc0, e, dist, nhit, &ntow);

        float deltaz = dist[0];
        float deltaphi = dist[1];
        int nSMDE = nhit[0];
        int nSMDP = nhit[1];
        float E0 = e[0];
        double p = gl_track->p().mag();
        double ratio = p/E0;
        double nsigma_primary = gl_track->nSigmaElectron();
        short charge1 = gl_track->charge();

	//..filling hists for single track rec. efficiency systematics (must be inclusive electrons)
	if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5){
	  if(TMath::Abs(DCA)<1.5 && pT>1.2 && nHits>=20 && hitsRatio>0.52 && dedxHits>=15 && TMath::Abs(pseudorapidity)<0.7){
	    firstTPC->Fill(pT, firstPoint);}
	  if(TMath::Abs(DCA)<1.5 && pT>1.2 && nHits>=20 && hitsRatio>0.52 && dedxHits>=15 && firstPoint<73){
	    eta_pt->Fill(pT, pseudorapidity);}	
	  if(TMath::Abs(DCA)<1.5 && pT>1.2 && nHits>=20 && hitsRatio>0.52 && TMath::Abs(pseudorapidity)<0.7 && firstPoint<73){
	    ndEdx_pt->Fill(pT, dedxHits);} 
	  if(TMath::Abs(DCA)<1.5 && pT>1.2 && dedxHits>=15 && TMath::Abs(pseudorapidity)<0.7 && firstPoint<73){
	    nHits_nMax_pt->Fill(pT, nHits, hitsRatio);}
	  if(pT>1.2 && nHits>=20 && hitsRatio>0.52 && dedxHits>=15 && TMath::Abs(pseudorapidity)<0.7 && firstPoint<73){
	    dca_pt->Fill(pT, DCA);}
	}
	
	//..good tracks
	if(TMath::Abs(DCA)<1.5 && pT>1.2 && nHits>=20 && hitsRatio>0.52 && dedxHits>=15 && TMath::Abs(pseudorapidity)<0.7 && firstPoint<73){
		
	  hdEdx->Fill(pT, gl_track->dEdx());
	  hnsigma->Fill(pT, nsigma_primary);

	  //..filling hists for possible investigation of EMC eff.
	  if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && -0.5<nsigma_primary && nsigma_primary<2.5){
	    hpE->Fill(pT, ratio);}
	  if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && 0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.){
	    hnSMDP->Fill(pT, nSMDP);}
	  if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDP>1 && 0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5){
	    hnSMDE->Fill(pT, nSMDE);}
	  if(TMath::Abs(deltaz)<3 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5){
	    hDeltaPhi->Fill(pT, deltaphi);}
	  if(TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5){
	    hDeltaZ->Fill(pT, deltaz);}

	  //..inclusive electrons
 	  if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5){
	    phi_pt->Fill(pT, phi);
	    electron_pt->Fill(pT);	
	  }
		
	  //..nsigma plots for purity
   	  if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2){
	    nsigma_purity_pt->Fill(pT, nsigma_primary);
	  }

	//..loop over global tracks
	for(int j=0;j<nGlobTracks;j++){

	  StMuTrack* gtrack = mCurrentMu->globalTracks(j);
	  if(gl_track->id() == gtrack->id()) continue;

  	  assert(gtrack);
	
	  double pT_global = gtrack->pt();
	  double nHits_global = gtrack->nHits();
          double possHits_global = gtrack->nHitsPoss();
          double hitsRatio_global = nHits_global/possHits_global;	
 	  double nsigma_global = gtrack->nSigmaElectron();
	  short charge2 = gtrack->charge();

	  //..good global tracks
	  if(pT_global>0.2 && nHits_global>=20 && hitsRatio_global>0.52 && TMath::Abs(nsigma_global)<5){
				
	     StPhysicalHelixD e1Helix = (gl_track->helix());
	     StPhysicalHelixD e2Helix = (gtrack->helix());

	     pair<double,double> s = e1Helix.pathLengths(e2Helix);
	     StThreeVectorD e1PosAtDca = e1Helix.at(s.first);
	     StThreeVectorD e2PosAtDca = e2Helix.at(s.second);
	     StThreeVectorD PosAtDca = (e1PosAtDca + e2PosAtDca)/2;
	     double pairDCA = (e1PosAtDca - e2PosAtDca).mag();//pairDCA

	     StThreeVectorD p1DCA=e1Helix.momentumAt(s.first,Bfield*tesla);
	     StThreeVectorD p2DCA=e2Helix.momentumAt(s.second,Bfield*tesla);

             StLorentzVectorD LV1,LV2;
             LV1.setPx(p1DCA.x());
             LV1.setPy(p1DCA.y());
             LV1.setPz(p1DCA.z());
             LV1.setE(p1DCA.massHypothesis(ELECMASS));

             LV2.setPx(p2DCA.x());
             LV2.setPy(p2DCA.y());
             LV2.setPz(p2DCA.z());
             LV2.setE(p2DCA.massHypothesis(ELECMASS));

             double massDCA=(LV1+LV2).m();

	     //..good pairs
	     if(pairDCA<1 && massDCA<0.3){
			
		//..EMC cuts QA		
		if(TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24){
		  if(charge1*charge2<0){hDeltaZ_unlike->Fill(pT,deltaz);}
                  else{hDeltaZ_like->Fill(pT,deltaz);}
		}
		if(TMath::Abs(deltaz)<3 && 0.3<ratio && ratio<2 && nSMDE>1 && nSMDP>1 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24){
		  if((charge1*charge2)<0){hDeltaPhi_unlike->Fill(pT,deltaphi);}
		  else{hDeltaPhi_like->Fill(pT,deltaphi);}
		}
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24){
		  if((charge1*charge2)<0){hpE_unlike->Fill(pT,ratio);}
		  else{hpE_like->Fill(pT,ratio);}
		}	
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDP>1 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24 && 0.3<ratio && ratio<2){
		  if((charge1*charge2)<0){smde_unlike->Fill(pT,nSMDE);}
                  else{smde_like->Fill(pT,nSMDE);}
		}
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24 && 0.3<ratio && ratio<2){
                  if((charge1*charge2)<0){smdp_unlike->Fill(pT,nSMDP);}
                  else{smdp_like->Fill(pT,nSMDP);}
                }
		//..invariant mass plot
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 
			0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3){
		  if(charge1*charge2<0){hinvmass_pt_unlike->Fill(pT,massDCA);}
		  else{hinvmass_pt_like->Fill(pT,massDCA);}
		}
		//..global partner nsigma
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 
			0.3<ratio && ratio<2 && -0.5<nsigma_primary && nsigma_primary<2.5 && massDCA<0.24){
		  if(charge1*charge2<0){nsigma_partner_unlike->Fill(pT,nsigma_global);}
		  else{nsigma_partner_like->Fill(pT,nsigma_global);}
		}
		//..photonic electron pT spectrum
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24){
		  if(charge1*charge2<0){photonic_pt_unlike->Fill(pT);}
                  else{photonic_pt_like->Fill(pT);}
		}
		//..calibration && nsigma efficiency
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 
			0.3<ratio && ratio<2 && massDCA<0.01 && -1<nsigma_global && nsigma_global<3){
		  if(charge1*charge2<0){nsigma_primary_unlike->Fill(pT,nsigma_primary);}
                  else{nsigma_primary_like->Fill(pT,nsigma_primary);}
		}
		//..EMC efficiency - no EMC cuts
		if(massDCA<0.01 && -0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3){
		  if(charge1*charge2<0){noEMC_unlike->Fill(nsigma_global,pT);}
                  else{noEMC_like->Fill(nsigma_global,pT);}
		}
		//..EMC efficiency - with EMC cuts
		if(massDCA<0.01 && -0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && 
			TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2){
		  if(charge1*charge2<0){EMC_unlike->Fill(nsigma_global,pT);}
                  else{EMC_like->Fill(nsigma_global,pT);}
		}
		//..plots for gamma, eta, pi0 embedding && systematics
		if(TMath::Abs(deltaz)<3 && TMath::Abs(deltaphi)<0.015 && nSMDE>1 && nSMDP>1 && 0.3<ratio && ratio<2 && 
			-0.5<nsigma_primary && nsigma_primary<2.5 && TMath::Abs(nsigma_global)<3 && massDCA<0.24){
		  if(charge1*charge2<0){
			eta_pt_unlike->Fill(pT, pseudorapidity);
			phi_pt_unlike->Fill(pT, phi);
			dca_pt_unlike->Fill(pT, DCA);
			nHits_nHitsRatio_pt_partner_unlike->Fill(pT, nHits_global, hitsRatio_global);
                        pT_partner_unlike->Fill(pT, pT_global);
                        pair_dca_unlike->Fill(pT, pairDCA);
		  }
		  else{
			eta_pt_like->Fill(pT, pseudorapidity);
			phi_pt_like->Fill(pT, phi);
			dca_pt_like->Fill(pT, DCA);
			nHits_nHitsRatio_pt_partner_like->Fill(pT, nHits_global, hitsRatio_global);
                        pT_partner_like->Fill(pT, pT_global);
                        pair_dca_like->Fill(pT, pairDCA);
		  }
		}
					
	     }//..good pairs
	  }//..good global track
	}//..loop over global tracks
     }//..good primary track
   }//..loop over primray tracks
 }//if passes refMult cut from StRefMultCorr
}//if passes event cuts

cout<<"KG4"<<endl;
*/
 return kStOK;

}


#if TOFILLTREE
//---------------------------------------------------------------

Int_t StMuJetAnalysisTreeMaker::doPrimTrks(){
  
  TStarJetPicoPrimaryTrack *MpTrack = new TStarJetPicoPrimaryTrack;
  
  //Float_t eta,  phi;
  Int_t flagval; 
  Int_t nTracks= mCurrentMu->primaryTracks()->GetEntries();
  cout<<"doPrimTrks() .. nTracks="<<nTracks<<endl;
  
  //muEvent = (StMuEvent*) mu->event();
  //StThreeVectorF vertex=muEvent->primaryVertexPosition();
 
  Int_t PrimTrks=0; 
  Int_t ncHit;
  //  Int_t keyC;
  // Int_t matchflag;

  fPrimIndexArray=new Int_t[nTracks];
     

  for (Int_t l=0; l<nTracks; l++) { 
    StMuTrack*  track =   mCurrentMu->primaryTracks(l); 
    StMuTrack* gl_track = (StMuTrack*)track->globalTrack();

    assert(track);
    //    if(track->dcaGlobal().mag()<1.5 && track->eta()<1.5 && track->eta()>-1.5 &&  track->flag() >= 0)
    //cout<<"track=  "<<l<<"  "<<fMatchTrArr[l]<<endl;    
    //    if(fMatchTrArr[l]==3)continue; //track is matched to EMC

    fPrimIndexArray[l]=-10;

    flagval=track->flag();
    ncHit=track->nHitsFit();
    //QAHist->hTrackNhits->Fill(ncHit);
//    const double m2=511e-6*511e-6; //electron mass

	Float_t pt=gl_track->pt();   

    // if(flagval >=0 && ncHit>fNHits && track->dcaGlobal().mag()<2 && track->eta()<1.5 && track->eta()>-1.5){
    if(flagval >=0 && flagval <1000 &&ncHit>fNHits && TMath::Abs(track->eta())<fEta && pt>1.2){    //fNHits=12 a fEta=1.5
      //cout<<"        "<<eta<<"  "<<phi<<"  "<<flagval<<"  "<<ncHit<<endl;
     
      MpTrack->Clear();

       //BEMC----------------------------------------------------------
       int id; int adc0; float e[2]; float dist[2]; int nhit[2]; int ntow;
      if( getBEMC(track, &id, &adc0, e, dist, nhit, &ntow))
	{
       FillBEMC(MpTrack, &adc0, e, dist, nhit, &ntow);}     
       //BEMC------------------------------------------------------------

      /*QAHist->hTrackNhitsCut->Fill(ncHit);
      QAHist->hTrackDCA->Fill(track->dcaGlobal().mag());
      QAHist->hTrackEta->Fill(eta);
      QAHist->hTrackPhi->Fill(phi);*/
 

      MpTrack->SetPx(static_cast<Float_t>(gl_track->momentum().x()));    
      MpTrack->SetPy(static_cast<Float_t>(gl_track->momentum().y()));      
      MpTrack->SetPz(static_cast<Float_t>(gl_track->momentum().z())); 
	
      MpTrack->SetDCA(static_cast<Float_t>(gl_track->dcaGlobal().mag()));
      MpTrack->SetdEdx(static_cast<Float_t>(gl_track->dEdx()));
      
      //******************************************************************************
//******************************************************KATKA*************************
/*      MpTrack->SetHelix(track->helix());
      MpTrack->SetHelixCurvature(static_cast<Float_t>(track->helix().curvature())); /
      MpTrack->SetHelixDipAngle(static_cast<Float_t>(track->helix().dipAngle()));    /
      MpTrack->SetHelixPhase(static_cast<Float_t>(track->helix().phase())); */         //*
//**************************KATKA*****************************************************
//***********************************************************************************
  
      MpTrack->SetBTofPid(gl_track->btofPidTraits()); //saving for beta calculation

      /*if (MpTrack->TOF_MatchFlag()){
	double px=track->momentum().x();
	double py=track->momentum().y();
	double pz=track->momentum().z();
	double p2=px*px+py*py+pz*pz;
	double p=TMath::Sqrt(px*px+py*py+pz*pz);
	double b=MpTrack->TOF_Beta();
	double dBeta=1-b*TMath::Sqrt(1-m2/p2);
	QAHist->hElectronPID->Fill(static_cast<Float_t>(track->nSigmaElectron()),dBeta,p);
	}*/

      MpTrack->SetNsigPion(static_cast<Float_t>(gl_track->nSigmaPion()));
      MpTrack->SetNsigKaon(static_cast<Float_t>(gl_track->nSigmaKaon()));
      MpTrack->SetNsigProton(static_cast<Float_t>(gl_track->nSigmaProton()));
      MpTrack->SetNsigElectron(static_cast<Float_t>(gl_track->nSigmaElectron()));
	
      MpTrack->SetCharge(static_cast<Char_t>(gl_track->charge()));
      MpTrack->SetFlag(static_cast<Short_t>(gl_track->flag()));
      MpTrack->SetNOfFittedHits(static_cast<UShort_t>(gl_track->nHitsFit()));
      MpTrack->SetNOfPossHits(static_cast<UShort_t>(gl_track->nHitsPoss()));
      //MpTrack->SetKey(static_cast<Int_t>(gl_track->id()));//KATKA: before <Int_t>(l) //track->id())); //Petr - pozor ugly shortcut

      //cout<<"l="<<l<<"id="<<track->id()<<endl;;
    //Petr    MpTrack.SetChi2(static_cast<Float_t>(track->chi2xy())); // chi2/ndf of the final fit
   //Petr      MpTrack.SetChi2PV(static_cast<Float_t>(track->chi2z())); //extra chi2 when Primary Vertex is included into the fit

      //Petr  MpTrack.SetFlag(flagval);
      MpTrack->SetNOfDedxHits(static_cast<Int_t>(gl_track->nHitsDedx()));

      //signed DCAxy
    

	//filling vectors of ids------------------------------------------------------------------
//cout<<"idJetTreePrimIndex before photo candidate cut "<<MEvent->GetHeader()->GetNOfPrimaryTracks()<<endl;
        if (PhotoCandidate(track,gl_track)){
	idMuDstPrimIndex.push_back(l);
	idJetTreePrimIndex.push_back(MEvent->GetHeader()->GetNOfPrimaryTracks());
//	cout<<"idJetTreePrimIndex"<<MEvent->GetHeader()->GetNOfPrimaryTracks()<<endl;
	}
	//fillinf vectors of ids----------------------------------------------------------------------
        
 
 /*   const StMuTrack* gltrack = track->globalTrack();
      StThreeVectorF prVtx = muEvent->primaryVertexPosition();
      if (gltrack) {
 //Petr 	MpTrack.SetsDCAxy(computeXY(prVtx,gltrack->helix()));
      }
      else {
 //Petr 	MpTrack.SetsDCAxy(-99);
	cout<<"Error: Track without associated global track!"<<endl;
	}*/
//	cout<<"at the end of loop"<<endl;
	fPrimIndexArray[PrimTrks]=l;
	MEvent->AddPrimaryTrack(MpTrack);
      PrimTrks++;
      
    }//end of "if"
    
  }//end of "for" through tracks
  delete MpTrack;
//cout<<"idJetTreePrimIndex size = "<<idJetTreePrimIndex.size()<<endl;
//cout<<"idMuDstPrimIndex size = "<<idMuDstPrimIndex.size()<<endl;
  cout<<PrimTrks<<endl;
  return PrimTrks;
}

//----------------------------------------------------------------------------
Int_t StMuJetAnalysisTreeMaker::doGlobTrks(){
  
  //TStarJetPicoGlobalTrack MpTrack;
  TStarJetPicoGlobalTrack *MpTrack = new TStarJetPicoGlobalTrack;
  
//  Float_t eta,  phi;// px, py, pT;
  Int_t flagval; 
  Int_t nTracks= mCurrentMu->globalTracks()->GetEntries();
  //StMuEvent* muEvent=mCurrentMu->event();

  // nMTracks=fmtracklist->GetEntries();
  
  //muEvent = (StMuEvent*) mu->event();
  //StThreeVectorF vertex=muEvent->primaryVertexPosition();
 
  Int_t GlobTrks=0; 
  Int_t ncHit;
  //  Int_t keyC;
  // Int_t matchflag;

  fGlobIndexArray=new Int_t[nTracks];


  for (Int_t l=0; l<nTracks; l++) { 
    StMuTrack*  track =   mCurrentMu->globalTracks(l); 

    assert(track);
    //    if(track->dcaGlobal().mag()<1.5 && track->eta()<1.5 && track->eta()>-1.5 &&  track->flag() >= 0)
    //cout<<"track=  "<<l<<"  "<<fMatchTrArr[l]<<endl;    
    //    if(fMatchTrArr[l]==3)continue; //track is matched to EMC

    fGlobIndexArray[l]=-10;


    flagval=track->flag();
    ncHit=track->nHitsFit();
    //QAHist->hTrackNhits->Fill(ncHit);
//    const double m2=511e-6*511e-6; //electron mass
    
//	Float_t pt=track->pt();
        Double_t nHits = track->nHitsFit();
        Double_t posshits = track->nHitsPoss();
        Double_t Ratio = nHits/posshits;
    // if(flagval >=0 && ncHit>fNHits && track->dcaGlobal().mag()<2 && track->eta()<1.5 && track->eta()>-1.5){
    if(flagval >=0 && flagval <1000 &&ncHit>fNHits && Ratio>0.52 && TMath::Abs(track->eta())<fEta){
      //cout<<"        "<<eta<<"  "<<phi<<"  "<<flagval<<"  "<<ncHit<<endl;


 
      MpTrack->Clear();
  
     
/*	//BEMC----------------------------------------------------------
        int id; int adc0; float e[2]; float dist[2]; int nhit[2]; int ntow;
        if(getBEMC(track, &id, &adc0, e, dist, nhit, &ntow)){
         FillBEMC(MpTrack, &id, &adc0, e, dist, nhit, &ntow);}               
       //BEMC------------------------------------------------------------*/
       



/*      QAHist->hTrackNhitsCut->Fill(ncHit);
      QAHist->hTrackDCA->Fill(track->dcaGlobal().mag());
      QAHist->hTrackEta->Fill(eta);
      QAHist->hTrackPhi->Fill(phi);*/
 
     

     
      MpTrack->SetPx(static_cast<Float_t>(track->momentum().x()));      
      MpTrack->SetPy(static_cast<Float_t>(track->momentum().y()));      
      MpTrack->SetPz(static_cast<Float_t>(track->momentum().z()));
	 
//      MpTrack->SetDCA(static_cast<Float_t>(track->dcaGlobal().mag()));
//      MpTrack->SetdEdx(static_cast<Float_t>(track->dEdx()));

   //******************************************************************************
//******************************************************KATKA*************************
//      MpTrack->SetHelix(track->helix());
//      MpTrack->SetHelixCurvature(static_cast<Float_t>(track->helix().curvature())); //*
//      MpTrack->SetHelixDipAngle(static_cast<Float_t>(track->helix().dipAngle()));    //*
//      MpTrack->SetHelixPhase(static_cast<Float_t>(track->helix().phase()));          //*
//**************************KATKA*****************************************************
//***********************************************************************************

      MpTrack->SetBTofPid(track->btofPidTraits());

      /*if (MpTrack->TOF_MatchFlag()){
	double px=track->momentum().x();
	double py=track->momentum().y();
	double pz=track->momentum().z();
	double p2=px*px+py*py+pz*pz;
	double p=TMath::Sqrt(px*px+py*py+pz*pz);
	double b=MpTrack->TOF_Beta();
	double dBeta=1-b*TMath::Sqrt(1-m2/p2);
	QAHist->hElectronPID->Fill(static_cast<Float_t>(track->nSigmaElectron()),dBeta,p);
	}*/

//      MpTrack->SetNsigPion(static_cast<Float_t>(track->nSigmaPion()));
//      MpTrack->SetNsigKaon(static_cast<Float_t>(track->nSigmaKaon()));
//      MpTrack->SetNsigProton(static_cast<Float_t>(track->nSigmaProton()));
      MpTrack->SetNsigElectron(static_cast<Float_t>(track->nSigmaElectron()));
	
      MpTrack->SetCharge(static_cast<Char_t>(track->charge()));
      MpTrack->SetFlag(static_cast<Short_t>(track->flag()));
      MpTrack->SetNOfFittedHits(static_cast<UShort_t>(track->nHitsFit()));
//      MpTrack->SetNOfPossHits(static_cast<UShort_t>(track->nHitsPoss()));
//      MpTrack->SetKey(static_cast<Int_t>(track->id()));

    //Petr    MpTrack.SetChi2(static_cast<Float_t>(track->chi2xy())); // chi2/ndf of the final fit
   //Petr      MpTrack.SetChi2PV(static_cast<Float_t>(track->chi2z())); //extra chi2 when Primary Vertex is included into the fit

      //Petr  MpTrack.SetFlag(flagval);
//      MpTrack->SetNOfDedxHits(static_cast<Int_t>(track->nHitsDedx()));

//cout<<"idJetTreeGlobIndex before global partner cut"<<MEvent->GetHeader()->GetNOfGlobalTracks()<<endl;
	//filling vectors of ids---------------------------------------------
	if (GlobalPartnerCandidate(track)){    
         idMuDstGlobIndex.push_back(l);
         idJetTreeGlobIndex.push_back(MEvent->GetHeader()->GetNOfGlobalTracks());
//cout<<"idJetTreeGlobIndex"<<MEvent->GetHeader()->GetNOfGlobalTracks()<<endl;
//	cout<<"idJetTreeGlobIndex = "<<MEvent->GetHeader()->GetNOfGlobalTracks()<<endl;
	 }
//cout<<"idJetTreeGlobIndex size = "<<idJetTreeGlobIndex.size()<<endl;
	//filling vectors of ids-----------------------------------------------


      fGlobIndexArray[GlobTrks]=l;
      MEvent->AddGlobalTrack(MpTrack);
      GlobTrks++;
     
    }//end of "if"
      
  }//end of "for" throught tracks
  delete MpTrack;
//cout<<"idJetTreeGlobIndex size = "<<idJetTreeGlobIndex.size()<<endl;
//cout<<"idMuDstGlobIndex size = "<<idMuDstGlobIndex.size()<<endl;
  cout<<GlobTrks<<endl;
  return GlobTrks;
}

//-----------------------------------------------------------------
Int_t StMuJetAnalysisTreeMaker::doPhotoElectrons(){

TStarJetPicoPairs *Pair = new TStarJetPicoPairs;
const Double_t ELECMASS=0.00051099906;

cout<<"idMuDstPrimIndex size = "<<idMuDstPrimIndex.size()<<endl;
cout<<"idMuDstGlobIndex size = "<<idMuDstGlobIndex.size()<<endl;
//Int_t pairCount = 0;
//Int_t photoCount = 0;
 for(Int_t j=0; j<idMuDstPrimIndex.size() ;j++){
	StMuTrack*  track1 =  mCurrentMu->primaryTracks(idMuDstPrimIndex[j]);
	StMuTrack* gl_track1 = (StMuTrack*)track1->globalTrack();

	for(Int_t k=0; k<idMuDstGlobIndex.size()  ;k++){
	  //if (idMuDstPrimIndex[j] == idMuDstGlobIndex[k]) continue;
	  StMuTrack* track2 = mCurrentMu->globalTracks(idMuDstGlobIndex[k]);
	  if (track1->id() == track2->id()) continue;

	  Pair->Clear();

   	  StPhysicalHelixD* e1Helix = &(gl_track1->helix());
          StPhysicalHelixD* e2Helix = &(track2->helix());

          pair<double,double> s = e1Helix->pathLengths(*e2Helix);
          StThreeVectorD e1PosAtDca = e1Helix->at(s.first);
          StThreeVectorD e2PosAtDca = e2Helix->at(s.second);
	  StThreeVectorD PosAtDca = (e1PosAtDca + e2PosAtDca)/2;
          Double_t pairDCA = (e1PosAtDca - e2PosAtDca).mag();//pairDCA

          StThreeVectorD p1DCA=e1Helix->momentumAt(s.first,Bfield*tesla);
          StThreeVectorD p2DCA=e2Helix->momentumAt(s.second,Bfield*tesla);

          StLorentzVectorD LV1,LV2;
          LV1.setPx(p1DCA.x());
          LV1.setPy(p1DCA.y());
          LV1.setPz(p1DCA.z());
          LV1.setE(p1DCA.massHypothesis(ELECMASS));

          LV2.setPx(p2DCA.x());
          LV2.setPy(p2DCA.y());
          LV2.setPz(p2DCA.z());
          LV2.setE(p2DCA.massHypothesis(ELECMASS));

          Double_t massDCA=(LV1+LV2).m();

	  //pairCount++;

          if (pairDCA < 1 && massDCA < 0.3){
		//photoCount++;
             Pair->SetmassDCA(massDCA);
	     Pair->SetpairDCA(pairDCA);
	     Pair->SetPrimId(idJetTreePrimIndex[j]);
	     Pair->SetGlobId(idJetTreeGlobIndex[k]);

	     //Pair->SetPositionDCA1(e1PosAtDca);
	     //Pair->SetPositionDCA2(e2PosAtDca);
	     Pair->SetPositionDCA(PosAtDca);

	     MEvent->AddPair(Pair);
	     }
        }//loop over global tracks
  }//loop over primary tracks
//cout<<"pair count = "<<pairCount<<endl;
//cout<<"photo count = "<<photoCount<<endl;
delete Pair;
}
//----------------------------------------------------------------

bool StMuJetAnalysisTreeMaker::PhotoCandidate(StMuTrack *t, StMuTrack *gl_t){

bool hits=false, dedxhits=false, ratio=false, nsigma=false, eta=false, pt=false, point=false;

  if (gl_t->nHitsFit()>20) {hits = true;}
 // if (t->nHitsFit()<20) return false;  
//couti<<"nhits "<<t->nHitsFit()<<endl;

  if(gl_t->nHitsDedx()>11){dedxhits = true;}

  Double_t nHits = gl_t->nHitsFit();
  Double_t posshits = gl_t->nHitsPoss();
//cout<<"posshits = "<<t->nHitsPoss()<<endl;
  Double_t Ratio = nHits/posshits;
//cout<<"ratio before condition = "<<Ratio<<endl;
  if (Ratio > 0.52) {ratio = true;}
//cout<<"ratio "<<Ratio<<endl;


  if (TMath::Abs(gl_t->nSigmaElectron())<3) {nsigma = true;}
//  if (TMath::Abs(t->nSigmaElectron())>3) return false; 
//cout<<"nsigma "<<t->nSigmaElectron()<<endl;


  if(TMath::Abs(gl_t->eta()) < 0.7) {eta = true;}
//  if(TMath::Abs(t->eta()) > 0.7) return false;
//cout<<"ate "<<t->eta()<<endl;

  if(gl_t->pt() > 0.2) {pt = true;}
//  if(gl_t->pt() < 0.2) return false;
//cout<<"pt "<<t->pt()<<endl;

  StThreeVectorF firstpoint = gl_t->firstPoint();
//cout<<"first point vector "<<firstpoint.x()<<" "<<firstpoint.y()<<" "<<firstpoint.z()<<endl;
  double firstPoint_x = firstpoint.x();
//cout<<"first point x = "<<firstPoint_x<<endl;
  double firstPoint_y = firstpoint.y();
//cout<<"first point y = "<<firstPoint_y<<endl;
  double firstPoint = TMath::Sqrt(firstPoint_x*firstPoint_x + firstPoint_y*firstPoint_y);
//cout<<"first point before condition "<<firstPoint<<endl;
  if(firstPoint < 73) {point = true;}
//	cout<<"first pint "<<firstPoint<<endl;

  if (hits==true && dedxhits == true && ratio==true && nsigma==true && eta==true && pt==true && point==true){
//cout<<"passed photo candidate cut"<<endl;
//cout<<"hits, ratio, nsigma, eta, pt = "<<t->nHitsFit()<<" , "<<Ratio<<" , "<<t->nSigmaElectron()<<" , "<<t->eta()<<" , "<<gl_t->pt()<<endl;
	return true;
  }
  else {//cout<<"didnt passed photo candidate cut"<<endl;
	return false;}

  //float nsigEl=t->NsigElectron();
  //if (nsigEl<-1.0 || nsigEl>3) return false; 
  
  /*float nsigPi=t->NsigPion();
  if (nsigPi<2.3 && nsigPi>-2.3) return false; //2.5*/


/*
  UChar_t betaMatchFlag=t->TOF_MatchFlag();
  Float_t tofYLocal=t->TOF_YLocal();

 
  if (TMath::Abs(t->NsigElectron())>3) return false; 
 
  bool TOFveto=false;
  bool TOFid=false;
  bool TOFhit=(betaMatchFlag==1 && fabs(tofYLocal)<1.8);
  if (TOFhit ){
    Float_t inverseBeta = 1/t->TOF_Beta();
    if (inverseBeta>1.02||inverseBeta<0.98) TOFveto=true;
    else TOFid=true;
  } 

  
  double p=t->P();
  if (p>1.4){ //use BEMC
    double bemcE = t->GetBEMC_e();
    if (bemcE<=0.) return false;
    double ratio=p/bemcE;
    if (ratio>2. || ratio<0.3) return false;
    if (TOFveto) return false;
  }
  else if (!TOFid) return false;
  return true;
*/

  }
//----------------------------------------
bool StMuJetAnalysisTreeMaker::GlobalPartnerCandidate(StMuTrack *t){
	bool pt=false, nsigma=false, nfit=false, ratio=false;
	if (t->pt() > 0.2) pt = true;
	//cout<<"pt = "<<t->pt()<<endl;}
	//cout<<"pt no matter the cut = "<<t->pt()<<endl;
	if (TMath::Abs(t->nSigmaElectron())<5) nsigma = true;
	//cout<<"nsigmaelectorn = "<<t->nSigmaElectron()<<endl;}
	//cout<<"nsigmaelectorn no matter the cut = "<<t->nSigmaElectron()<<endl;
	if (t->nHitsFit()>20){nfit = true;}
	
	Double_t nHits = t->nHitsFit();
	Double_t posshits = t->nHitsPoss();
	Double_t Ratio = nHits/posshits;
        if (Ratio > 0.52) {ratio = true;}

	if(pt==true && nsigma==true && nfit==true && ratio==true) {
//	cout<<"passed global cut"<<endl;
//cout<<"pt, sigma = "<<t->pt()<<" , "<<t->nSigmaElectron()<<endl;
	return true;}
	else {return false;}
}
//-----------------------------------------------------
/*int StMuJetAnalysisTreeMaker::FindPhotoElectrons(){
  bool flag[5000];
  UInt_t n=MEvent->GetHeader()->GetNOfPrimaryTracks();
  assert(n<=5000);
  cout<<"FindPhotoElectron ntracks="<<n<<endl; 
  if (n<=1) return 0;
  for (uint i=0;i<5000;i++)flag[i]=false;

  for (uint i=0;i<n-1;i++){
    TStarJetPicoPrimaryTrack *trk=MEvent->GetPrimaryTrack(i);
    if (! PhotoCandidate(trk) continue;//,fBEMC[i])) continue;

    for (uint j=0;j<n-1;j++){
      if (i==j) continue;
      TStarJetPicoPrimaryTrack *trk2=MEvent->GetPrimaryTrack(j);
      if (!PhotoCandidate(trk2) continue;//,fBEMC[j])) continue;

      double px1=trk->Px();
      double py1=trk->Py();
      double pz1=trk->Pz();
      double px2=trk2->Px();
      double py2=trk2->Py();
      double pz2=trk2->Pz();
      double E1=TMath::Sqrt(0.261121e-6 + px1*px1 + py1*py1 + pz1*pz1);  //0.26e-6 is a (mass of electron)^2
      double E2=TMath::Sqrt(0.261121e-6 + px2*px2 + py2*py2 + pz2*pz2);
      
      double px=px1+px2;
      double py=py1+py2;
      double pz=pz1+pz2;
      double E=E1+E2;

      double minv=E*E-px*px-py*py-pz*pz;
      if (minv<0.02) {
	flag[i]=true;
        flag[j]=true; 
      }
    }//j
 
  } //i

  for (uint i=0;i<n;i++)
    if (flag[i]){
      TStarJetPicoPrimaryTrack *trk=MEvent->GetPrimaryTrack(i);
      MEvent->AddPhotoElectron(trk);
    }
  return MEvent->GetNOfPhotoEl();
  }*/
//--------------------------------------------------------------------------------
/*Int_t StMuJetAnalysisTreeMaker::doCheckMatchedTracks(bool primary=false){
  
  if (mVerbose)cout << "----------- In doCheckMatchedTracks() -----------------" << endl;
  
  //before looping on all the towers: create an array with the index of tracks
  // which have a projection to the BEMC
  
  Int_t nhitcheck,flagval;
  Int_t nTr; 
  Int_t *indexArray;
  if (primary) {
    nTr=MEvent->GetHeader()->GetNOfPrimaryTracks();
    indexArray=fPrimIndexArray;
    cout<<"doCheckMatchedTracks for PrimTrks"<<endl;
  }
  else {
    nTr=MEvent->GetHeader()->GetNOfGlobalTracks();
    indexArray=fGlobIndexArray;
    cout<<"doCheckMatchedTracks for GlobTrks"<<endl;
  }
 
 //  nTracks= mCurrentMu->primaryTracks()->GetEntries(); 
 //  const Int_t nTr=nTracks;
  
  fMatchTrArr=new Int_t[nTr];
  fMatchTrEtaArr=new Float_t[nTr];
  fMatchTrPhiArr=new Float_t[nTr];
  fEtaArray=new Float_t[nTr];
  fPhiArray=new Float_t[nTr];

  for(Int_t i=0;i<nTr;i++){
    fMatchTrArr[i]=-10;
    fMatchTrEtaArr[i]=-10.;
    fMatchTrPhiArr[i]=-10.;
    fEtaArray[i]=-10.;
    fPhiArray[i]=-10.;
  }


  Int_t nMatched=0,n1=0,n2=0,n3=0;
  StEmcPosition mPosition;
  for(Int_t tr =0; tr<nTr; tr++) 
    {//loop over tracks
      Int_t l=indexArray[tr];
      StMuTrack *track;
      if (primary)track=mCurrentMu->primaryTracks(l); 
      else track=mCurrentMu->globalTracks(l); 
      assert(track);
      nhitcheck=track->nHitsFit();
      StThreeVectorD momentum,position;
      flagval=track->flag();

      nhitcheck=track->nHitsFit();
      n1++;   
      if(flagval<0 || nhitcheck<fNHits+1) continue;
      n2++;

      bool proj_ok= mPosition.projTrack(&position,&momentum,track,(Double_t) Bfield); 
      
      Float_t z,eta,phi;
      eta =position.pseudoRapidity(); 
      phi =position.phi();
      z   =position.z();
      Int_t m,e,s;
      Int_t okbin=mGeom->getBin(phi,eta,m,e,s);
      if (s <0 || okbin==1)continue;
      n3++;
      
      if (proj_ok){
	fMatchTrArr[tr]=1;
	fMatchTrEtaArr[tr]=eta;
	fMatchTrPhiArr[tr]=phi;
	fEtaArray[tr]=track->eta();
	fPhiArray[tr]=track->phi();
	nMatched++;
      }

    }


  if(mVerbose)printf("track # = %d and projected to BEMC # = %d \n",mCurrentMu->primaryTracks()->GetEntries(), nMatched);
  if(mVerbose)cout<<"n1,n2,n3="<<n1<<"  "<<n2<<"  "<<n3<<endl;

  return nMatched;
}*/
//--------------------------------------------------------------------
/*Int_t StMuJetAnalysisTreeMaker::doV0s(){

  TStarJetPicoV0 MpV0;
  StV0MuDst *v0MuDst;
  StMuTrack *globTrk;
  Int_t v0=0;
  Int_t posFound, negFound;

  Int_t nV0sMuDst = mCurrentMu->v0s()->GetEntries();
  
  for (Int_t l=0; l<nV0sMuDst; l++) { 
    v0MuDst =  (StV0MuDst*) mCurrentMu->v0s(l); 

    if( fTCutV0->CheckV0(v0MuDst,mCurrentMu ) == 1){
      MpV0.Clear();
      MpV0.SetPxPos(static_cast<Float_t>(v0MuDst->momPosX()));
      MpV0.SetPyPos(static_cast<Float_t>(v0MuDst->momPosY()));
      MpV0.SetPzPos(static_cast<Float_t>(v0MuDst->momPosZ()));
      MpV0.SetPxNeg(static_cast<Float_t>(v0MuDst->momNegX()));
      MpV0.SetPyNeg(static_cast<Float_t>(v0MuDst->momNegY()));
      MpV0.SetPzNeg(static_cast<Float_t>(v0MuDst->momNegZ()));
      MpV0.SetKeyPos(static_cast<Int_t>(v0MuDst->keyPos()));
      MpV0.SetKeyNeg(static_cast<Int_t>(v0MuDst->keyNeg()));
      MpV0.SetDcapn(static_cast<Int_t>(v0MuDst->dcaV0Daughters()));
      MpV0.SetDcaV0(static_cast<Float_t>(v0MuDst->dcaV0ToPrimVertex()));
      MpV0.SetDcaPos(static_cast<Float_t>(v0MuDst->dcaPosToPrimVertex()));
      MpV0.SetDcaNeg(static_cast<Float_t>(v0MuDst->dcaNegToPrimVertex()));
      MpV0.SetDecayLength(static_cast<Float_t>(v0MuDst->decayLengthV0()));
      MpV0.SetDedxPos(static_cast<Float_t>(v0MuDst->dedxPos()));
      MpV0.SetDedxNeg(static_cast<Float_t>(v0MuDst->dedxNeg()));

      posFound = 0;
      negFound = 0;
      for(Int_t glob=0; glob<mCurrentMu->globalTracks()->GetEntries(); glob++){
	globTrk =  mCurrentMu->globalTracks(glob);
	if( globTrk->id() == v0MuDst->keyPos()){
	  posFound = 1;
	  MpV0.SetNSigmaPionPos(static_cast<Float_t>(globTrk->nSigmaPion()));
	  MpV0.SetNSigmaProtonPos(static_cast<Float_t>(globTrk->nSigmaProton()));
	}
	if( globTrk->id() == v0MuDst->keyNeg()){
	  negFound = 1;
	  MpV0.SetNSigmaPionNeg(static_cast<Float_t>(globTrk->nSigmaPion()));
	  MpV0.SetNSigmaProtonNeg(static_cast<Float_t>(globTrk->nSigmaProton()));
	}
	if( negFound==1 && posFound ==1) break;
      }

      MEvent->AddV0(&MpV0);
      v0++;
    }
  }

  if(mVerbose)cout << "saving " << v0 << " out of " << nV0sMuDst << endl;
  return v0;  

}*/
//--------------------------------------------------------------------
Double_t StMuJetAnalysisTreeMaker::GetReactionPlane()
{ 
 if(mVerbose)cout << "----------- In GetReactionPlane() -----------------" << endl;
 TVector2 mQ;
 Double_t mQx=0., mQy=0.;
 Double_t order = 2; Int_t n,i;
 n=mCurrentMu->primaryTracks()->GetEntries();
 if (n==0) return 0.;
 for (i=0; i<n; i++) {
    StMuTrack* track =  mCurrentMu->primaryTracks(i);   // get pointer to primary track
   Float_t phi = track->phi();
   mQx += cos(phi * order);
   mQy += sin(phi * order);
 }
 
 mQ.Set(mQx, mQy);
 Double_t psi= mQ.Phi() / order;
 Float_t pi=TMath::Pi();
 return psi*180/pi;
}

//------------------------------------------------------------------
/*Int_t StMuJetAnalysisTreeMaker::doTowerMatching(bool primary=false){
  //true - match primary tracks, false - match secondary tracks

  if(mVerbose) cout << "----------- In doTowerMatching() -----------------" << endl;
  
  doCheckMatchedTracks(primary);

  Int_t nTrk; Int_t *indexArray;
  if (primary) {
    nTrk=MEvent->GetHeader()->GetNOfPrimaryTracks();
    indexArray=fPrimIndexArray;
  }
  else {
    nTrk=MEvent->GetHeader()->GetNOfGlobalTracks();
    indexArray=fGlobIndexArray;
  }
 
  StBemcTables* tables = mTables;
  assert(tables);
  Int_t flag=0;
  mEmcCol = (StEmcCollection*)mCurrentMu->emcCollection();
  StMuEmcCollection* muEmcCol = mCurrentMu->muEmcCollection();
  
  if(!mEmcCol || !muEmcCol)
    {
      printf("\n***-- no EMC Collection was found --***\n"); 
      return 1;
    }
    
 
  StSPtrVecEmcPoint& container =  mEmcCol->barrelPoints();
  
  if((mEmcCol->barrelPoints().size()==0 || container.size()==0)){
    cout<<"In doTrack Matching no container or no size"<<endl;
    return -1;
  }
  TStarJetPicoTower allTower;
  TStarJetPicoPrimaryTrack mTrack;
  //  allTower.SetObjArray();
 
  rBarrelPts = mEmcCol->barrelPoints().size();
    
  StEmcDetector* detector;
  
  detector=mEmcCol->detector(kBarrelEmcTowerId);
 

  Float_t etaTower, phiTower;
  if(detector) {
    Int_t mm = 0, ee = 0, ss = 0;	 
    Int_t statusAll  =-99;
    Int_t towerId = 0;
    fMatchedTow=0;

    // cout<<" counts   "<<fmtowerlist->GetEntries()<<endl;

    for (UInt_t i = 1; i < 121; i++){
      assert( detector->module(i));
      StSPtrVecEmcRawHit& emcTowerHits = detector->module(i)->hits();
      for (UInt_t j = 0; j < emcTowerHits.size(); j++){ 
	mm = (Int_t)emcTowerHits[j]->module();
	ee = (Int_t)emcTowerHits[j]->eta();
	ss = emcTowerHits[j]->sub();
	
	if(!(abs(mm)<=120 && abs(ee)<=20&&ss<=2)) continue;
	mGeom->getId(mm, ee, ss, towerId);
	tables->getStatus(BTOW, towerId, statusAll);
	
	if(statusAll!=1) continue;
	 
	//  Int_t matchflagtow;
	Float_t energyTowerAll=emcTowerHits[j]->energy();
	Float_t ADCTower=emcTowerHits[j]->adc();
	  
	if(energyTowerAll<0.15)continue;//added new Fri Oct31 08
	mGeom->getEtaPhi(towerId,etaTower,phiTower);
	
	  // 	    if(matchflagtow==1 || energyTowerAll<=0.) continue;
	  //   if(energyTowerAll<=0) continue;
	//Int_t didT=towerId;
	
	
	Int_t ehits,phits,clusterid;
	StEmcCluster* closeclust1 = findSMDCluster(etaTower,phiTower,2,clusterid);//SmdE
	if (closeclust1 !=NULL) {
	    ehits = SMDHits(2,clusterid);
	}
	StEmcCluster* closeclust2 = findSMDCluster(etaTower,phiTower,3,clusterid);//smdP
	if (closeclust2 !=NULL) {
	  phits = SMDHits(3,clusterid);
	}

	//Correction of PVZ shift
	
	Float_t etaC, T1;
	T1=2*atan(exp(-etaTower));
	Double_t zNew;
	if(etaTower!=0){zNew=231/tan(T1);} //231 cm = radius SMD
	if(etaTower==0){zNew=0;}
	Double_t zNom=zNew-PrimVertexZ;
	Double_t THETA=atan2(231,zNom);
	etaC=-log(tan(THETA/2));
	
	
	allTower.Clear();

	allTower.SetPhi(static_cast<Float_t>(phiTower));
	allTower.SetEta(static_cast<Float_t>(etaTower));
	allTower.SetId(static_cast<Int_t>(towerId));
	allTower.SetEnergy(static_cast<Float_t>(energyTowerAll));
	allTower.SetADC(static_cast<Int_t>(ADCTower));
	allTower.SetPhiCorrected(static_cast<Float_t>(phiTower)); //Not really corrected as no change is expected
	allTower.SetEtaCorrected(static_cast<Float_t>(etaC));
	
	allTower.SetSMDClusterP(static_cast<Int_t>(phits));
	allTower.SetSMDClusterE(static_cast<Int_t>(ehits));
	allTower.SetTowerStatus(static_cast<Int_t>(statusAll));
	Double_t ettower=energyTowerAll/TMath::CosH(etaC);
	  
	QAHist->hTowerEta->Fill(etaTower);
	QAHist->hTowerPhi->Fill(phiTower);
	
	QAHist->hTowerEtaC->Fill(etaC);
	QAHist->hTowerET->Fill(ettower);
	
	  
	Int_t count=0;
	assert(MEvent->GetHeader());
	
	for(Int_t tr =0; tr<nTrk; tr++) //this is nCand
	    {
	      if(fMatchTrArr[tr]!=1) continue;
	      
	      Float_t etaE,phiE,eta,phi;
	      etaE=fMatchTrEtaArr[tr];
	      phiE=fMatchTrPhiArr[tr];
	      eta=fEtaArray[tr];
	      phi=fPhiArray[tr];
	      
	      Float_t dphi = phiE - phiTower;
	      Float_t deta = etaE - etaTower;
	      
	      QAHist->hDeltaEtaDeltaPhiProj->Fill(deta,dphi);
	      
	      //cout<< " dphi  "<<dphi<<"  deta  "<<deta<<endl;
	      if (fabs(dphi) > 0.025) continue;
	      if (fabs(deta) > 0.025) continue;
	      QAHist->hDeltaEtaDeltaPhiProjAcc->Fill(deta,dphi);
	      
	      Int_t l=indexArray[tr];
	      StMuTrack* track;
	      if (primary)track=mCurrentMu->primaryTracks(l); 
	      else track=mCurrentMu->globalTracks(l); 
	      assert(track);
	      
	      Double_t trackorphi=track->phi();
	      Double_t trackoreta=track->eta();
	      Double_t diffeta=trackoreta-eta;
	      Double_t diffphi=trackorphi-phi;
	      
		
	      if(diffeta!=0 || diffphi!=0) {
		cout<<"==== DoTowerMatching ==== Sanity Check: Problem with the matched index "<<endl;
		QAHist->hEtaDiff->Fill(diffeta);
		QAHist->hPhiDiff->Fill(diffphi);
		  
		}
	      
		
	      Double_t pttrack=track->pt();
	      //Float_t POverE=-99;
	      
	      //POverE=TStarJetPicoUtils::GetTowerPoverE(MEvent,&allTower,tr);
	      //mTrack.Clear();
	      //mTrack.SetPx(track->momentum().x());
	      //mTrack.SetPy(track->momentum().y());
	      //mTrack.SetPz(track->momentum().z());
	      //POverE=TStarJetPicoUtils::GetTowerPoverE(&allTower,&mTrack);
	      
	      Float_t Pz=track->momentum().z();
	      Float_t POverE=TMath::Sqrt(pttrack*pttrack+Pz*Pz)/energyTowerAll;
	      //cout<<"POverE   "<<POverE<<"    "<<POverEMine<<endl;
	      QAHist->hTrackPoverE->Fill(POverE);
	      QAHist->hPtvsEtMatched->Fill(pttrack,ettower);
	      QAHist->hEtaMEtaT->Fill(eta,etaTower);
	      QAHist->hPhiMPhiT->Fill(phi,phiTower);
	      QAHist->hEtaDeltaEta->Fill(eta,deta);
	      QAHist->hPhiDeltaPhi->Fill(phi,dphi);
		  
		    
	      mTrack.SetEtaDiffHitProjected(static_cast<Float_t>(deta));
	      mTrack.SetPhiDiffHitProjected(static_cast<Float_t>(dphi));
	      
	      count++;
	      fNTOTMatchedTr++;
	      allTower.AddMatchedTrack(tr); 
	      fMatchTrArr[tr]=3;
	      //delete mPosition;
	      
	    }
	  
	  //Petr Chaloupka - hack to save only towers with associated tracks
	  if (allTower.GetNAssocTracks()>0){
	    MEvent->AddTower(&allTower);
	    flag++;
	  }
	  if(allTower.GetNAssocTracks()>0) fMatchedTow++;
      }
    }  
  }
  
  else {
    printf("no detector in EMC collection id \n");
  }
  
  if(mVerbose)cout<<"number of total associated tracks in the event = "<<fNTOTMatchedTr<<"  N of Towers with at least 1 matched track = "<<fMatchedTow<<endl;
  if(mVerbose)printf("Out StMuJetAnalysisTreeMaker:doMatchedTower \n");

   //POZOR tady je neco divneho - nekde asi tece pamet
  
  cout<<".. doTowerMatching DONE"<<endl;
  return flag;
}*/
//-----------------------------------------------------------------------------------
/*
StEmcCluster* StMuJetAnalysisTreeMaker::findSMDCluster(Float_t eta, Float_t phi,Int_t det, Int_t &clusterid)
{
  if(mVerbose==2)cout << "----------- In findSMDCluster() -----------------" << endl;
  //looking for corresponding cluster in SMD detector (det=2 SMDe, det=3 SMDp) with eta,phi of
  // TPC track
  StEmcCluster *clusterSMD = NULL;
  clusterid = -1;
  Float_t dmin = 9999;
  Float_t dmin_cut = 0.05; //this is cut on difference of SMD cluster eta,phi difference from TPC eta,phi 
  
  StDetectorId id = static_cast<StDetectorId>(det+kBarrelEmcTowerId);
  StEmcDetector* detector=mEmcCol->detector(id);
  if(detector) {
    StEmcClusterCollection* clusters = detector->cluster();
    if(clusters) {
      StSPtrVecEmcCluster& cl=clusters->clusters();
      //if(mVerbose) printf("number of SMD det %i clusters %i\n",det,cl.size());
      if(cl.size()>0) {
	for(Int_t i=0;i<(Int_t)cl.size();i++) if(cl[i]) { 
	  //going over all cluster
	  Float_t ETA = cl[i]->eta();
	  Float_t PHI = cl[i]->phi();
	  Float_t d = sqrt((eta-ETA)*(eta-ETA)+(phi-PHI)*(phi-PHI));
	  
	  if(d<dmin) {
	    //storing clusters with up to now best matching
	    dmin = d;
	    clusterid = i; 
	  }
	}//loop over clusters    
      }
      
      if(dmin<dmin_cut && clusterid!=-1) {
	clusterSMD = cl[clusterid];
	if(mVerbose==2)cout <<"Found SMD cluster for detector "<<det<<"\t distance="<< dmin << endl;
      
      }
    }//if clusters !NULL
  }//if detector !NULL
  
  return clusterSMD;
}

//------------------------------------------------------------
Int_t  StMuJetAnalysisTreeMaker::SMDHits(Int_t det, Int_t clustid)
{
  //looking for corresponding cluster in SMD detector (det=2 SMDe, det=3 SMDp) with eta,phi of
  Int_t nhits=0;
  StDetectorId id = static_cast<StDetectorId>(det+kBarrelEmcTowerId);
 StEmcDetector* detector=mEmcCol->detector(id);
  if(detector) {
   StEmcClusterCollection* clusters = detector->cluster();
    if(clusters) {
      StSPtrVecEmcCluster& cl=clusters->clusters();
      if (cl[clustid]) nhits = cl[clustid]->nHits();
     }    
    }
  return nhits;
}
*/
//------------------------------------------------------------------------
// same math as in $STAR/StRoot/StMiniMcMaker/StMiniMcMaker.cxx
/*Float_t StMuJetAnalysisTreeMaker::computeXY(const StThreeVectorF& pos, const StPhysicalHelixD &helix)
{
  // find the distance between the center of the circle and pos.                                                              
  // if this greater than the radius of curvature, then call                                                                                                           
  // it negative.                                                                                                                                                    
                                                                                                                                                                        
  double xCenter = helix.xcenter();
  double yCenter = helix.ycenter();
  double radius  = 1.0/helix.curvature();

  double dPosCenter
    = TMath::Sqrt( (pos.x()-xCenter) * (pos.x()-xCenter) +
                   (pos.y()-yCenter) * (pos.y()-yCenter));

  return (Float_t) radius - dPosCenter;
}*/
#endif
//----------------------------------------------------------------------------------
bool StMuJetAnalysisTreeMaker::getBEMC(StMuTrack *t, int *id, int *adc, float *ene, float *d, int *nep, int *towid)
{
  *id = -1; *adc = 0;
 
  for(int i=0;i<2;i++) { ene[i] = 0.; }
  for(int i=0;i<2;i++) { d[i] = 1.e9; }
  for(int i=0;i<2;i++) { nep[i] = 0; }
  *towid = -1;  
//for(int i=0;i<3;i++) { towid[i] = -1; }

bool assoc_ok = false;

  if(!mEmcCol) {
    LOG_WARN << " No Emc Collection for this event " << endm;
    return kFALSE;
  }

  StThreeVectorD position, momentum;
  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;
  //Double_t bFld = mBField*kilogauss/tesla; // bFld in Tesla
  bool ok       = false;
  bool okBSMDE  = false;
  bool okBSMDP  = false;
  if(mEmcPosition) {
    ok      = mEmcPosition->projTrack(&position, &momentum, t, Bfield, mEmcGeom[0]->Radius());
    okBSMDE = mEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, t, Bfield, mEmcGeom[2]->Radius());
    okBSMDP = mEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, t, Bfield, mEmcGeom[3]->Radius());
  }
//  if(!ok || !okBSMDE || !okBSMDP) 
  if(!ok) {
    LOG_WARN << " Projection failed for this track ... " << endm;
    return kFALSE;
  }
if(ok && okBSMDE && okBSMDP) {

  Int_t mod, eta, sub;
  StSPtrVecEmcPoint& bEmcPoints = mEmcCol->barrelPoints();
  int index=0;
  float mindist=1.e9;
  mEmcGeom[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub); //project on SMD plan
  for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++, index++) {
    bool associated=false;
    StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
    if(bEmcClusters.size()==0 ) continue;
    if(bEmcClusters[0]==NULL) continue;
    for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++) {
      StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
      for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) {
        if(mod == (Int_t)(*hIter)->module() && eta == (Int_t)(*hIter)->eta() && sub == (Int_t)(*hIter)->sub()) {
          associated=true;
          break;
        }
      }
      if(associated) {
        for(StPtrVecEmcRawHitIterator hitit=bEmcHits.begin(); hitit!=bEmcHits.end();hitit++) {
          if((*hitit)->energy()>ene[0]) {ene[0]=(*hitit)->energy();
	  				*towid = (*hitit)->softId(1);}
          if((int)(*hitit)->adc()>(*adc)) *adc=(*hitit)->adc();
        }
      }
    }
    StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
    StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);

    if(associated) {
      assoc_ok=true;
      *id = index;
      ene[1] = ene[1] + (*it)->energy(); //use point's energy, not tower cluster's energy

      float deltaphi=(*it)->position().phi()-positionBSMDP.phi();
      if(deltaphi>=TMath::Pi()) deltaphi=deltaphi-TMath::TwoPi();
      if(deltaphi<-TMath::Pi()) deltaphi=deltaphi+TMath::TwoPi();

      float rsmdp=mEmcGeom[3]->Radius();
      float pointz=(*it)->position().z();
      float deltaz=pointz-positionBSMDE.z();
      if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz)<mindist) {
        d[1]=deltaphi;
        d[0]=deltaz;
        if(smdeClusters.size()>=1) nep[0]=smdeClusters[0]->nHits();
        if(smdpClusters.size()>=1) nep[1]=smdpClusters[0]->nHits();
        mindist=sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz);
      }
    }//associated
  }

  } // end if (ok && okBSMDE && okBSMDP)

/* int towerId;
 mEmcGeom[0]->getId(position.phi(),position.pseudoRapidity(),towerId);
 towid = &towerId;*/

if(assoc_ok) return kTRUE;
else return kFALSE; //...added by Katka (new compilator claimed error because of not good return value)
}
//------------------------------------------------------------------------
void StMuJetAnalysisTreeMaker::FillBEMC(TStarJetPicoPrimaryTrack *track, int *adc, float *ene, float *d, int *nep, int *towid)
{
  
  track->SetBEMC_e0(ene[0]);
  track->SetBEMC_e(ene[1]);
//cout<<"e0 = "<<ene[0]<<endl;
//cout<<"e = "<<ene[1]<<endl;
  track->SetBEMC_adc(*adc);
  track->SetBEMC_towerid(*towid);
//  track->Set_emcPoint_id(*id);

  track->SetSMDE_Nhits(nep[0]);
  track->SetSMDP_Nhits(nep[1]);
  track->Set_emcdz(d[0]);
  track->Set_emcdphi(d[1]);
//cout<<"emcdz = "<<d[0]<<endl;
//cout<<"emcdphi = "<<d[1]<<endl;
}
//-----------------------------------------------------------------------------
void StMuJetAnalysisTreeMaker::initEmc()
{
  mEmcPosition = new StEmcPosition();
  for(int i=0;i<4;i++) {
	mEmcGeom[i] = StEmcGeom::getEmcGeom(detname[i].Data());
	}
}
//-------------------------------------------------------------------------
void StMuJetAnalysisTreeMaker::finishEmc()
{
  if(mEmcPosition) delete mEmcPosition;
  for(int i=0;i<4;i++) {
	mEmcGeom[i] = 0;
	}
}
//---------------------------------------------------------------------------
void StMuJetAnalysisTreeMaker::HistAllocation(){
Vz_before = new TH1F ("Vz_before","Vz_before",600,-300,300);
Vz_before->Sumw2();
DiffVz_before = new TH1F ("DiffVz_before","DiffVz_before",50,-5,5);
DiffVz_before->Sumw2();

hVz = new TH1F ("Vz","Vz",60,-30,30);
hVz->Sumw2();
hDiffVz = new TH1F("DiffVz","DiffVz",50,-3,3);
hDiffVz->Sumw2();

refMultiplicity = new TH1F("refMultiplicity","refMultiplicity;refMult;Counts",800,0,800);
refMultiplicity->Sumw2();
//centrality = new TH1F("centrality","centrality",9,0,9);
//centrality->Sumw2();
//mult_cent = new TH1F("mult_cent","mult_cent",800,0,800);
//mult_cent->Sumw2();

hdEdx = new TH2F("hdEdx","hdEdx;p_{T};dEdx",100,0,10,100,0,10);
hdEdx->Sumw2();
hnsigma = new TH2F("hnsigma","hnsigma;p_{T};n#sigma_{e}",100,0,10,200,-10,10);
hnsigma->Sumw2();

firstTPC = new TH2F("firstTPC","firstTPC;p_{T};1st TPC point",100,0,10,150,0,150);
ndEdx_pt = new TH2F("ndEdx_pt","ndEdx_pt;p_{T}; ndEdx hits",100,0,10,50,0,50);
nHits_nMax_pt = new TH3F("nHits_nMax_pt","nHits_nMax_pt;p_{T};nHits;nHits/nMax",100,0,10,50,0,50,50,0,1);

electron_pt = new TH1F ("electron_pt","electron_pt;p_{T} [GeV/c]",100,0,10);
electron_pt->Sumw2();
nsigma_purity_pt = new TH2F ("nsigma_purity","nsigma_purity;p_{T} [GeV/c];n#sigma_{e}",100,0,10,100,-10,10);
nsigma_purity_pt->Sumw2();

hpE = new TH2F("hpE","hpE;p_{T};p/E0",100,0,10,60,0,4);
hpE->Sumw2();
hpE_unlike = new TH2F ("pE_unlike","pE_unlike;p_{T} [GeV/c]; p/E0",100,0,10,60,0,4);
hpE_unlike->Sumw2();
hpE_like = new TH2F ("pE_like","pE_like;p_{T} [GeV/c]; p/E0",100,0,10,60,0,4);
hpE_like->Sumw2();
hDeltaPhi = new TH2F("hDeltaPhi","hDeltaPhi;p_{T};#Delta_{#varphi} [rad]",100,0,10,60,-0.1,0.1);
hDeltaPhi->Sumw2();
hDeltaPhi_unlike = new TH2F ("DeltaPhi_unlike","DeltaPhi_unlike;p_{T} [GeV/c];#Delta_{#varphi} [rad]",100,0,10,60,-0.1,0.1);
hDeltaPhi_unlike->Sumw2();
hDeltaZ = new TH2F("hDeltaZ","hDeltaZ;p_{T};#Delta_{z} [cm]",100,0,10,60,-8,8);
hDeltaZ->Sumw2();
hDeltaZ_unlike = new TH2F ("DeltaZ_unlike","DeltaZ_unlike;p_{T} [GeV/c];#Delta_{z} [cm]",100,0,10,60,-8,8);
hDeltaZ_unlike->Sumw2();
hDeltaPhi_like = new TH2F ("DeltaPhi_like","DeltaPhi_like;p_{T} [GeV/c];#Delta_{#varphi} [rad]",100,0,10,60,-0.1,0.1);
hDeltaPhi_like->Sumw2();
hDeltaZ_like = new TH2F ("DeltaZ_like","DeltaZ_like;p_{T} [GeV/c];#Delta_{z} [cm]",100,0,10,60,-8,8);
hDeltaZ_like->Sumw2();
hnSMDE = new TH2F("hnSMDE","hnSMDE;p_{T};nSMDE",100,0,10,15,0,5);
hnSMDE->Sumw2();
smde_unlike = new TH2F("smde_unlike","smde_unlike;p_{T} [GeV/c];SMDE hits",100,0,10,15,0,5);
smde_unlike->Sumw2();
smde_like = new TH2F("smde_like","smde_like;p_{T} [GeV/c];SMDE hits",100,0,10,15,0,5);
smde_like->Sumw2();
hnSMDP = new TH2F("hnSMDP","hnSMDP;p_{T};nSMDP",100,0,10,15,0,5);
hnSMDP->Sumw2();
smdp_unlike = new TH2F("smdp_unlike","smdp_unlike;p_{T} [GeV/c];SMDP hits",100,0,10,15,0,5);
smdp_unlike->Sumw2();
smdp_like = new TH2F("smdp_like","smdp_like;p_{T} [GeV/c];SMDP hits",100,0,10,15,0,5);
smdp_like->Sumw2();

hinvmass_pt_unlike = new TH2F ("invmass_pt_unlike","invmass_pt_unlike;p_{T} [GeV/c];m_{ee} [GeV/c^2]",100,0,10,60,0,0.3);
hinvmass_pt_unlike->Sumw2();
hinvmass_pt_like = new TH2F ("invmass_pt_like","invmass_pt_like;p_{T} [GeV/c];m_{ee} [GeV/c^2]",100,0,10,60,0,0.3);
hinvmass_pt_like->Sumw2();
nsigma_partner_unlike = new TH2F ("nsigma_partner_unlike","nsigma_partner_unlike;p_{T} [GeV/c];n#sigma_e",100,0,10,60,-5,5);
nsigma_partner_unlike->Sumw2();
nsigma_partner_like = new TH2F ("nsigma_partner_like","nsigma_partner_like;p_{T} [GeV/c];n#sigma_e",100,0,10,60,-5,5);
nsigma_partner_like->Sumw2();
        
nsigma_primary_unlike = new TH2F ("nsigma_primary_unlike","nsigma_primary_unlike;p_{T} [GeV/c];n#sigma_e",100,0,10,60,-5,5);
nsigma_primary_unlike->Sumw2();
nsigma_primary_like = new TH2F ("nsigma_primary_like","nsigma_primary_like;p_{T} [GeV/c];n#sigma_e",100,0,10,60,-5,5);
nsigma_primary_like->Sumw2();

photonic_pt_unlike = new TH1F("photonic_pt_unlike","photonic_pt_unlike;p_{T} [GeV/c]",100,0,10);
photonic_pt_unlike->Sumw2();
photonic_pt_like = new TH1F("photonic_pt_like","photonic_pt_like;p_{T} [GeV/c]",100,0,10);
photonic_pt_like->Sumw2();

EMC_unlike = new TH2F("EMC_unlike","EMC_unlike;n#sigma_{e};p_{T} [GeV/c]",60,-5,5,100,0,10);
EMC_unlike->Sumw2();
EMC_like = new TH2F("EMC_like","EMC_like;n#sigma_{e};p_{T} [GeV/c]",60,-5,5,100,0,10);
EMC_like->Sumw2();
noEMC_unlike = new TH2F("noEMC_unlike","noEMC_unlike;n#sigma_{e};p_{T} [GeV/c]",60,-5,5,100,0,10);
noEMC_unlike->Sumw2();
noEMC_like = new TH2F("noEMC_like","noEMC_like;n#sigma_{e};p_{T} [GeV/c]",60,-5,5,100,0,10);
noEMC_like->Sumw2();

eta_pt = new TH2F("eta","eta;p_{T} [GeV/c];eta",100,0,10,20,-1,1);
eta_pt->Sumw2();
eta_pt_unlike = new TH2F("eta_unlike","eta_unlike;p_{T} [GeV/c];eta",100,0,10,20,-1,1);
eta_pt_unlike->Sumw2();
eta_pt_like = new TH2F("eta_like","eta_like;p_{T} [GeV/c];eta",100,0,10,20,-1,1);
eta_pt_like->Sumw2();
phi_pt = new TH2F("phi","phi;p_{T} [GeV/c];phi [rad]",100,0,10,60,-4,4);
phi_pt->Sumw2();
phi_pt_unlike = new TH2F("phi_unlike","phi_unlike;p_{T} [GeV/c];phi [rad]",100,0,10,60,-4,4);
phi_pt_unlike->Sumw2();
phi_pt_like = new TH2F("phi_like","phi_like;p_{T} [GeV/c];phi [rad]",100,0,10,60,-4,4);
phi_pt_like->Sumw2();
dca_pt = new TH2F("dca","dca;p_{T} [GeV/c];DCA [cm]",100,0,10,20,0,4);
dca_pt->Sumw2();
dca_pt_unlike = new TH2F("dca_unlike","dca_unlike;p_{T} [GeV/c];DCA [cm]",100,0,10,20,0,4);
dca_pt_unlike->Sumw2();
dca_pt_like = new TH2F("dca_like","dca_like;p_{T} [GeV/c];DCA [cm]",100,0,10,20,0,4);
dca_pt_like->Sumw2();

nHits_nHitsRatio_pt_partner_unlike = new TH3F("nHits_nHitsRatio_pt_partner_unlike","nHits_nHitsRatio_pt_partner_unlike;p_{T} (GeV/c); nHits_partner, nHits/nMax partner",100, 0, 10,50,0,50,50,0,1);
nHits_nHitsRatio_pt_partner_unlike->Sumw2();
nHits_nHitsRatio_pt_partner_like = new TH3F("nHits_nHitsRatio_pt_partner_like","nHits_nHitsRatio_pt_partner_like;p_{T} (GeV/c); nHits_partner, nHits/nMax partner",100,0,10,50,0,50,50,0,1);
nHits_nHitsRatio_pt_partner_like->Sumw2();
pT_partner_unlike = new TH2F("pT_partner_unlike","pT_partner_unlike;p_{T} (primary) [GeV/c];p_{T} (GeV/c) (partner)",100,0,10,100,0,10);
pT_partner_unlike->Sumw2();
pT_partner_like = new TH2F("pT_partner_like","pT_partner_like;p_{T} (primary) [GeV/c]; p_{T} (partner) (GeV/c)",100,0,10,100,0,10);
pT_partner_like->Sumw2();
pair_dca_unlike = new TH2F("pair_dca_unlike","pair_dca_unlike;p_{T} [GeV/c];pairDCA",100,0,10,50,0,5);
pair_dca_unlike->Sumw2();
pair_dca_like = new TH2F("pair_dca_like","pair_dca_like;p_{T} [GeV/c];pairDCA",100,0,10,50,0,5);
pair_dca_like->Sumw2();



}
