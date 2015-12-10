/* Author Sevil Salur 2007
   Including Mark's PP study
   updated with TStarJetPico from Mateusz
   #StMuJetAnalysisTreeMaker.h  
*/

#ifndef StMuJetAnalysisTreeMaker_hh     
#define StMuJetAnalysisTreeMaker_hh
//
//  Include files
#include "StMaker.h"
#include <string>
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/filters/StEmcFilter.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StTriggerUtilities/StTriggerSimuMaker.h"
#include <TH1.h>
#include <TH2.h> 
#include <TH3.h>
#include "TObjArray.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "TVector3.h"

//  Forward declarations
class StMuTrack;
//class StMuV0I;
class TFile;
class TH1D;
//class TStMuEventAna;

class StEmcGeom;
class StEmcPosition;
class StEmcPoint;
class StSPtrVecEmcCluster;
class StPtrVecEmcCluster;
class StEmcCollection;
class StEmcFilter;
class StEmcCluster;
class StBemcTables; //v3.14
class StTriggerSimuMaker;
//class StV0MuDst;

class TStarJetPicoEvent;
//class TStMuCutV0Jet;
class TStMuCutEventJet;
class TStarJetPicoTriggerInfo;
//class TStarJetPicoQAHistograms;
class TStarJetPicoPrimaryTrack;
class TStarJetPicoGlobalTrack;
class TStarJetPicoPairs;
class StRefMultCorr;

#ifndef ST_NO_NAMESPACES
using std::string;
#endif 
//
//  The class declaration. It innherits from StMaker.
class StMuJetAnalysisTreeMaker : public StMaker {
public:

    StMuJetAnalysisTreeMaker(const Char_t *name="muAnalysis");   // constructor
    ~StMuJetAnalysisTreeMaker();                                 // destructor
    
    void     Clear(Option_t *option=""); // called after every event to cleanup
    void     setRootFile(const char* outpath,const char* filename);  
    Int_t    Init();                   // called once at the beginning of your job
    Int_t    Make();                   // invoked for every event
    Int_t    Finish();                 // called once at the end
    Int_t    doPrimTrks();
    Int_t    doGlobTrks();
    Int_t    doPhotoElectrons();
    //Int_t    doV0s();
    //Int_t    doCheckMatchedTracks(bool primary);
    //Int_t    doTowerMatching(bool primary); 
    //  Int_t    GetTowerE(int,int,int);
  
    //void     SetDoV0s(Bool_t doit){mDoV0 = doit;}
    void     SetDoPrimTrks(Bool_t doit){mDoPrimTrks = doit;}
    void     SetDoGlobTrks(Bool_t doit){mDoGlobTrks = doit;}
   
    //void SetNHits(Int_t nhits){fNHits=nhits;}
    //Int_t GetNHits()const{return fNHits;}
   
  
    //TStMuCutV0Jet*  GetV0Cuts(){return fTCutV0;};
    TStMuCutEventJet* GetEventCuts(){return  fTCutEvent;};

//    StEmcCluster* findSMDCluster(Float_t,Float_t,Int_t,Int_t&);
  //  Int_t  SMDHits(Int_t,Int_t);
    
/*    void fill_BEMC_information(TStarJetPicoPrimaryTrack *track, Int_t nTracks, StEmcCluster& assoc_bemc, StThreeVectorD& pos_bemc, StThreeVectorD& proj_bemc, Float_t emcdz_bemc, Float_t emcdphi_bemc);
    void fill_BEMC_information_global(TStarJetPicoGlobalTrack *track, Int_t nTracks, StEmcCluster& assoc_bemc, StThreeVectorD& pos_bemc, StThreeVectorD& proj_bemc, Float_t emcdz_bemc, Float_t emcdphi_bemc);
    void fill_SMDE_information(TStarJetPicoPrimaryTrack *track, Int_t nTracks, StEmcCluster& assoc_smde, StThreeVectorD& pos_smde, StThreeVectorD& proj_smde, Float_t emcdz_smde, Float_t emcdphi_smde);
    void fill_SMDE_information_global(TStarJetPicoGlobalTrack *track, Int_t nTracks, StEmcCluster& assoc_smde, StThreeVectorD& pos_smde, StThreeVectorD& proj_smde, Float_t emcdz_smde, Float_t emcdphi_smde);
    void fill_SMDP_information(TStarJetPicoPrimaryTrack *track, Int_t nTracks, StEmcCluster& assoc_smdp, StThreeVectorD& pos_smdp, StThreeVectorD& proj_smdp, Float_t emcdz_smdp, Float_t emcdphi_smdp);
    void fill_SMDP_information_global(TStarJetPicoGlobalTrack *track, Int_t nTracks, StEmcCluster& assoc_smdp, StThreeVectorD& pos_smdp, StThreeVectorD& proj_smdp, Float_t emcdz_smdp, Float_t emcdphi_smdp);
    void fill_PRS_information(TStarJetPicoPrimaryTrack *track, Int_t nTracks, StEmcCluster& assoc_prs, StThreeVectorD& pos_prs, StThreeVectorD& proj_prs, Float_t emcdz_prs, Float_t emcdphi_prs);
    void fill_PRS_information_global(TStarJetPicoGlobalTrack *track, Int_t nTracks, StEmcCluster& assoc_prs, StThreeVectorD& pos_prs, StThreeVectorD& proj_prs, Float_t emcdz_prs, Float_t emcdphi_prs);
     
    Bool_t getClusterInfo(StEmcCluster& cluster, Int_t n, Float_t& eta, Float_t& phi, Float_t& seed, Float_t& avgadd, Float_t& e);
    void getClusterEnergyProfile(StEmcCluster& cluster, Float_t& SeedE, Float_t& AvgAddE);
    Bool_t get_assoc_emc_cluster(StMuTrack* trk, StSPtrVecEmcCluster& clusters, StThreeVectorD& proj, Float_t& emcdz, Float_t& emcdphi, StThreeVectorD& pos, StEmcCluster& cluster_min, Float_t radius);
    Bool_t get_closest_assoc_clust(StSPtrVecEmcCluster& emccluster, StThreeVectorD& proj, StThreeVectorD& pos, Float_t radius, Int_t& i_min);
    void get_etower_adc(StEmcCluster& cluster, Int_t det, unsigned int *tower_adc, unsigned int *softId, Float_t *e);
    Bool_t get_proj_tower_adc_e_id(StMuTrack* track, unsigned int& tower_id, unsigned int& tower_adc, Float_t& tower_e);

*/
    bool getBEMC(StMuTrack *, int *, int *, float *, float *, int *, int *);
    void FillBEMC(TStarJetPicoPrimaryTrack *track, int *adc, float *ene, float *d, int *nep, int *towid);
void initEmc();
void finishEmc();


    void     SetVerbose(bool verb)           {mVerbose = verb;};
    Double_t GetReactionPlane();
   
   //inline  Int_t GetMatchedTracks()  {return rMatch;};

    void SetFlagData(Int_t Cut=1);

  
    Int_t GetEventCounter() {return mEventCounter;}
   Int_t GetInputEventCounter() {return mInputEventCounter;}

    
    virtual const char *GetCVS() const {
      static const char cvs[]="Tag $Name:  $ $Id: StMuJetAnalysisTreeMaker.h,v 1.1 2004/08/10 16:09:11 perev Exp $ built "__DATE__" "__TIME__ ; 
      return cvs;
    }
    
 protected:
    bool PhotoCandidate(StMuTrack *t, StMuTrack *glob_t);//(TStarJetPicoPrimaryTrack *t);//,double bemcE);
    bool GlobalPartnerCandidate(StMuTrack *t);

 private:
    //double fBEMC[5000]; //bemc energy for tracks
    void CleanBemc();
    //bool CheckForJpsiElectronCandidate();
    //Float_t computeXY(const StThreeVectorF&, const StPhysicalHelixD &);
    //int FindPhotoElectrons();   
 
    //Mateusz Class variables;
    TStarJetPicoEvent *MEvent;

    StMuDst* mCurrentMu;
    TTree *MTree;
    //TStarJetPicoQAHistograms *QAHist;
    //StRefMultCorr* fRefMultCorr;

    TString      fFilename;     //! Name of root file
    TString      fFilename_monitor;     //! Name of root file

    TString	file_name; //name of the muDst histograms file

    Int_t        mEventCounter;  //!
    Int_t        mInputEventCounter;
    TFile        *mFile;         //!
    TFile        *mFile_monitors;         //!
    
    // method (a simple track filter)
    //StEmcGeom       *mGeom;   

 /*   StEmcGeom* bemcGeom;
    StEmcGeom* smdeGeom;
    StEmcGeom* smdpGeom;
    StEmcGeom* prsGeom;*/

    StEmcGeom* mEmcGeom[4];
      
    StEmcCollection *mEmcCol;
    StEmcPosition* mEmcPosition;
    StEmcFilter     *mFilter;
    StBemcTables    *mTables;
    Int_t    mCounter;
    Float_t PrimVertexZ;
  
    
    Int_t *fPrimIndexArray;//!
    Int_t *fGlobIndexArray;//!
    //Int_t *fMatchTrArr;//!
    //Float_t *fMatchTrEtaArr;//!
    //Float_t *fMatchTrPhiArr;//!
    //Float_t *fEtaArray;//!
    //Float_t *fPhiArray;//!

    /*Int_t *idMuDstPrimIndexArray;
    Int_t *idMuDstGlobIndexArray;
    Int_t *idPrimIndexArray;
    Int_t *idGlobIndexArray;

    Int_t PrimIndex;
    Int_t GlobIndex;
*/

//commented to see if there will be Gexception error (trying to acces non existing STL, e.g. std::vector)
    std::vector<Int_t> idMuDstPrimIndex;
    std::vector<Int_t> idMuDstGlobIndex;
    std::vector<Int_t> idJetTreePrimIndex;
    std::vector<Int_t> idJetTreeGlobIndex;


    bool mVerbose;              

    Int_t   rBarrelPts; 
     
    Int_t    rMatch; //counter for matched
    Double_t Bfield; 
    Double_t mBField;           
    Float_t rplane;

    //primary track cuts    
    Int_t fNHits; //required number of hits
    Double_t  fEta;

    // Int_t nV0s;
    //Bool_t mDoV0;
    Bool_t mDoPrimTrks;
    Bool_t mDoGlobTrks;
   // Int_t fNTOTMatchedTr;
   // Int_t fMatchedTow;
    
    //this fills in events regardless of Primary Vertex existence/position: beyond cut/non-existing: only header&triggers are fille


    //TStMuCutV0Jet  *fTCutV0;        //! V0 cut class
     Int_t nFlagData;
    
    TStMuCutEventJet  *fTCutEvent;        //! event cut class

    //histogram filling functions
//    void FillHistos_beforeEvtCuts();
 //   void FillHistos_afterTrigCut();
 //   void FillHistos_afterEvtCuts();
 
    //Bool_t getClusterInfo(StEmcCluster& cluster, Float_t eta, Float_t phi, Float_t e);
    

void HistAllocation(); 

//******multiplicity histograms
//

StRefMultCorr *refmultCorrUtil;

TH1F *mult_mb;
   TH1F *mult_cent;
   TH1F *mult_cent_cut;
   TH1F *mult_cent_cut_5_10;
   TH1F *mult_cent_cut_10_20;
   TH1F *centrality;
   TH1F *mult_runid;
   TH1F *event_runid;
   TH1F *mult_runid_cent;
   TH1F *event_runid_cent;
   TH1F *mult_runid_cent_cut;
   TH1F *event_runid_cent_cut;
   
   TH1F *mult_runid_zoom;
   TH1F *mult_runid_zoom_corr;
   TH1F *event_runid_zoom;
   TH1F *mult_runid_cent_zoom;
   TH1F *mult_runid_cent_zoom_corr;
   TH1F *event_runid_cent_zoom;
   TH1F *mult_runid_cent_cut_zoom;
   TH1F *mult_runid_cent_cut_zoom_corr;
   TH1F *event_runid_cent_cut_zoom;
TH1F *multCorr;
TH1F *multCorr_cent;
TH1F *multCorr_cent_cut;
TH1F *hWeight;
//*****end of multiplicity histograms

TH1F *Vz_before;
TH1F *DiffVz_before;

TH1F *hVz;
TH1F *hDiffVz;

TH1F *refMultiplicity;
//TH1F *centrality;
//TH1F *mult_cent;

TH2F *hdEdx;
TH2F *hnsigma;

TH2F *eta_pt;
TH2F *phi_pt;
TH2F *dca_pt;
TH2F *firstTPC;
TH2F *ndEdx_pt;
TH3F *nHits_nMax_pt;

TH1F *electron_pt;
TH2F *nsigma_purity_pt;

TH2F *hDeltaZ;
TH2F *hDeltaZ_unlike;
TH2F *hDeltaZ_like;
TH2F *hDeltaPhi;
TH2F *hDeltaPhi_unlike;
TH2F *hDeltaPhi_like;
TH2F *hpE;
TH2F *hpE_unlike;
TH2F *hpE_like;
TH2F *hnSMDE;
TH2F *smde_unlike;
TH2F *smde_like;
TH2F *hnSMDP;
TH2F *smdp_unlike;
TH2F *smdp_like;

TH2F *hinvmass_pt_unlike;
TH2F *hinvmass_pt_like;
TH2F *nsigma_partner_unlike;
TH2F *nsigma_partner_like;

TH1F *photonic_pt_unlike;
TH1F *photonic_pt_like;

TH2F *nsigma_primary_unlike;
TH2F *nsigma_primary_like;

TH2F *noEMC_unlike;
TH2F *noEMC_like;
TH2F *EMC_unlike;
TH2F *EMC_like;

TH2F *eta_pt_unlike;
TH2F *phi_pt_unlike;
TH2F *dca_pt_unlike;
TH2F *eta_pt_like;
TH2F *phi_pt_like;
TH2F *dca_pt_like;

TH3F *nHits_nHitsRatio_pt_partner_like;
TH3F *nHits_nHitsRatio_pt_partner_unlike;
TH2F *pT_partner_unlike;
TH2F *pT_partner_like;
TH2F *pair_dca_unlike;
TH2F *pair_dca_like;


    ClassDef(StMuJetAnalysisTreeMaker,4)
};


inline void StMuJetAnalysisTreeMaker::SetFlagData(Int_t Cut){nFlagData=Cut;}   
#endif
