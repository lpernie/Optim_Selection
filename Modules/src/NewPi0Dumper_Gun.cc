// -*- C++ -*-
//
// Package:    NewPi0Dumper_Gun
// Class:      NewPi0Dumper_Gun
// 
/**\class NewPi0Dumper_Gun NewPi0Dumper_Gun.cc  ,Analysis/Modules/src/NewPi0Dumper_Gun.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Shahram Rahatlou
//         Created:  Wed Aug 25 10:44:55 CEST 2010
// $Id: NewPi0Dumper_Gun.cc,v 1.1 2013/02/04 16:16:13 lpernie Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <utility>
#include <cmath>
#include<algorithm>
#include<string>
#include<set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
//Geom
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "FWCore/Framework/interface/ESHandle.h"
//ES
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
//@#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "RecoEcal/EgammaClusterAlgos/interface/PreshowerClusterAlgo.h"
#include "Analysis/Modules/interface/PreshowerTools.h"
//@#include "Analysis/CalibTools/interface/EcalPreshowerHardcodedGeometry.h"
//@##include "Analysis/CalibTools/interface/EcalPreshowerHardcodedTopology.h"
//@##include "Analysis/CalibTools/interface/PreshowerCluster.h"
//@##include "Analysis/CalibTools/interface/PreshowerTools.h"
//@##include "Analysis/CalibTools/interface/GeometryService.h"
//@##include "Analysis/CalibTools/interface/EndcapTools.h"
//@#//EB Cont Correction
//@##include "Analysis/CalibTools/interface/EcalEnerCorr.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//PF
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapHardcodedTopology.h"
//@#include "Analysis/Pi0Calib/interface/ECALGeometry.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
//@#include "Analysis/Pi0Calib/interface/EcalCalibMap.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// Trigger
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"

// conversions
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

//Laser and transparency loss
#include "CalibCalorimetry/EcalLaserAnalyzer/interface/MEEBGeom.h"
#include "CalibCalorimetry/EcalLaserAnalyzer/interface/MEEEGeom.h"
#include "DataFormats/Provenance/interface/Timestamp.h"

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"

#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
//HLT
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include <DataFormats/Common/interface/TriggerResults.h>

//Json
#include "Analysis/Modules/interface/JSON.h"

// root
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <TH1F.h>
#include <TH2F.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRegexp.h"

using std::vector;
using std::cout;
using std::endl;
using std::pair;
using std::max_element;
using std::set;

//#define INDEXETA65

//Function
double max_array(double *A, int n);
double max(double x, double y);
void convxtalid(Int_t &nphi,Int_t &neta);
int diff_neta_s(Int_t neta1, Int_t neta2);
int diff_nphi_s(Int_t nphi1,Int_t nphi2);

//
// class declaration
//

struct TrackIsoVars {
    float   ptIso015;
    int   ntrkIso015;  
    float   ptIso035;
    int   ntrkIso035;  
    float   ptIso040;
    int   ntrkIso040;  
    TrackIsoVars() :
	  ptIso015(0.), ntrkIso015(0),
	  ptIso035(0.), ntrkIso035(0),
	  ptIso040(0.), ntrkIso040(0){ }
};

struct HCALIsoVars {
    float hcalIso005;
    float hcalIso010;
    float hcalIso040;
    HCALIsoVars() : hcalIso005(0.), hcalIso010(0.), hcalIso040(0.) { }
};

struct ClusterShape {
    float s1;
    float s4;
    float s9;
    float s25;
    float s4s9;
    float time;
    float ncry9;
    float ncry25;
    int flag;

    // major and minor cluster moments wrt principale axes:
    float sMaj9;
    float sMin9;
    float sMaj25;
    float sMin25;
    float sMaj49;
    float sMin49;

    float transparencyLoss;
    float laserCorr;
    float alpha;

    ClusterShape() : s1(-99.), s4(-99.), s9(-99.), s25(-99.), s4s9(-1.), 
    time(-9999.), ncry9(-1), ncry25(-1), flag(0.) ,
    sMaj9(-99.), sMin9(-99.), 
    sMaj25(-99.), sMin25(-99.),
    sMaj49(-99.), sMin49(-99.), 
    transparencyLoss(-99.),
    laserCorr(-99.), alpha(-99)
    { }
};

//Other class
class ecalRecHitPtrLess : public std::binary_function<EcalRecHit*, EcalRecHit*, bool>
{
    public:
	  bool operator()(EcalRecHit* x, EcalRecHit* y)
	  {
		return (x->energy() > y->energy());
	  }
};

class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
    public:
	  bool operator()(EcalRecHit x, EcalRecHit y)
	  {
		return (x.energy() > y.energy());
	  }
};

#define NCLUMAX 30000
#define NPI0MAX 15000

using namespace reco;

class NewPi0Dumper_Gun : public edm::EDAnalyzer {
    public:
	  explicit NewPi0Dumper_Gun(const edm::ParameterSet&);
	  ~NewPi0Dumper_Gun();


    private:
	  virtual void beginJob() ;
	  virtual void analyze(const edm::Event&, const edm::EventSetup&);
	  virtual void endJob() ;
          void FillTriggerInfo(const edm::Event&, const edm::EventSetup&);
	  template <class DetIdType, class RecHitCollectionType, class IteratorType> 
		void make3x3Clusters( const CaloGeometry* geometry, const RecHitCollectionType* hits, std::vector<CaloCluster>* clusters,
			  std::vector< ClusterShape >* shapes, edm::Timestamp const & iTime, vector<int>* Ncristal,vector<float>* S4S9,vector<float>* S1S9,vector<float>* S2S9, vector<float>* EtrSeed );
	  void fillPi0Tree(
		    math::PtEtaPhiMLorentzVector &pi0P4PV,  std::map<size_t,size_t> &savedCluEB,
		    math::PtEtaPhiMLorentzVector &g1P4PV, math::PtEtaPhiMLorentzVector &g2P4PV,
		    std::vector< ClusterShape > &shapes,
		    const CaloCluster* g1, const CaloCluster* g2, int &nClu, int i, int j, vector<int> &Ncristal, vector<float>& EtSeeds);
	  bool GetHLTResults(const edm::Event& iEvent, std::string s);

	  double DPhi(double phi1, double phi2);
	  double DR( double phi1, double phi2, double eta1, double eta2);
	  double min( double a, double b);

	  // ----------member data ---------------------------
	  CaloTopology *ebTopology;         // hardcoded topology
	  CaloTopology *eeTopology;         // hardcoded topology
	  EcalPreshowerTopology *estopology_;
	  const EcalPreshowerGeometry *esGeometry_;

	  //Json
	  JSON* myjson;

	  //General values
	  int             runn;
	  int             eventn;
	  int             ls;

	  // primary vertex
	  float chi2PV, ndofPV, ntrkPV, xPV, yPV, zPV;

	  //Paricles
	  std::vector<TLorentzVector> MyEBParicles;
	  std::vector<TLorentzVector> MyEEParicles;

	  //InputTag
	  edm::InputTag EBRecHitCollectionTag_;
	  edm::InputTag EERecHitCollectionTag_;
	  edm::InputTag ESRecHitCollectionTag_;
        edm::InputTag l1TriggerInput_;
	  edm::InputTag triggerTag_;

	  //Geometry value
	  bool  param_LogWeighted_ ;
	  float param_T0_barl_;
	  float param_T0_endc_;    
	  float param_T0_endcES_;
	  float param_W0_;
	  float param_X0_;         

	  // cuts
	  Float_t ptpi0Cut_, masspi0Cut_, ptCluCut_, s1CluCut_, s1CluCutEE_, s4s9CluCut_;
	  bool useES_E_;  

	  // trigger and conds foir GOODCOLL in MC
	  bool isBSC;
	  bool isGoodPrimaryVertex;
	  bool notScraping;
	  bool isHFMinBias;
	  bool isGOODCOLL;

	  //TFile
	  TFile* m_file;
	  std::string outfilename_;

          static const int MAXL1bits = 200;
          int nL1bits;
          int L1bits[MAXL1bits];
          int nL1bitsTech;
          int L1bitsTech[MAXL1bits];

	  //My Pi0 Branch
	  TTree* Tree_HLT; 
	  Int_t nPi0;
	  Int_t   STr2_NPi0_rec;
	  Int_t   STr2_Pi0recIsEB[NPI0MAX];
	  Float_t STr2_IsoPi0_rec[NPI0MAX];
	  Int_t   STr2_n1CrisPi0_rec[NPI0MAX];
	  Int_t   STr2_n2CrisPi0_rec[NPI0MAX];
	  Float_t STr2_mPi0_rec[NPI0MAX];
	  Float_t STr2_ptG1_rec[NPI0MAX];
	  Float_t STr2_ptG2_rec[NPI0MAX];
	  Float_t STr2_etaPi0_rec[NPI0MAX];
	  Float_t STr2_ptPi0_rec[NPI0MAX];
	  Float_t STr2_DeltaRG1G2[NPI0MAX];
	  Float_t STr2_Es_e1_1[NPI0MAX];
	  Float_t STr2_Es_e1_2[NPI0MAX];
	  Float_t STr2_Es_e2_1[NPI0MAX];
	  Float_t STr2_Es_e2_2[NPI0MAX];
	  Float_t STr2_S4S9_1[NPI0MAX];
	  Float_t STr2_S4S9_2[NPI0MAX];

	  TTree* Tree_HLT_clus;
	  Float_t ePi0[NPI0MAX];
	  Float_t massPi0[NPI0MAX];
	  Float_t ptPi0[NPI0MAX];
	  Float_t etaPi0[NPI0MAX];
	  Float_t phiPi0[NPI0MAX];
	  Float_t CristTot[NPI0MAX];
	  Float_t DRClus[NPI0MAX];
	  Float_t isoPi0[NPI0MAX];
	  Int_t ietaTTPi0[NPI0MAX];
	  Int_t iphiTTPi0[NPI0MAX];
	  Int_t indexClu1Pi0[NPI0MAX];
	  Int_t indexClu2Pi0[NPI0MAX];

	  //My Cluster Branch
	  Int_t nClu;
	  Int_t nCris[NCLUMAX];
	  Float_t EtSeed[NCLUMAX];
	  Float_t S1Clu[NCLUMAX];
	  Float_t S4Clu[NCLUMAX];
	  Float_t S9Clu[NCLUMAX];
	  Float_t S25Clu[NCLUMAX];
	  Float_t etaClu[NCLUMAX];
	  Float_t phiClu[NCLUMAX];
	  Float_t ptClu[NCLUMAX];
	  Float_t timeClu[NCLUMAX];
	  Int_t nCryClu[NCLUMAX];
	  Int_t indexCryClu[NCLUMAX][9];

	  Int_t flagClu[NCLUMAX];
	  Int_t ietaClu[NCLUMAX];
	  Int_t iphiClu[NCLUMAX];
	  Int_t iCryClu[NCLUMAX];
	  Int_t iSMClu[NCLUMAX];
	  Int_t imodClu[NCLUMAX];
	  Int_t iTTClu[NCLUMAX];
	  Int_t iTTetaClu[NCLUMAX];
	  Int_t iTTphiClu[NCLUMAX];

	  float sMaj9Clu[NCLUMAX];
	  float sMin9Clu[NCLUMAX];
	  float sMaj25Clu[NCLUMAX];
	  float sMin25Clu[NCLUMAX];
	  float sMaj49Clu[NCLUMAX];
	  float sMin49Clu[NCLUMAX];



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NewPi0Dumper_Gun::NewPi0Dumper_Gun(const edm::ParameterSet& ps)
{

    cout<<"Constructor"<<endl;

    //now do what ever initialization is needed

    //Json
    myjson = new JSON( "/afs/cern.ch/work/l/lpernie/pi0/Calibration/CMSSW_5_3_2_patch2/src/CalibCode/submit/common/goodrunlist_json2012.txt" );
    // collections to fetch
    useES_E_               = ps.getUntrackedParameter<bool>("useES_E",false);
    EBRecHitCollectionTag_= ps.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag");
    EERecHitCollectionTag_= ps.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag");
    ESRecHitCollectionTag_= ps.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag");
    triggerTag_           = ps.getUntrackedParameter<edm::InputTag>("triggerTag",edm::InputTag("TriggerResults"));

    param_LogWeighted_ = true;
    param_T0_barl_     = 5.7;
    param_T0_endc_     = 3.1;
    param_T0_endcES_   = 1.2;
    param_W0_          = 4.2;
    param_X0_          = 0.89;

    //cuts
    ptpi0Cut_             = ps.getUntrackedParameter<double>("ptpi0Cut", 0.700);
    masspi0Cut_           = ps.getUntrackedParameter<double>("masspi0Cut", 0.35);
    s1CluCut_             = ps.getUntrackedParameter<double>("s1CluCut", 0.350);
    s1CluCutEE_           = ps.getUntrackedParameter<double>("s1CluCutEE", 0.3);
    ptCluCut_             = ps.getUntrackedParameter<double>("ptCluCut", 0.3);
    s4s9CluCut_           = ps.getUntrackedParameter<double>("s4s9CluCut", 0.85);
    cout<<"Cut Used: ptpi0Cut: "<<ptpi0Cut_<<" s1CluCut "<<s1CluCut_<<" s1CluCutEE "<<s1CluCutEE_<<" ptCluCut: "<<ptCluCut_<<" s4s9CluCut: "<<s4s9CluCut_<<endl;
    //File
    outfilename_          = ps.getUntrackedParameter<std::string>("OutputFile","NewPi0Tuple.root");


    //Topology
    ebTopology = new CaloTopology();  
    eeTopology = new CaloTopology();  
    EcalBarrelHardcodedTopology* ebHCTopology=new EcalBarrelHardcodedTopology();
    EcalEndcapHardcodedTopology* eeHCTopology=new EcalEndcapHardcodedTopology();
    ebTopology->setSubdetTopology(DetId::Ecal,EcalBarrel,ebHCTopology);
    eeTopology->setSubdetTopology(DetId::Ecal,EcalEndcap,eeHCTopology);

    l1TriggerInput_ = ps.getUntrackedParameter<edm::InputTag>("L1TriggerTag", edm::InputTag("gtDigis"));

    //estopology_ = new EcalPreshowerHardcodedTopology();

    /*
    //@@endcapRawSuperClusterCollection_ = ps.getUntrackedParameter<edm::InputTag>("endcapRawSuperClusterCollection");
    conversionsTag_       = ps.getUntrackedParameter<edm::InputTag>("conversionsTag");
    outfilename_ContCorr_ = ps.getUntrackedParameter<std::string>("outfilename_ContCorr","ContCorr.root");
    externalGeometry_     = ps.getUntrackedParameter<std::string>("ExternalGeometry");
    storeMCTruth_         = ps.getUntrackedParameter<bool>("StoreMCTruth",false);
    storeConversions_     = ps.getUntrackedParameter<bool>("StoreConversions",true);
    useHCAL_ = ps.getUntrackedParameter<bool>("useHCAL",false);
    goodCollSelection_ = ps.getUntrackedParameter<bool>("goodCollSelection",false);
    useBeamSpotPosition_ = ps.getUntrackedParameter<bool>("useBeamSpotPosition",true);
    invertEELaserCorrections_ = ps.getUntrackedParameter<bool>("InvertEELaserCorrections",true);
    PVTag_                = ps.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag",
    edm::InputTag("hiSelectedVertex"));
    TracksTag_            = ps.getUntrackedParameter<edm::InputTag>("TracksTag",
    edm::InputTag("hiGlobalPrimTracks"));
    HBHETag_              = ps.getUntrackedParameter<edm::InputTag>("HBHETag",
    edm::InputTag("hbhereco"));
    HFTag_      = ps.getUntrackedParameter<edm::InputTag>( "HFTag",edm::InputTag("hfreco") );
    PFJetsTag_  = ps.getUntrackedParameter<edm::InputTag>("PFJetsTag", edm::InputTag("ak5PFJets"));
    PFMetTag_   = ps.getUntrackedParameter<edm::InputTag>("PFMetTag", edm::InputTag("pfMet"));
    DoOffGeom_         = ps.getUntrackedParameter<bool>("DoOffGeom",true);
    if(DoOffGeom_) cout<<"Using: Official geometry"<<endl;
    else           cout<<"Using: external geometry"<<endl;
    //@l1TriggerInput_ = ps.getUntrackedParameter<edm::InputTag>("L1TriggerTag", edm::InputTag("hltGtDigis"));
    //CaloMetTag_  = ps.getUntrackedParameter<edm::InputTag>("CaloMetTag", edm::InputTag("met"));

    //pfakt5JetCorrectionServiceTag_  
    //= ps.getUntrackedParameter<std::string>("pfakt5JetCorrectionService","ak5PFL2L3");

    mAPDPNRatiosRef_  = 0;
    mAPDPNRatios_ = 0;
    mAlphas_ = 0;

    //#TFile* f = TFile::Open("caloGeometry.root");
    TFile* f = TFile::Open( "/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/CalibCode/submit/common/caloGeometry.root"  );
    //TString path(pathExtGeo_);
    //TFile* f = TFile::Open( path.Data() );
    geom = ECALGeometry::getGeometry(f);
    //ES
    externalGeometryFile_ = TFile::Open(externalGeometry_.c_str());
    if(!externalGeometryFile_) cms::Exception("ExtGeom") << "External Geometry file (" << externalGeometry_ << ") not found" << endl;
    geom_ = ECALGeometry::getGeometry(externalGeometryFile_);
    GeometryService::setGeometryName(externalGeometry_);
    GeometryService::setGeometryPtr(geom_);

    hardcodedPreshowerGeometry_ = new EcalPreshowerHardcodedGeometry(geom_);
    hardcodedPreshowerGeometry_->initializeParms();
     */
}


NewPi0Dumper_Gun::~NewPi0Dumper_Gun()
{
    cout<<"Destructor"<<endl;
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    m_file->Write();
    m_file->Close();
    delete Tree_HLT;
    delete Tree_HLT_clus;
    delete myjson;

}

// ------------ method called once each job just before starting event loop  ------------
void NewPi0Dumper_Gun::beginJob()
{
    cout<<"Begin Job"<<endl;

    m_file = new TFile(outfilename_.c_str(),"RECREATE");
    if(!m_file) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << outfilename_.c_str() << "\n";
    m_file->cd();

    Tree_HLT = new TTree("Tree_HLT","Output TTree");
    Tree_HLT->Branch( "STr2_NPi0_rec",      &STr2_NPi0_rec,    "STr2_NPi0_rec/I");
    Tree_HLT->Branch( "STr2_Pi0recIsEB",    &STr2_Pi0recIsEB,    "STr2_Pi0recIsEB[STr2_NPi0_rec]/I");
    Tree_HLT->Branch( "STr2_IsoPi0_rec",    &STr2_IsoPi0_rec,     "STr2_IsoPi0_rec[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_n1CrisPi0_rec", &STr2_n1CrisPi0_rec,  "STr2_n1CrisPi0_rec[STr2_NPi0_rec]/I");
    Tree_HLT->Branch( "STr2_n2CrisPi0_rec", &STr2_n2CrisPi0_rec,  "STr2_n2CrisPi0_rec[STr2_NPi0_rec]/I");
    Tree_HLT->Branch( "STr2_mPi0_rec",      &STr2_mPi0_rec,       "STr2_mPi0_rec[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_ptG1_rec",      &STr2_ptG1_rec,       "STr2_ptG1_rec[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_ptG2_rec",      &STr2_ptG2_rec,       "STr2_ptG2_rec[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_etaPi0_rec",    &STr2_etaPi0_rec,     "STr2_etaPi0_rec[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_ptPi0_rec",     &STr2_ptPi0_rec,      "STr2_ptPi0_rec[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_DeltaRG1G2",    &STr2_DeltaRG1G2,     "STr2_DeltaRG1G2[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_Es_e1_1",       &STr2_Es_e1_1,        "STr2_Es_e1_1[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_Es_e1_2",       &STr2_Es_e1_2,        "STr2_Es_e1_2[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_Es_e2_1",       &STr2_Es_e2_1,        "STr2_Es_e2_1[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_Es_e2_2",       &STr2_Es_e2_2,        "STr2_Es_e2_2[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_S4S9_1",        &STr2_S4S9_1,         "STr2_S4S9_1[STr2_NPi0_rec]/F");
    Tree_HLT->Branch( "STr2_S4S9_2",        &STr2_S4S9_2,         "STr2_S4S9_2[STr2_NPi0_rec]/F");

    /// trigger variables
    Tree_HLT->Branch("nL1bits",&nL1bits,"nL1bits/I");
    Tree_HLT->Branch("L1bits",L1bits,"L1bits[nL1bits]/I");
    Tree_HLT->Branch("nL1bitsTech",&nL1bitsTech,"nL1bitsTech/I");
    Tree_HLT->Branch("L1bitsTech",L1bitsTech,"L1bitsTech[nL1bitsTech]/I");

    Tree_HLT_clus = new TTree("Tree_HLT_clus","Output TTree Cluster");
    Tree_HLT_clus->Branch("nClu",&nClu,"nClu/I");
    Tree_HLT_clus->Branch("nCris",&nCris,"nCris[nClu]/I");
    Tree_HLT_clus->Branch("EtSeed",&EtSeed,"EtSeed[nClu]/F");
    Tree_HLT_clus->Branch("ptClu",&ptClu,"ptClu[nClu]/F");
    Tree_HLT_clus->Branch("etaClu",&etaClu,"etaClu[nClu]/F");
    Tree_HLT_clus->Branch("phiClu",&phiClu,"phiClu[nClu]/F");
    Tree_HLT_clus->Branch("S1Clu",&S1Clu,"S1Clu[nClu]/F");
    Tree_HLT_clus->Branch("S4Clu",&S4Clu,"S4Clu[nClu]/F");
    Tree_HLT_clus->Branch("S9Clu",&S9Clu,"S9Clu[nClu]/F");
    Tree_HLT_clus->Branch("S25Clu",&S25Clu,"S25Clu[nClu]/F");
    Tree_HLT_clus->Branch("timeClu",&timeClu,"timeClu[nClu]/F");
    Tree_HLT_clus->Branch("nCryClu",&nCryClu,"nCryClu[nClu]/I");
    Tree_HLT_clus->Branch("flagClu",&flagClu,"flagClu[nClu]/I");

    Tree_HLT_clus->Branch("ietaClu",&ietaClu,"ietaClu[nClu]/I");
    Tree_HLT_clus->Branch("iphiClu",&iphiClu,"iphiClu[nClu]/I");
    Tree_HLT_clus->Branch("iCryClu",&iCryClu,"iCryClu[nClu]/I");
    Tree_HLT_clus->Branch("iSMClu",&iSMClu,"iSMClu[nClu]/I");
    Tree_HLT_clus->Branch("imodClu",&imodClu,"imodClu[nClu]/I");
    Tree_HLT_clus->Branch("iTTClu",&iTTClu,"iTTClu[nClu]/I");
    Tree_HLT_clus->Branch("iTTetaClu",&iTTetaClu,"iTTetaClu[nClu]/I");
    Tree_HLT_clus->Branch("iTTphiClu",&iTTphiClu,"iTTphiClu[nClu]/I");

    Tree_HLT_clus->Branch("sMaj9Clu", &sMaj9Clu, "sMaj9Clu[nClu]/F");
    Tree_HLT_clus->Branch("sMaj25Clu",&sMaj25Clu,"sMaj25Clu[nClu]/F");
    Tree_HLT_clus->Branch("sMaj49Clu",&sMaj49Clu,"sMaj49Clu[nClu]/F");
    Tree_HLT_clus->Branch("sMin9Clu", &sMin9Clu, "sMin9Clu[nClu]/F");
    Tree_HLT_clus->Branch("sMin25Clu",&sMin25Clu,"sMin25Clu[nClu]/F");
    Tree_HLT_clus->Branch("sMin49Clu",&sMin49Clu,"sMin49Clu[nClu]/F");

    Tree_HLT_clus->Branch("nPi0",&nPi0,"nPi0/I");
    Tree_HLT_clus->Branch("massPi0",&massPi0,"massPi0[nPi0]/F");
    Tree_HLT_clus->Branch("ePi0",&ePi0,"ePi0[nPi0]/F");
    Tree_HLT_clus->Branch("ptPi0",&ptPi0,"ptPi0[nPi0]/F");
    Tree_HLT_clus->Branch("etaPi0",&etaPi0,"etaPi0[nPi0]/F");
    Tree_HLT_clus->Branch("phiPi0",&phiPi0,"phiPi0[nPi0]/F");
    Tree_HLT_clus->Branch("CristTot",&CristTot,"CristTot[nPi0]/F");
    Tree_HLT_clus->Branch("DRClus",&DRClus,"DRClus[nPi0]/F");
    Tree_HLT_clus->Branch("indexClu1Pi0",&indexClu1Pi0,"indexClu1Pi0[nPi0]/I");
    Tree_HLT_clus->Branch("indexClu2Pi0",&indexClu2Pi0,"indexClu2Pi0[nPi0]/I");

}

// ------------ method called once each job just after ending the event loop  ------------
    void 
NewPi0Dumper_Gun::endJob()
{
    cout<<"End Job"<<endl;

    m_file->cd();

    Tree_HLT->Write();
    Tree_HLT_clus->Write();
}

// ------------ method called to for each event  ------------
    void
NewPi0Dumper_Gun::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    if ( !myjson->isGoodLS( iEvent.id().run() , iEvent.luminosityBlock() ) ) return;

    MyEBParicles.clear();
    MyEEParicles.clear();

    using namespace edm;
    using namespace reco;

    // tree variables
    runn = iEvent.id().run();
    eventn = iEvent.id().event();
    ls = iEvent.luminosityBlock();

    FillTriggerInfo(iEvent, iSetup);

    chi2PV = -1.; ndofPV= -1.; ntrkPV = -1.; ndofPV = -1.;
    xPV = 0.; yPV = 0.; zPV = 0.;
    TVector3 posPV(0.,0.,0.);
    isGoodPrimaryVertex = true;
    notScraping = true;

    //EB cluster
    Handle< EBRecHitCollection > ebHandle;
    iEvent.getByLabel ( EBRecHitCollectionTag_, ebHandle);
    const EBRecHitCollection* hits = ebHandle.product();

    // get calo geometry
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloGeometry* geometry = geoHandle.product();
    estopology_ = new EcalPreshowerTopology(geoHandle);
    const CaloSubdetectorGeometry *geometry_p = geometry->getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;
    esGeometry_ = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p));

    // 3x3 clusters
    std::vector< CaloCluster > clusters; // contains the output clusters
    std::vector< ClusterShape > shapes; // contains the output clusters
    clusters.clear();
    shapes.clear();
    std::vector<int> Ncristal_EB;
    std::vector<float> s4s9_EB; std::vector<float> s1s9_EB; std::vector<float> s2s9_EB;
    std::vector<float> EtSeed_EB;
    Ncristal_EB.clear(); s4s9_EB.clear(); s1s9_EB.clear(); s2s9_EB.clear();
    EtSeed_EB.clear();
 
    if( GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*") )
       make3x3Clusters<EBDetId,EBRecHitCollection,EBRecHitCollection::const_iterator>(geometry, hits, &clusters, &shapes, iEvent.time(), &Ncristal_EB, &s4s9_EB, &s1s9_EB, &s2s9_EB, &EtSeed_EB);
    
    //for(unsigned int d=0; d<Nclus_ieta_EB.size(); d++) cout<< Nclus_ieta_EB.size()-EtSeed_EB.size()<<"  "<<Nclus_ieta_EB[d]<<endl;
    //      std::cout << "+++ Clusters done" << std::endl;

    //cout << endl << "--- EB::clusters.size() = " << clusters.size() << endl;

    // keep track of clusters used to make pi0 candidates
    std::map<size_t,size_t> savedCluEB;
    nClu = 0;

    // loop over clusters to make Pi0
    nPi0 = 0;

    for(size_t i=0; i< clusters.size(); ++i) {

	  const CaloCluster* g1 = &(clusters[i]);

	  // track association
	  //float dclustrk1(9999.);
	  //const reco::Track* trk1 =0;

	  // 2nd loop over clusters
	  for(size_t j=i+1; j<clusters.size(); ++j) 
	  {
		// only store up to max number pi0 or clusters
		if(nPi0>= NPI0MAX/2 || nClu >= NCLUMAX/2) 
		{
		  cout << "No more space in the arrays: nPi0: " << nPi0
		       << "  nClu: " << nClu << endl;
		    continue;
		}

		const CaloCluster* g2 = &(clusters[j]);           
		// using the origin
		math::PtEtaPhiMLorentzVector g1P4( g1->energy()/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
		math::PtEtaPhiMLorentzVector g2P4( g2->energy()/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
		math::PtEtaPhiMLorentzVector pi0P4 = g1P4 + g2P4;

		// use PV 
		TVector3 tmpV( g1->x()-xPV, g1->y()-yPV, g1->z()-zPV);
		math::PtEtaPhiMLorentzVector g1P4PV( g1->energy()/cosh(tmpV.Eta()), 
			  tmpV.Eta(), tmpV.Phi(), 0. );

		tmpV = TVector3( g2->x()-xPV, g2->y()-yPV, g2->z()-zPV);
		math::PtEtaPhiMLorentzVector g2P4PV( g2->energy()/cosh(tmpV.Eta()), 
			  tmpV.Eta(), tmpV.Phi(), 0. );
		math::PtEtaPhiMLorentzVector pi0P4PV = g1P4PV + g2P4PV;

		// basic selection to reduce combinatorics
		if( pi0P4PV.mass() > masspi0Cut_) continue;
		if(pi0P4PV.pt() < ptpi0Cut_) continue;

		fillPi0Tree(
			  pi0P4PV, savedCluEB,
			  g1P4PV, g2P4PV,
			  shapes,
			  g1, g2, nClu, i, j, Ncristal_EB, EtSeed_EB);

		//new cut
		float nextClu = 999., Drtmp = 999.;
		for(size_t ind=0; ind<clusters.size(); ++ind){
		    const CaloCluster* Gtmp = &(clusters[ind]);
		    double deltaR1 = DR(Gtmp->phi(),g1P4.phi(),Gtmp->eta(),g1P4.eta());
		    double deltaR2 = DR(Gtmp->phi(),g2P4.phi(),Gtmp->eta(),g2P4.eta());
		    if( ind!=i && ind!=j && (deltaR1<Drtmp || deltaR2<Drtmp ) ){
			  nextClu = min(deltaR1,deltaR2);
			  Drtmp = nextClu;
		    }
		}
		if( clusters.size()<3 ) nextClu=999.;
		STr2_IsoPi0_rec[nPi0] = nextClu;
		if( g1->energy()>g2->energy() ){
		    STr2_n1CrisPi0_rec[nPi0] = Ncristal_EB[i];
		    STr2_n2CrisPi0_rec[nPi0] = Ncristal_EB[j];
		}
		else{
		    STr2_n1CrisPi0_rec[nPi0] = Ncristal_EB[j];
		    STr2_n2CrisPi0_rec[nPi0] = Ncristal_EB[i];
		}

		STr2_mPi0_rec[nPi0] = pi0P4.mass();
		STr2_ptG1_rec[nPi0] = g1P4.pt();
		STr2_ptG2_rec[nPi0] = g2P4.pt();
		STr2_etaPi0_rec[nPi0] = pi0P4.Eta();
		STr2_ptPi0_rec[nPi0] = pi0P4.Pt();
		STr2_Pi0recIsEB[nPi0] = 1;
		STr2_DeltaRG1G2[nPi0] = DR(g1P4PV.phi(),g2P4PV.phi(),g1P4PV.eta(),g2P4PV.eta() );
		STr2_Es_e1_1[nPi0] = -1.;
		STr2_Es_e1_2[nPi0] = -1.;
		STr2_Es_e2_1[nPi0] = -1.;
		STr2_Es_e2_2[nPi0] = -1.;
		STr2_S4S9_1[nPi0]  = s4s9_EB[i];
		STr2_S4S9_2[nPi0]  = s4s9_EB[j];

		nPi0++;
		TLorentzVector thispar; thispar.SetPtEtaPhiE(pi0P4.pt(),pi0P4.eta(),pi0P4.phi(),pi0P4.energy());
		//MyEBParicles.push_back(thispar);

	  } // loop over clusters (g2)

    } // loop over clusters to make pi0 

    // endcap 3x3 clusters

    // EE rechits
    Handle< EERecHitCollection > eeHandle;
    bool goodhandle = iEvent.getByLabel ( EERecHitCollectionTag_, eeHandle);
    if(!goodhandle) throw cms::Exception("Handle") << "Bad EE Handle\n";
    const EERecHitCollection* hitsEE = eeHandle.product();

    std::vector< CaloCluster > clustersEE; // contains the output clusters
    std::vector< ClusterShape > shapesEE; // contains the output clusters
    vector<int> Ncristal_EE;
    vector<float> s4s9_EE; vector<float> s1s9_EE; vector<float> s2s9_EE;
    vector<float> EtSeed_EE;
    Ncristal_EE.clear(); clustersEE.clear();
    EtSeed_EE.clear(); s4s9_EE.clear(); s1s9_EE.clear(); s2s9_EE.clear();

if( GetHLTResults(iEvent, "AlCa_EcalPi0EEonly.*")  )
    make3x3Clusters<EEDetId,EERecHitCollection,EERecHitCollection::const_iterator>(geometry, hitsEE, &clustersEE, &shapesEE, iEvent.time(), &Ncristal_EE, &s4s9_EE, &s1s9_EE, &s2s9_EE, &EtSeed_EE);

    //loop over eecluster to find matches with preshower
    std::vector< CaloCluster > eseeclusters_tot;
    std::vector< float > Ncristal_EE_TOT;
    std::vector< float > S4S9_TOT;
    std::vector< CaloCluster > eseeclusters;
    std::vector< float > Es_e1, Es_e2;
    std::vector< float > Es_e1_tot, Es_e2_tot;
    eseeclusters.clear(); eseeclusters_tot.clear(); Ncristal_EE_TOT.clear();  S4S9_TOT.clear();
    Es_e1.clear(); Es_e2.clear(); Es_e1_tot.clear(); Es_e2_tot.clear();

    bool Added_One=false;
    int i=0;

if( GetHLTResults(iEvent, "AlCa_EcalPi0EEonly.*") ){
    edm::Handle< ESRecHitCollection > esHandle;
    goodhandle = iEvent.getByLabel ( ESRecHitCollectionTag_, esHandle);
    if(!goodhandle) throw cms::Exception("Handle") << "Bad ES Handle\n";

    PreshowerTools esClusteringAlgo(geometry, estopology_, esHandle);


    for( std::vector<CaloCluster>::const_iterator eeclus_iter  = clustersEE.begin(); eeclus_iter != clustersEE.end(); ++eeclus_iter, i++ )
    {

	  Added_One=false;
	  double X = eeclus_iter->x();       double Y = eeclus_iter->y();       double Z = eeclus_iter->z();
	  const GlobalPoint point(X,Y,Z);
       
       DetId tmp1 = esGeometry_->getClosestCellInPlane(point,1);
       DetId tmp2 = esGeometry_->getClosestCellInPlane(point,2);

//cout<<tmp1.rawId()<<"   ID  "<<tmp2.rawId()<<endl;
	  if ((tmp1.rawId()!=0) && (tmp2.rawId()!=0) && fabs(point.eta())>1.66 && fabs(point.eta()<2.56 ) )//@ fiducial cut
	  {
		ESDetId tmp1_conversion (tmp1);
		ESDetId tmp2_conversion (tmp2);

		PreshowerCluster preshowerclusterp1 = esClusteringAlgo.makeOnePreshowerCluster( PreshowerTools::clusterwindowsize_, &tmp1_conversion);
		PreshowerCluster preshowerclusterp2 = esClusteringAlgo.makeOnePreshowerCluster( PreshowerTools::clusterwindowsize_, &tmp2_conversion);

		double e1 = preshowerclusterp1.energy();
		double e2 = preshowerclusterp2.energy();

//cout<<e1<<"   ES1  "<<e2<<" "<<PreshowerTools::clusterwindowsize_<<endl;
		// GeV to #MIPs
		e1 = e1 / PreshowerTools::mip_;
		e2 = e2 / PreshowerTools::mip_;

		double tempenergy = eeclus_iter->energy();

//cout<<e1<<"   ES2  "<<e2<<endl;
		if(e1 > 2.0 || e2 > 2.0) /// cut @ 2 MIPs as suggested by Ming @ DPG/EGamma Joint Meeting 19.03.2012 
		{
		    double deltaE = PreshowerTools::gamma_*(PreshowerTools::calib_planeX_*e1 + PreshowerTools::calib_planeY_*e2);

		    tempenergy = deltaE + eeclus_iter->energy();
		    eseeclusters.push_back( CaloCluster( tempenergy, eeclus_iter->position(), CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
		    Es_e1.push_back(PreshowerTools::gamma_*(PreshowerTools::calib_planeX_*e1)); Es_e2.push_back(PreshowerTools::gamma_*(PreshowerTools::calib_planeY_*e2));
		    Added_One=true;

		    double DZ2 = (preshowerclusterp2.z()-preshowerclusterp1.z() )/2.;
		    GlobalPoint posClu(preshowerclusterp1.x()*(1.+DZ2/preshowerclusterp1.z() ),preshowerclusterp2.y()*(1.-DZ2/preshowerclusterp2.z()),(preshowerclusterp1.z()+preshowerclusterp2.z() )/2. );
		    //The perfect cluster
		    if( fabs(preshowerclusterp1.z())>30  && fabs(preshowerclusterp2.z())>30 ){
			  math::XYZPoint posit(posClu.x(),posClu.y(),posClu.z());
			  eseeclusters_tot.push_back( CaloCluster( tempenergy, posit, CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
			  Ncristal_EE_TOT.push_back( Ncristal_EE[i] );
			  S4S9_TOT.push_back( s4s9_EE[i] );
			  Es_e1_tot.push_back(PreshowerTools::gamma_*(PreshowerTools::calib_planeX_*e1)); Es_e2_tot.push_back(PreshowerTools::gamma_*(PreshowerTools::calib_planeY_*e2)); 
		    }
		}

	  }
	  if( !Added_One ){
		eseeclusters.push_back( CaloCluster( eeclus_iter->energy(), eeclus_iter->position(), CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
		Es_e1.push_back(-1.); Es_e2.push_back(-1.);
		Added_One=false;
	  }

    }//end of the ES matching loop

}//HLT

    if(Ncristal_EE_TOT.size() != eseeclusters_tot.size() || Es_e1_tot.size()!=Ncristal_EE_TOT.size()) cout<<"WARNING Ncristal_EE_TOT.size() != eseeclusters_tot.size()"<<endl;
    if(Ncristal_EE.size() != eseeclusters.size() || Es_e1.size()!=Ncristal_EE.size()  ) cout<<"WARNING Ncristal_EE.size() != eseeclusters.size()"<<endl;

//    std::vector< CaloCluster > Allclusters; Allclusters.clear(); Allclusters=clusters;
//    std::vector< int > AllCristal; AllCristal.clear(); AllCristal=Ncristal_EB;
//    std::vector< float > s4s9_all; s4s9_all.clear(); s4s9_all=s4s9_EB;
//    std::vector< float > s1s9_all; s1s9_all.clear(); s1s9_all=s1s9_EB;
//    std::vector< float > s2s9_all; s2s9_all.clear(); s2s9_all=s2s9_EB;
//
//    for(unsigned int i=0; i<clustersEE.size(); i++){
//	  Allclusters.push_back(clustersEE[i]);
//	  AllCristal.push_back(Ncristal_EE[i]);
//	  s4s9_all.push_back( s4s9_EE[i] ); s1s9_all.push_back( s4s9_EE[i] ); s2s9_all.push_back( s4s9_EE[i] );
//    }

    std::map<size_t,size_t> savedCluEE; savedCluEE.clear();

    std::vector< CaloCluster > clustersEE_ES; clustersEE_ES.clear(); clustersEE_ES=clustersEE;
    if( !useES_E_ ){
	  clustersEE_ES.clear(); for(unsigned int i=0; i<eseeclusters_tot.size();i++ ) clustersEE_ES.push_back( eseeclusters_tot[i] );
	  Ncristal_EE.clear();   for(unsigned int i=0; i<Ncristal_EE_TOT.size();i++ )  Ncristal_EE.push_back( Ncristal_EE_TOT[i] );
	  s4s9_EE.clear();       for(unsigned int i=0; i<S4S9_TOT.size();i++ )         s4s9_EE.push_back( S4S9_TOT[i] );
	  Es_e1.clear();         for(unsigned int i=0; i<Es_e1_tot.size();i++ )        Es_e1.push_back( Es_e1_tot[i] );
	  Es_e2.clear();         for(unsigned int i=0; i<Es_e2_tot.size();i++ )        Es_e2.push_back( Es_e2_tot[i] );
    }
    else{
	  clustersEE_ES.clear(); for(unsigned int i=0; i<eseeclusters.size();i++ ) clustersEE_ES.push_back( eseeclusters[i] );
    }

    // loop over clusters to make Pi0
    for(size_t i=0; i< clustersEE_ES.size(); ++i) {

	  const CaloCluster* g1 = &(clustersEE_ES[i]);

	  // track association
	  math::XYZPoint impact1(999.,999.,999.);
	  //float dclustrk1(9999.);
	  //const reco::Track* trk1 = 0;     
	  // 2nd loop over clusters

	  for(size_t j=i+1; j<clustersEE_ES.size(); ++j) {

		// only store up to max number pi0 or clusters
	    if(nPi0>= NPI0MAX || nClu >= NCLUMAX) {

	      cout << "No more space in the arrays: nPi0: " << nPi0
		   << "  nClu: " << nClu << endl;
	      continue;
	    }

		const CaloCluster* g2 = &(clustersEE_ES[j]);           
		// using the origin
		math::PtEtaPhiMLorentzVector g1P4( g1->energy()/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
		math::PtEtaPhiMLorentzVector g2P4( g2->energy()/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
		math::PtEtaPhiMLorentzVector pi0P4 = g1P4 + g2P4;

		// use PV 
		TVector3 tmpV( g1->x()-xPV, g1->y()-yPV, g1->z()-zPV);
		math::PtEtaPhiMLorentzVector g1P4PV( g1->energy()/cosh(tmpV.Eta()), 
			  tmpV.Eta(), tmpV.Phi(), 0. );

		tmpV = TVector3( g2->x()-xPV, g2->y()-yPV, g2->z()-zPV);
		math::PtEtaPhiMLorentzVector g2P4PV( g2->energy()/cosh(tmpV.Eta()), 
			  tmpV.Eta(), tmpV.Phi(), 0. );
		math::PtEtaPhiMLorentzVector pi0P4PV = g1P4PV + g2P4PV;

		// basic selection to reduce combinatorics
		if( pi0P4PV.mass() > masspi0Cut_) continue;
		if( pi0P4PV.pt()   < ptpi0Cut_  ) continue;
		if(g1P4.eta()==g2P4.eta() && g1P4.phi()==g2P4.phi()) continue; 


		fillPi0Tree(
			  pi0P4PV, savedCluEE,
			  g1P4PV, g2P4PV,
			  shapesEE,
			  g1, g2, nClu, i, j, Ncristal_EE, EtSeed_EE);

		//new cut
		float nextClu = 999., Drtmp = 999.;
		for(size_t ind=0; ind<clustersEE_ES.size(); ++ind){
		    const CaloCluster* Gtmp = &(clustersEE_ES[ind]);
		    double deltaR1 = DR(Gtmp->phi(),g1P4.phi(),Gtmp->eta(),g1P4.eta());
		    double deltaR2 = DR(Gtmp->phi(),g2P4.phi(),Gtmp->eta(),g2P4.eta());
		    if( ind!=i && ind!=j && (deltaR1<Drtmp || deltaR2<Drtmp ) ){
			  nextClu = min(deltaR1,deltaR2);
			  Drtmp = nextClu;
		    }
		}
		if( clustersEE_ES.size()<3 ) nextClu=999.;
		STr2_IsoPi0_rec[nPi0] = nextClu;
		if( g1->energy()>g2->energy() ){
		    STr2_n1CrisPi0_rec[nPi0] = Ncristal_EE[i];
		    STr2_n2CrisPi0_rec[nPi0] = Ncristal_EE[j];
		}
		else{
		    STr2_n1CrisPi0_rec[nPi0] = Ncristal_EE[j];
		    STr2_n2CrisPi0_rec[nPi0] = Ncristal_EE[i];
		}
		STr2_mPi0_rec[nPi0] = pi0P4.mass();
		STr2_ptG1_rec[nPi0] = g1P4.pt();
		STr2_ptG2_rec[nPi0] = g2P4.pt();
		STr2_etaPi0_rec[nPi0] = pi0P4.Eta();
		STr2_ptPi0_rec[nPi0] = pi0P4.Pt();
		STr2_Pi0recIsEB[nPi0] = 2;
		STr2_DeltaRG1G2[nPi0] = DR(g1P4PV.phi(),g2P4PV.phi(),g1P4PV.eta(),g2P4PV.eta() );
		STr2_Es_e1_1[nPi0] = Es_e1[i];
		STr2_Es_e1_2[nPi0] = Es_e2[i];
		STr2_Es_e2_1[nPi0] = Es_e1[j];
		STr2_Es_e2_2[nPi0] = Es_e2[j];
		STr2_S4S9_1[nPi0]  = s4s9_EE[i];
		STr2_S4S9_2[nPi0]  = s4s9_EE[j];

		nPi0++;
		TLorentzVector thispar; thispar.SetPtEtaPhiE(pi0P4.pt(),pi0P4.eta(),pi0P4.phi(),pi0P4.energy());
		//MyEEParicles.push_back(thispar);

	  } // loop over clusters (g2)

    } // loop over clusters to make pi0 

    STr2_NPi0_rec = nPi0;

    Tree_HLT->Fill();
    Tree_HLT_clus->Fill();   


    //delete MyES;
    //delete MyES_clus;
}


//========================================================================================

template <class DetIdType, class RecHitCollectionType, class IteratorType> 
void NewPi0Dumper_Gun::make3x3Clusters(
	  const CaloGeometry* geometry, const RecHitCollectionType* ebHandle,
	  std::vector<CaloCluster>* clusters,
	  std::vector< ClusterShape >* shapes,
	  edm::Timestamp const & iTime, std::vector<int>* Ncristal, std::vector<float>* S4S9, std::vector<float>* S1S9, std::vector<float>* S2S9, std::vector<float>* EtrSeed )
{

    std::vector<EcalRecHit> seeds;
    seeds.clear();

    vector<DetId> usedXtals;
    usedXtals.clear();

    typedef std::map< DetId , bool > XtalInUse;
    XtalInUse isUsed; // map of which xtals have been used

    std::vector<DetId> detIdEBRecHits;
    std::vector<EcalRecHit> EBRecHits;

    // sort by energy and find the seeds
    for(IteratorType itb= ebHandle->begin(); itb != ebHandle->end(); ++itb) {
	  if(itb->energy() > 0.) {
		DetIdType EBoEE(itb->id());
		GlobalPoint posiz = geometry->getPosition( itb->id() );
		float EtSeeds = itb->energy()/cosh(posiz.eta());
//cout<<EtSeeds<<"  "<<itb->energy()<<" / "<<posiz.eta()<<endl;
		if( EBoEE.subdet() == EcalBarrel )     { if(itb->energy() > s1CluCut_ )         seeds.push_back( *itb ); }
		else if( EBoEE.subdet() == EcalEndcap) { if( EtSeeds > s1CluCutEE_ )        seeds.push_back( *itb ); }
		else cout<<"WARNING: cristal not in Barrel or in Endcap"<<endl;
	  }

    } // loop over xtals

    sort(seeds.begin(), seeds.end(), ecalRecHitLess());
    // loop over seeds and make clusters
    for (std::vector<EcalRecHit>::iterator itseed=seeds.begin();
		itseed!=seeds.end(); itseed++) {

	  DetId seed_id( itseed->id() );
	  DetIdType  seed_id_subdet(seed_id);

	  // check if seed already in use. If so go to next seed
	  std::map< DetId , bool >::const_iterator mapit = isUsed.find( seed_id );
	  if( mapit != isUsed.end() ) continue; // seed already in use

	  // find 3x3 matrix of xtals
	  int clusEtaSize_(3), clusPhiSize_(3);
	  std::vector<DetId> clus_v;

	  if( seed_id_subdet.subdet() == EcalBarrel ) 
		clus_v = ebTopology->getWindow(seed_id,clusEtaSize_,clusPhiSize_);       
	  else if( seed_id_subdet.subdet() == EcalEndcap ) 
		clus_v = eeTopology->getWindow(seed_id,clusEtaSize_,clusPhiSize_);       

	  // needed for position calculator
	  std::vector<std::pair<DetId,float> > clus_used;

	  // xtals actually used after removing those already used
	  vector<const EcalRecHit*> RecHitsInWindow; // 3x3 around seed. no sharing w/ other clusters. only unique xtals
	  vector<const EcalRecHit*> RecHitsInWindow3x3All; // 3x3 around seed. includes xtals from other clusters

	  float simple_energy = 0; 
	  float posTotalEnergy(0.); // need for position calculation

	  // make 3x3  cluster - reject overlaps
	  for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
		DetIdType thisId( *det );

		// find the rec hit
		IteratorType ixtal = ebHandle->find( thisId );
		if( ixtal == ebHandle->end() ) continue; // xtal not found
		RecHitsInWindow3x3All.push_back( &(*ixtal) );
		// skip this xtal if already used
		XtalInUse::const_iterator mapit = isUsed.find( thisId );
		if( mapit != isUsed.end() ) continue; // xtal already used

		RecHitsInWindow.push_back( &(*ixtal) );
		clus_used.push_back(std::make_pair(*det,1.));
		simple_energy +=  ixtal->energy();
		if(ixtal->energy()>0.) posTotalEnergy += ixtal->energy(); // use only pos energy for position
	  } // loop over xtals in the region

	  if(simple_energy <= 0) { 
		cout << "skipping cluster with negative energy " << simple_energy << endl; 
		continue;
	  }

	  float s4s9_tmp[4];
	  for(int i=0;i<4;i++){ 
		s4s9_tmp[i]= 0;
	  }

	  //        DetIdType  seed_id_subdet(seed_id);
	  int seed_ieta, seed_iphi;

	  if( seed_id_subdet.subdet() == EcalEndcap ) {
		seed_ieta = EEDetId(seed_id).ix();
		seed_iphi = EEDetId(seed_id).iy();

	  } else if( seed_id_subdet.subdet() == EcalBarrel ) {
		seed_ieta = EBDetId(seed_id).ieta();
		seed_iphi = EBDetId(seed_id).iphi();
		convxtalid( seed_iphi,seed_ieta);
	  }

	  // not sure works also for EE
#ifdef INDEXETA65
	  if( fabs(EBDetId(seed_id).ieta())>63 && fabs(EBDetId(seed_id).ieta())<67 )
		cout<<"iETA: "<<fabs(EBDetId(seed_id).ieta())<<" iPHI: "<<EBDetId(seed_id).iphi()<<" Index: "<<EBDetId(seed_id).hashedIndex()<<endl;
#endif
	  // energy of 3x3 cluster
	  float e3x3(0.);
	  int ncry9(0);
	  std::vector<std::pair<DetId,float> > enFracs;

	  // variables for position caculation
	  float xclu(0.), yclu(0.), zclu(0.); // temp var to compute weighted average
	  float trClu(0.), alphaClu(0.), laserCorrClu(0.);
	  //float trCluInvertZSide(0.);
	  float total_weight(0.);// to compute position
	  float total_TrWeight(0.);// to compute average transparency correction

	  double EnergyCristals[9] = {0.};
	  // loop over xtals and compute energy and position

	  for(unsigned int j=0; j<RecHitsInWindow.size();j++){

		DetIdType det(RecHitsInWindow[j]->id());

		int ieta, iphi;

		if( seed_id_subdet.subdet() == EcalEndcap ) {
		    //@@ ieta = EEDetId(seed_id).ix();
		    // iphi = EEDetId(seed_id).iy();
		    ieta = EEDetId(det).ix();
		    iphi = EEDetId(det).iy();

		} else if( seed_id_subdet.subdet() == EcalBarrel ) {
		    // ieta = EBDetId(seed_id).ieta();
		    // iphi = EBDetId(seed_id).iphi();
		    ieta = EBDetId(det).ieta();
		    iphi = EBDetId(det).iphi();
		    convxtalid(iphi,ieta);
		}
		else cout<<"WARNING: cristal not in Barrel or in Endcap"<<endl;

		// use calibration coeff for energy and position
		float en = RecHitsInWindow[j]->energy();
		int dx = diff_neta_s(seed_ieta,ieta);
		int dy = diff_nphi_s(seed_iphi,iphi);
		EnergyCristals[j] = en;

		if(abs(dx)<=1 && abs(dy)<=1) {
		    e3x3 += en;
		    ncry9++;
		    if(dx <= 0 && dy <=0){ s4s9_tmp[0] += en; }
		    if(dx >= 0 && dy <=0){ s4s9_tmp[1] += en; }
		    if(dx <= 0 && dy >=0){ s4s9_tmp[2] += en; }
		    if(dx >= 0 && dy >=0){ s4s9_tmp[3] += en; }
		    enFracs.push_back( std::make_pair( RecHitsInWindow[j]->id(), en ) );
		    // NOTA BENE: sto usando le frazioni per salvare energia rechit
		    isUsed[ RecHitsInWindow[j]->id() ] = true;
		}

		// compute position
		if(en>0.) {
		    float weight = std::max( float(0.), float(param_W0_ + log(en/posTotalEnergy)) );

		    GlobalPoint posThis = geometry->getPosition( det );

		    xclu += weight*posThis.x(); 
		    yclu += weight*posThis.y(); 
		    zclu += weight*posThis.z(); 
		    total_weight += weight;
		}

		//compute average transparency correction

	  } // loop over 3x3 rechits

	  float e2x2 = *max_element( s4s9_tmp,s4s9_tmp+4);
	  float s4s9 = e2x2/e3x3;
	  math::XYZPoint clusPos( xclu/total_weight, 
		    yclu/total_weight,
		    zclu/total_weight ); 

	  if(s4s9<s4s9CluCut_) continue;

	  //calculate e5x5 and fill energy fractions for 5x5 area - 3x3  already done
	  float e5x5 = 0;
	  int ncry25 = 0;
	  vector<const EcalRecHit*> RecHitsInWindow5x5; // 5x5 around seed. includes xtals from others
	  std::vector<DetId> clus_v5x5;
	  if( seed_id_subdet.subdet() == EcalBarrel ) 
		clus_v5x5 = ebTopology->getWindow(seed_id,5,5);       
	  else if( seed_id_subdet.subdet() == EcalEndcap ) 
		clus_v5x5 = eeTopology->getWindow(seed_id,5,5);       

	  for( std::vector<DetId>::const_iterator idItr = clus_v5x5.begin(); idItr != clus_v5x5.end(); idItr++){
		DetIdType det = *idItr;

		//inside collections
		IteratorType itdet = ebHandle->find(det);
		if(itdet == ebHandle->end()) continue;

		RecHitsInWindow5x5.push_back( &(*itdet) );
		e5x5 += itdet->energy();
		ncry25++;

		// check whether hit in 3x3 window - if not fraction = 0
		vector<const EcalRecHit*>::const_iterator in3x3 = find( RecHitsInWindow.begin(), RecHitsInWindow.end(), &(*itdet) ); 
		if( in3x3 == RecHitsInWindow.end() ) {
		    enFracs.push_back(std::make_pair(det,0.));
		}
	  }


	  vector<const EcalRecHit*> RecHitsInWindow7x7; // 7x7 around seed. includes xtals from others
	  std::vector<DetId> clus_v7x7;
	  if( seed_id_subdet.subdet() == EcalBarrel ) 
		clus_v7x7 = ebTopology->getWindow(seed_id,7,7);       
	  else if( seed_id_subdet.subdet() == EcalEndcap ) 
		clus_v7x7 = eeTopology->getWindow(seed_id,7,7);       

	  for( std::vector<DetId>::const_iterator idItr = clus_v7x7.begin(); idItr != clus_v7x7.end(); idItr++){
		DetIdType det = *idItr;

		//inside collections
		IteratorType itdet = ebHandle->find(det);
		if(itdet == ebHandle->end()) continue;
		RecHitsInWindow7x7.push_back( &(*itdet) );
	  }


	  // cluster shape
	  ClusterShape shape;
	  shape.s9 = e3x3;
	  shape.ncry9 = ncry9;
	  shape.s4 = e2x2;
	  shape.s4s9 = s4s9;
	  shape.s25 = e5x5;
	  shape.ncry25 = ncry25;
	  shape.s1 = itseed->energy();
	  shape.time = itseed->time();
	  shape.flag = itseed->recoFlag();
	  shape.transparencyLoss = float(trClu)/float(total_TrWeight);
	  shape.laserCorr = float(laserCorrClu)/float(total_TrWeight);
	  shape.alpha = float(alphaClu)/float(total_TrWeight);
	  //shape.transparencyLossInvertZSide = float(trCluInvertZSide)/float(total_TrWeight);

	  /// by now 2nd moments are implemented only for EB
	  if( seed_id_subdet.subdet() == EcalBarrel ) 
	  {
		// compute major and minor moments
		Cluster2ndMoments mom3x3 = EcalClusterTools::cluster2ndMoments( RecHitsInWindow3x3All );
		Cluster2ndMoments mom5x5 = EcalClusterTools::cluster2ndMoments( RecHitsInWindow5x5    );
		Cluster2ndMoments mom7x7 = EcalClusterTools::cluster2ndMoments( RecHitsInWindow7x7    );

		shape.sMaj9  = mom3x3.sMaj;
		shape.sMaj25 = mom5x5.sMaj;
		shape.sMaj49 = mom7x7.sMaj;
		shape.sMin9  = mom3x3.sMin;
		shape.sMin25 = mom5x5.sMin;
		shape.sMin49 = mom7x7.sMin;
	  }
	  else if( seed_id_subdet.subdet() == EcalEndcap )
	  {
		shape.sMaj9  = -999.999;
		shape.sMaj25 = -999.999;
		shape.sMaj49 = -999.999;
		shape.sMin9  = -999.999;
		shape.sMin25 = -999.999;
		shape.sMin49 = -999.999;
	  }

	  // compute pt of gamma and cut
	  float ptClus = e3x3/cosh(clusPos.eta());

	  if(ptClus<ptCluCut_) continue;

	  // energy corrections in bins of eta and energy
	  //if(!noCorrections_) e3x3 *= energyCorrection(e3x3,clusPos.eta());

	  // make calo clusters
	  Ncristal->push_back( RecHitsInWindow.size() ); //for the Branch
	  S4S9->push_back( s4s9 ); //for the Branch
	  S1S9->push_back( itseed->energy()/e3x3 ); //for the Branch
	  int IndMax=-1;
	  double maxEne = max_array( EnergyCristals, 9 );
	  for(int i=0; i<9; i++){ if(EnergyCristals[i]==maxEne) IndMax = i; }
	  for(int i=0; i<9; i++){ if(i == IndMax) EnergyCristals[i]=0.; }
	  double maxEne2 = max_array( EnergyCristals, 9);
	  S2S9->push_back( (maxEne+maxEne2)/e3x3 ); //for the Branch
	  //if(seed_id_subdet.subdet() == EcalBarrel) Nclus_ieta->push_back( seed_ieta );
	  //Nclus_id->push_back( seed_id );
	  //if(seed_id_subdet.subdet() == EcalEndcap ) cout<<"Seed "<<EndcapTools::getRingIndex( seed_id )<<endl;
	  clusters->push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_BARREL),
			  enFracs, CaloCluster::undefined, seed_id ) );
	  shapes->push_back( shape );

	  //recompute position just for Et
	  GlobalPoint posi;
	  if(itseed->energy() > 0.) posi = geometry->getPosition( itseed->id() );
	  EtrSeed->push_back( itseed->energy()/cosh(posi.eta())  ); //for the Branch

    } //loop over seeds to make clusters

}

//--------------------------------------------------------------------------------------------
void NewPi0Dumper_Gun::FillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //--------------------------------------------------------------------------------------------

  using namespace edm;
  using namespace reco;

  // trigger Handles
  Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel( l1TriggerInput_, gtReadoutRecord);

  ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);

  /// L1 trigger results for physics algorithms
  //const L1GtTriggerMenu* menu = menuRcd.product();
  const DecisionWord& gtDecisionWordBeforeMask = gtReadoutRecord->decisionWord();

  isBSC = false;

  nL1bits = int( gtDecisionWordBeforeMask.size() );

  for (int iBit = 0; iBit < nL1bits; ++iBit)
    L1bits[iBit] = gtDecisionWordBeforeMask[iBit] ;
  /// L1 technical triggers
  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
  nL1bitsTech = int(technicalTriggerWordBeforeMask.size());
  for(int iBit = 0; iBit < nL1bitsTech; ++iBit) {
    L1bitsTech[iBit] = technicalTriggerWordBeforeMask.at(iBit);
  }

  isBSC = L1bitsTech[0] && (L1bitsTech[40] || L1bitsTech[41]) 
    && ! (L1bitsTech[36] || L1bitsTech[37] || L1bitsTech[38] || L1bitsTech[39]) 
    && ! ((L1bitsTech[42] && ! L1bitsTech[43]) || (L1bitsTech[43] && ! L1bitsTech[42])) ;


} // end of FillTriggerInfo

/*
   float NewPi0Dumper_Gun::getLaserCorrection(DetId const & xid, edm::Timestamp const & iTime, bool invertZSide) const
   {
   float correctionFactor = 1.0;

   if (!mAPDPNRatios_ || !mAPDPNRatiosRef_)
   {
   edm::LogError("EcalLaserDbService") << "Laser DB info not found" << endl;
   return correctionFactor;
   }

   const EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap& laserRatiosMap =  mAPDPNRatios_->getLaserMap();
   const EcalLaserAPDPNRatios::EcalLaserTimeStampMap& laserTimeMap =  mAPDPNRatios_->getTimeMap();
   const EcalLaserAPDPNRatiosRefMap& laserRefMap =  mAPDPNRatiosRef_->getMap();
   const EcalLaserAlphaMap& laserAlphaMap =  mAlphas_->getMap();

   EcalLaserAPDPNRatios::EcalLaserAPDPNpair apdpnpair;
   EcalLaserAPDPNRatios::EcalLaserTimeStamp timestamp;
   EcalLaserAPDPNref apdpnref;
   EcalLaserAlpha alpha;

   if (xid.det()==DetId::Ecal) {
//    std::cout << " XID is in Ecal : ";
} else {
//    std::cout << " XID is NOT in Ecal : ";
edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL" << endl;
return correctionFactor;
} 

int hi = -1;
if (xid.subdetId()==EcalBarrel) {
//    std::cout << "EcalBarrel" << std::endl;
//    std::cout << "--> rawId() = " << xid.rawId() << "   id() = " << EBDetId( xid ).hashedIndex() << std::endl;
hi = EBDetId( xid ).hashedIndex();
} else if (xid.subdetId()==EcalEndcap) {
//    std::cout << "EcalEndcap" << std::endl;
hi = EEDetId( xid  ).hashedIndex() + EBDetId::MAX_HASH + 1;
} else {
//    std::cout << "NOT EcalBarrel or EcalEndCap" << std::endl;
edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL Barrel or Endcap" << endl;
return correctionFactor;
}

int iLM;
if (xid.subdetId()==EcalBarrel) {
EBDetId ebid( xid.rawId() );
iLM = MEEBGeom::lmr(ebid.ieta(), ebid.iphi());
} else if (xid.subdetId()==EcalEndcap) {
EEDetId eeid( xid.rawId() );
// SuperCrystal coordinates
MEEEGeom::SuperCrysCoord iX = (eeid.ix()-1)/5 + 1;
MEEEGeom::SuperCrysCoord iY = (eeid.iy()-1)/5 + 1;    
iLM = MEEEGeom::lmr(iX, iY, eeid.zside());    
} else {
edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL Barrel or Endcap" << endl;
return correctionFactor;
}
//  std::cout << " LM num ====> " << iLM << endl;

// get alpha, apd/pn ref, apd/pn pairs and timestamps for interpolation
DetId myId;
if (xid.subdetId()==EcalEndcap && invertZSide)
myId=EEDetId(EEDetId(xid).ix(),EEDetId(xid).iy(),-1*EEDetId(xid).zside());
else
myId=xid;

EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap::const_iterator itratio = laserRatiosMap.find(myId);
if (itratio != laserRatiosMap.end()) {
apdpnpair = (*itratio);
} else {
edm::LogError("EcalLaserDbService") << "error with laserRatiosMap!" << endl;     
return correctionFactor;
}

if (iLM-1< (int)laserTimeMap.size()) {
    timestamp = laserTimeMap[iLM-1];  
} else {
    edm::LogError("EcalLaserDbService") << "error with laserTimeMap!" << endl;     
    return correctionFactor;
}

EcalLaserAPDPNRatiosRefMap::const_iterator itref = laserRefMap.find(xid);
if ( itref != laserRefMap.end() ) {
    apdpnref = (*itref);
} else { 
    edm::LogError("EcalLaserDbService") << "error with laserRefMap!" << endl;     
    return correctionFactor;
}

EcalLaserAlphaMap::const_iterator italpha = laserAlphaMap.find(xid);
if ( italpha != laserAlphaMap.end() ) {
    alpha = (*italpha);
} else {
    edm::LogError("EcalLaserDbService") << "error with laserAlphaMap!" << endl;     
    return correctionFactor;
}

//   std::cout << " APDPN pair " << apdpnpair.p1 << " , " << apdpnpair.p2 << std::endl; 
//   std::cout << " TIME pair " << timestamp.t1.value() << " , " << timestamp.t2.value() << " iLM " << iLM << std::endl; 
//   std::cout << " LM module " << iLM << std::endl;
//   std::cout << " APDPN ref " << apdpnref << std::endl; 

//    std::cout << " ALPHA " << alpha << std::endl; 

// should implement some default in case of error...

// should do some quality checks first
// ...

// we will need to treat time differently...
// is time in DB same format as in MC?  probably not...

// interpolation

edm::TimeValue_t t = iTime.value();
edm::TimeValue_t t_i = 0, t_f = 0;
float p_i = 0, p_f = 0;

if ( t >= timestamp.t1.value() && t < timestamp.t2.value() ) {
    t_i = timestamp.t1.value();
    t_f = timestamp.t2.value();
    p_i = apdpnpair.p1;
    p_f = apdpnpair.p2;
} else if ( t >= timestamp.t2.value() && t <= timestamp.t3.value() ) {
    t_i = timestamp.t2.value();
    t_f = timestamp.t3.value();
    p_i = apdpnpair.p2;
    p_f = apdpnpair.p3;
} else if ( t < timestamp.t1.value() ) {
    t_i = timestamp.t1.value();
    t_f = timestamp.t2.value();
    p_i = apdpnpair.p1;
    p_f = apdpnpair.p2;
    //edm::LogWarning("EcalLaserDbService") << "The event timestamp t=" << t 
    //        << " is lower than t1=" << t_i << ". Extrapolating...";
} else if ( t > timestamp.t3.value() ) {
    t_i = timestamp.t2.value();
    t_f = timestamp.t3.value();
    p_i = apdpnpair.p2;
    p_f = apdpnpair.p3;
    //edm::LogWarning("EcalLaserDbService") << "The event timestamp t=" << t 
    //        << " is greater than t3=" << t_f << ". Extrapolating...";
}

if ( apdpnref != 0 && (t_i - t_f) != 0) {
    float interpolatedLaserResponse = p_i/apdpnref + (t-t_i)*(p_f-p_i)/apdpnref/(t_f-t_i);
    if ( interpolatedLaserResponse <= 0 ) {
	  edm::LogError("EcalLaserDbService") << "The interpolated laser correction is <= zero! (" 
		<< interpolatedLaserResponse << "). Using 1. as correction factor.";
	  return correctionFactor;
    } else {
	  //correctionFactor = interpolatedLaserResponse;
	  correctionFactor = 1/pow(interpolatedLaserResponse,alpha);
    }
    //  std::cout << "correction factor " << correctionFactor << std::endl;
} else {
    edm::LogError("EcalLaserDbService") 
	  << "apdpnref (" << apdpnref << ") "
	  << "or t_i-t_f (" << (t_i - t_f) << " is zero!";
    return correctionFactor;
}

return correctionFactor;
}


pair<float,float>
NewPi0Dumper_Gun::getTransparencyLoss(DetId const & xid, edm::Timestamp const & iTime, bool invertZSide) const
{
    //float correctionFactor = 1.0;
    pair<float,float> correctionFactor(1.,1.); //transparencyLoss, alpha

    if (!mAPDPNRatios_ || !mAPDPNRatiosRef_)
    {
	  edm::LogError("EcalLaserDbService") << "Laser DB info not found" << endl;
	  return correctionFactor;
    }

    const EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap& laserRatiosMap =  mAPDPNRatios_->getLaserMap();
    const EcalLaserAPDPNRatios::EcalLaserTimeStampMap& laserTimeMap =  mAPDPNRatios_->getTimeMap();
    const EcalLaserAPDPNRatiosRefMap& laserRefMap =  mAPDPNRatiosRef_->getMap();
    const EcalLaserAlphaMap& laserAlphaMap =  mAlphas_->getMap();

    EcalLaserAPDPNRatios::EcalLaserAPDPNpair apdpnpair;
    EcalLaserAPDPNRatios::EcalLaserTimeStamp timestamp;
    EcalLaserAPDPNref apdpnref;
    EcalLaserAlpha alpha;

    if (xid.det()==DetId::Ecal) {
	  //    std::cout << " XID is in Ecal : ";
    } else {
	  //    std::cout << " XID is NOT in Ecal : ";
	  edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL" << endl;
	  return correctionFactor;
    } 

    int hi = -1;
    if (xid.subdetId()==EcalBarrel) {
	  //    std::cout << "EcalBarrel" << std::endl;
	  //    std::cout << "--> rawId() = " << xid.rawId() << "   id() = " << EBDetId( xid ).hashedIndex() << std::endl;
	  hi = EBDetId( xid ).hashedIndex();
    } else if (xid.subdetId()==EcalEndcap) {
	  //    std::cout << "EcalEndcap" << std::endl;
	  hi = EEDetId( xid  ).hashedIndex() + EBDetId::MAX_HASH + 1;
    } else {
	  //    std::cout << "NOT EcalBarrel or EcalEndCap" << std::endl;
	  edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL Barrel or Endcap" << endl;
	  return correctionFactor;
    }

    int iLM;
    if (xid.subdetId()==EcalBarrel) {
	  EBDetId ebid( xid.rawId() );
	  iLM = MEEBGeom::lmr(ebid.ieta(), ebid.iphi());
    } else if (xid.subdetId()==EcalEndcap) {
	  EEDetId eeid( xid.rawId() );
	  // SuperCrystal coordinates
	  MEEEGeom::SuperCrysCoord iX = (eeid.ix()-1)/5 + 1;
	  MEEEGeom::SuperCrysCoord iY = (eeid.iy()-1)/5 + 1;    
	  iLM = MEEEGeom::lmr(iX, iY, eeid.zside());    
    } else {
	  edm::LogError("EcalLaserDbService") << " DetId is NOT in ECAL Barrel or Endcap" << endl;
	  return correctionFactor;
    }
    //  std::cout << " LM num ====> " << iLM << endl;

    // get alpha, apd/pn ref, apd/pn pairs and timestamps for interpolation
    DetId myId;
    if (xid.subdetId()==EcalEndcap && invertZSide)
	  myId=EEDetId(EEDetId(xid).ix(),EEDetId(xid).iy(),-1*EEDetId(xid).zside());
    else
	  myId=xid;

    EcalLaserAPDPNRatios::EcalLaserAPDPNRatiosMap::const_iterator itratio = laserRatiosMap.find(myId);
    if (itratio != laserRatiosMap.end()) {
	  apdpnpair = (*itratio);
    } else {
	  edm::LogError("EcalLaserDbService") << "error with laserRatiosMap!" << endl;     
	  return correctionFactor;
    }

    if (iLM-1< (int)laserTimeMap.size()) {
	  timestamp = laserTimeMap[iLM-1];  
    } else {
	  edm::LogError("EcalLaserDbService") << "error with laserTimeMap!" << endl;     
	  return correctionFactor;
    }

    EcalLaserAPDPNRatiosRefMap::const_iterator itref = laserRefMap.find(xid);
    if ( itref != laserRefMap.end() ) {
	  apdpnref = (*itref);
    } else { 
	  edm::LogError("EcalLaserDbService") << "error with laserRefMap!" << endl;     
	  return correctionFactor;
    }

    EcalLaserAlphaMap::const_iterator italpha = laserAlphaMap.find(xid);
    if ( italpha != laserAlphaMap.end() ) {
	  alpha = (*italpha);
    } else {
	  edm::LogError("EcalLaserDbService") << "error with laserAlphaMap!" << endl;     
	  return correctionFactor;
    }

    //   std::cout << " APDPN pair " << apdpnpair.p1 << " , " << apdpnpair.p2 << std::endl; 
    //   std::cout << " TIME pair " << timestamp.t1.value() << " , " << timestamp.t2.value() << " iLM " << iLM << std::endl; 
    //   std::cout << " LM module " << iLM << std::endl;
    //   std::cout << " APDPN ref " << apdpnref << std::endl; 

    //    std::cout << " ALPHA " << alpha << std::endl; 

    // should implement some default in case of error...

    // should do some quality checks first
    // ...

    // we will need to treat time differently...
    // is time in DB same format as in MC?  probably not...

    // interpolation

    edm::TimeValue_t t = iTime.value();
    edm::TimeValue_t t_i = 0, t_f = 0;
    float p_i = 0, p_f = 0;

    if ( t >= timestamp.t1.value() && t < timestamp.t2.value() ) {
	  t_i = timestamp.t1.value();
	  t_f = timestamp.t2.value();
	  p_i = apdpnpair.p1;
	  p_f = apdpnpair.p2;
    } else if ( t >= timestamp.t2.value() && t <= timestamp.t3.value() ) {
	  t_i = timestamp.t2.value();
	  t_f = timestamp.t3.value();
	  p_i = apdpnpair.p2;
	  p_f = apdpnpair.p3;
    } else if ( t < timestamp.t1.value() ) {
	  t_i = timestamp.t1.value();
	  t_f = timestamp.t2.value();
	  p_i = apdpnpair.p1;
	  p_f = apdpnpair.p2;
	  //edm::LogWarning("EcalLaserDbService") << "The event timestamp t=" << t 
	  //        << " is lower than t1=" << t_i << ". Extrapolating...";
    } else if ( t > timestamp.t3.value() ) {
	  t_i = timestamp.t2.value();
	  t_f = timestamp.t3.value();
	  p_i = apdpnpair.p2;
	  p_f = apdpnpair.p3;
	  //edm::LogWarning("EcalLaserDbService") << "The event timestamp t=" << t 
	  //        << " is greater than t3=" << t_f << ". Extrapolating...";
    }

    if ( apdpnref != 0 && (t_i - t_f) != 0) {
	  float interpolatedLaserResponse = p_i/apdpnref + (t-t_i)*(p_f-p_i)/apdpnref/(t_f-t_i);
	  if ( interpolatedLaserResponse <= 0 ) {
		edm::LogError("EcalLaserDbService") << "The interpolated laser correction is <= zero! (" 
		    << interpolatedLaserResponse << "). Using 1. as correction factor.";
		return correctionFactor;
	  } else {
		//correctionFactor = interpolatedLaserResponse;
		//correctionFactor = 1/pow(interpolatedLaserResponse,alpha);
		correctionFactor.first = interpolatedLaserResponse;
		correctionFactor.second = alpha;
	  }
	  //  std::cout << "correction factor " << correctionFactor << std::endl;
    } else {
	  edm::LogError("EcalLaserDbService") 
		<< "apdpnref (" << apdpnref << ") "
		<< "or t_i-t_f (" << (t_i - t_f) << " is zero!";
	  return correctionFactor;
    }

    return correctionFactor;
}
*/
void
NewPi0Dumper_Gun::fillPi0Tree(
	  math::PtEtaPhiMLorentzVector &pi0P4PV,  std::map<size_t,size_t> &savedCluEB,
	  math::PtEtaPhiMLorentzVector &g1P4PV, math::PtEtaPhiMLorentzVector &g2P4PV,
	  std::vector< ClusterShape > &shapes,
	  const CaloCluster* g1, const CaloCluster* g2, int &nClu, int i, int j, vector<int> &Ncristal, vector<float>& EtSeeds)
{

    massPi0[nPi0] =  pi0P4PV.mass();
    ePi0[nPi0] =  pi0P4PV.energy();
    ptPi0[nPi0] =  pi0P4PV.pt();
    etaPi0[nPi0] =  pi0P4PV.eta();
    phiPi0[nPi0] =  pi0P4PV.phi();
    CristTot[nPi0] =  Ncristal[i]+Ncristal[j];
    DRClus[nPi0] =  DR(g1P4PV.phi(),g2P4PV.phi(),g1P4PV.eta(),g2P4PV.eta());

    // 1st daugher index
    std::map<size_t,size_t>::const_iterator ind = savedCluEB.find( i );
    if( ind == savedCluEB.end() ) {
	  savedCluEB[ i ]  = nClu;
	  S9Clu[nClu] = g1P4PV.energy();
	  nCris[nClu] = Ncristal[i];
	  EtSeed[nClu] = EtSeeds[i];
	  ptClu[nClu] =  g1P4PV.pt();
	  etaClu[nClu] = g1P4PV.eta();
	  phiClu[nClu] = g1P4PV.phi();
	  S1Clu[nClu] =  shapes[i].s1;
	  S4Clu[nClu] =  shapes[i].s4;
	  S25Clu[nClu] = shapes[i].s25;
	  timeClu[nClu] = shapes[i].time;
	  nCryClu[nClu] = shapes[i].ncry9;
	  flagClu[nClu] = shapes[i].flag;
	  sMaj9Clu[nClu]  = shapes[i].sMaj9;
	  sMaj25Clu[nClu] = shapes[i].sMaj25;
	  sMaj49Clu[nClu] = shapes[i].sMaj49;
	  sMin9Clu[nClu]  = shapes[i].sMin9;
	  sMin25Clu[nClu] = shapes[i].sMin25;
	  sMin49Clu[nClu] = shapes[i].sMin49;

	  if( g1->seed().subdetId() == EcalBarrel ) 
	  {
		EBDetId sid( g1->seed() );
		ietaClu[nClu] = sid.ieta();
		iphiClu[nClu] = sid.iphi();
		iCryClu[nClu] = sid.ic();
		iSMClu[nClu] = sid.ism();
		imodClu[nClu] = sid.im();
		iTTClu[nClu] = sid.tower().iTT();
		iTTetaClu[nClu] = sid.tower_ieta();
		iTTphiClu[nClu] = sid.tower_iphi();
	  }
	  else if( g1->seed().subdetId() == EcalEndcap ) 
	  {  
		EEDetId sid( g1->seed() );
		ietaClu[nClu] = sid.ix();
		iphiClu[nClu] = sid.iy();
		iCryClu[nClu] = sid.ic();
		iSMClu[nClu] = sid.iquadrant();
		imodClu[nClu] = sid.isc();
		iTTClu[nClu] =    -999; //sid.tower().iTT();
		iTTetaClu[nClu] = -999; //sid.tower_ieta();
		iTTphiClu[nClu] = -999; //sid.tower_iphi();
	  }

	  nClu++;
    }

    // 2nd daugher index
    ind = savedCluEB.find( j );
    if( ind == savedCluEB.end() ) {
	  savedCluEB[ j ]  = nClu;
	  S9Clu[nClu] = g2P4PV.energy();
	  nCris[nClu] = Ncristal[j];
	  EtSeed[nClu] = EtSeeds[j]; 
	  ptClu[nClu] =  g2P4PV.pt();
	  etaClu[nClu] = g2P4PV.eta();
	  phiClu[nClu] = g2P4PV.phi();
	  S1Clu[nClu] =  shapes[j].s1;
	  S4Clu[nClu] =  shapes[j].s4;
	  S25Clu[nClu] = shapes[j].s25;
	  timeClu[nClu] = shapes[j].time;
	  nCryClu[nClu] = shapes[j].ncry9;
	  flagClu[nClu] = shapes[j].flag;
	  sMaj9Clu[nClu]  = shapes[j].sMaj9;
	  sMaj25Clu[nClu] = shapes[j].sMaj25;
	  sMaj49Clu[nClu] = shapes[j].sMaj49;
	  sMin9Clu[nClu]  = shapes[j].sMin9;
	  sMin25Clu[nClu] = shapes[j].sMin25;
	  sMin49Clu[nClu] = shapes[j].sMin49;

	  if( g2->seed().subdetId() == EcalBarrel ) 
	  {
		EBDetId sid( g2->seed() );
		ietaClu[nClu] = sid.ieta();
		iphiClu[nClu] = sid.iphi();
		iCryClu[nClu] = sid.ic();
		iSMClu[nClu] = sid.ism();
		imodClu[nClu] = sid.im();
		iTTClu[nClu] = sid.tower().iTT();
		iTTetaClu[nClu] = sid.tower_ieta();
		iTTphiClu[nClu] = sid.tower_iphi();
	  }
	  else if( g2->seed().subdetId() == EcalEndcap ) 
	  {
		EEDetId sid( g2->seed() );
		ietaClu[nClu] = sid.ix();
		iphiClu[nClu] = sid.iy();
		iCryClu[nClu] = sid.ic();
		iSMClu[nClu] = sid.iquadrant();
		imodClu[nClu] = sid.isc();
		iTTClu[nClu] =    -999; //sid.tower().iTT();
		iTTetaClu[nClu] = -999; //sid.tower_ieta();
		iTTphiClu[nClu] = -999; //sid.tower_iphi();
	  }
	  nClu++;
    }

    indexClu1Pi0[nPi0] = savedCluEB[i];
    indexClu2Pi0[nPi0] = savedCluEB[j];

    // 1st daughter always with higher energy
    if(g2P4PV.pt() > g1P4PV.pt()) {
	  indexClu1Pi0[nPi0] = savedCluEB[j];
	  indexClu2Pi0[nPi0] = savedCluEB[i];
    }

}//Fill TTree

bool
NewPi0Dumper_Gun::GetHLTResults(const edm::Event& iEvent, std::string s){

    edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
    iEvent.getByLabel(triggerTag_, hltTriggerResultHandle);

    edm::TriggerNames HLTNames;

    HLTNames = iEvent.triggerNames(*hltTriggerResultHandle);
    std::string tempnames;
    int hltCount = hltTriggerResultHandle->size();
    TRegexp reg(TString( s.c_str()) );

    for (int i = 0 ; i != hltCount; ++i) {
	  TString hltName_tstr(HLTNames.triggerName(i));
//cout<<"Trigger: "<<hltName_tstr<<endl;
	  std::string hltName_str(HLTNames.triggerName(i));
		if ( hltName_tstr.Contains(reg) ){
//		    cout<<hltTriggerResultHandle->accept(i)<<endl;
		    return hltTriggerResultHandle->accept(i);
            }
	  }
    return false;
} // HLT isValid



double NewPi0Dumper_Gun::DPhi(double phi1, double phi2)
{
    double deltaPhi = phi1 - phi2;
    if (deltaPhi > Geom::pi()) deltaPhi -= 2.*Geom::pi();
    if (deltaPhi < -Geom::pi()) deltaPhi += 2.*Geom::pi();

    return deltaPhi;
}
double NewPi0Dumper_Gun::DR( double phi1, double phi2, double eta1, double eta2)
{
    double deltaEta = eta1-eta2;
    double deltaPhi = DPhi(phi1,phi2);
    return std::sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
}
double NewPi0Dumper_Gun::min( double a, double b)
{
    if( a<=b ) return a;
    else       return b;
}

double max_array(double *A, int n){ 

    if(n==1)  return A[0]; 
    else      return max(A[0],max_array(A+1, n-1)); 
} 

double max (double x, double y) { 
    if(x>=y)  return x; 
    else      return y; 
}

void convxtalid(Int_t &nphi,Int_t &neta)
{

    if(neta > 0) neta -= 1;
    if(nphi > 359) nphi=nphi-360;

} //end of convxtalid

int diff_neta_s(Int_t neta1, Int_t neta2){
    Int_t mdiff;
    mdiff=(neta1-neta2);
    return mdiff;
}


int diff_nphi_s(Int_t nphi1,Int_t nphi2) {
    Int_t mdiff;
    if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
	  mdiff=nphi1-nphi2;
    }
    else {
	  mdiff=360-abs(nphi1-nphi2);
	  if(nphi1>nphi2) mdiff=-mdiff;
    }
    return mdiff;
}

//define this as a plug-in
DEFINE_FWK_MODULE(NewPi0Dumper_Gun);
