// system include files
#include <memory>
#include "DataFormats/GeometrySurface/interface/Line.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"

#include <RecoVertex/VertexPrimitives/interface/BasicVertexState.h>
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/SharedTracks.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/CandidatePtrTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "seedAnalyzer/seedAnalyzer/interface/trackVars2.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"


#include "RecoBTag/TrackProbability/interface/HistogramProbabilityEstimator.h"

class HistogramProbabilityEstimator;
#include <typeinfo>


#include "CondFormats/BTauObjects/interface/TrackProbabilityCalibration.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability2DRcd.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability3DRcd.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"

#include <stdio.h>   
#include <math.h>      


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

// 



class Analyzer_MINIAOD_new : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Analyzer_MINIAOD_new(const edm::ParameterSet&);
      ~Analyzer_MINIAOD_new();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      // ----------member data ---------------------------
      
      int jet_flavour_function(const pat::Jet& jet, bool usePhysForLightAndUndefined);
      double deltaPhi(double phi1, double phi2);
      void trackToTrack(std::vector<reco::TransientTrack> selTracks, std::vector<reco::TransientTrack>::const_iterator it, const pat::Jet &jet, const reco::Vertex &pv, std::vector<trackVars2>& nearTracks, std::vector<float> masses, int index);
      void checkEventSetup(const edm::EventSetup & iSetup);
      
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > CandidateToken;
      edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;
      edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
      edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
      edm::EDGetTokenT<double> rhoToken;
      
      edm::Service<TFileService> file;
      TTree *tree;
      
      std::auto_ptr<HistogramProbabilityEstimator> m_probabilityEstimator;
      bool m_computeProbabilities=1;
      unsigned long long  m_calibrationCacheId2D; 
      unsigned long long m_calibrationCacheId3D;
      
      double pv_x=-10;
      double pv_y=-10;
      double pv_z=-10;
      
      int evt=0;
      int lumi=0;
      int run=0;
      int nPVs=0;
      double rho=0.0;
      
      double jetpt;
      double jetptLOG;      
      double jeteta;
      double jetphi;
      double jetmass;
      
      double jetCSVv2;
      double jetCMVA;
      double jetTCHE;
      double jetSSVHE;
      double jetDeepCSV;
      double jetDeepFlavour;      
      
      int jetflavour;
      int jetNseeds;
      
      double seed_pt[10];
      double seed_eta[10];
      double seed_phi[10];
      double seed_mass[10];      
      double seed_dz[10];
      double seed_dxy[10];
      double seed_3D_ip[10];
      double seed_3D_sip[10];
      double seed_2D_ip[10];
      double seed_2D_sip[10];
      double seed_3D_signedIp[10];
      double seed_3D_signedSip[10];
      double seed_2D_signedIp[10];
      double seed_2D_signedSip[10];
      double seed_3D_TrackProbability[10];
      double seed_2D_TrackProbability[10];
      double seed_chi2reduced[10];
      double seed_nPixelHits[10];
      double seed_nHits[10];
      double seed_jetAxisDistance[10];
      double seed_jetAxisDlength[10];      
      
      //saving also some transformation of the variables
      double seed_ptLog[10];
      double seed_dzSymLog[10];
      double seed_dxySymLog[10];
      double seed_3D_ipLog[10];
      double seed_3D_sipLog[10];
      double seed_2D_ipLog[10];
      double seed_2D_sipLog[10];
      double seed_3D_signedIpSymLog[10];
      double seed_3D_signedSipSymLog[10];
      double seed_2D_signedIpSymLog[10];
      double seed_2D_signedSipSymLog[10];
      double seed_jetAxisDistanceLog[10];
      double seed_jetAxisDlengthLog[10];      
      
      double nearTracks_pt[200];
      double nearTracks_eta[200];
      double nearTracks_phi[200];
      double nearTracks_mass[200];
      double nearTracks_dz[200];
      double nearTracks_dxy[200];
      double nearTracks_3D_ip[200];
      double nearTracks_3D_sip[200];
      double nearTracks_2D_ip[200];
      double nearTracks_2D_sip[200];
      double nearTracks_PCAdist[200];
      double nearTracks_PCAdsig[200];      
      double nearTracks_PCAonSeed_x[200];
      double nearTracks_PCAonSeed_y[200];
      double nearTracks_PCAonSeed_z[200];      
      double nearTracks_PCAonSeed_xerr[200];
      double nearTracks_PCAonSeed_yerr[200];
      double nearTracks_PCAonSeed_zerr[200];      
      double nearTracks_PCAonTrack_x[200];
      double nearTracks_PCAonTrack_y[200];
      double nearTracks_PCAonTrack_z[200];      
      double nearTracks_PCAonTrack_xerr[200];
      double nearTracks_PCAonTrack_yerr[200];
      double nearTracks_PCAonTrack_zerr[200]; 
      double nearTracks_dotprodTrack[200];
      double nearTracks_dotprodSeed[200];
      double nearTracks_dotprodTrackSeed2D[200];
      double nearTracks_dotprodTrackSeed3D[200];
      double nearTracks_dotprodTrackSeedVectors2D[200];
      double nearTracks_dotprodTrackSeedVectors3D[200];      
      double nearTracks_PCAonSeed_pvd[200];
      double nearTracks_PCAonTrack_pvd[200];
      double nearTracks_PCAjetAxis_dist[200];
      double nearTracks_PCAjetMomenta_dotprod[200];
      double nearTracks_PCAjetDirs_DEta[200];
      double nearTracks_PCAjetDirs_DPhi[200];
      
      //saving also some transformation of the variables
      double nearTracks_ptLog[200];
      double nearTracks_dzLog[200];
      double nearTracks_dxyLog[200];
      double nearTracks_3D_ipLog[200];
      double nearTracks_3D_sipLog[200];
      double nearTracks_2D_ipLog[200];
      double nearTracks_2D_sipLog[200];
      double nearTracks_PCAdistLog[200];
      double nearTracks_PCAdsigLog[200];  
      double nearTracks_PCAonSeed_pvdLog[200];
      double nearTracks_PCAonTrack_pvdLog[200];
      double nearTracks_PCAjetAxis_distLog[200];
      double nearTracks_PCAjetDirs_DEtaLog[200];
      
      
      
      float min3DIPValue=0.005;
      float min3DIPSignificance=1.2;
      int max3DIPValue=9999.;
      int max3DIPSignificance=9999.;     

      
      std::vector<trackVars2> nearTracks;
      std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2>>> SortedSeedsMap;
    
      
};


Analyzer_MINIAOD_new::Analyzer_MINIAOD_new(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
    tree=file->make<TTree>("tree","tree");
    
    tree->Branch("lumi",&lumi, "lumi/I");
    tree->Branch("evt",&evt, "evt/I");
    tree->Branch("run",&run, "run/I");
    tree->Branch("nPVs",&nPVs, "nPVs/I");
    tree->Branch("rho",&rho, "rho/D");

    tree->Branch("pv_x",&pv_x, "pv_x/D");
    tree->Branch("pv_y",&pv_y, "pv_y/D");
    tree->Branch("pv_z",&pv_z, "pv_z/D");
    
    tree->Branch("jet_pt",&jetpt, "jet_pt/D");
    tree->Branch("jet_ptLOG",&jetptLOG, "jet_ptLOG/D");
    tree->Branch("jet_eta",&jeteta, "jet_eta/D");
    tree->Branch("jet_phi",&jetphi, "jet_phi/D");
    tree->Branch("jet_mass",&jetmass, "jet_mass/D");        
    tree->Branch("jet_CSVv2", &jetCSVv2,  "jet_CSVv2/D");
    tree->Branch("jet_CMVA", &jetCMVA,  "jet_CMVA/D");
    tree->Branch("jet_TCHE", &jetTCHE,  "jet_TCHE/D");
    tree->Branch("jet_SSVHE", &jetSSVHE,  "jet_SSVHE/D");    
    tree->Branch("jet_DeepCSV", &jetDeepCSV,  "jet_DeepCSV/D");
    tree->Branch("jet_DeepFlavour", &jetDeepFlavour,  "jet_DeepFlavour/D");
    tree->Branch("jet_flavour",&jetflavour, "jet_flavour/I");
    tree->Branch("jet_Nseeds",&jetNseeds, "jet_Nseeds/I"); 
   
    tree->Branch("seed_pt",&seed_pt, "seed_pt[10]/D");
    tree->Branch("seed_eta",&seed_eta, "seed_eta[10]/D");
    tree->Branch("seed_phi",&seed_phi, "seed_phi[10]/D");
    tree->Branch("seed_mass",&seed_mass, "seed_mass[10]/D");    
    tree->Branch("seed_dz", &seed_dz, "seed_dz[10]/D");
    tree->Branch("seed_dxy", &seed_dxy, "seed_dxy[10]/D");
    tree->Branch("seed_3D_ip", &seed_3D_ip, "seed_3D_ip[10]/D");
    tree->Branch("seed_3D_sip", &seed_3D_sip, "seed_3D_sip[10]/D");
    tree->Branch("seed_2D_ip", &seed_2D_ip, "seed_2D_ip[10]/D");
    tree->Branch("seed_2D_sip", &seed_2D_sip, "seed_2D_sip[10]/D");    
    tree->Branch("seed_3D_signedIp", &seed_3D_signedIp, "seed_3D_signedIp[10]/D");
    tree->Branch("seed_3D_signedSip", &seed_3D_signedSip, "seed_3D_signedSip[10]/D");
    tree->Branch("seed_2D_signedIp", &seed_2D_signedIp, "seed_2D_signedIp[10]/D");
    tree->Branch("seed_2D_signedSip", &seed_2D_signedSip, "seed_2D_signedSip[10]/D");
    tree->Branch("seed_3D_TrackProbability", &seed_3D_TrackProbability, "seed_3D_TrackProbability[10]/D");
    tree->Branch("seed_2D_TrackProbability", &seed_2D_TrackProbability, "seed_2D_TrackProbability[10]/D");    
    tree->Branch("seed_chi2reduced",&seed_chi2reduced, "seed_chi2reduced[10]/D");
    tree->Branch("seed_nPixelHits",&seed_nPixelHits, "seed_nPixelHits[10]/D");
    tree->Branch("seed_nHits",&seed_nHits, "seed_nHits[10]/D");
    tree->Branch("seed_jetAxisDistance",&seed_jetAxisDistance, "seed_jetAxisDistance[10]/D");
    tree->Branch("seed_jetAxisDlength",&seed_jetAxisDlength, "seed_jetAxisDlength[10]/D");    
    
    //saving also some transformation of the variables
    tree->Branch("seed_ptLog", &seed_ptLog, "seed_ptLog[10]/D");
    tree->Branch("seed_dzSymLog", &seed_dzSymLog, "seed_dzSymLog[10]/D");
    tree->Branch("seed_dxySymLog", &seed_dxySymLog, "seed_dxySymLog[10]/D");
    tree->Branch("seed_3D_ipLog", &seed_3D_ipLog, "seed_3D_ipLog[10]/D");
    tree->Branch("seed_3D_sipLog", &seed_3D_sipLog, "seed_3D_sipLog[10]/D");
    tree->Branch("seed_2D_ipLog", &seed_2D_ipLog, "seed_2D_ipLog[10]/D");
    tree->Branch("seed_2D_sipLog", &seed_2D_sipLog, "seed_2D_sipLog[10]/D");    
    tree->Branch("seed_3D_signedIpSymLog", &seed_3D_signedIpSymLog, "seed_3D_signedIpSymLog[10]/D");
    tree->Branch("seed_3D_signedSipSymLog", &seed_3D_signedSipSymLog, "seed_3D_signedSipSymLog[10]/D");
    tree->Branch("seed_2D_signedIpSymLog", &seed_2D_signedIpSymLog, "seed_2D_signedIpSymLog[10]/D");
    tree->Branch("seed_2D_signedSipSymLog", &seed_2D_signedSipSymLog, "seed_2D_signedSipSymLog[10]/D");
    tree->Branch("seed_jetAxisDistanceLog",&seed_jetAxisDistanceLog, "seed_jetAxisDistanceLog[10]/D");
    tree->Branch("seed_jetAxisDlengthLog",&seed_jetAxisDlengthLog, "seed_jetAxisDlengthLog[10]/D");    

    tree->Branch("nearTracks_pt", &nearTracks_pt, "nearTracks_pt[200]/D");
    tree->Branch("nearTracks_eta", &nearTracks_eta, "nearTracks_eta[200]/D");
    tree->Branch("nearTracks_phi", &nearTracks_phi, "nearTracks_phi[200]/D");
    tree->Branch("nearTracks_mass", &nearTracks_mass, "nearTracks_mass[200]/D");
    tree->Branch("nearTracks_dz", &nearTracks_dz, "nearTracks_dz[200]/D");
    tree->Branch("nearTracks_dxy", &nearTracks_dxy, "nearTracks_dxy[200]/D");
    tree->Branch("nearTracks_3D_ip", &nearTracks_3D_ip, "nearTracks_3D_ip[200]/D");
    tree->Branch("nearTracks_3D_sip", &nearTracks_3D_sip, "nearTracks_3D_sip[200]/D");
    tree->Branch("nearTracks_2D_ip", &nearTracks_2D_ip, "nearTracks_2D_ip[200]/D");
    tree->Branch("nearTracks_2D_sip", &nearTracks_2D_sip, "nearTracks_2D_sip[200]/D");
    tree->Branch("nearTracks_PCAdist", &nearTracks_PCAdist, "nearTracks_PCAdist[200]/D");
    tree->Branch("nearTracks_PCAdsig", &nearTracks_PCAdsig, "nearTracks_PCAdsig[200]/D");    
    tree->Branch("nearTracks_PCAonSeed_x", &nearTracks_PCAonSeed_x, "nearTracks_PCAonSeed_x[200]/D");
    tree->Branch("nearTracks_PCAonSeed_y", &nearTracks_PCAonSeed_y, "nearTracks_PCAonSeed_y[200]/D");
    tree->Branch("nearTracks_PCAonSeed_z", &nearTracks_PCAonSeed_z, "nearTracks_PCAonSeed_z[200]/D");
    tree->Branch("nearTracks_PCAonSeed_xerr", &nearTracks_PCAonSeed_xerr, "nearTracks_PCAonSeed_xerr[200]/D");
    tree->Branch("nearTracks_PCAonSeed_yerr", &nearTracks_PCAonSeed_yerr, "nearTracks_PCAonSeed_yerr[200]/D");
    tree->Branch("nearTracks_PCAonSeed_zerr", &nearTracks_PCAonSeed_zerr, "nearTracks_PCAonSeed_zerr[200]/D");
    tree->Branch("nearTracks_PCAonTrack_x", &nearTracks_PCAonTrack_x, "nearTracks_PCAonTrack_x[200]/D");
    tree->Branch("nearTracks_PCAonTrack_y", &nearTracks_PCAonTrack_y, "nearTracks_PCAonTrack_y[200]/D");
    tree->Branch("nearTracks_PCAonTrack_z", &nearTracks_PCAonTrack_z, "nearTracks_PCAonTrack_z[200]/D");
    tree->Branch("nearTracks_PCAonTrack_xerr", &nearTracks_PCAonTrack_xerr, "nearTracks_PCAonTrack_xerr[200]/D");
    tree->Branch("nearTracks_PCAonTrack_yerr", &nearTracks_PCAonTrack_yerr, "nearTracks_PCAonTrack_yerr[200]/D");
    tree->Branch("nearTracks_PCAonTrack_zerr", &nearTracks_PCAonTrack_zerr, "nearTracks_PCAonTrack_zerr[200]/D"); 
    tree->Branch("nearTracks_dotprodTrack", &nearTracks_dotprodTrack, "nearTracks_dotprodTrack[200]/D");
    tree->Branch("nearTracks_dotprodSeed", &nearTracks_dotprodSeed, "nearTracks_dotprodSeed[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeed2D", &nearTracks_dotprodTrackSeed2D, "nearTracks_dotprodTrackSeed2D[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeed3D", &nearTracks_dotprodTrackSeed3D, "nearTracks_dotprodTrackSeed3D[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeedVectors2D", &nearTracks_dotprodTrackSeedVectors2D, "nearTracks_dotprodTrackSeedVectors2D[200]/D");
    tree->Branch("nearTracks_dotprodTrackSeedVectors3D", &nearTracks_dotprodTrackSeedVectors3D, "nearTracks_dotprodTrackSeedVectors3D[200]/D");    
    tree->Branch("nearTracks_PCAonSeed_pvd", &nearTracks_PCAonSeed_pvd, "nearTracks_PCAonSeed_pvd[200]/D");
    tree->Branch("nearTracks_PCAonTrack_pvd", &nearTracks_PCAonTrack_pvd, "nearTracks_PCAonTrack_pvd[200]/D");
    tree->Branch("nearTracks_PCAjetAxis_dist",&nearTracks_PCAjetAxis_dist,"nearTracks_PCAjetAxis_dist[200]/D");
    tree->Branch("nearTracks_PCAjetMomenta_dotprod",&nearTracks_PCAjetMomenta_dotprod,"nearTracks_PCAjetMomenta_dotprod[200]/D");
    tree->Branch("nearTracks_PCAjetDirs_DEta",&nearTracks_PCAjetDirs_DEta,"nearTracks_PCAjetDirs_DEta[200]/D");
    tree->Branch("nearTracks_PCAjetDirs_DPhi",&nearTracks_PCAjetDirs_DPhi,"nearTracks_PCAjetDirs_DPhi[200]/D");
    
    
    //saving also some transformation of the variables
    tree->Branch("nearTracks_ptLog",&nearTracks_ptLog,"nearTracks_ptLog[200]/D");
    tree->Branch("nearTracks_dzLog",&nearTracks_dzLog,"nearTracks_dzLog[200]/D");
    tree->Branch("nearTracks_dxyLog",&nearTracks_dxyLog,"nearTracks_dxyLog[200]/D");
    tree->Branch("nearTracks_3D_ipLog",&nearTracks_3D_ipLog,"nearTracks_3D_ipLog[200]/D");
    tree->Branch("nearTracks_3D_sipLog",&nearTracks_3D_sipLog,"nearTracks_3D_sipLog[200]/D");
    tree->Branch("nearTracks_2D_ipLog",&nearTracks_2D_ipLog,"nearTracks_2D_ipLog[200]/D");
    tree->Branch("nearTracks_2D_sipLog",&nearTracks_2D_sipLog,"nearTracks_2D_sipLog[200]/D");
    tree->Branch("nearTracks_PCAdistLog",&nearTracks_PCAdistLog,"nearTracks_PCAdistLog[200]/D");
    tree->Branch("nearTracks_PCAdsigLog",&nearTracks_PCAdsigLog,"nearTracks_PCAdsigLog[200]/D");  
    tree->Branch("nearTracks_PCAonSeed_pvdLog",&nearTracks_PCAonSeed_pvdLog,"nearTracks_PCAonSeed_pvdLog[200]/D");
    tree->Branch("nearTracks_PCAonTrack_pvdLog",&nearTracks_PCAonTrack_pvdLog,"nearTracks_PCAonTrack_pvdLog[200]/D");
    tree->Branch("nearTracks_PCAjetAxis_distLog",&nearTracks_PCAjetAxis_distLog,"nearTracks_PCAjetAxis_distLog[200]/D");
    tree->Branch("nearTracks_PCAjetDirs_DEtaLog",&nearTracks_PCAjetDirs_DEtaLog,"nearTracks_PCAjetDirs_DEtaLog[200]/D");
    

    
    CandidateToken = consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    token_primaryVertex = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    token_beamSpot = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    jetsToken = consumes<std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
    rhoToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));

}

Analyzer_MINIAOD_new::~Analyzer_MINIAOD_new()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Analyzer_MINIAOD_new::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    
    evt=iEvent.id().event();
    lumi=iEvent.id().luminosityBlock();
    run=iEvent.id().run();    
    
    if (m_computeProbabilities) checkEventSetup(iSetup);
    
    //1st way to get objects
    edm::Handle<std::vector<pat::Jet> > jetsCollection;
    iEvent.getByToken(jetsToken, jetsCollection);
    std::vector<pat::Jet>  ak4jets = *jetsCollection.product();    
    
    //2nd way to get objects
    edm::Handle<edm::View<pat::PackedCandidate> > tracks;
    iEvent.getByToken(CandidateToken, tracks);

    edm::Handle<reco::VertexCollection > primaryVertices;
    iEvent.getByToken(token_primaryVertex, primaryVertices);    
    
    edm::Handle<double> rhoH;
    iEvent.getByToken(rhoToken,rhoH);
    rho = *rhoH;
    
    if(primaryVertices->size()!=0){
        
        nPVs=primaryVertices->size();
        
        const reco::Vertex &pv = (*primaryVertices)[0];
        GlobalPoint pvp(pv.x(),pv.y(),pv.z());
        //        //std::cout << pv.x() << std::endl;
        
        pv_x=pv.x();
        pv_y=pv.y();
        pv_z=pv.z();
        
        edm::ESHandle<TransientTrackBuilder> trackBuilder;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

        edm::Handle<reco::BeamSpot> beamSpot;
        iEvent.getByToken(token_beamSpot,beamSpot);
        
        
        //std::cout << "Transient tracks" << std::endl;
        std::vector<reco::TransientTrack> selectedTracks;
        
        //std::cout << "Masses" << std::endl;
        std::vector<float> masses;
        
        //if primary: build transient tracks form packedCandidates
        
        for(typename edm::View<pat::PackedCandidate>::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
        
            unsigned int k=track - tracks->begin();         
 
            if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5) {
                
                if (std::fabs(pv_z-trackBuilder->build(tracks->ptrAt(k)).track().vz())<0.5){
                    selectedTracks.push_back(trackBuilder->build(tracks->ptrAt(k)));
                    masses.push_back(tracks->ptrAt(k)->mass());
                    
                }                
            }            
        }
        
        
        
        for (std::vector<pat::Jet>::const_iterator iter = ak4jets.begin(); iter != ak4jets.end(); ++iter) {
            
            
            
            if (iter->pt()<20) continue;
           
            GlobalVector direction(iter->px(), iter->py(), iter->pz());
            
            for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
                
                float angular_distance=std::sqrt(std::pow(iter->eta()-it->track().eta(),2) + std::pow(deltaPhi(iter->phi(),it->track().phi()),2) );    
                
                if (angular_distance<0.4){
                
                    int index = it - selectedTracks.begin();
                    //std::cout<<"both simultaneously and independently?  "<<it - selectedTracks.begin()<<std::endl;
                    //std::cout<<"both simultaneously and independently?  "<<it->track().eta()<<" "<<it->track().pt()<<std::endl;
                    
                    std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*it, pv);        
                    std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*it, pv);
                    std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, pv);                   
                    TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
                    float length=999.99;
                    if (closest.isValid()) length=(closest.globalPosition() - pvp).mag(); 
                    
                    if(ip.first && ip.second.value() >= 0.0 && ip.second.significance() >= 1.0 && ip.second.value() <= max3DIPValue && ip.second.significance() <= max3DIPSignificance && 
                        it->track().normalizedChi2()<5. && std::fabs(it->track().dxy(pv.position())) < 2 && std::fabs(it->track().dz(pv.position())) < 17  && jet_dist.second.value()<0.07 && length<5. ) {
                    
                        std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*it,direction, pv); 
                    
                        nearTracks.clear();
                        //std::cout<< "nearTracks  " << nearTracks.size() << " after clear " << std::endl;   

                        trackToTrack(selectedTracks, it, *iter, pv, nearTracks, masses, index);
                    
                        std::sort (nearTracks.begin(), nearTracks.end(), sortfunction2());
                        nearTracks.resize(20);
                        
                        //std::cout<< "nearTracks " << nearTracks.size() << " after all " << std::endl;  
                        
                        SortedSeedsMap.insert(std::make_pair(-ipSigned.second.significance(), std::make_pair(&(*it), nearTracks)));
                        
                    }
                        
                }
                
            }
            
            //FILL EVERYTHING 
            std::fill_n(seed_pt, 10, 0.);
            std::fill_n(seed_eta, 10, 0.);
            std::fill_n(seed_phi, 10, 0.);
            std::fill_n(seed_mass, 10, 0.);
            std::fill_n(seed_dz, 10, 0.);
            std::fill_n(seed_dxy, 10, 0.);
            std::fill_n(seed_3D_ip, 10, 0.);
            std::fill_n(seed_3D_sip, 10, 0.);
            std::fill_n(seed_2D_ip, 10, 0.);
            std::fill_n(seed_2D_sip, 10, 0.);
            std::fill_n(seed_3D_signedIp, 10, 0.);
            std::fill_n(seed_3D_signedSip, 10, 0.);
            std::fill_n(seed_2D_signedIp, 10, 0.);
            std::fill_n(seed_2D_signedSip, 10, 0.);
            std::fill_n(seed_chi2reduced, 10, 0.);
            std::fill_n(seed_nPixelHits, 10, 0.);
            std::fill_n(seed_nHits, 10, 0.);
            std::fill_n(seed_jetAxisDistance, 10, 0.);
            std::fill_n(seed_jetAxisDlength, 10, 0.);   
            std::fill_n(seed_3D_TrackProbability, 10, 0.);
            std::fill_n(seed_2D_TrackProbability, 10, 0.);            
            
            std::fill_n(seed_ptLog, 10, 0.);
            std::fill_n(seed_dzSymLog, 10, 0.);
            std::fill_n(seed_dxySymLog, 10, 0.);
            std::fill_n(seed_3D_ipLog, 10, 0.);
            std::fill_n(seed_3D_sipLog, 10, 0.);
            std::fill_n(seed_2D_ipLog, 10, 0.);
            std::fill_n(seed_2D_sipLog, 10, 0.);
            std::fill_n(seed_3D_signedIpSymLog, 10, 0.);
            std::fill_n(seed_3D_signedSipSymLog, 10, 0.);
            std::fill_n(seed_2D_signedIpSymLog, 10, 0.);
            std::fill_n(seed_2D_signedSipSymLog, 10, 0.);
            std::fill_n(seed_jetAxisDistanceLog, 10, 0.);
            std::fill_n(seed_jetAxisDlengthLog, 10, 0.);
            
            std::fill_n(nearTracks_pt, 200, 0.);
            std::fill_n(nearTracks_eta, 200, 0.);
            std::fill_n(nearTracks_phi, 200, 0.);
            std::fill_n(nearTracks_dz, 200, 0.);
            std::fill_n(nearTracks_dxy, 200, 0.);
            std::fill_n(nearTracks_mass, 200, 0.);
            std::fill_n(nearTracks_3D_ip, 200, 0.);
            std::fill_n(nearTracks_3D_sip, 200, 0.);
            std::fill_n(nearTracks_2D_ip, 200, 0.);
            std::fill_n(nearTracks_2D_sip, 200, 0.);
            std::fill_n(nearTracks_PCAdist, 200, 0.);
            std::fill_n(nearTracks_PCAdsig, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_x, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_y, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_z, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_xerr, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_yerr, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_zerr, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_x, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_y, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_z, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_xerr, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_yerr, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_zerr, 200, 0.);
            std::fill_n(nearTracks_dotprodTrack, 200, 0.);
            std::fill_n(nearTracks_dotprodSeed, 200, 0.);
            std::fill_n(nearTracks_dotprodTrackSeed2D, 200, 0.);
            std::fill_n(nearTracks_dotprodTrackSeed3D, 200, 0.);
            std::fill_n(nearTracks_dotprodTrackSeedVectors2D, 200, 0.);
            std::fill_n(nearTracks_dotprodTrackSeedVectors3D, 200, 0.);       
            std::fill_n(nearTracks_PCAonSeed_pvd, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_pvd, 200, 0.);
            std::fill_n(nearTracks_PCAjetAxis_dist, 200, 0.);
            std::fill_n(nearTracks_PCAjetMomenta_dotprod, 200, 0.);
            std::fill_n(nearTracks_PCAjetDirs_DEta, 200, 0.);
            std::fill_n(nearTracks_PCAjetDirs_DPhi, 200, 0.);
            
            std::fill_n(nearTracks_ptLog, 200, 0.);
            std::fill_n(nearTracks_dzLog, 200, 0.);
            std::fill_n(nearTracks_dxyLog, 200, 0.);
            std::fill_n(nearTracks_3D_ipLog, 200, 0.);
            std::fill_n(nearTracks_3D_sipLog, 200, 0.);
            std::fill_n(nearTracks_2D_ipLog, 200, 0.);
            std::fill_n(nearTracks_2D_sipLog, 200, 0.);
            std::fill_n(nearTracks_PCAdistLog, 200, 0.);
            std::fill_n(nearTracks_PCAdsigLog, 200, 0.);
            std::fill_n(nearTracks_PCAonSeed_pvdLog, 200, 0.);
            std::fill_n(nearTracks_PCAonTrack_pvdLog, 200, 0.);
            std::fill_n(nearTracks_PCAjetAxis_distLog, 200, 0.);
            std::fill_n(nearTracks_PCAjetDirs_DEtaLog, 200, 0.);
            
            
            unsigned int seeds_max_counter=0;
            for(std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2>>>::const_iterator im = SortedSeedsMap.begin(); im != SortedSeedsMap.end(); im++){

                if(seeds_max_counter>=10) break;

                std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*im->second.first,direction, pv);        
                std::pair<bool,Measurement1D> ip2dSigned = IPTools::signedTransverseImpactParameter(*im->second.first,direction, pv);  
                std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*im->second.first, pv);        
                std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*im->second.first, pv);
                std::pair<double, Measurement1D> jet_distance =IPTools::jetTrackDistance(*im->second.first, direction, pv);         
                TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(im->second.first->impactPointState(),pv, direction,im->second.first->field());
                
                
                if (m_computeProbabilities) {
                    
                    std::pair<bool,double> probability;
                    
                    //probability with 3D ip
                    
                    probability = m_probabilityEstimator->probability(0, 0,ip.second.significance(),im->second.first->track(),*iter,pv);
                    double prob3D=(probability.first ? probability.second : -1.);
//                     std::cout<<prob3D<<std::endl;

                    //probability with 2D ip
                    
                    probability = m_probabilityEstimator->probability(0, 1,ip2d.second.significance(),im->second.first->track(),*iter,pv);
                    double prob2D=(probability.first ? probability.second : -1.);
//                     std::cout<<prob2D<<std::endl;
                                        
                    seed_3D_TrackProbability[seeds_max_counter]=prob3D;
                    seed_2D_TrackProbability[seeds_max_counter]=prob2D;

                  
                }

                seed_pt[seeds_max_counter]=im->second.first->track().pt();
                seed_eta[seeds_max_counter]=im->second.first->track().eta();
                seed_phi[seeds_max_counter]=im->second.first->track().phi();
                seed_mass[seeds_max_counter]=im->second.second.at(0).seedMass;
                seed_dz[seeds_max_counter]=im->second.first->track().dz(pv.position());
                seed_dxy[seeds_max_counter]=im->second.first->track().dxy(pv.position());
                seed_3D_ip[seeds_max_counter]=ip.second.value();
                seed_3D_sip[seeds_max_counter]=ip.second.significance();
                seed_2D_ip[seeds_max_counter]=ip2d.second.value();
                seed_2D_sip[seeds_max_counter]=ip2d.second.significance();
                seed_3D_signedIp[seeds_max_counter]=ipSigned.second.value();
                seed_3D_signedSip[seeds_max_counter]=ipSigned.second.significance();
                seed_2D_signedIp[seeds_max_counter]=ip2dSigned.second.value();
                seed_2D_signedSip[seeds_max_counter]=ip2dSigned.second.significance();
                seed_chi2reduced[seeds_max_counter]=im->second.first->track().normalizedChi2();
                seed_nPixelHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidPixelHits();
                seed_nHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidHits();
                
                seed_jetAxisDistance[seeds_max_counter]=std::fabs(jet_distance.second.value());                
                if (closest.isValid()) seed_jetAxisDlength[seeds_max_counter]=(closest.globalPosition() - pvp).mag(); 
                else seed_jetAxisDlength[seeds_max_counter]= -99;  
                
                seed_ptLog[seeds_max_counter]=log(im->second.first->track().pt());
                seed_dzSymLog[seeds_max_counter]=(15+log(fabs(im->second.first->track().dz(pv.position()))))*((im->second.first->track().dz(pv.position()) < 0) ? -1 : (im->second.first->track().dz(pv.position()) > 0));
                seed_dxySymLog[seeds_max_counter]=(15+log(fabs(im->second.first->track().dxy(pv.position()))))*((im->second.first->track().dxy(pv.position()) < 0) ? -1 : (im->second.first->track().dxy(pv.position()) > 0));
                seed_3D_ipLog[seeds_max_counter]=log(ip.second.value());
                seed_3D_sipLog[seeds_max_counter]=log(ip.second.significance());
                seed_2D_ipLog[seeds_max_counter]=log(ip2d.second.value());
                seed_2D_sipLog[seeds_max_counter]=log(ip2d.second.significance());
                seed_3D_signedIpSymLog[seeds_max_counter]=(8+log(fabs(ipSigned.second.value())))*((ipSigned.second.value() < 0) ? -1 : (ipSigned.second.value() > 0));
                seed_3D_signedSipSymLog[seeds_max_counter]=(2+log(fabs(ipSigned.second.significance())))*((ipSigned.second.significance() < 0) ? -1 : (ipSigned.second.significance() > 0));
                seed_2D_signedIpSymLog[seeds_max_counter]=(15+log(fabs(ip2dSigned.second.value())))*((ip2dSigned.second.value() < 0) ? -1 : (ip2dSigned.second.value() > 0));
                seed_2D_signedSipSymLog[seeds_max_counter]=(10+log(fabs(ip2dSigned.second.significance())))*((ip2dSigned.second.significance() < 0) ? -1 : (ip2dSigned.second.significance() > 0));
                seed_jetAxisDistanceLog[seeds_max_counter]=log(std::fabs(jet_distance.second.value()));
                seed_jetAxisDlengthLog[seeds_max_counter]=log(seed_jetAxisDlength[seeds_max_counter]);   
                
                
                for(unsigned int i=0; i< im->second.second.size(); i++) {
                    if((seeds_max_counter*20+i)>=200) break;
                    
                    //std::cout << im->second.second.at(i).dist << "  dist  "; 
                    //std::cout << im->second.second.at(i).dsig << "  dsig ";
                    nearTracks_pt[seeds_max_counter*20+i]=im->second.second.at(i).pt;
                    nearTracks_eta[seeds_max_counter*20+i]=im->second.second.at(i).eta;
                    nearTracks_phi[seeds_max_counter*20+i]=im->second.second.at(i).phi;
                    nearTracks_dz[seeds_max_counter*20+i]=im->second.second.at(i).dz;
                    nearTracks_dxy[seeds_max_counter*20+i]=im->second.second.at(i).dxy;
                    nearTracks_mass[seeds_max_counter*20+i]=im->second.second.at(i).mass;
                    nearTracks_3D_ip[seeds_max_counter*20+i]=im->second.second.at(i).t3Dip;
                    nearTracks_3D_sip[seeds_max_counter*20+i]=im->second.second.at(i).t3Dsip;
                    nearTracks_2D_ip[seeds_max_counter*20+i]=im->second.second.at(i).t2Dip;
                    nearTracks_2D_sip[seeds_max_counter*20+i]=im->second.second.at(i).t2Dsip;
                    nearTracks_PCAdist[seeds_max_counter*20+i]=im->second.second.at(i).dist;
                    nearTracks_PCAdsig[seeds_max_counter*20+i]=im->second.second.at(i).dsig;
                    nearTracks_PCAonSeed_x[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sx;
                    nearTracks_PCAonSeed_y[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sy;
                    nearTracks_PCAonSeed_z[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sz;
                    nearTracks_PCAonSeed_xerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_sxerr;
                    nearTracks_PCAonSeed_yerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_syerr;
                    nearTracks_PCAonSeed_zerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_szerr;
                    nearTracks_PCAonTrack_x[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tx;
                    nearTracks_PCAonTrack_y[seeds_max_counter*20+i]=im->second.second.at(i).PCA_ty;
                    nearTracks_PCAonTrack_z[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tz;
                    nearTracks_PCAonTrack_xerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_txerr;
                    nearTracks_PCAonTrack_yerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tyerr;
                    nearTracks_PCAonTrack_zerr[seeds_max_counter*20+i]=im->second.second.at(i).PCA_tzerr;
                    nearTracks_dotprodTrack[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrack;
                    nearTracks_dotprodSeed[seeds_max_counter*20+i]=im->second.second.at(i).dotprodSeed;
                    nearTracks_dotprodTrackSeed2D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeed2D;
                    nearTracks_dotprodTrackSeed3D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeed3D;
                    nearTracks_dotprodTrackSeedVectors2D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeedVectors2D;
                    nearTracks_dotprodTrackSeedVectors3D[seeds_max_counter*20+i]=im->second.second.at(i).dotprodTrackSeedVectors3D;

                    nearTracks_PCAonSeed_pvd[seeds_max_counter*20+i]=im->second.second.at(i).seedPCA_pv;
                    nearTracks_PCAonTrack_pvd[seeds_max_counter*20+i]=im->second.second.at(i).trackPCA_pv;
                    //std::cout <<  "after lambda  " << seeds_max_counter*20+i<< std::endl;
                    //std::cout << im->second.second.at(i).PCA_JetAxis_distance << "  " << im->second.second.at(i).PCAPair_Jet_dotprod << " "                    << im->second.second.at(i).PCAAxis_JetAxis_DEta <<" " <<im->second.second.at(i).PCAAxis_JetAxis_DPhi << std::endl;

                    nearTracks_PCAjetAxis_dist[seeds_max_counter*20+i]=im->second.second.at(i).PCA_JetAxis_distance;
                    nearTracks_PCAjetMomenta_dotprod[seeds_max_counter*20+i]=im->second.second.at(i).PCAPair_Jet_dotprod;

                    nearTracks_PCAjetDirs_DEta[seeds_max_counter*20+i]=im->second.second.at(i).PCAAxis_JetAxis_DEta;
                    nearTracks_PCAjetDirs_DPhi[seeds_max_counter*20+i]=im->second.second.at(i).PCAAxis_JetAxis_DPhi;
                    
                    nearTracks_ptLog[seeds_max_counter*20+i]=log(im->second.second.at(i).pt);
                    nearTracks_dxyLog[seeds_max_counter*20+i]=(16+log(fabs(im->second.second.at(i).dxy)))*((im->second.second.at(i).dxy < 0) ? -1 : (im->second.second.at(i).dxy > 0));
                    nearTracks_dzLog[seeds_max_counter*20+i]=(16+log(fabs(im->second.second.at(i).dz)))*((im->second.second.at(i).dz < 0) ? -1 : (im->second.second.at(i).dz > 0));
                    nearTracks_3D_ipLog[seeds_max_counter*20+i]=log(im->second.second.at(i).t3Dip);
                    nearTracks_3D_sipLog[seeds_max_counter*20+i]=log(im->second.second.at(i).t3Dsip);
                    nearTracks_2D_ipLog[seeds_max_counter*20+i]=log(im->second.second.at(i).t2Dip);
                    nearTracks_2D_sipLog[seeds_max_counter*20+i]=log(im->second.second.at(i).t2Dsip);
                    nearTracks_PCAdistLog[seeds_max_counter*20+i]=log(im->second.second.at(i).dist);
                    nearTracks_PCAdsigLog[seeds_max_counter*20+i]=log(im->second.second.at(i).dsig);  
                    nearTracks_PCAonSeed_pvdLog[seeds_max_counter*20+i]=log(im->second.second.at(i).seedPCA_pv);
                    nearTracks_PCAonTrack_pvdLog[seeds_max_counter*20+i]=log(im->second.second.at(i).trackPCA_pv);
                    nearTracks_PCAjetAxis_distLog[seeds_max_counter*20+i]=log(im->second.second.at(i).PCA_JetAxis_distance);
                    nearTracks_PCAjetDirs_DEtaLog[seeds_max_counter*20+i]=log(im->second.second.at(i).PCAAxis_JetAxis_DEta);
                    
                }
            
            seeds_max_counter++;
                
            }
            
            jetpt=iter->pt();
            jetptLOG=log(iter->pt());
            jeteta=iter->eta();
            jetphi=iter->phi();
            jetmass=iter->mass();
            jetCSVv2=iter->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            jetCMVA=iter->bDiscriminator("pfCombinedMVAV2BJetTags");
            jetTCHE=iter->bDiscriminator("pfTrackCountingHighEffBJetTags");
            jetSSVHE=iter->bDiscriminator("pfSimpleInclusiveSecondaryVertexHighEffBJetTags");
            
            jetDeepCSV=iter->bDiscriminator("pfDeepCSVJetTags:probbb")+iter->bDiscriminator("pfDeepCSVJetTags:probb");
            jetDeepFlavour=iter->bDiscriminator("pfDeepFlavourJetTags:probbb")+iter->bDiscriminator("pfDeepFlavourJetTags:probb");
            
            jetflavour=jet_flavour_function(*iter, 1);      
            jetNseeds=std::min((int)SortedSeedsMap.size(),10);

            
            tree->Fill();
            
            //CLEAR THE JET ATTRIBUTED TRACKS AND CLUSTERS
            SortedSeedsMap.clear(); 
            
            
        }
        
        
       
    }
    
    

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
//std::cout <<  " fine evento  " << std::endl;


}

double Analyzer_MINIAOD_new::deltaPhi(double phi1, double phi2) { 
    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI) result -= 2*PI;
    while (result <= -PI) result += 2*PI;
    return result;
    
}


// Jet Flavour
// BB 5
// B  5
// C  4
// UD 1
// S  3
// G  21
// NAN 0

int Analyzer_MINIAOD_new::jet_flavour_function(const pat::Jet& jet, bool usePhysForLightAndUndefined) {
    int hflav = abs(jet.hadronFlavour());
    int pflav = abs(jet.partonFlavour());
    int physflav = 0;
    if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
    std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
    std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();

    if(hflav == 5) { //B jet
        if(nbs > 1) return 5;
        else if(nbs == 1) {
            return 5;
        }
        else {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return 21;
                else if(physflav == 3) return 3;
                else if(physflav == 2 || physflav ==1) return 1;
                else return 0;
            }
            else return 0;
        }
    }
    else if(hflav == 4) { //C jet
        return 4;
    }
    else { //not a heavy jet
        if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
            if(usePhysForLightAndUndefined){
                if(physflav == 21) return 21;
                else if(physflav == 3) return 3;
                else if(physflav == 2 || physflav ==1) return 1;
                else return 0;
            }
            else return 0;
        }
        else if(usePhysForLightAndUndefined){
            if(physflav == 21) return 21;
            else if(physflav == 3) return 3;
            else if(physflav == 2 || physflav ==1) return 1;
            else return 0;
        }
        else {
            if(pflav == 21) return 21;
            else if(pflav == 3) return 3;
            else if(pflav == 2 || pflav ==1) return 1;
            else return 0;
        }
    }
}



void Analyzer_MINIAOD_new::trackToTrack(std::vector<reco::TransientTrack> selTracks, std::vector<reco::TransientTrack>::const_iterator it, 
                                        const pat::Jet &jet, const reco::Vertex &pv, std::vector<trackVars2>& nearTracks, std::vector<float> masses, int index ){
   
    
    
    //std::cout<<"in clustering function"<<std::endl;
    //std::cout<<"clustering?  "<<it->track().eta()<<" "<<it->track().pt()<<std::endl;
    //std::cout<<"in clustering function 31  "<<it-selTracks.begin()<<std::endl;
    //std::cout<<"in clustering function 31  "<<index<<std::endl;
    
    trackVars2 myTrack;
    
    GlobalPoint pvp(pv.x(),pv.y(),pv.z());    
    GlobalVector direction(jet.px(), jet.py(), jet.pz());
    
    for(std::vector<reco::TransientTrack>::const_iterator tt = selTracks.begin();tt!=selTracks.end(); ++tt )   {
    
        VertexDistance3D distanceComputer;
        TwoTrackMinimumDistance dist;
        
        if(*tt==*it) continue;
        if(std::fabs(pv.z()-tt->track().vz())>0.1) continue;
        
//         //std::cout<<"in clustering function tt"<<std::endl;
        
        if(dist.calculate(tt->impactPointState(),it->impactPointState())) {            
            
                       
            GlobalPoint ttPoint          = dist.points().first;
            GlobalError ttPointErr       = tt->impactPointState().cartesianError().position();
            GlobalPoint seedPosition     = dist.points().second;
            GlobalError seedPositionErr  = it->impactPointState().cartesianError().position();
            
            Measurement1D m = distanceComputer.distance(VertexState(seedPosition,seedPositionErr), VertexState(ttPoint, ttPointErr));

            GlobalPoint cp(dist.crossingPoint()); 
            GlobalVector PairMomentum(it->track().px()+tt->track().px(), it->track().py()+tt->track().py(), it->track().pz()+tt->track().pz());
            GlobalVector  PCA_pv(cp-pvp);

            float PCAseedFromPV =  (dist.points().second-pvp).mag();
            float PCAtrackFromPV =  (dist.points().first-pvp).mag();               
            float distance = dist.distance();            
            
            GlobalVector trackDir2D(tt->impactPointState().globalDirection().x(),tt->impactPointState().globalDirection().y(),0.); 
            GlobalVector seedDir2D(it->impactPointState().globalDirection().x(),it->impactPointState().globalDirection().y(),0.); 
            GlobalVector trackPCADir2D(dist.points().first.x()-pvp.x(),dist.points().first.y()-pvp.y(),0.); 
            GlobalVector seedPCADir2D(dist.points().second.x()-pvp.x(),dist.points().second.y()-pvp.y(),0.); 

            float dotprodTrack = (dist.points().first-pvp).unit().dot(tt->impactPointState().globalDirection().unit());
            float dotprodSeed = (dist.points().second-pvp).unit().dot(it->impactPointState().globalDirection().unit());                    
            float dotprodTrackSeed2D = trackDir2D.unit().dot(seedDir2D.unit());
            float dotprodTrackSeed3D = it->impactPointState().globalDirection().unit().dot(tt->impactPointState().globalDirection().unit());
            float dotprodTrackSeed2DV = trackPCADir2D.unit().dot(seedPCADir2D.unit());
            float dotprodTrackSeed3DV = (dist.points().second-pvp).unit().dot((dist.points().first-pvp).unit());

            std::pair<bool,Measurement1D> t_ip = IPTools::absoluteImpactParameter3D(*tt,pv);        
            std::pair<bool,Measurement1D> t_ip2d = IPTools::absoluteTransverseImpactParameter(*tt,pv);
            
            Line::PositionType pos(pvp);
            Line::DirectionType dir(direction);
            Line::DirectionType pairMomentumDir(PairMomentum);
            Line jetLine(pos,dir);   
            Line PCAMomentumLine(cp,pairMomentumDir);           
            
            float PCA_JetAxis_dist=jetLine.distance(cp).mag();
            float dotprodMomenta=PairMomentum.unit().dot(direction.unit());
            float dEta=std::fabs(PCA_pv.eta()-jet.eta());
            float dPhi=std::fabs(PCA_pv.phi()-jet.phi());   

            myTrack.set_values(tt->track().pt(), tt->track().eta(), tt->track().phi(),  tt->track().dz(pv.position()), tt->track().dxy(pv.position()), distance, m.significance(), 
                               seedPosition.x(), seedPosition.y(), seedPosition.z(), seedPositionErr.cxx(), seedPositionErr.cyy(), seedPositionErr.czz(), 
                               ttPoint.x(),  ttPoint.y(),  ttPoint.z(),  ttPointErr.cxx(),  ttPointErr.cyy(),  ttPointErr.czz(), 
                               dotprodTrack, dotprodSeed );    
            myTrack.set_index(-1);
            myTrack.set_distances(PCAseedFromPV, PCAtrackFromPV);
            myTrack.set_vars(masses[tt-selTracks.begin()],t_ip2d.second.value() , t_ip2d.second.significance(),
            t_ip.second.value() , t_ip.second.significance(), dotprodTrackSeed2D, dotprodTrackSeed3D, dotprodTrackSeed2DV, dotprodTrackSeed3DV );            
            myTrack.setSeedMass(masses[index]);   
            myTrack.set_JetAxisVars(PCA_JetAxis_dist,dotprodMomenta,dEta,dPhi);            
            nearTracks.push_back(myTrack);

            
        }
    }
}


void Analyzer_MINIAOD_new::checkEventSetup(const edm::EventSetup & iSetup) {
  
  using namespace edm;
  using namespace edm::eventsetup;

   const EventSetupRecord & re2D= iSetup.get<BTagTrackProbability2DRcd>();
   const EventSetupRecord & re3D= iSetup.get<BTagTrackProbability3DRcd>();
   unsigned long long cacheId2D= re2D.cacheIdentifier();
   unsigned long long cacheId3D= re3D.cacheIdentifier();

   if(cacheId2D!=m_calibrationCacheId2D || cacheId3D!=m_calibrationCacheId3D  )  //Calibration changed
   {
     //iSetup.get<BTagTrackProbabilityRcd>().get(calib);
     ESHandle<TrackProbabilityCalibration> calib2DHandle;
     iSetup.get<BTagTrackProbability2DRcd>().get(calib2DHandle);
     ESHandle<TrackProbabilityCalibration> calib3DHandle;
     iSetup.get<BTagTrackProbability3DRcd>().get(calib3DHandle);

     const TrackProbabilityCalibration *  ca2D= calib2DHandle.product();
     const TrackProbabilityCalibration *  ca3D= calib3DHandle.product();

     m_probabilityEstimator.reset(new HistogramProbabilityEstimator(ca3D,ca2D));

   }
   
   m_calibrationCacheId3D=cacheId3D;
   m_calibrationCacheId2D=cacheId2D;
   
}



// ------------ method called once each job just before starting event loop  ------------
void 
Analyzer_MINIAOD_new::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer_MINIAOD_new::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer_MINIAOD_new::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer_MINIAOD_new);
