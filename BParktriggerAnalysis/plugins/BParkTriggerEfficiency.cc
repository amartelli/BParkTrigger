#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TLorentzVector.h>



typedef edm::RefVector<std::vector<TrackingParticle>> TrackingParticleContainer;
typedef std::vector<TrackingParticle> TrackingParticleCollection;

typedef TrackingParticleRefVector::iterator tp_iterator;
typedef TrackingVertex::genv_iterator genv_iterator;
typedef TrackingVertex::g4v_iterator g4v_iterator;

typedef math::XYZPointD Point; 
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > XYZPointF;

edm::Service<TFileService> fs;

class BParkTriggerEfficiency : public edm::EDAnalyzer {
public:
  explicit BParkTriggerEfficiency(const edm::ParameterSet &conf);

  ~BParkTriggerEfficiency() override {}

  void analyze(const edm::Event &e, const edm::EventSetup &c) override;
  virtual void beginJob() override;
  virtual void endJob() override;

  int findGenBToKee(edm::Handle<edm::View<reco::GenParticle>> genPart,
		    int& B_bc, int& l1_bc, int& l2_bc, int& k_bc, 
		    TLorentzVector& Bgen, TLorentzVector& L1gen,
		    TLorentzVector& L2gen, TLorentzVector& Kgen,
		    Point& BdecayVtx);

private:
  edm::ParameterSet conf_;
  edm::EDGetToken trackingTruthTrk;
  edm::EDGetToken trackingTruthVtx;
  edm::EDGetToken genToken_;
  edm::EDGetTokenT<std::vector<int>> genParticleInts_;

  edm::EDGetTokenT<reco::TrackCollection> inputPixTrkToken_;
  edm::EDGetTokenT<reco::VertexCollection> inputPixVtxToken_;
  
  edm::EDGetTokenT<reco::SimToRecoCollection> trackSimToRecoAssociationToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;

  edm::EDGetTokenT<XYZPointF> genVtxPositionToken_;


  TH1F* h_dR_e1_pixel_Tracking;
  TH1F* h_dR_e2_pixel_Tracking;
  TH1F* h_dR_k_pixel_Tracking;


  TH1F* h_gen_pu[4];
  TH1F* h_gen_eta[4];
  TH1F* h_gen_phi[4];
  TH1F* h_gen_pt[4];

  TH1F* h_trkP_pu[4];
  TH1F* h_trkP_eta[4];
  TH1F* h_trkP_phi[4];
  TH1F* h_trkP_pt[4];
  //
  TH1F* h_trkP_Maxeta;
  TH1F* h_trkP_Minpt;

  TH1F* h_pixTrk_pu[4];
  TH1F* h_pixTrk_eta[4];
  TH1F* h_pixTrk_phi[4];
  TH1F* h_pixTrk_pt[4];

  TH1F* h_pixTrk_pR_pu[4];
  TH1F* h_pixTrk_pR_eta[4];
  TH1F* h_pixTrk_pR_phi[4];
  TH1F* h_pixTrk_pR_pt[4];
  //
  TH1F* h_pixTrk_pR_Maxeta;
  TH1F* h_pixTrk_pR_Minpt;


  TH1F* h_ratio_pixTrk_trkP_pu[4];
  TH1F* h_ratio_pixTrk_trkP_eta[4];
  TH1F* h_ratio_pixTrk_trkP_phi[4];
  TH1F* h_ratio_pixTrk_trkP_pt[4];
  //
  TH1F* h_ratio_pixTrk_trkP_Maxeta;
  TH1F* h_ratio_pixTrk_trkP_Minpt;


  TH1F* h_trkP_Bmass;
  TH1F* h_gen_Bmass;
  TH1F* h_pixTrk_Bmass;
  TH1F* h_pixTrk_pR_Bmass;
  TH1F* h_ratio_pixTrk_trkP_Bmass;

  TH1F* h_trkP_Bdxy;
  TH1F* h_pixTrk_pR_Bdxy;
  TH1F* h_ratio_pixTrk_trkP_Bdxy;


  int nTotEvents;
  int nEvents_Bdecay;
  int nEvents_trkP_matched;
  int nEvents_pixTrk_matched;
  int NVtxMatched_xy;
  int NVtxMatched;

  bool debug;
};


BParkTriggerEfficiency::BParkTriggerEfficiency(const edm::ParameterSet &conf):
  trackingTruthTrk(consumes<std::vector<TrackingParticle>>(conf.getParameter<edm::InputTag>("trackingTruth"))),
  trackingTruthVtx(consumes<std::vector<TrackingVertex>>(conf.getParameter<edm::InputTag>("trackingTruth"))),
  genToken_(consumes<edm::View<reco::GenParticle>>(conf.getParameter<edm::InputTag>("genParticles"))),
  genParticleInts_(consumes<std::vector<int>>(conf.getParameter<edm::InputTag>("genParticles"))),
  inputPixTrkToken_(consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("pixelTracks"))),
  inputPixVtxToken_(consumes<reco::VertexCollection>(conf.getParameter<edm::InputTag>("pixelVertex"))),
  trackSimToRecoAssociationToken_(consumes<reco::SimToRecoCollection>(conf.getParameter<edm::InputTag>("trackAssociation"))),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(conf.getUntrackedParameter<edm::InputTag>("puSummary"))),
  genVtxPositionToken_(consumes<XYZPointF>(conf.getParameter<edm::InputTag>("genPxyz0")))
{ 
  conf_ = conf; 
  debug = false;

  nTotEvents = 0;
  nEvents_Bdecay = 0;
  nEvents_trkP_matched = 0;  
  nEvents_pixTrk_matched = 0;
  
  NVtxMatched_xy = 0;
  NVtxMatched = 0;


  edm::Service<TFileService> fs;
  h_dR_e1_pixel_Tracking = fs->make<TH1F>("h_dR_e1_pixel_Tracking", "", 1000, 0., 15.);
  h_dR_e2_pixel_Tracking = fs->make<TH1F>("h_dR_e2_pixel_Tracking", "", 1000, 0., 15.);
  h_dR_k_pixel_Tracking = fs->make<TH1F>("h_dR_k_pixel_Tracking", "", 1000, 0., 15.);


  for(int ij=0; ij<4; ++ij){
    h_gen_pu[ij] = fs->make<TH1F>(Form("h_gen_pu_%d", ij), "", 100, 0, 100);
    h_gen_eta[ij] = fs->make<TH1F>(Form("h_gen_eta_%d", ij), "", 40, -5, 5);
    h_gen_phi[ij] = fs->make<TH1F>(Form("h_gen_phi_%d", ij), "", 40, -4, 4);
    h_gen_pt[ij] = fs->make<TH1F>(Form("h_gen_pt_%d", ij), "", 80, 0, 20);
    h_gen_pu[ij]->Sumw2();
    h_gen_eta[ij]->Sumw2();
    h_gen_phi[ij]->Sumw2();
    h_gen_pt[ij]->Sumw2();

    h_trkP_pu[ij] = fs->make<TH1F>(Form("h_trkP_pu_%d", ij), "", 100, 0, 100);
    h_trkP_eta[ij] = fs->make<TH1F>(Form("h_trkP_eta_%d", ij), "", 40, -5, 5);
    h_trkP_phi[ij] = fs->make<TH1F>(Form("h_trkP_phi_%d", ij), "", 40, -4, 4);
    h_trkP_pt[ij] = fs->make<TH1F>(Form("h_trkP_pt_%d", ij), "", 80, 0, 20);
    h_trkP_pu[ij]->Sumw2();
    h_trkP_eta[ij]->Sumw2();
    h_trkP_phi[ij]->Sumw2();
    h_trkP_pt[ij]->Sumw2();

    h_pixTrk_pu[ij] = fs->make<TH1F>(Form("h_pixTrk_pu_%d", ij), "", 100, 0, 100);
    h_pixTrk_eta[ij] = fs->make<TH1F>(Form("h_pixTrk_eta_%d", ij), "", 40, -5, 5);
    h_pixTrk_phi[ij] = fs->make<TH1F>(Form("h_pixTrk_phi_%d", ij), "", 40, -4, 4);
    h_pixTrk_pt[ij] = fs->make<TH1F>(Form("h_pixTrk_pt_%d", ij), "", 80, 0, 20);
    h_pixTrk_pu[ij]->Sumw2();
    h_pixTrk_eta[ij]->Sumw2();
    h_pixTrk_phi[ij]->Sumw2();
    h_pixTrk_pt[ij]->Sumw2();

    h_pixTrk_pR_pu[ij] = fs->make<TH1F>(Form("h_pixTrk_pR_pu_%d", ij), "", 100, 0, 100);
    h_pixTrk_pR_eta[ij] = fs->make<TH1F>(Form("h_pixTrk_pR_eta_%d", ij), "", 40, -5, 5);
    h_pixTrk_pR_phi[ij] = fs->make<TH1F>(Form("h_pixTrk_pR_phi_%d", ij), "", 40, -4, 4);
    h_pixTrk_pR_pt[ij] = fs->make<TH1F>(Form("h_pixTrk_pR_pt_%d", ij), "", 80, 0, 20);
    h_pixTrk_pR_pu[ij]->Sumw2();
    h_pixTrk_pR_eta[ij]->Sumw2();
    h_pixTrk_pR_phi[ij]->Sumw2();
    h_pixTrk_pR_pt[ij]->Sumw2();

    h_ratio_pixTrk_trkP_pu[ij] = fs->make<TH1F>(Form("h_ratio_pixTrk_trkP_pu_%d", ij), "", 100, 0, 100);
    h_ratio_pixTrk_trkP_eta[ij] = fs->make<TH1F>(Form("h_ratio_pixTrk_trkP_eta_%d", ij), "", 40, -5, 5);
    h_ratio_pixTrk_trkP_phi[ij] = fs->make<TH1F>(Form("h_ratio_pixTrk_trkP_phi_%d", ij), "", 40, -4, 4);
    h_ratio_pixTrk_trkP_pt[ij] = fs->make<TH1F>(Form("h_ratio_pixTrk_trkP_pt_%d", ij), "", 80, 0, 20);
    h_ratio_pixTrk_trkP_pu[ij]->Sumw2();
    h_ratio_pixTrk_trkP_eta[ij]->Sumw2();
    h_ratio_pixTrk_trkP_phi[ij]->Sumw2();
    h_ratio_pixTrk_trkP_pt[ij]->Sumw2();
  }
  //vs min pt max eta
  h_trkP_Maxeta = fs->make<TH1F>("h_trkP_Maxeta", "", 40, -5, 5);
  h_trkP_Minpt = fs->make<TH1F>("h_trkP_Minpt", "", 40, 0, 20);
  h_trkP_Maxeta->Sumw2();
  h_trkP_Minpt->Sumw2();
  //                
  h_pixTrk_pR_Maxeta = fs->make<TH1F>("h_pixTrk_pR_Maxeta", "", 40, -5, 5);
  h_pixTrk_pR_Minpt = fs->make<TH1F>("h_pixTrk_pR_Minpt", "", 40, 0, 20);
  h_pixTrk_pR_Maxeta->Sumw2();
  h_pixTrk_pR_Minpt->Sumw2();
  //
  h_ratio_pixTrk_trkP_Maxeta = fs->make<TH1F>("h_ratio_pixTrk_trkP_Maxeta", "", 40, -5, 5);
  h_ratio_pixTrk_trkP_Maxeta->Sumw2();
  h_ratio_pixTrk_trkP_Minpt = fs->make<TH1F>("h_ratio_pixTrk_trkP_Minpt", "", 40, 0, 20);
  h_ratio_pixTrk_trkP_Minpt->Sumw2();


  h_trkP_Bmass = fs->make<TH1F>("h_trkP_Bmass", "", 80, 0., 8.); 
  h_gen_Bmass = fs->make<TH1F>("h_gen_Bmass", "", 80, 0., 8.); 
  h_pixTrk_Bmass = fs->make<TH1F>("h_pixTrk_Bmass", "", 80, 0., 8.); 
  h_pixTrk_pR_Bmass = fs->make<TH1F>("h_pixTrk_pR_Bmass", "", 80, 0., 8.); 
  h_ratio_pixTrk_trkP_Bmass = fs->make<TH1F>("h_ratio_pixTrk_trkP_Bmass", "", 80, 0., 8.);

  h_trkP_Bmass->Sumw2();
  h_gen_Bmass->Sumw2();
  h_pixTrk_Bmass->Sumw2();
  h_pixTrk_pR_Bmass->Sumw2();
  h_ratio_pixTrk_trkP_Bmass->Sumw2();
  //                             
  h_trkP_Bdxy = fs->make<TH1F>("h_trkP_Bdxy", "", 50, -1., 1.);
  h_pixTrk_pR_Bdxy = fs->make<TH1F>("h_pixTrk_pR_Bdxy", "", 50, -1., 1.);
  h_ratio_pixTrk_trkP_Bdxy = fs->make<TH1F>("h_ratio_pixTrk_trkP_Bdxy", "", 50, -1., 1.);
  h_trkP_Bdxy->Sumw2();
  h_pixTrk_pR_Bdxy->Sumw2();
  h_ratio_pixTrk_trkP_Bdxy->Sumw2();

}

void BParkTriggerEfficiency::analyze(const edm::Event &event, const edm::EventSetup &c) {
  using namespace std;
  //  std::cout << " new event " << std::endl;

  ++nTotEvents; 

  edm::Handle<TrackingParticleCollection> trackingPH;
  edm::Handle<TrackingVertexCollection> trackingVH;

  event.getByToken(trackingTruthTrk, trackingPH);
  event.getByToken(trackingTruthVtx, trackingVH);


  edm::Handle<edm::View<reco::GenParticle>>  genPart;
  event.getByToken(genToken_,genPart);

  edm::Handle<std::vector<int>> barCodes;
  event.getByToken(genParticleInts_, barCodes);

  edm::Handle<reco::TrackCollection> pixelTrk;
  event.getByToken(inputPixTrkToken_, pixelTrk);

  edm::Handle<reco::VertexCollection> pixelVtx;
  event.getByToken(inputPixVtxToken_, pixelVtx);

  edm::Handle<reco::SimToRecoCollection > simToRecoCollection;
  event.getByToken(trackSimToRecoAssociationToken_, simToRecoCollection);
  reco::SimToRecoCollection simToRecoMap = *(simToRecoCollection.product());

  edm::Handle<std::vector<PileupSummaryInfo>> hPU;
  event.getByToken(puToken_, hPU);
  const auto& summary = *hPU;

  auto it = std::find_if(summary.begin(), summary.end(), [](const auto& s) { return s.getBunchCrossing() == 0; });
  if(it == summary.end()){
    std::cout << "Did not find PileupSummaryInfo in bunch crossing 0";
    return;
  }
  float trueNumInteractions = it->getTrueNumInteractions();

  const auto& genVtxPositionHandle = event.getHandle(genVtxPositionToken_);
  float PVx = genVtxPositionHandle->x();
  float PVy = genVtxPositionHandle->y();
  float PVz = genVtxPositionHandle->z();


  if(debug)  std::cout << " \n " << std::endl;


  //first get the B decay chain 
  int B_idx = -1;
  int l1_idx = -1;
  int l2_idx = -1;
  int k_idx = -1;

  TLorentzVector Bdecay[4];
  Point BdecayVtx;
  int BdecayOk = findGenBToKee(genPart, B_idx, l1_idx, l2_idx,  k_idx, Bdecay[3], Bdecay[0], Bdecay[1], Bdecay[2], BdecayVtx);
  
  if(!BdecayOk) return;
  ++nEvents_Bdecay;

  int etaMax = (std::abs(Bdecay[0].Eta()) > std::abs(Bdecay[1].Eta())) ? 0 : 1;
  etaMax = (std::abs(Bdecay[etaMax].Eta()) > std::abs(Bdecay[2].Eta())) ? etaMax : 2;
  float maxEta_Bdecay = Bdecay[etaMax].Eta();

  int ptMin = (Bdecay[0].Pt() < Bdecay[1].Pt()) ? 0 : 1;
  ptMin = (Bdecay[ptMin].Pt() < Bdecay[2].Pt()) ? ptMin : 2;
  float minPt_Bdecay = Bdecay[ptMin].Pt();


  if(debug)  
    std::cout << " >> B_idx = " << B_idx << " l1_idx = " << l1_idx 
	      << " l2_idx = " << l2_idx << " k_idx = " << k_idx << std::endl;

  const reco::GenParticlePtr BgenP = reco::GenParticlePtr(genPart, B_idx);
  const reco::Candidate* l1genP = BgenP->daughter(l1_idx);
  const reco::Candidate* l2genP = BgenP->daughter(l2_idx);  
  const reco::Candidate* kgenP = BgenP->daughter(k_idx);


  float BvtxX = BdecayVtx.x();
  float BvtxY = BdecayVtx.y();
  float BvtxZ = BdecayVtx.z();
  float Bdxy = sqrt(pow(BvtxX - PVx, 2) + pow(BvtxY - PVy, 2));


  // 0 = e1, 1 = e2, 3 = k, 4 = B
  bool matched_gen[4] = {0, 0, 0, 0};
  bool matched_pxlT[4] = {0, 0, 0, 0};
  TLorentzVector pixelTracksReco[4];
  
  if(debug)
    std::cout << " \n BgenP->pt() = " << BgenP->pt() 
	      << " l1genP->pt() = " << l1genP->pt()
	      << " l2genP->pt() = " << l2genP->pt() 
	      << " kgenP->pt() = " << kgenP->pt() << std::endl;
  
  //  std::cout << " Bgen.M = " << Bgen.M() << " L1gen.Pt() = " << L1gen.Pt() << std::endl;


  //then pick the trackingParticles associated to the B decay chain and the matched pixelTracks
  for(unsigned int i=0; i<trackingPH->size(); ++i){
    TrackingParticleRef tpr(trackingPH, i);
    //std::cout << " >>> tpr.pt() = " << tpr->pt() << " pdgId() = " << tpr->pdgId() << std::endl;
    // std::cout << " simToRecoCollection.isValid() = " << simToRecoCollection.isValid() 
    // 	      << " simToRecoMap.size() = " << simToRecoMap.size() << std::endl;

    /*
    const reco::GenParticleRefVector& genPV = tpr->genParticles();
    if(tpr->genParticles().size() > 0){
      reco::GenParticleRef genP = tpr->genParticles().at(0);
      std::cout << " >>> genP.id() = " << genP.id()  << " genP.key() = " << genP.key() 
		<< " genPart.id() = " << genPart.id() << " BgenP.id() = " << BgenP.id() << std::endl;
      std::cout << " associated genParticles = " << tpr->genParticles().size() << " genP[0]->pdgId() = " << genP->pdgId() 
		<< " genP[0]->pt() = " << genP->pt() << std::endl;
    }
    */


    int idxMatched = -1;
    
    if(std::abs(tpr->pt() - Bdecay[0].Pt()) < 1.e-4){
      idxMatched = 0;
      
      if(debug)
	std::cout << " tpr->pt() = " << tpr->pt() << " matched to l1 " 
		  << " pdgId() = " << tpr->pdgId() << std::endl;
    }
    if(std::abs(tpr->pt() - Bdecay[1].Pt()) < 1.e-4){
      idxMatched = 1;

      if(debug)
	std::cout << " tpr->pt() = " << tpr->pt() << " matched to l2 " 
		  << " pdgId() = " << tpr->pdgId() << std::endl;
    }
    if(std::abs(tpr->pt() - Bdecay[2].Pt()) < 1.e-4){
      idxMatched = 2;

      if(debug)
	std::cout << " tpr->pt() = " << tpr->pt() << " matched to k " 
		  << " pdgId() = " << tpr->pdgId() << std::endl;
    }

    if(idxMatched != -1 ) {
      matched_gen[idxMatched] = 1;
      h_trkP_eta[idxMatched]->Fill(Bdecay[idxMatched].Eta());
      h_trkP_phi[idxMatched]->Fill(Bdecay[idxMatched].Phi());
      h_trkP_pt[idxMatched]->Fill(Bdecay[idxMatched].Pt());
      h_trkP_pu[idxMatched]->Fill(trueNumInteractions);
    }
    else continue;

    //now look for trackingParticles-pixelTracks matching
    std::vector<std::pair<edm::RefToBase<reco::Track>, double>> trackV;

    if(simToRecoCollection.isValid() && simToRecoMap.find(tpr) != simToRecoMap.end()){
      trackV = simToRecoMap[tpr];
      //std::cout << " trackV.size() = " << trackV.size() << std::endl;
      if (trackV.size()!=0) {
	edm::RefToBase<reco::Track> track = trackV.begin()->first;
	//reco::TrackRef trkRef = trackV.begin().first;
	double associationQuality = trackV.begin()->second;
	// std::cout << " track->pt() = " << track->pt() << " quality = " << associationQuality 
	// 	  << " matched to " << idxMatched << std::endl;
	matched_pxlT[idxMatched] = 1;
	float massX = idxMatched == 2 ? Bdecay[2].M() : Bdecay[1].M();
	pixelTracksReco[idxMatched].SetPtEtaPhiM(track->pt(), track->eta(), track->phi(), massX);

	h_pixTrk_pu[idxMatched]->Fill(trueNumInteractions);
	h_pixTrk_eta[idxMatched]->Fill(pixelTracksReco[idxMatched].Eta());
	h_pixTrk_phi[idxMatched]->Fill(pixelTracksReco[idxMatched].Phi());
	h_pixTrk_pt[idxMatched]->Fill(pixelTracksReco[idxMatched].Pt());

	h_pixTrk_pR_pu[idxMatched]->Fill(trueNumInteractions);
	h_pixTrk_pR_eta[idxMatched]->Fill(Bdecay[idxMatched].Eta());
	h_pixTrk_pR_phi[idxMatched]->Fill(Bdecay[idxMatched].Phi());
	h_pixTrk_pR_pt[idxMatched]->Fill(Bdecay[idxMatched].Pt());
      } // trkP - pixelTrk match found
    }// trackingP-pixelTrak association map ok
  }//trackingP


  if(matched_gen[0] == matched_gen[1] && matched_gen[2] == matched_gen[1] && matched_gen[0] == 1){
    ++nEvents_trkP_matched;
    matched_gen[3] = 1;
    int idxMatched = 3;
    Bdecay[idxMatched] = TLorentzVector(Bdecay[0] + Bdecay[1] + Bdecay[2]);
    h_trkP_pu[idxMatched]->Fill(trueNumInteractions);
    h_trkP_eta[idxMatched]->Fill(Bdecay[idxMatched].Eta()); 
    h_trkP_phi[idxMatched]->Fill(Bdecay[idxMatched].Phi()); 
    h_trkP_pt[idxMatched]->Fill(Bdecay[idxMatched].Pt()); 
    h_trkP_Bmass->Fill(Bdecay[idxMatched].M()); 

    h_trkP_Bdxy->Fill(Bdxy);
    h_trkP_Maxeta->Fill(maxEta_Bdecay);
    h_trkP_Minpt->Fill(minPt_Bdecay);
    //    std::cout << " >> trackingParticles matched " << std::endl;
  }
  else{
    for(int ij=0; ij<4; ++ij){
      h_gen_pu[ij]->Fill(trueNumInteractions);
      h_gen_eta[ij]->Fill(Bdecay[ij].Eta());
      h_gen_phi[ij]->Fill(Bdecay[ij].Phi());
      h_gen_pt[ij]->Fill(Bdecay[ij].Pt());
      h_gen_Bmass->Fill(Bdecay[ij].M());
    }
  }
    
  if(matched_pxlT[0] == matched_pxlT[1] && matched_pxlT[2] == matched_pxlT[1] && matched_pxlT[0] == 1){
    ++nEvents_pixTrk_matched;
    matched_pxlT[3] = 1;
    int idxMatched = 3;
    pixelTracksReco[idxMatched] = TLorentzVector(pixelTracksReco[0] + pixelTracksReco[1] + pixelTracksReco[2]);
    h_pixTrk_pu[idxMatched]->Fill(trueNumInteractions);
    h_pixTrk_eta[idxMatched]->Fill(pixelTracksReco[idxMatched].Eta());
    h_pixTrk_phi[idxMatched]->Fill(pixelTracksReco[idxMatched].Phi());
    h_pixTrk_pt[idxMatched]->Fill(pixelTracksReco[idxMatched].Pt());
    h_pixTrk_Bmass->Fill(pixelTracksReco[idxMatched].M());

    h_pixTrk_pR_pu[idxMatched]->Fill(trueNumInteractions);
    h_pixTrk_pR_eta[idxMatched]->Fill(Bdecay[idxMatched].Eta());
    h_pixTrk_pR_phi[idxMatched]->Fill(Bdecay[idxMatched].Phi());
    h_pixTrk_pR_pt[idxMatched]->Fill(Bdecay[idxMatched].Pt());
    h_pixTrk_pR_Bmass->Fill(Bdecay[idxMatched].M());

    h_pixTrk_pR_Bdxy->Fill(Bdxy);
    h_pixTrk_pR_Maxeta->Fill(maxEta_Bdecay);
    h_pixTrk_pR_Minpt->Fill(minPt_Bdecay);
    //    std::cout << " >> pixelTracks matched " << std::endl;
  }

  int Nmatched = 0;
  int Nmatchedxy = 0;

  Point vtxPos;

  for(auto ij : *pixelVtx){
    if(!ij.isValid()) continue;
    Nmatched = 0;

    if(Nmatchedxy != 0 && 
       (std::abs(ij.position().x() - vtxPos.x()) > 1.e-4 || 
	std::abs(ij.position().y() - vtxPos.y()) > 1.e-4)) Nmatchedxy = 0;

    //    std::cout << " vtxPos = " << ij.position() << " trkSize = " << ij.tracksSize() << std::endl;
    for(reco::Vertex::trackRef_iterator kk = ij.tracks_begin(); kk != ij.tracks_end(); ++kk){
      //std::cout << " trkW = " << ij.trackWeight(*kk) << " trk_Pt = " << (*kk)->pt() << std::endl;

      if(std::abs((*kk)->pt() - pixelTracksReco[0].Pt()) < 1.e-4){
	//	std::cout << " tpr->pt() = " << (*kk)->pt() << " matched to l1 " << std::endl;
	++Nmatched;
	vtxPos = ij.position();
	++Nmatchedxy;
      }
      if(std::abs((*kk)->pt() - pixelTracksReco[1].Pt()) < 1.e-4){
	//	std::cout << " tpr->pt() = " << (*kk)->pt() << " matched to l2 " << std::endl;
	++Nmatched;
	vtxPos = ij.position();
	++Nmatchedxy;
      }
      if(std::abs((*kk)->pt() - pixelTracksReco[2].Pt()) < 1.e-4){
	//	std::cout << " tpr->pt() = " << (*kk)->pt() << " matched to k " << std::endl;
	++Nmatched;
	vtxPos = ij.position();
	++Nmatchedxy;
      }
    }
    if(Nmatched == 3) break;
  }

  if(Nmatched == 3) ++NVtxMatched; // std::cout << " vtxMatched " << std::endl;
  if(Nmatchedxy == 3) ++NVtxMatched_xy;

}



int BParkTriggerEfficiency::findGenBToKee(edm::Handle<edm::View<reco::GenParticle>> genPart, 
					  int& B_bc, int& l1_bc, int& l2_bc, int& k_bc, 
					  TLorentzVector& Bgen, TLorentzVector& L1gen,
					  TLorentzVector& L2gen, TLorentzVector& Kgen,
					  Point& BdecayVtx){

  int nGenPart = genPart->size();
  if(debug) std::cout << " \n\n nGenPart = " << nGenPart << std::endl;
  int genpart_B_index = -1;
  int genpart_lep1FromB_index = -1;
  int genpart_lep2FromB_index = -1;
  int genpart_kaonFromB_index = -1;
  //int genpart_JPsiFromB_index = -1;

  TLorentzVector ele1;
  TLorentzVector ele2;
  TLorentzVector kaon;

  int isNonResonant = -1;
  int leptonFlavour = 11;

  bool originVertexFound = false;
  Point originVertex;
  Point e1Vtx;
  Point e2Vtx;
  Point kVtx;

  for(int i_Bu=0; i_Bu<nGenPart; ++i_Bu){
    const reco::GenParticle Bpart = (*genPart)[i_Bu];
    int BpdgId = Bpart.pdgId();

    if(debug) std::cout << " particle pdg = " << BpdgId << std::endl;

    if(std::abs(BpdgId) == 521){

      if(debug) std::cout << " foundB pdg = " << BpdgId << std::endl;

      unsigned int B_nD = Bpart.numberOfDaughters();
      if(debug) std::cout << " B nDaug = " << B_nD << std::endl;

      if(B_nD > 0){
	if(genpart_lep1FromB_index == -1 ||
	   genpart_lep2FromB_index == -1 ||
	   genpart_kaonFromB_index == -1){
	  genpart_lep1FromB_index = -1;
	  genpart_lep2FromB_index = -1;
	  genpart_kaonFromB_index = -1;
	}
        for(unsigned int iD=0; iD < B_nD; ++iD){
          const reco::Candidate* B_daug = Bpart.daughter(iD);

          if(debug) std::cout << " B_daug->pdgId() = " << B_daug->pdgId() << std::endl;
	  if(B_daug->status() != 1) continue;
          if(B_daug->pdgId() == leptonFlavour){
	    if(originVertexFound && B_daug->vertex() != originVertex) continue;
	    ele1.SetPtEtaPhiM(B_daug->pt(), B_daug->eta(), B_daug->phi(), B_daug->mass());
	    genpart_lep1FromB_index = iD;
	    e1Vtx = B_daug->vertex();

	    if(e1Vtx == e2Vtx){
	      originVertex = e1Vtx;
	      originVertexFound = true;
	    }
	    else if(e1Vtx == kVtx){
	      originVertex = e1Vtx;
	      originVertexFound = true;
	    }
          }
          if(B_daug->pdgId() == -1 * leptonFlavour){
	    if(originVertexFound && B_daug->vertex() != originVertex) continue;
	    ele2.SetPtEtaPhiM(B_daug->pt(), B_daug->eta(), B_daug->phi(), B_daug->mass());
	    genpart_lep2FromB_index = iD;
	    e2Vtx = B_daug->vertex();

	    if(e2Vtx == e1Vtx){ 
	      originVertex = e2Vtx;
	      originVertexFound = true;
	    }
            else if(e2Vtx == kVtx){
	      originVertex = e2Vtx;
	      originVertexFound = true;
	    }
          }
          if(std::abs(B_daug->pdgId()) == 321){
	    if(originVertexFound && B_daug->vertex() != originVertex) continue;                                  
            kaon.SetPtEtaPhiM(B_daug->pt(), B_daug->eta(), B_daug->phi(), B_daug->mass());
            genpart_kaonFromB_index = iD;
	    kVtx = B_daug->vertex();

	    if(kVtx == e1Vtx){
              originVertex = kVtx;
              originVertexFound = true;
            }
            else if(e2Vtx == kVtx){
              originVertex = kVtx;
              originVertexFound = true;
            }
          }
          if(genpart_lep1FromB_index != -1 && originVertexFound == true &&
             genpart_lep2FromB_index != -1 && genpart_kaonFromB_index != -1){
            isNonResonant = 1;
	    if(kVtx != e1Vtx || e2Vtx != e1Vtx){
	      std::cout << " problem with Vtx " << " kVtx = " << kVtx << " e1Vtx = " << e1Vtx << " e2Vtx = " << e2Vtx << std::endl;
	      continue;
	    }
	    genpart_B_index = i_Bu;
	    BdecayVtx = e1Vtx;
	    //std::cout << " genpart_lep1FromB_index =  " << genpart_lep1FromB_index << " l1 vertex = " << e1Vtx << std::endl;
            if(debug) std::cout << " non resonant " << std::endl;
            //break;                                                             
          }
        }
      }
      /*
      else if(B_nD == 2){
	genpart_lep1FromB_index = -1;
	genpart_lep2FromB_index = -1;
	genpart_kaonFromB_index = -1;
	genpart_JPsiFromB_index = -1;
	for(unsigned int iD=0; iD < B_nD; ++iD){
	  const reco::Candidate* B_daug = Bpart.daughter(iD);
	  if(debug) std::cout << " B_daug->pdgId() = " << B_daug->pdgId() << std::endl;
	  if(B_daug->pdgId() == 443){
	    if(debug) std::cout << " JPsi found = " << std::endl;
	    genpart_JPsiFromB_index = iD;
	    unsigned int JPsi_nD = B_daug->numberOfDaughters();
	    if(JPsi_nD == 2){
	      for(unsigned int igD=0; igD < JPsi_nD; ++igD){
		const reco::Candidate* JPsi_daug = B_daug->daughter(igD);
		if(JPsi_daug->pdgId() == leptonFlavour){
		  ele1.SetPtEtaPhiM(JPsi_daug->pt(), JPsi_daug->eta(), JPsi_daug->phi(), JPsi_daug->mass());
		  genpart_lep1FromB_index = igD;
		}
		if(JPsi_daug->pdgId() == -1 * leptonFlavour){
		  ele2.SetPtEtaPhiM(JPsi_daug->pt(), JPsi_daug->eta(), JPsi_daug->phi(), JPsi_daug->mass());
		  genpart_lep2FromB_index = igD;
		}
	      }
	    }
	  }
	  if(std::abs(B_daug->pdgId()) == 321){
	    kaon.SetPtEtaPhiM(B_daug->pt(), B_daug->eta(), B_daug->phi(), B_daug->mass());
	    genpart_kaonFromB_index = iD;
	  }
	  if(genpart_lep1FromB_index != -1 && genpart_lep2FromB_index != -1 &&
	     genpart_kaonFromB_index != -1 && genpart_JPsiFromB_index != -1){
	    isNonResonant = 0;
	    genpart_B_index = i_Bu;                                                                    
	    if(debug) std::cout << " resonant " << std::endl;
	  }
	}//loop over B daughters                                                                               
      }//case JPsi                                                                                             
      */
    }//case B found                                                                                            
  }// loop over GenParticles  

  if(isNonResonant == 1){
    Bgen = TLorentzVector(ele1+ele2+kaon);
    if(debug){
      std::cout << " >>> non resonant " << std::endl;
      std::cout << " Bmass = " << Bgen.M() << std::endl;
    } 

    if(ele2.Pt() > ele1.Pt()){
      TLorentzVector dummy(ele1);
      int dummyIdx = genpart_lep1FromB_index;
      ele1 = ele2;
      ele2 = dummy;
      genpart_lep1FromB_index = genpart_lep2FromB_index;
      genpart_lep2FromB_index = dummyIdx;
    }

    L1gen = ele1;
    L2gen = ele2;
    Kgen = kaon;

    B_bc = genpart_B_index;
    l1_bc = genpart_lep1FromB_index;
    l2_bc = genpart_lep2FromB_index;
    k_bc = genpart_kaonFromB_index;
    if(debug) std::cout << " >> B_bc = " << B_bc << " l1_bc = " << l1_bc << " l2_bc = " << l2_bc << " k_bc = " << k_bc << std::endl;
    return 1;
  }
  else {
    B_bc = -1;
    l1_bc = -1;
    l2_bc = -1;
    k_bc = -1;
    return 0;
  }
} 

 
void BParkTriggerEfficiency::beginJob(){
}

void BParkTriggerEfficiency::endJob(){

  for(int ij=0; ij<4; ++ij){
    h_ratio_pixTrk_trkP_pu[ij]->Divide(h_pixTrk_pR_pu[ij], h_trkP_pu[ij]);
    h_ratio_pixTrk_trkP_eta[ij]->Divide(h_pixTrk_pR_eta[ij], h_trkP_eta[ij]);
    h_ratio_pixTrk_trkP_phi[ij]->Divide(h_pixTrk_pR_phi[ij], h_trkP_phi[ij]);
    h_ratio_pixTrk_trkP_pt[ij]->Divide(h_pixTrk_pR_pt[ij], h_trkP_pt[ij]);
  }

  h_ratio_pixTrk_trkP_Bmass->Divide(h_pixTrk_pR_Bmass, h_trkP_Bmass);

  h_ratio_pixTrk_trkP_Maxeta->Divide(h_pixTrk_pR_Maxeta, h_trkP_Maxeta);
  h_ratio_pixTrk_trkP_Minpt->Divide(h_pixTrk_pR_Minpt, h_trkP_Minpt);
  h_ratio_pixTrk_trkP_Bdxy->Divide(h_pixTrk_pR_Bdxy, h_trkP_Bdxy);

  std::cout << " \n nTotEvents = " << nTotEvents << " nEvents_Bdecay = " << nEvents_Bdecay 
	    << " nEvents_trkP_matched = " << nEvents_trkP_matched << " nEvents_pixTrk_matched = " << nEvents_pixTrk_matched << std::endl;
  std::cout << " eff_gen_toAll = " << 1.*nEvents_Bdecay/nTotEvents
	    << " eff_trkP_toGen = " << 1.*nEvents_trkP_matched/nEvents_Bdecay 
	    << " eff_pixTrk_toTrkP = " << 1.*nEvents_pixTrk_matched/nEvents_trkP_matched << std::endl;

  std::cout << " nRecoVtx with recoB = " << NVtxMatched << " nRecoVtx in xy with recoB = " << NVtxMatched_xy << std::endl;
}

//matching with cone - old
//{
  /*
  if(isNonResonant != -1){
    TLorentzVector Bcand(ele1+ele2+kaon);
    std::cout << " Bmass = " << Bcand.M() << std::endl;

    h_B_gen_eta->Fill(Bcand.Eta());
    h_B_gen_phi->Fill(Bcand.Phi());
    h_B_gen_pt->Fill(Bcand.Pt());

    int forRecoB = 0;
    TLorentzVector reco_e1;
    TLorentzVector reco_e2;
    TLorentzVector reco_k;


    for(auto ij: *pixelTrk){
      float trkEta = ij.eta();
      float trkPhi = ij.phi();
      float trkPt = ij.pt();

      float dR_e1 = reco::deltaR(trkEta, trkPhi, ele1.Eta(), ele1.Phi());
      float dR_e2 = reco::deltaR(trkEta, trkPhi, ele2.Eta(), ele2.Phi());
      float dR_k = reco::deltaR(trkEta, trkPhi, kaon.Eta(), kaon.Phi());
    
      bool is_e1Matched = (dR_e1 < dR_e2 && dR_e2 < dR_k);
      bool is_e2Matched = (dR_e2 < dR_e1 && dR_e2 < dR_k);
      bool is_kMatched = (dR_k < dR_e1 && dR_k < dR_e2);


      if(is_e1Matched){
	h_dR_e1_pixel_Tracking->Fill(dR_e1);
	if(dR_e1 < 0.05){
	  h_e1_pixTrk_eta->Fill(trkEta);
	  h_e1_pixTrk_phi->Fill(trkPhi);
	  h_e1_pixTrk_pt->Fill(trkPt);
	  reco_e1.SetPtEtaPhiM(trkPt, trkEta, trkPhi, ele1.M());
	  ++forRecoB;
	}
      }
      if(is_e2Matched){
	h_dR_e2_pixel_Tracking->Fill(dR_e2);
	if(dR_e2 < 0.05){
	  h_e2_pixTrk_eta->Fill(trkEta);
	  h_e2_pixTrk_phi->Fill(trkPhi);
	  h_e2_pixTrk_pt->Fill(trkPt);
	  reco_e2.SetPtEtaPhiM(trkPt, trkEta, trkPhi, ele2.M());
	  ++forRecoB;
	}
      }
      if(is_kMatched){
	h_dR_k_pixel_Tracking->Fill(dR_k);
	if(dR_k < 0.05){
	  h_k_pixTrk_eta->Fill(trkEta);
	  h_k_pixTrk_phi->Fill(trkPhi);
	  h_k_pixTrk_pt->Fill(trkPt);
	  reco_k.SetPtEtaPhiM(trkPt, trkEta, trkPhi, kaon.M());
	  ++forRecoB;
	}
      }
    }//loop over pixel tracks

    if(forRecoB == 3){
      TLorentzVector recoB(reco_k+reco_e1+reco_e2);

      h_B_pixTrk_eta->Fill(recoB.Eta());
      h_B_pixTrk_phi->Fill(recoB.Phi());
      h_B_pixTrk_pt->Fill(recoB.Pt());
    }
    
  }//found B gen
  */
//}


DEFINE_FWK_MODULE(BParkTriggerEfficiency);
