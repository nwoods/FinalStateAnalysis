
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///   MiniAODPhotonPFIsoCalculator.cc                                                                                   ///
///   Defines plugin PhotonPFIsoCalculator                                                                              ///
///   Copied from https://github.com/VBF-HZZ/UFHZZAnalysisRun2/blob/csa14/FSRPhotons/plugins/PhotonPFIsoCalculator.cc   ///
///       I take no credit                                                                                              ///
///           ... except for the variable caching which gives a giant speedup. I take credit for that.                  ///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "TMath.h"

class PhotonPFIsoCalculator : public edm::EDProducer {
public:
  explicit PhotonPFIsoCalculator(const edm::ParameterSet&);
  ~PhotonPFIsoCalculator();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::EDGetTokenT<reco::CandidateView> photonLabel_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfLabel_;
  StringCutObjectSelector<reco::Candidate> endcapDefinition_;
  StringCutObjectSelector<pat::PackedCandidate> pfSelection_;
  double deltaR_, deltaRself_, deltaZ_;
  double directional_;
  double vetoConeEndcaps_;
  bool   debug_;

};



PhotonPFIsoCalculator::PhotonPFIsoCalculator(const edm::ParameterSet& iConfig) :
  photonLabel_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("photonLabel"))),
  pfLabel_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  endcapDefinition_(iConfig.existsAs<std::string>("endcapDefinition") ? iConfig.getParameter<std::string>("endcapDefinition") : "abs(eta) > 1.479", true),
  pfSelection_(iConfig.getParameter<std::string>("pfSelection"), true),
  deltaR_(iConfig.getParameter<double>("deltaR")),
  deltaRself_(iConfig.getParameter<double>("deltaRself")),
  deltaZ_(iConfig.existsAs<double>("deltaZ") ? iConfig.getParameter<double>("deltaZ") : 0),
  directional_(iConfig.getParameter<bool>("directional")),
  vetoConeEndcaps_(iConfig.getParameter<double>("vetoConeEndcaps")),
  debug_(iConfig.getUntrackedParameter<bool>("debug",false))
{
  produces<edm::ValueMap<float> >().setBranchAlias("pfMuIso");
}

void PhotonPFIsoCalculator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  edm::Handle<reco::CandidateView> photons;
  iEvent.getByToken(photonLabel_,photons);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfLabel_,pfcands);

  std::vector<float> isoV;
  std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler isoF(*isoM);

  if (debug_) std::cout << "Run " << iEvent.id().run() << ", Event " << iEvent.id().event() << std::endl;

  // caching for hilarious speedup
  std::vector<double> pfPt = std::vector<double>(pfcands->size());
  std::vector<double> pfEta = std::vector<double>(pfcands->size());
  std::vector<double> pfPhi = std::vector<double>(pfcands->size());
  for(size_t iPF = 0; iPF < pfcands->size(); ++iPF)
    {
      pfPt.at(iPF) = pfcands->at(iPF).pt();
      pfEta.at(iPF) = pfcands->at(iPF).eta();
      pfPhi.at(iPF) = pfcands->at(iPF).phi();
    }

  for( reco::CandidateView::const_iterator ph = photons->begin(); ph != photons->end(); ++ph )
  {
    // more caching
    double phEta = ph->eta();
    double phPhi = ph->phi();
    bool phInEndcap = endcapDefinition_(*ph);

    if (debug_) std::cout << " Photon with pt = " << ph->pt() << ", eta = " << phEta << ", phi = " << phPhi << std::endl;

    Double_t ptSum =0.;  
    math::XYZVector isoAngleSum;
    std::vector<math::XYZVector> coneParticles;

    size_t iPF = 0;
    for( pat::PackedCandidateCollection::const_iterator pf = pfcands->begin(); pf!=pfcands->end(); ++pf )
    {
      double dr = deltaR(phEta, phPhi, pfEta[iPF], pfPhi[iPF]) ;
      if (dr>=deltaR_ || dr<deltaRself_ ) continue;
      if (deltaZ_>0 && fabs(pf->vz() - ph->vz()) > deltaZ_) continue;
      // dR Veto for Gamma: no-one in EB, dR > 0.08 in EE
      if (vetoConeEndcaps_ > 0 && phInEndcap && dr < vetoConeEndcaps_) continue;

      if (!pfSelection_(*pf)) continue; 
      if (debug_) std::cout << "   pfCandidate of pdgId " << pf->pdgId() 
          << ", pt = " << pfPt[iPF] << ", dr = " << dr << ", dz = " << (pf->vz()-ph->vz()) 
          << " is in cone... " << " from PV? " << pf->fromPV() << std::endl;
      if (debug_) std::cout << "          ...and passes all vetos, so it's added to the sum." << std::endl;

      // scalar sum
      ptSum += pfPt[iPF];

      if (directional_)
      {
        // directional sum
        math::XYZVector transverse(pfEta[iPF]-phEta, reco::deltaPhi(pfPhi[iPF],phPhi), 0);
        transverse *= pfPt[iPF] / transverse.rho();
        if (transverse.rho() > 0) 
        {
          isoAngleSum += transverse;
          coneParticles.push_back(transverse);
        }
      }
    }

    if (directional_) 
    {
      double directionalPT = 0;
      for (unsigned int iPtcl = 0; iPtcl < coneParticles.size(); ++iPtcl)
      {
        directionalPT += pow(TMath::ACos( coneParticles[iPtcl].Dot(isoAngleSum) / coneParticles[iPtcl].rho() / isoAngleSum.rho() ),2) * coneParticles[iPtcl].rho();
      }
      isoV.push_back(directionalPT);
    } 
    else 
    {
      isoV.push_back(ptSum);
    }

    iPF++;
  }

  isoF.insert(photons,isoV.begin(),isoV.end());
  isoF.fill();
  iEvent.put(isoM);
}

PhotonPFIsoCalculator::~PhotonPFIsoCalculator() { }
void PhotonPFIsoCalculator::beginJob() { }
void PhotonPFIsoCalculator::endJob() { }
DEFINE_FWK_MODULE(PhotonPFIsoCalculator);

