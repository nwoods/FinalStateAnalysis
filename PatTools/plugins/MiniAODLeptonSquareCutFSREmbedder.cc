//////////////////////////////////////////////////////////////////////////////
///                                                                        ///
///    MiniAODLeptonSquareCutFSREmbedder.cc                                ///
///                                                                        ///
///    From a collection of photons, a collection of muons, and a          ///
///        collection of electrons: make sure each photon is not in an     ///
///        electron supercluster, and pair it to its closest lepton.       ///
///        For each lepton, embed a reco::Candidate with the 4-momentum    ///
///        of all photons passing a dR and eT cut as a usercand. Cut       ///
///        strings may be supplied for all three typesof objects.          ///
///                                                                        ///
///    Author: Nate Woods, U. Wisconsin                                    ///
///                                                                        ///
//////////////////////////////////////////////////////////////////////////////


// system include files
#include <memory>
#include <iostream>
#include <math.h> // pow

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"


typedef reco::Candidate Cand;
typedef edm::Ptr<Cand> CandPtr;
typedef reco::PFCandidate PFCand;
typedef edm::Ptr<PFCand> PFCandPtr;
typedef reco::CandidateView CandView;
typedef pat::Electron Elec;
typedef edm::Ptr<pat::Electron> ElecPtr;
typedef edm::View<pat::Electron> ElecView;
typedef pat::Muon Muon;
typedef edm::Ptr<pat::Muon> MuonPtr;
typedef edm::View<pat::Muon> MuonView;

class MiniAODLeptonSquareCutFSREmbedder : public edm::EDProducer
{
public:
  explicit MiniAODLeptonSquareCutFSREmbedder(const edm::ParameterSet&);
  ~MiniAODLeptonSquareCutFSREmbedder();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  edm::EDGetTokenT<CandView> photons_;
  edm::EDGetTokenT<ElecView> electrons_;
  edm::EDGetTokenT<MuonView> muons_;

  bool isInSuperCluster(const CandPtr& cand, const std::vector<ElecPtr>& elecs) const;
  
  StringCutObjectSelector<Cand> phoSelection_;
  StringCutObjectSelector<Elec> eSelection_;
  StringCutObjectSelector<Muon> mSelection_;

  std::string fsrLabel_;

  const float scVetoDR_;
  const float scVetoDEta_;
  const float scVetoDPhi_;
  const float dRCut_;
  const float etCut_;
};


MiniAODLeptonSquareCutFSREmbedder::MiniAODLeptonSquareCutFSREmbedder(const edm::ParameterSet& iConfig):
  photons_(consumes<CandView>(iConfig.getParameter<edm::InputTag>("phoSrc"))),
  electrons_(consumes<ElecView>(iConfig.getParameter<edm::InputTag>("eSrc"))),
  muons_(consumes<MuonView>(iConfig.getParameter<edm::InputTag>("muSrc"))),
  phoSelection_(iConfig.exists("phoSelection") ? 
                iConfig.getParameter<std::string>("phoSelection") :
                ""),
  eSelection_(iConfig.exists("eSelection") ?
	      iConfig.getParameter<std::string>("eSelection") :
	      ""),
  mSelection_(iConfig.exists("muSelection") ? 
	      iConfig.getParameter<std::string>("muSelection") :
	      ""),
  fsrLabel_(iConfig.exists("fsrLabel") ?
            iConfig.getParameter<std::string>("fsrLabel") :
            "dREtFSRCand"),
  scVetoDR_(iConfig.exists("scVetoDR") ?
            float(iConfig.getParameter<double>("scVetoDR")) :
            0.15),
  scVetoDEta_(iConfig.exists("scVetoDEta") ?
              float(iConfig.getParameter<double>("scVetoDEta")) :
              0.05),
  scVetoDPhi_(iConfig.exists("scVetoDPhi") ?
              float(iConfig.getParameter<double>("scVetoDPhi")) :
              2.),
  dRCut_(iConfig.exists("dRCut") ?
         float(iConfig.getParameter<double>("dRCut")) :
         0.3),
  etCut_(iConfig.exists("etCut") ?
         float(iConfig.getParameter<double>("etCut")) :
         4.)
{
  produces<std::vector<Muon> >();
  produces<std::vector<Elec> >();
  produces<std::vector<PFCand> >();
}


MiniAODLeptonSquareCutFSREmbedder::~MiniAODLeptonSquareCutFSREmbedder()
{
}


void MiniAODLeptonSquareCutFSREmbedder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<Muon> > mOut( new std::vector<Muon> );
  std::auto_ptr<std::vector<Elec> > eOut( new std::vector<Elec> );
  edm::Handle<CandView> phos;
  iEvent.getByToken(photons_, phos);
  edm::Handle<edm::View<Elec> > elecs;
  iEvent.getByToken(electrons_, elecs);
  edm::Handle<edm::View<Muon> > mus;
  iEvent.getByToken(muons_, mus);

  // get a cleaned electron collection because we need it for the SC veto anyway
  std::vector<ElecPtr> cleanedElectrons;
  for(size_t iE = 0; iE < elecs->size(); ++iE)
    {
      ElecPtr elec = elecs->ptrAt(iE);
      if(eSelection_(*elec))
        cleanedElectrons.push_back(elec);
    }

  // associate photons to their closest leptons
  std::vector<std::vector<CandPtr> > phosByEle = std::vector<std::vector<CandPtr> >(elecs->size());
  std::vector<std::vector<CandPtr> > phosByMu = std::vector<std::vector<CandPtr> >(mus->size());

  for( size_t iPho = 0; iPho != phos->size(); ++iPho )
    {
      CandPtr pho = phos->ptrAt(iPho);
      
      // basic selection
      if (!phoSelection_(*pho))
        continue;

      // supercluster veto
      if(isInSuperCluster(pho, cleanedElectrons))
        continue;

      size_t iBestEle = 9999;
      size_t iBestMu = 9999;
      float dRBestEle = 9999.;
      float dRBestMu = 9999.;

      for(size_t iE = 0; iE < elecs->size(); ++iE)
        {
          float deltaR = reco::deltaR(pho->p4(), elecs->at(iE).p4());
          if(deltaR < dRBestEle)
            {
              iBestEle = iE;
              dRBestEle = deltaR;
            }
        }

      for(size_t iM = 0; iM < mus->size(); ++iM)
        {
          float deltaR = reco::deltaR(pho->p4(), mus->at(iM).p4());
          if(deltaR < dRBestMu)
            {
              iBestMu = iM;
              dRBestMu = deltaR;
            }
        }

      if(elecs->size() && dRBestEle < dRBestMu)
        phosByEle.at(iBestEle).push_back(pho);
      else if(mus->size())
        phosByMu.at(iBestMu).push_back(pho);
    }

  // because we need to embed a Ptr to the new photon, and that doesn't
  // exist do Event::put, we have to keep track of which new photon each
  // lepton is paired to. 9999 means no match.
  std::vector<size_t> newIndByEle = std::vector<size_t>();
  std::vector<size_t> newIndByMu = std::vector<size_t>();
  std::auto_ptr<std::vector<PFCand> > newPhos( new std::vector<PFCand> );
      
  for(size_t iE = 0; iE < elecs->size(); ++iE)
    {
      Elec e = elecs->at(iE);
      
      std::vector<CandPtr> passing = std::vector<CandPtr>();

      for(size_t iPho = 0; iPho < phosByEle[iE].size(); ++iPho)
        {
          CandPtr pho = phosByEle[iE][iPho];

          float dR = reco::deltaR(e.p4(), pho->p4());
          float et = pho->et();
          
          if(dR < dRCut_ && et > etCut_)
            passing.push_back(pho);
        }

      if(passing.size())
        {
          reco::Particle::LorentzVector p4 = passing[0]->p4();
          for(size_t i = 1; i < passing.size(); ++i)
            p4 += passing.at(i)->p4();
              
          newIndByEle.push_back(newPhos->size());
          newPhos->push_back(reco::PFCandidate(0, p4, reco::PFCandidate::gamma));
        }
      else
        newIndByEle.push_back(9999);
    }

  for(size_t iM = 0; iM < mus->size(); ++iM)
    {
      Muon m = mus->at(iM);
      
      std::vector<CandPtr> passing = std::vector<CandPtr>();
      
      for(size_t iPho = 0; iPho < phosByMu[iM].size(); ++iPho)
        {
          CandPtr pho = phosByMu[iM][iPho];

          float dR = reco::deltaR(m.p4(), pho->p4());
          float et = pho->et();

          if(dR < dRCut_ && et > etCut_)
            passing.push_back(pho);
        }

      if(passing.size())
        {
          reco::Particle::LorentzVector p4 = passing[0]->p4();
          for(size_t i = 1; i < passing.size(); ++i)
            p4 += passing.at(i)->p4();
              
          newIndByMu.push_back(newPhos->size());
          newPhos->push_back(reco::PFCandidate(0, p4, reco::PFCandidate::gamma));
        }
      else
        newIndByMu.push_back(9999);
    }

  edm::OrphanHandle<std::vector<PFCand> > newProd = iEvent.put(newPhos);

  for(size_t iE = 0; iE < elecs->size(); ++iE)
    {
      Elec e = elecs->at(iE);

      if(newIndByEle.at(iE) < newProd->size())
        {
          PFCandPtr phoPtr(newProd, newIndByEle.at(iE));

          e.addUserCand(fsrLabel_, phoPtr);
        }

      eOut->push_back(e);
    }

  for(size_t iM = 0; iM < mus->size(); ++iM)
    {
      Muon m = mus->at(iM);

      if(newIndByMu.at(iM) < newProd->size())
        {
          PFCandPtr phoPtr(newProd, newIndByMu.at(iM));
          m.addUserCand(fsrLabel_, phoPtr);
        }

      mOut->push_back(m);
    }

  iEvent.put( eOut );
  iEvent.put( mOut );

}


bool MiniAODLeptonSquareCutFSREmbedder::isInSuperCluster(const CandPtr& cand, 
                                                         const std::vector<ElecPtr>& elecs) const
{
  for(auto elec = elecs.begin(); elec != elecs.end(); ++elec)
    {
      float dR = reco::deltaR(cand->eta(), cand->phi(), (*elec)->superCluster()->eta(), (*elec)->superCluster()->phi());
      if(dR < scVetoDR_)
        return true;

      float dEta = fabs((*elec)->superCluster()->eta() - cand->eta());
      float dPhi = fabs(reco::deltaPhi((*elec)->superCluster()->phi(), cand->phi()));
      if(dEta < scVetoDEta_ && dPhi < scVetoDPhi_)
        return true;
    }

  return false;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODLeptonSquareCutFSREmbedder);

