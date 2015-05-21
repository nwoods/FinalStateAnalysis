//////////////////////////////////////////////////////////////////////////////
///                                                                        ///
///                                                                        ///
///    MiniAODObjectEmbedAKFSR.cc                                          ///
///                                                                        ///
///    Select FSR photons from lepton/photon pseudojets and embed them in  ///
///    the corresponding lepton.                                           ///
///                                                                        ///
///    Author: N. Woods, U. Wisconsin                                      ///
///                                                                        ///
///                                                                        ///
//////////////////////////////////////////////////////////////////////////////



#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/Math/interface/deltaR.h"

template<class T>
class MiniAODObjectEmbedAKFSR : public edm::EDProducer 
{
public:
  typedef std::vector<T> TCollection;
  typedef edm::Ptr<T> TPtr;
  typedef edm::View<T> TView;
  typedef edm::Ptr<reco::Jet> JetPtr;
  typedef edm::View<reco::Jet> JetView;
  typedef edm::Ptr<reco::Candidate> CandPtr;
  MiniAODObjectEmbedAKFSR(const edm::ParameterSet& pset);
  virtual ~MiniAODObjectEmbedAKFSR(){}
  void produce(edm::Event& evt, const edm::EventSetup& es);

private:
  edm::EDGetTokenT<edm::View<T> > src_;
  edm::EDGetTokenT<edm::View<reco::Jet> > jetSrc_;
  double maxDeltaR_;
  const std::string label_;
};

template<class T>
MiniAODObjectEmbedAKFSR<T>::MiniAODObjectEmbedAKFSR(const edm::ParameterSet& pset) :
  src_(consumes<TView>(pset.getParameter<edm::InputTag>("src"))),
  jetSrc_(consumes<JetView>(pset.getParameter<edm::InputTag>("jetSrc"))),
  maxDeltaR_(pset.exists("maxDeltaR") ?
             pset.getParameter<double>("maxDeltaR") :
             0.6),
  label_(pset.getParameter<std::string>("fsrLabel"))
{
  produces<TCollection>();
}

template<class T>
void MiniAODObjectEmbedAKFSR<T>::produce(
    edm::Event& evt, const edm::EventSetup& es) 
{
  std::auto_ptr<TCollection> output(new TCollection);

  edm::Handle<edm::View<T> > leptons;
  evt.getByToken(src_, leptons);
  output->reserve(leptons->size());

  unsigned int id = 0;
  if(leptons->size())
    id = abs(leptons->at(0).pdgId());

  edm::Handle<JetView> jets;
  evt.getByToken(jetSrc_, jets);

  unsigned int nKept = 0;

  for (size_t i = 0; i < leptons->size(); ++i) 
    {
      // Make a copy that we own
      T lep = leptons->at(i);
      reco::Candidate::LorentzVector lepP4 = lep.p4();
      auto lepSrcPtr = lep.sourceCandidatePtr(0);

      // Find the jet containing this lepton
      JetPtr thisJet = JetPtr();
      for (size_t j = 0; j < jets->size(); ++j) 
        {
          JetPtr jet = jets->ptrAt(j);
          
          reco::Candidate::LorentzVector jetP4 = jet->p4();
          double deltaR = reco::deltaR(lepP4, jetP4);
          if (deltaR > maxDeltaR_)
            continue;

          for(size_t iDau = 0; iDau < jet->numberOfDaughters(); ++iDau)
            {
              auto dau = jet->daughterPtr(iDau);

              if(abs(dau->pdgId()) == id)
                {
                  if(dau->masterClonePtr() == lepSrcPtr)
                      thisJet = jet;
                }
            }

          if(thisJet.isAvailable() && thisJet.isNonnull())
            break;
        }

      // If the lepton isn't in a jet, there's nothing more to do
      if(!(thisJet.isAvailable() && thisJet.isNonnull()))
        {
          output->push_back(lep);
          continue;
        }

      // get all photons, and all other leptons (which might be closer to a photon) from the jet
      std::vector<CandPtr> photons = std::vector<CandPtr>();
      std::vector<CandPtr> otherLeps = std::vector<CandPtr>();
      
      for(size_t iDau = 0; iDau < thisJet->numberOfDaughters(); ++iDau)
        {
          CandPtr dau = thisJet->daughterPtr(iDau);
          if(!(dau.isAvailable() && dau.isNonnull()))
            continue;

          if(abs(dau->pdgId()) == 11 || abs(dau->pdgId()) == 13)
            {
              if(dau->masterClonePtr() != lepSrcPtr)
                otherLeps.push_back(dau);
              continue;
            }
              
          photons.push_back(dau->masterClonePtr());
        }

      // pick the highest pt photon, ignoring photons with a closer lepton
      float bestPt = 0.;
      CandPtr bestPho = CandPtr();
      for(size_t iPho = 0; iPho < photons.size(); ++iPho)
        {
          float dR = reco::deltaR(photons[iPho]->p4(), lep.p4());

          bool foundCloserLep = false;
          for(size_t iLep = 0; iLep < otherLeps.size(); ++iLep)
            {
              if(reco::deltaR(photons[iPho]->p4(), otherLeps[iLep]->p4()) < dR)
                {
                  foundCloserLep = true;
                  break;
                }
            }
          if(foundCloserLep)
            break;

          if(photons[iPho]->pt() < bestPt)
            continue;

          bestPt = photons[iPho]->pt();
          bestPho = photons[iPho];
        }

      if(bestPt > 0.01 && bestPho.isAvailable() && bestPho.isNonnull())
        {
          nKept++;
          lep.addUserCand(label_, bestPho);
        }

      output->push_back(lep);
    }

  evt.put(output);
}


#include "FWCore/Framework/interface/MakerMacros.h"
typedef MiniAODObjectEmbedAKFSR<pat::Muon> MiniAODMuonEmbedAKFSR;
typedef MiniAODObjectEmbedAKFSR<pat::Electron> MiniAODElectronEmbedAKFSR;
DEFINE_FWK_MODULE(MiniAODMuonEmbedAKFSR);
DEFINE_FWK_MODULE(MiniAODElectronEmbedAKFSR);
