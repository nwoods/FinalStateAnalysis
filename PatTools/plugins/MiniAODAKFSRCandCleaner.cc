//////////////////////////////////////////////////////////////////////////////
///                                                                        ///
///    MiniAODAKFSRCandCleaner.cc                                          ///
///                                                                        ///
///    From a collection of lepton and photon reco Candidates (or          ///
///        similar), remove photons that aren't isolated or don't pass     ///
///        a string cut, and leptons that don't match PAT leptons passing  ///
///        some string cuts. Outputs as a collection of shallow clone Ptrs ///
///        to preserve provenance.                                         ///
///                                                                        ///
///    Author: Nate Woods, U. Wisconsin                                    ///
///                                                                        ///
//////////////////////////////////////////////////////////////////////////////


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
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
typedef reco::CandidateView CandView;
typedef pat::Electron Elec;
typedef edm::View<pat::Electron> ElecView;
typedef pat::Muon Muon;
typedef edm::View<pat::Muon> MuonView;
typedef reco::ShallowClonePtrCandidate CandClone;

class MiniAODAKFSRCandCleaner : public edm::EDProducer
{
public:
  explicit MiniAODAKFSRCandCleaner(const edm::ParameterSet&);
  ~MiniAODAKFSRCandCleaner();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  edm::EDGetTokenT<CandView> srcCands_;
  edm::EDGetTokenT<ElecView> electrons_;
  edm::EDGetTokenT<MuonView> muons_;
  std::vector<edm::EDGetTokenT<edm::ValueMap<float> > > isoMaps_;
  
  StringCutObjectSelector<Cand> selection_;
  StringCutObjectSelector<Elec> eSelection_;
  StringCutObjectSelector<Muon> mSelection_;

  const float isoCut_;
};


MiniAODAKFSRCandCleaner::MiniAODAKFSRCandCleaner(const edm::ParameterSet& iConfig):
  srcCands_(consumes<CandView>(iConfig.getParameter<edm::InputTag>("src"))),
  electrons_(consumes<ElecView>(iConfig.getParameter<edm::InputTag>("eSrc"))),
  muons_(consumes<MuonView>(iConfig.getParameter<edm::InputTag>("muSrc"))),
  selection_(iConfig.exists("phoSelection") ? 
	     iConfig.getParameter<std::string>("phoSelection") :
	     ""),
  eSelection_(iConfig.exists("eSelection") ?
	      iConfig.getParameter<std::string>("eSelection") :
	      ""),
  mSelection_(iConfig.exists("muSelection") ? 
	      iConfig.getParameter<std::string>("muSelection") :
	      ""),
  isoCut_(iConfig.exists("relIsoCut") ?
          float(iConfig.getParameter<double>("relIsoCut")) :
          1.)
{
  std::vector<edm::InputTag> isoTags = (iConfig.exists("isolations") ?
                                        iConfig.getParameter<std::vector<edm::InputTag> >("isoSrc") :
                                        std::vector<edm::InputTag>());
  isoMaps_ = std::vector<edm::EDGetTokenT<edm::ValueMap<float> > >();
  for(auto i = isoTags.begin(); i != isoTags.end(); ++i)
    {
      edm::EDGetTokenT<edm::ValueMap<float> > tok = consumes<edm::ValueMap<float> >(*i);
      isoMaps_.push_back(tok);
    }

  produces<std::vector<CandClone> >();
}


MiniAODAKFSRCandCleaner::~MiniAODAKFSRCandCleaner()
{
}


void MiniAODAKFSRCandCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<CandClone> > out( new std::vector<CandClone> );
  edm::Handle<CandView> cands;
  iEvent.getByToken(srcCands_, cands);
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electrons_, electrons);
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  std::vector<edm::Handle<edm::ValueMap<float> > > isoMaps = 
    std::vector<edm::Handle<edm::ValueMap<float> > >(isoMaps_.size());
  for(size_t i = 0; i < isoMaps.size(); ++i)
    iEvent.getByToken(isoMaps_.at(i), isoMaps.at(i));

  for( size_t iCand = 0; iCand != cands->size(); ++iCand )
    {
      CandPtr cand = cands->ptrAt(iCand);
      
      unsigned int id = abs(cand->pdgId());

      // if it's a photon, it must pass selection and be isolated
      if(id == 22)
        {
          if (!selection_(*cand))
            continue;

          // get total absolute isolation from valuemaps
          float iso = 0;
          for(auto iIso = isoMaps.begin(); iIso != isoMaps.end(); ++iIso)
            iso += (**iIso)[cand];

          // cut on relative isolation
          if((iso / cand->pt()) > isoCut_)
            continue;
        }
      else
        {
          // remove lepton candidates that don't correspond to good PAT leptons
          bool keep = false;
          if(id == 11)
            {
              for(auto iLep = electrons->begin(); iLep != electrons->end(); ++iLep)
                {
                  if(iLep->sourceCandidatePtr(0) == cand)
                    {
                      keep = eSelection_(*iLep);
                      break;
                    }
                }
            }
          else if(id == 13)
            {
              for(auto iLep = muons->begin(); iLep != muons->end(); ++iLep)
                {
                  if(iLep->sourceCandidatePtr(0) == cand)
                    {
                      keep = mSelection_(*iLep);
                      break;
                    }
                }
            }

          if(!keep)
            continue;
        }

      // store a reco::ShallowClonePtrCandidate to preserve provenance
      out->push_back(CandClone(cand));
    }

  iEvent.put( out );
}


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODAKFSRCandCleaner);
