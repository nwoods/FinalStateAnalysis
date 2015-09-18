/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// MiniAODFSRGenMatchProducer                                              //
//                                                                         //
// Tries to match FSR candidates attached to electrons and muons to        //
// packed gen photons. Embeds whether or not a match was found into the    //
// lepton as a userFloat (0 or 1).                                         //
// Also embeds userFloat "hasGenFSR" in the lepton to indicate (0 or 1)    //
// whether or not any gen FSR was found. Only gen FSR photons that have    //
// a matching reco photon are considered, to factor reco efficiency out    //
// of FSR efficiency studies.                                              //
//                                                                         //
// Nate Woods, U. Wisconsin                                                //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


// system include files
#include <memory>
#include <iostream>

// CMS include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef reco::Candidate Cand;
typedef edm::Ptr<Cand> CandPtr;
typedef reco::CandidateView CandView;
typedef pat::Electron Elec;
typedef edm::Ptr<pat::Electron> ElecPtr;
typedef edm::View<pat::Electron> ElecView;
typedef pat::Muon Muon;
typedef edm::Ptr<pat::Muon> MuonPtr;
typedef edm::View<pat::Muon> MuonView;
typedef pat::PackedGenParticle PGen;
typedef edm::Ptr<pat::PackedGenParticle> PGenPtr;
typedef edm::View<pat::PackedGenParticle> PGenView;


class MiniAODFSRGenMatchProducer : public edm::stream::EDProducer<>
{
public:
  explicit MiniAODFSRGenMatchProducer(const edm::ParameterSet& pset);
  ~MiniAODFSRGenMatchProducer() { ; }

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // Take a reco lepton, match it to a gen lepton, return gen FSR from that
  // gen lepton if there is any
  template<class P>
  std::pair<PGenPtr, std::vector<PGenPtr> > 
  getGenFSR(const P& lepton, const std::vector<PGenPtr>& possibleLeptons,
            const std::vector<PGenPtr>& pgens,
            const StringCutObjectSelector<PGen>& selection) const;
  // to be specialized for electron and muon to call the method above with 
  // correct selection
  template<class P>
  std::pair<PGenPtr, std::vector<PGenPtr> > 
  getGenFSR(const P& lepton, 
            const std::vector<PGenPtr>& possibleLeptons,
            const std::vector<PGenPtr>& pgens) const;

  // Embed gen match decisions in lepton
  template<class P>
  P& embedGenFSRMatch(P& lep, 
                      const std::vector<PGenPtr>& possibleLeptons,
                      const StringCutObjectSelector<P>& selection,
                      const std::vector<PGenPtr>& pgens) const;
  // to be specialized for electron and muon to call the method above with 
  // correct selection
  template<class P>
  P& embedGenFSRMatch(P& lep, 
                      const std::vector<PGenPtr>& possibleLeptons,
                      const std::vector<PGenPtr>& pgens) const;

  bool isAncestorOf(const PGenPtr& p, const reco::GenParticleRef& a) const;

  reco::GenParticleRef lastInterestingAncestor(const PGenPtr& p, 
                                               const unsigned int interestingId) const;
  bool hasInterestingAncestor(const PGenPtr& p, 
                              const unsigned int interestingId) const;
  const float maxDR_;
  const unsigned int ancestorId_;
  const std::vector<std::string> fsrLabels_;

  const edm::EDGetTokenT<ElecView> eSrc_;
  const edm::EDGetTokenT<MuonView> mSrc_;
  const edm::EDGetTokenT<CandView> gSrc_;
  const edm::EDGetTokenT<PGenView> gGenSrc_;

  const StringCutObjectSelector<Elec> eSelection_;
  const StringCutObjectSelector<Muon> mSelection_;
  const StringCutObjectSelector<Cand> gSelection_;
  const StringCutObjectSelector<PGen> eGenSelection_;
  const StringCutObjectSelector<PGen> mGenSelection_;
  const StringCutObjectSelector<PGen> gGenSelection_;
};


MiniAODFSRGenMatchProducer::MiniAODFSRGenMatchProducer(const edm::ParameterSet& pset) :
  maxDR_(pset.getParameter<double>("maxDR")),
  ancestorId_(pset.getParameter<unsigned int>("ancestorId")),
  fsrLabels_(pset.getParameter<std::vector<std::string> >("fsrLabels")),
  eSrc_(consumes<ElecView>(pset.getParameter<edm::InputTag>("eSrc"))),
  mSrc_(consumes<MuonView>(pset.getParameter<edm::InputTag>("mSrc"))),
  gSrc_(consumes<CandView>(pset.getParameter<edm::InputTag>("phoSrc"))),
  gGenSrc_(consumes<PGenView>(pset.getParameter<edm::InputTag>("genSrc"))),
  eSelection_(pset.getParameter<std::string>("eSelection")),
  mSelection_(pset.getParameter<std::string>("mSelection")),
  gSelection_(pset.getParameter<std::string>("phoSelection")),
  eGenSelection_(pset.getParameter<std::string>("eGenSelection")),
  mGenSelection_(pset.getParameter<std::string>("mGenSelection")),
  gGenSelection_(pset.getParameter<std::string>("phoGenSelection"))
{
  produces<std::vector<Elec> >();
  produces<std::vector<Muon> >();
}


template<>
Muon& MiniAODFSRGenMatchProducer::embedGenFSRMatch(Muon& lep,
                                                   const std::vector<PGenPtr>& possibleMatches,
                                                   const std::vector<PGenPtr>& pgens) const
{
  return embedGenFSRMatch(lep, possibleMatches, mSelection_, pgens);
}


template<>
Elec& MiniAODFSRGenMatchProducer::embedGenFSRMatch(Elec& lep,
                                                   const std::vector<PGenPtr>& possibleMatches,
                                                   const std::vector<PGenPtr>& pgens) const
{
  return embedGenFSRMatch(lep, possibleMatches, eSelection_, pgens);
}


template<class P>
P& MiniAODFSRGenMatchProducer::embedGenFSRMatch(P& lep, 
                                                const std::vector<PGenPtr>& possibleMatches,
                                                const StringCutObjectSelector<P>& selection,
                                                const std::vector<PGenPtr>& pgens) const
{
  std::pair<PGenPtr, std::vector<PGenPtr> > genFSR; // matching lepton and vector of gen FSR
  if(selection(lep))
    genFSR = getGenFSR(lep, possibleMatches, pgens);

  if(genFSR.second.size())
    {
      lep.addUserFloat("hasGenFSR", 1.);
      // Store some info about the highest pt gen FSR photon
      size_t indOfBest = 0;
      for(size_t i = 0; i < genFSR.second.size(); ++i)
        {
          if(genFSR.second.at(i)->pt() > genFSR.second.at(indOfBest)->pt())
            indOfBest = i;
        }
      lep.addUserCand("genFSR", genFSR.second.at(indOfBest));
      lep.addUserFloat("genFSRDR", reco::deltaR(genFSR.first->p4(),
                                                genFSR.second.at(indOfBest)->p4()));
    }
  else
    {
      lep.addUserFloat("hasGenFSR", 0.);
      return lep;
    }

  for(auto iLabel = fsrLabels_.begin(); 
      iLabel != fsrLabels_.end(); 
      ++iLabel)
    {
      if(!lep.hasUserCand(*iLabel))
        {
          lep.addUserFloat((*iLabel)+"GenMatch", 0.);
          continue;
        }

      const CandPtr fsr = lep.userCand(*iLabel);
      if((!fsr) || (!gSelection_(*fsr)))
        {
          lep.addUserFloat((*iLabel)+"GenMatch", 0.);
          continue;
        }          

      for(auto iGenPho = genFSR.second.begin(); iGenPho != genFSR.second.end(); ++iGenPho)
        {
          if(reco::deltaR(fsr->p4(), (*iGenPho)->p4()) < maxDR_)
            lep.addUserFloat((*iLabel)+"GenMatch", 1.);
          else
            lep.addUserFloat((*iLabel)+"GenMatch", 0.);
        }
    }

  return lep;
}


template<class P>
std::pair<PGenPtr, std::vector<PGenPtr> > 
MiniAODFSRGenMatchProducer::getGenFSR(const P& lep, 
                                      const std::vector<PGenPtr>& possibleMatches,
                                      const std::vector<PGenPtr>& pgens,
                                      const StringCutObjectSelector<PGen>& selection) const
{
  // find gen match for lepton
  PGenPtr match;
  for(size_t iGen = 0; iGen < possibleMatches.size(); ++iGen)
    {
      PGenPtr pg = possibleMatches.at(iGen);
      if(!pg)
        continue;

      if(reco::deltaR(pg->p4(), lep.p4()) >= maxDR_)
        continue;

      if(!hasInterestingAncestor(pg, ancestorId_))
        continue;

      if(!match || reco::deltaR(match->p4(), lep.p4()) >= reco::deltaR(pg->p4(), lep.p4()))
        match = pg;
    }

  std::vector<PGenPtr> genFSR;
  if(!match)
    return std::make_pair(match, genFSR);

  for(size_t iPho = 0; iPho < pgens.size(); ++iPho)
    {
      reco::GenParticleRef radiator = lastInterestingAncestor(pgens.at(iPho), 
                                                              abs(lep.pdgId()));
      if(!radiator)
        continue;

      if(isAncestorOf(match, radiator))
        genFSR.push_back(pgens.at(iPho));
    }

  return std::make_pair(match, genFSR);
}


template<>
std::pair<PGenPtr, std::vector<PGenPtr> > 
MiniAODFSRGenMatchProducer::getGenFSR(const Muon& lepton,
                                      const std::vector<PGenPtr>& possibleMatches,
                                      const std::vector<PGenPtr>& pgens) const
{
  return getGenFSR(lepton, possibleMatches, pgens, mGenSelection_);
}


template<>
std::pair<PGenPtr, std::vector<PGenPtr> > 
MiniAODFSRGenMatchProducer::getGenFSR(const Elec& lepton,
                                      const std::vector<PGenPtr>& possibleMatches,
                                      const std::vector<PGenPtr>& pgens) const
{
  return getGenFSR(lepton, possibleMatches, pgens, eGenSelection_);
} 


void MiniAODFSRGenMatchProducer::produce(edm::Event& ev, const edm::EventSetup& setup)
{
  edm::Handle<ElecView> elecs;
  ev.getByToken(eSrc_, elecs);
  edm::Handle<MuonView> muons;
  ev.getByToken(mSrc_, muons);
  edm::Handle<PGenView> pgens;
  ev.getByToken(gGenSrc_, pgens);
  edm::Handle<CandView> phocands;
  ev.getByToken(gSrc_, phocands);

  // light skims of gen and reco photons for a little speed boost
  std::vector<CandPtr> recPhos;
  for(size_t i = 0; i < phocands->size(); ++i)
    {
      CandPtr rg = phocands->ptrAt(i);
      if(abs(rg->pdgId()) == 22 && gSelection_(*rg))
        recPhos.push_back(rg);
    }

  //std::cout << recPhos.size() << " reco photons" << std::endl;

  std::vector<PGenPtr> genPhos;
  for(size_t i = 0; i < pgens->size(); ++i)
    {
      PGenPtr pg = pgens->ptrAt(i);
      if(abs(pg->pdgId()) != 22 || pg->status() != 1 || (!gGenSelection_(*pg)))
        continue;
      for(size_t j = 0; j < recPhos.size(); ++j)
        {
          if(reco::deltaR(pg->p4(), recPhos.at(j)->p4()) < maxDR_)
            {
              genPhos.push_back(pg);
              break;
            }
        }
    }

  //std::cout << genPhos.size() << " gen photons" << std::endl;

  // and some for leptons too
  std::vector<PGenPtr> genElecs;
  for(size_t i = 0; i < pgens->size(); ++i)
    {
      PGenPtr pg = pgens->ptrAt(i);
      if(abs(pg->pdgId()) != 11 || pg->status() != 1 || (!eGenSelection_(*pg)))
        continue;
      genElecs.push_back(pg);
    }

  //std::cout << genElecs.size() << " gen electrons" << std::endl;

  std::vector<PGenPtr> genMuons;
  for(size_t i = 0; i < pgens->size(); ++i)
    {
      PGenPtr pg = pgens->ptrAt(i);
      if(abs(pg->pdgId()) != 13 || pg->status() != 1 || (!mGenSelection_(*pg)))
        continue;
      genMuons.push_back(pg);
    }

  //std::cout << genMuons.size() << " gen muons" << std::endl;

  std::auto_ptr<std::vector<Elec> > eOut(new std::vector<Elec>);
  for(size_t i = 0; i < elecs->size(); ++i)
    {
      Elec lep = elecs->at(i);
      lep = embedGenFSRMatch(lep, genElecs, genPhos);
      eOut->push_back(lep);
    }

  std::auto_ptr<std::vector<Muon> > mOut(new std::vector<Muon>);
  for(size_t i = 0; i < muons->size(); ++i)
    {
      Muon lep = muons->at(i);
      lep = embedGenFSRMatch(lep, genMuons, genPhos);
      mOut->push_back(lep);
    }

  ev.put(eOut);
  ev.put(mOut);
}


bool
MiniAODFSRGenMatchProducer::isAncestorOf(const PGenPtr& p,
                                         const reco::GenParticleRef& a) const
{
  reco::GenParticleRef m = p->motherRef();

  while(m->pdgId() == p->pdgId())
    {
      if(m == a)
        return true;

      m = m->motherRef();
    }

  return false;
}


reco::GenParticleRef 
MiniAODFSRGenMatchProducer::lastInterestingAncestor(const PGenPtr& p, 
                                                    const unsigned int interestingId) const
{
  reco::GenParticleRef m = p->motherRef();

  while(m.isNonnull() && m->numberOfMothers() != 0)
    {
      if(abs(m->pdgId()) == interestingId)
        break;

      // nothing I care about should ever have more than 1 mother
      m = m->motherRef();
    }
  
  return m;
}


bool
MiniAODFSRGenMatchProducer::hasInterestingAncestor(const PGenPtr& p, 
                                                   const unsigned int interestingId) const
{
  return lastInterestingAncestor(p, interestingId).isNonnull();
}



//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODFSRGenMatchProducer);
