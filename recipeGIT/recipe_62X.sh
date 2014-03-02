#!/bin/bash
set -o errexit
set -o nounset

# Tags for 62X

pushd $CMSSW_BASE/src

echo "Checking out PAT dataformats"

# these 3 could be updated easily - well documented 
git cms-addpkg DataFormats/PatCandidates
git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg PhysicsTools/PatUtils

git cms-addpkg DataFormats/StdDictionaries
git cms-addpkg CommonTools/ParticleFlow

# if [ "$LIMITS" = "1" ]
# then
#    echo ""     
#    echo "======================================="
#    echo "You shouldnt run limits/fits from here." 
#    echo "LIMITS should be run from 6XX, and do not require the rest of the machinery. Please, change area."
#    echo "Check https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit for updates"
#    echo "======================================="
#    echo ""     
# fi


if [ "$PATPROD" = "1" ]
then
  echo "Checking out tuple production tags"
  git cms-addpkg DataFormats/CaloRecHit
  git cms-addpkg FWCore/GuiBrowsers
  #24/10/2012 LAG -- PF Isolation for Photons
  git cms-addpkg RecoParticleFlow/PFProducer
  #git cms-cvs-history import V00-00-12  CommonTools/RecoUtils
  #    V00-00-12, cvs up -r 1.4 CommonTools/RecoUtils/BuildFile.xml
  git cms-addpkg DataFormats/HLTReco
  git cms-addpkg JetMETCorrections/Type1MET
  git cms-addpkg RecoBTag/SecondaryVertex
  git cms-addpkg RecoVertex/AdaptiveVertexFinder

#  echo "Need to update recipe for Quark Gluon Jet ID - which is the correct tag?"
  #echo "Downloading Quark Gluon Jet ID"
  #cvs co -r v1-2-3 -d QuarkGluonTagger/EightTeV UserCode/tomc/QuarkGluonTagger/EightTeV
  # Quark-gluon tagging
#  git clone https://github.com/amarini/QuarkGluonTagger.git
#  pushd $CMSSW_BASE/src/QuarkGluonTagger
#  git checkout v1-2-6
#  pushd $CMSSW_BASE/src

  echo "Checking out Tau POG recipe"

  git cms-merge-topic -u cms-tau-pog:CMSSW_6_2_X_HighPt
#  git cms-addpkg RecoTauTag/RecoTau
#  git cms-addpkg RecoTauTag/Configuration
#  git cms-addpkg CondFormats/EgammaObjects

  #git cms-addpkg RecoTauTag/RecoTau  # recipe from christian, the merge topic complained in 539, it will probably work in 5314
  # to be checked
  #git cms-merge-topic -u cms-tau-pog:CMSSW_5_3_X


  echo "Checking out EGamma POG recipe for electron corrections"
  #cvs co -r V09-00-01 RecoEgamma/EgammaTools
  #cvs co -r FB_4Jun2013 EgammaAnalysis/ElectronTools
  git cms-addpkg RecoEgamma/EgammaTools
  git clone https://github.com/cms-analysis/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools
  pushd $CMSSW_BASE/src/EgammaAnalysis/ElectronTools
  git checkout EgammaAnalysis-ElectronTools-FB_4Jun2013
  pushd $CMSSW_BASE/src

  set +o errexit
  patch -N -p0 < FinalStateAnalysis/recipe/patches/Egamma_PassAll.patch
  set -o errexit

  # Add and patch to way speed up trigger matching
  # Don't crash if patch already applied.
  set +o errexit
  echo "Applying pat trigger matching speedup"
  patch -N -p0 < FinalStateAnalysis/recipe/patches/PassStrByRef_62X_test.patch
  set -o errexit

  #Get weight files
  pushd $CMSSW_BASE/src/EgammaAnalysis/ElectronTools/data
  cat download.url | xargs wget
  popd
  #apply some patches to make everything work
#  set +o errexit
#  patch -N -p0 < FinalStateAnalysis/recipe/patches/PATObject.h.patch
#  set -o errexit

#  echo "Applying Marias b-tag patch"   
#  #doubtful that we need it now... but just in case...
#  set +o errexit
#  patch -N -p0 < FinalStateAnalysis/recipe/patches/PhysicsToolsPatAlgos_fix_btags_52X.patch
#  set -o errexit
fi

popd

