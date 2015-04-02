#!/bin/bash
set -o errexit
set -o nounset

pushd $CMSSW_BASE/src

# Check CMSSW version
MAJOR_VERSION=`echo $CMSSW_VERSION | sed "s|CMSSW_\([0-9]\)_.*|\1|"`
MINOR_VERSION=`echo $CMSSW_VERSION | sed "s|CMSSW_\([0-9]\)_\([0-9]\)_.*|\2|"`

#for standalone version of svfit
# cvs co -r V00-01-04s TauAnalysis/CandidateTools
git clone https://github.com/cms-analysis/TauAnalysis-CandidateTools.git TauAnalysis/CandidateTools
pushd $CMSSW_BASE/src/TauAnalysis/CandidateTools
git checkout TauAnalysis-CandidateTools-V00-01-04s
popd

# Add and patch to way speed up trigger matching
# Don't crash if patch already applied.
set +o errexit
echo "Applying pat trigger matching speedup"
git cms-addpkg DataFormats/PatCandidates
git apply FinalStateAnalysis/recipe/patches/DataFormats_PatCandidates_TriggerEvent.cc.patch
set -o errexit

popd
