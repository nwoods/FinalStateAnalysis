############################################################################
##                                                                        ##
##   miniAODElectronMVAIDs_cfi.py                                         ##
##          configuration to embed MVA IDs in miniAOD pat::Electrons      ##
##                                                                        ##
##   Author: N. Woods, U. Wisconsin                                       ##
##                                                                        ##
############################################################################

import FWCore.ParameterSet.Config as cms



trigMVAWeights = [
    'EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_25ns_EB_BDT.weights.xml',
    'EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_25ns_EE_BDT.weights.xml',
    ]
nonTrigMVAWeights = [
    'EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml',
    'EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml',
    'EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml',
    'EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml',
    'EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml',
    'EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml',
    ]

miniAODElectronMVAID = cms.EDProducer(
    "MiniAODElectronMVAIDEmbedder",
    src=cms.InputTag("fixme!"),
    trigWeights = cms.vstring(*trigMVAWeights),
    trigLabel = cms.string('BDTIDTrig'), # triggering MVA ID userfloat key
    nonTrigWeights = cms.vstring(*nonTrigMVAWeights),
    nonTrigLabel = cms.string('BDTIDNonTrig') # nontriggering MVA ID userfloat key
    )
