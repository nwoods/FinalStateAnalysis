# Parameters to be used in production of ZZ ntuples
# Only parameters seen here are used. make_ntuples_cfg.py loads these first
# and then loads any modifications to these parameters from a custom param file
# passed via paramFile=/path/to/param/file.py

import FWCore.ParameterSet.Config as cms
from FinalStateAnalysis.Utilities.cfgtools import PSet
from collections import OrderedDict

zzEvVars = PSet()
zzObjVars = PSet()
zzDiObjVars = PSet()
eleVars = PSet()
muVars = PSet()

for fsr in ['dretFSR', 'dret15FSR', 'et4DR03FSR', 'et4DR01FSR']:
    brSuffix = fsr.replace("dret", "DREt").replace("et4DR", "Et4DR")
    for maybeDM in ['', "DM"]:
        for fsrVar in ['pt', 'eta', 'phi']:
            varCap = fsrVar[0].upper()+fsrVar[1:]
            setattr(zzObjVars, "object%s%s"%(brSuffix, varCap), 
                    cms.string(('? daughterHasUserCand({object_idx}, "%sCand") ? ' +
                                'daughterUserCand({object_idx}, "%sCand").%s() : -999.')%(fsr, fsr, fsrVar)))
                
            setattr(zzDiObjVars, "object1_object2_%s%s%s"%(varCap, brSuffix, maybeDM), 
                    cms.string(('diObjectP4WithUserCands%s({object1_idx}, {object2_idx}, "%sCand").%s')%(maybeDM, fsr, fsrVar))
                    )
        
            setattr(zzEvVars, '%s%s%s'%(varCap, brSuffix, maybeDM),
                    cms.string('p4WithUserCands%s("%sCand").%s'%(maybeDM, fsr, varCap)))
        
        setattr(eleVars, "objectRelPFIsoRho%s%s"%(brSuffix, maybeDM),
                cms.string(('({object}.chargedHadronIso()' +
                            '+max(0.0,{object}.neutralHadronIso()' +
                            '+{object}.photonIso()' +
                            '-daughterUserCandIsoContribution%s({object_idx}, "%sCand")' +
                            '-{object}.userFloat("rhoCSA14")*{object}.userFloat("EffectiveArea_HZZ4l2015")))' +
                            '/{object}.pt()')%(maybeDM, fsr))
                ),
        
        setattr(muVars, "objectRelPFIsoDB%s%s"%(brSuffix, maybeDM),
                cms.string(('({object}.chargedHadronIso()' +
                            '+max({object}.photonIso()' +
                            '-daughterUserCandIsoContribution%s({object_idx}, "%sCand")' +
                            '+{object}.neutralHadronIso()' +
                            '-0.5*{object}.puChargedHadronIso,0.0))' +
                            '/{object}.pt()')%(maybeDM, fsr))
                )
        
        setattr(zzDiObjVars, "object1_object2_Mass%s%s"%(brSuffix, maybeDM), 
                cms.string(('diObjectP4WithUserCands%s({object1_idx}, {object2_idx}, "%sCand").M')%(maybeDM, fsr))
                )
        
        setattr(zzEvVars, 'Mass%s%s'%(brSuffix, maybeDM),
                cms.string('p4WithUserCands%s("%sCand").M'%(maybeDM, fsr)))
        
        eleVars.objectDREt = cms.string(('? daughterHasUserCand({object_idx}, "%sCand") ? ' +
                                         'daughterAsElectron({object_idx}).userFloat("%sDREt") : ' +
                                         '-999.')%(fsr, fsr))
        
        muVars.objectDREt = cms.string(('? daughterHasUserCand({object_idx}, "%sCand") ? ' +
                                  'daughterAsMuon({object_idx}).userFloat("%sCandDREt") : ' +
                                  '-999.')%(fsr, fsr))

setattr(eleVars, "objectRelPFIsoRhoFSR",
        cms.string(('({object}.chargedHadronIso()' +
                    '+max(0.0,{object}.neutralHadronIso()' +
                    '+{object}.photonIso()' +
                    '-allFSRIsoContribution({object_idx}, "FSRCand")' +
                    '-{object}.userFloat("rhoCSA14")*{object}.userFloat("EffectiveArea_HZZ4l2015")))' +
                    '/{object}.pt()'))
        ),

setattr(muVars, "objectRelPFIsoDBFSR",
        cms.string(('({object}.chargedHadronIso()' +
                    '+max({object}.photonIso()' +
                    '-allFSRIsoContribution({object_idx}, "FSRCand")' +
                    '+{object}.neutralHadronIso()' +
                    '-0.5*{object}.puChargedHadronIso,0.0))' +
                    '/{object}.pt()'))
        )

zzEvVars.MassFSR = cms.string('p4fsr("FSRCand").M')
zzEvVars.PtFSR   = cms.string('p4fsr("FSRCand").pt')
zzEvVars.EtaFSR  = cms.string('p4fsr("FSRCand").eta')
zzEvVars.PhiFSR  = cms.string('p4fsr("FSRCand").phi')
zzEvVars.MtFSR   = cms.string('p4fsr("FSRCand").Mt')
zzEvVars.nJets   = cms.string('evt.jets.size')

zzDiObjVars.object1_object2_MassFSR  = cms.string('subcandPrimaryFSR({object1_idx}, {object2_idx}, "FSRCand").get.mass')
zzDiObjVars.object1_object2_PtFSR    = cms.string('subcandPrimaryFSR({object1_idx}, {object2_idx}, "FSRCand").get.pt')
zzDiObjVars.object1_object2_EtaFSR   = cms.string('subcandPrimaryFSR({object1_idx}, {object2_idx}, "FSRCand").get.eta')
zzDiObjVars.object1_object2_PhiFSR   = cms.string('subcandPrimaryFSR({object1_idx}, {object2_idx}, "FSRCand").get.phi')
zzDiObjVars.object1_object2_MtFSR    = cms.string('subcandPrimaryFSR({object1_idx}, {object2_idx}, "FSRCand").get.mt')
zzDiObjVars.object1_object2_FSRPt    = cms.string('? bestFSROfZ({object1_idx}, {object2_idx}, "FSRCand").isNonnull() ? '
                                                  'bestFSROfZ({object1_idx}, {object2_idx}, "FSRCand").pt() : -999.')
zzDiObjVars.object1_object2_FSREta   = cms.string('? bestFSROfZ({object1_idx}, {object2_idx}, "FSRCand").isNonnull() ? ' 
                                                  'bestFSROfZ({object1_idx}, {object2_idx}, "FSRCand").eta() : -999.')
zzDiObjVars.object1_object2_FSRPhi   = cms.string('? bestFSROfZ({object1_idx}, {object2_idx}, "FSRCand").isNonnull() ? '
                                                  'bestFSROfZ({object1_idx}, {object2_idx}, "FSRCand").phi() : -999.')


parameters = {
    # selections on all objects whether they're included in final states or not, done immediately after necessary variables are embedded
    'preselection' : OrderedDict(
        [
            # Veto electrons that are very close to muons
            ('e', {
                    'm' : {
                        'deltaR' : 0.05,
                        'selection' : 'userFloat("HZZ4lIDPassTight") > 0.5',
                        },
                    },
             ),
            # Remove jets that are near tight ID'd, no-FSR isolated leptons
            ('j', {
                    'selection' : 'pt > 30 && eta < 4.7 && eta > -4.7 && userFloat("puID") > 0.5 && userFloat("idLoose") > 0.5',
                    'e' : {
                        'deltaR' : 0.4,
                        'selection' : 'userFloat("HZZ4lIDPassTight") > 0.5 && userFloat("HZZ4lIsoPass") > 0.5',
                        },
                    'm' : {
                        'deltaR' : 0.4,
                        'selection' : 'userFloat("HZZ4lIDPassTight") > 0.5 && userFloat("HZZ4lIsoPass") > 0.5',
                        },
                    }
             )
            ]),
            
    # selections to include object in final state (should be looser than analysis selections)
    'finalSelection' : OrderedDict(
        [
            ('e', 'abs(superCluster().eta) < 3.0 && max(pt, userFloat("maxCorPt")) > 7'),
            ('m', 'max(pt, userFloat("maxCorPt")) > 4 && (isGlobalMuon || isTrackerMuon)'),
            ]
        ),
    
    # Don't automaticaly cross clean among FS objects
    'crossCleaning' : '',

    # additional variables for ntuple
    'eventVariables' : PSet(
        zzEvVars,
        HZZCategory = 'userFloat("HZZCategory")',
        D_bkg_kin = 'userFloat("p0plus_VAJHU") / (userFloat("p0plus_VAJHU") + userFloat("bkg_VAMCFM"))',
        D_bkg = 'userFloat("p0plus_VAJHU") * userFloat("p0plus_m4l") / '
            '(userFloat("p0plus_VAJHU") * userFloat("p0plus_m4l") + userFloat("bkg_VAMCFM") * userFloat("bkg_m4l"))',
        D_gg = 'userFloat("Dgg10_VAMCFM")',
        D_g4 = 'userFloat("p0plus_VAJHU") / (userFloat("p0plus_VAJHU") + userFloat("p0minus_VAJHU"))',
        Djet_VAJHU = '? evt.jets.size >= 2 ? userFloat("pvbf_VAJHU") / (userFloat("pvbf_VAJHU") + userFloat("phjj_VAJHU")) : -1',
        muVeto = 'vetoMuons(0.4, "isLooseMuon & pt > 10 & abs(eta) < 2.4").size()',
        muVetoIso = 'vetoMuons(0.4, "isLooseMuon & pt > 10 & abs(eta) < 2.4 & (chargedHadronIso()+max(0.0,neutralHadronIso()+photonIso()-userFloat(\'rhoCSA14\')*userFloat(\'EffectiveArea_HZZ4l2015\')))/pt()<0.2").size()',
        eVeto = 'vetoElectrons(0.4, "userFloat(\'CBIDLoose\')>0.5 & pt > 10 & abs(eta) < 2.5").size()',
        eVetoIso = 'vetoElectrons(0.4, "userFloat(\'CBIDLoose\')>0.5 & pt > 10 & abs(eta) < 2.5 & (chargedHadronIso()+max(0.0,neutralHadronIso()+photonIso()-userFloat(\'rhoCSA14\')*userFloat(\'EffectiveArea_HZZ4l2015\')))/pt() < 0.2").size()',
        
    ),
    # candidates of form: objectVarName = 'string expression for selection'
    'candidateVariables' : zzObjVars,
    'electronVariables' : eleVars,
    'muonVariables' : muVars,
    'tauVariables' : PSet(),
    'photonVariables' : PSet(),
    'jetVariables' : PSet(),
    # dicandidates of form: object1_object2_VarName = 'string expression for candidate'
    'dicandidateVariables' : zzDiObjVars,
}
