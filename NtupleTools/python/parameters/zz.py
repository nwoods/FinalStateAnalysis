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

for iso in ['', 'LooseIso', 'NIso']:
    for fsrVar in ['pt', 'eta', 'phi']:
        varCap = fsrVar[0].upper()+fsrVar[1:]
        setattr(zzObjVars, "objectAKFSR%s%s"%(iso, varCap), 
                cms.string(('? daughterHasUserCand({object_idx}, "akFSRCand%s") ? ' +
                            'daughterUserCand({object_idx}, "akFSRCand%s").%s() : -999.')%(iso, iso, fsrVar)))
            
        setattr(zzDiObjVars, "object1_object2_%sAKFSR%s"%(varCap, iso), 
                cms.string(('diObjectP4WithUserCands({object1_idx}, {object2_idx}, "akFSRCand%s").%s')%(iso, fsrVar))
                )

        setattr(zzEvVars, '%sAKFSR%s'%(varCap, iso),
                cms.string('p4WithUserCands("akFSRCand%s").%s'%(iso, varCap)))

    setattr(eleVars, "objectRelPFIsoRhoAKFSR%s"%iso,
            cms.string(('({object}.chargedHadronIso()' +
                        '+max(0.0,{object}.neutralHadronIso()' +
                        '+{object}.photonIso()' +
                        '-ptOfDaughterUserCand({object_idx}, "akFSR%s")' +
                        '-{object}.userFloat("rhoCSA14")*{object}.userFloat("EffectiveArea_HZZ4l2015")))' +
                        '/{object}.pt()')%(iso))
            ),
                        # '-(?daughterHasUserCand({object_idx}, "akFSRCand%s") ? ' +
                        # 'daughterUserCand({object_idx}, "akFSRCand%s").pt : ' +
                        # '0.)' +

    setattr(muVars, "objectRelPFIsoDBAKFSR%s"%iso,
            cms.string(('({object}.chargedHadronIso()' +
                        '+max({object}.photonIso()' +
                        '-ptOfDaughterUserCand({object_idx}, "akFSRCand%s")' +
                        '+{object}.neutralHadronIso()' +
                        '-0.5*{object}.puChargedHadronIso,0.0))' +
                        '/{object}.pt()')%(iso))
            )
                        # '-1*(?daughterHasUserCand({object_idx}, "akFSRCand%s") ? ' +
                        # 'daughterUserCand({object_idx}, "akFSRCand%s").pt : ' +
                        # '0.)' +


    setattr(zzDiObjVars, "object1_object2_MassAKFSR%s"%(iso), 
                cms.string(('diObjectP4WithUserCands({object1_idx}, {object2_idx}, "akFSRCand%s").M')%(iso))
                )

    setattr(zzEvVars, 'MassAKFSR%s'%(iso),
            cms.string('p4WithUserCands("akFSRCand%s").M'%(iso)))

eleVars.objectGenStatus = cms.string(('? (getDaughterGenParticle({object_idx}, 11, 0, 1).isAvailable && ' +
                                      'getDaughterGenParticle({object_idx}, 11, 0, 1).isNonnull) ? ' +
                                      'getDaughterGenParticle({object_idx}, 11, 0, 1).status : -999'))
muVars.objectGenStatus = cms.string(('? (getDaughterGenParticle({object_idx}, 13, 0, 1).isAvailable && ' +
                                     'getDaughterGenParticle({object_idx}, 13, 0, 1).isNonnull) ? ' +
                                     'getDaughterGenParticle({object_idx}, 13, 0, 1).status : -999'))

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
        HZZCategory = 'userFloat("HZZCategory")',
        D_bkg_kin = 'userFloat("p0plus_VAJHU") / (userFloat("p0plus_VAJHU") + userFloat("bkg_VAMCFM"))',
        D_bkg = 'userFloat("p0plus_VAJHU") * userFloat("p0plus_m4l") / '
            '(userFloat("p0plus_VAJHU") * userFloat("p0plus_m4l") + userFloat("bkg_VAMCFM") * userFloat("bkg_m4l"))',
        D_gg = 'userFloat("Dgg10_VAMCFM")',
        D_g4 = 'userFloat("p0plus_VAJHU") / (userFloat("p0plus_VAJHU") + userFloat("p0minus_VAJHU"))',
        Djet_VAJHU = '? evt.jets.size >= 2 ? userFloat("pvbf_VAJHU") / (userFloat("pvbf_VAJHU") + userFloat("phjj_VAJHU")) : -1',
        zzEvVars,
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

    






