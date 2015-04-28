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

doneYet = False
if not doneYet:
    for dR in [1, 3, 5]:
        for fsrVar in ['pt', 'eta', 'phi']:
            varCap = fsrVar[0].upper()+fsrVar[1:]
            setattr(zzObjVars, "objectAK%dFSR%s"%(dR, varCap), 
                    cms.string('({object}.hasUserCand("akFSRCand0p%d") ? {object}.userCand("akFSRCand0p%d").%s : -999.)'%(dR, dR, fsrVar)))
                
            setattr(zzDiObjVars, "object1_object2_%sAK%dFSR"%(varCap, dR), 
                    cms.string(('(daughterP4WithUserCand({object1_idx}, "akFSRCand0p%d") + ' +
                                'daughterP4WithUserCand({object2_idx}, "akFSRCand0p%d")).%s')%(dR, dR, varCap))
                    )
    
            setattr(zzEvVars, '%sAK%dFSR'%(varCap, dR),
                    cms.string('p4WithUserCands("akFSRCand0p%d").%s'%(dR, varCap)))
    
        setattr(zzDiObjVars, "object1_object2_MassAK%dFSR"%(dR), 
                cms.string(('(daughterP4WithUserCand({object1_idx}, "akFSRCand0p%d") + ' +
                            'daughterP4WithUserCand({object2_idx}, "akFSRCand0p%d")).M')%(dR, dR))
                )
    
        setattr(zzEvVars, 'MassAK%dFSR'%(dR),
                cms.string('p4WithUserCands("akFSRCand0p%d").M'%(dR)))

    doneYet = True


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
    'electronVariables' : PSet(),
    'muonVariables' : PSet(),
    'tauVariables' : PSet(),
    'photonVariables' : PSet(),
    'jetVariables' : PSet(),
    # dicandidates of form: object1_object2_VarName = 'string expression for candidate'
    'dicandidateVariables' : zzDiObjVars,
}

    






