import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from CommonTools.ParticleFlow.ParticleSelectors.pfSelectedMuons_cfi import pfSelectedMuons 
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import * 
from CommonTools.ParticleFlow.ParticleSelectors.pfSelectedElectrons_cfi import pfSelectedElectrons 
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

switchOnTrigger(process,sequence='patDefaultSequence',hltProcess='*')

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
		filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
		src=cms.InputTag('offlinePrimaryVertices')
)

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")

process.producePatPFMETCorrections.replace(
    process.pfCandMETcorr,
    process.type0PFMEtCorrection *
    process.patPFMETtype0Corr *
    process.pfCandMETcorr
)

########### Run PF2PAT

postfix = "PFlow"
usePF2PAT(process,
	  runPF2PAT=True,
          jetAlgo='AK5', 
	  runOnMC=True, 
	  postfix=postfix,
          jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
	  typeIMetCorrections=True,
	  pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
)

getattr(process,"pfPileUp"+postfix).checkClosestZVertex = False

########### Initialize lepton and PU removal from jets

usePFnoPU     = True
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

########## to turn on MET type-0+1 corrections

getattr(process,'patMETs'+postfix).metSource = cms.InputTag('patType1p2CorrectedPFMet'+postfix)

getattr(process,'patType1CorrectedPFMet'+postfix).srcType1Corrections = cms.VInputTag(
    cms.InputTag("patPFJetMETtype1p2Corr"+postfix,"type1"),
    cms.InputTag("patPFMETtype0Corr"+postfix)
)

########## to turn on MET x/y shift corrections

process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc

process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
		    cms.InputTag('pfJetMETcorr', 'type1'),
		        cms.InputTag('pfMEtSysShiftCorr')  
)
process.pfType1p2CorrectedMet.srcType1Corrections = cms.VInputTag(
		    cms.InputTag('pfJetMETcorr', 'type1'),
		    cms.InputTag('pfMEtSysShiftCorr')       
)

########### mu Trigger Matching

pathTriggerMu = 'path("HLT_Mu17_Mu8*")'

process.selectedTriggeredPatMuons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatMuons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerMu),   # selection of matches
                                                   maxDPtRel   = cms.double(0.5), 
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
                                                   )

process.selectedPatMuonsTriggerMatch = cms.EDProducer("PATTriggerMatchMuonEmbedder",
		     src     = cms.InputTag('selectedPatMuons'+postfix),
	             matches = cms.VInputTag('selectedTriggeredPatMuons')
)

############## Making Jets

process.goodJets = selectedPatJets.clone(
		src = cms.InputTag('selectedPatJets'+postfix),
                     cut = cms.string(
			'pt > 25. & abs(eta) < 2.5 &'
			'numberOfDaughters > 1 &'
			'neutralHadronEnergyFraction < 0.99 &'
			'neutralEmEnergyFraction < 0.99 &'
			'chargedEmEnergyFraction < 0.99 &'
			'chargedHadronEnergyFraction > 0 &'
			'chargedMultiplicity > 0'
		)
)

############## Making Z to mumu

process.matchedMuons0 = cms.EDProducer("MuScleFitPATMuonCorrector",
                         src = cms.InputTag("selectedPatMuonsTriggerMatch"),
                         debug = cms.bool(False),
                         identifier = cms.string("Summer12_DR53X_smearReReco"),
                         applySmearing = cms.bool(True),
                         fakeSmearing = cms.bool(False)
)

process.matchedMuons = selectedPatMuons.clone(
		src = cms.InputTag('matchedMuons0'),
		cut = cms.string(
			'pt > 20 & abs(eta) < 2.4 &'
			'isGlobalMuon & isPFMuon &'
			'globalTrack.normalizedChi2 < 10 &'
			'track.hitPattern.trackerLayersWithMeasurement > 5 &'
			'globalTrack.hitPattern.numberOfValidMuonHits > 0 &'
			'innerTrack.hitPattern.numberOfValidPixelHits > 0 &'
			'abs(dB) < 0.2 &'
			'numberOfMatchedStations > 1 &'
			'(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt,0.0))/pt < 0.2 &'
			'triggerObjectMatches.size > 0'
		)
)

process.zmuMatchedmuMatched = cms.EDProducer('CandViewShallowCloneCombiner',
			          decay = cms.string('matchedMuons@+ matchedMuons@-'),
				  cut   = cms.string('mass > 60.0 & mass < 120.0'),
				  name  = cms.string('Zmumatchedmumatched'),
				  roles = cms.vstring('matched1', 'matched2')
				  )

############## e Trigger Matching

pathTriggerEle = 'path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*")'

process.selectedTriggeredPatElectrons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatElectrons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerEle),  # selection of matches
                                                   maxDPtRel   = cms.double(0.5), 
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
                                                   )

process.selectedPatElectronsTriggerMatch = cms.EDProducer("PATTriggerMatchElectronEmbedder",
		     src     = cms.InputTag('selectedPatElectrons'+postfix),
	             matches = cms.VInputTag('selectedTriggeredPatElectrons')
		     )

##############

switchOnTriggerMatching(process,triggerMatchers = ['selectedTriggeredPatMuons','selectedTriggeredPatElectrons'],sequence ='patDefaultSequence',hltProcess = '*')

removeCleaningFromTriggerMatching(process)

############## Making Z to ee

process.matchedElectrons = selectedPatElectrons.clone(
		     src = cms.InputTag('selectedPatElectronsTriggerMatch'),
		     cut = cms.string(
			'pt > 10 & abs(eta) < 2.4 &'
			'(('
			 'abs(superCluster.eta) < 1.442 &'
			 'abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &'
			 'abs(deltaPhiSuperClusterTrackAtVtx) < 0.06 &'
			 'sigmaIetaIeta < 0.01 &'
			 'hadronicOverEm < 0.12'
			')|('
			 'abs(superCluster.eta) > 1.566 & abs(superCluster.eta) < 2.5 &'
			 'abs(deltaEtaSuperClusterTrackAtVtx) < 0.007 &'
			 'abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 &'
			 'sigmaIetaIeta < 0.03 &'
			 'hadronicOverEm < 0.10'
			')) &'
			'abs(dB) < 0.02 &'
			'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &'
			'(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et < 0.15 &'
			'passConversionVeto &'
			'gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 &'
			'triggerObjectMatches.size > 0'
		     )
)

process.zeleMatchedeleMatched = cms.EDProducer('CandViewShallowCloneCombiner',
			            decay = cms.string('matchedElectrons@+ matchedElectrons@-'),
				    cut   = cms.string('mass > 60.0 & mass < 120.0'),
				    name  = cms.string('Zelematchedelematched'),
				    roles = cms.vstring('matched1', 'matched2')
				    )

##############

process.GlobalTag.globaltag = 'START53_V23::All'

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)

process.source = cms.Source("PoolSource",
	#fileNames = cms.untracked.vstring('/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START50_V15-v1/0000/88AD6E87-E173-E111-9996-00E081791749.root')
	fileNames = cms.untracked.vstring('file:88AD6E87-E173-E111-9996-00E081791749.root')
			   )

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options.wantSummary = True
process.maxEvents.input = 1000

getattr(process,"pfElectronsFromVertex"+postfix).d0Cut = 0.02
getattr(process,"pfElectronsFromVertex"+postfix).dzCut = 0.1
getattr(process,"pfMuonsFromVertex"+postfix).d0Cut = 0.2
getattr(process,"pfMuonsFromVertex"+postfix).dzCut = 0.5

getattr(process,"pfIsolatedElectrons"+postfix).isolationValueMapsCharged = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFId'+postfix))
getattr(process,"pfIsolatedElectrons"+postfix).deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFId'+postfix)
getattr(process,"pfIsolatedElectrons"+postfix).isolationValueMapsNeutral = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFId'+postfix), cms.InputTag('elPFIsoValueGamma03PFId'+postfix))
getattr(process,"patElectrons"+postfix).isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag('elPFIsoValueCharged03PFId'+postfix),
        pfChargedAll = cms.InputTag('elPFIsoValueChargedAll03PFId'+postfix),
        pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFId'+postfix),
        pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral03PFId'+postfix),
        pfPhotons = cms.InputTag('elPFIsoValueGamma03PFId'+postfix)
        )

getattr(process,"pfIsolatedElectrons"+postfix).doDeltaBetaCorrection = True
getattr(process,"pfIsolatedElectrons"+postfix).isolationCut = 0.15
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection     = True
getattr(process,"pfIsolatedMuons"+postfix).isolationCut = 0.2

getattr(process,"patJets"+postfix).addTagInfos = True

########### top projections in PF2PAT:

getattr(process,"pfNoPileUp"+postfix).enable = usePFnoPU
getattr(process,"pfNoMuon"+postfix).enable = useNoMuon
getattr(process,"pfNoJet"+postfix).enable = useNoJet
getattr(process,"pfNoTau"+postfix).enable = useNoTau
getattr(process,"pfNoElectron"+postfix).enable = useNoElectron

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.demo3_ee = cms.EDAnalyzer('GenbAnalyzer',
	path = cms.untracked.string("."),
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("electron"),
)
process.demo3_mm = cms.EDAnalyzer('GenbAnalyzer',
	path = cms.untracked.string("."),
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("muon"),
)
process.MyProcess = cms.EDFilter('ZbFilter')
process.demo_ee = cms.EDAnalyzer('ZbAnalyzer',
	path = cms.untracked.string("."),
        pileup  = cms.untracked.string("S10"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(0),
        usePartonFlavour = cms.untracked.bool(False)
)
process.demo_mm = cms.EDAnalyzer('ZbAnalyzer',
	path = cms.untracked.string("."),
        pileup  = cms.untracked.string("S10"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(0),
        usePartonFlavour = cms.untracked.bool(False)
)
process.demo_ee_btag = cms.EDAnalyzer('ZbAnalyzer',
	path = cms.untracked.string("."),
        pileup  = cms.untracked.string("S10"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(0),
        usePartonFlavour = cms.untracked.bool(True)
)
process.demo_mm_btag = cms.EDAnalyzer('ZbAnalyzer',
	path = cms.untracked.string("."),
        pileup  = cms.untracked.string("S10"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(0),
        usePartonFlavour = cms.untracked.bool(True)
)
process.p = cms.Path(
   process.goodOfflinePrimaryVertices *
   getattr(process,"patPF2PATSequence"+postfix) *
   process.recoTauClassicHPSSequence *
   process.goodJets *
   process.pfMEtSysShiftCorrSequence *
   process.producePFMETCorrections *
   process.patTrigger *
   process.selectedTriggeredPatMuons *
   process.selectedPatMuonsTriggerMatch *
   process.matchedMuons0 *
   process.matchedMuons *
   process.zmuMatchedmuMatched *
   process.selectedTriggeredPatElectrons *
   process.selectedPatElectronsTriggerMatch *
   process.matchedElectrons *
   process.zeleMatchedeleMatched *
   process.demo3_ee * process.demo3_mm *
   process.MyProcess *
   process.demo_ee * process.demo_mm *
   process.demo_ee_btag * process.demo_mm_btag
   #process.dump
   )

process.out.fileName = 'patTuple.root'
process.outpath = cms.EndPath()

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('rootTuple.root')
)

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.out.outputCommands = patEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += [
	'keep *',
	'keep *_addPileupInfo_*_*',
	'keep *_matchedElectrons_*_*',
	'keep *_matchedMuons_*_*',
	'keep *_goodJets_*_*',
	'keep *_zeleMatchedeleMatched_*_*',
	'keep *_zmuMatchedmuMatched_*_*',
	'keep *_goodOfflinePrimaryVertices_*_*'
]

