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

switchOnTrigger(process,sequence='patDefaultSequence',hltProcess='*')

process.load("PhysicsTools.PatAlgos.patSequences_cff")

########### Run PF2PAT

postfix = "PFlow"
usePF2PAT(process,
	  runPF2PAT=True,
          jetAlgo='AK5', 
	  runOnMC=False, 
	  postfix=postfix,
          jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
	  )

########### Initialize lepton and PU removal from jets

usePFnoPU     = True
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = False # before tau top projection
useNoTau      = False # before MET top projection

########### mu Trigger Matching


pathTriggerMu = 'path("HLT_Mu17_Mu8*")'
#pathTriggerMu = 'path("HLT_Mu17_Mu9999999*")'

process.selectedTriggeredPatMuons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag( 'selectedPatMuonsPFlow') ,
                                                   matched = cms.InputTag( 'patTrigger' ),    # selections of trigger objects ,
                                                   matchedCuts = cms.string(pathTriggerMu),    # selection of matches ,
                                                   maxDPtRel   = cms.double( 0.5 ), 
                                                   maxDeltaR   = cms.double( 0.3 ) ,
                                                   resolveAmbiguities    = cms.bool( True ) ,
                                                   resolveByMatchQuality = cms.bool( True )
                                                   )

#process.patTriggerEvent.patTriggerMatches = [ "selectedTriggeredPatMuons" ]

process.selectedPatMuonsTriggerMatch = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
		     src     = cms.InputTag( "selectedPatMuonsPFlow" ),
	             matches = cms.VInputTag( "selectedTriggeredPatMuons" )
		     )

############## Making Jets

process.goodJets = selectedPatJets.clone(
                     src = cms.InputTag('selectedPatJetsPFlow'),
                     cut = cms.string('pt > 30. &'
                               'abs(eta) < 2.5'
                               )
                     )

############## Making Z to mumu

process.matchedMuons = selectedPatMuons.clone(
		src = cms.InputTag('selectedPatMuonsTriggerMatch'),
		cut = cms.string('isGlobalMuon & isTrackerMuon &'
			'innerTrack.hitPattern.trackerLayersWithMeasurement>8 &'  ## new requirement in 44X due to changes in tracking
			'userFloat("RelativePFIsolationDBetaCorr") < 0.2 &' # PF isolation   
			'abs(dB) < 0.02 &' 
			'normChi2 < 10 &'
			'innerTrack.hitPattern.numberOfValidPixelHits > 0 &'
			'numberOfMatchedStations>1 &'                                   # segments matched in at least two muon stations 
			'globalTrack.hitPattern.numberOfValidMuonHits > 0 &'    # one muon hit matched to the global fit
			'pt>20 &'
			'abs(eta) < 2.4 &'
		        ' triggerObjectMatches.size > 0'
			)
		)


process.zmuMatchedmuMatched = cms.EDProducer('CandViewShallowCloneCombiner',
			          decay = cms.string('matchedMuons@+ matchedMuons@-'),
				  cut   = cms.string('mass > 71.0 & mass < 111.0'),
				  name  = cms.string('Zmumatchedmumatched'),
				  roles = cms.vstring('matched1', 'matched2')
				  )

############## e Trigger Matching

pathTriggerEle = 'path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*")'

process.selectedTriggeredPatElectrons = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag( 'selectedPatElectronsPFlow') ,
                                                   matched = cms.InputTag( 'patTrigger' ),    # selections of trigger objects ,
                                                   matchedCuts = cms.string(pathTriggerEle),    # selection of matches ,
                                                   maxDPtRel   = cms.double( 0.5 ), 
                                                   maxDeltaR   = cms.double( 0.3 ) ,
                                                   resolveAmbiguities    = cms.bool( True ) ,
                                                   resolveByMatchQuality = cms.bool( True )
                                                   )
process.selectedPatElectronsTriggerMatch = cms.EDProducer( "PATTriggerMatchElectronEmbedder",
		     src     = cms.InputTag( "selectedPatElectronsPFlow" ),
	             matches = cms.VInputTag( "selectedTriggeredPatElectrons" )
		     )

##############

switchOnTriggerMatching(process,triggerMatchers = ['selectedTriggeredPatMuons','selectedTriggeredPatElectrons'],sequence ='patDefaultSequence',hltProcess = '*')

removeCleaningFromTriggerMatching(process)

############## Making Z to ee

process.matchedElectrons = selectedPatElectrons.clone(
		     src = cms.InputTag('selectedPatElectronsTriggerMatch'),
	             cut = cms.string(#'electronID("simpleEleId85relIso") == 5 &'
			              '((abs(superCluster.eta)< 1.442)||((1.566<(abs(superCluster.eta)))&&((abs(superCluster.eta))<2.50))) &'
				      'abs(dB) < 0.02 & '
				      'triggerObjectMatches.size > 0'                           
		     )		
)

process.zeleMatchedeleMatched = cms.EDProducer('CandViewShallowCloneCombiner',
			            decay = cms.string('matchedElectrons@+ matchedElectrons@-'),
				    cut   = cms.string('mass > 71.0 & mass < 111.0'),
				    name  = cms.string('Zelematchedelematched'),
				    roles = cms.vstring('matched1', 'matched2')
				    )




process.GlobalTag.globaltag = 'FT_R_53_V6::All'
process.source = cms.Source("PoolSource",
	#fileNames = cms.untracked.vstring('/store/data/Run2012A/DoubleElectron/AOD/13Jul2012-v1/00000/FA1B4710-F3D9-E111-858E-0024E876636C.root')
        fileNames = cms.untracked.vstring('file:FA1B4710-F3D9-E111-858E-0024E876636C.root')
			   )

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options.wantSummary = True
process.maxEvents.input = 100

getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection     = True    
getattr(process,"pfIsolatedElectrons"+postfix).doDeltaBetaCorrection = True


########### top projections in PF2PAT:

getattr(process,"pfNoPileUp"+postfix).enable = usePFnoPU
getattr(process,"pfNoMuon"+postfix).enable = useNoMuon
getattr(process,"pfNoJet"+postfix).enable = useNoJet
getattr(process,"pfNoTau"+postfix).enable = useNoTau
getattr(process,"pfNoElectron"+postfix).enable = useNoElectron


process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.MyProcess = cms.EDFilter('ZbFilter')

process.p = cms.Path(
   getattr(process,"patPF2PATSequence"+postfix) *
   process.goodJets *
   process.patTrigger *
   process.selectedTriggeredPatMuons *
   process.selectedPatMuonsTriggerMatch *
   process.matchedMuons *
   process.zmuMatchedmuMatched *
   process.selectedTriggeredPatElectrons *
   process.selectedPatElectronsTriggerMatch *
   process.matchedElectrons *
   process.zeleMatchedeleMatched *
   process.MyProcess 
   #process.dump
   )

process.out.fileName = 'patTuple.root'

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.out.outputCommands = patEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += [
        'keep *_addPileupInfo_*_*',
        'keep *_matchedElectrons_*_*',
        'keep *_matchedMuons_*_*',
        'keep *_goodJets_*_*',
        'keep *_zeleMatchedeleMatched_*_*',
        'keep *_zmuMatchedmuMatched_*_*'
        ]

