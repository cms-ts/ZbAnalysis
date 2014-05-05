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
pathTriggerEleMu = 'path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL*") || path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")'

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

########### mu Trigger Matching - em version

process.selectedTriggeredPatMuonsEM = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatMuons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerEleMu),   # selection of matches
                                                   maxDPtRel   = cms.double(0.5),
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
)

process.selectedPatMuonsTriggerMatchEM = cms.EDProducer("PATTriggerMatchMuonEmbedder",
                   src     = cms.InputTag('selectedPatMuons'+postfix),
                   matches = cms.VInputTag('selectedTriggeredPatMuonsEM')
)

############## Making Jets

process.goodJets = selectedPatJets.clone(
		src = cms.InputTag('selectedPatJets'+postfix),
                     cut = cms.string(
			'pt > 15. & abs(eta) < 2.55 &'
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
                               cut   = cms.string('mass > 50.0 & mass < 200.0'),
                               name  = cms.string('Zmumatchedmumatched'),
                               roles = cms.vstring('matched1', 'matched2')
)

############## Making Z to mumu - em version

process.matchedMuons0EM = cms.EDProducer("MuScleFitPATMuonCorrector",
                         src = cms.InputTag("selectedPatMuonsTriggerMatchEM"),
                         debug = cms.bool(False),
                         identifier = cms.string("Data2012_53X_ReReco"),
                         applySmearing = cms.bool(False),
                         fakeSmearing = cms.bool(False)
)

process.matchedMuonsEM = selectedPatMuons.clone(
             src = cms.InputTag('matchedMuons0EM'),
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

############## e Trigger Matching - em version

process.selectedTriggeredPatElectronsEM = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag('selectedPatElectrons'+postfix),
                                                   matched = cms.InputTag('patTrigger'),    # selections of trigger objects
                                                   matchedCuts = cms.string(pathTriggerEleMu),  # selection of matches
                                                   maxDPtRel   = cms.double(0.5),
                                                   maxDeltaR   = cms.double(0.3),
                                                   resolveAmbiguities    = cms.bool(True),
                                                   resolveByMatchQuality = cms.bool(True)
)

process.selectedPatElectronsTriggerMatchEM = cms.EDProducer("PATTriggerMatchElectronEmbedder",
                   src     = cms.InputTag('selectedPatElectrons'+postfix),
                   matches = cms.VInputTag('selectedTriggeredPatElectronsEM')
)

##############

switchOnTriggerMatching(process,triggerMatchers = ['selectedTriggeredPatMuons','selectedTriggeredPatElectrons','selectedTriggeredPatMuonsEM','selectedTriggeredPatElectronsEM'],sequence ='patDefaultSequence',hltProcess = '*')

removeCleaningFromTriggerMatching(process)

############## Making Z to ee

process.regressedElectrons = cms.EDProducer("RegressionEnergyPatElectronProducer",
           debug                = cms.untracked.bool(False),
           inputElectronsTag    = cms.InputTag('selectedPatElectronsTriggerMatch'),
           inputCollectionType  = cms.uint32(1), # 1:PATElectron
           useRecHitCollections = cms.bool(True), # True: RecHits have not been embedded into the PATElectron
           produceValueMaps     = cms.bool(False), # False for PAT
           regressionInputFile  = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root"),
           energyRegressionType = cms.uint32(2), # Regression type - 1: ECAL regression w/o subclusters 2: ECAL regression w/ subclusters
           rhoCollection        = cms.InputTag('kt6PFJets:rho:RECO'),
           vertexCollection     = cms.InputTag('offlinePrimaryVertices'),
           # Not used if inputCollectionType is set to 1                                         
           nameEnergyReg      = cms.string("eneRegForGsfEle"),
           nameEnergyErrorReg = cms.string("eneErrorRegForGsfEle"),
           # Used only if useRecHitCollections is set to true 
           recHitCollectionEB = cms.InputTag('reducedEcalRecHitsEB'),
           recHitCollectionEE = cms.InputTag('reducedEcalRecHitsEE')
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    ),
)

process.calibratedElectrons = cms.EDProducer("CalibratedPatElectronProducer",
                inputPatElectronsTag  = cms.InputTag('regressedElectrons'),
                #inputElectronsTag  = cms.InputTag('selectedParElectronsTriggerMatch'),
                # name of the ValueMaps containing the regression outputs                               
                #nameEnergyReg      = cms.InputTag('eleRegressionEnergy:eneRegForGsfEle'),        
                #nameEnergyErrorReg = cms.InputTag('eleRegressionEnergy:eneErrorRegForGsfEle'),   
                # The rechits are needed to compute r9                                     
                #recHitCollectionEB = cms.InputTag('reducedEcalRecHitsEB'),
                #recHitCollectionEE = cms.InputTag('reducedEcalRecHitsEE'),
                #outputGsfElectronCollectionLabel = cms.string('calibratedGsfElectrons'),         
                #nameNewEnergyReg       = cms.string('eneRegForGsfEle'), # the ValueMaps are re-created with the new collection as key.
                #nameNewEnergyErrorReg  = cms.string('eneErrorRegForGsfEle'),  
                isMC              = cms.bool(True),
                verbose           = cms.bool(False),
                synchronization   = cms.bool(False),
                updateEnergyError = cms.bool(True),
                correctionsType          = cms.int32(2),
                applyLinearityCorrection = cms.bool(True),
                combinationType          = cms.int32(3),
                lumiRatio                = cms.double(0.0),
                inputDataset                   = cms.string("Summer12_LegacyPaper"),
                combinationRegressionInputPath = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root"),
                scaleCorrectionsInputPath      = cms.string("EgammaAnalysis/ElectronTools/data/scalesNewReg-May2013.csv"),
                linearityCorrectionsInputPath  = cms.string("EgammaAnalysis/ElectronTools/data/linearityNewReg-May2013.csv")
)
 
process.matchedElectrons = selectedPatElectrons.clone(
		     src = cms.InputTag('calibratedElectrons'),
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
				    cut   = cms.string('mass > 50.0 & mass < 200.0'),
				    name  = cms.string('Zelematchedelematched'),
				    roles = cms.vstring('matched1', 'matched2')
)

############## Making Z to ee - em version

process.matchedElectronsEM = selectedPatElectrons.clone(
                   src = cms.InputTag('selectedPatElectronsTriggerMatchEM'),
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

############## Making Z to em

process.zeleMatchedmuMatched = cms.EDProducer('CandViewShallowCloneCombiner',
                                  decay = cms.string('matchedElectronsEM@+ matchedMuonsEM@-'),
                                  cut   = cms.string('mass > 50.0 & mass < 200.0'),
                                  name  = cms.string('ZelematchedmumatchedEM'),
                                  roles = cms.vstring('matched1', 'matched2')
)

##############

process.GlobalTag.globaltag = 'START53_V27::All'

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
	#fileNames = cms.untracked.vstring('store/mc/Summer12_DR53X/DYJets_0p0_1p2_2p10_3p15_4p15_CT10_8TeV-sherpa/AODSIM/PU_S10_START53_V7C-v2/20001/7ADC40F9-2C8B-E211-B250-0026189438FA.root')
	#fileNames = cms.untracked.vstring('file:7ADC40F9-2C8B-E211-B250-0026189438FA.root')
	#fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/1A92B137-32F0-E111-9426-E0CB4E553693.root')
	#fileNames = cms.untracked.vstring('file:1A92B137-32F0-E111-9426-E0CB4E553693.root')
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options.wantSummary = True
process.maxEvents.input = 2000

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
process.demoEleGen = cms.EDProducer('GenbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
)
process.demoMuoGen = cms.EDProducer('GenbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
)
process.MyProcess = cms.EDFilter('ZbFilter')
process.demoEle = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron")
)
process.demoElePum = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pum"),
        lepton  = cms.untracked.string("electron")
)
process.demoElePup = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pup"),
        lepton  = cms.untracked.string("electron")
)
process.demoEleUp = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC     = cms.untracked.string("S10"),
        pileupDT     = cms.untracked.string("ee"),
        lepton     = cms.untracked.string("electron"),
        JEC        = cms.untracked.double(1)
)
process.demoEleDown = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC       = cms.untracked.string("S10"),
        pileupDT       = cms.untracked.string("ee"),
        lepton       = cms.untracked.string("electron"),
        JEC          = cms.untracked.double(-1)
)
process.demoMuo = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon")
)
process.demoMuoPum = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pum"),
        lepton  = cms.untracked.string("muon")
)
process.demoMuoPup = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pup"),
        lepton  = cms.untracked.string("muon")
)
process.demoMuoUp = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC           = cms.untracked.string("S10"),
        pileupDT           = cms.untracked.string("mm"),
        lepton     = cms.untracked.string("muon"),
        JEC        = cms.untracked.double(1)
)
process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC             = cms.untracked.string("S10"),
        pileupDT             = cms.untracked.string("mm"),
        lepton       = cms.untracked.string("muon"),
        JEC          = cms.untracked.double(-1)
)
process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        usePartonFlavour = cms.untracked.bool(True)
)
process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        usePartonFlavour = cms.untracked.bool(True)
)
process.demoElePur = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True)
)
process.demoMuoPur = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True)
)
process.demoEleDR = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        useDeltaR = cms.untracked.bool(True)
)
process.demoMuoDR = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        useDeltaR = cms.untracked.bool(True)
)
process.demoEleJerUp = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1)
)       
process.demoEleJerDown = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1)
)       
process.demoMuoJerUp = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)       
process.demoMuoJerDown = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)       
process.demoEleDump = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron")
)
process.demoEleDumpPup = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pup")
)
process.demoEleDumpPum = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pum")
)
process.demoEleDumpUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(1)
)
process.demoEleDumpDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(-1)
)
process.demoEleDumpPur = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True)
)
process.demoEleDumpDR = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        useDeltaR = cms.untracked.bool(True)
)
process.demoEleDumpJerUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1)
)
process.demoEleDumpJerDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1)
)
process.demoMuoDump = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon")
)
process.demoMuoDumpPup = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pup")
)
process.demoMuoDumpPum = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pum")
)
process.demoMuoDumpUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(1)
)
process.demoMuoDumpDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(-1)
)
process.demoMuoDumpPur = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True)
)
process.demoMuoDumpDR = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        useDeltaR = cms.untracked.bool(True)
)
process.demoMuoDumpJerUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)
process.demoMuoDumpJerDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)

process.demoEle1b = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(1)
)
process.demoEle2b = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(2)
)
process.demoMuo1b = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        numB = cms.untracked.double(1)
)
process.demoMuo2b = cms.EDProducer('ZbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        numB = cms.untracked.double(2)
)
process.demoEleGen1b = cms.EDProducer('GenbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        numB  = cms.untracked.double(1)
)
process.demoEleGen2b = cms.EDProducer('GenbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        numB  = cms.untracked.double(2)
)
process.demoMuoGen1b = cms.EDProducer('GenbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        numB  = cms.untracked.double(1)
)
process.demoMuoGen2b = cms.EDProducer('GenbAnalyzer',
        path = cms.untracked.string("."),
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        numB  = cms.untracked.double(2)
)
process.demoEleDump1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        numB         = cms.untracked.double(1)
)
process.demoEleDump2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        numB         = cms.untracked.double(2)
)
process.demoMuoDump1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        numB         = cms.untracked.double(1)
)
process.demoMuoDump2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        numB         = cms.untracked.double(2)
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
   process.regressedElectrons *
   process.calibratedElectrons *
   process.matchedElectrons *
   process.zeleMatchedeleMatched *
   process.selectedTriggeredPatMuonsEM *
   process.selectedPatMuonsTriggerMatchEM *
   process.matchedMuons0EM *
   process.matchedMuonsEM *
   process.selectedTriggeredPatElectronsEM *
   process.selectedPatElectronsTriggerMatchEM *
   process.matchedElectronsEM *
   process.zeleMatchedmuMatched *
   process.demoEleGen * process.demoMuoGen *
#   process.MyProcess *
   process.demoEle *
   process.demoEleBtag *
   process.demoElePum * process.demoElePup *
   process.demoEleUp * process.demoEleDown *
   process.demoElePur *
   process.demoEleDR *
   process.demoEleJerUp * process.demoEleJerDown *
   process.demoMuo *
   process.demoMuoBtag *
   process.demoMuoPum * process.demoMuoPup *
   process.demoMuoUp * process.demoMuoDown *
   process.demoMuoPur *
   process.demoMuoDR *
   process.demoMuoJerUp * process.demoMuoJerDown *
   process.demoEleDump *
   process.demoEleDumpPup * process.demoEleDumpPum *
   process.demoEleDumpUp * process.demoEleDumpDown *
   process.demoEleDumpPur *
   process.demoEleDumpDR *
   process.demoEleDumpJerUp * process.demoEleDumpJerDown *
   process.demoMuoDump *
   process.demoMuoDumpPup * process.demoMuoDumpPum *
   process.demoMuoDumpUp * process.demoMuoDumpDown *
   process.demoMuoDumpPur *
   process.demoMuoDumpDR *
   process.demoMuoDumpJerUp * 
   process.demoMuoDumpJerDown *
   process.demoEle1b *
   process.demoEle2b *
   process.demoMuo1b *
   process.demoMuo2b *
   process.demoEleGen1b *
   process.demoEleGen2b *
   process.demoMuoGen1b *
   process.demoMuoGen2b *
   process.demoEleDump1b *
   process.demoEleDump2b *
   process.demoMuoDump1b *
   process.demoMuoDump2b
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
	'keep *_addPileupInfo_*_*',
	'keep *_matchedElectrons_*_*',
        'keep *_regressedElectrons_*_*',
        'keep *_calibratedElectrons_*_*',
	'keep *_matchedMuons_*_*',
	'keep *_matchedElectronsEM_*_*',
	'keep *_matchedMuonsEM_*_*',
	'keep *_goodJets_*_*',
	'keep *_zeleMatchedeleMatched_*_*',
	'keep *_zmuMatchedmuMatched_*_*',
	'keep *_zeleMatchedmuMatched_*_*',
	'keep *_goodOfflinePrimaryVertices_*_*'
]

