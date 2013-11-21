import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
		'file:patTuple_1_1_IPe.root'
        )
)

process.demoEle = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoElePum = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pum"),
	lepton  = cms.untracked.string("electron")
)

process.demoElePup = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pup"),
	lepton  = cms.untracked.string("electron")
)

process.demoEleUp = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(1)
)

process.demoEleDown = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton 	= cms.untracked.string("electron"),
	JEC     = cms.untracked.double(-1)
)

process.demoMuo = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoPum = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pum"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoPup = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pup"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton 	= cms.untracked.string("muon"),
	JEC    	= cms.untracked.double(1)
)

process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(-1)
)

process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoEle2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoMuo2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoEleMuo = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(0)
)

process.demoEleMuoUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(1)
)

process.demoEleMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC   = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuoPum = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pum"),
        lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuoPup = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pup"),
        lepton  = cms.untracked.string("electron+muon")
)

process.demoElePur = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	pcut = cms.untracked.bool(True)
)

process.demoMuoPur = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	pcut = cms.untracked.bool(True)
)

process.demoEleMuoPur = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	pcut = cms.untracked.bool(True)
)

process.demoEleDR = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoMuoDR = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoEleMuoDR = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	useDeltaR = cms.untracked.bool(True)
)

process.demoEleJerUp = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)

process.demoEleJerDown = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)

process.demoMuoJerUp = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)

process.demoMuoJerDown = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)

process.demoEleMuoJerUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JER     = cms.untracked.double(1)
)

process.demoEleMuoJerDown = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JER     = cms.untracked.double(-1)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ZbTree.root')
)
process.p = cms.Path(process.demoEle*process.demoElePum*process.demoElePup*process.demoEleUp*process.demoEleDown*process.demoMuo*process.demoMuoPum*process.demoMuoPup*process.demoMuoUp*process.demoMuoDown*process.demoEleBtag*process.demoMuoBtag*process.demoEle2*process.demoMuo2*process.demoEleMuo*process.demoEleMuoUp*process.demoEleMuoDown*process.demoEleMuoPum*process.demoEleMuoPup*process.demoElePur*process.demoMuoPur*process.demoEleMuoPur*process.demoEleDR*process.demoMuoDR*process.demoEleMuoDR*process.demoEleJerUp*process.demoEleJerDown*process.demoMuoJerUp*process.demoMuoJerDown*process.demoEleMuoJerUp*process.demoEleMuoJerDown)
