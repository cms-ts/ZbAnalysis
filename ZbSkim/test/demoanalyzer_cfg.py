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
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoElePum = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pum"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoElePup = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pup"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoEleUp = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoEleDown = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton 	= cms.untracked.string("electron"),
	JEC     = cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuo = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuoPum = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pum"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuoPup = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pup"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuoUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton 	= cms.untracked.string("muon"),
	JEC    	= cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(0),
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
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)
process.demoEleMuoUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)
process.demoEleMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC   = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC    	= cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)
process.demoEleMuoPum = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pum"),
        lepton  = cms.untracked.string("electron+muon"),
        JEC     = cms.untracked.double(0),
        usePartonFlavour = cms.untracked.bool(False)
)
process.demoEleMuoPup = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pup"),
        lepton  = cms.untracked.string("electron+muon"),
        JEC     = cms.untracked.double(0),
        usePartonFlavour = cms.untracked.bool(False)
)

process.demoElePur = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False),
	pcut = cms.untracked.bool(True)
)

process.demoMuoPur = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False),
	pcut = cms.untracked.bool(True)
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ZbTree.root')
)
process.p = cms.Path(process.demoEle*process.demoElePum*process.demoElePup*process.demoEleUp*process.demoEleDown*process.demoMuo*process.demoMuoPum*process.demoMuoPup*process.demoMuoUp*process.demoMuoDown*process.demoEleBtag*process.demoMuoBtag*process.demoEle2*process.demoMuo2*process.demoEleMuo*process.demoEleMuoUp*process.demoEleMuoDown*process.demoEleMuoPum*process.demoEleMuoPup*process.demoElePur*process.demoMuoPur)

