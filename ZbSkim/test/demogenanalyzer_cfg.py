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
	pileupMC     = cms.untracked.string("S10"),
	pileupDT     = cms.untracked.string("ee"),
	lepton     = cms.untracked.string("electron"),
	JEC        = cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoEleDown = cms.EDProducer('ZbAnalyzer',
	pileupMC       = cms.untracked.string("S10"),
	pileupDT       = cms.untracked.string("ee"),
	lepton 	     = cms.untracked.string("electron"),
	JEC    	     = cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuo = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("mm"),
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
	pileupMC 	   = cms.untracked.string("S10"),
	pileupDT 	   = cms.untracked.string("mm"),
	lepton 	   = cms.untracked.string("muon"),
	JEC    	   = cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC 	     = cms.untracked.string("S10"),
	pileupDT 	     = cms.untracked.string("mm"),
	lepton 	     = cms.untracked.string("muon"),
	JEC    	     = cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoEle2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("electron")
)

process.demoMuo2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("muon")
)

process.demoEleGen = cms.EDProducer('GenbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoMuoGen = cms.EDProducer('GenbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ZbTree.root')
)

process.demoEleDump = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron")
)

process.demoMuoDump = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon")
)

process.p = cms.Path(process.demoEle*process.demoElePum*process.demoElePup*process.demoEleUp*process.demoEleDown*process.demoMuo*process.demoMuoPum*process.demoMuoPup*process.demoMuoUp*process.demoMuoDown*process.demoEleBtag*process.demoMuoBtag*process.demoEle2*process.demoMuo2*process.demoEleGen*process.demoMuoGen*process.demoEleDump*process.demoMuoDump)

