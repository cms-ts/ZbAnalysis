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

process.filter = cms.EDFilter('ZbFilter'
)

process.demo_ee = cms.EDAnalyzer('ZbAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demo_ee_up = cms.EDAnalyzer('ZbAnalyzer',
	pileup     = cms.untracked.string("S10"),
	lepton     = cms.untracked.string("electron"),
	JEC        = cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demo_ee_down = cms.EDAnalyzer('ZbAnalyzer',
	pileup       = cms.untracked.string("S10"),
	lepton 	     = cms.untracked.string("electron"),
	JEC    	     = cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demo_mm = cms.EDAnalyzer('ZbAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demo_mm_up = cms.EDAnalyzer('ZbAnalyzer',
	pileup 	   = cms.untracked.string("S10"),
	lepton 	   = cms.untracked.string("muon"),
	JEC    	   = cms.untracked.double(1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demo_mm_down = cms.EDAnalyzer('ZbAnalyzer',
	pileup 	     = cms.untracked.string("S10"),
	lepton 	     = cms.untracked.string("muon"),
	JEC    	     = cms.untracked.double(-1),
	usePartonFlavour = cms.untracked.bool(False)
)

process.demo_ee_btag = cms.EDAnalyzer('ZbAnalyzer',
        pileup  = cms.untracked.string("S10"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demo_mm_btag = cms.EDAnalyzer('ZbAnalyzer',
        pileup  = cms.untracked.string("S10"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(0),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demo2_ee = cms.EDAnalyzer('ZJetsAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("electron"),
)

process.demo2_mm = cms.EDAnalyzer('ZJetsAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("muon"),
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ZbTree.root')
)

process.p = cms.Path(process.filter*process.demo_ee*process.demo_ee_up*process.demo_ee_down*process.demo_mm*process.demo_mm_up*process.demo_mm_down*process.demo_ee_btag*process.demo_mm_btag*process.demo2_ee*process.demo2_mm)

