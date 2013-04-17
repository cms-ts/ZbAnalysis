import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
		'file:patTuple.root'
        )
)

process.filter = cms.EDFilter('ZbFilter'
)

process.demo_ee = cms.EDAnalyzer('ZbAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(0)
)

process.demo_ee_up = cms.EDAnalyzer('ZbAnalyzer',
	pileup     = cms.untracked.string("S10"),
	lepton     = cms.untracked.string("electron"),
	JEC        = cms.untracked.double(1)
)

process.demo_ee_down = cms.EDAnalyzer('ZbAnalyzer',
	pileup       = cms.untracked.string("S10"),
	lepton 	     = cms.untracked.string("electron"),
	JEC    	     = cms.untracked.double(-1)
)

process.demo_mm = cms.EDAnalyzer('ZbAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(0)
)

process.demo_mm_up = cms.EDAnalyzer('ZbAnalyzer',
	pileup 	   = cms.untracked.string("S10"),
	lepton 	   = cms.untracked.string("muon"),
	JEC    	   = cms.untracked.double(1)
)

process.demo_mm_down = cms.EDAnalyzer('ZbAnalyzer',
	pileup 	     = cms.untracked.string("S10"),
	lepton 	     = cms.untracked.string("muon"),
	JEC    	     = cms.untracked.double(-1)
)

process.demo2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileup = cms.untracked.string("S10")
)

process.demo2_ee = cms.EDAnalyzer('ZbAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("electron"),
)

process.demo2_mm = cms.EDAnalyzer('ZbAnalyzer',
	pileup  = cms.untracked.string("S10"),
	lepton  = cms.untracked.string("muon"),
)

process.bjettag = cms.EDAnalyzer('PatBJetTagAnalyzer',
	jets = cms.InputTag( "goodJets" ),
	jetPtCut = cms.double( 30. ),
	jetEtaCut = cms.double( 2.5 )
)

process.bjetvtx = cms.EDAnalyzer('PatBJetVertexAnalyzer',
        jets = cms.InputTag( "goodJets" ),
        jetPtCut = cms.double( 30. ),
        jetEtaCut = cms.double( 2.5 )
)

process.TFileService = cms.Service("TFileService",
		                  fileName = cms.string('ZbTree.root')
			          )

process.p = cms.Path(process.filter*process.demo_ee*process.demo_mm*process.demo2*process.demo_mm_up*process.demo_mm_down*process.demo_ee_up*process.demo_ee_down*process.demo2_ee*process.demo2_mm)
#process.p = cms.Path(process.filter*process.demo*process.bjettag*process.bjetvtx)
