import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
		'file:patTuple.root'
        )
)

process.filter = cms.EDFilter('ZbFilter'
)

process.demo = cms.EDAnalyzer('ZbAnalyzer'
)

process.TFileService = cms.Service("TFileService",
		                  fileName = cms.string('ZbTree.root')
			          )

process.p = cms.Path(process.filter*process.demo)
