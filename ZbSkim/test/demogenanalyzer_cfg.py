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

process.demoEle1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(1)
)

process.demoEle2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(2)
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
	pileupMC     = cms.untracked.string("S10"),
	pileupDT     = cms.untracked.string("ee"),
	lepton     = cms.untracked.string("electron"),
	JEC        = cms.untracked.double(1)
)

process.demoEleDown = cms.EDProducer('ZbAnalyzer',
	pileupMC       = cms.untracked.string("S10"),
	pileupDT       = cms.untracked.string("ee"),
	lepton 	     = cms.untracked.string("electron"),
	JEC    	     = cms.untracked.double(-1)
)

process.demoMuo = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuo1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	numB = cms.untracked.double(1)	
)

process.demoMuo2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	numB = cms.untracked.double(2)
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
	pileupMC 	   = cms.untracked.string("S10"),
	pileupDT 	   = cms.untracked.string("mm"),
	lepton 	   = cms.untracked.string("muon"),
	JEC    	   = cms.untracked.double(1)
)

process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC 	     = cms.untracked.string("S10"),
	pileupDT 	     = cms.untracked.string("mm"),
	lepton 	     = cms.untracked.string("muon"),
	JEC    	     = cms.untracked.double(-1)
)

process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoEle2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoMuo2 = cms.EDAnalyzer('ZJetsAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoEleGen = cms.EDProducer('GenbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron")
)

process.demoEleGen1b = cms.EDProducer('GenbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	numB  = cms.untracked.double(1)
)

process.demoEleGen2b = cms.EDProducer('GenbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	numB  = cms.untracked.double(2)
)

process.demoMuoGen = cms.EDProducer('GenbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoGen1b = cms.EDProducer('GenbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	numB  = cms.untracked.double(1)
)

process.demoMuoGen2b = cms.EDProducer('GenbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	numB  = cms.untracked.double(2)
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

process.demoEleJerUp = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1)
)

process.demoEleJerDown = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
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

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ZbTree.root')
)

process.demoEleDump = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron")
)

process.demoEleDump1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        numB         = cms.untracked.double(1)
)

process.demoEleDump2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        numB         = cms.untracked.double(2)
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

process.demoMuoDump1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        numB         = cms.untracked.double(1)
)

process.demoMuoDump2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        numB         = cms.untracked.double(2)
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


process.p = cms.Path(process.demoEle*process.demoEle1b*process.demoEle2b*process.demoElePum*process.demoElePup*process.demoEleUp*process.demoEleDown*process.demoMuo*process.demoMuo1b*process.demoMuo2b*process.demoMuoPum*process.demoMuoPup*process.demoMuoUp*process.demoMuoDown*process.demoEleBtag*process.demoMuoBtag*process.demoEle2*process.demoMuo2*process.demoEleGen*process.demoEleGen1b*process.demoEleGen2b*process.demoMuoGen*process.demoMuoGen1b*process.demoMuoGen2b*process.demoElePur*process.demoMuoPur*process.demoEleDR*process.demoMuoDR*process.demoEleJerUp*process.demoEleJerDown*process.demoMuoJerUp*process.demoMuoJerDown*process.demoEleDump*process.demoEleDump1b*process.demoEleDump2b*process.demoMuoDump*process.demoMuoDump1b*process.demoMuoDump2b*process.demoEleDumpPup*process.demoEleDumpPum*process.demoEleDumpPum*process.demoEleDumpUp*process.demoEleDumpDown*process.demoEleDumpPur*process.demoEleDumpDR*process.demoEleDumpJerUp*process.demoEleDumpJerDown*process.demoMuoDump*process.demoMuoDumpPup*process.demoMuoDumpPum*process.demoMuoDumpUp*process.demoMuoDumpDown*process.demoMuoDumpPur*process.demoMuoDumpDR*process.demoMuoDumpJerUp*process.demoMuoDumpJerDown)
