import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = cms.untracked.string("WARNING")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
	# replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
		#'file:patTuple_1_1_IPe.root'
                'file:patTuple_2800_1_ma6.root'
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

process.demoElePum1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pum"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(1)
)

process.demoElePum2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pum"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(2)
)

process.demoElePup = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("ee_pup"),
	lepton  = cms.untracked.string("electron")
)

process.demoElePup1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pup"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(1)
)

process.demoElePup2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee_pup"),
        lepton  = cms.untracked.string("electron"),
        numB = cms.untracked.double(2)
)

process.demoEleUp = cms.EDProducer('ZbAnalyzer',
	pileupMC     = cms.untracked.string("S10"),
	pileupDT     = cms.untracked.string("ee"),
	lepton     = cms.untracked.string("electron"),
	JEC        = cms.untracked.double(1)
)

process.demoEleUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC     = cms.untracked.string("S10"),
        pileupDT     = cms.untracked.string("ee"),
        lepton     = cms.untracked.string("electron"),
        JEC        = cms.untracked.double(1),
        numB = cms.untracked.double(1)
)

process.demoEleUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC     = cms.untracked.string("S10"),
        pileupDT     = cms.untracked.string("ee"),
        lepton     = cms.untracked.string("electron"),
        JEC        = cms.untracked.double(1),
        numB = cms.untracked.double(2)
)

process.demoEleDown = cms.EDProducer('ZbAnalyzer',
	pileupMC       = cms.untracked.string("S10"),
	pileupDT       = cms.untracked.string("ee"),
	lepton 	     = cms.untracked.string("electron"),
	JEC    	     = cms.untracked.double(-1)
)

process.demoEleDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC       = cms.untracked.string("S10"),
        pileupDT       = cms.untracked.string("ee"),
        lepton       = cms.untracked.string("electron"),
        JEC          = cms.untracked.double(-1),
        numB = cms.untracked.double(1)
)

process.demoEleDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC       = cms.untracked.string("S10"),
        pileupDT       = cms.untracked.string("ee"),
        lepton       = cms.untracked.string("electron"),
        JEC          = cms.untracked.double(-1),
        numB = cms.untracked.double(2)
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

process.demoMuoPum1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pum"),
        lepton  = cms.untracked.string("muon"),
        numB = cms.untracked.double(1)
)

process.demoMuoPum2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pum"),
        lepton  = cms.untracked.string("muon"),
        numB = cms.untracked.double(2)
)

process.demoMuoPup = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm_pup"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuoPup1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pup"),
        lepton  = cms.untracked.string("muon"),
        numB = cms.untracked.double(1)
)

process.demoMuoPup2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm_pup"),
        lepton  = cms.untracked.string("muon"),
        numB = cms.untracked.double(2)
)

process.demoMuoUp = cms.EDProducer('ZbAnalyzer',
	pileupMC 	   = cms.untracked.string("S10"),
	pileupDT 	   = cms.untracked.string("mm"),
	lepton 	   = cms.untracked.string("muon"),
	JEC    	   = cms.untracked.double(1)
)

process.demoMuoUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC           = cms.untracked.string("S10"),
        pileupDT           = cms.untracked.string("mm"),
        lepton     = cms.untracked.string("muon"),
        JEC        = cms.untracked.double(1),
        numB = cms.untracked.double(1)
)

process.demoMuoUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC           = cms.untracked.string("S10"),
        pileupDT           = cms.untracked.string("mm"),
        lepton     = cms.untracked.string("muon"),
        JEC        = cms.untracked.double(1),
        numB = cms.untracked.double(2)
)

process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC 	     = cms.untracked.string("S10"),
	pileupDT 	     = cms.untracked.string("mm"),
	lepton 	     = cms.untracked.string("muon"),
	JEC    	     = cms.untracked.double(-1)
)

process.demoMuoDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC             = cms.untracked.string("S10"),
        pileupDT             = cms.untracked.string("mm"),
        lepton       = cms.untracked.string("muon"),
        JEC          = cms.untracked.double(-1),
        numB = cms.untracked.double(1)
)

process.demoMuoDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC             = cms.untracked.string("S10"),
        pileupDT             = cms.untracked.string("mm"),
        lepton       = cms.untracked.string("muon"),
        JEC          = cms.untracked.double(-1),
        numB = cms.untracked.double(2)
)

process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoEleBtag1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        usePartonFlavour = cms.untracked.bool(True),
        numB = cms.untracked.double(1)
)

process.demoEleBtag2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        usePartonFlavour = cms.untracked.bool(True),
        numB = cms.untracked.double(2)
)

process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        usePartonFlavour = cms.untracked.bool(True),
        numB = cms.untracked.double(1)
)

process.demoMuoBtag2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        usePartonFlavour = cms.untracked.bool(True),
        numB = cms.untracked.double(2)
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

process.demoElePur1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True),
        numB  = cms.untracked.double(1)
)

process.demoElePur2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True),
        numB  = cms.untracked.double(2)
)

process.demoMuoPur = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True)
)

process.demoMuoPur1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True),
        numB  = cms.untracked.double(1)
)

process.demoMuoPur2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True),
        numB  = cms.untracked.double(2)
)

process.demoEleJerUp = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1)
)

process.demoEleJerUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1),
        numB  = cms.untracked.double(1)
)

process.demoEleJerUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1),
        numB  = cms.untracked.double(2)
)

process.demoEleJerDown = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1)
)

process.demoEleJerDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1),
        numB  = cms.untracked.double(1)
)

process.demoEleJerDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1),
        numB  = cms.untracked.double(2)
)

process.demoMuoJerUp = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)

process.demoMuoJerUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1),
        numB  = cms.untracked.double(1)
)

process.demoMuoJerUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1),
        numB  = cms.untracked.double(2)
)

process.demoMuoJerDown = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)

process.demoMuoJerDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1),
        numB  = cms.untracked.double(1)
)

process.demoMuoJerDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1),
        numB  = cms.untracked.double(2)
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

process.demoEleDumpPup1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pup"),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpPup2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pup"),
        numB         = cms.untracked.double(2)
)

process.demoEleDumpPum = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pum")
)

process.demoEleDumpPum1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pum"),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpPum2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pileupDT = cms.untracked.string("ee_pum"),
        numB         = cms.untracked.double(2)
)

process.demoEleDumpUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(1)
)

process.demoEleDumpUp1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(1),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpUp2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(1),
        numB         = cms.untracked.double(2)
)

process.demoEleDumpDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(-1)
)

process.demoEleDumpDown1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(-1),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpDown2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(-1),
        numB         = cms.untracked.double(2)
)

process.demoEleDumpPur = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True)
)

process.demoEleDumpPur1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpPur2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True),
        numB         = cms.untracked.double(2)
)

process.demoEleDumpJerUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1)
)

process.demoEleDumpJerUp1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpJerUp2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1),
        numB         = cms.untracked.double(2)
)

process.demoEleDumpJerDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1)
)

process.demoEleDumpJerDown1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1),
        numB         = cms.untracked.double(1)
)

process.demoEleDumpJerDown2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1),
        numB         = cms.untracked.double(2)
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

process.demoMuoDumpPup1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pup"),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpPup2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pup"),
        numB         = cms.untracked.double(2)
)

process.demoMuoDumpPum = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pum")
)

process.demoMuoDumpPum1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pum"),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpPum2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pileupDT = cms.untracked.string("mm_pum"),
        numB         = cms.untracked.double(2)
)

process.demoMuoDumpUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(1)
)

process.demoMuoDumpUp1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(1),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpUp2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(1),
        numB         = cms.untracked.double(2)
)

process.demoMuoDumpDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(-1)
)

process.demoMuoDumpDown1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(-1),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpDown2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(-1),
        numB         = cms.untracked.double(2)
)

process.demoMuoDumpPur = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True)
)

process.demoMuoDumpPur1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpPur2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True),
        numB         = cms.untracked.double(2)
)

process.demoMuoDumpJerUp = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1)
)

process.demoMuoDumpJerUp1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpJerUp2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1),
        numB         = cms.untracked.double(2)
)

process.demoMuoDumpJerDown = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1)
)

process.demoMuoDumpJerDown1b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1),
        numB         = cms.untracked.double(1)
)

process.demoMuoDumpJerDown2b = cms.EDAnalyzer('ZbDumper',
        lepton       = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1),
        numB         = cms.untracked.double(2)
)
process.p = cms.Path(
process.demoEle*process.demoEle1b*process.demoEle2b*
process.demoElePum*process.demoElePum1b*process.demoElePum2b*
process.demoElePup*process.demoElePup1b*process.demoElePup2b*
process.demoEleUp*process.demoEleUp1b*process.demoEleUp2b*
process.demoEleDown*process.demoEleDown1b*process.demoEleDown2b*
process.demoMuo*process.demoMuo1b*process.demoMuo2b*
process.demoMuoPum*process.demoMuoPum1b*process.demoMuoPum2b*
process.demoMuoPup*process.demoMuoPup1b*process.demoMuoPup2b*
process.demoMuoUp*process.demoMuoUp1b*process.demoMuoUp2b*
process.demoMuoDown*process.demoMuoDown1b*process.demoMuoDown2b*
process.demoEleBtag*process.demoEleBtag1b*process.demoEleBtag2b*
process.demoMuoBtag*process.demoMuoBtag1b*process.demoMuoBtag2b*
process.demoEle2*process.demoMuo2*
process.demoEleGen*process.demoEleGen1b*process.demoEleGen2b*
process.demoMuoGen*process.demoMuoGen1b*process.demoMuoGen2b*
process.demoElePur*process.demoElePur1b*process.demoElePur2b*
process.demoMuoPur*process.demoMuoPur1b*process.demoMuoPur2b*
process.demoEleJerUp*process.demoEleJerUp1b*process.demoEleJerUp2b*
process.demoEleJerDown*process.demoEleJerDown1b*process.demoEleJerDown2b*
process.demoMuoJerUp*process.demoMuoJerUp1b*process.demoMuoJerUp2b*
process.demoMuoJerDown*process.demoMuoJerDown1b*process.demoMuoJerDown2b*
process.demoEleDump*process.demoEleDump1b*process.demoEleDump2b*
process.demoMuoDump*process.demoMuoDump1b*process.demoMuoDump2b*
process.demoEleDumpPup*process.demoEleDumpPup1b*process.demoEleDumpPup2b*
process.demoEleDumpPum*process.demoEleDumpPum1b*process.demoEleDumpPum2b*
process.demoEleDumpUp*process.demoEleDumpUp1b*process.demoEleDumpUp2b*
process.demoEleDumpDown*process.demoEleDumpDown1b*process.demoEleDumpDown2b*
process.demoEleDumpPur*process.demoEleDumpPur1b*process.demoEleDumpPur2b*
process.demoEleDumpJerUp*process.demoEleDumpJerUp1b*process.demoEleDumpJerUp2b*
process.demoEleDumpJerDown*process.demoEleDumpJerDown1b*process.demoEleDumpJerDown2b*
process.demoMuoDumpPup*process.demoMuoDumpPup1b*process.demoMuoDumpPup2b*
process.demoMuoDumpPum*process.demoMuoDumpPum1b*process.demoMuoDumpPum2b*
process.demoMuoDumpUp*process.demoMuoDumpUp1b*process.demoMuoDumpUp2b*
process.demoMuoDumpDown*process.demoMuoDumpDown1b*process.demoMuoDumpDown2b*
process.demoMuoDumpPur*process.demoMuoDumpPur1b*process.demoMuoDumpPur2b*
process.demoMuoDumpJerUp*process.demoMuoDumpJerUp1b*process.demoMuoDumpJerUp2b*
process.demoMuoDumpJerDown*process.demoMuoDumpJerDown1b*process.demoMuoDumpJerDown2b
)
