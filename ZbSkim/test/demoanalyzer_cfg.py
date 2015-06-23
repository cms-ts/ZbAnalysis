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

process.genBDW = cms.EDProducer('GenBWeightProducer',
	pprop = cms.FileInPath("ZbAnalysis/ZbSkim/data/bHProp.txt"),
	p411 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_411_decaytable.txt"),
	p421 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_421_decaytable.txt"),
	p431 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_431_decaytable.txt"),
	p441 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_441_decaytable.txt"),
	p511 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_511_decaytable.txt"),
	p521 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_521_decaytable.txt"),
	p531 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_531_decaytable.txt"),
	p541 = cms.FileInPath("ZbAnalysis/ZbSkim/data/fHadron_541_decaytable.txt")
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
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton  = cms.untracked.string("electron"),
	JEC     = cms.untracked.double(1)
)

process.demoEleUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(1),
	numB = cms.untracked.double(1)
)

process.demoEleUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(1),
	numB = cms.untracked.double(2)
)

process.demoEleDown = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
	lepton 	= cms.untracked.string("electron"),
	JEC     = cms.untracked.double(-1)
)

process.demoEleDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(-1),
	numB = cms.untracked.double(1)
)

process.demoEleDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JEC     = cms.untracked.double(-1),
	numB = cms.untracked.double(2)
)

process.demoMuo = cms.EDProducer('ZbAnalyzer',
	pileupMC = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon")
)

process.demoMuo1b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	numB = cms.untracked.double(1)
)

process.demoMuo2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
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
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton 	= cms.untracked.string("muon"),
	JEC    	= cms.untracked.double(1)
)

process.demoMuoUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(1),
	numB = cms.untracked.double(1)
)

process.demoMuoUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(1),
	numB = cms.untracked.double(2)
)

process.demoMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
	lepton  = cms.untracked.string("muon"),
	JEC     = cms.untracked.double(-1)
)

process.demoMuoDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(-1),
	numB = cms.untracked.double(1)
)

process.demoMuoDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JEC     = cms.untracked.double(-1),
	numB = cms.untracked.double(2)
)

process.demoEleBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoEleBtag1b = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        usePartonFlavour = cms.untracked.bool(True),
	numB = cms.untracked.double(1)
)

process.demoEleBtag2b = cms.EDProducer('ZbAnalyzer',
        pileupMC    = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        usePartonFlavour = cms.untracked.bool(True),
	numB = cms.untracked.double(2)
)

process.demoMuoBtag = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
	usePartonFlavour = cms.untracked.bool(True)
)

process.demoMuoBtag1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        usePartonFlavour = cms.untracked.bool(True),
	numB = cms.untracked.double(1)
)

process.demoMuoBtag2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        usePartonFlavour = cms.untracked.bool(True),
	numB = cms.untracked.double(2)
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
	lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuo1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
	numB = cms.untracked.double(1)
)

process.demoEleMuo2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
	numB = cms.untracked.double(2)
)

process.demoEleMuoUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(1)
)

process.demoEleMuoUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JEC     = cms.untracked.double(1),
	numB = cms.untracked.double(1)
)

process.demoEleMuoUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JEC     = cms.untracked.double(1),
	numB = cms.untracked.double(2)
)

process.demoEleMuoDown = cms.EDProducer('ZbAnalyzer',
	pileupMC   = cms.untracked.string("S10"),
	pileupDT = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JEC     = cms.untracked.double(-1)
)

process.demoEleMuoDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JEC     = cms.untracked.double(-1),
	numB = cms.untracked.double(1)
)

process.demoEleMuoDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JEC     = cms.untracked.double(-1),
	numB = cms.untracked.double(2)
)

process.demoEleMuoPum = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pum"),
        lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuoPum1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pum"),
        lepton  = cms.untracked.string("electron+muon"),
	numB = cms.untracked.double(1)
)

process.demoEleMuoPum2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pum"),
        lepton  = cms.untracked.string("electron+muon"),
	numB = cms.untracked.double(2)
)

process.demoEleMuoPup = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pup"),
        lepton  = cms.untracked.string("electron+muon")
)

process.demoEleMuoPup1b = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pup"),
        lepton  = cms.untracked.string("electron+muon"),
	numB = cms.untracked.double(1)
)

process.demoEleMuoPup2b = cms.EDProducer('ZbAnalyzer',
        pileupMC   = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("em_pup"),
        lepton  = cms.untracked.string("electron+muon"),
	numB = cms.untracked.double(2)
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
	numB = cms.untracked.double(1)
)

process.demoElePur2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        pcut = cms.untracked.bool(True),
	numB = cms.untracked.double(2)
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
	numB = cms.untracked.double(1)
)

process.demoMuoPur2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        pcut = cms.untracked.bool(True),
	numB = cms.untracked.double(2)
)

process.demoEleMuoPur = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	pcut = cms.untracked.bool(True)
)

process.demoEleMuoPur1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        pcut = cms.untracked.bool(True),
	numB = cms.untracked.double(1)
)

process.demoEleMuoPur2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        pcut = cms.untracked.bool(True),
	numB = cms.untracked.double(2)
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
	numB = cms.untracked.double(1)
)

process.demoEleJerUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(1),
	numB = cms.untracked.double(2)
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
	numB = cms.untracked.double(1)
)

process.demoEleJerDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("ee"),
        lepton  = cms.untracked.string("electron"),
        JER     = cms.untracked.double(-1),
	numB = cms.untracked.double(2)
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
	numB = cms.untracked.double(1)
)

process.demoMuoJerUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(1),
	numB = cms.untracked.double(2)
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
	numB = cms.untracked.double(1)
)

process.demoMuoJerDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC = cms.untracked.string("S10"),
        pileupDT = cms.untracked.string("mm"),
        lepton  = cms.untracked.string("muon"),
        JER     = cms.untracked.double(-1),
	numB = cms.untracked.double(2)
)

process.demoEleMuoJerUp = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JER     = cms.untracked.double(1)
)

process.demoEleMuoJerUp1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JER     = cms.untracked.double(1),
	numB = cms.untracked.double(1)
)

process.demoEleMuoJerUp2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JER     = cms.untracked.double(1),
	numB = cms.untracked.double(2)
)

process.demoEleMuoJerDown = cms.EDProducer('ZbAnalyzer',
	pileupMC  = cms.untracked.string("S10"),
	pileupDT  = cms.untracked.string("em"),
	lepton  = cms.untracked.string("electron+muon"),
	JER     = cms.untracked.double(-1)
)

process.demoEleMuoJerDown1b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JER     = cms.untracked.double(-1),
	numB = cms.untracked.double(1)
)

process.demoEleMuoJerDown2b = cms.EDProducer('ZbAnalyzer',
        pileupMC  = cms.untracked.string("S10"),
        pileupDT  = cms.untracked.string("em"),
        lepton  = cms.untracked.string("electron+muon"),
        JER     = cms.untracked.double(-1),
	numB = cms.untracked.double(2)
)
#process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ZbTree.root')
)
process.p = cms.Path(
process.genBDW*
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
process.demoEleMuo*process.demoEleMuo1b*process.demoEleMuo2b*
process.demoEleMuoUp*process.demoEleMuoUp1b*process.demoEleMuoUp2b*
process.demoEleMuoDown*process.demoEleMuoDown1b*process.demoEleMuoDown2b*
process.demoEleMuoPum*process.demoEleMuoPum1b*process.demoEleMuoPum2b*
process.demoEleMuoPup*process.demoEleMuoPup1b*process.demoEleMuoPup2b*
process.demoElePur*process.demoElePur1b*process.demoElePur2b*
process.demoMuoPur*process.demoMuoPur1b*process.demoMuoPur2b*
process.demoEleMuoPur*process.demoEleMuoPur1b*process.demoEleMuoPur2b*
process.demoEleJerUp*process.demoEleJerUp1b*process.demoEleJerUp2b*
process.demoEleJerDown*process.demoEleJerDown1b*process.demoEleJerDown2b*
process.demoMuoJerUp*process.demoMuoJerUp1b*process.demoMuoJerUp2b*
process.demoMuoJerDown*process.demoMuoJerDown1b*process.demoMuoJerDown2b*
process.demoEleMuoJerUp*process.demoEleMuoJerUp1b*process.demoEleMuoJerUp2b*
process.demoEleMuoJerDown*process.demoEleMuoJerDown1b*process.demoEleMuoJerDown2b
)
