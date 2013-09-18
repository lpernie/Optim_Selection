import FWCore.ParameterSet.Config as cms

isGun = True

if isGun:
   correctHits = False
   useHLTFilter = False
   process = cms.Process('PI0DUMPERGUN')
else:
   correctHits = True
   useHLTFilter = True
   process = cms.Process('PI0DUMPER')

if correctHits:
    print 'CORRECTING HITS'
    
if useHLTFilter:
    print 'FILTERING PI0 EVENTS'

if useHLTFilter:
    import copy
    from HLTrigger.HLTfilters.hltHighLevel_cfi import *
    process.AlcaP0Filter = copy.deepcopy(hltHighLevel)
    process.AlcaP0Filter.throw = cms.bool(False)
    process.AlcaP0Filter.HLTPaths = ['AlCa_EcalPi0_*']

import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi

from Geometry.CaloEventSetup.CaloTopology_cfi import *
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Geometry.TrackerGeometryBuilder.trackerGeometry_cfi')
process.load('RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi')
process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#load transparency loss
process.GlobalTag.globaltag = 'GR_R_53_V15::All'

process.options = cms.untracked.PSet( 
    wantSummary = cms.untracked.bool(False),
    SkipEvent   = cms.untracked.vstring('ProductNotFound')
) 

process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMenuConfig_cff')
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtStableParametersConfig_cff')
if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(False),
        doIntercalib = cms.bool(False),
        doLaserCorrections = cms.bool(True),
        EBRecHitCollection = cms.InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsEB', 'HLT'),
        EERecHitCollection = cms.InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsEE', 'HLT'),
        EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
        EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
    )

#NewPi0Dumper_Gun
process.load('Analysis.Modules.NewPi0Dumper_Gun_cfi')
if isGun:
   process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('/afs/cern.ch/work/l/lpernie/pi0/HLT_pi0/CMSSW_5_3_3_patch2/src/Analysis/Modules/test/HLT_test.root')
   process.newPi0Dumper_Gun.outfilename_ContCorr = cms.untracked.string('/tmp/ContCorr_0.root')
else:   
   process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('/afs/cern.ch/work/l/lpernie/pi0/HLT_pi0/CMSSW_5_3_3_patch2/src/Analysis/Modules/test/HLT_test.root')
process.newPi0Dumper_Gun.ExternalGeometry = cms.untracked.string('/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/CalibCode/submit/common/caloGeometry.root')
#Cuts for AlcaPi0
process.newPi0Dumper_Gun.useES = cms.untracked.bool(True)
if isGun:
   process.newPi0Dumper_Gun.StoreMCTruth = cms.untracked.bool(True)
else:
   process.newPi0Dumper_Gun.StoreMCTruth = cms.untracked.bool(False)
process.newPi0Dumper_Gun.OnlyContCorr = cms.untracked.bool(False)
process.newPi0Dumper_Gun.ptpi0Cut = 0.7
process.newPi0Dumper_Gun.s1CluCutEE = 0.5
process.newPi0Dumper_Gun.s1CluCut = 0.35
process.newPi0Dumper_Gun.ptCluCut = 0.35
process.newPi0Dumper_Gun.s4s9CluCut = 0.8
process.newPi0Dumper_Gun.DoOffGeom = cms.untracked.bool(False)
process.newPi0Dumper_Gun.StoreTransparencyCorrection = cms.untracked.bool(True)
### choosing proper input tag (recalibration module changes the collection names)
if( not(isGun) and correctHits) :
    process.newPi0Dumper_Gun.EBRecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEB')
    process.newPi0Dumper_Gun.EERecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEE')
    process.newPi0Dumper_Gun.ESRecHitCollectionTag = cms.untracked.InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES', 'HLT')
if isGun:
    process.newPi0Dumper_Gun.EBRecHitCollectionTag = cms.untracked.InputTag('ecalRecHit','EcalRecHitsEB', 'RECO')
    process.newPi0Dumper_Gun.EERecHitCollectionTag = cms.untracked.InputTag('ecalRecHit','EcalRecHitsEE', 'RECO')
    process.newPi0Dumper_Gun.ESRecHitCollectionTag = cms.untracked.InputTag('ecalPreshowerRecHit','EcalRecHitsES','RECO')

process.newPi0Dumper_Gun.preshRecHitProducer = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES")
process.newPi0Dumper_Gun.PFRecHitCollectionTag = cms.untracked.InputTag('particleFlowClusterECAL','','RECO')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.p = cms.Path()

#use HLT
if useHLTFilter:
    process.p *= process.AlcaP0Filter

if correctHits :
    print 'ADDING RECALIB RECHIT MODULE WITH PARAMETERS'
    print 'ENERGY SCALE '+str(process.ecalPi0ReCorrected.doEnergyScale)
    print 'INTERCALIBRATION '+str(process.ecalPi0ReCorrected.doIntercalib)
    print 'LASER '+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected

#use beam spot
if (process.newPi0Dumper_Gun.useBeamSpotPosition == True):
    process.p *= process.offlineBeamSpot

#build ntuple
process.p *= process.newPi0Dumper_Gun

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/data/Run2012A/MinimumBias/RECO/PromptReco-v1/000/190/913/3C2E8833-A385-E111-BBEE-E0CB4E553651.root',
        'root://eoscms//eos/cms/store/data/Run2012A/MinimumBias/RECO/PromptReco-v1/000/190/456/18309A2F-F880-E111-A9D5-001D09F24353.root', 
        'root://eoscms//eos/cms/store/data/Run2012A/MinimumBias/RECO/PromptReco-v1/000/190/456/36AD3832-FD80-E111-82AD-BCAEC518FF8F.root',
        'root://eoscms//eos/cms/store/data/Run2012A/MinimumBias/RECO/PromptReco-v1/000/191/201/760DC119-2887-E111-8FE6-5404A63886CC.root'
        #'root://eoscms//eos/cms', 
   )
)
