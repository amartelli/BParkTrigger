import FWCore.ParameterSet.Config as cms

process = cms.Process('testBParkHLT')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('FWCore/MessageService/MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(
        ## 0PU
        #'file:/eos/cms/store/user/amartell/BPrk/MC_forHLT2021/BPH_RAW-RECO_forHLT_noTagMu_part0_noPU_valid.root',
        #'file:/eos/cms/store/user/amartell/BPrk/MC_forHLT2021/BPH_RAW-RECO_forHLT_noTagMu_part1_noPU_valid.root',

        ## withPU
        'file:/eos/cms/store/user/amartell/BPrk/MC_forHLT2021/BPH_RAW-RECO_forHLT_noTagMu_part0_pAll_wPU_valid.root',
        'file:/eos/cms/store/user/amartell/BPrk/MC_forHLT2021/BPH_RAW-RECO_forHLT_noTagMu_part1_pAll_wPU_valid.root',
    )
)


process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('BPark_HLTEff_out.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.myOutputTest = cms.EDAnalyzer(
    'BParkTriggerEfficiency',
    trackingTruth = cms.InputTag('mix:MergedTrackTruth:DIGI2RAW'),
    pixelTracks = cms.InputTag('pixelTracks::RECO'),
    pixelVertex = cms.InputTag('pixelVertices::RECO'),
    genParticles = cms.InputTag("genParticles::DIGI2RAW"),
    trackAssociation = cms.InputTag("trackingParticlePixelTrackAsssociation"),
    puSummary = cms.untracked.InputTag("addPileupInfo::DIGI2RAW"),
    dumpVertexes = cms.untracked.bool(True),
    dumpOnlyBremsstrahlung = cms.untracked.bool(False)  
)


process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("file:BPark_HLTEff_valid_withPU.root")
                                   #fileName = cms.string("file:BPark_HLTEff_valid_noPU.root")
                                   fileName = cms.string("file:testVtx.root")
                               )

process.p = cms.EndPath(process.myOutputTest)
