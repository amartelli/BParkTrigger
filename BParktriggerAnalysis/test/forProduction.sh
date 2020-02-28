#### RAW with Run3 configuration from DAS relval  (release CMSSW_11_0_0)

# 0PU
#cmsDriver.py step1 --filein file:/eos/cms/store/user/fiorendi/parking/MC_forHLT2021/BPH-RunIIFall18GS-00214_forHLT_noTagMu_part0.root --fileout file:/tmp/amartell/BPH-RunIIFall18_DIGIRAW-00214_forHLT_noTagMu_part0_noPU.root --mc --eventcontent FEVTDEBUGHLT  --datatier GEN-SIM-DIGI-RAW --conditions auto:phase1_2021_realistic --step DIGI,L1,DIGI2RAW --nThreads 8 --geometry DB:Extended --era Run3 --python_filename step1_DIGI_RAW_noPU_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1

# withPU  (55 to 75)
#cmsDriver.py step1 --filein file:/eos/cms/store/user/fiorendi/parking/MC_forHLT2021/BPH-RunIIFall18GS-00214_forHLT_noTagMu_part0.root --fileout file:/tmp/amartell/BPH-RunIIFall18_DIGIRAW-00214_forHLT_noTagMu_part0_wPU.root --pileup_input das:/RelValMinBias_14TeV/CMSSW_11_0_0_pre13-110X_mcRun3_2021_realistic_v6-v1/GEN-SIM --mc --eventcontent FEVTDEBUGHLT  --pileup Run3_Flat55To75_PoissonOOTPU --datatier GEN-SIM-DIGI-RAW --conditions auto:phase1_2021_realistic --step DIGI,L1,DIGI2RAW --nThreads 8 --geometry DB:Extended --era Run3 --python_filename step1_DIGI_RAW_wPU_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1                                                                                                                             


#### HLT with pixelTracks  (managed to run on @patatrack02 => need access rights  & release CMSSW_11_1_0_pre2_Patatrack)

#cmsDriver.py step3 --conditions auto:phase1_2021_realistic -n 100 --era Run3 --eventcontent RECOSIM,DQM --runUnscheduled --customise_unsch=RecoPixelVertexing/Configuration/customizePixelTracksSoAonCPU.customizePixelTracksSoAonCPU -s RAW2DIGI:RawToDigi_pixelOnly,RECO:reconstruction_pixelTrackingOnly,VALIDATION:@pixelTrackingOnlyValidation,DQM:@pixelTrackingOnlyDQM --datatier GEN-SIM-RECO,DQMIO --geometry DB:Extended --fileout file:step3.root --dirin=/eos/cms/store/user/amartell/BPrk/MC_forHLT2021/ --filein BPH-RunIIFall18_DIGIRAW-00214_forHLT_noTagMu_part0_noPU.root --no_exec --python_filename=step3_RAW-RECOwHLT_Valid.py

# relax cuts for tracking particles (add missing lines)
##from Configuration.Eras.Era_Run3_cff import Run3
##process = cms.Process('RECO',Run3)
#from SimGeneral.MixingModule.trackingTruthProducer_cfi import trackingParticles
#trackingParticles.vertexDistanceCut = cms.double(100.)
## import of standard configurations

# keep trackingParticles and association map (trackingP - pixelTracks) to output
## Additional output definition
#process.RECOSIMoutput.outputCommands += cms.untracked.vstring('keep *_*_MergedTrackTruth_*', 
#                                                              'keep *_trackingParticlePixelTrackAsssociation_*_*')


##input GEN-SIM in    /eos/cms/store/user/fiorendi/parking/MC_forHLT2021/
##input RAW in        /eos/cms/store/user/amartell/BPrk/MC_forHLT2021/
##with HLT same dir   /eos/cms/store/user/amartell/BPrk/MC_forHLT2021/
