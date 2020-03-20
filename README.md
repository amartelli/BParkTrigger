# BParkTrigger

## production
-   RAW with Run3 configuration from DAS relval  (release CMSSW_11_0_0), 0PU and 55-75PU

-   HLT with pixelTracks: on @patatrack02 => need access rights  & release CMSSW_11_1_0_pre2_Patatrack

-   cmsDrivers in BParkTrigger/BParktriggerAnalysis/test/forProduction.sh


## analysis
```shell
cmsrel CMSSW_11_1_0_pre2_Patatrack
cd CMSSW_11_1_0_pre2_Patatrack/src
cmsenv
git clone git@github.com:amartelli/BParkTrigger.git
```

- efficiency of reco wrt gen:
   - *pixelTracks wrt trackingParticles*	
   - *and pixelTracks with common pixelVertex*

- selections at reco-level MISSING => next steps
- estimate of efficiency vs rates MISSING => next steps
