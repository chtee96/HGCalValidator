#!/bin/sh

# This is the shell script that will generate all the skeletons using cmsDriver commands.
# The commands included have been taken from runTheMatrix with the following command:
#
# runTheMatrix.py -w upgrade -l 20488.0 --command="--no_exec" --dryRun
#
# For all commands remove --filein and --fileout options.
# Add python_filename option
#
# The first command combines step1 and step2 (GSD):
# - run up to DIGI...HLT:@fake2
# The following changes are implemented on top:
# --beamspot HLLHC ➜ --beamspot NoSmear
# --eventcontent FEVTDEBUG ➜ --eventcontent FEVTDEBUGHLT
# Removed --relval option.
#
# The second command is step3 removing overlap with step2 (RECO):
# - remove MINIAODSIM from event content and data tier
# - remove PAT from steps (-s)
# - remove @miniAODValidation from VALIDATION step
# - remove @miniAODDQM from DQM step
#
# The third command is a copy of the second only re-running RECO (for NTUP):
# - remove DQM from event content
# - remove DQMIO from data tier
# - add --processName=NTUP option
#
# Those commands should be regularly checked and, in case of changes, propagated into this script!


  cmsDriver.py CE_E_Front_200um_cfi \
    --conditions auto:phase2_realistic_T21 \
    -n 10 \
    --era Phase2C11I13M9 \
    --eventcontent FEVTDEBUGHLT \
    -s GEN,SIM,DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,HLT:@fake2 \
    --datatier GEN-SIM \
    --beamspot NoSmear \
    --geometry Extended2026D88 \
    --no_exec \
    --python_filename=GSD_fragment.py


  cmsDriver.py step3 \
    --conditions auto:phase2_realistic_T21 \
    -n 10 \
    --era Phase2C11I13M9 \
    --eventcontent FEVTDEBUGHLT,DQM \
    -s RAW2DIGI,L1Reco,RECO,RECOSIM,VALIDATION:@phase2Validation,DQM:@phase2 \
    --datatier GEN-SIM-RECO,DQMIO \
    --geometry Extended2026D88 \
    --no_exec \
    --python_filename=RECO_fragment.py

    # define what to inject where
    hook="# customisation of the process."
    content="process.TFileService = cms.Service(\n    \"TFileService\",\n    fileName = cms.string(\"hgc.root\"),\n    closeFileFast = cms.untracked.bool(True)\n)"

    # do the injection
    sed "/$hook/a $content" -i RECO_fragment.py



  cmsDriver.py step3 \
    --conditions auto:phase2_realistic_T21 \
    -n 10 \
    --era Phase2C11I13M9 \
    --eventcontent FEVTDEBUGHLT \
    -s RAW2DIGI,L1Reco,RECO,RECOSIM \
    --datatier GEN-SIM-RECO \
    --geometry Extended2026D86 \
    --no_exec \
    --processName=NTUP \
    --python_filename=NTUP_fragment.py




