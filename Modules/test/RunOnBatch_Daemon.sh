#!/bin/bash
cd /afs/cern.ch/work/l/lpernie/pi0/HLT_pi0/CMSSW_5_3_3_patch2/src/Analysis/Modules/test
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`
echo 'python /afs/cern.ch/work/l/lpernie/pi0/HLT_pi0/CMSSW_5_3_3_patch2/src/Analysis/Modules/test/RunOnBatch.py'
python /afs/cern.ch/work/l/lpernie/pi0/HLT_pi0/CMSSW_5_3_3_patch2/src/Analysis/Modules/test/RunOnBatch.py
