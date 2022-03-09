#!/bin/bash

for NX in 64  #64 #128
do
for SUBX in 2 4 8  #4 8 #16 32 64  
do
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_v1.csh > run_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_2D.csh > run2D_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_dev.csh > run_${SUBX}_${NX}.csh
sed -e "s/SUBXX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_auto_dev_2D_nodebug.csh > run_2Ddev_nodebug_${SUBX}_${NX}.csh
chmod 700 run_2Ddev_nodebug_${SUBX}_${NX}.csh
./run_2Ddev_nodebug_${SUBX}_${NX}.csh
done
done
