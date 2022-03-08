#!/bin/bash

for NX in 128
do
for SUBX in 2 4 8 16 #4 8 #16 32 64  
do
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_v1.csh > run_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_2D.csh > run2D_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_normal.csh > run_${SUBX}_${NX}.csh
sed -e "s/SUBXX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_auto_normal_2D.csh > run_2Dnormal_${SUBX}_${NX}.csh
chmod 700 run_2Dnormal_${SUBX}_${NX}.csh
./run_2Dnormal_${SUBX}_${NX}.csh
done
done
