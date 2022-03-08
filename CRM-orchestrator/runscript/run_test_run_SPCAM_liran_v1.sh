#!/bin/bash

for NX in 32 64 128
do
for SUBX in  6 12
do
sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_v1.csh > run_${SUBX}_${NX}.csh
chmod 700 run_${SUBX}_${NX}.csh
./run_${SUBX}_${NX}.csh
done
done
