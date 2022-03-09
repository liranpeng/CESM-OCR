#!/bin/bash

for NX in 32  #64 #128
do
for SUBXX in 2 4 8 #4 8 #16 32 64  
do
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_v1.csh > run_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_3D.csh > run3D_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_dev.csh > run_${SUBX}_${NX}.csh
sed -e "s/SUBXXXX/${SUBXX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_auto_dev_3D.csh > run_3Ddev_${SUBXX}_${NX}.csh
chmod 700 run_3Ddev_${SUBXX}_${NX}.csh
./run_3Ddev_${SUBXX}_${NX}.csh
done
done