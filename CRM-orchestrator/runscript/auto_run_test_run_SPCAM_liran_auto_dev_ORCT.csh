#!/bin/bash

for NX in 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380  #64 #128
do
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_v1.csh > run_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_2D.csh > run2D_${SUBX}_${NX}.csh
#sed -e "s/SBX/${SUBX}/g; s/NNNX/${NX}/g" test_run_SPCAM_liran_dev.csh > run_${SUBX}_${NX}.csh
sed -e "s/ORCTTT/${NX}/g" run_2Ddev_2_64_sample.csh > run_2Ddev_ORCT_${NX}.csh
chmod 700 run_2Ddev_ORCT_${NX}.csh
#./run_2Ddev_ORCT_${NX}.csh
done
