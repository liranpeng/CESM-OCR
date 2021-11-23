#!/bin/csh
# This script automatically download cesm2 
# Once the job is finished, it helps to submit the job

set run_time       = 01:00:00
#set queue          = skx-normal
set queue          = development
set account        = ATM20009
set run_start_date = "0001-01-01"
set pcount         = 50

## ====================================================================
#   define case
## ====================================================================

setenv CCSMTAG     CESM-2021-01
setenv CASE        SPdev_liran_1mom_v124
#bloss(2021-01-22): Revert to basic case for testing
#setenv CASESET     HIST_CAM%SPCAMS_CLM50%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
#setenv CASERES     f19_f19_mg17
#CESM standard case
#setenv CASESET     B1850
#setenv CASERES     f19_g17
# Two moment microphysics SPCAM for testing
#   and movement toward aerosol-cloud interactions
setenv CASESET     FSPCAMS
setenv CASERES     f45_g37
setenv PROJECT     ATM20009
setenv JOB_QUEUE   $queue
setenv SCRIPTDIR   $HOME/CleanVersion/$CCSMTAG/CRM-orchestrator/runscript
## ====================================================================
#   define directories <Please make sure the directories are correct>
## ====================================================================

#setenv MACH      stampede2-skx
#setenv MACH      stampede2-knl
setenv MACH      frontera
setenv CCSMROOT  $HOME/CleanVersion/$CCSMTAG
setenv CASEROOT  $SCRATCH/CESM2_case/$CASE
setenv PTMP      $SCRATCH/
setenv RUNDIR    $PTMP/$CASE/run
setenv ARCHDIR   $PTMP/archive/$CASE
setenv DATADIR   /scratch1/07088/tg863871/inputdata # pritch, link to bloss' downloaded input data.
# Note: my input folder is /scratch1/07088/tg863871/inputdata
## ====================================================================
#   Download model source code <This part only need to do once>
## ====================================================================
# A reference for CESM run https://escomp.github.io/CESM/release-cesm2/downloading_cesm.html
# Run svn ls https://svn-ccsm-models.cgd.ucar.edu/ww3/release_tags, permanently accepting the certificate 
# when prompted, then retry the CESM download (starting over at the top of these instructions).
# ===================================================================== 	
#git clone -b release-cesm2.2.0 https://github.com/ESCOMP/CESM.git $CCSMROOT
#cd $CCSMROOT
#./manage_externals/checkout_externals
#./manage_externals/checkout_externals --logging
#svn ls https://svn-ccsm-models.cgd.ucar.edu/ww3/release_tags
## ====================================================================
#   create new case, configure, compile and run
## ====================================================================

rm -rf $CASEROOT
rm -rf $PTMP/$CASE
cp /home1/07088/tg863871/CleanVersion/CESM-2021-01/CRM-orchestrator/runscript/CRM/TwoExecutableDriver.F90 /home1/07088/tg863871/CleanVersion/CESM-2021-01/components/cam/src/physics/spcam/crm/
#------------------
## create new case
#------------------

cd $CCSMROOT/cime/scripts

#bloss(2021-01-22): Revert to basic pelayout for testing
#./create_newcase --case $CASEROOT --pecount $pcount --pesfile ./pelayout_frontera01.xml --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT --run-unsupported
./create_newcase --case $CASEROOT  --pecount $pcount --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT  --run-unsupported

cd  $CASEROOT

#cat <<EOF >> user_nl_cam
#nhtfrq =   1
#mfilt  = 1
#avgflag_pertape = 'A'
#fincl1 = 'CRM_U','CRM_W','CRM_T'
#EOF

xmlchange --file env_batch.xml --id JOB_QUEUE --val $queue
xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val $run_time
#xmlchange --file env_run.xml --id RUN_STARTDATE --val $run_start_date
xmlchange --file env_run.xml --id STOP_OPTION --val nhour
xmlchange --file env_run.xml --id STOP_N --val 10
xmlchange --file env_run.xml --id ATM_NCPL --val 432  


#xmlchange --file env_run.xml --id run_data_archive --val "FALSE"
#xmlchange --file env_run.xml --id RESUBMIT --val 4

./case.setup

cp $SCRIPTDIR/src.cam/*F90 SourceMods/src.cam/.

./case.build
pwd
cp -R $SCRIPTDIR/CRM/ CRM/
cd CRM

#mpif90 -I../$CASE/bld/lib/include -I../$CASE/bld/intel/impi/nodebug/nothreads/mct/include -I../$CASE/bld/intel/impi/nodebug/nothreads/mct/bld/intel/impi/nodebug/nothreads/mct/mct/mct -I../$CASE/bld/atm/obj -O2 -debug minimal -o ../$CASE/bld/crm.exe TwoExecutableDriver.F90

cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/TwoExecutableDriver.o . 
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_crm_module_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmdims.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/ppgrid.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/phys_grid.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/mct/noesmf/c1a1l1i1o1r1g1w1i1e1/csm_share/shr_kind_mod.o .

#mpif90 -I../$CASE/bld/lib/include -I../$CASE/bld/atm/obj shr_kind_mod.o spmd_utils.o crmx_crm_module.o ppgrid.o crmdims.o TwoExecutableDriver.o -O2 -debug minimal -o ../$CASE/bld/crm.exe
#mpif90 -I../$CASE/bld/lib/include -I../$CASE/bld/atm/obj -O2 -debug minimal -o ../$CASE/bld/crm.exe TwoExecutableDriver.F90

mpif90  -o /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/crm.exe crmdims.o ppgrid.o crmx_crm_module_ORC.o  shr_kind_mod.o phys_grid.o TwoExecutableDriver.o -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -latm -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lice  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/mct/noesmf/lib/ -lclm  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -locn  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lrof  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lglc  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lwav  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lesp  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -liac -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/mct/noesmf/c1a1l1i1o1r1g1w1i1e1/lib -lcsm_share -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/lib -lpiof -lpioc -lgptl -lmct -lmpeu  -lpthread   -mkl=cluster   -L/opt/apps/intel19/impi19_0/parallel-netcdf/4.7.4/x86_64 -lnetcdff -Wl,--as-needed,-L/opt/apps/intel19/impi19_0/parallel-netcdf/4.7.4/x86_64/lib -lnetcdff -lnetcdf  -L/opt/apps/intel19/impi19_0/pnetcdf/1.11.2/lib -lpnetcdf



cd ..

#./case.submit

