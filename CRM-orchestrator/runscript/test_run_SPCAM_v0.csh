#!/bin/csh
# This script automatically download cesm2 
# Once the job is finished, it helps to submit the job

set run_time       = 01:00:00
#set queue          = skx-normal
set queue          = development
set account        = ATM20009
set run_start_date = "2003-01-02"
set pcount         = 50

## ====================================================================
#   define case
## ====================================================================

setenv CCSMTAG     CESM-2021-01
setenv CASE        SPdevel_al
#bloss(2021-01-22): Revert to basic case for testing
#setenv CASESET     HIST_CAM%SPCAMS_CLM50%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
#setenv CASERES     f19_f19_mg17
#CESM standard case
#setenv CASESET     B1850
#setenv CASERES     f19_g17
# Two moment microphysics SPCAM for testing
#   and movement toward aerosol-cloud interactions
setenv CASESET     FSPCAMM
setenv CASERES     f45_g37
setenv PROJECT     ATM20009
setenv JOB_QUEUE   $queue
setenv SCRIPTDIR   $PWD
## ====================================================================
#   define directories <Please make sure the directories are correct>
## ====================================================================

#setenv MACH      stampede2-skx
#setenv MACH      stampede2-knl
setenv MACH      frontera
setenv CCSMROOT  $HOME/$CCSMTAG
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

#------------------
## create new case
#------------------

cd $CCSMROOT/cime/scripts

#bloss(2021-01-22): Revert to basic pelayout for testing
#./create_newcase --case $CASEROOT --pecount $pcount --pesfile ./pelayout_frontera01.xml --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT --run-unsupported
./create_newcase --case $CASEROOT  --pecount $pcount --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT  --run-unsupported

cd  $CASEROOT

xmlchange --file env_batch.xml --id JOB_QUEUE --val $queue
xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val $run_time
#xmlchange --file env_run.xml --id RUN_STARTDATE --val $run_start_date
xmlchange --file env_run.xml --id STOP_OPTION --val nhour
xmlchange --file env_run.xml --id STOP_N --val 1
#xmlchange --file env_run.xml --id run_data_archive --val "FALSE"
#xmlchange --file env_run.xml --id RESUBMIT --val 4

./case.setup

cp $SCRIPTDIR/src.cam/*F90 SourceMods/src.cam/.

./case.build

cp -R $SCRIPTDIR/CRM/ CRM/
cd CRM
mpif90 -I../$CASE/bld/lib/include -I../$CASE/bld/atm/obj -o ../$CASE/bld/crm.exe TwoExecutableDriver.f90
cd ..

#./case.submit

