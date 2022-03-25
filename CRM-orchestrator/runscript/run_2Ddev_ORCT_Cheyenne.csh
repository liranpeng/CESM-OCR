#!/bin/csh
# This script automatically download cesm2 
# Once the job is finished, it helps to submit the job
set run_time       = 01:00:00
#set queue          = skx-normal
#set queue          = economy
set queue          = regular
set account        = UWAS0096
set run_start_date = "0001-01-01"
set pcount         = 36
set taskPernode    = 36
set emailaddress   = liranp@uci.edu
## ====================================================================
#   define case
## ====================================================================
setenv CCSMTAG     CESM-OCR
#setenv CASE        SPdev_liran_1mom_v864
#bloss(2021-01-22): Revert to basic case for testing
#setenv CASESET     HIST_CAM%SPCAMS_CLM50%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
#setenv CASERES     f19_f19_mg17
#CESM standard case
#setenv CASESET     B1850
#setenv CASERES     f19_g17
# Two moment microphysics SPCAM for testing
#   and movement toward aerosol-cloud interactions
setenv CASESET     FSPCAMS
setenv CASERES     f10_f10_mg37
setenv PROJECT     UWAS0096
setenv UNAME       lpeng
setenv JOB_QUEUE   $queue
setenv SCRATCH     /glade/scratch/lpeng 
setenv SCRIPTDIR   $HOME/repositories/$CCSMTAG/CRM-orchestrator/runscript
### GRID OPTIONS <Liran>
set crm_nx_in         = 1024         # <<< change this one!
set crm_ny_in         = 1
set crm_dx_in         = 200
set crm_dt_in         = 0.4
set crm_nz_in         = 24
set spcam_subx_in     = 4
set spcam_suby_in     = 1
set spcam_orctotal_in = 456
@ CRM_pcount       = $spcam_orctotal_in * $spcam_subx_in * $spcam_suby_in
@ NPNN = $pcount +  $CRM_pcount
@ NNODE = $NPNN / $taskPernode + 1
setenv CASE       scalling_timming_GCMRes_${CASERES}_GCMTask${pcount}_crmnx${crm_nx_in}_crmdt04_crmny${crm_ny_in}_subx${spcam_subx_in}_suby${spcam_suby_in}_${spcam_orctotal_in}orc_${NNODE}nodes_${queue}
## ====================================================================
#   define directories <Please make sure the directories are correct>
## ====================================================================
#setenv MACH      stampede2-skx
#setenv MACH      stampede2-knl
setenv MACH      cheyenne
setenv CCSMROOT  $HOME/repositories/$CCSMTAG
setenv CASEROOT  $SCRATCH/CESM2_case/$CASE
setenv PTMP      $SCRATCH/
setenv RUNDIR    $PTMP/$CASE/run
setenv ARCHDIR   $PTMP/archive/$CASE
# Note: my input folder is $SCRATCH/inputdata
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
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/TwoExecutableDriver.F90 $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crm/
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/src.cam/crmx_crm_module_ORC.F90 $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crm/
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/src.cam/crmx_crm_module.F90 $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crm/
#------------------
## create new case
#------------------
cd $CCSMROOT/cime/scripts
#bloss(2021-01-22): Revert to basic pelayout for testing
#./create_newcase --case $CASEROOT --pecount $pcount --pesfile ./pelayout_frontera01.xml --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT --run-unsupported
./create_newcase --case $CASEROOT  --pecount $pcount --res $CASERES --machine $MACH --compset $CASESET  --output-root $CASEROOT  --run-unsupported
cd  $CASEROOT
sed -e "s/SUBXdim/$spcam_subx_in/g; s/SUBYdim/$spcam_suby_in/g; s/ORCT/$spcam_orctotal_in/g;" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/crmdims.sample > $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crmdims.F90 
sed -e "s/GCM_pcount/$pcount/g; s/CRM_pcount/$CRM_pcount/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/env_mach_specific.sample  > $CASEROOT/env_mach_specific.xml
sed -e "s/NXX/$crm_nx_in/g; s/NYY/$crm_ny_in/g; s/DXX/$crm_dx_in/g; s/DTT/$crm_dt_in/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/config_component.sample > $HOME/repositories/CESM-OCR/components/cam/cime_config/config_component.xml
sed -e "s/NPN/$taskPernode/g; s/USERNAME/$UNAME/g; s/NNNODE/$NNODE/g; s/CCCASE/$CASE/g; s/PPRO/$PROJECT/g; s/QQUE/$queue/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.sample > $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.send
#cat <<EOF >> user_nl_drv
#atm_cpl_dt = 10
#EOF
cat <<EOF >> user_nl_cam
fincl2 = 'PS:I'
nhtfrq = 0,1
mfilt = 0,1
EOF
mkdir $SCRATCH/CESM2_case/$CASE/archive
cd  $CASEROOT
./xmlchange --file env_batch.xml --id JOB_QUEUE --val $queue
./xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val $run_time
./xmlchange --file env_run.xml --id STOP_OPTION --val ndays
./xmlchange --file env_run.xml --id STOP_N --val 2
./xmlchange --file env_run.xml --id DOUT_S_ROOT --val $SCRATCH/CESM2_case/$CASE/archive
./xmlchange --file env_run.xml --id DOUT_S --val TRUE
./case.setup
./xmlchange --file env_run.xml --id run_data_archive --val "FALSE"
#xmlchange --file env_run.xml --id RESUBMIT --val 4
cp $SCRIPTDIR/src.cam/*F90 SourceMods/src.cam/.
qcmd -A ${PROJECT} -- ./case.build --skip-provenance-check
pwd
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.send .case.run
cp -R $SCRIPTDIR/CRM/ CRM/
cd CRM
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_task_init_ORC.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/TwoExecutableDriver.o . 
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_crm_module_ORC.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmdims.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/ppgrid.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/phys_grid.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_mpi.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/task_assign_bnd.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/task_dispatch.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/task_exchange.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_task_util_ORC.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_kurant_ORC.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_pressure_ORC.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_press_rhs_ORC.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/intel/mpt/nodebug/nothreads/mct/mct/noesmf/c1a1l1i1o1r1g1w1i1e1/csm_share/shr_kind_mod.o .
cp $SCRATCH/CESM2_case/$CASE/$CASE/bld/intel/mpt/nodebug/nothreads/mct/gptl/perf_mod.o .
mpif90  -o $SCRATCH/CESM2_case/$CASE/$CASE/bld/crm.exe perf_mod.o task_exchange.o crmx_pressure_ORC.o crmx_press_rhs_ORC.o crmx_kurant_ORC.o task_dispatch.o task_assign_bnd.o crmx_mpi.o crmx_task_util_ORC.o crmx_task_init_ORC.o crmdims.o ppgrid.o crmx_crm_module_ORC.o  shr_kind_mod.o phys_grid.o TwoExecutableDriver.o -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -latm -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -lice  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/intel/mpt/nodebug/nothreads/mct/mct/noesmf/lib/ -lclm  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -locn  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -lrof  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -lglc  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -lwav  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -lesp  -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/lib/ -liac -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/intel/mpt/nodebug/nothreads/mct/mct/noesmf/c1a1l1i1o1r1g1w1i1e1/lib -lcsm_share -L$SCRATCH/CESM2_case/$CASE/$CASE/bld/intel/mpt/nodebug/nothreads/mct/lib -lpiof -lpioc -lgptl -lmct -lmpeu   -mkl=cluster  -L/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.21/intel/19.0.5//lib -lpnetcdf -L/glade/u/apps/ch/opt/netcdf/4.7.3/intel/19.0.5//lib -lnetcdff -lnetcdf
cd  $CASEROOT
./case.submit
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.send .case.run

echo "=================================================================================================="
echo "=================================================================================================="
echo "The model is currently run using interactive jobs on Cheyenne"
echo "Case Directory is"
echo "cd" $CASEROOT
echo "qsub -I -l select="${NNODE}":ncpus=36:mpiprocs=36 -l walltime=01:00:00 -q regular -A "${account}
echo "./.case.run"
echo "=================================================================================================="
echo "=================================================================================================="

#qsub -l select=${NNODE}:ncpus=${taskPernode}:mpiprocs=${taskPernode}:ompthreads=2 -l -walltime ${run_time} -q ${queue}  -A ${account} .case.run -resubmit
#cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.send .case.run
#cd ..
#./case.submit

