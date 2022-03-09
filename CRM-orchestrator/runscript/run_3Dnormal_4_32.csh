#!/bin/csh
# This script automatically download cesm2 
# Once the job is finished, it helps to submit the job
set run_time       = 02:00:00
#set queue          = skx-normal
set queue          = normal
set account        = ATM20009
set run_start_date = "0001-01-01"
set pcount         = 50
set CRM_pcount     = 12          
@ NPNN = $pcount +  $CRM_pcount
set NNODE = 2
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
setenv CASERES     f45_g37
setenv PROJECT     ATM20009
setenv JOB_QUEUE   $queue
setenv SCRIPTDIR   $HOME/repositories/$CCSMTAG/CRM-orchestrator/runscript
### GRID OPTIONS <Liran>
set crm_nx_in         = 32         # <<< change this one!
set crm_ny_in         = 32
set crm_dx_in         = 4000
set crm_dt_in         = 1
set crm_nz_in         = 24
set spcam_subx_in     = 4
set spcam_suby_in     = 4
set spcam_orctotal_in = 6
@ CRM_pcount       = $spcam_orctotal_in * $spcam_subx_in * $spcam_suby_in
@ NPNN = $pcount +  $CRM_pcount
@ NNODE = $NPNN / 50 + 1
setenv CASE       scalling2_crmnx${crm_nx_in}_crmny${crm_ny_in}_subx${spcam_subx_in}_suby${spcam_suby_in}_${spcam_orctotal_in}orc_${NNODE}nodes_${queue}
## ====================================================================
#   define directories <Please make sure the directories are correct>
## ====================================================================
#setenv MACH      stampede2-skx
#setenv MACH      stampede2-knl
setenv MACH      frontera
setenv CCSMROOT  $HOME/repositories/$CCSMTAG
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
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/TwoExecutableDriver.F90 $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crm/
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/src.cam/crmx_crm_module_ORC.F90 $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crm/
cp $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/src.cam/crmx_crm_module.F90 $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crm/
#------------------
## create new case
#------------------
cd $CCSMROOT/cime/scripts
#bloss(2021-01-22): Revert to basic pelayout for testing
#./create_newcase --case $CASEROOT --pecount $pcount --pesfile ./pelayout_frontera01.xml --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT --run-unsupported
./create_newcase --case $CASEROOT  --pecount $pcount --res $CASERES --machine $MACH --compset $CASESET --input-dir $DATADIR --output-root $CASEROOT  --run-unsupported
cd  $CASEROOT
sed -e "s/SUBXdim/$spcam_subx_in/g; s/SUBYdim/$spcam_suby_in/g; s/ORCT/$spcam_orctotal_in/g;" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/crmdims.sample > $HOME/repositories/CESM-OCR/components/cam/src/physics/spcam/crmdims.F90 
sed -e "s/NPN/$NPNN/g; s/NNODE/$NNODE/g; s/CASE_FOLDER/$CASE/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.sample > $CASEROOT/.case.run
sed -e "s/GCM_pcount/$pcount/g; s/CRM_pcount/$CRM_pcount/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/env_mach_specific.sample  > $CASEROOT/env_mach_specific.xml
sed -e "s/NXX/$crm_nx_in/g; s/NYY/$crm_ny_in/g; s/DXX/$crm_dx_in/g; s/DTT/$crm_dt_in/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/config_component.sample > $HOME/repositories/CESM-OCR/components/cam/cime_config/config_component.xml
#cat <<EOF >> user_nl_drv
#atm_cpl_dt = 10
#EOF
xmlchange --file env_batch.xml --id JOB_QUEUE --val $queue
xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val $run_time
#xmlchange --file env_run.xml --id RUN_STARTDATE --val $run_start_date
xmlchange --file env_run.xml --id STOP_OPTION --val nhour
xmlchange --file env_run.xml --id STOP_N --val 2
xmlchange --file env_run.xml --id ATM_NCPL --val 432
./case.setup
#xmlchange --file env_run.xml --id run_data_archive --val "FALSE"
#xmlchange --file env_run.xml --id RESUBMIT --val 4
cp $SCRIPTDIR/src.cam/*F90 SourceMods/src.cam/.
./case.build
pwd
cp -R $SCRIPTDIR/CRM/ CRM/
cd CRM
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_task_init_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/TwoExecutableDriver.o . 
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_crm_module_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmdims.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/ppgrid.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/phys_grid.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_mpi.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/task_assign_bnd.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/task_dispatch.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/task_exchange.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_task_util_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_kurant_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_pressure_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/atm/obj/crmx_press_rhs_ORC.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/mct/noesmf/c1a1l1i1o1r1g1w1i1e1/csm_share/shr_kind_mod.o .
cp /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/gptl/perf_mod.o .
mpif90  -o /scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/crm.exe perf_mod.o task_exchange.o crmx_pressure_ORC.o crmx_press_rhs_ORC.o crmx_kurant_ORC.o task_dispatch.o task_assign_bnd.o crmx_mpi.o crmx_task_util_ORC.o crmx_task_init_ORC.o crmdims.o ppgrid.o crmx_crm_module_ORC.o  shr_kind_mod.o phys_grid.o TwoExecutableDriver.o -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -latm -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lice  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/mct/noesmf/lib/ -lclm  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -locn  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lrof  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lglc  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lwav  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -lesp  -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/lib/ -liac -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/mct/noesmf/c1a1l1i1o1r1g1w1i1e1/lib -lcsm_share -L/scratch1/07088/tg863871/CESM2_case/$CASE/$CASE/bld/intel/impi/nodebug/nothreads/mct/lib -lpiof -lpioc -lgptl -lmct -lmpeu  -lpthread   -mkl=cluster   -L/opt/apps/intel19/impi19_0/parallel-netcdf/4.7.4/x86_64 -lnetcdff -Wl,--as-needed,-L/opt/apps/intel19/impi19_0/parallel-netcdf/4.7.4/x86_64/lib -lnetcdff -lnetcdf  -L/opt/apps/intel19/impi19_0/pnetcdf/1.11.2/lib -lpnetcdf
cd  $CASEROOT
./case.submit
sbatch --time ${run_time} -p ${queue} --nodes=${NNODE} --account ATM20009 .case.run --resubmit
cd ..
#./case.submit

