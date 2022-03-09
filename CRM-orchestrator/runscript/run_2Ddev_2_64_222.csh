#!/bin/csh
# This script automatically download cesm2 
# Once the job is finished, it helps to submit the job
set run_time       = 02:00:00
#set queue          = skx-normal
set queue          = development
set account        = ATM20009
set run_start_date = "0001-01-01"
set pcount         = 50
set CRM_pcount     = 12          
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
set crm_nx_in         = 64         # <<< change this one!
set crm_ny_in         = 1
set crm_dx_in         = 200
set crm_dt_in         = 0.5
set crm_nz_in         = 24
set spcam_subx_in     = 2
set spcam_suby_in     = 1
set spcam_orctotal_in = 384
@ CRM_pcount       = $spcam_orctotal_in * $spcam_subx_in * $spcam_suby_in
@ NPNN = $pcount +  $CRM_pcount
@ NPNNX = 50
@ NNODE = $NPNN / $NPNNX + 1
setenv CASE       scalli_check_${NPNNX}_crmnx${crm_nx_in}_crmny${crm_ny_in}_subx${spcam_subx_in}_suby${spcam_suby_in}_${spcam_orctotal_in}orc_${NNODE}nodes_${queue}

setenv CASEROOT  $SCRATCH/CESM2_case/$CASE

setenv JOB_QUEUE   $queue
setenv SCRIPTDIR   $HOME/repositories/$CCSMTAG/CRM-orchestrator/runscript
### GRID OPTIONS <Liran>
set crm_nx_in         = 64         # <<< change this one!
set crm_ny_in         = 1
set crm_dx_in         = 200
set crm_dt_in         = 0.5
set crm_nz_in         = 24
set spcam_subx_in     = 2
set spcam_suby_in     = 1
set spcam_orctotal_in = 384
@ CRM_pcount       = $spcam_orctotal_in * $spcam_subx_in * $spcam_suby_in
@ NPNN = $pcount +  $CRM_pcount
@ NPNNX = 50
@ NNODE = $NPNN / $NPNNX + 1
setenv CASE       scalli_check_${NPNNX}_crmnx${crm_nx_in}_crmny${crm_ny_in}_subx${spcam_subx_in}_suby${spcam_suby_in}_${spcam_orctotal_in}orc_${NNODE}nodes_${queue}
cd  $CASEROOT
set NODETEST     = 17
./case.submit
sed -e "s/NPN/$NPNNX/g; s/NNODEEE/$NODETEST/g; s/CASE_FOLDER/$CASE/g" $HOME/repositories/CESM-OCR/CRM-orchestrator/runscript/CRM/case.run.sample > $CASEROOT/.case.run
sbatch --time ${run_time} -p ${queue} --nodes=${NNODE} --account ATM20009 .case.run --resubmit
cd ..
#./case.submit

