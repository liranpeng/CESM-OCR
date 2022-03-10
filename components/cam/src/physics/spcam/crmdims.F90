module crmdims
#define ORCHESTRATOR
#ifdef CRM
    use shr_kind_mod,    only: r8 => shr_kind_r8
    implicit none

       integer, parameter :: nclubbvars = 17

       integer, parameter ::  crm_nx=SPCAM_NX, crm_ny=SPCAM_NY, crm_nz=SPCAM_NZ
       real(r8), parameter :: crm_dx=SPCAM_DX, crm_dy=SPCAM_DX, crm_dt=SPCAM_DT


#ifdef ORCHESTRATOR
       integer, parameter ::  orc_total = 30
       integer, parameter ::  SPCAM_ORC_NSUBDOMAINS_X = 2
       integer, parameter ::  SPCAM_ORC_NSUBDOMAINS_Y = 1
       integer, parameter ::  orc_nsubdomains_x=SPCAM_ORC_NSUBDOMAINS_X
       integer, parameter ::  orc_nsubdomains_y=SPCAM_ORC_NSUBDOMAINS_Y
       integer, parameter ::  orc_nsubdomains = orc_nsubdomains_x*orc_nsubdomains_y
       integer, parameter ::  orc_rank_total = orc_total*orc_nsubdomains
       integer, parameter ::  orc_nx = int(crm_nx/orc_nsubdomains_x)
       integer, parameter ::  orc_ny = int(crm_ny/orc_nsubdomains_y)
#endif

#endif

end module crmdims
