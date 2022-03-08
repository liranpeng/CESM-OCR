#define ORCHESTRATOR
module crm_physics
!-----------------------------------------------------------------------
! Purpose: 
! 
!    Provides the CAM interface to the crm code.  
!
! Revision history: 
! June, 2009, Minghuai Wang:  
!          crm_physics_tend 
! July, 2009, Minghuai Wang: m2005_effradius
!
!---------------------------------------------------------------------------
   use crmx_grid
   use crmx_mpi
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, pver, pverp
#ifdef CRM
   use cam_abortutils,  only: endrun
   use physics_types,   only: physics_state, physics_tend
   use constituents,    only: cnst_add, cnst_get_ind, cnst_set_spec_class, cnst_spec_class_cldphysics, &
                              cnst_spec_class_gas, cnst_name, cnst_longname, sflxnam, apcnst, bpcnst, pcnst
#ifdef m2005
   use module_ecpp_ppdriver2,   only: papampollu_init
   use crmx_ecppvars,   only: NCLASS_CL,ncls_ecpp_in,NCLASS_PR
#endif
#ifdef ORCHESTRATOR
   use spmd_utils, only: npes
   use cam_logfile,      only: iulog ! pritch
   use mpi
#endif                    
   implicit none 
#ifdef ORCHESTRATOR
   !include 'mpif.h' 
   !use mpi
#endif
   public
   save

   character(len=2)  :: spcam_direction='NS'     ! SPCAM 2D orientation

   public :: crm_physics_tend, crm_physics_register, crm_physics_init
   public :: crm_implements_cnst, crm_init_cnst
   public :: m2005_effradius

   integer :: crm_u_idx, crm_v_idx, crm_w_idx, crm_t_idx
   integer :: crm_qt_idx, crm_nc_idx, crm_qr_idx, crm_nr_idx, crm_qi_idx, crm_ni_idx 
   integer :: crm_qs_idx, crm_ns_idx, crm_qg_idx, crm_ng_idx, crm_qc_idx, crm_qp_idx, crm_qn_idx
   integer :: crm_t_rad_idx, crm_qv_rad_idx, crm_qc_rad_idx, crm_qi_rad_idx, crm_cld_rad_idx
   integer :: crm_nc_rad_idx, crm_ni_rad_idx, crm_qs_rad_idx, crm_ns_rad_idx, crm_qrad_idx
   integer :: crm_qaerwat_idx, crm_dgnumwet_idx
   integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx, snow_sed_idx, prec_pcw_idx, snow_pcw_idx
   integer :: cldo_idx, cld_idx, cldtop_idx
   integer :: rei_idx, rel_idx, rprdtot_idx, nevapr_idx, prain_idx
   integer :: wsedl_idx, dei_idx, des_idx, mu_idx, lambdac_idx
   integer :: rate1_cw2pr_st_idx
   integer :: qme_idx, icwmrdp_idx, rprddp_idx, icwmrsh_idx, rprdsh_idx
   integer :: nevapr_shcu_idx, nevapr_dpcu_idx, ast_idx
   integer :: fice_idx,acldy_cen_idx, cmfmc_sh_idx
   integer :: clubb_buffer_idx, tk_crm_idx, tke_idx, kvm_idx, kvh_idx, pblh_idx, tpert_idx
   integer :: sh_frac_idx, dp_frac_idx

   integer :: &
              ixcldliq,      &! cloud liquid amount index
              ixcldice,      &! cloud ice amount index
              ixnumliq,      &! cloud liquid number index
              ixnumice        ! cloud ice water index

   integer            :: nmodes

   integer, parameter :: ncnst = 4       ! Number of constituents
   integer            :: ncnst_use
   character(len=8), parameter :: &      ! Constituent names
                                 cnst_names(ncnst) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)

   logical           :: use_spcam, prog_modal_aero, do_clubb_sgs
   logical           :: is_spcam_m2005, is_spcam_sam1mom

   integer          :: crm_nx_ny

#endif

!========================================================================================================
contains
!========================================================================================================

!---------------------------------------------------------------------------------------------------------
subroutine crm_physics_register()
#ifdef CRM
!-------------------------------------------------------------------------------------------------------
! 
! Purpose:  add necessary fileds into physics buffer
!
!--------------------------------------------------------------------------------------------------------
  use spmd_utils,      only: masterproc
  use physconst,       only: mwdry, cpair
  use physics_buffer,  only: dyn_time_lvls, pbuf_add_field, dtype_r8
  use phys_control,    only: phys_getopts, cam_physpkg_is
  use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy, crm_dt, nclubbvars,orc_nsubdomains_x,orc_nsubdomains_y,orc_nsubdomains,orc_total
  use cam_history_support,only: add_hist_coord
  use crmx_setparm_mod,     only: setparm
  use rad_constituents, only: rad_cnst_get_info

  is_spcam_m2005   = cam_physpkg_is('spcam_m2005')
  is_spcam_sam1mom = cam_physpkg_is('spcam_sam1mom')

  call phys_getopts( use_spcam_out           = use_spcam)
  call phys_getopts( prog_modal_aero_out     = prog_modal_aero)
  call phys_getopts( do_clubb_sgs_out        = do_clubb_sgs)

  call rad_cnst_get_info(0, nmodes=nmodes)
  
  ! Register microphysics constituents and save indices.

  ncnst_use = 2
  call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
  call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
        longname='Grid box averaged cloud ice amount', is_convtran1=.true.)
  if (is_spcam_m2005) then
     call cnst_add(cnst_names(3), mwdry, cpair, 0._r8, ixnumliq, &
          longname='Grid box averaged cloud liquid number', is_convtran1=.false.)
     call cnst_add(cnst_names(4), mwdry, cpair, 0._r8, ixnumice, &
          longname='Grid box averaged cloud ice number', is_convtran1=.false.)
     ncnst_use = 4
  end if

  if(masterproc) then
      print*,'_________________________________________'
      print*,'_ Super-parameterization run ____________'
      print*,'crm_nx=',crm_nx,'   crm_ny=',crm_ny,'   crm_nz=',crm_nz
      print*,'crm_dx=',crm_dx,'   crm_dy=',crm_dy,'   crm_dt=',crm_dt
      print*,'crm_subx=',orc_nsubdomains_x,'   crm_suby=',orc_nsubdomains_y
      print*,'orc subdomains=',orc_nsubdomains,'    orc_total=',orc_total
      if (is_spcam_sam1mom) print*,'Microphysics: SAM1MOM'
      if (is_spcam_m2005)   print*,'Microphysics: M2005'
      print*,'_________________________________________'
  end if

  if (do_clubb_sgs) then
      call pbuf_add_field('CLUBB_BUFFER','global', dtype_r8, (/pcols,crm_nx,crm_ny,crm_nz+1,nclubbvars/), clubb_buffer_idx)
      call pbuf_add_field('tke',         'global', dtype_r8, (/pcols, pverp/),                            tke_idx)
      call pbuf_add_field('kvm',         'global', dtype_r8, (/pcols, pverp/),                            kvm_idx)
      call pbuf_add_field('kvh',         'global', dtype_r8, (/pcols, pverp/),                            kvh_idx)
      call pbuf_add_field('pblh',        'global', dtype_r8, (/pcols, pverp/),                            pblh_idx)
      call pbuf_add_field('tpert',       'global', dtype_r8, (/pcols, pverp/),                            tpert_idx)
  end if

  call setparm()

  call pbuf_add_field('CRM_U',       'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_u_idx)
  call pbuf_add_field('CRM_V',       'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_v_idx)
  call pbuf_add_field('CRM_W',       'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_w_idx)
  call pbuf_add_field('CRM_T',       'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_t_idx)
  call pbuf_add_field('CLDO',        'global',  dtype_r8, (/pcols, pver, dyn_time_lvls/),                cldo_idx)
  call pbuf_add_field('CLD',         'global',  dtype_r8, (/pcols, pver, dyn_time_lvls/),                cld_idx)
  call pbuf_add_field('AST',         'global',  dtype_r8, (/pcols, pver, dyn_time_lvls/),                ast_idx)

  call pbuf_add_field('CRM_T_RAD',   'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),             crm_t_rad_idx)
  call pbuf_add_field('CRM_QV_RAD',  'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),             crm_qv_rad_idx)
  call pbuf_add_field('CRM_QC_RAD',  'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),             crm_qc_rad_idx)
  call pbuf_add_field('CRM_QI_RAD',  'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),             crm_qi_rad_idx)
  call pbuf_add_field('CRM_CLD_RAD', 'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),             crm_cld_rad_idx)
  call pbuf_add_field('CRM_QRAD',    'global',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz/),              crm_qrad_idx)

  call pbuf_add_field('PREC_DP',     'physpkg', dtype_r8, (/pcols/),                                     prec_dp_idx)
  call pbuf_add_field('SNOW_DP',     'physpkg', dtype_r8, (/pcols/),                                     snow_dp_idx)
  call pbuf_add_field('PREC_SH',     'physpkg', dtype_r8, (/pcols/),                                     prec_sh_idx)
  call pbuf_add_field('SNOW_SH',     'physpkg', dtype_r8, (/pcols/),                                     snow_sh_idx)
  call pbuf_add_field('PREC_SED',    'physpkg', dtype_r8, (/pcols/),                                     prec_sed_idx)
  call pbuf_add_field('SNOW_SED',    'physpkg', dtype_r8, (/pcols/),                                     snow_sed_idx)
  call pbuf_add_field('PREC_PCW',    'physpkg', dtype_r8, (/pcols/),                                     prec_pcw_idx)
  call pbuf_add_field('SNOW_PCW',    'physpkg', dtype_r8, (/pcols/),                                     snow_pcw_idx)
  call pbuf_add_field('CLDTOP',      'physpkg', dtype_r8, (/pcols,1/),                                   cldtop_idx )
  call pbuf_add_field('RPRDTOT',     'physpkg' ,dtype_r8, (/pcols,pver/),                                rprdtot_idx )
  call pbuf_add_field('ICWMRSH',     'physpkg' ,dtype_r8, (/pcols,pver/),                                icwmrsh_idx )
  call pbuf_add_field('RPRDSH',      'physpkg' ,dtype_r8, (/pcols,pver/),                                rprdsh_idx )
  call pbuf_add_field('NEVAPR_SHCU', 'physpkg' ,dtype_r8, (/pcols,pver/),                                nevapr_shcu_idx )
  call pbuf_add_field('ICWMRDP',     'physpkg', dtype_r8, (/pcols,pver/),                                icwmrdp_idx)
  call pbuf_add_field('RPRDDP',      'physpkg', dtype_r8, (/pcols,pver/),                                rprddp_idx)
  call pbuf_add_field('NEVAPR_DPCU', 'physpkg', dtype_r8, (/pcols,pver/),                                nevapr_dpcu_idx)
  call pbuf_add_field('REI',         'physpkg', dtype_r8, (/pcols,pver/),                                rei_idx)
  call pbuf_add_field('REL',         'physpkg', dtype_r8, (/pcols,pver/),                                rel_idx)
  call pbuf_add_field('NEVAPR',      'physpkg', dtype_r8, (/pcols,pver/),                                nevapr_idx)
  call pbuf_add_field('PRAIN',       'physpkg', dtype_r8, (/pcols,pver/),                                prain_idx)
  call pbuf_add_field('WSEDL',       'physpkg', dtype_r8, (/pcols,pver/),                                wsedl_idx)
  call pbuf_add_field('QME',         'physpkg', dtype_r8, (/pcols,pver/),                                qme_idx)
  call pbuf_add_field('DEI',         'physpkg', dtype_r8, (/pcols,pver/),                                dei_idx)
  call pbuf_add_field('DES',         'physpkg', dtype_r8, (/pcols,pver/),                                des_idx)
  call pbuf_add_field('MU',          'physpkg', dtype_r8, (/pcols,pver/),                                mu_idx)
  call pbuf_add_field('LAMBDAC',     'physpkg', dtype_r8, (/pcols,pver/),                                lambdac_idx)
  call pbuf_add_field('CMFMC_SH',    'physpkg' ,dtype_r8, (/pcols,pverp/),                               cmfmc_sh_idx )

  call pbuf_add_field('FICE',        'physpkg', dtype_r8, (/pcols,pver/),                                fice_idx)

  if (prog_modal_aero) then
     call pbuf_add_field('RATE1_CW2PR_ST','physpkg', dtype_r8, (/pcols,pver/),                            rate1_cw2pr_st_idx)
     call pbuf_add_field('CRM_QAERWAT',   'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz, nmodes/), crm_qaerwat_idx)
     call pbuf_add_field('CRM_DGNUMWET',  'physpkg',  dtype_r8, (/pcols,crm_nx, crm_ny, crm_nz, nmodes/), crm_dgnumwet_idx)
  endif

  if (is_spcam_m2005) then
    call pbuf_add_field('CRM_NC_RAD',     'physpkg',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),        crm_nc_rad_idx)
    call pbuf_add_field('CRM_NI_RAD',     'physpkg',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),        crm_ni_rad_idx)
    call pbuf_add_field('CRM_QS_RAD',     'physpkg',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),        crm_qs_rad_idx)
    call pbuf_add_field('CRM_NS_RAD',     'physpkg',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),        crm_ns_rad_idx)

    ! Fields for crm_micro array
    call pbuf_add_field('CRM_QT',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qt_idx)
    call pbuf_add_field('CRM_NC',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_nc_idx)
    call pbuf_add_field('CRM_QR',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qr_idx)
    call pbuf_add_field('CRM_NR',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_nr_idx)
    call pbuf_add_field('CRM_QI',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qi_idx)
    call pbuf_add_field('CRM_NI',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_ni_idx)
    call pbuf_add_field('CRM_QS',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qs_idx)
    call pbuf_add_field('CRM_NS',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_ns_idx)
    call pbuf_add_field('CRM_QG',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qg_idx)
    call pbuf_add_field('CRM_NG',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_ng_idx)
    call pbuf_add_field('CRM_QC',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qc_idx)
  else
    call pbuf_add_field('CRM_QT',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qt_idx)
    call pbuf_add_field('CRM_QP',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qp_idx)
    call pbuf_add_field('CRM_QN',         'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/),         crm_qn_idx)
  endif

   
  if (is_spcam_m2005) then
     call pbuf_add_field('TK_CRM',        'global',  dtype_r8, (/pcols, pver/),                           tk_crm_idx)
     ! total (all sub-classes) cloudy fractional area in previous time step 
     call pbuf_add_field('ACLDY_CEN',     'global',  dtype_r8, (/pcols,pver/),                            acldy_cen_idx) 
  endif 

! Adding crm dimensions to cam history
  call add_hist_coord('crm_nx'       ,crm_nx,  'CRM NX')
  call add_hist_coord('crm_ny'       ,crm_ny,  'CRM NY')
  call add_hist_coord('crm_nz'       ,crm_nz,  'CRM NZ')
  call add_hist_coord('crm_z1'       ,crm_nz+1,'CRM_Z1')

  call add_hist_coord('pverp'        ,pverp,     'pverp ')
  call add_hist_coord('pver'         ,pver,      'pver  ')

! ifdef needed because of NCLASS_CL
#ifdef m2005
  call add_hist_coord('NCLASS_CL'    ,NCLASS_CL,'NCLASS_CL')
  call add_hist_coord('ncls_ecpp_in' ,ncls_ecpp_in,'ncls_ecpp_in')
  call add_hist_coord('NCLASS_PR'    ,NCLASS_PR,'NCLASS_PR')
#endif

#endif

end subroutine crm_physics_register
!=========================================================================================================

subroutine crm_physics_init(pbuf2d)
!-------------------------------------------------------------------------------------------------------
! 
! Purpose: initialize some variables, and add necessary fileds into output fields 
!
!--------------------------------------------------------------------------------------------------------
  use physics_buffer,  only: physics_buffer_desc, pbuf_set_field, pbuf_get_index
#ifdef CRM
  use physconst,       only: tmelt, rair, cpair, rh2o, latvap, latice
  use constituents,    only: pcnst, cnst_species_class, cnst_spec_class_gas
  use cam_history,     only: addfld, add_default, horiz_only
  use crmdims,         only: crm_nx, crm_ny, crm_nz
  use ndrop,           only: ndrop_init
  use gas_wetdep_opts, only: gas_wetdep_method
  use micro_mg_utils,  only: micro_mg_utils_init
  use time_manager,    only: is_first_step

    use cam_history,     only: fieldname_len
#ifdef MODAL_AERO
    use modal_aero_data, only: cnst_name_cw, ntot_amode, &
                               lmassptr_amode, lmassptrcw_amode, &
                               nspec_amode, numptr_amode, numptrcw_amode
#endif

#endif

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
       
#ifdef CRM
   integer :: l, lphase, lspec
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit

! local variables
  integer :: i, m, mm
  integer :: icldphy                 ! index for cloud physic species (water vapor and cloud hydrometers)

  character(len=128):: errstring     ! return status (non-blank for error return)

  crm_nx_ny      = crm_nx*crm_ny

  !------------------------- 
  ! Make sure gas_wetdep_method is set to 'MOZ' as 'NEU' is not currently supported by SPCAM
  !  'MOZ' for spcam_sam1mom
  !  'OFF' for spcam_m2005
  if (is_spcam_sam1mom) then
     if (gas_wetdep_method /= 'MOZ') call endrun( "crm_physics: gas_wetdep_method must be set to 'MOZ' ")
  elseif (is_spcam_m2005) then
     if (gas_wetdep_method /= 'OFF') call endrun( "crm_physics: gas_wetdep_method must be set to 'OFF' ")
  else
     call endrun( "crm_physics: don't know how gas_wetdep_method should be set")
  endif

  !------------------------- 
  ! Initialize the micro_mg_utils
  ! Value of dcs in MG 1.0 is 400.e-6_r8
  call micro_mg_utils_init(r8, rair, rh2o, cpair, tmelt, latvap, latice, 400.e-6_r8,  errstring)

  !------------------------- 
  ! Register general history fields
  do m = 1, ncnst_use
     call cnst_get_ind(cnst_names(m), mm)
     if ( any(mm == (/ ixcldliq, ixcldice /)) ) then
        ! mass mixing ratios
        call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg   ', cnst_longname(mm))
        call addfld(sflxnam(mm),   horiz_only,  'A', 'kg/m2/s ', trim(cnst_name(mm))//' surface flux')
     else if ( any(mm == (/ ixnumliq, ixnumice /)) ) then
        ! number concentrations
        call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg    ', cnst_longname(mm))
        call addfld(sflxnam(mm),   horiz_only,  'A', '1/m2/s  ', trim(cnst_name(mm))//' surface flux')
     else
        call endrun( "crm_physics: Could not call addfld for constituent with unknown units.")
     endif
  end do

  do m=1, pcnst 
    if(cnst_name(m) == 'DMS') then 
       call addfld('DMSCONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'DMS tendency from ZM convection')
    end if
    if(cnst_name(m) == 'SO2') then 
       call addfld('SO2CONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'SO2 tendency from ZM convection')
     end if
  end do

  call addfld ('CRM_TK',    (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm^2/s',  'Eddy viscosity from CRM')
  call addfld ('CRM_TKH',   (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm^2/s',  'Eddy viscosity from CRM')

  call addfld ('SPCLD3D  ', (/ 'lev' /), 'A', 'fraction', 'cloud fraction on GCM grids')
  call addfld ('MU_CRM   ', (/ 'lev' /), 'A', 'Pa/s',     'mass flux up from CRM')
  call addfld ('MD_CRM   ', (/ 'lev' /), 'A', 'Pa/s',     'mass flux down from CRM')
  call addfld ('DU_CRM   ', (/ 'lev' /), 'A', '/s',       'detrainment from updraft from CRM')
  call addfld ('EU_CRM   ', (/ 'lev' /), 'A', '/s',       'entraiment rate from updraft')
  call addfld ('ED_CRM   ', (/ 'lev' /), 'A', '/s',       'entraiment rate from downdraft')
  call addfld ('SPQRL    ', (/ 'lev' /), 'A', 'K/s',      'long-wave heating rate')
  call addfld ('SPQRS    ', (/ 'lev' /), 'A', 'K/s',      'short-wave heating rate')
  call addfld ('LENGC    ', (/ 'ilev' /), 'A', 'm  ',      'Mixing length scale for the calcuation of vertical difusivity')

  call addfld ('SPKVH     ',(/ 'ilev' /), 'A', 'm2/s    ', 'Vertical diffusivity used in dropmixnuc in the MMF call')
  call addfld ('SPLCLOUD  ',(/ 'lev' /), 'A', '        ', 'Liquid cloud fraction')
  call add_default ('SPKVH     ', 1, ' ')
  call add_default ('SPLCLOUD  ', 1, ' ')

  call addfld ('SPCLDTOT', horiz_only, 'A', 'fraction',    'Vertically-integrated total cloud from CRM'     )
  call addfld ('SPCLDLOW', horiz_only, 'A', 'fraction',    'Vertically-integrated low cloud from CRM'       )
  call addfld ('SPCLDMED', horiz_only, 'A', 'fraction',    'Vertically-integrated mid-level cloud from CRM' )
  call addfld ('SPCLDHGH', horiz_only, 'A', 'fraction',    'Vertically-integrated high cloud from CRM'      )
  call add_default ('SPCLDTOT', 1, ' ')
  call add_default ('SPCLDLOW', 1, ' ')
  call add_default ('SPCLDMED', 1, ' ')
  call add_default ('SPCLDHGH', 1, ' ')

  call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg   ', trim(cnst_name(ixcldliq))//' after physics'  )
  call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg   ', trim(cnst_name(ixcldliq))//' before physics' )
  call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg   ', trim(cnst_name(ixcldice))//' after physics'  )
  call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg   ', trim(cnst_name(ixcldice))//' before physics' )

  call addfld ('PRES    ',(/ 'lev' /), 'A', 'Pa      ','Pressure'                                )
  call addfld ('DPRES   ',(/ 'lev' /), 'A', 'Pa      ','Pressure thickness of layer'             )
  call addfld ('SPDT    ',(/ 'lev' /), 'A', 'K/s     ','T tendency due to CRM'                   )
  call addfld ('SPDQ    ',(/ 'lev' /), 'A', 'kg/kg/s ','Q tendency due to CRM'                   )
  call addfld ('SPDQC   ',(/ 'lev' /), 'A', 'kg/kg/s ','QC tendency due to CRM'                  )
  call addfld ('SPDQI   ',(/ 'lev' /), 'A', 'kg/kg/s ','QI tendency due to CRM'                  )
  call addfld ('SPMC    ',(/ 'lev' /), 'A', 'kg/m2/s ','Total mass flux from CRM'                )
  call addfld ('SPMCUP  ',(/ 'lev' /), 'A', 'kg/m2/s ','Updraft mass flux from CRM'              )
  call addfld ('SPMCDN  ',(/ 'lev' /), 'A', 'kg/m2/s ','Downdraft mass flux from CRM'            )
  call addfld ('SPMCUUP ',(/ 'lev' /), 'A', 'kg/m2/s ','Unsaturated updraft mass flux from CRM'  )
  call addfld ('SPMCUDN ',(/ 'lev' /), 'A', 'kg/m2/s ','Unsaturated downdraft mass flux from CRM')
  call addfld ('SPQC    ',(/ 'lev' /), 'A', 'kg/kg   ','Cloud water from CRM'                    )
  call addfld ('SPQI    ',(/ 'lev' /), 'A', 'kg/kg   ','Cloud ice from CRM'                      )
  call addfld ('SPQS    ',(/ 'lev' /), 'A', 'kg/kg   ','Snow from CRM'                           )
  call addfld ('SPQG    ',(/ 'lev' /), 'A', 'kg/kg   ','Graupel from CRM'                        )
  call addfld ('SPQR    ',(/ 'lev' /), 'A', 'kg/kg   ','Rain from CRM'                           )
  call addfld ('SPQTFLX ',(/ 'lev' /), 'A', 'kg/m2/s ','Nonprecip. water flux from CRM'          )
  call addfld ('SPUFLX  ',(/ 'lev' /), 'A', 'm2/s2   ','x-momentum flux from CRM'                )
  call addfld ('SPVFLX  ',(/ 'lev' /), 'A', 'm2/s2   ','y-momentum flux from CRM'                )
  call addfld ('SPQTFLXS',(/ 'lev' /), 'A', 'kg/m2/s ','SGS Nonprecip. water flux from CRM'      )
  call addfld ('SPTKE   ',(/ 'lev' /), 'A', 'kg/m/s2 ','Total TKE in CRM'                        )
  call addfld ('SPTKES  ',(/ 'lev' /), 'A', 'kg/m/s2 ','SGS TKE in CRM'                          )
  call addfld ('SPTK    ',(/ 'lev' /), 'A', 'm2/s    ','SGS TK in CRM'                           )
  call addfld ('SPQPFLX ',(/ 'lev' /), 'A', 'kg/m2/s ','Precip. water flux from CRM'             )
  call addfld ('SPPFLX  ',(/ 'lev' /), 'A', 'm/s     ','Precipitation flux from CRM'             )
  call addfld ('SPQTLS  ',(/ 'lev' /), 'A', 'kg/kg/s ','L.S. Vapor Tendency from CRM'            )
  call addfld ('SPQTTR  ',(/ 'lev' /), 'A', 'kg/kg/s ','Nonprec. water transport from CRM'       )
  call addfld ('SPQPTR  ',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water transport from CRM'          )
  call addfld ('SPQPEVP ',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water evaporation from CRM'        )
  call addfld ('SPQPFALL',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water fall-out from CRM'           )
  call addfld ('SPQPSRC ',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water source from CRM'             )
  call addfld ('SPTLS   ',(/ 'lev' /), 'A', 'kg/kg/s ','L.S. LIWSE Tendency from CRM'            )
  call addfld ('TIMINGF ', horiz_only, 'A', '        ','CRM CPU usage efficiency: 1 - ideal'     )
  call addfld ('CLOUDTOP',(/ 'lev' /), 'A', '        ','Cloud Top PDF'                           )

  !------------------------- 
  ! Register m2005 history fields
  if (is_spcam_m2005) then
     call addfld ('SPNC    ',(/ 'lev' /), 'A', '/kg   ','Cloud water dropet number from CRM')
     call addfld ('SPNI    ',(/ 'lev' /), 'A', '/kg   ','Cloud ice crystal number from CRM')
     call addfld ('SPNS    ',(/ 'lev' /), 'A', '/kg   ','Snow particle number from CRM')
     call addfld ('SPNG    ',(/ 'lev' /), 'A', '/kg   ','Graupel particle number from CRM')
     call addfld ('SPNR    ',(/ 'lev' /), 'A', '/kg   ','Rain particle number from CRM')
     call add_default ('SPNC    ', 1, ' ')
     call add_default ('SPNI    ', 1, ' ')
     call add_default ('SPNS    ', 1, ' ')
     call add_default ('SPNG    ', 1, ' ')
     call add_default ('SPNR    ', 1, ' ')

     call addfld ('CRM_FLIQ ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Liquid'      )
     call addfld ('CRM_FICE ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Ice'         )
     call addfld ('CRM_FRAIN',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Rain'        )
     call addfld ('CRM_FSNOW',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Snow'        )
     call addfld ('CRM_FGRAP',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Graupel'     )
     call addfld ('CRM_QS  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Snow mixing ratio from CRM'             )
     call addfld ('CRM_QG  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Graupel mixing ratio from CRM'          )
     call addfld ('CRM_QR  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Rain mixing ratio from CRM'             )

     call addfld ('CRM_NC  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Cloud water dropet number from CRM'     )
     call addfld ('CRM_NI  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Cloud ice crystal number from CRM'      )
     call addfld ('CRM_NS  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Snow particle number from CRM'          )
     call addfld ('CRM_NG  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Graupel particle number from CRM'       )
     call addfld ('CRM_NR  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Rain particle number from CRM'          )

     ! below is for *instantaneous* crm output
     call addfld ('CRM_AUT  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Autoconversion cloud waterfrom CRM'     )
     call addfld ('CRM_ACC  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Accretion cloud water from CRM'         )
     call addfld ('CRM_EVPC ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Evaporation cloud water from CRM'       )
     call addfld ('CRM_EVPR ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Evaporation rain from CRM'              )
     call addfld ('CRM_MLT  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Melting ice snow graupel from CRM'      )
     call addfld ('CRM_SUB  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Sublimation ice snow graupel from CRM'  )
     call addfld ('CRM_DEP  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Deposition ice snow graupel from CRM'   )
     call addfld ('CRM_CON  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Condensation cloud water from CRM'      )

     ! below is for *gcm-grid and time-step-avg* process output
     call addfld ('A_AUT  ',(/ 'lev' /), 'A', '/s   ','Avg autoconversion cloud water from CRM'            )
     call addfld ('A_ACC  ',(/ 'lev' /), 'A', '/s   ','Avg accretion cloud water from CRM'                 )
     call addfld ('A_EVPC ',(/ 'lev' /), 'A', '/s   ','Avg evaporation cloud water from CRM'               )
     call addfld ('A_EVPR ',(/ 'lev' /), 'A', '/s   ','Avg evaporation rain from CRM'                      )
     call addfld ('A_MLT  ',(/ 'lev' /), 'A', '/s   ','Avg melting ice snow graupel from CRM'              )
     call addfld ('A_SUB  ',(/ 'lev' /), 'A', '/s   ','Avg sublimation ice snow graupel from CRM'          )
     call addfld ('A_DEP  ',(/ 'lev' /), 'A', '/s   ','Avg deposition ice snow graupel from CRM'           )
     call addfld ('A_CON  ',(/ 'lev' /), 'A', '/s   ','Avg condensation cloud water from CRM'              )

     call addfld ('CRM_REL  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale droplet effective radius')
     call addfld ('CRM_REI  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale ice crystal effective radius')
     call addfld ('CRM_DEI  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale Mitchell ice effective diameter')
     call addfld ('CRM_DES  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale snow effective diameter')
     call addfld ('CRM_MU   ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', &
                                                 'cloud scale droplet size distribution shape parameter for radiation')
     call addfld ('CRM_LAMBDA',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers',  &
                                                 'cloud scale slope of droplet distribution for radiation')
     call addfld ('CRM_TAU  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', '1',           'cloud scale cloud optical depth'  )
     call addfld ('CRM_WVAR' , (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm/s',         'vertical velocity variance from CRM')

     call addfld ('CRM_FSNT',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA shortwave fluxes at CRM grids')
     call addfld ('CRM_FSNTC', (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA clear-sky shortwave fluxes at CRM grids')
     call addfld ('CRM_FSNS',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net surface shortwave fluxes at CRM grids')
     call addfld ('CRM_FSNSC', (/'crm_nx','crm_ny'/),          'A',  'unitless',  &
                                                 'net surface clear-sky shortwave fluxes at CRM grids')
     call addfld ('CRM_FLNT',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA longwave fluxes at CRM grids')
     call addfld ('CRM_FLNTC', (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net TOA clear-sky longwave fluxes at CRM grids')
     call addfld ('CRM_FLNS',  (/'crm_nx','crm_ny'/),          'A',  'unitless', 'net surface longwave fluxes at CRM grids')
     call addfld ('CRM_FLNSC', (/'crm_nx','crm_ny'/),          'A',  'unitless',  &
                                                 'net surface clear-sky longwave fluxes at CRM grids')

     call addfld ('CRM_AODVIS', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 550nm in CRM grids',&
                                                                                                               flag_xyfill=.true.)
     call addfld ('CRM_AOD400', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 400nm in CRM grids',&
                                                                                                               flag_xyfill=.true.)
     call addfld ('CRM_AOD700', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 700nm in CRM grids', &
                                                                                                               flag_xyfill=.true.)
     call addfld ('CRM_AODVISZ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'unitless',  &
                                                  'Aerosol optical depth at each layer at 500nm in CRM grids', flag_xyfill=.true.)
     call addfld ('AOD400',      horiz_only,                    'A', 'unitless', 'Aerosol optical depth at 400nm', &
                                                                                                               flag_xyfill=.true.)
     call addfld ('AOD700',      horiz_only,                    'A', 'unitless', 'Aerosol optical depth at 700nm', &
                                                                                                               flag_xyfill=.true.)
     call add_default ('AOD400',  1, ' ')
     call add_default ('AOD700',  1, ' ')
  endif

  !------------------------- 
  ! Register CLUBB history fields
  if (do_clubb_sgs) then
     call addfld ('UP2     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'u prime ^2 from clubb')
     call addfld ('VP2     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'v prime ^2 from clubb')
     call addfld ('WPRTP   ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'mkg/skg',    'w prime * rt prime from clubb')
     call addfld ('WPTHLP  ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'mK/s',       'w prime * th_l prime from clubb')
     call addfld ('WP2     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'w prime ^2 from clubb')
     call addfld ('WP3     ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^3/s^3',    'w prime ^3 from clubb')
     call addfld ('RTP2    ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', '(kg/kg)2',   'r_t prime ^2 from clubb')
     call addfld ('THLP2   ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'K^2',        'th_l_prime ^2 from clubb')
     call addfld ('RTPTHLP ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'kgK/kg',     'r_t prime * th_l prime  from clubb')
     call addfld ('UPWP    ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'u prime * w prime from clubb')
     call addfld ('VPWP    ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'm^2/s^2',    'v prime * w prime from clubb')
     call addfld ('CRM_CLD ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'fraction',   'cloud fraction from clubb')
     call addfld ('T_TNDCY ',    (/'crm_nx','crm_ny','crm_z1'/), 'A', 'K/s',        't tendency from clubb')
     call addfld ('QV_TNDCY ',   (/'crm_nx','crm_ny','crm_z1'/), 'A', 'kg/kg/s',    'water vapor tendency from clubb')
     call addfld ('QC_TNDCY ',   (/'crm_nx','crm_ny','crm_z1'/), 'A', 'kg/kg/s',    'liquid condensate tendency from clubb')
     call addfld ('CLUBB_TK',    (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm^2/s',      'Eddy viscosity from clubb')
     call addfld ('CLUBB_TKH',   (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm^2/s',      'Eddy viscosity from clubb')
     call addfld ('CRM_RELVAR',  (/'crm_nx','crm_ny','crm_nz'/), 'A',   '',         'cloud water relative variance from clubb')
     call addfld ('ACCRE_ENHAN', (/'crm_nx','crm_ny','crm_nz'/), 'A', '',           'Accretion enhancment from clubb')
     call addfld ('QCLVAR',      (/'crm_nx','crm_ny','crm_nz'/), 'A', '(kg/kg)^2',  'cloud water variance from clubb')
     ! add GCM-scale output
     call addfld ('SPUP2',        (/ 'lev' /), 'A', 'm^2/s^2',  'u prime ^2 from clubb on GCM grids')
     call addfld ('SPVP2',        (/ 'lev' /), 'A', 'm^2/s^2',  'v prime ^2 from clubb on GCM grids')
     call addfld ('SPWPRTP',      (/ 'lev' /), 'A', 'mkg/skg',  'w prime * rt prime from clubb on GCM grids')
     call addfld ('SPWPTHLP',     (/ 'lev' /), 'A', 'mK/s',     'w prime * th_l prime from clubb on GCM grids')
     call addfld ('SPWP2',        (/ 'lev' /), 'A', 'm^2/s^2',  'w prime ^2 from clubb on GCM grids')
     call addfld ('SPWP3',        (/ 'lev' /), 'A', 'm^3/s^3',  'w prime ^3 from clubb on GCM grids')
     call addfld ('SPRTP2',       (/ 'lev' /), 'A', '(kg/kg)2', 'r_t prime ^2 from clubb on GCM grids')
     call addfld ('SPTHLP2',      (/ 'lev' /), 'A', 'K^2',      'th_l_prime ^2 from clubb on GCM grids')
     call addfld ('SPRTPTHLP',    (/ 'lev' /), 'A', 'kgK/kg',   'r_t prime * th_l prime  from clubb on GCM grids')
     call addfld ('SPUPWP',       (/ 'lev' /), 'A', 'm^2/s^2',  'u prime * w prime from clubb on GCM grids')
     call addfld ('SPVPWP',       (/ 'lev' /), 'A', 'm^2/s^2',  'v prime * w prime from clubb on GCM grids')
     call addfld ('SPCRM_CLD ',   (/ 'lev' /), 'A', 'fraction', 'cloud fraction from clubb on GCM grids')
     call addfld ('SPT_TNDCY ',   (/ 'lev' /), 'A', 'K/s',      't tendency from clubb on GCM grids')
     call addfld ('SPQV_TNDCY ',  (/ 'lev' /), 'A', 'kg/kg/s',  'water vapor tendency from clubb on GCM grids')
     call addfld ('SPQC_TNDCY ',  (/ 'lev' /), 'A', 'kg/kg/s',  'liquid condensate tendency from clubb on GCM grids')
     call addfld ('SPCLUBB_TK',   (/ 'lev' /), 'A', 'm^2/s',    'Eddy viscosity from clubb on GCM grids')
     call addfld ('SPCLUBB_TKH',  (/ 'lev' /), 'A', 'm^2/s',    'Eddy viscosity from clubb on GCM grids')
     call addfld ('SPRELVAR',     (/ 'lev' /), 'A', '',         'cloud water relative variance from clubb on GCM grids')
     call addfld ('SPACCRE_ENHAN',(/ 'lev' /), 'A', '',         'Accretion enhancment from clubb on GCM grids')
     call addfld ('SPQCLVAR',     (/ 'lev' /), 'A', '',         'cloud water variance from clubb on GCM grids')
  endif


  !------------------------- 
  ! Register ECPP history fields
  ! ifdef needed because of ECPP parameters such as NCLASS_CL and ncls_ecpp_in and papampollu_init 
#ifdef m2005
  if (is_spcam_m2005) then

     call papampollu_init ()

     call addfld ('ABND    ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for full time period at layer boundary')
     call addfld ('ABND_TF ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for end-portion of time period at layer boundary')
     call addfld ('MASFBND ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                  'sub-class vertical mass flux (kg/m2/s) at layer boundary')
     call addfld ('ACEN    ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for full time period at layer center')
     call addfld ('ACEN_TF ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'cloud fraction for each sub-sub class for end-portion of time period at layer center')
     call addfld ('RHCEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                  'relative humidity for each sub-sub calss at layer center')
     call addfld ('QCCEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg',    &
                  'cloud water for each sub-sub class at layer center')
     call addfld ('QICEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg',    &
                  'cloud ice for each sub-sub class at layer center')
     call addfld ('QSINK_AFCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',    &
                  'cloud water loss rate from precip. using cloud water after precip. for each sub-sub class at layer center')
     call addfld ('QSINK_BFCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',    &
                 'cloud water loss rate from precip. using cloud water before precip. for each sub-sub class at layer center')
     call addfld ('QSINK_AVGCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',   &
      'cloud water loss rate from precip. using averaged cloud water and precip. rate for each sub-sub class at layer center')
     call addfld ('PRAINCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg/s',  &
                  ' cloud water loss rate from precipitation (kg/kg/s) for each sub-sub class at layer center')
     call addfld ('PRECRCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                  'liquid (rain) precipitation rate for each sub-sub class at layer center')
     call addfld ('PRECSCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                  'solid (snow, graupel,...) precipitation rate for each sub-sub class at layer center')
     call addfld ('WUPTHRES',      (/ 'ilev' /), 'A', 'm/s',    'vertical velocity threshold for updraft')
     call addfld ('WDNTHRES',      (/ 'ilev' /), 'A', 'm/s',    'vertical velocity threshold for dndraft')
     call addfld ('WWQUI_CEN',     (/ 'lev' /),  'A', 'm2/s2',  'vertical velocity variance in the quiescent class, layer center')
     call addfld ('WWQUI_CLD_CEN', (/ 'lev' /),  'A', 'm2/s2',  &
                                                      'vertical velocity variance in the cloudy quiescent class, layer center')
     call addfld ('WWQUI_BND',     (/ 'ilev' /), 'A', 'm2/s2',  &
                                                      'vertical velocity variance in the quiescent class, layer boundary')
     call addfld ('WWQUI_CLD_BND', (/ 'ilev' /), 'A', 'm2/s2',  &
                                                      'vertical velocity variance in the cloudy quiescent class, layer boundary')
  endif
#endif

  !------------------------- 
  ! Register modal aerosol history fields
  ! ifdef needed because of use of cnst_name_cw which not defined if not modal aerosols
#ifdef MODAL_AERO
  if (prog_modal_aero) then

    call ndrop_init()

    do m=1, pcnst
       if(cnst_species_class(m).eq.cnst_spec_class_gas) then
          fieldname = trim(cnst_name(m)) // '_mixnuc1sp'
          long_name = trim(cnst_name(m)) // ' dropmixnuc mixnuc column tendency in the mmf one '
          call addfld( fieldname,  horiz_only, 'A', unit, long_name)
          call add_default( fieldname, 1, ' ' )
       end if
    end do

  endif

#endif
  
  ! These variables do not vary in CRM
  call pbuf_set_field (pbuf2d, prec_dp_idx,   0.0_r8)
  call pbuf_set_field (pbuf2d, prec_sh_idx,   0.0_r8)
  call pbuf_set_field (pbuf2d, snow_sh_idx,   0.0_r8)
  call pbuf_set_field (pbuf2d, snow_dp_idx,   0.0_r8)
  call pbuf_set_field (pbuf2d, prec_sed_idx,  0.0_r8)
  call pbuf_set_field (pbuf2d, snow_sed_idx,  0.0_r8)
  call pbuf_set_field (pbuf2d, prec_pcw_idx,  0.0_r8)
  call pbuf_set_field (pbuf2d, snow_pcw_idx,  0.0_r8)

  
  call addfld ('CRM_U   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'm/s     ', 'CRM x-wind'                          )
  call addfld ('CRM_V   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'm/s     ', 'CRM y-wind'                          )
  call addfld ('CRM_W   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'm/s     ', 'CRM z-wind'                          )
  call addfld ('CRM_T   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'K       ', 'CRM Temperature'                     )
  call addfld ('CRM_QV  ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Water Vapor'                     )
  call addfld ('CRM_QC  ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Cloud Water'                     )
  call addfld ('CRM_QI  ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Cloud Ice'                       )
  call addfld ('CRM_QPC ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Precipitating Water'             )
  call addfld ('CRM_QPI ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Precipitating Ice'               )
  call addfld ('CRM_PREC',(/'crm_nx','crm_ny'/),           'I', 'm/s     ', 'CRM Precipitation Rate'              )
  call addfld ('CRM_QRS ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'K/s     ', 'CRM Shortwave radiative heating rate')
  call addfld ('CRM_QRL ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'K/s     ', 'CRM Longwave radiative heating rate' )

  call add_default ('SPDT    ', 1, ' ')
  call add_default ('SPDQ    ', 1, ' ')
  call add_default ('SPDQC   ', 1, ' ')
  call add_default ('SPDQI   ', 1, ' ')
  call add_default ('SPMC    ', 1, ' ')
  call add_default ('SPMCUP  ', 1, ' ')
  call add_default ('SPMCDN  ', 1, ' ')
  call add_default ('SPMCUUP ', 1, ' ')
  call add_default ('SPMCUDN ', 1, ' ')
  call add_default ('SPQC    ', 1, ' ')
  call add_default ('SPQI    ', 1, ' ')
  call add_default ('SPQS    ', 1, ' ')
  call add_default ('SPQG    ', 1, ' ')
  call add_default ('SPQR    ', 1, ' ')
  call add_default ('SPQTFLX ', 1, ' ')
  call add_default ('SPQTFLXS', 1, ' ')
  call add_default ('SPTKE   ', 1, ' ')
  call add_default ('SPTKES  ', 1, ' ')
  call add_default ('SPTK    ', 1, ' ')
  call add_default ('SPQPFLX ', 1, ' ')
  call add_default ('SPPFLX  ', 1, ' ')
  call add_default ('SPQTLS  ', 1, ' ')
  call add_default ('SPQTTR  ', 1, ' ')
  call add_default ('SPQPTR  ', 1, ' ')
  call add_default ('SPQPEVP ', 1, ' ')
  call add_default ('SPQPFALL', 1, ' ')
  call add_default ('SPQPSRC ', 1, ' ')
  call add_default ('SPTLS   ', 1, ' ')
  call add_default ('CLOUDTOP', 1, ' ')
  call add_default ('TIMINGF ', 1, ' ')

  sh_frac_idx = pbuf_get_index('SH_FRAC')
  dp_frac_idx = pbuf_get_index('DP_FRAC')
  call pbuf_set_field (pbuf2d, sh_frac_idx,   0.0_r8)
  call pbuf_set_field (pbuf2d, dp_frac_idx,   0.0_r8)

  call pbuf_set_field (pbuf2d, cmfmc_sh_idx,    0.0_r8)
  call pbuf_set_field (pbuf2d, rprdsh_idx,      0.0_r8)
  call pbuf_set_field (pbuf2d, icwmrsh_idx,     0.0_r8)
  call pbuf_set_field (pbuf2d, nevapr_shcu_idx, 0.0_r8)

  call pbuf_set_field (pbuf2d, icwmrdp_idx,     0.0_r8)
  call pbuf_set_field (pbuf2d, fice_idx,        0.0_r8)
 
  call pbuf_set_field (pbuf2d, prain_idx,       0.0_r8)
  call pbuf_set_field (pbuf2d, rprdtot_idx,     0.0_r8)
  call pbuf_set_field (pbuf2d, nevapr_idx,      0.0_r8)

  if (is_first_step()) then
     call pbuf_set_field (pbuf2d, ast_idx,         0.0_r8)
  end if
#endif
end subroutine crm_physics_init

!=========================================================================================================

function crm_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: crm_implements_cnst    ! return value

#ifdef CRM
   !-----------------------------------------------------------------------

   crm_implements_cnst = any(name == cnst_names)

#endif
end function crm_implements_cnst

!===============================================================================

subroutine crm_init_cnst(name, q)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   !-----------------------------------------------------------------------

#ifdef CRM
   if (crm_implements_cnst(name)) q = 0.0_r8
#endif

end subroutine crm_init_cnst

!===============================================================================

!---------------------------------------------------------------------------------------------------------
   subroutine crm_physics_tend(ztodt, state, tend, ptend, pbuf, cam_in)

!------------------------------------------------------------------------------------------
!  Purpose: to update state from CRM physics. 
! 
! Revision history: 
!
! June, 2009, Minghuai Wang: 
!          These codes are taken out from tphysbc.F90 
!       in the spcam3.5, developed by Marat Khairoutdinov 
!       (mkhairoutdin@ms.cc.sunysb.edu). Here we try to follow the procedure 
!      in 'Interface to Column Physics and Chemistry packages' to implement 
!      the CRM physics.
! July, 13, 2009, Minghuai Wang: 
!      Hydrometer numbers are outputed from SAM when Morrison's microphysics is used, 
!      and will be used in the radiative transfer code to calculate radius. 
! July, 15, 2009, Minghuai Wang: 
!      Get modal aerosol, and use it in the SAM. 
! 
!-------------------------------------------------------------------------------------------
#ifdef CRM
   use crmx_micro_params
   use crmx_params,         only: ggr,cp,fac_cond,fac_sub
   use spmd_utils,          only: masterproc,iam
   use shr_spfn_mod,        only: gamma => shr_spfn_gamma
   use time_manager,        only: is_first_step, get_nstep
   use cam_history,         only: outfld
   use perf_mod
   use crmdims,             only: crm_nx, crm_ny, crm_nz
   use physconst,           only: cpair, latvap, gravit
   use constituents,        only: pcnst, cnst_get_ind
   use crmx_crm_module,     only: crm
   use crmx_microphysics,   only: nmicro_fields
   use physconst,           only: latvap
   use check_energy,        only: check_energy_chng
   use phys_grid,           only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use modal_aero_calcsize, only: modal_aero_calcsize_sub
   use micro_mg_utils,      only: size_dist_param_liq, mg_liq_props, mincld, qsmall
   use phys_grid, only: get_rlon_p, get_rlat_p, get_gcol_all_p

#ifdef MODAL_AERO
     use crmclouds_camaerosols, only: crmclouds_mixnuc_tend, spcam_modal_aero_wateruptake_dr
     use ndrop,                 only: loadaer
#endif
#ifdef m2005
     use module_ecpp_ppdriver2, only: parampollu_driver2
     use crmx_ecppvars,         only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
     use module_data_ecpp1,     only: dtstep_pp_input
#endif
#ifdef SPCAM_CLUBB_SGS
     use cloud_cover_diags, only: cloud_cover_diags_out 
     use pkg_cldoptics,     only: cldovrlap
#endif

#endif 

   use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, dyn_time_lvls, pbuf_get_field
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init, &
                              physics_state_copy, physics_ptend_sum, physics_ptend_scale
   use camsrfexch,      only: cam_in_t

#ifdef ORCHESTRATOR
    use crmdims,         only:orc_rank_total,orc_total,orc_nsubdomains,orc_nx,orc_ny,orc_nsubdomains_x,orc_nsubdomains_y
#endif

   integer, parameter :: nmicro_fields_total = 3
   real(r8), intent(in)              :: ztodt                          ! 2 delta t (model time increment)
   type(physics_state), intent(in)   :: state   
   type(physics_tend), intent(in)    :: tend
   type(physics_ptend ), intent(out) :: ptend
   type(physics_buffer_desc),pointer :: pbuf(:)
   type (cam_in_t), intent(in)       :: cam_in
#ifdef CRM

   type(physics_state) :: state_loc   ! local copy of state   
   type(physics_tend)  :: tend_loc    ! local copy of tend   
   type(physics_ptend) :: ptend_loc   ! local copy of ptend   

   ! convective precipitation variables
   real(r8), pointer :: prec_dp(:)                ! total precipitation from ZM convection [m/s]
   real(r8), pointer :: snow_dp(:)                ! snow from ZM convection                [m/s]

   real(r8), pointer ::  nc_rad(:,:,:,:) ! rad cloud water droplet number [#/kg]
   real(r8), pointer ::  ni_rad(:,:,:,:) ! rad cloud ice crystal number [#/kg]
   real(r8), pointer ::  qs_rad(:,:,:,:) ! rad cloud snow mass [kg/kg]
   real(r8), pointer ::  ns_rad(:,:,:,:) ! rad cloud snow crystal number [#/kg]
   real(r8), pointer ::  cld_rad(:,:,:,:)  ! cloud fraction

   real(r8), pointer ::  t_rad (:,:,:,:) ! rad temperuture
   real(r8), pointer ::  qv_rad(:,:,:,:) ! rad vapor
   real(r8), pointer ::  qc_rad(:,:,:,:) ! rad cloud water
   real(r8), pointer ::  qi_rad(:,:,:,:) ! rad cloud ice
   real(r8), pointer ::  crm_qrad(:,:,:,:)
   real(r8), pointer ::  clubb_buffer  (:,:,:,:,:)

   real(r8),pointer  :: cldtop_pbuf(:)  ! cloudtop location for pbuf

   real(r8),pointer  :: tk_crm_ecpp(:,:)
   real(r8),pointer  :: acldy_cen_tbeg(:,:)        ! cloud fraction
   real(r8), pointer, dimension(:,:) :: cldo
 
!
!--------------------------- Local variables -----------------------------------------------------------
!
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer nstep     ! time steps

   real(r8) qc_crm (pcols,crm_nx, crm_ny, crm_nz)
   real(r8) qi_crm (pcols,crm_nx, crm_ny, crm_nz)
   real(r8) qpc_crm(pcols,crm_nx, crm_ny, crm_nz)
   real(r8) qpi_crm(pcols,crm_nx, crm_ny, crm_nz)

   real(r8),allocatable :: crm_cld(:,:,:,:)
   real(r8),allocatable :: clubb_tk(:,:,:,:)
   real(r8),allocatable :: clubb_tkh(:,:,:,:)
   real(r8),allocatable :: relvar(:,:,:,:)
   real(r8),allocatable :: accre_enhan(:,:,:,:)
   real(r8),allocatable :: qclvar(:,:,:,:)

   real(r8) crm_tk(pcols,crm_nx, crm_ny, crm_nz)
   real(r8) crm_tkh(pcols,crm_nx, crm_ny, crm_nz)
   real(r8) cld3d_crm(pcols, crm_nx, crm_ny, crm_nz)   ! 3D instaneous cloud fraction 
   real(r8) prec_crm(pcols,crm_nx, crm_ny)
   real(r8) mctot(pcols,pver)        ! total cloud mass flux
   real(r8) mcup(pcols,pver)         ! cloud updraft mass flux
   real(r8) mcdn(pcols,pver)         ! cloud downdraft mass flux
   real(r8) mcuup(pcols,pver)        ! unsaturated updraft mass flux
   real(r8) mcudn(pcols,pver)        ! unsaturated downdraft mass flux
   real(r8) spqc(pcols,pver)         ! cloud water
   real(r8) spqi(pcols,pver)         ! cloud ice
   real(r8) spqs(pcols,pver)         ! snow
   real(r8) spqg(pcols,pver)         ! graupel
   real(r8) spqr(pcols,pver)         ! rain
   real(r8) spnc(pcols,pver)         ! cloud water droplet (#/kg)
   real(r8) spni(pcols,pver)         ! cloud ice crystal number (#/kg)
   real(r8) spns(pcols,pver)         ! snow particle number (#/kg)
   real(r8) spng(pcols,pver)         ! graupel particle number (#/kg)
   real(r8) spnr(pcols,pver)         ! rain particle number (#/kg)
   real(r8) wvar_crm (pcols,crm_nx, crm_ny, crm_nz)   ! vertical velocity variance (m/s)

   real(r8) aut_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Cloud water autoconversion (1/s)
   real(r8) acc_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Cloud water accretion by rain (1/s)
   real(r8) evpc_crm (pcols,crm_nx, crm_ny, crm_nz)  ! Cloud water evaporation (1/s)
   real(r8) evpr_crm (pcols,crm_nx, crm_ny, crm_nz)  ! Rain evaporation (1/s)
   real(r8) mlt_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Ice, snow, graupel melting (1/s)
   real(r8) sub_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Ice, snow, graupel sublimation (1/s)
   real(r8) dep_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Ice, snow, graupel deposition (1/s)
   real(r8) con_crm (pcols,crm_nx, crm_ny, crm_nz)   ! Cloud water condensation (1/s)
   real(r8) aut_crm_a (pcols,pver)   ! Cloud water autoconversion (1/s)
   real(r8) acc_crm_a (pcols,pver)   ! Cloud water accretion by rain (1/s)
   real(r8) evpc_crm_a (pcols,pver)  ! Cloud water evaporation (1/s)
   real(r8) evpr_crm_a (pcols,pver)  ! Rain evaporation (1/s)
   real(r8) mlt_crm_a (pcols,pver)   ! Ice, snow, graupel melting (1/s)
   real(r8) sub_crm_a (pcols,pver)   ! Ice, snow, graupel sublimation (1/s)
   real(r8) dep_crm_a (pcols,pver)   ! Ice, snow, graupel deposition (1/s)
   real(r8) con_crm_a (pcols,pver)   ! Cloud water condensation (1/s)
 
   real(r8) flux_qt(pcols,pver)      ! nonprecipitating water flux
   real(r8) flux_u(pcols,pver)       ! x-momentum flux
   real(r8) flux_v(pcols,pver)       ! y-momentum flux
   real(r8) fluxsgs_qt(pcols,pver)   ! sgs nonprecipitating water flux
   real(r8) tkez(pcols,pver)         ! tke profile [kg/m/s2]
   real(r8) tkesgsz(pcols,pver)      ! sgs tke profile  [kg/m/s2]
   real(r8) flux_qp(pcols,pver)      ! precipitating water flux
   real(r8) precflux(pcols,pver)     ! precipitation flux
   real(r8) qt_ls(pcols,pver)        ! water tendency due to large-scale
   real(r8) qt_trans(pcols,pver)     ! nonprecip water tendency due to transport
   real(r8) qp_trans(pcols,pver)     ! precip water tendency due to transport
   real(r8) qp_fall(pcols,pver)      ! precip water tendency due to fall-out
   real(r8) qp_evp(pcols,pver)       ! precip water tendency due to evap
   real(r8) qp_src(pcols,pver)       ! precip water tendency due to conversion
   real(r8) t_ls(pcols,pver)         ! tendency of crm's liwse due to large-scale
   real(r8) cldtop(pcols,pver)
   real(r8) cwp   (pcols,pver)       ! in-cloud cloud (total) water path (kg/m2)
   real(r8) gicewp(pcols,pver)       ! grid-box cloud ice water path  (g/m2)
   real(r8) gliqwp(pcols,pver)       ! grid-box cloud liquid water path (g/m2)
   real(r8) gwp   (pcols,pver)       ! grid-box cloud (total) water path (kg/m2)
   real(r8) tgicewp(pcols)           ! Vertically integrated ice water path (kg/m2
   real(r8) tgliqwp(pcols)           ! Vertically integrated liquid water path (kg/m2)
   real(r8) cicewp(pcols,pver)       ! in-cloud cloud ice water path (kg/m2)
   real(r8) cliqwp(pcols,pver)       ! in-cloud cloud liquid water path (kg/m2)
   real(r8) tgwp   (pcols)           ! Vertically integrated (total) cloud water path  (kg/m2)
   real(r8) precc(pcols)             ! convective precipitation [m/s]
   real(r8) precl(pcols)             ! large scale precipitation [m/s]
   real(r8) precsc(pcols)            ! convecitve snow   [m/s]
   real(r8) precsl(pcols)            ! convective snow   [m/s]
   real(r8) cltot(pcols)             ! Diagnostic total cloud cover
   real(r8) cllow(pcols)             ! Diagnostic low cloud cover
   real(r8) clmed(pcols)             ! Diagnostic mid  cloud cover
   real(r8) clhgh(pcols)             ! Diagnostic hgh  cloud cover
   real(r8) :: ftem(pcols,pver)      ! Temporary workspace for outfld variables
   real(r8) ul(pver)
   real(r8) vl(pver)

   real(r8) :: mu_crm(pcols,pver)   
   real(r8) :: md_crm(pcols,pver) 
   real(r8) :: du_crm(pcols,pver)
   real(r8) :: eu_crm(pcols,pver)
   real(r8) :: ed_crm(pcols,pver)
   real(r8) :: tk_crm(pcols,pver)
   real(r8) :: jt_crm(pcols)
   real(r8) :: mx_crm(pcols)
   real(r8) :: ideep_crm(pcols)


   integer itim 
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction

   real(r8),allocatable ::  na(:)              ! aerosol number concentration [/m3]
   real(r8),allocatable ::  va(:)              ! aerosol voume concentration [m3/m3]
   real(r8),allocatable ::  hy(:)              ! aerosol bulk hygroscopicity
   real(r8),allocatable ::  naermod(:,:)       ! Aerosol number concentration [/m3]
   real(r8),allocatable ::  vaerosol(:,:)      ! aerosol volume concentration [m3/m3]
   real(r8),allocatable ::  hygro(:,:)         ! hygroscopicity of aerosol mode 
   integer   phase                             ! phase to determine whether it is interstitial, cloud-borne, or the sum. 

   real(r8) cs(pcols, pver)                     ! air density  [kg/m3]

   real(r8),allocatable :: qicecen(:,:,:,:,:)   ! cloud ice (kg/kg)
   real(r8),allocatable :: qlsink_afcen(:,:,:,:,:)  ! cloud water loss rate from precipitation calculated
                                                                                      ! cloud water before precipitatinog (/s)
   real(r8),allocatable :: qlsink_bfcen(:,:,:,:,:)  ! cloud water loss rate from precipitation calculated
                                                                                      ! cloud water before precipitatinog (/s)
   real(r8),allocatable :: qlsink_avgcen(:,:,:,:,:)  ! cloud water loss rate from precipitation calculated
                                                                                      ! from praincen and qlcoudcen averaged over
                                                                                      ! ntavg1_ss time step (/s)
   real(r8),allocatable :: praincen(:,:,:,:,:)  ! cloud water loss rate from precipitation (kg/kg/s)
   real(r8),allocatable :: wupthresh_bnd(:,:) 
   real(r8),allocatable :: wdownthresh_bnd(:,:)

   ! CRM column radiation stuff:
   real(r8) prectend(pcols) ! tendency in precipitating water and ice
   real(r8) precstend(pcols) ! tendency in precipitating ice
   real(r8) icesink(pcols) ! sink of
   real(r8) tau00  ! surface stress
   real(r8) wnd  ! surface wnd
   real(r8) bflx   ! surface buoyancy flux (Km/s)
   real(r8) :: temp_diff_crm_u,temp_diff_crm_v,temp_diff_crm_w,temp_diff_crm_t,temp_diff_crm_micro1,temp_diff_crm_micro2,temp_diff_crm_micro3
   real(r8) :: temp_diff_crm_qrad,temp_diff_qc_crm,temp_diff_qi_crm,temp_diff_qpc_crm,temp_diff_qpi_crm,temp_diff_t_rad,temp_diff_qv_rad,temp_diff_qi_rad,temp_diff_cld_rad,temp_diff_cld3d_crm
   real(r8) :: out_precc,out_precl,out_precsc,out_precsl,out_cltot,out_clhgh,out_clmed,out_cllow,out_prectend, out_precstend,out_timing_factor,out_global_rank 
   real(r8) :: out_ocnfrac,out_wnd,out_tau00,out_bflx,out_fluxu0,out_fluxv0,out_fluxt0,out_fluxq0,out_taux_crm,out_tauy_crm,out_z0m,out_it,out_jt,out_jt_crm,out_mx_crm
   real(r8), dimension(pver) :: out_qltend,out_qcltend,out_qiltend,out_sltend,out_cld,out_cldtop,out_gicewp,out_gliqwp,out_mctot,out_mcup
   real(r8), dimension(pver) :: out_mcdn,out_mcuup,out_mcudn,out_spqc,out_spqi,out_spqs,out_spqg,out_spqr,out_mu_crm,out_md_crm,out_du_crm
   real(r8), dimension(pver) :: out_eu_crm,out_ed_crm,out_tkez,out_tkesgsz,out_tk_crm,out_flux_u,out_precflux,out_flux_v,out_flux_qt,out_fluxsgs_qt,out_flux_qp,out_qt_ls
   real(r8), dimension(pver) :: out_qt_trans,out_qp_trans,out_qp_fall,out_qp_evp,out_qp_src, out_t_ls
   real(r8), dimension(pver) :: work_crm_u,work_crm_v,work_crm_w,work_crm_t,work_crm_micro1,work_crm_micro2,work_crm_micro3
   real(r8), dimension(pver) :: work_org_crm_u,work_org_crm_v,work_org_crm_w,work_org_crm_t,work_org_crm_micro1,work_org_crm_micro2,work_org_crm_micro3
   real(r8), dimension(pver) :: work_crm_qrad,work_qc_crm,work_qi_crm,work_qpc_crm,work_qpi_crm,work_t_rad
   real(r8), dimension(pver) :: work_org_crm_qrad,work_org_qc_crm,work_org_qi_crm,work_org_qpc_crm,work_org_qpi_crm,work_org_t_rad
   real(r8), dimension(pver) :: work_qv_rad,work_qc_rad,work_qi_rad,work_cld_rad,work_cld3d_crm
   real(r8), dimension(pver) :: work_org_qv_rad,work_org_qc_rad,work_org_qi_rad,work_org_cld_rad,work_org_cld3d_crm
   real(r8), dimension(pver) :: diff_crm_u,maxdiff_crm_u,sqrt_crm_u,diff_crm_v,maxdiff_crm_v,sqrt_crm_v,diff_crm_w,maxdiff_crm_w,sqrt_crm_w,diff_crm_t,maxdiff_crm_t,sqrt_crm_t
   real(r8), dimension(pver) :: diff_crm_micro3,maxdiff_crm_micro3,sqrt_crm_micro3,diff_crm_micro2,maxdiff_crm_micro2,sqrt_crm_micro2,diff_crm_micro1,maxdiff_crm_micro1,sqrt_crm_micro1
   real(r8), dimension(pver) :: diff_crm_qrad,maxdiff_crm_qrad,sqrt_crm_qrad,diff_qc_rad,maxdiff_qc_rad,sqrt_qc_rad,diff_qi_crm,maxdiff_qi_crm,sqrt_qi_crm,diff_qc_crm,maxdiff_qc_crm,sqrt_qc_crm
   real(r8), dimension(pver) :: diff_qpc_crm,maxdiff_qpc_crm,sqrt_qpc_crm,diff_qpi_crm,maxdiff_qpi_crm,sqrt_qpi_crm,diff_t_rad,maxdiff_t_rad,sqrt_t_rad,diff_qv_rad,maxdiff_qv_rad,sqrt_qv_rad
   real(r8), dimension(pver) :: diff_qi_rad,maxdiff_qi_rad,sqrt_qi_rad,diff_cld_rad,maxdiff_cld_rad,sqrt_cld_rad,diff_cld3d_crm,maxdiff_cld3d_crm,sqrt_cld3d_crm

   real(r8), dimension(20)   :: out_qtotcrm
   real(r8), dimension(crm_nx, crm_ny, crm_nz) :: out_crm_u,out_crm_v,out_crm_w,out_crm_t,out_crm_qrad,out_qc_crm,out_qi_crm,out_qpc_crm
   real(r8), dimension(crm_nx, crm_ny, crm_nz) :: out_qpi_crm,out_t_rad,out_qv_rad,out_qc_rad,out_qi_rad,out_cld_rad,out_cld3d_crm,out_crm_tk,out_crm_tkh
   real(r8), dimension(crm_nx, crm_ny, crm_nz,nmicro_fields_total) ::  out_crm_micro
   real(r8), dimension(crm_nx, crm_ny) :: out_prec_crm
   real(r8), dimension(pcols, crm_nx, crm_ny, crm_nz) :: orc_crm_u,orc_crm_v,orc_crm_w,orc_crm_t,orc_crm_qrad,orc_qc_crm,orc_qi_crm,orc_qpc_crm
   real(r8), dimension(pcols, crm_nx, crm_ny, crm_nz) :: orc_qpi_crm,orc_t_rad,orc_qv_rad,orc_qc_rad,orc_qi_rad,orc_cld_rad,orc_cld3d_crm
   real(r8), dimension(pcols, crm_nx, crm_ny, crm_nz,nmicro_fields_total) :: orc_crm_micro 
   real(r8), dimension(pcols, crm_nx, crm_ny) :: orc_prec_crm
   real(r8) taux_crm(pcols)  ! zonal CRM surface stress perturbation
   real(r8) tauy_crm(pcols)  ! merid CRM surface stress perturbation
   real(r8) z0m(pcols)  ! surface momentum roughness length
   real(r8), pointer, dimension(:,:) :: qrs, qrl        ! rad heating rates
   real(r8), pointer, dimension(:,:,:,:)   :: crm_u
   real(r8), pointer, dimension(:,:,:,:)   :: crm_v
   real(r8), pointer, dimension(:,:,:,:)   :: crm_w
   real(r8), pointer, dimension(:,:,:,:)   :: crm_t
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qt
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qp
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qn
   real(r8), pointer, dimension(:,:,:,:)   :: crm_nc
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qr
   real(r8), pointer, dimension(:,:,:,:)   :: crm_nr
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qi
   real(r8), pointer, dimension(:,:,:,:)   :: crm_ni
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qs
   real(r8), pointer, dimension(:,:,:,:)   :: crm_ns
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qg
   real(r8), pointer, dimension(:,:,:,:)   :: crm_ng
   real(r8), pointer, dimension(:,:,:,:)   :: crm_qc

   real(r8), allocatable, dimension(:,:,:,:,:) :: crm_micro

   integer                         :: pblh_idx
   real(r8), pointer, dimension(:) :: pblh
   
   real(r8), pointer, dimension(:,:) :: wsedl

   real(r8),allocatable :: acen(:,:,:,:,:)   ! cloud fraction for each sub-sub class for full time period
   real(r8),allocatable :: acen_tf(:,:,:,:,:) ! cloud fraction for end-portion of time period
   real(r8),allocatable :: rhcen(:,:,:,:,:)  ! relative humidity (0-1)
   real(r8),allocatable :: qcloudcen(:,:,:,:,:)  ! cloud water (kg/kg)
   real(r8),allocatable :: qlsinkcen(:,:,:,:,:)  ! cloud water loss rate from precipitation (/s??)
   real(r8),allocatable :: precrcen(:,:,:,:,:)   ! liquid (rain) precipitation rate (kg/m2/s)
   real(r8),allocatable :: precsolidcen(:,:,:,:,:)   ! solid (rain) precipitation rate (kg/m2/s)
   real(r8),allocatable :: wwqui_cen(:,:)            ! vertical velocity variance in quiescent class (m2/s2)
   real(r8),allocatable :: wwqui_cloudy_cen(:,:)     ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
   ! at layer boundary
   real(r8),allocatable :: abnd(:,:,:,:,:)       ! cloud fraction for each sub-sub class for full time period
   real(r8),allocatable :: abnd_tf(:,:,:,:,:)    ! cloud fraction for end-portion of time period
   real(r8),allocatable :: massflxbnd(:,:,:,:,:) ! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
   real(r8),allocatable :: wwqui_bnd(:,:)        ! vertical velocity variance in quiescent class (m2/s2)
   real(r8),allocatable :: wwqui_cloudy_bnd(:,:) ! vertical velocity variance in quiescent, and cloudy class (m2/s2)

   integer,  pointer :: nmxrgn(:)      ! Number of maximally overlapped regions
   real(r8), pointer :: pmxrgn(:,:)    ! Maximum values of pressure for each

   real(r8), allocatable :: spup2(:,:) 
   real(r8), allocatable :: spvp2(:,:)
   real(r8), allocatable :: spwprtp(:,:)
   real(r8), allocatable :: spwpthlp(:,:)
   real(r8), allocatable :: spwp2(:,:)
   real(r8), allocatable :: spwp3(:,:)
   real(r8), allocatable :: sprtp2(:,:)
   real(r8), allocatable :: spthlp2(:,:)
   real(r8), allocatable :: sprtpthlp(:,:)
   real(r8), allocatable :: spupwp(:,:)
   real(r8), allocatable :: spvpwp(:,:)
   real(r8), allocatable :: spcrm_cld(:,:)
   real(r8), allocatable :: spt_tndcy(:,:)
   real(r8), allocatable :: spqv_tndcy(:,:)
   real(r8), allocatable :: spqc_tndcy(:,:)
   real(r8), allocatable :: spclubb_tk(:,:)
   real(r8), allocatable :: spclubb_tkh(:,:)
   real(r8), allocatable :: sprelvar(:,:)
   real(r8), allocatable :: spaccre_enhan(:,:)
   real(r8), allocatable :: spqclvar(:,:)

   real(r8)              :: spcld3d (pcols,pver)
   real(r8) :: work_u(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_v(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_w(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_tabs(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qcl(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qci(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qpl(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qpi(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qv(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_q(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qn(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_t(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_qp(crm_nx, crm_ny, crm_nz)
   real(r8) :: work_u0(crm_nz)
   real(r8) :: work_v0(crm_nz)
   real(r8) :: work_t0(crm_nz)
   real(r8) :: work_t00(crm_nz)
   real(r8) :: work_tabs0(crm_nz)
   real(r8) :: work_q0(crm_nz)
   real(r8) :: work_qv0(crm_nz)
   real(r8) :: work_qn0(crm_nz)
   real(r8) :: work_qp0(crm_nz)
   real(r8) :: work_tke0(crm_nz)
   real(r8) :: work_gamaz(crm_nz)
   real(r8) :: work_omn,work_omp

   real(r8) :: tmp4d(pcols,crm_nx, crm_ny, crm_nz)
   real(r8) :: tmp2d(pcols,pver)

   ! Surface fluxes
   real(r8) ::  fluxu0           ! surface momenment fluxes
   real(r8) ::  fluxv0           ! surface momenment fluxes
   real(r8) ::  fluxt0           ! surface sensible heat fluxes
   real(r8) ::  fluxq0           ! surface latent heat fluxes
   real(r8) ::  dtstep_pp        ! time step for the ECPP (seconds)
   integer  ::  necpp            ! the number of GCM time step in which ECPP is called once.
   real(r8) ::  longitude0
   real(r8) ::  latitude0
   real(r8) :: cldd(pver)
   real(r8) radflux(pcols)   ! radiative fluxes from radiation calculation (qrs + qrl)
   CHARACTER(LEN=6) :: crm_number
   real(r8) qtot(pcols, 3)      ! total water
   real(r8) qt_hydro(pcols, 2)  ! total hydrometer 
   real(r8) qt_cloud(pcols, 3)  ! total cloud water
   real(r8) qtv(pcols, 3)       ! total water vapor
   real(r8) qli_hydro(pcols, 2) ! column-integraetd rain + snow + graupel 
   real(r8) qi_hydro(pcols, 2)  ! column-integrated snow water + graupel water
   real(r8) sfactor 

   real(r8) zero(pcols)          ! zero
   real(r8) timing_factor(pcols) ! factor for crm cpu-usage: 1 means no subcycling

   real(r8) qtotcrm(pcols, 20)   ! the toal water calculated in crm.F90

   real(r8), parameter :: rhow = 1000._r8
   real(r8), parameter :: bc   = 2._r8
   real(r8) :: t, mu, acn, dumc, dunc, pgam, lamc
   real(r8) :: dunc_arr(pcols,pver)
   integer gcolindex(pcols)  ! array of global latitude indices
   integer ii, jj,status
   integer iii,FlagEnd,it,jt
   integer i, k, m,lchnk_save,i_save
   integer ifld
   logical :: ls, lu, lv, lq(pcnst)
   integer,parameter :: structlen = 49
   integer,parameter :: singlelen = 30
   integer,parameter :: flen      = structlen*pver+singlelen+1+20+pcols
   integer,parameter :: flen2     = 17*orc_nx*orc_ny*crm_nz + orc_nx*orc_ny*crm_nz*nmicro_fields_total+orc_nx*orc_ny
   integer,dimension(structlen) :: blocklengths
   INTEGER,dimension(structlen) :: types
   integer :: crm_comm
   integer :: myrank_global, numproc_global
   integer :: myrank_crm, numproc_crm
#ifdef ORCHESTRATOR
   integer ::dest,ierr,ierr1,ierr2,ierr3,fcount,chnksz,kk,ll,iorc,crm_step,crm_start_ind,crm_end_ind,icheck_column,icheck_iorc
   real(r8), allocatable :: flattened_crm_inout(:)
   real(r8),dimension(flen) :: Var_Flat
   real(r8),dimension(flen2) :: Var_Flat2
   integer,parameter :: structleno = 39
   integer,parameter :: singleleno = 27
   integer,parameter :: fleno      = structleno*pver+singleleno+1+20
   integer,parameter :: nflat =  17*(orc_nx*orc_ny*crm_nz) + orc_nx*orc_ny*crm_nz*nmicro_fields_total+orc_nx*orc_ny
   real(r8),dimension(fleno+nflat) :: CRM_Var_Flat
   double precision :: wall(6), sys(6), usr(6)
#endif
   zero = 0.0_r8

!========================================================
!========================================================
!  CRM (Superparameterization).
! Author: Marat Khairoutdinov (mkhairoutdin@ms.cc.sunysb.edu)
!========================================================
   call t_stampf(wall(1), usr(1), sys(1))
   call t_startf ('crm')

   allocate(crm_micro(pcols,crm_nx,crm_ny,crm_nz,nmicro_fields_total))

   ! Initialize stuff:
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   ls           = .TRUE.
   lq(:)        = .FALSE.
   lq(1)        = .TRUE.
   lq(ixcldliq) = .TRUE.
   lq(ixcldice) = .TRUE.
   lu           = .FALSE.
   lv           = .FALSE.
   call physics_ptend_init(ptend,     state%psetcols, 'crm', lu=lu, lv=lv, ls=ls, lq=lq)  ! Initialize output physics_ptend object
   call physics_ptend_init(ptend_loc, state%psetcols, 'crm', lu=lu, lv=lv, ls=ls, lq=lq)  ! Initialize local physics_ptend object

   nstep = get_nstep()

   lchnk = state%lchnk
   ncol  = state%ncol
   lchnk_save = lchnk
   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cld_idx,  cld,  start=(/1,1,itim/), kount=(/pcols,pver,1/))

   call physics_state_copy(state, state_loc)
   tend_loc  = tend
 
   !------------------------- 
   ! Set up general fields
   call pbuf_get_field (pbuf, crm_u_idx,      crm_u)
   call pbuf_get_field (pbuf, crm_v_idx,      crm_v)
   call pbuf_get_field (pbuf, crm_w_idx,      crm_w)
   call pbuf_get_field (pbuf, crm_t_idx,      crm_t)
   call pbuf_get_field (pbuf, crm_qrad_idx,   crm_qrad)
   call pbuf_get_field (pbuf, crm_t_rad_idx,  t_rad)
   call pbuf_get_field (pbuf, crm_qv_rad_idx, qv_rad)
   call pbuf_get_field (pbuf, crm_qc_rad_idx, qc_rad)
   call pbuf_get_field (pbuf, crm_qi_rad_idx, qi_rad)
   call pbuf_get_field (pbuf, crm_cld_rad_idx, cld_rad)

   call pbuf_get_field (pbuf, prec_dp_idx,   prec_dp)
   call pbuf_get_field (pbuf, snow_dp_idx,   snow_dp)


   !------------------------- 
   ! setup CLUBB fields
   if (do_clubb_sgs) then
      allocate(nmxrgn        (pcols))
      allocate(pmxrgn        (pcols,pverp))
      allocate(spup2         (pcols, pver))
      allocate(spvp2         (pcols, pver))
      allocate(spwprtp       (pcols, pver))
      allocate(spwpthlp      (pcols, pver))
      allocate(spwp2         (pcols, pver))
      allocate(spwp3         (pcols, pver))
      allocate(sprtp2        (pcols, pver))
      allocate(spthlp2       (pcols, pver))
      allocate(sprtpthlp     (pcols, pver))
      allocate(spupwp        (pcols, pver))
      allocate(spvpwp        (pcols, pver))
      allocate(spcrm_cld     (pcols, pver))
      allocate(spt_tndcy     (pcols, pver))
      allocate(spqv_tndcy    (pcols, pver))
      allocate(spqc_tndcy    (pcols, pver))
      allocate(spclubb_tk    (pcols, pver))
      allocate(spclubb_tkh   (pcols, pver))
      allocate(sprelvar      (pcols, pver))
      allocate(spaccre_enhan (pcols, pver))
      allocate(spqclvar      (pcols, pver))
      allocate(crm_cld       (pcols,crm_nx, crm_ny, crm_nz+1))
      allocate(clubb_tk      (pcols,crm_nx, crm_ny, crm_nz))
      allocate(clubb_tkh     (pcols,crm_nx, crm_ny, crm_nz))
      allocate(relvar        (pcols,crm_nx, crm_ny, crm_nz))
      allocate(accre_enhan   (pcols,crm_nx, crm_ny, crm_nz))
      allocate(qclvar        (pcols,crm_nx, crm_ny, crm_nz))

      call pbuf_get_field (pbuf, clubb_buffer_idx,  clubb_buffer)

   endif

   !------------------------- 
   ! Setup m2005 fields
   if (is_spcam_m2005) then
     allocate(na (pcols))
     allocate(va (pcols))
     allocate(hy (pcols))
     allocate(naermod (pver, nmodes))
     allocate(vaerosol (pver, nmodes))
     allocate(hygro (pver, nmodes))

     call pbuf_get_field(pbuf, crm_nc_rad_idx, nc_rad)
     call pbuf_get_field(pbuf, crm_ni_rad_idx, ni_rad)
     call pbuf_get_field(pbuf, crm_qs_rad_idx, qs_rad)
     call pbuf_get_field(pbuf, crm_ns_rad_idx, ns_rad)
     call pbuf_get_field(pbuf, crm_qt_idx,     crm_qt)
     call pbuf_get_field(pbuf, crm_nc_idx,     crm_nc)
     call pbuf_get_field(pbuf, crm_qr_idx,     crm_qr)
     call pbuf_get_field(pbuf, crm_nr_idx,     crm_nr)
     call pbuf_get_field(pbuf, crm_qi_idx,     crm_qi)
     call pbuf_get_field(pbuf, crm_ni_idx,     crm_ni)
     call pbuf_get_field(pbuf, crm_qs_idx,     crm_qs)
     call pbuf_get_field(pbuf, crm_ns_idx,     crm_ns)
     call pbuf_get_field(pbuf, crm_qg_idx,     crm_qg)
     call pbuf_get_field(pbuf, crm_ng_idx,     crm_ng)
     call pbuf_get_field(pbuf, crm_qc_idx,     crm_qc)

   !------------------------- 
   ! Setup sam1mom fields
   else if (is_spcam_sam1mom) then
     call pbuf_get_field(pbuf, crm_qt_idx,     crm_qt)
     call pbuf_get_field(pbuf, crm_qp_idx,     crm_qp)
     call pbuf_get_field(pbuf, crm_qn_idx,     crm_qn)
   endif


   !------------------------- 
   ! Setup ECPP fields
   ! ifdef needed because of use of NCLASS_CL
#ifdef m2005 
   if (is_spcam_m2005) then
      allocate(acen (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(acen_tf (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(rhcen (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(qcloudcen (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(qlsinkcen (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(precrcen (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(precsolidcen (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(wwqui_cen (pcols, pver))
      allocate(wwqui_cloudy_cen (pcols, pver))

      ! at layer boundary
      allocate(abnd (pcols,pver+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(abnd_tf (pcols,pver+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(massflxbnd (pcols,pver+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(wwqui_bnd (pcols, pver+1))
      allocate(wwqui_cloudy_bnd (pcols, pver+1))

      allocate(qicecen         (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(qlsink_afcen    (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(qlsink_bfcen    (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(qlsink_avgcen   (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(praincen        (pcols,pver,NCLASS_CL,ncls_ecpp_in,NCLASS_PR))
      allocate(wupthresh_bnd   (pcols, pverp)) 
      allocate(wdownthresh_bnd (pcols, pverp))

      call pbuf_get_field(pbuf, tk_crm_idx,    tk_crm_ecpp)
      call pbuf_get_field(pbuf, acldy_cen_idx, acldy_cen_tbeg)

      if(is_first_step())then
          acldy_cen_tbeg(:ncol,:) = cld(:ncol, :)
      end if

   end if
#endif

   !------------------------- 
   ! Initialize all aerosol and gas species
   ! When ECPP is used, dropmixnuc and all transport(deep and shallow) are done in ECPP.
   if (is_spcam_sam1mom) then
         state_loc%q(:ncol, :pver, :pcnst)   = 1.e-36_r8
         ! set the values which SPCAM uses back to state
         state_loc%q(:ncol, :pver, 1)        = state%q(:ncol, :pver, 1)
         state_loc%q(:ncol, :pver, ixcldice) = state%q(:ncol, :pver, ixcldice)
         state_loc%q(:ncol, :pver, ixcldliq) = state%q(:ncol, :pver, ixcldliq)
   endif

call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)

   !------------------------- 
   !------------------------- 
   ! On the first_step, initialize values only and do not call CRM
   !------------------------- 
   !------------------------- 
   if(is_first_step()) then
       do k=1,crm_nz
          m = pver-k+1
          do i=1,ncol

             if (spcam_direction == 'NS') then
                if(crm_ny.eq.1) then ! change domain orientation only for 2D CRM
                  crm_u(i,:,:,k) = state_loc%v(i,m)
                  crm_v(i,:,:,k) = state_loc%u(i,m)
                else
                  crm_u(i,:,:,k) = state_loc%u(i,m)
                  crm_v(i,:,:,k) = state_loc%v(i,m)
                end if
             else if( spcam_direction == 'WE') then
                crm_u(i,:,:,k) = state_loc%u(i,m)
                crm_v(i,:,:,k) = state_loc%v(i,m)
             endif
             crm_w(i,:,:,k) = 0._r8
             !write (13,*),'MPI init w',myrank_global,i,k,crm_w(i,:,:,k)
             crm_t(i,:,:,k) = state_loc%t(i,m)

             if (is_spcam_sam1mom)  then
                crm_qt(i,:,:,k) = state_loc%q(i,m,1)+state_loc%q(i,m,ixcldliq)+state_loc%q(i,m,ixcldice)
                crm_qp(i,:,:,k) = 0.0_r8
                crm_qn(i,:,:,k) = state_loc%q(i,m,ixcldliq)+state_loc%q(i,m,ixcldice)

            else if (is_spcam_m2005) then
                crm_qt(i,:,:,k) = state_loc%q(i,m,1)+state_loc%q(i,m,ixcldliq)
                crm_nc(i,:,:,k) = 0.0_r8
                crm_qr(i,:,:,k) = 0.0_r8
                crm_nr(i,:,:,k) = 0.0_r8
                crm_qi(i,:,:,k) = state_loc%q(i,m,ixcldice) 
                crm_ni(i,:,:,k) = 0.0_r8
                crm_qs(i,:,:,k) = 0.0_r8
                crm_ns(i,:,:,k) = 0.0_r8
                crm_qg(i,:,:,k) = 0.0_r8
                crm_ng(i,:,:,k) = 0.0_r8
                crm_qc(i,:,:,k) = state_loc%q(i,m,ixcldliq)


                nc_rad(i,:,:,k)   = 0._r8
                ni_rad(i,:,:,k)   = 0._r8
                qs_rad(i,:,:,k)   = 0.0_r8
                ns_rad(i,:,:,k)   = 0.0_r8
                wvar_crm(i,:,:,k) = 0.0_r8
                aut_crm(i,:,:,k)  = 0.0_r8
                acc_crm(i,:,:,k)  = 0.0_r8
                evpc_crm(i,:,:,k) = 0.0_r8
                evpr_crm(i,:,:,k) = 0.0_r8
                mlt_crm(i,:,:,k)  = 0.0_r8
                sub_crm(i,:,:,k)  = 0.0_r8
                dep_crm(i,:,:,k)  = 0.0_r8
                con_crm(i,:,:,k)  = 0.0_r8
             endif

             if (do_clubb_sgs) then
                ! In the inital run, variables are set in clubb_sgs_setup at the first time step
                clubb_buffer(i,:,:,k,:) = 0.0_r8  
             endif

             crm_qrad (i,:,:,k)    = 0._r8
             qc_crm (i,:,:,k)      = 0._r8
             qi_crm (i,:,:,k)      = 0._r8
             qpc_crm(i,:,:,k)      = 0._r8
             qpi_crm(i,:,:,k)      = 0._r8
             t_rad  (i,:,:,k)      = state_loc%t(i,m)
             qv_rad (i,:,:,k)      = state_loc%q(i,m,1)
             qc_rad (i,:,:,k)      = 0._r8
             qi_rad (i,:,:,k)      = 0._r8
             cld_rad(i,:,:,k)      = 0._r8
         end do
       end do
       !do i=1,ncol
       !  write(iulog,*) 'Check Init u: ',i,crm_u(i,:,:,:)
       !end do
       ! use radiation from grid-cell mean radctl on first time step
       prec_crm (:,:,:)          = 0._r8
       ptend_loc%q(:,:,1)        = 0._r8
       ptend_loc%q(:,:,ixcldliq) = 0._r8
       ptend_loc%q(:,:,ixcldice) = 0._r8
       ptend_loc%s(:,:)          = 0._r8
       precc(:)                  = 0._r8
       precl(:)                  = 0._r8
       precsc(:)                 = 0._r8
       precsl(:)                 = 0._r8
       cltot(:)                  = 0._r8
       clhgh(:)                  = 0._r8
       clmed(:)                  = 0._r8
       cllow(:)                  = 0._r8
       cld(:,:)                  = 0._r8
       cldtop(:,:)               = 0._r8
       gicewp(:,:)               = 0._r8
       gliqwp(:,:)               = 0._r8
       mctot(:,:)                = 0._r8
       mcup(:,:)                 = 0._r8
       mcdn(:,:)                 = 0._r8
       mcuup(:,:)                = 0._r8
       mcudn(:,:)                = 0._r8
       spqc(:,:)                 = 0._r8
       spqi(:,:)                 = 0._r8
       spqs(:,:)                 = 0._r8
       spqg(:,:)                 = 0._r8
       spqr(:,:)                 = 0._r8
       cld3d_crm (:,:,:,:)       = 0._r8
       flux_qt(:,:)              = 0._r8
       flux_u(:,:)               = 0._r8
       flux_v(:,:)               = 0._r8
       fluxsgs_qt(:,:)           = 0._r8
       tkez(:,:)                 = 0._r8
       tkesgsz(:,:)              = 0._r8
       flux_qp(:,:)              = 0._r8
       precflux(:,:)             = 0._r8
       qt_ls(:,:)                = 0._r8
       qt_trans(:,:)             = 0._r8
       qp_trans(:,:)             = 0._r8
       qp_fall(:,:)              = 0._r8
       qp_evp(:,:)               = 0._r8
       qp_src(:,:)               = 0._r8
       z0m(:)                    = 0._r8
       taux_crm(:)               = 0._r8
       tauy_crm(:)               = 0._r8
       t_ls(:,:)                 = 0._r8


       if (is_spcam_m2005) then
         spnc(:,:)          = 0._r8
         spni(:,:)          = 0._r8
         spns(:,:)          = 0._r8
         spng(:,:)          = 0._r8
         spnr(:,:)          = 0._r8
         aut_crm_a(:,:)     = 0._r8
         acc_crm_a(:,:)     = 0._r8
         evpc_crm_a(:,:)    = 0._r8
         evpr_crm_a(:,:)    = 0._r8
         mlt_crm_a(:,:)     = 0._r8
         sub_crm_a(:,:)     = 0._r8
         dep_crm_a(:,:)     = 0._r8
         con_crm_a(:,:)     = 0._r8
         abnd               = 0.0_r8
         abnd_tf            = 0.0_r8
         massflxbnd         = 0.0_r8
         acen               = 0.0_r8
         acen_tf            = 0.0_r8
         rhcen              = 0.0_r8
         qcloudcen          = 0.0_r8
         qicecen            = 0.0_r8
         qlsinkcen          = 0.0_r8
         precrcen           = 0.0_r8
         precsolidcen       = 0.0_r8
         wupthresh_bnd      = 0.0_r8
         wdownthresh_bnd    = 0.0_r8
         qlsink_afcen       = 0.0_r8
         qlsink_bfcen       = 0.0_r8
         qlsink_avgcen      = 0.0_r8
         praincen           = 0.0_r8

         ! default is clear, non-precipitating, and quiescent class
         abnd(:,:,1,1,1)    = 1.0_r8
         abnd_tf(:,:,1,1,1) = 1.0_r8
         acen(:,:,1,1,1)    = 1.0_r8
         acen_tf(:,:,1,1,1) = 1.0_r8
         wwqui_cen          = 0.0_r8
         wwqui_bnd          = 0.0_r8
         wwqui_cloudy_cen   = 0.0_r8
         wwqui_cloudy_bnd   = 0.0_r8
         tk_crm             = 0.0_r8

         ! turbulence
         cs(:ncol, 1:pver) = state_loc%pmid(:ncol, 1:pver)/(287.15_r8*state_loc%t(:ncol, 1:pver))

       endif 
       FlagEnd = 0
   !------------------------- 
   !------------------------- 
   ! not is_first_step
   !------------------------- 
   !------------------------- 

   else

       call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)

       !do i=1,ncol
       !  write(iulog,*) 'Check step2 u: ',myrank_global,i,crm_u(i,:,:,:)
       !end do


       !write (iulog,*),'MPI init w1',myrank_global,crm_w(1,:,:,:)
       ptend_loc%q(:,:,1)        = 0._r8
       ptend_loc%q(:,:,ixcldliq) = 0._r8
       ptend_loc%q(:,:,ixcldice) = 0._r8
       ptend_loc%s(:,:)          = 0._r8
       cwp      = 0._r8
       gicewp   = 0._r8
       gliqwp   = 0._r8
       cltot    = 0._r8
       clhgh    = 0._r8
       clmed    = 0._r8
       cllow    = 0._r8

       qc_crm   = 0._r8
       qi_crm   = 0._r8
       qpc_crm  = 0._r8
       qpi_crm  = 0._r8
       prec_crm = 0._r8

       ! Populate the internal crm_micro array
       if (is_spcam_sam1mom) then
          crm_micro(:,:,:,:,1)  = crm_qt(:,:,:,:)
          crm_micro(:,:,:,:,2)  = crm_qp(:,:,:,:)
          crm_micro(:,:,:,:,3)  = crm_qn(:,:,:,:)
       else if (is_spcam_m2005) then
          crm_micro(:,:,:,:,1)  = crm_qt(:,:,:,:)
          crm_micro(:,:,:,:,2)  = crm_nc(:,:,:,:)
          crm_micro(:,:,:,:,3)  = crm_qr(:,:,:,:)
          crm_micro(:,:,:,:,4)  = crm_nr(:,:,:,:)
          crm_micro(:,:,:,:,5)  = crm_qi(:,:,:,:)
          crm_micro(:,:,:,:,6)  = crm_ni(:,:,:,:)
          crm_micro(:,:,:,:,7)  = crm_qs(:,:,:,:)
          crm_micro(:,:,:,:,8)  = crm_ns(:,:,:,:)
          crm_micro(:,:,:,:,9)  = crm_qg(:,:,:,:)
          crm_micro(:,:,:,:,10) = crm_ng(:,:,:,:)
          crm_micro(:,:,:,:,11) = crm_qc(:,:,:,:)

          ! initialize gcm-time-step-avg output at start of each time step
          aut_crm_a  = 0.0_r8
          acc_crm_a  = 0.0_r8
          evpc_crm_a = 0.0_r8
          evpr_crm_a = 0.0_r8
          mlt_crm_a  = 0.0_r8
          sub_crm_a  = 0.0_r8
          dep_crm_a  = 0.0_r8
          con_crm_a  = 0.0_r8
       endif

       call t_startf ('crm_call')

       do m=1,crm_nz
         k = pver-m+1
         do i = 1,ncol
             crm_qrad(i,:,:,m) = crm_qrad(i,:,:,m) / state_loc%pdel(i,k) ! for energy conservation
         end do
       end do

      if (is_spcam_m2005) then
          cs(1:ncol, 1:pver) = state_loc%pmid(1:ncol, 1:pver)/(287.15_r8*state_loc%t(1:ncol, 1:pver))
      end if
      do i = 1,ncol
         tau00  = sqrt(cam_in%wsx(i)**2 + cam_in%wsy(i)**2)
         wnd    = sqrt(state_loc%u(i,pver)**2 + state_loc%v(i,pver)**2)
         bflx   = cam_in%shf(i)/cpair + 0.61_r8*state_loc%t(i,pver)*cam_in%lhf(i)/latvap
         fluxu0 = cam_in%wsx(i)     !N/m2
         fluxv0 = cam_in%wsy(i)     !N/m2
         fluxt0 = cam_in%shf(i)/cpair  ! K Kg/ (m2 s)
         fluxq0 = cam_in%lhf(i)/latvap ! Kg/(m2 s)

         !
         ! calculate total water before calling crm
         ! total hydrometer water (rain, snow, and graupel)
         if (is_spcam_m2005) then
            qt_hydro(i, 1)  = 0.0_r8
            qli_hydro(i, 1) = 0.0_r8
            qi_hydro(i, 1)  = 0.0_r8
            do m=1, crm_nz
              k=pver-m+1
              do ii=1, crm_nx
              do jj=1, crm_ny
                qt_hydro(i,1)  = qt_hydro(i,1)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit 
                qli_hydro(i,1) = qli_hydro(i,1)+(crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit
                qi_hydro(i,1)  = qi_hydro(i,1)+(crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * state_loc%pdel(i,k)/gravit
              end do
              end do 
            end do
            qt_hydro(i,1)  = qt_hydro(i,1)  / (crm_nx_ny)
            qli_hydro(i,1) = qli_hydro(i,1) / (crm_nx_ny)
            qi_hydro(i,1)  = qi_hydro(i,1)  / (crm_nx_ny)

            ! total cloud water and total water vapor
            qt_cloud(i,1) = 0._r8
            qtv(i,1)      = 0._r8
            do k=1, pver
               qt_cloud(i,1) = qt_cloud(i,1) + (state_loc%q(i,k,ixcldliq)+state_loc%q(i,k,ixcldice)) * state_loc%pdel(i,k)/gravit
               qtv(i,1)      = qtv(i,1) + state_loc%q(i,k,1) * state_loc%pdel(i,k)/gravit
            end do 

            ! total water 
            qtot(i,1) = qt_hydro(i,1) + qt_cloud(i,1) + qtv(i,1)

         else if (is_spcam_sam1mom) then
            qli_hydro(i, 1) = 0.0_r8
            qi_hydro(i, 1)  = 0.0_r8
            do m=1, crm_nz
              k=pver-m+1
              do ii=1, crm_nx
              do jj=1, crm_ny
                sfactor        = max(0._r8,min(1._r8,(crm_t(i,ii,jj,m)-268.16_r8)*1._r8/(283.16_r8-268.16_r8)))
                qli_hydro(i,1) = qli_hydro(i,1)+crm_qp(i,ii,jj,m) * state_loc%pdel(i,k)/gravit
                qi_hydro(i,1)  = qi_hydro(i,1)+crm_qp(i,ii,jj,m) * (1-sfactor) * state_loc%pdel(i,k)/gravit
              end do
              end do
            end do
            qli_hydro(i,1) = qli_hydro(i,1) / (crm_nx_ny)
            qi_hydro(i,1)  = qi_hydro(i,1)  / (crm_nx_ny)    

            ! total cloud water and total water vapor, and energy
            qt_cloud(i,1) = 0._r8
            qtv(i,1)      = 0._r8
            do k=1, pver
               qt_cloud(i,1) = qt_cloud(i,1) + (state_loc%q(i,k,ixcldliq)+state_loc%q(i,k,ixcldice)) * state_loc%pdel(i,k)/gravit
               qtv(i,1)      = qtv(i,1) + state_loc%q(i,k,1) * state_loc%pdel(i,k)/gravit
            end do
         endif

! ifdef required because of loadaer
#ifdef MODAL_AERO
         if (prog_modal_aero) then
            do k=1, pver
              phase = 1  ! interstital aerosols only
              do m=1, nmodes
               call loadaer( &
                     state_loc, pbuf, i, i, k, &
                     m, cs, phase, na, va, &
                     hy)
                naermod(k, m)  = na(i)
                vaerosol(k, m) = va(i)
                hygro(k, m)    = hy(i)
              end do    
            end do
         endif
#endif

         if (spcam_direction == 'NS') then
            if(crm_ny.eq.1) then
              ul(:) = state_loc%v(i,:)  ! change orientation only if 2D CRM
              vl(:) = state_loc%u(i,:)
            else
              ul(:) = state_loc%u(i,:)
              vl(:) = state_loc%v(i,:)
            end if
         else if (spcam_direction == 'WE') then
              ul(:) = state_loc%u(i,:)
              vl(:) = state_loc%v(i,:)
         endif
       ! Outsource the call to CRM to the separate Orchestrator executable by
       ! sending required inputs...
       end do !do i = 1,ncol
#ifdef ORCHESTRATOR
       ! =================================================
       !
       !
       ! =================================================
       call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
       do i = 1,ncol
         if (state%isorchestrated(i)) then
          do iorc=1,orc_nsubdomains
         !if (state%crmrank0(i) .ne. -2) then ! is this column coupled to an external CRM?
          ! This part is to test the idear of creating a new data type for MPI
          ! transfer data. 
          ! =================================================
          ! Start to flatten the arry and send it to CRM. 
          ! =================================================
          dest = state%crmrank(i,iorc)
          !i_save  = i+state%crmrank1(i)
          i_save  = i
          if (mod(crm_nx,orc_nsubdomains).eq.0) then
            crm_step = int(crm_nx/orc_nsubdomains)
          else
            write(iulog,*) 'WARNING!! crm_nx is not dividable by orc_nsubdomains'
          end if 
          crm_start_ind = (iorc-1)*crm_step+1
          crm_end_ind   = (iorc)*crm_step-1+1
          call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
          chnksz = orc_nx*orc_ny*crm_nz
          call get_gcol_all_p(lchnk, pcols, gcolindex)
          latitude0 = get_rlat_p(lchnk, i_save)*57.296_r8
          longitude0 = get_rlon_p(lchnk, i_save)*57.296_r8 
          Var_Flat(        1)                             = real(lchnk,8)
          Var_Flat(        2)                             = real(i_save,8)
          Var_Flat(        3)                             = state_loc%ps(i_save)
          Var_Flat(        4)                             = state_loc%phis(i_save)
          Var_Flat(        5)                             = ztodt
          Var_Flat(        6)                             = precc(i_save)
          Var_Flat(        7)                             = precl(i_save)
          Var_Flat(        8)                             = precsc(i_save)
          Var_Flat(        9)                             = precsl(i_save)
          Var_Flat(       10)                             = cltot(i_save)
          Var_Flat(       11)                             = clhgh(i_save)
          Var_Flat(       12)                             = clmed(i_save)
          Var_Flat(       13)                             = cllow(i_save)
          Var_Flat(       14)                             = prectend(i_save)
          Var_Flat(       15)                             = precstend(i_save)
          Var_Flat(       16)                             = cam_in%ocnfrac(i_save)
          Var_Flat(       17)                             = wnd
          Var_Flat(       18)                             = tau00
          Var_Flat(       19)                             = bflx
          Var_Flat(       20)                             = fluxu0
          Var_Flat(       21)                             = fluxv0
          Var_Flat(       22)                             = fluxt0
          Var_Flat(       23)                             = fluxq0
          Var_Flat(       24)                             = taux_crm(i_save)
          Var_Flat(       25)                             = tauy_crm(i_save)
          Var_Flat(       26)                             = z0m(i_save)
          Var_Flat(       27)                             = timing_factor(i_save)
          Var_Flat(       28)                             = latitude0
          Var_Flat(       29)                             = longitude0
          Var_Flat(       30)                             = FlagEnd
          Var_Flat(        1+singlelen:   pver+singlelen) = state_loc%t(i_save,:)
          Var_Flat( 1*pver+1+singlelen: 2*pver+singlelen) = state_loc%q(i_save,:,1)
          Var_Flat( 2*pver+1+singlelen: 3*pver+singlelen) = state_loc%q(i_save,:,ixcldliq)
          Var_Flat( 3*pver+1+singlelen: 4*pver+singlelen) = state_loc%q(i_save,:,ixcldice)
          Var_Flat( 4*pver+1+singlelen: 5*pver+singlelen) = ul(:)
          Var_Flat( 5*pver+1+singlelen: 6*pver+singlelen) = vl(:) 
          Var_Flat( 6*pver+1+singlelen: 7*pver+singlelen) = state_loc%pmid(i_save,:)
          Var_Flat( 7*pver+1+singlelen: 8*pver+singlelen) = state_loc%pdel(i_save,:) 
          Var_Flat( 8*pver+1+singlelen: 9*pver+singlelen) = state_loc%zm(i_save,:)
          Var_Flat( 9*pver+1+singlelen:10*pver+singlelen+1) = state_loc%zi(i_save,:)
          Var_Flat(10*pver+2+singlelen:11*pver+singlelen+1) = cld(i_save,1:pver)
          Var_Flat(11*pver+2+singlelen:12*pver+singlelen+1) = cldtop(i_save,1:pver)
          Var_Flat(12*pver+2+singlelen:13*pver+singlelen+1) = gicewp(i_save,1:pver)
          Var_Flat(13*pver+2+singlelen:14*pver+singlelen+1) = gliqwp(i_save,1:pver)
          Var_Flat(14*pver+2+singlelen:15*pver+singlelen+1) = mctot(i_save,1:pver)
          Var_Flat(15*pver+2+singlelen:16*pver+singlelen+1) = mcup(i_save,1:pver)
          Var_Flat(16*pver+2+singlelen:17*pver+singlelen+1) = mcdn(i_save,1:pver)
          Var_Flat(17*pver+2+singlelen:18*pver+singlelen+1) = mcuup(i_save,1:pver)
          Var_Flat(18*pver+2+singlelen:19*pver+singlelen+1) = mcudn(i_save,1:pver)
          Var_Flat(19*pver+2+singlelen:20*pver+singlelen+1) = spqc(i_save,1:pver)
          Var_Flat(20*pver+2+singlelen:21*pver+singlelen+1) = spqi(i_save,1:pver)
          Var_Flat(21*pver+2+singlelen:22*pver+singlelen+1) = spqs(i_save,1:pver)
          Var_Flat(22*pver+2+singlelen:23*pver+singlelen+1) = spqg(i_save,1:pver)
          Var_Flat(23*pver+2+singlelen:24*pver+singlelen+1) = spqr(i_save,1:pver) 
          Var_Flat(24*pver+2+singlelen:25*pver+singlelen+1) = mu_crm(i_save,1:pver)
          Var_Flat(25*pver+2+singlelen:26*pver+singlelen+1) = md_crm(i_save,1:pver)
          Var_Flat(26*pver+2+singlelen:27*pver+singlelen+1) = du_crm(i_save,1:pver)
          Var_Flat(27*pver+2+singlelen:28*pver+singlelen+1) = eu_crm(i_save,1:pver)
          Var_Flat(28*pver+2+singlelen:29*pver+singlelen+1) = ed_crm(i_save,1:pver)
          Var_Flat(29*pver+2+singlelen:30*pver+singlelen+1) = tkez(i_save,1:pver)
          Var_Flat(30*pver+2+singlelen:31*pver+singlelen+1) = tkesgsz(i_save,1:pver)
          Var_Flat(31*pver+2+singlelen:32*pver+singlelen+1) = tk_crm(i_save,1:pver)
          Var_Flat(32*pver+2+singlelen:33*pver+singlelen+1) = flux_u(i_save,1:pver)
          Var_Flat(33*pver+2+singlelen:34*pver+singlelen+1) = flux_v(i_save,1:pver)
          Var_Flat(34*pver+2+singlelen:35*pver+singlelen+1) = flux_qt(i_save,1:pver)
          Var_Flat(35*pver+2+singlelen:36*pver+singlelen+1) = fluxsgs_qt(i_save,1:pver)
          Var_Flat(36*pver+2+singlelen:37*pver+singlelen+1) = flux_qp(i_save,1:pver)
          Var_Flat(37*pver+2+singlelen:38*pver+singlelen+1) = precflux(i_save,1:pver)
          Var_Flat(38*pver+2+singlelen:39*pver+singlelen+1) = qt_ls(i_save,1:pver)
          Var_Flat(39*pver+2+singlelen:40*pver+singlelen+1) = qt_trans(i_save,1:pver)
          Var_Flat(40*pver+2+singlelen:41*pver+singlelen+1) = qp_trans(i_save,1:pver)
          Var_Flat(41*pver+2+singlelen:42*pver+singlelen+1) = qp_fall(i_save,1:pver)
          Var_Flat(42*pver+2+singlelen:43*pver+singlelen+1) = qp_evp(i_save,1:pver)
          Var_Flat(43*pver+2+singlelen:44*pver+singlelen+1) = qp_src(i_save,1:pver)
          Var_Flat(44*pver+2+singlelen:45*pver+singlelen+1) = t_ls(i_save,1:pver)
          Var_Flat(45*pver+2+singlelen:46*pver+singlelen+1) = ptend_loc%q(i_save,:,1)
          Var_Flat(46*pver+2+singlelen:47*pver+singlelen+1) = ptend_loc%q(i_save,:,ixcldliq)
          Var_Flat(47*pver+2+singlelen:48*pver+singlelen+1) = ptend_loc%q(i_save,:,ixcldice)
          Var_Flat(48*pver+2+singlelen:49*pver+singlelen+1) = ptend_loc%s(i_save,:)
          Var_Flat(49*pver+2+singlelen:49*pver+21+singlelen)= qtotcrm(i_save,1:20)
          Var_Flat(49*pver+22+singlelen:49*pver+21+singlelen+pcols)=gcolindex(1:pcols) 
          fcount = 0
          do kk=1,crm_nz
            do jj=1,crm_ny
              do ii=crm_start_ind,crm_end_ind
                fcount = fcount + 1
                Var_Flat2(fcount)              = crm_u(i_save,ii,jj,kk)
                Var_Flat2(fcount + 1 * chnksz) = crm_v(i_save,ii,jj,kk)
                Var_Flat2(fcount + 2 * chnksz) = crm_w(i_save,ii,jj,kk)
                Var_Flat2(fcount + 3 * chnksz) = crm_t(i_save,ii,jj,kk)
                Var_Flat2(fcount + 4 * chnksz) = crm_qrad(i_save,ii,jj,kk)
                Var_Flat2(fcount + 5 * chnksz) = qc_crm(i_save,ii,jj,kk)
                Var_Flat2(fcount + 6 * chnksz) = qi_crm(i_save,ii,jj,kk)
                Var_Flat2(fcount + 7 * chnksz) = qpc_crm(i_save,ii,jj,kk)
                Var_Flat2(fcount + 8 * chnksz) = qpi_crm(i_save,ii,jj,kk)
                Var_Flat2(fcount + 9 * chnksz) = t_rad(i_save,ii,jj,kk)
                Var_Flat2(fcount + 10* chnksz) = qv_rad(i_save,ii,jj,kk)
                Var_Flat2(fcount + 11* chnksz) = qc_rad(i_save,ii,jj,kk)
                Var_Flat2(fcount + 12* chnksz) = qi_rad(i_save,ii,jj,kk)
                Var_Flat2(fcount + 13* chnksz) = cld_rad(i_save,ii,jj,kk)
                Var_Flat2(fcount + 14* chnksz) = cld3d_crm(i_save,ii,jj,kk)
                Var_Flat2(fcount + 15* chnksz) = crm_tk(i_save,ii,jj,kk)
                Var_Flat2(fcount + 16* chnksz) = crm_tkh(i_save,ii,jj,kk)
              end do
            end do
          end do

          fcount = 17*chnksz
          do ll=1,nmicro_fields_total
            do kk=1,crm_nz
              do jj=1,crm_ny
                do ii=crm_start_ind,crm_end_ind
                  fcount=fcount+1
                  Var_Flat2(fcount) = crm_micro(i_save,ii,jj,kk,ll)
                end do
              end do
            end do
          end do
        !write(iulog,*) 'check na:',na,nb,nc
        !write(iulog,*) 'Send 1:',i_save,lchnk,crm_micro(i_save,1:10,:,1:5,1)
        !write(iulog,*) 'Send 3:',i_save,lchnk,crm_micro(i_save,1:10,:,1:5,3)
          !write (iulog,*),'Finish send: ', crm_micro(i_save,:,:,:,3) 
          fcount = 17*chnksz + orc_nx*orc_ny*crm_nz*nmicro_fields_total
          do jj=1,crm_ny
            do ii=crm_start_ind,crm_end_ind
              fcount = fcount + 1
              Var_Flat2(fcount) = prec_crm(i_save,ii,jj)
            end do
          end do

          fcount = 17*chnksz + orc_nx*orc_ny*crm_nz*nmicro_fields_total+orc_nx*orc_ny
          do kk=1,crm_nz
            fcount = fcount + 1
            Var_Flat2(fcount + 0 * crm_nz) = work_u0(kk)
            Var_Flat2(fcount + 1 * crm_nz) = work_v0(kk)
            Var_Flat2(fcount + 2 * crm_nz) = work_t0(kk)
            Var_Flat2(fcount + 3 * crm_nz) = work_t00(kk)
            Var_Flat2(fcount + 4 * crm_nz) = work_tabs0(kk)
            Var_Flat2(fcount + 5 * crm_nz) = work_q0(kk)
            Var_Flat2(fcount + 6 * crm_nz) = work_qv0(kk)
            Var_Flat2(fcount + 7 * crm_nz) = work_qn0(kk)
            Var_Flat2(fcount + 8 * crm_nz) = work_qp0(kk)
            Var_Flat2(fcount + 9 * crm_nz) = work_tke0(kk)
          end do
          call MPI_Send(Var_Flat(:)         ,flen,MPI_REAL8 ,dest,9018,MPI_COMM_WORLD,ierr)
          call MPI_Send(Var_Flat2(:)         ,flen2,MPI_REAL8,dest,9019,MPI_COMM_WORLD,ierr)
          end do !do iorc=1,orc_nsubdomains
        end if ! if (state%isorchestrated(i)) then
      end do !do i = 1,ncol
#endif

      do i = 1,ncol
         i_save  = i
         nsubdomains_x  = orc_nsubdomains_x
         nsubdomains_y  = orc_nsubdomains_y
         call crm_define_grid()
 
        if (state%isofflinecrm(i))then
         call crm (lchnk,      i_save,                                                                                            &
             state_loc%t(i_save,:),   state_loc%q(i_save,:,1),    state_loc%q(i_save,:,ixcldliq), state_loc%q(i_save,:,ixcldice),                &
             ul(:),              vl(:),                                                                                      &
             state_loc%ps(i_save),    state_loc%pmid(i_save,:),   state_loc%pdel(i_save,:),  state_loc%phis(i_save),                             &
             state_loc%zm(i_save,:),  state_loc%zi(i_save,:),     ztodt,                pver,                                          &
             ptend_loc%q(i_save,:,1), ptend_loc%q(i_save,:,ixcldliq),ptend_loc%q(i_save,:,ixcldice), ptend_loc%s(i_save,:),                      &
             crm_u(i_save,:,:,:),     crm_v(i_save,:,:,:),       crm_w(i_save,:,:,:),        crm_t(i_save,:,:,:),          crm_micro(i_save,:,:,:,:), &
             crm_qrad(i_save,:,:,:),                                                                                              &
             qc_crm(i_save,:,:,:),    qi_crm(i_save,:,:,:),      qpc_crm(i_save,:,:,:),      qpi_crm(i_save,:,:,:),                              &
             prec_crm(i_save,:,:),    t_rad(i_save,:,:,:),       qv_rad(i_save,:,:,:),                                                      &
             qc_rad(i_save,:,:,:),    qi_rad(i_save,:,:,:),      cld_rad(i_save,:,:,:),      cld3d_crm(i_save, :, :, :),                         &
#ifdef m2005
             nc_rad(i_save,:,:,:),    ni_rad(i_save,:,:,:),      qs_rad(i_save,:,:,:),       ns_rad(i_save,:,:,:),         wvar_crm(i_save,:,:,:),    &
             aut_crm(i_save,:,:,:),   acc_crm(i_save,:,:,:),     evpc_crm(i_save,:,:,:),     evpr_crm(i_save,:,:,:),       mlt_crm(i_save,:,:,:),     &
             sub_crm(i_save,:,:,:),   dep_crm(i_save,:,:,:),     con_crm(i_save,:,:,:),                                                     &
             aut_crm_a(i_save,:),     acc_crm_a(i_save,:),       evpc_crm_a(i_save,:),       evpr_crm_a(i_save,:),         mlt_crm_a(i_save,:),       &
             sub_crm_a(i_save,:),     dep_crm_a(i_save,:),       con_crm_a(i_save,:),                                                       &
#endif
             precc(i_save),           precl(i_save),             precsc(i_save),             precsl(i_save),                                     &
             cltot(i_save),           clhgh(i_save),             clmed(i_save),              cllow(i_save),                cld(i_save,:),  cldtop(i_save,:), &
             gicewp(i_save,:),        gliqwp(i_save,:),                                                                                &
             mctot(i_save,:),         mcup(i_save,:),            mcdn(i_save,:),             mcuup(i_save,:),              mcudn(i_save,:),           &
             spqc(i_save,:),          spqi(i_save,:),            spqs(i_save,:),             spqg(i_save,:),               spqr(i_save,:),            &
#ifdef m2005
             spnc(i_save,:),          spni(i_save,:),            spns(i_save,:),             spng(i_save,:),               spnr(i_save,:),            &
#ifdef MODAL_AERO
             naermod,            vaerosol,             hygro,                                                                &
#endif 
#endif
#ifdef SPCAM_CLUBB_SGS
             clubb_buffer(i_save,:,:,:,:),                                                                                        &
             crm_cld(i_save,:, :, :),                                                                                             &
             clubb_tk(i_save, :, :, :), clubb_tkh(i_save, :, :, :),                                                                    &
             relvar(i_save,:, :, :),  accre_enhan(i_save, :, :, :),  qclvar(i_save, :, :, :),                                               &
#endif
             crm_tk(i_save, :, :, :), crm_tkh(i_save, :, :, :),                                                                        &
             mu_crm(i_save,:),        md_crm(i_save,:),          du_crm(i_save,:),           eu_crm(i_save,:),                                   & 
             ed_crm(i_save,:),        jt_crm(i_save),            mx_crm(i_save),                                                            &
#ifdef m2005
             abnd(i_save,:,:,:,:),    abnd_tf(i_save,:,:,:,:),   massflxbnd(i_save,:,:,:,:), acen(i_save,:,:,:,:),         acen_tf(i_save,:,:,:,:),   &
             rhcen(i_save,:,:,:,:),   qcloudcen(i_save,:,:,:,:), qicecen(i_save,:,:,:,:),    qlsink_afcen(i_save,:,:,:,:),                       &
             precrcen(i_save,:,:,:,:),                      precsolidcen(i_save,:,:,:,:),                                              &
             qlsink_bfcen(i_save,:,:,:,:),                  qlsink_avgcen(i_save,:,:,:,:),                       praincen(i_save,:,:,:,:),  &
             wupthresh_bnd(i_save,:), wdownthresh_bnd(i_save,:),                                                                       &
             wwqui_cen(i_save,:),     wwqui_bnd(i_save,:),       wwqui_cloudy_cen(i_save,:), wwqui_cloudy_bnd(i_save,:),                         &
#endif
             tkez(i_save,:),          tkesgsz(i_save,:),         tk_crm(i_save, :),                                                         &
             flux_u(i_save,:),        flux_v(i_save,:),          flux_qt(i_save,:),          fluxsgs_qt(i_save,:),         flux_qp(i_save,:),         &
             precflux(i_save,:),      qt_ls(i_save,:),           qt_trans(i_save,:),         qp_trans(i_save,:),           qp_fall(i_save,:),         &
             qp_evp(i_save,:),        qp_src(i_save,:),          t_ls(i_save,:),             prectend(i_save),             precstend(i_save),         &
             cam_in%ocnfrac(i_save),  wnd,                  tau00,                 bflx,                                          & 
             fluxu0,             fluxv0,               fluxt0,                fluxq0,                                        & 
             taux_crm(i_save),        tauy_crm(i_save),          z0m(i_save),                timing_factor(i_save),        qtotcrm(i_save, :)         )   

        end if !if (state%isofflinecrm(i))then
        !if (state%isorchestrated(i))then
        !write(iulog,*) 'After default crm:',i_save,crm_micro(i_save,:,:,:,3)
        !end if
        ! end if   ! if (state%isorchestrated(i))then
       end do  !do i = 1,ncol
       write (iulog,*),'Finish call crm module at rank: ',myrank_global
       do i = 1,ncol
        ! re: next args in call to crm () are ptend* which despite being declared inout are actually output.

          ! ****************
          ! pritch left off here, need to deal with GCM containers for crm-scale variables.
          ! PROBLEM: *can't* just make a flattened buffer vector for all crm_* variables 
          ! that are intent(inout) because orchestrator means CRM dims can be different on 
          ! different GCM columns. 

          ! If we want orchestrator to allow that maybe we need to (1) disable typical CRM_* pbuf / CRM-scale output
          ! stuff using #ifndef ORCHESTRATOR preprocessor directives.
          ! and likewise (2) redefine links between GCM containers like "crm_u" and
          ! CRM-intrinsic prognostics e.g. "u" 
          ! inside crmx_crm_module.F90 etc with analogous #ifdef shunts.
          ! We would still need to pull out the info for radiative transfer that
          ! happens outside the CRM but would need to fix the # of columns used
          ! for RT across all CRMs for that to happen.

          ! This is all adding up to a pretty overwhelming amount of data
          ! structure work / debugging.
          ! A gentler thought then is to keep CRM dims fixed and allow the MPI-decomp of orchestrator to
          ! offset the smaller time step of UP but not have UP domains contain
          ! more columns. Then again this may not make sense since with UP
          ! vertical grids, the CFL is controlled by dz more than dx.

          ! Maybe a middle ground is to keep the current crm_* variables
          ! including in pbuf but re-imagine them (for orchestrated big-CRMs)
          ! to be undersamplings (e.g. for the radiation) of the full domain.
          ! will need to be careful that doesn't re-initialize corners of the
          ! domain interior to crm_module.

          ! @ Peter what are your thoughts on this?
          ! *** MEANWHILE pritch is forging ahead under the middle ground above.
          ! **********************

#ifdef m2005
          call endrun('crm_physics: ORCHESTRATOR not compatible with m2005 microphysics') ! more stuff would be needed to pass to/fro CRM for this.
#endif
          ! ----------------------- receive outputs from orchestrator
          ! left off here, need to _Recv all the intent(inout) stuff.
         ! left off: need to add some ability to cross-check received output for
         ! a state variable or two (e.g. crm-scale temperature)
         ! against standard inline solution below to build credibility / as a
         ! unit test.

! Liran Start Here   ======================================================
         ! write (iulog,*) 'Liran Check Start Receving data from CRM'
          ! ---------- flattened multi-dim CRM inputs, prepackaged in
          ! crm_physics -------
#ifdef ORCHESTRATOR
         if (state%isorchestrated(i)) then
         do iorc=1,orc_nsubdomains
         !if (state%crmrank0(i) .ne. -2) then ! is this column coupled to an external CRM?
          call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
          call MPI_Recv(CRM_Var_Flat,fleno+nflat,MPI_REAL8,MPI_ANY_SOURCE,8006,MPI_COMM_WORLD,status,ierr)
          out_precc               = CRM_Var_Flat(1 ) 
          out_precl               = CRM_Var_Flat(2 )
          out_precsc              = CRM_Var_Flat(3 )
          out_precsl              = CRM_Var_Flat(4 )
          out_cltot               = CRM_Var_Flat(5 )
          out_clhgh               = CRM_Var_Flat(6 )
          out_clmed               = CRM_Var_Flat(7 )
          out_cllow               = CRM_Var_Flat(8 )
          out_prectend            = CRM_Var_Flat(9 )
          out_precstend           = CRM_Var_Flat(10)
          out_ocnfrac             = CRM_Var_Flat(11)
          out_wnd                 = CRM_Var_Flat(12)
          out_tau00               = CRM_Var_Flat(13)
          out_bflx                = CRM_Var_Flat(14)
          out_fluxu0              = CRM_Var_Flat(15)
          out_fluxv0              = CRM_Var_Flat(16)
          out_fluxt0              = CRM_Var_Flat(17)
          out_fluxq0              = CRM_Var_Flat(18)
          out_taux_crm            = CRM_Var_Flat(19)
          out_tauy_crm            = CRM_Var_Flat(20)
          out_z0m                 = CRM_Var_Flat(21)
          out_timing_factor       = CRM_Var_Flat(22)
          out_global_rank         = CRM_Var_Flat(23)
          out_it                  = CRM_Var_Flat(24)
          out_jt                  = CRM_Var_Flat(25)
          out_jt_crm              = CRM_Var_Flat(26)
          out_mx_crm              = CRM_Var_Flat(27)
!print*,'test_new add',out_jt_crm,out_mx_crm
          out_qltend              = CRM_Var_Flat(        1+singleleno:  pver+singleleno)
          out_qcltend             = CRM_Var_Flat( 1*pver+1+singleleno:2*pver+singleleno) 
          out_qiltend             = CRM_Var_Flat( 2*pver+1+singleleno:3*pver+singleleno) 
          out_sltend              = CRM_Var_Flat( 3*pver+1+singleleno:4*pver+singleleno)
          out_cld(:)              = CRM_Var_Flat( 4*pver+1+singleleno:5*pver+singleleno)
          out_cldtop(:)           = CRM_Var_Flat( 5*pver+1+singleleno:6*pver+singleleno)
          out_gicewp(:)           = CRM_Var_Flat( 6*pver+1+singleleno:7*pver+singleleno) 
          out_gliqwp(:)           = CRM_Var_Flat( 7*pver+1+singleleno:8*pver+singleleno)
          out_mctot(:)            = CRM_Var_Flat( 8*pver+1+singleleno:9*pver+singleleno)
          out_mcup(:)             = CRM_Var_Flat( 9*pver+1+singleleno:10*pver+singleleno) 
          out_mcdn(:)             = CRM_Var_Flat(10*pver+1+singleleno:11*pver+singleleno)
          out_mcuup(:)            = CRM_Var_Flat(11*pver+1+singleleno:12*pver+singleleno)
          out_mcudn(:)            = CRM_Var_Flat(12*pver+1+singleleno:13*pver+singleleno)
          out_spqc(:)             = CRM_Var_Flat(13*pver+1+singleleno:14*pver+singleleno) 
          out_spqi(:)             = CRM_Var_Flat(14*pver+1+singleleno:15*pver+singleleno)
          out_spqs(:)             = CRM_Var_Flat(15*pver+1+singleleno:16*pver+singleleno) 
          out_spqg(:)             = CRM_Var_Flat(16*pver+1+singleleno:17*pver+singleleno) 
          out_spqr(:)             = CRM_Var_Flat(17*pver+1+singleleno:18*pver+singleleno) 
          out_mu_crm(:)           = CRM_Var_Flat(18*pver+1+singleleno:19*pver+singleleno)
          out_md_crm(:)           = CRM_Var_Flat(19*pver+1+singleleno:20*pver+singleleno)
          out_du_crm(:)           = CRM_Var_Flat(20*pver+1+singleleno:21*pver+singleleno)
          out_eu_crm(:)           = CRM_Var_Flat(21*pver+1+singleleno:22*pver+singleleno)
          out_ed_crm(:)           = CRM_Var_Flat(22*pver+1+singleleno:23*pver+singleleno) 
          out_tkez(:)             = CRM_Var_Flat(23*pver+1+singleleno:24*pver+singleleno)
          out_tkesgsz(:)          = CRM_Var_Flat(24*pver+1+singleleno:25*pver+singleleno)
          out_tk_crm(:)           = CRM_Var_Flat(25*pver+1+singleleno:26*pver+singleleno)
          out_flux_u(:)           = CRM_Var_Flat(26*pver+1+singleleno:27*pver+singleleno)
          out_flux_v(:)           = CRM_Var_Flat(27*pver+1+singleleno:28*pver+singleleno)
          out_flux_qt(:)          = CRM_Var_Flat(28*pver+1+singleleno:29*pver+singleleno) 
          out_fluxsgs_qt(:)       = CRM_Var_Flat(29*pver+1+singleleno:30*pver+singleleno) 
          out_flux_qp(:)          = CRM_Var_Flat(30*pver+1+singleleno:31*pver+singleleno) 
          out_precflux(:)         = CRM_Var_Flat(31*pver+1+singleleno:32*pver+singleleno)
          out_qt_ls(:)            = CRM_Var_Flat(32*pver+1+singleleno:33*pver+singleleno) 
          out_qt_trans(:)         = CRM_Var_Flat(33*pver+1+singleleno:34*pver+singleleno) 
          out_qp_trans(:)         = CRM_Var_Flat(34*pver+1+singleleno:35*pver+singleleno)
          out_qp_fall(:)          = CRM_Var_Flat(35*pver+1+singleleno:36*pver+singleleno)
          out_qp_evp(:)           = CRM_Var_Flat(36*pver+1+singleleno:37*pver+singleleno)
          out_qp_src(:)           = CRM_Var_Flat(37*pver+1+singleleno:38*pver+singleleno)
          out_t_ls(:)             = CRM_Var_Flat(38*pver+1+singleleno:39*pver+singleleno)
          out_qtotcrm(1:20)       = CRM_Var_Flat(39*pver+1+singleleno:39*pver+20+singleleno)

          it = int(out_it)
          jt = int(out_jt)
          icheck_iorc = mod(int(out_global_rank-npes),orc_nsubdomains)
          icheck_column = (int(out_global_rank)-npes-icheck_iorc+1)/orc_nsubdomains+1   
          fcount = structleno*pver+20+singleleno+1
          do kk=1,crm_nz
            do jj=1,orc_ny
              do ii=1,orc_nx
                fcount = fcount + 1
                out_crm_u(ii,jj,kk)    = CRM_Var_Flat(fcount)             
                out_crm_v(ii,jj,kk)    = CRM_Var_Flat(fcount + 1 * chnksz) 
                out_crm_w(ii,jj,kk)    = CRM_Var_Flat(fcount + 2 * chnksz) 
                out_crm_t(ii,jj,kk)    = CRM_Var_Flat(fcount + 3 * chnksz) 
                out_crm_qrad(ii,jj,kk) = CRM_Var_Flat(fcount + 4 * chnksz) 
                out_qc_crm(ii,jj,kk)   = CRM_Var_Flat(fcount + 5 * chnksz) 
                out_qi_crm(ii,jj,kk)   = CRM_Var_Flat(fcount + 6 * chnksz) 
                out_qpc_crm(ii,jj,kk)  = CRM_Var_Flat(fcount + 7 * chnksz) 
                out_qpi_crm(ii,jj,kk)  = CRM_Var_Flat(fcount + 8 * chnksz) 
                out_t_rad(ii,jj,kk)    = CRM_Var_Flat(fcount + 9 * chnksz) 
                out_qv_rad(ii,jj,kk)   = CRM_Var_Flat(fcount + 10* chnksz)
                out_qc_rad(ii,jj,kk)   = CRM_Var_Flat(fcount + 11* chnksz)
                out_qi_rad(ii,jj,kk)   = CRM_Var_Flat(fcount + 12* chnksz) 
                out_cld_rad(ii,jj,kk)  = CRM_Var_Flat(fcount + 13* chnksz) 
                out_cld3d_crm(ii,jj,kk)= CRM_Var_Flat(fcount + 14* chnksz)
                out_crm_tk(ii,jj,kk)   = CRM_Var_Flat(fcount + 15* chnksz)
                out_crm_tkh(ii,jj,kk)  = CRM_Var_Flat(fcount + 16* chnksz)
              end do
            end do
          end do

          fcount = structleno*pver+20+singleleno + 1 + 17*chnksz
          do ll=1,nmicro_fields_total
            do kk=1,crm_nz
              do jj=1,orc_ny
                do ii=1,orc_nx
                  fcount=fcount+1
                  out_crm_micro(ii,jj,kk,ll) = CRM_Var_Flat(fcount)
                end do
              end do
            end do
          end do

          fcount = structleno*pver+20+singleleno + 17*chnksz + 1 + orc_nx*orc_ny*crm_nz*nmicro_fields_total
          do jj=1,orc_ny
            do ii=1,orc_nx
              fcount = fcount + 1
              out_prec_crm(ii,jj) = CRM_Var_Flat(fcount)
            end do
          end do
         call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
         lchnk = lchnk_save
         if (state%isofflinecrm(i)) then
          do kk=1,crm_nz
            work_crm_u(kk)      = 0.0
            work_crm_v(kk)      = 0.0
            work_crm_w(kk)      = 0.0
            work_crm_t(kk)      = 0.0
            work_crm_micro1(kk) = 0.0
            work_crm_micro2(kk) = 0.0
            work_crm_micro3(kk) = 0.0
            work_crm_qrad(kk)   = 0.0
            work_qc_crm(kk)     = 0.0
            work_qi_crm(kk)     = 0.0
            work_qpc_crm(kk)    = 0.0
            work_qpi_crm(kk)    = 0.0
            work_t_rad(kk)      = 0.0
            work_qv_rad(kk)     = 0.0
            work_qc_rad(kk)     = 0.0
            work_qi_rad(kk)     = 0.0
            work_cld_rad(kk)    = 0.0
            work_cld3d_crm(kk)  = 0.0
            work_org_crm_u(kk)      = 0.0
            work_org_crm_v(kk)      = 0.0
            work_org_crm_w(kk)      = 0.0
            work_org_crm_t(kk)      = 0.0
            work_org_crm_micro1(kk) = 0.0
            work_org_crm_micro2(kk) = 0.0
            work_org_crm_micro3(kk) = 0.0
            work_org_crm_qrad(kk)   = 0.0
            work_org_qc_crm(kk)     = 0.0
            work_org_qi_crm(kk)     = 0.0
            work_org_qpc_crm(kk)    = 0.0
            work_org_qpi_crm(kk)    = 0.0
            work_org_t_rad(kk)      = 0.0
            work_org_qv_rad(kk)     = 0.0
            work_org_qc_rad(kk)     = 0.0
            work_org_qi_rad(kk)     = 0.0
            work_org_cld_rad(kk)    = 0.0
            work_org_cld3d_crm(kk)  = 0.0
            temp_diff_crm_u         = -99999.9 
            temp_diff_crm_v         = -99999.9
            temp_diff_crm_w         = -99999.9
            temp_diff_crm_t         = -99999.9
            temp_diff_crm_micro1         = -99999.9
            temp_diff_crm_micro2         = -99999.9
            temp_diff_crm_micro3         = -99999.9
            temp_diff_crm_qrad         = -99999.9
            temp_diff_qc_crm         = -99999.9
            temp_diff_qi_crm         = -99999.9
            temp_diff_qpc_crm         = -99999.9
            temp_diff_qpi_crm         = -99999.9
            temp_diff_t_rad         = -99999.9
            temp_diff_qv_rad         = -99999.9
            temp_diff_qi_rad         = -99999.9
            temp_diff_cld_rad         = -99999.9
            temp_diff_cld3d_crm         = -99999.9
            do jj=1,orc_ny
              do ii=1,orc_nx
                orc_crm_u    (icheck_column,ii+it,jj+jt,kk)    = out_crm_u    (ii,jj,kk)
                orc_crm_v    (icheck_column,ii+it,jj+jt,kk)    = out_crm_v    (ii,jj,kk)
                orc_crm_w    (icheck_column,ii+it,jj+jt,kk)    = out_crm_w    (ii,jj,kk)
                orc_crm_t    (icheck_column,ii+it,jj+jt,kk)    = out_crm_t    (ii,jj,kk)
                orc_crm_micro(icheck_column,ii+it,jj+jt,kk,:)  = out_crm_micro(ii,jj,kk,:) 
                orc_crm_qrad (icheck_column,ii+it,jj+jt,kk)    = out_crm_qrad (ii,jj,kk)
                orc_qc_crm   (icheck_column,ii+it,jj+jt,kk)    = out_qc_crm   (ii,jj,kk)
                orc_qi_crm   (icheck_column,ii+it,jj+jt,kk)    = out_qi_crm   (ii,jj,kk)
                orc_qpc_crm  (icheck_column,ii+it,jj+jt,kk)    = out_qpc_crm  (ii,jj,kk)
                orc_qpi_crm  (icheck_column,ii+it,jj+jt,kk)    = out_qpi_crm  (ii,jj,kk)
                orc_t_rad    (icheck_column,ii+it,jj+jt,kk)    = out_t_rad    (ii,jj,kk)
                orc_qv_rad   (icheck_column,ii+it,jj+jt,kk)    = out_qv_rad   (ii,jj,kk)
                orc_qc_rad   (icheck_column,ii+it,jj+jt,kk)    = out_qc_rad   (ii,jj,kk)
                orc_qi_rad   (icheck_column,ii+it,jj+jt,kk)    = out_qi_crm   (ii,jj,kk)
                orc_cld_rad  (icheck_column,ii+it,jj+jt,kk)    = out_cld_rad  (ii,jj,kk) 
                orc_cld3d_crm(icheck_column,ii+it,jj+jt,kk)    = out_cld3d_crm(ii,jj,kk)
                work_org_crm_u(kk)       = work_org_crm_u(kk)     + crm_u(icheck_column,ii+it,jj+jt,kk)
                work_org_crm_v(kk)       = work_org_crm_v(kk)     + crm_v(icheck_column,ii+it,jj+jt,kk)
                work_org_crm_w(kk)       = work_org_crm_w(kk)     + crm_w(icheck_column,ii+it,jj+jt,kk)
                work_org_crm_t(kk)       = work_org_crm_t(kk)     + crm_t(icheck_column,ii+it,jj+jt,kk)
                work_org_crm_micro1(kk)  = work_org_crm_micro1(kk) + crm_micro(icheck_column,ii+it,jj+jt,kk,1)
                work_org_crm_micro2(kk)  = work_org_crm_micro2(kk) + crm_micro(icheck_column,ii+it,jj+jt,kk,2)
                work_org_crm_micro3(kk)  = work_org_crm_micro3(kk) + crm_micro(icheck_column,ii+it,jj+jt,kk,3)
                work_org_crm_qrad(kk)    = work_org_crm_qrad(kk)  + crm_qrad(icheck_column,ii+it,jj+jt,kk)
                work_org_qc_crm(kk)      = work_org_qc_crm(kk)    + qc_crm(icheck_column,ii+it,jj+jt,kk)
                work_org_qi_crm(kk)      = work_org_qi_crm(kk)    + qi_crm(icheck_column,ii+it,jj+jt,kk)
                work_org_qpc_crm(kk)     = work_org_qpc_crm(kk)   + qpc_crm(icheck_column,ii+it,jj+jt,kk)
                work_org_qpi_crm(kk)     = work_org_qpi_crm(kk)   + qpi_crm(icheck_column,ii+it,jj+jt,kk)
                work_org_t_rad(kk)       = work_org_t_rad(kk)     + t_rad(icheck_column,ii+it,jj+jt,kk)
                work_org_qv_rad(kk)      = work_org_qv_rad(kk)    + qv_rad(icheck_column,ii+it,jj+jt,kk)
                work_org_qc_rad(kk)      = work_org_qc_rad(kk)    + qc_rad(icheck_column,ii+it,jj+jt,kk)
                work_org_qi_rad(kk)      = work_org_qi_rad(kk)    + qi_crm(icheck_column,ii+it,jj+jt,kk)
                work_org_cld_rad(kk)     = work_org_cld_rad(kk)   + cld_rad(icheck_column,ii+it,jj+jt,kk)
                work_org_cld3d_crm(kk)   = work_org_cld3d_crm(kk) + cld3d_crm(icheck_column,ii+it,jj+jt,kk)
                work_crm_u(kk)       = work_crm_u(kk)     + out_crm_u(ii,jj,kk)
                work_crm_v(kk)       = work_crm_v(kk)     + out_crm_v(ii,jj,kk)
                work_crm_w(kk)       = work_crm_w(kk)     + out_crm_w(ii,jj,kk)
                work_crm_t(kk)       = work_crm_t(kk)     + out_crm_t(ii,jj,kk)
                work_crm_micro1(kk)  = work_crm_micro1(kk) + out_crm_micro(ii,jj,kk,1)
                work_crm_micro2(kk)  = work_crm_micro2(kk) + out_crm_micro(ii,jj,kk,2)
                work_crm_micro3(kk)  = work_crm_micro3(kk) + out_crm_micro(ii,jj,kk,3) 
                work_crm_qrad(kk)    = work_crm_qrad(kk)  + out_crm_qrad(ii,jj,kk)
                work_qc_crm(kk)      = work_qc_crm(kk)    + out_qc_crm(ii,jj,kk)
                work_qi_crm(kk)      = work_qi_crm(kk)    + out_qi_crm(ii,jj,kk)
                work_qpc_crm(kk)     = work_qpc_crm(kk)   + out_qpc_crm(ii,jj,kk)
                work_qpi_crm(kk)     = work_qpi_crm(kk)   + out_qpi_crm(ii,jj,kk)
                work_t_rad(kk)       = work_t_rad(kk)     + out_t_rad(ii,jj,kk)
                work_qv_rad(kk)      = work_qv_rad(kk)    + out_qv_rad(ii,jj,kk)
                work_qc_rad(kk)      = work_qc_rad(kk)    + out_qc_rad(ii,jj,kk)
                work_qi_rad(kk)      = work_qi_rad(kk)    + out_qi_crm(ii,jj,kk)
                work_cld_rad(kk)     = work_cld_rad(kk)   + out_cld_rad(ii,jj,kk)
                work_cld3d_crm(kk)   = work_cld3d_crm(kk) + out_cld3d_crm(ii,jj,kk)
                diff_crm_u(kk)       = diff_crm_u(kk) + abs(out_crm_u(ii,jj,kk)-crm_u(icheck_column,ii+it,jj+jt,kk))
                maxdiff_crm_u(kk)    = max(temp_diff_crm_u,abs(out_crm_u(ii,jj,kk)-crm_u(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_crm_u      = maxdiff_crm_u(kk)
                diff_crm_v(kk)       = diff_crm_v(kk) + abs(out_crm_v(ii,jj,kk)-crm_v(icheck_column,ii+it,jj+jt,kk))
                maxdiff_crm_v(kk)    = max(temp_diff_crm_v,abs(out_crm_v(ii,jj,kk)-crm_v(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_crm_v      = maxdiff_crm_v(kk)
                diff_crm_w(kk)       = diff_crm_w(kk) + abs(out_crm_w(ii,jj,kk)-crm_w(icheck_column,ii+it,jj+jt,kk))
                maxdiff_crm_w(kk)    = max(temp_diff_crm_w,abs(out_crm_w(ii,jj,kk)-crm_w(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_crm_w      = maxdiff_crm_w(kk)
                diff_crm_t(kk)       = diff_crm_t(kk) + abs(out_crm_t(ii,jj,kk)-crm_t(icheck_column,ii+it,jj+jt,kk))
                maxdiff_crm_t(kk)    = max(temp_diff_crm_t,abs(out_crm_t(ii,jj,kk)-crm_t(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_crm_t      = maxdiff_crm_t(kk)
                diff_crm_micro1(kk)       = diff_crm_micro1(kk) + abs(out_crm_micro(ii,jj,kk,1)-crm_micro(icheck_column,ii+it,jj+jt,kk,1))
                maxdiff_crm_micro1(kk)    = max(temp_diff_crm_micro1,abs(out_crm_micro(ii,jj,kk,1)-crm_micro(icheck_column,ii+it,jj+jt,kk,1)))
                temp_diff_crm_micro1      = maxdiff_crm_micro1(kk)
                diff_crm_micro2(kk)       = diff_crm_micro2(kk) + abs(out_crm_micro(ii,jj,kk,2)-crm_micro(icheck_column,ii+it,jj+jt,kk,2))
                maxdiff_crm_micro2(kk)    = max(temp_diff_crm_micro2,abs(out_crm_micro(ii,jj,kk,2)-crm_micro(icheck_column,ii+it,jj+jt,kk,2)))
                temp_diff_crm_micro2      = maxdiff_crm_micro2(kk)
                diff_crm_micro3(kk)       = diff_crm_micro3(kk) + abs(out_crm_micro(ii,jj,kk,3)-crm_micro(icheck_column,ii+it,jj+jt,kk,3))
                maxdiff_crm_micro3(kk)    = max(temp_diff_crm_micro3,abs(out_crm_micro(ii,jj,kk,3)-crm_micro(icheck_column,ii+it,jj+jt,kk,3)))
                temp_diff_crm_micro3      = maxdiff_crm_micro3(kk)
                diff_crm_qrad(kk)       = diff_crm_qrad(kk) + abs(out_crm_qrad(ii,jj,kk)-crm_qrad(icheck_column,ii+it,jj+jt,kk))
                maxdiff_crm_qrad(kk)    = max(temp_diff_crm_qrad,abs(out_crm_qrad(ii,jj,kk)-crm_qrad(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_crm_qrad      = maxdiff_crm_qrad(kk)
                diff_qc_crm(kk)       = diff_qc_crm(kk) + abs(out_qc_crm(ii,jj,kk)-qc_crm(icheck_column,ii+it,jj+jt,kk))
                maxdiff_qc_crm(kk)    = max(temp_diff_qc_crm,abs(out_qc_crm(ii,jj,kk)-qc_crm(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_qc_crm      = maxdiff_qc_crm(kk)
                diff_qi_crm(kk)       = diff_qi_crm(kk) + abs(out_qi_crm(ii,jj,kk)-qi_crm(icheck_column,ii+it,jj+jt,kk))
                maxdiff_qi_crm(kk)    = max(temp_diff_qi_crm,abs(out_qi_crm(ii,jj,kk)-qi_crm(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_qi_crm      = maxdiff_qi_crm(kk)   
                diff_qpc_crm(kk)       = diff_qpc_crm(kk) + abs(out_qpc_crm(ii,jj,kk)-qpc_crm(icheck_column,ii+it,jj+jt,kk))
                maxdiff_qpc_crm(kk)    = max(temp_diff_qpc_crm,abs(out_qpc_crm(ii,jj,kk)-qpc_crm(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_qpc_crm      = maxdiff_qpc_crm(kk)
                diff_qpi_crm(kk)       = diff_qpi_crm(kk) + abs(out_qpi_crm(ii,jj,kk)-qpi_crm(icheck_column,ii+it,jj+jt,kk))
                maxdiff_qpi_crm(kk)    = max(temp_diff_qpi_crm,abs(out_qpi_crm(ii,jj,kk)-qpi_crm(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_qpi_crm      = maxdiff_qpi_crm(kk)
                diff_t_rad(kk)       = diff_t_rad(kk) + abs(out_t_rad(ii,jj,kk)-t_rad(icheck_column,ii+it,jj+jt,kk))
                maxdiff_t_rad(kk)    = max(temp_diff_t_rad,abs(out_t_rad(ii,jj,kk)-t_rad(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_t_rad      = maxdiff_t_rad(kk)
                diff_qv_rad(kk)       = diff_qv_rad(kk) + abs(out_qv_rad(ii,jj,kk)-qv_rad(icheck_column,ii+it,jj+jt,kk))
                maxdiff_qv_rad(kk)    = max(temp_diff_qv_rad,abs(out_qv_rad(ii,jj,kk)-qv_rad(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_qv_rad      = maxdiff_qv_rad(kk)
                diff_qi_rad(kk)       = diff_qi_rad(kk) + abs(out_qi_rad(ii,jj,kk)-qi_rad(icheck_column,ii+it,jj+jt,kk))
                maxdiff_qi_rad(kk)    = max(temp_diff_qi_rad,abs(out_qi_rad(ii,jj,kk)-qi_rad(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_qi_rad      = maxdiff_qi_rad(kk)
                diff_cld_rad(kk)       = diff_cld_rad(kk) + abs(out_cld_rad(ii,jj,kk)-cld_rad(icheck_column,ii+it,jj+jt,kk))
                maxdiff_cld_rad(kk)    = max(temp_diff_cld_rad,abs(out_cld_rad(ii,jj,kk)-cld_rad(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_cld_rad      = maxdiff_cld_rad(kk)
                diff_cld3d_crm(kk)       = diff_cld3d_crm(kk) + abs(out_cld3d_crm(ii,jj,kk)-cld3d_crm(icheck_column,ii+it,jj+jt,kk))
                maxdiff_cld3d_crm(kk)    = max(temp_diff_cld3d_crm,abs(out_cld3d_crm(ii,jj,kk)-cld3d_crm(icheck_column,ii+it,jj+jt,kk)))
                temp_diff_cld3d_crm      = maxdiff_cld3d_crm(kk)
              end do
            end do
            work_crm_u(kk)       = work_crm_u(kk) / (orc_nx*orc_ny)
            work_crm_v(kk)       = work_crm_v(kk) / (orc_nx*orc_ny)
            work_crm_w(kk)       = work_crm_w(kk) / (orc_nx*orc_ny)
            work_crm_t(kk)       = work_crm_t(kk) / (orc_nx*orc_ny)
            work_crm_micro1(kk)  = work_crm_micro1(kk) / (orc_nx*orc_ny)
            work_crm_micro2(kk)  = work_crm_micro2(kk) / (orc_nx*orc_ny)
            work_crm_micro3(kk)  = work_crm_micro3(kk) / (orc_nx*orc_ny)
            work_crm_qrad(kk)    = work_crm_qrad(kk) / (orc_nx*orc_ny)
            work_qc_crm(kk)      = work_qc_crm(kk) / (orc_nx*orc_ny)
            work_qi_crm(kk)      = work_qi_crm(kk) / (orc_nx*orc_ny)
            work_qpc_crm(kk)     = work_qpc_crm(kk) / (orc_nx*orc_ny)
            work_qpi_crm(kk)     = work_qpi_crm(kk) / (orc_nx*orc_ny)
            work_t_rad(kk)       = work_t_rad(kk) / (orc_nx*orc_ny)
            work_qv_rad(kk)      = work_qv_rad(kk) / (orc_nx*orc_ny)
            work_qc_rad(kk)      = work_qc_rad(kk) / (orc_nx*orc_ny)
            work_qi_rad(kk)      = work_qi_rad(kk) / (orc_nx*orc_ny)
            work_cld_rad(kk)     = work_cld_rad(kk)  / (orc_nx*orc_ny)
            work_cld3d_crm(kk)   = work_cld3d_crm(kk)  / (orc_nx*orc_ny)
            work_org_crm_u(kk)       = work_org_crm_u(kk) / (orc_nx*orc_ny)
            work_org_crm_v(kk)       = work_org_crm_v(kk) / (orc_nx*orc_ny)
            work_org_crm_w(kk)       = work_org_crm_w(kk) / (orc_nx*orc_ny)
            work_org_crm_t(kk)       = work_org_crm_t(kk) / (orc_nx*orc_ny)
            work_org_crm_micro1(kk)  = work_org_crm_micro1(kk) / (orc_nx*orc_ny)
            work_org_crm_micro2(kk)  = work_org_crm_micro2(kk) / (orc_nx*orc_ny)
            work_org_crm_micro3(kk)  = work_org_crm_micro3(kk) / (orc_nx*orc_ny)
            work_org_crm_qrad(kk)    = work_org_crm_qrad(kk) / (orc_nx*orc_ny)
            work_org_qc_crm(kk)      = work_org_qc_crm(kk) / (orc_nx*orc_ny)
            work_org_qi_crm(kk)      = work_org_qi_crm(kk) / (orc_nx*orc_ny)
            work_org_qpc_crm(kk)     = work_org_qpc_crm(kk) / (orc_nx*orc_ny)
            work_org_qpi_crm(kk)     = work_org_qpi_crm(kk) / (orc_nx*orc_ny)
            work_org_t_rad(kk)       = work_org_t_rad(kk) / (orc_nx*orc_ny)
            work_org_qv_rad(kk)      = work_org_qv_rad(kk) / (orc_nx*orc_ny)
            work_org_qc_rad(kk)      = work_org_qc_rad(kk) / (orc_nx*orc_ny)
            work_org_qi_rad(kk)      = work_org_qi_rad(kk) / (orc_nx*orc_ny)
            work_org_cld_rad(kk)     = work_org_cld_rad(kk)  / (orc_nx*orc_ny)
            work_org_cld3d_crm(kk)   = work_org_cld3d_crm(kk)  / (orc_nx*orc_ny)
            diff_crm_u(kk)           = diff_crm_u(kk)  / (orc_nx*orc_ny)
            sqrt_crm_u(kk)           = sqrt(diff_crm_u(kk)) 
            diff_crm_v(kk)           = diff_crm_v(kk)  / (orc_nx*orc_ny)
            sqrt_crm_v(kk)           = sqrt(diff_crm_v(kk))
            diff_crm_w(kk)           = diff_crm_w(kk)  / (orc_nx*orc_ny)
            sqrt_crm_w(kk)           = sqrt(diff_crm_w(kk))
            diff_crm_t(kk)           = diff_crm_t(kk)  / (orc_nx*orc_ny)
            sqrt_crm_t(kk)           = sqrt(diff_crm_t(kk))
            diff_crm_micro1(kk)           = diff_crm_micro1(kk)  / (orc_nx*orc_ny)
            sqrt_crm_micro1(kk)           = sqrt(diff_crm_micro1(kk))
            diff_crm_micro2(kk)           = diff_crm_micro2(kk)  / (orc_nx*orc_ny)
            sqrt_crm_micro2(kk)           = sqrt(diff_crm_micro2(kk))
            diff_crm_micro3(kk)           = diff_crm_micro3(kk)  / (orc_nx*orc_ny)
            sqrt_crm_micro3(kk)           = sqrt(diff_crm_micro3(kk))
            diff_crm_qrad(kk)           = diff_crm_qrad(kk)  / (orc_nx*orc_ny)
            sqrt_crm_qrad(kk)           = sqrt(diff_crm_qrad(kk))
            diff_qc_crm(kk)           = diff_qc_crm(kk)  / (orc_nx*orc_ny)
            sqrt_qc_crm(kk)           = sqrt(diff_qc_crm(kk))
            diff_qi_crm(kk)           = diff_qi_crm(kk)  / (orc_nx*orc_ny)
            sqrt_qi_crm(kk)           = sqrt(diff_qi_crm(kk))
            diff_qpc_crm(kk)           = diff_qpc_crm(kk)  / (orc_nx*orc_ny)
            sqrt_qpc_crm(kk)           = sqrt(diff_qpc_crm(kk))
            diff_qpi_crm(kk)           = diff_qpi_crm(kk)  / (orc_nx*orc_ny)
            sqrt_qpi_crm(kk)           = sqrt(diff_qpi_crm(kk))
            diff_t_rad(kk)           = diff_t_rad(kk)  / (orc_nx*orc_ny)
            sqrt_t_rad(kk)           = sqrt(diff_t_rad(kk))
            diff_qv_rad(kk)           = diff_qv_rad(kk)  / (orc_nx*orc_ny)
            sqrt_qv_rad(kk)           = sqrt(diff_qv_rad(kk))
            diff_qc_rad(kk)           = diff_qc_rad(kk)  / (orc_nx*orc_ny)
            sqrt_qc_rad(kk)           = sqrt(diff_qc_rad(kk))
            diff_qi_rad(kk)           = diff_qi_rad(kk)  / (orc_nx*orc_ny)
            sqrt_qi_rad(kk)           = sqrt(diff_qi_rad(kk))
            diff_cld_rad(kk)           = diff_cld_rad(kk)  / (orc_nx*orc_ny)
            sqrt_cld_rad(kk)           = sqrt(diff_cld_rad(kk))
            diff_cld3d_crm(kk)           = diff_cld3d_crm(kk)  / (orc_nx*orc_ny)
            sqrt_cld3d_crm(kk)           = sqrt(diff_cld3d_crm(kk))
            write(iulog,1011),int(out_global_rank),int(icheck_column),kk,work_crm_u(kk),work_org_crm_u(kk),diff_crm_u(kk),sqrt_crm_u(kk),maxdiff_crm_u(kk)
            write(iulog,1012),int(out_global_rank),int(icheck_column),kk,work_crm_v(kk),work_org_crm_v(kk),diff_crm_v(kk),sqrt_crm_v(kk),maxdiff_crm_v(kk)
            write(iulog,1013),int(out_global_rank),int(icheck_column),kk,work_crm_w(kk),work_org_crm_w(kk),diff_crm_w(kk),sqrt_crm_w(kk),maxdiff_crm_w(kk)
            write(iulog,1014),int(out_global_rank),int(icheck_column),kk,work_crm_t(kk),work_org_crm_t(kk),diff_crm_t(kk),sqrt_crm_t(kk),maxdiff_crm_t(kk)
            write(iulog,1015),int(out_global_rank),int(icheck_column),kk,work_crm_micro1(kk),work_org_crm_micro1(kk),diff_crm_micro1(kk),sqrt_crm_micro1(kk),maxdiff_crm_micro1(kk)
            write(iulog,1016),int(out_global_rank),int(icheck_column),kk,work_crm_micro2(kk),work_org_crm_micro2(kk),diff_crm_micro2(kk),sqrt_crm_micro2(kk),maxdiff_crm_micro2(kk)
            write(iulog,1017),int(out_global_rank),int(icheck_column),kk,work_crm_micro3(kk),work_org_crm_micro3(kk),diff_crm_micro3(kk),sqrt_crm_micro3(kk),maxdiff_crm_micro3(kk)
            write(iulog,1018),int(out_global_rank),int(icheck_column),kk,work_crm_qrad(kk),work_org_crm_qrad(kk),diff_crm_qrad(kk),sqrt_crm_qrad(kk),maxdiff_crm_qrad(kk)
            write(iulog,1019),int(out_global_rank),int(icheck_column),kk,work_qc_crm(kk),work_org_qc_crm(kk),diff_qc_crm(kk),sqrt_qc_crm(kk),maxdiff_qc_crm(kk)
            write(iulog,1020),int(out_global_rank),int(icheck_column),kk,work_qi_crm(kk),work_org_qi_crm(kk),diff_qi_crm(kk),sqrt_qi_crm(kk),maxdiff_qi_crm(kk)
            write(iulog,1021),int(out_global_rank),int(icheck_column),kk,work_qpc_crm(kk),work_org_qpc_crm(kk),diff_qpc_crm(kk),sqrt_qpc_crm(kk),maxdiff_qpc_crm(kk)
            write(iulog,1022),int(out_global_rank),int(icheck_column),kk,work_qpi_crm(kk),work_org_qpi_crm(kk),diff_qpi_crm(kk),sqrt_qpi_crm(kk),maxdiff_qpi_crm(kk)
            write(iulog,1023),int(out_global_rank),int(icheck_column),kk,work_t_rad(kk),work_org_t_rad(kk),diff_t_rad(kk),sqrt_t_rad(kk),maxdiff_t_rad(kk)
            write(iulog,1024),int(out_global_rank),int(icheck_column),kk,work_qv_rad(kk),work_org_qv_rad(kk),diff_qv_rad(kk),sqrt_qv_rad(kk),maxdiff_qv_rad(kk)
            write(iulog,1025),int(out_global_rank),int(icheck_column),kk,work_qc_rad(kk),work_org_qc_rad(kk),diff_qc_rad(kk),sqrt_qc_rad(kk),maxdiff_qc_rad(kk)
            write(iulog,1026),int(out_global_rank),int(icheck_column),kk,work_qi_rad(kk),work_org_qi_rad(kk),diff_qi_rad(kk),sqrt_qi_rad(kk),maxdiff_qi_rad(kk)
            write(iulog,1027),int(out_global_rank),int(icheck_column),kk,work_cld_rad(kk),work_org_cld_rad(kk),diff_cld_rad(kk),sqrt_cld_rad(kk),maxdiff_cld_rad(kk) 
            write(iulog,1028),int(out_global_rank),int(icheck_column),kk,work_cld3d_crm(kk),work_org_cld3d_crm(kk),diff_cld3d_crm(kk),sqrt_cld3d_crm(kk),maxdiff_cld3d_crm(kk)
          end do
1011  format ('=>crm_u',3I4,5E15.6)
1012  format ('=>crm_v',3I4,5E15.6)
1013  format ('=>crm_w',3I4,5E15.6)
1014  format ('=>crm_t',3I4,5E15.6)
1015  format ('=>crm_micro1',3I4,5E15.6)
1016  format ('=>crm_micro2',3I4,5E15.6)
1017  format ('=>crm_micro3',3I4,5E15.6)
1018  format ('=>crm_qrad',3I4,5E15.6)
1019  format ('=>qc_crm',3I4,5E15.6)
1020  format ('=>qi_crm',3I4,5E15.6)
1021  format ('=>qpc_crm',3I4,5E15.6)
1022  format ('=>qpi_crm',3I4,5E15.6)
1023  format ('=>t_rad',3I4,5E15.6)
1024  format ('=>qv_rad',3I4,5E15.6)
1025  format ('=>qc_rad',3I4,5E15.6)
1026  format ('=>qi_rad',3I4,5E15.6)
1027  format ('=>cld_rad',3I4,5E15.6)
1028  format ('=>cld3d_crm',3I4,5E15.6)
            do jj=1,orc_ny
              do ii=1,orc_nx
                orc_prec_crm(i,ii+it,jj+jt)  = out_prec_crm(ii,jj)
              end do
            end do
          !write (iulog,*),'Check print 17',int(out_global_rank),it,jt,icheck_column,icheck_iorc,out_crm_u
          
          do kk=1,crm_nz
            do jj=1,orc_ny
              do ii=1,orc_nx
                write(iulog,1111),int(out_global_rank),int(icheck_column),ii,kk,int(ztodt),crm_u(icheck_column,ii+it,jj+jt,kk),orc_crm_u(icheck_column,ii+it,jj+jt,kk)
                write(iulog,1112),int(out_global_rank),int(icheck_column),ii,kk,int(ztodt),crm_w(icheck_column,ii+it,jj+jt,kk),orc_crm_w(icheck_column,ii+it,jj+jt,kk)
                write(iulog,1113),int(out_global_rank),int(icheck_column),ii,kk,int(ztodt),crm_t(icheck_column,ii+it,jj+jt,kk),orc_crm_t(icheck_column,ii+it,jj+jt,kk)
                write(iulog,1114),int(out_global_rank),int(icheck_column),ii,kk,int(ztodt),crm_micro(icheck_column,ii+it,jj+jt,kk,1),orc_crm_micro(icheck_column,ii+it,jj+jt,kk,1)
                write(iulog,1115),int(out_global_rank),int(icheck_column),ii,kk,int(ztodt),crm_micro(icheck_column,ii+it,jj+jt,kk,2),orc_crm_micro(icheck_column,ii+it,jj+jt,kk,2)
                write(iulog,1116),int(out_global_rank),int(icheck_column),ii,kk,int(ztodt),crm_micro(icheck_column,ii+it,jj+jt,kk,3),orc_crm_micro(icheck_column,ii+it,jj+jt,kk,3)
              end do
            end do
          end do
1111  format ('C_u',5I4,2E15.6)
1112  format ('C_w',5I4,2E15.6)
1113  format ('C_t',5I4,2E15.6)
1114  format ('C_m1',5I4,2E15.6)
1115  format ('C_m2',5I4,2E15.6)
1116  format ('C_m3',5I4,2E15.6)
        else
          cltot  (icheck_column)                           = out_cltot 
          clhgh  (icheck_column)                           = out_clhgh
          clmed  (icheck_column)                           = out_clmed
          cllow  (icheck_column)                           = out_cllow
          precc  (icheck_column)                           = out_precc
          precl  (icheck_column)                           = out_precl
          precsc (icheck_column)                           = out_precsc
          precsl (icheck_column)                           = out_precsl
          wnd                                              = out_wnd
          tau00                                            = out_tau00
          bflx                                             = out_bflx
          fluxu0                                           = out_fluxu0
          fluxv0                                           = out_fluxv0
          fluxt0                                           = out_fluxt0
          fluxq0                                           = out_fluxq0
          taux_crm(icheck_column)                          = out_taux_crm
          tauy_crm(icheck_column)                          = out_tauy_crm
          z0m(icheck_column)                               = out_z0m
          timing_factor(icheck_column)                     = out_timing_factor
          jt_crm(icheck_column)                            = out_jt_crm
          mx_crm(icheck_column)                            = out_mx_crm

          do kk=1,crm_nz
            do jj=1,orc_ny
              do ii=1,orc_nx
               crm_tk   (icheck_column,ii+it,jj+jt,kk)            = out_crm_tk   (ii, jj, kk)
               crm_tkh  (icheck_column,ii+it,jj+jt,kk)            = out_crm_tkh  (ii, jj, kk) 
               crm_u    (icheck_column,ii+it,jj+jt,kk)            = out_crm_u    (ii, jj, kk)
               crm_v    (icheck_column,ii+it,jj+jt,kk)            = out_crm_v    (ii, jj, kk)
               crm_w    (icheck_column,ii+it,jj+jt,kk)            = out_crm_w    (ii, jj, kk)
               crm_t    (icheck_column,ii+it,jj+jt,kk)            = out_crm_t    (ii, jj, kk)
               t_rad    (icheck_column,ii+it,jj+jt,kk)            = out_t_rad    (ii, jj, kk)
               qv_rad   (icheck_column,ii+it,jj+jt,kk)            = out_qv_rad   (ii, jj, kk)
               qc_rad   (icheck_column,ii+it,jj+jt,kk)            = out_qc_rad   (ii, jj, kk)
               qi_rad   (icheck_column,ii+it,jj+jt,kk)            = out_qi_rad   (ii, jj, kk)
               cld_rad  (icheck_column,ii+it,jj+jt,kk)            = out_cld_rad  (ii, jj, kk)
               cld3d_crm(icheck_column,ii+it,jj+jt,kk)            = out_cld3d_crm(ii, jj, kk)
              end do
            end do
          end do

            do jj=1,orc_ny
              do ii=1,orc_nx
               prec_crm(icheck_column,ii+it,jj+jt)              = out_prec_crm (ii, jj)
              end do
            end do


          do ll=1,nmicro_fields_total
            do kk=1,crm_nz
              do jj=1,orc_ny
                do ii=1,orc_nx
                 crm_micro(icheck_column,ii+it,jj+jt,kk,ll)       = out_crm_micro(ii,jj,kk,ll)
                end do
              end do
            end do
          end do

        do kk=1,crm_nz
          ptend_loc%q(icheck_column,kk,1)                  = out_qltend(kk)
          ptend_loc%q(icheck_column,kk,ixcldliq)           = out_qcltend(kk)
          ptend_loc%q(icheck_column,kk,ixcldice)           = out_qiltend(kk)
          ptend_loc%s(icheck_column,kk)                    = out_sltend(kk)
          cld(icheck_column,kk)                            = out_cld(kk)
          cldtop(icheck_column,kk)                         = out_cldtop(kk)
          gicewp(icheck_column,kk)                         = out_gicewp(kk)
          gliqwp(icheck_column,kk)                         = out_gliqwp(kk)
          mctot(icheck_column,kk)                          = out_mctot (kk)
          mcup(icheck_column,kk)                           = out_mcup(kk)
          mcdn(icheck_column,kk)                           = out_mcdn(kk)
          mcuup(icheck_column,kk)                          = out_mcuup(kk)
          mcudn(icheck_column,kk)                          = out_mcudn (kk)
          spqc(icheck_column,kk)                           = out_spqc(kk)
          spqi(icheck_column,kk)                           = out_spqi(kk)
          spqs(icheck_column,kk)                           = out_spqs(kk)
          spqg(icheck_column,kk)                           = out_spqg(kk)
          spqr(icheck_column,kk)                           = out_spqr  (kk)
          mu_crm(icheck_column,kk)                         = out_mu_crm(kk)
          md_crm(icheck_column,kk)                         = out_md_crm(kk)
          du_crm(icheck_column,kk)                         = out_du_crm(kk)
          eu_crm(icheck_column,kk)                         = out_eu_crm(kk)
          ed_crm(icheck_column,kk)                         = out_ed_crm(kk)
          tkez(icheck_column,kk)                           = out_tkez(kk)
          tkesgsz(icheck_column,kk)                        = out_tkesgsz(kk)
          tk_crm(icheck_column,kk)                         = out_tk_crm(kk)
          flux_u(icheck_column,kk)                         = out_flux_u(kk)
          flux_v(icheck_column,kk)                         = out_flux_v(kk)
          flux_qt(icheck_column,kk)                        = out_flux_qt(kk)
          fluxsgs_qt(icheck_column,kk)                     = out_fluxsgs_qt(kk)
          flux_qp(icheck_column,kk)                        = out_flux_qp(kk)
          precflux(icheck_column,kk)                       = out_precflux(kk) 
          qt_ls(icheck_column,kk)                          = out_qt_ls(kk)
          qt_trans(icheck_column,kk)                       = out_qt_trans(kk)
          qp_trans(icheck_column,kk)                       = out_qp_trans(kk)
          qp_fall(icheck_column,kk)                        = out_qp_fall(kk)
          qp_evp(icheck_column,kk)                         = out_qp_evp(kk)
          qp_src(icheck_column,kk)                         = out_qp_src(kk)
          t_ls(icheck_column,kk)                           = out_t_ls(kk)
          prectend(icheck_column)                          = out_prectend
          precstend(icheck_column)                         = out_precstend
       end do

        do kk=1,20
          qtotcrm(icheck_column,kk)                        = out_qtotcrm(kk) 
        end do

        end if !if (state%isofflinecrm(i)) then
       end do ! do iorc=1,orc_nsubdomains
        endif ! is this column column coupled to an external (orchestrated) CRM?
       end do !do i = 1,ncol

#endif

       do i = 1,ncol
           ! Retrieve the values back out of the internal crm array structure
           if (is_spcam_sam1mom) then
              crm_qt(i,:,:,:) = crm_micro(i,:,:,:,1)
              crm_qp(i,:,:,:) = crm_micro(i,:,:,:,2)
              crm_qn(i,:,:,:) = crm_micro(i,:,:,:,3)
           else if (is_spcam_m2005) then 
              crm_qt(i,:,:,:) = crm_micro(i,:,:,:,1)
              crm_nc(i,:,:,:) = crm_micro(i,:,:,:,2)
              crm_qr(i,:,:,:) = crm_micro(i,:,:,:,3)
              crm_nr(i,:,:,:) = crm_micro(i,:,:,:,4)
              crm_qi(i,:,:,:) = crm_micro(i,:,:,:,5)
              crm_ni(i,:,:,:) = crm_micro(i,:,:,:,6)
              crm_qs(i,:,:,:) = crm_micro(i,:,:,:,7)
              crm_ns(i,:,:,:) = crm_micro(i,:,:,:,8)
              crm_qg(i,:,:,:) = crm_micro(i,:,:,:,9)
              crm_ng(i,:,:,:) = crm_micro(i,:,:,:,10)
              crm_qc(i,:,:,:) = crm_micro(i,:,:,:,11)
          endif
       if (i.eq.ncol) then
         FlagEnd = 1
       endif
       end do ! i (loop over ncol)
       call t_stopf('crm_call')
       call t_stampf(wall(2), usr(2), sys(2))
       wall(1) = wall(2)-wall(1)
       call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
       write(iulog,*),'Rank,time: ',myrank_global,wall(1)
       ! There is no separate convective and stratiform precip for CRM:
       precc(:ncol)  = precc(:ncol) + precl(:ncol)
       precl(:ncol)  = 0._r8
       precsc(:ncol) = precsc(:ncol) + precsl(:ncol)
       precsl(:ncol) = 0._r8

       prec_dp(:ncol)= precc(:ncol)
       snow_dp(:ncol)= precsc(:ncol)

       do m=1,crm_nz
         k = pver-m+1
         do i = 1,ncol
           crm_qrad(i,:,:,m) = crm_qrad(i,:,:,m) * state_loc%pdel(i,k) ! for energy conservation
         end do
       end do

       call outfld('PRES    ',state_loc%pmid ,pcols   ,lchnk   )
       call outfld('DPRES   ',state_loc%pdel ,pcols   ,lchnk   )
       call outfld('CRM_U   ',crm_u          ,pcols   ,lchnk   )
       call outfld('CRM_V   ',crm_v          ,pcols   ,lchnk   )
       call outfld('CRM_W   ',crm_w          ,pcols   ,lchnk   )
       call outfld('CRM_T   ',crm_t          ,pcols   ,lchnk   )
       call outfld('CRM_QC  ',qc_crm         ,pcols   ,lchnk   )
       call outfld('CRM_QI  ',qi_crm         ,pcols   ,lchnk   )
       call outfld('CRM_QPC ',qpc_crm        ,pcols   ,lchnk   )
       call outfld('CRM_QPI ',qpi_crm        ,pcols   ,lchnk   )
       call outfld('CRM_PREC',prec_crm       ,pcols   ,lchnk   )
       call outfld('CRM_TK ', crm_tk(:, :, :, :)   ,pcols   ,lchnk   )
       call outfld('CRM_TKH', crm_tkh(:, :, :, :)  ,pcols   ,lchnk   )

       if (is_spcam_sam1mom) then
          tmp4d(:ncol,:,:,:) = crm_qt(:ncol,:,:,:)-qc_crm(:ncol,:,:,:)-qi_crm(:ncol,:,:,:)
          call outfld('CRM_QV  ',tmp4d,pcols   ,lchnk   )
       else if (is_spcam_m2005) then 
          tmp4d(:ncol,:,:,:) = crm_qt(:ncol,:,:,:)-qc_crm(:ncol,:,:,:)
          call outfld('CRM_QV  ',tmp4d, pcols   ,lchnk   )
       endif


       if (is_spcam_m2005) then
          call outfld('CRM_NC ',  crm_nc       ,pcols   ,lchnk)
          call outfld('CRM_NI ',  crm_ni       ,pcols   ,lchnk)
          call outfld('CRM_NR ',  crm_nr       ,pcols   ,lchnk)
          call outfld('CRM_NS ',  crm_ns       ,pcols   ,lchnk)
          call outfld('CRM_NG ',  crm_ng       ,pcols   ,lchnk)
          call outfld('CRM_WVAR', wvar_crm     ,pcols   ,lchnk)
          call outfld('CRM_QR ',  crm_qr       ,pcols   ,lchnk)
          call outfld('CRM_QS ',  crm_qs       ,pcols   ,lchnk)
          call outfld('CRM_QG ',  crm_qg       ,pcols   ,lchnk)
          call outfld('CRM_AUT',  aut_crm      ,pcols   ,lchnk)
          call outfld('CRM_ACC',  acc_crm      ,pcols   ,lchnk)
          call outfld('CRM_EVPC', evpc_crm     ,pcols   ,lchnk)
          call outfld('CRM_EVPR', evpr_crm     ,pcols   ,lchnk)
          call outfld('CRM_MLT',  mlt_crm      ,pcols   ,lchnk)
          call outfld('CRM_SUB',  sub_crm      ,pcols   ,lchnk)
          call outfld('CRM_DEP',  dep_crm      ,pcols   ,lchnk)
          call outfld('CRM_CON',  con_crm      ,pcols   ,lchnk)

          ! output for time-mean-avg
          call outfld('A_AUT',    aut_crm_a    , pcols  ,lchnk)
          call outfld('A_ACC',    acc_crm_a    , pcols  ,lchnk)
          call outfld('A_EVPC',   evpc_crm_a   , pcols  ,lchnk)
          call outfld('A_EVPR',   evpr_crm_a   , pcols  ,lchnk)
          call outfld('A_MLT',    mlt_crm_a    , pcols  ,lchnk)
          call outfld('A_SUB',    sub_crm_a    , pcols  ,lchnk)
          call outfld('A_DEP',    dep_crm_a    , pcols  ,lchnk)
          call outfld('A_CON',    con_crm_a    , pcols  ,lchnk)
       endif

       if(do_clubb_sgs) then
          call outfld('UP2     ' , clubb_buffer(:, :, :, :, 1)   ,pcols   ,lchnk)
          call outfld('VP2     ' , clubb_buffer(:, :, :, :, 2)   ,pcols   ,lchnk)
          call outfld('WPRTP   ' , clubb_buffer(:, :, :, :, 3)   ,pcols   ,lchnk)
          call outfld('WPTHLP  ' , clubb_buffer(:, :, :, :, 4)   ,pcols   ,lchnk)
          call outfld('WP2     ' , clubb_buffer(:, :, :, :, 5)   ,pcols   ,lchnk)
          call outfld('WP3     ' , clubb_buffer(:, :, :, :, 6)   ,pcols   ,lchnk)
          call outfld('RTP2    ' , clubb_buffer(:, :, :, :, 7)   ,pcols   ,lchnk)
          call outfld('THLP2   ' , clubb_buffer(:, :, :, :, 8)   ,pcols   ,lchnk)
          call outfld('RTPTHLP ' , clubb_buffer(:, :, :, :, 9)   ,pcols   ,lchnk)
          call outfld('UPWP    ' , clubb_buffer(:, :, :, :, 10)  ,pcols   ,lchnk)
          call outfld('VPWP    ' , clubb_buffer(:, :, :, :, 11)  ,pcols   ,lchnk)
          call outfld('CRM_CLD ' , clubb_buffer(:, :, :, :, 12)  ,pcols   ,lchnk)
          call outfld('T_TNDCY ' , clubb_buffer(:, :, :, :, 13)  ,pcols   ,lchnk)
          call outfld('QC_TNDCY' , clubb_buffer(:, :, :, :, 14)  ,pcols   ,lchnk)
          call outfld('QV_TNDCY' , clubb_buffer(:, :, :, :, 15)  ,pcols   ,lchnk)
          call outfld('CLUBB_TK ', clubb_tk(:, :, :, :)          ,pcols   ,lchnk)
          call outfld('CLUBB_TKH', clubb_tkh(:, :, :, :)         ,pcols   ,lchnk)
          call outfld('CRM_RELVAR', relvar(:, :, :, :)           ,pcols   ,lchnk)
          call outfld('QCLVAR'   , qclvar(:, :, :, :)            ,pcols   ,lchnk)
          call outfld('ACCRE_ENHAN', accre_enhan(:, :, :, :)     ,pcols   ,lchnk)
       
          spup2      = 0.0_r8; spvp2          = 0.0_r8; spwprtp    = 0.0_r8; spwpthlp  = 0.0_r8
          spwp2      = 0.0_r8; spwp3          = 0.0_r8; sprtp2     = 0.0_r8; spthlp2   = 0.0_r8
          sprtpthlp  = 0.0_r8; spupwp         = 0.0_r8; spvpwp     = 0.0_r8; spcrm_cld = 0.0_r8
          spt_tndcy  = 0.0_r8; spqc_tndcy     = 0.0_r8; spqv_tndcy = 0.0_r8
          spclubb_tk = 0.0_r8; spclubb_tkh    = 0.0_r8
          sprelvar   = 0.0_r8; spaccre_enhan  = 0.0_r8; spqclvar   = 0.0_r8

          do i=1, ncol
           do jj=1, crm_ny
           do ii=1, crm_nx
             do m=1, crm_nz+1
               k = pver-m+1
               spup2(i,k)      = spup2(i,k)      + clubb_buffer(i, ii, jj, m, 1)  / (crm_nx_ny)
               spvp2(i,k)      = spvp2(i,k)      + clubb_buffer(i, ii, jj, m, 2)  / (crm_nx_ny)
               spwprtp(i,k)    = spwprtp(i,k)    + clubb_buffer(i, ii, jj, m, 3)  / (crm_nx_ny)
               spwpthlp(i,k)   = spwpthlp(i,k)   + clubb_buffer(i, ii, jj, m, 4)  / (crm_nx_ny)
               spwp2(i,k)      = spwp2(i,k)      + clubb_buffer(i, ii, jj, m, 5)  / (crm_nx_ny)
               spwp3(i,k)      = spwp3(i,k)      + clubb_buffer(i, ii, jj, m, 6)  / (crm_nx_ny)
               sprtp2(i,k)     = sprtp2(i,k)     + clubb_buffer(i, ii, jj, m, 7)  / (crm_nx_ny)
               spthlp2(i,k)    = spthlp2(i,k)    + clubb_buffer(i, ii, jj, m, 8)  / (crm_nx_ny)  
               sprtpthlp(i,k)  = sprtpthlp(i,k)  + clubb_buffer(i, ii, jj, m, 9)  / (crm_nx_ny) 
               spupwp(i,k)     = spupwp(i,k)     + clubb_buffer(i, ii, jj, m, 10) / (crm_nx_ny)
               spupwp(i,k)     = spupwp(i,k)     + clubb_buffer(i, ii, jj, m, 11) / (crm_nx_ny)
               spcrm_cld(i,k)  = spcrm_cld(i,k)  + clubb_buffer(i, ii, jj, m, 12) / (crm_nx_ny)
               spt_tndcy(i,k)  = spt_tndcy(i,k)  + clubb_buffer(i, ii, jj, m, 13) / (crm_nx_ny)
               spqc_tndcy(i,k) = spqc_tndcy(i,k) + clubb_buffer(i, ii, jj, m, 14) / (crm_nx_ny)
               spqv_tndcy(i,k) = spqv_tndcy(i,k) + clubb_buffer(i, ii, jj, m, 15) / (crm_nx_ny)
             end do
             do m=1, crm_nz
               k = pver-m+1
               spclubb_tk(i,k)    = spclubb_tk(i,k)    + clubb_tk(i, ii, jj, m)    / (crm_nx_ny)
               spclubb_tkh(i,k)   = spclubb_tkh(i,k)   + clubb_tkh(i, ii, jj, m)   / (crm_nx_ny)
               sprelvar(i,k)      = sprelvar(i,k)      + relvar(i, ii, jj, m)      / (crm_nx_ny)
               spaccre_enhan(i,k) = spaccre_enhan(i,k) + accre_enhan(i, ii, jj, m) / (crm_nx_ny)
               spqclvar(i,k)      = spqclvar(i,k)      + qclvar(i, ii, jj, m)      / (crm_nx_ny)
             end do
           end do
           end do
          end do

          call outfld('SPUP2',         spup2         ,pcols   ,lchnk)
          call outfld('SPVP2',         spvp2         ,pcols   ,lchnk)
          call outfld('SPWPRTP',       spwprtp       ,pcols   ,lchnk)
          call outfld('SPWPTHLP',      spwpthlp      ,pcols   ,lchnk)
          call outfld('SPWP2',         spwp2         ,pcols   ,lchnk)
          call outfld('SPWP3',         spwp3         ,pcols   ,lchnk)
          call outfld('SPRTP2',        sprtp2        ,pcols   ,lchnk)
          call outfld('SPTHLP2',       spthlp2       ,pcols   ,lchnk)
          call outfld('SPRTPTHLP',     sprtpthlp     ,pcols   ,lchnk)
          call outfld('SPUPWP',        spupwp        ,pcols   ,lchnk)
          call outfld('SPVPWP',        spvpwp        ,pcols   ,lchnk)
          call outfld('SPCRM_CLD',     spcrm_cld     ,pcols   ,lchnk)
          call outfld('SPT_TNDCY',     spt_tndcy     ,pcols   ,lchnk)
          call outfld('SPQC_TNDCY',    spqc_tndcy    ,pcols   ,lchnk)
          call outfld('SPQV_TNDCY',    spqv_tndcy    ,pcols   ,lchnk)
          call outfld('SPCLUBB_TK ',   spclubb_tk    ,pcols   ,lchnk)
          call outfld('SPCLUBB_TKH',   spclubb_tkh   ,pcols   ,lchnk)
          call outfld('SPRELVAR',      sprelvar      ,pcols, lchnk)
          call outfld('SPACCRE_ENHAN', spaccre_enhan ,pcols, lchnk)
          call outfld('SPQCLVAR',      spqclvar      ,pcols, lchnk)
       endif ! if do_clubb_sgs

       spcld3d = 0.0_r8
       do i=1, ncol
         do jj=1, crm_ny
         do ii=1, crm_nx
           do m=1, crm_nz
             k = pver-m+1
             spcld3d(i,k) = spcld3d(i,k) + cld3d_crm(i,ii,jj,m) / (crm_nx_ny) 
           end do
         end do
         end do
       end do
       call outfld('SPCLD3D', spcld3d, pcols, lchnk)

       ifld = pbuf_get_index('QRL')
       call pbuf_get_field(pbuf, ifld, qrl)
       ifld = pbuf_get_index('QRS')
       call pbuf_get_field(pbuf, ifld, qrs)
       do k =1 , pver
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)/state_loc%pdel(i,k)
             qrl(i,k) = qrl(i,k)/state_loc%pdel(i,k)
          end do
       end do

       !
       ! add radiation tendencies to levels above CRM domain and 2 top CRM levels
       ! The radiation tendencies in the top 4 GCM levels are set to be zero in the CRM
       ptend_loc%s(:ncol, :pver-crm_nz+2) = qrs(:ncol,:pver-crm_nz+2)+qrl(:ncol,:pver-crm_nz+2)
   

       ! calculate the radiative fluxes from the radiation calculation
       ! This will be used to check energe conservations
       radflux(:) = 0.0_r8
       do k=1, pver
          do i=1, ncol
            radflux(i) = radflux(i) + (qrs(i,k)+qrl(i,k)) * state_loc%pdel(i,k)/gravit
          end do
       end do

       ftem(:ncol,:pver) = (ptend_loc%s(:ncol,:pver)-qrs(:ncol,:pver)-qrl(:ncol,:pver))/cpair

       tmp2d(:ncol,:) = qrl(:ncol,:)/cpair
       call outfld('SPQRL   ',tmp2d                     ,pcols   ,lchnk)

       tmp2d(:ncol,:) = qrs(:ncol,:)/cpair
       call outfld('SPQRS   ',tmp2d                     ,pcols   ,lchnk)

       call outfld('SPDT    ',ftem                      ,pcols   ,lchnk)
       call outfld('SPDQ    ',ptend_loc%q(1,1,1)        ,pcols   ,lchnk)
       call outfld('SPDQC   ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk)
       call outfld('SPDQI   ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk)
       call outfld('SPMC    ',mctot                     ,pcols   ,lchnk)
       call outfld('SPMCUP  ',mcup                      ,pcols   ,lchnk)
       call outfld('SPMCDN  ',mcdn                      ,pcols   ,lchnk)
       call outfld('SPMCUUP ',mcuup                     ,pcols   ,lchnk)
       call outfld('SPMCUDN ',mcudn                     ,pcols   ,lchnk)
       call outfld('SPQC    ',spqc                      ,pcols   ,lchnk)
       call outfld('SPQI    ',spqi                      ,pcols   ,lchnk)
       call outfld('SPQS    ',spqs                      ,pcols   ,lchnk)
       call outfld('SPQG    ',spqg                      ,pcols   ,lchnk)
       call outfld('SPQR    ',spqr                      ,pcols   ,lchnk)
       call outfld('SPQTFLX ',flux_qt                   ,pcols   ,lchnk)
       call outfld('SPUFLX  ',flux_u                    ,pcols   ,lchnk)
       call outfld('SPVFLX  ',flux_v                    ,pcols   ,lchnk)
       call outfld('SPTKE   ',tkez                      ,pcols   ,lchnk)
       call outfld('SPTKES  ',tkesgsz                   ,pcols   ,lchnk)
       call outfld('SPTK    ',tk_crm                    ,pcols   ,lchnk)
       call outfld('SPQTFLXS',fluxsgs_qt                ,pcols   ,lchnk)
       call outfld('SPQPFLX ',flux_qp                   ,pcols   ,lchnk)
       call outfld('SPPFLX  ',precflux                  ,pcols   ,lchnk)
       call outfld('SPQTLS  ',qt_ls                     ,pcols   ,lchnk)
       call outfld('SPQTTR  ',qt_trans                  ,pcols   ,lchnk)
       call outfld('SPQPTR  ',qp_trans                  ,pcols   ,lchnk)
       call outfld('SPQPEVP ',qp_evp                    ,pcols   ,lchnk)
       call outfld('SPQPFALL',qp_fall                   ,pcols   ,lchnk)
       call outfld('SPQPSRC ',qp_src                    ,pcols   ,lchnk)
       call outfld('SPTLS   ',t_ls                      ,pcols   ,lchnk)
       call outfld('CLOUDTOP',cldtop                    ,pcols   ,lchnk)
       call outfld('TIMINGF ',timing_factor             ,pcols   ,lchnk)

       if (is_spcam_m2005) then
          call outfld('SPNC    ',spnc         ,pcols   ,lchnk)
          call outfld('SPNI    ',spni         ,pcols   ,lchnk)
          call outfld('SPNS    ',spns         ,pcols   ,lchnk)
          call outfld('SPNG    ',spng         ,pcols   ,lchnk)
          call outfld('SPNR    ',spnr         ,pcols   ,lchnk)
       endif

       if (.not. do_clubb_sgs) then
          call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
          call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
          call outfld('CLDMED  ',clmed  ,pcols,lchnk)
          call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
          call outfld('CLOUD   ',cld,  pcols,lchnk)
       end if

       !
       ! Compute liquid water paths (for diagnostics only)
       tgicewp(:ncol) = 0._r8
       tgliqwp(:ncol) = 0._r8
       do k=1,pver
          do i = 1,ncol
             cicewp(i,k) = gicewp(i,k) * 1.0e-3_r8 / max(0.01_r8,cld(i,k)) ! In-cloud ice water path.  g/m2 --> kg/m2
             cliqwp(i,k) = gliqwp(i,k) * 1.0e-3_r8 / max(0.01_r8,cld(i,k)) ! In-cloud liquid water path. g/m2 --> kg/m2
             tgicewp(i)  = tgicewp(i) + gicewp(i,k) *1.0e-3_r8             ! grid cell mean ice water path.  g/m2 --> kg/m2
             tgliqwp(i)  = tgliqwp(i) + gliqwp(i,k) *1.0e-3_r8             ! grid cell mean ice water path.  g/m2 --> kg/m2
          end do
       end do
       tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
       gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver)
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)


       call outfld('SPCLDTOT',cltot  ,pcols,lchnk)
       call outfld('SPCLDHGH',clhgh  ,pcols,lchnk)
       call outfld('SPCLDMED',clmed  ,pcols,lchnk)
       call outfld('SPCLDLOW',cllow  ,pcols,lchnk)

       if(do_clubb_sgs) then
       ! Determine parameters for maximum/random overlap
#ifdef SPCAM_CLUBB_SGS
          call cldovrlap(lchnk, ncol, state%pint, cld, nmxrgn, pmxrgn)
          call cloud_cover_diags_out(lchnk, ncol, cld, state%pmid, nmxrgn, pmxrgn )
#endif
          deallocate(pmxrgn)
          deallocate(nmxrgn)
          deallocate(spup2)
          deallocate(spvp2)
          deallocate(spwprtp)
          deallocate(spwpthlp)
          deallocate(spwp2)
          deallocate(spwp3)
          deallocate(sprtp2)
          deallocate(spthlp2)
          deallocate(sprtpthlp)
          deallocate(spupwp)
          deallocate(spvpwp)
          deallocate(spcrm_cld)
          deallocate(spt_tndcy)
          deallocate(spqv_tndcy)
          deallocate(spqc_tndcy)
          deallocate(spclubb_tk)
          deallocate(spclubb_tkh)
          deallocate(sprelvar)
          deallocate(spaccre_enhan)
          deallocate(spqclvar)
          deallocate(crm_cld)
          deallocate(clubb_tk)
          deallocate(clubb_tkh)
          deallocate(relvar)
          deallocate(accre_enhan)
          deallocate(qclvar)
       endif

       call outfld('CLOUDTOP',cldtop,  pcols,lchnk)
       call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
       call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
       call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
       call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
       call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
       call outfld('ICLDIWP' ,cicewp , pcols,lchnk)

       ! Calculate fields which are needed elsewhere in CAM
       call pbuf_get_Field(pbuf, ast_idx, cld)  ! AST gets values in cld

       ! Find the cldtop for the physics buffer looking for the first location that has a value in the CRM cldtop field
       call pbuf_get_field(pbuf, cldtop_idx, cldtop_pbuf)
       cldtop_pbuf = pver
       do i=1,ncol
          do k=1,pver
             if (cldtop(i,k) > 1._r8/(crm_nx_ny)) then
                cldtop_pbuf(i)=k
                exit
             end if
          end do
       end do
       
       cs(:ncol, 1:pver) = state_loc%pmid(:ncol, 1:pver) / (287.15_r8*state_loc%t(:ncol, 1:pver))

       call pbuf_get_Field(pbuf, wsedl_idx, wsedl)
       if (is_spcam_m2005) then
          dunc_arr(:,:) = state_loc%q(:,:,ixnumliq)/ max(mincld,cld(:,:))
       else
          dunc_arr(:ncol,1:pver) = 100.e6_r8 / cs(:ncol,1:pver)
       end if
       do i=1,ncol
          do k=1,pver
             t    = state_loc%t(i,k)
             mu   = 1.496e-6_r8 * t**1.5_r8/(t+120._r8)
             acn  = gravit*rhow/(18._r8*mu)
             dumc = min( state_loc%q(i,k,ixcldliq) / max(mincld,cld(i,k)),0.005_r8 )
             dunc = dunc_arr(i,k)
             call size_dist_param_liq(mg_liq_props, dumc,dunc,cs(i,k),pgam,lamc)
             if (dumc >= qsmall) then
                wsedl(i,k)=acn*gamma(4._r8+bc+pgam)/(lamc**bc*gamma(pgam+4._r8))
             else
                wsedl(i,k)=0._r8
             endif
          end do
       end do

       if (is_spcam_m2005) then

          ! For convective transport
          do i=1, ncol
            ideep_crm(i) = i*1.0_r8
          end do
       endif
       call outfld('MU_CRM  ', mu_crm, pcols, lchnk)
       call outfld('MD_CRM  ', md_crm, pcols, lchnk)
       call outfld('EU_CRM  ', eu_crm, pcols, lchnk)
       call outfld('DU_CRM  ', du_crm, pcols, lchnk)
       call outfld('ED_CRM  ', ed_crm, pcols, lchnk)

! NAG requires ifdef because tk_crm_ecpp dereferened when not allocated
#ifdef m2005 
       if (is_spcam_m2005) then

         qlsinkcen = qlsink_avgcen

         ! copy local tk_crm into pbuf copy
         tk_crm_ecpp = tk_crm

         call outfld('ACEN    '     , acen               , pcols, lchnk)
         call outfld('ABND    '     , abnd               , pcols, lchnk)
         call outfld('ACEN_TF '     , acen_tf            , pcols, lchnk)
         call outfld('ABND_TF '     , abnd_tf            , pcols, lchnk)
         call outfld('MASFBND '     , massflxbnd         , pcols, lchnk)
         call outfld('RHCEN   '     , rhcen              , pcols, lchnk)
         call outfld('QCCEN   '     , qcloudcen          , pcols, lchnk)
         call outfld('QICEN   '     , qicecen            , pcols, lchnk)
         call outfld('QSINK_AFCEN'  , qlsink_afcen       , pcols, lchnk)
         call outfld('PRECRCEN'     , precrcen           , pcols, lchnk)
         call outfld('PRECSCEN'     , precsolidcen       , pcols, lchnk)
         call outfld('WUPTHRES'     , wupthresh_bnd      , pcols, lchnk)
         call outfld('WDNTHRES'     , wdownthresh_bnd    , pcols, lchnk)
         call outfld('WWQUI_CEN'    , wwqui_cen          , pcols, lchnk)
         call outfld('WWQUI_CLD_CEN', wwqui_cloudy_cen   , pcols, lchnk)
         call outfld('WWQUI_BND'    , wwqui_bnd          , pcols, lchnk)
         call outfld('WWQUI_CLD_BND', wwqui_cloudy_bnd   , pcols, lchnk)
         call outfld('QSINK_BFCEN'  , qlsink_bfcen       , pcols, lchnk)
         call outfld('QSINK_AVGCEN' , qlsink_avgcen      , pcols, lchnk)
         call outfld('PRAINCEN'     , praincen           , pcols, lchnk)
       endif
#endif

       if (is_spcam_m2005) then
            call cnst_get_ind('NUMLIQ', ixnumliq)
            call cnst_get_ind('NUMICE', ixnumice)
            ptend_loc%lq(ixnumliq)      = .TRUE.
            ptend_loc%lq(ixnumice)      = .TRUE.
            ptend_loc%q(:, :, ixnumliq) = 0._r8
            ptend_loc%q(:, :, ixnumice) = 0._r8

            do i = 1, ncol
             do k=1, crm_nz 
               m= pver-k+1
               do ii=1, crm_nx
               do jj=1, crm_ny
                 ptend_loc%q(i,m,ixnumliq) = ptend_loc%q(i,m,ixnumliq) + crm_nc(i,ii,jj,k) 
                 ptend_loc%q(i,m,ixnumice) = ptend_loc%q(i,m,ixnumice) + crm_ni(i,ii,jj,k)
               end do
               end do
               ptend_loc%q(i,m,ixnumliq) = (ptend_loc%q(i,m,ixnumliq)/(crm_nx_ny) - state_loc%q(i,m,ixnumliq))/ztodt
               ptend_loc%q(i,m,ixnumice) = (ptend_loc%q(i,m,ixnumice)/(crm_nx_ny) - state_loc%q(i,m,ixnumice))/ztodt
             end do
            end do
       end if

       ! Sum into overall ptend
       call physics_ptend_sum(ptend_loc, ptend, ncol)

       call physics_update(state_loc, ptend_loc, ztodt, tend_loc)

       ! calculate column water of rain, snow and graupel 
       if(is_spcam_m2005) then
         do i=1, ncol
            qt_hydro(i, 2)  = 0.0_r8
            qli_hydro(i, 2) = 0.0_r8
            qi_hydro(i, 2)  = 0.0_r8
            qtot(i, 3)      = 0.0_r8
            qt_cloud(i, 3)  = 0.0_r8
            qtv(i, 3)       = 0.0_r8
            do m=1, crm_nz
              k=pver-m+1
              do ii=1, crm_nx
              do jj=1, crm_ny
                qt_hydro(i,2)  = qt_hydro(i,2)  + (crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit 
                qli_hydro(i,2) = qli_hydro(i,2) + (crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit
                qi_hydro(i,2)  = qi_hydro(i,2)  + (crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit
                qtot(i, 3)     = qtot(i,3)      + (crm_qr(i,ii,jj,m)+crm_qs(i,ii,jj,m)+crm_qg(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit + (crm_qt(i,ii,jj,m)+crm_qi(i,ii,jj,m)) * state_loc%pdel(i,k)/gravit
                qt_cloud(i, 3) = qt_cloud(i, 3) + (crm_qt(i,ii,jj,m)+crm_qi(i,ii,jj,m)) * &
                                 state_loc%pdel(i,k)/gravit 
              end do
              end do
            end do
            qt_hydro(i,2)  = qt_hydro(i,2)  / (crm_nx_ny)
            qli_hydro(i,2) = qli_hydro(i,2) / (crm_nx_ny)
            qi_hydro(i,2)  = qi_hydro(i,2)  / (crm_nx_ny)
            qtot(i, 3)     = qtot(i, 3)     / (crm_nx_ny)
            qt_cloud(i, 3) = qt_cloud(i, 3) / (crm_nx_ny)
         end do 
       else if(is_spcam_sam1mom) then 
         do i=1, ncol
            qli_hydro(i, 2) = 0.0_r8
            qi_hydro(i, 2) = 0.0_r8
            do m=1, crm_nz
              k=pver-m+1
              do ii=1, crm_nx
              do jj=1, crm_ny
                sfactor        = max(0._r8,min(1._r8,(crm_t(i,ii,jj,m)-268.16_r8)*1._r8/(283.16_r8-268.16_r8)))
                qli_hydro(i,2) = qli_hydro(i,2)+crm_qp(i,ii,jj,m) * state_loc%pdel(i,k)/gravit
                qi_hydro(i,2)  = qi_hydro(i,2) +crm_qp(i,ii,jj,m) * (1-sfactor) * state_loc%pdel(i,k)/gravit
              end do
              end do
            end do
            qli_hydro(i,2) = qli_hydro(i,2) / (crm_nx_ny)
            qi_hydro(i,2)  = qi_hydro(i,2)  / (crm_nx_ny)

            ! total cloud water and total water vapor, and energy
            qt_cloud(i,2) = 0._r8
            qtv(i,2)      = 0._r8
            do k=1, pver
              qt_cloud(i,2) = qt_cloud(i,2) + (state_loc%q(i,k,ixcldliq)+state_loc%q(i,k,ixcldice)) * state_loc%pdel(i,k)/gravit
              qtv(i,2)      = qtv(i,2)      + state_loc%q(i,k,1)                                    * state_loc%pdel(i,k)/gravit
            end do
         end do
       end if

       ! check water and energy conservation
       call check_energy_chng(state_loc, tend_loc, "crm_tend", nstep, ztodt, zero, &
                prec_dp(:ncol)+(qli_hydro(:ncol,2)-qli_hydro(:ncol,1))/ztodt/1000._r8,  &
                snow_dp(:ncol)+(qi_hydro(:ncol,2)-qi_hydro(:ncol,1))/ztodt/1000._r8, radflux)

       !
       ! calculate total water after crm update
       ! total hydrometer water (rain, snow, and graupel)
       if (is_spcam_m2005) then
          do i=1, ncol

             ! total cloud water and total water vapor
             qt_cloud(i,2) = 0._r8
             qtv(i,2)      = 0._r8
             do k=1, pver
               qt_cloud(i,2) = qt_cloud(i,2) + (state_loc%q(i,k,ixcldliq)+state_loc%q(i,k,ixcldice)) * state_loc%pdel(i,k)/gravit
               qtv(i,2)      = qtv(i,2)      + state_loc%q(i,k,1) * state_loc%pdel(i,k)/gravit
             end do
             ! total water
             qtot(i,2) = qt_hydro(i,2) + qt_cloud(i,2) + qtv(i,2)

             ! to check water conservations
             if(abs((qtot(i,2)+(precc(i)+precl(i))*1000_r8*ztodt)-qtot(i,1))/qtot(i,1).gt.1.0e-5_r8) then 
                write(0, *) 'water before crm call', i, lchnk, qtot(i,1), qtv(i,1), qt_cloud(i,1), qt_hydro(i,1)
                write(0, *) 'water after crm call', i, lchnk, qtot(i,2)+(precc(i)+precl(i))*1000*ztodt, &
                             qtv(i,2), qt_cloud(i,2), qt_hydro(i,2), (precc(i)+precl(i))*1000*ztodt
                write(0, *) 'water, nstep, crm call2', nstep, i, lchnk, &
                             ((qtot(i,2)+(precc(i)+precl(i))*1000_r8*ztodt)-qtot(i,1))/qtot(i,1)
                write(0, *) 'water, calcualted in crm.F90', i, lchnk, qtotcrm(i, 1), qtotcrm(i, 9), &
                             qtot(i, 3)+(precc(i)+precl(i))*1000_r8*ztodt, qt_cloud(i, 3), qtv(i,2)+qt_cloud(i,2)
                write(0, *) 'water, temperature', i, lchnk, state_loc%t(i,pver)
             end if
          end do ! end i
       endif

   end if ! (is_first_step())

   call t_stopf('crm')

! ifdef needed because of use of dtstep_pp_input and spcam_modal_aero_wateruptake_dr
#ifdef m2005
   if (is_spcam_m2005) then
       call t_startf('bc_aerosols_mmf')

       where(qc_rad(:ncol,:,:,:crm_nz)+qi_rad(:ncol,:,:,:crm_nz) > 1.0e-10_r8) 
          cld_rad(:ncol,:,:,:crm_nz) = cld_rad(:ncol,:,:,:crm_nz)
       elsewhere
          cld_rad(:ncol,:,:,:crm_nz) = 0.0_r8
       endwhere

       ! temporarily turn on all lq, so it is allocated
       lq(:) = .true.
       call physics_ptend_init(ptend_loc, state_loc%psetcols, 'crm_physics', lq=lq)

       ! set all ptend%lq to false as they will be set in modal_aero_calcsize_sub
       ptend%lq(:) = .false.
       call modal_aero_calcsize_sub (state_loc, ptend_loc, ztodt, pbuf)
       call spcam_modal_aero_wateruptake_dr(state_loc, pbuf)

      ! Wet deposition is done in ECPP,
      ! So tendency from wet depostion is not updated in mz_aero_wet_intr (mz_aerosols_intr.F90)
      ! tendency from other parts of crmclouds_aerosol_wet_intr are still updated here.

      ! Sum into overall ptend
      call physics_ptend_sum(ptend_loc, ptend, ncol)
      call physics_update(state_loc, ptend_loc, ztodt, tend_loc)


      pblh_idx  = pbuf_get_index('pblh')
      call pbuf_get_field(pbuf, pblh_idx, pblh)

      !
      !   ECPP is called at every 3rd GCM time step.
      !   GCM time step is 10 minutes, and ECPP time step is 30 minutes.
      !
      dtstep_pp = dtstep_pp_input
      necpp = dtstep_pp/ztodt

      ! Only call ECPP every necpp th time step
      ! !!!BE CAUTIOUS (Minghuai Wang, 2017-02)!!!!: 
      ! ptend_loc from crmclouds_mixnuc_tend and parampollu_driver2 has
      ! to be multiplied by necpp, as the updates in state occure in tphysbc_spcam, 
      ! and the normal time step used in tphysbc_spcam is short 
      ! and ECPP time step is longer (by a facotr of ncecpp). 
      ! Otherwise, this will lead to underestimation in wet scavenging.
      ! 
      if(nstep.ne.0 .and. mod(nstep, necpp).eq.0) then
        call t_startf('crmclouds_mixnuc')

        call crmclouds_mixnuc_tend (state_loc, ptend_loc, dtstep_pp, cam_in%cflx, pblh, pbuf,  &
                                    wwqui_cen, wwqui_cloudy_cen, wwqui_bnd, wwqui_cloudy_bnd)

        ! scale ptend_loc by necpp 
        call physics_ptend_scale(ptend_loc, necpp*1.0_r8, ncol)
        ! Sum into overall ptend
        call physics_ptend_sum(ptend_loc, ptend, ncol)  
        call physics_update(state_loc, ptend_loc, ztodt, tend_loc)
        call t_stopf('crmclouds_mixnuc')

        call t_startf('ecpp')
        call parampollu_driver2(state_loc, ptend_loc, pbuf, dtstep_pp, dtstep_pp,  &
                                acen, abnd, acen_tf, abnd_tf, massflxbnd,   &
                                rhcen, qcloudcen, qlsinkcen, precrcen, precsolidcen, acldy_cen_tbeg )
        ! scale ptend_loc by necpp 
        call physics_ptend_scale(ptend_loc, necpp*1.0_r8, ncol)
        ! Sum into overall ptend
        call physics_ptend_sum(ptend_loc, ptend, ncol) 
        call physics_update(state_loc, ptend_loc, ztodt, tend_loc)
        call t_stopf ('ecpp')
      end if


      call t_stopf('bc_aerosols_mmf')
   endif ! /*m2005*/
#endif

   ! save for old cloud fraction in the MMF simulations
   cldo(:ncol, :) = cld(:ncol, :)

   deallocate(crm_micro)
   
   if (is_spcam_m2005) then
      deallocate(acen)
      deallocate(acen_tf)
      deallocate(rhcen)
      deallocate(qcloudcen)
      deallocate(qlsinkcen)
      deallocate(precrcen)
      deallocate(precsolidcen)
      deallocate(wwqui_cen)
      deallocate(wwqui_cloudy_cen)
      deallocate(abnd)
      deallocate(abnd_tf)
      deallocate(massflxbnd)
      deallocate(wwqui_bnd)
      deallocate(wwqui_cloudy_bnd)
      deallocate(qicecen)
      deallocate(qlsink_afcen)
      deallocate(qlsink_bfcen)
      deallocate(qlsink_avgcen)
      deallocate(praincen)
      deallocate(wupthresh_bnd)
      deallocate(wdownthresh_bnd)

      deallocate(na)
      deallocate(va)
      deallocate(hy)
      deallocate(naermod)
      deallocate(vaerosol)
      deallocate(hygro)
   end if

#endif

end subroutine crm_physics_tend

!=====================================================================================================

subroutine m2005_effradius(ql, nl,qi,ni,qs, ns, cld, pres, tk, effl, effi, effl_fn, deffi, lamcrad, pgamrad, des)
!-----------------------------------------------------------------------------------------------------
! 
! This subroutine is used to calculate droplet and ice crystal effective radius, which will be used 
! in the CAM radiation code. The method to calcualte effective radius is taken out of the Morrision's
! two momenent scheme from M2005MICRO_GRAUPEL. It is also very similar with the subroutine of effradius in 
! the module of cldwat2m in the CAM source codes. 
!
! Adopted by Minghuai Wang (Minghuai.Wang@pnl.gov). 
!
!-----------------------------------------------------------------------------------------------------
   ! ----------------------------------------------------------- !
   ! Calculate effective radius for pass to radiation code       !
   ! If no cloud water, default value is 10 micron for droplets, !
   ! 25 micron for cloud ice.                                    !
   ! Be careful of the unit of effective radius : [micro meter]  !
   ! ----------------------------------------------------------- !
  use shr_spfn_mod, only: gamma => shr_spfn_gamma
  implicit none

     real(r8), intent(in)    :: ql                ! Mean LWC of pixels [ kg/kg ]
     real(r8), intent(in)    :: nl                ! Grid-mean number concentration of cloud liquid droplet [#/kg]
     real(r8), intent(in)    :: qi                ! Mean IWC of pixels [ kg/kg ]
     real(r8), intent(in)    :: ni                ! Grid-mean number concentration of cloud ice    droplet [#/kg]
     real(r8), intent(in)    :: qs                ! mean snow water content [kg/kg]
     real(r8), intent(in)    :: ns                ! Mean snow crystal number concnetration [#/kg]
     real(r8), intent(in)    :: cld               ! Physical stratus fraction
     real(r8), intent(in)    :: pres               ! Air pressure [Pa] 
     real(r8), intent(in)    :: tk                ! air temperature [K]

     real(r8), intent(out)   :: effl              ! Effective radius of cloud liquid droplet [micro-meter]
     real(r8), intent(out)   :: effi              ! Effective radius of cloud ice    droplet [micro-meter]
     real(r8), intent(out)   :: effl_fn           ! effl for fixed number concentration of nlic = 1.e8
     real(r8), intent(out)   :: deffi             ! ice effective diameter for optics (radiation)
     real(r8), intent(out)   :: pgamrad           ! gamma parameter for optics (radiation)
     real(r8), intent(out)   :: lamcrad           ! slope of droplet distribution for optics (radiation)
     real(r8), intent(out)   :: des               ! snow effective diameter for optics (radiation) [micro-meter]

#ifdef CRM
     real(r8)  qlic                               ! In-cloud LWC [kg/m3]
     real(r8)  qiic                               ! In-cloud IWC [kg/m3]
     real(r8)  nlic                               ! In-cloud liquid number concentration [#/kg]
     real(r8)  niic                               ! In-cloud ice    number concentration [#/kg]

     real(r8)  cldm                               ! Constrained stratus fraction [no]
     real(r8)  mincld                             ! Minimum stratus fraction [no]

     real(r8)  lami, laml, lammax, lammin, pgam, lams, lammaxs, lammins

     real(r8)  dcs    !autoconversion size threshold   [meter]
     real(r8)  di, ci     ! cloud ice mass-diameter relationship
     real(r8)  ds, cs     ! snow crystal mass-diameter relationship 
     real(r8)  qsmall
     real(r8)  rho       ! air density [kg/m3]
     real(r8)  rhow      ! liquid water density [kg/m3]
     real(r8)  rhoi      ! ice density [kg/m3]
     real(r8)  rhos      ! snow density [kg/m3]
     real(r8)  res       ! effective snow diameters
     real(r8)  pi

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     pi = 3.1415926535897932384626434_r8
     qsmall = 1.0e-14_r8  ! in the SAM source code (module_mp_graupel)
     rhow = 997._r8       ! in module_mp_graupel, SAM
     rhoi = 500._r8       ! in both CAM and SAM

     dcs = 125.e-6_r8   ! in module_mp_graupel, SAM 
     ci = rhoi * pi/6._r8
     di = 3._r8

     !   for snow water
     rhos = 100._r8     ! in both SAM and CAM5 
     cs = rhos*pi/6._r8
     ds = 3._r8


     rho = pres / (287.15_r8*tk)    ! air density [kg/m3]

     mincld  = 0.0001_r8
     cldm    = max(cld,mincld)
     qlic    = min(5.e-3_r8,max(0._r8,ql/cldm))
     qiic    = min(5.e-3_r8,max(0._r8,qi/cldm))
     nlic    = max(nl,0._r8)/cldm
     niic    = max(ni,0._r8)/cldm

!------------------------------------------------------
!     Effective diameters of snow crystals
!------------------------------------------------------
     if(qs.gt.1.0e-7_r8) then 
       lammaxs=1._r8/10.e-6_r8
       lammins=1._r8/2000.e-6_r8
       lams = (gamma(1._r8+ds)*cs * ns/qs)**(1._r8/ds)
       lams = min(lammaxs,max(lams,lammins))
       res = 1.5_r8/lams*1.0e6_r8
     else
       res = 500._r8 
     end if 

     !
     ! from Hugh Morrision: rhos/917 accouts for assumptions about 
     ! ice density in the Mitchell optics. 
     !
     des = res * rhos/917._r8 *2._r8

   ! ------------------------------------- !
   ! Effective radius of cloud ice droplet !
   ! ------------------------------------- !

     if( qiic.ge.qsmall ) then
         niic   = min(niic,qiic*1.e20_r8)
         lammax = 1._r8/1.e-6_r8      ! in module_mp_graupel, SAM 
         lammin = 1._r8/(2._r8*dcs+100.e-6_r8)    ! in module_mp_graupel, SAM 
         lami   = (gamma(1._r8+di)*ci*niic/qiic)**(1._r8/di)
         lami   = min(lammax,max(lami,lammin))
         effi   = 1.5_r8/lami*1.e6_r8
     else
         effi   = 25._r8
     endif

     !--hm ice effective radius for david mitchell's optics
     !--ac morrison indicates that this is effective diameter
     !--ac morrison indicates 917 (for the density of pure ice..)
     deffi  = effi *rhoi/917._r8*2._r8

   ! ---------------------------------------- !
   ! Effective radius of cloud liquid droplet !
   ! ---------------------------------------- !

     if( qlic.ge.qsmall ) then
         !        Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM).
         !        See also Morrison and Grabowski (2007, JAS, Eq. (2))
         nlic   = min(nlic,qlic*1.e20_r8)

         ! set the minimum droplet number as 20/cm3.

         pgam   = 0.0005714_r8*(nlic*rho/1.e6_r8) + 0.2714_r8
         pgam   = 1._r8/(pgam**2)-1._r8
         pgam   = min(10._r8,max(pgam,2._r8))   ! in module_mp_graupel, SAM
         laml   = (pi/6._r8*rhow*nlic*gamma(pgam+4._r8)/(qlic*gamma(pgam+1._r8)))**(1._r8/3._r8)
         lammin = (pgam+1._r8)/50.e-6_r8    ! in cldwat2m, CAM
         lammax = (pgam+1._r8)/2.e-6_r8     ! in cldwat2m, CAM   ! cldwat2m should be used, 
                                                                 ! if lammax is too large, this will lead to crash in 
                                                                 ! src/physics/rrtmg/cloud_rad_props.F90 because 
                                                                 ! klambda-1 can be zero in gam_liquid_lw and gam_liquid_sw
                                                                 !  and g_lambda(kmu,klambda-1) will not be defined. 

         laml   = min(max(laml,lammin),lammax)
         effl   =  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8  ! in module_mp_graupel, SAM
         lamcrad  = laml 
         pgamrad  = pgam
     else
         ! chose 10. over 25, since 10 is a more reasonable value for liquid droplet
         effl   = 10._r8     ! in cldwat2m, CAM
         lamcrad  = 0.0_r8
         pgamrad  = 0.0_r8
     endif

   ! ---------------------------------------------------------------------- !
   ! Recalculate effective radius for constant number, in order to separate !
   ! first and second indirect effects. Assume constant number of 10^8 kg-1 !
   ! ---------------------------------------------------------------------- !

     nlic = 1.e8_r8
     if( qlic.ge.qsmall ) then
         !        Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM). 
         !        See also Morrison and Grabowski (2007, JAS, Eq. (2))  
         nlic   = min(nlic,qlic*1.e20_r8)
         pgam   = 0.0005714_r8*(nlic/1.e6_r8/rho) + 0.2714_r8
         pgam   = 1._r8/(pgam**2)-1._r8
         pgam   = min(10._r8,max(pgam,2._r8))   ! in module_mp_graupel, SAM
         laml   = (pi/6._r8*rhow*nlic*gamma(pgam+4._r8)/(qlic*gamma(pgam+1._r8)))**(1._r8/3._r8)
         lammin = (pgam+1._r8)/60.e-6_r8    ! in module_mp_graupel, SAM
         lammax = (pgam+1._r8)/1.e-6_r8     ! in module_mp_graupel, SAM

         laml   = min(max(laml,lammin),lammax)
         effl_fn   =  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8  ! in module_mp_graupel, SAM
     else
         ! chose 10. over 25, since 10 is a more reasonable value for liquid droplet.
         effl_fn   = 10._r8     ! in cldwat2m, CAM
     endif

   return
#endif
end subroutine m2005_effradius

end module crm_physics
