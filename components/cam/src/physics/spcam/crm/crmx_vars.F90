module crmx_vars

use crmx_grid
#ifdef CRM 
#ifdef MODAL_AERO 
use modal_aero_data,   only: ntot_amode
#endif
#endif

implicit none
!--------------------------------------------------------------------
! prognostic variables:
real, allocatable, dimension(:,:,:)  :: u
real, allocatable, dimension(:,:,:)  :: v
real, allocatable, dimension(:,:,:)  :: w
real, allocatable, dimension(:,:,:)  :: t

allocate (u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)) ! x-wind
allocate (v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)) ! y-wind
allocate (w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )) ! z-wind
allocate (t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)) ! liquid/ice water static energy 

!--------------------------------------------------------------------
! diagnostic variables:
real, allocatable, dimension(:,:,:)  :: p
real, allocatable, dimension(:,:,:)  :: tabs
real, allocatable, dimension(:,:,:)  :: qv
real, allocatable, dimension(:,:,:)  :: qcl
real, allocatable, dimension(:,:,:)  :: qpl
real, allocatable, dimension(:,:,:)  :: qci
real, allocatable, dimension(:,:,:)  :: qpi
real, allocatable, dimension(:,:,:)  :: tke2
real, allocatable, dimension(:,:,:)  :: tk2

allocate (p      (0:nx, (1-YES3D):ny, nzm))     ! perturbation pressure (from Poison eq)
allocate (tabs   (nx, ny, nzm)            )     ! temperature
allocate (qv      (nx, ny, nzm)           )     ! water vapor
allocate (qcl     (nx, ny, nzm)           )     ! liquid water  (condensate)
allocate (qpl     (nx, ny, nzm)           )     ! liquid water  (precipitation)
allocate (qci     (nx, ny, nzm)           )     ! ice water  (condensate)
allocate (qpi     (nx, ny, nzm)           )     ! ice water  (precipitation)

allocate (tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))   ! SGS TKE
allocate (tk2     (0:nxp1, (1-YES3D):nyp1, nzm)) ! SGS eddyviscosity
        
!--------------------------------------------------------------------
! time-tendencies for prognostic variables
real, allocatable, dimension(:,:,:,:)  :: dudt
real, allocatable, dimension(:,:,:,:)  :: dvdt
real, allocatable, dimension(:,:,:,:)  :: dwdt

allocate ( dudt   (nxp1, ny, nzm, 3))
allocate ( dvdt   (nx, nyp1, nzm, 3))
allocate ( dwdt   (nx, ny, nz,  3))

!----------------------------------------------------------------
! Temporary storage array:
real, allocatable, dimension(:,:,:)  :: misc
allocate ( misc(nx, ny, nz))
!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:
real, allocatable, dimension(:,:)  :: fluxbu
real, allocatable, dimension(:,:)  :: fluxbv
real, allocatable, dimension(:,:)  :: fluxbt
real, allocatable, dimension(:,:)  :: fluxbq
real, allocatable, dimension(:,:)  :: fluxtu
real, allocatable, dimension(:,:)  :: fluxtv
real, allocatable, dimension(:,:)  :: fluxtt
real, allocatable, dimension(:,:)  :: fluxtq
real, allocatable, dimension(:,:)  :: fzero
real, allocatable, dimension(:,:)  :: precsfc
real, allocatable, dimension(:,:)  :: precssfc

allocate ( fluxbu (nx, ny))
allocate ( fluxbv (nx, ny))
allocate ( fluxbt (nx, ny))
allocate ( fluxbq (nx, ny))
allocate ( fluxtu (nx, ny))
allocate ( fluxtv (nx, ny))
allocate ( fluxtt (nx, ny))
allocate ( fluxtq (nx, ny))
allocate ( fzero  (nx, ny))
allocate ( precsfc(nx,ny) ) ! surface precip. rate
allocate ( precssfc(nx,ny)) ! surface ice precip. rate
                
!-----------------------------------------------------------------
! profiles 
real, allocatable, dimension(:)  :: t0
real, allocatable, dimension(:)  :: q0
real, allocatable, dimension(:)  :: qv0
real, allocatable, dimension(:)  :: tabs0
real, allocatable, dimension(:)  :: tl0
real, allocatable, dimension(:)  :: tv0         
real, allocatable, dimension(:)  :: u0
real, allocatable, dimension(:)  :: v0
real, allocatable, dimension(:)  :: tg0
real, allocatable, dimension(:)  :: qg0
real, allocatable, dimension(:)  :: ug0
real, allocatable, dimension(:)  :: vg0
real, allocatable, dimension(:)  :: p0
real, allocatable, dimension(:)  :: tke0
real, allocatable, dimension(:)  :: t01
real, allocatable, dimension(:)  :: q01
real, allocatable, dimension(:)  :: qp0
real, allocatable, dimension(:)  :: qn0

allocate ( t0(nzm))
allocate ( q0(nzm))
allocate ( qv0(nzm))
allocate ( tabs0(nzm))
allocate ( tl0(nzm))
allocate ( tv0(nzm))
allocate ( u0(nzm))
allocate ( v0(nzm))
allocate ( tg0(nzm))
allocate ( qg0(nzm))
allocate ( ug0(nzm))
allocate ( vg0(nzm))
allocate ( p0(nzm))
allocate (tke0(nzm))
allocate ( t01(nzm))
allocate ( q01(nzm))
allocate ( qp0(nzm))
allocate ( qn0(nzm))
!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics 

real  sstobs, lhobs, shobs
!----------------------------------------------------------------
!  Domain top stuff:

real   gamt0    ! gradient of t() at the top,K/m
real   gamq0    ! gradient of q() at the top,g/g/m

!-----------------------------------------------------------------
! reference vertical profiles:
real, allocatable, dimension(:)  :: prespot
real, allocatable, dimension(:)  :: rho
real, allocatable, dimension(:)  :: rhow
real, allocatable, dimension(:)  :: bet
real, allocatable, dimension(:)  :: gamaz
real, allocatable, dimension(:)  :: wsub
real, allocatable, dimension(:)  :: qtend
real, allocatable, dimension(:)  :: ttend
real, allocatable, dimension(:)  :: utend
real, allocatable, dimension(:)  :: vtend
 
allocate ( prespot(nzm))  ! (1000./pres)**R/cp
allocate (  rho(nzm))	  ! air density at pressure levels,kg/m3 
allocate (  rhow(nz))   ! air density at vertical velocity levels,kg/m3
allocate (  bet(nzm))	  ! = ggr/tv0
allocate (  gamaz(nzm)) ! ggr/cp*z
allocate (  wsub(nz))   ! Large-scale subsidence velocity,m/s
allocate (  qtend(nzm)) ! Large-scale tendency for total water
allocate (  ttend(nzm)) ! Large-scale tendency for temp.
allocate (  utend(nzm)) ! Large-scale tendency for u
allocate (  vtend(nzm)) ! Large-scale tendency for v


!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nlsf	! number of large-scale forcing profiles
integer nrfc	! number of radiative forcing profiles
integer nsfc	! number of surface forcing profiles
integer nsnd	! number of observed soundings
integer nzlsf	! number of large-scale forcing profiles
integer nzrfc	! number of radiative forcing profiles
integer nzsnd	! number of observed soundings

real, allocatable :: dqls(:,:) ! Large-scale tendency for total water
real, allocatable :: dtls(:,:) ! Large-scale tendency for temp.
real, allocatable :: ugls(:,:) ! Large-scale wind in X-direction
real, allocatable :: vgls(:,:) ! Large-scale wind in Y-direction
real, allocatable :: wgls(:,:) ! Large-scale subsidence velocity,m/s
real, allocatable :: pres0ls(:)! Surface pressure, mb
real, allocatable :: zls(:,:)  ! Height
real, allocatable :: pls(:,:)  ! Pressure
real, allocatable :: dayls(:)  ! Large-scale forcing arrays time (days) 
real, allocatable :: dtrfc(:,:)! Radiative tendency for pot. temp.
real, allocatable :: dayrfc(:) ! Radiative forcing arrays time (days) 
real, allocatable :: prfc(:,:) ! Pressure/Height
real, allocatable :: sstsfc(:) ! SSTs
real, allocatable :: shsfc(:)   ! Sensible heat flux,W/m2
real, allocatable :: lhsfc(:)  ! Latent heat flux,W/m2
real, allocatable :: tausfc(:) ! Surface drag,m2/s2
real, allocatable :: daysfc(:) ! Surface forcing arrays time (days) 
real, allocatable :: usnd(:,:) ! Observed zonal wind
real, allocatable :: vsnd(:,:) ! Observed meriod wind
real, allocatable :: tsnd(:,:) ! Observed Abs. temperature
real, allocatable :: qsnd(:,:) ! Observed Moisture
real, allocatable :: zsnd(:,:) ! Height
real, allocatable :: psnd(:,:) ! Pressure
real, allocatable :: daysnd(:) ! number of sounding samples
 
!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!
real, allocatable ::  sstxy(:,:)	!  surface temperature xy-distribution
real, allocatable ::  fcory(:)      !  Coriolis parameter xy-distribution
real, allocatable ::  fcorzy(:)      !  z-Coriolis parameter xy-distribution
real, allocatable ::  latitude(:,:)	     ! latitude (degrees)
real, allocatable ::  longitude(:,:)	     ! longitude(degrees)
real, allocatable ::  prec_xy(:,:) ! mean precip. rate for outout
real, allocatable ::  shf_xy(:,:) ! mean precip. rate for outout
real, allocatable ::  lhf_xy(:,:) ! mean precip. rate for outout
real, allocatable ::  lwns_xy(:,:) ! mean net lw at SFC
real, allocatable ::  swns_xy(:,:) ! mean net sw at SFC
real, allocatable ::  lwnsc_xy(:,:) ! clear-sky mean net lw at SFC
real, allocatable ::  swnsc_xy(:,:) ! clear-sky mean net sw at SFC
real, allocatable ::  lwnt_xy(:,:) ! mean net lw at TOA
real, allocatable ::  swnt_xy(:,:) ! mean net sw at TOA
real, allocatable ::  lwntc_xy(:,:) ! clear-sky mean net lw at TOA
real, allocatable ::  swntc_xy(:,:) ! clear-sky mean net sw at TOA
real, allocatable ::  solin_xy(:,:) ! solar TOA insolation
real, allocatable ::  pw_xy(:,:)   ! precipitable water
real, allocatable ::  cw_xy(:,:)   ! cloud water path
real, allocatable ::  iw_xy(:,:)   ! ice water path
real, allocatable ::  cld_xy(:,:)   ! cloud frequency
real, allocatable ::  u200_xy(:,:) ! u-wind at 200 mb
real, allocatable ::  usfc_xy(:,:) ! u-wind at at the surface
real, allocatable ::  v200_xy(:,:) ! v-wind at 200 mb
real, allocatable ::  vsfc_xy(:,:) ! v-wind at the surface
real, allocatable ::  w500_xy(:,:) ! w at 500 mb
real, allocatable ::  qocean_xy(:,:) ! ocean cooling in W/m2

allocate ( sstxy(0:nx,(1-YES3D):ny))
allocate ( fcory(0:ny))
allocate ( fcorzy(ny) )
allocate ( latitude(nx,ny))
allocate ( longitude(nx,ny)) 
allocate ( prec_xy(nx,ny))
allocate ( shf_xy(nx,ny) )
allocate ( lhf_xy(nx,ny))
allocate ( lwns_xy(nx,ny))
allocate ( swns_xy(nx,ny))
allocate ( lwnsc_xy(nx,ny))
allocate ( swnsc_xy(nx,ny))
allocate ( lwnt_xy(nx,ny))
allocate ( swnt_xy(nx,ny))
allocate ( lwntc_xy(nx,ny))
allocate ( swntc_xy(nx,ny))
allocate ( solin_xy(nx,ny))
allocate ( pw_xy(nx,ny))
allocate ( cw_xy(nx,ny))
allocate ( iw_xy(nx,ny))
allocate ( cld_xy(nx,ny))
allocate ( u200_xy(nx,ny)) 
allocate ( usfc_xy(nx,ny))
allocate ( v200_xy(nx,ny))
allocate ( vsfc_xy(nx,ny))
allocate ( w500_xy(nx,ny))
allocate ( qocean_xy(nx,ny))

!----------------------------------------------------------------------
!	Vertical profiles of quantities sampled for statitistics purposes:
real, allocatable ::  twle(:)
real, allocatable ::  twsb(:)
real, allocatable ::  precflux(:)
real, allocatable ::  uwle(:)
real, allocatable ::  uwsb(:)
real, allocatable ::  vwle(:)
real, allocatable ::  vwsb(:)
real, allocatable ::  radlwup(:)
real, allocatable ::  radlwdn(:)
real, allocatable ::  radswup(:)
real, allocatable ::  radswdn(:)
real, allocatable ::  radqrlw(:)
real, allocatable ::  radqrsw(:)  
real, allocatable ::  tkeleadv(:)
real, allocatable ::  tkelepress(:)
real, allocatable ::  tkelediss(:)
real, allocatable ::  tkelediff(:)
real, allocatable ::  tkelebuoy(:)
real, allocatable ::  t2leadv(:)
real, allocatable ::  t2legrad(:)
real, allocatable ::  t2lediff(:)
real, allocatable ::  t2leprec(:)
real, allocatable ::  t2lediss(:)
real, allocatable ::  q2leadv(:)
real, allocatable ::  q2legrad(:)
real, allocatable ::  q2lediff(:)
real, allocatable ::  q2leprec(:)
real, allocatable ::  q2lediss(:)
real, allocatable ::  twleadv(:)
real, allocatable ::  twlediff(:)
real, allocatable ::  twlepres(:)
real, allocatable ::  twlebuoy(:)
real, allocatable ::  twleprec(:)
real, allocatable ::  qwleadv(:)
real, allocatable ::  qwlediff(:)
real, allocatable ::  qwlepres(:)
real, allocatable ::  qwlebuoy(:)
real, allocatable ::  qwleprec(:)
real, allocatable ::  momleadv(:,:)
real, allocatable ::  momlepress(:,:)
real, allocatable ::  momlebuoy(:,:)
real, allocatable ::  momlediff(:,:)
real, allocatable ::  tadv(:)
real, allocatable ::  tdiff(:)
real, allocatable ::  tlat(:)
real, allocatable ::  tlatqi(:)
real, allocatable ::  qifall(:)
real, allocatable ::  qpfall(:)
real, allocatable ::  tdiff_xy(:)
real, allocatable ::  tdiff_z(:)
real, allocatable ::  ttest0(:)
real, allocatable ::  ttest1(:)
real, allocatable ::  ttest2(:,:)

allocate ( twle(nz))
allocate ( twsb(nz))
allocate ( precflux(nz))
allocate ( uwle(nz))
allocate ( uwsb(nz))
allocate ( vwle(nz))
allocate ( vwsb(nz))
allocate ( radlwup(nz))
allocate ( radlwdn(nz))
allocate ( radswup(nz))
allocate ( radswdn(nz))
allocate ( radqrlw(nz))
allocate ( radqrsw(nz))
allocate ( tkeleadv(nz))
allocate ( tkelepress(nz))
allocate ( tkelediss(nz))
allocate ( tkelediff(nz))
allocate (tkelebuoy(nz))
allocate ( t2leadv(nz))
allocate (t2legrad(nz))
allocate (t2lediff(nz))
allocate (t2leprec(nz))
allocate (t2lediss(nz))
allocate ( q2leadv(nz))
allocate (q2legrad(nz))
allocate (q2lediff(nz))
allocate (q2leprec(nz))
allocate (q2lediss(nz))
allocate ( twleadv(nz))
allocate (twlediff(nz))
allocate (twlepres(nz))
allocate (twlebuoy(nz))
allocate (twleprec(nz))
allocate ( qwleadv(nz))
allocate (qwlediff(nz))
allocate (qwlepres(nz))
allocate (qwlebuoy(nz))
allocate (qwleprec(nz))
allocate (  momleadv(nz,3))
allocate (momlepress(nz,3))
allocate (momlebuoy(nz,3))
allocate (  momlediff(nz,3))
allocate (tadv(nz))
allocate (tdiff(nz))
allocate (tlat(nz))
allocate ( tlatqi(nz))
allocate (qifall(nz))
allocate (qpfall(nz))
allocate (tdiff_xy(nz))
allocate ( tdiff_z(nz))
allocate ( ttest0(nzm))
allocate ( ttest1(nz))
allocate ( ttest2(nz, 10))  !+++mhwang test


! register functions:


real, external :: esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
real, external :: qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm
integer, external :: lenstr, bytes_in_rec

! energy conservation diagnostics:
 
  real(kind=selected_real_kind(12)) total_water_before, total_water_after
  real(kind=selected_real_kind(12)) total_water_evap, total_water_prec, total_water_ls
!#ifdef CLUBB_CRM
  real(kind=selected_real_kind(12)) total_water_clubb
  real(kind=selected_real_kind(12)) total_energy_before, total_energy_after
  real(kind=selected_real_kind(12)) total_energy_evap, total_energy_prec, total_energy_ls
  real(kind=selected_real_kind(12)) total_energy_clubb, total_energy_rad
!#endif
  real(kind=selected_real_kind(12))  qtotmicro(5)  ! total water for water conservation test in microphysics +++mhwang

!===========================================================================
! UW ADDITIONS

! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
     icondavg_satdn, icondavg_satup, icondavg_env
real, allocatable :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
real, allocatable :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
character(LEN=8), allocatable :: condavgname(:) ! array of short names
character(LEN=25), allocatable :: condavglongname(:) ! array of long names

real, allocatable :: qlsvadv(:)
real, allocatable :: tlsvadv(:)
real, allocatable :: ulsvadv(:)
real, allocatable :: vlsvadv(:)

real, allocatable :: qnudge(:)
real, allocatable :: tnudge(:)
real, allocatable :: unudge(:)
real, allocatable :: vnudge(:)

real, allocatable :: qstor(:)
real, allocatable :: tstor(:)
real, allocatable :: ustor(:)
real, allocatable :: vstor(:)
real, allocatable :: qtostor(:)

real, allocatable :: utendcor(:)
real, allocatable :: vtendcor(:)

real, allocatable :: CF3D(:,:,:)

allocate (  qlsvadv(nzm)) ! Large-scale vertical advection tendency for total water
allocate (  tlsvadv(nzm)) ! Large-scale vertical advection tendency for temperature
allocate (  ulsvadv(nzm)) ! Large-scale vertical advection tendency for zonal velocity
allocate (  vlsvadv(nzm)) ! Large-scale vertical advection tendency for meridional velocity

allocate (  qnudge(nzm)) ! Nudging of horiz.-averaged total water profile
allocate (  tnudge(nzm)) ! Nudging of horiz.-averaged temperature profile
allocate (  unudge(nzm)) ! Nudging of horiz.-averaged zonal velocity
allocate (  vnudge(nzm)) ! Nudging of horiz.-averaged meridional velocity

allocate (  qstor(nzm)) ! Storage of horiz.-averaged total water profile
allocate (  tstor(nzm)) ! Storage of horiz.-averaged temperature profile
allocate (  ustor(nzm)) ! Storage of horiz.-averaged zonal velocity
allocate (  vstor(nzm)) ! Storage of horiz.-averaged meridional velocity
allocate (  qtostor(nzm)) ! Storage of horiz.-averaged total water profile (vapor + liquid)

allocate (  utendcor(nzm)) ! coriolis acceleration of zonal velocity
allocate (  vtendcor(nzm)) ! coriolis acceleration of meridional velocity

allocate (  CF3D(1:nx, 1:ny, 1:nzm))  ! Cloud fraction 
                                ! =1.0 when there is no fractional cloudiness scheme
                                ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

! 850 mbar horizontal winds
real, allocatable :: u850_xy(:,:)
real, allocatable :: v850_xy(:,:)
allocate ( u850_xy(nx,ny)) ! zonal velocity at 850 mb
allocate ( v850_xy(nx,ny)) ! meridional velocity at 850 mb

! Surface pressure
real, allocatable :: psfc_xy(:,:)
allocate ( psfc_xy(nx,ny)) ! pressure (in millibar) at lowest grid point

! Saturated water vapor path, useful for computing column relative humidity
real, allocatable :: swvp_xy(:,:)
allocate ( swvp_xy(nx,ny))  ! saturated water vapor path (wrt water)

! Cloud and echo top heights, and cloud top temperature (instantaneous)
real, allocatable :: cloudtopheight(:,:)
real, allocatable :: echotopheight(:,:)
real, allocatable :: cloudtoptemp(:,:)
allocate ( cloudtopheight(nx,ny))
allocate (  echotopheight(nx,ny))
allocate (   cloudtoptemp(nx,ny))

! END UW ADDITIONS
!===========================================================================
! Initial bubble parameters. Activated when perturb_type = 2
  real bubble_x0  
  real bubble_y0 
  real bubble_z0 
  real bubble_radius_hor 
  real bubble_radius_ver 
  real bubble_dtemp 
  real bubble_dq 
  real, allocatable ::  naer(:,:)     ! Aerosol number concentration [/m3]
  real, allocatable ::  vaer(:,:)     ! aerosol volume concentration [m3/m3]
  real, allocatable ::  hgaer(:,:)    ! hygroscopicity of aerosol mode

end module crmx_vars
