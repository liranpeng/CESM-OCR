
subroutine crm_deallocate
	
use crmx_vars

implicit none
 deallocate ( z      )! height of the pressure levels above surface,m
 deallocate ( pres  )! pressure,mb at scalar levels
 deallocate ( zi     )! height of the interface levels
 deallocate ( presi  )! pressure,mb at interface levels
 deallocate ( adz   )! ratio of the thickness of scalar levels to dz 
 deallocate ( adzw  )! ratio of the thinckness of w levels to dz
	
 deallocate (u ) ! x-wind
 deallocate (v ) ! y-wind
 deallocate (w ) ! z-wind
 deallocate (t ) ! liquid/ice water static energy 
 deallocate (micro_field)
 deallocate( qn)
 deallocate (tke)   ! SGS TKE 
 deallocate (tk ) ! SGS eddy viscosity 
 deallocate (tkh)! SGS eddy conductivity
	
 deallocate (p   )     ! perturbation pressure (from Poison eq)
 deallocate (tabs)     ! temperature
 deallocate (qv  )     ! water vapor
 deallocate (qcl )     ! liquid water  (condensate)
 deallocate (qpl )     ! liquid water  (precipitation)
 deallocate (qci )     ! ice water  (condensate)
 deallocate (qpi )     ! ice water  (precipitation)
 deallocate (qp)
 deallocate (f )


 deallocate (tke2)   ! SGS TKE
 deallocate (tk2 ) ! SGS eddyviscosity

 deallocate ( dudt)
 deallocate ( dvdt)
 deallocate ( dwdt)

 deallocate ( misc)

 deallocate ( fluxb)
 deallocate ( fluxt)
 deallocate ( fluxbu )
 deallocate ( fluxbv )
 deallocate ( fluxbt )
 deallocate ( fluxbq )
 deallocate ( fluxtu )
 deallocate ( fluxtv )
 deallocate ( fluxtt )
 deallocate ( fluxtq )
 deallocate ( fzero  )
 deallocate ( precsfc) ! surface precip. rate
 deallocate ( precssfc) ! surface ice precip. rate

 deallocate ( t0)
 deallocate ( q0)
 deallocate ( qv0)
 deallocate ( tabs0)
 deallocate ( tl0)
 deallocate ( tv0)
 deallocate ( u0)
 deallocate ( v0)
 deallocate ( tg0)
 deallocate ( qg0)
 deallocate ( ug0)
 deallocate ( vg0)
 deallocate ( p0)
 deallocate (tke0)
 deallocate ( t01)
 deallocate ( q01)
 deallocate ( qp0)
 deallocate ( qn0)

 deallocate ( prespot)  ! (1000./pres)**R/cp
 deallocate (  rho)    ! air density at pressure levels,kg/m3 
 deallocate (  rhow)   ! air density at vertical velocity levels,kg/m3
 deallocate (  bet)    ! = ggr/tv0
 deallocate (  gamaz) ! ggr/cp*z
 deallocate (  wsub)   ! Large-scale subsidence velocity,m/s
 deallocate (  qtend) ! Large-scale tendency for total water
 deallocate (  ttend) ! Large-scale tendency for temp.
 deallocate (  utend) ! Large-scale tendency for u
 deallocate (  vtend) ! Large-scale tendency for v

 deallocate ( sstxy)
 deallocate ( fcory)
 deallocate ( fcorzy )
 deallocate ( latitude)
 deallocate ( longitude) 
 deallocate ( prec_xy)
 deallocate ( shf_xy )
 deallocate ( lhf_xy)
 deallocate ( lwns_xy)
 deallocate ( swns_xy)
 deallocate ( lwnsc_xy)
 deallocate ( swnsc_xy)
 deallocate ( lwnt_xy)
 deallocate ( swnt_xy)
 deallocate ( lwntc_xy)
 deallocate ( swntc_xy)
 deallocate ( solin_xy)
 deallocate ( pw_xy)
 deallocate ( cw_xy)
 deallocate ( iw_xy)
 deallocate ( cld_xy)
 deallocate ( u200_xy) 
 deallocate ( usfc_xy)
 deallocate ( v200_xy)
 deallocate ( vsfc_xy)
 deallocate ( w500_xy)
 deallocate ( qocean_xy)

 deallocate ( twle)
 deallocate ( twsb)
 deallocate ( precflux)
 deallocate ( uwle)
 deallocate ( uwsb)
 deallocate ( vwle)
 deallocate ( vwsb)
 deallocate ( radlwup)
 deallocate ( radlwdn)
 deallocate ( radswup)
 deallocate ( radswdn)
 deallocate ( radqrlw)
 deallocate ( radqrsw)
 deallocate ( tkeleadv)
 deallocate ( tkelepress)
 deallocate ( tkelediss)
 deallocate ( tkelediff)
 deallocate (tkelebuoy)
 deallocate ( t2leadv)
 deallocate (t2legrad)
 deallocate (t2lediff)
 deallocate (t2leprec)
 deallocate (t2lediss)
 deallocate ( q2leadv)
 deallocate (q2legrad)
 deallocate (q2lediff)
 deallocate (q2leprec)
 deallocate (q2lediss)
 deallocate ( twleadv)
 deallocate (twlediff)
 deallocate (twlepres)
 deallocate (twlebuoy)
 deallocate (twleprec)
 deallocate ( qwleadv)
 deallocate (qwlediff)
 deallocate (qwlepres)
 deallocate (qwlebuoy)
 deallocate (qwleprec)
 deallocate (  momleadv)
 deallocate (momlepress)
 deallocate (momlebuoy)
 deallocate (  momlediff)
 deallocate (tadv)
 deallocate (tdiff)
 deallocate (tlat)
 deallocate ( tlatqi)
 deallocate (qifall)
 deallocate (qpfall)
 deallocate (tdiff_xy)
 deallocate ( tdiff_z)
 deallocate ( ttest0)
 deallocate ( ttest1)
 deallocate ( ttest2)  !+++mhwang test

 deallocate (  qlsvadv) ! Large-scale vertical advection tendency for total water
 deallocate (  tlsvadv) ! Large-scale vertical advection tendency for temperature
 deallocate (  ulsvadv) ! Large-scale vertical advection tendency for zonal velocity
 deallocate (  vlsvadv) ! Large-scale vertical advection tendency for meridional velocity

 deallocate (  qnudge) ! Nudging of horiz.-averaged total water profile
 deallocate (  tnudge) ! Nudging of horiz.-averaged temperature profile
 deallocate (  unudge) ! Nudging of horiz.-averaged zonal velocity
 deallocate (  vnudge) ! Nudging of horiz.-averaged meridional velocity

 deallocate (  qstor) ! Storage of horiz.-averaged total water profile
 deallocate (  tstor) ! Storage of horiz.-averaged temperature profile
 deallocate (  ustor) ! Storage of horiz.-averaged zonal velocity
 deallocate (  vstor) ! Storage of horiz.-averaged meridional velocity
 deallocate (  qtostor) ! Storage of horiz.-averaged total water profile (vapor + liquid)

 deallocate (  utendcor) ! coriolis acceleration of zonal velocity
 deallocate (  vtendcor) ! coriolis acceleration of meridional velocity

 deallocate (  CF3D)  ! Cloud fraction 
                                 ! =1.0 when there is no fractional cloudiness
                                 ! scheme
                                 ! = cloud fraction produced by fractioal
                                 ! cloudiness scheme when 


 deallocate ( u850_xy) ! zonal velocity at 850 mb
 deallocate ( v850_xy) ! meridional velocity at 850 mb

 deallocate ( psfc_xy) ! pressure (in millibar) at lowest grid point

 deallocate ( swvp_xy) ! saturated water vapor path (wrt water)

 deallocate ( cloudtopheight)
 deallocate (  echotopheight)
 deallocate (   cloudtoptemp)


end
