
subroutine crm_deallocate
	
use crmx_vars

implicit none
  if (allocated(z)) deallocate ( z      )! height of the pressure levels above surface,m
  if (allocated(pres)) deallocate ( pres  )! pressure,mb at scalar levels
  if (allocated(zi)) deallocate ( zi     )! height of the interface levels
  if (allocated(presi)) deallocate ( presi  )! pressure,mb at interface levels
  if (allocated(adz)) deallocate ( adz   )! ratio of the thickness of scalar levels to dz 
  if (allocated(adzw)) deallocate ( adzw  )! ratio of the thinckness of w levels to dz
	
  if (allocated(u)) deallocate (u ) ! x-wind
  if (allocated(v)) deallocate (v ) ! y-wind
  if (allocated(w)) deallocate (w ) ! z-wind
  if (allocated(t)) deallocate (t ) ! liquid/ice water static energy 
  if (allocated(micro_field)) deallocate (micro_field)
  if (allocated(qn)) deallocate( qn)
  if (allocated(tke)) deallocate (tke)   ! SGS TKE 
  if (allocated(tk)) deallocate (tk ) ! SGS eddy viscosity 
  if (allocated(tkh)) deallocate (tkh)! SGS eddy conductivity
	
  if (allocated(p)) deallocate (p   )     ! perturbation pressure (from Poison eq)
  if (allocated(tabs)) deallocate (tabs)     ! temperature
  if (allocated(qv)) deallocate (qv  )     ! water vapor
  if (allocated(qcl)) deallocate (qcl )     ! liquid water  (condensate)
  if (allocated(qpl)) deallocate (qpl )     ! liquid water  (precipitation)
  if (allocated(qci)) deallocate (qci )     ! ice water  (condensate)
  if (allocated(qpi)) deallocate (qpi )     ! ice water  (precipitation)
  if (allocated(qp)) deallocate (qp)
  if (allocated(f)) deallocate (f )


  if (allocated(tke2)) deallocate (tke2)   ! SGS TKE
  if (allocated(tk2)) deallocate (tk2 ) ! SGS eddyviscosity

  if (allocated(dudt)) deallocate ( dudt)
  if (allocated(dvdt)) deallocate ( dvdt)
  if (allocated(dwdt)) deallocate ( dwdt)

  if (allocated(misc)) deallocate ( misc)

  if (allocated(fluxb)) deallocate ( fluxb)
  if (allocated(fluxt)) deallocate ( fluxt)
  if (allocated(fluxbu)) deallocate ( fluxbu )
  if (allocated(fluxbv)) deallocate ( fluxbv )
  if (allocated(fluxbt)) deallocate ( fluxbt )
  if (allocated(fluxbq)) deallocate ( fluxbq )
  if (allocated(fluxtu)) deallocate ( fluxtu )
  if (allocated(fluxtv)) deallocate ( fluxtv )
  if (allocated(fluxtt)) deallocate ( fluxtt )
  if (allocated(fluxtq)) deallocate ( fluxtq )
  if (allocated(fzero)) deallocate ( fzero  )
  if (allocated(precsfc)) deallocate ( precsfc) ! surface precip. rate
  if (allocated(precssfc)) deallocate ( precssfc) ! surface ice precip. rate

  if (allocated(t0)) deallocate ( t0)
  if (allocated(q0)) deallocate ( q0)
  if (allocated(qv0)) deallocate ( qv0)
  if (allocated(tabs0)) deallocate ( tabs0)
  if (allocated(tl0)) deallocate ( tl0)
  if (allocated(tv0)) deallocate ( tv0)
  if (allocated(u0)) deallocate ( u0)
  if (allocated(v0)) deallocate ( v0)
  if (allocated(tg0)) deallocate ( tg0)
  if (allocated(qg0)) deallocate ( qg0)
  if (allocated(ug0)) deallocate ( ug0)
  if (allocated(vg0)) deallocate ( vg0)
  if (allocated(p0)) deallocate ( p0)
  if (allocated(tke0)) deallocate (tke0)
  if (allocated(t01)) deallocate ( t01)
  if (allocated(q01)) deallocate ( q01)
  if (allocated(qp0)) deallocate ( qp0)
  if (allocated(qn0)) deallocate ( qn0)

  if (allocated(prespot)) deallocate ( prespot)  ! (1000./pres)**R/cp
  if (allocated(rho)) deallocate (  rho)    ! air density at pressure levels,kg/m3 
  if (allocated(rhow)) deallocate (  rhow)   ! air density at vertical velocity levels,kg/m3
  if (allocated(bet)) deallocate (  bet)    ! = ggr/tv0
  if (allocated(gamaz)) deallocate (  gamaz) ! ggr/cp*z
  if (allocated(wsub)) deallocate (  wsub)   ! Large-scale subsidence velocity,m/s
  if (allocated(qtend)) deallocate (  qtend) ! Large-scale tendency for total water
  if (allocated(ttend)) deallocate (  ttend) ! Large-scale tendency for temp.
  if (allocated(utend)) deallocate (  utend) ! Large-scale tendency for u
  if (allocated(vtend)) deallocate (  vtend) ! Large-scale tendency for v

  if (allocated(sstxy)) deallocate ( sstxy)
  if (allocated(fcory)) deallocate ( fcory)
  if (allocated(fcorzy)) deallocate ( fcorzy )
  if (allocated(latitude)) deallocate ( latitude)
  if (allocated(longitude)) deallocate ( longitude) 
  if (allocated(prec_xy)) deallocate ( prec_xy)
  if (allocated(shf_xy)) deallocate ( shf_xy )
  if (allocated(lhf_xy)) deallocate ( lhf_xy)
  if (allocated(lwns_xy)) deallocate ( lwns_xy)
  if (allocated(swns_xy)) deallocate ( swns_xy)
  if (allocated(lwnsc_xy)) deallocate ( lwnsc_xy)
  if (allocated(swnsc_xy)) deallocate ( swnsc_xy)
  if (allocated(lwnt_xy)) deallocate ( lwnt_xy)
  if (allocated(swnt_xy)) deallocate ( swnt_xy)
  if (allocated(lwntc_xy)) deallocate ( lwntc_xy)
  if (allocated(swntc_xy)) deallocate ( swntc_xy)
  if (allocated(solin_xy)) deallocate ( solin_xy)
  if (allocated(pw_xy)) deallocate ( pw_xy)
  if (allocated(cw_xy)) deallocate ( cw_xy)
  if (allocated(iw_xy)) deallocate ( iw_xy)
  if (allocated(cld_xy)) deallocate ( cld_xy)
  if (allocated(u200_xy)) deallocate ( u200_xy) 
  if (allocated(usfc_xy)) deallocate ( usfc_xy)
  if (allocated(v200_xy)) deallocate ( v200_xy)
  if (allocated(vsfc_xy)) deallocate ( vsfc_xy)
  if (allocated(w500_xy)) deallocate ( w500_xy)
  if (allocated(qocean_xy)) deallocate ( qocean_xy)

  if (allocated(twle)) deallocate ( twle)
  if (allocated(twsb)) deallocate ( twsb)
  if (allocated(precflux)) deallocate ( precflux)
  if (allocated(uwle)) deallocate ( uwle)
  if (allocated(uwsb)) deallocate ( uwsb)
  if (allocated(vwle)) deallocate ( vwle)
  if (allocated(vwsb)) deallocate ( vwsb)
  if (allocated(radlwup)) deallocate ( radlwup)
  if (allocated(radlwdn)) deallocate ( radlwdn)
  if (allocated(radswup)) deallocate ( radswup)
  if (allocated(radswdn)) deallocate ( radswdn)
  if (allocated(radqrlw)) deallocate ( radqrlw)
  if (allocated(radqrsw)) deallocate ( radqrsw)
  if (allocated(tkeleadv)) deallocate ( tkeleadv)
  if (allocated(tkelepress)) deallocate ( tkelepress)
  if (allocated(tkelediss)) deallocate ( tkelediss)
  if (allocated(tkelediff)) deallocate ( tkelediff)
  if (allocated(tkelebuoy)) deallocate (tkelebuoy)
  if (allocated(t2leadv)) deallocate ( t2leadv)
  if (allocated(t2legrad)) deallocate (t2legrad)
  if (allocated(t2lediff)) deallocate (t2lediff)
  if (allocated(t2leprec)) deallocate (t2leprec)
  if (allocated(t2lediss)) deallocate (t2lediss)
  if (allocated(q2leadv)) deallocate ( q2leadv)
  if (allocated(q2legrad)) deallocate (q2legrad)
  if (allocated(q2lediff)) deallocate (q2lediff)
  if (allocated(q2leprec)) deallocate (q2leprec)
  if (allocated(q2lediss)) deallocate (q2lediss)
  if (allocated(twleadv)) deallocate ( twleadv)
  if (allocated(twlediff)) deallocate (twlediff)
  if (allocated(twlepres)) deallocate (twlepres)
  if (allocated(twlebuoy)) deallocate (twlebuoy)
  if (allocated(twleprec)) deallocate (twleprec)
  if (allocated(qwleadv)) deallocate ( qwleadv)
  if (allocated(qwlediff)) deallocate (qwlediff)
  if (allocated(qwlepres)) deallocate (qwlepres)
  if (allocated(qwlebuoy)) deallocate (qwlebuoy)
  if (allocated(qwleprec)) deallocate (qwleprec)
  if (allocated(momleadv)) deallocate (  momleadv)
  if (allocated(momlepress)) deallocate (momlepress)
  if (allocated(momlebuoy)) deallocate (momlebuoy)
  if (allocated(momlediff)) deallocate (  momlediff)
  if (allocated(tadv)) deallocate (tadv)
  if (allocated(tdiff)) deallocate (tdiff)
  if (allocated(tlat)) deallocate (tlat)
  if (allocated(tlatqi)) deallocate ( tlatqi)
  if (allocated(qifall)) deallocate (qifall)
  if (allocated(qpfall)) deallocate (qpfall)
  if (allocated(tdiff_xy)) deallocate (tdiff_xy)
  if (allocated(tdiff_z)) deallocate ( tdiff_z)
  if (allocated(ttest0)) deallocate ( ttest0)
  if (allocated(ttest1)) deallocate ( ttest1)
  if (allocated(ttest2)) deallocate ( ttest2)  !+++mhwang test

  if (allocated(qlsvadv)) deallocate (  qlsvadv) ! Large-scale vertical advection tendency for total water
  if (allocated(tlsvadv)) deallocate (  tlsvadv) ! Large-scale vertical advection tendency for temperature
  if (allocated(ulsvadv)) deallocate (  ulsvadv) ! Large-scale vertical advection tendency for zonal velocity
  if (allocated(vlsvadv)) deallocate (  vlsvadv) ! Large-scale vertical advection tendency for meridional velocity

  if (allocated(qnudge)) deallocate (  qnudge) ! Nudging of horiz.-averaged total water profile
  if (allocated(tnudge)) deallocate (  tnudge) ! Nudging of horiz.-averaged temperature profile
  if (allocated(unudge)) deallocate (  unudge) ! Nudging of horiz.-averaged zonal velocity
  if (allocated(vnudge)) deallocate (  vnudge) ! Nudging of horiz.-averaged meridional velocity

  if (allocated(qstor)) deallocate (  qstor) ! Storage of horiz.-averaged total water profile
  if (allocated(tstor)) deallocate (  tstor) ! Storage of horiz.-averaged temperature profile
  if (allocated(ustor)) deallocate (  ustor) ! Storage of horiz.-averaged zonal velocity
  if (allocated(vstor)) deallocate (  vstor) ! Storage of horiz.-averaged meridional velocity
  if (allocated(qtostor)) deallocate (  qtostor) ! Storage of horiz.-averaged total water profile (vapor + liquid)

  if (allocated(utendcor)) deallocate (  utendcor) ! coriolis acceleration of zonal velocity
  if (allocated(vtendcor)) deallocate (  vtendcor) ! coriolis acceleration of meridional velocity

  if (allocated(CF3D)) deallocate (  CF3D)  ! Cloud fraction 
                                 ! =1.0 when there is no fractional cloudiness
                                 ! scheme
                                 ! = cloud fraction produced by fractioal
                                 ! cloudiness scheme when 


  if (allocated(u850_xy)) deallocate ( u850_xy) ! zonal velocity at 850 mb
  if (allocated(v850_xy)) deallocate ( v850_xy) ! meridional velocity at 850 mb

  if (allocated(psfc_xy)) deallocate ( psfc_xy) ! pressure (in millibar) at lowest grid point

  if (allocated(swvp_xy)) deallocate ( swvp_xy) ! saturated water vapor path (wrt water)

  if (allocated(cloudtopheight)) deallocate ( cloudtopheight)
  if (allocated(echotopheight)) deallocate (  echotopheight)
  if (allocated(cloudtoptemp)) deallocate (   cloudtoptemp)


end
