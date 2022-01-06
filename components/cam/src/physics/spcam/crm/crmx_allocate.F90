
subroutine crm_allocate
	
use crmx_vars

implicit none
 if (.not. allocated(z   )) allocate ( z(nz)      )! height of the pressure levels above surface,m
 if (.not. allocated(pres)) allocate ( pres(nzm)  )! pressure,mb at scalar levels
 if (.not. allocated(zi  )) allocate ( zi(nz)     )! height of the interface levels
 if (.not. allocated(presi)) allocate ( presi(nz)  )! pressure,mb at interface levels
 if (.not. allocated(adz)) allocate ( adz(nzm)   )! ratio of the thickness of scalar levels to dz 
 if (.not. allocated(adzw)) allocate ( adzw(nz)   )! ratio of the thinckness of w levels to dz
z=0.
pres=0.
zi=0.
presi=0.
adz=0.
adzw=0.
	
 if (.not. allocated(u)) allocate (u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)) ! x-wind
 if (.not. allocated(v)) allocate (v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)) ! y-wind
 if (.not. allocated(w)) allocate (w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )) ! z-wind
 if (.not. allocated(t)) allocate (t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)) ! liquid/ice water static energy 
 if (.not. allocated(micro_field)) allocate (micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields))
 if (.not. allocated(qn)) allocate( qn(nx,ny,nzm))
 if (.not. allocated(tke)) allocate (tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )   ! SGS TKE 
 if (.not. allocated(tk)) allocate (tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ) ! SGS eddy viscosity 
 if (.not. allocated(tkh)) allocate (tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) )! SGS eddy conductivity
 if (.not. allocated(def2)) allocate (def2(nx,ny,nzm))
u=0.
v=0.
w=0.
t=0.
micro_field=0.
qn=0.
tke=0.
tk=0.
tkh=0.
def2=0.
	
 if (.not. allocated(p)) allocate (p      (0:nx, (1-YES3D):ny, nzm))     ! perturbation pressure (from Poison eq)
 if (.not. allocated(tabs)) allocate (tabs   (nx, ny, nzm)            )     ! temperature
 if (.not. allocated(qv)) allocate (qv      (nx, ny, nzm)           )     ! water vapor
 if (.not. allocated(qcl)) allocate (qcl     (nx, ny, nzm)           )     ! liquid water  (condensate)
 if (.not. allocated(qpl)) allocate (qpl     (nx, ny, nzm)           )     ! liquid water  (precipitation)
 if (.not. allocated(qci)) allocate (qci     (nx, ny, nzm)           )     ! ice water  (condensate)
 if (.not. allocated(qpi)) allocate (qpi     (nx, ny, nzm)           )     ! ice water  (precipitation)
 if (.not. allocated(qp)) allocate (qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
 if (.not. allocated(f)) allocate (f (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
p=0.
tabs=0.
qv=0.
qcl=0.
qpl=0.
qci=0.
qpi=0.
qp=0.
f=0.
 if (.not. allocated(tke2)) allocate (tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))   ! SGS TKE
 if (.not. allocated(tk2)) allocate (tk2     (0:nxp1, (1-YES3D):nyp1, nzm)) ! SGS eddyviscosity
tke2=0.
tk2=0.
 if (.not. allocated(dudt)) allocate ( dudt   (nxp1, ny, nzm, 3))
 if (.not. allocated(dvdt)) allocate ( dvdt   (nx, nyp1, nzm, 3))
 if (.not. allocated(dwdt)) allocate ( dwdt   (nx, ny, nz,  3))
dudt=0.
dvdt=0.
dwdt=0.

 if (.not. allocated(misc)) allocate ( misc(nx, ny, nz))
 if (.not. allocated(fluxb)) allocate ( fluxb  (nx, ny))
 if (.not. allocated(fluxt)) allocate ( fluxt  (nx, ny))
 if (.not. allocated(fluxbu)) allocate ( fluxbu (nx, ny))
 if (.not. allocated(fluxbv)) allocate ( fluxbv (nx, ny))
 if (.not. allocated(fluxbt)) allocate ( fluxbt (nx, ny))
 if (.not. allocated(fluxbq)) allocate ( fluxbq (nx, ny))
 if (.not. allocated(fluxtu)) allocate ( fluxtu (nx, ny))
 if (.not. allocated(fluxtv)) allocate ( fluxtv (nx, ny))
 if (.not. allocated(fluxtt)) allocate ( fluxtt (nx, ny))
 if (.not. allocated(fluxtq)) allocate ( fluxtq (nx, ny))
 if (.not. allocated(fzero)) allocate ( fzero  (nx, ny))
 if (.not. allocated(precsfc)) allocate ( precsfc(nx,ny) ) ! surface precip. rate
 if (.not. allocated(precssfc)) allocate ( precssfc(nx,ny)) ! surface ice precip. rate
 if (.not. allocated(fluxbtmp)) allocate ( fluxbtmp(nx,ny))
 if (.not. allocated(fluxttmp)) allocate ( fluxttmp(nx,ny))
misc=0.
fluxb=0.
fluxt=0.
fluxbu=0.
fluxbv=0.
fluxbt=0.
fluxbq=0.
fluxtu=0.
fluxtv=0.
fluxtt=0.
fluxtq=0.
fzero=0.
precsfc=0.
precssfc=0.
fluxbtmp=0.
fluxttmp=0.

 if (.not. allocated(t0)) allocate ( t0(nzm))
 if (.not. allocated(q0)) allocate ( q0(nzm))
 if (.not. allocated(qv0)) allocate ( qv0(nzm))
 if (.not. allocated(tabs0)) allocate ( tabs0(nzm))
 if (.not. allocated(tl0)) allocate ( tl0(nzm))
 if (.not. allocated(tv0)) allocate ( tv0(nzm))
 if (.not. allocated(u0)) allocate ( u0(nzm))
 if (.not. allocated(v0)) allocate ( v0(nzm))
 if (.not. allocated(tg0)) allocate ( tg0(nzm))
 if (.not. allocated(qg0)) allocate ( qg0(nzm))
 if (.not. allocated(ug0)) allocate ( ug0(nzm))
 if (.not. allocated(vg0)) allocate ( vg0(nzm))
 if (.not. allocated(p0)) allocate ( p0(nzm))
 if (.not. allocated(tke0)) allocate (tke0(nzm))
 if (.not. allocated(t01)) allocate ( t01(nzm))
 if (.not. allocated(q01)) allocate ( q01(nzm))
 if (.not. allocated(qp0)) allocate ( qp0(nzm))
 if (.not. allocated(qn0)) allocate ( qn0(nzm))
t0=0.
q0=0.
qv0=0.
tabs0=0.
tl0=0.
tv0=0.
u0=0.
v0=0.
tg0=0.
qg0=0.
ug0=0.
vg0=0.
p0=0.
tke0=0.
t01=0.
q01=0.
qp0=0.
qn0=0.


 if (.not. allocated(prespot)) allocate ( prespot(nzm))  ! (1000./pres)**R/cp
 if (.not. allocated(rho)) allocate (  rho(nzm))    ! air density at pressure levels,kg/m3 
 if (.not. allocated(rhow)) allocate (  rhow(nz))   ! air density at vertical velocity levels,kg/m3
 if (.not. allocated(bet)) allocate (  bet(nzm))    ! = ggr/tv0
 if (.not. allocated(gamaz)) allocate (  gamaz(nzm)) ! ggr/cp*z
 if (.not. allocated(wsub)) allocate (  wsub(nz))   ! Large-scale subsidence velocity,m/s
 if (.not. allocated(qtend)) allocate (  qtend(nzm)) ! Large-scale tendency for total water
 if (.not. allocated(ttend)) allocate (  ttend(nzm)) ! Large-scale tendency for temp.
 if (.not. allocated(utend)) allocate (  utend(nzm)) ! Large-scale tendency for u
 if (.not. allocated(vtend)) allocate (  vtend(nzm)) ! Large-scale tendency for v
prespot=0.
rho=0.
rhow=0.
bet=0.
gamaz=0.
wsub=0.
qtend=0.
ttend=0.
utend=0.
vtend=0.

 if (.not. allocated(sstxy)) allocate ( sstxy(0:nx,(1-YES3D):ny))
 if (.not. allocated(fcory)) allocate ( fcory(0:ny))
 if (.not. allocated(fcorzy)) allocate ( fcorzy(ny) )
 if (.not. allocated(latitude)) allocate ( latitude(nx,ny))
 if (.not. allocated(longitude)) allocate ( longitude(nx,ny)) 
 if (.not. allocated(prec_xy)) allocate ( prec_xy(nx,ny))
 if (.not. allocated(shf_xy)) allocate ( shf_xy(nx,ny) )
 if (.not. allocated(lhf_xy)) allocate ( lhf_xy(nx,ny))
 if (.not. allocated(lwns_xy)) allocate ( lwns_xy(nx,ny))
 if (.not. allocated(swns_xy)) allocate ( swns_xy(nx,ny))
 if (.not. allocated(lwnsc_xy)) allocate ( lwnsc_xy(nx,ny))
 if (.not. allocated(swnsc_xy)) allocate ( swnsc_xy(nx,ny))
 if (.not. allocated(lwnt_xy)) allocate ( lwnt_xy(nx,ny))
 if (.not. allocated(swnt_xy)) allocate ( swnt_xy(nx,ny))
 if (.not. allocated(lwntc_xy)) allocate ( lwntc_xy(nx,ny))
 if (.not. allocated(swntc_xy)) allocate ( swntc_xy(nx,ny))
 if (.not. allocated(solin_xy)) allocate ( solin_xy(nx,ny))
 if (.not. allocated(pw_xy)) allocate ( pw_xy(nx,ny))
 if (.not. allocated(cw_xy)) allocate ( cw_xy(nx,ny))
 if (.not. allocated(iw_xy)) allocate ( iw_xy(nx,ny))
 if (.not. allocated(cld_xy)) allocate ( cld_xy(nx,ny))
 if (.not. allocated(u200_xy)) allocate ( u200_xy(nx,ny)) 
 if (.not. allocated(usfc_xy)) allocate ( usfc_xy(nx,ny))
 if (.not. allocated(v200_xy)) allocate ( v200_xy(nx,ny))
 if (.not. allocated(vsfc_xy)) allocate ( vsfc_xy(nx,ny))
 if (.not. allocated(w500_xy)) allocate ( w500_xy(nx,ny))
 if (.not. allocated(qocean_xy)) allocate ( qocean_xy(nx,ny))
sstxy=0.
fcory=0.
fcorzy=0.
latitude=0.
longitude=0.
prec_xy=0.
shf_xy=0.
lhf_xy=0.
lwns_xy=0.
swns_xy=0.
lwnsc_xy=0.
swnsc_xy=0.
lwnt_xy=0.
swnt_xy=0.
lwntc_xy=0.
swntc_xy=0.
solin_xy=0.
pw_xy=0.
cw_xy=0.
iw_xy=0.
cld_xy=0.
u200_xy=0.
usfc_xy=0.
v200_xy=0.
vsfc_xy=0.
w500_xy=0.
qocean_xy=0.

 if (.not. allocated(twle)) allocate ( twle(nz))
 if (.not. allocated(twsb)) allocate ( twsb(nz))
 if (.not. allocated(precflux)) allocate ( precflux(nz))
 if (.not. allocated(uwle)) allocate ( uwle(nz))
 if (.not. allocated(uwsb)) allocate ( uwsb(nz))
 if (.not. allocated(vwle)) allocate ( vwle(nz))
 if (.not. allocated(vwsb)) allocate ( vwsb(nz))
 if (.not. allocated(radlwup)) allocate ( radlwup(nz))
 if (.not. allocated(radlwdn)) allocate ( radlwdn(nz))
 if (.not. allocated(radswup)) allocate ( radswup(nz))
 if (.not. allocated(radswdn)) allocate ( radswdn(nz))
 if (.not. allocated(radqrlw)) allocate ( radqrlw(nz))
 if (.not. allocated(radqrsw)) allocate ( radqrsw(nz))
 if (.not. allocated(tkeleadv)) allocate ( tkeleadv(nz))
 if (.not. allocated(tkelepress)) allocate ( tkelepress(nz))
 if (.not. allocated(tkelediss)) allocate ( tkelediss(nz))
 if (.not. allocated(tkelediff)) allocate ( tkelediff(nz))
 if (.not. allocated(tkelebuoy)) allocate (tkelebuoy(nz))
 if (.not. allocated(t2leadv)) allocate ( t2leadv(nz))
 if (.not. allocated(t2legrad)) allocate (t2legrad(nz))
 if (.not. allocated(t2lediff)) allocate (t2lediff(nz))
 if (.not. allocated(t2leprec)) allocate (t2leprec(nz))
 if (.not. allocated(t2lediss)) allocate (t2lediss(nz))
twle=0.
twsb=0.
precflux=0.
uwle=0.
uwsb=0.
vwle=0.
vwsb=0.
radlwup=0.
radlwdn=0.
radswup=0.
radswdn=0.
radqrlw=0.
radqrsw=0.
tkeleadv=0.
tkelepress=0.
tkelediss=0.
tkelediff=0.
tkelebuoy=0.
t2leadv=0.
t2legrad=0.
t2lediff=0.
t2leprec=0.
t2lediss=0.

 if (.not. allocated(q2leadv)) allocate ( q2leadv(nz))
 if (.not. allocated(q2legrad)) allocate (q2legrad(nz))
 if (.not. allocated(q2lediff)) allocate (q2lediff(nz))
 if (.not. allocated(q2leprec)) allocate (q2leprec(nz))
 if (.not. allocated(q2lediss)) allocate (q2lediss(nz))
 if (.not. allocated(twleadv)) allocate ( twleadv(nz))
 if (.not. allocated(twlediff)) allocate (twlediff(nz))
 if (.not. allocated(twlepres)) allocate (twlepres(nz))
 if (.not. allocated(twlebuoy)) allocate (twlebuoy(nz))
 if (.not. allocated(twleprec)) allocate (twleprec(nz))
 if (.not. allocated(qwleadv)) allocate ( qwleadv(nz))
 if (.not. allocated(qwlediff)) allocate (qwlediff(nz))
 if (.not. allocated(qwlepres)) allocate (qwlepres(nz))
 if (.not. allocated(qwlebuoy)) allocate (qwlebuoy(nz))
 if (.not. allocated(qwleprec)) allocate (qwleprec(nz))
 if (.not. allocated(momleadv)) allocate (  momleadv(nz,3))
 if (.not. allocated(momlepress)) allocate (momlepress(nz,3))
 if (.not. allocated(momlebuoy)) allocate (momlebuoy(nz,3))
 if (.not. allocated(momlediff)) allocate (  momlediff(nz,3))
 if (.not. allocated(tadv)) allocate (tadv(nz))
 if (.not. allocated(tdiff)) allocate (tdiff(nz))
 if (.not. allocated(tlat)) allocate (tlat(nz))
 if (.not. allocated(tlatqi)) allocate ( tlatqi(nz))
 if (.not. allocated(qifall)) allocate (qifall(nz))
 if (.not. allocated(qpfall)) allocate (qpfall(nz))
 if (.not. allocated(tdiff_xy)) allocate (tdiff_xy(nz))
 if (.not. allocated(tdiff_z)) allocate ( tdiff_z(nz))
 if (.not. allocated(ttest0)) allocate ( ttest0(nzm))
 if (.not. allocated(ttest1)) allocate ( ttest1(nz))
 if (.not. allocated(ttest2)) allocate ( ttest2(nz, 10))  !+++mhwang test
q2leadv=0.
q2legrad=0.
q2lediff=0.
q2leprec=0.
q2lediss=0.
twleadv=0.
twlediff=0.
twlepres=0.
twlebuoy=0.
twleprec=0.
qwleadv=0.
qwlediff=0.
qwlepres=0.
qwlebuoy=0.
qwleprec=0.
momleadv=0.
momlepress=0.
momlebuoy=0.
momlediff=0.
tadv=0.
tdiff=0.
tlat=0.
tlatqi=0.
qifall=0.
qpfall=0.
tdiff_xy=0.
tdiff_z=0.
ttest0=0.
ttest1=0.
ttest2=0.


 if (.not. allocated(qlsvadv)) allocate (  qlsvadv(nzm)) ! Large-scale vertical advection tendency for total water
 if (.not. allocated(tlsvadv)) allocate (  tlsvadv(nzm)) ! Large-scale vertical advection tendency for temperature
 if (.not. allocated(ulsvadv)) allocate (  ulsvadv(nzm)) ! Large-scale vertical advection tendency for zonal velocity
 if (.not. allocated(vlsvadv)) allocate (  vlsvadv(nzm)) ! Large-scale vertical advection tendency for meridional velocity
qlsvadv=0.
tlsvadv=0.
ulsvadv=0.
vlsvadv=0.

 if (.not. allocated(qnudge)) allocate (  qnudge(nzm)) ! Nudging of horiz.-averaged total water profile
 if (.not. allocated(tnudge)) allocate (  tnudge(nzm)) ! Nudging of horiz.-averaged temperature profile
 if (.not. allocated(unudge)) allocate (  unudge(nzm)) ! Nudging of horiz.-averaged zonal velocity
 if (.not. allocated(vnudge)) allocate (  vnudge(nzm)) ! Nudging of horiz.-averaged meridional velocity
qnudge=0.
tnudge=0.
unudge=0.
vnudge=0.
 if (.not. allocated(qstor)) allocate (  qstor(nzm)) ! Storage of horiz.-averaged total water profile
 if (.not. allocated(tstor)) allocate (  tstor(nzm)) ! Storage of horiz.-averaged temperature profile
 if (.not. allocated(ustor)) allocate (  ustor(nzm)) ! Storage of horiz.-averaged zonal velocity
 if (.not. allocated(vstor)) allocate (  vstor(nzm)) ! Storage of horiz.-averaged meridional velocity
 if (.not. allocated(qtostor)) allocate (  qtostor(nzm)) ! Storage of horiz.-averaged total water profile (vapor + liquid)
qstor=0.
tstor=0.
ustor=0.
vstor=0.
qtostor=0.

 if (.not. allocated(utendcor)) allocate (  utendcor(nzm)) ! coriolis acceleration of zonal velocity
 if (.not. allocated(vtendcor)) allocate (  vtendcor(nzm)) ! coriolis acceleration of meridional velocity
utendcor=0.
vtendcor=0.
 if (.not. allocated(CF3D)) allocate (  CF3D(1:nx, 1:ny, 1:nzm))  ! Cloud fraction 
                                 ! =1.0 when there is no fractional cloudiness
                                 ! scheme
                                 ! = cloud fraction produced by fractioal
                                 ! cloudiness scheme when 

CF3D=0.
 if (.not. allocated(u850_xy)) allocate ( u850_xy(nx,ny)) ! zonal velocity at 850 mb
 if (.not. allocated(v850_xy)) allocate ( v850_xy(nx,ny)) ! meridional velocity at 850 mb

 if (.not. allocated(psfc_xy)) allocate ( psfc_xy(nx,ny)) ! pressure (in millibar) at lowest grid point

 if (.not. allocated(swvp_xy)) allocate ( swvp_xy(nx,ny)) ! saturated water vapor path (wrt water)

 if (.not. allocated(cloudtopheight)) allocate ( cloudtopheight(nx,ny))
 if (.not. allocated(echotopheight)) allocate (  echotopheight(nx,ny))
 if (.not. allocated(cloudtoptemp)) allocate (   cloudtoptemp(nx,ny))

u850_xy=0.
v850_xy=0.
psfc_xy=0.
swvp_xy=0.
cloudtopheight=0.
echotopheight=0.
cloudtoptemp=0.



end
