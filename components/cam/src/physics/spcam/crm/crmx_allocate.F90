
subroutine crm_allocate
	
use crmx_vars

implicit none
 allocate ( z(nz)      )! height of the pressure levels above surface,m
 allocate ( pres(nzm)  )! pressure,mb at scalar levels
 allocate ( zi(nz)     )! height of the interface levels
 allocate ( presi(nz)  )! pressure,mb at interface levels
 allocate ( adz(nzm)   )! ratio of the thickness of scalar levels to dz 
 allocate ( adzw(nz)   )! ratio of the thinckness of w levels to dz
z=0.
pres=0.
zi=0.
presi=0.
adz=0.
adzw=0.
	
 allocate (u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)) ! x-wind
 allocate (v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)) ! y-wind
 allocate (w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )) ! z-wind
 allocate (t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)) ! liquid/ice water static energy 
 allocate (micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields))
 allocate( qn(nx,ny,nzm))
 allocate (tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )   ! SGS TKE 
 allocate (tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ) ! SGS eddy viscosity 
 allocate (tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) )! SGS eddy conductivity
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
	
 allocate (p      (0:nx, (1-YES3D):ny, nzm))     ! perturbation pressure (from Poison eq)
 allocate (tabs   (nx, ny, nzm)            )     ! temperature
 allocate (qv      (nx, ny, nzm)           )     ! water vapor
 allocate (qcl     (nx, ny, nzm)           )     ! liquid water  (condensate)
 allocate (qpl     (nx, ny, nzm)           )     ! liquid water  (precipitation)
 allocate (qci     (nx, ny, nzm)           )     ! ice water  (condensate)
 allocate (qpi     (nx, ny, nzm)           )     ! ice water  (precipitation)
 allocate (qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
 allocate (f (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
p=0.
tabs=0.
qv=0.
qcl=0.
qpl=0.
qci=0.
qpi=0.
qp=0.
f=0.
 allocate (tke2(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))   ! SGS TKE
 allocate (tk2     (0:nxp1, (1-YES3D):nyp1, nzm)) ! SGS eddyviscosity
tke2=0.
tk2=0.
 allocate ( dudt   (nxp1, ny, nzm, 3))
 allocate ( dvdt   (nx, nyp1, nzm, 3))
 allocate ( dwdt   (nx, ny, nz,  3))
dudt=0.
dvdt=0.
dwdt=0.

 allocate ( misc(nx, ny, nz))
 allocate ( fluxb  (nx, ny))
 allocate ( fluxt  (nx, ny))
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


 allocate ( prespot(nzm))  ! (1000./pres)**R/cp
 allocate (  rho(nzm))    ! air density at pressure levels,kg/m3 
 allocate (  rhow(nz))   ! air density at vertical velocity levels,kg/m3
 allocate (  bet(nzm))    ! = ggr/tv0
 allocate (  gamaz(nzm)) ! ggr/cp*z
 allocate (  wsub(nz))   ! Large-scale subsidence velocity,m/s
 allocate (  qtend(nzm)) ! Large-scale tendency for total water
 allocate (  ttend(nzm)) ! Large-scale tendency for temp.
 allocate (  utend(nzm)) ! Large-scale tendency for u
 allocate (  vtend(nzm)) ! Large-scale tendency for v
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


 allocate (  qlsvadv(nzm)) ! Large-scale vertical advection tendency for total water
 allocate (  tlsvadv(nzm)) ! Large-scale vertical advection tendency for temperature
 allocate (  ulsvadv(nzm)) ! Large-scale vertical advection tendency for zonal velocity
 allocate (  vlsvadv(nzm)) ! Large-scale vertical advection tendency for meridional velocity
qlsvadv=0.
tlsvadv=0.
ulsvadv=0.
vlsvadv=0.

 allocate (  qnudge(nzm)) ! Nudging of horiz.-averaged total water profile
 allocate (  tnudge(nzm)) ! Nudging of horiz.-averaged temperature profile
 allocate (  unudge(nzm)) ! Nudging of horiz.-averaged zonal velocity
 allocate (  vnudge(nzm)) ! Nudging of horiz.-averaged meridional velocity
qnudge=0.
tnudge=0.
unudge=0.
vnudge=0.
 allocate (  qstor(nzm)) ! Storage of horiz.-averaged total water profile
 allocate (  tstor(nzm)) ! Storage of horiz.-averaged temperature profile
 allocate (  ustor(nzm)) ! Storage of horiz.-averaged zonal velocity
 allocate (  vstor(nzm)) ! Storage of horiz.-averaged meridional velocity
 allocate (  qtostor(nzm)) ! Storage of horiz.-averaged total water profile (vapor + liquid)
qstor=0.
tstor=0.
ustor=0.
vstor=0.
qtostor=0.

 allocate (  utendcor(nzm)) ! coriolis acceleration of zonal velocity
 allocate (  vtendcor(nzm)) ! coriolis acceleration of meridional velocity
utendcor=0.
vtendcor=0.
 allocate (  CF3D(1:nx, 1:ny, 1:nzm))  ! Cloud fraction 
                                 ! =1.0 when there is no fractional cloudiness
                                 ! scheme
                                 ! = cloud fraction produced by fractioal
                                 ! cloudiness scheme when 

CF3D=0.
 allocate ( u850_xy(nx,ny)) ! zonal velocity at 850 mb
 allocate ( v850_xy(nx,ny)) ! meridional velocity at 850 mb

 allocate ( psfc_xy(nx,ny)) ! pressure (in millibar) at lowest grid point

 allocate ( swvp_xy(nx,ny)) ! saturated water vapor path (wrt water)

 allocate ( cloudtopheight(nx,ny))
 allocate (  echotopheight(nx,ny))
 allocate (   cloudtoptemp(nx,ny))

u850_xy=0.
v850_xy=0.
psfc_xy=0.
swvp_xy=0.
cloudtopheight=0.
echotopheight=0.
cloudtoptemp=0.
end
