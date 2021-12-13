program TwoExecutableDriver
  use shr_kind_mod,      only: r8 => SHR_KIND_R8
  use shr_kind_mod,      only: i8 => SHR_KIND_I8
!  use ppgrid,           only: pver
!  use crmx_saturation, only:sat_mixrat_liq
  use crmdims,         only: crm_nx, crm_ny, crm_nz
  use crmx_crm_module_orc,     only: crm_orc
!  use ppgrid,          only: pcols
!  use spmd_utils, only: npes
  implicit none

  !--------------------------------------------------------------------------
  !bloss: MPI-related variables
  !--------------------------------------------------------------------------
  integer :: crm_comm,crm_comm_in
  integer :: myrank_global, numproc_global,crm_comm_color
  integer :: myrank_crm, numproc_crm,myrank_crm_in, numproc_crm_in
  integer :: ierr,status
  integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars
  integer, parameter :: pcols = 16
  integer, parameter :: pver = 26
  integer :: i,from,EndFlag
  integer, dimension(1) :: destGCM0 = -999
  CHARACTER(LEN=6) :: crm_number
  integer gcolindex(pcols)  ! array of global latitude indices
  ! local copies of input ingredients to CRM ( to be received)
  integer :: inp01_lchnk,inp02_i ! pritch note important for integers to
 ! be MPI_received that they masquerade as 1-element arrays here.
  real(r8), dimension(pver) :: inp03_tl,inp04_ql,inp05_qccl,inp06_qiil,inp07_ul,inp08_vl
  real(r8), dimension(pver) :: inp10_pmid,inp11_pdel,inp13_zm,inp14_zi
  real(r8) :: inp09_ps,inp12_phis,inp15_ztodt

! output ingrdients from CRM (to be sent)
  real(r8),dimension(pver) :: out01_qltend,out02_qcltend,out03_qiltend,out04_sltend
! inout ingredients:
  real(r8), dimension(crm_nx,crm_ny,crm_nz) :: outin01_crm_u,outin02_crm_v,outin03_crm_w,outin04_crm_t,&
        outin06_crm_qrad,outin07_qc_crm,outin08_qi_crm,outin09_qpc_crm,outin10_qpi_crm,outin12_t_rad,&
        outin13_qv_rad,outin14_qc_rad,outin15_qi_rad,outin16_cld_rad,outin17_cld3d_crm,crm_tk,crm_tkh,out_cld_rad
  real(r8), dimension(crm_nx,crm_ny) :: outin11_prec_crm
  real(r8), allocatable :: flattened_crm_inout(:),out_Var_Flat(:)
  real(r8),dimension(crm_nx, crm_ny, crm_nz,nmicro_fields+1) :: outin05_crm_micro
  real (r8) :: precc,precl,precsc,precsl,cltot,clhgh,clmed,cllow,lon,lat
  real (r8), dimension(pver) :: cld,cldtop,gicewp,gliqwp,mctot,mcup,mcdn,mcuup,mcudn
  real (r8), dimension(pver) :: spqc,spqi,spqs,spqg,spqr
  real (r8), dimension(pver) :: mu_crm,md_crm,du_crm,eu_crm,ed_crm,&
                tkez,tkesgsz,tk_crm,flux_u,flux_v,flux_qt,fluxsgs_qt,flux_qp
  real(r8), dimension(pver) :: precflux,qt_ls,qt_trans,qp_trans,qp_fall,qp_evp,qp_src,t_ls
  real(r8), dimension(20) :: qtotcrm
  real(r8) :: prectend,precstend,ocnfrac,wnd,tau00,bflx,fluxu0,fluxv0,fluxt0,fluxq0
  real(r8) :: taux_crm,tauy_crm,z0m,timing_factor,jt_crm,mx_crm
  integer :: chnksz,nflat,fcount,ii,jj,kk,ll
  integer,parameter :: structlen  = 49
  integer,parameter :: singlelen  = 30
  integer,parameter :: flen       = structlen*pver+singlelen+1+20+pcols
  integer,parameter :: flen2      = 17*crm_nx*crm_nz*crm_nz + crm_nx*crm_ny*crm_nz*nmicro_fields+crm_nx*crm_ny
  integer,parameter :: structleno = 37
  integer,parameter :: singleleno = 22
  integer,parameter :: rank_offset=2
  integer,parameter :: fleno      = structleno*pver+singleleno+1+20
  real(r8),dimension(flen ) :: inp_Var_Flat
  real(r8),dimension(flen2) :: inp_Var_Flat2

  include 'mpif.h'

  ! initialize MPI alongside CESM call in cime_comp_mod.F90
  call mpi_init(ierr)
  !call MPI_Init_thread(ierr)
  if(ierr.eq.0) then
!bloss     write(13,*) 'CRM: successful call to mpi_init'
  else
     STOP 'TwoExecutableDriver.F90: MPI_Init Failed'
  end if

  ! Split MPI_COMM_WORLD into two communicators:
  !  - global_comm: used by CIME for CESM-related communication
  !  - crm_comm: used by CRM routines to receive data from CESM
  !         for separate-executable CRM computations and
  !         to send the resulting tendencies back to CESM.
  call mpi_comm_split(MPI_COMM_WORLD, 1, 0, crm_comm, ierr)
  if(ierr.eq.0) then
!bloss     write(13,*) 'CRM: successful call to mpi_split'
  else
     STOP 'TwoExecutableDriver.F90: MPI_Split Failed'
  end if

  ! get information on MPI_COMM_WORLD
  call mpi_comm_size(MPI_COMM_WORLD, numproc_global, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)

  ! get information on MPI_COMM_WORLD
  call mpi_comm_size(crm_comm, numproc_crm, ierr)
  call mpi_comm_rank(crm_comm, myrank_crm, ierr)
 
 ! ---------- basic comm pool size sanity checks by bloss -------------
  923 format(I6.6)
  write(crm_number,923) myrank_crm
  open(unit=13,file='crm.log.'//TRIM(crm_number),form='formatted')
  ! ----------- GCM handshake from spcam_drivers --------------
  EndFlag  = 1
  do i = 1,numproc_crm-1
    if((myrank_global.lt.(50+i*rank_offset)).and.(myrank_global.ge.(50+(i-1)*rank_offset))) then
       crm_comm_color = int(myrank_crm/2)+2
       write(13,924) numproc_global, myrank_global,numproc_crm,myrank_crm,crm_comm_color
924        format('MPI_COMM_WORLD: size/myrank = ',2i5,', crm_comm: size/myrank= ',3i5)
       call mpi_comm_split(crm_comm, crm_comm_color, 0, crm_comm_in, ierr)
       call mpi_comm_size(crm_comm_in, numproc_crm_in, ierr)
       call mpi_comm_rank(crm_comm_in, myrank_crm_in, ierr)
!      open(unit=13,file='crm.log.'//TRIM(crm_number),form='formatted')    
       write(13,*) 'Liran Check:', myrank_crm, myrank_global,numproc_crm_in,myrank_crm_in
       call MPI_Barrier(crm_comm,ierr)
    end if
    if(myrank_global.eq.(50+(i-1)*rank_offset)) then
      ! Recieve rank of host GCM column linked to this CRM, for eventual MPI_Send
      call MPI_Recv(destGCM0,1,MPI_INTEGER,MPI_ANY_SOURCE,54321,MPI_COMM_WORLD,status,ierr)
      if (ierr.eq.0) then
        write(13,*) 'CRM rank',myrank_crm,crm_comm_in,' got handshake; its GCM dest rank=',destGCM0,myrank_global
        call mpi_comm_size(crm_comm_in, numproc_crm_in, ierr)
        call mpi_comm_rank(crm_comm_in, myrank_crm_in, ierr)    
      else 
        write (13,*) 'MPI_Recv from spcam_driver handshake failed for CRM rank ',myrank_crm,',ierr=',ierr
      end if
      EndFlag  = 0
      allocate(out_Var_Flat(fleno+nflat))
     end if
  end do
  EndFlag  = 0
  do while (EndFlag.eq.0)
  ! ----------- Receive from cesm.exe/crm_physics inputs to CRM right before call to crm() ----------------
  chnksz = crm_nx*crm_nz*crm_nz
  nflat =  17*chnksz + crm_nx*crm_ny*crm_nz*nmicro_fields +crm_nx*crm_ny
! =====================================================================================
! Receive from GCM, in order, the input ingredients expected by crm subroutine:
! =====================================================================================
  call MPI_Recv(inp_Var_Flat(:) ,flen   ,MPI_REAL8,MPI_ANY_SOURCE,9018,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp_Var_Flat2(:) ,flen2,MPI_REAL8,MPI_ANY_SOURCE,9019,MPI_COMM_WORLD,status,ierr)
  inp01_lchnk         = int(inp_Var_Flat(1))
  inp02_i             = int(inp_Var_Flat(2))
  inp09_ps            =     inp_Var_Flat(3)
  inp12_phis          =     inp_Var_Flat(4)
  inp15_ztodt         =     inp_Var_Flat(5) 
  precc               =     inp_Var_Flat(6)
  precl               =     inp_Var_Flat(7)
  precsc              =     inp_Var_Flat(8)
  precsl              =     inp_Var_Flat(9)
  cltot               =     inp_Var_Flat(10)
  clhgh               =     inp_Var_Flat(11)
  clmed               =     inp_Var_Flat(12)
  cllow               =     inp_Var_Flat(13)
  prectend            =     inp_Var_Flat(14)
  precstend           =     inp_Var_Flat(15)
  ocnfrac             =     inp_Var_Flat(16)
  wnd                 =     inp_Var_Flat(17)
  tau00               =     inp_Var_Flat(18)
  bflx                =     inp_Var_Flat(19)
  fluxu0              =     inp_Var_Flat(20)
  fluxv0              =     inp_Var_Flat(21)
  fluxt0              =     inp_Var_Flat(22)
  fluxq0              =     inp_Var_Flat(23)
  taux_crm            =     inp_Var_Flat(24)
  tauy_crm            =     inp_Var_Flat(25)
  z0m                 =     inp_Var_Flat(26)
  timing_factor       =     inp_Var_Flat(27)
  lat                 =     inp_Var_Flat(28)
  lon                 =     inp_Var_Flat(29)
  EndFlag             =     inp_Var_Flat(30)
  inp03_tl(1:pver)    = inp_Var_Flat(        1+singlelen:   pver+singlelen)
  inp04_ql(1:pver)    = inp_Var_Flat( 1*pver+1+singlelen: 2*pver+singlelen)
  inp05_qccl(1:pver)  = inp_Var_Flat( 2*pver+1+singlelen: 3*pver+singlelen)
  inp06_qiil(1:pver)  = inp_Var_Flat( 3*pver+1+singlelen: 4*pver+singlelen) 
  inp07_ul(1:pver)    = inp_Var_Flat( 4*pver+1+singlelen: 5*pver+singlelen)
  inp08_vl(1:pver)    = inp_Var_Flat( 5*pver+1+singlelen: 6*pver+singlelen)   
  inp10_pmid(1:pver)  = inp_Var_Flat( 6*pver+1+singlelen: 7*pver+singlelen)
  inp11_pdel(1:pver)  = inp_Var_Flat( 7*pver+1+singlelen: 8*pver+singlelen)
  inp13_zm(1:pver)    = inp_Var_Flat( 8*pver+1+singlelen: 9*pver+singlelen)
  inp14_zi(1:pver+1)  = inp_Var_Flat( 9*pver+1+singlelen:10*pver+singlelen+1)
  cld(1:pver)         = inp_Var_Flat(10*pver+2+singlelen:11*pver+singlelen+1)
  cldtop(1:pver)      = inp_Var_Flat(11*pver+2+singlelen:12*pver+singlelen+1)
  gicewp(1:pver)      = inp_Var_Flat(12*pver+2+singlelen:13*pver+singlelen+1)
  gliqwp(1:pver)      = inp_Var_Flat(13*pver+2+singlelen:14*pver+singlelen+1)
  mctot(1:pver)       = inp_Var_Flat(14*pver+2+singlelen:15*pver+singlelen+1)
  mcup(1:pver)        = inp_Var_Flat(15*pver+2+singlelen:16*pver+singlelen+1)
  mcdn(1:pver)        = inp_Var_Flat(16*pver+2+singlelen:17*pver+singlelen+1)
  mcuup(1:pver)       = inp_Var_Flat(17*pver+2+singlelen:18*pver+singlelen+1)
  mcudn(1:pver)       = inp_Var_Flat(18*pver+2+singlelen:19*pver+singlelen+1)
  spqc(1:pver)        = inp_Var_Flat(19*pver+2+singlelen:20*pver+singlelen+1)
  spqi(1:pver)        = inp_Var_Flat(20*pver+2+singlelen:21*pver+singlelen+1)
  spqs(1:pver)        = inp_Var_Flat(21*pver+2+singlelen:22*pver+singlelen+1)
  spqg(1:pver)        = inp_Var_Flat(22*pver+2+singlelen:23*pver+singlelen+1)
  spqr(1:pver)        = inp_Var_Flat(23*pver+2+singlelen:24*pver+singlelen+1)
  mu_crm(1:pver)      = inp_Var_Flat(24*pver+2+singlelen:25*pver+singlelen+1)
  md_crm(1:pver)      = inp_Var_Flat(25*pver+2+singlelen:26*pver+singlelen+1)
  du_crm(1:pver)      = inp_Var_Flat(26*pver+2+singlelen:27*pver+singlelen+1)
  eu_crm(1:pver)      = inp_Var_Flat(27*pver+2+singlelen:28*pver+singlelen+1)
  ed_crm(1:pver)      = inp_Var_Flat(28*pver+2+singlelen:29*pver+singlelen+1)
  tkez(1:pver)        = inp_Var_Flat(29*pver+2+singlelen:30*pver+singlelen+1)
  tkesgsz(1:pver)     = inp_Var_Flat(30*pver+2+singlelen:31*pver+singlelen+1)
  tk_crm(1:pver)      = inp_Var_Flat(31*pver+2+singlelen:32*pver+singlelen+1)
  flux_u(1:pver)      = inp_Var_Flat(32*pver+2+singlelen:33*pver+singlelen+1)
  flux_v(1:pver)      = inp_Var_Flat(33*pver+2+singlelen:34*pver+singlelen+1)
  flux_qt(1:pver)     = inp_Var_Flat(34*pver+2+singlelen:35*pver+singlelen+1)
  fluxsgs_qt(1:pver)  = inp_Var_Flat(35*pver+2+singlelen:36*pver+singlelen+1)
  flux_qp(1:pver)     = inp_Var_Flat(36*pver+2+singlelen:37*pver+singlelen+1)
  precflux(1:pver)    = inp_Var_Flat(37*pver+2+singlelen:38*pver+singlelen+1)
  qt_ls(1:pver)       = inp_Var_Flat(38*pver+2+singlelen:39*pver+singlelen+1)
  qt_trans(1:pver)    = inp_Var_Flat(39*pver+2+singlelen:40*pver+singlelen+1)
  qp_trans(1:pver)    = inp_Var_Flat(40*pver+2+singlelen:41*pver+singlelen+1)
  qp_fall(1:pver)     = inp_Var_Flat(41*pver+2+singlelen:42*pver+singlelen+1)
  qp_evp(1:pver)      = inp_Var_Flat(42*pver+2+singlelen:43*pver+singlelen+1)
  qp_src(1:pver)      = inp_Var_Flat(43*pver+2+singlelen:44*pver+singlelen+1)
  t_ls(1:pver)        = inp_Var_Flat(44*pver+2+singlelen:45*pver+singlelen+1)
  out01_qltend(1:pver)= inp_Var_Flat(45*pver+2+singlelen:46*pver+singlelen+1)
  out02_qcltend(1:pver)=inp_Var_Flat(46*pver+2+singlelen:47*pver+singlelen+1)
  out03_qiltend(1:pver)=inp_Var_Flat(47*pver+2+singlelen:48*pver+singlelen+1)
  out04_sltend(1:pver)= inp_Var_Flat(48*pver+2+singlelen:49*pver+singlelen+1)
  qtotcrm(1:20)       = inp_Var_Flat(49*pver+2+singlelen:49*pver+21+singlelen)
  gcolindex(1:pcols)  = inp_Var_Flat(49*pver+22+singlelen:49*pver+21+singlelen+pcols)
  fcount = 0
  do ii=1,crm_nx
    do jj=1,crm_ny
      do kk=1,crm_nz
        fcount = fcount + 1
        outin01_crm_u(ii,jj,kk)    = inp_Var_Flat2(fcount)
        outin02_crm_v(ii,jj,kk)    = inp_Var_Flat2(fcount + 1 * chnksz)
        outin03_crm_w(ii,jj,kk)    = inp_Var_Flat2(fcount + 2 * chnksz)
        outin04_crm_t(ii,jj,kk)    = inp_Var_Flat2(fcount + 3 * chnksz)
        outin06_crm_qrad(ii,jj,kk) = inp_Var_Flat2(fcount + 4 * chnksz)
        outin07_qc_crm(ii,jj,kk)   = inp_Var_Flat2(fcount + 5 * chnksz)
        outin08_qi_crm(ii,jj,kk)   = inp_Var_Flat2(fcount + 6 * chnksz)
        outin09_qpc_crm(ii,jj,kk)  = inp_Var_Flat2(fcount + 7 * chnksz)
        outin10_qpi_crm(ii,jj,kk)  = inp_Var_Flat2(fcount + 8 * chnksz)
        outin12_t_rad(ii,jj,kk)    = inp_Var_Flat2(fcount + 9 * chnksz)
        outin13_qv_rad(ii,jj,kk)   = inp_Var_Flat2(fcount + 10* chnksz)
        outin14_qc_rad(ii,jj,kk)   = inp_Var_Flat2(fcount + 11* chnksz)
        outin15_qi_rad(ii,jj,kk)   = inp_Var_Flat2(fcount + 12* chnksz)
        out_cld_rad(ii,jj,kk)      = inp_Var_Flat2(fcount + 13* chnksz)
        outin16_cld_rad(ii,jj,kk)  = inp_Var_Flat2(fcount + 14* chnksz)
        crm_tk(ii,jj,kk)           = inp_Var_Flat2(fcount + 15* chnksz)
        crm_tkh(ii,jj,kk)          = inp_Var_Flat2(fcount + 16* chnksz)
      end do
    end do
  end do

  fcount = 17*chnksz
  do ii=1,crm_nx
    do jj=1,crm_ny
      do kk=1,crm_nz
        do ll=1,nmicro_fields
          fcount=fcount+1
          outin05_crm_micro(ii,jj,kk,ll) = inp_Var_Flat2(fcount)
        end do
      end do
    end do
  end do

  fcount = 17*chnksz + crm_nx*crm_ny*crm_nz*nmicro_fields
  do ii=1,crm_nx
    do jj=1,crm_ny
      fcount = fcount + 1
      outin11_prec_crm(ii,jj) = inp_Var_Flat2(fcount)
    end do
  end do
! Preparing to call the CRM

  !call mpi_comm_size(crm_comm_in, numproc_crm_in, ierr)
  !call mpi_comm_rank(crm_comm_in, myrank_crm_in, ierr)


  do i = 1,numproc_crm-1
    if(myrank_global.eq.(50+(i-1)*rank_offset)) then
        call mpi_comm_size(crm_comm_in, numproc_crm_in, ierr)
        call mpi_comm_rank(crm_comm_in, myrank_crm_in, ierr)
    end if
  end do


  write(13,*) 'Liran Check again:', myrank_crm, myrank_global, numproc_crm_in,myrank_crm_in 
  call crm_orc (numproc_crm_in,myrank_crm_in,lon,lat,gcolindex,inp01_lchnk, inp02_i,                            &
            inp03_tl(:),inp04_ql(:),inp05_qccl(:),inp06_qiil(:), &
            inp07_ul(:),inp08_vl(:),inp09_ps,inp10_pmid(:),inp11_pdel(:), &
            inp12_phis,inp13_zm(:),inp14_zi(:),inp15_ztodt,pver, &
            out01_qltend(:),out02_qcltend(:),out03_qiltend(:),out04_sltend(:),&
            outin01_crm_u(:,:,:),outin02_crm_v(:,:,:),outin03_crm_w(:,:,:),&
            outin04_crm_t(:,:,:),outin05_crm_micro(:,:,:,:),outin06_crm_qrad(:,:,:),&
            outin07_qc_crm(:,:,:),outin08_qi_crm(:,:,:),outin09_qpc_crm(:,:,:),&
            outin10_qpi_crm(:,:,:),outin11_prec_crm(:,:),outin12_t_rad(:,:,:),&
            outin13_qv_rad(:,:,:),outin14_qc_rad(:,:,:),outin15_qi_rad(:,:,:),&
            outin16_cld_rad(:,:,:),outin17_cld3d_crm(:,:,:),&
#ifdef m2005 ! to be added.
#endif            
            precc,precl,precsc,precsl,&
            cltot,clhgh,clmed,&
            cllow,                cld(:),  cldtop(:), &
            gicewp(:),        gliqwp(:),&
            mctot(:),         mcup(:),            mcdn(:),&
            mcuup(:),         mcudn(:),           &
            spqc(:),          spqi(:),            spqs(:),&
            spqg(:),               spqr(:),            &
#ifdef m2005
#endif
#ifdef SPCAM_CLUBB_SGS
#endif
            crm_tk( :, :, :), crm_tkh( :, :, :),&
          mu_crm(:),        md_crm(:),          du_crm(:), eu_crm(:),&
             ed_crm(:),        jt_crm,            mx_crm,&
             tkez(:),          tkesgsz(:),         tk_crm( :),&
             flux_u(:),        flux_v(:),          flux_qt(:), fluxsgs_qt(:),         flux_qp(:),         &
             precflux(:),      qt_ls(:),           qt_trans(:), qp_trans(:),           qp_fall(:),         &
             qp_evp(:),        qp_src(:),          t_ls(:), prectend,             precstend,         &
             ocnfrac,  wnd,                  tau00, bflx,                                          &
             fluxu0,             fluxv0,               fluxt0, fluxq0,                                        &
             taux_crm,        tauy_crm,          z0m, timing_factor,        qtotcrm( :)         &
            )
! ====================== DONE CALLING CRM -- TIME TO SEND OUTPUTS TO GCM

  out_Var_Flat(        1)                               = precc
  out_Var_Flat(        2)                               = precl
  out_Var_Flat(        3)                               = precsc
  out_Var_Flat(        4)                               = precsl
  out_Var_Flat(        5)                               = cltot
  out_Var_Flat(        6)                               = clhgh
  out_Var_Flat(        7)                               = clmed
  out_Var_Flat(        8)                               = cllow
  out_Var_Flat(        9)                               = prectend
  out_Var_Flat(       10)                               = precstend
  out_Var_Flat(       11)                               = ocnfrac
  out_Var_Flat(       12)                               = wnd
  out_Var_Flat(       13)                               = tau00
  out_Var_Flat(       14)                               = bflx
  out_Var_Flat(       15)                               = fluxu0
  out_Var_Flat(       16)                               = fluxv0
  out_Var_Flat(       17)                               = fluxt0
  out_Var_Flat(       18)                               = fluxq0
  out_Var_Flat(       19)                               = taux_crm
  out_Var_Flat(       20)                               = tauy_crm
  out_Var_Flat(       21)                               = z0m
  out_Var_Flat(       22)                               = timing_factor
  out_Var_Flat(        1+singleleno:   pver+singleleno) = out01_qltend
  out_Var_Flat( 1*pver+1+singleleno: 2*pver+singleleno) = out02_qcltend
  out_Var_Flat( 2*pver+1+singleleno: 3*pver+singleleno) = out03_qiltend
  out_Var_Flat( 3*pver+1+singleleno: 4*pver+singleleno) = out04_sltend
  out_Var_Flat( 4*pver+1+singleleno: 5*pver+singleleno) = cld(:)
  out_Var_Flat( 5*pver+1+singleleno: 6*pver+singleleno) = cldtop(:)
  out_Var_Flat( 6*pver+1+singleleno: 7*pver+singleleno) = gicewp(:)
  out_Var_Flat( 7*pver+1+singleleno: 8*pver+singleleno) = gliqwp(:)
  out_Var_Flat( 8*pver+1+singleleno: 9*pver+singleleno) = mctot(:)
  out_Var_Flat( 9*pver+1+singleleno:10*pver+singleleno) = mcup(:)
  out_Var_Flat(10*pver+1+singleleno:11*pver+singleleno) = mcdn(:)
  out_Var_Flat(11*pver+1+singleleno:12*pver+singleleno) = mcuup(:)
  out_Var_Flat(12*pver+1+singleleno:13*pver+singleleno) = mcudn(:)
  out_Var_Flat(13*pver+1+singleleno:14*pver+singleleno) = spqc(:)
  out_Var_Flat(14*pver+1+singleleno:15*pver+singleleno) = spqi(:)
  out_Var_Flat(15*pver+1+singleleno:16*pver+singleleno) = spqs(:)
  out_Var_Flat(16*pver+1+singleleno:17*pver+singleleno) = spqg(:)
  out_Var_Flat(17*pver+1+singleleno:18*pver+singleleno) = spqr(:)
  out_Var_Flat(18*pver+1+singleleno:19*pver+singleleno) = mu_crm(:)
  out_Var_Flat(19*pver+1+singleleno:20*pver+singleleno) = md_crm(:)
  out_Var_Flat(20*pver+1+singleleno:21*pver+singleleno) = du_crm(:)
  out_Var_Flat(21*pver+1+singleleno:22*pver+singleleno) = eu_crm(:)
  out_Var_Flat(22*pver+1+singleleno:23*pver+singleleno) = ed_crm(:)
  out_Var_Flat(23*pver+1+singleleno:24*pver+singleleno) = tkez(:)
  out_Var_Flat(24*pver+1+singleleno:25*pver+singleleno) = tkesgsz(:)
  out_Var_Flat(25*pver+1+singleleno:26*pver+singleleno) = tk_crm(:)
  out_Var_Flat(26*pver+1+singleleno:27*pver+singleleno) = flux_u(:)
  out_Var_Flat(27*pver+1+singleleno:28*pver+singleleno) = flux_qt(:)
  out_Var_Flat(28*pver+1+singleleno:29*pver+singleleno) = fluxsgs_qt(:)
  out_Var_Flat(29*pver+1+singleleno:30*pver+singleleno) = flux_qp(:)
  out_Var_Flat(30*pver+1+singleleno:31*pver+singleleno) = qt_ls(:)
  out_Var_Flat(31*pver+1+singleleno:32*pver+singleleno) = qt_trans(:)
  out_Var_Flat(32*pver+1+singleleno:33*pver+singleleno) = qp_trans(:)
  out_Var_Flat(33*pver+1+singleleno:34*pver+singleleno) = qp_fall(:)
  out_Var_Flat(34*pver+1+singleleno:35*pver+singleleno) = qp_evp(:)
  out_Var_Flat(35*pver+1+singleleno:36*pver+singleleno) = qp_src(:)
  out_Var_Flat(36*pver+1+singleleno:37*pver+singleleno) = t_ls(:)
  out_Var_Flat(37*pver+1+singleleno:37*pver+20+singleleno) = qtotcrm(1:20)

  fcount = 37*pver+20+singleleno
  do ii=1,crm_nx
    do jj=1,crm_ny
      do kk=1,crm_nz
        fcount = fcount + 1
        out_Var_Flat(fcount)              = outin01_crm_u(ii,jj,kk)
        out_Var_Flat(fcount + 1 * chnksz) = outin02_crm_v(ii,jj,kk)
        out_Var_Flat(fcount + 2 * chnksz) = outin03_crm_w(ii,jj,kk)
        out_Var_Flat(fcount + 3 * chnksz) = outin04_crm_t(ii,jj,kk)
        out_Var_Flat(fcount + 4 * chnksz) = outin06_crm_qrad(ii,jj,kk)
        out_Var_Flat(fcount + 5 * chnksz) = outin07_qc_crm(ii,jj,kk)
        out_Var_Flat(fcount + 6 * chnksz) = outin08_qi_crm(ii,jj,kk)
        out_Var_Flat(fcount + 7 * chnksz) = outin09_qpc_crm(ii,jj,kk)
        out_Var_Flat(fcount + 8 * chnksz) = outin10_qpi_crm(ii,jj,kk)
        out_Var_Flat(fcount + 9 * chnksz) = outin12_t_rad(ii,jj,kk)
        out_Var_Flat(fcount + 10* chnksz) = outin13_qv_rad(ii,jj,kk)
        out_Var_Flat(fcount + 11* chnksz) = outin14_qc_rad(ii,jj,kk)
        out_Var_Flat(fcount + 12* chnksz) = outin15_qi_rad(ii,jj,kk) 
        out_Var_Flat(fcount + 13* chnksz) = outin16_cld_rad(ii,jj,kk)
        out_Var_Flat(fcount + 14* chnksz) = outin17_cld3d_crm(ii,jj,kk)
        out_Var_Flat(fcount + 15* chnksz) = crm_tk(ii,jj,kk)
        out_Var_Flat(fcount + 16* chnksz) = crm_tkh(ii,jj,kk)
      end do
    end do
  end do
  fcount = 37*pver+20+singleleno + 16* chnksz
  do ii=1,crm_nx
    do jj=1,crm_ny
      do kk=1,crm_nz
        do ll=1,nmicro_fields
          fcount=fcount+1
          out_Var_Flat(fcount) = outin05_crm_micro(ii,jj,kk,ll) 
        end do
      end do
    end do
  end do
  fcount = 37*pver+20+singleleno + 16* chnksz +  crm_nx*crm_ny*crm_nz*nmicro_fields
  do ii=1,crm_nx
    do jj=1,crm_ny
      fcount = fcount + 1
      out_Var_Flat(fcount) = outin11_prec_crm(ii,jj)
    end do
  end do
  write(13,*) 'CRM chunk,i',myrank_crm, destGCM0,inp01_lchnk,inp02_i
  call MPI_Send(out_Var_Flat,fleno+nflat,MPI_REAL8,destGCM0,8006,MPI_COMM_WORLD,ierr)
 end do
  call MPI_comm_free(crm_comm, ierr)
  call MPI_barrier(MPI_COMM_WORLD, ierr)
  call MPI_finalize(ierr)
end program TwoExecutableDriver
