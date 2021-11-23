program TwoExecutableDriver
  use shr_kind_mod,      only: r8 => SHR_KIND_R8
  use shr_kind_mod,      only: i8 => SHR_KIND_I8
  use ppgrid,           only: pver
  use crmdims,         only: crm_nx, crm_ny, crm_nz
  use crmx_crm_module,     only: crm
  use ppgrid,          only: pcols, begchunk, endchunk
  use spmd_utils, only: npes
  implicit none

  !--------------------------------------------------------------------------
  !bloss: MPI-related variables
  !--------------------------------------------------------------------------
  integer :: crm_comm
  integer :: myrank_global, numproc_global
  integer :: myrank_crm, numproc_crm
  integer :: ierr,status
  integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

  integer :: i,from
  integer, dimension(1) :: destGCM = -999
  CHARACTER(LEN=6) :: crm_number

  ! local copies of input ingredients to CRM ( to be received)
  integer :: inp01_lchnk(1),inp02_i(1) ! pritch note important for integers to
 ! be MPI_received that they masquerade as 1-element arrays here.
  real(r8), dimension(pver) :: inp03_tl,inp04_ql,inp05_qccl,inp06_qiil,inp07_ul,inp08_vl
  real(r8), dimension(pver) :: inp10_pmid,inp11_pdel,inp13_zm,inp14_zi
  real(r8), dimension(1) :: inp09_ps,inp12_phis,inp15_ztodt

! output ingrdients from CRM (to be sent)
  real(r8),dimension(pver) :: out01_qltend,out02_qcltend,out03_qiltend,out04_sltend
! inout ingredients:
  real(r8), dimension(crm_nx,crm_ny,crm_nz) :: outin01_crm_u,outin02_crm_v,outin03_crm_w,outin04_crm_t,&
        outin06_crm_qrad,outin07_qc_crm,outin08_qi_crm,outin09_qpc_crm,outin10_qpi_crm,outin12_t_rad,&
        outin13_qv_rad,outin14_qc_rad,outin15_qi_rad,outin16_cld_rad,outin17_cld3d_crm,crm_tk,crm_tkh
  real(r8), dimension(crm_nx,crm_ny) :: outin11_prec_crm
  real(r8), allocatable :: flattened_crm_inout(:)
  real(r8),allocatable :: outin05_crm_micro(:,:,:,:)
  real (r8), dimension(1) :: precc,precl,precsc,precsl,cltot,clhgh,clmed,cllow
  real (r8), dimension(pver) :: cld,cldtop,gicewp,gliqwp,mctot,mcup,mcdn,mcuup,mcudn
  real (r8), dimension(pver) :: spqc,spqi,spqs,spqg,spqr
  real (r8), dimension(pver) :: mu_crm,md_crm,du_crm,eu_crm,ed_crm,&
                tkez,tkesgsz,tk_crm,flux_u,flux_v,flux_qt,fluxsgs_qt,flux_qp
  real(r8), dimension(pver) :: precflux,qt_ls,qt_trans,qp_trans,qp_fall,qp_evp,qp_src,t_ls,qtotcrm
  real(r8), dimension(1) :: prectend,precstend,ocnfrac,wnd,tau00,bflx,fluxu0,fluxv0,fluxt0,fluxq0
  real(r8), dimension(1) :: taux_crm,z0m,timing_factor,jt_crm,mx_crm
  integer :: chnksz,nflat,fcount,ii,jj,kk,ll,g2ctype
  integer :: g2cMPItype
  type gcm2crm
     real(r8) :: cld(pver)
     real(r8) :: cldtop(pver)
     real(r8) :: gicewp(pver)
     real(r8) :: gliqwp(pver)
     real(r8) :: mctot(pver)
     real(r8) :: mcup(pver)
     real(r8) :: mcdn(pver)
     real(r8) :: mcuup(pver)
     real(r8) :: mcudn(pver)
     real(r8) :: spqc(pver)
     real(r8) :: spqi(pver)
     real(r8) :: spqs(pver)
     real(r8) :: spqg(pver)
     real(r8) :: spqr(pver)
     real(r8) :: mu_crm(pver)
     real(r8) :: md_crm(pver)
     real(r8) :: du_crm(pver)
     real(r8) :: eu_crm(pver)
     real(r8) :: ed_crm(pver)
     real(r8) :: tkez(pver)
     real(r8) :: tkesgsz(pver)
     real(r8) :: tk_crm(pver)
     real(r8) :: flux_u(pver)
     real(r8) :: flux_v(pver)
     real(r8) :: flux_qt(pver)
     real(r8) :: fluxsgs_qt(pver)
     real(r8) :: flux_qp(pver)
     real(r8) :: precflux(pver)
     real(r8) :: qt_ls(pver)
     real(r8) :: qt_trans(pver)
     real(r8) :: qp_trans(pver)
     real(r8) :: qp_fall(pver)
     real(r8) :: qp_evp(pver)
     real(r8) :: qp_src(pver)
     real(r8) :: t_ls(pver)
  end type gcm2crm
  
  type(gcm2crm) :: g2c 

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

     do i = 0,numproc_crm-1
        if(i.eq.myrank_crm) then
           write(13,924) numproc_global, myrank_global, numproc_crm, myrank_crm
924        format('MPI_COMM_WORLD: size/myrank = ',2i5,', crm_comm: size/myrank = ',2i5) 
        end if
        call MPI_Barrier(crm_comm,ierr)
     end do

  ! ----------- GCM handshake from spcam_drivers --------------

  ! Recieve rank of host GCM column linked to this CRM, for eventual MPI_Send
  call MPI_Recv(destGCM,1,MPI_INTEGER,MPI_ANY_SOURCE,54321,MPI_COMM_WORLD,status,ierr)
  if (ierr.eq.0) then
    write(13,*) 'CRM rank',myrank_crm,' got handshake; its GCM dest rank=',destGCM
  else 
    write (13,*) 'MPI_Recv from spcam_driver handshake failed for CRM rank ',myrank_crm,',ierr=',ierr
  end if

  ! INSERT a time step loop here?
  ! ----------- Receive from cesm.exe/crm_physics inputs to CRM right before call to crm() ----------------

! Receive from GCM, in order, the input ingredients expected by crm subroutine:
  call MPI_Recv(inp01_lchnk  ,1   ,MPI_INTEGER,MPI_ANY_SOURCE,9001,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp02_i      ,1   ,MPI_INTEGER,MPI_ANY_SOURCE,9002,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp03_tl(:)  ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9003,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp04_ql(:)  ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9004,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp05_qccl(:),pver,MPI_REAL8  ,MPI_ANY_SOURCE,9005,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp06_qiil(:),pver,MPI_REAL8  ,MPI_ANY_SOURCE,9006,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp07_ul(:)  ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9007,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp08_vl(:)  ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9008,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp09_ps     ,1   ,MPI_REAL8  ,MPI_ANY_SOURCE,9009,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp10_pmid   ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9010,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp11_pdel   ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9011,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp12_phis   ,1   ,MPI_REAL8  ,MPI_ANY_SOURCE,9012,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp13_zm     ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9013,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp14_zi     ,pver,MPI_REAL8  ,MPI_ANY_SOURCE,9014,MPI_COMM_WORLD,status,ierr)
  call MPI_Recv(inp15_ztodt  ,1   ,MPI_REAL8  ,MPI_ANY_SOURCE,9015,MPI_COMM_WORLD,status,ierr)
  call MPI_RECV(g2cMPItype   ,1   ,MPI_INTEGER,MPI_ANY_SOURCE,9016,MPI_COMM_WORLD,status,ierr)
  write (13,*) 'inp01_lchnk=',inp01_lchnk
  write (13,*) 'inp02_i=',inp02_i
  write (13,*) 'inp03_tl=',inp03_tl
  write (13,*) 'inp04_ql=',inp04_ql
  write (13,*) 'inp05_qccl=',inp05_qccl
  write (13,*) 'inp06_qiil=',inp06_qiil
  write (13,*) 'inp07_ul=',inp07_ul
  write (13,*) 'inp08_vl=',inp08_vl
  write (13,*) 'inp09_ps=',inp09_ps
  write (13,*) 'inp10_pmid=',inp10_pmid
  write (13,*) 'inp11_pdel=',inp11_pdel
  write (13,*) 'inp12_phis=',inp12_phis
  write (13,*) 'inp13_zm=',inp13_zm
  write (13,*) 'inp14_zi=',inp14_zi
  write (13,*) 'inp15_ztodt=',inp15_ztodt
  write (13,*) 'g2ctype=',g2cMPItype
  ! call MPI_Send(g2cMPItype,1,MPI_INTEGER,dest,9026,MPI_COMM_WORLD,ierr)
  ! call MPI_Send(g2c,1,g2cMPItype,dest,9025,MPI_COMM_WORLD,ierr)
  call MPI_RECV(g2c%cld(1)   ,1   ,g2cMPItype,MPI_ANY_SOURCE,9017,MPI_COMM_WORLD,status,ierr)
  cld(1:pver) = g2c%cld(1:pver)
  write (13,*) 'cld=',cld(1:pver)
  cldtop(1:pver) = g2c%cldtop(1:pver)
  gicewp(1:pver) = g2c%gicewp(1:pver)
  gliqwp(1:pver) = g2c%gliqwp(1:pver)
  mctot(1:pver) = g2c%mctot(1:pver)
  mcup(1:pver) = g2c%mcup(1:pver)
  mcdn(1:pver) = g2c%mcdn(1:pver)
  mcuup(1:pver) = g2c%mcuup(1:pver)
  mcudn(1:pver) = g2c%mcudn(1:pver)
  spqc(1:pver) = g2c%spqc(1:pver)
  spqi(1:pver) = g2c%spqi(1:pver)
  spqs(1:pver) = g2c%spqs(1:pver)
  spqg(1:pver) = g2c%spqg(1:pver)
  spqr(1:pver) = g2c%spqr(1:pver)
  mu_crm(1:pver) = g2c%mu_crm(1:pver)
  md_crm(1:pver) = g2c%md_crm(1:pver)
  du_crm(1:pver) = g2c%du_crm(1:pver)
  eu_crm(1:pver) = g2c%eu_crm(1:pver)
  ed_crm(1:pver) = g2c%ed_crm(1:pver)
  tkez(1:pver) = g2c%tkez(1:pver)
  tkesgsz(1:pver) = g2c%tkesgsz(1:pver)
  tk_crm(1:pver) = g2c%tk_crm(1:pver)
  flux_u(1:pver) = g2c%flux_u(1:pver)
  flux_v(1:pver) = g2c%flux_v(1:pver)
  flux_qt(1:pver) = g2c%flux_qt(1:pver)
  fluxsgs_qt(1:pver) = g2c%fluxsgs_qt(1:pver)
  flux_qp(1:pver) = g2c%flux_qp(1:pver)
  precflux(1:pver) = g2c%precflux(1:pver)
  qt_ls(1:pver) = g2c%qt_ls(1:pver)
  qt_trans(1:pver) = g2c%qt_trans(1:pver)
  qp_trans(1:pver) = g2c%qp_trans(1:pver)
  qp_fall(1:pver) = g2c%qp_fall(1:pver)
  qp_evp(1:pver) = g2c%qp_evp(1:pver)
  qp_src(1:pver) = g2c%qp_src(1:pver)
  t_ls(1:pver) = g2c%t_ls(1:pver)

  !call MPI_Recv(precc,1,MPI_REAL8,MPI_ANY_SOURCE,9017,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(precl,1,MPI_REAL8,MPI_ANY_SOURCE,9018,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(precsc,1,MPI_REAL8,MPI_ANY_SOURCE,9019,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(precsl,1,MPI_REAL8,MPI_ANY_SOURCE,9020,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(cltot,1,MPI_REAL8,MPI_ANY_SOURCE,9021,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(clhgh,1,MPI_REAL8,MPI_ANY_SOURCE,9022,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(clmed,1,MPI_REAL8,MPI_ANY_SOURCE,9023,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(cllow,1,MPI_REAL8,MPI_ANY_SOURCE,9024,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(jt_crm,pver,MPI_REAL8,MPI_ANY_SOURCE,9044,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(mx_crm,pver,MPI_REAL8,MPI_ANY_SOURCE,9045,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(prectend,1,MPI_REAL8,MPI_ANY_SOURCE,9062,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(precstend,1,MPI_REAL8,MPI_ANY_SOURCE,9063,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(ocnfrac,1,MPI_REAL8,MPI_ANY_SOURCE,9064,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(wnd,1,MPI_REAL8,MPI_ANY_SOURCE,9065,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(tau00,1,MPI_REAL8,MPI_ANY_SOURCE,9066,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(bflx,1,MPI_REAL8,MPI_ANY_SOURCE,9067,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(fluxu0,1,MPI_REAL8,MPI_ANY_SOURCE,9068,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(fluxv0,1,MPI_REAL8,MPI_ANY_SOURCE,9069,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(fluxt0,1,MPI_REAL8,MPI_ANY_SOURCE,9070,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(fluxq0,1,MPI_REAL8,MPI_ANY_SOURCE,9071,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(taux_crm,1,MPI_REAL8,MPI_ANY_SOURCE,9072,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(z0m,1,MPI_REAL8,MPI_ANY_SOURCE,9073,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(timing_factor,1,MPI_REAL8,MPI_ANY_SOURCE,9074,MPI_COMM_WORLD,status,ierr)
  !call MPI_Recv(qtotcrm(:),pver,MPI_REAL8,MPI_ANY_SOURCE,9075,MPI_COMM_WORLD,status,ierr)

#ifdef ITSREADY
! Preparing to call the CRM (assuming we can link it in here!)
  call crm (inp01_lchnk, inp02_i,                            &
            inp03_tl(:),inp04_ql(:),inp05_qccl(:),inp06_qiil(:), &
            inp07_ul(:),inp08_vl(:),inp09_ps,inp10_pmid(:),inp11_pmid(:), &
            inp12_phis(:),inp13_zm(:),inp14_zi(:),inp15_ztodt, &
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
             cam_in%ocnfrac(i),  wnd,                  tau00, bflx,                                          &
             fluxu0,             fluxv0,               fluxt0, fluxq0,                                        &
             taux_crm,        tauy_crm,          z0m, timing_factor,        qtotcrm( :)         &
            )

! ====================== DONE CALLING CRM -- TIME TO SEND OUTPUTS TO GCM

  call MPI_Send(out01_qltend,pver,MPI_REAL8,destGCM,8001,MPI_COMM_WORLD,ierr)
  call MPI_Send(out02_qcltend,pver,MPI_REAL8,destGCM,8002,MPI_COMM_WORLD,ierr)
  call MPI_Send(out03_qiltend,pver,MPI_REAL8,destGCM,8003,MPI_COMM_WORLD,ierr)
  call MPI_Send(out04_sltend,pver,MPI_REAL8,destGCM,8004,MPI_COMM_WORLD,ierr)
  
! pritch, left off here. Need to MPI_Send that huge list of things that was
! intent(inout)

! Liran Start Here ======================================================
! Pack data here
  fcount = 0
  do ii=1,crm_nx
    do jj=1,crm_ny
      do kk=1,crm_nz
        fcount = fcount + 1
        flattened_crm_inout(fcount)              = outin01_crm_u(ii,jj,kk)
        flattened_crm_inout(fcount + 1 * chnksz) = outin02_crm_v(ii,jj,kk)
        flattened_crm_inout(fcount + 2 * chnksz) = outin03_crm_w(ii,jj,kk)
        flattened_crm_inout(fcount + 3 * chnksz) = outin04_crm_t(ii,jj,kk)
        flattened_crm_inout(fcount + 4 * chnksz) = outin06_crm_qrad(ii,jj,kk)
        flattened_crm_inout(fcount + 5 * chnksz) = outin07_qc_crm(ii,jj,kk)
        flattened_crm_inout(fcount + 6 * chnksz) = outin08_qi_crm(ii,jj,kk)
        flattened_crm_inout(fcount + 7 * chnksz) = outin09_qpc_crm(ii,jj,kk)
        flattened_crm_inout(fcount + 8 * chnksz) = outin10_qpi_crm(ii,jj,kk)
        flattened_crm_inout(fcount + 9 * chnksz) = outin12_t_rad(ii,jj,kk)
        flattened_crm_inout(fcount + 10* chnksz) = outin13_qv_rad(ii,jj,kk)
        flattened_crm_inout(fcount + 11* chnksz) = outin14_qc_rad(ii,jj,kk)
        flattened_crm_inout(fcount + 12* chnksz) = outin15_qi_rad(ii,jj,kk) 
        flattened_crm_inout(fcount + 13* chnksz) = outin16_cld_rad(ii,jj,kk)
        flattened_crm_inout(fcount + 14* chnksz) = outin17_cld3d_crm(ii,jj,kk)
        flattened_crm_inout(fcount + 15* chnksz) = crm_tk(ii,jj,kk)
        flattened_crm_inout(fcount + 16* chnksz) = crm_tkh(ii,jj,kk)
      end do
    end do
  end do
  fcount = 0 + 17*chnksz
  do ii=1,crm_nx
    do jj=1,crm_ny
      do kk=1,crm_nz
        do ll=1,nmicro_fields
          fcount=fcount+1
          flattened_crm_inout(fcount) = outin05_crm_micro(ii,jj,kk,ll) 
        end do
      end do
    end do
  end do
   ! next, unpack the 2D precip array.
  fcount = 0 + 17*chnksz + crm_nx*crm_ny*crm_nz*nmicro_fields
  do ii=1,crm_nx
    do jj=1,crm_ny
      fcount = fcount + 1
      flattened_crm_inout(fcount) = outin11_prec_crm(ii,jj)
    end do
  end do
  call MPI_Send(flattened_crm_inout,nflat,MPI_REAL8,destGCM,8005,MPI_COMM_WORLD,ierr)
  call MPI_Send(cld(:),pver,MPI_REAL8,destGCM,8006,MPI_COMM_WORLD,ierr)
  call MPI_Send(cldtop(:),pver,MPI_REAL8,destGCM,8007,MPI_COMM_WORLD,ierr)
  call MPI_Send(gicewp(:),pver,MPI_REAL8,destGCM,8008,MPI_COMM_WORLD,ierr)
  call MPI_Send(gliqwp(:),pver,MPI_REAL8,destGCM,8009,MPI_COMM_WORLD,ierr)
  call MPI_Send(mctot(:),pver,MPI_REAL8,destGCM,8010,MPI_COMM_WORLD,ierr)
  call MPI_Send(mcup(:),pver,MPI_REAL8,destGCM,8011,MPI_COMM_WORLD,ierr)
  call MPI_Send(mcdn(:),pver,MPI_REAL8,destGCM,8012,MPI_COMM_WORLD,ierr)
  call MPI_Send(mcuup(:),pver,MPI_REAL8,destGCM,8013,MPI_COMM_WORLD,ierr)
  call MPI_Send(mcudn(:),pver,MPI_REAL8,destGCM,8014,MPI_COMM_WORLD,ierr)
  call MPI_Send(spqc(:),pver,MPI_REAL8,destGCM,8015,MPI_COMM_WORLD,ierr)
  call MPI_Send(spqi(:),pver,MPI_REAL8,destGCM,8016,MPI_COMM_WORLD,ierr)
  call MPI_Send(spqs(:),pver,MPI_REAL8,destGCM,8017,MPI_COMM_WORLD,ierr)
  call MPI_Send(spqg(:),pver,MPI_REAL8,destGCM,8018,MPI_COMM_WORLD,ierr)
  call MPI_Send(spqr(:),pver,MPI_REAL8,destGCM,8019,MPI_COMM_WORLD,ierr)
  call MPI_Send(mu_crm(:),pver,MPI_REAL8,destGCM,8020,MPI_COMM_WORLD,ierr)
  call MPI_Send(md_crm(:),pver,MPI_REAL8,destGCM,8021,MPI_COMM_WORLD,ierr)
  call MPI_Send(du_crm(:),pver,MPI_REAL8,destGCM,8022,MPI_COMM_WORLD,ierr)
  call MPI_Send(eu_crm(:),pver,MPI_REAL8,destGCM,8023,MPI_COMM_WORLD,ierr)
  call MPI_Send(ed_crm(:),pver,MPI_REAL8,destGCM,8024,MPI_COMM_WORLD,ierr)
  call MPI_Send(tkez(:),pver,MPI_REAL8,destGCM,8025,MPI_COMM_WORLD,ierr)
  call MPI_Send(tkesgsz(:),pver,MPI_REAL8,destGCM,8026,MPI_COMM_WORLD,ierr)
  call MPI_Send(tk_crm(:),pver,MPI_REAL8,destGCM,8027,MPI_COMM_WORLD,ierr)
  call MPI_Send(flux_u(:),pver,MPI_REAL8,destGCM,8028,MPI_COMM_WORLD,ierr)
  call MPI_Send(flux_v(:),pver,MPI_REAL8,destGCM,8029,MPI_COMM_WORLD,ierr)
  call MPI_Send(flux_qt(:),pver,MPI_REAL8,destGCM,8030,MPI_COMM_WORLD,ierr)
  call MPI_Send(fluxsgs_qt(:),pver,MPI_REAL8,destGCM,8031,MPI_COMM_WORLD,ierr)
  call MPI_Send(flux_qp(:),pver,MPI_REAL8,destGCM,8032,MPI_COMM_WORLD,ierr)
  call MPI_Send(qt_ls(:),pver,MPI_REAL8,destGCM,8033,MPI_COMM_WORLD,ierr)
  call MPI_Send(qt_trans(:),pver,MPI_REAL8,destGCM,8034,MPI_COMM_WORLD,ierr)
  call MPI_Send(qp_trans(:),pver,MPI_REAL8,destGCM,8035,MPI_COMM_WORLD,ierr)
  call MPI_Send(qp_fall(:),pver,MPI_REAL8,destGCM,8036,MPI_COMM_WORLD,ierr)
  call MPI_Send(qp_evp(:),pver,MPI_REAL8,destGCM,8037,MPI_COMM_WORLD,ierr)
  call MPI_Send(qp_src(:),pver,MPI_REAL8,destGCM,8038,MPI_COMM_WORLD,ierr)
  call MPI_Send(t_ls(:),pver,MPI_REAL8,destGCM,8039,MPI_COMM_WORLD,ierr)
  call MPI_Send(qtotcrm(:),pver,MPI_REAL8,destGCM,8040,MPI_COMM_WORLD,ierr)
  call MPI_Send(precc,1,MPI_REAL8,destGCM,8041,MPI_COMM_WORLD,ierr)
  call MPI_Send(precl,1,MPI_REAL8,destGCM,8042,MPI_COMM_WORLD,ierr)
  call MPI_Send(precsc,1,MPI_REAL8,destGCM,8043,MPI_COMM_WORLD,ierr)
  call MPI_Send(precsl,1,MPI_REAL8,destGCM,8044,MPI_COMM_WORLD,ierr)
  call MPI_Send(cltot,1,MPI_REAL8,destGCM,8045,MPI_COMM_WORLD,ierr)
  call MPI_Send(clhgh,1,MPI_REAL8,destGCM,8046,MPI_COMM_WORLD,ierr)
  call MPI_Send(clmed,1,MPI_REAL8,destGCM,8047,MPI_COMM_WORLD,ierr)
  call MPI_Send(cllow,1,MPI_REAL8,destGCM,8048,MPI_COMM_WORLD,ierr)
  call MPI_Send(prectend,1,MPI_REAL8,destGCM,8049,MPI_COMM_WORLD,ierr)
  call MPI_Send(precstend,1,MPI_REAL8,destGCM,8050,MPI_COMM_WORLD,ierr)
  call MPI_Send(ocnfrac,1,MPI_REAL8,destGCM,8051,MPI_COMM_WORLD,ierr)
  call MPI_Send(wnd,1,MPI_REAL8,destGCM,8052,MPI_COMM_WORLD,ierr)
  call MPI_Send(tau00,1,MPI_REAL8,destGCM,8053,MPI_COMM_WORLD,ierr)
  call MPI_Send(bflx,1,MPI_REAL8,destGCM,8054,MPI_COMM_WORLD,ierr)
  call MPI_Send(fluxu0,1,MPI_REAL8,destGCM,8055,MPI_COMM_WORLD,ierr)
  call MPI_Send(fluxv0,1,MPI_REAL8,destGCM,8056,MPI_COMM_WORLD,ierr)
  call MPI_Send(fluxt0,1,MPI_REAL8,destGCM,8057,MPI_COMM_WORLD,ierr)
  call MPI_Send(fluxq0,1,MPI_REAL8,destGCM,8058,MPI_COMM_WORLD,ierr)
  call MPI_Send(taux_crm,1,MPI_REAL8,destGCM,8059,MPI_COMM_WORLD,ierr)
  call MPI_Send(z0m,1,MPI_REAL8,destGCM,8060,MPI_COMM_WORLD,ierr)
  call MPI_Send(timing_factor,1,MPI_REAL8,destGCM,8061,MPI_COMM_WORLD,ierr)
! Liran End Here   ======================================================

#endif

  call MPI_comm_free(crm_comm, ierr)

  call MPI_barrier(MPI_COMM_WORLD, ierr)

  call MPI_finalize(ierr)

end program TwoExecutableDriver
