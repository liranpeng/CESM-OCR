subroutine diffuse_scalar (fdiff,flux,f2lediff,f2lediss,fwlediff,doit)

use crmx_grid
use crmx_vars, only: rho, rhow,fluxb,fluxt,f,tkh
use crmx_sgs
implicit none

! input:	
!real,allocatable, intent(in) :: f(:,:,:)	! scalar
!real,allocatable, intent(in) :: fluxb(:,:)		! bottom flux
!real,allocatable, intent(in) :: fluxt(:,:)		! top flux
!real, allocatable, dimension(:,:,:) :: f
!real, allocatable, dimension(:,:)  :: fluxb,fluxt
!real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! scalar
!real fluxb(nx,ny)               ! bottom flux
!real fluxt(nx,ny)               ! top flux
real flux(nz)
real f2lediff(nz),f2lediss(nz),fwlediff(nz)
real fdiff(nz)
logical doit
! Local
!real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real, allocatable, dimension(:,:,:) :: df
integer i,j,k

!call t_startf ('diffuse_scalars')

!if (.not. allocated(fluxb)) allocate ( fluxb(nx,ny))
!if (.not. allocated(fluxt)) allocate ( fluxt(nx,ny))
!if (.not. allocated(f)) allocate ( f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
allocate ( df(nx, ny, nzm))

do k=1,nzm
   do j=1,ny
    do i=1,nx
     df(i,j,k)=f(i,j,k)
    end do
   end do
end do

!df(:,:,:) = f(:,:,:)

if(RUN3D) then
  call diffuse_scalar3D (rho,rhow,flux)
else  
  call diffuse_scalar2D (rho,rhow,flux)
endif

do k=1,nzm
   fdiff(k)=0.
   do j=1,ny
    do i=1,nx
     fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
    end do
   end do
end do

!call t_stopf ('diffuse_scalars')
!deallocate(df)
end subroutine diffuse_scalar 
