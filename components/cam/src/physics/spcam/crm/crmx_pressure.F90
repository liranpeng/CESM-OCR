! Non-blocking receives before blocking sends

subroutine pressure
	
!       Original pressure solver based on horizontal slabs
!       (C) 1998, 2002 Marat Khairoutdinov
!       Works only when the number of slabs is equal to the number of processors.
!       Therefore, the number of processors shouldn't exceed the number of levels nzm
!       Also, used for a 2D version 
!       For more processors for the given number of levels and 3D, use pressure_big

use crmx_vars
use crmx_params, only: dowallx, dowally, docolumn
use crmx_grid, only: dompi
use crmx_task_util_mpi

implicit none
	
integer :: npressureslabs,nzslab,nx2,ny2,n3i,n3j
real, allocatable, dimension(:,:,:) :: fp
real, allocatable, dimension(:,:,:) :: ff
real, allocatable, dimension(:,:,:,:) :: buff_slabs
real, allocatable, dimension(:,:,:,:) :: buff_subs
real, allocatable, dimension(:,:,:,:) :: bufp_slabs
real, allocatable, dimension(:,:,:,:) :: bufp_subs
real, allocatable, dimension(:,:) :: work
real, allocatable, dimension(:) :: trigxi
real, allocatable, dimension(:) :: trigxj
real(kind=selected_real_kind(12)), allocatable, dimension(:) :: a
real(kind=selected_real_kind(12)), allocatable, dimension(:) :: c
real(kind=selected_real_kind(12)), allocatable, dimension(:) :: fff
real(kind=selected_real_kind(12)), allocatable, dimension(:) :: alfa
real(kind=selected_real_kind(12)), allocatable, dimension(:) :: beta
integer, allocatable, dimension(:) :: reqs_in
integer, allocatable, dimension(:) :: iii
integer, allocatable, dimension(:) :: jjj
logical flag(nsubdomains)

integer ifaxj(100),ifaxi(100)

real(kind=selected_real_kind(12)) b,e
real(kind=selected_real_kind(12)) xi,xj,xnx,xny,ddx2,ddy2,pii,factx,facty,eign

integer i, j, k, id, jd, m, n, it, jt, ii, jj, tag, rf
integer nyp22, n_in, count
integer iwall,jwall,myrank_global,ierr
integer,parameter :: DBL = selected_real_kind(12)

include 'mpif.h'

npressureslabs = nsubdomains
nzslab = max(1,nzm / npressureslabs)
nx2=nx_gl+2
ny2=ny_gl+2*YES3D
n3i=3*nx_gl/2+1
n3j=3*ny_gl/2+1

 allocate (fp(nx2,ny2,nzslab))  ! global rhs and array for FTP coefficeients
 allocate (ff(nx+1,ny+2*YES3D,nzm))  ! local (subdomain's) version of f
 allocate (buff_slabs(nxp1,nyp2,nzslab,npressureslabs))
 allocate (buff_subs(nxp1,nyp2,nzslab,nsubdomains))
 allocate (bufp_slabs(nx, (1-YES3D):ny, nzslab,npressureslabs))
 allocate (bufp_subs(nx, (1-YES3D):ny, nzslab,nsubdomains))
 allocate (work(nx2,ny2))
 allocate (trigxi(n3i))
 allocate (trigxj(n3j))
 allocate (a(nzm))
 allocate (c(nzm))
 allocate (fff(nzm))
 allocate (alfa(nzm-1))
 allocate (beta(nzm-1))
 allocate (reqs_in(nsubdomains))
 allocate (iii(0:nx_gl))
 allocate (jjj(0:ny_gl))

 fp=0.
 ff=0.
 buff_slabs=0.
 buff_subs=0.
 bufp_slabs=0.
 bufp_subs=0.
 work=0.
 trigxi=0.
 trigxj=0.
 a=0.
 c=0.
 fff=0.
 alfa=0.
 beta=0.
 reqs_in=0
 iii=0
 jjj=0
 flag=0

!common/tmpstack/f,ff,buff_slabs,buff_subs

! check if the grid size allows the computation:


if(nsubdomains.gt.nzm) then
  if(masterproc) print*,'pressure_orig: nzm < nsubdomains. STOP'
  call task_abort
endif

if(mod(nzm,npressureslabs).ne.0) then
  if(masterproc) print*,'pressure_orig: nzm/npressureslabs is not round number. STOP'
  call task_abort
endif

!-----------------------------------------------------------------

if(docolumn) return

if(dowallx) then
  iwall=1
else
  iwall=0
end if
if(RUN2D) then  
  nyp22=1
  jwall=0
else
  nyp22=nyp2
  if(dowally) then
    jwall=2
  else
    jwall=0
  end if
endif
	
!-----------------------------------------------------------------
!  Compute the r.h.s. of the Poisson equation for pressure
call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
if(myrank_global.eq.0) then
  !print*, 'Liran Check u0 org ',u
  print*, 'Liran Check p0 org ',p
end if
call press_rhs()
if(myrank_global.eq.0) then
  print*, 'Liran Check p org ',p
end if 
!-----------------------------------------------------------------	 
!   Form the horizontal slabs of right-hand-sides of Poisson equation 
!   for the global domain. Request sending and receiving tasks.

!  Non-blocking receive first:

n_in = 0
do m = 0,nsubdomains-1

  if(rank.lt.npressureslabs.and.m.ne.nsubdomains-1) then

    n_in = n_in + 1
    if(dompi) then
      call task_receive_float_ORC(bufp_subs(0,1-YES3D,1,n_in), &
                             nzslab*nxp1*nyp1,reqs_in(n_in))
    else
      call task_receive_float(bufp_subs(0,1-YES3D,1,n_in), &
                             nzslab*nxp1*nyp1,reqs_in(n_in))
    end if
     do k = 1,nzslab 
      do j = 1,nyp2 
       do i = 1,nxp1 
         buff_subs(i,j,k,n_in) = bufp_subs(i,j,k,n_in) 
       end do 
      end do 
     end do    
 
    flag(n_in) = .false.
! if(myrank_global.eq.0) then
!  print*, 'Liran p org ',p
!  end if 
  endif
 ! if(myrank_global.eq.0) then
 ! print*, 'Liran buff_subs org ',buff_subs
 ! end if
  if(rank.lt.npressureslabs.and.m.eq.nsubdomains-1) then

    if(dompi) then
      call task_rank_to_index_ORC(rank,it,jt)	  
    else
      call task_rank_to_index(rank,it,jt)
    end if
    n = rank*nzslab
    do k = 1,nzslab
     do j = 1,ny
       do i = 1,nx
         fp(i+it,j+jt,k) = p(i,j,k+n)
       end do
     end do
    end do
  endif

end do ! m
  !if(myrank_global.eq.0) then
  !print*, 'Liran fp org ',it,fp
  !print*, 'Liran p001 org ',it,p
  !end if

! Blocking send now:


do m = 0,nsubdomains-1

  if(m.lt.npressureslabs.and.m.ne.rank) then

    n = m*nzslab + 1
    if(dompi) then
      call task_bsend_float_ORC(m,p(0,1-YES3D,n),nzslab*nxp1*nyp1, 33)
    else
      call task_bsend_float(m,p(0,1-YES3D,n),nzslab*nxp1*nyp1, 33)
    end if
  endif

end do ! m
  !if(myrank_global.eq.0) then
  !print*, 'Liran p 2org ',p
  !end if

! Fill slabs when receive buffers are full:

count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
        if(dompi) then
          call task_test_ORC(reqs_in(m), flag(m), rf, tag)
        else
	  call task_test(reqs_in(m), flag(m), rf, tag)
        end if
        if(flag(m)) then 
	   count=count-1
           if(dompi) then
             call task_rank_to_index_ORC(rf,it,jt)
           else
             call task_rank_to_index(rf,it,jt)	  
           end if
           do k = 1,nzslab
            do j = 1,ny
             do i = 1,nx
               fp(i+it,j+jt,k) = bufp_subs(i,j,k,m)
             end do
            end do
           end do
	endif   
   endif
  end do
end do
  !if(myrank_global.eq.0) then
  !print*, 'Liran fp 2 org ',fp
  !end if
!-------------------------------------------------
! Perform Fourier transformation for a slab:

if(rank.lt.npressureslabs) then

 call fftfax_crm(nx_gl,ifaxi,trigxi)
 if(RUN3D) call fftfax_crm(ny_gl,ifaxj,trigxj)
print*,'org fp',rank,fp
 do k=1,nzslab
print*,'org fp k',k,fp(1,1,k)
   call fft991_crm(fp(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,-1)

  if(RUN3D) then
     call fft991_crm(fp(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,-1)
  end if

 end do 

endif


! Synchronize all slabs:
if(dompi) then
  call task_barrier_ORC()
else
  call task_barrier()
end if

!-------------------------------------------------
!   Send Fourier coeffiecients back to subdomains:

! Non-blocking receive first:

n_in = 0
do m = 0, nsubdomains-1
		
   if(dompi) then
     call task_rank_to_index_ORC(m,it,jt)
  else
     call task_rank_to_index(m,it,jt)
  end if

   if(rank.lt.npressureslabs.and.m.eq.rank) then

     n = rank*nzslab
     do k = 1,nzslab
      do j = 1,nyp22-jwall
        do i = 1,nxp1-iwall
          ff(i,j,k+n) = fp(i+it,j+jt,k) 
        end do
      end do
     end do 

   end if

   if(m.lt.npressureslabs-1.or.m.eq.npressureslabs-1 &
                            .and.rank.ge.npressureslabs) then

     n_in = n_in + 1
   if(dompi) then
     call task_receive_float_ORC(buff_slabs(1,1,1,n_in), &
                                nzslab*nxp1*nyp22,reqs_in(n_in))
   else
     call task_receive_float(buff_slabs(1,1,1,n_in), &
                                nzslab*nxp1*nyp22,reqs_in(n_in))
   end if


     do k = 1,nzslab
      do j = 1-YES3D,ny
       do i = 1,nx
         bufp_slabs(i,j,k,n_in) = buff_slabs(i,j,k,n_in)
       end do
      end do
     end do

     flag(n_in) = .false.	    
   endif

end do ! m

! Blocking send now:

do m = 0, nsubdomains-1

   call task_rank_to_index(m,it,jt)

   if(rank.lt.npressureslabs.and.m.ne.rank) then

     do k = 1,nzslab
      do j = 1,nyp22
       do i = 1,nxp1
         buff_subs(i,j,k,1) = fp(i+it,j+jt,k)
         bufp_subs(i,j,k,1) = fp(i+it,j+jt,k)
       end do
      end do
     end do
   if(dompi) then
     call task_bsend_float_ORC(m, buff_subs(1,1,1,1),nzslab*nxp1*nyp22,44)
   else
     call task_bsend_float(m, buff_subs(1,1,1,1),nzslab*nxp1*nyp22,44)
   endif

   endif

end do ! m



! Fill slabs when receive buffers are complete:


count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
     if(dompi) then
        call task_test_ORC(reqs_in(m), flag(m), rf, tag)
     else
	call task_test(reqs_in(m), flag(m), rf, tag)
     end if
        if(flag(m)) then 
	   count=count-1
           n = rf*nzslab           
           do k = 1,nzslab
             do j=1,nyp22
               do i=1,nxp1
                 ff(i,j,k+n) = buff_slabs(i,j,k,m)
               end do
             end do
           end do
	endif   
   endif
  end do
end do
!-------------------------------------------------
!   Solve the tri-diagonal system for Fourier coeffiecients 
!   in the vertical for each subdomain:

do k=1,nzm
    a(k)=rhow(k)/(adz(k)*adzw(k)*dz*dz)
    c(k)=rhow(k+1)/(adz(k)*adzw(k+1)*dz*dz)	 
end do 
if(dompi) then
  call task_rank_to_index_ORC(rank,it,jt)
else
  call task_rank_to_index(rank,it,jt)
end if
	
ddx2=1._DBL/(dx*dx)
ddy2=1._DBL/(dy*dy)
pii = acos(-1._DBL)
xnx=pii/nx_gl
xny=pii/ny_gl 	 
do j=1,nyp22-jwall
   if(dowally) then
      jd=j+jt-1
      facty = 1.d0
   else
      jd=(j+jt-0.1)/2.
      facty = 2.d0
   end if
   xj=jd
   do i=1,nxp1-iwall
      if(dowallx) then
        id=i+it-1
        factx = 1.d0
      else
        id=(i+it-0.1)/2.
        factx = 2.d0
      end if
      fff(1:nzm) = ff(i,j,1:nzm)
      xi=id
      eign=(2._DBL*cos(factx*xnx*xi)-2._DBL)*ddx2+ & 
            (2._DBL*cos(facty*xny*xj)-2._DBL)*ddy2
      if(id+jd.eq.0) then               
         b=1._DBL/(eign*rho(1)-a(1)-c(1))
         alfa(1)=-c(1)*b
         beta(1)=fff(1)*b
      else
         b=1._DBL/(eign*rho(1)-c(1))
         alfa(1)=-c(1)*b
         beta(1)=fff(1)*b
      end if
      do k=2,nzm-1
        e=1._DBL/(eign*rho(k)-a(k)-c(k)+a(k)*alfa(k-1))
        alfa(k)=-c(k)*e
        beta(k)=(fff(k)-a(k)*beta(k-1))*e
      end do

      fff(nzm)=(fff(nzm)-a(nzm)*beta(nzm-1))/ &
	        (eign*rho(nzm)-a(nzm)+a(nzm)*alfa(nzm-1))
	  
      do k=nzm-1,1,-1
       fff(k)=alfa(k)*fff(k+1)+beta(k)
      end do
      ff(i,j,1:nzm) = fff(1:nzm)

   end do  
end do 
if(dompi) then
  call task_barrier_ORC()
else
  call task_barrier()
end if
!-----------------------------------------------------------------	 
!   Send the Fourier coefficient to the tasks performing
!   the inverse Fourier transformation:

! Non-blocking receive first:

n_in = 0
do m = 0,nsubdomains-1

  if(rank.lt.npressureslabs.and.m.ne.nsubdomains-1) then
    n_in = n_in + 1
    if(dompi) then
      call task_receive_float_ORC(buff_subs(1,1,1,n_in), &
                                nzslab*nxp1*nyp22, reqs_in(n_in))
    else
      call task_receive_float(buff_subs(1,1,1,n_in), &
                                nzslab*nxp1*nyp22, reqs_in(n_in))
    end if
    flag(n_in) = .false.	    
  endif

  if(rank.lt.npressureslabs.and.m.eq.nsubdomains-1) then

    if(dompi) then
      call task_rank_to_index_ORC(rank,it,jt)
    else
      call task_rank_to_index(rank,it,jt)	  
    end if
    n = rank*nzslab
    do k = 1,nzslab
     do j = 1,nyp22-jwall
       do i = 1,nxp1-iwall
         fp(i+it,j+jt,k) = ff(i,j,k+n)
       end do
     end do
    end do

  endif

end do ! m

! Blocking send now:

do m = 0,nsubdomains-1

  if(m.lt.npressureslabs.and.m.ne.rank) then
    n = m*nzslab+1
    if(dompi) then
      call task_bsend_float_ORC(m,ff(1,1,n),nzslab*nxp1*nyp22, 33)
    else
      call task_bsend_float(m,ff(1,1,n),nzslab*nxp1*nyp22, 33)
    end if
  endif

end do ! m


! Fill slabs when receive buffers are full:


count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
        if(dompi) then
        call task_test_ORC(reqs_in(m), flag(m), rf, tag)
        else
	  call task_test(reqs_in(m), flag(m), rf, tag)
        end if
        if(flag(m)) then 
	   count=count-1
           if(dompi) then
             call task_rank_to_index_ORC(rf,it,jt)
           else
             call task_rank_to_index(rf,it,jt)	  
           end if
           do k = 1,nzslab
            do j = 1,nyp22-jwall
             do i = 1,nxp1-iwall
                fp(i+it,j+jt,k) = buff_subs(i,j,k,m)
             end do
            end do
           end do
	endif   
   endif
  end do
end do
!if(myrank_global.eq.0) then
!  print*, 'Liran Check fp org ',fp
!end if
!-------------------------------------------------
!   Perform inverse Fourier transformation:

if(rank.lt.npressureslabs) then

 do k=1,nzslab

  if(RUN3D) then
     call fft991_crm(fp(1,1,k),work,trigxj,ifaxj,nx2,1,ny_gl,nx_gl+1,+1)
  end if
	 
   call fft991_crm(fp(1,1,k),work,trigxi,ifaxi,1,nx2,nx_gl,ny_gl,+1)

 end do 

endif
if(dompi) then
  call task_barrier_ORC()
else
  call task_barrier()
end if
!-----------------------------------------------------------------	 
!   Fill the pressure field for each subdomain: 

do i=1,nx_gl
 iii(i)=i
end do
iii(0)=nx_gl
do j=1,ny_gl
 jjj(j)=j
end do
jjj(0)=ny_gl

! Non-blocking receive first:

n_in = 0
do m = 0, nsubdomains-1
		
   if(dompi) then
     call task_rank_to_index_ORC(m,it,jt)
   else
     call task_rank_to_index(m,it,jt)
   end if

   if(m.lt.npressureslabs-1.or.  &
		m.eq.npressureslabs-1.and.rank.ge.npressureslabs) then

     n_in = n_in + 1
     if(dompi) then
       call task_receive_float_ORC(bufp_slabs(0,1-YES3D,1,n_in), &
                                  nzslab*nxp1*nyp1, reqs_in(n_in))
     else
       call task_receive_float(bufp_slabs(0,1-YES3D,1,n_in), &
                                  nzslab*nxp1*nyp1, reqs_in(n_in))
     end if
     do k = 1,nzslab 
      do j = 1-YES3D,ny 
       do i = 0,nx 
         buff_slabs(i,j,k,n_in) = bufp_slabs(i,j,k,n_in) 
       end do 
      end do 
     end do 
     flag(n_in) = .false.    

   endif

   if(rank.lt.npressureslabs.and.m.eq.rank) then

     n = rank*nzslab
     do k = 1,nzslab
      do j = 1-YES3D,ny
       jj=jjj(j+jt)
        do i = 0,nx
	 ii=iii(i+it)
          p(i,j,k+n) = fp(ii,jj,k) 
        end do
      end do
     end do 

   end if

end do ! m


! Blocking send now:

do m = 0, nsubdomains-1
   if(dompi) then
     call task_rank_to_index_ORC(m,it,jt)
   else
     call task_rank_to_index(m,it,jt)
   end if

   if(rank.lt.npressureslabs.and.m.ne.rank) then

     do k = 1,nzslab
      do j = 1-YES3D,ny
       jj=jjj(j+jt)
       do i = 0,nx
         ii=iii(i+it)
         bufp_subs(i,j,k,1) = fp(ii,jj,k)
         buff_subs(i,j,k,1) = bufp_subs(i,j,k,1)
       end do
      end do
     end do
     if(dompi) then
       call task_bsend_float_ORC(m, bufp_subs(0,1-YES3D,1,1), nzslab*nxp1*nyp1,44)
     else
       call task_bsend_float(m, bufp_subs(0,1-YES3D,1,1), nzslab*nxp1*nyp1,44)
     end if
   endif

end do ! m

! Fill the receive buffers:

count = n_in
do while (count .gt. 0)
  do m = 1,n_in
   if(.not.flag(m)) then
        if(dompi) then
          call task_test_ORC(reqs_in(m), flag(m), rf, tag)
        else
	  call task_test(reqs_in(m), flag(m), rf, tag)
        end if
        if(flag(m)) then 
	   count=count-1
           n = rf*nzslab           
           do k = 1,nzslab
            do j=1-YES3D,ny
             do i=0,nx
               p(i,j,k+n) = bufp_slabs(i,j,k,m)
             end do
            end do
           end do
        endif   
   endif
  end do
end do

call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
!if(myrank_global.eq.0) then
!  print*, 'Liran Check p org 2 ',p
!end if

if(dompi) then
  call task_barrier_ORC()
else
  call task_barrier()
end if
!  Add pressure gradient term to the rhs of the momentum equation:

call press_grad()


 deallocate (fp)  ! global rhs and array for FTP coefficeients
 deallocate (ff)  ! local (subdomain's) version of f
 deallocate (buff_slabs)
 deallocate (buff_subs)
 deallocate (bufp_slabs)
 deallocate (bufp_subs)
 deallocate (work)
 deallocate (trigxi)
 deallocate (trigxj)
 deallocate (a)
 deallocate (c)
 deallocate (fff)
 deallocate (alfa)
 deallocate (beta)
 deallocate (reqs_in)
 deallocate (iii)
 deallocate (jjj)

end 



