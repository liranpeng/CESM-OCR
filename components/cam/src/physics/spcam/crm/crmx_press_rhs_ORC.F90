module crmx_crm_press_rhs_ORC
use crmx_mpi
use crmx_vars
use crmx_params, only: dowallx, dowally
use crmx_task_util_mpi

implicit none
private
save

public :: press_rhs_ORC

contains

subroutine press_rhs_ORC

!       right-hand-side of the Poisson equation for pressure
         
	
real *8 dta,rdx,rdy,rdz,btat,ctat,rup,rdn
integer i,j,k,ic,jc,kc,ierr
integer rank,numtasks

include 'mpif.h'

call task_start_ORC(rank,numtasks)

if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then

    do k=1,nzm
     do j=1,ny
      dudt(1,j,k,na) = 0.
     end do
    end do	

end if

if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

    do k=1,nzm
     do i=1,nx
      dvdt(i,1,k,na) = 0.
     end do
    end do	

end if

call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)

!write(0, *) 'Liran CRM_ORC1 dudt',myrank_global,nx,ny,nsubdomains_x,dudt

if(dompi) then
   call task_bound_duvdt_ORC()
else
   call bound_duvdt()	   
endif

!write(0, *) 'Liran CRM_ORC2 dudt',myrank_global,dudt

dta=1./dt3(na)/at
rdx=1./dx
rdy=1./dy
btat=bt/at
ctat=ct/at

if(RUN3D) then

do k=1,nzm
 kc=k+1 
 rdz=1./(adz(k)*dz)
 rup = rhow(kc)/rho(k)*rdz
 rdn = rhow(k)/rho(k)*rdz
 do j=1,ny
  jc=j+1 
  do i=1,nx
   ic=i+1
   p(i,j,k)=(rdx*(u(ic,j,k)-u(i,j,k))+ &
             rdy*(v(i,jc,k)-v(i,j,k))+ &
             (w(i,j,kc)*rup-w(i,j,k)*rdn) )*dta + &	
            (rdx*(dudt(ic,j,k,na)-dudt(i,j,k,na))+ &
             rdy*(dvdt(i,jc,k,na)-dvdt(i,j,k,na))+ &
             (dwdt(i,j,kc,na)*rup-dwdt(i,j,k,na)*rdn) ) + &
       btat*(rdx*(dudt(ic,j,k,nb)-dudt(i,j,k,nb))+ &
             rdy*(dvdt(i,jc,k,nb)-dvdt(i,j,k,nb))+ &
             (dwdt(i,j,kc,nb)*rup-dwdt(i,j,k,nb)*rdn) ) + &
       ctat*(rdx*(dudt(ic,j,k,nc)-dudt(i,j,k,nc))+ &
             rdy*(dvdt(i,jc,k,nc)-dvdt(i,j,k,nc))+ &
             (dwdt(i,j,kc,nc)*rup-dwdt(i,j,k,nc)*rdn) )
   p(i,j,k)=p(i,j,k)*rho(k)
  end do
 end do
end do


else

j=1

do k=1,nzm
 kc=k+1 
 rdz=1./(adz(k)*dz)
 rup = rhow(kc)/rho(k)*rdz
 rdn = rhow(k)/rho(k)*rdz
 do i=1,nx
  ic=i+1
  p(i,j,k)=(rdx*(u(ic,j,k)-u(i,j,k))+ &
                (w(i,j,kc)*rup-w(i,j,k)*rdn) )*dta + &
                (rdx*(dudt(ic,j,k,na)-dudt(i,j,k,na))+ &
                (dwdt(i,j,kc,na)*rup-dwdt(i,j,k,na)*rdn) ) + &
           btat*(rdx*(dudt(ic,j,k,nb)-dudt(i,j,k,nb))+ &
                 (dwdt(i,j,kc,nb)*rup-dwdt(i,j,k,nb)*rdn) ) + &
           ctat*(rdx*(dudt(ic,j,k,nc)-dudt(i,j,k,nc))+ &
                 (dwdt(i,j,kc,nc)*rup-dwdt(i,j,k,nc)*rdn) )
  p(i,j,k)=p(i,j,k)*rho(k)
  if(p(i,j,k) .gt. 100) then
  print*,'Liran check p rhs0',myrank_global,i,k,p(i,j,k),u(ic,j,k),u(i,j,k)
  print*,'Liran check p rh0q',myrank_global,i,k,w(i,j,kc),w(i,j,k)
  print*,'Liran check p rh1q',myrank_global,i,k,dudt(ic,j,k,na),dudt(i,j,k,na),dudt(i,j,k,nb)
  print*,'Liran check p rh2q',myrank_global,i,k,dwdt(i,j,kc,na),dwdt(i,j,k,na)
  print*,'Liran check p rhs1',myrank_global,i,k,(rdx*(u(ic,j,k)-u(i,j,k))+(w(i,j,kc)*rup-w(i,j,k)*rdn))*dta
  print*,'Liran check p rhs2',myrank_global,i,k,(dwdt(i,j,kc,na)*rup-dwdt(i,j,k,na)*rdn) 
  print*,'Liran check p rhs3',myrank_global,i,k,btat*(rdx*(dudt(ic,j,k,nb)-dudt(i,j,k,nb))+(dwdt(i,j,kc,nb)*rup-dwdt(i,j,k,nb)*rdn))
  print*,'Liran check p rhs4',myrank_global,i,k,ctat*(rdx*(dudt(ic,j,k,nc)-dudt(i,j,k,nc))+(dwdt(i,j,kc,nc)*rup-dwdt(i,j,k,nc)*rdn))
  print*,'Liran check p rhs5',myrank_global,i,k,(w(i,j,kc)*rup-w(i,j,k)*rdn)*dta
  end if
 end do
end do


endif
if(dompi) then
  call task_barrier_ORC()
else 
  call task_barrier()
end if
end subroutine press_rhs_ORC

end module crmx_crm_press_rhs_ORC
