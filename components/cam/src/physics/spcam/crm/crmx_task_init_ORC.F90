module crmx_task_init_mpi
use crmx_mpi
use crmx_grid
use crmx_task_util_mpi

implicit none
private
save

public :: task_init_ORC

contains

subroutine task_init_ORC(ncomm,npro,ntask)
		
! crm_comm_in,numproc_crm_in,myrank_crm_in
!   Check things, initialize multitasking:
!use crmx_mpi
!use crmx_grid
implicit none
integer, intent(in) :: ncomm,npro,ntask
integer itasks,ntasks,numproc_in_check,myrank_in_check,ierr,crm_comm_in_dup
include 'mpif.h'
print *,'Enter ORC'
if(YES3D .ne. 1 .and. YES3D .ne. 0) then
  print*,'YES3D is not 1 or 0. STOP'
  stop
endif

if(YES3D .eq. 1 .and. ny_gl .lt. 4) then
  print*,'ny_gl is too small for a 3D case.STOP'
  stop
endif

if(YES3D .eq. 0 .and. ny_gl .ne. 1) then
  print*,'ny_gl should be 1 for a 2D case. STOP'
  stop
endif
if(nsubdomains.eq.1) then

  rank =0
  ntasks = 1
  dompi = .false.
  print *,'Something is wrong here',npro,ntask
else
  print *,'Enter task start ORC',crm_comm_in, numproc_crm_in,myrank_crm_in
  call mpi_comm_size(crm_comm_in, numproc_in_check, ierr)
  call mpi_comm_rank(crm_comm_in, myrank_in_check, ierr) 
  print *,'Liran check888',numproc_in_check,myrank_in_check
  ntasks = myrank_in_check
  call task_start_ORC(rank,ntasks)
  print *,'start finish',dompi,rank,ntasks

  !call systemf('hostname')

  if(numproc_in_check.ne.nsubdomains) then
    if(masterproc) print *,'number of processors is not equal to nsubdomains!',&
             '  ntasks=',ntasks,'   nsubdomains=',nsubdomains
    call task_abort_ORC() 
  endif
  print *,'barrier start',ntasks,nsubdomains      
  call task_barrier_ORC()
print *,'barrier finish'
  call task_ranks_ORC()
print *,'rank finish'        
end if ! nsubdomains.eq.1

#ifndef CRM
print *,'CaseName start'
do itasks=0,nsubdomains-1
   !call task_barrier()
   if(itasks.eq.rank) then
    open(8,file='./CaseName',status='old',form='formatted')
    read(8,'(a)') case
    close (8)
   endif
end do
print *,'CaseName end'
#endif  /*CRM*/

masterproc = rank.eq.0

#ifndef CRM
print *,'MPI tasks start'
if(masterproc) print *,'number of MPI tasks:',ntasks
print *,'MPI tasks end'
#endif /*CRM*/
	

end

end module crmx_task_init_mpi
