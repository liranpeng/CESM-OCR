subroutine task_init_ORC(ncomm,npro,ntask)
		
! crm_comm_in,numproc_crm_in,myrank_crm_in
!   Check things, initialize multitasking:

use crmx_grid
implicit none
integer, intent(in) :: ncomm,npro,ntask
integer itasks,ntasks
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
ntasks = ntask
rank=npro
if(nsubdomains.eq.1) then

  rank =0
  ntasks = 1
  dompi = .false.
  print *,'Something is wrong here',npro,ntask
else
  print *,'Enter task start ORC',ncomm,npro,ntask
  call task_start_ORC(rank, ntasks)

  dompi = .true.

  !call systemf('hostname')

  if(ntasks.ne.nsubdomains) then
    if(masterproc) print *,'number of processors is not equal to nsubdomains!',&
             '  ntasks=',ntasks,'   nsubdomains=',nsubdomains
    call task_abort_ORC() 
  endif
        
  call task_barrier_ORC()

  call task_ranks_ORC()
        
end if ! nsubdomains.eq.1

#ifndef CRM
do itasks=0,nsubdomains-1
   call task_barrier()
   if(itasks.eq.rank) then
    open(8,file='./CaseName',status='old',form='formatted')
    read(8,'(a)') case
    close (8)
   endif
end do
#endif  /*CRM*/

masterproc = rank.eq.0

#ifndef CRM
if(masterproc) print *,'number of MPI tasks:',ntasks
#endif /*CRM*/
	

end
