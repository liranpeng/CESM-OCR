
subroutine setperturb_ORC(iseed)

!  Random noise
!  This surboutine has been updated for SPCAM5 (Minghuai.Wang@pnnl.gov, April, 2012). 
!  Now the random generator is seeded based on the global column id, which gets rid
!  of the dependence of the SPCAM reulst on pcols. 

use crmx_vars
use crmx_sgs, only: setperturb_sgs
use crmx_task_util_mpi, only: task_sum_real_ORC, task_rank_to_index_ORC

implicit none

integer, intent(in) :: iseed

integer i,j,k
real rrr,ranf_
integer, allocatable :: rndm_seed(:)
integer :: rndm_seed_sz
real :: t02(nzm), t02b(nzm)
real :: tke02(nzm)

!bloss
integer :: it, jt
integer :: i_global, j_global

!call ranset_(30*rank)
call random_seed(size=rndm_seed_sz)
allocate(rndm_seed(rndm_seed_sz))

rndm_seed = iseed
call random_seed(put=rndm_seed)

!bloss(ORC): SGS perturbations are deterministic.  No special ORC treatment needed.
call setperturb_sgs(0)  ! set sgs fields

call task_rank_to_index_ORC(rank,it,jt)

t02 = 0.0
tke02 = 0.0
do k=1,nzm
 do j_global=1,ny_gl  ! Loop over all points on global CRM domain
  do i_global=1,nx_gl
    rrr=1.-2.*ranf_()

    if(k.le.5) then
      i = i_global - it  ! get local values of i and j
      j = j_global - jt
      if((j.ge.1).AND.(j.le.ny).AND.(i.ge.1).AND.(i.le.nx)) then
        t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k) ! only add perturbation within this subdomain
        t02(k) = t02(k) + t(i,j,k)/(nx*ny)
      end if
    endif
  end do ! i_global
 end do ! j_global
end do !bloss(ORC): Separate perturbation and energy conservation

if(dompi) then
  call task_sum_real_ORC(t02,t02b,nzm) ! sum across subdomains
  t02(:) = t02b(:)/float(nsubdomains) 
end if

! energy conservation +++mhwang (2012-06)
do k=1, nzm
 if(k.le.5) then
   do j=1, ny
    do i=1, nx
      t(i,j,k) = t(i,j,k) * t0(k)/t02(k)
    end do
  end do
 end if
end do

deallocate(rndm_seed)

end

