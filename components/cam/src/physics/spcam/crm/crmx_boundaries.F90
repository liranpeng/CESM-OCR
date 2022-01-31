
subroutine boundaries(flag)

use crmx_grid, only: dompi
use crmx_task_util_mpi
        	
implicit none
integer flag

!call t_startf ('boundaries')

if(dompi) then
  !call periodic(flag)
  call task_boundaries_ORC(flag)
else
  call periodic(flag)
end if

!call t_stopf ('boundaries')

end subroutine boundaries
