
	subroutine task_start_ORC(ncomm_in,numproc_in,myrank_in)
	integer rank,numtasks
        integer numproc_global,myrank_global,ierr
        integer myrank_crm,numproc_crm,crm_comm	

        include 'mpif.h'
        ! get information on MPI_COMM_WORLD
        call mpi_comm_size(ncomm_in, numproc_in_check, ierr)
        call mpi_comm_rank(ncomm_in, myrank_in_check, ierr)
! Split MPI_COMM_WORLD into two communicators:
!  - global_comm: used by CIME for CESM-related communication
!  - crm_comm: used by CRM routines to receive data from CESM
!         for separate-executable CRM computations and
!         to send the resulting tendencies back to CESM.
!call mpi_comm_split(MPI_COMM_WORLD, 1, 0, crm_comm, ierr)

! get information on MPI_COMM_WORLD

	print*, 'Liran Check MPI start1',numproc_in,numproc_in_check
        print*, 'Liran Check MPI start2',myrank_in,myrank_in_check
	!stop
	end
	
!----------------------------------------------------------------------
	
	subroutine task_abort_ORC()
        use crmx_grid, only: dompi, nstep,nstop
        include 'mpif.h'
        integer ierr, rc

        if(dompi) then
!          call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
!bloss          call MPI_FINALIZE(ierr)
        endif

!bloss: call task_stop instead
!bloss        call exit(999) ! to avolid resubmission when finished
        call task_stop_ORC()
	end

!----------------------------------------------------------------------

        subroutine task_stop_ORC()

        use crmx_grid, only: dompi,nstep,nstop,nelapse
        include 'mpif.h'
        integer ierr

        if(dompi) then
          call MPI_FINALIZE(ierr)
        endif

        if(nstep.ge.nstop) then
           call exit(9) ! avoid resubmission when finished
        elseif(nelapse.eq.0) then
           call exit(0) !bloss: clean exit condition for restart
        else
           call exit(1) !bloss: avoid resubmission if ending in error
        end if

        end

!----------------------------------------------------------------------
	
	subroutine task_finish_ORC()
	print*,'program is finished...'
	stop
	end

!----------------------------------------------------------------------
        subroutine task_barrier_ORC()

        use crmx_grid, only: dompi
        implicit none
        include 'mpif.h'
        integer ierr

        if(dompi) then
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end if

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_float_ORC(rank_from,buffer,length)
        implicit none
        integer rank_from       ! broadcasting task's rank
        real buffer(*)          ! buffer of data
        integer length          ! buffers' length
	print*, 'MPIsndf call from a single task program! Exiting...'
        stop
        end

!----------------------------------------------------------------------

	subroutine task_send_float_ORC(rank_to,buffer,length,tag,request)
	implicit none
	integer rank_to		! receiving task's rank
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request		! request id
	print*, 'MPIsndf call from a single task program! Exiting...'
	stop
	end

!----------------------------------------------------------------------

	subroutine task_send_integer_ORC(rank_to,buffer,length,tag,request)
	implicit none
	integer rank_to		! receiving task's rank
	integer buffer(*)	! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request
	print*, 'MPIsndi call from a single task program! Exiting...'
	stop
	end
	
!----------------------------------------------------------------------

	subroutine task_send_character_ORC(rank_to,buffer,length,tag,request)
	implicit none
	integer rank_to		! receiving task's rank
	character*1 buffer(*)	! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request
	print*, 'MPIsndi call from a single task program! Exiting...'
	stop
	end
	
!----------------------------------------------------------------------

        subroutine task_receive_float_ORC(buffer,length,request)
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer request
	print*, 'MPIrcvf call from a single task program! Exiting...'
	stop
	end

!----------------------------------------------------------------------

        subroutine task_receive_charcater_ORC(buffer,length,request)
	character*1 buffer(*)	! buffer of data
	integer length		! buffers' length
	integer request
	print*, 'MPIrcvi call from a single task program! Exiting...'
	stop
	end

!----------------------------------------------------------------------

        subroutine task_receive_integer_ORC(buffer,length,request)
	integer buffer(*)	! buffer of data
	integer length		! buffers' length
	integer request
	print*, 'MPIrcvi call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_bsend_float_ORC(rank_to,buffer,length,tag)
        integer rank_to         ! receiving task's rank
        real buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        print*, 'MPI call from a single task program! Exiting...'  
	stop
        return
        end

!----------------------------------------------------------------------
        subroutine task_wait_ORC(request,rank,tag)
	integer request
	integer rank, tag
	return
	end

!----------------------------------------------------------------------
        
        subroutine task_waitall_ORC(count,reqs,ranks,tags)
 	integer count,reqs(count)
	integer ranks(count),tags(count)
	return
	end

!----------------------------------------------------------------------
        subroutine task_test_ORC(request,flag,rank,tag)
	integer request
	integer rank, tag
	logical flag
	print*, 'MPItst call from a single task program! Exiting...'
	stop
	end

!----------------------------------------------------------------------

        subroutine task_sum_real_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end

!----------------------------------------------------------------------

        subroutine task_sum_real8_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_sum_integer_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_max_real_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	return
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_max_integer_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_min_real_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_min_integer_ORC(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
!----------------------------------------------------------------------

        subroutine task_receive_character_ORC(buffer,length,request)
        character*1 buffer(*)   ! buffer of data
        integer length          ! buffers' length
        integer request
	print*, 'MPI call from a single task program! Exiting...'
        stop
        end
!----------------------------------------------------------------------
       subroutine task_rank_to_index_ORC(rank,i,j)
        integer rank, i, j
        i=0
        j=0
        end
!----------------------------------------------------------------------
       subroutine task_bound_duvdt_ORC()
       return
       end
!----------------------------------------------------------------------
       subroutine task_boundaries_ORC(flag)
       integer flag
       end


