module crmx_task_util_mpi
use crmx_mpi
use crmx_grid

implicit none
private
save

public :: task_start_ORC
public :: task_abort_ORC
public :: task_stop_ORC
public :: task_finish_ORC
public :: task_bcast_float_ORC
public :: task_send_float_ORC
public :: task_send_integer_ORC
public :: task_send_character_ORC
public :: task_receive_float_ORC
public :: task_receive_charcater_ORC
public :: task_receive_integer_ORC
public :: task_bsend_float_ORC
public :: task_wait_ORC
public :: task_waitall_ORC
public :: task_test_ORC
public :: task_sum_real_ORC
public :: task_sum_real8_ORC
public :: task_sum_integer_ORC
public :: task_max_real_ORC
public :: task_max_integer_ORC
public :: task_min_real_ORC
public :: task_min_integer_ORC
public :: task_receive_character_ORC
public :: task_rank_to_index_ORC
public :: task_bound_duvdt_ORC
public :: task_boundaries_ORC
public :: task_barrier_ORC

contains

	subroutine task_start_ORC(rank,numtasks)
	integer numproc_in_check,myrank_in_check,crm_comm_check
        integer numproc_global_check,myrank_global_check
        integer crm_comm_in_check,numproc_crm_check,myrank_crm_check,ierr
        integer rank,numtasks
        include 'mpif.h'
        print *,'Liran check984'
        call mpi_comm_size(MPI_COMM_WORLD, numproc_global_check, ierr)
        call mpi_comm_rank(MPI_COMM_WORLD, myrank_global_check, ierr)
        print *,'Liran check985',numproc_global_check,numproc_global
        print *,'Liran check986',myrank_global_check,myrank_global
! get information on MPI_COMM_WORLD
        call mpi_comm_size(crm_comm, numproc_crm_check, ierr)
        call mpi_comm_rank(crm_comm, myrank_crm_check, ierr)
        print *,'Liran check987',numproc_crm_check,numproc_crm
        print *,'Liran check988',myrank_crm_check,myrank_crm
        !crm_comm_color = int((myrank_global_check-50)/2)
!print *,'Liran check987',crm_comm_color
        !call mpi_comm_split(, crm_comm_color, 0, crm_comm_in_check, ierr)
        print *,'Liran check989',crm_comm_color,ierr,crm_comm_in
        ! get information on MPI_COMM_WORLD
        call mpi_comm_size(crm_comm_in, numproc_in_check, ierr)
        call mpi_comm_rank(crm_comm_in, myrank_in_check, ierr)
        print *,'Liran check999',numproc_in_check,myrank_in_check
        rank = myrank_in_check
        numtasks = numproc_in_check
! Split MPI_COMM_WORLD into two communicators:
!  - global_comm: used by CIME for CESM-related communication
!  - crm_comm: used by CRM routines to receive data from CESM
!         for separate-executable CRM computations and
!         to send the resulting tendencies back to CESM.
!call mpi_comm_split(MPI_COMM_WORLD, 1, 0, crm_comm, ierr)

! get information on MPI_COMM_WORLD

        print*, 'Liran Check MPI start',numproc_in_check,myrank_in_check
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
        use crmx_mpi
        use crmx_grid, only: dompi
        implicit none
        include 'mpif.h'
        integer ierr

        if(dompi) then
          call MPI_BARRIER(crm_comm_in,ierr)
        end if

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_float_ORC(rank_from,buffer,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'

        integer rank_from       ! broadcasting task's rank
        real buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer ierr, real_size

        if(sizeof(buffer(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_BCAST(buffer,length,real_size,rank_from,crm_comm_in,ierr)

        return  
        end

!----------------------------------------------------------------------

	subroutine task_send_float_ORC(rank_to,buffer,length,tag,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        integer rank_to         ! receiving task's rank
        real buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer request         ! request id
        integer ierr, real_size

        if(sizeof(buffer(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_ISEND(buffer,length,real_size,rank_to,tag,crm_comm_in,request,ierr)

        
        return
	end

!----------------------------------------------------------------------

	subroutine task_send_integer_ORC(rank_to,buffer,length,tag,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        integer rank_to         ! receiving task's rank
        integer buffer(*)       ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer request
        integer ierr

        call MPI_ISEND(buffer,length,MPI_INTEGER,rank_to,tag, &
                                        crm_comm_in,request,ierr)

        return
	end
	
!----------------------------------------------------------------------

	subroutine task_send_character_ORC(rank_to,buffer,length,tag,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        integer rank_to         ! receiving task's rank
        character*1 buffer(*)   ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer request
        integer ierr

        call MPI_ISEND(buffer,length,MPI_CHARACTER,rank_to,tag, &
                                        crm_comm_in,request,ierr)

        return
	end
	
!----------------------------------------------------------------------

        subroutine task_receive_float_ORC(buffer,length,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        real buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer request
        integer ierr, real_size

        if(sizeof(buffer(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_IRECV(buffer,length,real_size,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,crm_comm_in,request,ierr)

        return
	end

!----------------------------------------------------------------------

        subroutine task_receive_charcater_ORC(buffer,length,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        character*1 buffer(*)   ! buffer of data
        integer length          ! buffers' length
        integer request
        integer ierr

        call MPI_IRECV(buffer,length,MPI_CHARACTER,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,crm_comm_in,request,ierr)

        return
	end

!----------------------------------------------------------------------

        subroutine task_receive_integer_ORC(buffer,length,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        integer buffer(*)       ! buffer of data
        integer length          ! buffers' length
        integer request
        integer ierr

        call MPI_IRECV(buffer,length,MPI_INTEGER,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,crm_comm_in,request,ierr)

        return
	end
!----------------------------------------------------------------------

        subroutine task_bsend_float_ORC(rank_to,buffer,length,tag)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        integer rank_to         ! receiving task's rank
        real buffer(*)          ! buffer of data
        integer length          ! buffers' length
        integer tag             ! tag of the message
        integer ierr, real_size

        if(sizeof(buffer(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if
         call MPI_SEND(buffer,length,real_size,rank_to,tag,crm_comm_in,ierr)
        
        return
        end

!----------------------------------------------------------------------
        subroutine task_wait_ORC(request,rank,tag)
        use crmx_mpi
        implicit none
        include 'mpif.h'
        integer status(MPI_STATUS_SIZE),request
        integer rank, tag
        integer ierr
        call MPI_WAIT(request,status,ierr)
        rank = status(MPI_SOURCE)
        tag = status(MPI_TAG)

        return
	end

!----------------------------------------------------------------------
        
        subroutine task_waitall_ORC(count,reqs,ranks,tags)
        use crmx_mpi
        use crmx_grid, only: dompi
        implicit none
        include 'mpif.h'
        integer count,reqs(count)
        integer stats(MPI_STATUS_SIZE,1000),ranks(count),tags(count)
        integer ierr, i
        if(dompi) then
        call MPI_WAITALL(count,reqs,stats,ierr)
        if(count.gt.1000) then
            print*,'task_waitall: count > 1000 !'
            call task_abort()
        end if
        do i = 1,count
          ranks(i) = stats(MPI_SOURCE,i)
          tags(i) = stats(MPI_TAG,i)
        end do
        end if

        return
	end

!----------------------------------------------------------------------

        subroutine task_test_ORC(request,flag,rank,tag)
        use crmx_mpi
        implicit none
        include 'mpif.h'
        integer rank, tag
        logical flag
        integer ierr
        integer status(MPI_STATUS_SIZE),request
        call MPI_TEST(request,flag,status,ierr)
        if(flag) then
          rank = status(MPI_SOURCE)
          tag = status(MPI_TAG)
        endif

        return
	end

!----------------------------------------------------------------------

        subroutine task_sum_real_ORC(buffer_in,buffer_out,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        real buffer_in(*)       ! buffer of data
        real buffer_out(*)      ! buffer of data
        integer length          ! buffers' length
        integer ierr, real_size

        if(sizeof(buffer_in(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                           real_size,MPI_SUM,crm_comm_in,ierr)

        return
	end

!----------------------------------------------------------------------

        subroutine task_sum_real8_ORC(buffer_in,buffer_out,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        real(8) buffer_in(*)    ! buffer of data
        real(8) buffer_out(*)   ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                         MPI_REAL8,MPI_SUM, crm_comm_in,ierr)

        return
	end
!----------------------------------------------------------------------

        subroutine task_sum_integer_ORC(buffer_in,buffer_out,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        integer buffer_in(*)    ! buffer of data
        integer buffer_out(*)   ! buffer of data
        integer length          ! buffers' length
        integer ierr

        call MPI_ALLREDUCE(buffer_in,buffer_out,length, &
                        MPI_INTEGER,MPI_SUM, crm_comm_in,ierr)

        return
	end
!----------------------------------------------------------------------

        subroutine task_max_real_ORC(buffer_in,buffer_out,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        real buffer_in(*)       ! buffer of data
        real buffer_out(*)      ! buffer of data
        integer length          ! buffers' length
        integer ierr, real_size

        if(sizeof(buffer_in(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_ALLREDUCE(buffer_in,buffer_out, &
                          length,real_size,MPI_MAX,crm_comm_in,ierr)

        return
	end
!----------------------------------------------------------------------

        subroutine task_max_integer_ORC(buffer_in,buffer_out,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
       
        integer buffer_in(*)    ! buffer of data
        integer buffer_out(*)   ! buffer of data
        integer length          ! buffers' length
        integer ierr,numproc_in_check,myrank_in_check
        integer myrank_global
        call mpi_comm_rank(MPI_COMM_WORLD, myrank_global, ierr)
        call mpi_comm_size(crm_comm_in, numproc_in_check, ierr)
        call mpi_comm_rank(crm_comm_in, myrank_in_check, ierr)
        print *,'Liran check9790',numproc_in_check,myrank_in_check,myrank_global
        !if (myrank_in_check.eq.1) then
          call MPI_ALLREDUCE(buffer_in,buffer_out, &
                          length,MPI_INTEGER,MPI_MAX,crm_comm_in,ierr)
        !end if
        print *,'Liran check9991'
        return
	end
!----------------------------------------------------------------------

        subroutine task_min_real_ORC(buffer_in,buffer_out,length)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        real buffer_in(*)       ! buffer of data
        real buffer_out(*)      ! buffer of data
        integer length          ! buffers' length
        integer ierr, real_size

        if(sizeof(buffer_in(1)).eq.4) then
         real_size=MPI_REAL
        else
         real_size=MPI_REAL8
        end if

        call MPI_ALLREDUCE(buffer_in,buffer_out, &
                            length,real_size,MPI_MIN,crm_comm_in,ierr)
        return
	end
!----------------------------------------------------------------------

        subroutine task_min_integer_ORC(buffer_in,buffer_out,length)
        use crmx_mpi	
        implicit none
        include 'mpif.h'        
        
        integer buffer_in(*)    ! buffer of data
        integer buffer_out(*)   ! buffer of data
        integer length          ! buffers' length
        integer ierr
        call MPI_ALLREDUCE(buffer_in,buffer_out, &
                  length,MPI_INTEGER,MPI_MIN,crm_comm_in,ierr)

        return
	end
!----------------------------------------------------------------------

        subroutine task_receive_character_ORC(buffer,length,request)
        use crmx_mpi
        implicit none
        include 'mpif.h'        
        
        character*1 buffer(*)   ! buffer of data
        integer length          ! buffers' length
        integer request
        integer ierr

        call MPI_IRECV(buffer,length,MPI_CHARACTER,MPI_ANY_SOURCE, &
                MPI_ANY_TAG,crm_comm_in,request,ierr)

        return
        end
!----------------------------------------------------------------------
       subroutine task_rank_to_index_ORC(rank,i,j)
        !   returns the pair of  beginning indeces for the subdomain on the  
        !   global grid given the subdomain's rank.

        use crmx_domain

        implicit none

        integer rank, i, j
                
        j = rank/nsubdomains_x
        i = rank - j*nsubdomains_x
        
        i = i * (nx_gl/nsubdomains_x)
        j = j * (ny_gl/nsubdomains_y)
        end
!----------------------------------------------------------------------
       subroutine task_bound_duvdt_ORC()

        !  These routine exchanges subdomain overlaping information 
        !  for horizontal velocity tendencies

        use crmx_vars
        implicit none

        integer i,j,k,n,m,tag(2),rf(2),bufflen,nsent
        real, allocatable, dimension(:) :: buff_send
        real, allocatable, dimension(:,:) :: buff_recv
        integer reqs_in(2)

         allocate ( buff_send(max(nx,ny)*nzm))      
         allocate ( buff_recv(max(nx,ny)*nzm,2)      )

        bufflen = max(nx,ny)*nzm
        ! Any messages to send?

        nsent=0
        if(rank.ne.rankww) nsent=nsent+1
        if(rank.ne.rankss) nsent=nsent+1

        ! Non-blocking receive first:

         do m =  1,nsent
            call task_receive_float_ORC(buff_recv(1,m),bufflen,reqs_in(m))
         end do

        ! Blocking send second:

         n=0
         do k=1,nzm
            do j=1,ny
             n=n+1
             buff_send(n) = dudt(1,j,k,na)
            end do
         end do

         if(rank.ne.rankww) then
               call task_bsend_float_ORC(rankww, buff_send,n,54)
         else
               do i=1,n
                  buff_recv(i,2) = buff_send(i)
               end do
               tag(2) = 54
               rf(2) = rankee
         end if

         if(RUN3D) then

          n=0
          do k=1,nzm
            do i=1,nx
              n=n+1
              buff_send(n) = dvdt(i,1,k,na)
            end do
          end do


          if(rank.ne.rankss) then
             call task_bsend_float_ORC(rankss, buff_send, n, 54)
          else
             do i=1,n
               buff_recv(i,2) = buff_send(i)
             end do
             tag(2) = 54
             rf(2) = ranknn
          end if

         endif


! Wait until messages are received::

          call task_waitall_ORC(nsent,reqs_in,rf,tag)

! Fill data:

          do m =  1,1+YES3D
            if(tag(m).ne.54) then
               print*,'MPI:Wrong message tag in task_bound_duvdt.'
               print*,'    expected 54  Received:',tag(m)
               call task_abort()
            endif
            if(rf(m).eq.rankee) then
               n=0
               do k=1,nzm
                 do j=1,ny
                   n = n+1
                   dudt(nxp1,j,k,na) = buff_recv(n,m)
                 end do
               end do
            elseif (rf(m).eq.ranknn) then
               n=0
               do k=1,nzm
                 do i=1,nx
                   n = n+1
                   dvdt(i,nyp1,k,na) = buff_recv(n,m)
                 end do
               end do
            end if
          end do
          end
!----------------------------------------------------------------------
       subroutine task_boundaries_ORC(flag)
!  These routines exchanges overlaping information for various
!  subdomains to
        !  be able to calculate various quantities in common /com3d/
        !  near the subdomain boundaries.

        use crmx_vars
        use crmx_microphysics
        use crmx_sgs
        use crmx_params, only: dotracers, dosgs
        use crmx_crmtracers
        implicit none

        integer flag,i

        if(flag.eq.0) then

         call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
         call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
         ! use w at the top level  - 0s anyway - to exchange the sst
         ! boundaries (for
         ! surface fluxes call
         w(1:nx,1:ny,nz) = sstxy(1:nx,1:ny)
         call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)   
         sstxy(0:nx,1-YES3D:ny) = w(0:nx,1-YES3D:ny,nz)
         w(0:nx+1,1-YES3D:ny+YES3D,nz) = 0. ! fill it back with 0s

        endif

        if(flag.eq.2) then

         call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV,1)
         call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3,2)
         call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV,3)       

         call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4)
         do i = 1,nsgs_fields
            if(dosgs.and.advect_sgs) &
             call task_exchange(sgs_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm, &
                                                                   3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+i)
         end do
         do i = 1,nmicro_fields
        !!$    if(   i.eq.index_water_vapor             &
        !!$     .or. docloud.and.flag_precip(i).ne.1    &
        !!$     .or. doprecip.and.flag_precip(i).eq.1 ) &
             call task_exchange(micro_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,&
                                        3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+i)
         end do
         if(dotracers) then
           do i=1,ntracers
             call task_exchange(tracer(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm, &
                           3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+i)
           end do
         end if


        endif


        if(flag.eq.3) then

         call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
         do i = 1,nsgs_fields
            if(dosgs.and.advect_sgs) &
             call task_exchange(sgs_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+i)
         end do
         do i = 1,nmicro_fields
        !!$    if(   i.eq.index_water_vapor             &
        !!$     .or. docloud.and.flag_precip(i).ne.1    &
        !!$     .or. doprecip.and.flag_precip(i).eq.1 ) &
             call task_exchange(micro_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,&
                                                                     1,1,1,1,4+nsgs_fields+nsgs_fields_diag+i)
         end do
         if(dotracers) then
           do i=1,ntracers
             call task_exchange(tracer(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm, &
                                      1,1,1,1,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+i)
           end do
         end if

        end if

        if(flag.eq.4) then

         do i = 1,nsgs_fields_diag
            if(dosgs) &
             call task_exchange(sgs_field_diag(:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nzm,&
                           1+dimx1_d,dimx2_d-nx,YES3D+dimy1_d,1-YES3D+dimy2_d-ny,4+nsgs_fields+i)
         end do

        end if



       end

end module crmx_task_util_mpi
