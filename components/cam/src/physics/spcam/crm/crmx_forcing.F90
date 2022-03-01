
subroutine forcing(printflag)
	
        use crmx_grid, only: dompi 
        use crmx_vars
        use crmx_params
        use crmx_microphysics, only: micro_field, index_water_vapor, total_water
        use crmx_task_util_mpi, only: task_sum_real_ORC
        implicit none

        real coef,qneg(nzm),qpoz(nzm), factor,buffer1(nzm,2),buffer2(nzm,2)
        integer i,j,k,nneg(nzm),printflag

        coef = 1./3600.

        do k=1,nzm

            qpoz(k) = 0.
            qneg(k) = 0.
            nneg(k) = 0

if (printflag.eq.1) then
print*,'Liran ttend org',k,ttend(k),utend(k),qtend(k)
end if
if (printflag.eq.2) then
print*,'Liran ttend',k,ttend(k),utend(k),qtend(k)
end if

            do j=1,ny
             do i=1,nx

              t(i,j,k)=t(i,j,k) + ttend(k) * dtn
              micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor) + qtend(k) * dtn
              if(micro_field(i,j,k,index_water_vapor).lt.0.) then
                   nneg(k) = nneg(k) + 1
                   qneg(k) = qneg(k) + micro_field(i,j,k,index_water_vapor)
              else
                   qpoz(k) = qpoz(k) + micro_field(i,j,k,index_water_vapor)
              end if
              dudt(i,j,k,na)=dudt(i,j,k,na) + utend(k)
              dvdt(i,j,k,na)=dvdt(i,j,k,na) + vtend(k)
             end do
            end do
        end do
        if(dompi) then
          buffer1(:,1) = qneg
          buffer1(:,2) = qpoz
          call task_sum_real_ORC(buffer1,buffer2,2*nzm)
          qneg(:)    = buffer2(:, 1)
          qpoz(:)    = buffer2(:, 2)
        end if
        do k=1,nzm    
            if(nneg(k).gt.0.and.qpoz(k)+qneg(k).gt.0.) then
             factor = 1. + qneg(k)/qpoz(k)
             do j=1,ny
              do i=1,nx
               micro_field(i,j,k,index_water_vapor) = max(0.,micro_field(i,j,k,index_water_vapor)*factor)
              end do
             end do
            end if

        end do

end 

