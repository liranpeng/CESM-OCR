
subroutine buoyancy(printflag)

use crmx_vars
use crmx_params
implicit none
	
integer i,j,k,kb,printflag
real betu, betd

if(docolumn) return

do k=2,nzm	
 kb=k-1
 betu=adz(kb)/(adz(k)+adz(kb))
 betd=adz(k)/(adz(k)+adz(kb))
 do j=1,ny
  do i=1,nx
if (printflag.eq.1) then
print*,'Liran dwdt0',i,j,nx,ny,k,na,dwdt(i,j,k,na)
end if
if (printflag.eq.2) then
print*,'Liran dwdt0 org',i,j,nx,ny,k,na,dwdt(i,j,k,na)
end if

   dwdt(i,j,k,na)=dwdt(i,j,k,na) +  &
      bet(k)*betu* &
     ( tabs0(k)*(epsv*(qv(i,j,k)-qv0(k))-(qcl(i,j,k)+qci(i,j,k)-qn0(k)+qpl(i,j,k)+qpi(i,j,k)-qp0(k))) &
       +(tabs(i,j,k)-tabs0(k))*(1.+epsv*qv0(k)-qn0(k)-qp0(k)) ) &
    + bet(kb)*betd* &
     ( tabs0(kb)*(epsv*(qv(i,j,kb)-qv0(kb))-(qcl(i,j,kb)+qci(i,j,kb)-qn0(kb)+qpl(i,j,kb)+qpi(i,j,kb)-qp0(kb))) &
       +(tabs(i,j,kb)-tabs0(kb))*(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb)) )  
if (printflag.eq.1) then
print*,'Liran dwdt1',dwdt(i,j,k,na),tabs0(k)*(epsv*(qv(i,j,k)-qv0(k))-(qcl(i,j,k)+qci(i,j,k)-qn0(k)+qpl(i,j,k)+qpi(i,j,k)-qp0(k))) 
print*,'Liran dwdt2',dwdt(i,j,k,na),(tabs(i,j,k)-tabs0(k))*(1.+epsv*qv0(k)-qn0(k)-qp0(k))
print*,'Liran dwdt21',dwdt(i,j,k,na),tabs(i,j,k),tabs0(k)
print*,'Liran dwdt22',dwdt(i,j,k,na),qv0(k),qn0(k),qp0(k)
print*,'Liran dwdt3',dwdt(i,j,k,na),tabs0(kb)*(epsv*(qv(i,j,kb)-qv0(kb))-(qcl(i,j,kb)+qci(i,j,kb)-qn0(kb)+qpl(i,j,kb)+qpi(i,j,kb)-qp0(kb))) 
print*,'Liran dwdt4',dwdt(i,j,k,na),(tabs(i,j,kb)-tabs0(kb))*(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb))
print*,'Liran dwdt5',dwdt(i,j,k,na),(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb))
print*,'Liran dwdt6',dwdt(i,j,k,na),tabs(i,j,kb),tabs0(kb),(tabs(i,j,kb)-tabs0(kb))
end if
if (printflag.eq.2) then
print*,'Liran dwdt1 org',dwdt(i,j,k,na),tabs0(k)*(epsv*(qv(i,j,k)-qv0(k))-(qcl(i,j,k)+qci(i,j,k)-qn0(k)+qpl(i,j,k)+qpi(i,j,k)-qp0(k)))
print*,'Liran dwdt2 org',dwdt(i,j,k,na),(tabs(i,j,k)-tabs0(k))*(1.+epsv*qv0(k)-qn0(k)-qp0(k))
print*,'Liran dwdt21 org',dwdt(i,j,k,na),tabs(i,j,k),tabs0(k)
print*,'Liran dwdt22 org',dwdt(i,j,k,na),qv0(k),qn0(k),qp0(k)
print*,'Liran dwdt3 org',dwdt(i,j,k,na),tabs0(kb)*(epsv*(qv(i,j,kb)-qv0(kb))-(qcl(i,j,kb)+qci(i,j,kb)-qn0(kb)+qpl(i,j,kb)+qpi(i,j,kb)-qp0(kb)))
print*,'Liran dwdt4 org',dwdt(i,j,k,na),(tabs(i,j,kb)-tabs0(kb))*(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb))
print*,'Liran dwdt5 org',dwdt(i,j,k,na),(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb))
print*,'Liran dwdt6 org',dwdt(i,j,k,na),tabs(i,j,kb),tabs0(kb),(tabs(i,j,kb)-tabs0(kb))
end if

  end do ! i
 end do ! j
end do ! k

end subroutine buoyancy


