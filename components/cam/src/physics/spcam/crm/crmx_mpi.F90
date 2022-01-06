module crmx_mpi

implicit none

integer :: numproc_global,myrank_global
integer :: crm_comm, numproc_crm, myrank_crm
integer :: crm_comm_in,crm_comm_in0, numproc_crm_in,myrank_crm_in
integer :: crm_comm_color,ORC_count
end module crmx_mpi
