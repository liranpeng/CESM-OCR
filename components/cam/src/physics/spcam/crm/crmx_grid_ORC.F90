module crmx_grid_ORC

use crmx_domain
use crmx_advection, only: NADV, NADVS

implicit none

character(6) :: version = '6.10.4'
character(8) :: version_date = 'Feb 2013'
        
integer :: nx ,ny ,nz ,nzm 

integer :: nsubdomains 

logical :: RUN3D 
logical :: RUN2D 

integer :: nxp1 ,nyp1 ,nxp2 ,nyp2 ,nxp3 ,nyp3 ,nxp4 ,nyp4 

integer :: dimx1_u ,dimx2_u ,dimy1_u ,dimy2_u 
integer :: dimx1_v ,dimx2_v ,dimy1_v ,dimy2_v
integer :: dimx1_w ,dimx2_w ,dimy1_w ,dimy2_w 
integer :: dimx1_s ,dimx2_s ,dimy1_s ,dimy2_s 
integer :: dimx1_d, dimx2_d, dimy1_d, dimy2_d

integer :: ncols 
integer, parameter :: nadams = 3

! Vertical grid parameters:
real, allocatable, dimension(:)  :: z
real, allocatable, dimension(:)  :: pres
real, allocatable, dimension(:)  :: zi
real, allocatable, dimension(:)  :: presi
real, allocatable, dimension(:)  :: adz
real, allocatable, dimension(:)  :: adzw

real pres0      ! Reference surface pressure, Pa

integer:: nstep =0! current number of performed time steps 
integer  ncycle  ! number of subcycles over the dynamical timestep
integer icycle  ! current subcycle 
integer:: na=1, nb=2, nc=3 ! indeces for swapping the rhs arrays for AB scheme
real at, bt, ct ! coefficients for the Adams-Bashforth scheme 
real dtn	! current dynamical timestep (can be smaller than dt)
real dt3(3) 	! dynamical timesteps for three most recent time steps
real(kind=selected_real_kind(12)):: time=0.	! current time in sec.
real day	! current day (including fraction)
real dtfactor   ! dtn/dt
        
!  MPI staff:     
integer rank   ! rank of the current subdomain task (default 0) 
integer ranknn ! rank of the "northern" subdomain task
integer rankss ! rank of the "southern" subdomain task
integer rankee ! rank of the "eastern"  subdomain task
integer rankww ! rank of the "western"  subdomain task
integer rankne ! rank of the "north-eastern" subdomain task
integer ranknw ! rank of the "north-western" subdomain task
integer rankse ! rank of the "south-eastern" subdomain task
integer ranksw ! rank of the "south-western" subdomain task
logical dompi  ! logical switch to do multitasking
logical masterproc ! .true. if rank.eq.0 
	
character(80) case   ! id-string to identify a case-name(set in CaseName file)

logical dostatis     ! flag to permit the gathering of statistics
logical dostatisrad  ! flag to permit the gathering of radiation statistics
integer nstatis	! the interval between substeps to compute statistics

logical :: compute_reffc = .false. 
logical :: compute_reffi = .false. 

logical notopened2D  ! flag to see if the 2D output datafile is opened	
logical notopened3D  ! flag to see if the 3D output datafile is opened	
logical notopenedmom ! flag to see if the statistical moment file is opened

!-----------------------------------------
! Parameters controled by namelist PARAMETERS

real:: dx =0. 	! grid spacing in x direction
real:: dy =0.	! grid spacing in y direction
real:: dz =0.	! constant grid spacing in z direction (when dz_constant=.true.)
logical:: doconstdz = .false.  ! do constant vertical grid spacing set by dz

integer:: nstop =0   ! time step number to stop the integration
integer:: nelapse =999999999! time step number to elapse before stoping

real:: dt=0.	! dynamical timestep
real:: day0=0.	! starting day (including fraction)

integer:: nrad =1  ! frequency of calling the radiation routines
integer:: nprint =1000	! frequency of printing a listing (steps)
integer:: nrestart =0 ! switch to control starting/restarting of the model
integer:: nstat =1000	! the interval in time steps to compute statistics
integer:: nstatfrq =50 ! frequency of computing statistics 

logical:: restart_sep =.false.  ! write separate restart files for sub-domains
integer:: nrestart_skip =0 ! number of skips of writing restart (default 0)
logical:: output_sep =.false.   ! write separate 3D and 2D files for sub-domains

character(80):: caseid =''! id-string to identify a run	
character(80):: caseid_restart =''! id-string for branch restart file 
character(80):: case_restart =''! id-string for branch restart file 

logical:: doisccp = .false.
logical:: domodis = .false.
logical:: domisr = .false.
logical:: dosimfilesout = .false.

logical:: doSAMconditionals = .false. !core updraft,downdraft conditional statistics
logical:: dosatupdnconditionals = .false.!cloudy updrafts,downdrafts and cloud-free
logical:: doscamiopdata = .false.! initialize the case from a SCAM IOP netcdf input file
logical:: dozero_out_day0 = .false.
character(len=120):: iopfile=''
character(256):: rundatadir ='./RUNDATA' ! path to data directory

integer:: nsave3D =1000     ! frequency of writting 3D fields (steps)
integer:: nsave3Dstart =99999999! timestep to start writting 3D fields
integer:: nsave3Dend  =99999999 ! timestep to end writting 3D fields
logical:: save3Dbin =.false.   ! save 3D data in binary format(no 2-byte compression)
logical:: save3Dsep =.false.   ! use separate file for each time point for2-model
real   :: qnsave3D =0.    !threshold manimum cloud water(kg/kg) to save 3D fields
logical:: dogzip3D =.false.    ! gzip compress a 3D output file   
logical:: rad3Dout = .false. ! output additional 3D radiation foelds (like reff)

integer:: nsave2D =1000     ! frequency of writting 2D fields (steps)
integer:: nsave2Dstart =99999999! timestep to start writting 2D fields
integer:: nsave2Dend =99999999  ! timestep to end writting 2D fields
logical:: save2Dbin =.false.   ! save 2D data in binary format, rather than compressed
logical:: save2Dsep =.false.   ! write separate file for each time point for 2D output
logical:: save2Davg =.false.   ! flag to time-average 2D output fields (default .false.)
logical:: dogzip2D =.false.    ! gzip compress a 2D output file if save2Dsep=.true.   

integer:: nstatmom =1000! frequency of writting statistical moment fields (steps)
integer:: nstatmomstart =99999999! timestep to start writting statistical moment fields
integer:: nstatmomend =99999999  ! timestep to end writting statistical moment fields
logical:: savemomsep =.false.! use one file with stat moments for each time point
logical:: savemombin =.false.! save statistical moment data in binary format

integer:: nmovie =1000! frequency of writting movie fields (steps)
integer:: nmoviestart =99999999! timestep to start writting statistical moment fields
integer:: nmovieend =99999999  ! timestep to end writting statistical moment fields

logical :: isInitialized_scamiopdata = .false.
logical :: wgls_holds_omega = .false.

!-----------------------------------------
end module crmx_grid_ORC
