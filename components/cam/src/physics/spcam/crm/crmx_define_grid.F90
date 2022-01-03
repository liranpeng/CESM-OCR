
subroutine crm_define_grid
	
use crmx_vars

implicit none
	
nx = nx_gl/nsubdomains_x
ny = ny_gl/nsubdomains_y
nz = nz_gl+1
nzm = nz-1

nsubdomains = nsubdomains_x * nsubdomains_y

if(nsubdomains.eq.1) then
  dompi = .false.
else
  dompi = .true.
endif

RUN3D = ny_gl.gt.1
RUN2D = .not.RUN3D

nxp1 = nx + 1
nyp1 = ny + 1 * YES3D
nxp2 = nx + 2
nyp2 = ny + 2 * YES3D
nxp3 = nx + 3
nyp3 = ny + 3 * YES3D
nxp4 = nx + 4
nyp4 = ny + 4 * YES3D

dimx1_u = -1                !!-1        -1        -1        -1
dimx2_u = nxp3              !!nxp3      nxp3      nxp3      nxp3
dimy1_u = 1-(2+NADV)*YES3D  !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-2*YES3D
dimy2_u = nyp2+NADV         !!nyp5      nyp4      nyp3      nyp2
dimx1_v = -1-NADV           !!-4        -3        -2        -1
dimx2_v = nxp2+NADV         !!nxp5      nxp4      nxp3      nxp2
dimy1_v = 1-2*YES3D         !!1-2*YES3D 1-2*YES3D 1-2*YES3D 1-2*YES3D
dimy2_v = nyp3              !!nyp3       nyp3      nyp3      nyp3
dimx1_w = -1-NADV           !!-4        -3        -2        -1
dimx2_w = nxp2+NADV         !!nxp5      nxp4      nxp3      nxp2
dimy1_w = 1-(2+NADV)*YES3D  !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-2*YES3D
dimy2_w = nyp2+NADV         !!nyp5      nyp4      nyp3      nyp2
dimx1_s = -2-NADVS          !!-4        -3        -2        -2
dimx2_s = nxp3+NADVS        !!nxp5      nxp4      nxp3      nxp3
dimy1_s = 1-(3+NADVS)*YES3D !!1-5*YES3D 1-4*YES3D 1-3*YES3D 1-3*YES3D
dimy2_s = nyp3+NADVS        !!nyp5      nyp4      nyp3      nyp3

dimx1_d=0
dimx2_d=nxp1
dimy1_d=1-YES3D
dimy2_d=nyp1

ncols = nx*ny

end
