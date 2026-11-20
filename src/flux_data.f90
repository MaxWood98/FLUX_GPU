!flux 2d data module 
!max wood
!version : 2.0.0
!updated : 15-11-26

!module 
module flux_data

!data types
use ISO_FORTRAN_ENV, only: in16=>int16
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64
use ISO_FORTRAN_ENV, only: sp=>real32 
use ISO_FORTRAN_ENV, only: dp=>real64

!imports 
use io_utilities

!options type
type options_type
    logical :: cdisplay
    logical, dimension(:), allocatable :: rk_dissipation
    integer(in64) :: niter_max,rk_niter,num_threads,flux_method_idx
    real(dp) :: cfl,aoadeg,aoarad,gamma,R,k2,k4,residual_convtol,outflow_pratio
    real(dp) :: machinf,tinf,rhoinf,cinf,pinf,velinf,uinf,vinf,p0inf,rho0inf,t0inf
    character(len=:), allocatable :: meshpath,meshname,mode,flux_method,adjoint_objective
end type options_type 

!mesh type
type mesh_type
    integer(in32) :: nvertex,nedge,ncell 
    integer(in32) :: cell_max_nedge,cell_min_nedge,ncell_edges
    integer(in32), dimension(:), allocatable :: edge_order,cell_order
    integer(in32), dimension(:), allocatable :: edges_v1,edges_v2,edges_c1,edges_c2
    integer(in32), dimension(:), allocatable :: cells_edges_start,cells_edges
    real(dp), dimension(:), allocatable :: cells_edges_sign
    real(dp), dimension(:), allocatable :: edges_length,edges_nx,edges_ny,edges_dx,edges_dy,cells_volume
    real(dp), dimension(:,:), allocatable :: vertices
end type mesh_type

!methods
contains

!set default options ===============
subroutine set_default_options(options)

!variables - inout 
type(options_type) :: options 

!general
options%cdisplay = .true.
options%meshpath = 'grid'

!flow conditions
options%aoadeg = 0.0d0 
options%machinf = 0.5d0 
options%gamma = 1.4d0 
options%R = 287.058d0
options%tinf = 288.0d0
options%rhoinf = 1.225d0 
options%outflow_pratio = 0.8d0

!solve options
options%niter_max = 10000
options%cfl = 2.0
options%k2 = 0.5d0 
options%k4 = 0.1d0
options%num_threads = 1
options%residual_convtol = -8.0 
options%flux_method = 'vanleer'
options%flux_method_idx = 0

!adjoint options
options%adjoint_objective = 'cd'

!rk options
options%rk_niter = 4
allocate(options%rk_dissipation(options%rk_niter))
options%rk_dissipation(1) = .false.
options%rk_dissipation(2) = .false.
options%rk_dissipation(3) = .false.
options%rk_dissipation(4) = .false.
end subroutine set_default_options



end module flux_data