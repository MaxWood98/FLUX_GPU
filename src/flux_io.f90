!flux 2d io module 
!max wood
!version : 2.0.0
!updated : 15-11-26

!module 
module flux_io
use flux_data
contains 

!read command arguments subroutine ===============
subroutine get_command_arguments(options)
implicit none

!variables - import
type(options_type) :: options 

!variables - local 
integer(in64) :: nargs

!check and process supplied command arguments 
nargs = command_argument_count()

!set mode
if (nargs == 0) then !default to solve mode
    options%mode = 'solve'
else !set mode from first argument 
    options%mode = get_command_argument_n_str(1)
end if 
end subroutine get_command_arguments

!read options ===============
subroutine read_options(options)
implicit none

!variables - inout
type(options_type) :: options 

!variables - Local
integer(in32) :: fh
integer(in32) :: ii
integer(in64) :: itemtemp
integer(in64), dimension(:), allocatable :: dissflags

!check if file exists 
if (.NOT. file_exists('flux_options')) then 
    write(*,'(A)') '** cannot locate options file: '//trim('flux_options')
    stop 
end if 

!open file 
open(newunit=fh,file='flux_options')

!set options 
call set_log_opt(options%cdisplay,fh,'console_display')
call set_str_opt(options%meshpath,fh,'meshpath')
call set_real_opt(options%aoadeg,fh,'aoadeg')
call set_real_opt(options%machinf,fh,'machinf')
call set_real_opt(options%gamma,fh,'gamma')
call set_real_opt(options%R,fh,'R')
call set_real_opt(options%tinf,fh,'tinf')
call set_real_opt(options%rhoinf,fh,'rhoinf')
call set_real_opt(options%outflow_pratio,fh,'outflow_pratio')
itemtemp = -1 
call set_int_opt(itemtemp,fh,'niter_max')
if (itemtemp .GE. 0) then 
    options%niter_max = int(itemtemp,in64)
end if 
call set_real_opt(options%cfl,fh,'cfl')
call set_real_opt(options%k2,fh,'k2')
call set_real_opt(options%k4,fh,'k4')
call set_str_opt(options%flux_method,fh,'flux_method')
if (options%flux_method == 'vanleer') then 
    options%flux_method_idx = 0
elseif (options%flux_method == 'ausm') then 
    options%flux_method_idx = 1
elseif (options%flux_method == 'ausm+up') then 
    options%flux_method_idx = 2
elseif (options%flux_method == 'roe') then 
    options%flux_method_idx = 3
elseif (options%flux_method == 'jameson') then 
    options%flux_method_idx = 4
end if 
call set_str_opt(options%adjoint_objective,fh,'adjoint_objective')
itemtemp = -1 
call set_int_opt(itemtemp,fh,'num_threads')
if (itemtemp .GE. 0) then 
    options%num_threads = int(itemtemp,in64)
end if 
call set_real_opt(options%residual_convtol,fh,'residual_convtol')
call set_int_opt_arr1(dissflags,fh,'RKstagediss_eval')
if (allocated(dissflags)) then 
    if (allocated(options%rk_dissipation)) then 
        deallocate(options%rk_dissipation)
    end if 
    allocate(options%rk_dissipation(size(dissflags)))
    do ii=1,size(dissflags)
        if (dissflags(ii) == 1) then 
            options%rk_dissipation(ii) = .true.
        else
            options%rk_dissipation(ii) = .false.
        end if 
    end do 
    deallocate(dissflags)
end if 

!close file 
close(fh)
end subroutine read_options

!read mesh ===============
subroutine read_mesh(mesh,filename)
implicit none

!variables - inout
character(*) :: filename
type(mesh_type) :: mesh 

!variables - local
integer(in64) :: ii
integer(in64) :: iostatus
character(len=100) :: rtemp 

!check if file exists 
if (.NOT. file_exists(trim(filename))) then 
    write(*,'(A)') '** cannot locate mesh file: '//trim(filename)
    stop 
end if 

!open mesh file
open(11,file=filename) 

!read number of cells
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'ncell') then 
        read(rtemp(9:len_trim(rtemp)),*) mesh%ncell
    end if 
end do 
rewind(11)

!read edges
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'nedge') then 
        read(rtemp(8:len_trim(rtemp)),*) mesh%nedge
        allocate(mesh%edges_v1(mesh%nedge))
        allocate(mesh%edges_v2(mesh%nedge))
        allocate(mesh%edges_c1(mesh%nedge))
        allocate(mesh%edges_c2(mesh%nedge))
        do ii=1,mesh%nedge
            read(11,*) mesh%edges_v1(ii),mesh%edges_v2(ii),mesh%edges_c1(ii),mesh%edges_c2(ii)
        end do 
        exit
    end if 
end do 
rewind(11)

!read vertices
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:7) == 'nvertex') then 
        read(rtemp(11:len_trim(rtemp)),*) mesh%nvertex
        allocate(mesh%vertices(mesh%nvertex,2))
        do ii=1,mesh%nvertex
            read(11,*) mesh%vertices(ii,:)
        end do 
        exit 
    end if 
end do 

!close mesh file 
close(11)
end subroutine read_mesh

end module flux_io
