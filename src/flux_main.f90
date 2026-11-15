!flux 2d - a cell based euler solver
!max wood
!version : 2.0.0
!updated : 15-11-26

!program 
program flux2d
use flux_io
use flux_mesh
implicit none 

!variables
type(mesh_type) :: mesh
type(options_type) :: options

!set default options 
call set_default_options(options)

!read command arguments 
call get_command_arguments(options)

!read options 
call read_options(options)

!display splash 
if (options%cdisplay) then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                flux GPU 2d                 |'
    write(*,'(A)')'|      2d unstructured euler flow solver     |'
    write(*,'(A)')'|       Version 2.0.0 || 15/11/2026          |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    if (options%mode == 'solve') then !primal solve
        write(*,'(A)') '                {primal solve}                 '
    elseif (options%mode == 'adjoint') then !adjoint solve
        write(*,'(A)') '                {adjoint solve}                 '
    end if 
    write(*,'(A)') ' '
end if

!import and preprocess mesh 
if (options%cdisplay) then
    write(*,'(A)') '--> reading mesh: '//trim(options%meshpath)//trim(options%meshname)
end if 
call read_mesh(mesh,options%meshpath)
if (options%cdisplay) then
    write(*,'(A,I0,A)') '    {ncell = ',mesh%ncell,'}' 
    write(*,'(A,I0,A)') '    {nedge = ',mesh%nedge,'}' 
    write(*,'(A,I0,A)') '    {nvertex = ',mesh%nvertex,'}' 
end if 
if (options%cdisplay) then
    write(*,'(A)') '--> preprocessing mesh'
end if 
! if (options%cdisplay) then
!     write(*,'(A,E12.6,A,E12.6,A)') '    {cell volume (max/min) = ',maxval(mesh%cells_volume),' / ',&
!     minval(mesh%cells_volume),'}' 
! end if 






end program flux2d