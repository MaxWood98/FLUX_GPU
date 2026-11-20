!flux 2d geometry module 
!max wood
!version : 1.0.0
!updated : 19-11-26

!module 
module flux_geometry
use flux_data
contains

!edge volume ===============
function edge_volume(v1,v2) result(aseg) !+ve area for CCW oriented shapes
implicit none 

!variables - import
real(dp) :: v1(2),v2(2)

!result
real(dp) :: aseg

!segment area
aseg = 0.5d0*(v1(1)*v2(2) - v2(1)*v1(2))
return 
end function edge_volume

end module flux_geometry