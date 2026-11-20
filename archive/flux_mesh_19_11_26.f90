!flux 2d mesh module 
!max wood
!version : 2.0.0
!updated : 15-11-26

!module 
module flux_mesh

!imports    
use flux_data
use flux_geometry

!local types
type cell_edges_accum
    integer(in32) :: nedge
    integer(in32), allocatable, dimension(:) :: edges
end type cell_edges_accum

!subrotuines 
contains

!initialise mesh ===============
subroutine initialise_mesh(mesh)
implicit none

!variables - inout 
type(mesh_type) :: mesh

!variables - local
integer(in32) :: ii,jj 
integer(in32) :: ceentry
! type(cell_edges_accum) :: cell_edges(mesh%ncell)

!sort edges 
call sort_mesh_edges(mesh)

! !count and accumulate edges per cell
! do ii = 1,mesh%ncell
!     cell_edges(ii)%nedge = 0
! end do 
! do ii = 1,mesh%nedge
!     cell_edges(mesh%edges_c1(ii))%nedge = cell_edges(mesh%edges_c1(ii))%nedge + 1
!     if (mesh%edges_c2(ii) > 0) then
!         cell_edges(mesh%edges_c2(ii))%nedge = cell_edges(mesh%edges_c2(ii))%nedge + 1
!     end if 
! end do
! mesh%cell_max_nedge = maxval(cell_edges(:)%nedge)
! mesh%cell_min_nedge = minval(cell_edges(:)%nedge)
! do ii = 1,mesh%ncell
!     allocate(cell_edges(ii)%edges(cell_edges(ii)%nedge))
!     cell_edges(ii)%nedge = 0
! end do
! do ii = 1,mesh%nedge
!     cell_edges(mesh%edges_c1(ii))%nedge = cell_edges(mesh%edges_c1(ii))%nedge + 1
!     cell_edges(mesh%edges_c1(ii))%edges(cell_edges(mesh%edges_c1(ii))%nedge) = ii
!     if (mesh%edges_c2(ii) > 0) then
!         cell_edges(mesh%edges_c2(ii))%nedge = cell_edges(mesh%edges_c2(ii))%nedge + 1
!         cell_edges(mesh%edges_c2(ii))%edges(cell_edges(mesh%edges_c2(ii))%nedge) = ii
!     end if 
! end do 

! !build cell edges csr arrays
! mesh%ncell_edges = 0
! do ii = 1,mesh%ncell
!     mesh%ncell_edges = mesh%ncell_edges + cell_edges(ii)%nedge
! end do
! allocate(mesh%cells_edges_start(mesh%ncell+1))
! allocate(mesh%cells_edges(mesh%ncell_edges))
! allocate(mesh%cells_edges_sign(mesh%ncell_edges))
! ceentry = 1
! do ii = 1,mesh%ncell
!     mesh%cells_edges_start(ii) = ceentry
!     do jj=1,cell_edges(ii)%nedge
!         mesh%cells_edges(ceentry) = cell_edges(ii)%edges(jj)
!         if (ii == mesh%edges_c1(cell_edges(ii)%edges(jj))) then 
!             mesh%cells_edges_sign(ceentry) = 1
!         else
!             mesh%cells_edges_sign(ceentry) = -1
!         end if 
!         ceentry = ceentry + 1
!     end do 
! end do
! mesh%cells_edges_start(mesh%ncell+1) = ceentry

!build cell-edge linking csr arrays
call link_cells2edges(mesh)

!debug ===
do ii = 1,4!mesh%ncell
    print *, '====== cell: ',ii
    do jj=mesh%cells_edges_start(ii),mesh%cells_edges_start(ii+1)-1
        ceentry = mesh%cells_edges(jj)
        print *, 'edge = ',ceentry,' : ',mesh%edges_c1(ceentry),mesh%edges_c2(ceentry),' sign = ',mesh%cells_edges_sign(jj)
    end do 
end do 
!debug ===

!evaluate cell volumes 
call get_cell_volumes(mesh)

!get edge geometries 
call get_edge_geometries(mesh)
end subroutine initialise_mesh

!build cell-edge csr ===============
subroutine link_cells2edges(mesh)
implicit none 

!variables - inout
type(mesh_type) :: mesh

!variables - local
integer(in32) :: ii 
integer(in32) :: cidx,ceidx
integer(in32) :: cell_nedge(mesh%ncell)

!count edges in each cell
cell_nedge(:) = 0
do ii=1,mesh%nedge
    cell_nedge(mesh%edges_c1(ii)) = cell_nedge(mesh%edges_c1(ii)) + 1
    if (mesh%edges_c2(ii) > 0) then
        cell_nedge(mesh%edges_c2(ii)) = cell_nedge(mesh%edges_c2(ii)) + 1
    end if 
end do
mesh%cell_max_nedge = maxval(cell_nedge)
mesh%cell_min_nedge = minval(cell_nedge)

!get cell start indices 
allocate(mesh%cells_edges_start(mesh%ncell+1))
mesh%cells_edges_start(:) = 0 
mesh%cells_edges_start(1) = 1
do ii=2,mesh%ncell+1 
    mesh%cells_edges_start(ii) = mesh%cells_edges_start(ii-1) + cell_nedge(ii-1)
end do 

!accumulate edges to cells and set the edge sign 
mesh%ncell_edges = mesh%cells_edges_start(mesh%ncell+1) - 1
allocate(mesh%cells_edges(mesh%ncell_edges))
allocate(mesh%cells_edges_sign(mesh%ncell_edges))

! print *, mesh%ncell_edges

! mesh%cells_edges(:) = 0 
! mesh%cells_edges_sign(:) = 0 

cell_nedge(:) = 0 
do ii=1,mesh%nedge

    cidx = mesh%edges_c1(ii)
    ceidx = mesh%cells_edges_start(cidx) + cell_nedge(cidx)
    mesh%cells_edges(ceidx) = ii
    mesh%cells_edges_sign(ceidx) = 1
    cell_nedge(cidx) = cell_nedge(cidx) + 1

    cidx = mesh%edges_c2(ii)
    if (cidx > 0) then
        ceidx = mesh%cells_edges_start(cidx) + cell_nedge(cidx)
        mesh%cells_edges(ceidx) = ii
        mesh%cells_edges_sign(ceidx) = -1
        cell_nedge(cidx) = cell_nedge(cidx) + 1
    end if 

end do 


end subroutine link_cells2edges

!get edge geometries ===============
subroutine get_edge_geometries(mesh)
implicit none

!variables - inout
type(mesh_type) :: mesh

!variables - local
integer(in32) :: ii

!allocate
allocate(mesh%edges_length(mesh%nedge))
allocate(mesh%edges_nx(mesh%nedge))
allocate(mesh%edges_ny(mesh%nedge))
allocate(mesh%edges_dx(mesh%nedge))
allocate(mesh%edges_dy(mesh%nedge))

!evaluate 
do ii=1,mesh%nedge
    mesh%edges_dx(ii) = mesh%vertices(mesh%edges_v2(ii),1) - mesh%vertices(mesh%edges_v1(ii),1)
    mesh%edges_dy(ii) = mesh%vertices(mesh%edges_v2(ii),2) - mesh%vertices(mesh%edges_v1(ii),2)
    mesh%edges_length(ii) = sqrt(mesh%edges_dx(ii)**2 + mesh%edges_dy(ii)**2)
    mesh%edges_nx(ii) = mesh%edges_dy(ii)/mesh%edges_length(ii)
    mesh%edges_ny(ii) = -mesh%edges_dx(ii)/mesh%edges_length(ii)
end do 
end subroutine get_edge_geometries

!get cell volumes ===============
subroutine get_cell_volumes(mesh)
implicit none

!variables - inout
type(mesh_type) :: mesh

!variables - local
integer(in32) :: ii
real(dp) :: evolume

!evaluate 
allocate(mesh%cells_volume(mesh%ncell))
mesh%cells_volume(:) = 0.0d0
do ii = 1,mesh%nedge
    evolume = edge_volume(mesh%vertices(mesh%edges_v1(ii),:),mesh%vertices(mesh%edges_v2(ii),:))
    mesh%cells_volume(mesh%edges_c1(ii)) = mesh%cells_volume(mesh%edges_c1(ii)) + evolume
    if (mesh%edges_c2(ii) > 0) then
        mesh%cells_volume(mesh%edges_c2(ii)) = mesh%cells_volume(mesh%edges_c2(ii)) - evolume
    end if
end do
end subroutine get_cell_volumes

!sort cells ===============


!sort edges ===============
subroutine sort_mesh_edges(mesh)
implicit none

!variables - inout
type(mesh_type) :: mesh

!variables - local
integer(in32) :: ii
integer(in32) :: edge_order(mesh%nedge)
integer(in32) :: edge_v1_tmp(mesh%nedge),edge_v2_tmp(mesh%nedge),edge_c1_tmp(mesh%nedge),edge_c2_tmp(mesh%nedge)

!sort edges
call sort_edges(mesh%edges_c1,mesh%edges_c2,edge_order)

!sort arrays
do ii=1,mesh%nedge
    edge_v1_tmp(ii) = mesh%edges_v1(edge_order(ii))
    edge_v2_tmp(ii) = mesh%edges_v2(edge_order(ii))
    edge_c1_tmp(ii) = mesh%edges_c1(edge_order(ii))
    edge_c2_tmp(ii) = mesh%edges_c2(edge_order(ii))
end do
mesh%edges_v1(:) = edge_v1_tmp(:)
mesh%edges_v2(:) = edge_v2_tmp(:)
mesh%edges_c1(:) = edge_c1_tmp(:)
mesh%edges_c2(:) = edge_c2_tmp(:)
end subroutine sort_mesh_edges

!sort edges ===============
subroutine sort_edges(cell1, cell2, order)
    implicit none 
    integer(in32), intent(in) :: cell1(:), cell2(:)
    integer(in32), intent(out) :: order(:)
    integer(in32), allocatable :: temp(:)
    integer(in32) :: n, i

    !sort
    n = size(cell1)
    allocate(temp(n))
    do i = 1, n
        order(i) = i
    end do
    call merge_sort(1,n)
    deallocate(temp)

contains

    ! Recursive merge sort
    recursive subroutine merge_sort(lo,hi)
        implicit none 
        integer(in32), intent(in) :: lo,hi
        integer(in32) :: mid
        if (lo >= hi) return
        mid = (lo+hi)/2
        call merge_sort(lo,mid)
        call merge_sort(mid+1,hi)
        call merge(lo,mid,hi)
    end subroutine merge_sort

    !Compare two edges
    logical function edge_less(e1,e2)
        implicit none
        integer(in32), intent(in) :: e1,e2
        integer(in32) :: a1,b1,a2,b2
        integer(in32) :: t

        !Exract cell indices 
        a1 = cell1(e1)
        b1 = cell2(e1)
        a2 = cell1(e2)
        b2 = cell2(e2)

        !Canonicalise edge 1
        if (a1 < 0) then
            t  = a1
            a1 = b1
            b1 = t
        end if
        if (a1 > 0 .and. b1 > 0) then
            if (a1 > b1) then
                t  = a1
                a1 = b1
                b1 = t
            end if
        end if

        !Canonicalise edge 2
        if (a2 < 0) then
            t  = a2
            a2 = b2
            b2 = t
        end if
        if (a2 > 0 .and. b2 > 0) then
            if (a2 > b2) then
                t  = a2
                a2 = b2
                b2 = t
            end if
        end if

        !Compare first cell
        if (a1 < a2) then
            edge_less = .true.
            return
        else if (a1 > a2) then
            edge_less = .false.
            return
        end if

        !Compare second cell / boundary identifier
        if (b1 < b2) then
            edge_less = .true.
        else if (b1 > b2) then
            edge_less = .false.
        else
            edge_less = e1 < e2
        end if
    end function edge_less

    !Merge two sorted sections
    subroutine merge(lo,mid,hi)
        implicit none 
        integer(in32), intent(in) :: lo,mid,hi
        integer(in32) :: i,j,k
        i = lo
        j = mid+1
        k = lo
        do while (i <= mid .and. j <= hi)
            if (edge_less(order(i),order(j))) then
                temp(k) = order(i)
                i = i+1
            else
                temp(k) = order(j)
                j = j+1
            end if
            k = k+1
        end do
        do while (i <= mid)
            temp(k) = order(i)
            i = i+1
            k = k+1
        end do
        do while (j <= hi)
            temp(k) = order(j)
            j = j+1
            k = k+1
        end do
        order(lo:hi) = temp(lo:hi)
    end subroutine merge
end subroutine sort_edges

end module flux_mesh