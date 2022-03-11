module mpi_util

    use mpi
    use config, only : nmember

    implicit none

    integer, protected          :: myid, nproc
    logical, protected          :: is_root 
    integer, parameter, private :: nxb = 1   !block size along x  
    integer, parameter, private :: nyb = 1   !block size along y  
    integer,            private :: nproc_x, nproc_y
    real*8,             private :: begin

    type local_domain
        integer                            :: loc_nx, loc_ny
        integer, dimension(:), allocatable :: xloc, yloc

        integer                            :: loc_nx_u, loc_ny_v
        integer, dimension(:), allocatable :: xloc_u, yloc_v
    end type local_domain

    type(local_domain), dimension(:), allocatable, target :: cpu

    !====================================================

    type request_list
        integer                     :: req
        integer                     :: idx
        type(request_list), pointer :: next => null()
    end type request_list

    type(request_list), pointer, private     :: req_list => null()
    type(request_list), pointer              :: req_ptr  => null()

contains

    subroutine letkf_init
        implicit none
        integer               :: mpierr
        integer, dimension(2) :: dims

        call mpi_init(mpierr)
        call mpi_comm_rank(mpi_comm_world, myid,  mpierr)
        call mpi_comm_size(mpi_comm_world, nproc, mpierr)
        begin   = mpi_wtime()
        is_root = myid == 0

        dims = 0
        call mpi_dims_create(nproc, 2, dims, mpierr)

        allocate(cpu(0:nproc-1))

        nproc_x = dims(1)
        nproc_y = dims(2)
    end subroutine letkf_init

    subroutine letkf_finalize
        implicit none
        integer :: mpierr

        deallocate(cpu)
        call mpi_finalize(mpierr)
    end subroutine letkf_finalize

    function timer() result(elapsed_time)
        implicit none
        real :: elapsed_time

        elapsed_time = mpi_wtime() - begin
    end function timer

    subroutine letkf_local_info(nx, ny)
        implicit none
        integer, intent(in) ::  nx, ny
        integer             ::  id, id_x, id_y, start, stride
        integer             ::  np_x, np_y
        integer             ::   i, j 

        do id_y = 0, nproc_y-1
        do id_x = 0, nproc_x-1
            id = id_x + id_y * nproc_x

            !split along x
            np_x   = 0
            start  = id_x * nxb+1
            stride = nproc_x * nxb

            do i = start, nx, stride
            do j = i, min(nx, i+nxb-1)
                np_x = np_x + 1
            end do
            end do
            cpu(id)%loc_nx = np_x

            allocate(cpu(id)%xloc(np_x))

            np_x = 0   !as counter
            do i = start, nx, stride
            do j = i, min(nx, i+nxb-1)
                np_x = np_x + 1
                cpu(id)%xloc(np_x) = j
            end do
            end do


            !split along y
            np_y   = 0
            start  = id_y * nyb+1
            stride = nproc_y * nyb

            do i = start, ny, stride
            do j = i, min(ny, i+nyb-1)
                np_y = np_y + 1
            end do
            end do
            cpu(id)%loc_ny = np_y

            allocate(cpu(id)%yloc(np_y))

            np_y = 0   !as counter
            do i = start, ny, stride
            do j = i, min(ny, i+nyb-1)
                np_y = np_y + 1
                cpu(id)%yloc(np_y) = j
            end do
            end do

            !split along x for u
            np_x   = 0
            start  = id_x * nxb+1
            stride = nproc_x * nxb

            do i = start, nx+1, stride
            do j = i, min(nx+1, i+nxb-1)
                np_x = np_x + 1
            end do
            end do
            cpu(id)%loc_nx_u = np_x

            allocate(cpu(id)%xloc_u(np_x))

            np_x = 0   !as counter
            do i = start, nx+1, stride
            do j = i, min(nx+1, i+nxb-1)
                np_x = np_x + 1
                cpu(id)%xloc_u(np_x) = j
            end do
            end do


            !split along y for v
            np_y   = 0
            start  = id_y * nyb+1
            stride = nproc_y * nyb

            do i = start, ny+1, stride
            do j = i, min(ny+1, i+nyb-1)
                np_y = np_y + 1
            end do
            end do
            cpu(id)%loc_ny_v = np_y

            allocate(cpu(id)%yloc_v(np_y))

            np_y = 0   !as counter
            do i = start, ny+1, stride
            do j = i, min(ny+1, i+nyb-1)
                np_y = np_y + 1
                cpu(id)%yloc_v(np_y) = j
            end do
            end do
        end do
        end do


        if(is_root) then
            open(99, file='rsl.out.0000', form='formatted')
            do i = 1, nproc
                write(99, *) '==========================================='
                write(99, 2001) i, cpu(i-1)%loc_nx, cpu(i-1)%loc_ny
                write(99, *) cpu(i-1)%xloc
                write(99, *) cpu(i-1)%yloc
            end do
            close(99)
        end if
 2001   format("process(",i3,")",2x,"dims:",1x,"(",i3,",",i3,")")
    end subroutine letkf_local_info

    subroutine letkf_scatter_grid(global, local, stagger)

        implicit none

        real, dimension(:,:,:),   intent(in)  :: global
        real, dimension(:,:,:,:), intent(out) :: local
        integer,                  intent(in)  :: stagger

        !local
        integer                               :: i, nx, ny, nz
        integer                               :: loc_nx, loc_ny
        integer                               :: st, ed
        integer, dimension(4)                 :: dims
        integer, dimension(:), pointer        :: idxx => null()
        integer, dimension(:), pointer        :: idxy => null()
        real,    dimension(:), allocatable    :: tmp

        !local mpi
        integer                               :: mpierr
        integer, dimension(0:nproc-1)         :: sendcnts, sdispls
        integer, dimension(0:nproc-1)         :: recvcnts, rdispls

        dims = shape(local)
        nx   = dims(1)
        ny   = dims(2)
        nz   = dims(3)

        sdispls  = 0
        rdispls  = 0
        sendcnts = 0
        recvcnts = 0

        rdispls (0:nmember-1) = [ (i, i = 0, nmember-1) ] * nx * ny * nz 
        recvcnts(0:nmember-1) = nx * ny * nz
        !==============================================================

        if(myid < nmember) then
            allocate(tmp(size(global)))

            do i = 0, nproc-1
                loc_nx = cpu(i) % loc_nx
                loc_ny = cpu(i) % loc_ny
                if(stagger == 1) loc_nx = cpu(i) % loc_nx_u
                if(stagger == 2) loc_ny = cpu(i) % loc_ny_v

                if(stagger == 1) then
                    idxx => cpu(i) % xloc_u
                else
                    idxx => cpu(i) % xloc
                end if

                if(stagger == 2) then
                    idxy => cpu(i) % yloc_v
                else
                    idxy => cpu(i) % yloc
                end if

                if(i == 0) then
                    sdispls(i) = 0
                else
                    sdispls(i) = sdispls(i-1) + sendcnts(i-1)
                end if
                sendcnts(i) = loc_nx * loc_ny * nz

                st = sdispls(i)+1
                ed = sdispls(i)+sendcnts(i)
                tmp(st:ed) = reshape(global(idxx,idxy,:), [sendcnts(i)])

                nullify(idxx, idxy)
            end do
        end if

        call mpi_alltoallv(  tmp, sendcnts, sdispls, mpi_real, &
                           local, recvcnts, rdispls, mpi_real, &
                           mpi_comm_world, mpierr)
        if(allocated(tmp)) deallocate(tmp)

    end subroutine letkf_scatter_grid

    subroutine letkf_gather_grid(local, global, stagger)

        implicit none

        real, dimension(:,:,:,:), intent(in)    :: local
        real, dimension(:,:,:),   intent(inout) :: global
        integer,                  intent(in)    :: stagger

        !local
        integer                               :: i, nx, ny, nz
        integer                               :: loc_nx, loc_ny
        integer                               :: st, ed
        integer, dimension(4)                 :: dims
        integer, dimension(:), pointer        :: idxx => null()
        integer, dimension(:), pointer        :: idxy => null()
        real,    dimension(:), allocatable    :: tmp

        !local mpi
        integer                               :: mpierr
        integer, dimension(0:nproc-1)         :: sendcnts, sdispls
        integer, dimension(0:nproc-1)         :: recvcnts, rdispls


        !determine shape of local array
        dims = shape(local)
        nx   = dims(1)
        ny   = dims(2)
        nz   = dims(3)

        sdispls  = 0
        rdispls  = 0
        sendcnts = 0
        recvcnts = 0

        sdispls (0:nmember-1) = [ (i, i = 0, nmember-1) ] * nx * ny * nz 
        sendcnts(0:nmember-1) = nx * ny * nz
        !==============================================================

        if(myid < nmember) then
            allocate(tmp(size(global)))

            do i = 0, nproc-1
                loc_nx = cpu(i) % loc_nx
                loc_ny = cpu(i) % loc_ny
                if(stagger == 1) loc_nx = cpu(i) % loc_nx_u
                if(stagger == 2) loc_ny = cpu(i) % loc_ny_v

                if(i == 0) then
                    rdispls(i) = 0
                else
                    rdispls(i) = rdispls(i-1) + recvcnts(i-1)
                end if
                recvcnts(i) = loc_nx * loc_ny * nz
            end do
        end if

        call mpi_alltoallv(local, sendcnts, sdispls, mpi_real, &
                             tmp, recvcnts, rdispls, mpi_real, &
                           mpi_comm_world, mpierr)

        if(allocated(tmp)) then
            do i = 0, nproc-1
                loc_nx = cpu(i) % loc_nx
                loc_ny = cpu(i) % loc_ny
                if(stagger == 1) loc_nx = cpu(i) % loc_nx_u
                if(stagger == 2) loc_ny = cpu(i) % loc_ny_v

                if(stagger == 1) then
                    idxx => cpu(i) % xloc_u
                else
                    idxx => cpu(i) % xloc
                end if

                if(stagger == 2) then
                    idxy => cpu(i) % yloc_v
                else
                    idxy => cpu(i) % yloc
                end if

                st = rdispls(i)+1
                ed = rdispls(i)+recvcnts(i)
                global(idxx,idxy,:) = reshape(tmp(st:ed), [loc_nx, loc_ny, nz])

                nullify(idxx, idxy)
            end do

            deallocate(tmp)
        end if
 
    end subroutine letkf_gather_grid

    subroutine letkf_scatter_hcoord(global_lat, local_lat, &
                                    global_lon, local_lon, stagger)

        implicit none

        real, dimension(:,:), intent(in)      :: global_lat, global_lon
        real, dimension(:,:), intent(out)     ::  local_lat,  local_lon
        integer,              intent(in)      :: stagger  !0: no   stagger
                                                          !1: u    stagger
                                                          !2: v    stagger
        !local
        integer                               :: i, nx, ny
        integer                               :: loc_nx, loc_ny
        integer                               :: st, ed
        integer, dimension(2)                 :: dims
        integer, dimension(:), pointer        :: idxx => null()
        integer, dimension(:), pointer        :: idxy => null()
        real,    dimension(:), allocatable    :: tmp_lat, tmp_lon

        !local mpi
        integer                               :: mpierr
        integer, dimension(0:nproc-1)         :: sendcnts, sdispls
        integer                               :: recvcnts


        !determine shape of local array
        dims = shape(local_lat)
        nx   = dims(1)
        ny   = dims(2)

        sdispls  = 0
        sendcnts = 0
        recvcnts = nx * ny

        !==============================================================

        if(is_root) then
            allocate(tmp_lat(size(global_lat)), &
                     tmp_lon(size(global_lon)))

            do i = 0, nproc-1
                loc_nx = cpu(i) % loc_nx
                loc_ny = cpu(i) % loc_ny
                if(stagger == 1) loc_nx = cpu(i) % loc_nx_u
                if(stagger == 2) loc_ny = cpu(i) % loc_ny_v

                if(stagger == 1) then
                    idxx => cpu(i) % xloc_u
                else
                    idxx => cpu(i) % xloc
                end if

                if(stagger == 2) then
                    idxy => cpu(i) % yloc_v
                else
                    idxy => cpu(i) % yloc
                end if

                if(i == 0) then
                    sdispls(i) = 0
                else
                    sdispls(i) = sdispls(i-1) + sendcnts(i-1)
                end if
                sendcnts(i) = loc_nx * loc_ny

                st = sdispls(i)+1
                ed = sdispls(i)+sendcnts(i)
                tmp_lat(st:ed) = reshape(global_lat(idxx,idxy), [sendcnts(i)])
                tmp_lon(st:ed) = reshape(global_lon(idxx,idxy), [sendcnts(i)])

                nullify(idxx, idxy)
            end do
        end if

        call mpi_scatterv(   tmp_lat, sendcnts, sdispls, mpi_real, &
                           local_lat, recvcnts,          mpi_real, &
                           0, mpi_comm_world, mpierr)
        call mpi_scatterv(   tmp_lon, sendcnts, sdispls, mpi_real, &
                           local_lon, recvcnts,          mpi_real, &
                           0, mpi_comm_world, mpierr)
        if(is_root) deallocate(tmp_lat, tmp_lon)

    end subroutine letkf_scatter_hcoord


    subroutine letkf_scatter_vcoord(global, local, stagger)

        use param, only : g
        implicit none

        real, dimension(:,:,:), intent(in)       :: global
        real, dimension(:,:,:), intent(out)      :: local
        integer,                intent(in)       :: stagger  !0: no   stagger
                                                             !1: w/ph stagger
        !local
        integer                                  :: i, nx, ny, nz, nz_ph
        integer                                  :: loc_nx, loc_ny
        integer, dimension(3)                    :: dims
        integer                                  :: st, ed
        integer, dimension(:),       pointer     :: idxx => null()
        integer, dimension(:),       pointer     :: idxy => null()
        real,    dimension(:),       allocatable :: tmp
        real,    dimension(:,:,:),   allocatable :: tmp3d
        real,    dimension(:,:,:,:), allocatable :: tmp4d
        real,    dimension(nmember)              :: x

        !local mpi
        integer                               :: mpierr
        integer, dimension(0:nproc-1)         :: sendcnts, sdispls
        integer, dimension(0:nproc-1)         :: recvcnts, rdispls


        !determine shape of local array
        dims = shape(local)
        nx = dims(1)
        ny = dims(2)
        nz = dims(3)

        select case (stagger)
        case (-1, 1)
            nz_ph = nz
        case default
            nz_ph = nz + 1
        end select

        if(stagger /= -1) then
            allocate(tmp3d(nx, ny, nz_ph), &
                     tmp4d(nx, ny, nz_ph, 0:nmember-1))

            sdispls  = 0
            rdispls  = 0
            sendcnts = 0
            recvcnts = 0

            rdispls (0:nmember-1) = [ (i, i = 0, nmember-1) ] * nx * ny * nz_ph 
            recvcnts(0:nmember-1) = nx * ny * nz_ph
            !==============================================================

            if(myid < nmember) then
                allocate(tmp(size(global)))

                do i = 0, nproc-1
                    loc_nx = cpu(i) % loc_nx
                    loc_ny = cpu(i) % loc_ny

                    idxx => cpu(i) % xloc
                    idxy => cpu(i) % yloc

                    if(i == 0) then
                        sdispls(i) = 0
                    else
                        sdispls(i) = sdispls(i-1) + sendcnts(i-1)
                    end if
                    sendcnts(i) = loc_nx * loc_ny * nz_ph

                    st = sdispls(i)+1
                    ed = sdispls(i)+sendcnts(i)
                    tmp(st:ed) = reshape(global(idxx,idxy,:), [sendcnts(i)])

                    nullify(idxx, idxy)
                end do
            end if

            call mpi_alltoallv(  tmp, sendcnts, sdispls, mpi_real, &
                               tmp4d, recvcnts, rdispls, mpi_real, &
                               mpi_comm_world, mpierr)
            if(allocated(tmp)) deallocate(tmp)

            !do ensemble mean of tmp4d and convert unit from m^2/s^2 to m
            x = 1.0
            call sgemv('n', nx*ny*nz_ph, nmember, 1.0/(g*nmember), tmp4d, nx*ny*nz_ph, x, 1, 0.0, tmp3d, 1)

            deallocate(tmp4d)

            select case (stagger)
            case(1)
                local = tmp3d
            case default
                local = (tmp3d(:,:,2:nz_ph)+tmp3d(:,:,1:nz)) * 0.5
            end select

            deallocate(tmp3d)
        else
            sdispls  = 0
            sendcnts = 0

            !==============================================================

            if(is_root) then
                allocate(tmp(size(global)))

                do i = 0, nproc-1
                    loc_nx = cpu(i) % loc_nx
                    loc_ny = cpu(i) % loc_ny

                    idxx => cpu(i) % xloc
                    idxy => cpu(i) % yloc

                    if(i == 0) then
                        sdispls(i) = 0
                    else
                        sdispls(i) = sdispls(i-1) + sendcnts(i-1)
                    end if
                    sendcnts(i) = loc_nx * loc_ny

                    st = sdispls(i)+1
                    ed = sdispls(i)+sendcnts(i)
                    tmp(st:ed) = reshape(global(idxx,idxy,1), [sendcnts(i)])

                    nullify(idxx, idxy)
                end do
            end if

            call mpi_scatterv(   tmp, sendcnts, sdispls, mpi_real, &
                               local,    nx*ny,          mpi_real, &
                               0, mpi_comm_world, mpierr)

            if(is_root) deallocate(tmp)
        end if

    end subroutine letkf_scatter_vcoord

    subroutine request_append
        implicit none
        integer                             :: n

        if(associated(req_list)) then
            n = req_ptr%idx
            allocate(req_ptr%next)
            req_ptr => req_ptr%next
            req_ptr%idx = n+1
        else
            allocate(req_list)
            req_ptr => req_list
            req_ptr%idx = 1
        end if
    end subroutine request_append

    subroutine wait_jobs(stlt, edlt)
        implicit none
        integer, value, intent(in)         :: stlt, edlt
        integer                            :: n
        integer, dimension(:), allocatable :: req
        integer                            :: mpierr
        type(request_list), pointer        :: next

        n = edlt-stlt+1
        allocate(req(n))

        n = 0
        req_ptr => req_list
        do
            if(associated(req_ptr)) then
                if(req_ptr%idx >= stlt .and. req_ptr%idx <= edlt) then
                    n       =  n+1
                    req(n)  =  req_ptr%req

                    next => req_ptr%next
                    nullify(req_ptr%next)
                    deallocate(req_ptr)
                else if(req_ptr%idx == stlt-1) then
                    next => req_ptr%next
                    nullify(req_ptr%next)
                else
                    next => req_ptr%next
                end if

                req_ptr => next
            else
                exit
            end if
        end do

        call mpi_waitall(n, req(1:n), mpi_statuses_ignore, mpierr)

        deallocate(req)
        if(stlt == 1) nullify(req_list)

    end subroutine wait_jobs

end module mpi_util
