module localization
!   use kdtree2_module
    use gts_omboma
    use simulated_radar
    use config
    use param
    use kdtree
    implicit none

    private
    public   :: lz_structure, build_tree, destroy_tree, get_hlz, get_vlz, destroy_vlz

    type kdtree_type
        integer                :: mytype
!       type(kdtree2), pointer :: tree
        type(kd_root)          :: tree
    end type kdtree_type

    type vlz_structure
        integer, dimension(:), allocatable  :: idx
        real,    dimension(:), allocatable  :: hdistance
        integer, dimension(:), allocatable  :: idk
        real,    dimension(:), allocatable  :: vdistance
    end type vlz_structure

    type lz_structure
        integer                             :: mytype
        integer, dimension(:), allocatable  :: idx
        real,    dimension(:), allocatable  :: hdistance
        type(vlz_structure),   allocatable  :: vert
    end type lz_structure

    type list_type
        integer                       :: idx
        real                          :: hdistance
        integer                       :: idk
        real                          :: vdistance
        type(list_type), pointer      :: next => null()
    end type list_type

    type(kdtree_type), dimension(:), allocatable, target ::   gts_tree
    type(kdtree_type), dimension(:), allocatable, target :: radar_tree

contains

    function build_tree(obs) result(succeed)
        implicit none
        class(obs_structure), dimension(:), intent(in)  :: obs
        logical                                         :: succeed

        !local
        type(kdtree_type), dimension(:), pointer :: tree    => null()
        type(list_type),                 pointer :: head    => null()
        type(list_type),                 pointer :: current => null()
        integer                                  :: i, obs_type, ntype
        logical                                  :: use_this

        ntype    = 0
        succeed  = .false.

        select type (obs)
        type is (gts_structure)

            do obs_type = 1, num_gts_indexes
                if(obs(obs_type)%nobs > 0) then
                    select case (obs_type)
                    case (synop)
                        use_this = synop_nml%use_it
                    case (metar)
                        use_this = metar_nml%use_it
                    case (ships)
                        use_this = ships_nml%use_it
                    case (sound)
                        use_this = sound_nml%use_it
                    case (gpspw)
                        use_this = gpspw_nml%use_it
                    case default
                        use_this = .false.
                    end select

                    if(use_this) then
                        ntype = ntype + 1
                        call append(head, current, obs_type)
                    end if
                end if
            end do

            if(ntype == 0) return

            allocate(gts_tree(ntype))
            tree => gts_tree

        type is (radar_structure)

            do obs_type = 1, num_radar_indexes
                if(obs(obs_type)%nobs > 0) then
                    select case (obs_type)
                    case (dbz)
                        use_this = radar_nml%dbz%use_it
                    case (vr)
                        use_this = radar_nml%vr%use_it
                    case (zdr)
                        use_this = radar_nml%zdr%use_it
                    case (kdp)
                        use_this = radar_nml%kdp%use_it
                    case default
                        use_this = .false.
                    end select

                    if(use_this) then
                        ntype = ntype + 1
                        call append(head, current, obs_type)
                    end if
                end if
            end do

            if(ntype == 0) return

            allocate(radar_tree(ntype))
            tree => radar_tree

        class default
            stop "Unknown type!"
        end select

        succeed = .true.

        current => head
        do i = 1, ntype
            if(i > 1) current => current % next
            obs_type       =  current % idx
            tree(i)%mytype =  obs_type
            call kd_init(tree(i)%tree, obs(obs_type)%lon, obs(obs_type)%lat)
        end do

        nullify(tree, current)
        call destroy(head)

    end function build_tree

    subroutine destroy_tree
        implicit none
        integer :: i

        if(allocated(radar_tree)) then
            do i = 1, size(radar_tree)
                call kd_free(radar_tree(i)%tree)
            end do
            deallocate(radar_tree)
        end if

        if(allocated(gts_tree)) then
            do i = 1, size(gts_tree)
                call kd_free(gts_tree(i)%tree)
            end do
            deallocate(  gts_tree)
        end if
    end subroutine destroy_tree

    function get_hlz(name, obs_lz, ivar, lat, lon) result(fail)
        implicit none
        character(len=*),                              intent(in)     :: name
        type(lz_structure), dimension(:), allocatable, intent(in out) :: obs_lz
        integer,                                       intent(in)     :: ivar
        real,                                          intent(in)     :: lat, lon
        logical                                                       :: fail

        !local
        integer, dimension(max_lz_pts) :: idx
        real,    dimension(max_lz_pts) :: distance
        integer                        :: i, nlz
        real                           :: hclr

        fail = .true.

        select case (name)
        case ('gts')
            if(.not. allocated(gts_tree)) return

            allocate(obs_lz(size(gts_tree)))

            do i = 1, size(gts_tree)
                obs_lz(i) % mytype    = gts_tree(i) % mytype

                select case(gts_tree(i) % mytype)
                case (synop)
                    hclr =  synop_nml%hclr(ivar)  * 1e3
                case (metar)
                    hclr =  metar_nml%hclr(ivar)  * 1e3
                case (gpspw)
                    hclr =  gpspw_nml%hclr(ivar)  * 1e3
                case (ships)
                    hclr =  ships_nml%hclr(ivar)  * 1e3
                case (sound)
                    hclr =  sound_nml%hclr(ivar)  * 1e3
                case default
                    stop "Unknown obs type"
                end select

                call kd_search_radius(gts_tree(i) % tree, lon, lat, &
                                      hclr, idx, distance, nlz)

                if(nlz > 0) then
                    fail = .false.

                    allocate(obs_lz(i) %       idx(nlz), &
                             obs_lz(i) % hdistance(nlz))

                    obs_lz(i) % idx       = idx(1:nlz)
                    obs_lz(i) % hdistance = distance(1:nlz)
                end if
            end do
        case ('radar')
            if(.not. allocated(radar_tree)) return

            allocate(obs_lz(size(radar_tree)))

            do i = 1, size(radar_tree)
                obs_lz(i) % mytype    = radar_tree(i) % mytype

                select case(radar_tree(i) % mytype)
                case (dbz)
                    hclr = radar_nml%dbz%hclr(ivar) * 1e3
                case (vr)
                    hclr = radar_nml% vr%hclr(ivar) * 1e3
                case (zdr)
                    hclr = radar_nml%zdr%hclr(ivar) * 1e3
                case (kdp)
                    hclr = radar_nml%kdp%hclr(ivar) * 1e3
                case default
                    stop "Unknown obs type"
                end select

                call kd_search_radius(radar_tree(i) % tree, lon, lat, &
                                      hclr, idx, distance, nlz)

                if(nlz > 0) then
                    fail = .false.

                    allocate(obs_lz(i) %       idx(nlz), &
                             obs_lz(i) % hdistance(nlz))

                    obs_lz(i) % idx       = idx(1:nlz)
                    obs_lz(i) % hdistance = distance(1:nlz)
                end if
            end do
        case default
            stop "Error: Unknown name!"
        end select

        if(fail) deallocate(obs_lz)
 
    end function get_hlz

    function get_vlz(obs, obs_lz, ivar, alt) result(fail)
        implicit none
        class(obs_structure), dimension(:), intent(in)     :: obs
        type(  lz_structure), dimension(:), intent(in out) :: obs_lz
        integer,                            intent(in)     :: ivar
        real,                               intent(in)     :: alt
        logical                                            :: fail

        !local
        type(list_type), pointer :: head    => null()
        type(list_type), pointer :: current => null()
        integer                  :: i, j, k, idx, obs_type, nlz
        real                     :: vdistance, vclr

        fail = .true.

        select type (obs)
        type is (radar_structure)
            do i = 1, size(obs_lz)
                if(.not. allocated(obs_lz(i) % idx)) cycle

                nlz      = 0
                obs_type = obs_lz(i) % mytype

                select case (obs_type)
                case (dbz)
                    vclr = radar_nml%dbz%vclr(ivar) * 1e3
                case (vr)
                    vclr = radar_nml% vr%vclr(ivar) * 1e3
                case (zdr)
                    vclr = radar_nml%zdr%vclr(ivar) * 1e3
                case (kdp)
                    vclr = radar_nml%kdp%vclr(ivar) * 1e3
                end select

                do j = 1, size(obs_lz(i)%idx)
                    idx = obs_lz(i)%idx(j)
                    vdistance = abs(obs(obs_type)%alt(idx) - alt)
                    if(vdistance < vclr) then
                        nlz = nlz+1
                        call append(head, current, idx, obs_lz(i)%hdistance(j), 0, vdistance)
                    end if
                end do

                if(nlz > 0) then
                    fail = .false.

                    allocate(obs_lz(i)%vert)
                    allocate(obs_lz(i)%vert%idx      (nlz), &
                             obs_lz(i)%vert%hdistance(nlz), &
                             obs_lz(i)%vert%idk      (nlz), &
                             obs_lz(i)%vert%vdistance(nlz))

                    current => head
                    do j = 1, nlz
                        if(j > 1) current => current%next

                        obs_lz(i)%vert%idx(j)       = current%idx
                        obs_lz(i)%vert%hdistance(j) = current%hdistance
                        obs_lz(i)%vert%idk(j)       = current%idk
                        obs_lz(i)%vert%vdistance(j) = current%vdistance
                    end do


                    nullify(current)
                    call destroy(head)
                end if
            end do

        type is (gts_structure)
            do i = 1, size(obs_lz)
                if(.not. allocated(obs_lz(i) % idx)) cycle

                nlz      = 0
                obs_type = obs_lz(i) % mytype

                select case (obs_type)
                case (synop, metar, gpspw, ships)
                    select case (obs_type)
                    case (synop)
                        vclr =  synop_nml%vclr(ivar)  * 1e3
                    case (metar)
                        vclr =  metar_nml%vclr(ivar)  * 1e3
                    case (gpspw)
                        vclr =  gpspw_nml%vclr(ivar)  * 1e3
                    case (ships)
                        vclr =  ships_nml%vclr(ivar)  * 1e3
                    end select
     
                    do j = 1, size(obs_lz(i)%idx)
                        idx = obs_lz(i)%idx(j)
                        vdistance = abs(obs(obs_type)%surf%h(idx) - alt)
                        if(vdistance < vclr) then
                            nlz = nlz+1
                            call append(head, current, idx, obs_lz(i)%hdistance(j), 0, vdistance)
                        end if
                    end do
                case (sound)
                    vclr =  sound_nml%vclr(ivar) * 1E3
                    do j = 1, size(obs_lz(i)%idx)
                        idx  = obs_lz(i)%idx(j)
                        do k = 1, obs(obs_type)%vert(idx)%nlev
                            vdistance = abs(obs(obs_type)%vert(idx)%h(k) - alt)
                            if(vdistance < vclr) then
                                nlz = nlz+1
                                call append(head, current, idx, obs_lz(i)%hdistance(j), k, vdistance)
                            end if
                        end do
                    end do
                end select

                if(nlz > 0) then
                    fail = .false.

                    allocate(obs_lz(i)%vert)
                    allocate(obs_lz(i)%vert%idx      (nlz), &
                             obs_lz(i)%vert%hdistance(nlz), &
                             obs_lz(i)%vert%idk      (nlz), &
                             obs_lz(i)%vert%vdistance(nlz))

                    current => head
                    do j = 1, nlz
                        if(j > 1) current => current%next

                        obs_lz(i)%vert%idx(j)       = current%idx
                        obs_lz(i)%vert%hdistance(j) = current%hdistance
                        obs_lz(i)%vert%idk(j)       = current%idk
                        obs_lz(i)%vert%vdistance(j) = current%vdistance
                    end do


                    nullify(current)
                    call destroy(head)
                end if
            end do
        class default
            stop "Unknown type!"
        end select

    end function get_vlz

    subroutine destroy_vlz(obs_lz)
        implicit none
        type(lz_structure), dimension(:), allocatable, intent(in out) :: obs_lz

        !local
        integer :: i

        if(allocated(obs_lz)) then
            do i = 1, size(obs_lz)
                if(allocated(obs_lz(i)%vert)) deallocate(obs_lz(i)%vert)
            end do
        end if
    end subroutine destroy_vlz

    !=============================================================

    subroutine append(head, current, idx, hdistance, idk, vdistance)
        implicit none
        type(list_type), pointer, intent(in out) :: head, current
        integer,                  intent(in)     :: idx
        integer,        optional, intent(in)     :: idk
        real,           optional, intent(in)     :: hdistance, vdistance


        if(associated(head)) then
            allocate(current % next)
            current => current % next
        else
            allocate(head)
            current => head
        end if

        current % idx = idx
        if(present(hdistance)) current % hdistance = hdistance
        if(present(idk))       current % idk       = idk
        if(present(vdistance)) current % vdistance = vdistance
    end subroutine append

    subroutine destroy(head)
        implicit none
        type(list_type), pointer, intent(in out) :: head
        type(list_type), pointer                 :: current, next

        current => head
        do while(associated(current))
            next => current % next
            deallocate(current)
            current => next
        end do

        nullify(head)
    end subroutine destroy


end module localization
