module localization
    use kdtree2_module
    use gts_omboma
    use simulated_radar
    use config
    use param
    implicit none

    private
    public   :: lz_structure, build_tree, destroy_tree, get_lz, Gaspari_Cohn_1999

    type kdtree_type
        integer                :: mytype
        type(kdtree2), pointer :: tree
    end type kdtree_type

    type lz_structure
        integer                             :: mytype
        integer, dimension(:), allocatable  :: idx
        real,    dimension(:), allocatable  :: r2
    end type lz_structure

    type list_type
        integer                       :: obs_type
        real                          :: hclr_inv
        real                          :: vclr_inv
        type(list_type), pointer      :: next => null()
    end type list_type

    type(kdtree_type), dimension(:), allocatable, target ::   gts_tree
    type(kdtree_type), dimension(:), allocatable, target :: radar_tree

contains

    function build_tree(obs, ivar) result(succeed)
        implicit none
        class(obs_structure), dimension(:), intent(in)  :: obs
        logical                                         :: succeed

        !local
        type(kdtree_type), dimension(:), pointer :: tree    => null()
        type(list_type),                 pointer :: head    => null()
        type(list_type),                 pointer :: current => null()
        type(gts_config),                pointer :: gts_nml => null()
        type(radar_variable_config),     pointer :: rad_nml => null()
        integer                                  :: i, ivar, obs_type, ntype
        integer                                  :: dim
        real                                     :: hclr_inv, vclr_inv
        real, dimension(:,:), allocatable        :: xyz

        ntype    = 0
        succeed  = .false.

        select type (obs)
        type is (gts_structure)

            do obs_type = 1, num_gts_indexes     !loop all gts data type
                if(obs(obs_type)%nobs > 0) then
                    select case (obs_type)
                    case (synop)
                        gts_nml => synop_nml
                    case (metar)
                        gts_nml => metar_nml
                    case (ships)
                        gts_nml => ships_nml
                    case (sound)
                        gts_nml => sound_nml
                    case (gpspw)
                        gts_nml => gpspw_nml
                    case default
                        cycle
                    end select

                    if(gts_nml % use_it .and. gts_nml%hclr(ivar) > 0.) then
                        ntype    = ntype + 1
                        hclr_inv = 1.0 / (gts_nml%hclr(ivar)  * 1e3)

                        if(gts_nml%vclr(ivar) > 0.) then
                            vclr_inv = 1.0 / (gts_nml%vclr(ivar)  * 1e3)
                        else
                            vclr_inv = -1. !no vertical localization
                        end if

                        call append(head, current, obs_type, hclr_inv, vclr_inv)
                    end if

                    if(associated(gts_nml)) nullify(gts_nml)
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
                        rad_nml => radar_nml%dbz
                    case (vr)
                        rad_nml => radar_nml%vr
                    case (zdr)
                        rad_nml => radar_nml%zdr
                    case (kdp)
                        rad_nml => radar_nml%kdp
                    case default
                        cycle
                    end select

                    if(rad_nml % use_it .and. rad_nml%hclr(ivar) > 0.) then
                        ntype    = ntype + 1
                        hclr_inv = 1.0 / (rad_nml%hclr(ivar)  * 1e3)

                        if(rad_nml%vclr(ivar) > 0.) then
                            vclr_inv = 1.0 / (rad_nml%vclr(ivar)  * 1e3)
                        else
                            vclr_inv = -1. !no vertical localization
                        end if

                        call append(head, current, obs_type, hclr_inv, vclr_inv)
                    end if

                    if(associated(rad_nml)) nullify(rad_nml)
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
            obs_type       =  current % obs_type
            tree(i)%mytype =  obs_type

            !normalize by localization length scale
            allocate(xyz, source=obs(obs_type)%xyz)
            xyz(1:2,:) = xyz(1:2,:) * current % hclr_inv

            if(vclr_inv > 0.) then
                dim      = 3
                xyz(3,:) = xyz(3,:) * current % vclr_inv
            else
                dim      = 2
                xyz(3,:) = -1.
            end if

            !build kdtree for this variable
            tree(i)%tree => kdtree2_create(xyz, dim=dim)
            deallocate(xyz)
        end do

        nullify(tree, current)
        call destroy(head)

    end function build_tree

    subroutine destroy_tree
        implicit none
        integer :: i

        if(allocated(radar_tree)) then
            do i = 1, size(radar_tree)
                call kdtree2_destroy(radar_tree(i)%tree)
            end do
            deallocate(radar_tree)
        end if

        if(allocated(gts_tree)) then
            do i = 1, size(gts_tree)
                call kdtree2_destroy(gts_tree(i)%tree)
            end do
            deallocate(  gts_tree)
        end if
    end subroutine destroy_tree

    function get_lz(name, obs_lz, ivar, xyz) result(fail)
        implicit none
        character(len=*),                              intent(in)     :: name
        type(lz_structure), dimension(:), allocatable, intent(in out) :: obs_lz
        real,               dimension(3),              intent(in)     :: xyz
        integer,                                       intent(in)     :: ivar
        logical                                                       :: fail

        !local
        type(kdtree2_result), dimension(:), allocatable :: results
        type(gts_config),                   pointer     :: gts_nml => null()
        type(radar_variable_config),        pointer     :: rad_nml => null()

        !normalized searching radius
        real, parameter                                 :: r2 = gc1999 * gc1999

        integer                                         :: i, nlz
        real                                            :: hclr_inv, vclr_inv
        real,                 dimension(3)              :: tmp

        fail = .true.

        select case (name)
        case ('gts')
            if(.not. allocated(gts_tree)) return

            allocate(obs_lz(size(gts_tree)))

            do i = 1, size(gts_tree)
                obs_lz(i) % mytype = gts_tree(i) % mytype

                select case(gts_tree(i) % mytype)
                case (synop)
                    gts_nml => synop_nml
                case (metar)
                    gts_nml => metar_nml
                case (ships)
                    gts_nml => ships_nml
                case (sound)
                    gts_nml => sound_nml
                case (gpspw)
                    gts_nml => gpspw_nml
                end select

                allocate(results(gts_nml%max_lz_pts))

                hclr_inv = 1.0 / (gts_nml%hclr(ivar)  * 1e3)

                if(gts_nml%vclr(ivar) > 0.) then
                    vclr_inv = 1.0 / (gts_nml%vclr(ivar)  * 1e3)
                else
                    vclr_inv = -1.0
                end if

                !normalize by localization length scale
                tmp(1:2) = xyz(1:2) * hclr_inv

                if(vclr_inv > 0.) then  !3D localization
                    tmp(3) = xyz(3) * vclr_inv

                    !NOTE: the radius of kdtree2 is r2 not r
                    call kdtree2_r_nearest(gts_tree(i) % tree, tmp(1:3), r2, nlz, gts_nml%max_lz_pts, results)
                else !2D localization
                    !NOTE: the radius of kdtree2 is r2 not r
                    call kdtree2_r_nearest(gts_tree(i) % tree, tmp(1:2), r2, nlz, gts_nml%max_lz_pts, results)
                end if

                nullify(gts_nml)

                if(nlz > 0) then
                    fail = .false.

                    allocate(obs_lz(i) % idx(nlz), &
                             obs_lz(i) % r2(nlz))

                    obs_lz(i) % idx = results(1:nlz) % idx
                    obs_lz(i) % r2  = results(1:nlz) % dis
                end if

                deallocate(results)
            end do
        case ('radar')
            if(.not. allocated(radar_tree)) return

            allocate(obs_lz(size(radar_tree)))

            do i = 1, size(radar_tree)
                obs_lz(i) % mytype = radar_tree(i) % mytype

                select case (radar_tree(i) % mytype)
                case (dbz)
                    rad_nml => radar_nml%dbz
                case (vr)
                    rad_nml => radar_nml%vr
                case (zdr)
                    rad_nml => radar_nml%zdr
                case (kdp)
                    rad_nml => radar_nml%kdp
                end select

                allocate(results(rad_nml%max_lz_pts))

                hclr_inv = 1.0 / (rad_nml%hclr(ivar)  * 1e3)

                if(rad_nml%vclr(ivar) > 0.) then
                    vclr_inv = 1.0 / (rad_nml%vclr(ivar)  * 1e3)
                else
                    vclr_inv = -1.0
                end if

                !normalize by localization length scale
                tmp(1:2) = xyz(1:2) * hclr_inv

                if(vclr_inv > 0.) then  !3D localization
                    tmp(3) = xyz(3) * vclr_inv

                    !NOTE: the radius of kdtree2 is r2 not r
                    call kdtree2_r_nearest(radar_tree(i) % tree, tmp(1:3), r2, nlz, rad_nml%max_lz_pts, results)
                else !2D localization
                    !NOTE: the radius of kdtree2 is r2 not r
                    call kdtree2_r_nearest(radar_tree(i) % tree, tmp(1:2), r2, nlz, rad_nml%max_lz_pts, results)
                end if

                nullify(rad_nml)

                if(nlz > 0) then
                    fail = .false.

                    allocate(obs_lz(i) % idx(nlz), &
                             obs_lz(i) % r2(nlz))

                    obs_lz(i) % idx = results(1:nlz) % idx
                    obs_lz(i) % r2  = results(1:nlz) % dis  !return is also r2
                end if

                deallocate(results)
            end do
        case default
            stop "Error: Unknown name!"
        end select

        if(fail) deallocate(obs_lz)
 
    end function get_lz

    pure function Gaspari_Cohn_1999(x) result(func)
        implicit none
        real, intent(in) :: x
        real             :: func

        !local
        real, parameter  :: a  = sqrt(10. / 3.)
        real, parameter  :: a1 = -0.25 
        real, parameter  :: a2 = 0.5
        real, parameter  :: a3 = 0.625
        real, parameter  :: a4 = -5./3.
        real, parameter  :: a5 = 1.
        real, parameter  :: b1 = 1./12.
        real, parameter  :: b2 = -0.5
        real, parameter  :: b3 = 0.625
        real, parameter  :: b4 = 5./3.
        real, parameter  :: b5 = -5.
        real, parameter  :: b6 = 4.
        real, parameter  :: b7 = -2./3.
        real             :: z

        !x is already normalized by length scale
        z = x / a

        if( z <= 1. ) then
            func = z * z * ( z * ( z * ( a1 * z + a2 ) + a3 ) + a4 ) + a5
        else if( z <= 2. ) then
            func = z * ( z * ( z * ( z * ( b1 * z + b2 ) + b3 ) + b4 ) + b5 ) + b6 + b7 / z
        else
            func = 0.
        end if
    end function Gaspari_Cohn_1999


    !=============================================================

    subroutine append(head, current, obs_type, hclr_inv, vclr_inv)
        implicit none
        type(list_type), pointer, intent(in out) :: head, current
        integer,                  intent(in)     :: obs_type
        real,           optional, intent(in)     :: hclr_inv, vclr_inv


        if(associated(head)) then
            allocate(current % next)
            current => current % next
        else
            allocate(head)
            current => head
        end if

        current % obs_type = obs_type
        if(present(hclr_inv)) current % hclr_inv = hclr_inv
        if(present(vclr_inv)) current % vclr_inv = vclr_inv
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
