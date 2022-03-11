module letkf_core

    use localization,    only : lz_structure, get_lz, build_tree, destroy_tree, &
                                Gaspari_Cohn_1999
    use grid,            only : grid_structure
    use gts_omboma,      only : wrfda_gts
    use simulated_radar, only : cwb_radar
    use projection,      only : proj_type
    use param
    use mpi_util
    use eigen
    use config

    implicit none

    private
    public  :: letkf_driver

contains

    subroutine letkf_driver(wrf, gts, rad, proj)
        implicit none

        type(grid_structure),         intent(inout)   :: wrf
        type(wrfda_gts),              intent(in   )   :: gts
        type(cwb_radar),              intent(in   )   :: rad
        type(proj_type),              intent(in   )   :: proj

        !localization
        type(lz_structure), dimension(:), allocatable ::   gts_lz
        type(lz_structure), dimension(:), allocatable :: radar_lz

        !local
        integer                                       :: ivar
        integer                                       :: i, j, k
        integer                                       :: nx, ny, nz
        integer                                       :: ids, ide, jds, jde
        integer                                       :: loc_nx, loc_ny
        integer                                       :: Hstag,  Vstag
        logical                                       :: Hreset, Vreset
        real                                          :: inflat
        real, dimension(3)                            :: xyz
        real, dimension(nmember)                      :: xb
        real, dimension(:),       allocatable         :: yo
        real, dimension(:,:),     allocatable         :: yb, lat, lon
        real, dimension(:,:,:),   allocatable         :: alt, hgt, mu
        real, dimension(:,:,:,:), allocatable         :: var
        logical, dimension(2)                         :: succeed, fail

        call wait_jobs(1, req_ptr%idx)
        nx = wrf % nx
        ny = wrf % ny

        call letkf_local_info(nx, ny)

        Hstag  = 0
        Vstag  = 0

        update : do ivar = 1, max_vars
            if(len_trim(var_update(ivar)) == 0) exit
            if(is_root) print '(f7.3,1x,a,1x,10a)', timer(), "sec ==========> update", var_update(ivar)

            succeed(1) = build_tree(gts % platform, ivar)
            succeed(2) = build_tree(rad % radarobs, ivar)

            if(all(.not. succeed)) cycle

            inflat = (nmember-1) / multi_infl(ivar)
            nz     = wrf % nz

            loc_nx = cpu(myid)%loc_nx
            loc_ny = cpu(myid)%loc_ny

            select case(trim(var_update(ivar)))
            case ('U')
                loc_nx = cpu(myid)%loc_nx_u
            case ('V')
                loc_ny = cpu(myid)%loc_ny_v
            case ('W', 'PH')
                nz = nz + 1
            case ('MU')
                nz = 1
            end select

            allocate(var( loc_nx, loc_ny, nz, 0:nmember-1 ))

            !=================================================================================


            select case(trim(var_update(ivar)))
            case ('U')
                Hreset = check_coordinate(Hstag, 1)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % u, var, Hstag)
            case ('V')
                Hreset = check_coordinate(Hstag, 2)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % v, var, Hstag)
            case ('W')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 1)
                call letkf_scatter_grid  (wrf % w, var, Hstag)
            case ('T')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % t, var, Hstag)
            case ('QVAPOR')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % qv, var, Hstag)
            case ('QRAIN')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % qr, var, Hstag)
            case ('QSNOW')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % qs, var, Hstag)
            case ('QGRAUP')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % qg, var, Hstag)
            case ('QHAIL')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % qh, var, Hstag)
            case ('QNRAIN')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % nqr, var, Hstag)
            case ('QNSNOW')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % nqs, var, Hstag)
            case ('QNGRAUPEL')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % nqg, var, Hstag)
            case ('QNHAIL')
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % nqh, var, Hstag)
            case ('P')   !full pressure
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 0)
                call letkf_scatter_grid  (wrf % p, var, Hstag)
            case ('MU') !full MU
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag,-1)
                if(allocated(wrf % mu)) then
                    allocate(mu(nx, ny, 1))
                    mu(:,:,1) = wrf % mu
                end if
                call letkf_scatter_grid  (mu, var, Hstag)
            case ('PH')   !full geopotential height
                Hreset = check_coordinate(Hstag, 0)
                Vreset = check_coordinate(Vstag, 1)
                call letkf_scatter_grid  (wrf % ph, var, Hstag)
            case default
                print *, var_update(ivar)
                stop "Need to code for unknown variable"
            end select


            if((.not. allocated(lat)) .or. &
               (.not. allocated(lon)) .or. Hreset) then

                if( allocated(lat) ) deallocate(lat)
                if( allocated(lon) ) deallocate(lon)

                allocate(lat( loc_nx, loc_ny ), &
                         lon( loc_nx, loc_ny ))

                select case (Hstag)
                case(0)
                    call letkf_scatter_hcoord(wrf % xlat,   lat, &
                                              wrf % xlon,   lon, Hstag)
                case(1)
                    call letkf_scatter_hcoord(wrf % xlat_u, lat, &
                                              wrf % xlon_u, lon, Hstag)
                case(2)
                    call letkf_scatter_hcoord(wrf % xlat_v, lat, &
                                              wrf % xlon_v, lon, Hstag)
                end select
            end if


            if(.not. allocated(alt) .or. Vreset) then

                if( allocated(alt) ) deallocate(alt)

                allocate(alt( cpu(myid)%loc_nx, &
                              cpu(myid)%loc_ny, nz))

                select case (Vstag)
                case(-1)
                    if(allocated(wrf % hgt)) then
                        allocate(hgt(nx, ny, 1))
                        hgt(:,:,1) = wrf % hgt
                    end if
                    call letkf_scatter_vcoord(      hgt, alt, Vstag)
                    if(allocated(hgt)) deallocate(hgt)
                case(0, 1)
                    call letkf_scatter_vcoord(wrf % ph,  alt, Vstag)
                end select
            end if


            do j = 1, cpu(myid)%loc_ny
            do i = 1, cpu(myid)%loc_nx
                xyz(1:2) = proj % lonlat_to_xy(lon(i,j), lat(i,j))

                do k = 1, nz 
                    xyz(3)  = alt(i,j,k)

                    !localization
                    fail(1) = get_lz(  'gts',   gts_lz, ivar, xyz)
                    fail(2) = get_lz('radar', radar_lz, ivar, xyz)

                    if(all(fail)) cycle

                    !collect data that we need
                    call letkf_yoyb(ivar,  gts,   gts_lz, &
                                           rad, radar_lz, yo, yb)

                    if( allocated(yo) ) then
                        !main letkf process
                        xb           = var(i,j,k,:)
                        var(i,j,k,:) = letkf_solve(xb, yo, yb, inflat,           &
                                            use_RTPP(ivar), RTPP_Alpha(ivar),    &
                                            use_RTPS(ivar), RTPS_Alpha(ivar))

                        deallocate(yo, yb)
                    end if

                    if(allocated(  gts_lz)) deallocate(  gts_lz)
                    if(allocated(radar_lz)) deallocate(radar_lz)
                end do
            end do
            end do


            select case(trim(var_update(ivar)))
            case ('U')
                call letkf_gather_grid(var, wrf % u, Hstag)
            case ('V')                             
                call letkf_gather_grid(var, wrf % v, Hstag)
            case ('W')                             
                call letkf_gather_grid(var, wrf % w, Hstag)
            case ('T')                             
                call letkf_gather_grid(var, wrf % t, Hstag)
            case ('QVAPOR')
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf %  qv, Hstag)
            case ('QRAIN')                
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf %  qr, Hstag)
            case ('QSNOW')               
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf %  qs, Hstag)
            case ('QGRAUP')               
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf %  qg, Hstag)
            case ('QHAIL')                
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf %  qh, Hstag)
            case ('QNRAIN')               
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf % nqr, Hstag)
            case ('QNSNOW')               
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf % nqs, Hstag)
            case ('QNGRAUPEL')              
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf % nqg, Hstag)
            case ('QNHAIL')               
                call letkf_tune_q(nz, var)
                call letkf_gather_grid(var, wrf % nqh, Hstag)
            case ('P')               
                call letkf_gather_grid(var, wrf % p,   Hstag)
            case ('MU')               
                call letkf_gather_grid(var,      mu,   Hstag)
                if(allocated(mu)) then
                    wrf % mu = mu(:,:,1)
                    deallocate(mu)
                end if
            case ('PH')               
                call letkf_gather_grid(var, wrf % ph,  Hstag)
            case default
                stop "Need to code for unknown variable"
            end select


            deallocate(var)
            call destroy_tree

        end do update
    end subroutine letkf_driver

    subroutine letkf_yoyb(ivar,  gts,   gts_lz, &
                                 rad, radar_lz, yo, yb)
        implicit none
        integer,                            intent(in)                :: ivar
        type(wrfda_gts),                    intent(in),  target       :: gts
        type(cwb_radar),                    intent(in),  target       :: rad
        type(lz_structure), dimension(:),   intent(in),  allocatable  :: gts_lz
        type(lz_structure), dimension(:),   intent(in),  allocatable  :: radar_lz
        real,               dimension(:),   intent(out), allocatable  :: yo
        real,               dimension(:,:), intent(out), allocatable  :: yb

        type data_list
            real                            :: yo_prime, error
            real, dimension(:), allocatable :: yb_prime
            type(data_list), pointer        :: next => null()
        end type data_list

        !local
        type(data_list),     pointer         :: head      => null()
        type(data_list),     pointer         :: current   => null()
        type(gts_config),    pointer         :: obs_nml   => null()
        type(radar_variable_config), pointer :: radar_var => null()

        real,    dimension(nmember)         :: bg
        real                                :: mean, omm, std, err
        real                                :: error_inv
        integer                             :: i, j, k, obs_type, total
        integer                             :: nvar
        real,    dimension(:), allocatable  :: err_muti, err_rej
        logical, dimension(:), allocatable  :: is_assim

        total = 0

        if(allocated(gts_lz)) then
            do i = 1, size(gts_lz)
                obs_type = gts_lz(i) % mytype

                if(allocated(gts_lz(i)%idx)) then
                    select case (obs_type)
                    case (synop, ships, metar)
                        select case(obs_type)
                        case (synop)
                            obs_nml  => synop_nml
                        case (ships)
                            obs_nml  => ships_nml
                        case (metar)
                            obs_nml  => metar_nml
                        end select

                        nvar = 5

                        allocate(err_muti (nvar), &
                                 err_rej  (nvar), &
                                 is_assim (nvar))

                        if(obs_nml%hclr(ivar) > 0.) then
                            is_assim  = [ obs_nml%u%is_assim(ivar), &
                                          obs_nml%v%is_assim(ivar), &  
                                          obs_nml%t%is_assim(ivar), &  
                                          obs_nml%p%is_assim(ivar), &  
                                          obs_nml%q%is_assim(ivar) ]
                        else
                            is_assim  = .false.
                        end if

                        err_muti  = [ obs_nml%u%err_muti, &
                                      obs_nml%v%err_muti, &
                                      obs_nml%t%err_muti, &
                                      obs_nml%p%err_muti, &
                                      obs_nml%q%err_muti ]
                        err_rej   = [ obs_nml%u%err_rej, &
                                      obs_nml%v%err_rej, &
                                      obs_nml%t%err_rej, &
                                      obs_nml%p%err_rej, &
                                      obs_nml%q%err_rej ]
                    case (sound)
                        obs_nml  => sound_nml

                        nvar = 4

                        allocate(err_muti (nvar), &
                                 err_rej  (nvar), &
                                 is_assim (nvar))

                        if(obs_nml%hclr(ivar) > 0.) then
                            is_assim  = [ obs_nml%u%is_assim(ivar), &
                                          obs_nml%v%is_assim(ivar), &  
                                          obs_nml%t%is_assim(ivar), &  
                                          obs_nml%q%is_assim(ivar) ]
                        else
                            is_assim  = .false.
                        end if

                        err_muti  = [ obs_nml%u%err_muti, &
                                      obs_nml%v%err_muti, &
                                      obs_nml%t%err_muti, &
                                      obs_nml%q%err_muti ]
                        err_rej   = [ obs_nml%u%err_rej, &
                                      obs_nml%v%err_rej, &
                                      obs_nml%t%err_rej, &
                                      obs_nml%q%err_rej ]
                    case (gpspw)
                        obs_nml  => gpspw_nml

                        nvar = 1

                        allocate(err_muti (nvar), &
                                 err_rej  (nvar), &
                                 is_assim (nvar))

                        if(obs_nml%hclr(ivar) > 0.) then
                            is_assim  = obs_nml%tpw%is_assim(ivar)
                        else
                            is_assim  = .false.
                        end if

                        err_muti  = obs_nml%tpw%err_muti
                        err_rej   = obs_nml%tpw%err_rej
                    end select


                    do j = 1, size(gts_lz(i)%idx)  !loop over available data
                        associate( idx     => gts_lz(i) % idx(j),               &
                                   r2      => gts_lz(i) % r2(j),                & 
                                   qc      => gts % platform(obs_type) % qc,    &
                                   error   => gts % platform(obs_type) % error, &
                                   obs     => gts % platform(obs_type) % obs,   &
                                   hdxb    => gts % platform(obs_type) % hdxb )
                            do k = 1, nvar  !number of observation variables
                                if(is_assim(k) .and. any(qc(k,idx,:) >= 0)) then
                                    bg    = hdxb(k,idx,:)
                                    mean  = sum(bg) * nmember_inv
                                    bg    = bg - mean
                                    omm   = obs(k,idx) - mean
                                    std   = sqrt(dot_product(bg,bg) * nmember_1_inv)
                                    err   = error(k,idx) * err_muti(k)

                                    if(abs(omm) > sqrt(std * std + err * err) * err_rej(k)) cycle

                                    !variance = variance * exp(r^2 / (2 * rloc^2))
                                    !note! We multiply localization funcion  on error not variance
                                    !so the weighting function should be square root
                                    !that's why it is 0.25 not 0.5 here
                                    if(weight_function /= 1) then
                                        error_inv = 1.0 / (err * exp( 0.25 * r2 ))
                                    else
                                        !Gaspari_Cohn_1999's input must be r
                                        !Gaspari_Cohn_1999 is applied on variance, but here is error,
                                        !so you need square root 
                                        error_inv = sqrt(Gaspari_Cohn_1999(sqrt(r2))) / err
                                    end if
                                    omm       = omm * error_inv
                                    bg        = bg  * error_inv
                                    call append(head, current, omm, bg)
 
                                    total = total+1
                                end if
                            end do
                        end associate
                    end do

                    nullify(obs_nml)
                    deallocate(err_muti, err_rej, is_assim)
                end if
            end do
        end if

        !=======================================================================================

        if(allocated(radar_lz)) then
            allocate(err_muti(1), err_rej(1), is_assim(1))

            do i = 1, size(radar_lz)
                obs_type = radar_lz(i) % mytype

                if(allocated(radar_lz(i)%idx)) then
                    select case (obs_type)
                    case (dbz)
                        radar_var => radar_nml%dbz
                    case (vr)
                        radar_var => radar_nml%vr
                    case (zdr)
                        radar_var => radar_nml%zdr
                    case (kdp)
                        radar_var => radar_nml%kdp
                    end select

                    is_assim(1)  =  radar_var%hclr(ivar) > 0.
                    err_muti(1)  =  radar_var%error
                    err_rej(1)   =  radar_var%err_rej
 
                    if(is_assim(1)) then
                        do j = 1, size(radar_lz(i)%idx)
                            associate( idx     => radar_lz(i) % idx(j),            &
                                       r2      => radar_lz(i) % r2(j),             & 
                                       hdxb    => rad % radarobs(obs_type) % hdxb, &
                                       obs     => rad % radarobs(obs_type) % obs )
                                bg    = hdxb(idx,:)
                                mean  = sum(bg) * nmember_inv
                                bg    = bg - mean
                                omm   = obs(idx) - mean
                                std   = sqrt(dot_product(bg,bg) * nmember_1_inv)
                                err   = err_muti(1)

                                if(obs_type == dbz) then
                                    if(abs(omm) > sqrt(std * std + err * err) * err_rej(1) .and. &
                                       obs(idx) /= norain_value) cycle
                                    if(obs(idx) == norain_value .and. mean == norain_value) cycle
                                else
                                    if(abs(omm) > sqrt(std * std + err * err) * err_rej(1)) cycle
                                end if

                                !variance = variance * exp(r^2 / (2 * rloc^2))
                                !note! We multiply localization funcion  on error not variance
                                !so the weighting function should be square root
                                !that's why it is 0.25 not 0.5 here
                                if(weight_function /= 1) then
                                    error_inv = 1.0 / (err * exp( 0.25 * r2 ))
                                else
                                    !Gaspari_Cohn_1999's input must be r
                                    !Gaspari_Cohn_1999 is applied on variance, but here is error,
                                    !so you need square root 
                                    error_inv = sqrt(Gaspari_Cohn_1999(sqrt(r2))) / err
                                end if
                                omm       = omm * error_inv
                                bg        = bg  * error_inv
                                call append(head, current, omm, bg)
                            end associate
 
                            total = total+1
                        end do
                    end if

                    nullify(radar_var)
                end if
            end do

            deallocate(err_muti, err_rej, is_assim)
        end if

        !==================================================================================

        if(total == 0) return
        allocate(yo(total), yb(nmember, total))

        current => head
        do i = 1, total
            if(i > 1) current => current % next

            yo(i)     = current % yo_prime
            yb(:,i)   = current % yb_prime
        end do

        nullify(current)
        call destroy(head)

    contains
 
        subroutine append(mylist, ptr, yo_prime, yb_prime)
            implicit none
            type(data_list), pointer, intent(in out) :: mylist, ptr

            !local
            real, dimension(:),       intent(in)     :: yb_prime
            real,                     intent(in)     :: yo_prime

            if(associated(mylist)) then
                allocate(ptr%next)
                ptr => ptr%next
            else
                allocate(mylist)
                ptr => mylist
            end if

            allocate(ptr % yb_prime(nmember))

            ptr % yo_prime = yo_prime
            ptr % yb_prime = yb_prime
        end subroutine append

        subroutine destroy(mylist)
            implicit none
            type(data_list), pointer, intent(in out) :: mylist
            type(data_list), pointer                 :: current, next

            current => mylist
            do while(associated(current))
                next => current % next
                deallocate(current)
                current => next
            end do

            nullify(mylist)
        end subroutine destroy

    end subroutine letkf_yoyb
 

    function letkf_solve(xb, yo, yb, inflat, use_RTPP, RTPP_Alpha, use_RTPS, RTPS_Alpha) result(xa)
        implicit none
        real,   dimension(nmember), intent(in) :: xb
        real,   dimension(:),       intent(in) :: yo
        real,   dimension(:,:),     intent(in) :: yb
        real,                       intent(in) :: inflat, RTPP_Alpha, RTPS_Alpha
        logical,                    intent(in) :: use_RTPP, use_RTPS
        real,   dimension(nmember)             :: xa

        !local
        integer                                :: i, nobs
#ifdef REAL64
        real*8                                 :: xb_mean, inflat_r8
        real*8, dimension(nmember)             :: xb_prime, wbar
        real*8, dimension(nmember, nmember)    :: identity, wbar2d, w
        real*8, dimension(:),   allocatable    :: yo_r8
        real*8, dimension(:,:), allocatable    :: yb_r8
#else
        real                                   :: xb_mean
        real,   dimension(nmember)             :: xb_prime, wbar
        real,   dimension(nmember, nmember)    :: identity, wbar2d, w
#endif

        !RTPP, RTPS
        real                                   :: xa_std, xb_std
        real                                   :: xa_mean
        real,   dimension(nmember)             :: xa_prime


#ifdef REAL64
        identity = 0d0
        do concurrent(i=1:nmember)
            identity(i,i) = 1d0
        end do
#else
        identity = 0.0
        do concurrent(i=1:nmember)
            identity(i,i) = 1.0
        end do
#endif

        nobs = size(yo)

#ifdef REAL64
        allocate(yb_r8(nmember, nobs), yo_r8(nobs))

        !real to double
        inflat_r8  = dble(inflat)
        yb_r8      = dble(yb)
        yo_r8      = dble(yo)

        call dsyrk('l', 'n', nmember, nobs, 1d0, yb_r8, nmember, inflat_r8, identity, nmember)
        call inverse_matrix(identity, nmember)
        call dgemv('n', nmember, nobs, 1d0, yb_r8, nmember, yo_r8, 1, 0d0, wbar, 1)
        call dsymv('l', nmember, 1d0, identity, nmember, wbar, 1, 0d0, xb_prime, 1)

        deallocate(yb_r8, yo_r8)
#else
        call ssyrk('l', 'n', nmember, nobs, 1.0, yb, nmember, inflat, identity, nmember)
        call inverse_matrix(identity, nmember)
        call sgemv('n', nmember, nobs, 1.0, yb, nmember, yo, 1, 0.0, wbar, 1)
        call ssymv('l', nmember, 1.0, identity, nmember, wbar, 1, 0.0, xb_prime, 1)
#endif

        wbar2d = spread(xb_prime, 2, nmember)

        call sqrt_matrix(w, nmember)
#ifdef REAL64
        call daxpy(nmember*nmember, sqrt(dble(nmember-1)), w, 1, wbar2d, 1)
#else
        call saxpy(nmember*nmember, sqrt(float(nmember-1)), w, 1, wbar2d, 1)
#endif

        xb_mean  = sum(xb) * nmember_inv
        xb_prime = xb - xb_mean
        wbar     = xb_mean
#ifdef REAL64
        call dgemv('t', nmember, nmember, 1d0, wbar2d, nmember, xb_prime, 1, 1d0, wbar, 1)
#else
        call sgemv('t', nmember, nmember, 1.0, wbar2d, nmember, xb_prime, 1, 1.0, wbar, 1)
#endif
        xa       = wbar


        !=====================================================================================

        if(use_RTPP .or. use_RTPS) then
            xa_mean  = sum(xa) * nmember_inv
            xa_prime = xa - xa_mean

            if(use_RTPP) then
                xa_prime = (1.0 - RTPP_alpha) * xa_prime + RTPP_alpha * xb_prime
            end if
            if(use_RTPS) then
                xb_std   = dot_product(xb_prime, xb_prime)
                xa_std   = dot_product(xa_prime, xa_prime)
                xa_prime = xa_prime * (RTPS_alpha * sqrt(xb_std / xa_std) - RTPS_alpha + 1.0)
            end if

            xa = xa_mean + xa_prime
        end if

    end function letkf_solve

    subroutine letkf_tune_q(nz, q)
        implicit none
        integer,                  intent(in)                :: nz
        real, dimension(:,:,:,:), intent(inout), contiguous :: q

        !local
        integer                     :: i, j, k
        real                        :: ratio
        real,    dimension(nmember) :: var
        logical, dimension(nmember) :: is_positive

        do k = 1, nz
        do j = 1, cpu(myid) % loc_ny
        do i = 1, cpu(myid) % loc_nx

            var              = q(i, j, k, :)
            is_positive      = var > 0.
            ratio            = sum(var) / sum(var, mask=is_positive)

            !tune q to eliminate negative value but preserve mean
            where( var < 0. )
                var = 0.
            else where
                var = ratio * var
            end where

            q(i, j, k, :) = var

        end do
        end do
        end do
    end subroutine letkf_tune_q

    function check_coordinate(previous, now) result(reset)
        implicit none
        integer, intent(inout) :: previous
        integer, intent(in)    :: now
        logical                :: reset

        if(previous /= now) then
            previous = now
            reset = .true.
        else
            reset = .false.
        end if
    end function check_coordinate

end module letkf_core
