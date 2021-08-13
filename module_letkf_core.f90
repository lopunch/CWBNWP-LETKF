module letkf_core

    use localization,    only : lz_structure, get_hlz, get_vlz, destroy_vlz
    use grid,            only : grid_structure
    use gts_omboma,      only : wrfda_gts
    use simulated_radar, only : cwb_radar
    use param
    use mpi_util
    use eigen
    use config

    implicit none

    private
    public  :: letkf_driver

contains

    subroutine letkf_driver(wrf, gts, rad)
        implicit none

        type(grid_structure),         intent(inout)   :: wrf
        type(wrfda_gts),              intent(in   )   :: gts
        type(cwb_radar),              intent(in   )   :: rad

        !localization
        type(lz_structure), dimension(:), allocatable ::   gts_lz
        type(lz_structure), dimension(:), allocatable :: radar_lz

        !local
        integer                                       :: ivar
        integer                                       :: i, j, k
        integer                                       :: nx, ny, nz
        real                                          :: inflat
        real, dimension(nmember)                      :: xb
        real, dimension(:),       allocatable         :: yo
        real, dimension(:,:),     allocatable         :: yb, lat, lon
        real, dimension(:,:,:),   allocatable         :: alt
        real, dimension(:,:,:,:), allocatable         :: var
        logical, dimension(2)                         :: fail_h, fail_v

        call wait_jobs(1, 3)
        nx = wrf % nx
        ny = wrf % ny

        call letkf_local_info(nx, ny)

        update : do ivar = 1, max_vars
            if(len_trim(var_update(ivar)) == 0) exit
            if(is_root) print '(f7.3,1x,a,1x,10a)', timer(), "sec ==========> update", var_update(ivar)

            inflat = (nmember-1) / multi_infl(ivar)
            nz     = wrf % nz

            select case(trim(var_update(ivar)))
            case ('U')
                call letkf_scatter_hcoord(wrf % xlat_u, lat, &
                                          wrf % xlon_u, lon, 1)
                call letkf_scatter_vcoord(wrf % ph, nz, alt   )
                call letkf_scatter_grid  (wrf % u,  nz, var, 1)
            case ('V')
                call letkf_scatter_hcoord(wrf % xlat_v, lat, &
                                          wrf % xlon_v, lon, 2)
                call letkf_scatter_vcoord(wrf % ph, nz, alt   )
                call letkf_scatter_grid(  wrf % v,  nz, var, 2)
            case ('W')
                nz = nz+1
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz, alt, 1)
                call letkf_scatter_grid(  wrf % w,  nz, var)
            case ('T')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz, alt)
                call letkf_scatter_grid(  wrf % t,  nz, var)
            case ('QVAPOR')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz, alt)
                call letkf_scatter_grid(  wrf % qv, nz, var)
            case ('QRAIN')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz,  alt)
                call letkf_scatter_grid(  wrf % qr, nz,  var)
            case ('QSNOW')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz,  alt)
                call letkf_scatter_grid(  wrf % qs, nz,  var)
            case ('QGRAUP')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz,  alt)
                call letkf_scatter_grid(  wrf % qg, nz,  var)
            case ('QHAIL')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph, nz,  alt)
                call letkf_scatter_grid(  wrf % qh, nz,  var)
            case ('QNRAIN')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph,  nz,  alt)
                call letkf_scatter_grid(  wrf % nqr, nz,  var)
            case ('QNSNOW')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph,  nz,  alt)
                call letkf_scatter_grid(  wrf % nqs, nz,  var)
            case ('QNGRAUPEL')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph,  nz,  alt)
                call letkf_scatter_grid(  wrf % nqg, nz,  var)
            case ('QNHAIL')
                call letkf_scatter_hcoord(wrf % xlat, lat, &
                                          wrf % xlon, lon)
                call letkf_scatter_vcoord(wrf % ph,  nz,  alt)
                call letkf_scatter_grid(  wrf % nqh, nz,  var)
            case default
                print *, var_update(ivar)
                stop "Need to code for unknown variable"
            end select

            do j = 1, cpu(myid) % loc_ny
            do i = 1, cpu(myid) % loc_nx

                !horizontal localization
                fail_h(1) = get_hlz(  'gts',   gts_lz, ivar, lat(i,j), lon(i,j))
                fail_h(2) = get_hlz('radar', radar_lz, ivar, lat(i,j), lon(i,j))

                if(all(fail_h)) cycle

                do k = 1, nz 
                    !vertical localization
                    if(fail_h(1)) then
                        fail_v(1) = .true.
                    else
                        fail_v(1) = get_vlz(gts%platform,   gts_lz, ivar, alt(i,j,k))
                    end if

                    if(fail_h(2)) then
                        fail_v(2) = .true.
                    else
                        fail_v(2) = get_vlz(rad%radarobs, radar_lz, ivar, alt(i,j,k))
                    end if

                    if(all(fail_v)) cycle

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

                    call destroy_vlz(  gts_lz)
                    call destroy_vlz(radar_lz)
                end do

                if(allocated(  gts_lz)) deallocate(  gts_lz)
                if(allocated(radar_lz)) deallocate(radar_lz)
            end do
            end do

            select case(trim(var_update(ivar)))
            case ('U')
                call letkf_gather_grid(var, wrf % u, 1)
            case ('V')                             
                call letkf_gather_grid(var, wrf % v, 2)
            case ('W')                             
                call letkf_gather_grid(var, wrf % w)
            case ('T')                             
                call letkf_gather_grid(var, wrf % t)
            case ('QVAPOR')
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf %  qv)
            case ('QRAIN')                
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf %  qr)
            case ('QSNOW')               
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf %  qs)
            case ('QGRAUP')               
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf %  qg)
            case ('QHAIL')                
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf %  qh)
            case ('QNRAIN')               
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf % nqr)
            case ('QNSNOW')               
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf % nqs)
            case ('QNGRAUPEL')              
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf % nqg)
            case ('QNHAIL')               
                where(var < 0.) var = 0.
                call letkf_gather_grid(var, wrf % nqh)
            case default
                stop "Need to code for unknown variable"
            end select

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
        real                                :: mean, omm, std, err, f1, f2
        real                                :: hclr_inv, vclr_inv,  hclr_rej, vclr_rej, hdist, vdist
        real                                :: error_inv
        integer                             :: i, j, k, obs_type, idx, total
        real,    dimension(:), allocatable  :: err_muti, err_rej
        logical, dimension(:), allocatable  :: is_assim

        total = 0

        if(allocated(gts_lz)) then
            do i = 1, size(gts_lz)
                obs_type = gts_lz(i) % mytype

                if(allocated(gts_lz(i)%vert)) then
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

                        allocate(err_muti (5), &
                                 err_rej  (5), &
                                 is_assim (5))

                        is_assim  = [ obs_nml%u%is_assim(ivar), &
                                      obs_nml%v%is_assim(ivar), &  
                                      obs_nml%t%is_assim(ivar), &  
                                      obs_nml%p%is_assim(ivar), &  
                                      obs_nml%q%is_assim(ivar) ]
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

                        allocate(err_muti (4), &
                                 err_rej  (4), &
                                 is_assim (4))

                        is_assim  = [ obs_nml%u%is_assim(ivar), &
                                      obs_nml%v%is_assim(ivar), &  
                                      obs_nml%t%is_assim(ivar), &  
                                      obs_nml%q%is_assim(ivar) ]
                        err_muti  = [ obs_nml%u%err_muti, &
                                      obs_nml%v%err_muti, &
                                      obs_nml%t%err_muti, &
                                      obs_nml%q%err_muti ]
                        err_rej   = [ obs_nml%u%err_rej, &
                                      obs_nml%v%err_rej, &
                                      obs_nml%t%err_rej, &
                                      obs_nml%q%err_rej ]
                    end select

                    hclr_inv  =   1.0 / (obs_nml%hclr(ivar) * 1e3)
                    hclr_rej  =   obs_nml%hclr_rej
                    vclr_inv  =   1.0 / (obs_nml%vclr(ivar) * 1e3)
                    vclr_rej  =   obs_nml%vclr_rej

                    do j = 1, size(gts_lz(i)%vert%idx)

                        idx   = gts_lz(i)%vert%idx(j)
                        hdist = gts_lz(i)%vert%hdistance(j)
                        vdist = gts_lz(i)%vert%vdistance(j)

                        associate( qc    => gts % platform(obs_type) % qc,    &
                                   error => gts % platform(obs_type) % error, &
                                   obs   => gts % platform(obs_type) % obs,   &
                                   hdxb  => gts % platform(obs_type) % hdxb )
                            do k = 1, size(obs, 1)
                                if(is_assim(k) .and. any(qc(k,idx,:) >= 0)) then
                                    bg    = hdxb(k,idx,:)
                                    mean  = sum(bg) * nmember_inv
                                    bg    = bg - mean
                                    omm   = obs(k,idx) - mean
                                    std   = sqrt(dot_product(bg,bg) * nmember_1_inv)
                                    err   = error(k,idx) * err_muti(k)

                                    if(abs(omm) > sqrt(std * std + err * err) * err_rej(k)) cycle

                                    f1    = hdist * hclr_rej * hclr_inv
                                    f2    = vdist * vclr_rej * vclr_inv

                                    !variance = variance * exp(r^2 / (2 * rloc^2))
                                    !note! We multiply localization funcion  on error not variance
                                    !so the weighting function should be square root
                                    !that's why it is 0.25 not 0.5 here
                                    error_inv = 1.0 / (err * exp( 0.25 * (f1*f1 + f2*f2) ))
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

                if(allocated(radar_lz(i)%vert)) then
                    hclr_rej  =   radar_nml%hclr_rej
                    vclr_rej  =   radar_nml%vclr_rej

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

                    hclr_inv     =  1.0 / (radar_var%hclr(ivar) * 1e3)
                    vclr_inv     =  1.0 / (radar_var%vclr(ivar) * 1e3)
                    is_assim(1)  =  radar_var%is_assim(ivar)
                    err_muti(1)  =  radar_var%error
                    err_rej(1)   =  radar_var%err_rej
 
                    do j = 1, size(radar_lz(i)%vert%idx)
                        idx   = radar_lz(i)%vert%idx(j)
                        hdist = radar_lz(i)%vert%hdistance(j)
                        vdist = radar_lz(i)%vert%vdistance(j)

                        if(is_assim(1)) then
                            associate( hdxb => rad % radarobs(obs_type) % hdxb, &
                                       obs  => rad % radarobs(obs_type) % obs )
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
                            end associate

                            f1    = hdist * hclr_rej * hclr_inv
                            f2    = vdist * vclr_rej * vclr_inv
                            !variance = variance * exp(r^2 / (2 * rloc^2))
                            !note! We multiply localization funcion  on error not variance
                            !so the weighting function should be square root
                            !that's why it is 0.25 not 0.5 here
                            error_inv = 1.0 / (err * exp( 0.25 * (f1*f1 + f2*f2) ))
                            omm       = omm * error_inv
                            bg        = bg  * error_inv
                            call append(head, current, omm, bg)
 
                            total = total+1
                        end if
                    end do

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
        integer                                :: i, nobs, shp(2)
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

        shp  = shape(yb)
        nobs = shp(2)

#ifdef REAL64
        allocate(yb_r8(shp(1), shp(2)), yo_r8(nobs))

        !real to double
        inflat_r8  = inflat
        yb_r8      = yb
        yo_r8      = yo

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

        if(use_RTPP .or. use_RTPP) then
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

end module letkf_core
