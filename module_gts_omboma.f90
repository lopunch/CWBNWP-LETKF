module gts_omboma

    use config, only : nmember
    use param
    use mpi_util

    implicit none

    private
    public  :: variables, gts_structure, wrfda_gts

    type variables
        integer                                 :: nlev
        real, dimension(:), allocatable         :: h, pre

        !===============================================
        real,     dimension(:,:),   allocatable :: obs       ![nvar, nlev]
        real,     dimension(:,:),   allocatable :: error     ![nvar, nlev]
        real,     dimension(:,:,:), allocatable :: omb       ![nvar, nlev, nmember]      
        integer,  dimension(:,:,:), allocatable :: qc        ![nvar, nlev, nmember]
    end type variables

    type, extends(obs_structure) :: gts_structure
        character(len=5), dimension(:), allocatable  :: id

        !obs on surface
        type(variables),                allocatable  :: surf          !dim is number of data
        !obs have vertical profile
        type(variables),  dimension(:), allocatable  :: vert          !dim is number of data
    end type gts_structure

    !========================================================================================

    type wrfda_gts
        type(gts_structure), dimension(:), allocatable :: platform
    contains
        procedure, public, pass(self) :: read_data  => read_gts_omboma
        procedure, public, pass(self) :: write_data => write_gts_omboma
        procedure, public, pass(self) :: distribute => gts_distribute
    end type wrfda_gts

    !========================================================================================

    type alt_info
        real, dimension(:), allocatable              :: alt
    end type alt_info

    type alt_structure
        integer                                       :: nobs  = 0
        character(len=20), dimension(:), allocatable  :: id
        type(alt_info),    dimension(:), allocatable  :: info        !dim is number of data
    end type alt_structure
 
contains

    subroutine read_gts_omboma(self, filename, obascii)
        implicit none
        class(wrfda_gts), intent(in out) :: self
        character(len=*), intent(in)     :: filename, obascii

        !local
        type(alt_structure), dimension(num_gts_indexes)   :: gtsalt

        integer                                           :: ierr
        integer                                           :: nobs, nlev, obs_type
        character(len=20)                                 :: iv_type, obs_name

        integer                                           :: n, k
        integer                                           :: kk, l, nreq
        real, dimension(5)                                :: oma
        logical, parameter                                :: write_gts = .false.
        
        !read altitude info
        call read_alt_info(gtsalt, obascii)

        open(30, file=filename, iostat=ierr)
        if(ierr /= 0) stop "open gts_omboma error"

        if(write_gts) open(40, file='gts_out_'//filename(20:22), iostat=ierr)
        associate( platform => self % platform )
            report: do
                read(30, '(a20,i8)', iostat=ierr) iv_type, nobs
                if(ierr < 0) exit
                if(ierr > 0) stop "read gts_omboma error"
                obs_name = trim(adjustl(iv_type))

                if(write_gts) write(40, '(a20,i8)', iostat=ierr) iv_type, nobs

                select case (trim(obs_name))
                case ('synop', 'ships', 'buoy', 'metar', 'sonde_sfc', 'tamdar_sfc')
                    if(nobs > 0) then
                        select case (trim(obs_name))    
                        case ('synop') 
                            obs_type = synop 
                        case ('ships') 
                            obs_type = ships 
                        case ('buoy') 
                            obs_type = buoy
                        case ('metar') 
                            obs_type = metar
                        case ('sonde_sfc') 
                            obs_type = sonde_sfc
                        case ('tamdar_sfc') 
                            obs_type = tamdar_sfc
                        end select

                        platform(obs_type)%nobs = nobs 
                        allocate(platform(obs_type)%id (nobs), &
                                 platform(obs_type)%lat(nobs), &
                                 platform(obs_type)%lon(nobs))

                        allocate(platform(obs_type)%surf)
                        allocate(platform(obs_type)%surf%h    (nobs), &
                                 platform(obs_type)%surf%pre  (nobs), &
                                 platform(obs_type)%surf%obs  (5,nobs),         &
                                 platform(obs_type)%surf%error(5,nobs),         &
                                 platform(obs_type)%surf%omb  (5,nobs,0:nmember-1),         &
                                 platform(obs_type)%surf%qc   (5,nobs,0:nmember-1))

                        do n = 1, nobs 
                            read(30, '(2i8)') nlev, nreq
                            if(write_gts) write(40, '(2i8)') nlev, nreq

                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                kk, l, platform(obs_type)%id(n), &
                                platform(obs_type)%lat(n),  &
                                platform(obs_type)%lon(n),  &   
                                platform(obs_type)%surf%pre(n),                                             &
                                platform(obs_type)%surf%obs(1,n),   platform(obs_type)%surf%omb  (1,n,0),       &
                                platform(obs_type)%surf%qc (1,n,0), platform(obs_type)%surf%error(1,n),       &
                                oma(1),                        &
                                platform(obs_type)%surf%obs(2,n),   platform(obs_type)%surf%omb  (2,n,0),       &
                                platform(obs_type)%surf%qc (2,n,0), platform(obs_type)%surf%error(2,n),       &
                                oma(2),                        &
                                platform(obs_type)%surf%obs(3,n),   platform(obs_type)%surf%omb  (3,n,0),       &
                                platform(obs_type)%surf%qc (3,n,0), platform(obs_type)%surf%error(3,n),       &
                                oma(3),                      &
                                platform(obs_type)%surf%obs(4,n),   platform(obs_type)%surf%omb  (4,n,0),       &
                                platform(obs_type)%surf%qc (4,n,0), platform(obs_type)%surf%error(4,n),       &
                                oma(4),                      &
                                platform(obs_type)%surf%obs(5,n),   platform(obs_type)%surf%omb  (5,n,0),       &
                                platform(obs_type)%surf%qc (5,n,0), platform(obs_type)%surf%error(5,n),       &
                                oma(5)

                            platform(obs_type)%surf%h(n) = get_alt(gtsalt(obs_type), platform(obs_type)%id(n), 1)

                            if(write_gts) &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    kk, l, platform(obs_type)%id(n), &
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%surf%pre(n),                                             &
                                    platform(obs_type)%surf%obs(1,n),   platform(obs_type)%surf%omb  (1,n,0),       &
                                    platform(obs_type)%surf%qc (1,n,0), platform(obs_type)%surf%error(1,n),       &
                                    oma(1),                        &
                                    platform(obs_type)%surf%obs(2,n),   platform(obs_type)%surf%omb  (2,n,0),       &
                                    platform(obs_type)%surf%qc (2,n,0), platform(obs_type)%surf%error(2,n),       &
                                    oma(2),                        &
                                    platform(obs_type)%surf%obs(3,n),   platform(obs_type)%surf%omb  (3,n,0),       &
                                    platform(obs_type)%surf%qc (3,n,0), platform(obs_type)%surf%error(3,n),       &
                                    oma(3),                      &
                                    platform(obs_type)%surf%obs(4,n),   platform(obs_type)%surf%omb  (4,n,0),       &
                                    platform(obs_type)%surf%qc (4,n,0), platform(obs_type)%surf%error(4,n),       &
                                    oma(4),                      &
                                    platform(obs_type)%surf%obs(5,n),   platform(obs_type)%surf%omb  (5,n,0),       &
                                    platform(obs_type)%surf%qc (5,n,0), platform(obs_type)%surf%error(5,n),       &
                                    oma(5)
                        end do  !nobs

                        !bg = obs - omb
                        platform(obs_type)%surf%omb(:,:,0) = platform(obs_type)%surf%obs - platform(obs_type)%surf%omb(:,:,0)
                    end if
                case ('pilot', 'profiler', 'geoamv', 'qscat', 'polaramv')    
                    if(nobs > 0) then
                        select case (trim(obs_name))    
                        case ('pilot') 
                            obs_type = pilot
                        case ('profiler') 
                            obs_type = profiler
                        case ('geoamv') 
                            obs_type = geoamv
                        case ('qscat') 
                            obs_type = qscat
                        case ('polaramv') 
                            obs_type = polaramv
                        end select

                        platform(obs_type)%nobs  = nobs 
                        allocate(platform(obs_type)%id  (nobs), &
                                 platform(obs_type)%lat (nobs), &
                                 platform(obs_type)%lon (nobs), &
                                 platform(obs_type)%vert(nobs))
                        
                        do n = 1, nobs 
                            read(30,'(2i8)') nlev, nreq
                            platform(obs_type)%vert(n)%nlev = nlev

                            if(write_gts) write(40,'(2i8)') nlev, nreq

                            allocate(platform(obs_type)%vert(n)%pre(nlev), &
                                     platform(obs_type)%vert(n)%h  (nlev), &
                                     platform(obs_type)%vert(n)%obs  (2,nlev), &
                                     platform(obs_type)%vert(n)%error(2,nlev), &
                                     platform(obs_type)%vert(n)%omb  (2,nlev,0:nmember-1), &
                                     platform(obs_type)%vert(n)%qc   (2,nlev,0:nmember-1))

                            do k = 1, nlev
                                read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    kk, l, platform(obs_type)%id(n), &
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%vert(n)%pre(k),                                          &
                                    platform(obs_type)%vert(n)%obs(1,k),    platform(obs_type)%vert(n)%omb  (1,k,0),       &
                                    platform(obs_type)%vert(n)%qc (1,k,0),  platform(obs_type)%vert(n)%error(1,k),       &
                                    oma(1),                                     &
                                    platform(obs_type)%vert(n)%obs(2,k),    platform(obs_type)%vert(n)%omb  (2,k,0),       &
                                    platform(obs_type)%vert(n)%qc (2,k,0),  platform(obs_type)%vert(n)%error(2,k),       &
                                    oma(2)

                                platform(obs_type)%vert(n)%h(k) = get_alt(gtsalt(obs_type), platform(obs_type)%id(n), k)
                               !if(k == 1) platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                                if(write_gts)  &
                                    write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                        kk, l, platform(obs_type)%id(n), &
                                        platform(obs_type)%lat(n),  &
                                        platform(obs_type)%lon(n),  &   
                                        platform(obs_type)%vert(n)%pre(k),                                          &
                                        platform(obs_type)%vert(n)%obs(1,k),    platform(obs_type)%vert(n)%omb  (1,k,0),       &
                                        platform(obs_type)%vert(n)%qc (1,k,0),  platform(obs_type)%vert(n)%error(1,k),       &
                                        oma(1),                                     &
                                        platform(obs_type)%vert(n)%obs(2,k),    platform(obs_type)%vert(n)%omb  (2,k,0),       &
                                        platform(obs_type)%vert(n)%qc (2,k,0),  platform(obs_type)%vert(n)%error(2,k),       &
                                        oma(2)

                            end do   !nlev

                            !bg = obs - omb
                            platform(obs_type)%vert(n)%omb(:,:,0) = platform(obs_type)%vert(n)%obs - platform(obs_type)%vert(n)%omb(:,:,0)
                        end do  !nobs
                    end if
                case ('gpspw')
                    if(nobs > 0) then
                        obs_type = gpspw

                        platform(obs_type)%nobs  = nobs 
                        allocate(platform(obs_type)%id (nobs), &
                                 platform(obs_type)%lat(nobs), &
                                 platform(obs_type)%lon(nobs))

                        allocate(platform(obs_type)%surf)
                        allocate(platform(obs_type)%surf%h(nobs),   &
                                 platform(obs_type)%surf%obs  (1,nobs), &
                                 platform(obs_type)%surf%error(1,nobs), &
                                 platform(obs_type)%surf%omb  (1,nobs,0:nmember-1), &
                                 platform(obs_type)%surf%qc   (1,nobs,0:nmember-1))
                        
                        do n = 1, nobs 
                            read(30,'(2i8)') nlev, nreq
                            if(write_gts) write(40,'(2i8)') nlev, nreq

                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                kk, l, platform(obs_type)%id(n), &
                                platform(obs_type)%lat(n),  &
                                platform(obs_type)%lon(n),  &   
                                platform(obs_type)%surf%h  (n),                                              &
                                platform(obs_type)%surf%obs(1,n),    platform(obs_type)%surf%omb  (1,n,0),     &
                                platform(obs_type)%surf%qc (1,n,0),  platform(obs_type)%surf%error(1,n),     &
                                oma(1)

                           !platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                            if(write_gts) &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    kk, l, platform(obs_type)%id(n), &
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%surf%h  (n),                                              &
                                    platform(obs_type)%surf%obs(1,n),    platform(obs_type)%surf%omb  (1,n,0),     &
                                    platform(obs_type)%surf%qc (1,n,0),  platform(obs_type)%surf%error(1,n),     &
                                    oma(1)

                        end do  !nobs

                        !bg = obs - omb
                        platform(obs_type)%surf%omb(:,:,0) = platform(obs_type)%surf%obs - platform(obs_type)%surf%omb(:,:,0)
                    end if
                case ('sound', 'tamdar', 'airep')
                    if(nobs > 0) then
                        select case (trim(obs_name))    
                        case ('sound') 
                            obs_type = sound
                        case ('tamdar') 
                            obs_type = tamdar
                        case ('airep') 
                            obs_type = airep
                        end select

                        platform(obs_type)%nobs  = nobs 
                        allocate(platform(obs_type)%id  (nobs), &
                                 platform(obs_type)%lat (nobs), &
                                 platform(obs_type)%lon (nobs), &
                                 platform(obs_type)%vert(nobs))
                        
                        do n = 1, nobs 
                            read(30,'(2i8)') nlev, nreq
                            platform(obs_type)%vert(n)%nlev = nlev

                            if(write_gts) write(40,'(2i8)') nlev, nreq

                            allocate(platform(obs_type)%vert(n)%pre(nlev), &
                                     platform(obs_type)%vert(n)%h  (nlev), &
                                     platform(obs_type)%vert(n)%obs  (4,nlev), &
                                     platform(obs_type)%vert(n)%error(4,nlev), &
                                     platform(obs_type)%vert(n)%omb  (4,nlev,0:nmember-1), &
                                     platform(obs_type)%vert(n)%qc   (4,nlev,0:nmember-1))

                            do k = 1, nlev
                                read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    kk, l, platform(obs_type)%id(n), &
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%vert(n)%pre(k),                                          &
                                    platform(obs_type)%vert(n)%obs(1,k),    platform(obs_type)%vert(n)%omb  (1,k,0),       &
                                    platform(obs_type)%vert(n)%qc (1,k,0),  platform(obs_type)%vert(n)%error(1,k),       &
                                    oma(1),                                     &
                                    platform(obs_type)%vert(n)%obs(2,k),    platform(obs_type)%vert(n)%omb  (2,k,0),       &
                                    platform(obs_type)%vert(n)%qc (2,k,0),  platform(obs_type)%vert(n)%error(2,k),       &
                                    oma(2),                                     &
                                    platform(obs_type)%vert(n)%obs(3,k),    platform(obs_type)%vert(n)%omb  (3,k,0),       &
                                    platform(obs_type)%vert(n)%qc (3,k,0),  platform(obs_type)%vert(n)%error(3,k),       &
                                    oma(1),                                     &
                                    platform(obs_type)%vert(n)%obs(4,k),    platform(obs_type)%vert(n)%omb  (4,k,0),       &
                                    platform(obs_type)%vert(n)%qc (4,k,0),  platform(obs_type)%vert(n)%error(4,k),       &
                                    oma(2)

                                platform(obs_type)%vert(n)%h(k) = get_alt(gtsalt(obs_type), platform(obs_type)%id(n), k)
                               !if(k == 1) platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                                if(write_gts) &
                                    write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                        kk, l, platform(obs_type)%id(n), &
                                        platform(obs_type)%lat(n),  &
                                        platform(obs_type)%lon(n),  &   
                                        platform(obs_type)%vert(n)%pre(k),                                          &
                                        platform(obs_type)%vert(n)%obs(1,k),    platform(obs_type)%vert(n)%omb  (1,k,0),       &
                                        platform(obs_type)%vert(n)%qc (1,k,0),  platform(obs_type)%vert(n)%error(1,k),       &
                                        oma(1),                                     &
                                        platform(obs_type)%vert(n)%obs(2,k),    platform(obs_type)%vert(n)%omb  (2,k,0),       &
                                        platform(obs_type)%vert(n)%qc (2,k,0),  platform(obs_type)%vert(n)%error(2,k),       &
                                        oma(2),                                     &
                                        platform(obs_type)%vert(n)%obs(3,k),    platform(obs_type)%vert(n)%omb  (3,k,0),       &
                                        platform(obs_type)%vert(n)%qc (3,k,0),  platform(obs_type)%vert(n)%error(3,k),       &
                                        oma(1),                                     &
                                        platform(obs_type)%vert(n)%obs(4,k),    platform(obs_type)%vert(n)%omb  (4,k,0),       &
                                        platform(obs_type)%vert(n)%qc (4,k,0),  platform(obs_type)%vert(n)%error(4,k),       &
                                        oma(2)

                            end do   !nlev

                            !bg = obs - omb
                            platform(obs_type)%vert(n)%omb(:,:,0) = platform(obs_type)%vert(n)%obs - platform(obs_type)%vert(n)%omb(:,:,0)
                        end do  !nobs
                    end if
                case ('gpsref')
                    if(nobs > 0) then
                        obs_type = gpsref

                        platform(obs_type)%nobs  = nobs 
                        allocate(platform(obs_type)%id  (nobs), &
                                 platform(obs_type)%lat (nobs), &
                                 platform(obs_type)%lon (nobs), &
                                 platform(obs_type)%vert(nobs))

                        do n = 1, nobs 
                            read(30,'(2i8)') nlev, nreq
                            platform(obs_type)%vert(n)%nlev = nlev
                        
                            if(write_gts) write(40,'(2i8)') nlev, nreq

                            allocate(platform(obs_type)%vert(n)%h(nlev), &
                                     platform(obs_type)%vert(n)%obs  (1,nlev), &
                                     platform(obs_type)%vert(n)%error(1,nlev), &
                                     platform(obs_type)%vert(n)%omb  (1,nlev,0:nmember-1), &
                                     platform(obs_type)%vert(n)%qc   (1,nlev,0:nmember-1))

                            do k = 1, nlev
                                read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    kk, l, platform(obs_type)%id(n), &
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%vert(n)%h  (k),                                              &
                                    platform(obs_type)%vert(n)%obs(1,k),    platform(obs_type)%vert(n)%omb  (1,k,0),     &
                                    platform(obs_type)%vert(n)%qc (1,k,0),  platform(obs_type)%vert(n)%error(1,k),     &
                                    oma(1)

                               !platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                                if(write_gts) &
                                    write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                        kk, l, platform(obs_type)%id(n), &
                                        platform(obs_type)%lat(n),  &
                                        platform(obs_type)%lon(n),  &   
                                        platform(obs_type)%vert(n)%h  (k),                                              &
                                        platform(obs_type)%vert(n)%obs(1,k),    platform(obs_type)%vert(n)%omb  (1,k,0),     &
                                        platform(obs_type)%vert(n)%qc (1,k,0),  platform(obs_type)%vert(n)%error(1,k),     &
                                        oma(1)
                            end do

                            !bg = obs - omb
                            platform(obs_type)%vert(n)%omb(:,:,0) = platform(obs_type)%vert(n)%obs - platform(obs_type)%vert(n)%omb(:,:,0)
                        end do  !nobs
                    end if
                end select
            end do report
        end associate

        close(30)
        if(write_gts) close(40)
    end subroutine read_gts_omboma

    subroutine gts_distribute(self, root)

        implicit none
        class(wrfda_gts), intent(in out)               :: self
        integer,          intent(in)                   :: root

        !local
        integer                                        :: obs_type, nobs, nlev
        integer                                        :: n, k, m

        !local mpi 
        integer,  dimension(0:nproc-1)                 :: recvcount, displs
        integer                                        :: mpierr, sendcount
        integer,  dimension(:),     allocatable        :: i4_nobs, i4_nlev

        real,     dimension(:,:),   allocatable        :: local_omb
        integer,  dimension(:,:),   allocatable        :: local_qc

        !broadcast number of obs of each obs type
        allocate(i4_nobs(num_gts_indexes))
        if(myid == root) then
            do obs_type = 1, num_gts_indexes
                i4_nobs(obs_type) = self%platform(obs_type)%nobs
            end do
        end if

        call mpi_bcast(i4_nobs, num_gts_indexes, mpi_integer, root, mpi_comm_world, mpierr)

        !start to exchange obs info
        associate( platform => self % platform )
            do obs_type = 1, num_gts_indexes
                nobs = i4_nobs(obs_type)
                platform(obs_type)%nobs = nobs
                if(nobs > 0) then

                    if(.not. allocated(platform(obs_type)%id)) then
                        allocate(platform(obs_type)%id (nobs), &
                                 platform(obs_type)%lat(nobs), &
                                 platform(obs_type)%lon(nobs))
                    end if

                    call request_append
                    call mpi_ibcast(platform(obs_type)%id, 5*nobs, mpi_character, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(platform(obs_type)%lat,  nobs, mpi_real,      root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(platform(obs_type)%lon,  nobs, mpi_real,      root, mpi_comm_world, req_ptr%req, mpierr)

                    select case (obs_type)
                    case (synop, ships, buoy, metar, sonde_sfc, tamdar_sfc)
                        if(allocated(platform(obs_type)%surf)) then
                            allocate(local_omb(5,nobs), &
                                     local_qc (5,nobs))

                            local_omb = platform(obs_type)%surf%omb(:,:,0)
                            local_qc  = platform(obs_type)%surf%qc (:,:,0)

                            sendcount = 5 * nobs
                        else
                            allocate(platform(obs_type)%surf)
                            allocate(platform(obs_type)%surf%h  (nobs), &
                                     platform(obs_type)%surf%pre(nobs), &
                                     platform(obs_type)%surf%obs  (5,nobs), &
                                     platform(obs_type)%surf%error(5,nobs), &
                                     platform(obs_type)%surf%omb  (5,nobs,0:nmember-1), &
                                     platform(obs_type)%surf%qc   (5,nobs,0:nmember-1))

                            sendcount = 0
                        end if

                        recvcount = 0
                        displs    = 5 * nobs * nmember
                        do k = 0, nmember-1
                            m            = k+root
                            recvcount(m) = 5 * nobs
                            displs(m)    = k * 5 * nobs
                        end do


                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%h,       nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%pre,     nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%obs,   5*nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%error, 5*nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)

                        call request_append
                        call mpi_iallgatherv(                local_omb, sendcount,         mpi_real, &
                                           platform(obs_type)%surf%omb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_iallgatherv(                local_qc, sendcount,         mpi_integer, &
                                           platform(obs_type)%surf%qc, recvcount, displs, mpi_integer, mpi_comm_world, req_ptr%req, mpierr)

                        if(allocated(local_omb)) deallocate(local_omb,  local_qc)

                    case (pilot, profiler, geoamv, qscat, polaramv)
                        allocate(i4_nlev(nobs))
                        if( allocated(platform(obs_type)%vert) ) then
                            do n = 1, nobs
                                i4_nlev(n) = platform(obs_type)%vert(n)%nlev
                            end do
                        else
                            allocate(platform(obs_type)%vert(nobs))
                        end if

                        call mpi_bcast(i4_nlev, nobs, mpi_integer,     root, mpi_comm_world, mpierr)

                        do n = 1, nobs
                            nlev = i4_nlev(n)
                            platform(obs_type)%vert(n)%nlev = nlev

                            if(allocated(platform(obs_type)%vert(n)%pre)) then
                                allocate(local_omb(2,nlev), &
                                         local_qc (2,nlev))

                                local_omb = platform(obs_type)%vert(n)%omb(:,:,0)
                                local_qc  = platform(obs_type)%vert(n)%qc (:,:,0)

                                sendcount = 2 * nlev
                            else
                                allocate(platform(obs_type)%vert(n)%h  (nlev), &
                                         platform(obs_type)%vert(n)%pre(nlev), &
                                         platform(obs_type)%vert(n)%obs  (2,nlev), &
                                         platform(obs_type)%vert(n)%error(2,nlev), &
                                         platform(obs_type)%vert(n)%omb  (2,nlev,0:nmember-1), &
                                         platform(obs_type)%vert(n)%qc   (2,nlev,0:nmember-1))

                                sendcount = 0
                            end if

                            recvcount = 0
                            displs    = 2 * nlev * nmember
                            do k = 0, nmember-1
                                m            = k+root
                                recvcount(m) = 2 * nlev
                                displs(m)    = k * 2 * nlev
                            end do

                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%h,       nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%pre,     nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%obs,   2*nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)

                            call request_append
                            call mpi_iallgatherv(                   local_omb, sendcount,         mpi_real, &
                                               platform(obs_type)%vert(n)%omb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_iallgatherv(                   local_qc, sendcount,         mpi_integer, &
                                               platform(obs_type)%vert(n)%qc, recvcount, displs, mpi_integer, mpi_comm_world, req_ptr%req, mpierr)

                            if(allocated(local_omb)) deallocate(local_omb,  local_qc)

                        end do

                        deallocate(i4_nlev)
                    case (gpspw)
                        if( allocated(platform(obs_type)%surf) ) then
                            allocate(local_omb  (1,nobs), &
                                     local_qc   (1,nobs))

                            local_omb = platform(obs_type)%surf%omb(:,:,0)
                            local_qc  = platform(obs_type)%surf%qc (:,:,0)

                            sendcount = nobs
                        else
                            allocate(platform(obs_type)%surf)
                            allocate(platform(obs_type)%surf%h(nobs), &
                                     platform(obs_type)%surf%obs  (1,nobs), &
                                     platform(obs_type)%surf%error(1,nobs), &
                                     platform(obs_type)%surf%omb  (1,nobs,0:nmember-1), &
                                     platform(obs_type)%surf%qc   (1,nobs,0:nmember-1))

                             sendcount = 0
                        end if

                        recvcount = 0
                        displs    = nobs * nmember
                        do k = 0, nmember-1
                            m            = k+root
                            recvcount(m) = nobs
                            displs(m)    = k * nobs
                        end do

                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%h,     nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%obs,   nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_ibcast(platform(obs_type)%surf%error, nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)

                        call request_append
                        call mpi_iallgatherv(                local_omb, sendcount,         mpi_real, &
                                           platform(obs_type)%surf%omb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr)
                        call request_append
                        call mpi_iallgatherv(                local_qc, sendcount,         mpi_integer, &
                                           platform(obs_type)%surf%qc, recvcount, displs, mpi_integer, mpi_comm_world, req_ptr%req, mpierr)

                        if(allocated(local_omb)) deallocate(local_omb,  local_qc)

                    case (sound, tamdar, airep)
                        allocate(i4_nlev(nobs))
                        if( allocated(platform(obs_type)%vert) ) then
                            do n = 1, nobs
                                i4_nlev(n) = platform(obs_type)%vert(n)%nlev
                            end do
                        else
                            allocate(platform(obs_type)%vert(nobs))
                        end if

                        call mpi_bcast(i4_nlev, nobs, mpi_integer,     root, mpi_comm_world, mpierr)

                        do n = 1, nobs
                            nlev = i4_nlev(n)
                            platform(obs_type)%vert(n)%nlev = nlev

                            if(allocated(platform(obs_type)%vert(n)%pre)) then
                                allocate(local_omb  (4,nlev), &
                                         local_qc   (4,nlev))

                                local_omb = platform(obs_type)%vert(n)%omb(:,:,0)
                                local_qc  = platform(obs_type)%vert(n)%qc (:,:,0)

                                sendcount = 4 * nlev
                            else
                                allocate(platform(obs_type)%vert(n)%h  (nlev), &
                                         platform(obs_type)%vert(n)%pre(nlev), &
                                         platform(obs_type)%vert(n)%obs  (4,nlev), &
                                         platform(obs_type)%vert(n)%error(4,nlev), &
                                         platform(obs_type)%vert(n)%omb  (4,nlev,0:nmember-1), &
                                         platform(obs_type)%vert(n)%qc   (4,nlev,0:nmember-1))

                                sendcount = 0
                            end if

                            recvcount = 0
                            displs    = 4 * nlev * nmember
                            do k = 0, nmember-1
                                m            = k+root
                                recvcount(m) = 4 * nlev
                                displs(m)    = k * 4 * nlev
                            end do

                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%h,       nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%pre,     nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%obs,   4*nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%error, 4*nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)

                            call request_append
                            call mpi_iallgatherv(                   local_omb, sendcount,         mpi_real, &
                                               platform(obs_type)%vert(n)%omb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_iallgatherv(                   local_qc, sendcount,         mpi_integer, &
                                               platform(obs_type)%vert(n)%qc, recvcount, displs, mpi_integer, mpi_comm_world, req_ptr%req, mpierr)

                            if(allocated(local_omb)) deallocate(local_omb,  local_qc)

                        end do

                        deallocate(i4_nlev)

                    case (gpsref)
                        allocate(i4_nlev(nobs))
                        if( allocated(platform(obs_type)%vert) ) then
                            do n = 1, nobs
                                i4_nlev(n) = platform(obs_type)%vert(n)%nlev
                            end do
                        else
                            allocate(platform(obs_type)%vert(nobs))
                        end if

                        call mpi_bcast(i4_nlev, nobs, mpi_integer,     root, mpi_comm_world, mpierr)

                        do n = 1, nobs
                            nlev = i4_nlev(n)
                            platform(obs_type)%vert(n)%nlev = nlev
 
                            if(allocated(platform(obs_type)%vert(n)%h)) then
                                allocate(local_omb  (1,nlev), &
                                         local_qc   (1,nlev))

                                local_omb  = platform(obs_type)%vert(n)%omb(:,:,0)
                                local_qc   = platform(obs_type)%vert(n)%qc (:,:,0)

                                sendcount = nlev
                            else
                                allocate(platform(obs_type)%vert(n)%h(nlev), &
                                         platform(obs_type)%vert(n)%obs  (1,nlev), &
                                         platform(obs_type)%vert(n)%error(1,nlev), &
                                         platform(obs_type)%vert(n)%omb  (1,nlev,0:nmember-1), &
                                         platform(obs_type)%vert(n)%qc   (1,nlev,0:nmember-1))

                                 sendcount = 0
                            end if

                            recvcount = 0
                            displs    = nlev * nmember
                            do k = 0, nmember-1
                                m            = k+root
                                recvcount(m) = nlev
                                displs(m)    = k * nlev
                            end do

                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%h,     nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%obs,   nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_ibcast(platform(obs_type)%vert(n)%error, nlev, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)

                            call request_append
                            call mpi_iallgatherv(                    local_omb, sendcount,         mpi_real, &
                                                platform(obs_type)%vert(n)%omb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr)
                            call request_append
                            call mpi_iallgatherv(                    local_qc, sendcount,         mpi_integer, &
                                                platform(obs_type)%vert(n)%qc, recvcount, displs, mpi_integer, mpi_comm_world, req_ptr%req, mpierr)

                            if(allocated(local_omb)) deallocate(local_omb,  local_qc)

                        end do

                        deallocate(i4_nlev)

                    end select
                end if
            end do
        end associate

        deallocate(i4_nobs)
    end subroutine gts_distribute

    subroutine write_gts_omboma(self, id)
        implicit none
        class(wrfda_gts), intent(in) :: self
        integer,          intent(in) :: id
        character(len=3)             :: proc

        !local
        integer                      :: ierr, ifile
        integer                      :: nobs, nlev, obs_type
        integer                      :: n, k

        if(myid == id) then

        associate( platform => self % platform )
            do ifile = 0, nmember-1
                write(proc, '(i3.3)') ifile+1
                open(40, file='gts_out_'//proc, iostat=ierr)

                do obs_type = 1, num_gts_indexes
                    nobs = platform(obs_type)%nobs
                    if(nobs == 0) cycle
                    write(40,'(a20,i8)') trim(gts_names(obs_type)), nobs

                    select case (obs_type)
                    case (synop, ships, buoy, metar, sonde_sfc, tamdar_sfc)
                        do n = 1, nobs 
                            write(40, '(2i8)') 1, 1

                            write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                n, 1, platform(obs_type)%id(n),  &
                                platform(obs_type)%lat(n),  &
                                platform(obs_type)%lon(n),  &   
                                platform(obs_type)%surf%pre(n),                                             &
                                platform(obs_type)%surf%obs(1,n),       platform(obs_type)%surf%omb  (1,n,ifile),       &
                                platform(obs_type)%surf%qc (1,n,ifile), platform(obs_type)%surf%error(1,n),       &
                                0.0,                      &
                                platform(obs_type)%surf%obs(2,n),       platform(obs_type)%surf%omb  (2,n,ifile),       &
                                platform(obs_type)%surf%qc (2,n,ifile), platform(obs_type)%surf%error(2,n),       &
                                0.0,                      &
                                platform(obs_type)%surf%obs(3,n),       platform(obs_type)%surf%omb  (3,n,ifile),       &
                                platform(obs_type)%surf%qc (3,n,ifile), platform(obs_type)%surf%error(3,n),       &
                                0.0,                      &
                                platform(obs_type)%surf%obs(4,n),       platform(obs_type)%surf%omb  (4,n,ifile),       &
                                platform(obs_type)%surf%qc (4,n,ifile), platform(obs_type)%surf%error(4,n),       &
                                0.0,                      &
                                platform(obs_type)%surf%obs(5,n),       platform(obs_type)%surf%omb  (5,n,ifile),       &
                                platform(obs_type)%surf%qc (5,n,ifile), platform(obs_type)%surf%error(5,n),       &
                                0.0
                        end do  !nobs
                    case (pilot, profiler, geoamv, qscat, polaramv)
                        do n = 1, nobs 
                            nlev = platform(obs_type)%vert(n)%nlev

                            write(40,'(2i8)') nlev, 1

                            do k = 1, nlev
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    n, k, platform(obs_type)%id(n), &   
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%vert(n)%pre(k),                                          &
                                    platform(obs_type)%vert(n)%obs(1,k),        platform(obs_type)%vert(n)%omb  (1,k,ifile),       &
                                    platform(obs_type)%vert(n)%qc (1,k,ifile),  platform(obs_type)%vert(n)%error(1,k),       &
                                    0.0,                                     &
                                    platform(obs_type)%vert(n)%obs(2,k),        platform(obs_type)%vert(n)%omb  (2,k,ifile),       &
                                    platform(obs_type)%vert(n)%qc (2,k,ifile),  platform(obs_type)%vert(n)%error(2,k),       &
                                    0.0

                            end do !nlev
                        end do  !nobs
                    case (gpspw)
                        do n = 1, nobs 
                            write(40,'(2i8)') 1, 1

                            write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                n, 1, platform(obs_type)%id(n),  &   
                                platform(obs_type)%lat(n),  &
                                platform(obs_type)%lon(n),  &   
                                platform(obs_type)%surf%h  (n),                                                            &
                                platform(obs_type)%surf%obs(1,n),        platform(obs_type)%surf%omb  (1,n,ifile),     &
                                platform(obs_type)%surf%qc (1,n,ifile),  platform(obs_type)%surf%error(1,n),     &
                                0.0

                        end do  !nobs
                    case (sound, tamdar, airep)
                        do n = 1, nobs 
                            nlev = platform(obs_type)%vert(n)%nlev

                            write(40,'(2i8)') nlev, 1

                            do k = 1, nlev
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    n, k, platform(obs_type)%id(n), &   
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%vert(n)%pre(k),                                          &
                                    platform(obs_type)%vert(n)%obs(1,k),        platform(obs_type)%vert(n)%omb  (1,k,ifile),       &
                                    platform(obs_type)%vert(n)%qc (1,k,ifile),  platform(obs_type)%vert(n)%error(1,k),       &
                                    0.0,                                     &
                                    platform(obs_type)%vert(n)%obs(2,k),        platform(obs_type)%vert(n)%omb  (2,k,ifile),       &
                                    platform(obs_type)%vert(n)%qc (2,k,ifile),  platform(obs_type)%vert(n)%error(2,k),       &
                                    0.0,                                     &
                                    platform(obs_type)%vert(n)%obs(3,k),        platform(obs_type)%vert(n)%omb  (3,k,ifile),       &
                                    platform(obs_type)%vert(n)%qc (3,k,ifile),  platform(obs_type)%vert(n)%error(3,k),       &
                                    0.0,                                     &
                                    platform(obs_type)%vert(n)%obs(4,k),        platform(obs_type)%vert(n)%omb  (4,k,ifile),       &
                                    platform(obs_type)%vert(n)%qc (4,k,ifile),  platform(obs_type)%vert(n)%error(4,k),       &
                                    0.0

                            end do   !nlev
                        end do  !nobs
                    case (gpsref)
                        do n = 1, nobs 
                            nlev = platform(obs_type)%vert(n)%nlev

                            write(40,'(2i8)') nlev, 1

                            do k = 1, nlev
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr)    &
                                    n, k, platform(obs_type)%id(n), &   
                                    platform(obs_type)%lat(n),  &
                                    platform(obs_type)%lon(n),  &   
                                    platform(obs_type)%vert(n)%h  (k),                                              &
                                    platform(obs_type)%vert(n)%obs(1,k),        platform(obs_type)%vert(n)%omb  (1,k,ifile),     &
                                    platform(obs_type)%vert(n)%qc (1,k,ifile),  platform(obs_type)%vert(n)%error(1,k),     &
                                    0.0
                            end do   !nlev

                        end do  !nobs
                    end select
                end do

                close(40)
            end do
        end associate

        end if

        call mpi_barrier(mpi_comm_world, n)
        call mpi_abort(mpi_comm_world, k, n)
    end subroutine write_gts_omboma

    subroutine read_alt_info(gtsalt, filename)
        type(alt_structure), dimension(:), intent(in out) :: gtsalt
        character(len=*),                  intent(in)     :: filename
        real                :: lat, lon
        character(len=10)   :: fmt_name
        character(len=19)   :: date
        character(len=40)   :: source
        character(len=45)   :: info_fmt
        character(len=29)   :: srfc_fmt
        character(len=60)   :: each_fmt


        integer             :: gpszd, amdar
        real                :: missing, alt, dummyr
        integer             :: total, level, fm, obs_type
        integer             :: dummyi, iost, k, n
        character(len=40)   :: id
        character(len=160)  :: dummyc
        integer, dimension(num_gts_indexes) :: reports

        reports = 0

        open(11,file=filename,status='old')
        read(11,fmt='(A6,1X,I7,2X,A6,1X,F8.0)') dummyc, total, dummyc, missing
        read(11,fmt='(6(A6,1X,I7,2X))') dummyc, gtsalt( synop)%nobs, &
                                        dummyc, gtsalt( metar)%nobs, &
                                        dummyc, gtsalt( ships)%nobs, &
                                        dummyc, gtsalt( buoy) %nobs, &
                                        dummyc, gtsalt( bogus)%nobs, &
                                        dummyc, gtsalt( sound)%nobs
        read(11,fmt='(6(A6,1X,I7,2X))') dummyc, amdar,             &
                                        dummyc, gtsalt( airep)%nobs, &
                                        dummyc, gtsalt(tamdar)%nobs,&
                                        dummyc, gtsalt( pilot)%nobs, &
                                        dummyc, gtsalt( satem)%nobs, &
                                        dummyc, dummyi 
        read(11,fmt='(6(A6,1X,I7,2X))') dummyc, gtsalt( gpspw)%nobs, &
                                        dummyc, gpszd,              &
                                        dummyc, gtsalt(gpsref)%nobs,&
                                        dummyc, gtsalt(gpseph)%nobs,&
                                        dummyc, gtsalt( ssmt1)%nobs, &
                                        dummyc, gtsalt( ssmt2)%nobs
        read(11,fmt='(6(A6,1X,I7,2X))') dummyc, dummyi,  &
                                        dummyc, gtsalt( qscat)%nobs, &
                                        dummyc, gtsalt(profiler)%nobs, &
                                        dummyc, gtsalt( airsr)%nobs

        gtsalt( gpspw)%nobs = gtsalt( gpspw)%nobs + gpszd
        gtsalt(tamdar)%nobs = gtsalt(tamdar)%nobs + amdar

        do obs_type = 1, num_gts_indexes
            if(gtsalt(obs_type)%nobs > 0) then
                n = gtsalt(obs_type)%nobs
                allocate(gtsalt(obs_type)%id  (n), &
                         gtsalt(obs_type)%info(n))
            end if
        end do

        do 
            read(11,'(A)') dummyc
            if(dummyc(1:6) == 'EACH  ') exit
        end do

        read(11,'(A,1X,A)') &
            fmt_name, info_fmt, &
            fmt_name, srfc_fmt, &
            fmt_name, each_fmt

        read(11,'(A)')

        do 
            read(11, fmt=info_fmt, iostat=iost) &
                dummyc, date, source, level, lat, lon, alt, id
            if(iost > 0) stop "Problem"
            if(iost < 0) exit


            if (dummyc(6:6) == ' ') then
                read(dummyc(4:5), '(I2)') fm
            else
                read(dummyc(4:6), '(I3)') fm
            end if


            select case (fm)
            case (12)
                reports(synop) = reports(synop)+1
                n              = reports(synop)
                gtsalt(synop)%id(n) = trim(id) 
                allocate(gtsalt(synop)%info(n)%alt(1))

                read(11,'(A)')
                read(11,each_fmt) &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    gtsalt(synop)%info(n)%alt(1)
            case (13, 17)
                reports(ships) = reports(ships)+1
                n              = reports(ships)
                gtsalt(ships)%id(n) = trim(id) 
                allocate(gtsalt(ships)%info(n)%alt(1))

                read(11,'(A)')
                read(11,each_fmt) &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    gtsalt(ships)%info(n)%alt(1)
            case (15:16)
                reports(metar) = reports(metar)+1
                n              = reports(metar)
                gtsalt(metar)%id(n) = trim(id) 
                allocate(gtsalt(metar)%info(n)%alt(1))

                read(11,'(A)')
                read(11,each_fmt) &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    gtsalt(metar)%info(n)%alt(1)
            case (32:34)
                reports(pilot) = reports(pilot)+1
                n              = reports(pilot)
                gtsalt(pilot)%id(n) = trim(id) 
                allocate(gtsalt(pilot)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(pilot)%info(n)%alt(k)
                end do
            case (35:38)
                reports(sound) = reports(sound)+1
                n              = reports(sound)
                gtsalt(sound)%id(n) = trim(id) 
                allocate(gtsalt(sound)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(sound)%info(n)%alt(k)
                end do
            case (101)
                reports(tamdar) = reports(tamdar)+1
                n               = reports(tamdar)
                gtsalt(tamdar)%id(n) = trim(id) 
                allocate(gtsalt(tamdar)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(tamdar)%info(n)%alt(k)
                end do
            case (161)
                reports(mtgirs) = reports(mtgirs)+1
                n               = reports(mtgirs)
                gtsalt(mtgirs)%id(n) = trim(id) 
                allocate(gtsalt(mtgirs)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(mtgirs)%info(n)%alt(k)
                end do
            case (86)
                reports(satem) = reports(satem)+1
                n              = reports(satem)
                gtsalt(satem)%id(n) = trim(id) 
                allocate(gtsalt(satem)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(satem)%info(n)%alt(k)
                end do
            case (42, 96:97)
                reports(airep) = reports(airep)+1
                n              = reports(airep)
                gtsalt(airep)%id(n) = trim(id) 
                allocate(gtsalt(airep)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(airep)%info(n)%alt(k)
                end do
            case (111, 114)
                reports(gpspw) = reports(gpspw)+1
                n              = reports(gpspw)
                gtsalt(gpspw)%id(n) = trim(id) 

                allocate(gtsalt(gpspw)%info(n)%alt(1))
                gtsalt(gpspw)%info(n)%alt(1) = alt

                read(11,'(A)')
            case (116)
                reports(gpsref) = reports(gpsref)+1
                n               = reports(gpsref)
                gtsalt(gpsref)%id(n) = trim(id) 
                allocate(gtsalt(gpsref)%info(n)%alt(1))

                read(11,'(A)')
                read(11,each_fmt) &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    gtsalt(gpsref)%info(n)%alt(1)
            case (121)
                reports(ssmt1) = reports(ssmt1)+1
                n              = reports(ssmt1)
                gtsalt(ssmt1)%id(n) = trim(id) 
                allocate(gtsalt(ssmt1)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(ssmt1)%info(n)%alt(k)
                end do
            case (122)
                reports(ssmt2) = reports(ssmt2)+1
                n              = reports(ssmt2)
                gtsalt(ssmt2)%id(n) = trim(id) 
                allocate(gtsalt(ssmt2)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(ssmt2)%info(n)%alt(k)
                end do
            case (281)
                reports(qscat) = reports(qscat)+1
                n              = reports(qscat)
                gtsalt(qscat)%id(n) = trim(id) 
                allocate(gtsalt(qscat)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(qscat)%info(n)%alt(k)
                end do
            case (132)
                reports(profiler) = reports(profiler)+1
                n                 = reports(profiler)
                gtsalt(profiler)%id(n) = trim(id) 
                allocate(gtsalt(profiler)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(profiler)%info(n)%alt(k)
                end do
            case (135)
                reports(bogus) = reports(bogus)+1
                n              = reports(bogus)
                gtsalt(bogus)%id(n) = trim(id) 
                allocate(gtsalt(bogus)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(bogus)%info(n)%alt(k)
                end do
            case (18, 19)
                reports(buoy) = reports(buoy)+1
                n             = reports(buoy)
                gtsalt(buoy)%id(n) = trim(id) 
                allocate(gtsalt(buoy)%info(n)%alt(1))

                read(11,'(A)')
                read(11,each_fmt) &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    dummyr, dummyi, dummyr, &
                    gtsalt(buoy)%info(n)%alt(1)
            case (133)
                reports(airsr) = reports(airsr)+1
                n              = reports(airsr)
                gtsalt(airsr)%id(n) = trim(id)
                allocate(gtsalt(airsr)%info(n)%alt(level))

                read(11,'(A)')
                do k = 1, level
                    read(11,each_fmt) &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        dummyr, dummyi, dummyr, &
                        gtsalt(airsr)%info(n)%alt(k)
                end do
            end select
        end do
        close(11)

    end subroutine read_alt_info

    real function get_alt(gtsalt, id, level)
        implicit none
        type(alt_structure), intent(in) :: gtsalt
        character(len=*),    intent(in) :: id
        integer,             intent(in) :: level

        !local
        integer                         :: i

        do i = 1, gtsalt%nobs
            if(trim(gtsalt%id(i)) == trim(id)) then
                get_alt = gtsalt%info(i)%alt(level)
                return
            end if
        end do

        stop "ID not found!!"
    end function get_alt


end module gts_omboma
