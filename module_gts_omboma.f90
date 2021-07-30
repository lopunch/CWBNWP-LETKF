module gts_omboma

    use config, only : nmember
    use param
    use mpi_util

    implicit none

    private
    public  :: gts_structure, wrfda_gts

    type, extends(obs_structure) :: gts_structure
        character(len=5), dimension(:),     allocatable :: id
        real,             dimension(:),     allocatable :: pre

        !===============================================
        real,             dimension(:,:),   allocatable :: obs       ![nvar, nobs]
        real,             dimension(:,:),   allocatable :: error     ![nvar, nobs]
        real,             dimension(:,:,:), allocatable :: hdxb      ![nvar, nobs, nmember]      
        integer,          dimension(:,:,:), allocatable :: qc        ![nvar, nobs, nmember]
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
        type, extends(obs_structure) :: vert_structure
            character(len=5)                        :: id
            real,     dimension(:), allocatable     :: pre

            !===============================================
            real,     dimension(:,:),   allocatable :: obs       ![nvar, nobs]
            real,     dimension(:,:),   allocatable :: error     ![nvar, nobs]
            real,     dimension(:,:,:), allocatable :: hdxb      ![nvar, nobs, nmember]      
            integer,  dimension(:,:,:), allocatable :: qc        ![nvar, nobs, nmember]
        end type vert_structure

        type(vert_structure), dimension(:), allocatable   :: vert
        type(alt_structure),  dimension(num_gts_indexes)  :: gtsalt

        integer                                           :: ierr
        integer                                           :: nobs, nlev, obs_type
        character(len=20)                                 :: iv_type, obs_name

        integer                                           :: n, k, total
        integer                                           :: kk, l, nreq
        real, dimension(5)                                :: oma
        logical, parameter                                :: write_gts = .false.
        
        !read altitude info
        call read_alt_info(gtsalt, obascii)

        open(30, file=filename, iostat=ierr)
        if(ierr /= 0) stop "open gts_omboma error"

        if(write_gts) open(40, file='gts_out_'//filename(20:22), iostat=ierr)
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

                    associate( platform => self % platform(obs_type) )
                        platform % nobs = nobs 
                        allocate(platform % id   (nobs), &
                                 platform % lat  (nobs), &
                                 platform % lon  (nobs), &
                                 platform % alt  (nobs), &
                                 platform % pre  (nobs), &
                                 platform % obs  (5,nobs),         &
                                 platform % error(5,nobs),         &
                                 platform % hdxb (5,nobs,0:nmember-1),         &
                                 platform % qc   (5,nobs,0:nmember-1))

                        do n = 1, nobs 
                            read(30, '(2i8)') nlev, nreq
                            if(write_gts) write(40, '(2i8)') nlev, nreq

                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                platform % id (n),     platform % lat(n),  &
                                platform % lon(n),     platform % pre(n),  &
                                platform % obs(1,n),   platform % hdxb (1,n,0),       &
                                platform % qc (1,n,0), platform % error(1,n), oma(1),         &
                                platform % obs(2,n),   platform % hdxb (2,n,0),       &
                                platform % qc (2,n,0), platform % error(2,n), oma(2),         &
                                platform % obs(3,n),   platform % hdxb (3,n,0),       &
                                platform % qc (3,n,0), platform % error(3,n), oma(3),         &
                                platform % obs(4,n),   platform % hdxb (4,n,0),       &
                                platform % qc (4,n,0), platform % error(4,n), oma(4),         &
                                platform % obs(5,n),   platform % hdxb (5,n,0),       &
                                platform % qc (5,n,0), platform % error(5,n), oma(5)

                            platform % alt(n) = get_alt(gtsalt(obs_type), platform % id(n), 1)

                            if(write_gts) &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                    platform % id (n),     platform % lat(n),  &
                                    platform % lon(n),     platform % pre(n),  &
                                    platform % obs(1,n),   platform % hdxb (1,n,0),       &
                                    platform % qc (1,n,0), platform % error(1,n), oma(1),         &
                                    platform % obs(2,n),   platform % hdxb (2,n,0),       &
                                    platform % qc (2,n,0), platform % error(2,n), oma(2),         &
                                    platform % obs(3,n),   platform % hdxb (3,n,0),       &
                                    platform % qc (3,n,0), platform % error(3,n), oma(3),         &
                                    platform % obs(4,n),   platform % hdxb (4,n,0),       &
                                    platform % qc (4,n,0), platform % error(4,n), oma(4),         &
                                    platform % obs(5,n),   platform % hdxb (5,n,0),       &
                                    platform % qc (5,n,0), platform % error(5,n), oma(5)
                        end do  !nobs

                        !hdxb = obs - omb
                        platform % hdxb(:,:,0) = platform % obs - platform % hdxb(:,:,0)
                    end associate
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

                    total = 0
                    allocate(vert(nobs))
                    
                    do n = 1, nobs 
                        read(30,'(2i8)') nlev, nreq
                        vert(n) % nobs = nlev

                        if(write_gts) write(40,'(2i8)') nlev, nreq

                        allocate(vert(n) % lat  (nlev), &
                                 vert(n) % lon  (nlev), &  
                                 vert(n) % alt  (nlev), &
                                 vert(n) % pre  (nlev), &
                                 vert(n) % obs  (2,nlev), &
                                 vert(n) % error(2,nlev), &
                                 vert(n) % hdxb (2,nlev,0:nmember-1), &
                                 vert(n) % qc   (2,nlev,0:nmember-1))

                        do k = 1, nlev
                            total = total + 1

                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                vert(n) % id,         vert(n) % lat(k),    &
                                vert(n) % lon(k),     vert(n) % pre(k),    &
                                vert(n) % obs(1,k),   vert(n) % hdxb (1,k,0),       &
                                vert(n) % qc (1,k,0), vert(n) % error(1,k), oma(1),                  &
                                vert(n) % obs(2,k),   vert(n) % hdxb (2,k,0),       &
                                vert(n) % qc (2,k,0), vert(n) % error(2,k), oma(2)

                            vert(n) % alt(k) = get_alt(gtsalt(obs_type), vert(n) % id, k)
                           !if(k == 1) platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                            if(write_gts)  &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                    vert(n) % id,         vert(n) % lat(k),    &
                                    vert(n) % lon(k),     vert(n) % pre(k),    &
                                    vert(n) % obs(1,k),   vert(n) % hdxb (1,k,0),       &
                                    vert(n) % qc (1,k,0), vert(n) % error(1,k), oma(1),                  &
                                    vert(n) % obs(2,k),   vert(n) % hdxb (2,k,0),       &
                                    vert(n) % qc (2,k,0), vert(n) % error(2,k), oma(2)

                        end do   !nlev

                        !hdxb = obs - omb
                        vert(n) % hdxb(:,:,0) = vert(n) % obs - vert(n) % hdxb(:,:,0)
                    end do  !nobs

                    !=======================================================================

                    !put data from vert to platform
                    associate( platform => self % platform(obs_type) )
                        platform % nobs = total
                        allocate(platform % id   (total), &
                                 platform % lat  (total), &
                                 platform % lon  (total), &
                                 platform % alt  (total), &
                                 platform % pre  (total), &
                                 platform % obs  (2,total),         &
                                 platform % error(2,total),         &
                                 platform % hdxb (2,total,0:nmember-1),         &
                                 platform % qc   (2,total,0:nmember-1))

                        total = 0

                        do n = 1, nobs 
                        do k = 1, vert(n) % nobs
                            total = total + 1
                            platform % id   (total)     = vert(n) % id
                            platform % lat  (total)     = vert(n) % lat(k)
                            platform % lon  (total)     = vert(n) % lon(k)
                            platform % alt  (total)     = vert(n) % alt(k)
                            platform % pre  (total)     = vert(n) % pre(k)
                            platform % obs  (:,total)   = vert(n) % obs  (:,k)
                            platform % error(:,total)   = vert(n) % error(:,k)
                            platform % hdxb (:,total,0) = vert(n) % hdxb (:,k,0)
                            platform % qc   (:,total,0) = vert(n) % qc   (:,k,0)
                        end do
                        end do
                    end associate

                    deallocate(vert)

                end if
            case ('gpspw')
                if(nobs > 0) then
                    obs_type = gpspw

                    associate( platform => self % platform(obs_type) )
                        platform % nobs  = nobs 
                        allocate(platform % id (nobs), &
                                 platform % lat(nobs), &
                                 platform % lon(nobs), &
                                 platform % alt(nobs), &
                                 platform % obs  (1,nobs), &
                                 platform % error(1,nobs), &
                                 platform % hdxb (1,nobs,0:nmember-1), &
                                 platform % qc   (1,nobs,0:nmember-1))
                        
                        do n = 1, nobs 
                            read(30,'(2i8)') nlev, nreq
                            if(write_gts) write(40,'(2i8)') nlev, nreq

                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                platform % id (n),     platform % lat(n),  &
                                platform % lon(n),     platform % alt(n),  &
                                platform % obs(1,n),   platform % hdxb (1,n,0),     &
                                platform % qc (1,n,0), platform % error(1,n), oma(1)

                           !platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                            if(write_gts) &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                    platform % id (n),     platform % lat(n),  &
                                    platform % lon(n),     platform % alt(n),  &
                                    platform % obs(1,n),   platform % hdxb (1,n,0),     &
                                    platform % qc (1,n,0), platform % error(1,n), oma(1)

                        end do  !nobs

                        !hdxb = obs - omb
                        platform % hdxb(:,:,0) = platform % obs - platform % hdxb(:,:,0)
                    end associate
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

                    total = 0
                    allocate(vert(nobs))
                    
                    do n = 1, nobs 
                        read(30,'(2i8)') nlev, nreq
                        vert(n) % nobs = nlev

                        if(write_gts) write(40,'(2i8)') nlev, nreq

                        allocate(vert(n) % lat  (nlev), &
                                 vert(n) % lon  (nlev), &  
                                 vert(n) % alt  (nlev), &
                                 vert(n) % pre  (nlev), &
                                 vert(n) % obs  (4,nlev), &
                                 vert(n) % error(4,nlev), &
                                 vert(n) % hdxb (4,nlev,0:nmember-1), &
                                 vert(n) % qc   (4,nlev,0:nmember-1))

                        do k = 1, nlev
                            total = total + 1
                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                vert(n) % id,         vert(n) % lat(k),    &
                                vert(n) % lon(k),     vert(n) % pre(k),    &
                                vert(n) % obs(1,k),   vert(n) % hdxb (1,k,0),       &
                                vert(n) % qc (1,k,0), vert(n) % error(1,k), oma(1),                  &
                                vert(n) % obs(2,k),   vert(n) % hdxb (2,k,0),       &
                                vert(n) % qc (2,k,0), vert(n) % error(2,k), oma(2),                  &
                                vert(n) % obs(3,k),   vert(n) % hdxb (3,k,0),       &
                                vert(n) % qc (3,k,0), vert(n) % error(3,k), oma(3),                  &
                                vert(n) % obs(4,k),   vert(n) % hdxb (4,k,0),       &
                                vert(n) % qc (4,k,0), vert(n) % error(4,k), oma(4)

                            vert(n) % alt(k) = get_alt(gtsalt(obs_type), vert(n) % id, k)
                           !if(k == 1) platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                            if(write_gts) &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                    vert(n) % id,         vert(n) % lat(k),    &
                                    vert(n) % lon(k),     vert(n) % pre(k),    &
                                    vert(n) % obs(1,k),   vert(n) % hdxb (1,k,0),       &
                                    vert(n) % qc (1,k,0), vert(n) % error(1,k), oma(1),                  &
                                    vert(n) % obs(2,k),   vert(n) % hdxb (2,k,0),       &
                                    vert(n) % qc (2,k,0), vert(n) % error(2,k), oma(2),                  &
                                    vert(n) % obs(3,k),   vert(n) % hdxb (3,k,0),       &
                                    vert(n) % qc (3,k,0), vert(n) % error(3,k), oma(3),                  &
                                    vert(n) % obs(4,k),   vert(n) % hdxb (4,k,0),       &
                                    vert(n) % qc (4,k,0), vert(n) % error(4,k), oma(4)
                        end do   !nlev

                        !hdxb = obs - omb
                        vert(n) % hdxb(:,:,0) = vert(n) % obs - vert(n) % hdxb(:,:,0)
                    end do  !nobs

                    !=======================================================================

                    !put data from vert to platform
                    associate( platform => self % platform(obs_type) )
                        platform % nobs = total
                        allocate(platform % id   (total), &
                                 platform % lat  (total), &
                                 platform % lon  (total), &
                                 platform % alt  (total), &
                                 platform % pre  (total), &
                                 platform % obs  (4,total),         &
                                 platform % error(4,total),         &
                                 platform % hdxb (4,total,0:nmember-1),         &
                                 platform % qc   (4,total,0:nmember-1))

                        total = 0

                        do n = 1, nobs 
                        do k = 1, vert(n) % nobs
                            total = total + 1
                            platform % id   (total)     = vert(n) % id
                            platform % lat  (total)     = vert(n) % lat(k)
                            platform % lon  (total)     = vert(n) % lon(k)
                            platform % alt  (total)     = vert(n) % alt(k)
                            platform % pre  (total)     = vert(n) % pre(k)
                            platform % obs  (:,total)   = vert(n) % obs  (:,k)
                            platform % error(:,total)   = vert(n) % error(:,k)
                            platform % hdxb (:,total,0) = vert(n) % hdxb (:,k,0)
                            platform % qc   (:,total,0) = vert(n) % qc   (:,k,0)
                        end do
                        end do

                        deallocate(vert)
                    end associate

                end if
            case ('gpsref')
                if(nobs > 0) then
                    obs_type = gpsref

                    total = 0
                    allocate(vert(nobs))

                    do n = 1, nobs 
                        read(30,'(2i8)') nlev, nreq
                        vert(n) % nobs = nlev
                    
                        if(write_gts) write(40,'(2i8)') nlev, nreq

                        allocate(vert(n) % lat  (nlev), &
                                 vert(n) % lon  (nlev), &  
                                 vert(n) % alt  (nlev), &
                                 vert(n) % obs  (1,nlev), &
                                 vert(n) % error(1,nlev), &
                                 vert(n) % hdxb (1,nlev,0:nmember-1), &
                                 vert(n) % qc   (1,nlev,0:nmember-1))

                        do k = 1, nlev
                            total = total + 1
                            read(30,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                vert(n) % id,         vert(n) % lat(k),    &
                                vert(n) % lon(k),     vert(n) % alt(k),    &
                                vert(n) % obs(1,k),   vert(n) % hdxb (1,k,0),     &
                                vert(n) % qc (1,k,0), vert(n) % error(1,k), oma(1)

                           !platform(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                            if(write_gts) &
                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) kk, l, &
                                    vert(n) % id,         vert(n) % lat(k),    &
                                    vert(n) % lon(k),     vert(n) % alt(k),    &
                                    vert(n) % obs(1,k),   vert(n) % hdxb (1,k,0),     &
                                    vert(n) % qc (1,k,0), vert(n) % error(1,k), oma(1)
                        end do

                        !hdxb = obs - omb
                        vert(n) % hdxb(:,:,0) = vert(n) % obs - vert(n) % hdxb(:,:,0)
                    end do  !nobs

                    !=======================================================================

                    !put data from vert to platform
                    associate( platform => self % platform(obs_type) )
                        platform % nobs = total
                        allocate(platform % id   (total), &
                                 platform % lat  (total), &
                                 platform % lon  (total), &
                                 platform % alt  (total), &
                                 platform % obs  (1,total),         &
                                 platform % error(1,total),         &
                                 platform % hdxb (1,total,0:nmember-1),         &
                                 platform % qc   (1,total,0:nmember-1))

                        total = 0

                        do n = 1, nobs 
                        do k = 1, vert(n) % nobs
                            total = total + 1
                            platform % id   (total)     = vert(n) % id
                            platform % lat  (total)     = vert(n) % lat(k)
                            platform % lon  (total)     = vert(n) % lon(k)
                            platform % alt  (total)     = vert(n) % alt(k)
                            platform % obs  (:,total)   = vert(n) % obs  (:,k)
                            platform % error(:,total)   = vert(n) % error(:,k)
                            platform % hdxb (:,total,0) = vert(n) % hdxb (:,k,0)
                            platform % qc   (:,total,0) = vert(n) % qc   (:,k,0)
                        end do
                        end do
                    end associate

                    deallocate(vert)

                end if
            end select
        end do report

        close(30)
        if(write_gts) close(40)
    end subroutine read_gts_omboma

    subroutine gts_distribute(self, root)

        implicit none
        class(wrfda_gts), intent(in out)               :: self
        integer,          intent(in)                   :: root

        !local
        integer                                        :: obs_type, nobs, nvar
        integer                                        :: k, m

        !local mpi 
        integer,  dimension(0:nproc-1)                 :: recvcount, displs
        integer                                        :: mpierr, sendcount
        integer,  dimension(:),     allocatable        :: i4_nobs

        real,     dimension(:,:),   allocatable        :: local_hdxb
        integer,  dimension(:,:),   allocatable        :: local_qc

        !broadcast number of obs of each obs type
        allocate(i4_nobs(num_gts_indexes))
        if(myid == root) then
            do obs_type = 1, num_gts_indexes
                i4_nobs(obs_type) = self % platform(obs_type) % nobs
            end do
        end if

        call mpi_bcast(i4_nobs, num_gts_indexes, mpi_integer, root, mpi_comm_world, mpierr)

        !start to exchange obs info
        do obs_type = 1, num_gts_indexes
            nobs = i4_nobs(obs_type)
            associate( platform => self % platform(obs_type) )
                platform % nobs = nobs
                if(nobs > 0) then

                    if(.not. allocated(platform % id)) then
                        allocate(platform % id (nobs), &
                                 platform % lat(nobs), &
                                 platform % lon(nobs), &
                                 platform % alt(nobs))
                    end if

                    call request_append
                    call mpi_ibcast(platform % id, 5*nobs, mpi_character, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(platform % lat,  nobs, mpi_real,      root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(platform % lon,  nobs, mpi_real,      root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(platform % alt,  nobs, mpi_real,      root, mpi_comm_world, req_ptr%req, mpierr)

                    select case (obs_type)
                    case (synop, ships, buoy, metar, sonde_sfc, tamdar_sfc)
                        nvar = 5
                    case (pilot, profiler, geoamv, qscat, polaramv)
                        nvar = 2
                    case (sound, tamdar, airep)
                        nvar = 4
                    case (gpspw, gpsref)
                        nvar = 1
                    end select

                    if(allocated(platform % hdxb)) then
                        allocate(local_hdxb(nvar,nobs), &
                                 local_qc  (nvar,nobs))

                        local_hdxb = platform % hdxb(:,:,0)
                        local_qc   = platform % qc  (:,:,0)

                        sendcount = nvar * nobs
                    else
                        allocate(platform % obs  (nvar,nobs), &
                                 platform % error(nvar,nobs), &
                                 platform % hdxb (nvar,nobs,0:nmember-1), &
                                 platform % qc   (nvar,nobs,0:nmember-1))

                        if(obs_type /= gpspw .and. obs_type /= gpsref) then
                            allocate(platform % pre(nobs))
                        end if

                        sendcount = 0
                    end if

                    recvcount = 0
                    displs    = nvar * nobs * nmember
                    do k = 0, nmember-1
                        m            = k+root
                        recvcount(m) = nvar * nobs
                        displs(m)    = k * recvcount(m)
                    end do

                    call request_append
                    call mpi_ibcast(platform % obs,   nvar*nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(platform % error, nvar*nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)

                    if(obs_type /= gpspw .and. obs_type /= gpsref) then
                        call request_append
                        call mpi_ibcast(platform % pre,    nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                    end if

                    call request_append
                    call mpi_iallgatherv(     local_hdxb, sendcount,         mpi_real, &
                                         platform % hdxb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_iallgatherv(      local_qc, sendcount,         mpi_integer, &
                                          platform % qc, recvcount, displs, mpi_integer, mpi_comm_world, req_ptr%req, mpierr)

                    if(allocated(local_hdxb)) deallocate(local_hdxb,  local_qc)
                end if
            end associate
        end do

        deallocate(i4_nobs)
    end subroutine gts_distribute

    subroutine write_gts_omboma(self, id)
        implicit none
        class(wrfda_gts), intent(in) :: self
        integer,          intent(in) :: id
        character(len=3)             :: proc

        !local
        integer                      :: ierr, ifile
        integer                      :: nobs, obs_type
        integer                      :: n, k

        if(myid == id) then
            do ifile = 0, nmember-1
                write(proc, '(i3.3)') ifile+1
                open(40, file='gts_out_'//proc, iostat=ierr)

                do obs_type = 1, num_gts_indexes
                    associate( platform => self % platform(obs_type) )
                        nobs = platform % nobs
                        if(nobs == 0) cycle
                        write(40,'(a20,i8)') trim(gts_names(obs_type)), nobs

                        select case (obs_type)
                        case (synop, ships, buoy, metar, sonde_sfc, tamdar_sfc)
                            do n = 1, nobs 
                                write(40, '(2i8)') 1, 1

                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) n, 1, &
                                    platform % id (n),         platform % lat(n),  &
                                    platform % lon(n),         platform % pre(n),  &
                                    platform % obs(1,n),       platform % hdxb (1,n,ifile),      &
                                    platform % qc (1,n,ifile), platform % error(1,n), 0.0,       &
                                    platform % obs(2,n),       platform % hdxb (2,n,ifile),      &
                                    platform % qc (2,n,ifile), platform % error(2,n), 0.0,       &
                                    platform % obs(3,n),       platform % hdxb (3,n,ifile),      &
                                    platform % qc (3,n,ifile), platform % error(3,n), 0.0,       &
                                    platform % obs(4,n),       platform % hdxb (4,n,ifile),      &
                                    platform % qc (4,n,ifile), platform % error(4,n), 0.0,       &
                                    platform % obs(5,n),       platform % hdxb (5,n,ifile),      &
                                    platform % qc (5,n,ifile), platform % error(5,n), 0.0
                            end do  !nobs
                        case (pilot, profiler, geoamv, qscat, polaramv)
                            do n = 1, nobs 
                                write(40, '(2i8)') 1, 1

                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) n, 1, &
                                    platform % id (n),         platform % lat(n),  &
                                    platform % lon(n),         platform % pre(n),  &
                                    platform % obs(1,n),       platform % hdxb (1,n,ifile),      &
                                    platform % qc (1,n,ifile), platform % error(1,n), 0.0,       &
                                    platform % obs(2,n),       platform % hdxb (2,n,ifile),      &
                                    platform % qc (2,n,ifile), platform % error(2,n), 0.0
                            end do  !nobs
                        case (gpspw, gpsref)
                            do n = 1, nobs 
                                write(40, '(2i8)') 1, 1

                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) n, 1, &
                                    platform % id (n),         platform % lat(n),  &
                                    platform % lon(n),         platform % alt(n),  &
                                    platform % obs(1,n),       platform % hdxb (1,n,ifile),      &
                                    platform % qc (1,n,ifile), platform % error(1,n), 0.0
                            end do  !nobs
                        case (sound, tamdar, airep)
                            do n = 1, nobs 
                                write(40, '(2i8)') 1, 1

                                write(40,'(2i8,a5,2f9.2,f17.7,5(2f17.7,i8,2f17.7))', iostat=ierr) n, 1, &
                                    platform % id (n),         platform % lat(n),  &
                                    platform % lon(n),         platform % pre(n),  &
                                    platform % obs(1,n),       platform % hdxb (1,n,ifile),      &
                                    platform % qc (1,n,ifile), platform % error(1,n), 0.0,       &
                                    platform % obs(2,n),       platform % hdxb (2,n,ifile),      &
                                    platform % qc (2,n,ifile), platform % error(2,n), 0.0,       &
                                    platform % obs(3,n),       platform % hdxb (3,n,ifile),      &
                                    platform % qc (3,n,ifile), platform % error(3,n), 0.0,       &
                                    platform % obs(4,n),       platform % hdxb (4,n,ifile),      &
                                    platform % qc (4,n,ifile), platform % error(4,n), 0.0
                            end do  !nobs
                        end select
                    end associate
                end do

                close(40)
            end do
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
