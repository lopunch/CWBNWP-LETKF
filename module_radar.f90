module simulated_radar

    use config, only : nmember
    use param
    use mpi_util

    implicit none

    private
    public  :: radar_structure, cwb_radar

    type, extends(obs_structure)           :: radar_structure
        real, dimension(:),   allocatable  :: obs, alt
        real, dimension(:,:), allocatable  :: omb
    end type radar_structure

    !========================================================================================

    type cwb_radar
        type(radar_structure), dimension(:), allocatable :: radarobs
    contains
        procedure, public, pass(self) :: read_data  => read_radar
        procedure, public, pass(self) :: write_data => write_radar
        procedure, public, pass(self) :: distribute => radar_distribute
    end type cwb_radar

contains

    subroutine read_radar(self, filename, varname)
        implicit none
        class(cwb_radar), intent(in out)  :: self
        character(len=*), intent(in)      :: filename, varname

        !local
        integer              :: ierr
        integer              :: n, nobs, obs_type

        logical, parameter   :: write_test = .false.


        open(50, file=filename, iostat=ierr)
        if(ierr /= 0) then
            print '(A,i3.3,A)', "CPU: ", myid, " open "//varname//"_letkf  error"
            stop
        end if

        read(50, '(i10)', iostat=ierr) nobs
        if(ierr < 0) then
            close(50)
            return
        end if
        if(ierr > 0) then
            print '(A,i3.3,A)', "CPU: ", myid, "read "//varname//"_letkf nobs error"
            stop
        end if

        if(nobs <= 0) return

        if(write_test) then
            open(30, file=varname//'_letkf_out')
            write(30, '(i10)') nobs
        end if

        select case(varname)
        case ('VR')
            obs_type = vr
        case ('MR')
            obs_type = dbz
        case ('MD')
            obs_type = zdr
        case ('MK')
            obs_type = kdp
        end select

        associate( radarobs => self % radarobs )
            radarobs(obs_type)%nobs  = nobs
            allocate(radarobs(obs_type)%obs(nobs),  &
                     radarobs(obs_type)%lat(nobs),  &
                     radarobs(obs_type)%lon(nobs),  &
                     radarobs(obs_type)%alt(nobs),  &
                     radarobs(obs_type)%omb(nobs,0:nmember-1))

            do n = 1, nobs
                read(50, '(5(f10.4,1x))', iostat=ierr)    &
                             radarobs(obs_type)%obs(n),   &
                             radarobs(obs_type)%omb(n,0), &    !read background value, not omb yet
                             radarobs(obs_type)%lon(n),   &
                             radarobs(obs_type)%lat(n),   &
                             radarobs(obs_type)%alt(n)
                if(ierr < 0) exit
                if(ierr > 0) then
                    print *, "read "//varname//"_letkf data error"
                    stop
                end if

               !radarobs(obs_type)%xy(:,n) = latlon_to_dist( [lat, lon] )

                if(write_test) then
                    write(30,'(5(f10.4,1x))')           &
                             radarobs(obs_type)%obs(n), &
                             radarobs(obs_type)%omb(n,0), &
                             radarobs(obs_type)%lon(n),   &
                             radarobs(obs_type)%lat(n),   &
                             radarobs(obs_type)%alt(n)
                end if
            end do 
        end associate

        close(50)
        if(write_test) close(30)

    end subroutine read_radar

    subroutine radar_distribute(self, root)
        implicit none
        class(cwb_radar), intent(in out)        :: self
        integer,          intent(in)            :: root

        !local
        integer                                 :: obs_type, nobs
        integer                                 :: k, m

        !local mpi 
        integer,  dimension(0:nproc-1)          :: recvcount, displs
        integer                                 :: mpierr, sendcount
        integer,  dimension(:),     allocatable :: i4_nobs
        real,     dimension(:),     allocatable :: local_omb


        !broadcast number of obs of each obs type
        allocate(i4_nobs(num_radar_indexes))
        if(myid == root) then
            do obs_type = 1, num_radar_indexes
                i4_nobs(obs_type) = self%radarobs(obs_type)%nobs
            end do
        end if

        call mpi_bcast(i4_nobs, num_radar_indexes, mpi_integer, root, mpi_comm_world, mpierr)

        !start exchange obs info
        associate( radarobs => self%radarobs )
            do obs_type = 1, num_radar_indexes
                nobs = i4_nobs(obs_type)
                radarobs(obs_type)%nobs = nobs
                if(nobs > 0) then
                    if(allocated(radarobs(obs_type)%obs)) then
                        allocate(local_omb(nobs))
                        local_omb = radarobs(obs_type)%omb(:,0)
                        sendcount = nobs
                    else
                        allocate(radarobs(obs_type)%obs(nobs),  &
                                 radarobs(obs_type)%lat(nobs),  &
                                 radarobs(obs_type)%lon(nobs),  &
                                 radarobs(obs_type)%alt(nobs),  &
                                 radarobs(obs_type)%omb(nobs,0:nmember-1))
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
                    call mpi_ibcast(radarobs(obs_type)%obs, nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(radarobs(obs_type)%lat, nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(radarobs(obs_type)%lon, nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_ibcast(radarobs(obs_type)%alt, nobs, mpi_real, root, mpi_comm_world, req_ptr%req, mpierr)
                    call request_append
                    call mpi_iallgatherv(             local_omb, sendcount,         mpi_real, &
                                         radarobs(obs_type)%omb, recvcount, displs, mpi_real, mpi_comm_world, req_ptr%req, mpierr )

                    if(allocated(local_omb)) deallocate(local_omb)
                end if
            end do
        end associate

        deallocate(i4_nobs)
    end subroutine radar_distribute

    subroutine write_radar(self, id)
        implicit none
        class(cwb_radar), intent(in) :: self
        integer,          intent(in) :: id
        character(len=3)             :: proc
        character(len=2), parameter  :: obs_names(num_radar_indexes) = &
            (/ "MR", "VR", "ZD", "KP" /)

        !local
        integer                      :: ierr, ifile
        integer                      :: nobs, obs_type
        integer                      :: n

        if(myid == id) then

        associate( radarobs => self % radarobs )
            do ifile = 0, nmember-1
                write(proc, '(i3.3)') ifile+1
                do obs_type = 1, num_radar_indexes
                    nobs = radarobs(obs_type)%nobs
                    if(nobs > 0) then
                        open(40, file=obs_names(obs_type)//'_out_'//proc, iostat=ierr)

                        write(40, '(i10)', iostat=ierr) nobs
                        do n = 1, nobs
                            write(40,'(5(f10.4,1x))')           &
                                     radarobs(obs_type)%obs(n), &
                                     radarobs(obs_type)%omb(n,ifile), &
                                     radarobs(obs_type)%lon(n),   &
                                     radarobs(obs_type)%lat(n),   &
                                     radarobs(obs_type)%alt(n)
                        end do

                        close(40)
                    end if
                end do
            end do
        end associate

        end if

    end subroutine write_radar

!   function handle_error(ierr, message) result(nodata)
!       implicit none
!       integer,          intent(in) :: ierr
!       character(len=*), intent(in) :: message
!       logical                      :: nodata

!       if(ierr < 0) then
!           nodata = .true.
!       else if(ierr == 0) then
!           nodata = .false.
!       else 
!           print *, message
!           stop
!       end if
!   end function handle_error

end module simulated_radar
