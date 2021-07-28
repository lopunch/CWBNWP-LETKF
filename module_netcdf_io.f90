module netcdf_io

    use netcdf

    implicit none

    private
    public  :: read_nc, write_nc, open_wrf_ncdf, create_wrf_ncdf


    type read_nc

        integer :: ncid

    contains

        procedure, public, pass(self)  :: get_variable_0d
        procedure, public, pass(self)  :: get_variable_1d
        procedure, public, pass(self)  :: get_variable_2d
        procedure, public, pass(self)  :: get_variable_3d
        procedure, public, pass(self)  :: get_attribute
        procedure, public, pass(self)  :: get_dimension
        procedure, public, pass(self)  :: close_reading

        generic,   public              :: get_variable   => get_variable_0d, &
                                                            get_variable_1d, &
                                                            get_variable_2d, &
                                                            get_variable_3d
    end type read_nc

    type write_nc

        integer :: ncid

    contains

        procedure, public, pass(self)  :: copy_header_from
        procedure, public, pass(self)  :: write_variable_1d
        procedure, public, pass(self)  :: write_variable_2d
        procedure, public, pass(self)  :: write_variable_3d
        procedure, public, pass(self)  :: write_variable_others
        procedure, public, pass(self)  :: close_writing

        generic,   public              :: write_variable => write_variable_1d, &
                                                            write_variable_2d, &
                                                            write_variable_3d, &
                                                            write_variable_others
    end type write_nc


    interface read_nc
        module procedure  ::   open_wrf_ncdf
    end interface read_nc

    interface write_nc
        module procedure  :: create_wrf_ncdf
    end interface write_nc

    type output_structure
        character(len=50)     :: name
        integer               :: xtype
        integer               :: ndims
        integer, dimension(4) :: dims = 1
        integer, dimension(:), allocatable :: dimids
    end type
    type(output_structure), dimension(:), allocatable  :: vars
    logical,                dimension(:), allocatable  :: I_am_not_write

contains

    type(read_nc)  function   open_wrf_ncdf(filename) result(self)
        implicit none
        character(len=*), intent(in) :: filename

        call check( nf90_open  (filename, nf90_nowrite, self % ncid) )
    end function open_wrf_ncdf

    type(write_nc) function create_wrf_ncdf(filename) result(self)
        implicit none
        character(len=*), intent(in) :: filename

        call check( nf90_create(filename, nf90_clobber, self % ncid) )
    end function create_wrf_ncdf

    subroutine close_reading(self)
        implicit none
        class(read_nc),  intent(in)  :: self

        call check( nf90_close(self % ncid) )
    end subroutine close_reading

    subroutine close_writing(self)
        implicit none
        class(write_nc), intent(in)  :: self

        call check( nf90_close(self % ncid) )
    end subroutine close_writing

    function get_dimension(self, varname) result(dsize)
        implicit none
        class(read_nc),     intent(in)  :: self
        character(len=*),   intent(in)  :: varname
        integer                         :: dsize

        integer                         :: dimid

        call check( nf90_inq_dimid(self % ncid, varname, dimid) )
        call check( nf90_inquire_dimension(self % ncid, dimid, len=dsize))
    end function get_dimension

    function get_attribute(self, varname) result(att)
        implicit none
        class(read_nc),     intent(in)  :: self
        character(len=*),   intent(in)  :: varname
        real                            :: att

        call check( nf90_get_att(self % ncid, nf90_global, varname, att) )
    end function get_attribute

    subroutine get_variable_0d(self, varname, variable)
        implicit none
        class(read_nc),     intent(in)  :: self
        character(len=*),   intent(in)  :: varname
        real,               intent(out) :: variable

        integer                         :: varid

        call check( nf90_inq_varid(self % ncid, varname, varid)    )
        call check( nf90_get_var  (self % ncid,   varid, variable) )
    end subroutine get_variable_0d

    subroutine get_variable_1d(self, varname, variable)
        implicit none
        class(read_nc),     intent(in)  :: self
        character(len=*),   intent(in)  :: varname
        real, dimension(:), intent(out) :: variable

        integer                         :: varid

        call check( nf90_inq_varid(self % ncid, varname, varid)    )
        call check( nf90_get_var  (self % ncid,   varid, variable) )
    end subroutine get_variable_1d

    subroutine get_variable_2d(self, varname, variable)
        implicit none
        class(read_nc),       intent(in)  :: self
        character(len=*),     intent(in)  :: varname
        real, dimension(:,:), intent(out) :: variable

        integer                           :: varid

        call check( nf90_inq_varid(self % ncid, varname, varid)    )
        call check( nf90_get_var  (self % ncid,   varid, variable) )
    end subroutine get_variable_2d

    subroutine get_variable_3d(self, varname, variable)
        implicit none
        class(read_nc),         intent(in)  :: self
        character(len=*),       intent(in)  :: varname
        real, dimension(:,:,:), intent(out) :: variable

        integer                             :: varid

        call check( nf90_inq_varid(self % ncid, varname, varid)    )
        call check( nf90_get_var  (self % ncid,   varid, variable) )
    end subroutine get_variable_3d

    !==========================================================================

    subroutine copy_header_from(self, you)
        implicit none
        class(write_nc), intent(in)   :: self
        integer,         intent(in)   :: you

        !nc file headers
        integer           :: ndims, nvars, natts, unlimitid
        integer           :: varnatts
        character(len=50) :: dimname, attname

        !local              
        integer           :: dummyid, dimval
        integer           :: dimid, attid, varid

        call check( nf90_inquire(you, ndims, nvars, natts, unlimitid) )

        do dimid = 1, ndims
            call check( nf90_inquire_dimension(you, dimid, dimname, dimval) )
            if(dimid == unlimitid) then
                call check( nf90_def_dim(self % ncid, dimname, nf90_unlimited, dummyid) )
            else
                call check( nf90_def_dim(self % ncid, dimname, dimval, dummyid) )
            end if
        end do

        do attid = 1, natts
            call check( nf90_inq_attname(you, nf90_global, attid, attname) )
            call check( nf90_copy_att(you, nf90_global, attname, &
                                      self % ncid, nf90_global) )
        end do

        allocate(vars(nvars), I_am_not_write(nvars))
        I_am_not_write = .true.

        do varid = 1, nvars
            call check( nf90_inquire_variable(you, varid, name=vars(varid)%name, xtype=vars(varid)%xtype, &
                ndims=vars(varid)%ndims) )   !get number of dimension then allocate array

            allocate(vars(varid)%dimids(vars(varid)%ndims))

            call check( nf90_inquire_variable(you, varid, dimids=vars(varid)%dimids, natts=varnatts) )
            call check( nf90_def_var(self % ncid, vars(varid)%name, vars(varid)%xtype, vars(varid)%dimids, dummyid) )

            do attid = 1, varnatts
                call check( nf90_inq_attname(you, dummyid, attid, attname) )
                call check( nf90_copy_att(you, dummyid, attname, self % ncid, dummyid) )
            end do

            vars(varid)%dims = 1

            do dimid = 1, vars(varid)%ndims
                call check( nf90_inquire_dimension(you, vars(varid)%dimids(dimid), len=vars(varid)%dims(dimid)) )
            end do
        end do

        call check( nf90_enddef(self % ncid) )
    end subroutine copy_header_from

    subroutine write_variable_1d(self, varname, grid)
        implicit none
        class(write_nc),    intent(in) :: self
        character(len=*),   intent(in) :: varname
        real, dimension(:), intent(in) :: grid
        
        integer                        :: varid
        logical                        :: find_var

        find_var = .false.
        do varid = 1, size(vars)
            if(varname == trim(vars(varid)%name)) then
                find_var = .true.
                I_am_not_write(varid) = .false.
                exit
            end if
        end do
        if(find_var) call check( nf90_put_var(self % ncid, varid, grid) )
    end subroutine write_variable_1d

    subroutine write_variable_2d(self, varname, grid)
        implicit none
        class(write_nc),      intent(in) :: self
        character(len=*),     intent(in) :: varname
        real, dimension(:,:), intent(in) :: grid
        
        integer                          :: varid
        logical                          :: find_var

        find_var = .false.
        do varid = 1, size(vars)
            if(varname == trim(vars(varid)%name)) then
                find_var = .true.
                I_am_not_write(varid) = .false.
                exit
            end if
        end do
        if(find_var) call check( nf90_put_var(self % ncid, varid, grid) )
    end subroutine write_variable_2d

    subroutine write_variable_3d(self, varname, grid)
        implicit none
        class(write_nc),        intent(in) :: self
        character(len=*),       intent(in) :: varname
        real, dimension(:,:,:), intent(in) :: grid
        
        integer                            :: varid
        logical                            :: find_var

        find_var = .false.
        do varid = 1, size(vars)
            if(varname == trim(vars(varid)%name)) then
                find_var = .true.
                I_am_not_write(varid) = .false.
                exit
            end if
        end do
        if(find_var) call check( nf90_put_var(self % ncid, varid, grid) )
    end subroutine write_variable_3d

    subroutine write_variable_others(self, varname, you)
        implicit none
        class(write_nc),  intent(in) :: self
        character(len=*), intent(in) :: varname
        integer,          intent(in) :: you
        
        integer                                            :: varid
        real,              dimension(:,:,:,:), allocatable :: tmp_r4
        double precision,  dimension(:,:,:,:), allocatable :: tmp_r8
        integer,           dimension(:,:,:,:), allocatable :: tmp_i4
        character(len=50)                                  :: string

        if(varname == 'OTHERS') then
            do varid = 1, size(vars)
                if(I_am_not_write(varid)) then
                    select case(vars(varid)%xtype)
                    case (nf90_char)
                        call check( nf90_get_var(        you, varid, string, count=[vars(varid)%dims(1)]) )
                        call check( nf90_put_var(self % ncid, varid, string, count=[vars(varid)%dims(1)]) )
                    case (nf90_int)
                        allocate(tmp_i4(vars(varid)%dims(1), &
                                        vars(varid)%dims(2), &
                                        vars(varid)%dims(3), &
                                        vars(varid)%dims(4)))
                        call check( nf90_get_var(        you, varid, tmp_i4) )
                        call check( nf90_put_var(self % ncid, varid, tmp_i4) )
                        deallocate(tmp_i4)
                    case (nf90_double)
                        allocate(tmp_r8(vars(varid)%dims(1), &
                                        vars(varid)%dims(2), &
                                        vars(varid)%dims(3), &
                                        vars(varid)%dims(4)))
                        call check( nf90_get_var(        you, varid, tmp_r8) )
                        call check( nf90_put_var(self % ncid, varid, tmp_r8) )
                        deallocate(tmp_r8)
                    case (nf90_float)
                        allocate(tmp_r4(vars(varid)%dims(1), &
                                        vars(varid)%dims(2), &
                                        vars(varid)%dims(3), &
                                        vars(varid)%dims(4)))
                        call check( nf90_get_var(        you, varid, tmp_r4) )
                        call check( nf90_put_var(self % ncid, varid, tmp_r4) )
                        deallocate(tmp_r4)
                    end select
                end if
            end do
        else
            stop "Unknown variable name"
        end if
    end subroutine write_variable_others

    !=============================================================================

    subroutine check(ierr)
        implicit none
        integer, intent(in) :: ierr
        
        if(ierr /= nf90_noerr) then
!           call letkf_abort(trim(nf90_strerror(ierr)))
            stop "Netcdf error"
        end if
    end subroutine check

end module netcdf_io
