module grid
    
    use param
    use netcdf_io, only : read_nc, write_nc, open_wrf_ncdf, create_wrf_ncdf
    use config,    only : wrf_mp_physics,  &
                          wrf_mp_hail_opt, &
                          wrf_hypsometric_opt, &
                          NT2LOG, NT2Dm, NT2D0, NT2D6, NT2De

    implicit none

    private
    public  :: grid_structure!, latlon_to_dist

    !reference information
!   real, protected    :: lon0, lat0, lat1, lat2, n, f, rh0
!   real, protected    :: swxy11(2)

    type grid_structure

        integer :: nx,   ny,   nz,                     &
                   model_MP_momentR, model_MP_momentS, &
                   model_MP_momentG, model_MP_momentH

        real    :: dx, dy, cen_lon, cen_lat, truelat1, truelat2, sta_lon

        logical :: model_MP_graupel, model_MP_hail

        real    :: model_uR_value, model_uS_value,     &
                   model_uG_value, model_uH_value,     &
                   model_rhoR_value, model_rhoS_value, &
                   model_rhoG_value, model_rhoH_value, &
                   MASSaR, MASSaS, MASSaG, MASSaH,     &
                   MASSbR, MASSbS, MASSbG, MASSbH,     &
                   LDARmax,LDASmax,LDAGmax,LDAHmax,    &
                   LDARmin,LDASmin,LDAGmin,LDAHmin

        real, dimension(:,:),   allocatable :: xlon,   xlat,   &
                                               xlon_u, xlat_u, &
                                               xlon_v, xlat_v
        real, dimension(:,:),   allocatable :: mu, psfc

        real, dimension(:,:,:), allocatable ::   u,   v,   w,   t,  p, ph, qv
        real, dimension(:,:,:), allocatable ::  qr,  qs,  qg,  qh, &
                                               nqr, nqs, nqg, nqh
        real, dimension(:,:,:), allocatable :: rhoa!, DmR, DmS, DmG, DmH

    contains

        procedure, public, pass(self) :: define_wrf_mp_physics
        procedure, public, pass(self) :: read_model
        procedure, public, pass(self) :: write_model
        procedure, public, pass(self) :: distribute_geoinfo

    end type grid_structure

    character(len=:), allocatable     :: input_filename

contains

    subroutine define_wrf_mp_physics(self)
        implicit none
        class(grid_structure), intent(in out)   :: self

        !local
        character(len=255)                      :: message

        self % model_MP_momentR = 1
        self % model_MP_momentS = 1
        self % model_MP_momentG = 1
        self % model_MP_momentH = 1

        select case ( wrf_mp_physics )
        case ( WRFMPLIN )
!           if(is_root) print *, 'wrf microphysics : Lin et al. scheme'
            self % model_MP_graupel = .true.
            self % model_MP_hail    = .false.

        case ( WRFMPWSM5 )
!           if(is_root) print *, 'wrf microphysics : WSM5 scheme'
            self % model_MP_graupel = .false.
            self % model_MP_hail    = .false.

        case ( WRFMPWSM6 )
!           if(is_root) print *, 'wrf microphysics : WSM6 scheme'
            if ( wrf_mp_hail_opt == 0 ) then
                self % model_MP_graupel = .true.
                self % model_MP_hail    = .false.
            else
                self % model_MP_graupel = .false.
                self % model_MP_hail    = .true.
            end if

        case ( WRFMPGSFCGCE )
!           if(is_root) print *, 'wrf microphysics : Goddard scheme'
            if ( wrf_mp_hail_opt == 0 ) then
                self % model_MP_graupel = .true.
                self % model_MP_hail    = .false.
            else
                self % model_MP_graupel = .false.
                self % model_MP_hail    = .true.
            end if

        case ( WRFMPTHOMPSON )
!           if(is_root) print *, 'wrf microphysics : New Thompson et al. scheme'
            self % model_MP_graupel = .true.
            self % model_MP_hail    = .false.
            self % model_MP_momentR = 2
            self % model_uR_value   = 0.
            self % model_rhoR_value = 1000.
            self % MASSaR  = pi/6. * self % model_rhoR_value
            self % MASSbR  = 3.
            self % MASSaS  = 0.069
            self % MASSbS  = 2.
            self % LDARmax = 1.E6
            self % LDARmin = 1.

        case ( WRFMPMILBRANDT )
!           if(is_root) print *, 'wrf microphysics : Milbrandt-Yau scheme'
            self % model_MP_graupel = .true.
            self % model_MP_hail    = .true.
            self % model_MP_momentR = 2
            self % model_MP_momentS = 2
            self % model_MP_momentG = 2
            self % model_MP_momentH = 2
            self % model_uR_value   = 0.
            self % model_uS_value   = 0.
            self % model_uG_value   = 0.
            self % model_uH_value   = 0.
            self % model_rhoR_value = 1000.
            self % model_rhoS_value =  100.
            self % model_rhoG_value =  400.
            self % model_rhoH_value =  900.
            self % MASSaR  = pi/6. * self % model_rhoR_value
            self % MASSbR  = 3.
            self % MASSaS  = 0.1597
            self % MASSbS  = 2.078
            self % MASSaG  = pi/6. * self % model_rhoG_value
            self % MASSbG  = 3.
            self % MASSaH  = pi/6. * self % model_rhoH_value
            self % MASSbH  = 3.
            self % LDARmax = 1.E6
            self % LDARmin = 1.
            self % LDASmax = 1.E10
            self % LDASmin = 1.
            self % LDAGmax = 1.E6
            self % LDAGmin = 1.
            self % LDAHmax = 1.E6
            self % LDAHmin = 1.

        case ( WRFMPMORR )
!           if(is_root) print *, 'wrf microphysics : Morrison scheme'
            if ( wrf_mp_hail_opt == 0 ) then
                self % model_MP_graupel = .true.
                self % model_MP_hail    = .false.
            else
                self % model_MP_graupel = .false.
                self % model_MP_hail    = .true.
            end if
            self % model_MP_momentR = 2
            self % model_MP_momentS = 2
            self % model_MP_momentG = 2
            self % model_MP_momentH = 2
            self % model_uR_value   = 0.
            self % model_uS_value   = 0.
            self % model_uG_value   = 0.
            self % model_uH_value   = 0.
            self % model_rhoR_value = 997.
            self % model_rhoS_value = 100.
            self % model_rhoG_value = 400.
            self % model_rhoH_value = 900.
            self % MASSaR  = pi/6. * self % model_rhoR_value
            self % MASSbR  = 3.
            self % MASSaS  = pi/6. * self % model_rhoS_value
            self % MASSbS  = 3.
            self % MASSaG  = pi/6. * self % model_rhoG_value
            self % MASSbG  = 3.
            self % MASSaH  = pi/6. * self % model_rhoH_value
            self % MASSbH  = 3.
            self % LDARmax = 1./20.E-6 
            self % LDARmin = 1./2800.E-6 
            self % LDASmax = 1./10.E-6
            self % LDASmin = 1./2000.E-6
            self % LDAGmax = 1./20.E-6 
            self % LDAGmin = 1./2000.E-6
            self % LDAHmax = 1./20.E-6
            self % LDAHmin = 1./2000.E-6

        case ( WRFMPWDM5 )
!           if(is_root) print *, 'wrf microphysics : WDM5 scheme'
            self % model_MP_graupel = .false.
            self % model_MP_hail    = .false.
            self % model_MP_momentR = 2
            self % model_uR_value   = 0.
            self % model_rhoR_value = 1000.
            self % MASSaR  = pi/6. * self % model_rhoR_value
            self % MASSbR  = 3.
            self % LDARmax = 5.E4
            self % LDARmin = 2.E3

        case ( WRFMPWDM6 )
!           if(is_root) print *, 'wrf microphysics : WDM6 scheme'
            if ( wrf_mp_hail_opt .eq. 0 ) then
                self % model_MP_graupel = .true.
                self % model_MP_hail    = .false.
            else
                self % model_MP_graupel = .false.
                self % model_MP_hail    = .true.
            end if
            self % model_MP_momentR = 2
            self % model_uR_value   = 0.
            self % model_rhoR_value = 1000.
            self % MASSaR  = pi/6. * self % model_rhoR_value
            self % MASSbR  = 3.
            self % LDARmax = 5.E4
            self % LDARmin = 2.E3

        case default
            write(message, '(A)') 'please set wrf_mp_physics = ',WRFMPLIN,WRFMPWSM5,WRFMPWSM6,&
                        WRFMPGSFCGCE,WRFMPMILBRANDT,WRFMPMORR,WRFMPWDM5,WRFMPWDM6         
!           call letkf_abort(message)
        end select

    end subroutine define_wrf_mp_physics

    subroutine read_model(self, filename)

        implicit none

        class(grid_structure), intent(in out)  :: self
        character(len=*),      intent(in)      :: filename

        !local
        type(read_nc)                          :: nc
        integer                                :: k
        real                                   :: t00, p00, tlp, tiso, p_strat, tlp_strat, p_top
        real                                   :: p00_inv, p_strat_inv, p1000mb_inv
        real, dimension(:),     allocatable    :: rdnw, znw, znu
        real, dimension(:,:),   allocatable    :: mub, alb, al, pfu, pfd, phm
        real, dimension(:,:,:), allocatable    :: pb, phb
        real, dimension(:,:,:), allocatable    :: temp, t_init


        !open file
        nc = open_wrf_ncdf(filename)
        
        !save filename for output
        allocate(character(len=len_trim(filename)) :: input_filename)
        input_filename = filename

        !global attributes
        self % nx       = nc % get_dimension('west_east')
        self % ny       = nc % get_dimension('south_north')
        self % nz       = nc % get_dimension('bottom_top')
        self % dx       = nc % get_attribute('DX')
        self % dy       = nc % get_attribute('DY')
        self % cen_lon  = nc % get_attribute('CEN_LON')
        self % cen_lat  = nc % get_attribute('CEN_LAT')
        self % truelat1 = nc % get_attribute('TRUELAT1')
        self % truelat2 = nc % get_attribute('TRUELAT2')
        self % sta_lon  = nc % get_attribute('STAND_LON')


        !allocate variables
        associate( nx => self % nx, ny => self % ny, nz => self % nz )
            allocate(self % xlon  (nx,  ny),  &
                     self % xlat  (nx,  ny),  & 
                     self % xlon_u(nx+1,ny),  &
                     self % xlat_u(nx+1,ny),  &
                     self % xlon_v(nx,ny+1),  &
                     self % xlat_v(nx,ny+1),  &
                     self % mu  (nx, ny),     &
                     self % psfc(nx, ny),     &
                     self % u (nx+1, ny, nz), &
                     self % v (nx, ny+1, nz), &
                     self % w (nx, ny, nz+1), &
                     self % ph(nx, ny, nz+1), &
                     self % t (nx, ny, nz),   &
                     self % p (nx, ny, nz),   &
                     self % qv(nx, ny, nz),   &
                     self % qr(nx, ny, nz),   &
                     self % qs(nx, ny, nz))

            if ( self % model_MP_graupel ) allocate( self % qg(nx, ny, nz) )
            if ( self % model_MP_hail    ) allocate( self % qh(nx, ny, nz) )

            if( self % model_MP_momentR >= 2 .or. self % model_MP_momentS >= 2 .or. &
                self % model_MP_momentG >= 2 .or. self % model_MP_momentH >= 2 ) then
                allocate( self % rhoa(nx, ny, nz) )

                if ( self % model_MP_momentR >= 2 ) then
                    allocate( self % nqr(nx, ny, nz) )
                   !if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) allocate( self % DmR(nx, ny, nz) )
                end if

                if ( self % model_MP_momentS >= 2 ) then
                    allocate( self % nqs(nx, ny, nz) )
                   !if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) allocate( self % DmS(nx, ny, nz) )
                end if

                if ( self % model_MP_momentG >= 2 ) then
                    allocate( self % nqg(nx, ny, nz) )
                   !if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) allocate( self % DmG(nx, ny, nz) )
                end if

                if ( self % model_MP_momentH >= 2 ) then
                    allocate( self % nqh(nx, ny, nz) )
                   !if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) allocate( self % DmH(nx, ny, nz) )
                end if
            end if

            allocate( phb(nx, ny, nz+1), pb (nx, ny, nz) )
        end associate

        !read variables
        call nc % get_variable('XLONG',   self % xlon)
        call nc % get_variable('XLAT',    self % xlat)
        call nc % get_variable('XLONG_U', self % xlon_u)
        call nc % get_variable('XLAT_U',  self % xlat_u)
        call nc % get_variable('XLONG_V', self % xlon_v)
        call nc % get_variable('XLAT_V',  self % xlat_v)
        call nc % get_variable('PSFC',    self % psfc)
        call nc % get_variable('MU',      self % mu)
        call nc % get_variable('U',       self % u)
        call nc % get_variable('V',       self % v)
        call nc % get_variable('W',       self % w)
        call nc % get_variable('PH',      self % ph)
        call nc % get_variable('T',       self % t)
        call nc % get_variable('P',       self % p)
        call nc % get_variable('QVAPOR',  self % qv)
        call nc % get_variable('QRAIN',   self % qr)
        call nc % get_variable('QSNOW',   self % qs)
        call nc % get_variable('PHB',            phb)
        call nc % get_variable('PB',             pb)


        if ( self % model_MP_graupel ) then
            call nc % get_variable('QGRAUP',  self % qg)
            if ( self % model_MP_momentG >= 2 ) then
                call nc % get_variable('QNGRAUPEL', self % nqg)
            end if
        end if

        if ( self % model_MP_hail ) then
            call nc % get_variable('QHAIL',   self % qh)
            if ( self % model_MP_momentH >= 2 ) then
                call nc % get_variable('QNHAIL',    self % nqh)
            end if
        end if

        if ( self % model_MP_momentR >= 2 ) then
            call nc % get_variable('QNRAIN',  self % nqr)
        end if

        if ( self % model_MP_momentS >= 2 ) then
            call nc % get_variable('QNSNOW',  self % nqs)
        end if

        !check hydrometers

                                       where( self % qr < 0.0 ) self % qr = 0.0
                                       where( self % qs < 0.0 ) self % qs = 0.0
        if ( self % model_MP_graupel ) where( self % qg < 0.0 ) self % qg = 0.0
        if ( self % model_MP_hail    ) where( self % qh < 0.0 ) self % qh = 0.0


        !Calculate reference information
!       lon0   = self % sta_lon  * d2r
!       lat0   = self % cen_lat  * d2r
!       lat1   = self % truelat1 * d2r
!       lat2   = self % truelat2 * d2r
!       n      = log(cos(lat1)/cos(lat2)) / &
!                log(tan(0.25*pi+0.5*lat2)*cotan(0.25*pi+0.5*lat1))
!      !f      = cos(lat1)*tan(0.25*pi+0.5*lat1)**n / n
!      !rh0    = earthradius * f * cotan(0.25*pi+0.5*lat0)**n
!       f      = cos(lat1)*exp(n*log(tan(0.25*pi+0.5*lat1))) / n
!       rh0    = earthradius * f * exp(n*log(cotan(0.25*pi+0.5*lat0)))
!       swxy11 = latlon_to_dist( [self % xlat(1,1), self % xlon(1,1)] )


        !if use double moment scheme, do some calculations
        if( self % model_MP_momentR >= 2 .or. self % model_MP_momentS >= 2 .or. &
            self % model_MP_momentG >= 2 .or. self % model_MP_momentH >= 2 ) then

            call nc % get_variable('T00',       t00)
            call nc % get_variable('P00',       p00)
            call nc % get_variable('TLP',       tlp)
            call nc % get_variable('TISO',      tiso)
            call nc % get_variable('P_STRAT',   p_strat)
            call nc % get_variable('TLP_STRAT', tlp_strat)

            p00_inv     = 1.0 / p00
            p_strat_inv = 1.0 / p_strat
            p1000mb_inv = 1.0 / p1000mb

            associate( nx => self % nx, ny => self % ny, nz => self % nz )
                allocate( temp(nx, ny, nz), t_init(nx, ny, nz), mub(nx, ny) )
            end associate

            temp   = max( tiso, t00 + tlp * log(pb * p00_inv) )
            where ( pb < p_strat ) temp   =  tiso + tlp_strat * log(pb * p_strat_inv)

           !t_init = temp * (p00 / pb) ** (r_d/cp)
            t_init = exp(r_d/cp * log(temp * (p00 / pb)))

            !get full mu in mub
            call nc % get_variable('MUB',   mub)
            call saxpy(size(mub), 1.0,  self % mu, 1,  mub, 1)

            if (wrf_hypsometric_opt == 1) then

                associate( nx => self % nx, ny => self % ny, nz => self % nz )
                    allocate(rdnw(nz), alb(nx, ny), al(nx, ny))
                end associate

                call nc % get_variable('RDNW', rdnw)
                do k = 1, self % nz
                   !alb = (r_d*p1000mb_inv) * t_init(:,:,k) * (pb(:,:,k)*p1000mb_inv) ** cvpm
                    alb = (r_d*p1000mb_inv) * t_init(:,:,k) * exp(cvpm * log(pb(:,:,k)*p1000mb_inv))
                    al  = -1./mub * (alb*self%mu+rdnw(k) * (self%ph(:,:,k+1)-self%ph(:,:,k)))
                    self % rhoa(:,:,k) = 1./(alb+al)
                end do

                deallocate(rdnw,alb,al)

            elseif (wrf_hypsometric_opt == 2) then

                associate( nx => self % nx, ny => self % ny, nz => self % nz )
                    allocate(znw(nz+1), znu(nz),                    &
                             pfu(nx, ny), pfd(nx, ny), phm(nx, ny), &
                             alb(nx, ny),  al(nx, ny))
                end associate

                call nc % get_variable('P_TOP', p_top)
                call nc % get_variable('ZNW',   znw)
                call nc % get_variable('ZNU',   znu)

                do k = 1, self % nz
                    pfu = mub * znw(k+1) + p_top
                    pfd = mub * znw(k)   + p_top
                    phm = mub * znu(k)   + p_top
                   !alb = (r_d*p1000mb_inv) * t_init(:,:,k) * (pb(:,:,k)*p1000mb_inv) ** cvpm
                    alb = (r_d*p1000mb_inv) * t_init(:,:,k) * exp(cvpm * log(pb(:,:,k)*p1000mb_inv))
                    al  = (self%ph (:,:,k+1) - self%ph (:,:,k)  + &
                                phb(:,:,k+1) -      phb(:,:,k)) / (phm * log(pfd/pfu)) - alb
                    self % rhoa(:,:,k) = 1./(alb+al)
                end do

                deallocate(alb,al,pfu,pfd,phm,znw,znu)
            end if

            deallocate(mub, temp, t_init)

!           if ( self % model_MP_momentR >= 2 ) then
!               call nqx_to_DmX(self % qr,      &
!                               self % rhoa,    &
!                               self % LDARmax, &
!                               self % LDARmin, &
!                               self % MASSaR,  &
!                               self % MASSbR,  &
!                               self % model_uR_value, &
!                               .true.,         &
!                               self % nqr,     &
!                               self % DmR)
!           end if

!           if ( self % model_MP_momentS >= 2 ) then
!               call nqx_to_DmX(self % qs,      &
!                               self % rhoa,    &
!                               self % LDASmax, &
!                               self % LDASmin, &
!                               self % MASSaS,  &
!                               self % MASSbS,  &
!                               self % model_uS_value, &
!                               .true.,         &
!                               self % nqs,     &
!                               self % DmS)
!           end if

!           if ( self % model_MP_momentG >= 2 ) then
!               call nqx_to_DmX(self % qg,      &
!                               self % rhoa,    &
!                               self % LDAGmax, &
!                               self % LDAGmin, &
!                               self % MASSaG,  &
!                               self % MASSbG,  &
!                               self % model_uG_value, &
!                               .true.,         &
!                               self % nqg,     &
!                               self % DmG)
!           end if

!           if ( self % model_MP_momentH >= 2 ) then
!               call nqx_to_DmX(self % qh,      &
!                               self % rhoa,    &
!                               self % LDAHmax, &
!                               self % LDAHmin, &
!                               self % MASSaH,  &
!                               self % MASSbH,  &
!                               self % model_uH_value, &
!                               .true.,         &
!                               self % nqh,     &
!                               self % DmH)
!           end if
        end if

        !end of reading file, close it
        call nc % close_reading

        !sum based-state and perturbed pressure to get full pressure
        call saxpy(size(self % ph), 1.0, phb, 1, self % ph, 1)
        call saxpy(size(self %  p), 1.0,  pb, 1, self %  p, 1)

        deallocate(phb, pb)

    end subroutine read_model

    subroutine write_model(self, filename)

        implicit none
        class(grid_structure), intent(in out)  :: self
        character(len=*),      intent(in)      :: filename

        type( read_nc) :: nc1
        type(write_nc) :: nc2

        nc1 =   open_wrf_ncdf(input_filename)
        nc2 = create_wrf_ncdf(      filename)

        call nc2 % copy_header_from( nc1 % ncid )
        call nc2 % write_variable('U',      self % u)
        call nc2 % write_variable('V',      self % v)
        call nc2 % write_variable('W',      self % w)
        call nc2 % write_variable('T',      self % t)
        call nc2 % write_variable('QVAPOR', self % qv)
        call nc2 % write_variable('QRAIN',  self % qr)
        call nc2 % write_variable('QSNOW',  self % qs)
        if ( self % model_MP_graupel ) then
            call nc2 % write_variable('QGRAUP', self % qg)
        end if
        if ( self % model_MP_hail ) then
            call nc2 % write_variable('QHAIL',  self % qh)
        end if

        if( self % model_MP_momentR >= 2 ) then
!           call nqx_to_DmX(self % qr,      &
!                           self % rhoa,    &
!                           self % LDARmax, &
!                           self % LDARmin, &
!                           self % MASSaR,  &
!                           self % MASSbR,  &
!                           self % model_uR_value, &
!                           .false.,        &
!                           self % nqr,     &
!                           self % DmR)
            call nc2 % write_variable('QNRAIN',  self % nqr)
        end if

        if( self % model_MP_momentS >= 2 ) then
!           call nqx_to_DmX(self % qs,      &
!                           self % rhoa,    &
!                           self % LDASmax, &
!                           self % LDASmin, &
!                           self % MASSaS,  &
!                           self % MASSbS,  &
!                           self % model_uS_value, &
!                           .false.,        &
!                           self % nqs,     &
!                           self % DmS)
            call nc2 % write_variable('QNSNOW',  self % nqs)
        end if

        if( self % model_MP_momentG >= 2 ) then
!           call nqx_to_DmX(self % qg,      &
!                           self % rhoa,    &
!                           self % LDAGmax, &
!                           self % LDAGmin, &
!                           self % MASSaG,  &
!                           self % MASSbG,  &
!                           self % model_uG_value, &
!                           .false.,        &
!                           self % nqg,     &
!                           self % DmG)
            call nc2 % write_variable('QNGRAUPEL',  self % nqg)
        end if

        if( self % model_MP_momentH >= 2 ) then
!           call nqx_to_DmX(self % qh,      &
!                           self % rhoa,    &
!                           self % LDAHmax, &
!                           self % LDAHmin, &
!                           self % MASSaH,  &
!                           self % MASSbH,  &
!                           self % model_uH_value, &
!                           .false.,        &
!                           self % nqh,     &
!                           self % DmH)
            call nc2 % write_variable('QNHAIL',  self % nqh)
        end if

        call nc2 % write_variable('OTHERS',  nc1 % ncid)
        call nc2 % close_writing
        call nc1 % close_reading


        !release allocated arrays
        deallocate(self % xlon  ,  &
                   self % xlat  ,  & 
                   self % xlon_u,  &
                   self % xlat_u,  &
                   self % xlon_v,  &
                   self % xlat_v,  &
                   self % mu  ,    &
                   self % psfc,    &
                   self % u ,      &
                   self % v ,      &
                   self % w ,      &
                   self % ph,      &
                   self % t ,      &
                   self % p ,      &
                   self % qv,      &
                   self % qr,      &
                   self % qs)

        if ( self % model_MP_graupel ) deallocate( self % qg )
        if ( self % model_MP_hail    ) deallocate( self % qh )

        if( self % model_MP_momentR >= 2 .or. self % model_MP_momentS >= 2 .or. &
            self % model_MP_momentG >= 2 .or. self % model_MP_momentH >= 2 ) then
            deallocate( self % rhoa )

            if ( self % model_MP_momentR >= 2 ) then
                deallocate( self % nqr )
!               if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) deallocate( self % DmR )
            end if

            if ( self % model_MP_momentS >= 2 ) then
                deallocate( self % nqs )
!               if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) deallocate( self % DmS )
            end if

            if ( self % model_MP_momentG >= 2 ) then
                deallocate( self % nqg )
!               if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) deallocate( self % DmG )
            end if

            if ( self % model_MP_momentH >= 2 ) then
                deallocate( self % nqh )
!               if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) deallocate( self % DmH )
            end if
        end if

    end subroutine write_model

    subroutine distribute_geoinfo(self, root)
        use mpi_util
        implicit none

        class(grid_structure), intent(in out) :: self
        integer,               intent(in)     :: root
        integer                               :: mpierr

        call request_append
        call mpi_ibcast(self%nx, 1, mpi_integer, root, mpi_comm_world, req_ptr%req, mpierr)
        call request_append
        call mpi_ibcast(self%ny, 1, mpi_integer, root, mpi_comm_world, req_ptr%req, mpierr)
        call request_append
        call mpi_ibcast(self%nz, 1, mpi_integer, root, mpi_comm_world, req_ptr%req, mpierr)
    end subroutine distribute_geoinfo

    !=================================================================================================

!   pure subroutine nqx_to_DmX(qx, rhoa, LDAmax, LDAmin, a, b, u, forward, nqx, DmX)
!       implicit none
!       real, dimension(:,:,:), intent(in)     :: qx, rhoa
!       real,                   intent(in)     :: LDAmax, LDAmin, a, b, u
!       logical,                intent(in)     :: forward
!       real, dimension(:,:,:), intent(in out) :: nqx, DmX

!       if(forward) then
!           nqx = nqx * rhoa
!           nqx = check_Nt(LDAmax, LDAmin, a, b, u, rhoa, qx, nqx)
!           where( qx <= 1.E-6 .or. nqx < 0. ) nqx = 0.

!           if ( NT2LOG ) then
!               nqx = log(nqx)
!               where( qx <= 1.E-6 .or. nqx < 0. ) nqx = -5.
!           end if

!           if ( NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) then
!               if ( NT2Dm ) DmX = calc_Dm(rhoa,    qx, nqx, a, b, forward)
!               if ( NT2D0 ) DmX = calc_D0(rhoa, u, qx, nqx, a, b, forward)
!               if ( NT2De ) DmX = calc_De(rhoa, u, qx, nqx, a, b, forward)
!               !if ( NT2D6 ) then ; add in the future
!               !end if
!               !-----------------
!               where( qx <= 1.E-6 .or. nqx < 0. ) DmX = 0.
!               where( DmX > 10.E-3 ) DmX  = 10.E-3 
!               DmX = DmX * 1.E3 ! m to mm
!           end if
!       else
!           if(NT2LOG) then

!               where ( nqx /= nqx ) nqx = -5.
!               nqx = exp(nqx)
!               nqx = check_Nt(LDAmax, LDAmin, a, b, u, rhoa, qx, nqx)
!               nqx = nqx / rhoa
!               where ( qx <= 1.E-6 ) nqx = 0.

!           else if(NT2Dm .or. NT2D0 .or. NT2D6 .or. NT2De ) then

!               DmX = DmX * 1.E-3
!               nqx = DmX ** b
!               if ( NT2Dm ) nqx = calc_Dm(rhoa,    qx, nqx, a, b, forward)
!               if ( NT2D0 ) nqx = calc_D0(rhoa, u, qx, nqx, a, b, forward)
!               if ( NT2De ) nqx = calc_De(rhoa, u, qx, nqx, a, b, forward)

!               nqx = check_Nt(LDAmax, LDAmin, a, b, u, rhoa, qx, nqx)
!               nqx = nqx / rhoa

!               where ( qx <= 1.E-6 .or. DmX <= 0. ) nqx  =   0.
!               where ( DmX < 1.E-5 ) nqx  =   0.

!           else

!               nqx = check_Nt(LDAmax, LDAmin, a, b, u, rhoa, qx, nqx)
!               nqx = nqx / rhoa
!               where ( qx <= 1.E-6 ) nqx  =   0.

!           end if
!       end if
!   end subroutine nqx_to_DmX

!   elemental function check_Nt(LDAMAX,LDAMIN,MASSa,MASSb,U,RHOa3D,Qx3D,NT3D) result (NT_OUT)
!       implicit none
!       real, intent(in) :: LDAMAX, LDAMIN
!       real, intent(in) :: MASSa,  MASSb
!       real, intent(in) :: U
!       real, intent(in) :: RHOa3D 
!       real, intent(in) :: Qx3D 
!       real, intent(in) :: NT3D
!       real             :: NT_OUT

!       real :: LDA3D
!       real :: tmp1, tmp1LN
!       real :: tmp2, tmp2LN


!       tmp1   = 1.+U
!       tmp2   = tmp1+MASSb
!       tmp1LN = lgamma(tmp1)
!       tmp2LN = lgamma(tmp2)

!       LDA3D  = exp( (1./MASSb) * ( tmp2LN - tmp1LN + &
!                   log(MASSa)+log(NT3D)-log(RHOa3D*Qx3D)) ) 

!       if( LDA3D > LDAMAX ) then
!           LDA3D  =  LDAMAX
!           NT_OUT = exp( (tmp1LN - tmp2LN) + MASSb*log(LDA3D) +&
!                    log(RHOa3D*Qx3D) - log(MASSa) )
!       else if( LDA3D < LDAMIN ) then
!           LDA3D  =  LDAMIN
!           NT_OUT = exp( (tmp1LN - tmp2LN) + MASSb*log(LDA3D) +&
!                    log(RHOa3D*Qx3D) - log(MASSa) )
!       else
!           NT_OUT = NT3D
!       end if
!   end function check_Nt

!   elemental function calc_Dm(rhoa, qx, nqx, a, b, forward) result(Dm)
!       implicit none
!       real,    intent(in)  :: rhoa, qx, nqx, a, b
!       logical, intent(in)  :: forward
!       real                 :: Dm

!       if(forward) then
!           Dm = ( (rhoa * qx) / (a * nqx) )**(1.0 / b)
!       else
!           Dm =   (rhoa * qx) / (a * nqx)
!       end if
!   end function calc_Dm

!   elemental function calc_D0(rhoa, u, qx, nqx, a, b, forward) result(D0)
!       implicit none
!       real,    intent(in)  :: rhoa, u, qx, nqx, a, b
!       logical, intent(in)  :: forward
!       real                 :: D0
!       real, dimension(4)   :: gamma_u

!       gamma_u(1) = gamma(u+1.0)
!       gamma_u(2) = (u+1.0) * gamma_u(1)
!       gamma_u(3) = (u+2.0) * gamma_u(2)
!       gamma_u(4) = (u+3.0) * gamma_u(3)

!       if(forward) then
!          !D0 = ( ( (rhoa * qx) / (a * nqx) )**(1.0 / b) ) * &
!          !     ( (gamma(u+1.0)/gamma(u+4.0))**(1.0 / b) ) * &
!          !        gamma(u+2.0)/gamma(u+1.0)  
!           D0 = ( ( (rhoa * qx) / (a * nqx) )**(1.0 / b) ) * &
!                ( (gamma_u(1)/gamma_u(4))**(1.0 / b) ) * &
!                   gamma_u(2)/gamma_u(1)
!       else
!          !D0 = ( (rhoa * qx) / (a * nqx) ) * &
!          !     ( (gamma(u+1.0)/gamma(u+4.0)) ) * &
!          !       (gamma(u+2.0)/gamma(u+1.0))**b
!           D0 = ( (rhoa * qx) / (a * nqx) ) * &
!                ( (gamma_u(1)/gamma_u(4)) ) * &
!                  (gamma_u(2)/gamma_u(1))**b
!       end if
!   end function calc_D0

!   elemental function calc_De(rhoa, u, qx, nqx, a, b, forward) result(De)
!       implicit none
!       real,    intent(in)  :: rhoa, u, qx, nqx, a, b
!       logical, intent(in)  :: forward
!       real                 :: De
!       real, dimension(4)   :: gamma_u

!       gamma_u(1) = gamma(u+1.0)
!       gamma_u(2) = (u+1.0) * gamma_u(1)
!       gamma_u(3) = (u+2.0) * gamma_u(2)
!       gamma_u(4) = (u+3.0) * gamma_u(3)

!       if(forward) then
!          !De = ( ( (rhoa * qx) / (a * nqx) )**(1.0 / b) ) * &
!          !     ( (gamma(u+1.0)/gamma(u+4.0))**(1.0 / b) ) * &
!          !        gamma(u+4.0)/gamma(u+3.0)  
!           De = ( ( (rhoa * qx) / (a * nqx) )**(1.0 / b) ) * &
!                ( (gamma_u(1)/gamma_u(4))**(1.0 / b) ) * &
!                   gamma_u(4)/gamma_u(3)  
!       else
!          !De = ( (rhoa * qx) / (a * nqx) ) * &
!          !     ( (gamma(u+1.0)/gamma(u+4.0)) ) * &
!          !       (gamma(u+4.0)/gamma(u+3.0))**b
!           De = ( (rhoa * qx) / (a * nqx) ) * &
!                ( (gamma_u(1)/gamma_u(4)) ) * &
!                  (gamma_u(4)/gamma_u(3))**b
!       end if
!   end function calc_De

    subroutine print_cpu_time(message)
        implicit none
        character(len=*), intent(in) :: message
        real                         :: time

        call cpu_time(time)
        print '(f6.2,2x,a)', time, message
    end subroutine print_cpu_time

!   pure function latlon_to_dist(latlon) result(dist)
!       implicit none
!       real, dimension(2), intent(in) :: latlon  ![ lat,  lon]
!       real, dimension(2)             :: dist    ![disx, disy]
!       real                           :: rh, dlon

!      !rh      = earthradius * f * cotan(0.25*pi+0.5*latlon(1)*d2r)**n
!       rh      = earthradius * f * exp(n*log(cotan(0.25*pi+0.5*latlon(1)*d2r)))
!       dlon    = n*(latlon(2)*d2r-lon0)
!       dist(1) =      rh*sin(dlon)
!       dist(2) =  rh0-rh*cos(dlon)
!   end function latlon_to_dist

end module grid
