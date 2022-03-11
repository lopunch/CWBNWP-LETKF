module projection
    use config, only : sta_lon, cen_lat, truelat1, truelat2, truelat2
    use param,  only : pi, d2r, earthradius

    implicit none

    private
    public  :: proj_type

    type proj_type
        private

        real :: lon0, n, f, rh0
    contains
        procedure, public, pass(self) :: init => proj_init
        procedure, public, pass(self) :: lonlat_to_xy
    end type proj_type

contains

    pure subroutine proj_init(self)
        implicit none
        class(proj_type), intent(in out) :: self
        real                             :: lat0, lat1, lat2


        lat0          = cen_lat  * d2r
        lat1          = truelat1 * d2r
        lat2          = truelat2 * d2r
        self % lon0   = sta_lon  * d2r
        self % n      = log( cos(lat1) / cos(lat2) ) / &
                        log( tan(0.5*(0.5*pi+lat2)) * cotan(0.5*(0.5*pi+lat1)) )
        self % f      = cos(lat1) * exp(self % n * log( tan(0.5*(0.5*pi+lat1)) )) / self % n
        self % rh0    = earthradius * self % f * exp(self % n * log(cotan(0.5*(0.5*pi+lat0))))
    end subroutine proj_init

    pure function lonlat_to_xy(self, lon, lat) result(xy)
        implicit none
        class(proj_type),   intent(in)           :: self
        real,               intent(in)           :: lon, lat
        real, dimension(2)                       :: xy

        real                                     :: rh, dlon

        rh    = earthradius * self % f * exp(self % n * log(cotan(0.5*(0.5*pi+lat*d2r))))
        dlon  = self % n * (lon*d2r - self % lon0)

        xy(1) =              rh*sin(dlon) 
        xy(2) = self % rh0 - rh*cos(dlon)
    end function lonlat_to_xy

end module projection
