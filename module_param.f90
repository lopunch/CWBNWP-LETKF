module param
    implicit none

    type, abstract :: obs_structure
        integer                           :: nobs = 0
        real, dimension(:),   allocatable :: lat, lon, alt
        real, dimension(:,:), allocatable :: xyz          ![   3, nobs]
    end type obs_structure

!--------------------------------------------------------------
! for wrf_mp_physics
!--------------------------------------------------------------
    integer, parameter :: WRFMPLIN  = 2
    integer, parameter :: WRFMPWSM5 = 4
    integer, parameter :: WRFMPWSM6 = 6
    integer, parameter :: WRFMPGSFCGCE   = 7
    integer, parameter :: WRFMPTHOMPSON  = 8
    integer, parameter :: WRFMPMILBRANDT = 9
    integer, parameter :: WRFMPMORR = 10
    integer, parameter :: WRFMPWDM5 = 14
    integer, parameter :: WRFMPWDM6 = 16
    integer, parameter :: WRFMPNSSL2MOM  = 17
    integer, parameter :: WRFMPNSSL1MOM  = 19
    integer, parameter :: WRFMPNSSL2MOMG = 22
!--------------------------------------------------------------
! for WRF gts
!--------------------------------------------------------------
    integer, parameter :: sound     = 1
    integer, parameter :: synop     = 2
    integer, parameter :: pilot     = 3
    integer, parameter :: satem     = 4
    integer, parameter :: geoamv    = 5
    integer, parameter :: polaramv  = 6
    integer, parameter :: airep     = 7
    integer, parameter :: gpspw     = 8
    integer, parameter :: gpsref    = 9
    integer, parameter :: metar     = 10
    integer, parameter :: ships     = 11
    integer, parameter :: ssmi_rv   = 12
    integer, parameter :: ssmi_tb   = 13
    integer, parameter :: ssmt1     = 14
    integer, parameter :: ssmt2     = 15
    integer, parameter :: qscat     = 16
    integer, parameter :: profiler  = 17
    integer, parameter :: buoy      = 18
    integer, parameter :: bogus     = 19
    integer, parameter :: pseudo    = 20
    integer, parameter :: radar     = 21
    integer, parameter :: radiance  = 22
    integer, parameter :: airsr     = 23
    integer, parameter :: sonde_sfc = 24
    integer, parameter :: mtgirs    = 25
    integer, parameter :: tamdar    = 26
    integer, parameter :: tamdar_sfc = 27
    integer, parameter :: rain      = 28
    integer, parameter :: gpseph    = 29
    integer, parameter :: num_gts_indexes = 29

    character(len=14), parameter :: gts_names(num_gts_indexes) = (/ &
          "sound         ", &
          "synop         ", &
          "pilot         ", &
          "satem         ", &
          "geoamv        ", &
          "polaramv      ", &
          "airep         ", &
          "gpspw         ", &
          "gpsrf         ", &
          "metar         ", &
          "ships         ", &
          "ssmi_rv       ", &
          "ssmi_tb       ", &
          "ssmt1         ", &
          "ssmt2         ", &
          "qscat         ", &
          "profiler      ", &
          "buoy          ", &
          "bogus         ", &
          "pseudo        ", &
          "radar         ", &
          "radiance      ", &
          "airs retrieval", &
          "sonde_sfc     ", &
          "mtgirs        ", &
          "tamdar        ", &
          "tamdar_sfc    ", &
          "rain          ", &
          "gpseph        " /)

!-----------------------------------------------------------------
! for radar
!-----------------------------------------------------------------
    integer, parameter  :: dbz = 1
    integer, parameter  :: vr  = 2
    integer, parameter  :: zdr = 3
    integer, parameter  :: kdp = 4
    integer, parameter  :: num_radar_indexes = 4

    character(len=3), parameter :: radar_names(num_radar_indexes) = &
        (/ "MR ", "VR ", "ZDR", "KDP" /)

!-----------------------------------------------------------------
! Constants
!-----------------------------------------------------------------
    real,    parameter :: pi  = acos(-1.)
    real,    parameter :: d2r = pi / 180. 
    real,    parameter :: r2d = 180. / pi
    real,    parameter :: earthradius = 6.37122e6
    real,    parameter :: g = 9.81
    real,    parameter :: p1000mb = 100000.
    real,    parameter :: t0      = 300.
    real,    parameter :: r_d     = 287.
    real,    parameter :: cp      = 7.*r_d*0.5
    real,    parameter :: cv      = cp-r_d
    real,    parameter :: cvpm    = -cv/cp
    real,    parameter :: gc1999  = 2. * sqrt(10./3.)   !Gaspari and Cohn 1999 

!-----------------------------------------------------------------
! Computation
!-----------------------------------------------------------------
    real               :: nmember_inv
    real               :: nmember_1_inv

contains

    subroutine set_ensemble_constants(nmember)
        implicit none
        integer, intent(in) :: nmember

        nmember_inv   = 1.0 /  nmember
        nmember_1_inv = 1.0 / (nmember-1)
    end subroutine set_ensemble_constants

end module param
