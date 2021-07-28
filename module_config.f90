module config
    implicit none

    integer, parameter               :: max_vars = 15

    !observations_nml
    type radar_variable_config
        logical                      :: use_it        = .false.
        real                         :: error         = 1.
        real                         :: err_rej       = 5.
        real                         :: threshold_val = 0.
        real                         :: threshold_alt = 50.
        real,    dimension(max_vars) :: hclr          = 0.
        real,    dimension(max_vars) :: vclr          = 0.
        logical, dimension(max_vars) :: is_assim      = .false.
    end type radar_variable_config

    type gts_variable_config
        real                         :: err_muti = 1.
        real                         :: err_rej  = 5.
        logical, dimension(max_vars) :: is_assim = .false.
    end type gts_variable_config

    !==========================================================

    type radar_config
        real                         :: hclr_rej     = 3.
        real                         :: vclr_rej     = 3.
        type(radar_variable_config)  :: dbz, vr, zdr, kdp
    end type radar_config

    type gts_config
        logical                      :: use_it    = .false.
        real                         :: hclr_rej  = 3.
        real                         :: vclr_rej  = 3.
        real, dimension(max_vars)    :: hclr      = 0.
        real, dimension(max_vars)    :: vclr      = 0.
        type(gts_variable_config)    :: u, v, t, p, q, tpw, ref
    end type gts_config

    !==========================================================

    type(radar_config), target       :: radar_nml
    type(  gts_config), target       :: synop_nml, &
                                        ships_nml, &
                                        metar_nml, &
                                        sound_nml, &
                                        gpspw_nml

    !control_nml
    real,    protected  :: norain_value          = -5.
    logical, protected  :: write_analy_mean      = .true. , &
                           deterministic_update  = .false., &
                           NT2LOG                = .false., &
                           NT2Dm                 = .false., &
                           NT2D0                 = .false., &
                           NT2De                 = .false., &
                           NT2D6                 = .false.
    integer, protected  :: wrf_mp_physics        = -1,   &
                           wrf_mp_hail_opt       = -1,   &
                           wrf_hypsometric_opt   =  2,   &
                           nmember               = -1
    character(len=10), dimension(max_vars), protected :: var_update = ''


    !inflation_nml
    real,    dimension(max_vars), protected :: multi_infl =   1., &
                                               RTPS_Alpha = 0.85, &
                                               RTPP_Alpha = 0.85

    logical, dimension(max_vars), protected :: use_RTPS  = .false., &
                                               use_RTPP  = .false.

contains
    
    subroutine read_namelist(filename)
        implicit none
        character(len=*), intent(in)     :: filename

        namelist /control/  norain_value,     &
                            write_analy_mean, &
                            deterministic_update, &
                            NT2LOG,           &
                            NT2Dm,            &
                            NT2D0,            &
                            NT2De,            &
                            NT2D6,            &
                            wrf_mp_physics,   &
                            wrf_mp_hail_opt,  &
                            wrf_hypsometric_opt, &
                            nmember,          &
                            var_update

        namelist /observations/ radar_nml, &
                                synop_nml, &
                                ships_nml, &
                                metar_nml, &
                                sound_nml, &
                                gpspw_nml


        namelist /inflation/   multi_infl, use_RTPS, RTPS_Alpha, &
                                           use_RTPP, RTPP_Alpha

        !local
        integer       :: ierr


        open(20, file=filename, iostat=ierr)
        if(ierr /= 0) stop "No input.nml file"

        read(20, nml=control, iostat=ierr)
        if(ierr /= 0) stop "read namelist of control_nml fail!"

        read(20, nml=observations, iostat=ierr)
        if(ierr /= 0) stop "read namelist of observations_nml fail!"

        read(20, nml=inflation, iostat=ierr)
        if(ierr /= 0) stop "read namelist of inflation_nml fail!"

        close(20)

        if( nmember == -1 ) stop "Please input ensemble size in init_nml: nmember"

    end subroutine read_namelist

end module config
