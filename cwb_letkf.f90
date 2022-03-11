program main
    use mpi_util,         only : letkf_init, letkf_finalize, timer, &
                                 myid, nproc, is_root
    use param,            only : num_gts_indexes, num_radar_indexes, set_ensemble_constants
    use config,           only : nmember, read_namelist, write_analy_mean
    use grid,             only : grid_structure
    use eigen,            only : set_optimal_workspace_for_eigen, destroy_eigen_array
    use gts_omboma,       only : wrfda_gts
    use simulated_radar,  only : cwb_radar
    use letkf_core,       only : letkf_driver
    use projection,       only : proj_type

    implicit none
    type(grid_structure)  :: wrf
    type(wrfda_gts)       :: gts
    type(cwb_radar)       :: rad
    type(proj_type)       :: proj
    character(len=3)      :: proc

    call letkf_init

    allocate( gts % platform(num_gts_indexes), &
              rad % radarobs(num_radar_indexes) )

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> reading namelist"
    call read_namelist("../input/input.nml")

    if(nproc < nmember) stop "Number of process should be greater than ensemble sizes!!"

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> Initialize map factors"
    call proj % init()

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> Initialize MP constants and LAPACK working space"
    call wrf  % define_wrf_mp_physics()
    call set_optimal_workspace_for_eigen(nmember)
    call set_ensemble_constants(nmember)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> reading model data"
    if(myid < nmember) then 
        write(proc, '(i3.3)') myid+1

        call wrf % read_model("../input/wrfinput_nc_"//proc)
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> read obs data"
    if(nproc-myid-1 < nmember) then 
        write(proc, '(i3.3)') nmember-(nproc-myid-1)

        call gts % read_data (proj, "../input/gts_letkf_"//proc, "../input/obs_gts")
        call rad % read_data (proj, "../input/VR_letkf_"//proc, "VR") 
        call rad % read_data (proj, "../input/MR_letkf_"//proc, "MR") 
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> exchanging obs data"
    call wrf % distribute_geoinfo(root=0)
    call gts % distribute(root=nproc-nmember)
    call rad % distribute(root=nproc-nmember)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> get into letkf core"
    call letkf_driver(wrf, gts, rad, proj)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> finish letkf core"
    call destroy_eigen_array

    deallocate(gts % platform, &
               rad % radarobs)

    if(write_analy_mean) then
        if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> write analysis mean"
        call wrf % write_mean(filename="../output/wrfout_nc_mean", root=nproc-1)
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> write analysis ensemble"
    if(myid < nmember) then 
        write(proc, '(i3.3)') myid+1
        call wrf % write_model("../output/wrfout_nc_"//proc)
    end if


    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> finish all steps"
    call letkf_finalize
end program main
