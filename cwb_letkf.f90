program main
    use mpi_util,         only : letkf_init, letkf_finalize, timer, wait_jobs, &
                                 myid, nproc, req_ptr, is_root
    use param,            only : num_gts_indexes, num_radar_indexes, set_ensemble_constants
    use config,           only : nmember, read_namelist, write_analy_mean
    use grid,             only : grid_structure
    use eigen,            only : set_optimal_workspace_for_eigen, destroy_eigen_array
    use localization,     only : build_tree, destroy_tree
    use gts_omboma,       only : wrfda_gts
    use simulated_radar,  only : cwb_radar
    use letkf_core,       only : letkf_driver

    implicit none
    type(grid_structure)  :: wrf
    type(wrfda_gts)       :: gts
    type(cwb_radar)       :: rad
    logical, dimension(2) :: succeed
    character(len=3)      :: proc

    call letkf_init

    allocate( gts % platform(num_gts_indexes), &
              rad % radarobs(num_radar_indexes) )

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> reading namelist"
    call read_namelist("../input/input.nml")
    if(nproc < nmember) stop "Number of process should be greater than ensemble sizes!!"

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> set optimal workspace for eigen"
    call set_optimal_workspace_for_eigen(nmember)
    call set_ensemble_constants(nmember)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> define wrf mp_physics"
    call wrf % define_wrf_mp_physics()

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> read model data"
    if(myid < nmember) then 
        write(proc, '(i3.3)') myid+1

        call wrf % read_model("../input/wrfinput_nc_"//proc)
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> read obs data"
    if(nproc-myid-1 < nmember) then 
        write(proc, '(i3.3)') nmember-(nproc-myid-1)

        call gts % read_data ("../input/gts_letkf_"//proc, "../input/obs_gts")
        call rad % read_data ("../input/VR_letkf_"//proc, "VR") 
        call rad % read_data ("../input/MR_letkf_"//proc, "MR") 
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> exchanging obs data"
    call wrf % distribute_geoinfo(root=0)
    call gts % distribute(root=nproc-nmember)
    call rad % distribute(root=nproc-nmember)
    call wait_jobs(4, req_ptr%idx)

!   call rad % write_data(id=40)
!   call gts % write_data(id=40)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> building gts tree"
    succeed(1) = build_tree(gts % platform)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> building radar tree"
    succeed(2) = build_tree(rad % radarobs)

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> get into letkf core"
    if(any(succeed)) then
        call letkf_driver(wrf, gts, rad)
        call destroy_tree
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> finish letkf core"
    call destroy_eigen_array

    deallocate(gts % platform, &
               rad % radarobs)

    if(write_analy_mean) then
        if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> write analysis mean"
        call wrf % write_mean(filename="output/wrfout_nc_mean", root=nproc-1)
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> write analysis ensemble"
    if(myid < nmember) then 
        call wrf % write_model("output/wrfout_nc_"//proc)
    end if

    if(is_root) print '(f7.3,1x,a)', timer(), "sec ==========> finish all steps"
    call letkf_finalize
end program main
