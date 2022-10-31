!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module stats
  use param, only : mytype
  use decomp_2d, only : DECOMP_INFO
  implicit none

  character(len=*), parameter :: io_statistics = "statistics-io", &
       stat_dir = "statistics"

  integer :: stats_time

  real(mytype), dimension(:,:,:), allocatable ::  pmean
  real(mytype), dimension(:,:,:), allocatable ::  umean, uumean
  real(mytype), dimension(:,:,:), allocatable ::  vmean, vvmean
  real(mytype), dimension(:,:,:), allocatable ::  wmean, wwmean
  real(mytype), dimension(:,:,:), allocatable ::  uvmean, uwmean
  real(mytype), dimension(:,:,:), allocatable ::  vwmean
  real(mytype), dimension(:,:,:), allocatable ::  phimean, phiphimean

  real(mytype), dimension(:,:,:), allocatable ::  uuumean, uuvmean
  real(mytype), dimension(:,:,:), allocatable ::  uuwmean, uvvmean
  real(mytype), dimension(:,:,:), allocatable ::  uvwmean, uwwmean
  real(mytype), dimension(:,:,:), allocatable ::  vvwmean, vwwmean
  real(mytype), dimension(:,:,:), allocatable ::  vvvmean, wwwmean

  real(mytype), dimension(:,:,:), allocatable ::  pdudxmean, pdvdymean, pdwdzmean
  real(mytype), dimension(:,:,:), allocatable ::  pumean, pvmean, pwmean

  real(mytype), dimension(:,:,:), allocatable :: dudxmean, dudymean, dudzmean
  real(mytype), dimension(:,:,:), allocatable :: dvdxmean, dvdymean, dvdzmean
  real(mytype), dimension(:,:,:), allocatable :: dwdxmean, dwdymean, dwdzmean

  real(mytype), dimension(:,:,:), allocatable :: dudxdudxmean, dudxdvdxmean
  real(mytype), dimension(:,:,:), allocatable :: dudxdwdxmean, dvdxdvdxmean
  real(mytype), dimension(:,:,:), allocatable :: dvdxdwdxmean, dwdxdwdxmean

  real(mytype), dimension(:,:,:), allocatable :: dudydudymean, dudydvdymean
  real(mytype), dimension(:,:,:), allocatable :: dudydwdymean, dvdydvdymean
  real(mytype), dimension(:,:,:), allocatable :: dvdydwdymean, dwdydwdymean

  real(mytype), dimension(:,:,:), allocatable :: pdvdy_q1mean, pdvdy_q2mean
  real(mytype), dimension(:,:,:), allocatable :: pdvdy_q3mean, pdvdy_q4mean

  real(mytype), dimension(:,:,:), allocatable :: lambda2mean, lambda22mean

  
  real(mytype), dimension(:,:), allocatable :: coefy
  real(mytype) :: coefy_bc1(3), coefy_bcn(3)
  real(mytype) :: coefx_bc1(3), coefx_bcn(3)
  real(mytype) :: coefz_bc1(3), coefz_bcn(3)

  real(mytype) :: coefx_b, coefx_f
  real(mytype) :: coefz_b, coefz_f

  type(DECOMP_INFO) :: dstat_plane

  private
  public overall_statistic

contains

  subroutine init_statistic_adios2

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use param, only : istatbudget, istatpstrain, istatlambda2
    use var, only : numscalar
    
    implicit none

    integer :: ierror
    character(len=30) :: varname
    integer :: is
    logical, save :: initialised = .false.

    if (.not. initialised) then
       call decomp_2d_init_io(io_statistics)
    
       call decomp_2d_register_variable(io_statistics, "umean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "vmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "wmean", 3, 1, 3, mytype, opt_nplanes=1)
       
       call decomp_2d_register_variable(io_statistics, "pmean", 3, 1, 3, mytype, opt_nplanes=1)
       
       call decomp_2d_register_variable(io_statistics, "uumean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "vvmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "wwmean", 3, 1, 3, mytype, opt_nplanes=1)

       call decomp_2d_register_variable(io_statistics, "uvmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "uwmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "vwmean", 3, 1, 3, mytype, opt_nplanes=1)

       if (istatbudget) then
        call decomp_2d_register_variable(io_statistics, "uuumean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "uuvmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "uuwmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "uvvmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "uvwmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "uwwmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "vvvmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "vvwmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "vwwmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "wwwmean", 3, 1, 3, mytype, opt_nplanes=1)

        call decomp_2d_register_variable(io_statistics, "pdudxmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pdvdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pdwdzmean", 3, 1, 3, mytype, opt_nplanes=1)

        call decomp_2d_register_variable(io_statistics, "pumean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pvmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pwmean", 3, 1, 3, mytype, opt_nplanes=1)

        call decomp_2d_register_variable(io_statistics, "dudxmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dudymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dudzmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdxmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdzmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dwdxmean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dwdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dwdzmean", 3, 1, 3, mytype, opt_nplanes=1)

        call decomp_2d_register_variable(io_statistics, "dudydudymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dudydvdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dudydwdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdydvdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdydwdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dwdydwdymean", 3, 1, 3, mytype, opt_nplanes=1)

        call decomp_2d_register_variable(io_statistics, "dudydudymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dudydvdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dudydwdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdydvdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dvdydwdymean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "dwdydwdymean", 3, 1, 3, mytype, opt_nplanes=1)
       endif

       if (istatpstrain) then
        call decomp_2d_register_variable(io_statistics, "pdvdy_q1mean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pdvdy_q2mean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pdvdy_q3mean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "pdvdy_q4mean", 3, 1, 3, mytype, opt_nplanes=1)
       endif

       if (istatlambda2) then
        call decomp_2d_register_variable(io_statistics, "lambda2mean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "lambda22mean", 3, 1, 3, mytype, opt_nplanes=1)
       endif
       do is=1, numscalar
          write(varname,"('phi',I2.2)") is
          call decomp_2d_register_variable(io_statistics, varname, 3, 1, 3, mytype, opt_nplanes=1)
          write(varname,"('phiphi',I2.2)") is
          call decomp_2d_register_variable(io_statistics, varname, 3, 1, 3, mytype, opt_nplanes=1)
       enddo

       initialised = .true.
    endif
       
  end subroutine init_statistic_adios2
  
  !
  ! Initialize to zero all statistics
  !
  subroutine init_statistic

    use param, only : zero, iscalar, istatbudget, istatpstrain, istatlambda2
    use variables, only : nx, ny, nz
    use decomp_2d, only : zsize, decomp_info_init
    use MPI

    implicit none
    integer :: code
    allocate(umean(zsize(1),zsize(2),1))
    allocate(vmean(zsize(1),zsize(2),1))
    allocate(wmean(zsize(1),zsize(2),1))
    allocate(pmean(zsize(1),zsize(2),1))
    allocate(uumean(zsize(1),zsize(2),1))
    allocate(vvmean(zsize(1),zsize(2),1))
    allocate(wwmean(zsize(1),zsize(2),1))
    allocate(uvmean(zsize(1),zsize(2),1))
    allocate(uwmean(zsize(1),zsize(2),1))
    allocate(vwmean(zsize(1),zsize(2),1))

    pmean = zero
    umean = zero
    uumean = zero
    vmean = zero
    vvmean = zero
    wmean = zero
    wwmean = zero
    uvmean = zero
    uwmean = zero
    vwmean = zero

    if (istatbudget) then
      allocate(uuumean(zsize(1),zsize(2),1))
      allocate(uuvmean(zsize(1),zsize(2),1))
      allocate(uuwmean(zsize(1),zsize(2),1))
      allocate(uvvmean(zsize(1),zsize(2),1))
      allocate(uvwmean(zsize(1),zsize(2),1))
      allocate(uwwmean(zsize(1),zsize(2),1))
      allocate(vvvmean(zsize(1),zsize(2),1))
      allocate(vvwmean(zsize(1),zsize(2),1))
      allocate(vwwmean(zsize(1),zsize(2),1))
      allocate(wwwmean(zsize(1),zsize(2),1))

      allocate(pdudxmean(zsize(1),zsize(2),1))
      allocate(pdvdymean(zsize(1),zsize(2),1))
      allocate(pdwdzmean(zsize(1),zsize(2),1))

      allocate(phimean(zsize(1),zsize(2),1))
      allocate(phiphimean(zsize(1),zsize(2),1))

      allocate(pumean(zsize(1),zsize(2),1))
      allocate(pvmean(zsize(1),zsize(2),1))
      allocate(pwmean(zsize(1),zsize(2),1))

      allocate(dudxmean(zsize(1),zsize(2),1))
      allocate(dudymean(zsize(1),zsize(2),1))
      allocate(dudzmean(zsize(1),zsize(2),1))
      allocate(dvdxmean(zsize(1),zsize(2),1))
      allocate(dvdymean(zsize(1),zsize(2),1))
      allocate(dvdzmean(zsize(1),zsize(2),1))
      allocate(dwdxmean(zsize(1),zsize(2),1))
      allocate(dwdymean(zsize(1),zsize(2),1))
      allocate(dwdzmean(zsize(1),zsize(2),1))

      allocate(dudxdudxmean(zsize(1),zsize(2),1))
      allocate(dudxdvdxmean(zsize(1),zsize(2),1))
      allocate(dudxdwdxmean(zsize(1),zsize(2),1))
      allocate(dvdxdvdxmean(zsize(1),zsize(2),1))
      allocate(dvdxdwdxmean(zsize(1),zsize(2),1))
      allocate(dwdxdwdxmean(zsize(1),zsize(2),1))

      allocate(dudydudymean(zsize(1),zsize(2),1))
      allocate(dudydvdymean(zsize(1),zsize(2),1))
      allocate(dudydwdymean(zsize(1),zsize(2),1))
      allocate(dvdydvdymean(zsize(1),zsize(2),1))
      allocate(dvdydwdymean(zsize(1),zsize(2),1))
      allocate(dwdydwdymean(zsize(1),zsize(2),1))


      uuumean = zero
      uuvmean = zero
      uuwmean = zero
      uvvmean = zero
      uvwmean = zero
      uwwmean = zero
      vvvmean = zero
      vvwmean = zero
      vwwmean = zero
      wwwmean = zero

      pdudxmean = zero
      pdvdymean = zero
      pdwdzmean = zero

      pumean = zero
      pvmean = zero
      pwmean = zero

      dudxmean = zero
      dudymean = zero
      dudzmean = zero
      dvdxmean = zero
      dvdymean = zero
      dvdzmean = zero
      dwdxmean = zero
      dwdymean = zero
      dwdzmean = zero

      dudxdudxmean = zero
      dudxdvdxmean = zero
      dudxdwdxmean = zero
      dvdxdvdxmean = zero
      dvdxdwdxmean = zero
      dwdxdwdxmean = zero

      dudydudymean = zero
      dudydvdymean = zero
      dudydwdymean = zero
      dvdydvdymean = zero
      dvdydwdymean = zero
      dwdydwdymean = zero
    endif

    if (istatpstrain) then
      if (.not.istatbudget) then
        write(*,*) "If p-strain conditional averaging used, istatbudget must be T"
        call MPI_Abort(MPI_COMM_WORLD,1,code)
      endif

      allocate(pdvdy_q1mean(zsize(1),zsize(2),1))
      allocate(pdvdy_q2mean(zsize(1),zsize(2),1))
      allocate(pdvdy_q3mean(zsize(1),zsize(2),1))
      allocate(pdvdy_q4mean(zsize(1),zsize(2),1))

      pdvdy_q1mean = zero
      pdvdy_q2mean = zero
      pdvdy_q3mean = zero
      pdvdy_q4mean = zero
    endif

    if (istatlambda2) then
#ifndef HAVE_LAPACK 
      write(*,*) "Cannot calaculate mean lambda2 without LAPACK"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
#endif  

      allocate(lambda2mean(zsize(1),zsize(2),1))
      allocate(lambda22mean(zsize(1),zsize(2),1))

      lambda2mean = zero
      lambda22mean = zero
    endif
  
    call decomp_info_init(nx, ny, 1,dstat_plane)

    call init_statistic_adios2

    call grad_init

  end subroutine init_statistic

  !
  ! Read all statistics if possible
  !
  subroutine restart_statistic

    use param, only : initstat, irestart, ifirst, zero
    use variables, only : nstat
    use var, only : tmean

    implicit none

    ! No reading for statistics when nstat > 1 or no restart
    call init_statistic()

    if (nstat.gt.1 .or. irestart.eq.0) then
       initstat = ifirst
       return
    else
       call init_statistic_adios2
    endif


    ! Temporary array
    tmean = zero

    ! Read all statistics
    call read_or_write_all_stats(.true.)

  end subroutine restart_statistic

  function gen_statname(stat) result(newname)

    implicit none
    
    character(len=*), intent(in) :: stat
    character(len=30) :: newname
    
#ifndef ADIOS2
    write(newname, "(A,'.dat',I7.7)") stat, stats_time
#else
    write(newname, *) stat
#endif
    
  end function gen_statname

  !
  ! Statistics: perform all IO
  !
  subroutine read_or_write_all_stats(flag_read)

    use param, only : iscalar, itime, istatbudget, istatpstrain, initstat2, istatlambda2
    use variables, only : numscalar
    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_write_mode, decomp_2d_read_mode, &
         decomp_2d_open_io, decomp_2d_close_io, decomp_2d_start_io, decomp_2d_end_io

    implicit none

    ! Argument
    logical, intent(in) :: flag_read

    ! Local variables
    integer :: is, it
    character(len=30) :: filename
    integer :: io_mode
    integer :: ierror
    

    ! File ID to read or write
    if (flag_read) then
        it = itime - 1
    else
        it = itime
     endif
     stats_time = it

    if (nrank==0) then
      print *,'==========================================================='
      if (flag_read) then
        print *,'Reading stat file', stats_time
      else
        print *,'Writing stat file', stats_time
      endif
    endif

    if (flag_read) then
       io_mode = decomp_2d_read_mode
    else
       io_mode = decomp_2d_write_mode
    endif
#ifdef ADIOS2
    call decomp_2d_open_io(io_statistics, stat_dir, io_mode)
    call decomp_2d_start_io(io_statistics, stat_dir)
#endif
    
    call read_or_write_one_stat(flag_read, gen_statname("pmean"), pmean)
    call read_or_write_one_stat(flag_read, gen_statname("umean"), umean)
    call read_or_write_one_stat(flag_read, gen_statname("vmean"), vmean)
    call read_or_write_one_stat(flag_read, gen_statname("wmean"), wmean)

    call read_or_write_one_stat(flag_read, gen_statname("uumean"), uumean)
    call read_or_write_one_stat(flag_read, gen_statname("vvmean"), vvmean)
    call read_or_write_one_stat(flag_read, gen_statname("wwmean"), wwmean)

    call read_or_write_one_stat(flag_read, gen_statname("uvmean"), uvmean)
    call read_or_write_one_stat(flag_read, gen_statname("uwmean"), uwmean)
    call read_or_write_one_stat(flag_read, gen_statname("vwmean"), vwmean)

    if (istatbudget) then
      call read_or_write_one_stat(flag_read, gen_statname("uuumean"), uuumean)
      call read_or_write_one_stat(flag_read, gen_statname("uuvmean"), uuvmean)
      call read_or_write_one_stat(flag_read, gen_statname("uuwmean"), uuwmean)
      call read_or_write_one_stat(flag_read, gen_statname("uvvmean"), uvvmean)
      call read_or_write_one_stat(flag_read, gen_statname("uvwmean"), uvwmean)
      call read_or_write_one_stat(flag_read, gen_statname("uwwmean"), uwwmean)
      call read_or_write_one_stat(flag_read, gen_statname("vvvmean"), vvvmean)
      call read_or_write_one_stat(flag_read, gen_statname("vvwmean"), vvwmean)
      call read_or_write_one_stat(flag_read, gen_statname("vwwmean"), vwwmean)
      call read_or_write_one_stat(flag_read, gen_statname("wwwmean"), wwwmean)
      
      call read_or_write_one_stat(flag_read, gen_statname("pdudxmean"), pdudxmean)
      call read_or_write_one_stat(flag_read, gen_statname("pdvdymean"), pdvdymean)
      call read_or_write_one_stat(flag_read, gen_statname("pdwdzmean"), pdwdzmean)

      call read_or_write_one_stat(flag_read, gen_statname("pumean"), pumean)
      call read_or_write_one_stat(flag_read, gen_statname("pvmean"), pvmean)
      call read_or_write_one_stat(flag_read, gen_statname("pwmean"), pwmean)

      call read_or_write_one_stat(flag_read, gen_statname("dudxmean"), dudxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dudymean"), dudymean)
      call read_or_write_one_stat(flag_read, gen_statname("dudzmean"), dudzmean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdxmean"), dvdxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdymean"), dvdymean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdzmean"), dvdzmean)
      call read_or_write_one_stat(flag_read, gen_statname("dwdxmean"), dwdxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dwdymean"), dwdymean)
      call read_or_write_one_stat(flag_read, gen_statname("dwdzmean"), dwdzmean)

      call read_or_write_one_stat(flag_read, gen_statname("dudxdudxmean"), dudxdudxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dudxdvdxmean"), dudxdvdxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dudxdwdxmean"), dudxdwdxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdxdvdxmean"), dvdxdvdxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdxdwdxmean"), dvdxdwdxmean)
      call read_or_write_one_stat(flag_read, gen_statname("dwdxdwdxmean"), dwdxdwdxmean)

      call read_or_write_one_stat(flag_read, gen_statname("dudydudymean"), dudydudymean)
      call read_or_write_one_stat(flag_read, gen_statname("dudydvdymean"), dudydvdymean)
      call read_or_write_one_stat(flag_read, gen_statname("dudydwdymean"), dudydwdymean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdydvdymean"), dvdydvdymean)
      call read_or_write_one_stat(flag_read, gen_statname("dvdydwdymean"), dvdydwdymean)
      call read_or_write_one_stat(flag_read, gen_statname("dwdydwdymean"), dwdydwdymean)
    endif

    if (istatpstrain .and. itime > initstat2) then
      call read_or_write_one_stat(flag_read, gen_statname("pdvdy_q1mean"), pdvdy_q1mean)
      call read_or_write_one_stat(flag_read, gen_statname("pdvdy_q2mean"), pdvdy_q2mean)
      call read_or_write_one_stat(flag_read, gen_statname("pdvdy_q3mean"), pdvdy_q3mean)
      call read_or_write_one_stat(flag_read, gen_statname("pdvdy_q4mean"), pdvdy_q4mean)
    endif

    if (istatlambda2 .and. itime > initstat2) then
      call read_or_write_one_stat(flag_read, gen_statname("lambda2mean"), lambda2mean)
      call read_or_write_one_stat(flag_read, gen_statname("lambda22mean"), lambda22mean)
    endif
#ifdef ADIOS2
    call decomp_2d_end_io(io_statistics, stat_dir)
    call decomp_2d_close_io(io_statistics, stat_dir)
#endif
    
    if (nrank==0) then
      if (flag_read) then
        print *,'Read stat done!'
      else
        print *,'Write stat done!'
      endif
      print *,'==========================================================='
    endif

  end subroutine read_or_write_all_stats

  !
  ! Statistics: perform one IO
  !
  subroutine read_or_write_one_stat(flag_read, filename, array)

    use decomp_2d, only : mytype, xstS, xenS, zsize
    use decomp_2d_io, only : decomp_2d_read_one, decomp_2d_write_one

    implicit none

    ! Arguments
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: array
    integer :: ierror

    if (flag_read) then
       ! There was a check for nvisu = 1 before
       call decomp_2d_read_one(3, array, stat_dir, filename, io_statistics, reduce_prec=.false.,opt_decomp=dstat_plane)
    else
       call decomp_2d_write_one(3, array, stat_dir, filename, 0, io_statistics, reduce_prec=.false.,opt_decomp=dstat_plane)
    endif

  end subroutine read_or_write_one_stat

  !
  ! Statistics : Intialize, update and perform IO
  !
  subroutine overall_statistic(ux1,uy1,uz1,phi1,pp3,ep1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io
    use tools, only : rescale_pressure
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use ibm_param, only : ubcx, ubcy, ubcz
    use var, only : ux2, ux3, uy2, uy3, uz2, uz3
    use var, only : ta2
    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ppi3, dip3
    use var, only : pp2, ppi2, dip2
    use var, only : pp1, ta1, di1
    implicit none

    !! Inputs
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar),intent(in) :: phi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    !! Locals

    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dudy, dvdy, dwdy
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dudx, dvdx, dwdx

    integer :: is
    character(len=30) :: filename

    if (itime.lt.initstat) then
       return
    elseif (itime.eq.initstat) then
       call init_statistic()

       if (istatout < 1) istatout = icheckpoint
    elseif (itime.eq.ifirst) then
       if (itempaccel.eq.1) then
          call init_statistic()
       else
          call restart_statistic()
       endif
       if (istatout < 1) istatout = icheckpoint
    endif
    
    if (mod(itime,istatcalc) /= 0) return
    if (itempaccel == 1 .and. mod(itime,istatout)/=0) return 

    !! Mean pressure
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(td1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    ! Convert to physical pressure
    call rescale_pressure(td1)

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_x_to_y(td1,td2)

    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
    call transpose_y_to_z(td2,td3)

    call update_average_scalar(pmean, td3, ep1)

    !! Mean velocity
    call update_average_vector(umean, vmean, wmean, &
                               ux3, uy3, uz3, td3)

    !! Second-order velocity moments
    call update_variance_vector(uumean, vvmean, wwmean, uvmean, uwmean, vwmean, &
                                ux3, uy3, uz3, td3)

    if (istatbudget) then
      call grad_x(ux1, ta1)
      call grad_x(uy1, tb1)
      call grad_x(uz1, tc1)
      
      !y-derivatives
      call grad_y(ux2,ta2)
      call grad_y(uy2,tb2)
      call grad_y(uz2,tc2)

      !!z-derivatives

      call grad_z(ux3,ta3)
      call grad_z(uy3,tb3)
      call grad_z(uz3,tc3)

      call transpose_y_to_z(ta2,dudy)
      call transpose_y_to_z(tb2,dvdy)
      call transpose_y_to_z(tc2,dwdy)

      call transpose_x_to_y(ta1,ta2)
      call transpose_x_to_y(tb1,tb2)
      call transpose_x_to_y(tc1,tc2)

      call transpose_y_to_z(ta2,dudx)
      call transpose_y_to_z(tb2,dvdx)
      call transpose_y_to_z(tc2,dwdx)

    call update_skewness_tensor(uuumean,uuvmean,uuwmean,uvvmean,uvwmean,uwwmean,&
                               vvvmean,vvwmean,vwwmean,wwwmean,&
                               ux3,uy3,uz3,td3)  

    call update_pvelograd_vector(pdudxmean,pdvdymean,pdwdzmean,dudx,dvdy,tc3,td3)

    call update_pvelo_vector(pumean, pvmean, pwmean, ux3, uy3, uz3,td3)

    call update_pvelograd_vector(pdudxmean,pdvdymean,pdwdzmean,dudx,dvdy,tc3,td3)

    call update_velograd_tensor(dudxmean,dudymean,dudzmean,dvdxmean,dvdymean,&
                                dvdzmean,dwdxmean,dwdymean,dwdzmean,&
                                dudx,dudy,ta3,dvdx,dvdy,tb3,dwdx,dwdy,tc3,td3)

    call update_velograd2_tensor(dudxdudxmean,dudxdvdxmean,dudxdwdxmean,&
                                 dvdxdvdxmean, dvdxdwdxmean,dwdxdwdxmean,&
                                 dudx,dvdx,dwdx,td3)

    call update_velograd2_tensor(dudydudymean,dudydvdymean,dudydwdymean,&
                                 dvdydvdymean, dvdydwdymean,dwdydwdymean,&
                                 dudy,dvdy,dwdy,td3)                         
    endif                             

    if (istatpstrain .and. itime>initstat2) then
      call update_pstrain_cond_avg(pdvdy_q1mean,pdvdy_q2mean,&
                                  pdvdy_q3mean,pdvdy_q4mean,&
                                  td3,dvdy,pmean,dvdymean,td3)
    endif

    if (istatlambda2) then
      call update_lambda2_avg(lambda2mean,lambda22mean,&
                                dudx,dudy,ta3,dvdx,dvdy,&
                                tb3,dwdx,dwdy,tc3,td3)
    endif

    ! Write all statistics
    if (mod(itime,istatout)==0) then
       call read_or_write_all_stats(.false.)
    endif

  end subroutine overall_statistic

  !
  ! Basic function, can be applied to arrays
  !
  elemental real(mytype) function one_minus_ep1(var, ep1)

    use decomp_2d, only : mytype
    use param, only : iibm, one

    implicit none

    ! inputs
    real(mytype), intent(in) :: var, ep1

    if (iibm==2) then
      one_minus_ep1 = (one - ep1) * var
    else
      one_minus_ep1 = var
    endif

  end function one_minus_ep1

  !
  ! Update um, the average of ux
  !
  subroutine update_average_scalar(um, ux, ep,mask,istat2)

    use decomp_2d, only : mytype, xsize, zsize
    use param, only : itime, initstat, istatcalc, itempaccel, initstat2
    use var, only : di1, tmean

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: um
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux
    logical, dimension(zsize(1),zsize(2),zsize(3)), optional, intent(in) :: mask
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ep
    logical, optional, intent(in) :: istat2

    real(mytype), dimension(zsize(1), zsize(2)) :: stat_z
    real(mytype) :: stat_inc
    integer :: i, j
    if (present(mask)) then
      stat_z = sum(ux,dim=3,mask=mask) / real(zsize(3),kind=mytype)

    else
      stat_z = sum(ux,dim=3) / real(zsize(3),kind=mytype)
    endif

    if (itempaccel==1) then
      um(:,:,1) = stat_z(:,:)
      return
    endif

    if (present(istat2)) then
      if (istat2) then
        stat_inc = real((itime-initstat2)/istatcalc+1, kind=mytype)
      else
        stat_inc = real((itime-initstat)/istatcalc+1, kind=mytype)
      endif
    else
      stat_inc = real((itime-initstat)/istatcalc+1, kind=mytype)
    endif

    do j = 1, zsize(2)
      do i = 1, zsize(1)
        um(i,j,1) = um(i,j,1) + (stat_z(i,j) - um(i,j,1))/ stat_inc
      enddo
    enddo

  end subroutine update_average_scalar

  ! subroutine average_z(array, stat_xz,mask)
  !   use decomp_2d
  !   use param, only : zero

  !   real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: array
  !   real(mytype), dimension(zsize(1),zsize(2),1) :: stat_xz
  !   logical, dimension(zsize(1),zsize(2),zsize(3)), optional, intent(in) :: mask
  !   real(mytype) :: avg_z
  !   integer :: i, j, k

    
  !   stat_xz = sum(array,dim=3,mask=mask) / real(zsize(3),kind=mytype)


  !   do i = 1, zsize(1)
  !     do j = 1, zsize(2)
  !       avg_z = sum(array)
  !       stat_xz(i,j,1) = stat_xz(i,j,1) + array(i,j,k)/zsize(3)
  !     enddo
  !   enddo


  !  end subroutine
  !
  ! Update (um, vm, wm), the average of (ux, uy, uz)
  !
  subroutine update_average_vector(um, vm, wm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, zsize

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2), 1), intent(inout) :: um, vm, wm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(um, ux, ep)
    call update_average_scalar(vm, uy, ep)
    call update_average_scalar(wm, uz, ep)

  end subroutine update_average_vector

  !
  ! Update (uum, vvm, wwm, uvm, uwm, vwm), the variance of the vector (ux, uy, uz)
  !
  subroutine update_variance_vector(uum, vvm, wwm, uvm, uwm, vwm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, zsize

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: uum, vvm, wwm, uvm, uwm, vwm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(uum, ux*ux, ep)
    call update_average_scalar(vvm, uy*uy, ep)
    call update_average_scalar(wwm, uz*uz, ep)
    call update_average_scalar(uvm, ux*uy, ep)
    call update_average_scalar(uwm, ux*uz, ep)
    call update_average_scalar(vwm, uy*uz, ep)

  end subroutine update_variance_vector

  subroutine update_skewness_tensor(uuum,uuvm,uuwm,uvvm,uvwm,&
                                    uwwm,vvvm,vvwm,vwwm,wwwm,&
                                    ux, uy, uz, ep)
  use decomp_2d, only : mytype, xsize, zsize
  real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: uuum,uuvm,uuwm,uvvm,uvwm,&
                                                                 uwwm,vvvm,vvwm,vwwm,wwwm
  real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

  call update_average_scalar(uuum, ux*ux*ux, ep)
  call update_average_scalar(uuvm, ux*ux*uy, ep)
  call update_average_scalar(uuwm, ux*ux*uz, ep)
  call update_average_scalar(uvvm, ux*uy*uy, ep)
  call update_average_scalar(uvwm, ux*uy*uz, ep)
  call update_average_scalar(uwwm, ux*uz*uz, ep)
  call update_average_scalar(vvvm, uy*uy*uy, ep)
  call update_average_scalar(vvwm, uy*uy*uz, ep)
  call update_average_scalar(vwwm, uy*uz*uz, ep)
  call update_average_scalar(wwwm, uz*uz*uz, ep)
                                    
  end subroutine update_skewness_tensor

  subroutine update_pvelograd_vector(pdudxm,pdvdym,pdwdzm,dudx,dvdy,dwdz,ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: pdudxm,pdvdym,pdwdzm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: dudx, dvdy, dwdz, ep

    call update_average_scalar(pdudxm, ep*dudx, ep)
    call update_average_scalar(pdvdym, ep*dvdy, ep)
    call update_average_scalar(pdwdzm, ep*dwdz, ep)
  end subroutine

  subroutine update_pvelo_vector(pum, pvm, pwm, ux, uy, uz, ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: pum, pvm, pwm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(pum, ep*ux, ep)
    call update_average_scalar(pvm, ep*uy, ep)
    call update_average_scalar(pwm, ep*uz, ep)

  end subroutine update_pvelo_vector

  subroutine update_velograd_tensor(dudxm, dudym, dudzm, dvdxm,&
                                    dvdym, dvdzm, dwdxm, dwdym,&
                                    dwdzm, dudx, dudy, dudz, dvdx,&
                                    dvdy, dvdz, dwdx, dwdy, dwdz, ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: dudxm, dudym, dudzm
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: dvdxm, dvdym, dvdzm
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: dwdxm, dwdym, dwdzm


    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dudy, dudz, ep
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dvdx, dvdy, dvdz   
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dwdx, dwdy, dwdz

    call update_average_scalar(dudxm, dudx, ep)
    call update_average_scalar(dudym, dudy, ep)
    call update_average_scalar(dudzm, dudz, ep)
    call update_average_scalar(dvdxm, dvdx, ep)
    call update_average_scalar(dvdym, dvdy, ep)
    call update_average_scalar(dvdzm, dvdz, ep)
    call update_average_scalar(dwdxm, dudx, ep)
    call update_average_scalar(dwdym, dwdy, ep)
    call update_average_scalar(dwdzm, dwdz, ep)


  end subroutine update_velograd_tensor

  subroutine update_velograd2_tensor(dudxdudxm, dudxdvdxm, dudxdwdxm,&
                                     dvdxdvdxm, dvdxdwdxm, dwdxdwdxm,&
                                     dudx,dvdx,dwdx,ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: dudxdudxm, dudxdvdxm
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: dudxdwdxm, dvdxdvdxm
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: dvdxdwdxm, dwdxdwdxm

    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dvdx, dwdx, ep

    call update_average_scalar(dudxdudxm, dudx*dudx, ep)
    call update_average_scalar(dudxdvdxm, dudx*dvdx, ep)
    call update_average_scalar(dudxdwdxm, dudx*dwdx, ep)
    call update_average_scalar(dvdxdvdxm, dvdx*dvdx, ep)
    call update_average_scalar(dvdxdwdxm, dvdx*dwdx, ep)
    call update_average_scalar(dwdxdwdxm, dwdx*dwdx, ep)

  end subroutine update_velograd2_tensor

  subroutine update_pstrain_cond_avg(pdvdy_q1m,pdvdy_q2m,&
                                     pdvdy_q3m,pdvdy_q4m,&
                                     p,dvdy, pm, dvdym,ep)
    use decomp_2d, only : mytype, xsize, zsize
    use param, only : zero, zpfive
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: pdvdy_q1m, pdvdy_q2m
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: pdvdy_q3m, pdvdy_q4m
    real(mytype), dimension(zsize(1),zsize(2),1), intent(in) :: pm, dvdym
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: p,dvdy, ep

    logical, dimension(zsize(1),zsize(2),zsize(3)) :: mask1, mask2, mask3, mask4
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: p_fluct, dvdy_fluct
    integer :: i,j,k

    do k = 1, zsize(3)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          p_fluct(i,j,k) = p(i,j,k) - pm(i,j,1)
          dvdy_fluct(i,j,k) = dvdy(i,j,k) - dvdym(i,j,1)
          mask1(i,j,k) = p_fluct(i,j,k) .gt. zero .and.  dvdy_fluct(i,j,k) .gt. zero
          mask2(i,j,k) = p_fluct(i,j,k) .lt. zero .and.  dvdy_fluct(i,j,k) .gt. zero
          mask3(i,j,k) = p_fluct(i,j,k) .lt. zero .and.  dvdy_fluct(i,j,k) .lt. zero
          mask4(i,j,k) = p_fluct(i,j,k) .gt. zero .and.  dvdy_fluct(i,j,k) .lt. zero
        enddo
      enddo
    enddo

    call update_average_scalar(pdvdy_q1m, p_fluct*dvdy_fluct, ep, mask=mask1,istat2=.true.)
    call update_average_scalar(pdvdy_q2m, p_fluct*dvdy_fluct, ep, mask=mask2,istat2=.true.)
    call update_average_scalar(pdvdy_q3m, p_fluct*dvdy_fluct, ep, mask=mask3,istat2=.true.)
    call update_average_scalar(pdvdy_q4m, p_fluct*dvdy_fluct, ep, mask=mask4,istat2=.true.)

  end subroutine update_pstrain_cond_avg

  subroutine update_lambda2_avg(lambda2m,lambda22m,&
                              dudx,dudy,dudz,dvdx,dvdy,dvdz,&
                              dwdx,dwdy,dwdz,ep)
    use decomp_2d, only : mytype, xsize, zsize
    use MPI
    use param

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: lambda2m,lambda22m
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dudy, dudz
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dvdx, dvdy, dvdz
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dwdx, dwdy, dwdz, ep

    real(mytype), dimension(3,3) :: S, Omega, S2O2, U_mat
    integer :: lwork, info, info_sum, nb
    real(mytype), dimension(3) :: lambda
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: lambda2
    real(mytype), dimension(:), allocatable :: work
    integer :: i,j,k,code

#ifdef HAVE_LAPACK
    interface
      integer function ilaenv(ISPEC,NAME,OPTS, N1, N2, N3, N4)
        integer :: ispec
        character(*) :: NAME, OPTS
        integer :: N1, N2, N3, N4
      end function ilaenv
    end interface

    nb = ilaenv(1,'DSYTRD','NU',3,3,-1,-1)
    lwork = (nb+2)*3
    allocate(work(lwork))

    info_sum = 0
    do k = 1, zsize(3)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          S(1,1) = dudx(i,j,k)
          S(2,2) = dvdy(i,j,k)
          S(3,3) = dwdz(i,j,k)

          S(1,2) = zpfive* ( dudy(i,j,k) + dvdx(i,j,k))
          S(2,1) = S(1,2)

          S(1,3) = zpfive* ( dudz(i,j,k) + dwdx(i,j,k))
          S(3,1) = S(1,3)

          S(3,2) = zpfive* ( dwdy(i,j,k) + dvdz(i,j,k))
          S(2,3) = S(3,2)

          Omega(1,1) = zero
          Omega(2,2) = zero
          Omega(3,3) = zero

          Omega(1,2) = zpfive* ( dudy(i,j,k) - dvdx(i,j,k))
          Omega(2,1) = -Omega(1,2)

          Omega(1,3) = zpfive* ( dudz(i,j,k) - dwdx(i,j,k))
          Omega(3,1) = -Omega(1,3)

          Omega(3,2) = zpfive* ( dwdy(i,j,k) - dvdz(i,j,k))
          Omega(2,3) = -Omega(3,2)

          S2O2 = matmul(S,S) + matmul(Omega,Omega)

          call dsyev('N','U',3,S2O2,U_mat,lambda,work,lwork,info)
          lambda2(i,j,k) = lambda(2)
          info_sum = info_sum + info
        enddo
      enddo
    enddo

    if (info_sum/=0) then
      write(*,*) "Problems calculating eigenvalules"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif

    call update_average_scalar(lambda2m, lambda2, ep)
    call update_average_scalar(lambda22m, lambda2*lambda2, ep)

#endif    
  end subroutine                                                         
  subroutine grad_init
    use param, only : zpfive, one, three,onepfive, two
    use var, only : yp, ny, dx, dz

    real(mytype) :: ddy1, ddy2
    integer :: i

    allocate(coefy(2:ny-1,3))

    do i = 2, ny-1
      ddy1 = yp(i) - yp(i-1)
      ddy2 = yp(i+1) - yp(i)

      coefy(i,1) = -(ddy2)/(ddy1 * (ddy1 + ddy2))
      coefy(i,2) = (ddy2 - ddy1) / (ddy1 * ddy2)
      coefy(i,3) = ddy1 / (ddy2 * (ddy1 + ddy2))
    enddo

    coefx_b = - one / dx
    coefx_f = one /dx

    coefz_b = - one/dz
    coefz_f = one/dz

    ddy1 = yp(2) - yp(1)
    ddy2 = yp(3) - yp(2)

    coefy_bc1(1) = -( two*ddy1 + ddy2)/(ddy1*(ddy1 + ddy2))
    coefy_bc1(2) = (ddy1 + ddy2)/(ddy1*ddy2)
    coefy_bc1(3) = -ddy1/(ddy2*(ddy1 + ddy2))

    ddy1 = yp(ny-1) - yp(ny-2)
    ddy2 = yp(ny-2) - yp(ny-3)

    coefy_bcn(1) = ddy2/(ddy1*(ddy1 + ddy2))
    coefy_bcn(2) = -(ddy1 + ddy2)/(ddy1*ddy2)
    coefy_bcn(3) = (two*ddy2 + ddy1)/(ddy2*(ddy1 + ddy2))

    coefx_bc1(1) = - onepfive/ dx
    coefx_bc1(2) = two / dx
    coefx_bc1(3) = - zpfive / dx

    coefx_bcn(1) = zpfive / dx
    coefx_bcn(2) = -two / dx
    coefx_bcn(3) = onepfive / dx

    coefz_bc1(1) = - onepfive/ dz
    coefz_bc1(2) = two / dz
    coefz_bc1(3) = - zpfive / dz

    coefz_bcn(1) = zpfive / dz
    coefz_bcn(2) = -two / dz
    coefz_bcn(3) = onepfive / dz
  end subroutine

  subroutine grad_x(array, grad)
    use decomp_2d, only : xsize
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: array
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: grad

    integer :: i, j, k


    do k = 1, xsize(3)
      do j = 1, xsize(2)

        grad(1,j,k) = coefx_bc1(1)*array(1,j,k) + coefx_bc1(2)*array(2,j,k) &
                       + coefx_bc1(3)*array(3,j,k)

        do i = 2,xsize(1)-1
          grad(i,j,k) = coefx_b*array(i-1,j,k) + coefx_f*array(i+1,j,k)
        enddo

        grad(xsize(1),j,k) = coefx_bcn(1)*array(xsize(1)-2,j,k) + &
                      coefx_bcn(2)*array(xsize(1)-1,j,k) + &
                      coefx_bcn(3)*array(xsize(1),j,k)
      enddo
    enddo

  end subroutine grad_x

  subroutine grad_y(array, grad)
    use decomp_2d, only : ysize
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)), intent(in) :: array
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)), intent(out) :: grad
    integer :: i, j, k

    do k = 1, ysize(3)
      do j = 2, ysize(2)-1
        do i = 1,ysize(1)
          grad(i,j,k) = coefy(j,1)*array(i,j-1,k) &
                        + coefy(j,2)*array(i,j,k) &
                        + coefy(j,3)*array(i,j+1,k)
        enddo
      enddo

      do i = 1, ysize(1)
        grad(i,1,k) = coefy_bc1(1)*array(i,1,k) + coefy_bc1(2)*array(i,2,k) &
                        + coefy_bc1(3)*array(i,3,k)

        grad(i,ysize(2),k) = coefy_bcn(1)*array(i,ysize(2)-2,k) + &
                              coefy_bcn(2)*array(i,ysize(2)-1,k) + &
                              coefy_bcn(3)*array(i,ysize(2),k)
      enddo
    enddo
  end subroutine grad_y

  subroutine grad_z(array, grad)
    use decomp_2d, only : zsize
    
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: array
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(out) :: grad
    integer :: i, j, k

    do k = 2, zsize(3)-1
      do j = 1, zsize(2)
        do i = 1,zsize(1)
          grad(i,j,k) = coefz_b*array(i,j,k-1) + coefz_f*array(i,j,k+1)
        enddo
      enddo
    enddo

    do j = 1, zsize(2)
      do i = 1,zsize(1)
        grad(i,j,1) = coefz_bc1(1)*array(i,j,1) + coefz_bc1(2)*array(i,j,2) &
                        + coefz_bc1(3)*array(i,j,3)       

        grad(i,j,zsize(3)) = coefz_bcn(1)*array(i,j,zsize(3)-2) + &
                        coefz_bcn(2)*array(i,j,zsize(3)-1) + &
                        coefz_bcn(3)*array(i,j,zsize(3))
      enddo
    enddo
  end subroutine grad_z
endmodule stats

