!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module stats
  use param, only : mytype
  use decomp_2d, only : DECOMP_INFO
  implicit none

#ifdef HAVE_FFTW
  include "fftw3.f"
#endif

  character(len=*), parameter :: io_statistics = "statistics-io", &
       stat_dir = "statistics"

  integer :: stats_time
  real(mytype), dimension(:,:,:), allocatable ::  uvwp_mean
  real(mytype), dimension(:,:,:), allocatable ::  uu_mean
  real(mytype), dimension(:,:,:), allocatable ::  pu_mean
  real(mytype), dimension(:,:,:), allocatable ::  uuu_mean
  real(mytype), dimension(:,:,:), allocatable ::  pdudx_mean
  real(mytype), dimension(:,:,:), allocatable ::  dudx_mean
  real(mytype), dimension(:,:,:), allocatable ::  dudxdudx_mean

  real(mytype), dimension(:,:,:), allocatable ::  pdvdy_q_mean

  real(mytype), dimension(:,:,:), allocatable :: uv_quadrant_mean
  real(mytype), dimension(:,:,:), allocatable :: uuuu_mean
  real(mytype), dimension(:,:,:), allocatable :: autocorr_mean

  complex(mytype), dimension(:,:,:,:), allocatable :: spectra_2d_mean
  complex(mytype), dimension(:,:,:,:), allocatable :: spectra_z_z_mean

  real(mytype), dimension(:,:,:), allocatable :: lambda2

  type(DECOMP_INFO) :: uvwp_info, uu_info, pu_info
  type(DECOMP_INFO) :: uuu_info, pdudx_info, dudx_info
  type(DECOMP_INFO) :: dudxdudx_info, pdvdy_q_info
  type(DECOMP_INFO) :: uv_quadrant_info, uuuu_info
  real(mytype), dimension(:,:,:), allocatable :: lambda2mean, lambda22mean

  
  real(mytype), dimension(:,:), allocatable :: coefy
  real(mytype) :: coefy_bc1(3), coefy_bcn(3)
  real(mytype) :: coefx_bc1(3), coefx_bcn(3)
  real(mytype) :: coefz_bc1(3), coefz_bcn(3)

  real(mytype) :: coefx_b, coefx_f
  real(mytype) :: coefz_b, coefz_f
  real(mytype), dimension(:), allocatable :: h_quads
  integer :: spectra_level, nspectra, spectra_nlocs
  real(mytype), dimension(:), allocatable :: spectra_xlocs
  integer, allocatable, dimension(:) :: spectra_x_indices
  type(DECOMP_INFO) :: dft_info
  type(DECOMP_INFO) :: dstat_plane
  integer(8) :: plan_x, plan_z
  integer :: zdft_size,xdft_size

  real(mytype) :: autocorr_xlocs
  real(mytype) :: autocorr_max_sep
  integer, dimension(:), allocatable :: autocorr_x_indices
  integer :: spectra_corr_nlocs, spectra_corr_thisrank, spectra_corr_comm
  real(mytype), dimension(:), allocatable :: spectra_corr_ylocs
  integer, dimension(:), allocatable :: spectra_corr_yinds, spectra_corr_ranks

  private
  public :: overall_statistic, h_quads, spectra_level, lambda2,&
            autocorr_max_sep, autocorr_xlocs, spectra_corr_nlocs, spectra_corr_ylocs

contains

  subroutine init_statistic_adios2

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use param, only : istatbudget, istatpstrain, istatlambda2, istatquadrant
    use param, only : istatflatness
    use var, only : numscalar
    
    implicit none

    integer :: ierror
    character(len=30) :: varname
    integer :: is
    logical, save :: initialised = .false.

    if (.not. initialised) then
       call decomp_2d_init_io(io_statistics)
    
       call decomp_2d_register_variable(io_statistics, "uvwp_mean", 3, 0, 0, mytype,opt_decomp=uvwp_info)
       call decomp_2d_register_variable(io_statistics, "uu_mean", 3, 0, 0, mytype,opt_decomp=uu_info)
              
       if (istatbudget) then
        call decomp_2d_register_variable(io_statistics, "pu_mean", 3, 0, 0, mytype,opt_decomp=pu_info)
        call decomp_2d_register_variable(io_statistics, "uuu_mean", 3, 0, 0, mytype,opt_decomp=uuu_info)
        call decomp_2d_register_variable(io_statistics, "pdudx_mean", 3, 0, 0, mytype,opt_decomp=pdudx_info)
        call decomp_2d_register_variable(io_statistics, "dudx_mean", 3, 0, 0, mytype,opt_decomp=pdudx_info)
        call decomp_2d_register_variable(io_statistics, "dudxdudx_mean", 3, 0, 0, mytype,opt_decomp=dudxdudx_info)

      endif

       if (istatpstrain) then
        call decomp_2d_register_variable(io_statistics, "pdvdy_q_mean", 3, 0, 0, mytype, opt_decomp=pdvdy_q_info)
       endif

       if (istatlambda2) then
        call decomp_2d_register_variable(io_statistics, "lambda2mean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "lambda22mean", 3, 1, 3, mytype, opt_nplanes=1)
       endif

       if (istatflatness) then
        call decomp_2d_register_variable(io_statistics, "uuuu_mean", 3, 0, 0, mytype, opt_decomp=uuuu_info)
       endif
       if (istatquadrant) then
        call decomp_2d_register_variable(io_statistics, "uv_quadrant_mean", 3, 0, 0, mytype,opt_decomp=uv_quadrant_info)
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

    use param, only : zero, iscalar, istatbudget, istatpstrain
    use param, only : istatquadrant, istatlambda2, nquads, istatspectra
    use param, only : istatflatness, istatautocorr
    use variables, only : nx, ny, nz
    use decomp_2d, only : zsize, decomp_info_init
    use MPI

    implicit none
    integer :: code

    allocate(uvwp_mean(zsize(1),zsize(2),4))
    allocate(uu_mean(zsize(1),zsize(2),6))

    uvwp_mean = zero
    uu_mean = zero

    if (istatbudget) then
      allocate(pu_mean(zsize(1),zsize(2),3))
      allocate(uuu_mean(zsize(1),zsize(2),10))
      allocate(pdudx_mean(zsize(1),zsize(2),3))
      allocate(dudx_mean(zsize(1),zsize(2),9))
      allocate(dudxdudx_mean(zsize(1),zsize(2),12))

      pu_mean = zero
      uuu_mean = zero
      pdudx_mean = zero
      dudx_mean = zero
      dudxdudx_mean = zero
      
    endif

    if (istatflatness) then
      allocate(uuuu_mean(zsize(1),zsize(2),3))
      uuuu_mean = zero
    endif

    if (istatpstrain) then
      if (.not.istatbudget) then
        write(*,*) "If p-strain conditional averaging used, istatbudget must be T"
        call MPI_Abort(MPI_COMM_WORLD,1,code)
      endif

      allocate(pdvdy_q_mean(zsize(1),zsize(2),4))

      pdvdy_q_mean = zero
    endif

    if (istatquadrant) then
      allocate(uv_quadrant_mean(zsize(1),zsize(2),4*nquads))

      uv_quadrant_mean = zero
    endif
    if (istatlambda2) then

      allocate(lambda2mean(zsize(1),zsize(2),1))
      allocate(lambda22mean(zsize(1),zsize(2),1))
      allocate(lambda2(zsize(1),zsize(2),zsize(3)))

      lambda2mean = zero
      lambda22mean = zero
    endif

    if (istatspectra) then
#ifndef HAVE_FFTW
      write(*,*) "Cannot use spectra unless built with FFTw"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
#endif      
      call init_spectra
    endif
    
    if (istatautocorr) then
      call init_autocorrelation
    endif

    call decomp_info_init(nx, ny, 4, uvwp_info)
    call decomp_info_init(nx, ny, 6, uu_info)
    call decomp_info_init(nx, ny, 3, pu_info)
    call decomp_info_init(nx, ny, 10, uuu_info)
    call decomp_info_init(nx, ny, 3, pdudx_info)
    call decomp_info_init(nx, ny, 9, dudx_info)
    call decomp_info_init(nx, ny, 12, dudxdudx_info)

    call decomp_info_init(nx, ny, 4*nquads, uv_quadrant_info)

    call decomp_info_init(nx, ny, 4, pdvdy_q_info)
    call decomp_info_init(nx, ny, 3, uuuu_info)
    call decomp_info_init(nx, ny, 1, dstat_plane)
    
    call init_statistic_adios2

    call grad_init

  end subroutine init_statistic

  subroutine init_spectra
    use param
    use decomp_2d
    use dbg_schemes, only : abs_prec
    use MPI
    use variables, only: nx, nz, ny, yp
    real(mytype), dimension(:,:,:), allocatable :: spectra_z_in
    complex(mytype), dimension(:,:,:), allocatable :: spectra_x_out, spectra_z_out, spectra_x_in
    real(mytype) :: x
    integer :: code, i, j, check, fl

    integer :: color, key
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 
    integer, allocatable, dimension(:) :: code_list
    character(len=80) :: xfmt

#ifdef HAVE_FFTW
    integer :: plan_type = FFTW_MEASURE

    zdft_size = nz/2+1
    xdft_size = nx/2+1

    call decomp_info_init(nx,ny,zdft_size, dft_info)

    if (spectra_level == 1) then
      nspectra = 4
    else if (spectra_level == 2) then
      nspectra = 8
    else if (spectra_level == 3) then
      nspectra = 8 + spectra_corr_nlocs*3
      allocate(spectra_corr_yinds(spectra_corr_nlocs))
      allocate(spectra_corr_ranks(spectra_corr_nlocs))
      
      ! create communicators


      if (nclx) then
        call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, code)
        key = coords(2)
        color = coords(1)

        call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key,color, spectra_corr_comm,code)
        allocate(code_list(dims(1)))

      else
        call MPI_CART_GET(DECOMP_2D_COMM_CART_Z, 2, dims, periods, coords, code)
        key = coords(1)
        color = coords(2)

        call MPI_Comm_split(DECOMP_2D_COMM_CART_Z, key,color, spectra_corr_comm,code)
        allocate(code_list(dims(2)))

      endif

      do i = 1, spectra_corr_nlocs
        if (spectra_corr_ylocs(i)>yp(ny).or.spectra_corr_ylocs(i)<yp(1)) then
          if (nrank==0) write(*,*) "Invalid y location for spectra"
          call MPI_Abort(MPI_COMM_WORLD,1,code)
        endif
        do j = 1, ny
          if (spectra_corr_ylocs(i)<yp(j)) then
            if (abs_prec(spectra_corr_ylocs(i)-yp(j))>abs_prec(spectra_corr_ylocs(i)-yp(j-1))) then
              spectra_corr_yinds(i) = j-1
            else
              spectra_corr_yinds(i) = j
            endif
            exit
          endif
        enddo
        ! send rank for data later

        if (nclx) then
          if (spectra_corr_yinds(i)>=xstart(2).and.spectra_corr_yinds(i)<=xend(2)) then
            check = 1
          else
            check = 0
          endif
        else
          if (spectra_corr_yinds(i)>=zstart(2).and.spectra_corr_yinds(i)<=zend(2)) then
            check = 1
          else
            check = 0
          endif
        endif

        call MPI_Allgather(check,1,MPI_INTEGER,code_list,1,MPI_INTEGER,spectra_corr_comm,code)

        do j =0, size(code_list) -1
          if (code_list(j+1)==1) spectra_corr_ranks(i) = j
        enddo
      enddo

      spectra_corr_thisrank = color
      deallocate(code_list)

      if (nrank ==0) then
        open(newunit=fl,file='parameters.json',status='old',action='readwrite',position='append')
        backspace(fl)
        write(fl,"(A ,' : {')") '  "spectra_corr"'
        if (spectra_corr_nlocs>1) then
          write(xfmt,'(A,I0,A)') "( A, ': [',g0,",spectra_corr_nlocs-1,"(',',g0),']')"
        else
          write(xfmt,'(A,I0,A)') "( A, ': [',g0,']')"
        endif
        write(fl,xfmt) '    "y_locs"',spectra_corr_ylocs
        if (istatautocorr) write(fl,'(A)') "  },"
        if (.not.istatautocorr) write(fl,'(A)') "  }"
          write(fl,'(A)') "}"
        close(fl)
      endif
    else
      write(*,*) "Invalid spectra level"
    call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif

    if (nclx) then
      allocate(spectra_2d_mean(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3),nspectra))
      spectra_2d_mean = zero
    else
      allocate(spectra_z_z_mean(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3),nspectra))
      spectra_z_z_mean = zero
    endif

    allocate(spectra_x_out(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)))
    allocate(spectra_z_out(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)))

    allocate(spectra_x_in(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)))
    allocate(spectra_z_in(zsize(1),zsize(2),zsize(3)))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_r2c(plan_z, 1, zsize(3), &
         zsize(1)*zsize(2), spectra_z_in, zsize(3), &
         zsize(1)*zsize(2), 1, spectra_z_out, zdft_size, &
         zsize(1)*zsize(2), 1, plan_type)
#else
    call sfftw_plan_many_dft_r2c(plan_z, 1, zsize(3), &
         zsize(1)*zsize(2), spectra_z_in, zsize(3), &
         zsize(1)*zsize(2), 1, spectra_z_out, zdft_size, &
         zsize(1)*zsize(2), 1, plan_type)
#endif
   
#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(plan_x, 1, xsize(1), &
        dft_info%xsz(2)*dft_info%xsz(3), spectra_x_in, xsize(1), 1, &
        xsize(1), spectra_x_out, dft_info%xsz(1), 1, dft_info%xsz(1), &
        FFTW_FORWARD, plan_type)
#else
    call sfftw_plan_many_dft(plan_x, 1, xsize(1), &
        dft_info%xsz(2)*dft_info%xsz(3), spectra_x_in, xsize(1), 1, &
        xsize(1), spectra_x_out, dft_info%xsz(1), 1, dft_info%xsz(1), &
        FFTW_FORWARD, plan_type)
#endif
    deallocate(spectra_x_in,spectra_x_out)
    deallocate(spectra_z_in,spectra_z_out)
#endif    
  end subroutine

  subroutine init_autocorrelation
    use param
    use decomp_2d
    use var, only : ny

    integer :: nlocs, xlocs, x_ref
    integer :: i, j, k, fl
    character(len=80) :: xfmt
    nlocs = (xlx - two*autocorr_max_sep) / autocorr_xlocs +1
    xlocs = int(two*autocorr_max_sep/real(dx)) + 1

    allocate(autocorr_mean(xlocs,xsize(2),nlocs))
    autocorr_mean = zero
    allocate(autocorr_x_indices(nlocs))
    do i = 1, nlocs
      autocorr_x_indices(i) = int(autocorr_max_sep/dx) + int((i-1)*autocorr_xlocs/dx,kind=mytype)
    enddo

    if (nrank ==0) then
      open(newunit=fl,file='parameters.json',status='old',action='readwrite',position='append')
      backspace(fl)
      write(fl,"(A ,' : {')") '  "autocorrelation"'
      if (nlocs>1) then
        write(xfmt,'(A,I0,A)') "( A, ': [',g0,",nlocs-1,"(',',g0),']')"
      else
        write(xfmt,'(A,I0,A)') "( A, ': [',g0,']')"
      endif
      write(fl,"(A ,': [',I0,',',I0,',',I0,'],')") '    "shape"',xlocs, ny, nlocs
      write(fl,xfmt) '    "x_locs"',(autocorr_x_indices(:)-1)*dx
      write(fl,'(A)') "  }"
      write(fl,'(A)') "}"
      close(fl)
    endif
  end subroutine
  !
  ! Read all statistics if possible
  !
  subroutine restart_statistic

    use param, only : initstat, irestart, ifirst, itempaccel
    use param, only: zero,itime, istatspectra,initstat2, istatautocorr
    use variables, only : nstat
    use var, only : tmean, nclx

    implicit none
    logical :: read_spectra, read_autocorrelation
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
    read_spectra = istatspectra .and. (itime>=initstat2.or.itempaccel==1)
    read_autocorrelation = istatautocorr .and. (itime>=initstat2.or.itempaccel==1)

    call put_write_read_start(.true.,.true.,read_spectra,read_autocorrelation)
    ! Read all statistics
    call read_or_write_all_stats(.true.)

    if (read_spectra) then
      if(nclx) then
        call read_write_spectra(.true.,spectra_2d_mean,'spectra_2d',1)
      else
        call read_write_spectra(.true.,spectra_z_z_mean,'spectra_z_z',3)
      endif
    endif

    if (read_autocorrelation) then
      call read_autocorr('statistics/'//gen_statname("autocorr_mean"))
    endif

    call put_write_read_end(.true.)

  end subroutine restart_statistic

  function gen_statname(stat) result(newname)

    implicit none
    
    character(len=*), intent(in) :: stat
    character(len=50) :: newname
    
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

    use param, only : iscalar, itime, istatbudget, istatpstrain, nclx, itempaccel
    use param, only : initstat2, istatlambda2, istatquadrant, istatflatness, istatspectra
    use param, only : zpfive, zptwofive, two
    use variables, only : numscalar
    use decomp_2d, only : nrank, zsize, zstart,zend
    use decomp_2d_io, only : decomp_2d_write_mode, decomp_2d_read_mode, &
         decomp_2d_open_io, decomp_2d_close_io, decomp_2d_start_io, decomp_2d_end_io
    use var ,only: nz, nx
    implicit none

    ! Argument
    logical, intent(in) :: flag_read

    ! Local variables
    integer :: is, it
    character(len=30) :: filename
    integer :: io_mode
    integer :: i, j
    real(mytype) :: factor



    if (flag_read) then
       io_mode = decomp_2d_read_mode
    else
       io_mode = decomp_2d_write_mode
    endif
#ifdef ADIOS2
    call decomp_2d_open_io(io_statistics, stat_dir, io_mode)
    call decomp_2d_start_io(io_statistics, stat_dir)
#endif
    
    call read_or_write_one_stat(flag_read, gen_statname("uvwp_mean"), uvwp_mean,uvwp_info)
    call read_or_write_one_stat(flag_read, gen_statname("uu_mean"), uu_mean,uu_info)
   
    if (istatbudget) then
      call read_or_write_one_stat(flag_read, gen_statname("uuu_mean"), uuu_mean, uuu_info)
      call read_or_write_one_stat(flag_read, gen_statname("pdudx_mean"), pdudx_mean, pdudx_info)
      call read_or_write_one_stat(flag_read, gen_statname("dudx_mean"), dudx_mean,dudx_info)
      call read_or_write_one_stat(flag_read, gen_statname("pu_mean"), pu_mean,pu_info)
      call read_or_write_one_stat(flag_read, gen_statname("dudxdudx_mean"), dudxdudx_mean,dudxdudx_info)

    endif

    if (istatflatness) then
      call read_or_write_one_stat(flag_read, gen_statname("uuuu_mean"), uuuu_mean, uuuu_info)
    endif

    if (istatpstrain .and. (itime>=initstat2.or.itempaccel==1)) then
      call read_or_write_one_stat(flag_read, gen_statname("pdvdy_q_mean"), pdvdy_q_mean, pdvdy_q_info)
    endif

    if (istatquadrant .and. (itime>=initstat2.or.itempaccel==1)) then
      call read_or_write_one_stat(flag_read, gen_statname("uv_quadrant_mean"), uv_quadrant_mean, uv_quadrant_info)
    endif

    if (istatlambda2.and. (itime>=initstat2.or.itempaccel==1)) then
      call read_or_write_one_stat(flag_read, gen_statname("lambda2mean"), lambda2mean,dstat_plane)
      call read_or_write_one_stat(flag_read, gen_statname("lambda22mean"), lambda22mean,dstat_plane)
    endif

#ifdef ADIOS2
    call decomp_2d_end_io(io_statistics, stat_dir)
    call decomp_2d_close_io(io_statistics, stat_dir)
#endif
    
  if (flag_read) then

    if (check_stat_correct(gen_statname("pdudx_mean"))) then
      if (nrank ==0 ) write(*,*) "Correcting gradient means on read!!!!"
      factor = zpfive*real(nz,kind=mytype)/(real(nz)-two)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          pdudx_mean(i,j,3) = pdudx_mean(i,j,3)*factor

          dudx_mean(i,j,3) = dudx_mean(i,j,3)*factor
          dudx_mean(i,j,6) = dudx_mean(i,j,6)*factor
          dudx_mean(i,j,9) = dudx_mean(i,j,9)*factor

          if (zstart(1)==1.and.i==1) cycle
          if (zend(1)==nx.and.i==zsize(1)) cycle

          pdudx_mean(i,j,1) = pdudx_mean(i,j,1)*zpfive

          dudx_mean(i,j,1) = dudx_mean(i,j,1)*zpfive
          dudx_mean(i,j,4) = dudx_mean(i,j,4)*zpfive
          dudx_mean(i,j,7) = dudx_mean(i,j,7)*zpfive

          dudxdudx_mean(i,j,1) = dudxdudx_mean(i,j,1)*zptwofive
          dudxdudx_mean(i,j,2) = dudxdudx_mean(i,j,2)*zptwofive
          dudxdudx_mean(i,j,3) = dudxdudx_mean(i,j,3)*zptwofive
          dudxdudx_mean(i,j,4) = dudxdudx_mean(i,j,4)*zptwofive
          dudxdudx_mean(i,j,5) = dudxdudx_mean(i,j,5)*zptwofive
          dudxdudx_mean(i,j,6) = dudxdudx_mean(i,j,6)*zptwofive
                    
        enddo
      enddo
    endif
  endif

  end subroutine read_or_write_all_stats

  function check_stat_correct(fname) result(update)
    use decomp_2d
    use MPI

    character(len=*), intent(in) :: fname
    integer :: values(13), check
    integer :: date(9), ref_date(6)
    integer :: code, i
    logical :: update

    call stat('statistics/'//trim(fname),values,check)
    if (check /=0 .and. nrank==0) then
      write(*,*) "There has been a problem"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif

    ref_date = [0,0,15,13,0,123]
    call ltime(values(10),date)

    do i = 6, 1,-1
      if (ref_date(i)>date(i)) then
        update = .true.
        return
      endif
      if (ref_date(i)<date(i)) then
        update = .false.
        return
      endif
    enddo

  end function
  !
  ! Statistics: perform one IO
  !
  subroutine read_or_write_one_stat(flag_read, filename, array,opt_decomp)

    use decomp_2d, only : mytype, xstS, xenS, zsize
    use decomp_2d_io, only : decomp_2d_read_one, decomp_2d_write_one

    implicit none

    ! Arguments
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename
    real(mytype), dimension(:,:,:), intent(inout) :: array
    type(DECOMP_INFO), intent(in) :: opt_decomp
    integer :: ierror

    if (flag_read) then
       ! There was a check for nvisu = 1 before
       call decomp_2d_read_one(3, array, stat_dir, filename, io_statistics, reduce_prec=.false.,opt_decomp=opt_decomp)
    else
       call decomp_2d_write_one(3, array, stat_dir, filename, 0, io_statistics, reduce_prec=.false.,opt_decomp=opt_decomp)
    endif

  end subroutine read_or_write_one_stat

  subroutine read_write_spectra(flag_read, array, file_base,ipencil)
    use decomp_2d
    use MPI
    use variables, only: nx, ny
    use decomp_2d_io

    logical, intent(in) :: flag_read
    complex(mytype), dimension(:,:,:,:) :: array
    character(len=*) :: file_base
    integer,intent(in) :: ipencil
    character(len=80) :: fn, suffix
    integer :: i

    if (spectra_level.ge.1) then
      if (flag_read) then

        fn = gen_statname('statistics/'//trim(file_base)//'_uu')
        call decomp_2d_read_one(ipencil,array(:,:,:,1),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_vv')
        call decomp_2d_read_one(ipencil,array(:,:,:,2),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_ww')
        call decomp_2d_read_one(ipencil,array(:,:,:,3),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_uv')
        call decomp_2d_read_one(ipencil,array(:,:,:,4),fn,opt_decomp=dft_info)
      else
        fn = gen_statname('statistics/'//trim(file_base)//'_uu')
        call decomp_2d_write_one(ipencil,array(:,:,:,1),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_vv')
        call decomp_2d_write_one(ipencil,array(:,:,:,2),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_ww')
        call decomp_2d_write_one(ipencil,array(:,:,:,3),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_uv')
        call decomp_2d_write_one(ipencil,array(:,:,:,4),fn,opt_decomp=dft_info)

      endif 
    endif

    if (spectra_level.ge.2) then
      if (flag_read) then

        fn = gen_statname('statistics/'//trim(file_base)//'_omega_y')
        call decomp_2d_read_one(ipencil,array(:,:,:,5),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_pdudx')
        call decomp_2d_read_one(ipencil,array(:,:,:,6),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_pdvdy')
        call decomp_2d_read_one(ipencil,array(:,:,:,7),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_pdwdz')
        call decomp_2d_read_one(ipencil,array(:,:,:,8),fn,opt_decomp=dft_info)
      else
        fn = gen_statname('statistics/'//trim(file_base)//'_omega_y')
        call decomp_2d_write_one(ipencil,array(:,:,:,5),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_pdudx')
        call decomp_2d_write_one(ipencil,array(:,:,:,6),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_pdvdy')
        call decomp_2d_write_one(ipencil,array(:,:,:,7),fn,opt_decomp=dft_info)
        fn = gen_statname('statistics/'//trim(file_base)//'_pdwdz')
        call decomp_2d_write_one(ipencil,array(:,:,:,8),fn,opt_decomp=dft_info)

      endif 
    endif
    if (spectra_level.ge.3) then
      if (flag_read) then
        do i = 1, spectra_corr_nlocs
          write(suffix,'("_corr",I0)') i
          fn = gen_statname('statistics/'//trim(file_base)//'_'//'uu'//trim(suffix))
          call decomp_2d_read_one(ipencil,array(:,:,:,8+3*(i-1)+1),fn,opt_decomp=dft_info)

          fn = gen_statname('statistics/'//trim(file_base)//'_'//'vv'//trim(suffix))
          call decomp_2d_read_one(ipencil,array(:,:,:,8+3*(i-1)+2),fn,opt_decomp=dft_info)          
          
          fn = gen_statname('statistics/'//trim(file_base)//'_'//'ww'//trim(suffix))
          call decomp_2d_read_one(ipencil,array(:,:,:,8+3*(i-1)+3),fn,opt_decomp=dft_info)
        enddo
      else
        do i = 1, spectra_corr_nlocs
          write(suffix,'("_corr",I0)') i
          fn = gen_statname('statistics/'//trim(file_base)//'_'//'uu'//trim(suffix))
          call decomp_2d_write_one(ipencil,array(:,:,:,8+3*(i-1)+1),fn,opt_decomp=dft_info)

          fn = gen_statname('statistics/'//trim(file_base)//'_'//'vv'//trim(suffix))
          call decomp_2d_write_one(ipencil,array(:,:,:,8+3*(i-1)+2),fn,opt_decomp=dft_info)          
          
          fn = gen_statname('statistics/'//trim(file_base)//'_'//'ww'//trim(suffix))
          call decomp_2d_write_one(ipencil,array(:,:,:,8+3*(i-1)+3),fn,opt_decomp=dft_info)
        enddo
      endif 
    endif
  end subroutine

  subroutine read_autocorr(filename)
    use decomp_2d
    use MPI
    use var
    use param

    character(len=*) :: filename
    integer :: sizes(3), subsizes(3), starts(3)
    integer :: fh, split_comm_z, key, color, newtype, code
    integer :: nlocs, xlocs
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 

    nlocs = (xlx - two*autocorr_max_sep) / autocorr_xlocs + 1
    xlocs = int(two*autocorr_max_sep/real(dx)) + 1

    sizes(:) = [xlocs, ny, nlocs]
    subsizes(:) = [xlocs, xsize(2), nlocs]
    starts(:) = [0, xstart(2)-1, 0]

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, code)
    key = coords(2)
    color = coords(1)

    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key,color, split_comm_z,code)

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
           MPI_ORDER_FORTRAN, real_type, newtype, code)
    call MPI_TYPE_COMMIT(newtype,code)

    call MPI_File_open(split_comm_z,filename,MPI_MODE_RDONLY, MPI_INFO_NULL, &
                      fh, code)
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,real_type, &
                      newtype,'native',MPI_INFO_NULL,code)
    call MPI_FILE_READ_ALL(fh, autocorr_mean, &
                      subsizes(1)*subsizes(2)*subsizes(3), &
                      real_type, MPI_STATUS_IGNORE, code)
    call MPI_FILE_CLOSE(fh,code)
    call MPI_TYPE_FREE(newtype,code)

    call MPI_Comm_free(split_comm_z,code)

  end subroutine read_autocorr

  subroutine write_autocorr(filename)
    use decomp_2d
    use MPI
    use var
    use param

    character(len=*) :: filename
    integer :: sizes(3), subsizes(3), starts(3)
    integer :: fh, split_comm_z, key, color, newtype, code
    integer :: nlocs, xlocs
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 
    real(mytype), dimension(:), allocatable :: autocorr_g
    integer, dimension(:), allocatable :: recvcounts, displs

    nlocs = (xlx - two*autocorr_max_sep) / autocorr_xlocs + 1
    xlocs = int(two*autocorr_max_sep/real(dx)) + 1

    sizes(:) = [xlocs, ny, nlocs]
    subsizes(:) = [xlocs, xsize(2), nlocs]
    starts(:) = [0, xstart(2)-1, 0]

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, code)
    key = coords(2)
    color = coords(1)

    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key,color, split_comm_z,code)

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
           MPI_ORDER_FORTRAN, real_type, newtype, code)
    call MPI_TYPE_COMMIT(newtype,code)

    if (key == 0) then
      call MPI_File_open(split_comm_z,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,&
                        MPI_INFO_NULL, fh, code)
      call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,code)                        
      call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,real_type, &
                        newtype,'native',MPI_INFO_NULL,code)
      call MPI_FILE_WRITE_ALL(fh, autocorr_mean, &
                        subsizes(1)*subsizes(2)*subsizes(3), &
                        real_type, MPI_STATUS_IGNORE, code)
      call MPI_FILE_CLOSE(fh,code)
    endif
    call MPI_TYPE_FREE(newtype,code)
    call MPI_Comm_free(split_comm_z,code)
  
  end subroutine write_autocorr
  !
  ! Statistics : Intialize, update and perform IO
  !
  subroutine put_write_read_start(flag_read,write,write_spectra,write_autocorrelation)
    use decomp_2d
    use var, only: itime
    logical, intent(in) :: flag_read,write,write_spectra,write_autocorrelation

    integer :: it
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
        print *,'Reading stat files', stats_time
      else
        print *,'Writing stat files', stats_time
      endif
      if (write) print *, '    basic statistics'
      if (write_autocorrelation) print *, '    autocorrelation'
      if (write_spectra) print *, '    spectra'
  
    endif
    
  end subroutine

  subroutine put_write_read_end(flag_read)
    use decomp_2d
    logical, intent(in) :: flag_read

    if (nrank==0) then
      if (flag_read) then
        print *,'Read stat done!'
      else
        print *,'Write stat done!'
      endif
      print *,'==========================================================='
    endif

  end subroutine
  subroutine overall_statistic(ux1,uy1,uz1,phi1,pp3,ep1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io
    use tools, only : rescale_pressure
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,di3
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
    logical :: any, depend,spectra, write, write_spectra

    if (itime.lt.initstat) then
       return
    elseif (itime.eq.initstat) then
       call init_statistic()
       
    elseif (itime.eq.ifirst) then
       if (itempaccel.eq.1) then
          call init_statistic()
       else
          call restart_statistic()
       endif
    endif
    
    call process_statistics(any, depend,spectra, write, write_spectra)
    if (any) then

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

      call update_average_scalar(uvwp_mean(:,:,4), td3, ep1)

      !! Mean velocity
      call update_average_vector(uvwp_mean(:,:,1),uvwp_mean(:,:,2), uvwp_mean(:,:,3), &
                                ux3, uy3, uz3, td3)

      !! Second-order velocity moments
      call update_variance_vector(uu_mean(:,:,1),uu_mean(:,:,2), uu_mean(:,:,3),&
                                  uu_mean(:,:,4), uu_mean(:,:,5), uu_mean(:,:,6), &
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

      call update_skewness_tensor(uuu_mean(:,:,1), uuu_mean(:,:,2), uuu_mean(:,:,3),&
                                uuu_mean(:,:,4),uuu_mean(:,:,5),uuu_mean(:,:,6),&
                                uuu_mean(:,:,7),uuu_mean(:,:,8),uuu_mean(:,:,9),&
                                uuu_mean(:,:,10), ux3, uy3, uz3, td3)  

      call update_pvelograd_vector(pdudx_mean(:,:,1), pdudx_mean(:,:,2), pdudx_mean(:,:,3),&
                                  dudx, dvdy, tc3, td3)

      call update_pvelo_vector(pu_mean(:,:,1), pu_mean(:,:,2), pu_mean(:,:,3), ux3, uy3, uz3,td3)

      call update_velograd_tensor(dudx_mean, dudx,dudy,ta3,dvdx,dvdy,&
                                  tb3,dwdx,dwdy,tc3,td3)

      call update_velograd2_tensor(dudxdudx_mean,&
                                  dudx,dvdx,dwdx,dudy,dvdy,dwdy,td3)
                        
      endif                             

      if (istatflatness) then
        call update_flatness(uuuu_mean,ux3, uy3, uz3,td3)
      endif
      if (istatpstrain .and. depend) then
        call update_pstrain_cond_avg(pdvdy_q_mean(:,:,1),pdvdy_q_mean(:,:,2),&
                                    pdvdy_q_mean(:,:,3),pdvdy_q_mean(:,:,4),&
                                    td3,dvdy,uvwp_mean(:,:,4),dudx_mean(:,:,5),td3)
      endif

      if (istatquadrant .and. depend) then
        call update_uv_quad_avg(uv_quadrant_mean, uu_mean(:,:,1), uu_mean(:,:,2),&
                                uvwp_mean(:,:,1), uvwp_mean(:,:,2),ux3,uy3, td3)
      endif

      if (istatspectra.and. spectra) then
        if (nclx) then
          call update_spectra_avg_xz(spectra_2d_mean,ux3,uy3,uz3,td3,dvdy, ta3, dwdx)
        else
          call update_spectra_avg_z(spectra_z_z_mean,ux3,uy3,uz3,td3,dudx, dvdy, ta3,dwdx)
        endif
      endif

      if (istatautocorr.and. spectra) then
        call update_autocorr_x(autocorr_mean,ux3)
      endif

      if (istatlambda2.and. depend) then
        call update_lambda2_avg(lambda2mean,lambda22mean,&
                                  dudx,dudy,ta3,dvdx,dvdy,&
                                  tb3,dwdx,dwdy,tc3,td3)
      endif
    endif

    ! Write all statistics
    if (write .or. write_spectra) call put_write_read_start(.false., write,&
                                                           write_spectra.and.istatspectra,&
                                                           write_spectra.and.istatautocorr)
    if (write) then
       call read_or_write_all_stats(.false.)
    endif

    if (istatautocorr .and. write_spectra) then
      call write_autocorr('statistics/'//gen_statname("autocorr_mean"))
   endif

    if (istatspectra .and. write_spectra) then
      if (mod(itime,ispectout)==0) then

        if(nclx) then
          call read_write_spectra(.false.,spectra_2d_mean,'spectra_2d',1)
        else
          call read_write_spectra(.false.,spectra_z_z_mean,'spectra_z_z',3)
        endif
      endif
    endif
    if (write .or. write_spectra) call put_write_read_end(.false.)

  end subroutine overall_statistic

  subroutine process_statistics(any,depend,spectra,write,write_spectra)
    use param, only: ispectstart, itempaccel, istatcalc, initstat2, ispectout, istatout
    use var, only: itime
    logical, intent(out) :: any, depend, spectra, write, write_spectra

    
    if (itempaccel == 1) then
      any = mod(itime,istatout)==0
      depend = any
      spectra = mod(itime,ispectout)==0 .and. itime >=ispectstart
      write = any 
      write_spectra =  spectra
    else
      any = mod(itime,istatcalc) == 0
      depend = itime >= initstat2
      spectra = depend
      write = mod(itime,istatout)==0
      write_spectra = spectra .and. mod(itime,ispectout)==0
    endif

  end subroutine process_statistics
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
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: um
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
      um(:,:) = stat_z(:,:)
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
        um(i,j) = um(i,j) + (stat_z(i,j) - um(i,j))/ stat_inc
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
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: um, vm, wm
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
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: uum, vvm, wwm, uvm, uwm, vwm
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
  real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: uuum,uuvm,uuwm,uvvm,uvwm,&
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

    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pdudxm,pdvdym,pdwdzm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: dudx, dvdy, dwdz, ep

    call update_average_scalar(pdudxm, ep*dudx, ep)
    call update_average_scalar(pdvdym, ep*dvdy, ep)
    call update_average_scalar(pdwdzm, ep*dwdz, ep)
  end subroutine

  subroutine update_pvelo_vector(pum, pvm, pwm, ux, uy, uz, ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pum, pvm, pwm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(pum, ep*ux, ep)
    call update_average_scalar(pvm, ep*uy, ep)
    call update_average_scalar(pwm, ep*uz, ep)

  end subroutine update_pvelo_vector

  subroutine update_velograd_tensor(dudxm, dudx, dudy, dudz, dvdx,&
                                    dvdy, dvdz, dwdx, dwdy, dwdz, ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),9), intent(inout) :: dudxm


    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dudy, dudz, ep
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dvdx, dvdy, dvdz   
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dwdx, dwdy, dwdz

    call update_average_scalar(dudxm(:,:,1), dudx, ep)
    call update_average_scalar(dudxm(:,:,2), dudy, ep)
    call update_average_scalar(dudxm(:,:,3), dudz, ep)
    call update_average_scalar(dudxm(:,:,4), dvdx, ep)
    call update_average_scalar(dudxm(:,:,5), dvdy, ep)
    call update_average_scalar(dudxm(:,:,6), dvdz, ep)
    call update_average_scalar(dudxm(:,:,7), dudx, ep)
    call update_average_scalar(dudxm(:,:,8), dwdy, ep)
    call update_average_scalar(dudxm(:,:,9), dwdz, ep)


  end subroutine update_velograd_tensor

  subroutine update_velograd2_tensor(dudxdudxm,&
                                     dudx,dvdx,dwdx,dudy,dvdy,dwdy,ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),12), intent(inout) :: dudxdudxm

    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dvdx, dwdx, ep
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudy, dvdy, dwdy

    call update_average_scalar(dudxdudxm(:,:,1), dudx*dudx, ep)
    call update_average_scalar(dudxdudxm(:,:,2), dudx*dvdx, ep)
    call update_average_scalar(dudxdudxm(:,:,3), dudx*dwdx, ep)
    call update_average_scalar(dudxdudxm(:,:,4), dvdx*dvdx, ep)
    call update_average_scalar(dudxdudxm(:,:,5), dvdx*dwdx, ep)
    call update_average_scalar(dudxdudxm(:,:,6), dwdx*dwdx, ep)

    call update_average_scalar(dudxdudxm(:,:,7), dudy*dudy, ep)
    call update_average_scalar(dudxdudxm(:,:,8), dudy*dvdy, ep)
    call update_average_scalar(dudxdudxm(:,:,9), dudy*dwdy, ep)
    call update_average_scalar(dudxdudxm(:,:,10), dvdy*dvdy, ep)
    call update_average_scalar(dudxdudxm(:,:,11), dvdy*dwdy, ep)
    call update_average_scalar(dudxdudxm(:,:,12), dwdy*dwdy, ep)

  end subroutine update_velograd2_tensor
  
  subroutine update_flatness(uuuum, ux, uy, uz, ep)
    use decomp_2d, only : mytype, xsize, zsize
    real(mytype), dimension(zsize(1),zsize(2),3), intent(inout) :: uuuum
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep
  
    call update_average_scalar(uuuum(:,:,1), ux*ux*ux*ux, ep)
    call update_average_scalar(uuuum(:,:,2), uy*uy*uy*uy, ep)
    call update_average_scalar(uuuum(:,:,3), uz*uz*uz*uz, ep)

  end subroutine update_flatness

  subroutine update_uv_quad_avg(uvqm, uum, vvm, um, vm, ux, uy, ep)
    use dbg_schemes, only : sqrt_prec, abs_prec
    use decomp_2d, only : mytype, xsize, zsize
    use param, only: nquads, zero
    use MPI
    real(mytype), dimension(zsize(1),zsize(2),4*nquads), intent(inout) :: uvqm
    real(mytype), dimension(zsize(1),zsize(2)), intent(in) :: um, vm
    real(mytype), dimension(zsize(1),zsize(2)), intent(in) :: uum, vvm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, ep

    logical, dimension(:,:,:,:,:) , allocatable :: mask
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: uv_fluct

    logical :: quad(4), thresh
    integer :: i, j, k, l, ierr
    real(mytype) :: u_rms, v_rms, u_fluct, v_fluct

    allocate(mask(zsize(1),zsize(2),zsize(3),nquads,4))
    uv_fluct = zero
      do k = 1, zsize(3)
        do j = 1, zsize(2)
          do i = 1, zsize(1)
            u_rms = sqrt_prec(uum(i,j) - um(i,j)*um(i,j))
            v_rms = sqrt_prec(vvm(i,j) - vm(i,j)*vm(i,j))

            u_fluct = ux(i,j,k) - um(i,j)
            v_fluct = uy(i,j,k) - vm(i,j)

            uv_fluct(i,j,k) = u_fluct*v_fluct

            quad(1) =  u_fluct .gt. zero .and. v_fluct .gt. zero
            quad(2) =  u_fluct .lt. zero .and. v_fluct .gt. zero
            quad(3) =  u_fluct .lt. zero .and. v_fluct .lt. zero
            quad(4) =  u_fluct .gt. zero .and. v_fluct .lt. zero

            do l = 1, nquads

              thresh = abs_prec(uv_fluct(i,j,k)) > h_quads(l)*u_rms*v_rms
              

              mask(i,j,k,l,1) =  quad(1) .and. thresh
              mask(i,j,k,l,2) =  quad(2) .and. thresh
              mask(i,j,k,l,3) =  quad(3) .and. thresh
              mask(i,j,k,l,4) =  quad(4) .and. thresh

          enddo
        enddo
      enddo
    enddo

    if (any(isnan(uv_fluct))) then
      write(*,*) "NaN detected in quadrant analysis"
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif 
    do j = 1, nquads
      do i = 1, 4
        call update_average_scalar(uvqm(:,:,(j-1)*4+i),uv_fluct,&
                                  ep,mask=mask(:,:,:,j,i),istat2=.true.)
      enddo
    enddo

    deallocate(mask)

  end subroutine update_uv_quad_avg
  subroutine update_pstrain_cond_avg(pdvdy_q1m,pdvdy_q2m,&
                                     pdvdy_q3m,pdvdy_q4m,&
                                     p,dvdy, pm, dvdym,ep)
    use decomp_2d, only : mytype, xsize, zsize
    use param, only : zero, zpfive
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pdvdy_q1m, pdvdy_q2m
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pdvdy_q3m, pdvdy_q4m
    real(mytype), dimension(zsize(1),zsize(2)), intent(in) :: pm, dvdym
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: p,dvdy, ep

    logical, dimension(zsize(1),zsize(2),zsize(3)) :: mask1, mask2, mask3, mask4
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: p_fluct, dvdy_fluct
    integer :: i,j,k

    do k = 1, zsize(3)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          p_fluct(i,j,k) = p(i,j,k) - pm(i,j)
          dvdy_fluct(i,j,k) = dvdy(i,j,k) - dvdym(i,j)
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
    use param
    USE dbg_schemes, only : sqrt_prec, acos_prec, cos_prec

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: lambda2m,lambda22m
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dudy, dudz
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dvdx, dvdy, dvdz
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dwdx, dwdy, dwdz, ep

    real(mytype), dimension(3,3) :: S, Omega, S2O2, U_mat
    real(mytype) :: p1, p2, p, q, det_B, phi
    integer :: i,j,k

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

          p1 = S2O2(1,1)*S2O2(1,1) + S2O2(2,2)*S2O2(2,2) + S2O2(3,3)*S2O2(3,3)
          p2 = S2O2(1,2)*S2O2(1,2) + S2O2(1,3)*S2O2(1,3) + S2O2(2,3)*S2O2(2,3)

          p = sqrt_prec((p1 + two*p2)/six)
          q = (S2O2(1,1) + S2O2(2,2) + S2O2(3,3))/three

          U_mat(1,1) = (S2O2(1,1) - q)/p
          U_mat(2,2) = (S2O2(2,2) - q)/p
          U_mat(3,3) = (S2O2(3,3) - q)/p
          U_mat(1,2) = S2O2(1,2)/p
          U_mat(2,1) = U_mat(1,2)
          U_mat(1,3) = S2O2(1,3)/p
          U_mat(3,1) = U_mat(1,3)
          U_mat(2,3) = S2O2(2,3)/p
          U_mat(3,2) = U_mat(2,3)
          
          det_B = U_mat(1,1)*(U_mat(2,2)*U_mat(3,3) -U_mat(2,3)*U_mat(3,2) ) &
                  - U_mat(1,2)*(U_mat(2,1)*U_mat(3,3) -U_mat(2,3)*U_mat(3,1) ) &
                   + U_mat(1,3)*(U_mat(2,1)*U_mat(3,2) -U_mat(2,2)*U_mat(3,1) )

          phi = acos_prec(zpfive*det_B)/3 + four*pi/three
          lambda2(i,j,k) = q + two*p*cos_prec(phi)
        enddo
      enddo
    enddo

    call update_average_scalar(lambda2m, lambda2, ep,istat2=.true.)
    call update_average_scalar(lambda22m, lambda2*lambda2, ep,istat2=.true.)
  end subroutine

  subroutine fft_2d_calc(fft_2d,val)
    use variables, only: nx,nz
    use param
    use decomp_2d
    use MPI
    use param, only: zero

    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: val
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(out) :: fft_2d

    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)) :: spec_z
    complex(mytype), dimension(dft_info%ysz(1),dft_info%ysz(2),dft_info%ysz(3)) :: spec_y
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)) :: spec_x

#ifdef HAVE_FFTW
    call dfftw_execute_dft_r2c(plan_z,val,spec_z)
    call transpose_z_to_y(spec_z,spec_y,dft_info)
    call transpose_y_to_x(spec_y,spec_x,dft_info)

    call dfftw_execute_dft(plan_x,spec_x,fft_2d)
#endif
  end subroutine
  subroutine spectra_2d_calc(spec_2d,spec1, spec2)
    use variables, only: nx,nz
    use param
    use decomp_2d
    use MPI
    use param, only: zero

    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(in) :: spec1
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(in) :: spec2
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(out) :: spec_2d
    real(mytype) :: ar, ai, br, bi
    integer :: i,j,k
    real(mytype) :: norm

    norm = real(nz,kind=mytype)*real(nz,kind=mytype)*real(nx,kind=mytype)*real(nx,kind=mytype)
    spec_2d(:,:,:) = zero
    do k = 1, dft_info%xsz(3)
      do j = 1, dft_info%xsz(2)
        do i = 1, dft_info%xsz(1)
          
          ar = spec1(i,j,k)%re
          ai = spec1(i,j,k)%im
          br = spec2(i,j,k)%re
          bi = spec2(i,j,k)%im
          spec_2d(i,j,k) = cmplx(ar*br + ai*bi,-ai*br + bi*ar)/norm
          ! spec_2d(i,j,k) = conjg(spec1(i,j,k))*spec2(i,j,k)/norm

        enddo
      enddo
    enddo

  end subroutine
  subroutine spectra_z_calc(spec_zm,spec1,spec2)
    use variables, only: nx,nz
    use param
    use decomp_2d
    use MPI
    use param, only: zero

    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)) :: spec_zm
    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)), intent(in) :: spec1, spec2

    ! local args
    real(mytype) :: norm
    integer :: i, j, k
    real(mytype) :: ar, ai, br, bi

    ! calculate spectra of local rank
    norm = real(nz,kind=mytype)*real(nz,kind=mytype)

    do k = 1, zdft_size
      do j = 1,zsize(2)
        do i = 1, zsize(1)
          ar = spec1(i,j,k)%re
          ai = spec1(i,j,k)%im
          br = spec2(i,j,k)%re
          bi = spec2(i,j,k)%im
          spec_zm(i,j,k) = cmplx(ar*br + ai*bi,-ai*br + bi*ar)/norm
          ! spec_zm(i,j,k) = conjg(spec1(i,j,k))*spec2(i,j,k)/norm
        enddo
      enddo
    enddo
    
  end subroutine

  subroutine spectra_corr_2d_calc(spec_2d,spec1, spec2,index)
    use variables, only: nx,nz, yp
    use param
    use decomp_2d
    use MPI
    use param, only: zero

    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(in) :: spec1
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(in) :: spec2
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)), intent(out) :: spec_2d
    integer, intent(in) :: index
    ! local args
    integer :: i,j,k
    real(mytype) :: norm
    integer :: local_index, spec_size, code
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(3)) :: spec_plane
    real(mytype) :: ar, br, ai, bi
    ! calculate spectra of local rank
    norm = real(nz,kind=mytype)*real(nz,kind=mytype)*real(nx,kind=mytype)*real(nx,kind=mytype)
    spec_size = dft_info%xsz(1)*dft_info%xsz(3)

    ! broadcast data
    if (spectra_corr_ranks(index)==spectra_corr_thisrank) then
      local_index = spectra_corr_yinds(index) - xstart(2) + 1

      do k = 1, dft_info%xsz(3)
        do i = 1, dft_info%xsz(1)
          spec_plane(i,k) = spec1(i,local_index,k)
        enddo
      enddo
    endif

    call MPI_Bcast(spec_plane,spec_size,complex_type,&
                  spectra_corr_ranks(index),spectra_corr_comm,&
                  code)

    spec_2d(:,:,:) = zero
    do k = 1, dft_info%xsz(3)
      do j = 1, dft_info%xsz(2)
        do i = 1, dft_info%xsz(1)
          ar = spec_plane(i,k)%re
          ai = spec_plane(i,k)%im
          br = spec2(i,j,k)%re
          bi = spec2(i,j,k)%im
          spec_2d(i,j,k) = cmplx(ar*br + ai*bi,-ai*br + bi*ar)/norm

          ! spec_2d(i,j,k) = conjg(spec_plane(i,k))*spec2(i,j,k)/norm

        enddo
      enddo
    enddo
  end subroutine

  subroutine spectra_corr_z_calc(spec_zm,spec1,spec2, index)
    use variables, only: nx,nz
    use param
    use decomp_2d
    use MPI
    use param, only: zero

    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)) :: spec_zm
    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)), intent(in) :: spec1, spec2
    integer :: index

    ! local args
    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(3)) :: spec_plane
    real(mytype) :: norm
    integer :: i, j, k
    integer :: local_index, spec_size, code
    real(mytype) :: ar, ai, br, bi


    ! calculate spectra of local rank
    norm = real(nz,kind=mytype)*real(nz,kind=mytype)


    ! calculate spectra of local rank
    spec_size = dft_info%zsz(1)*dft_info%zsz(3)

    ! broadcast data
    if (spectra_corr_ranks(index)==spectra_corr_thisrank) then
      local_index = spectra_corr_yinds(index) - zstart(2) + 1

      do k = 1, dft_info%zsz(3)
        do i = 1, dft_info%zsz(1)
          spec_plane(i,k) = spec1(i,local_index,k)
        enddo
      enddo
    endif

    call MPI_Bcast(spec_plane,spec_size,complex_type,&
                  spectra_corr_ranks(index),spectra_corr_comm,&
                  code)

    do k = 1, zdft_size
      do j = 1,zsize(2)
        do i = 1, zsize(1)
          ar = spec_plane(i,k)%re
          ai = spec_plane(i,k)%im
          br = spec2(i,j,k)%re
          bi = spec2(i,j,k)%im
          spec_zm(i,j,k) = cmplx(ar*br + ai*bi,-ai*br + bi*ar)/norm

          ! spec_zm(i,j,k) = conjg(spec_plane(i,k))*spec2(i,j,k)/norm
        enddo
      enddo
    enddo
  end subroutine
  subroutine update_spectra_avg_xz(spec_2d_m, ux3, uy3, uz3, p3,dvdy3, dudz3, dwdx3)
    use MPI
    use decomp_2d
    use variables, only : nx, nz
    use param, only : itime, initstat,initstat2, istatcalc, itempaccel, initstat2
    use var, only: te3, tf3, tg3

    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3),nspectra) :: spec_2d_m

    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ux3, uy3, uz3, p3
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: dudz3,dvdy3,dwdx3

    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3),nspectra) :: spec_2d_ml
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)) :: u_spec, v_spec, w_spec, p_spec
    complex(mytype), dimension(dft_info%xsz(1),dft_info%xsz(2),dft_info%xsz(3)) :: u_spec1, w_spec1,dvdy_spec, omega_spec

    real(mytype), dimension(zsize(2),6) :: uvw_m
    real(mytype) :: stat_inc
    integer :: i,j,k

    call compute_avg_xz_u(uvw_m)
    if (spectra_level.ge.1) then

      do k =1, zsize(3)
        do j =1, zsize(2)
          do i =1, zsize(1)
            te3(i,j,k) = ux3(i,j,k) - uvw_m(j,1)
            tf3(i,j,k) = uy3(i,j,k) - uvw_m(j,2)
            tg3(i,j,k) = uz3(i,j,k) - uvw_m(j,3)
          enddo
        enddo
      enddo

      call fft_2d_calc(u_spec,te3)
      call fft_2d_calc(v_spec,tf3)
      call fft_2d_calc(w_spec,tg3)

      call spectra_2d_calc(spec_2d_ml(:,:,:,1),u_spec,u_spec)
      call spectra_2d_calc(spec_2d_ml(:,:,:,2),v_spec,v_spec)
      call spectra_2d_calc(spec_2d_ml(:,:,:,3),w_spec,w_spec)
      call spectra_2d_calc(spec_2d_ml(:,:,:,4),u_spec,v_spec)
    
    endif

    if (spectra_level .ge.2) then
      do k =1, zsize(3)
        do j =1, zsize(2)
          do i =1, zsize(1)
            te3(i,j,k) = p3(i,j,k) - uvw_m(j,4)
            tf3(i,j,k) = dvdy3(i,j,k) - uvw_m(j,5)
            tg3(i,j,k) = dudz3(i,j,k) - dwdx3(i,j,k) - uvw_m(j,6)
          enddo
        enddo
      enddo

      call fft_2d_calc(p_spec,te3)
      call fft_2d_calc(dvdy_spec,tf3)
      call fft_2d_calc(omega_spec,tg3)

      do k = 1, dft_info%xsz(3)
        do j = 1, dft_info%xsz(2)
          do i = 1, dft_info%xsz(1)
            u_spec1(i,j,k) = u_spec(i,j,k)*cmplx(0,1,kind=mytype)
            w_spec1(i,j,k) = w_spec(i,j,k)*cmplx(0,1,kind=mytype)
          enddo
        enddo
      enddo

      call spectra_2d_calc(spec_2d_ml(:,:,:,5),omega_spec,v_spec)
      call spectra_2d_calc(spec_2d_ml(:,:,:,6),p_spec,u_spec1)
      call spectra_2d_calc(spec_2d_ml(:,:,:,7),p_spec,dvdy_spec)
      call spectra_2d_calc(spec_2d_ml(:,:,:,8),p_spec,w_spec1)

    endif

    if (spectra_level.ge.3) then
      do i = 1, spectra_corr_nlocs
        call spectra_corr_2d_calc(spec_2d_ml(:,:,:,8+3*(i-1)+1),u_spec,u_spec,i)
        call spectra_corr_2d_calc(spec_2d_ml(:,:,:,8+3*(i-1)+2),v_spec,v_spec,i)
        call spectra_corr_2d_calc(spec_2d_ml(:,:,:,8+3*(i-1)+3),w_spec,w_spec,i)
      enddo
    endif
    if (itempaccel==1) then
      spec_2d_m = spec_2d_ml

    else
      stat_inc = real((itime-initstat2)/istatcalc+1, kind=mytype)

      spec_2d_m = spec_2d_m + (spec_2d_ml - spec_2d_m)/ stat_inc
    endif

  end subroutine

  subroutine compute_avg_xz_u(uvw_m)
    use MPI
    use decomp_2d
    use var, only: nx
    use param, only: zero
    real(mytype), dimension(zsize(2),6), intent(out) :: uvw_m
    real(mytype), dimension(zsize(2),6) :: uvw_l
    
    integer :: color, key
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 
    integer :: split_comm_y
    integer :: i,j,k, code

    uvw_l = zero
    do j =1, zsize(2)
      do i =1,zsize(1)
        uvw_l(j,1) = uvw_l(j,1) + uvwp_mean(i,j,1)
        uvw_l(j,2) = uvw_l(j,2) + uvwp_mean(i,j,2)
        uvw_l(j,3) = uvw_l(j,3) + uvwp_mean(i,j,3)
      enddo
    enddo

    if (spectra_level.ge.2) then
      do j =1, zsize(2)
        do i =1,zsize(1)
          uvw_l(j,4) = uvw_l(j,4) + uvwp_mean(i,j,4)
          uvw_l(j,5) = uvw_l(j,5) + dudx_mean(i,j,5)
          uvw_l(j,6) = uvw_l(j,6) + dudx_mean(i,j,3) - dudx_mean(i,j,7)
        enddo
      enddo
    endif
    uvw_l(:,:) = uvw_l(:,:)/real(nx,kind=mytype)

    call MPI_CART_GET(DECOMP_2D_COMM_CART_Z, 2, dims, periods, coords, code)
    key = coords(2)
    color = coords(1)

    call MPI_Comm_split(DECOMP_2D_COMM_CART_Z, key,color, split_comm_y,code)

    call MPI_Allreduce(uvw_l,uvw_m,size(uvw_m),&
                       real_type,MPI_SUM,split_comm_y,code)
    call MPI_Comm_free(split_comm_y,code)

  end subroutine
  
  subroutine update_spectra_avg_z(spec_zm, ux3, uy3, uz3,p3,dudx3, dvdy3, dudz3,dwdx3)
    use MPI
    use decomp_2d
    use variables, only : nx, nz
    use param, only : itime, initstat,initstat2, istatcalc, itempaccel, initstat2
    use var, only : te3, tf3, tg3, th3
    implicit none

    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3),nspectra) :: spec_zm

    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ux3, uy3, uz3, p3
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: dudx3,dvdy3,dudz3,dwdx3
    
    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3),nspectra) :: spec_z_ml
    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)) :: u_spec, v_spec, w_spec, w_spec1
    complex(mytype), dimension(dft_info%zsz(1),dft_info%zsz(2),dft_info%zsz(3)) :: p_spec, omega_spec,dudx_spec, dvdy_spec

    real(mytype) :: stat_inc
    integer :: i,j,k

    if (spectra_level.ge.1) then

      do k =1, zsize(3)
        do j =1, zsize(2)
          do i =1, zsize(1)
            te3(i,j,k) = ux3(i,j,k) - uvwp_mean(i,j,1)
            tf3(i,j,k) = uy3(i,j,k) - uvwp_mean(i,j,2)
            tg3(i,j,k) = uz3(i,j,k) - uvwp_mean(i,j,3)
          enddo
        enddo
      enddo

#ifdef HAVE_FFTW
      call dfftw_execute_dft_r2c(plan_z,te3,u_spec)
      call dfftw_execute_dft_r2c(plan_z,tf3,v_spec)
      call dfftw_execute_dft_r2c(plan_z,tg3,w_spec)
#endif

      call spectra_z_calc(spec_z_ml(:,:,:,1),u_spec,u_spec)
      call spectra_z_calc(spec_z_ml(:,:,:,2),v_spec,v_spec)
      call spectra_z_calc(spec_z_ml(:,:,:,3),w_spec,w_spec)
      call spectra_z_calc(spec_z_ml(:,:,:,4),u_spec,v_spec)
    
    endif

    if (spectra_level .ge. 2) then

      do k =1, zsize(3)
        do j =1, zsize(2)
          do i =1, zsize(1)
            te3(i,j,k) = p3(i,j,k) - uvwp_mean(i,j,4)
            tf3(i,j,k) = dudx3(i,j,k) - dudx_mean(i,j,1)
            tg3(i,j,k) = dvdy3(i,j,k) - dudx_mean(i,j,5)
            th3(i,j,k) = dudz3(i,j,k) - dwdx3(i,j,k) - dudx_mean(i,j,3) + dudx_mean(i,j,7)
          enddo
        enddo
      enddo
      
      w_spec1 = cmplx(-aimag(w_spec),real(w_spec))
#ifdef HAVE_FFTW
      call dfftw_execute_dft_r2c(plan_z,te3,p_spec)
      call dfftw_execute_dft_r2c(plan_z,tf3,dudx_spec)
      call dfftw_execute_dft_r2c(plan_z,tg3,dvdy_spec)
      call dfftw_execute_dft_r2c(plan_z,th3,omega_spec)
#endif
      call spectra_z_calc(spec_z_ml(:,:,:,5),omega_spec,v_spec)
      call spectra_z_calc(spec_z_ml(:,:,:,6),p_spec,dudx_spec)
      call spectra_z_calc(spec_z_ml(:,:,:,7),p_spec,dvdy_spec)
      call spectra_z_calc(spec_z_ml(:,:,:,8),p_spec,w_spec1)

    endif

    if (spectra_level.ge.3) then
      do i = 1, spectra_corr_nlocs
        call spectra_corr_z_calc(spec_z_ml(:,:,:,8+3*(i-1)+1),u_spec,u_spec,i)
        call spectra_corr_z_calc(spec_z_ml(:,:,:,8+3*(i-1)+2),v_spec,v_spec,i)
        call spectra_corr_z_calc(spec_z_ml(:,:,:,8+3*(i-1)+3),w_spec,w_spec,i)
      enddo
    endif

    if (itempaccel==1) then
      spec_zm = spec_z_ml

    else
      stat_inc = real((itime-initstat2)/istatcalc+1, kind=mytype)

      spec_zm = spec_zm + (spec_z_ml - spec_zm)/ stat_inc
    endif    
  end subroutine

  subroutine update_autocorr_x(autocorr_m,ux3)
    use decomp_2d
    use param
    use var, only : nz
    use MPI
    real(mytype), dimension(:,:,:) :: autocorr_m
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ux3

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: u_fluct1
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: u_fluct2
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: u_fluct3
    real(mytype), dimension(:,:,:), allocatable :: autocorr_ml, autocorr_ml2
    integer :: i, j, k, l, stat_inc
    integer :: nlocs, xlocs,autocorr_min

    integer :: split_comm_y, key, color, code

    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 

    do k = 1, zsize(3)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          u_fluct3(i,j,k) = ux3(i,j,k) - uvwp_mean(i,j,1)
        enddo
      enddo
    enddo

    call transpose_z_to_y(u_fluct3,u_fluct2)
    call transpose_y_to_x(u_fluct2,u_fluct1)

    nlocs = (xlx - two*autocorr_max_sep) / autocorr_xlocs + 1
    xlocs = int(two*autocorr_max_sep/real(dx)) + 1
    
    allocate(autocorr_ml(xlocs,xsize(2),nlocs))
    allocate(autocorr_ml2(xlocs,xsize(2),nlocs))

    autocorr_ml(:,:,:) = zero

    do l = 1,xsize(3)
      do k = 1, nlocs
        autocorr_min = autocorr_x_indices(k) - int(autocorr_max_sep/real(dx))
        do j = 1, xsize(2)
          do i = 1, xlocs
            autocorr_ml(i,j,k) = autocorr_ml(i,j,k) &
                                   + u_fluct1(autocorr_x_indices(k),j,l)&
                                    *u_fluct1(autocorr_min +i,j,l)
          enddo
        enddo
      enddo
    enddo
    autocorr_ml(:,:,:) = autocorr_ml(:,:,:)/real(nz,kind=mytype)

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, code)
    key = coords(1)
    color = coords(2)

    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key,color, split_comm_y,code)

    call MPI_Allreduce(autocorr_ml,autocorr_ml2,size(autocorr_ml),&
                       real_type,MPI_SUM,split_comm_y,code)

    call MPI_Comm_free(split_comm_y,code)

    if (itempaccel==1) then
      autocorr_m = autocorr_ml2

    else
      stat_inc = real((itime-initstat2)/istatcalc+1, kind=mytype)

      autocorr_m = autocorr_m + (autocorr_ml2 - autocorr_m)/ stat_inc
    endif                           

    deallocate(autocorr_ml,autocorr_ml2)
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

    coefx_b = - zpfive / dx
    coefx_f = zpfive /dx

    coefz_b = - zpfive/dz
    coefz_f = zpfive/dz

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

