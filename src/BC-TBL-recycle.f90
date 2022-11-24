!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module tbl_recy

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: fs
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character
  real(mytype), allocatable, dimension(:,:,:) :: source
  real(mytype), dimension(:,:), allocatable :: recy_mean_z, recy_mean_t
  real(mytype), dimension(:,:), allocatable :: inlt_mean_z, inlt_mean_t
  real(mytype), dimension(:), allocatable :: u_infty_file
  real(mytype) :: delta_inlt_old
  real(mytype) :: dv_inlt, dv_recy
  real(mytype) :: delta_inlt, delta_recy

  integer :: plane_index, tbl_recy_log
  logical, parameter :: write_logs = .true.
  integer :: test_unit
  abstract interface 
  subroutine u_infty_interface(index,u_infty, u_infty_grad)
     import mytype
     integer, intent(in) :: index
     real(mytype), intent(out) :: u_infty, u_infty_grad

  end subroutine
   end interface

   procedure(u_infty_interface), pointer :: u_infty_calc


   real(mytype), dimension(:,:), allocatable :: avg_fstream
#ifdef BL_DEBG
   real(mytype), allocatable, dimension(:) :: dbg_u_inner, dbg_u_outer
   real(mytype), allocatable, dimension(:) :: dbg_v_inner, dbg_v_outer
   real(mytype), allocatable, dimension(:) :: dbg_w_inner, dbg_w_outer
   real(mytype), allocatable, dimension(:,:) :: dbg_w_fluct, dbg_v_fluct, dbg_u_fluct

   real(mytype), allocatable, dimension(:) :: dbg_u, dbg_v, dbg_w
   real(mytype), allocatable, dimension(:) :: dbg_y_plus_inlt, dbg_y_plus_recy
   real(mytype), allocatable, dimension(:) :: dbg_eta_inlt, dbg_eta_recy
   real(mytype), allocatable, dimension(:) :: dbg_u_fluct_in, dbg_v_fluct_in, dbg_w_fluct_in

   real(mytype) :: dbg_gamma
#endif

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tbl_recy, boundary_conditions_tbl_recy,&
            postprocess_tbl_recy, visu_tbl_recy, &
            visu_tbl_recy_init, momentum_forcing_tbl_recy,&
            u_infty_calc, restart_tbl_recy

contains

  subroutine init_tbl_recy (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use param , only : zptwofive
    use MPI
    use dbg_schemes


    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype), dimension(zsize(1),zsize(2)) :: u_avg_z
    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg, c, xdiff
    integer :: k,j,i,ierror,ii,is,it,code, jdx

    integer, dimension (:), allocatable :: seed
    real(mytype), dimension(nx,ny) :: u_mean, v_mean

#ifdef BL_DEBG
    real(mytype), allocatable, dimension(:,:,:) :: ux_dbg, uy_dbg, uz_dbg
    real(mytype) :: dbg_t_recy1, dbg_t_recy2
    real(mytype), allocatable, dimension(:) :: y_vals
    real(mytype), allocatable, dimension(:) :: y_in, y_out
    real(mytype), allocatable, dimension(:) :: u_out, u_in
    real(mytype), allocatable, dimension(:,:) :: fluct_out, fluct_in

    real(mytype) :: dbg_dv_inlt, dbg_dv_recy
    real(mytype) :: dbg_delta_inlt, dbg_delta_recy
    integer :: dbg_unit
    character(len=256) dbg_name

    allocate(ux_dbg(xsize(1),xsize(2),xsize(3)))
    allocate(uy_dbg(xsize(1),xsize(2),xsize(3)))
    allocate(uz_dbg(xsize(1),xsize(2),xsize(3)))

    allocate(fluct_in(xsize(2), xsize(3)))
    allocate(fluct_out(xsize(2), xsize(3)))
    allocate(y_vals(ny))

    allocate(y_in(ny), y_out(ny))
    allocate(u_out(xsize(2)), u_in(ny))
#endif

    call setup_tbl_recy

    if (iscalar==1) then

       phi1(:,:,:,:) = zptwofive !change as much as you want
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi1(:,xsize(2),:,:) = zptwofive
          endif

    endif
    ux1=zero;uy1=zero;uz1=zero 

#ifdef BL_DEBG
    if (nrank.eq.0) write(*,*) "# TBL recy debug:"
    if (nrank.eq.0) write(*,*) "# Creating debug velocities"
    do j = 1, ny
      y_vals(j) = real(j,mytype)
    enddo

    do j = 1, xsize(2)
      jdx  = xstart(2) + j -1
      ux_dbg(:,j,:) = y_vals(jdx)
      uy_dbg(:,j,:) = y_vals(jdx)
      uz_dbg(:,j,:) = y_vals(jdx)
    enddo

    if (nrank.eq.0) write(*,*) "# Resetting debug means"
   call compute_recycle_mean(ux_dbg,&
                              uy_dbg,&
                              uz_dbg,&
                              reset=.true.)

   if (nrank.eq.0) write(*,*) "# Checking means"
   call check_mean_avg("recy_t",y_vals, recy_mean_t)
   call check_mean_avg("inlt_t",y_vals, inlt_mean_t)
   call check_mean_avg("recy_z",y_vals, recy_mean_z)
   call check_mean_avg("inlt_z",y_vals, inlt_mean_z)

    if (nrank.eq.0) write(*,*) "# Writing debug scalings"
    call GetScalings(dbg_dv_inlt, dbg_dv_recy, dbg_delta_inlt, dbg_delta_recy, reset=.true.)

    call output_scalings("reset",dbg_dv_inlt, dbg_dv_recy, dbg_delta_inlt, dbg_delta_recy)

    if(nrank.eq.0) write(*,*) "# Test recycle on known field - should not change"
    
    dbg_t_recy1 = t_recy1
    dbg_t_recy2 = t_recy2
 
    t_recy1 = - one
    t_recy2 = - one

    call mean_flow_inlt_calc(dbg_u, dbg_v)
    call check_mean_flow_recy('interp')
    call output_mean_flow_recy('interp-check')


    if(nrank.eq.0) write(*,*) "# Repeating for mean update"
    call compute_recycle_mean(ux_dbg,&
                              uy_dbg,&
                              uz_dbg)

   call check_mean_avg("recy_t",y_vals, recy_mean_t)
   call check_mean_avg("inlt_t",y_vals, inlt_mean_t)
   call check_mean_avg("recy_z",y_vals, recy_mean_z)
   call check_mean_avg("inlt_z",y_vals, inlt_mean_z)
                           

        ! report average quantities
   call GetScalings(dbg_dv_inlt, dbg_dv_recy, dbg_delta_inlt, dbg_delta_recy, reset=.true.)
   call output_scalings("update",dbg_dv_inlt, dbg_dv_recy, dbg_delta_inlt, dbg_delta_recy)

   do i = 1, ny
      y_in(i)  = real(i - one    ,mytype)
      y_out(i) = real(i - zpfive ,mytype)
      u_in(i)  = real(i -one     , mytype)
   enddo

   do j = 1, xsize(2)
      jdx = xstart(2) + j -1
      fluct_in(j,:) = real(jdx -one     , mytype)
   enddo

   call mean_interp(y_in, y_out, u_in, u_out )
   call check_mean_interp('mean_interp',u_out, y_out)

   call fluct_interp(y_in, y_out, fluct_in, fluct_out)
   call check_fluct_interp('fluct_interp',fluct_out, y_out)

   ! check fluctuation interp
   do j = 1, xsize(2)
      jdx = xstart(2) + j -1
      ux1(:,j,:) = real(jdx -one     , mytype)
      uy1(:,j,:) = real(jdx -one     , mytype)
      uz1(:,j,:) = real(jdx -one     , mytype)
   enddo

   recy_mean_z(:,:) = zero

   call fluct_flow_inlt_calc(ux1, two*uy1, three*uz1,&
             dbg_u_fluct, dbg_v_fluct, dbg_w_fluct )
   
   call check_fluct_interp('u_fluct',dbg_u_fluct, y_in)
   call check_fluct_interp('v_fluct',dbg_v_fluct, two*y_in)
   call check_fluct_interp('w_fluct',dbg_w_fluct, three*y_in)

#endif

      call nickels_init(u_mean, v_mean)


      do k=1,xsize(3)
         do j=1,xsize(2)
            jdx  = xstart(2) + j -1
            do i=1,xsize(1)
               ux1(i,j,k)= u_mean(i,jdx)
               uy1(i,j,k)= v_mean(i,jdx)
               uz1(i,j,k) = zero
            enddo
         enddo
      enddo

      call compute_recycle_mean(ux1, uy1, uz1, reset=.true.)

#ifdef BL_DEBG

      call GetScalings(dbg_dv_inlt, dbg_dv_recy, dbg_delta_inlt, dbg_delta_recy, reset=.true.)
      call output_scalings("real",dbg_dv_inlt, dbg_dv_recy, dbg_delta_inlt, dbg_delta_recy)

      
      call mean_flow_inlt_calc(dbg_u, dbg_v)
      call output_mean_flow_recy('interp_real')
      
      t_recy1 = dbg_t_recy1
      t_recy2 = dbg_t_recy2


    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_tbl_recy

  subroutine setup_tbl_recy
   use MPI
   use dbg_schemes, only : abs_prec

   real(mytype) :: x, xdiff, x_read
   integer :: i, ierror, unit, ios

   allocate(recy_mean_t(3,ny))
   allocate(recy_mean_z(3,ny))
   allocate(inlt_mean_t(3,ny))
   allocate(inlt_mean_z(3,ny))

   if (plane_location.gt.xlx.and. nrank.eq.0) then
      write(*,*) "Plane location must be less than domain size"
      call MPI_Abort(MPI_COMM_WORLD,1,ierror)
    endif

    if (plane_location.le.0.and. nrank.eq.0) then
      write(*,*) "Plane location must be more than 0"
      call MPI_Abort(MPI_COMM_WORLD,1,ierror)
    endif

    do i = 1, nx
      x = real(i-1,mytype)*dx
      if (x.gt.plane_location) then
         xdiff = x-plane_location
         if (abs_prec(xdiff).gt.zpfive*dx) then
            plane_index = i - 1

         else
            plane_index = i
         endif
         exit
      endif
    enddo

    if (iaccel.eq.0) then
      u_infty_calc => zpg_BL
   else if (iaccel.eq.1) then
      u_infty_calc => tanh_BL
   else if (iaccel .eq. 2) then
      u_infty_calc => tanh_cubic_BL
   else if (iaccel .eq. 3) then
      open(newunit=unit,file=accel_file,action='read',status='old',iostat=ios)
      if (ios /= 0) then
         write(*,*) "Error reading file "//trim(adjustl(accel_file))
         call MPI_Abort(MPI_COMM_WORLD, 1,ierror)
      endif

      allocate(u_infty_file(nx))

      do i = 1, nx
         x = real(i-1,kind=mytype)*dx
         read(unit,fmt=*) &
               x_read, u_infty_file(i)

         if (abs_prec(x_read-x) > real(1e-8,kind=mytype)) then
            write(*,*) "x coordinate in file does not match case setup"
            call MPI_Abort(MPI_COMM_WORLD, 1,ierror)
         endif

      enddo
      u_infty_calc => file_BL
   else
      write(*,*) "Invalid iaccel value"
      call MPI_Abort(MPI_COMM_WORLD, 1,ierror)
   endif

   if (nrank == 0) then
      open(unit=test_unit,file='tbl_info.txt',status='replace',action='write')
      close(test_unit)
   endif
   call MPI_Barrier(MPI_COMM_WORLD,ierror)
   open(unit=test_unit,file='tbl_info.txt',status='old',action='write')

   call accel_source

#ifdef BL_DEBG   
  
   allocate(dbg_u_inner(xsize(2)), dbg_v_inner(xsize(2)))
   allocate(dbg_u_outer(xsize(2)), dbg_v_outer(xsize(2)))
   allocate(dbg_w_inner(xsize(2)), dbg_w_outer(xsize(2)))

   allocate(dbg_u(xsize(2)), dbg_v(xsize(2)), dbg_w(xsize(2)))

   allocate(dbg_u_fluct_in(xsize(2)))
   allocate(dbg_v_fluct_in(xsize(2)))
   allocate(dbg_w_fluct_in(xsize(2)))

   allocate(dbg_u_fluct(xsize(2), xsize(3)))
   allocate(dbg_v_fluct(xsize(2), xsize(3)))
   allocate(dbg_w_fluct(xsize(2), xsize(3)))

   allocate(dbg_y_plus_inlt(ny))
   allocate(dbg_y_plus_recy(ny))
   allocate(dbg_eta_inlt(ny))
   allocate(dbg_eta_recy(ny))
#endif   
  end subroutine

  subroutine restart_tbl_recy(it)
   integer, intent(in) :: it
   integer :: unit
   character(len=80) :: fn

   real(mytype), dimension(ny) :: u_inlt, v_inlt
   call setup_tbl_recy

   write(fn,'(A,I7.7)') 'tbl_recy/restart-',it

   open(newunit=unit,file=fn,status='old',action='read',access='stream')
   read(unit) recy_mean_t, inlt_mean_t, u_inlt, v_inlt,delta_inlt_old
   close(unit)
   
   if(nrank.eq.0) write(*,*) "Reading tbl_recy restart information"
  end subroutine

  subroutine write_restart_tbl_recy(it, force_write)
   use MPI
   integer, intent(in) :: it
   logical, intent(in), optional :: force_write


   integer :: unit
   integer, dimension(2) :: dims, coords
   logical, dimension(2) :: periods 
   integer :: split_comm, color, key, ierr
   character(len=80) :: fn
   integer, allocatable, dimension(:) :: ldispl,displs, recvcounts
   real(mytype), dimension(ny) :: inlt_v, inlt_u
   logical :: force_local

   if (.not.present(force_write)) then
      force_local = .false.
   else
      force_local = force_write
   endif

   if (mod(it,icheckpoint).ne.0 .and..not.force_local) return

   inlt_u = zero
   inlt_v = zero
                  
   if(nrank .eq.0) then
      write(*,*) "Writing tbl_recy restart file"
      write(fn,'(A,I7.7)') 'tbl_recy/restart-',it

      call system('mkdir -p tbl_recy')

      open(newunit=unit,file=fn,status='replace',action='write',access='stream')
      write(unit) recy_mean_t, inlt_mean_t, inlt_u, inlt_v, delta_inlt_old
      close(unit)
   endif

  end subroutine

  subroutine zpg_BL(index,u_infty, u_infty_grad)
   use param
   use variables
   use dbg_schemes, only : tanh_prec, cosh_prec

   integer, intent(in) :: index
   real(mytype), intent(out) :: u_infty, u_infty_grad

   u_infty = one
   u_infty_grad = zero

   end subroutine

   subroutine tanh_BL(index,u_infty, u_infty_grad)
      use param
      use variables
      use dbg_schemes, only : tanh_prec, cosh_prec
   
      integer, intent(in) :: index
      real(mytype), intent(out) :: u_infty, u_infty_grad
      
      real(mytype) :: x_coord, eps

      if (t<50) then
         eps = t/fifty
      else
         eps = one
      endif
      x_coord = real(index - 1, mytype) * dx
   
      u_infty = one +  eps*half *(U_ratio - one)*(&
                  tanh_prec( alpha_accel*(x_coord - accel_centre )) &
                     + one )
   
      u_infty_grad = eps*half*alpha_accel*(U_ratio-one)*cosh_prec( alpha_accel*(x_coord &
                      - accel_centre ) )**(-two)
   end subroutine

   subroutine tanh_cubic_BL(index,u_infty, u_infty_grad)
      use param
      use variables
      use var, only : t
      use dbg_schemes, only : tanh_prec, cosh_prec
      use MPI
      integer, intent(in) :: index
      real(mytype), intent(out) :: u_infty, u_infty_grad
      real(mytype) :: x, x_1, x_2, eps, a, b, c, d, inflection
      integer :: code

      if (t<50) then
         eps = t/fifty
      else
         eps = one
      endif

      x = real(index - 1, mytype) * dx

      x_1 = accel_centre + atanh(two*iaccel_thresh/(U_ratio-one)-one)/alpha_accel
      x_2 = accel_centre

      a = (U_ratio*alpha_accel*x_1 - U_ratio*alpha_accel*x_2 + two*U_ratio - alpha_accel*x_1 + alpha_accel*x_2 - two)&
         /(two*(x_1**three - three*x_1**two*x_2 + three*x_1*x_2**two - x_2**three))

      b = (-two*U_ratio*alpha_accel*x_1**two + U_ratio*alpha_accel*x_1*x_2 &
         + U_ratio*alpha_accel*x_2**two - three*U_ratio*x_1 - three*U_ratio*x_2 &
         + two*alpha_accel*x_1**two - alpha_accel*x_1*x_2 - alpha_accel*x_2**two &
            + three*x_1 + three*x_2)/(2*(x_1**three - three*x_1**two*x_2 + three*x_1*x_2**two - x_2**three))

      c = x_1*(U_ratio*alpha_accel*x_1**two + U_ratio*alpha_accel*x_1*x_2 &
         - two*U_ratio*alpha_accel*x_2**two + six*U_ratio*x_2 &
         - alpha_accel*x_1**two - alpha_accel*x_1*x_2 + 2*alpha_accel*x_2**two &
         - six*x_2)/(2*(x_1**three - three*x_1**two*x_2 + three*x_1*x_2**two - x_2**three))

      d = (-U_ratio*alpha_accel*x_1**three*x_2 + U_ratio*alpha_accel*x_1**two*x_2**two &
         + U_ratio*x_1**three - three*U_ratio*x_1**two*x_2 + alpha_accel*x_1**three*x_2 &
         - alpha_accel*x_1**two*x_2**two + x_1**three - three*x_1**two*x_2 &
         + six*x_1*x_2**2 - two*x_2**three)/(2*(x_1**three - three*x_1**two*x_2 &
         + three*x_1*x_2**two - x_2**three))

      inflection = -b / (three*a)
      if (inflection > x_1  .and. inflection < x_2) then
           if (nrank.eq.0) then
            write(*,*) "Incorrect parameters"
            call MPI_Abort(MPI_COMM_WORLD,1,code)
           endif
           
      endif
      if (x>accel_centre) then
         call tanh_BL(index,u_infty, u_infty_grad)
      
      else if (x>x_1) then
         U_infty = one - eps + eps*(a*x**three + b*x**two + c*x + d)
         U_infty_grad = eps*(three*a*x**two + two*b*x + c)
      else
         u_infty = one
         u_infty_grad = zero
      endif
   end subroutine

   subroutine file_BL(index,u_infty, u_infty_grad)
      use param, only : dx
      use var, only : t
      integer, intent(in) :: index
      real(mytype), intent(out) :: u_infty, u_infty_grad

      real(mytype) :: eps
      if (t<50) then
         eps = t/fifty
      else
         eps = one
      endif

      u_infty = u_infty_file(1) + eps*(u_infty_file(index) - u_infty_file(1))

      if (index == 1) then
         u_infty_grad = eps*(u_infty_file(2)-u_infty_file(1))/dx
      else if (index == nx) then
         u_infty_grad = eps*(u_infty_file(nx)-u_infty_file(nx-1))/dx
      else
         u_infty_grad = eps*0.5*(u_infty_file(index+1)-u_infty_file(index-1))/dx
      endif

   end subroutine

  subroutine accel_source
   use param
   use variables
   use dbg_schemes
   implicit none

   real(mytype) :: x_coord, u_infty, u_infty_grad, t_tmp
   integer :: i

   t_tmp = t
   t = 50.0
   call alloc_x(source)
   do i = 1, nx
      call u_infty_calc(i, u_infty, u_infty_grad)

      source(i,:,:) = u_infty*u_infty_grad
   enddo

   if (nrank .eq.0) then
      open(unit=13,file='u_infty.csv',action='write', status='replace')

      write(13,*) "i, x, u_infty, u_infty_grad, source"
      do i = 1, nx
         x_coord = real(i-1, mytype)*dx
         call u_infty_calc(i, u_infty, u_infty_grad)

         write(13,"(I0,*(',',g0))") i, x_coord, u_infty, u_infty_grad, source(i,1,1)
      enddo
      close(13)
   endif
   t = t_tmp

   end subroutine

  subroutine momentum_forcing_tbl_recy(dux1, duy1, duz1, ux1, uy1, uz1)
   real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
   real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

  end subroutine momentum_forcing_tbl_recy
  !********************************************************************
  subroutine boundary_conditions_tbl_recy (ux,uy,uz,phi)
   use MPI
    use navier, only : tbl_flrt
    use param , only : zero, zptwofive
    use dbg_schemes, only: cos_prec
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype), dimension(xsize(2)) :: u_mean, v_mean
    real(mytype), dimension(xsize(2),xsize(3)) :: u_fluct,&
                                                 v_fluct,&
                                                 w_fluct

    real(mytype), dimension(nx) :: v_infty

    real(mytype) :: x, y, z, u_infty, max_u, max_v, max_w, max_work
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx, dudx
    integer :: i, j, k, ierr
    

    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet
    max_work = maxval(abs(ux))
    call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    max_work = maxval(abs(uy))
    call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    max_work = maxval(abs(uz))
    call MPI_Allreduce(max_work,max_w,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (nrank ==0) &
       write(test_unit,'("Iteration: ",I0," max u v w ",g0," ",g0," ",g0," ")') itime, max_u, max_v, max_w

    call compute_recycle_mean(ux, uy, uz)
    call GetScalings(dv_inlt, dv_recy, delta_inlt, delta_recy)

    call mean_flow_inlt_calc(u_mean, v_mean)
    
    max_work = maxval(abs(u_mean))
    call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    max_work = maxval(abs(v_mean))
    call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)

    if (nrank ==0) &
       write(test_unit,'("Iteration: ",I0," max inlet mean u v ",g0," ",g0," ")') itime,max_u, max_v
    call fluct_flow_inlt_calc(ux, uy, uz, u_fluct, v_fluct, w_fluct)

    max_work = maxval(abs(u_fluct))
    call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    max_work = maxval(abs(v_fluct))
    call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    max_work = maxval(abs(w_fluct))
    call MPI_Allreduce(max_work,max_w,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    
    if (nrank ==0) &
       write(test_unit,'("Iteration: ",I0," max fluct u v w ",g0," ",g0," ",g0," ")') itime, max_u, max_v, max_w

    do k = 1,xsize(3)
      do j = 1, xsize(2)
         bxx1(j, k) = u_mean(j) + u_fluct(j,k)
         bxy1(j, k) = v_mean(j) + v_fluct(j,k)
         bxz1(j, k) = w_fluct(j,k)
      enddo
   enddo
    !INLET FOR SCALAR, TO BE CONSISTENT WITH INITIAL CONDITION
    if (iscalar==1) then
       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(1,:,:,:)=zptwofive
             if ((xstart(2)==1)) then
                phi(:,1,:,:) = one
             endif
             if ((xend(2)==ny)) THEN
                phi(:,xsize(2),:,:) = zptwofive
             endif
          enddo
       enddo
    endif

    !OUTFLOW based on a 1D convection equation

    udx=one/dx
    udy=one/dy
    udz=one/dz
    uddx=half/dx
    uddy=half/dy
    uddz=half/dz

    do k=1,xsize(3)
       do j=1,xsize(2)

          cx=ux(nx,j,k)*gdt(itr)*udx

          if (cx<zero) cx=zero
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
          if (iscalar==1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
          enddo
    enddo

    !! Bottom Boundary
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i, k) = zero
          byy1(i, k) = zero
          byz1(i, k) = zero
        enddo
      enddo
    endif

    !! Top Boundary
    if (nclyn == 2) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             call u_infty_calc(i, u_infty, dudx)

             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1,k) - dudx*(yp(ny) - yp(ny-1))
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
    endif

    !SCALAR   
    if (itimescheme/=7) then
    if (iscalar/=0) then
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi(1,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi(1,xsize(2),:,:) = phi(1,xsize(2)-1,:,:)
          endif
    endif
    endif

    !update of the flow rate (what is coming in the domain is getting out)
    call tbl_flrt(ux,uy,uz)

    return
  end subroutine boundary_conditions_tbl_recy
  subroutine v_infty_calc(uy,byy)
   real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: uy
   real(mytype), dimension(xsize(1),xsize(3)) :: byy
   real(mytype) :: a, b, c,d_inv,c_inv, u_infty, dudx, h1,h2,h3,h4
   integer :: i, k

   if (xsize(2)>2) then ! second order
      h2 = yp(ny-1) - yp(ny-2)
      h3 = yp(ny) - yp(ny-1)

      a = h3/(h2*(h2 + h3))
      b = -(h2 + h3)/(h2*h3)
      c_inv     = h3*(h2 + h3)/(h2 + two*h3)

      do i = 1, xsize(1)
         call u_infty_calc(i, u_infty, dudx)

         do k = 1, xsize(3)
            byy(i,k) = c_inv*(-dudx -b*uy(i,xsize(2)-1,k) &
                                       -a*uy(i,xsize(2)-2,k))
         enddo
      enddo
   else ! first order
      h3 = yp(ny) - yp(ny-1)

      do i = 1, xsize(1)
         call u_infty_calc(i, u_infty, dudx)

         do k = 1, xsize(3)
            byy(i,k) = uy(i,xsize(2)-1,k) - dudx*h3
         enddo
      enddo
   endif
  end subroutine
  subroutine mean_flow_inlt_calc(um_inlt, vm_inlt)
   real(mytype), dimension(:), intent(out) :: um_inlt
   real(mytype), dimension(:), intent(out) :: vm_inlt

   ! local variable declaration
   real(mytype), dimension(:), allocatable :: u_inner
   real(mytype), dimension(:), allocatable :: u_outer

   real(mytype), dimension(:), allocatable :: v_inner
   real(mytype), dimension(:), allocatable :: v_outer

   real(mytype), dimension(:), allocatable :: y_plus_inlt
   real(mytype), dimension(:), allocatable :: y_plus_recy
   real(mytype), dimension(:), allocatable :: eta_inlt
   real(mytype), dimension(:), allocatable :: eta_recy

   real(mytype) :: u_tmp, eps, y, u_inf, dudx
   real(mytype) :: gamma, w, u_infty
   integer :: j, jdx
#ifdef BL_DEBG
   character(20) :: fname
#endif
   allocate(y_plus_inlt(ny))
   allocate(y_plus_recy(ny))
   allocate(eta_inlt(ny))
   allocate(eta_recy(ny))

   allocate(u_inner(xsize(2)))
   allocate(u_outer(xsize(2)))

   allocate(v_inner(xsize(2)))
   allocate(v_outer(xsize(2)))

   gamma = dv_recy / dv_inlt ! u_tau_inlt / u_tau_recy
   u_infty = recy_mean_t(1,ny)



   eps=one

   do j = 1, ny
      if (istret==0) y=real(j-1,mytype)*dy
      if (istret/=0) y=yp(j)

      y_plus_inlt(j) = y / dv_inlt
      y_plus_recy(j) = y / dv_recy
      eta_inlt(j) = y / delta_inlt
      eta_recy(j) = y / delta_recy      
   enddo

#ifdef BL_DEBG
   dbg_gamma          = gamma
   dbg_eta_inlt(:)    = eta_inlt(:)
   dbg_eta_recy(:)    = eta_recy(:)
   dbg_y_plus_inlt(:) = y_plus_inlt(:)
   dbg_y_plus_recy(:) = y_plus_recy(:)

#endif

   ! interpolating u for inner and outer
   call mean_interp(y_plus_recy, y_plus_inlt,&
               recy_mean_t(1,:), u_inner)

   call mean_interp(eta_recy, eta_inlt,&
               recy_mean_t(1,:), u_outer)

   ! interpolating v for inner and outer
   call mean_interp(y_plus_recy, y_plus_inlt,&
               recy_mean_t(2,:), v_inner)

   call mean_interp(eta_recy, eta_inlt,&
               recy_mean_t(2,:), v_outer)                

   
   do j = 1, xsize(2)
      jdx = xstart(2) + j -1
      ! u mean inlet calculation
      w = recycleWeighting(eta_inlt(jdx))

      ! inner u
      u_tmp = gamma*u_inner(j)*(one - w) 
#ifdef BL_DEBG
      dbg_u_inner(j) = u_inner(j)
#endif
      ! outer u
      um_inlt(j) = u_tmp + w*(gamma*u_outer(j) +&
                        (one - gamma)*u_infty)
#ifdef BL_DEBG

      dbg_u_outer(j) = u_outer(j)

      dbg_v_inner(j) = v_inner(j)
      dbg_v_outer(j) = v_outer(j)
#endif
      ! weighting v with weighting function
      vm_inlt(j) = v_inner(j)*(one - w) + &
               w*v_outer(j)


   enddo

   call u_infty_inlt_calc(um_inlt,u_inf)
   call u_infty_calc(1,u_infty,dudx)
   um_inlt(:) = um_inlt(:)*u_infty/u_inf

   if (nrank ==0) &
      write(test_unit,'("Iteration: ",I0," u_infty ",g0," ",g0)') itime, u_infty, u_inf
#ifdef BL_DEBG
   write(fname,"('real-',I0)") itime

   do j = 1, xsize(2)
      dbg_u(j) = um_inlt(j)
      dbg_v(j) = vm_inlt(j)
   enddo
   call output_mean_flow_recy(trim(fname))
#endif



  end subroutine
  subroutine u_infty_inlt_calc(u_inlt,u_inf)
   use MPI
   use decomp_2d

   real(mytype), dimension(xsize(2)), intent(in) :: u_inlt
   real(mytype), intent(out) :: u_inf

   integer :: ierr, brank, rank
   integer, dimension(2) :: dims, coords
   logical, dimension(2) :: periods 

   call MPI_cart_get(DECOMP_2D_COMM_CART_X,2,dims,&
                     periods,coords,ierr)
   call MPI_cart_rank(DECOMP_2D_COMM_CART_X,[dims(1)-1,dims(2)-1],&
                     brank,ierr)
   call MPI_Comm_rank(DECOMP_2D_COMM_CART_X,rank,ierr)
   
   if (brank == rank) then
      u_inf = u_inlt(xsize(2))
   endif

   call MPI_Bcast(u_inf,1,real_type,brank,&
                  DECOMP_2D_COMM_CART_X,ierr)

  end subroutine
  
  subroutine fluct_flow_inlt_calc(ux, uy, uz, u_fluct, v_fluct, w_fluct)
   use MPI
   real(mytype), dimension(:,:,:), intent(in) :: ux, uy, uz
   real(mytype), dimension(:,:), intent(out) :: u_fluct
   real(mytype), dimension(:,:), intent(out) :: v_fluct
   real(mytype), dimension(:,:), intent(out) :: w_fluct

   real(mytype), dimension(:,:), allocatable :: u_fluct_in
   real(mytype), dimension(:,:), allocatable :: v_fluct_in
   real(mytype), dimension(:,:), allocatable :: w_fluct_in

   real(mytype), dimension(:,:), allocatable :: u_fluct_inner, u_fluct_outer
   real(mytype), dimension(:,:), allocatable :: v_fluct_inner, v_fluct_outer 
   real(mytype), dimension(:,:), allocatable :: w_fluct_inner, w_fluct_outer

   real(mytype), dimension(:), allocatable :: y_plus_inlt
   real(mytype), dimension(:), allocatable :: y_plus_recy
   real(mytype), dimension(:), allocatable :: eta_inlt
   real(mytype), dimension(:), allocatable :: eta_recy
   real(mytype) :: y, w, u_tmp, gamma, factor
   real(mytype) :: max_u, max_v, max_w, max_work

   integer :: j,k, jdx, ierr
#ifdef BL_DEBG
   character(20) :: fname
#endif
   allocate(u_fluct_in(xsize(2),xsize(3)))
   allocate(v_fluct_in(xsize(2),xsize(3)))
   allocate(w_fluct_in(xsize(2),xsize(3)))

   allocate(u_fluct_inner(xsize(2),xsize(3)))
   allocate(v_fluct_inner(xsize(2),xsize(3)))
   allocate(w_fluct_inner(xsize(2),xsize(3)))

   allocate(u_fluct_outer(xsize(2),xsize(3)))
   allocate(v_fluct_outer(xsize(2),xsize(3)))
   allocate(w_fluct_outer(xsize(2),xsize(3)))

   max_work = maxval(abs(ux))
   call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(uy))
   call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(uz))
   call MPI_Allreduce(max_work,max_w,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   if (nrank ==0) &
       write(test_unit,'("Iteration: ",I0," (fluct) max u v w ",g0," ",g0," ",g0," ")') itime, max_u, max_v, max_w

   max_work = maxval(abs(recy_mean_t(1,:)))
   call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(recy_mean_t(2,:)))
   call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)

   if (nrank ==0) &
      write(test_unit,'("Iteration: ",I0," max mean t u v ",g0," ",g0," ")') itime,max_u, max_v
   
   max_work = maxval(abs(recy_mean_z(1,:)))
   call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(recy_mean_z(2,:)))
   call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(recy_mean_z(3,:)))
   call MPI_Allreduce(max_work,max_w,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)

   if (nrank ==0) &
      write(test_unit,'("Iteration: ",I0," max mean z u v w ",g0," ",g0," ",g0," ")') itime,max_u, max_v, max_w
   
   do k = 1, xsize(3)
      do j = 1, xsize(2)
         jdx = xstart(2) + j -1
         u_fluct_in(j,k) =  ux(plane_index,j,k) - recy_mean_z(1,jdx)
         v_fluct_in(j,k) =  uy(plane_index,j,k) - recy_mean_z(2,jdx)
         w_fluct_in(j,k) =  uz(plane_index,j,k) - recy_mean_z(3,jdx)
      enddo
   enddo

   max_work = maxval(abs(u_fluct_in))
   call MPI_Allreduce(max_work,max_u,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(v_fluct_in))
   call MPI_Allreduce(max_work,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   max_work = maxval(abs(w_fluct_in))
   call MPI_Allreduce(max_work,max_w,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
   
   if (nrank ==0) &
      write(test_unit,'("Iteration: ",I0," max fluct in u v w ",g0," ",g0," ",g0," ")') itime, max_u, max_v, max_w


   allocate(y_plus_inlt(ny))
   allocate(y_plus_recy(ny))
   allocate(eta_inlt(ny))
   allocate(eta_recy(ny))

   gamma = dv_recy / dv_inlt
   
   do j = 1, ny
      if (istret==0) y=real(j-1,mytype)*dy
      if (istret/=0) y=yp(j)

      y_plus_inlt(j) = y / dv_inlt
      y_plus_recy(j) = y / dv_recy
      eta_inlt(j) = y / delta_inlt
      eta_recy(j) = y / delta_recy      
   enddo

#ifdef BL_DEBG
   dbg_gamma          = gamma
   dbg_eta_inlt(:)    = eta_inlt(:)
   dbg_eta_recy(:)    = eta_recy(:)
   dbg_y_plus_inlt(:) = y_plus_inlt(:)
   dbg_y_plus_recy(:) = y_plus_recy(:)

#endif

   ! interpolating u for inner and outer
   call fluct_interp(y_plus_recy, y_plus_inlt,&
                     u_fluct_in, u_fluct_inner)

   call fluct_interp(eta_recy, eta_inlt,&
                     u_fluct_in, u_fluct_outer)

   ! interpolating v for inner and outer
   call fluct_interp(y_plus_recy, y_plus_inlt,&
                     v_fluct_in, v_fluct_inner)

   call fluct_interp(eta_recy, eta_inlt,&
                     v_fluct_in, v_fluct_outer)

   ! interpolating w for inner and outer
   call fluct_interp(y_plus_recy, y_plus_inlt,&
                     w_fluct_in, w_fluct_inner)

   call fluct_interp(eta_recy, eta_inlt,&
                     w_fluct_in, w_fluct_outer)               

#ifdef BL_DEBG
   do j = 1, xsize(2)
      dbg_u_fluct_in(j) = u_fluct_in(j,1)
      dbg_v_fluct_in(j) = v_fluct_in(j,1)
      dbg_w_fluct_in(j) = w_fluct_in(j,1)

      dbg_u_inner(j) = u_fluct_inner(j,1)
      dbg_u_outer(j) = u_fluct_outer(j,1)
      dbg_v_inner(j) = v_fluct_inner(j,1)
      dbg_v_outer(j) = v_fluct_outer(j,1)
      dbg_w_inner(j) = w_fluct_inner(j,1)
      dbg_w_outer(j) = w_fluct_outer(j,1)
   enddo
#endif       
   
   ! if (t < spinup_time) then
   !    factor = fluct_multiply(gamma*u_fluct_inner)
   ! else
   !    factor = one
   ! endif
   factor = one
   do k = 1, xsize(3)
      do j = 1, xsize(2)
         jdx = xstart(2) + j -1
         ! u mean inlet calculation
         w = recycleWeighting(eta_inlt(jdx))

         ! u
         u_tmp = gamma*factor*u_fluct_inner(j,k)*(one - w) 
         u_fluct(j,k) = u_tmp + w*gamma*u_fluct_outer(j,k)

         ! v
         u_tmp = gamma*factor*v_fluct_inner(j,k)*(one - w) 
         v_fluct(j,k) = u_tmp + w*gamma*v_fluct_outer(j,k)

         ! w
         u_tmp = gamma*factor*w_fluct_inner(j,k)*(one - w) 
         w_fluct(j,k) = u_tmp + w*gamma*w_fluct_outer(j,k)

      enddo
   enddo

#ifdef BL_DEBG
   write(fname,"('fluct-',I0)") itime

   do j = 1, xsize(2)
      dbg_u(j) = u_fluct(j,1)
      dbg_v(j) = v_fluct(j,1)
      dbg_w(j) = w_fluct(j,1)
   enddo
   call output_fluct_flow_recy(trim(fname))
#endif   
   end subroutine

   ! subroutine v_infty_calc(v_infty)
   !    use var, only: uy1
   !    use param
   !    real(mytype), dimension(:,:) :: v_infty
   !    integer :: i, k
   !    real(mytype) :: dudx, u_infty
            
   !    if (iaccel.eq.0) then
   !       do i = 1, xsize(1)
   !          call u_infty_calc(i, u_infty, dudx)
   !          do k= 1, xsize(3)
   !             v_infty(i,k) = uy1(i, xsize(2) - 1, k) - dudx*(yp(ny)- yp(ny - 1))
   !          enddo
   !       enddo
         
   !    else
   !       call compute_disp_thick
   !       do i = 1, xsize(1)
   !          call u_infty_calc(i, u_infty, dudx)
   !          do k= 1, xsize(3)
   !             v_infty(i,k) = u_infty*disp_thick_grad(i) +&
   !                            (disp_thick(i) - yly)*dudx
   !          enddo
   !       enddo
   !    endif

   ! end subroutine

   function recycleWeighting(eta) result(w)
      use dbg_schemes, only : tanh_prec
      real(mytype), intent(in) :: eta
      real(mytype) :: w

      real(mytype), parameter :: a = four,&
                                 b = zpthree

      w = half * ( one + tanh_prec(a*(eta - b)/&
            ((one - two*b)*eta + b))/tanh_prec(a) )

   end function

   real(mytype) function fluct_multiply(u_fluct_inner)
      use MPI
      real(mytype), dimension(:,:), intent(in) :: u_fluct_inner

      integer :: j, k, ierr
      real(mytype) :: local_max, max_val

      local_max = maxval(dabs(u_fluct_inner))

      call MPI_Allreduce(local_max, max_val, 1, real_type,&
                         MPI_MAX, DECOMP_2D_COMM_CART_X, ierr)

      fluct_multiply = max(init_noise/max_val, one)                         

   end function
  !********************************************************************
  !********************************************************************
  
  subroutine nickels_init(u_mean, v_mean)
      real(mytype), dimension(:,:), intent(out) :: u_mean
      real(mytype), dimension(:, :), intent(out) :: v_mean

      real(mytype), dimension(:), allocatable :: Cf, delta_plus
      real(mytype), dimension(:), allocatable :: y_plus, u_plus

      real(mytype) :: u_tau , delta_v, u_infty, dudx
      integer(4) :: i,j

      allocate(Cf(nx), delta_plus(nx))
      allocate(y_plus(ny), u_plus(ny))

      call delta_plus_correlation(re, delta_plus)
      call Cf_correlation(re, Cf)

      do i = 1, nx
          call u_infty_calc(i,u_infty,dudx)
          call initWallUnitCalc(Cf(i), u_tau, delta_v)

          y_plus(:) = yp(:)/delta_v 
      
          call uplus_nickels(delta_plus(i), u_tau, y_plus, u_plus)
          
          do j =1, ny
              u_mean(i,j) = u_plus(j)*u_infty*u_tau
          enddo
          
      enddo

      call compute_v_mean(u_mean, v_mean)

      deallocate(y_plus, u_plus)
      deallocate(Cf, delta_plus)
   end subroutine

   subroutine uplus_nickels(delta_plus, u_tau, y_plus, u_plus)
      real(mytype), intent(in) :: delta_plus, u_tau
      real(mytype), dimension(:), intent(in) ::  y_plus
      real(mytype), dimension(:), intent(out) :: u_plus

      real(mytype) :: b, U_infty_plus, a, y, eta
      real(mytype), parameter :: yc_plus = twelve
      real(mytype), parameter :: kappa = 0.384_mytype
      real(mytype) :: u_plus_o, u_plus_s, u_plus_w
      integer :: i


      U_infty_plus = one / u_tau
      
      a = dlog(zpsix**six * delta_plus**six / yc_plus**six)/ (six*kappa)

      b = U_infty_plus - yc_plus - a

      do i = 1, ny

         y = y_plus(i)/yc_plus
         eta = y_plus(i)/delta_plus

         u_plus_o = dlog((one + (zpsix * y)**six)/(one + eta**six))/ (six*kappa)

         u_plus_s = yc_plus*( one - (one + two*y + onepfive*y*y)*dexp(-three*y))

         u_plus_w = b*(one - dexp(-five*(eta**four + eta**eight) &
                     / (one + five * eta**three)))

         u_plus(i) = u_plus_o + u_plus_s + u_plus_w
      enddo

   end subroutine

   subroutine compute_v_mean(u_mean, v_mean)
      real(mytype), dimension(:,:), intent(in) :: u_mean
      real(mytype), dimension(:,:), intent(out) :: v_mean

      real(mytype), dimension(:,:), allocatable :: dvdy
      real(mytype) :: dys, v_tmp
      integer :: i, j, k

      allocate(dvdy(nx,ny))

      do j = 1, ny
         dvdy(1, j) = -(u_mean(2,j) - u_mean(1,j)) / dx
         do i = 2, nx -1       
              dvdy(i, j) =-zpfive * (u_mean(i+1,j) - u_mean(i-1,j)) / dx
          enddo   
          dvdy(nx, j) = -(u_mean(nx,j) - u_mean(nx-1,j)) / dx
      enddo

      ! integrating dvdy using trapezoid rule
      do j = 2, ny
         do i = 1, nx 
            v_tmp = zero
            do k = 2, j
               dys = yp(k) - yp(k-1)
               v_tmp = v_tmp + zpfive * ( dvdy(i,k-1) + dvdy(i,k) )*dys
            enddo
            v_mean(i, j) = v_tmp
          enddo   
      enddo

      deallocate(dvdy)
  end subroutine

   subroutine Cf_correlation(Re_theta_target,C_f)
      real(mytype), intent(in) :: Re_theta_target
      real(mytype), dimension(:), intent(out) :: C_f
      real(mytype) :: Re_theta, x

      integer :: i

      do i = 1, nx
          x = real(i-1,mytype)*dx
          Re_theta = ( 0.015_mytype*x*Re_theta_target  + Re_theta_target**1.25_mytype)**0.8_mytype
          C_f(i) = 0.024_mytype*Re_theta**(-0.25_mytype)
      enddo

  end subroutine Cf_correlation

  subroutine delta_plus_correlation(Re_theta_target, delta_plus)
      real(mytype), intent(in) :: Re_theta_target
      real(mytype), dimension(:), intent(out) :: delta_plus
      real(mytype) :: Re_theta, x

      integer(4) :: i

      do i = 1, nx
         x = real(i-1,mytype)*dx
         Re_theta = ( 0.015_mytype*x*Re_theta_target  + Re_theta_target**1.25_mytype)**0.8_mytype
         delta_plus(i) = 1.13_mytype*Re_theta**(0.843_mytype)
      enddo
  end subroutine
  elemental subroutine initWallUnitCalc(Cf, u_tau, delta_v)
      real(mytype), intent(in) :: Cf
      real(mytype), intent(out) :: u_tau, delta_v

      u_tau = dsqrt(zpfive*Cf)
      delta_v = one/(u_tau*re)

  end subroutine

  subroutine compute_recycle_mean( ux, uy, uz, reset)
   use decomp_2d
   use MPI
   use iso_fortran_env, only: output_unit

   real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uy, uz
   logical, optional, intent(in) :: reset

   real(mytype), allocatable, dimension(:,:) :: recy_mean_z_local
   real(mytype), allocatable, dimension(:,:) :: inlt_mean_z_local

   real(mytype), allocatable, dimension(:,:) :: recy_mean_z_ls
   real(mytype), allocatable, dimension(:,:) :: inlt_mean_z_ls


   integer :: k, j, ierr,  color, key, split_comm_y, split_comm_z, asize
   integer, dimension(2) :: dims, coords
   logical, dimension(2) :: periods 
   integer, allocatable, dimension(:) :: ldispl,displs, recvcounts
   logical :: reset_local
   real(mytype) :: T_period, dtdivT, u_infty, dudx

   allocate(recy_mean_z_local(3, xsize(2)))
   allocate(inlt_mean_z_local(3, xsize(2)))

   allocate(recy_mean_z_ls(3, xsize(2)))
   allocate(inlt_mean_z_ls(3, xsize(2)))

   recy_mean_z_local(:,:) = zero
   inlt_mean_z_local(:,:) = zero

   call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierr)

   allocate(displs(dims(1)), recvcounts(dims(1)))
   allocate(ldispl(dims(1)))
   
   color = coords(1)
   key = coords(2)

   call MPI_Comm_split(DECOMP_2D_COMM_CART_X, color, key, split_comm_z,ierr)
   call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key, color, split_comm_y,ierr)

   do k = 1, xsize(3)
      do j = 1, xsize(2)

         inlt_mean_z_local(1,j) = inlt_mean_z_local(1,j) + ux(1,j, k)/nz
         inlt_mean_z_local(2,j) = inlt_mean_z_local(2,j) + uy(1,j, k)/nz
         inlt_mean_z_local(3,j) = inlt_mean_z_local(3,j) + uz(1,j, k)/nz

         recy_mean_z_local(1,j) = recy_mean_z_local(1,j) + ux(plane_index,j, k)/nz
         recy_mean_z_local(2,j) = recy_mean_z_local(2,j) + uy(plane_index,j, k)/nz
         recy_mean_z_local(3,j) = recy_mean_z_local(3,j) + uz(plane_index,j, k)/nz

      enddo
   enddo

   asize = xsize(2)*3
   call MPI_Allreduce(inlt_mean_z_local, inlt_mean_z_ls, asize, real_type,&
                      MPI_SUM, split_comm_z, ierr)

   call MPI_Allreduce(recy_mean_z_local, recy_mean_z_ls, asize, real_type,&
                      MPI_SUM, split_comm_z, ierr) 

    
   call MPI_Allgather(xstart(2), 1, MPI_INTEGER, ldispl,&
                     1, MPI_INTEGER,split_comm_y,ierr)

   displs = 3*(ldispl - 1)
   call MPI_Allgather(xsize(2)*3, 1, MPI_INTEGER, recvcounts,&
               1, MPI_INTEGER,split_comm_y,ierr)                     

   call MPI_Allgatherv(inlt_mean_z_ls, xsize(2)*3, real_type, &
                       inlt_mean_z,recvcounts,displs,real_type,&
                       split_comm_y,ierr)

   call MPI_Allgatherv(recy_mean_z_ls, xsize(2)*3, real_type, &
                       recy_mean_z,recvcounts,displs,real_type,&
                       split_comm_y,ierr)                       

   if (.not. present(reset)) reset_local = .false.
   if (present(reset)) reset_local = reset

   if (reset_local) then
      inlt_mean_t(:,:) = inlt_mean_z(:,:)
      recy_mean_t(:,:) = recy_mean_z(:,:)

   else

      if (itime .lt. t_avg1.and.itime.lt. t_avg2) then
         T_period = 1

      elseif (itime.lt. t_avg2) then
         T_period = 50

      else
         T_period =50 + itime - t_avg2

      endif

      dtdivT = dt / T_period
      inlt_mean_t(:,:) = dtdivT *  inlt_mean_z(:,:) + &
                           ( one - dtdivT )*inlt_mean_t(:,:)

      recy_mean_t(:,:) = (dtdivT) *  recy_mean_z(:,:) + &
                           ( one - dtdivT )*recy_mean_t(:,:)                           
   endif

   call u_infty_calc(1,u_infty,dudx)
   inlt_mean_t(1,:) = inlt_mean_t(1,:)*u_infty/inlt_mean_t(1,ny)


   call MPI_Comm_free(split_comm_y,ierr)
   call MPI_Comm_free(split_comm_z,ierr)
   deallocate(inlt_mean_z_local)
   deallocate(recy_mean_z_local)
   deallocate(recy_mean_z_ls)
   deallocate(inlt_mean_z_ls)

   deallocate(displs,recvcounts,ldispl)
  end subroutine

  subroutine GetScalings(delta_v_inlt, delta_v_recy, delta_i, delta_r, reset)
   use dbg_schemes
   use param, only : irestart
   real(mytype), intent(out) :: delta_v_inlt, delta_v_recy, delta_i, delta_r
   logical, intent(in), optional :: reset
   real(mytype) :: delta_meas

   logical, save :: first_call = .true.
   real(mytype) :: dyy, u_infty_recy, u_infty_inlt, u_thresh
   real(mytype) :: theta_inlt, theta_recy, int_inlt
   real(mytype) :: int_recy, mid_u, u_tau
   real(mytype), parameter :: alp = 0.3
   integer :: i,j, unit, pos, nreads
   logical :: reset_local

   ! compute theta
   theta_inlt = zero
   theta_recy = zero

   u_infty_inlt = inlt_mean_t(1,ny)
   u_infty_recy = recy_mean_t(1,ny)

   do j = 2, ny
      mid_u = zpfive * (inlt_mean_t(1,j) + inlt_mean_t(1,j-1))/u_infty_inlt
      int_inlt = (mid_u)*( one - mid_u)

      mid_u = zpfive * (recy_mean_t(1,j) + recy_mean_t(1,j-1))/u_infty_recy
      int_recy = mid_u*( one - mid_u)

      theta_inlt = theta_inlt + int_inlt*(yp(j) - yp(j-1))
      theta_recy = theta_recy + int_recy*(yp(j) - yp(j-1))
   enddo

   ! compute delta
   do j = 1, ny
      u_thresh = 0.99_mytype*u_infty_recy
      if (recy_mean_t(1,j) .gt. u_thresh) then
         delta_r = yp(j-1) + (yp(j) - yp(j-1))*(u_thresh - recy_mean_t(1,j-1))
         exit
      endif
   enddo

   do j = 1, ny
      u_thresh = 0.99_mytype*u_infty_inlt
      if (inlt_mean_t(1,j) .gt. u_thresh) then
         delta_meas = yp(j-1) + (yp(j) - yp(j-1))*(u_thresh - inlt_mean_t(1,j-1))
         exit
      endif
   enddo

   if (irestart.ne.0) first_call = .false.
   if (first_call) then
      delta_i = delta_meas
      
      if (.not. present(reset)) reset_local = .false.
      if (present(reset)) reset_local = reset

      if (reset_local) then
         first_call = .true.
      else
         first_call = .false.
      endif
   else if (itr .eq. 1) then
      delta_i = delta_inlt_old + alp*( one - theta_inlt)*delta_inlt_old
      if (abs_prec(delta_i - delta_meas) > one) then
         delta_i = delta_meas + sign(one,delta_i - delta_meas)
      endif
   else 
      delta_i = delta_inlt_old
   endif
   
   ! compute friction velocity
   if (istret==0) dyy = dy
   if (istret/=0) dyy = yp(2) - yp(1)

   u_tau = sqrt_prec((recy_mean_t(1,2) - recy_mean_t(1,1))/dyy/re)
   delta_v_recy = one / (re*u_tau)

   u_tau = u_tau * (theta_recy / theta_inlt)** 0.125_mytype
   delta_v_inlt = one / (re*u_tau)

   if (write_logs.and.nrank.eq.0) then
      if (itime==1) then
         open(newunit=tbl_recy_log,file='tbl_recy.log',status='replace',action='write')
         write(tbl_recy_log,"(A,*(',',A))") "itime","theta_inlt", "theta_recy",&
                                 "delta_r","delta_meas","delta_old","delta_i",&
                                 "u_tau_inlt","delta_v_recy","delta_v_inlt"
      
      else if (itime == ifirst) then
         open(newunit=tbl_recy_log,file='tbl_recy.log',status='old',action='readwrite')
         rewind(tbl_recy_log)
         nreads = itime/ilist + 1
         if (ilist==1) nreads = itime/ilist - 1
         do i = 1, nreads
            read(tbl_recy_log,*)
         enddo
      endif

      if ((mod(itime,ilist)==0.or.itime==1) .and. itr.eq.1) then
            write(tbl_recy_log,"(I0,*(',',g0))") itime,theta_inlt,theta_recy,delta_r,&
            delta_meas,delta_inlt_old, delta_i, u_tau, delta_v_recy,&
            delta_v_inlt
      endif
      if (itime == ilast) close(tbl_recy_log)
   endif

   delta_inlt_old = delta_i
  end subroutine

  subroutine mean_interp(y_in, y_out, u_in, u_out)
   real(mytype), dimension(:),intent(in) :: y_in, y_out, u_in
   real(mytype), dimension(:), intent(out) :: u_out

   integer :: j, k, j_lower, j_upper, idx

   do j = 1, xsize(2)
      idx = xstart(2) + j -1
      do k = 1, ny
         if (y_in(k).gt.y_out(idx )) then
            j_lower = k - 1
            j_upper = k

            u_out(j) = u_in(j_lower) + &
                     (y_out(idx) - y_in(j_lower))*&
                     (u_in(j_upper) - u_in(j_lower))
            exit
         elseif (k == ny) then
            u_out(j) = u_in(ny)
            exit

         endif

      enddo
   enddo

  end subroutine mean_interp

  subroutine fluct_interp(y_in, y_out, u_in, u_out)
   use iso_fortran_env, only: output_unit
   use decomp_2d
   use MPI

   real(mytype), dimension(:),intent(in) :: y_in, y_out
   real(mytype), dimension(:,:),intent(in) :: u_in
   real(mytype), dimension(:,:), intent(out) :: u_out
   integer :: i, j, k

   integer :: color, key
   integer, dimension(2) :: dims, coords
   logical, dimension(2) :: periods 

   integer :: sendcount, ierr, split_comm_y
   real(mytype), allocatable, dimension(:,:) :: u_in_local
   integer, allocatable, dimension(:) :: ldispl,displs,recvcounts

   integer :: j_lower, j_upper, idx, send_glb
   real(mytype), allocatable, dimension(:,:) :: u_in_switched

   call MPI_Cart_get(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierr)
   
   key = coords(1)
   color = coords(2)

   call MPI_Comm_split(DECOMP_2D_COMM_CART_X, color, key, split_comm_y,ierr)

   allocate(u_in_local(xsize(3),ny))
   allocate(displs(dims(1)))
   allocate(recvcounts(dims(1)))
   allocate(u_in_switched(xsize(3),xsize(2)))
   allocate(ldispl(dims(1)))

   call MPI_Allgather(xstart(2), 1, MPI_INTEGER, ldispl,&
                     1, MPI_INTEGER,split_comm_y,ierr)  
   call MPI_Allgather(xsize(2)*xsize(3), 1, MPI_INTEGER, recvcounts,&
                     1, MPI_INTEGER,split_comm_y,ierr)                        
   
   sendcount = xsize(2)*xsize(3)
   do j = 1, xsize(2)
      do k = 1, xsize(3)
         u_in_switched(k,j) = u_in(j,k)
      enddo
   enddo

   displs = xsize(3)*(ldispl - 1)
   call MPI_Allgatherv(u_in_switched, sendcount, real_type,u_in_local,&
                        recvcounts,displs,real_type,split_comm_y,&
                        ierr)

   do i = 1, xsize(3)
      do j = 1, xsize(2)
         idx = xstart(2) + j -1
         do k = 1, ny
            if (y_in(k).gt.y_out(idx )) then
               j_lower = k - 1
               j_upper = k

               u_out(j,i) = u_in_local(i,j_lower) + &
                        (y_out(idx) - y_in(j_lower)) * &
                        (u_in_local(i,j_upper) - u_in_local(i,j_lower))
               exit
            elseif (k == ny) then
               u_out(j,i) = zero
               exit

            endif

         enddo
      enddo
   enddo

   call MPI_Comm_free(split_comm_y,ierr)

   deallocate(u_in_switched)
   deallocate(u_in_local)
  end subroutine fluct_interp

   ! subroutine initialise_avg
   !    use var, only : ux1
   !    real(mytype), dimension(:,:), allocatable :: avg_z

   !    allocate(avg_z(xsize(1),xsize(2)))
   !    allocate(avg_fstream(nx,ny))
   !    allocate(disp_thick(nx),disp_thick_grad(nx))

   !    call compute_avg_z(ux1,avg_z)
   !    call gather_avg(avg_z,avg_fstream)
   ! end subroutine

   ! subroutine gather_avg(avg_local, avg)
   !    use decomp_2d
   !    use MPI
   !    real(mytype), dimension(:,:), intent(in) :: avg_local
   !    real(mytype), dimension(:,:), intent(out) :: avg 

   !    integer :: split_comm_y
   !    integer :: sendcount, ierr, i
   !    integer :: color, key
   !    integer,dimension(:), allocatable :: recv_counts
   !    integer,dimension(:), allocatable :: displs
   !    logical :: periods(2)
   !    integer :: coords(2), dims(2)

   !    call MPI_Cart_get(DECOMP_2D_COMM_CART_X,2,dims,periods,coords,ierr)

   !    color = coords(2)
   !    key = coords(1)
   !    call MPI_comm_split(DECOMP_2D_COMM_CART_X,color, key, split_comm_y,ierr)

   !    allocate(recv_counts(dims(1)))
   !    allocate(displs(dims(1)))

   !    sendcount = xsize(1)*xsize(2)

   !    call MPI_Allgather(sendcount,1,MPI_INTEGER,recv_counts,&
   !                   1,MPI_INTEGER,split_comm_y,ierr)
   !    do i = 2, dims(1)
   !       displs(i) = sum(recv_counts(1:i-1))
   !    enddo
   !    displs(1) = 0

   !    call MPI_Allgatherv(avg_local,sendcount,real_type,avg,&
   !                         recv_counts,displs,real_type,&
   !                         split_comm_y,ierr)

   !    call MPI_comm_free(split_comm_y,ierr)
   ! end subroutine
   ! subroutine compute_avg_z(ux,avg_z)
   !    use decomp_2d
   !    use MPI

   !    real(mytype), dimension(:,:,:), intent(in) :: ux
   !    real(mytype), dimension(:,:), intent(out) :: avg_z
   !    real(mytype), dimension(:,:), allocatable :: avg_z_local
   !    integer :: i,j,k
   !    integer :: color, key, ierr
   !    integer :: split_comm_z, size
   !    logical :: periods(2)
   !    integer :: coords(2), dims(2)

   !    allocate(avg_z_local(xsize(1), xsize(2)))

   !    do i = 1, xsize(1)
   !       do j = 1, xsize(2)
   !          avg_z_local(i,j) = zero
   !          do k = 1, xsize(3)
   !             avg_z_local(i,j) = avg_z_local(i,j)  + ux(i,j,k)/nz
   !          enddo
   !       enddo
   !    enddo
      
   !    call MPI_Cart_get(DECOMP_2D_COMM_CART_X,2,dims,periods,coords,ierr)
   !    color = coords(1)
   !    key = coords(2)
   !    call MPI_comm_split(DECOMP_2D_COMM_CART_X,color, key, split_comm_z,ierr)
   !    size = xsize(1)*xsize(2)
   !    call MPI_Allreduce(avg_z_local, avg_z, size, real_type,&
   !                       MPI_SUM,split_comm_z,ierr)

   !    call MPI_Comm_free(split_comm_z,ierr)

   ! end subroutine

   ! subroutine compute_disp_thick
   !    use var, only : ux1
   !    use param

   !    real(mytype), dimension(:,:), allocatable :: avg_z
   !    real(mytype), dimension(:,:), allocatable :: avg_z_full

   !    real(mytype), dimension(:), allocatable :: u_infty
   !    real(mytype) :: t_period, integrand, midp, eps
   !    integer :: i, j

   !    allocate(avg_z_full(nx, ny))
   !    allocate(avg_z(xsize(1),xsize(2)))
   !    allocate(u_infty(xsize(1)))

   !    call compute_avg_z(ux1, avg_z)
   !    call gather_avg(avg_z, avg_z_full)

   !    ! compute time average
   !    if (t < t_avg_fstream) then
   !       t_period = ten
   !    else
   !       t_period = 100.0_mytype
   !    endif

   !    eps = dt/t_period
   !    avg_fstream(:,:) = eps*avg_z_full(:,:) + (one - eps)*avg_fstream(:,:)

   !    ! compute displacement thickness global array

   !    do i = 1, nx
   !       disp_thick(i) = zero
   !       do j = 1, ny -1
   !          midp = zpfive*(avg_fstream(i,j) + avg_fstream(i,j+1))
   !          integrand =  one - midp/avg_fstream(i,ny)
   !          dy = yp(j+1) - yp(j)
   !          disp_thick(i) = disp_thick(i) + integrand*dy
   !       enddo
   !    enddo

   !    disp_thick_grad(1) = (disp_thick(2) - disp_thick(1))/dx
   !    do i = 2, nx -1
   !       disp_thick_grad(i) = zpfive*(disp_thick(i-1) - disp_thick(i+1))/dx
   !    enddo
   !    disp_thick_grad(nx) = (disp_thick(nx) - disp_thick(nx-1))/dx
   !    if (nrank .eq.0) write(*,*) disp_thick
   ! end subroutine 
  !############################################################################
  subroutine postprocess_tbl_recy(ux1,uy1,uz1,ep1)

    USE MPI
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    character(len=30) :: filename
    integer :: nfil, i

   call write_restart_tbl_recy(itime)

  end subroutine postprocess_tbl_recy

  subroutine visu_tbl_recy_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "mpi_rank", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_tbl_recy_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_tbl
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs TBL-specific visualization
  !!
  !############################################################################
  subroutine visu_tbl_recy(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)

    call write_field(di1, ".", "vort", trim(num), flush=.true.) ! Reusing temporary array, force flush

    ! write mpi ranks
    di1(:,:,:) = real(nrank,kind=mytype)
    call write_field(di1, ".", "mpi_rank", trim(num), flush=.true.) ! Reusing temporary array, force flush

  end subroutine visu_tbl_recy

#ifdef BL_DEBG
   subroutine check_mean_avg(name,y_vals, mean_array)
      use MPI
      use decomp_2d
      use dbg_schemes, only: abs_prec

      character(len=*) :: name
      real(mytype), dimension(:) :: y_vals
      real(mytype), dimension(:,:) :: mean_array

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      real(mytype) :: dbg_val
      integer :: ierror, i,j,k

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         close(13)
      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      
      do i = 0 , nproc
         if (i == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')

            do j = 1, ny
               do k = 1, 3
                  dbg_val = abs_prec(mean_array(k,j) - y_vals(j))
                  
                  if (( dbg_val/y_vals(j) > 1e-10_mytype)) then
                     write(13,"('Rank:',I3,' Failure at position (reset) ',I3,I3,' Value:',F17.8, F17.8, ' X COMM coords:',I3,I3)") nrank, j, k, mean_array(k,j), y_vals(j),coords(1), coords(2)
                  endif
               enddo
            enddo
            close(13)
         endif
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo

   end subroutine

   subroutine output_scalings(name, dbg_dv_inlt1, dbg_dv_recy1, dbg_delta_inlt1, dbg_delta_recy1)
      use MPI
      use decomp_2d
      character(len=*) :: name
      real(mytype), intent(in) :: dbg_dv_inlt1, dbg_dv_recy1,&
                                  dbg_delta_inlt1, dbg_delta_recy1

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      integer :: ierror, j

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         close(13)
      endif
      
      do j = 0, nproc
         if (j == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')
            write(13,"('Reset Rank: ',I5,' dv_inlt: ',F17.8,' dv_recy',F17.8, ' delta_inlt',F17.8,' delta_recy',F17.8)")&
                        nrank, dbg_dv_inlt1, dbg_dv_recy1, dbg_delta_inlt1, dbg_delta_recy1
            close(13)
         endif
   
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo

   end subroutine

   subroutine check_mean_interp(name, u, u_exp)
      use MPI
      use decomp_2d
      use dbg_schemes, only: abs_prec

      character(len=*), intent(in) :: name
      real(mytype), dimension(:), intent(in) :: u, u_exp

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      integer :: ierror, i,j
      real(mytype) :: dbg_val

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         close(13)
      endif
      
      do j = 0, nproc
         if (j == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')
            do  i = 1, xsize(2)
               dbg_val = u_exp(xstart(2)+ i -1) - u(i)
               if (abs_prec(dbg_val) > 1e-10_mytype) then
                  write(13,"('Rank: ',I5,' J (global): ',I5,' J (local): ',I5,' u_exp',F17.8, ' u',F17.8)")&
                        nrank, xstart(2)+ i -1, i, u_exp(xstart(2)+ i -1), u(i)
               endif
            enddo
            
            close(13)
         endif
   
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo

   end subroutine

   subroutine check_mean_flow_recy(name)
      use MPI
      use decomp_2d
      use dbg_schemes, only: abs_prec

      character(len=*), intent(in) :: name

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      integer :: ierror, i,j, idx
      real(mytype) :: dbg_val

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         close(13)
      endif
      
      do j = 0, nproc
         if (j == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')
            do  i = 1, xsize(2)
               idx = xstart(2) + i -1
               dbg_val = dbg_u(i) - recy_mean_t(1,idx)

               if (abs_prec(dbg_val/recy_mean_t(1,idx)) > 1e-7_mytype) then
                  write(13,"('Rank: ',I5,' J (global): ',I5,' J (local): ',I5,' interp',F17.8, ' mean',F17.8, ' diff',ES17.8)")&
                        nrank, idx, i, dbg_u(i), recy_mean_t(1,idx), abs_prec(dbg_val/recy_mean_t(1,idx))
               endif
            enddo
            
            close(13)
         endif
   
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo
   end subroutine
   subroutine output_mean_flow_recy(name)
      use MPI
      use decomp_2d

      character(len=*), intent(in) :: name

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      integer :: ierror, i,j , idx
      real(mytype) :: w

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         write(13, *) "j, y, gamma, w, y_plus_inlt, y_plus_recy, eta_inlt, eta_recy, u_mean, u_inner,"//&
                  " u_outer, u, v_mean, v_inner, v_outer, v"
         close(13)
      endif
      
      do j = 0, nproc, dims(2)
         if (j == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')
            do  i = 1, xsize(2)
               idx = xstart(2) + i - 1
               w = recycleWeighting(dbg_eta_inlt(idx))
               write(13,"(I0, *(',',g0))") idx, yp(idx), dbg_gamma, w, dbg_y_plus_inlt(idx), dbg_y_plus_recy(idx),&
                                       dbg_eta_inlt(idx), dbg_eta_recy(idx), recy_mean_t(1,idx), &
                                       dbg_u_inner(i), dbg_u_outer(i), dbg_u(i), recy_mean_t(2,idx),&
                                       dbg_v_inner(i), dbg_v_outer(i), dbg_v(i)
            enddo
            
            close(13)
         endif
   
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo
   end subroutine

   subroutine output_fluct_flow_recy(name)
      use MPI
      use decomp_2d

      character(len=*), intent(in) :: name

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      integer :: ierror, i,j, idx
      real(mytype) :: w

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         write(13, *) "jdx, j, y, gamma, w, y_plus_inlt, y_plus_recy, eta_inlt, eta_recy,u_in,  u_mean, u_inner,"//&
                  " u_outer, u, v_in, v_mean, v_inner, v_outer, v, w_in,  w_mean , w_inner, w_outer, w"
         close(13)
      endif
      
      do j = 0, nproc, dims(2)
         if (j == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')
            do  i = 1, xsize(2)
               idx = xstart(2) + i - 1
               w = recycleWeighting(dbg_eta_inlt(idx))
               write(13,"(I0,',', I0, *(',',g0))") idx, i, yp(idx), dbg_gamma, w, dbg_y_plus_inlt(idx), dbg_y_plus_recy(idx),&
                                       dbg_eta_inlt(idx), dbg_eta_recy(idx),  &
                                       dbg_u_fluct_in(i), recy_mean_z(1,idx), dbg_u_inner(i), dbg_u_outer(i), dbg_u(i), &
                                       dbg_v_fluct_in(i), recy_mean_z(2,idx), dbg_v_inner(i), dbg_v_outer(i), dbg_v(i),&
                                       dbg_w_fluct_in(i), recy_mean_z(3,idx), dbg_w_inner(i), dbg_w_outer(i), dbg_w(i)
            enddo
            
            close(13)
         endif
   
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo
   end subroutine   

   subroutine check_fluct_interp(name, u, u_exp)
      use MPI
      use decomp_2d
      use dbg_schemes, only: abs_prec

      character(len=*), intent(in) :: name
      real(mytype), dimension(:,:), intent(in) :: u
      real(mytype), dimension(:), intent(in) :: u_exp

      integer, dimension(2) :: dims, coords
      logical, dimension(2) :: periods 
      integer :: ierror, i,j,k
      real(mytype) :: dbg_val

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierror)

      if (nrank.eq.0) then
         open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',status='replace')
         close(13)
      endif
      
      do j = 0, nproc
         if (j == nrank) then
            open(13, file='tbl_recy-dbg-'//trim(name)//'.log',action='write',access='append',status='old')
            do k = 1, xsize(3)
               do  i = 1, xsize(2)
                  dbg_val = u_exp(xstart(2)+ i -1) - u(i,k)
                  if (abs_prec(dbg_val) > 1e-10_mytype) then
                     write(13,"('Rank: ',I5,' J (global): ',I5,' J (local): ',I5,' K (global): ',I5,' K (local): ',I5,' u_exp',F17.8, ' u',F17.8)")&
                           nrank, xstart(2)+ i -1, i, xstart(3)+ k -1, k, u_exp(xstart(2)+ i -1), u(i,k)
                  endif
               enddo
            enddo
            
            close(13)
         endif
   
         call MPI_Barrier(MPI_COMM_WORLD,ierror)
      enddo

   end subroutine

#endif
end module tbl_recy
