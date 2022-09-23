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
  real(mytype), allocatable, dimension(:) :: um_inlt_ini, vm_inlt_ini
  real(mytype), allocatable, dimension(:,:,:) :: source
  real(mytype), dimension(:,:), allocatable :: recy_mean_z, recy_mean_t
  real(mytype), dimension(:,:), allocatable :: inlt_mean_z, inlt_mean_t
  real(mytype), dimension(:), allocatable :: disp_thick, disp_thick_grad
  real(mytype) :: delta_inlt_old
  real(mytype) :: dv_inlt, dv_recy
  real(mytype) :: delta_inlt, delta_recy

  integer :: plane_index
  logical, parameter :: write_logs = .true.
  
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

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg, c, xdiff
    integer :: k,j,i,ierror,ii,is,it,code, jdx

    integer, dimension (:), allocatable :: seed
    real(mytype), dimension(:,:), allocatable :: u_mean, v_mean

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

    allocate(u_mean(nx,ny), v_mean(nx,ny))

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

   if (irestart.ne.1) then       
      call nickels_init(u_mean, v_mean)
      
      do j = 1, xsize(2)
         jdx  = xstart(2) + j -1
         um_inlt_ini(j) = u_mean(1,jdx)
         vm_inlt_ini(j) = v_mean(1,jdx)
      enddo

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

      ! reseting inly_ini variables
      do j = 1, xsize(2)
         jdx  = xstart(2) + j -1
         um_inlt_ini(j) = u_mean(1,jdx)
         vm_inlt_ini(j) = v_mean(1,jdx)
      enddo

#endif
    else
    endif



#ifdef BL_DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    deallocate(u_mean, v_mean)
    return
  end subroutine init_tbl_recy

  subroutine setup_tbl_recy
   use MPI
   use dbg_schemes, only : abs_prec

   real(mytype) :: x, xdiff
   integer :: i, ierror

   allocate(recy_mean_t(3,ny))
   allocate(recy_mean_z(3,ny))
   allocate(inlt_mean_t(3,ny))
   allocate(inlt_mean_z(3,ny))

   allocate(um_inlt_ini(xsize(2)), vm_inlt_ini(xsize(2)))

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
   else
      write(*,*) "Invalid iaccel value"
      call MPI_Abort(MPI_COMM_WORLD, 1,ierror)
   endif

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

   um_inlt_ini(:) = u_inlt(xstart(2):xend(2))
   vm_inlt_ini(:) = v_inlt(xstart(2):xend(2))
   
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

   call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierr)

   allocate(displs(dims(1)), recvcounts(dims(1)))
   allocate(ldispl(dims(1)))
   
   key = coords(1)
   color = coords(2)


   call MPI_Comm_split(DECOMP_2D_COMM_CART_X,color, key,split_comm,ierr)
   if (color.ne.0) then
      call MPI_Comm_free(split_comm,ierr)
      return
   endif

   call MPI_Allgather(xstart(2), 1, MPI_INTEGER, ldispl,&
   1, MPI_INTEGER,split_comm,ierr)

   displs = ldispl - 1
   call MPI_Allgather(xsize(2), 1, MPI_INTEGER, recvcounts,&
   1, MPI_INTEGER,split_comm,ierr) 

   call MPI_Gatherv(um_inlt_ini,xsize(2),real_type,inlt_u,recvcounts,&
                  displs,real_type,0,split_comm,ierr)

   call MPI_Gatherv(vm_inlt_ini,xsize(2),real_type,inlt_v,recvcounts,&
                  displs,real_type,0,split_comm,ierr)
                  
   if(nrank .eq.0) then
      write(*,*) "Writing tbl_recy restart file"
      write(fn,'(A,I7.7)') 'tbl_recy/restart-',it

      call system('mkdir -p tbl_recy')

      open(newunit=unit,file=fn,status='replace',action='write',access='stream')
      write(unit) recy_mean_t, inlt_mean_t, inlt_u, inlt_v, delta_inlt_old
      close(unit)
   endif
   call MPI_Comm_free(split_comm,ierr)

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
  subroutine accel_source
   use param
   use variables
   use dbg_schemes
   implicit none

   real(mytype) :: x_coord, u_infty, u_infty_grad
   integer :: i

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
   end subroutine

  subroutine momentum_forcing_tbl_recy(dux1, duy1, duz1, ux1, uy1, uz1)
   real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
   real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

   if (iaccel.eq.1) then
      dux1(:,:,:,1) = dux1(:,:,:,1) + source(:,:,:)
   endif

  end subroutine momentum_forcing_tbl_recy
  !********************************************************************
  subroutine boundary_conditions_tbl_recy (ux,uy,uz,phi)
   use MPI
    use navier, only : tbl_flrt
    use param , only : zero, zptwofive
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype), dimension(:), allocatable :: u_mean, v_mean
    real(mytype), dimension(:,:), allocatable :: u_fluct,&
                                                 v_fluct,&
                                                 w_fluct

    real(mytype), dimension(:,:), allocatable :: v_infty

    real(mytype) :: x, y, z, u_infty
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx, dudx
    integer :: i, j, k
    
    allocate(u_mean(xsize(2)))
    allocate(v_mean(xsize(2)))

    allocate(u_fluct(xsize(2),xsize(3)))
    allocate(v_fluct(xsize(2),xsize(3)))
    allocate(w_fluct(xsize(2),xsize(3)))
    allocate(v_infty(xsize(1), xsize(3)))

    !INFLOW with an update of bxx1, byy1 and bzz1 at the inlet
    
    call compute_recycle_mean(ux, uy, uz)
    call GetScalings(dv_inlt, dv_recy, delta_inlt, delta_recy)

    call mean_flow_inlt_calc(u_mean, v_mean)
    call fluct_flow_inlt_calc(ux, uy, uz, u_fluct, v_fluct, w_fluct)

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

   !  call v_infty_calc(v_infty)

    !! Top Boundary
    if (nclyn == 2) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
            call u_infty_calc(i, u_infty, dudx)

             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1, k) - dudx*(yp(ny)- yp(ny - 1))
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

   real(mytype) :: u_tmp, eps, y
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



   if (t < t_recy1) then
      eps = zero
   else if (t < t_recy2) then
      eps = (t - t_recy1)/ (t_recy2 - t_recy1)
   else
      eps = one
   endif

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
      u_tmp = u_tmp + w*(gamma*u_outer(j) +&
                        (one - gamma)*u_infty)
#ifdef BL_DEBG

      dbg_u_outer(j) = u_outer(j)

      dbg_v_inner(j) = v_inner(j)
      dbg_v_outer(j) = v_outer(j)
#endif
      ! weighted with value from previous time step
      um_inlt(j) = eps*u_tmp + (one - eps)*um_inlt_ini(j)
      um_inlt_ini(j) = um_inlt(j)
      ! weighting v with weighting function
      u_tmp = v_inner(j)*(one - w) + &
               w*v_outer(j)

      ! weighted with value from previous time step
      vm_inlt(j) = eps*u_tmp + (one - eps)*vm_inlt_ini(j)
      vm_inlt_ini(j) = vm_inlt(j)

   enddo

#ifdef BL_DEBG
   write(fname,"('real-',I0)") itime

   do j = 1, xsize(2)
      dbg_u(j) = um_inlt(j)
      dbg_v(j) = vm_inlt(j)
   enddo
   call output_mean_flow_recy(trim(fname))
#endif



  end subroutine

  subroutine fluct_flow_inlt_calc(ux, uy, uz, u_fluct, v_fluct, w_fluct)
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
   integer :: j,k, jdx
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

   do k = 1, xsize(3)
      do j = 1, xsize(2)
         jdx = xstart(2) + j -1
         u_fluct_in(j,k) =  ux(plane_index,j,k) - recy_mean_z(1,jdx)
         v_fluct_in(j,k) =  uy(plane_index,j,k) - recy_mean_z(2,jdx)
         w_fluct_in(j,k) =  uz(plane_index,j,k) - recy_mean_z(3,jdx)
      enddo
   enddo

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

      real(mytype) :: u_tau , delta_v
      integer(4) :: i,j

      allocate(Cf(nx), delta_plus(nx))
      allocate(y_plus(ny), u_plus(ny))

      call delta_plus_correlation(re, delta_plus)
      call Cf_correlation(re, Cf)

      do i = 1, nx

          call initWallUnitCalc(Cf(i), u_tau, delta_v)

          y_plus(:) = yp(:)/delta_v 
      
          call uplus_nickels(delta_plus(i), u_tau, y_plus, u_plus)
          
          do j =1, ny
              u_mean(i,j) = u_plus(j)*u_tau
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
   real(mytype) :: T_period, dtdivT

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
         T_period = 100.0

      elseif (itime.lt. t_avg2) then
         T_period = 1000.0

      else
         T_period = 1000.0 + itime - t_avg2

      endif

      dtdivT = dt / T_period
      inlt_mean_t(:,:) = dtdivT *  inlt_mean_z(:,:) + &
                           ( one - dtdivT )*inlt_mean_t(:,:)

      recy_mean_t(:,:) = (dtdivT) *  recy_mean_z(:,:) + &
                           ( one - dtdivT )*recy_mean_t(:,:)                           
   endif

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
      if ((mod(itime,ilist)==0.or.itime==1) .and. itr.eq.1) then
         if (itime==1) then
            open(newunit=unit,file='tbl_recy.log',status='replace',action='readwrite')
            write(unit,"(A,*(',',A))") "itime","theta_inlt", "theta_recy",&
                                    "delta_r","delta_meas","delta_old","delta_i",&
                                    "u_tau_inlt","delta_v_recy","delta_v_inlt"
         else
            open(newunit=unit,file='tbl_recy.log',status='old',action='readwrite')
            rewind(unit)
            nreads = itime/ilist + 1
            if (ilist==1) nreads = itime/ilist - 1
            do i = 1, nreads
               read(unit,*)
            enddo
         endif
         
         write(unit,"(I0,*(',',g0))") itime,theta_inlt,theta_recy,delta_r,&
                                  delta_meas,delta_inlt_old, delta_i, u_tau, delta_v_recy,&
                                  delta_v_inlt

         close(unit)
          
      endif
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
