!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module channel

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character
  real(mytype), dimension(:), allocatable :: body_force, base_force
  PRIVATE ! All functions/subroutines private by default
  procedure(real(mytype)), pointer :: temp_accel_calc
  real(mytype), dimension(:), allocatable :: ub_temp_accel
  real(mytype) :: source_val, ub_old

  PUBLIC :: init_channel, boundary_conditions_channel, postprocess_channel, &
            visu_channel, visu_channel_init, momentum_forcing_channel, &
            geomcomplex_channel, body_force, temp_accel_init, body_forces_init,&
            temp_accel_calc, write_params_channel

contains
  !############################################################################
  subroutine init_channel (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use dbg_schemes, only: exp_prec, abs_prec, sqrt_prec
#ifdef DEBG 
    use tools, only : avg3d
#endif
    

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    real(mytype) :: ftent, u_max
    integer :: k,j,i,fh,ierror,ii,is,it,code, jj

    !integer, dimension (:), allocatable :: seed
    !real(mytype), dimension(:,:,:), allocatable :: urand

    integer :: xshift, yshift, zshift

    integer ( kind = 4 ) :: seed1, seed2, seed3, seed11, seed22, seed33
    integer ( kind = 4 ) :: return_30k
    integer ( kind = 4 ), parameter :: nsemini = 1000 ! For the moment we fix it but after this can go in the input file
    real(mytype), dimension(3,nsemini) :: eddy, posvor
    real(mytype)     :: volsemini, rrand, ddx, ddy, ddz, lsem, upr, vpr, wpr, rrand1
    real(mytype), dimension(3) :: dim_min, dim_max
    real( kind = 8 ) :: r8_random
    external r8_random, return_30k
#ifdef DEBG 
    real(mytype) avg_param
#endif


    if (idir_stream /= 1 .and. idir_stream /= 3) then
       if (nrank == 0) then
          write(*,*) '!! ERROR in imposing sorce term for momentum !!'
          write(*,*) '!! idir_stream ', idir_stream
          write(*,*) '!! idir_stream has to be:'
          write(*,*) '!! - 1 for streamwise direction in X'
          write(*,*) '!! - 3 for streamwise direction in Z'
          write(*,*) '!! Y is not supported and other values do not make sense'
          write(*,*) '!! Calculation will be now stop'
        endif
        call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
    endif

    if (iscalar==1) then
       if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) then
          write(*,*) 'Imposing linear temperature profile'
       end if
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
             if (istret/=0) y=yp(j+xstart(2)-1)
             do i=1,xsize(1)
                phi1(i,j,k,:) = one - y/yly
             enddo
          enddo
       enddo

       phi1(:,:,:,:) = zero !change as much as you want
       if ((nclyS1 == 2).and.(xstart(2) == 1)) then
         !! Generate a hot patch on bottom boundary
         phi1(:,1,:,:) = one
       endif
       if ((nclySn == 2).and.(xend(2) == ny)) then
         phi1(:,xsize(2),:,:) = zero
       endif
    endif
!
    ux1=zero
    uy1=zero
    uz1=zero
    byx1=zero;byy1=zero;byz1=zero
    ! if to decide type of initialization to apply 
    if (iin == 0) then ! laminar flow
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
             if (istret/=0) y=yp(j+xstart(2)-1)-yly*half
             um=exp_prec(-zptwo*y*y)
             do i=1,xsize(1)
                if (idir_stream == 1) then
                   ux1(i,j,k)=one-y*y
                   uy1(i,j,k)=zero
                   uz1(i,j,k)=sin(real(i-1,mytype)*dx)+cos(real(k-1,mytype)*dz)
                else
                   uz1(i,j,k)=one-y*y
                   uy1(i,j,k)=zero
                   ux1(i,j,k)=zero
                endif
             enddo
          enddo
       enddo     
    elseif (iin <= 2) then ! Traditional init to turbulent flows using random numbers + lam profile
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
       !modulation of the random noise + initial velocity profile
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
             if (istret/=0) y=yp(j+xstart(2)-1)-yly*half
             um=exp_prec(-zptwo*y*y)
             do i=1,xsize(1)
                if (use_center) u_max = one
                if (.not.use_center) u_max = onepfive
                if (idir_stream == 1) then
                   ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+u_max*(one-y*y)
                   uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
                   uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
                else
                   uz1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+u_max*(one-y*y)
                   uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
                   ux1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
                endif
             enddo
          enddo
       enddo
    ! iin = 3 is for inlet-outlet files
    elseif (iin == 4) then ! Simplified version of SEM 
       dim_min(1) = zero
       dim_min(2) = zero
       dim_min(3) = zero
       dim_max(1) = xlx
       dim_max(2) = yly
       dim_max(3) = zlz
       volsemini = xlx * yly * zlz
       ! 3 int to get different random numbers
       seed1 =  2345
       seed2 = 13456
       seed3 = 24567
       do jj=1,nsemini
          ! Vortex Position
          do ii=1,3
             seed11 = return_30k(seed1+jj*2+ii*379)
             seed22 = return_30k(seed2+jj*5+ii*5250)
             seed33 = return_30k(seed3+jj*3+ii*8170)
             rrand1  = real(r8_random(seed11, seed22, seed33),mytype)
             call random_number(rrand)
             !write(*,*) ' rr r1 ', rrand, rrand1
             posvor(ii,jj) = dim_min(ii)+(dim_max(ii)-dim_min(ii))*rrand
          enddo
          ! Eddy intensity
          do ii=1,3
             seed11 = return_30k(seed1+jj*7+ii*7924)
             seed22 = return_30k(seed2+jj*11+ii*999)
             seed33 = return_30k(seed3+jj*5+ii*5054)
             rrand1  = real(r8_random(seed11, seed22, seed33),mytype)
             call random_number(rrand)
             !write(*,*) ' rr r1 ', rrand, rrand1
             if (rrand <= zpfive) then
                eddy(ii,jj) = -one
             else
                eddy(ii,jj) = +one
             endif 
          enddo
       enddo
       !do jj=1,nsemini
       !   write(*,*) 'posvor ', posvor(1,jj), posvor(2,jj), posvor(3,jj)
       !   write(*,*) 'eddy   ', eddy(1,jj)  , eddy(2,jj)  , eddy(3,jj)
       !   write(*,*) '  '
       !enddo
       ! Loops to apply the fluctuations 
       do k=1,xsize(3)
          z=real((k+xstart(3)-1-1),mytype)*dz
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
             if (istret/=0) y=yp(j+xstart(2)-1)
             do i=1,xsize(1)
                x=real(i-1,mytype)*dx
                lsem = 0.15_mytype ! For the moment we keep it constant
                upr = zero
                vpr = zero
                wpr = zero
                do jj=1,nsemini
                   ddx = abs_prec(x-posvor(1,jj))
                   ddy = abs_prec(y-posvor(2,jj))
                   ddz = abs_prec(z-posvor(3,jj))
                   if (ddx < lsem .and. ddy < lsem .and. ddz < lsem) then
                      ! coefficients for the intensity of the fluctuation
                      ftent = (one-ddx/lsem)*(one-ddy/lsem)*(one-ddz/lsem)
                      ftent = ftent / (sqrt_prec(two/three*lsem))**3
                      upr = upr + eddy(1,jj) * ftent
                      vpr = vpr + eddy(2,jj) * ftent
                      wpr = wpr + eddy(3,jj) * ftent
                   endif
                enddo
                upr = upr * sqrt_prec(volsemini/nsemini)
                vpr = vpr * sqrt_prec(volsemini/nsemini)
                wpr = wpr * sqrt_prec(volsemini/nsemini)
                ! 
                um  = one-(y-yly*half)**2 ! we can use a better arroximation 
                if (idir_stream == 1) then
                   ux1(i,j,k)=upr*sqrt_prec(two/three*init_noise*um) + um
                   uy1(i,j,k)=vpr*sqrt_prec(two/three*init_noise*um)
                   uz1(i,j,k)=wpr*sqrt_prec(two/three*init_noise*um)
                else
                   uz1(i,j,k)=upr*sqrt_prec(two/three*init_noise*um) + um
                   uy1(i,j,k)=vpr*sqrt_prec(two/three*init_noise*um)
                   ux1(i,j,k)=wpr*sqrt_prec(two/three*init_noise*um)
                endif
             enddo
          enddo
       enddo
    endif
   
    !INIT FOR G AND U=MEAN FLOW + NOISE 
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    avg_param = zero
    call avg3d (ux1, avg_param)
    if (nrank == 0) write(*,*)'## SUB Channel Init ux_avg ', avg_param
    avg_param = zero
    call avg3d (uy1, avg_param)
    if (nrank == 0) write(*,*)'## SUB Channel Init uy_avg ', avg_param
    avg_param = zero
    call avg3d (uz1, avg_param)
    if (nrank == 0) write(*,*)'## SUB Channel Init uz_avg ', avg_param
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_channel
  !############################################################################
  subroutine write_params_channel
   use param
   use variables
   real(mytype), dimension(:), allocatable :: u_b, t_b
   character(80) :: xfmt, yfmt
   integer :: fl, i
   if (nrank .ne. 0) return
   open(newunit=fl,file='parameters.json',status='old',action='write',position='append')

   if (ibodyforces.eq.1) then
      write(yfmt,'(A,I0,A)') "( A, ': [',g0,",ny-1,"(',',g0),'],')"
      
      write(fl,"(A ,': {')") '  "bodyforces"'
      write(fl,yfmt) '    "bf_array"', body_force
      write(fl,"(A,':'I0,',')") '    "ibftype"', ibftype
      write(fl,"(A,':'I0)")'    "itempbf"', itempbf

      if (itempaccel /= 0) write(fl,'(A)') "  },"
      if (itempaccel == 0) write(fl,'(A)') "  }"
   endif

   if(itempaccel == 1) then
      write(fl,"(A ,': {')") '  "temp_accel"'
      if (iacceltype== 1) then
         write(fl,"(A,': ',A,',')") '    "profile"','"linear"'
         write(fl,"(A,': ',g0,',')") '    "t_start"',t_start
         write(fl,"(A,': ',g0,',')") '    "t_end"',t_start
         write(fl,"(A,': ',g0,',')") '    "Re_ratio"',Re_ratio
      else if (iacceltype==2) then
         write(fl,"(A,': ',A,',')") '    "profile"','"spatial equiv"'
         write(fl,"(A,': ',g0,',')") '    "U_ratio"',U_ratio
         write(fl,"(A,': ',g0,',')") '    "x0"',accel_centre
         write(fl,"(A,': ',g0,',')") '    "alpha_accel"',alpha_accel
      endif

      allocate(u_b(ilast/ilist))
      allocate(t_b(ilast/ilist))
      
      do i = 1,ilast/ilist
         t_b(i) = real(i,kind=mytype)*dt*ilist
         u_b(i) = temp_accel_calc( t_b(i))
      enddo
      write(yfmt,'(A,I0,A)') "( A, ': [',g0,",ilast/ilist-1,"(',',g0),'],')"
      write(fl,yfmt) '    "t"', t_b
      write(yfmt,'(A,I0,A)') "( A, ': [',g0,",ilast/ilist-1,"(',',g0),']')"
      write(fl,yfmt) '    "U_b"', u_b

      deallocate(u_b,t_b)
      write(fl,'(A)') "  }"
   endif
   write(fl,'(A)') "}"
   close(fl)
  end subroutine write_params_channel

  !############################################################################
  subroutine boundary_conditions_channel (ux,uy,uz,phi)

    use param
    use var, only : di2
    use variables
    use decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: u_b

    if (use_center) then
      u_b = two/three
    else
      u_b = one
    endif
    if (.not. cpg) then ! if not constant pressure gradient
       
      if (.not. use_center) then
         if (itempaccel==1) u_b = temp_accel_calc(t) 
       
         if (idir_stream == 1) then
            call channel_cfr(ux,u_b)
         else
            call channel_cfr(uz,u_b)
         endif
      else
         if (itempaccel==1) u_b = temp_accel_calc(t) 
       
         if (idir_stream == 1) then
            call channel_c_center(ux,u_b)
         else
            call channel_c_center(uz,u_b)
         endif
      endif
      if (itr==1) source_val = compute_amp(u_b)
    end if

    if (iscalar /= 0) then
       if (iimplicit <= 0) then
          if ((nclyS1 == 2).and.(xstart(2) == 1)) then
             !! Generate a hot patch on bottom boundary
             phi(:,1,:,:) = one
          endif
          if ((nclySn == 2).and.(xend(2) == ny)) THEN
             phi(:,xsize(2),:,:) = zero
          endif
       else
          !
          ! Implicit boundary conditions are usually given in input file
          ! It is possible to modify g_sc here
          ! It is not possible to modify alpha_sc and beta_sc here
          !
          ! Bottom temperature if alpha_sc(:,1)=1 and beta_sc(:,1)=0 (default)
          !if (nclyS1.eq.2) g_sc(:,1) = one
          ! Top temperature if alpha_sc(:,2)=1 and beta_sc(:,2)=0 (default)
          !if (nclySn.eq.2) g_sc(:,2) = zero
       endif
    endif

  end subroutine boundary_conditions_channel
  !############################################################################
  !!
  !!  SUBROUTINE: channel_cfr
  !!      AUTHOR: Kay Schäfer
  !! DESCRIPTION: Inforces constant flow rate without need of data transposition
  !!
  !############################################################################
  subroutine channel_cfr (ux, constant)

    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux
    real(mytype), intent(in) :: constant

    integer :: code, i, j, k, jloc
    real(mytype) :: can, ub, uball, coeff

    ub = zero
    uball = zero
    coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

    do k = 1, xsize(3)
       do jloc = 1, xsize(2)
          j = jloc + xstart(2) - 1
          do i = 1, xsize(1)
            ub = ub + ux(i,jloc,k) / ppy(j)
          enddo
       enddo
    enddo

    ub = ub * coeff

    call MPI_ALLREDUCE(ub,uball,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can = - (constant - uball)

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) 'UT', uball, can

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux(i,j,k) = ux(i,j,k) - can
        enddo
      enddo
    enddo

  end subroutine channel_cfr

  subroutine channel_c_center(ux, constant)
   use MPI
   real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux
   real(mytype), intent(in) :: constant
   integer, dimension(:),allocatable :: ranks_incl,displs, recvcounts

   integer :: y_mid, y_ind
   integer :: ierr, i, j, k
   real(mytype) :: u_c, u_c_global, can

   y_mid = (ny+1) / 2
   if (y_mid>=xstart(2).and.y_mid<=xend(2)) then
      y_ind = y_mid - xstart(2) +1
      u_c = sum(ux(:,y_ind,:))/nx/nz
   else
      u_c = zero
   endif
   
   call MPI_Allreduce(u_c,u_c_global,1,real_type,MPI_SUM,DECOMP_2D_COMM_CART_X,ierr)
   can = -(constant - u_c_global)

   if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
   write(*,*) 'UT', u_c_global, can

   do k=1,xsize(3)
      do j=1,xsize(2)
         do i=1,xsize(1)
            ux(i,j,k) = ux(i,j,k) - can
         enddo
      enddo
   enddo
  end subroutine
  !! calculation of temporal acceleration

   function compute_amp(u_b) result(dudt)
      real(mytype),intent(in) :: u_b
      real(mytype) :: dudt

      integer :: unit

      if (itempbf ==0) then
         dudt = one
         return
      endif
      dudt = (u_b - ub_old)/gdt(itr)
      ub_old = u_b
      if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) then
         write(*,*) 'itempbf: dudt', dudt
         if (itime==1) then
            open(newunit=unit,file='dudt.csv',status='replace',action='write',position='append')
            write(unit,*) '# itime, t, dudt'
         else 
            open(newunit=unit,file='dudt.csv',status='old',action='write',position='append')
         endif
         write(unit,'(I0,",", g0,",", g0)') itime, t, dudt
         close(unit)
      endif
   end function
  real(mytype) function linear_temp(temp)
   real(mytype), intent(in) :: temp

   if (temp > t_end) then
      linear_temp = Re_ratio
   else if (temp > t_start) then
      linear_temp = one + (temp - t_start)&
                        * (Re_ratio - one)&
                        / (t_end - t_start)
   else
      linear_temp = one
   endif

  end function linear_temp

  real(mytype) function spatial_equiv(temp)
     real(mytype), intent(in) :: temp

      integer :: i 

      i = int(temp/dt)

      spatial_equiv = ub_temp_accel(i)
  end function spatial_equiv
  !############################################################################
  !############################################################################
  subroutine postprocess_channel(ux1,uy1,uz1,pp3,phi1,ep1)

    use var, ONLY : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    
    logical :: exists
    integer :: unit
    character(len=128) :: fname

    if (mod(itime,istatout)==0 .and. itime>=initstat) then
      if (nrank==0) then
         inquire(file="body_force",exist=exists)
         if (.not. exists) call system("mkdir -p body_force")
         write(fname,'(A,"/",A,"-",I7.7)') "body_force", "bodyf", itime
         open(newunit=unit,file=fname,status='replace',action='write',access='stream')
         write(unit) body_force
         close(unit)
      endif

    endif
  end subroutine postprocess_channel
  subroutine visu_channel_init(visu_initialised)

   use decomp_2d, only : mytype
   use decomp_2d_io, only : decomp_2d_register_variable
   use visu, only : io_name, output2D
   use param, only : istatlambda2
   
   implicit none

   logical, intent(out) :: visu_initialised

   if (istatlambda2) then
      call decomp_2d_register_variable(io_name, "lambda2", 1, 0, output2D, mytype)
   endif
   visu_initialised = .true.
    
  end subroutine visu_channel_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_channel
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs channel-specific visualization
  !!
  !############################################################################
  subroutine visu_channel(ux1, uy1, uz1, pp3, phi1, ep1, num)

   use visu, only : write_field
   use stats, only : lambda2
   use param, only : istatlambda2, initstat
   use var, ONLY : nxmsize, nymsize, nzmsize, ta1, ta2, itime
   use decomp_2d
   implicit none

   real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
   real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
   real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
   real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
   character(len=32), intent(in) :: num


  if (istatlambda2.and.(itime>initstat2.or.(itempaccel.eq.1.and.itime>initstat))) then
     call transpose_z_to_y(lambda2,ta2)
     call transpose_y_to_x(ta2,ta1)
     call write_field(ta1, ".", "lambda2", trim(num), flush=.true.) ! Reusing temporary array, force flush
  endif
  end subroutine visu_channel
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies rotation for t < spinup_time.
  !!
  !############################################################################
  subroutine momentum_forcing_channel(dux1, duy1, duz1, ux1, uy1, uz1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
    integer :: j, jloc

    if (cpg) then
        !! fcpg: add constant pressure gradient in streamwise direction
        if (idir_stream == 1) then
           dux1(:,:,:,1) = dux1(:,:,:,1) + fcpg !* (re/re_cent)**2
        else
           duz1(:,:,:,1) = duz1(:,:,:,1) + fcpg !* (re/re_cent)**2
        endif
    endif

    ! To update to take into account possible flow in z dir
    if (itime < spinup_time .and. iin <= 2) then
       if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
          write(*,*) 'Rotating turbulent channel at speed ',wrotation
       dux1(:,:,:,1) = dux1(:,:,:,1) - wrotation*uy1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + wrotation*ux1(:,:,:)
    endif

    if (itempbf==2) call compute_tempbf(ux1)

    if (ibodyforces.eq.1) then
      if (idir_stream == 1) then
         do j = 1,xsize(2)
            jloc = j + xstart(2) -1
            dux1(:,j,:,1) = dux1(:,j,:,1) + source_val*body_force(jloc)
         enddo

      else
         do j = 1,xsize(2)
            jloc = j + xstart(2) -1
            duz1(:,j,:,1) = duz1(:,j,:,1) + source_val*body_force(jloc)
         enddo
      endif


    endif
  end subroutine momentum_forcing_channel
  subroutine compute_tempbf(ux)
   use dbg_schemes, only: abs_prec
   use MPI
   real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux
   real(mytype), dimension(xsize(2)) :: u_l, u
   real(mytype), dimension(ny) :: u_g, dudy
   integer, dimension(:), allocatable :: recvcounts, displs
   integer :: code, split_comm_y, split_comm_z, key, color, i, j, k
   
   integer, dimension(2) :: dims, coords
   logical, dimension(2) :: periods
   real(mytype) :: a, b, c, ddy1, ddy2

   if (itempbf==2) then
      u_l(:) = zero
      do k = 1, xsize(3)
         do j = 1, xsize(2)
            do i = 1, xsize(1)
               u_l(j) = u_l(j) + ux(i,j,k)
            enddo
         enddo
      enddo

      u_l(:) = u_l(:)/real(nx*nz,kind=mytype)

      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, code)

      key = coords(1)
      color = coords(2)
  
      call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key,color, split_comm_y,code)
      call MPI_Comm_split(DECOMP_2D_COMM_CART_X, color,key, split_comm_z,code)
  
      call MPI_Allreduce(u_l,u,xsize(2),&
                         real_type,MPI_SUM,split_comm_y,code)

      allocate(recvcounts(dims(1)),displs(dims(1)))
      call MPI_Allgather(xsize(2),1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,split_comm_z,code)
      call MPI_Allgather(xstart(2)-1,1,MPI_INTEGER,displs,1,MPI_INTEGER,split_comm_z,code)
      call MPI_AllGatherv(u,xsize(2),real_type,u_g, recvcounts,displs,&
                           real_type,split_comm_z,code)

      do j = 2, ny-1
         ddy1 = yp(j) - yp(j-1)
         ddy2 = yp(j+1) - yp(j)

         a = -(ddy2)/(ddy1 * (ddy1 + ddy2))
         b = (ddy2 - ddy1) / (ddy1 * ddy2)
         c = ddy1 / (ddy2 * (ddy1 + ddy2))

         dudy(j) = a*u_g(j-1) + b*u_g(j) + c*u_g(j+1)
      enddo
      
      ddy1 = yp(2) - yp(1)
      ddy2 = yp(3) - yp(2)

      a = -( two*ddy1 + ddy2)/(ddy1*(ddy1 + ddy2))
      b = (ddy1 + ddy2)/(ddy1*ddy2)
      c = -ddy1/(ddy2*(ddy1 + ddy2))

      dudy(1) = a*u_g(1) + b*u_g(2) + c*u_g(3)
      
      ddy1 = yp(ny-1) - yp(ny-2)
      ddy2 = yp(ny) - yp(ny-1)

      a = ddy2/(ddy1*(ddy1 + ddy2))
      b = -(ddy1 + ddy2)/(ddy1*ddy2)
      c = (two*ddy2 + ddy1)/(ddy2*(ddy1 + ddy2))
  
      dudy(ny) = a*u_g(ny-2) + b*u_g(ny-1) + c*u_g(ny)

      do j = 1, ny
         body_force(j) = zpfive*base_force(j)*(dudy(j) - dudy(ny-j+1))*sign(one,-(yp(j)-half*yly))
      enddo
      call MPI_Comm_free(split_comm_y,code)
      call MPI_Comm_free(split_comm_z,code)
   endif
  end subroutine
  subroutine mass_flow_conserve(du1)
   use MPI
   real(mytype), intent(inout), dimension(xsize(1), xsize(2), xsize(3),ntime) :: du1
   real(mytype) :: dvol, int_du, coeff, int_du_all
   integer :: jloc, i, j, k, code

   coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

   int_du = zero
   do k = 1, xsize(3)
      do jloc = 1, xsize(2)
         j = jloc + xstart(2) - 1
         do i = 1, xsize(1)
            int_du = int_du + du1(i,jloc,k,1) / ppy(j)
         enddo
      enddo
   enddo

   int_du = int_du*coeff
   call MPI_ALLREDUCE(int_du,int_du_all,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

   du1(:,:,:,1) = du1(:,:,:,1) - int_du_all

  end subroutine
  subroutine temp_accel_init

   use param
   use dbg_schemes, only : tanh_prec

   integer :: i
   real(mytype) :: dudx, x=0

   if (itempaccel == 1) then
      if (iacceltype == 1) then
         temp_accel_calc => linear_temp
      else if (iacceltype == 2) then
         allocate(ub_temp_accel(ilast))
         ub_temp_accel(1) = one
         do i = 2, ilast
            x = x + ub_temp_accel(i-1)*dt
            ub_temp_accel(i) = one +   half *(U_ratio - one)*(&
                                       tanh_prec( alpha_accel*(x &
                                       - accel_centre ))  + one )
         enddo

         temp_accel_calc => spatial_equiv
      else
         write(*,*) "Invalid temporal acceleration profile"
      endif
   endif
  end subroutine
  subroutine body_forces_init
   use dbg_schemes, only : abs_prec, sin_prec, tanh_prec
   integer :: j
   real(mytype) :: lim, y
   if (ibodyforces.eq.1) then

      allocate(body_force(ny))
      body_force = zero
      ub_old = one

      if (ibftype.eq.1) then
         lim = one - bf_ext
         do j = 1, ny
            if (istret==0) y=real(j,mytype)*dy-yly*half
            if (istret/=0) y=yp(j)-yly*half
            
            if (abs_prec(y)>lim) then
               body_force(j) = bf_amp*(abs_prec(y)-lim)/bf_ext
            endif
         enddo
      elseif (ibftype.eq.2) then
         lim = one - bf_ext
         do j = 1, ny
            if (istret==0) y=real(j,mytype)*dy-yly*half
            if (istret/=0) y=yp(j)-yly*half
            
            if (abs_prec(y)>lim) then
               body_force(j) = bf_amp*sin_prec(pi*(one - abs_prec(y))/bf_ext)**2
            endif
         enddo
      elseif (ibftype.eq.3) then
         allocate(base_force(ny))
         do j = 1, ny
            if (istret==0) y=real(j,mytype)*dy
            if (istret/=0) y=yp(j)
            
            base_force(j) = zpfive*bf_amp*(tanh_prec(bf_alp*(y-bf_ext)) &
                                             + tanh_prec(-bf_alp*(y-(yly-bf_ext))))
         enddo
         body_force(:) = base_force(:)
      endif
   endif

   
  end subroutine
  !############################################################################
  !############################################################################
  subroutine geomcomplex_channel(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,yp,remp)

    use decomp_2d, only : mytype
    use param, only : zero, one, two, ten
    use ibm

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: remp
    integer                    :: j
    real(mytype)               :: ym
    real(mytype)               :: zeromach
    real(mytype)               :: h

    epsi(:,:,:) = zero
    h = (yly - two) / two

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = ten*zeromach

    do j=nyi,nyf
       ym=yp(j)
       if ((ym.le.h).or.(ym.ge.(h+two))) then
          epsi(:,j,:)=remp
       endif
    enddo

    return
  end subroutine geomcomplex_channel
  !############################################################################
end module channel
