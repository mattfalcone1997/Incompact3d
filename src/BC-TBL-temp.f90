module tbl_temp
    use decomp_2d
    implicit none
    private

    procedure(real(mytype)), pointer :: wall_velocity
    real(mytype), allocatable, dimension(:) :: uw_temp_accel
    integer :: itime_ref
    logical :: equiv_initialised = .false.

    public :: setup_tbl_temp, init_tbl_temp, momentum_forcing_tbl_temp
    public :: boundary_conditions_tbl_temp, postprocess_tbl_temp
    public :: visu_tbl_temp, visu_tbl_temp_init, write_params_tbl_temp
    public :: wall_velocity

    contains
    subroutine write_params_tbl_temp
        use param
        use variables
        real(mytype), dimension(:), allocatable :: u_w, t_b
        character(80) :: yfmt
        integer :: fl, i

        allocate(u_w(ilast/ilist))
        allocate(t_b(ilast/ilist))
        
        if (iaccel .eq.2) then
            do i = 1,ilast/ilist
                t_b(i) = real(i,kind=mytype)*dt*ilist
                u_w(i) = one
            enddo
        else
            do i = 1,ilast/ilist
                t_b(i) = real(i,kind=mytype)*dt*ilist
                u_w(i) = wall_velocity( t_b(i))
            enddo
        endif

        if (nrank.ne.0) then
            deallocate(u_w,t_b)
            return
        endif

        open(newunit=fl,file='parameters.json',status='old',action='write',position='append')
        write(fl,"(A ,': {')") '  "tbl_temp_accel"'
        if(iaccel==1) then
            write(fl,"(A,': ',A,',')") '    "profile"','"linear"'
            write(fl,"(A,': ',g0,',')") '    "U_ratio"',U_ratio
            write(fl,"(A,': ',g0,',')") '    "t_start"',t_start
            write(fl,"(A,': ',g0,',')") '    "t_end"',t_end
        else if (iaccel==2) then
            write(fl,"(A,': ',A,',')") '    "profile"','"spatial equiv"'
            write(fl,"(A,': ',g0,',')") '    "U_ratio"',U_ratio
            write(fl,"(A,': ',g0,',')") '    "x0"',accel_centre
            write(fl,"(A,': ',g0,',')") '    "alpha_accel"',alpha_accel
            write(fl,"(A,': ',g0,',')") '    "re_ref"',re_ref
        endif
        
        write(yfmt,'(A,I0,A)') "( A, ': [',g0,",ilast/ilist-1,"(',',g0),'],')"
        write(fl,yfmt) '    "t"', t_b
        write(yfmt,'(A,I0,A)') "( A, ': [',g0,",ilast/ilist-1,"(',',g0),']')"
        write(fl,yfmt) '    "U_w"', u_w

        write(fl,'(A)') "  }"
        write(fl,'(A)') "}"

        close(fl)
        deallocate(u_w,t_b)
    end subroutine
    subroutine setup_tbl_temp()
        use param, only : itype, itype_tbl_temp, iaccel, t_start, t_end
        use param, only : one, half, U_ratio, accel_centre, re_ref, ilast
        use param, only: alpha_accel
        use var, only : dt
        use dbg_schemes, only : tanh_prec
        use MPI
        integer :: code, i
        real(mytype) :: x

        if (itype.ne.itype_tbl_temp) return


        if (iaccel == 0) then
            wall_velocity => const_velo
        else if (iaccel == 1) then
            if (t_start > t_end.and.nrank==0) then
                write(*,*) "t_start cannot be greater than t_end" 
                call MPI_Abort(MPI_COMM_WORLD,1,code)
            endif
            wall_velocity => linear_accel
        else if (iaccel == 2) then
            if (re_ref<0.and.nrank==0) then
                write(*,*) "re_ref must be greater than 0"
                call MPI_Abort(MPI_COMM_WORLD,1,code)
            endif
            allocate(uw_temp_accel(ilast))
            uw_temp_accel(1)= one
            do i = 2, ilast
                x = x + uw_temp_accel(i-1)*dt
                uw_temp_accel(i) = one +   half *(U_ratio - one)*(&
                                           tanh_prec( alpha_accel*(x &
                                           - accel_centre ))  + one )
             enddo
             equiv_initialised = .false.
             wall_velocity => spatial_equiv
        endif


    end subroutine
    subroutine init_tbl_temp(ux1,uy1,uz1,ep1,phi1)

        use decomp_2d_io
        use param , only : one, zptwofive, zpfive, Re_D, re
        use param, only : init_noise, istret, two, zero
        use MPI
        use dbg_schemes
        use var, only: numscalar, dy, byx1, byy1, byz1, byxn, byyn, byzn, ny, yp
    
    
        implicit none
    
        real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
        real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

        real(mytype), dimension(ny) :: u_mean
        real(mytype) :: D, theta_sl, arg, y, y_thresh, um
        integer :: i, j, k, jloc, unit, code, ii

        D = Re_D/re
        theta_sl = 54.0_mytype/re
        do j = 1, ny
            if (istret==0) y = real(j-1,mytype)*dy
            if (istret/=0) y = yp(j)

            arg = zpfive*D*(one- y/D)/theta_sl
            u_mean(j) = zpfive + zpfive*tanh_prec(arg)
        enddo

    y_thresh = D - two*theta_sl*atanh(0.98_mytype)
    ! variation of the channel flow fluctuation inlet

    ux1=zero
    uy1=zero
    uz1=zero

    byx1=one;byy1=zero;byz1=zero
    byxn=zero;byyn=zero;byzn=zero

    ! if to decide type of initialization to apply 
    
    call system_clock(count=code)
    call random_seed(size = ii)
    call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)

    !modulation of the random noise + initial velocity profile
    do k=1,xsize(3)
        do j=1,xsize(2)
            um = one !y/y_thresh
            do i = 1,xsize(1)
                jloc = xstart(2) + j -1
                if (istret==0) y = real(jloc-1,mytype)*dy
                if (istret/=0) y = yp(jloc)
                if (y<=y_thresh) then
                    ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+u_mean(jloc)
                    uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
                    uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
                else
                    ux1(i,j,k)=u_mean(jloc)
                    uy1(i,j,k)=zero
                    uz1(i,j,k)=zero
                endif
            enddo
        enddo
    enddo

    if (nrank==0) then
        open(newunit=unit,file="tbl_temp.csv",action='write',status='new')
        write(unit,'("# ",A,",",A,",",A)') "J", "y", "u_mean"
        do j = 1, ny
            if (istret==0) y = real(j-1,mytype)*dy
            if (istret/=0) y = yp(j)
            write(unit,'(I0,",",g0,","g0)') j, y, u_mean(j)
        enddo
        close(unit)
    endif

    end subroutine init_tbl_temp
    subroutine momentum_forcing_tbl_temp(dux1, duy1, duz1, ux1, uy1, uz1)
        use var, only: ntime, t
        use param, only: ten, Re_D, re, x0_tr_tbl, y0_tr_tbl, ys_tr_tbl,A_tr, zero, one
        use var, only : ta1
        real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
        real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
        
        ! if (t>ten.or.t<one) return

        ! A_tr = 0.1_mytype
        ! x0_tr_tbl = 100.0_mytype
        ! y0_tr_tbl= Re_D/re
        ! ys_tr_tbl=Re_D/re

        ! call tbl_tripping(dux1,ta1)

    end subroutine momentum_forcing_tbl_temp

    subroutine boundary_conditions_tbl_temp (ux,uy,uz,phi)
        use param, only : ncly1, nclyn, one, zero, re
        use param, only: ilast, ifirst, iaccel, re_ref
        use var, only: numscalar, byx1, byy1, byz1, byxn, byyn
        use var, only:  byzn, yp, t, itime, itr
        use variables, only: ilist

        implicit none
    
        real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
        real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
        real(mytype) :: wall_velo, theta
        integer :: unit

        wall_velo = wall_velocity(t)
        theta = compute_theta(ux)
        if ((mod(itime,ilist)==0.or. itime == ifirst .or. itime == ilast) .and. nrank ==0.and.itr==1) then
             if (iaccel .ne. 2)   write(*,*) 'U_w', wall_velo
             if (iaccel .eq.2) then
                write(*,*) 'U_w', wall_velo, 're_theta', re*theta, 're_ref', re_ref
                if (itime==ifirst) then
                    open(newunit=unit,file='theta.csv',status='replace',action='write')
                    write(unit,'("#",A,",",A,",",A)') "itime", "re_theta", "U_w"
                else
                    open(newunit=unit,file='theta.csv',status='old',action='write',position='append')
                endif

                write(unit,'(I0,",",g0,",",g0)') itime, re*theta, wall_velo
                close(unit)
             endif
        endif

        if (ncly1==2) then
            byx1=wall_velo
            byy1=zero;byz1=zero
        endif
        if (nclyn==2) then
            byxn=zero;byyn=zero;byzn=zero
        endif

    end subroutine boundary_conditions_tbl_temp

    function const_velo(time) result(wall_velo)
        use param, only: one
        real(mytype), intent(in) :: time
        real(mytype) :: wall_velo
        wall_velo = one
    end function const_velo

    function linear_accel(time) result(wall_velo)
        use param, only: t_start, t_end, U_ratio, one
        real(mytype), intent(in) :: time
        real(mytype) :: wall_velo
        if (time < t_start) then
            wall_velo = one
        elseif (time > t_end) then
            wall_velo = U_ratio
        else
            wall_velo = one + (U_ratio - one)*(time - t_start)/(t_end - t_start)
        endif
    end function linear_accel

    function spatial_equiv(time) result(wall_velo)
        use param, only: U_ratio, re, one, re_ref
        use var, only : ux1, itime
        real(mytype), intent(in) :: time
        real(mytype) :: wall_velo
        real(mytype) :: theta

        if (equiv_initialised) then
            wall_velo = uw_temp_accel(itime-itime_ref)
        else
            theta = compute_theta(ux1)
            if (theta*re > re_ref) then
                equiv_initialised = .true.
                itime_ref=itime-1
                wall_velo = uw_temp_accel(itime-itime_ref)
            else
                wall_velo = one
            endif
        endif

        

    end function spatial_equiv

    function compute_theta(ux) result(theta)
        use param, only: one, zero
        use var, only : ppy,dy, yly
        use MPI
        use decomp_2d, only : xsize
        real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux
        real(mytype) :: theta, theta_local, coeff
        integer :: i, j, k, jloc, code

        coeff = dy / (real(xsize(1) * zsize(3), kind=mytype))

        theta_local = zero
        theta = zero
        do k = 1, xsize(3)
            do j = 1, xsize(2)
                jloc = j + xstart(2) - 1
                do i = 1, xsize(1)
                    theta_local = theta_local + ux(i,j,k)*(one - ux(i,j,k))/ppy(jloc)
                enddo
            enddo
        enddo

        theta_local = theta_local*coeff
        call MPI_ALLREDUCE(theta_local,theta,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    end function

    subroutine postprocess_tbl_temp(ux1,uy1,uz1,ep1)
        real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1

    end subroutine postprocess_tbl_temp

    subroutine visu_tbl_temp_init (visu_initialised)
        logical, intent(out) :: visu_initialised


    end subroutine visu_tbl_temp_init

    subroutine visu_tbl_temp(ux1, uy1, uz1, pp3, phi1, ep1, num)
        use var, only: numscalar
        use var, only: nzmsize,npress
        real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
        real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
        real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
        real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
        character(len=32), intent(in) :: num

    end subroutine visu_tbl_temp
end module