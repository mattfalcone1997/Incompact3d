module tbl_temp
    use decomp_2d
    implicit none
    private

    procedure(real(mytype)), pointer :: wall_velocity

    public :: setup_tbl_temp, init_tbl_temp, momentum_forcing_tbl_temp
    public :: boundary_conditions_tbl_temp, postprocess_tbl_temp
    public :: visu_tbl_temp, visu_tbl_temp_init

    contains
    subroutine setup_tbl_temp()
        use param, only : itype, itype_tbl_temp, iaccel, t_start, t_end
        use MPI
        integer :: code

        if (itype.ne.itype_tbl_temp) return


        if (iaccel == 0) then
            wall_velocity => const_velo
        else if (iaccel == 1) then
            if (t_start > t_end) then
                write(*,*) "t_start cannot be greater than t_end" 
                call MPI_Abort(MPI_COMM_WORLD,1,code)
            endif
            wall_velocity => linear_accel
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
        use param, only : ncly1, nclyn, one, zero
        use param, only: ilast, ifirst
        use var, only: numscalar, byx1, byy1, byz1, byxn, byyn
        use var, only:  byzn, yp, t, itime
        use variables, only: ilist

        implicit none
    
        real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
        real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

        if ((mod(itime,ilist)==0.or. itime == ifirst .or. itime == ilast) .and. nrank ==0) &
                write(*,*) 'U_w', wall_velocity(t)

        if (ncly1==2) then
            byx1=wall_velocity(t)
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