!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use var
  use case

  use transeq, only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
       calc_divu_constraint, solve_poisson, cor_vel
  use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
  use turbine, only : compute_turbines
  use ibm_param
  use ibm, only : body
  use genepsi, only : genepsi3d

  implicit none


  call init_xcompact3d()

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     call simu_stats(2)

     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype.eq.itype_abl.or.iturbine.ne.0).and.(ifilter.ne.0).and.(ilesmod.ne.0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

        if (imove.eq.1) then ! update epsi for moving objects
          if ((iibm.eq.2).or.(iibm.eq.3)) then
             call genepsi3d(ep1)
          else if (iibm.eq.1) then
             call body(ux1,uy1,uz1,ep1)
          endif
        endif
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)
#ifdef DEBG
        call check_transients()
#endif
        
        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
           call velocity_to_momentum(rho1,ux1,uy1,uz1)
        endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

     enddo !! End sub timesteps

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
     call simu_stats(5)
  enddo !! End time loop

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_init
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use sandbox, only : init_sandbox
  use forces

  use var

  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, compute_cfldiff, &
       init_inflow_outflow

  use param, only : ilesmod, jles,itype, iaccel, itype_tbl_temp
  use param, only : irestart

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : nstat, nvisu, nprobe, ilist

  use les, only: init_explicit_les
  use turbine, only: init_turbines

  use visu, only : visu_init, visu_ready

  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body

  use probes, only : init_probes, init_line_probes
   use tbl_recy, only : restart_tbl_recy

   use channel, only : body_forces_init, temp_accel_init
   use tbl_temp, only : setup_tbl_temp
  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd, i
  logical :: back
  character(len=80) :: InputFN, FNBase
    
  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'Xcompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'Xcompact3d is run with the provided file -->', trim(InputFN)
  endif

#ifdef ADIOS2
  if (nrank .eq. 0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_io_init()
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (iforces.eq.1) then
     call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif

  if (itype == itype_channel) then
   call body_forces_init
   call temp_accel_init
  endif

  if (itype == itype_tbl_temp) then
   call setup_tbl_temp
  endif
  !####################################################################
  ! initialise visu
  if (ivisu.ne.0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     end if

     if (itype == itype_tbl_recy) then
      call restart_tbl_recy(ifirst-1)
     endif
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif
  
!   call postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)
!   call MPI_Barrier(MPI_COMM_WORLD,ierr)
!   if (nrank  == 0) write(*,*) "End."
!    call MPI_ABORT(MPI_COMM_WORLD,1,ierr); stop

  if ((ioutflow.eq.1).or.(iin.eq.3)) then
     call init_inflow_outflow()
  end if

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if ((iibm.eq.1).or.(iibm.eq.3)) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()
  call init_line_probes()

  if (iturbine.ne.0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank.eq.0)then
        open(42,file='time_evol.dat',form='formatted')
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif

  call write_params_json
  call write_run_info

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise

  use tools, only : simu_stats
  use param, only : itype, jles, ilesmod
  use probes, only : finalize_probes
  use visu, only : visu_finalise
  use les, only: finalise_explicit_les

  implicit none

  integer :: ierr
  
  if (itype==2) then
     if(nrank.eq.0)then
        close(42)
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        close(38)
     endif
  endif
  
  call simu_stats(4)
  call finalize_probes()
  call visu_finalise()
  if (ilesmod.ne.0) then
     if (jles.gt.0) call finalise_explicit_les()
  endif
  call decomp_2d_io_finalise()
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d

subroutine check_transients()

  use decomp_2d, only : mytype

  use var
  use tools, only : avg3d
  
  implicit none

  real(mytype) avg_param
  
  avg_param = zero
  call avg3d (dux1, avg_param)
  if (nrank == 0) write(*,*)'## Main dux1 ', avg_param
  avg_param = zero
  call avg3d (duy1, avg_param)
  if (nrank == 0) write(*,*)'## Main duy1 ', avg_param
  avg_param = zero
  call avg3d (duz1, avg_param)
  if (nrank == 0) write(*,*)'## Main duz1 ', avg_param
  
end subroutine check_transients



subroutine write_params_json
   use param
   use variables
   use decomp_2d, only : nrank
   use tbl_recy, only : write_params_tbl_recy
   use channel, only : write_params_channel
   use tbl_temp, only: write_params_tbl_temp
   use stats, only : h_quads
   implicit none

   real(mytype), dimension(nx) :: xcoords, u_infty
   real(mytype), dimension(nz) :: zcoords
   real(mytype), dimension(:), allocatable :: u_b, t_b
   real(mytype) :: u_infty_grad, t_tmp
   character(80) :: xfmt, zfmt,yfmt
   integer :: fl, i
   ! itype
   ! mesh
   ! geometry
   ! re
   ! optionals
   ! acceleration
   ! tbl_recy

   if (nrank .ne. 0) return

   do i =1, nx
      xcoords(i) = real(i-1,kind=mytype)*dx
   enddo
   if (istret.eq.0) then
      do i = 1, ny
         yp(i) = real(i-1,kind=mytype)*dy
      enddo
   endif
   do i =1, nz
      zcoords(i) = real(i-1,kind=mytype)*dz
   enddo

   open(newunit=fl,file='parameters.json',status='replace',action='write')

   write(fl,'(A)') "{"
   write(fl,"(A ,':',I0,',')") '  "itype"',itype
   write(fl,"(A ,':',I0,',')") '  "istatcalc"',istatcalc
   write(fl,"(A ,':',I0,',')") '  "initstat"',initstat
   write(fl,"(A ,':',I0,',')") '  "initstat2"',initstat2

   write(fl,"(A ,' : {')") '  "geometry"'
   write(fl,"(A ,': ',g0,',')") '    "xlx"',xlx
   write(fl,"(A ,' : ',g0,',')") '    "yly"', yly
   write(fl,"(A ,': ',g0)") '    "zlz"', zlz
   write(fl,*) "  },"

   write(xfmt,'(A,I0,A)') "( A, ': [',g0,",nx-1,"(',',g0),'],')"
   write(yfmt,'(A,I0,A)') "( A, ': [',g0,",ny-1,"(',',g0),'],')"
   write(zfmt,'(A,I0,A)') "( A, ': [',g0,",nz-1,"(',',g0),']')"

   write(fl,"(A ,': {')") '  "mesh"'
   write(fl,"(A ,': [',I0,',',I0,',',I0,'],')") '    "sizes"',nx, ny, nz
   write(fl,xfmt) '    "xcoords"',xcoords
   write(fl,yfmt) '    "ycoords"',yp
   write(fl,zfmt) '    "zcoords"',zcoords
   write(fl,'(A)') "  },"   
   if (istatquadrant) then
      write(fl,"(A ,': {')") '  "uv_quadrant"'
      if (nquads > 1) then
         write(xfmt,'(A,I0,A)') "( A, ': [',g0,",nquads-1,"(',',g0),']')"
      else
         xfmt =  "( A, ': [',g0,']')"
      endif
      write(fl,xfmt) '    "h_quads"',h_quads
      write(fl,'(A)') "  },"
   endif
   write(fl,"(A ,':',g0,',')") '  "re"',re
   if (istatautocorr) then
      write(fl,'(A," : ",g0,",")') '    "dt"',dt
   else
      write(fl,'(A," : ",g0)') '    "dt"',dt
   endif   

   close(fl)

   if (itype .eq. itype_tbl_recy) then
      call write_params_tbl_recy
   elseif (itype.eq. itype_channel) then
      call write_params_channel
   elseif(itype.eq.itype_tbl_temp) then
      call write_params_tbl_temp
   endif
end subroutine

subroutine write_run_info
   use param
   use decomp_2d, only : nrank

   implicit none
   integer :: fl, recl
   character(80) :: run_str

   if (nrank .ne. 0) return

   if (irestart==1) then
      open(newunit=fl,file='run_log.json',status='old',action='readwrite',position='rewind')
      recl=0
      do
         read(fl,*,end=10)
         recl = recl + 1
      enddo
      10 write(run_str,'("  ",A," ",I0)') "run",(recl -2)/5 +1
      backspace(fl); backspace(fl); backspace(fl)
      write(fl,'(A)') "  },"
      write(fl,"(A ,': {')") '  "'//trim(adjustl(run_str))//'"'
   else
      open(newunit=fl,file='run_log.json',status='replace',action='write')
      write(fl,'(A)') "{"
      write(fl,"(A ,': {')") '  "run 1"'
   endif

   write(fl,'(A," : ",I0,",")') '    "itime0"',ifirst
   write(fl,'(A," : ",g0,",")') '    "t0"', t0
   write(fl,'(A," : ",g0)') '    "dt"',dt
   write(fl,'(A)') "  }"
   write(fl,'(A)') "}"

   close(fl)
end subroutine