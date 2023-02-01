!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!###########################################################################
!
!  SUBROUTINE: parameter
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the
!              simulation.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!###########################################################################
subroutine parameter(input_i3d)

  use mpi
  
  use iso_fortran_env

  use param
  use variables
  use complex_geometry
  use decomp_2d
  use ibm_param
  use stats, only: 
  use dbg_schemes, only: sin_prec, cos_prec

  use var, only : dphi1

  use lockexch, only : pfront
  use stats, only : h_quads, spectra_level, autocorr_max_sep,autocorr_xlocs,&
                    spectra_corr_nlocs,spectra_corr_ylocs
  use probes, only : nprobes, nlineprobes, setup_probes, flag_all_digits,&
                     flag_extra_probes, xyzprobes, zlineprobes, probe_freq
  use visu, only : output2D
  use forces, only : iforces, nvol, xld, xrd, yld, yud!, zld, zrd

  implicit none

  character(len=80), intent(in) :: input_i3d
  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total, iostat, ierr

  NAMELIST /BasicParam/ p_row, p_col, nx, ny, nz, istret, beta, xlx, yly, zlz, &
       itype, iin, re, u1, u2, init_noise, inflow_noise, &
       dt, ifirst, ilast, &
       numscalar, iibm, ilmn, &
       ilesmod, iscalar, &
       nclx1, nclxn, ncly1, nclyn, nclz1, nclzn, &
       ivisu, ipost, &
       gravx, gravy, gravz, &
       cpg, idir_stream, &
       ifilter, C_filter, iturbine, use_center
  NAMELIST /NumOptions/ ifirstder, isecondder, itimescheme, iimplicit, &
       nu0nu, cnu, ipinter
  NAMELIST /InOutParam/ irestart, icheckpoint, ioutput, nvisu, ilist, iprocessing, &
       ninflows, ntimesteps, inflowpath, ioutflow, output2D, nprobes, log_cputime,&
       nlineprobes
  NAMELIST /Statistics/ wrotation,spinup_time, nstat, initstat, &
            istatcalc, istatbudget,istatpstrain,istatlambda2, initstat2,&
            istatquadrant, nquads, istatflatness, istatspectra,spectra_level,&
            istatautocorr, autocorr_max_sep,autocorr_xlocs, spectra_corr_nlocs
  NAMELIST /ProbesParam/ flag_all_digits, flag_extra_probes, xyzprobes, zlineprobes, probe_freq
  NAMELIST /ScalarParam/ sc, ri, uset, cp, &
       nclxS1, nclxSn, nclyS1, nclySn, nclzS1, nclzSn, &
       scalar_lbound, scalar_ubound, sc_even, sc_skew, &
       alpha_sc, beta_sc, g_sc, Tref
  NAMELIST /LESModel/ jles, smagcst, smagwalldamp, nSmag, walecst, maxdsmagcst, iwall
  NAMELIST /WallModel/ smagwalldamp
  NAMELIST /Tripping/ itrip,A_tr,xs_tr_tbl,ys_tr_tbl,ts_tr_tbl,x0_tr_tbl, y0_tr_tbl,x_tr_stop
  NAMELIST /ibmstuff/ cex,cey,cez,ra,nobjmax,nraf,nvol,iforces, npif, izap, ianal, imove, thickness, chord, omega ,ubcx,ubcy,ubcz,rads, c_air
  NAMELIST /ForceCVs/ xld, xrd, yld, yud!, zld, zrd
  NAMELIST /LMN/ dens1, dens2, prandtl, ilmn_bound, ivarcoeff, ilmn_solve_temp, &
       massfrac, mol_weight, imultispecies, primary_species, &
       Fr, ibirman_eos
  NAMELIST /ABL/ z_zero, iwallmodel, k_roughness, ustar, dBL, &
       imassconserve, ibuoyancy, iPressureGradient, iCoriolis, CoriolisFreq, &
       istrat, idamping, iheight, TempRate, TempFlux, itherm, gravv, UG, T_wall, T_top, ishiftedper, iconcprec, pdl 
  NAMELIST /CASE/ tgv_twod, pfront
  NAMELIST/ALMParam/iturboutput,NTurbines,TurbinesPath,NActuatorlines,ActuatorlinesPath,eps_factor,rho_air
  NAMELIST/ADMParam/Ndiscs,ADMcoords,C_T,aind,iturboutput,rho_air
  NAMELIST/TBLRecy/plane_location, t_avg1,t_avg2, t_recy1, t_recy2, iaccel,t_trip
  NAMELIST/tanhAccelTBL/U_ratio, accel_centre, alpha_accel
  NAMELIST/tanhAccelCubicTBL/U_ratio, accel_centre, alpha_accel, iaccel_thresh
  NAMELIST/fileAccel/accel_file
  NAMELIST/bodyForce/ibodyforces,ibftype
  NAMELIST/linearBodyF/linear_amp,linear_ext
  NAMELIST/tempAccel/itempaccel, iacceltype, istatout,ispectout,ispectstart
  NAMELIST/linear_prof/Re_ratio, t_start, t_end
  NAMELIST/spatial_equiv/U_ratio, accel_centre, alpha_accel
  NAMELIST/hquadrant/h_quads
  NAMELIST/spectra_corr/spectra_corr_ylocs

#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter start'
#endif

  if (nrank==0) then
     write(*,*) '==========================================================='
     write(*,*) '======================Xcompact3D==========================='
     write(*,*) '===Copyright (c) 2018 Eric Lamballais and Sylvain Laizet==='
     write(*,*) '===Modified by Felipe Schuch and Ricardo Frantz============'
     write(*,*) '===Modified by Paul Bartholomew, Georgios Deskos and======='
     write(*,*) '===Sylvain Laizet -- 2018- ================================'
     write(*,*) '==========================================================='
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif

  call parameter_defaults()

  !! Read parameters
  open(10, file=input_i3d)

  !! These are the 'essential' parameters
  read(10, nml=BasicParam); rewind(10)
  read(10, nml=NumOptions); rewind(10)
  read(10, nml=InOutParam); rewind(10)
  read(10, nml=Statistics); rewind(10)

  if (istatquadrant) then
   allocate(h_quads(nquads))
   read(10,nml=hquadrant); rewind(10)
  endif

  if (istatspectra) then
   if (spectra_level.lt.1) then
      write(*,*) "spectra_level must be set"
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif
   if (spectra_level.eq.3) then
      if (spectra_corr_nlocs<1) then
         if (nrank ==0 )write(*,*) "spectra_level must be set"
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
      allocate(spectra_corr_ylocs(spectra_corr_nlocs))
      read(10,nml=spectra_corr); rewind(10)
   endif
  endif

  if (istatautocorr) then
   if(autocorr_max_sep<0) then
      write(*,*) "max sep must be set"
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif
   if(autocorr_xlocs<0) then
      write(*,*) "xlocs must be set"
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif
  endif

  if (iibm.ne.0) then
     read(10, nml=ibmstuff); rewind(10)
  endif
  if (nprobes.gt.0.or.nlineprobes.gt.0) then
     call setup_probes()
     read(10, nml=ProbesParam); rewind(10)
  endif
  if (iforces.eq.1) then
     allocate(xld(nvol), xrd(nvol), yld(nvol), yud(nvol))!, zld(nvol), zrd(nvol))
     read(10, nml=ForceCVs); rewind(10)
  endif
  
  !! Set Scalar BCs same as fluid (may be overridden) [DEFAULT]
  nclxS1 = nclx1; nclxSn = nclxn
  nclyS1 = ncly1; nclySn = nclyn
  nclzS1 = nclz1; nclzSn = nclzn
  
  if (numscalar.ne.0) then
     iscalar = 1

     !! Allocate scalar arrays and set sensible defaults
     allocate(massfrac(numscalar))
     allocate(mol_weight(numscalar))
     massfrac(:) = .FALSE.
     mol_weight(:) = one
     allocate(sc(numscalar), ri(numscalar), uset(numscalar), cp(numscalar))
     ri(:) = zero
     uset(:) = zero
     cp(:) = zero
     if (iimplicit.gt.0) then
        allocate(xcst_sc(numscalar))
        xcst_sc(:) = zero
        allocate(alpha_sc(numscalar,2), beta_sc(numscalar,2), g_sc(numscalar,2))
        ! Default scalar BC : dirichlet BC, zero value
        alpha_sc = one
        beta_sc = zero
        g_sc = zero
     endif

     ! In case of symmetry, scalars are even by default
     allocate(sc_even(numscalar))
     sc_even(:) = .true.

     ! Skew-symmetric convection of scalars, off by default
     allocate(sc_skew(numscalar))
     sc_skew(:) = .false.

     allocate(scalar_lbound(numscalar), scalar_ubound(numscalar))
     scalar_lbound(:) = -huge(one)
     scalar_ubound(:) = huge(one)
  endif

  if (ilmn) then
     if (istret.ne.0) then
        if (nrank.eq.0) then
           print *, "WARNING: LMN solver does not currently support stretching!"
           stop
        endif
     endif
     read(10, nml=LMN); rewind(10)

     do is = 1, numscalar
        if (massfrac(is)) then
           imultispecies = .TRUE.
        endif
     enddo

     if (imultispecies) then
        if (primary_species.lt.1) then
           if (nrank==0) then
              write(*,*)  "Error: you must set a primary species for multispecies flow"
              write(*,*)  "       solver will enforce Y_p = 1 - sum_s Y_s, s != p."
              stop
           endif
        else if (.not.massfrac(primary_species)) then
           if (nrank==0) then
              write(*,*)  "Error: primary species must be a massfraction!"
           endif
        endif
     endif
  endif
  if (numscalar.ne.0) then
     read(10, nml=ScalarParam); rewind(10)
  endif
  ! !! These are the 'optional'/model parameters
  ! read(10, nml=ScalarParam)
  if(ilesmod==0) then
     nu0nu=four
     cnu=0.44_mytype
  endif
  if(ilesmod.ne.0) then
     read(10, nml=LESModel); rewind(10)
  endif
  
  if (itype.eq.itype_tbl) then
     read(10, nml=Tripping); rewind(10)
  endif

  if (itype.eq.itype_tbl_recy) then
      read(10, nml=TBLRecy); rewind(10)
      read(10, nml=Tripping); rewind(10)
      if(iaccel.eq.1) then
         read(10, nml=tanhAccelTBL); rewind(10)
      else if(iaccel.eq.2) then
         read(10, nml=tanhAccelCubicTBL); rewind(10)
      else if (iaccel == 3) then
         read(10, nml=fileAccel); rewind(10)
      endif

   endif
  
   if (x_tr_stop .le. zero) then
    x_tr_stop = xlx
  endif

  if(itype==itype_channel) then
   read(10,nml=bodyForce,iostat=iostat); rewind(10)
   if (iostat ==0.and.ibodyforces.eq.1) then
      if (ibftype.eq.1) then
         read(10,nml=linearBodyF); rewind(10)
      endif
   else
      ibodyforces = 0
   endif

   read(10,nml=tempAccel,iostat=iostat); rewind(10)
   if (iostat==0.and.itempaccel==1) then
      if (iacceltype==1) then
         read(10,nml=linear_prof); rewind(10)
      else if (iacceltype==2) then
         read(10,nml=spatial_equiv); rewind(10)
      endif

      if (istatout<1) then
         write(*,*) "If temporal acceleration is used istatout "//&
                     "must be specified and valid"
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
      if (ispectout<1.and.istatspectra) then
         write(*,*) "If temporal acceleration is used "//&
                     "ispectout must be specified and valid"
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif

      if (mod(ispectout,istatout)/= 0.and.istatspectra) then
         write(*,*) "ispectout must be multiple of istatout"
            call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif

      if (ispectstart<1.and.istatspectra) then
         write(*,*) "If temporal acceleration is used "//&
                     "ispectstart must be specified and valid"
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif
      
   endif
  endif

   if (itempaccel == 0) then
      istatout = icheckpoint
      ispectout = icheckpoint
   endif

  if (itype.eq.itype_abl) then
     read(10, nml=ABL); rewind(10)
  endif
  if (iturbine.eq.1) then
     read(10, nml=ALMParam); rewind(10)
  else if (iturbine.eq.2) then
     read(10, nml=ADMParam); rewind(10)
  endif
  ! read(10, nml=TurbulenceWallModel)
  read(10, nml=CASE); rewind(10) !! Read case-specific variables
  close(10)

  ! allocate(sc(numscalar),cp(numscalar),ri(numscalar),group(numscalar))

  if (nclx1.eq.0.and.nclxn.eq.0) then
     nclx=.true.
     nxm=nx
  else
     nclx=.false.
     nxm=nx-1
  endif
  if (ncly1.eq.0.and.nclyn.eq.0) then
     ncly=.true.
     nym=ny
  else
     ncly=.false.
     nym=ny-1
  endif
  if (nclz1.eq.0.and.nclzn.eq.0) then
     nclz=.true.
     nzm=nz
  else
     nclz=.false.
     nzm=nz-1
  endif

  dx=xlx/real(nxm,mytype)
  dy=yly/real(nym,mytype)
  dz=zlz/real(nzm,mytype)

  dx2 = dx * dx
  dy2 = dy * dy
  dz2 = dz * dz

  xnu=one/re
  !! Constant pressure gradient, re = Re_tau -> use to compute Re_centerline
  if (cpg) then
    re_cent = (re/0.116_mytype)**(1.0_mytype/0.88_mytype)
    xnu = one/re_cent ! viscosity based on Re_cent to keep same scaling as CFR
    !
    fcpg = two/yly * (re/re_cent)**2
  end if

  if (ilmn) then
     if (ivarcoeff) then
        npress = 2 !! Need current pressure and previous iterate
     else
        npress = 1
     endif
  endif

  ! 2D snapshot is not compatible with coarse visualization
  if (output2D.ne.0) nvisu = 1
#ifdef ADIOS2
  if (nvisu .ne. 1) then
     if (nrank .eq. 0) then
        print *, "ADIOS2 output is not compatible with coarse visualisation"
        print *, "disabling coarse visualisation"
        print *, "To compress the IO, see ADIOS2 options"
     endif
     nvisu = 1
  endif
#if defined(DOUBLE_PREC) && defined(SAVE_SINGLE)
  print *, "ADIOS2 does not support mixing the simulation and output precision"
  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
#endif
#endif

  if (iimplicit.ne.0) then
     if ((itimescheme==5).or.(itimescheme==6)) then
        if (nrank==0) write(*,*) 'Error: implicit Y diffusion not yet compatible with RK time schemes'
        stop
     endif
     if (isecondder==5) then
        if (nrank==0) write(*,*)  "Warning : support for implicit Y diffusion and isecondder=5 is experimental"
     endif
     if (iimplicit==1) then
        xcst = dt * xnu
     else if (iimplicit==2) then
        xcst = dt * xnu * half
     else
        if (nrank==0) write(*,*)  'Error: wrong value for iimplicit ', iimplicit
        stop
     endif
     if (iscalar.eq.1) xcst_sc = xcst / sc
  endif

  if (itype==itype_tbl.and.A_tr .gt. zero.and.nrank==0)  write(*,*)  "TBL tripping is active"

  anglex = sin_prec(pi*angle/onehundredeighty)
  angley = cos_prec(pi*angle/onehundredeighty)
  !###########################################################################
  ! Log-output
  !###########################################################################
  if (nrank==0) call system('mkdir data out probes 2> /dev/null')

#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter input.i3d done'
#endif
  if (nrank==0) then
     print *,'==========================================================='
     if (itype.eq.itype_user) then
        print *,'User-defined simulation'
     elseif (itype.eq.itype_lockexch) then
        print *,'Simulating lock-exchange'
     elseif (itype.eq.itype_tgv) then
        print *,'Simulating TGV'
     elseif (itype.eq.itype_channel) then
        print *,'Simulating channel'
     elseif (itype.eq.itype_hill) then
        print *,'Simulating periodic hill'
     elseif (itype.eq.itype_cyl) then
        print *,'Simulating cylinder'
     elseif (itype.eq.itype_dbg) then
        print *,'Debug schemes'
     elseif (itype.eq.itype_mixlayer) then
        print *,'Mixing layer'
     elseif (itype.eq.itype_jet) then
        print *,'Jet'
     elseif (itype.eq.itype_tbl) then
        print *,'Turbulent boundary layer'
      elseif (itype.eq.itype_tbl_recy) then
      print *,'Turbulent boundary layer using recycling method'
     elseif (itype.eq.itype_abl) then
        print *,'Atmospheric boundary layer'
     elseif (itype.eq.itype_uniform) then
        print *,'Uniform flow'
     elseif (itype.eq.itype_sandbox) then
           print *,'Sandbox'
     else
        print *,'Unknown itype: ', itype
        stop
     endif
     print *,'==========================================================='
     if (itype.eq.itype_channel) then
       if (.not.cpg) then
         write(*,*) 'Channel forcing with constant flow rate (CFR)'
         write(*,"(' Re_cl                  : ',F17.3)") re
       else 
         write(*,*) 'Channel forcing with constant pressure gradient (CPG)'
         write(*,"(' Re_tau                 : ',F17.3)") re
         write(*,"(' Re_cl (estimated)      : ',F17.3)") re_cent
         write(*,"(' fcpg                   : ',F17.8)") fcpg
       end if
     else
       write(*,"(' Reynolds number Re     : ',F17.3)") re
     endif
     write(*,"(' xnu                    : ',F17.8)") xnu
     write(*,*) '==========================================================='
     write(*,"(' p_row, p_col           : ',I9, I8)") p_row, p_col
     write(*,*) '==========================================================='
     write(*,"(' Time step dt           : ',F17.8)") dt
     !
     if (itimescheme.eq.1) then
       !print *,'Temporal scheme        : Forwards Euler'
       write(*,"(' Temporal scheme        : ',A20)") "Forwards Euler"
     elseif (itimescheme.eq.2) then
       !print *,'Temporal scheme        : Adams-bashforth 2'
       write(*,"(' Temporal scheme        : ',A20)") "Adams-bashforth 2"
     elseif (itimescheme.eq.3) then
       !print *,'Temporal scheme        : Adams-bashforth 3'
       write(*,"(' Temporal scheme        : ',A20)") "Adams-bashforth 3"
     elseif (itimescheme.eq.4) then
       !print *,'Temporal scheme        : Adams-bashforth 4'
       write(*,"(' Temporal scheme        : ',A20)") "Adams-bashforth 4"
       print *,'Error: Adams-bashforth 4 not implemented!'
       stop
     elseif (itimescheme.eq.5) then
       !print *,'Temporal scheme        : Runge-kutta 3'
       write(*,"(' Temporal scheme        : ',A20)") "Runge-kutta 3"
     elseif (itimescheme.eq.6) then
       !print *,'Temporal scheme        : Runge-kutta 4'
       write(*,"(' Temporal scheme        : ',A20)") "Runge-kutta 4"
       print *,'Error: Runge-kutta 4 not implemented!'
       stop
     else
       print *,'Error: itimescheme must be specified as 1-6'
       stop
     endif
     !
     if (iimplicit.ne.0) then
       if (iimplicit.eq.1) then
         write(*,"('                          ',A40)") "With backward Euler for Y diffusion"
       else if (iimplicit.eq.2) then
         write(*,"('                          ',A40)") "With CN for Y diffusion"
       endif
     endif
     !
     if (ilesmod.ne.0) then
       write(*,*) '                   : DNS'
     else
       if (jles==1) then
          write(*,*) '                   : Phys Smag'
       else if (jles==2) then
          write(*,*) '                   : Phys WALE'
       else if (jles==3) then
          write(*,*) '                   : Phys dyn. Smag'
       else if (jles==4) then
          write(*,*) '                   : iSVV'
       else
       endif
     endif
     write(*,*) '==========================================================='
     write(*,"(' ifirst                 : ',I17)") ifirst
     write(*,"(' ilast                  : ',I17)") ilast
     write(*,*) '==========================================================='
     write(*,"(' Lx                     : ',F17.8)") xlx
     write(*,"(' Ly                     : ',F17.8)") yly
     write(*,"(' Lz                     : ',F17.8)") zlz
     write(*,"(' nx                     : ',I17)") nx
     write(*,"(' ny                     : ',I17)") ny
     write(*,"(' nz                     : ',I17)") nz
     write(*,*) '==========================================================='
     write(*,"(' istret                 : ',I17)") istret
     write(*,"(' beta                   : ',F17.8)") beta
     write(*,*) '==========================================================='
     write(*,"(' nu0nu                  : ',F17.8)") nu0nu
     write(*,"(' cnu                    : ',F17.8)") cnu
     write(*,*) '==========================================================='
     if (iscalar==0) write(*,"(' Scalar                 : ',A17)") "off"
     if (iscalar==1) write(*,"(' Scalar                 : ',A17)") "on"
     write(*,"(' numscalar              : ',I17)") numscalar
     if (iscalar.eq.1) then
       do is=1, numscalar
          write(*,"(' Schmidt number sc(',I2,')  : ',F17.8)") is, sc(is)
          write(*,"(' Richardson n.  ri(',I2,')  : ',F17.8)") is, ri(is)
          if (scalar_lbound(is).gt.-huge(one)) then
             write(*,"(' Lower bound      (',I2,')  : ',F17.8)") is, scalar_lbound(is)
          else
             ! This is the default option, no information printed in the listing
          endif
          if (scalar_ubound(is).lt.huge(one)) then
             write(*,"(' Upper bound      (',I2,')  : ',F17.8)") is, scalar_ubound(is)
          else
             ! This is the default option, no information printed in the listing
          endif
          if (iscalar.eq.1) then
             if (nclxS1.eq.1 .or. nclxSn.eq.1 .or. &
                 nclyS1.eq.1 .or. nclySn.eq.1 .or. &
                 nclzS1.eq.1 .or. nclzSn.eq.1) then
                if (sc_even(is)) then
                   ! This is the default option, no information printed in the listing
                else
                   write(*,"(' Scalar ',I2,' is odd')") is
                endif
             endif
             if (sc_skew(is)) then
                write(*,"(' Scalar ',I2,' with skew-symmetric convection')") is
             else
                ! This is the default option, no information printed in the listing
             endif
          endif
       end do
     endif
     write(*,*) '==========================================================='
     write(*,"(' spinup_time            : ',I17)") spinup_time
     write(*,"(' wrotation              : ',F17.8)") wrotation
     write(*,*) '==========================================================='
     if (iibm==0) write(*,"(' Immersed boundary      : ',A17)") "off"
     if (iibm.gt.1) then
      write(*,"(' Immersed boundary      : ',A17)") "on"
      write(*,"(' iibm                   : ',I17)") iibm
     end if
     if (iibm==1) write(*,*) 'Simple immersed boundary method'
     if (iibm==2) then
       write(*,*) 'Lagrangian polynomial reconstruction'
       write(*,*) '==========================================================='
       write(*,"(' npif                   : ',I17)") npif
       write(*,"(' izap                   : ',I17)") izap
       write(*,"(' nraf                   : ',I17)") nraf
       write(*,"(' nobjmax                : ',I17)") nobjmax
     end if
     write(*,*) '==========================================================='
     write(*,"(' Boundary condition velocity field: ')")
     write(*,"(' nclx1, nclxn           : ',I15,',',I1 )") nclx1,nclxn
     write(*,"(' ncly1, nclyn           : ',I15,',',I1 )") ncly1,nclyn
     write(*,"(' nclz1, nclzn           : ',I15,',',I1 )") nclz1,nclzn
     write(*,*) '==========================================================='
     if ((iscalar==1).or.(ilmn)) then
       write(*,"(' Boundary condition scalar field: ')")
       write(*,"(' nclxS1, nclxSn         : ',I15,',',I1 )") nclxS1,nclxSn
       write(*,"(' nclyS1, nclySn         : ',I15,',',I1 )") nclyS1,nclySn
       write(*,"(' nclzS1, nclzSn         : ',I15,',',I1 )") nclzS1,nclzSn
       write(*,*) '==========================================================='
     endif

#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
     write(*,*) 'Numerical precision: Double, saving in single'
#else
     print *,'Numerical precision: Double'
#endif
#else
     write(*,*) 'Numerical precision: Single'
#endif
     write(*,*) '==========================================================='
     write(*,"(' High and low speed : u1=',F6.2,' and u2=',F6.2)") u1,u2
     write(*,"(' Gravity vector     : (gx, gy, gz)=(',F15.8,',',F15.8,',',F15.8,')')") gravx, gravy, gravz
     if (ilmn) then
        write(*,*)  "LMN                : Enabled"
        if (ivarcoeff) then
           write(*,*)  "LMN-Poisson solver : Variable-coefficient"
        else
           write(*,*)  "LMN-Poisson solver : Constant-coefficient"
        endif
        if (ilmn_bound) then
           write(*,*)  "LMN boundedness    : Enforced"
        else
           write(*,*)  "LMN boundedness    : Not enforced"
        endif
        write(*,"(' dens1 and dens2    : ',F6.2' ',F6.2)") dens1, dens2
        write(*,"(' Prandtl number Re  : ',F15.8)") prandtl
     endif
     if (angle.ne.0.) write(*,"(' Solid rotation     : ',F6.2)") angle
     write(*,*) ' '
     !! Print case-specific information
     if (itype==itype_lockexch) then
        write(*,*)  "Initial front location: ", pfront
     elseif (itype==itype_tgv) then
        write(*,*)  "TGV 2D: ", tgv_twod
     endif
     write(*,*) '==========================================================='
     if (itype.eq.itype_tbl_recy) then
      write(*,*) "Reycling-rescaling TBL parameters:"
      write(*,"(' Recycling plane loc    : ',F17.8)") plane_location
      write(*,"(' time average change 1  : ',F17.8)") t_avg1
      write(*,"(' time average change 2  : ',F17.8)") t_avg2
      write(*,"(' time recy change 1     : ',F17.8)") t_recy1
      write(*,"(' time recy change 2     : ',F17.8)") t_recy2

      if (iaccel.eq.1) then
         write(*,*) "TBL has spatial acceleration. Parameters:"
         write(*,"(' U_ratio                : ',F17.8)") U_ratio
         write(*,"(' Accel alpha            : ',F17.8)") alpha_accel
         write(*,"(' Accel centre           : ',F17.8)") accel_centre
         write(*,*) '==========================================================='

      else
         write(*,*) "TBL doesn't have spatial acceleration."

      endif
     endif
     if (ibodyforces.eq.1) then
      write(*,*) "Streamwise body forces activated:"
      if (ibftype.eq.1) then
         write(*,*) "Linear body force selected"
         write(*,"(' Body force amplitude   : ',F17.8)") linear_amp
         write(*,"(' Body force extent   : ',F17.8)") linear_ext
      endif
      write(*,*) '==========================================================='
   else if (itype==itype_channel) then
      write(*,*) "No streamwise body forces"
   endif
  endif
  
  if (iibm.eq.3) then ! This is only for the Cubic Spline Reconstruction
     npif=npif+1
  endif

#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter done'
#endif

  return
end subroutine parameter

!###########################################################################
!
!  SUBROUTINE: parameter_defaults
! DESCRIPTION: Sets the default simulation parameters.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!###########################################################################
subroutine parameter_defaults()

  use param
  use variables
  use decomp_2d
  use complex_geometry

  use probes, only : nprobes, flag_all_digits, flag_extra_probes, probe_freq, nlineprobes
  use visu, only : output2D
  use forces, only : iforces, nvol
  use stats, only : spectra_level, autocorr_max_sep
  use stats, only : autocorr_xlocs,spectra_corr_nlocs

  implicit none

  integer :: i

  ro = 99999999._mytype
  angle = zero
  u1 = 2
  u2 = 1
  init_noise = zero
  inflow_noise = zero
  iin = 0
  itimescheme = 4
  iimplicit = 0
  istret = 0
  ipinter=3
  beta = 0
  iscalar = 0
  cont_phi = 0
  filepath = './data/'
  irestart = 0
  itime0 = 0
  t0 = zero
  datapath = './data/'

  !! LES stuff
  SmagWallDamp=0
  nSmag=1

  !! IBM stuff
  nraf = 0
  nobjmax = 0

  nvol = 0
  iforces = 0
  itrip = 0
  wrotation = zero
  irotation = 0
  itest=1

  !! Gravity field
  gravx = zero
  gravy = zero
  gravz = zero

  !! LMN stuff
  ilmn = .FALSE.
  ilmn_bound = .TRUE.
  pressure0 = one
  prandtl = one
  dens1 = one
  dens2 = one
  ivarcoeff = .FALSE.
  npress = 1 !! By default people only need one pressure field
  ilmn_solve_temp = .FALSE.
  imultispecies = .FALSE.
  Fr = zero
  ibirman_eos = .FALSE.

  primary_species = -1

  !! Channel
  cpg = .false.
  idir_stream = 1

  !! Filter
  ifilter=0
  C_filter=0.49_mytype

  !! ABL
  z_zero=zpone
  k_roughness=zpfour
  ustar=0.45_mytype
  dBL=250._mytype
  iPressureGradient=1
  iwallmodel=1
  imassconserve=0
  ibuoyancy=1
  iheight=0
  itherm=1
  idamping=0
  gravv=9.81_mytype
  TempRate=-zptwofive/3600_mytype
  TempFlux=0.24_mytype
  UG=[zero,zero,zero]
  ishiftedper=0
  iconcprec=0
  pdl=zero
  !! Turbine modelling
  iturbine=0
  rho_air=one

  !! IO
  ivisu = 1
  ipost = 0
  iprocessing = huge(i)
  initstat = huge(i)
  ninflows=1
  ntimesteps=1
  inflowpath='./'
  ioutflow=0
  output2D = 0
  nprobes=0
  nlineprobes=0
  probe_freq = 0
  log_cputime=.false.

  !! PROBES
  flag_all_digits = .false.
  flag_extra_probes = .false.


  ipost = 0
  iibm=0
  npif = 2
  izap = 1

  imodulo2 = 1

  !! CASE specific variables
  tgv_twod = .FALSE.

  !! TRIPPING
  A_tr=zero
  xs_tr_tbl=1.402033_mytype
  ys_tr_tbl=0.350508_mytype
  ts_tr_tbl=1.402033_mytype
  x0_tr_tbl=3.505082_mytype
  y0_tr_tbl=0.05_mytype
  x_tr_stop=- one
  ! Accel
  iaccel = 0
  U_ratio = one
  t_trip = 600.0_mytype

   ! body forces
  ibodyforces = 0
  linear_amp = zero
  linear_ext = zero

  ! temporal acceleration
  itempaccel = 0
  Re_ratio = one
  t_start = zero
  t_end = zero
  use_center = .false.

  ! statistics
  istatcalc = 1
  istatbudget = .true.
  istatpstrain = .false.
  istatlambda2 = .false.
  istatquadrant = .false.
  istatspectra = .false.
  istatautocorr = .false.
  istatflatness = .false.
  istatout = -1
  ispectout = -1
  spectra_level=-1
  autocorr_max_sep = -1
  autocorr_xlocs = -1
  spectra_corr_nlocs = -1
  ispectstart = -1
end subroutine parameter_defaults
