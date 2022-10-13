module poisson_mod
    use decomp_2d
    integer :: isource = 0
    real(mytype) :: source_cnst,a,b
    real(mytype), dimension(:,:,:), allocatable :: source
    real(mytype), dimension(:,:,:), allocatable :: rhs

    character(len=*), parameter :: io_name = "solution-io"
    integer :: ioxdmf
    logical, save :: filenamedigits = .true.
    character(len=9) :: ifilenameformat = '(I7.7)'

    contains
    subroutine initialise_poisson
        use MPI
        use variables, only : nx, ny, nz, nxm, nym, nzm
        use variables, only : p_row, p_col, nstat,nvisu
        use var
        use decomp_2d_poisson
        use decomp_2d_io
        use dbg_schemes, only : abs_prec, sqrt_prec
        implicit none

        integer :: ierr, i, i_loc, j
        integer :: x_ind1, y_ind1, x_ind2, y_ind2
        real(mytype) :: eps1x, eps1y, eps2x, eps2y, eps1, eps2
        real(mytype) :: p1(2), p2(2)

        real(mytype) :: x,y
    
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    
        call read_inputs()

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
        call init_variables
        call schemes()

        call decomp_2d_poisson_init()
        call decomp_info_init(nxm,nym,nzm,phG)

        call system("mkdir -p data")

        allocate(source(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize))
        allocate(rhs(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize))

        if (isource.eq.1) then
            ! do j = ph1%zst(1),ph1%zen(1)
            !     x = (real(j,kind=mytype) - zpfive)*dx
            !     source(j,:,:) = a*x + b
            ! enddo
            do j = ph1%zst(2),ph1%zen(2)
                x = (real(j,kind=mytype) - zpfive)*dy
                source(:,j,:) = a*x + b
            enddo
        elseif (isource.eq.2) then
            p1 = [0.4,0.7]
            p2 = [0.6,0.7]

            source = zero
            eps1x=1e20
            eps1y=1e20
            eps2x=1e20
            eps2y=1e20
            
            do i = ph1%zst(1),ph1%zen(1)
                x = (real(i,kind=mytype) - zpfive)*dx
                eps1 = abs_prec(x - p1(1))
                if (eps1<eps1x) then
                    eps1x = eps1
                    x_ind1 = i
                endif
                eps2 = abs_prec(x - p2(1))
                if (eps2<eps2x) then
                    eps2x = eps2
                    x_ind2 = i
                endif
            enddo

            do j = ph1%zst(2),ph1%zen(2)
                y = 0.5*(yp(j) + yp(j+1))
                eps1 = abs_prec(y - p1(2))
                if (eps1<eps1y) then
                    eps1y = eps1
                    y_ind1 = j
                endif
                eps2 = abs_prec(y - p2(2))
                if (eps2<eps2y) then
                    eps2y = eps2
                    y_ind2 = j
                endif
            enddo

            eps1 = sqrt_prec(eps1x*eps1x + eps1y*eps1y)
            call MPI_Allreduce(eps1,eps1x,1,real_type,MPI_MIN,DECOMP_2D_COMM_CART_Z,ierr)
            if (eps1x == eps1) source(x_ind1,y_ind1,:) = -one

            eps2 = sqrt_prec(eps2x*eps2x + eps2y*eps2y)
            call MPI_Allreduce(eps2,eps2x,1,real_type,MPI_MIN,DECOMP_2D_COMM_CART_Z,ierr)
            if (eps2x == eps2) source(x_ind2,y_ind2,:) = one
        else
            write(*,*) "problems pal"
            call MPI_Abort(MPI_COMM_WORLD,1,ierr)
        endif

        rhs(:,:,:) = source(:,:,:)
        call post_init

    end subroutine initialise_poisson

    subroutine finalize_poisson
        use MPI
        use decomp_2d
        use decomp_2d_io, only : decomp_2d_io_finalise
        implicit none
        integer :: ierr

        call decomp_2d_io_finalise()
        call decomp_2d_finalize
        CALL MPI_FINALIZE(ierr)
    end subroutine

    subroutine read_inputs
        use decomp_2d
        use variables
        use param
        use decomp_2d_poisson

        implicit none
        integer :: unit
        namelist/basic_params/p_row, p_col, nx, ny, nz,&
                             istret, beta, xlx, yly, zlz,&
                             isource
        namelist/boundaries/nclx1, nclxn, ncly1, nclyn, nclz1, nclzn
        namelist/poissonlist/ a,b,bcx,bcy,bcz
        NAMELIST/NumOptions/ ifirstder, isecondder, itimescheme, iimplicit, &
                                nu0nu, cnu, ipinter
        
                                open(newunit=unit,file='poisson.i3d',action='read',status='old')
        
        read(unit,nml=basic_params) ; rewind(unit)
        read(unit,nml=boundaries) ; rewind(unit)
        read(unit,nml=NumOptions) ; rewind(unit)
        if (isource .eq.1.or.isource .eq.2) &
            read(unit,nml=poissonlist) ; rewind(unit)

        close(unit)

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
        
        nstat = 1
        nvisu = 1
        nprobe = 1

        dx=xlx/real(nxm,mytype)
        dy=yly/real(nym,mytype)
        dz=zlz/real(nzm,mytype)

        dx2 = dx * dx
        dy2 = dy * dy
        dz2 = dz * dz

        xnu=one/re

    end subroutine
    subroutine post_init
        use decomp_2d_io
        implicit none
        call decomp_2d_init_io(io_name)
        call decomp_2d_register_variable(io_name, "source", 3, 0, 0, mytype,opt_decomp=ph1)
        call decomp_2d_register_variable(io_name, "rhs", 3, 0, 0, mytype,opt_decomp=ph1)

    end subroutine

    
  !
  ! Write a snapshot
  !
  subroutine write_snapshot()

    use decomp_2d, only : transpose_z_to_y, transpose_y_to_x
    use decomp_2d, only : mytype, xsize, ysize, zsize
    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_start_io

    use param, only : nrhotime, ilmn, iscalar, ioutput, irestart

    use variables, only : sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6
    use variables, only : sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6
    use variables, only : sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6
    use variables, only : numscalar

    use var, only : pp1, ta1, di1, nxmsize
    use var, only : pp2, ppi2, dip2, ph2, nymsize
    use var, only : ppi3, dip3, ph3, nzmsize
    use var, only : npress

    use tools, only : rescale_pressure

    implicit none

    !! inputs
    ! Local variables
    integer :: is
    integer :: ierr
    character(len=30) :: scname
    integer :: mode
    logical, save :: outloc_init = .false.
    logical :: dir_exists

#ifdef ADIOS2
    call decomp_2d_start_io(io_name, "data")
#endif
    
    ! Write XDMF header
    call write_xdmf_header(".", "snapshot")

    ! Write velocity
    call write_field(source, ".",'source')
    call write_field(rhs, ".",'rhs')

    call write_xdmf_footer()
  end subroutine write_snapshot

  subroutine end_snapshot(itime)

    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_end_io
    use param, only : istret, xlx, yly, zlz
    use variables, only : nx, ny, nz, beta
    use var, only : dt,t

    implicit none

    integer, intent(in) :: itime

    character(len=32) :: fmt2, fmt3, fmt4
    integer :: is
    integer :: ierr
    
    ! Write XDMF footer
    call write_xdmf_footer()

#ifdef ADIOS2
    call decomp_2d_end_io(io_name, "data")
#endif
    
  end subroutine end_snapshot

  !
  ! Write the header of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_header(pathname, filename)

    use variables, only : nvisu, yp
    use param, only : dx,dy,dz,istret,zpfive
    use decomp_2d, only : mytype, nrank, xszV, yszV, zszV, ystV
    use variables, only : nx, ny, nz, nxm, nym, nzm

    implicit none

    ! Arguments
    character(len=*), intent(in) :: pathname, filename

    ! Local variables
    integer :: i,k,j
    real(mytype) :: xc(nxm), zc(nzm), yc(nym)

    if (nrank.eq.0) then
      OPEN(newunit=ioxdmf,file="./data/"//gen_snapshotname(pathname, filename, "xdmf"))

      write(ioxdmf,'(A22)')'<?xml version="1.0" ?>'
      write(ioxdmf,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(ioxdmf,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(ioxdmf,*)'<Domain>'
      if (istret.ne.0) then
        do i=1,nxm
          xc(i) = (real(i-1,mytype) + zpfive)*dx*nvisu
        enddo
        do k=1,nzm
          zc(k) = (real(k-1,mytype) + zpfive)*dz*nvisu
        enddo
        do j=1,nym,nvisu
            yc(k) = zpfive*(yp(j) + yp(j+1))
        enddo

        write(ioxdmf,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
        write(ioxdmf,*)'        Dimensions="',nzm,nym,nxm,'">'
        write(ioxdmf,*)'    </Topology>'
        write(ioxdmf,*)'    <Geometry name="geo" Type="VXVYVZ">'
          write(ioxdmf,*)'        <DataItem Dimensions="',nxm,'" NumberType="Float" Precision="4" Format="XML">'
          write(ioxdmf,*)'        ',xc(:)

        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'        <DataItem Dimensions="',nym,'" NumberType="Float" Precision="4" Format="XML">'
        write(ioxdmf,*)'        ',yc(:)

        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'        <DataItem Dimensions="',nzm,'" NumberType="Float" Precision="4" Format="XML">'
        write(ioxdmf,*)'        ',zc(:)

        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      else
        write(ioxdmf,*)'    <Topology name="topo" TopologyType="3DCoRectMesh"'
        write(ioxdmf,*)'        Dimensions="',nzm,nym,nxm,'">'
        
        write(ioxdmf,*)'    </Topology>'
        write(ioxdmf,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
        write(ioxdmf,*)'        <!-- Origin -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        write(ioxdmf,*)'        0.0 0.0 0.0'
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'        <!-- DxDyDz -->'
        write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
        write(ioxdmf,*)'        ',nvisu*dz,nvisu*dy,nvisu*dx
        write(ioxdmf,*)'        </DataItem>'
        write(ioxdmf,*)'    </Geometry>'
      endif
      write(ioxdmf,*)'    <Grid Name="1" GridType="Uniform">'
      write(ioxdmf,*)'        <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(ioxdmf,*)'        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    endif
  end subroutine write_xdmf_header

  !
  ! Write the footer of the XDMF file
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_xdmf_footer()

    use decomp_2d, only : nrank

    implicit none

    if (nrank.eq.0) then
      write(ioxdmf,'(/)')
      write(ioxdmf,*)'    </Grid>'
      write(ioxdmf,*)'</Domain>'
      write(ioxdmf,'(A7)')'</Xdmf>'
      close(ioxdmf)
    endif

  end subroutine write_xdmf_footer

  !
  ! Write the given field for visualization
  ! Adapted from https://github.com/fschuch/Xcompact3d/blob/master/src/visu.f90
  !
  subroutine write_field(f1, pathname, filename, skip_ibm, flush)

    use mpi
    
    use var, only : ep1
    use var, only : zero, one
    use var, only : uvisu
    use param, only : iibm
    use decomp_2d, only : mytype, xsize, xszV, yszV, zszV, ph1
    use decomp_2d, only : nrank, fine_to_coarseV
    use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_write_plane
    use variables, only : nx, ny, nz, nxm, nym, nzm

    implicit none

    real(mytype), intent(in), dimension(ph1%zsz(1),ph1%zsz(2),ph1%zsz(3)) :: f1
    character(len=*), intent(in) :: pathname, filename
    logical, optional, intent(in) :: skip_ibm, flush

    logical :: mpiio, force_flush
    
    integer :: ierr

#ifndef ADIOS2
    mpiio = .true.
#else
    mpiio = .false.
#endif

    if (present(flush)) then
       force_flush = flush
    else
       force_flush = .false.
    end if

    if (nrank.eq.0) then
        write(ioxdmf,*)'        <Attribute Name="'//filename//'" Center="Node">'
#ifndef ADIOS2
        write(ioxdmf,*)'           <DataItem Format="Binary"'
#else
        write(ioxdmf,*)'           <DataItem Format="HDF"'
#endif
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
        if (output2D.eq.0) then
            write(ioxdmf,*)'            DataType="Float" Precision="4" Endian="little" Seek="0"'
        else
            write(ioxdmf,*)'            DataType="Float" Precision="8" Endian="little" Seek="0"'
        endif
#else
        write(ioxdmf,*)'            DataType="Float" Precision="8" Endian="little" Seek="0"'
#endif
#else
        write(ioxdmf,*)'            DataType="Float" Precision="4" Endian="little" Seek="0"'
#endif
        write(ioxdmf,*)'            Dimensions="',nzm,nym,nxm,'">'
        write(ioxdmf,*)'              '//gen_h5path(gen_filename(pathname, filename, 'bin'), '1')
        write(ioxdmf,*)'           </DataItem>'
        write(ioxdmf,*)'        </Attribute>'
    endif


    if (mpiio .or. (iibm == 2) .or. force_flush) then
        !! XXX: This (re)uses a temporary array for data - need to force synchronous writes.
        uvisu = zero
        
        call decomp_2d_write_one(3,f1,"data",gen_filename(pathname, filename, 'bin'),0,io_name,&
            opt_deferred_writes=.false.,opt_decomp=ph1)
    else
        call decomp_2d_write_one(3,f1,"data",gen_filename(pathname, filename, 'bin'),0,io_name)
    end if

  end subroutine write_field

  function gen_snapshotname(pathname, varname, ext)
    character(len=*), intent(in) :: pathname, varname, ext
#ifndef ADIOS2
    character(len=(len(pathname) + 1 + len(varname) + 1  + len(ext))) :: gen_snapshotname
    gen_snapshotname = gen_filename(pathname, varname, ext)
#else
    character(len=(len(varname) + 1 + len(num) + 1 + len(ext))) :: gen_snapshotname
    write(gen_snapshotname, "(A)") varname//'.'//ext
#endif
  end function gen_snapshotname
  
  function gen_filename(pathname, varname, ext)

    character(len=*), intent(in) :: pathname, varname, ext
#ifndef ADIOS2
    character(len=(len(pathname) + 1 + len(varname) + 1 + len(ext)+1)) :: gen_filename
    gen_filename = adjustl(pathname//'/'//varname//'.'//ext)
#else
    character(len=len(varname)) :: gen_filename
    write(gen_filename, "(A)") varname
#endif
    
  end function gen_filename
  function gen_h5path(filename, num)

    character(len=*), intent(in) :: filename, num
#ifndef ADIOS2
    character(len=*), parameter :: path_to_h5file = "./"
    character(len=(len(path_to_h5file) + len(filename))) :: gen_h5path
    write(gen_h5path, "(A)") path_to_h5file//filename
#else
    character(len=*), parameter :: path_to_h5file = "../data.hdf5:/Step"
    character(len=(len(path_to_h5file) + len(num) + 1+ len(filename))) :: gen_h5path
    write(gen_h5path, "(A)") path_to_h5file//num//"/"//filename
#endif
    
  end function gen_h5path
  
end module

program poisson_test
    use poisson_mod
    use decomp_2d_poisson

    implicit none

    call initialise_poisson
    
    call poisson(rhs)

    call write_snapshot

    call finalize_poisson

end program