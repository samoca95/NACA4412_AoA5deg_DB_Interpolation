C-----------------------------------------------------------------------
C  nek5000 user-file template
C
C  user specified routines:
C     - userbc : boundary conditions
C     - useric : initial conditions
C     - uservp : variable properties
C     - userf  : local acceleration term for fluid
C     - userq  : local source term for scalars
C     - userchk: general purpose routine for checking errors etc.
C
C-----------------------------------------------------------------------

      subroutine uservp(ix,iy,iz,eg) ! set variable properties
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      udiff  = 0.0
      utrans = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userf(ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userq(ix,iy,iz,eg) ! set source term
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      qvol   = 0.0
      source = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userbc(ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux = 0.0
      uy = 0.0
      uz = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine useric(ix,iy,iz,ieg) ! set up initial conditions
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux = 0.0
      uy = 0.0
      uz = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userchk()

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
cc      include 'INPUT'                   ! PARAM(12) (DT)
      include 'TOTAL'
cc            include 'SOLN'                    ! T
cc MA:      include 'MASS_DEF'
cc            include 'MASS'                    !BM1 for lambda2
c      include 'NEKUSE'
c      include 'CHKPOINT'
cc      include 'TSTEP'                   ! ISTEP
c      include 'PARALLEL'
      include 'MPIF_DEF'
      include 'mpif.h'

      character*132 inputname1,hdr,temp_n_file

cc MA: files with interpolated fields:
      character*132 interp_out_name
      character*132 ensemble_out_name

      integer*8 nfiles,ssf,ntot,ifld,i,j,k
      real nu,rho,mtimes,ftime
      integer*8 lt
      parameter (lt=lx1*ly1*lz1*lelv)

      integer nelx, nely, nelz
      integer n_g, ngx, ngy, ngz

      real x_min, x_max, y_min, y_max, z_min, z_max

cc MA: interepolating mesh
      integer*8 i_field,i_point,i_buff
      real pts(ldim,lhis)
      real x_pts(IntSize), y_pts(IntSize), z_pts(IntSize)

cc MA: output fields
      integer n_out_fields
      parameter(n_out_fields=7)
      parameter(ntot=lx1*ly1*lz1*lelt)
      real struct(lx1*ly1*lz1*lelt,n_out_fields)
      real vorticity(ntot,3), work1(ntot), work2(ntot)

cc MA: for parallelization
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer rnk,ierr
      integer status(mpi_status_size)
      integer output_counter, counter_b
      integer nnhis,npr      !!! MA: nhis already defined

cc MA: info header
      real mReynolds
      real mtime
      real mdomain_x, mdomain_y, mdomain_z
      integer mnelx, mnely, mnelz
      integer polyx, polyy, polyz

      real    buff_fieldout(lhis,n_out_fields)

cc MA: buffers for mpi_send / mpi_recv
      real    b_vx(lhis), b_vy(lhis), b_vz(lhis), b_la(lhis)
cc MA: omega
      real    b_ox(lhis), b_oy(lhis), b_oz(lhis)

      real    dist(lhis)

      real    rst(lhis*ldim)
      integer rcode(lhis),elid(lhis),proc(lhis)
      common /hpts_r/ rst
      common /hpts_i/ rcode,elid,proc

      integer icalld
      save    icalld
      data    icalld  /0/
      character*80 val1, val2, val3, val4, val5, val6, val7

      integer inth_hpts
      save    inth_hpts

      integer nfail

cc MA: flag to interpolate
      logical ifintp, ifensemble, ifderiv_out

      nfiles=param(68) 		! Number of input files
      nu=param(2)		! Kinematic viscosity
      rho=1.0			! Fluid density

cc MA: interpolation:
      ifintp = .TRUE.

cc MA: number of elements in x and y directions:
      nelx = lelx
      nely = lely
      nelz = lelz

      if (ifintp) then

      if (nid.eq.0) write(*,*) 'IntSize',IntSize,'lhis*np',lhis*np

cc MA: read the grid points in the x direction
      open(unit=111,form='unformatted',file='ZSTRUCT/x.fort')
      read(111) ngx
      if (nid.eq.0) write(*,*) 'ngx',ngx
      read(111) (x_pts(i),i=1,ngx)
      close(111)

cc MA: read the grid points in the y direction
      open(unit=112,form='unformatted',file='ZSTRUCT/y.fort')
      read(112) ngy
      if (nid.eq.0) write(*,*) 'ngy',ngy
      read(112) (y_pts(i),i=1,ngy)
      close(112)

cc MA: read the grid points in the z direction
      open(unit=113,form='unformatted',file='ZSTRUCT/z.fort')
      read(113) ngz
      if (nid.eq.0) write(*,*) 'ngz',ngz
      read(113) (z_pts(i),i=1,ngz)
      close(113)

cc MA: check input mesh
      if (ngx.ne.ngy) then
         if (nid.eq.0) write(*,*) 'ABORT: ngx!=ngy'
         call exitt
      endif
      if (ngx.ne.ngz) then
         if (nid.eq.0) write(*,*) 'ABORT: ngx!=ngz'
         call exitt
      endif
      if (ngy.ne.ngz) then
         if (nid.eq.0) write(*,*) 'ABORT: ngy!=ngz'
         call exitt
      endif

      n_g = ngx

      x_min=1E10
      y_min=1E10
      z_min=1E10
      x_max=-1E10
      y_max=-1E10
      z_max=-1E10

      if(nid.eq.0) then
        do i=1,n_g
          x_min = min(x_min,x_pts(i))
          x_max = max(x_max,x_pts(i))
          y_min = min(y_min,y_pts(i))
          y_max = max(y_max,y_pts(i))
          z_min = min(z_min,z_pts(i))
          z_max = max(z_max,z_pts(i))
        enddo
        write(*,*) 'x_min,x_max',x_min,x_max
        write(*,*) 'y_min,y_max',y_min,y_max
        write(*,*) 'z_min,z_max',z_min,z_max
      endif


cc      write(*,*) 'nid, np, lp', nid, np, lp

      if (mod(n_g,np).eq.0) then
        nnhis=n_g/np
        npr=0
      else
        npr=np-mod(n_g,np)
        if (nid.lt.np-npr) then
          nnhis = n_g/np+1
        else
          nnhis = n_g/np
        endif
      endif

cc MA: for each id,
cc MA: output local number and total number of interpolating points
cc      write(*,*) 'id, local n intpt (total)',nid,nnhis,n_g

cc MA: assign to local processor

      do i_point=1,nnhis
        if (mod(n_g,np).eq.0) then
          i_buff = nnhis*nid+i_point
        else
          if (nid.lt.np-npr) then
            i_buff = nnhis*nid+i_point
          else
            i_buff = nnhis*nid+np-npr+i_point
          endif
        endif
          pts(1,i_point) = x_pts(i_buff)
          pts(2,i_point) = y_pts(i_buff)
          pts(3,i_point) = z_pts(i_buff)
      enddo

      call rzero(struct,ntot*n_out_fields)

      endif ! ifintp

!!!!!! DO FOR EACH INPUT FILE

      do ssf =1,nfiles

      write(temp_n_file,'(i5.5)') ssf

      inputname1 = 'la2naca_wing?.f'//trim(temp_n_file)
cc      call read_hdr(inputname1,ftime) ! We read header to get the times

      if (nid.eq.0) write(*,*) '**FIELD, file name',ssf,inputname1


cc      call load_field(inputname1)
      call full_restart(inputname1,1)

      if(nid.eq.0) write(*,*) 'vx(1), vx(2), vx(3)',
     &vx(1,1,1,1), vx(2,1,1,1), vx(3,1,1,1)

c      call lambda2(T(1,1,1,1,1))
c      call col2  (T(1,1,1,1,1),bm1,lt)
c      call dssum (T(1,1,1,1,1),nx1,ny1,nz1)
c      call col2  (T(1,1,1,1,1),binvm1,lt)
c
      call comp_vort3(vorticity,work1,work2,vx,vy,vz)

      do j=1,ntot
        struct(j,1) = vx(j,1,1,1)
        struct(j,2) = vy(j,1,1,1)
        struct(j,3) = vz(j,1,1,1)
c        struct(j,4) = T(j,1,1,1,1)
        struct(j,5) = vorticity(j,1)
        struct(j,6) = vorticity(j,2)
        struct(j,7) = vorticity(j,3)
      enddo

      if(ifintp) then

cc MA: INTERPOLATION:
! Do necessary checks before interpolation
      if(icalld.eq.0) then
         if(nid.eq.0) then
            ! Compare number of points in interpolating mesh npoints
            ! with maximum defined in SIZE file lhis
            if(nnhis.gt.lhis) then
               write(6,*) 'Increase lhis to npoints',lhis,nnhis
               call exitt
            endif
            write(6,*) 'found ', nnhis, ' points in interp. mesh'
         endif
         ! Setup for interpolation tool. Use default tolerance of -1
         call intpts_setup(-1.0,inth_hpts)
      endif

!     Compare number of points in interpolating mesh npoints
!     with maximum defined in SIZE file lhis
      if(nnhis.gt.lhis) then
         if(nid.eq.0) write(6,*)
     &        'ABORT: lhis too low, increase in SIZE', nnhis, lhis
         call exitt
      endif

!     Start interpolation
      if(icalld.eq.0) then

!     Conceptually, locate npt points. The data corresponding to each point
!     is whether it is inside an element, closesst to a border, not found.
!     Also identify the processor where the point was found, the element
!     where the point was found, parametric coordinates of the point and
!     the distance squared from found to sought point

         call findpts(inth_hpts,rcode,1,
     &        proc,1,
     &        elid,1,
     &        rst,ndim,
     &        dist,1,
     &        pts(1,1),ndim,
     &        pts(2,1),ndim,
     &        pts(3,1),ndim,nnhis)


!     Check the return code
         do i=1,nnhis
            ! Interpolating point is on boundary or outside the SEM mesh
            if(rcode(i).eq.1) then
               if(dist(i).gt.1e-12) then
                  write(6,'(A,4E15.7)')
     &                 ' WARNING: point on boundary or outside
     &the mesh xy[z]d^2:',
     &                 (pts(k,i),k=1,ndim),dist(i)
               endif
            ! Interpolating ponit is not in the SEM mesh
            elseif(rcode(i).eq.2) then
               nfail = nfail + 1
c               write(6,'(A,3E15.7)')
c     &              ' WARNING: point not within mesh xy[z]: !',
c     &              (pts(k,i),k=1,ndim)
            endif
         enddo

         icalld = 1

      endif ! icalld

      do i_field=1,n_out_fields

       call findpts_eval(inth_hpts,buff_fieldout(1,i_field),1,
     &        rcode,1,
     &        proc,1,
     &        elid,1,
     &        rst,ndim,nnhis,
     &        struct(1,i_field))         ! write NEK calculated fields

      enddo ! each field to be interpolated

!      Write out the interpolated field
      if(nid.eq.0) then

cc MA: generate the file(s) to store the interpolated fields

        interp_out_name = 'intVX.f'//trim(temp_n_file)
        open(unit=137,form='unformatted',file=interp_out_name)
        interp_out_name = 'intVY.f'//trim(temp_n_file)
        open(unit=138,form='unformatted',file=interp_out_name)
        interp_out_name = 'intVZ.f'//trim(temp_n_file)
        open(unit=139,form='unformatted',file=interp_out_name)
c        interp_out_name = 'intLA.f'//trim(temp_n_file)
c        open(unit=140,form='unformatted',file=interp_out_name)
        interp_out_name = 'intOX.f'//trim(temp_n_file)
        open(unit=141,form='unformatted',file=interp_out_name)
        interp_out_name = 'intOY.f'//trim(temp_n_file)
        open(unit=142,form='unformatted',file=interp_out_name)
        interp_out_name = 'intOZ.f'//trim(temp_n_file)
        open(unit=143,form='unformatted',file=interp_out_name)

        mReynolds = 2500  ! Reynolds number
        mdomain_x = 2     ! Domain size (D X)
        mdomain_y = 2     ! Domain size (D Y)
        mdomain_z = 4     ! Domain size (D Z)
        mtime = ftime     ! Time of the snapshot
        mnelx = lx1       ! number of elements
        mnely = ly1
        mnelz = lz1
        polyx = mnelx-1   ! polynomial order
        polyy = mnely-1
        polyz = mnelz-1

!     Parameters to write on the header of the int_fld file
        write(val1,'(1p15e17.9)') mReynolds                     ! Reynolds number
        write(val2,'(1p15e17.9)') mdomain_x,mdomain_y,mdomain_z ! domain size
        write(val3,'(9i9)') mnelx,mnely,mnelz                   ! number of elements
        write(val4,'(9i9)') polyx,polyy,polyz                   ! polynomial order
        write(val5,'(1p15e17.9)') mtime                         ! time of the snapshot
        write(val6,'(9i9)')       nnhis                   ! number of points in int. mesh
        write(val7,'(9i9)')       n_out_fields                  ! number of fields

!     Write header
        write(137) '(Re ='//trim(val1)
     &    //') (Lx, Ly, Lz ='//trim(val2)
     &    //') (nelx, nely, nelz ='//trim(val3)
     &    //') (Polynomial order ='//trim(val4)
     &    //') (time ='//trim(val5)
     &    //') (nnhis ='//trim(val6)
     &    //') (n_out_fields ='//trim(val7)
     &    //')'

!     Write values corresponding to the header
        write(137) mReynolds,
     &    mdomain_x, mdomain_y, mdomain_z,
     &    mnelx    , mnely    , mnelz,
     &    polyx    , polyy    , polyz,
     &    mtime,
     &    nnhis,
     &    n_out_fields

!     Write header
        write(138) '(Re ='//trim(val1)
     &    //') (Lx, Ly, Lz ='//trim(val2)
     &    //') (nelx, nely, nelz ='//trim(val3)
     &    //') (Polynomial order ='//trim(val4)
     &    //') (time ='//trim(val5)
     &    //') (nnhis ='//trim(val6)
     &    //') (n_out_fields ='//trim(val7)
     &    //')'

!     Write values corresponding to the header
        write(138) mReynolds,
     &    mdomain_x, mdomain_y, mdomain_z,
     &    mnelx    , mnely    , mnelz,
     &    polyx    , polyy    , polyz,
     &    mtime,
     &    nnhis,
     &    n_out_fields

!     Write header
        write(139) '(Re ='//trim(val1)
     &    //') (Lx, Ly, Lz ='//trim(val2)
     &    //') (nelx, nely, nelz ='//trim(val3)
     &    //') (Polynomial order ='//trim(val4)
     &    //') (time ='//trim(val5)
     &    //') (nnhis ='//trim(val6)
     &    //') (n_out_fields ='//trim(val7)
     &    //')'

!     Write values corresponding to the header
        write(139) mReynolds,
     &    mdomain_x, mdomain_y, mdomain_z,
     &    mnelx    , mnely    , mnelz,
     &    polyx    , polyy    , polyz,
     &    mtime,
     &    nnhis,
     &    n_out_fields

c!     Write header
c        write(140) '(Re ='//trim(val1)
c     &    //') (Lx, Ly, Lz ='//trim(val2)
c     &    //') (nelx, nely, nelz ='//trim(val3)
c     &    //') (Polynomial order ='//trim(val4)
c     &    //') (time ='//trim(val5)
c     &    //') (nnhis ='//trim(val6)
c     &    //') (n_out_fields ='//trim(val7)
c     &    //')'

c!     Write values corresponding to the header
c        write(140) mReynolds,
c     &    mdomain_x, mdomain_y, mdomain_z,
c     &    mnelx    , mnely    , mnelz,
c     &    polyx    , polyy    , polyz,
c     &    mtime,
c     &    nnhis,
c     &    n_out_fields

!     Write header
        write(141) '(Re ='//trim(val1)
     &    //') (Lx, Ly, Lz ='//trim(val2)
     &    //') (nelx, nely, nelz ='//trim(val3)
     &    //') (Polynomial order ='//trim(val4)
     &    //') (time ='//trim(val5)
     &    //') (nnhis ='//trim(val6)
     &    //') (n_out_fields ='//trim(val7)
     &    //')'

!     Write values corresponding to the header
        write(141) mReynolds,
     &    mdomain_x, mdomain_y, mdomain_z,
     &    mnelx    , mnely    , mnelz,
     &    polyx    , polyy    , polyz,
     &    mtime,
     &    nnhis,
     &    n_out_fields
 
!     Write header
         write(142) '(Re ='//trim(val1)
     &    //') (Lx, Ly, Lz ='//trim(val2)
     &    //') (nelx, nely, nelz ='//trim(val3)
     &    //') (Polynomial order ='//trim(val4)
     &    //') (time ='//trim(val5)
     &    //') (nnhis ='//trim(val6)
     &    //') (n_out_fields ='//trim(val7)
     &    //')'
 
c!     Write values corresponding to the header
         write(142) mReynolds,
     &    mdomain_x, mdomain_y, mdomain_z,
     &    mnelx    , mnely    , mnelz,
     &    polyx    , polyy    , polyz,
     &    mtime,
     &    nnhis,
     &    n_out_fields
 
c!     Write header
         write(143) '(Re ='//trim(val1)
     &    //') (Lx, Ly, Lz ='//trim(val2)
     &    //') (nelx, nely, nelz ='//trim(val3)
     &    //') (Polynomial order ='//trim(val4)
     &    //') (time ='//trim(val5)
     &    //') (nnhis ='//trim(val6)
     &    //') (n_out_fields ='//trim(val7)
     &    //')'
 
c!     Write values corresponding to the header
         write(143) mReynolds,
     &    mdomain_x, mdomain_y, mdomain_z,
     &    mnelx    , mnely    , mnelz,
     &    polyx    , polyy    , polyz,
     &    mtime,
     &    nnhis,
     &    n_out_fields

        interp_out_name = 'intHistory.f'//trim(temp_n_file)
        open(unit=133,access='append',file=interp_out_name)

      endif ! just with the first nid!!

! ! ! TO THINK: Write the "long" header just in History
        output_counter=nnhis

      do i_buff=1,nnhis
        b_vx(i_buff)=0
        b_vy(i_buff)=0
        b_vz(i_buff)=0
      enddo

cc MA: buffer for MPI
      do i_buff=1,nnhis
cc        b_vx(i_buff)=pts(1,i_buff)
cc        b_vy(i_buff)=pts(2,i_buff)
cc        b_vz(i_buff)=pts(3,i_buff)
        b_vx(i_buff)=buff_fieldout(i_buff,1)
        b_vy(i_buff)=buff_fieldout(i_buff,2)
        b_vz(i_buff)=buff_fieldout(i_buff,3)
c        b_la(i_buff)=buff_fieldout(i_buff,4)
         b_ox(i_buff)=buff_fieldout(i_buff,5)
         b_oy(i_buff)=buff_fieldout(i_buff,6)
         b_oz(i_buff)=buff_fieldout(i_buff,7)
      enddo

      if (nid.ne.0) then

        call mpi_send(output_counter,1,nekreal,
     &           0,0,nekcomm,ierr)
        call mpi_send(b_vx,output_counter,nekreal,
     &           0,1,nekcomm,ierr)
        call mpi_send(b_vy,output_counter,nekreal,
     &           0,2,nekcomm,ierr)
        call mpi_send(b_vz,output_counter,nekreal,
     &           0,3,nekcomm,ierr)
c        call mpi_send(b_la,output_counter,nekreal,
c     &           0,4,nekcomm,ierr)
         call mpi_send(b_ox,output_counter,nekreal,
     &           0,5,nekcomm,ierr)
         call mpi_send(b_oy,output_counter,nekreal,
     &           0,6,nekcomm,ierr)
         call mpi_send(b_oz,output_counter,nekreal,
     &           0,7,nekcomm,ierr)

      else

         x_min=minval(b_vx)
         x_max=maxval(b_vx)
         y_min=minval(b_vy)
         y_max=maxval(b_vy)
         z_min=minval(b_vz)
         z_max=maxval(b_vz)

         write(6,*) 'rnk, n pts',nid,output_counter
         write(6,*) 'rnk, minvx, maxvx',nid,x_min,x_max
         write(6,*) 'rnk, minvx, maxvy',nid,y_min,y_max
         write(6,*) 'rnk, minvx, maxvz',nid,z_min,z_max
         do i=1,10
           write(6,*) 'i, b_vx(i)',i,b_vx(i)
         enddo
         do i=output_counter-10,output_counter
           write(6,*) 'i, b_vx(i)',i,b_vx(i)
         enddo

         write(133,'(9i9)') nid,output_counter
         write(137) (b_vx(i),i=1,output_counter)
         write(138) (b_vy(i),i=1,output_counter)
         write(139) (b_vz(i),i=1,output_counter)
c         write(140) (b_la(i),i=1,output_counter)
         write(141) (b_ox(i),i=1,output_counter)
         write(142) (b_oy(i),i=1,output_counter)
         write(143) (b_oz(i),i=1,output_counter)

         do rnk=1,lp-1

              do i_buff=1,nnhis
                b_vx(i_buff)=0
                b_vy(i_buff)=0
                b_vz(i_buff)=0
              enddo

              call mpi_recv(counter_b,1,nekreal,
     &           rnk,0,nekcomm,status,ierr)

              call mpi_recv(b_vx,lhis,nekreal,
     &              rnk,1,nekcomm,status,ierr)
              call mpi_recv(b_vy,lhis,nekreal,
     &              rnk,2,nekcomm,status,ierr)
              call mpi_recv(b_vz,lhis,nekreal,
     &              rnk,3,nekcomm,status,ierr)
c              call mpi_recv(b_la,lhis,nekreal,
c     &              rnk,4,nekcomm,status,ierr)
              call mpi_recv(b_ox,lhis,nekreal,
     &              rnk,5,nekcomm,status,ierr)
              call mpi_recv(b_oy,lhis,nekreal,
     &              rnk,6,nekcomm,status,ierr)
              call mpi_recv(b_oz,lhis,nekreal,
     &              rnk,7,nekcomm,status,ierr)


         x_min=minval(b_vx)
         x_max=maxval(b_vx)
         y_min=minval(b_vy)
         y_max=maxval(b_vy)
         z_min=minval(b_vz)
         z_max=maxval(b_vz)

         write(6,*) 'rnk, n pts',rnk,output_counter
         write(6,*) 'rnk, minvx, maxvx',rnk,x_min,x_max
         write(6,*) 'rnk, minvx, maxvy',rnk,y_min,y_max
         write(6,*) 'rnk, minvx, maxvz',rnk,z_min,z_max
         do i=1,10
           write(6,*) 'i, b_vx(i)',i,b_vx(i)
         enddo
         do i=output_counter-10,output_counter
           write(6,*) 'i, b_vx(i)',i,b_vx(i)
         enddo

         write(133,'(9i9)') rnk,counter_b
         write(137) (b_vx(i),i=1,counter_b)
         write(138) (b_vy(i),i=1,counter_b)
         write(139) (b_vz(i),i=1,counter_b)
c         write(140) (b_la(i),i=1,counter_b)
         write(141) (b_ox(i),i=1,counter_b)
         write(142) (b_oy(i),i=1,counter_b)
         write(143) (b_oz(i),i=1,counter_b)

         enddo  ! all the nid!=0

      endif



      endif ! ifintp

      close(133)

      close(137)
      close(138)
      close(139)
c      close(140)
      close(141)
      close(142)
      close(143)
c      write(6,*) 'Output interporlated field in int_fld, nid:', nid

      enddo ! ssf (input files)


      call exitt

      return
      end


c-----------------------------------------------------------------------

      subroutine load_field(field)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'STATS'           ! Statistics specific variables
cc MA      include 'INPUT_DEF'
      include 'INPUT'           ! if3d
cc MA      include 'SOLN_DEF'
      include 'SOLN'
cc MA      include 'TSTEP_DEF'
      include 'TSTEP'

      character*132 field

      call load_fld(field)

c      call outpost(vx,vy,vz,pr,t,'new')

      return
      end

c-----------------------------------------------------------------------
      subroutine read_hdr(field,mtimee)

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'STATS'           ! Statistics specific variables
cc MA      include 'INPUT_DEF'
      include 'INPUT'           ! if3d
cc MA      include 'SOLN_DEF'
      include 'SOLN'
cc MA      include 'TSTEP_DEF'
      include 'TSTEP'

      character*132 field,hdr,fmt1
      character*10 tmpchar
      integer twdsize,mnelx,mnely,mnelz,nelo,isteps,fid0,nfileoo
      integer stat_gnum
      real mtimee

      open(unit=33,file=field,form='unformatted')
      read(33) hdr
      close(33)

      fmt1 = '(1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,
     &1x,i9,1x,i6,1x,i6,1x,10a)'

      read(hdr,fmt1) twdsize,mnelx,mnely,mnelz,nelo,
     &     stat_gnum,mtimee,isteps,fid0,nfileoo,tmpchar

      return
      end
c-----------------------------------------------------------------------

      subroutine usrdat()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'ZPER'            ! For nelx,nely,nelz - needed for z_average

      return
      end

c-----------------------------------------------------------------------

      subroutine usrdat2()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end

c-----------------------------------------------------------------------

      subroutine usrdat3()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
c----------------------------------------------------------------------

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl(flag)
      logical flag

      call qthermal_ig(flag)

      return
      end
