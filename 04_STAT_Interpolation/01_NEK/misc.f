!!======================================================================
!!
!!   Bunch of miscelineous routines.
!!   misc.f made by Prabal Negi
!!   Routine are not mine...
!!
!!======================================================================
 
      real function step(x)
c     nima 
c     belongs to bla version 2.2
c     for more info see the bla.f file
c     
      real x
      if(x.le.0.02) then
         step=0.
      else
         if(x.le.0.98) then
            step=1./(1.+exp(1./(x-1.)+1./x))
         else
            step=1.
         endif
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine fix_mygll()
c     Set the GLL points as defined in the CASENAME.grid file

      include 'SIZE'
      include 'TOTAL'

      real dum(3)


      CHARACTER CB*3
      integer Iel,NFACES,npts
      parameter(npts=1000001)
      real*8 xnew,ynew
      real*8 xexact(npts),yexact(npts)
      common /nacap/ xexact,yexact


c      param(66) = 6.   ! These give the std nek binary i/o and are 
c      param(67) = 6.   ! good default values

      NFACES=2*NDIM

      ifxyo = .true.
c      call opcopy  (vx,vy,vz,xm1,ym1,zm1)
      call outpost (vx,vy,vz,pr,t,'gri')
      n = nx1*ny1*nz1*nelv

      open(unit=19,file='naca4412.dat')
      do i=1,npts
         read(19,*) xexact(i),yexact(i)
      enddo
      close (19)

      IFIELD = 1
c      if (nid.eq.0) then
c      if (istep.eq.2) then
      do  Iel=1,NELV            !do ieg=1,nelgt

!     Iel=GLLEL(ieg) 

c     Read GLL points (2D) from file
         do IFACE = 1,NFACES
            CB = CBC(IFACE,Iel,IFIELD)
            if (CB.EQ.'W  ') then
               CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1
     $              ,IFACE)
               do iz=KZ1,KZ2
                  do iy=KY1,KY2
                     do  ix=KX1,KX2
c     call fix_naca(0.4,0.04,0.12,xm1(ix,iy
c     $                       ,iz,Iel),ym1(ix,iy,iz,Iel),xnew,ynew)
                        call fix_naca_pts(xm1(ix,iy,iz,Iel),ym1(ix,iy
     $                        ,iz,Iel),xnew,ynew)
                        xm1(ix,iy,iz,Iel) = xnew
                        ym1(ix,iy,iz,Iel) = ynew


                     enddo
                  enddo
               enddo
            endif
         enddo

      enddo
c      endif
c      endif
!      call opcopy  (vx,vy,vz,xm1,ym1,zm1)
      call outpost (vx,vy,vz,pr,t,'gri ')
c      call exitt


      return
      end

!---------------------------------------------------------------------- 

      subroutine fix_naca_pts(xfoil,yfoil,xnew,ynew)

      real*8 xfoil,yfoil,xnew,ynew,dist_min
      integer npts,counter,i
      parameter(npts=1000001)

      real*8 dist(npts)

      real*8 xexact(npts),yexact(npts)
      common /nacap/ xexact,yexact

      dist_min = 1.
c     Computed the distance for each point

      do i=1,npts
         dist(i) = sqrt(abs(xexact(i)-xfoil)**2+abs(yexact(i)-yfoil)**2)
         if (dist(i).lt.dist_min) then
            counter = i
            dist_min = dist(i)
         endif
      enddo

      xnew = xexact(counter)
      ynew = yexact(counter)

      return
      end
c-----------------------------------------------------------------------

