c
c
c  START OF EFIELD
c
c
      subroutine efield
     >  (nx,ny,nz
     >  ,id1,jd1,kd1,istag,jstag,kstag
     >  ,dx,dy,dz,elec
     :  ,gx,gy,gz
     :  ,iunit)  
     
      USE INDEX_MODULE
      USE GRID_MODULE
#ifdef MPI
      USE COMMASMPI_MODULE
#endif
c
c
c 7.23.2002: replaced 1st order centered difference for Ez at k=1
c            with a 2nd order forward difference.
c
c 5.13.2002: update for vertical stretching
c
      implicit none
c

      TYPE(VARIABLE)     :: elec(neelec)

      real gx(nx)
      real gy(ny)
      real gz(nz)

      integer ix,jy,kz
      integer nx,ny,nz
      integer id1,jd1,kd1,istag,jstag,kstag
      
      integer iunit
      
      integer ibg  ! =1 for phi=0 at k=1, -1 for phi=0 at k=0
      
      integer imx,jmx,kmx
c
      real emax
      real emin
      real dx,dy,dz
      real dxi,dyi,dzi
      real dxi2,dyi2,dzi2
      
      integer im1,ip1
      integer jm1,jp1
      integer km1,kp1
c
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 
#ifdef MPI
      integer ixb, jyb, kzb
      integer ixe, jye, kze
#endif
c
c  compute electric field - three components
c
c
c   
      dxi = 1./dx
      dxi2 = .5*dxi
      dyi = 1./dy
      dyi2 = .5*dyi
      dzi = 1./dz
      dzi2 = .5*dzi
c
c  x-component
c
      if ( nx .gt. 2 ) then
c
      elec(iex)%flt3d(nx,:,:) = 0.0
      elec(iex)%flt3d(:,ny,:) = 0.0
      elec(iex)%flt3d(:,:,nz) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 2
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1

      do kz = kzb,kze
       do jy = jyb,jye
        do ix = ixb,ixe
#else
      do kz = 1, nz-1
       do jy = 1, ny-1
        do ix = 2, nx-2
#endif
         im1 = Max(1,ix-1)
         ip1 = Min(nx-1,ix+1)

         elec(iex)%flt3d(ix,jy,kz) = 
     >       -(elec(ipot)%flt3d(ip1,jy,kz)-elec(ipot)%flt3d(im1,jy,kz))*0.5*gx(ix)

        enddo !! ix

#ifdef MPI
        if (ixbeg .eq. nxbeg) then
          elec(iex)%flt3d(1,jy,kz) = 
     >        -(elec(ipot)%flt3d(2,jy,kz)-elec(ipot)%flt3d(1,jy,kz))*gx(1)
        endif
         
        if (ixend .eq. nxend) then
          elec(iex)%flt3d(nx-2,jy,kz) = 
     >        -(elec(ipot)%flt3d(nx-1,jy,kz)-elec(ipot)%flt3d(nx-2,jy,kz))*gx(nx-2)

          elec(iex)%flt3d(nx-1,jy,kz) = 
     >        -(elec(ipot)%flt3d(nx-1,jy,kz)-elec(ipot)%flt3d(nx-3,jy,kz))*gx(nx-1)
        endif
#else
        do ix = 1, nx-1, nx-2

         im1 = Max(1,ix-1)
         ip1 = Min(nx-1,ix+1)
      
         elec(iex)%flt3d(ix,jy,kz) = 
     >       -(elec(ipot)%flt3d(ip1,jy,kz)-elec(ipot)%flt3d(im1,jy,kz))*gx(ix)
      
        end do  !! ix
#endif
       end do  !! jy
      end do  !! kz
c
c
      end if
c
c  y-component
c
      if ( ny .gt. 2 ) then
c
      elec(iey)%flt3d(nx,:,:) = 0.0
      elec(iey)%flt3d(:,ny,:) = 0.0
      elec(iey)%flt3d(:,:,nz) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 2
      if (jyend .eq. nyend) jye = jyend-jybeg-1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
       do jy = jyb,jye

        jm1 = Max(1,jy-1)
        jp1 = Min(ny-1,jy+1)

        do ix = ixb,ixe
#else
      do kz = 1, nz-1
       do jy = 2, ny-2

        jm1 = Max(1,jy-1)
        jp1 = Min(ny-1,jy+1)

        do ix = 1, nx-1
#endif
         elec(iey)%flt3d(ix,jy,kz) = 
     >    -(elec(ipot)%flt3d(ix,jp1,kz)-elec(ipot)%flt3d(ix,jm1,kz))*0.5*gy(jy)

        end do !! ix
       end do !! jy

#ifdef MPI
       ixb = 1
       ixe = itile
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do ix = 1, nx-1
        if (jybeg .eq. nybeg) then

         elec(iey)%flt3d(ix,1,kz) = 
     >       -(elec(ipot)%flt3d(ix,2,kz)-elec(ipot)%flt3d(ix,1,kz))*gy(1)

        endif
        
        if (jyend .eq. nyend) then

         elec(iey)%flt3d(ix,ny-2,kz) = 
     >       -(elec(ipot)%flt3d(ix,ny-1,kz)-elec(ipot)%flt3d(ix,ny-3,kz))*gy(ny-2)

         elec(iey)%flt3d(ix,ny-1,kz) = 
     >       -(elec(ipot)%flt3d(ix,ny-1,kz)-elec(ipot)%flt3d(ix,ny-2,kz))*gy(ny-1)

        endif
        
       enddo !!ix
#else
       do jy = 1, ny-1, ny-2

        jm1 = Max(1,jy-1)
        jp1 = Min(ny-1,jy+1)

        do ix = 1, nx-1

         elec(iey)%flt3d(ix,jy,kz) = 
     >       -(elec(ipot)%flt3d(ix,jp1,kz)-elec(ipot)%flt3d(ix,jm1,kz))*gy(jy)
      
        end do !! ix
       end do !! jy
#endif
      end do !! kz
c
c
      end if
c
c  z-component
c
      if ( nz .gt. 2 ) then
c

      elec(iez)%flt3d(nx,:,:) = 0.0
      elec(iez)%flt3d(:,ny,:) = 0.0
      elec(iez)%flt3d(:,:,nz) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg-1

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do jy = jyb,jye
       do kz = kzb,kze

        km1 = Max(1,kz-1)
        kp1 = Min(nz-1,kz+1)

        do ix = ixb,ixe
#else
      do jy = 1, ny-1
       do kz = 2, nz-2

        km1 = Max(1,kz-1)
        kp1 = Min(nz-1,kz+1)

        do ix = 1, nx-1
#endif
         elec(iez)%flt3d(ix,jy,kz) = 
     >    -(elec(ipot)%flt3d(ix,jy,kp1)-elec(ipot)%flt3d(ix,jy,km1))*0.5*gz(kz)
      
        end do !! ix
       end do !! kz

!       do kz = 1, nz-1, nz-1

        km1 = Max(1,kz-1)
        kp1 = Min(nz,kz+1)

#ifdef MPI
        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do ix = ixb,ixe
#else
        do ix = 1, nx-1
#endif
         elec(iez)%flt3d(ix,jy,1) = 
     >   -(elec(ipot)%flt3d(ix,jy,1))*gz(1)


         elec(iez)%flt3d(ix,jy,nz-1) = 
     >   -(elec(ipot)%flt3d(ix,jy,nz)-elec(ipot)%flt3d(ix,jy,nz-2))*0.5*gz(nz-1)
      
        end do !! ix
!       end do

      end do !! jy
c
c
      endif



c
c 1st order difference
c      elec()%flt3d(ix,jy,01,iez) = -elec()%flt3d(ix,jy,01,ipot)*dzi*
c     :   gt(ix,jy,01,imapz)

c 2nd order difference
c      elec()%flt3d(ix,jy,01,iez) = -elec()%flt3d(ix,jy,02,ipot)*dzi2*
c     :   gt(ix,jy,01,imapz)


c 2nd order difference -- not working right?
c      elec()%flt3d(ix,jy,01,iez) = 
c     :    -0.5*gt(ix,jy,2,imapz)/dz* 
c     :          (  -3.0*elec()%flt3d(ix,jy,1,ipot)  
c     :              + 4*elec()%flt3d(ix,jy,2,ipot) 
c     :              -   elec()%flt3d(ix,jy,3,ipot)   )

c      elec()%flt3d(ix,jy,nz-kstag,iez) = 
c     > -(elec()%flt3d(ix,jy,nz-kstag,ipot)-elec()%flt3d(ix,jy,nz-2*kstag,ipot))*dzi*
c     :   gt(ix,jy,nz-kstag,imapz)

c
c  e-magnitude
c
#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
       do jy = jyb,jye
        do ix = ixb,ixe
#else
      do kz = 1,nz-kstag
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
#endif
        elec(iemag)%flt3d(ix,jy,kz) = 
     >      Sqrt( elec(iex)%flt3d(ix,jy,kz)**2
     >          + elec(iey)%flt3d(ix,jy,kz)**2
     >          + elec(iez)%flt3d(ix,jy,kz)**2 ) 
        end do
       end do
      end do
c
      emax = 0.0
      emin = 0.0
      imx = 1
      jmx = 1
      kmx = 1
#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg+1-jstag

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg+1-istag

      do kz = kzb,kze
       do jy = jyb,jye
        do ix = ixb,ixe
#else
      do kz = 1,nz-kstag
       do jy = 1,ny-jstag
        do ix = 1,nx-istag
#endif
         IF ( emax .lt. elec(iemag)%flt3d(ix,jy,kz) ) THEN
          imx = ix
          jmx = jy
          kmx = kz
         ENDIF
      
         emax = max(emax,elec(iemag)%flt3d(ix,jy,kz))
         emin = min(emin,elec(iemag)%flt3d(ix,jy,kz))
      
        end do
       end do
      end do

#ifdef MPI
      write(0,*) 'EFIELD: E-MAX=',emax,'for my_rank=',my_rank
      write(0,*) 'EFIELD: E-MIN=',emin,'for my_rank=',my_rank
      write(0,*) 'EFIELD: ijk,Ez at emax: ',imx,jmx,kmx,elec(iez)%flt3d(imx,jmx,kmx),'for my_rank=',my_rank
#else
      write(iunit,*) 'EFIELD--E-MAX=',emax
      write(iunit,*) 'EFIELD--E-MIN=',emin
      write(iunit,*) 'ijk,Ez at emax: ',imx,jmx,kmx,elec(iez)%flt3d(imx,jmx,kmx)
#endif
c
      return
      end
c
c
c  END OF EFIELD
c
c

c
c
c  START OF EFIELD2 : For filling EX2, EY2, EZ2 with unstaggered values (when the arrays are turned on)
c
c
      subroutine efield2
     >  (nx,ny,nz
     >  ,id1,jd1,kd1,istag,jstag,kstag
     >  ,dx,dy,dz,elec
     :  ,gx,gy,gz
     :  ,iunit)  
     
      USE INDEX_MODULE
      USE GRID_MODULE
#ifdef MPI
      USE COMMASMPI_MODULE
#endif
c
c
c 7.23.2002: replaced 1st order centered difference for Ez at k=1
c            with a 2nd order forward difference.
c
c 5.13.2002: update for vertical stretching
c
      implicit none
c

      TYPE(VARIABLE)     :: elec(neelec)

      real gx(nx)
      real gy(ny)
      real gz(nz)

      integer ix,jy,kz
      integer nx,ny,nz
      integer id1,jd1,kd1,istag,jstag,kstag
      
      integer iunit
      
      integer ibg  ! =1 for phi=0 at k=1, -1 for phi=0 at k=0
      
      integer imx,jmx,kmx
c
      real emax
      real emin
      real dx,dy,dz
      real dxi,dyi,dzi
      real dxi2,dyi2,dzi2
      
      integer im1,ip1
      integer jm1,jp1
      integer km1,kp1
c
!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 
#ifdef MPI
      integer ixb, jyb, kzb
      integer ixe, jye, kze
#endif
c
c  compute electric field - three components
c
c
c   
      dxi = 1./dx
      dxi2 = .5*dxi
      dyi = 1./dy
      dyi2 = .5*dyi
      dzi = 1./dz
      dzi2 = .5*dzi
c
c  x-component
c
      if ( nx .gt. 2 ) then
c
      elec(iex2)%flt3d(nx,:,:) = 0.0
      elec(iex2)%flt3d(:,ny,:) = 0.0
      elec(iex2)%flt3d(:,:,nz) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixbeg .eq. nxbeg) ixb = 2
      if (ixend .eq. nxend) ixe = ixend-ixbeg-1

      do kz = kzb,kze
       do jy = jyb,jye
        do ix = ixb,ixe
#else
      do kz = 1, nz-1
       do jy = 1, ny-1
        do ix = 2, nx-2
#endif
         im1 = Max(1,ix-1)
         ip1 = Min(nx-1,ix+1)

         elec(iex2)%flt3d(ix,jy,kz) = 
     >       -(elec(ipot)%flt3d(ip1,jy,kz)-elec(ipot)%flt3d(im1,jy,kz))*0.5*gx(ix)

        enddo !! ix

#ifdef MPI
        if (ixbeg .eq. nxbeg) then
          elec(iex2)%flt3d(1,jy,kz) = 
     >        -(elec(ipot)%flt3d(2,jy,kz)-elec(ipot)%flt3d(1,jy,kz))*gx(1)
        endif
         
        if (ixend .eq. nxend) then
          elec(iex2)%flt3d(nx-2,jy,kz) = 
     >        -(elec(ipot)%flt3d(nx-1,jy,kz)-elec(ipot)%flt3d(nx-2,jy,kz))*gx(nx-2)

          elec(iex2)%flt3d(nx-1,jy,kz) = 
     >        -(elec(ipot)%flt3d(nx-1,jy,kz)-elec(ipot)%flt3d(nx-3,jy,kz))*gx(nx-1)
        endif
#else
        do ix = 1, nx-1, nx-2

         im1 = Max(1,ix-1)
         ip1 = Min(nx-1,ix+1)
      
         elec(iex2)%flt3d(ix,jy,kz) = 
     >       -(elec(ipot)%flt3d(ip1,jy,kz)-elec(ipot)%flt3d(im1,jy,kz))*gx(ix)
      
        end do  !! ix
#endif
       end do  !! jy
      end do  !! kz
c
c
      end if
c
c  y-component
c
      if ( ny .gt. 2 ) then
c
      elec(iey2)%flt3d(nx,:,:) = 0.0
      elec(iey2)%flt3d(:,ny,:) = 0.0
      elec(iey2)%flt3d(:,:,nz) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jyb = 1
      jye = jtile
      if (jybeg .eq. nybeg) jyb = 2
      if (jyend .eq. nyend) jye = jyend-jybeg-1

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do kz = kzb,kze
       do jy = jyb,jye

        jm1 = Max(1,jy-1)
        jp1 = Min(ny-1,jy+1)

        do ix = ixb,ixe
#else
      do kz = 1, nz-1
       do jy = 2, ny-2

        jm1 = Max(1,jy-1)
        jp1 = Min(ny-1,jy+1)

        do ix = 1, nx-1
#endif
         elec(iey2)%flt3d(ix,jy,kz) = 
     >    -(elec(ipot)%flt3d(ix,jp1,kz)-elec(ipot)%flt3d(ix,jm1,kz))*0.5*gy(jy)

        end do !! ix
       end do !! jy

#ifdef MPI
       ixb = 1
       ixe = itile
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do ix = 1, nx-1
        if (jybeg .eq. nybeg) then

         elec(iey2)%flt3d(ix,1,kz) = 
     >       -(elec(ipot)%flt3d(ix,2,kz)-elec(ipot)%flt3d(ix,1,kz))*gy(1)

        endif
        
        if (jyend .eq. nyend) then

         elec(iey2)%flt3d(ix,ny-2,kz) = 
     >       -(elec(ipot)%flt3d(ix,ny-1,kz)-elec(ipot)%flt3d(ix,ny-3,kz))*gy(ny-2)

         elec(iey2)%flt3d(ix,ny-1,kz) = 
     >       -(elec(ipot)%flt3d(ix,ny-1,kz)-elec(ipot)%flt3d(ix,ny-2,kz))*gy(ny-1)

        endif
        
       enddo !!ix
#else
       do jy = 1, ny-1, ny-2

        jm1 = Max(1,jy-1)
        jp1 = Min(ny-1,jy+1)

        do ix = 1, nx-1

         elec(iey2)%flt3d(ix,jy,kz) = 
     >       -(elec(ipot)%flt3d(ix,jp1,kz)-elec(ipot)%flt3d(ix,jm1,kz))*gy(jy)
      
        end do !! ix
       end do !! jy
#endif
      end do !! kz
c
c
      end if
c
c  z-component
c
      if ( nz .gt. 2 ) then
c

      elec(iez2)%flt3d(nx,:,:) = 0.0
      elec(iez2)%flt3d(:,ny,:) = 0.0
      elec(iez2)%flt3d(:,:,nz) = 0.0

#ifdef MPI
      kzb = 1
      kze = ktile
      if (kzbeg .eq. nzbeg) kzb = 2
      if (kzend .eq. nzend) kze = kzend-kzbeg-1

      jyb = 1
      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixb = 1
      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg

      do jy = jyb,jye
       do kz = kzb,kze

        km1 = Max(1,kz-1)
        kp1 = Min(nz-1,kz+1)

        do ix = ixb,ixe
#else
      do jy = 1, ny-1
       do kz = 2, nz-2

        km1 = Max(1,kz-1)
        kp1 = Min(nz-1,kz+1)

        do ix = 1, nx-1
#endif
         elec(iez2)%flt3d(ix,jy,kz) = 
     >    -(elec(ipot)%flt3d(ix,jy,kp1)-elec(ipot)%flt3d(ix,jy,km1))*0.5*gz(kz)
      
        end do !! ix
       end do !! kz

!       do kz = 1, nz-1, nz-1

        km1 = Max(1,kz-1)
        kp1 = Min(nz,kz+1)

#ifdef MPI
        ixb = 1
        ixe = itile
        if (ixend .eq. nxend) ixe = ixend-ixbeg

        do ix = ixb,ixe
#else
        do ix = 1, nx-1
#endif
         elec(iez2)%flt3d(ix,jy,1) = 
     >   -(elec(ipot)%flt3d(ix,jy,1))*gz(1)


         elec(iez2)%flt3d(ix,jy,nz-1) = 
     >   -(elec(ipot)%flt3d(ix,jy,nz)-elec(ipot)%flt3d(ix,jy,nz-2))*0.5*gz(nz-1)
      
        end do !! ix
!       end do

      end do !! jy
c
c
      endif



c
      return
      end
c
c
c  END OF EFIELD2
c
c
