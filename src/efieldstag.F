c
c
c  START OF EFIELDSTAG: Calculates efield on staggered C-grid
c
c
      subroutine efieldstag
     >  (nx,ny,nz,nor
     >  ,elec
     :  ,gxt,gyt,gzt
     :  ,nbw,nbe,nbs,nbn,nnz
     :  ,phi0
     :  ,iunit)  
     
      USE INDEX_MODULE
      USE GRID_MODULE
      USE COMMASMPI_MODULE
c
c
c 7.23.2002: replaced 1st order centered difference for Ez at k=1
c            with a 2nd order forward difference.
c
c 5.13.2002: update for vertical stretching
c
      implicit none
c
      integer nx,ny,nz,nor

      TYPE(VARIABLE)     :: elec(neelec)

!      real gxt(nx,4)
!      real gyt(ny,4)
!      real gzt(nz,4)

      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer    nbw,nbe,nbs,nbn,nnz
      real phi0(-nbw:nx+nbe,-nbs:ny+nbn,0:nnz-1)

      integer ix,jy,kz
      integer, parameter :: istag = 1, jstag = 1, kstag = 1
      
      integer iunit
      
      integer ibg  ! =1 for phi=0 at k=1, -1 for phi=0 at k=0
      
      integer imx,jmx,kmx, imx2, jmx2, kmx2
c
      real emag
      real emax, emax2
      real emin
c      real dx,dy,dz
c      real dxi,dyi,dzi
c      real dxi2,dyi2,dzi2
      
      integer im1,ip1
      integer jm1,jp1
      integer km1,kp1


      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: ixep, jyep, kzep

#ifdef MPI
      INCLUDE "mpif.h"

      integer, parameter :: ntot = 50
      real  mpitotin(ntot), mpitotout(ntot)
      integer n

#endif
c
c
c
c
c  compute electric field - three components
c
c
c   
c      dxi = 1./dx
c      dxi2 = .5*dxi
c      dyi = 1./dy
c      dyi2 = .5*dyi
c      dzi = 1./dz
c      dzi2 = .5*dzi
c
c  x-component
c
      if ( nx .gt. 2 ) then
c
      elec(iex)%flt3d(:,ny,:) = 0.0
      elec(iex)%flt3d(:,:,nz) = 0.0

c
      kzb = 1
      jyb = 1
      ixb = 1
#ifdef MPI
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixe = itile+1
      ixep = itile
      if (ixend .eq. nxend) ixep = ixend-ixbeg
#else

      kze = nz-1
      jye = ny-1
      ixe = nx
      ixep = nx-1
#endif

      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe
      
      im1 = ix-1
      ip1 = Min(ix,nx-1)
      
      elec(iex)%flt3d(ix,jy,kz) = -(phi0(ix,jy,kz) - phi0(im1,jy,kz))*gxt(ip1,4)
      
      ENDDO
      
      end do
      end do
c
c
      end if
c
c  y-component
c
      if ( ny .gt. 2 ) then
c
      elec(iey)%flt3d(nx,:,:) = 0.0
      elec(iey)%flt3d(:,:,nz) = 0.0

      kzb = 1
      jyb = 1
      ixb = 1
#ifdef MPI
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jye = jtile+1
      jyep = jtile
      if (jyend .eq. nyend) jyep = jyend-jybeg

      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
#else

      kze = nz-1
      jye = ny
      ixe = nx-1
      jyep = ny-1
#endif

      do kz = kzb,kze
       do jy = jyb,jye
!      do kz = 1, nz-1
!       do jy = 1, ny

        jm1 = jy-1
        jp1 = Min(jy,jyep)
        
        do ix = ixb,ixe
!        do ix = 1, nx-1
      
         elec(iey)%flt3d(ix,jy,kz) = -(phi0(ix,jy,kz) - phi0(ix,jm1,kz))*gyt(jp1,4)

      
        end do
       end do
      end do
c
c
      ELSE
       elec(iey)%flt3d(:,:,:) = 0.0
      end if
c
c  z-component
c
      if ( nz .gt. 2 ) then
c

      elec(iez)%flt3d(nx,:,:) = 0.0
      elec(iez)%flt3d(:,ny,:) = 0.0

      kzb = 1
      jyb = 1
      ixb = 1
#ifdef MPI
      kze = ktile
!      if (kzend .eq. nzend) kze = kzend-kzbeg

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
#else

      kze = nz
      jye = ny-1
      ixe = nx-1
      kzep = nz
#endif

      do kz = kzb,kze
!      do kz = 1, nz
        km1 = kz-1
        kp1 = kz ! Max(1,Min(kz,nz))
       do jy = jyb,jye
       do ix = ixb,ixe
!       do jy = 1, ny-1

!        do ix = 1, nx-1
      
         elec(iez)%flt3d(ix,jy,kz) = -(phi0(ix,jy,kz) - phi0(ix,jy,km1))*gzt(kp1,4)

      
        end do
       end do
      
!      print*,'ez : ',kz,elec(iez)%flt3d(nx/2,ny/2,kz),phi0(nx/2,ny/2,kz), phi0(nx/2,ny/2,kz-1),gzt(kp1,4)
      
      end do
c
c
!      stop
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

      kzb = 1
      jyb = 1
      ixb = 1
#ifdef MPI
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
#else

      kze = nz-1
      jye = ny-1
      ixe = nx-1
      kzep = nz
#endif
      do kz = kzb,kze
      do jy = jyb,jye
        jp1 = jy+1
      do ix = ixb,ixe
      
       ip1 = ix+1

c      do kz = 1,nz-kstag
c      do jy = 1,ny-jstag
c      do ix = 1,nx-istag 
        elec(iemag)%flt3d(ix,jy,kz) = 
     > Sqrt( 0.25*(( elec(iex)%flt3d(ix,jy,kz) + elec(iex)%flt3d(ip1,jy ,kz  ) )**2
     >     +       ( elec(iey)%flt3d(ix,jy,kz) + elec(iey)%flt3d(ix ,jp1,kz  ) )**2
     >     +       ( elec(iez)%flt3d(ix,jy,kz) + elec(iez)%flt3d(ix ,jy ,kz+1) )**2 )) 
     
       IF ( iemagstag > 1 ) THEN
        elec(iemagstag)%flt3d(ix,jy,kz) = 
     >  Sqrt( (( elec(iex)%flt3d(ix,jy,kz) )**2
     >     +  ( elec(iey)%flt3d(ix,jy,kz) )**2
     >     +  ( elec(iez)%flt3d(ix,jy,kz) )**2 )) 
       
       ENDIF
     
      end do
      end do
      end do
c
      emax  = 0.0
      emax2 = 0.0
      emin  = 0.0
        imx = 1
        jmx = 1
        kmx = 1
        imx2 = 1
        jmx2 = 1
        kmx2 = 1
      kzb = 1
      jyb = 1
      ixb = 1
#ifdef MPI
      kze = ktile
      if (kzend .eq. nzend) kze = kzend-kzbeg

      jye = jtile
      if (jyend .eq. nyend) jye = jyend-jybeg

      ixe = itile
      if (ixend .eq. nxend) ixe = ixend-ixbeg
#else

      kze = nz-1
      jye = ny-1
      ixe = nx-1
#endif
      do kz = kzb,kze
      do jy = jyb,jye
      do ix = ixb,ixe

c      do kz = 1,nz-kstag
c      do jy = 1,ny-jstag
c      do ix = 1,nx-istag
      IF ( emax .lt. elec(iemag)%flt3d(ix,jy,kz) ) THEN
        imx = ix
        jmx = jy
        kmx = kz
      ENDIF
      emax = max(emax,elec(iemag)%flt3d(ix,jy,kz))
      emin = min(emin,elec(iemag)%flt3d(ix,jy,kz))

       emag = 
     > Sqrt( (( elec(iex)%flt3d(ix,jy,kz) )**2
     >     +  ( elec(iey)%flt3d(ix,jy,kz) )**2
     >     +  ( elec(iez)%flt3d(ix,jy,kz) )**2 )) 

      IF ( emax2 .lt. emag ) THEN
        emax2 = emag
        imx2 = ix
        jmx2 = jy
        kmx2 = kz
      ENDIF

!        elec(iemag)%flt3d(ix,jy,kz) = Max( elec(iemag)%flt3d(ix,jy,kz), emag )

      end do
      end do
      end do
c

#ifdef MPI
! find global integrated rate max
       mpitotin(1)  =  emax
       mpitotin(2)  =  emax2
       n = 2

      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MAX, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
       emax = mpitotout(1)  
       emax2 = mpitotout(2)  
      ENDIF

! find global integrated rate max
       mpitotin(1)  = emin
       n = 1

      CALL MPI_Reduce(mpitotin, mpitotout, n, MPI_REAL, MPI_MIN, 0, my_comm, mpi_error_code)

       
      IF ( my_rank == 0 ) THEN
         emin = mpitotout(1)  
      ENDIF

#endif
      
      IF ( my_rank == 0 ) THEN
        write(iunit,*) 'EFIELD--E-MAX=',emax
        write(iunit,*) 'EFIELD--E-MIN=',emin
        write(iunit,*) 'ijk,Ez at emax: ',imx,jmx,kmx,elec(iez)%flt3d(imx,jmx,kmx),
     :               0.5*( elec(iez)%flt3d(imx,jmx,kmx) + elec(iez)%flt3d(imx,jmx,kmx+1) )
        write(iunit,*) 'EFIELD2--E-MAX=',emax2
        write(iunit,*) 'ijk,Ez at emax2: ',imx2,jmx2,kmx2,elec(iez)%flt3d(imx2,jmx2,kmx2)
      ENDIF
c
c
      return
      end
c
c
c  END OF EFIELD
c
c
