c
c     mudtst5: test code for vertical stretching and boundary below ground
c
c      sets up a simple dipole with two spheres of charge and
c       solves for the potential everywhere and vertical electric field
c       along the symmetry axis
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 5.0                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose
c
c     test program for the MUDPACK solver mud3sp
c
c ... required MUDPACK files
c
c     mud3sp.f, mudcom.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud3sp
c
c **********************************************************
c **********************************************************
c
c     a sample program/test driver for mud3sp is below. it can be
c     executed as an initial test.  the output is listed for the
c     test case described.
c
c     test the driver below by solving the separable elliptic pde
c
c        d(x**2*dp/dx)/dx + d(y**2*dp/dy)/dy + d(z**2*dp/dz)/dz
c
c        - (x+y+z)*pe(x,y,z) = r(x,y,z)
c
c
c     on the cube [0.5,1.0]**3 with specified boundary conditions at x,y,z=1.0
c     and derivative boundary conditions of the form
c
c        dp/dx - pe(0.5,y,z) = f(0.5,y,z) at x = 0.5
c
c        dp/dy - pe(x,0.5,z) = g(x,0.5,z) at y = 0.5
c
c        dp/dz - pe(x,y,0.5) = h(x,y,0.5) at z = 0.5
c
c     generate the approximation on a  65 x 49 x 41 grid using the exact
c     solution
c
c          pe(x,y,z) = (x*y*z)**4
c
c     for testing.  notice the solver mud3 would require far more work space
c     than mud3sp if used to solve this pde.
c
c ******************************************************
c     output (32 bit floating point arithmetic)
c *******************************************************
c
c     mud3sp test
c
c     input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
c     nze =  2 nzf =  1
c     ixp =  2 jyq =  3 kzr =  5
c     iex =  6 jey =  5 kez =  4
c     nx =  65 ny =  49 nz =  41 iguess =  0 maxcy =  3
c     method =  0 work space length input =  514258
c     xa =  0.50 xb =  1.00
c     yc =  0.50 yd =  1.00
c     ze =  0.50 zf =  1.00
c     tolmax =  0.000E+00
c
c     multigrid options
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  1
c
c     discretization call to mud3sp intl =  0
c     ierror =  0 minimum work space =  472557
c
c     approximation call to mud3sp
c     intl =  1 method =  0 iguess =  0 maxcy =  3
c     ierror =  0
c     maximum error  =   0.719E-05
c
c ************************************************************
c     end of output
c ************************************************************
c
c 7/10/2006: (ERM) Removed some items from call list that are now in elec_module.
c                  Pass in just the potential from the elec array as 'pot'
c
      subroutine mud(nx,ny,nz,nor,
     :   dx,dy,dz,gxt,gyt,gzt,
     :   nnx,nny,nnz,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,
     :   llwork, istretch,
     :   t0,ibg,
     :   elec,iestag,iunit,lcalce,bcx,bcy)
      
      USE GRID_MODULE
      USE ELEC_MODULE
      USE INDEX_MODULE
      
      implicit none

      integer nx,ny,nz,nor

      TYPE(VARIABLE)     :: elec(neelec)
      
!      real               :: gxt(nx,4)
!      real               :: gyt(ny,4)
!      real               :: gzt(nz,4)
      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer iixp,jjyq,kkzr,iiex,jjey,kkez
      integer*8 llwork
      integer nnx,nny,nnz
      integer nbw,nbe,nbs,nbn
      real dx,dy,dz
      integer ibg
      logical istretch
      integer iunit
      integer iestag
      logical lcalce   ! whether to calculate electric field components when iestag=1
      
c input arrays: t0 holds net charge
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
      
      integer bcx, bcy
c      real pot(nx,ny,nz)
      
c   local variables/arrays
      
      real, allocatable ::  pot2(:,:,:) 
      real, allocatable ::  rhs(:,:,:)
      real, allocatable ::  work(:)
c      real pot2(nnx,nny,nnz) 
c      real rhs(nnx,nny,nnz)
c      real work(llwork)
      integer ix,jy,kz
      
c ####
c   begin execute
c ####      
       
      IF ( ny .gt. 2 ) THEN
      
       allocate (pot2(nnx,nny,nnz))
       allocate (rhs(nnx,nny,nnz))
       allocate (work(llwork))
      
       call mud1(nx,ny,nz,nor,
     :   dx,dy,dz,
     :   nnx,nny,nnz,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,
     :   llwork,istretch,
     :   t0,ibg,ibdcnd,jbdcnd,kbdcnd,kdbcdn,dlz,
     :   elec(ipot)%flt3d,potfair,rhofair,z1d2,nzmax1,
     :   pot2,pot2,rhs,rhs,work,iunit, bcx, bcy)
      
       IF ( iestag .eq. 1 .and. lcalce ) THEN

        call efieldstag
     >  (nx,ny,nz,nor
     >  ,elec
     :  ,gxt,gyt,gzt
     :  ,nbw,nbe,nbs,nbn,nnz
     :  ,pot2
     :  ,iunit)  
       
       ENDIF
       
       deallocate (pot2)
       deallocate (rhs)
       deallocate (work)
       
      ELSE
c      
!       write(0,*) 'call mud2d'
       call mud2d(nx,ny,nz,nor,
     :   dx,dy,dz,gxt,gyt,gzt,
     :   nnx,nny,nnz,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,
     :   llwork, istretch,
     :   t0,ibg,
     :   elec,iestag,iunit,lcalce, bcx, bcy)
      
!      write(iunit,*) 'done with mud2d'
      ENDIF
      
      RETURN
      END
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mud1(nx,ny,nz,nor,
     :   dx,dy,dz,
     :   nnx,nny,nnz,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,
     :   llwork,istretch,
     :   t0,ibg,ibdcnd,jbdcnd,kbdcnd,kdbcdn,dlz,
     :   pot,potfair,rhofair,z1d2,nzmax1,
     :   phi,phi0,rhs,rhs0,work,iunit, bcx, bcy)
c
c
      implicit none
      
      integer ndebug
      parameter ( ndebug = 0 )

      integer nx,ny,nz,nor
      integer iunit
      integer iixp,jjyq,kkzr,iiex,jjey,kkez
      integer*8 llwork
      integer nnx,nny,nnz
      integer nbw,nbe,nbs,nbn
      real dx,dy,dz
      integer ibg,ibdcnd,jbdcnd,kbdcnd,kdbcdn

      logical istretch

      real phi(nnx,nny,nnz),rhs(nnx,nny,nnz),work(llwork)
      real phi0(-nbw:nx+nbe,-nbs:ny+nbn,0:nnz-1)
      real rhs0(-nbw:nx+nbe,-nbs:ny+nbn,0:nnz-1)
      
c input arrays: t0 holds net charge
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      real pot(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)

      real potfair(nnz),rhofair(nnz)
      
      integer bcx, bcy
      
      integer ix,jy,kz
c
c     set grid sizes with parameter statements
c
c      integer iixp,jjyq,kkzr,iiex,jjey,kkez,llwork,nnx,nny,nnz
c      parameter(iixp=3,jjyq=3,kkzr=2)
c      parameter(iiex=5,jjey=5,kkez=6)
c      parameter (nnx = iixp*2**(iiex-1)+1)
c      parameter (nny = jjyq*2**(jjey-1)+1)
c      parameter (nnz = kkzr*2**(kkez-1)+1)
      
c      integer nxr,nyr,nzr
c      parameter(nxr=nnx,nyr=nny,nzr=41)
c
c     set work space length approximation (see mud3sp.d)
c
c      parameter(llwork = 7*(nnx+2)*(nny+2)*(nnz+2)/2 )
c
c
c     dimension solution,right hand side, and work arrays
c
c
c pretend that the 'real nz' is nnz-1, so that phi0 and rhs0
c  will correspond to the 'regular' grid (add buffer zones next)
c
      integer*8 iprm(22)
      integer mgopt(4)
      real fprm(8)

      real eperao
      parameter (eperao  = 8.8592e-12 )
      
c      real rho(nnx,nny,nzr),bn,bp
      
c      real z1d(nzr,4)
      real z0
      real dzbottom, dzmax, zbot,zstr
      integer nbndlyr
      real zfind
      integer nzmax,nzmax1,nz1d
      real dlz
!      parameter(nzmax1=400)
      real z1d2(nzmax1,4)
      real dlz1,dlznew
c
c timing vars
c
      real tim1,tim2,timarr(2)
      real second,commas_dtime

c
c     put integer and floating point arguments names in contiguous
c     storeage for labelling
c
      integer*8 intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
! NOTE: need these dummy common blocks to align the variable names
! in memory so that the equivalence with iprm etc. will work correctly.
      common/itmsp3/intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz2,x,y,z,cxx,cyy,czz,cx,cy,cz,cex,cey,cez,ce,errm
      real pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k,ierror
      equivalence(intl,iprm)
      equivalence(xa,fprm)
c      equivalence(phi(1,1,2),phi0(1,1,1))
c      equivalence(rhs(1,1,2),rhs0(1,1,1))
c
c     declare coefficient and boundary condition input subroutines external
c
      external cfx,cfy,cfz,cfzu,bndc
!      common/strz/z1d2,dlz,nz1d
c
c ************************************************************
!      nzmax = nzmax1
c
c     set input integer parameters
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 2
      nxb1 = 2
      IF ( bcx .eq. 2 ) THEN
!       nxa = 0
!       nxb1 = 0
      ENDIF
      nyc = 2
      nyd = 2
      IF ( bcy .eq. 2 ) THEN
!       nyc = 0
!       nyd = 0
      ENDIF
      nze = 1 
      nzf = 1  
      IF ( kbdcnd .eq. 2 ) THEN
      nze = 1 
      nzf = 1
      ELSEIF ( kbdcnd .eq. 1 ) THEN
      nze = 1
      nzf = 1  ! changed from 2 to 1 on 12/31/2004 by erm 
               ! solution seems to be more stable as a result, especially
               ! when there is an upper level grid stretch.
               ! Also seems to preserve net charge better than nzf = 2
      ELSEIF ( kbdcnd .eq. 3 ) THEN
      nze = 2 
      nzf = 2  
      ENDIF
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      kzr = kkzr
      iex = iiex
      jey = jjey
      kez = kkez
      nxm = nnx
      nym = nny
      nzm = nnz
c      nz1d = nnz
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set method of relaxation--point gauss/seidel only for mud3sp
c
      method = 0
c
c     set default multigrid options
c
      mgopt(1) = 2  ! kcycle (default 2)
      mgopt(2) = 1  ! iprer  (default 2)
      mgopt(3) = 1  ! ipost  (default 1)
      mgopt(4) = 3  ! intpol (default 3)
c
c     set full multigrid cycling by flagging no initial guess
c
      iguess = 0
c
c     set a limit of three w(2,1) cycles from the finest level
c
      maxcy = 2
c
c     set mesh increments
c
c      dlx = (xb-xa)/float(nxm-1)
c      dly = (yd-yc)/float(ny-1)
c      dlz = (zf-ze)/float(nz-1)
c      dlx = 500.
c      dly = 500.
c      dlz = 500.
      
c      dlx = 250.
c      dly = 250.
c      dlz = 250.
c
c     set end points of solution cube in (x,y,z) space
c
      xa = 0.
      xb = (nxm-1)*dx
      yc = 0.
      yd = (nym-1)*dy
      ze = -dlz*0.5
      zf = ze + (nzm-1)*dlz
      
      zbot = 0.0

c      STOP
c
c     set for no error control (this allows phi to be equivalenced with work)
c
      tolmax = 0.0
c
c     set right hand side in rhs 
c     and initialize phi to zero
c
c      real phi0(-nbw:nx+nbe,-nbs:ny+nbn,0:nnz-1)
c      real rhs0(-nbw:nx+nbe,-nbs:ny+nbn,0:nnz-1)

      do kz = 0, nnz-1
      do jy = -nbs, ny+nbn
      do ix = -nbw, nx+nbe
      
      IF ( kz .ge. 1 .and. kz .le. nz-1 .and. 
     :     jy .ge. 1 .and. jy .le. ny-1 .and.
     :     ix .ge. 1 .and. ix .le. nx-1 ) THEN
        rhs0(ix,jy,kz) = -t0(ix,jy,kz)/eperao
        phi0(ix,jy,kz) = pot(ix,jy,kz)
      ELSEIF ( kz .eq. 0 ) THEN
        rhs0(ix,jy,kz) = 0.0
        phi0(ix,jy,kz) = 0.0
      ELSE
        rhs0(ix,jy,kz) = -rhofair(kz)/eperao
        phi0(ix,jy,kz) = potfair(kz)
      ENDIF

      end do
      end do
      end do
      
      IF ( bcx .eq. 2 ) rhs0(nx,1:ny-1,1:nz-1) = rhs0(1,1:ny-1,1:nz-1)
      IF ( bcy .eq. 2 ) rhs0(1:nx-1,ny,1:nz-1) = rhs0(1:nx-1,1,1:nz-1)

c      do k=1,nz-1
c      do j=1,nym
c        do i=1,nxm
c
c        end do
c      end do
c      end do

c      DO k=nz,nnz
c      do j=1,nym
c      do i=1,nxm
c        phi(i,j,k) = 0.0
c        rhs(i,j,k) = 0.0
c      end do
c      end do
c      ENDDO
c
c     set specified values at lower boundary
c
c      do j=1,nym
c      do i=1,nxm
c        phi(i,j,1) = 0.0
c        rhs(i,j,1) = 0.0
c      end do
c      end do
      
      IF ( ndebug .ge. 1 ) THEN
      write(6,50)
   50 format(//' mud3sp debug test ')
c
c     print input arguments
c
      write(6,100)intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +          nxm,nym,nzm,iguess,maxcy,method,nwork,xa,xb,yc,yd,ze,zf,
     +          tolmax,mgopt(1),mgopt(2),mgopt(3),mgopt(4)
  100 format(/' input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' nze = ',i2, ' nzf = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' kzr = ',i2,
     +/' iex = ',i2, ' jey = ',i2, ' kez = ',i2,
     +/' nxm = ',i3,' nym = ',i3,' nzm = ',i3, ' iguess = ',i2,
     +/' maxcy = ',i2,
     +/' method = ',i2, ' work space length input = ',i7,
     +/' xa = ',f7.2,' xb = ',f12.2,
     +/' yc = ',f7.2,' yd = ',f12.2,
     +/' ze = ',f7.2,' zf = ',f12.2,
     +/' tolmax = ' ,e10.3
     +//' multigrid options '
     +/' kcycle = ',i2
     +/' iprer = ',i2
     +/' ipost = ',i2
     +/' intpol = ',i2 )
     
      write(6,*) '(zf-ze)/float(nzm-1) = ',(zf-ze)/float(nzm-1)
      
      ENDIF ! (ndebug)
      tim1 = commas_dtime(timarr) 
c
c     discretize pde
c
      IF (ndebug .ge. 1) THEN
        write(*,104) intl
      ENDIF
  104 format(/' discretization call to mud3sp', ' intl = ', i2)
      IF ( istretch ) THEN
      call mud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
      ELSE
      call mud3sp(iprm,fprm,work,cfx,cfy,cfzu,bndc,rhs,phi,mgopt,ierror)
      ENDIF
      IF (ndebug .ge. 1) THEN
        write (*,105) ierror,iprm(21)
      ENDIF

  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      
      IF (ierror.gt.0) THEN
       call exit(0)
      ENDIF
c
c     approximate pde
c
      intl = 1
      IF (ndebug .ge. 1) THEN
        write(*,106) intl,method,iguess,maxcy
      ENDIF
  106 format(/' approximation call to mud3sp ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      IF ( istretch ) THEN
      call mud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
      ELSE
      call mud3sp(iprm,fprm,work,cfx,cfy,cfzu,bndc,rhs,phi,mgopt,ierror)
      ENDIF

      IF (ndebug .ge. 1) THEN
        write (*,107) ierror
      ENDIF
  107 format(' ierror = ',i2)

      if (ierror.gt.0) call exit(0)

      tim2 = 0.0
      tim2 = commas_dtime(timarr)
      write(iunit,*) 'time1 for mud = ',tim2,timarr(1)

c      tim1 = commas_dtime(timarr) 
c      intl = 1
c      iguess = 1
c      write(*,106) intl,method,iguess,maxcy
c      call mud3sp(iprm,fprm,work,cfx,cfy,cfz,bndc,rhs,phi,mgopt,ierror)
c      write (*,107) ierror
c      if (ierror.gt.0) call exit(0)
c
c      tim2 = 0.0
c      tim2 = commas_dtime(timarr)
c      write(6,*) 'time2 for mud = ',tim2,timarr(1)
c
c     compute and print maximum error
c
      IF ( ndebug .ge. 1 ) THEN
      write(6,*) 'k,z,phi,rho through nx/2,ny/2'
      DO k=nz-1,2,-1
c        cz = 1./z1d2(k,3)*(1./z1d2(k+1,4)-1./z1d2(k,4) )/dlz
c        cz = (dlz*(k-1)+dlz/2.)/z1d(k,1)
        write(6,*) k,z1d2(k,1),phi0(nx/2,ny/2,k),
     :   -0.5/dlz*(phi0(nx/2,ny/2,k+1)-phi0(nx/2,ny/2,k-1))*
     :     z1d2(k,3),
     :   t0(nx/2,ny/2,k)
c     :   ,cz
      ENDDO
      k=1
        write(6,*) k,z1d2(k,1),phi0(nxm/2,nym/2,k),
     :   -0.5/dlz*(phi0(nx/2,ny/2,k+1)-0.0)*z1d2(k,3),
     :   t0(nx/2,ny/2,k)

      k=1
        write(6,*) k-1,-z1d2(k,1),phi(nxm/2,nym/2,k)
      
      ENDIF ! (ndebug .ge. 1)
      
c      write(*,108) errm
c  108 format(' maximum error  =  ',e10.3)
c      end
c
c  put pot into elec
c
      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx
         pot(ix,jy,kz) = phi0(ix,jy,kz)
      end do
      end do
      end do


      
      RETURN
      END
      
c #####################################################################
c
c     coefficient subroutines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfx(x,cxx,cx,cex)
      implicit none
      real x,cxx,cx,cex
      cxx = 1.0
      cx = 0.0
      cex = 0.0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfy(y,cyy,cy,cey)
      implicit none
      real y,cyy,cy,cey
      cyy = 1.0
      cy = 0.0
      cey = 0.0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfzu(z,czz,cz,cez)
      implicit none
      real z,czz,cz,cez
      czz = 1.0
      cz = 0.0
      cez = 0.0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfz(z,czz,cz,cez)
      USE ELEC_MODULE
      implicit none
      integer kz
      real z,czz,cz,cez
!      real dlz
!      integer nzmax,nz1d
!      parameter (nzmax=400)
!      real z1d2(nzmax,4)
!      common/strz/z1d2,dlz,nz1d

c      kz = Int(z/dlz) 
      kz = Max(0,NInt(z/dlz - 1.5)) + 1
c      write(6,*) 'z,kz=',z,kz
      IF ( kz .gt. 0 ) THEN
        czz = z1d2(Min(kz,nz1d),3)**(2)
      ELSE
        czz = z1d2(1,3)**(2)
      ENDIF
      IF ( kz .lt. nz1d .and. kz .gt. 1 ) THEN
         cz = z1d2(kz,3)*(z1d2(kz+1,3)-z1d2(kz-1,3))/dlz*0.5
      ELSEIF ( kz .le. 1 ) THEN
         cz = z1d2(1,3)*(z1d2(2,3)-z1d2(1,3))/dlz
      ELSE ! ( kz .ge. nz1d ) THEN
         cz = z1d2(nz1d,3)*(z1d2(nz1d,3)-z1d2(nz1d-1,3))/dlz
      ENDIF

c      write(6,*) 'z,kz,czz,cz=',z,kz,czz,cz
c      czz = 1.0
c      cz = 0.0
      
      cez = 0.0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cfzlgt(z,czz,cz,cez)
      USE ELEC_MODULE
      implicit none
      integer kz
      real z,czz,cz,cez
!      real dlz
!      integer nzmax,nz1d
!      parameter (nzmax=400)
!      real z1d2(nzmax,4)
!      common/strz/z1d2,dlz,nz1d

c      kz = Int(z/dlz) 
      kz = Max(0,NInt(z/dlzlgt - 1.5)) + 1
c      write(6,*) 'z,kz=',z,kz
      IF ( kz .gt. 0 ) THEN
        czz = z1d2lgt(Min(kz,nz1dlgt),3)**(2)
      ELSE
        czz = z1d2lgt(1,3)**(2)
      ENDIF
      IF ( kz .lt. nz1dlgt .and. kz .gt. 1 ) THEN
         cz = z1d2lgt(kz,3)*(z1d2lgt(kz+1,3)-z1d2lgt(kz-1,3))/dlzlgt*0.5
      ELSEIF ( kz .le. 1 ) THEN
         cz = z1d2lgt(1,3)*(z1d2lgt(2,3)-z1d2lgt(1,3))/dlzlgt
      ELSE ! ( kz .ge. nz1d ) THEN
         cz = z1d2lgt(nz1dlgt,3)*(z1d2lgt(nz1dlgt,3)-z1d2lgt(nz1dlgt-1,3))/dlzlgt
      ENDIF

c      write(6,*) 'z,kz,czz,cz=',z,kz,czz,cz
c      czz = 1.0
c      cz = 0.0
      
      cez = 0.0
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bndc(kbdy,xory,yorz,alfa,gbdy)
c
c     input derivative boundary conditions to mud3sp
c     at lower x,y,z boundaries
c
      implicit none
      real alfa
      integer kbdy
      real xory,yorz,cons,gbdy,x,y,z,pxx,pyy,pzz,px,py,pz,pe
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/ftmsp3/xa,xb,yc,yd,ze,zf,tolmax,relmax

      real eztop, pottop
      
      common/bndcz/eztop,pottop
      
      
      IF ( (kbdy .ge. 1 .and. kbdy .le. 4)  ) THEN
        alfa = 0.0
        gbdy = 0.0
      ELSEIF ( kbdy.eq.6 ) THEN 
        alfa = 0.0
        gbdy = 0.0 ! eztop
        gbdy = eztop
      ENDIF
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE mudpacksetup(nx,ny,nz,dx,dy,dz,
     :           nyslm,nxslm,nzslm,
     :           dxx,dyy,dzz,dzmax,
     :           nxb2s,nyb2s,nzbt2s,nzb2s,nb2s,
     :           iixps,jjyqs,kkzrs,
     :           iiexs,jjeys,kkezs,
     :           llworks,
     :           nnxs,nnys,nnzs,
     :           bcx,bcy,
     :           istretch,imudstop,iunit)

      USE ELEC_MODULE
      USE COMMASMPI_MODULE
      
      implicit none
      
c input variables      
      integer nx,ny,nz
      integer nyslm,nxslm,nzslm
      real    dx,dy,dz,dzmax
      real dzz(nz),dxx(nx),dyy(ny)     ! dz(k),dx(i),dy(j)
      logical istretch
      integer imudstop
      integer nxb2s(2),nyb2s(2),nzbt2s(2),nzb2s(2),nb2s(2)
      integer iixps(2),jjyqs(2),kkzrs(2)
      integer iiexs(2),jjeys(2),kkezs(2)
      integer*8 llworks(2)
      integer nnxs(2),nnys(2),nnzs(2)
      integer iunit
      integer bcx,bcy
      
      
c local variables
      integer n
      integer*8 n1,n2,n3
      

         IF ( lightintx .ge. 1 ) dslight = dx/Float(lightintx)
         
#ifdef MPI
         nxslm = NInt((nxend-1)*dx/dslight) + 1
         nyslm = NInt((nyend-1)*dy/dslight) + 1
#else
         nxslm = NInt((nx-1)*dx/dslight) + 1
         nyslm = NInt((ny-1)*dy/dslight) + 1
#endif

         IF ( Abs(dx/dslight - Float(NInt(dx/dslight) ))
     :        .gt. 0.001 )  Then
          write(0,*) 'dslight is not an integer multiple of dx!'
          write(iunit,*) 'dslight is not an integer multiple of dx!'
          imudstop = 1
         ENDIF

#ifdef MPI
c      if (debug_mpi) write(0,*) "DRIVER: DEBUG 1"
#else
!      if (ndebug .gt. 0) write(iunit,*) "DRIVER: DEBUG 1"
#endif

c        dzmax = gt(1,1,nz,mzdist) - gt(1,1,nz-1,mzdist)
        dzmax = dzz(nz-1)

        write(iunit,*) 'dzmax = ',dzmax

c
c IF vertical spacing is constant, then can have dslight < dz
c

         IF ( dzz(1) .eq. dzz(nz-1) ) THEN

        !   dslightz = dz ! Min(dz,dslight)
        !   IF ( dslight .le. 0.5*dz ) 
           dslightz = dz/Max( 1, NInt(dz/dslight) )

         IF ( lightintz .eq. -1 ) lightintz = Max( 1, NInt( dz/dslight ) )
         IF ( lightintz .ge. 1 ) dslightz = dz/Float(lightintz)

           nzslm = (nz-1)*lightintz + 1
           
           istretch = .false.
         ELSE
           istretch = .true.
           IF ( lightintz .eq. 2 ) THEN
!            lightintz = 1
            nzslm = lightintz*nz
           dslightz = dzmax/lightintz
!             write(0,*) 'WARNING: lightintz .eq. 2 not completed yet for stretched vertical grid!'
!             STOP
           
           ELSEIF ( lightintz .gt. 2 ) THEN
             write(0,*) 'lightintz .gt. 2 not supported yet for stretched vertical grid!'
             STOP
           ELSE
           lightintz = 1
           dslightz = dzmax
           nzslm = nz
           ENDIF
         ENDIF

        write(iunit,*) 'dslightz = ',dslightz
         
         nb2s(1) = Min( maxbzone, Max(NInt(4000./dx) , 8) )
         nb2s(2) = nb2s(1)*Nint(dx/dslight) ! Min( maxbzone, Max(NInt(4000./dslight) , 8) )
         

c         ktopmin = Min(8, Nint(8.0*500./dzmax) )
         nxb2s(1) = nx + 2*nb2s(1) + 1
         nyb2s(1) = ny + 2*nb2s(1) + 1
         nzbt2s(1) = nz + Max(8,NInt(4000./dzmax + 0.5) ) ! nb2s(2)
         nzb2s(1) = nzbt2s(1) + 1
         

         IF ( bcx .eq. 2 ) THEN
           nxb2s(1) = nx
         ENDIF
         IF ( bcy .eq. 2 ) THEN
           nyb2s(1) = ny
         ENDIF

         nxb2s(2) = nxslm + 2*nb2s(2) + 1
         nyb2s(2) = nyslm + 2*nb2s(2) + 1
         
!         IF ( istretch .and. lightintz .ge. 1 ) THEN
!         
!         ELSE
         nzbt2s(2) = nzslm + Max(8,NInt(4000./dslightz + 0.5) ) ! nb2s(2)
         nzb2s(2) = nzbt2s(2) + 1
!         ENDIF

         IF ( bcx .eq. 2 ) THEN
           nxb2s(2) = nxslm
         ENDIF
         IF ( bcy .eq. 2 ) THEN
           nyb2s(2) = nyslm
         ENDIF

         IF ( ny .eq. 2 ) THEN ! 2D case
           nyb2s(1) = ny
           nyb2s(2) = ny
           nyslm = ny
         ENDIF
         
         DO n=1,2

         write(iunit,*) 'mudsetup: n,nxb2s,nyb2s,nzb2s = ',n,nxb2s(n),nyb2s(n),nzb2s(n)


         IF ( n .eq. 3 .and. ( lightintx .le. 2 .or. Mod(lightintx,2) .eq. 0) ) THEN
c            nxb2s(2) = lightintx*(nxb2s(1) - 1)
c            nxb2s(2) = nxb2s(2) + 1
           iixps(2) = iixps(1)
           iiexs(2) = lightintx/2 + iiexs(1)
           
           write(iunit,*) 'lightintx,iixps,iiexs(2) = ',lightintx,iixps(2),iiexs(2)
           
!         ELSEIF ( n .eq. 2 .and. lightintz .eq. 1 ) THEN
            
         ELSE 

         IF ( nxb2s(n) .le. 33 ) THEN
           iixps(n) = 2
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 41 ) THEN
           iixps(n) = 5
           iiexs(n) = 4
         ELSEIF ( nxb2s(n) .le. 49 ) THEN
           iixps(n) = 3
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 65 ) THEN
           iixps(n) = 2
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 81 ) THEN
           iixps(n) = 5
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 97 ) THEN
           iixps(n) = 3
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 113 ) THEN
           iixps(n) = 7
           iiexs(n) = 5
         ELSEIF ( nxb2s(n) .le. 129 ) THEN
           iixps(n) = 2
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 161 ) THEN
           iixps(n) = 5
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 193 ) THEN
           iixps(n) = 3
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 225 ) THEN
           iixps(n) = 7
           iiexs(n) = 6
         ELSEIF ( nxb2s(n) .le. 257 ) THEN
           iixps(n) = 2
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 321 ) THEN
           iixps(n) = 5
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 385 ) THEN
           iixps(n) = 3
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 449 ) THEN
           iixps(n) = 7
           iiexs(n) = 7
         ELSEIF ( nxb2s(n) .le. 513 ) THEN
           iixps(n) = 2
           iiexs(n) = 9
         ELSEIF ( nxb2s(n) .le. 641 ) THEN
           iixps(n) = 5
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 768 ) THEN
           iixps(n) = 3
           iiexs(n) = 9
         ELSEIF ( nxb2s(n) .le. 896 ) THEN
           iixps(n) = 7
           iiexs(n) = 8
         ELSEIF ( nxb2s(n) .le. 1025 ) THEN
           iixps(n) = 2
           iiexs(n) = 10
         ELSEIF ( nxb2s(n) .le. 1281 ) THEN
           iixps(n) = 5
           iiexs(n) = 9
         ELSEIF ( nxb2s(n) .le. 1535 ) THEN
           iixps(n) = 3
           iiexs(n) = 10
         ELSEIF ( nxb2s(n) .le. 2561 ) THEN
           iixps(n) = 5
           iiexs(n) = 10
         ELSE
           imudstop = 1
           write(iunit,*) 
     :      'WARNING! nxb2 too large for mudpack, will stop',
     :       ' if ipelec .gt. 0'
         ENDIF
         
         ENDIF
         
         IF ( nyb2s(n) .eq. nxb2s(n) ) THEN
           jjyqs(n) = iixps(n)
           jjeys(n) = iiexs(n)
         ELSEIF ( ny .eq. 2 ) THEN
           jjyqs(n) = 1
           jjeys(n) = 1
         ELSEIF ( nyb2s(n) .le. 33 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 41 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 4
         ELSEIF ( nyb2s(n) .le. 49 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 65 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 81 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 97 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 113 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 5
         ELSEIF ( nyb2s(n) .le. 129 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 161 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 193 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 225 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 6
         ELSEIF ( nyb2s(n) .le. 257 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 321 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 385 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 449 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 7
         ELSEIF ( nyb2s(n) .le. 513 ) THEN
           jjyqs(n) = 2
           jjeys(n) = 9
         ELSEIF ( nyb2s(n) .le. 641 ) THEN
           jjyqs(n) = 5
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 769 ) THEN
           jjyqs(n) = 3
           jjeys(n) = 9
         ELSEIF ( nyb2s(n) .le. 897 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 8
         ELSEIF ( nyb2s(n) .le. 1793 ) THEN
           jjyqs(n) = 7
           jjeys(n) = 9
         ELSE
           imudstop = 1
           write(iunit,*) 
     :      'WARNING! nyb2 too large for mudpack, will stop',
     :       ' if ipelec .gt. 0'
         ENDIF
         
         IF ( n .eq. 3 .and. ( lightintz .le. 2 .or. Mod(lightintz,2) .eq. 0) ) THEN
            nzbt2s(2) = lightintz*(nzbt2s(1) - 1)
            nzb2s(2) = nzbt2s(2) + 1
           kkzrs(2) = kkzrs(1)
           kkezs(2) = lightintz/2 + kkezs(1)
           
!         ELSEIF ( n .eq. 2 .and. lightintz .eq. 1 ) THEN
            
         ELSE 
         
         IF ( nzb2s(n) .le. 41 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 4
         ELSEIF ( nzb2s(n) .le. 49 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 5
         ELSEIF ( nzb2s(n) .le. 65 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 81 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 5
         ELSEIF ( nzb2s(n) .le. 97 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 113 ) THEN
           kkzrs(n) = 7
           kkezs(n) = 5
         ELSEIF ( nzb2s(n) .le. 129 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 7
         ELSEIF ( nzb2s(n) .le. 161 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 193 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 7
         ELSEIF ( nzb2s(n) .le. 225 ) THEN
           kkzrs(n) = 7
           kkezs(n) = 6
         ELSEIF ( nzb2s(n) .le. 257 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 8
         ELSEIF ( nzb2s(n) .le. 321 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 7
         ELSEIF ( nzb2s(n) .le. 385 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 8
         ELSEIF ( nzb2s(n) .le. 449 ) THEN
           kkzrs(n) = 7
           kkezs(n) = 7
         ELSEIF ( nzb2s(n) .le. 513 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 9
         ELSEIF ( nzb2s(n) .le. 641 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 8
         ELSEIF ( nzb2s(n) .le. 769 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 9
         ELSEIF ( nzb2s(n) .le. 897 ) THEN
           kkzrs(n) = 7
           kkezs(n) = 8
         ELSEIF ( nzb2s(n) .le. 1025 ) THEN
           kkzrs(n) = 2
           kkezs(n) = 10
         ELSEIF ( nzb2s(n) .le. 1281 ) THEN
           kkzrs(n) = 5
           kkezs(n) = 9
         ELSEIF ( nzb2s(n) .le. 1535 ) THEN
           kkzrs(n) = 3
           kkezs(n) = 10
         ELSEIF ( nzb2s(n) .le. 1793 ) THEN
           kkzrs(n) = 7
           kkezs(n) = 9
         ELSE
           imudstop = 1
           write(iunit,*) 
     :      'WARNING! nzb2 too large for mudpack, will stop',
     :       ' if ipelec .gt. 0',nzb2s(n)
         ENDIF
         
         ENDIF
         
         nnxs(n) = iixps(n)*2**(iiexs(n)-1)+1
         nnys(n) = jjyqs(n)*2**(jjeys(n)-1)+1
         nnzs(n) = kkzrs(n)*2**(kkezs(n)-1)+1

         IF ( bcx .eq. 2 ) THEN
           IF ( nnxs(n) .ne. nxb2s(n) ) THEN
             write(0,*) 'STOP! must have nnxs(n) = nxb2s(n) ',
     :                  n,nnxs(n),nxb2s(n)
           STOP
           ENDIF
         ENDIF
         IF ( bcy .eq. 2 ) THEN
           IF ( nnys(n) .ne. nyb2s(n) ) THEN
             write(0,*) 'STOP! must have nnys(n) = nyb2s(n) ',
     :                  n,nnys(n),nyb2s(n)
           STOP
           ENDIF
         ENDIF

#ifdef MPI
         IF ( nyend .gt. 2 ) THEN
#else
         IF ( ny .gt. 2 ) THEN
#endif
           n1 = nnxs(n)
           n2 = nnys(n)
           n3 = nnzs(n)
!           llworks(n) = 7*(nnxs(n)+2)*(nnys(n)+2)*(nnzs(n)+2)/2
           llworks(n) = 7*((n1+2)*(n2+2)*(n3+2)/2)
           write(iunit,*) 'llworks(',n,') = ',llworks(n)
         ELSE
           llworks(n) = 5*((nnxs(n))*(nnzs(n)) + 2*(nnxs(n)+nnzs(n)) )
         ENDIF

         write(iunit,*) 'mudpack dimensions for grid ',n,' are ',
     :      nnxs(n),nnys(n),nnzs(n)

         nxb2s(n) = nnxs(n)
         nyb2s(n) = nnys(n)

         ENDDO ! n

         nz1d = nnzs(1)
         nz1dlgt = nnzs(2)

         
         IF ( nz1d .gt. nzmax1 ) THEN
           write(iunit,*) 
     :     'Increase size of nzmax1 in elec_module to at least',
     :      nnzs(1)
           write(0,*) 'Increase size of nzmax1 in driver to at least',
     :      nnzs(1)
           STOP
         ENDIF

         IF ( nz1dlgt .gt. nzmax1lgt ) THEN
           write(iunit,*) 
     :     'Increase size of nzmax1lgt in elec_module to at least',
     :      nnzs(2)
           write(0,*) 'Increase size of nzmax1lgt in driver to at least',
     :      nnzs(2)
           STOP
         ENDIF
         

      RETURN
      END
      
      
