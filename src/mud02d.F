c
c     mud2d: Calls 2-d Poisson solver for 2-d mode of the model
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
      subroutine mud2d(nx,ny,nz,nor,
     :   dx,dy,dz,gxt,gyt,gzt,
     :   nnx,nny,nnz,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,
     :   llwork, istretch,
     :   t0,ibg,
     :   elec,iestag,iunit,lcalce, bcx, bcy)
      
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
      integer bcx, bcy
      
c input arrays: t0 holds net charge
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)
c      real pot(nx,ny,nz)
      
c   local variables/arrays
      
      real, allocatable ::  pot2(:,:) 
      real, allocatable ::  pot3d(:,:,:) 
      real, allocatable ::  rhs(:,:)
      real, allocatable ::  work(:)
c      real pot2(nnx,nny,nnz) 
c      real rhs(nnx,nny,nnz)
c      real work(llwork)
      integer ix,jy,kz, i,j,k
      
c ####
c   begin execute
c ####      
       allocate (pot2(nnx,nnz))
       allocate (rhs(nnx,nnz))
       allocate (work(llwork))
      
!       write(0,*) 'call mud2d1'
       call mud2d1(nx,ny,nz,nor,
     :   dx,dy,dz,
     :   nnx,nny,nnz,
     :   iixp,jjyq,kkzr,iiex,jjey,kkez,
     :   nbw,nbe,nbs,nbn,
     :   llwork,istretch,
     :   t0,ibg,ibdcnd,jbdcnd,kbdcnd,kdbcdn,dlz,
     :   elec(ipot)%flt3d,potfair,rhofair,z1d2,nzmax1,
     :   pot2,pot2,rhs,rhs,work,iunit, bcx, bcy)
      
       IF ( iestag .eq. 1 .and. lcalce ) THEN

       allocate (pot3d(nnx,ny,nnz))

        DO k = 1,nnz
          DO i = 1,nnx
            pot3d(i,1,k) = pot2(i,k)
          ENDDO
        ENDDO
        call efieldstag
     >  (nx,ny,nz,nor
     >  ,elec
     :  ,gxt,gyt,gzt
     :  ,nbw,nbe,-1,0,nnz
     :  ,pot3d
     :  ,iunit)  
       
       deallocate  (pot3d)

       ENDIF
       
       deallocate (pot2)
       deallocate (rhs)
       deallocate (work)
      RETURN
      END
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mud2d1(nx,ny,nz,nor,
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

      real phi(nnx,nnz),rhs(nnx,nnz),work(llwork)
      real phi0(-nbw:nx+nbe,0:nnz-1)
      real rhs0(-nbw:nx+nbe,0:nnz-1)
      
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
c     put integer and floating point argument names in contiguous
c     storage for labelling in vectors iprm,fprm
c
      integer*8 iprm(16)
      integer mgopt(4)
      real fprm(6)
      integer*8 intl,nxa,nxb1,nyc,nyd,ixp,jyq,iex,jey,nxm,nym,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
! NOTE: need these dummy common blocks to align the variable names
! in memory so that the equivalence with iprm etc. will work correctly.
      common/itmud2sp/intl,nxa,nxb1,nyc,nyd,ixp,jyq,iex,jey,nxm,nym,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)

c
c     put integer and floating point arguments names in contiguous
c     storeage for labelling
c
c      integer intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
c     +         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
c      common/itmsp3/intl,nxa,nxb1,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
c     +         kez,nxm,nym,nzm,iguess,maxcy,method,nwork,lwrkqd,itero
c      real xa,xb,yc,yd,ze,zf,tolmax,relmax
c      common/ftmsp3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz2,x,y,z,cxx,cyy,czz,cx,cy,cz,cex,cey,cez,ce,errm
      real pxx,pyy,pzz,px,py,pz,pe
      integer i,j,k,ierror
c      equivalence(intl,iprm)
c      equivalence(xa,fprm)
c      equivalence(phi(1,1,2),phi0(1,1,1))
c      equivalence(rhs(1,1,2),rhs0(1,1,1))
c
c     declare coefficient and boundary condition input subroutines external
c
      external cfx,cfy,cfz,cfzu,bndc2d
      
      logical :: debug = .false.
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
      IF ( debug ) write(0,*) 'start of mud2d1'
      
      nxa = 2
      nxb1 = 2
      IF ( bcx .eq. 2 ) THEN
       nxa = 0
       nxb1 = 0
      ENDIF
      nyc = 1
      nyd = 1
c      nze = 1 
c      nzf = 1  
      IF ( kbdcnd .eq. 2 ) THEN
      nyc = 1
      nyd = 1
      ELSEIF ( kbdcnd .eq. 1 ) THEN
      nyc = 1
      nyd = 1  ! changed from 2 to 1 on 12/31/2004 by erm 
               ! solution seems to be more stable as a result, especially
               ! when there is an upper level grid stretch.
               ! Also seems to preserve net charge better than nzf = 2
      ELSEIF ( kbdcnd .eq. 3 ) THEN
      nyc = 2
      nyd = 2
      ENDIF
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = kkzr! jjyq
      iex = iiex
      jey = kkez ! jjey
      nxm = nnx
      nym = nnz
c      nzm = nnz
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
      mgopt(4) = 1  ! intpol (default 3)
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
c      yc = 0.
c      yd = (nym-1)*dy
      yc = -dlz*0.5
      yd = yc + (nym-1)*dlz
      
      IF ( debug ) print*, 'xa,xb,yc,yd,dx,dlz = ',xa,xb,yc,yd,dx,dlz
      
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
c      do jy = -nbs, ny+nbn
      do ix = -nbw, nx+nbe
      
      IF ( kz .ge. 1 .and. kz .le. nz-1 .and. 
     :     ix .ge. 1 .and. ix .le. nx-1 ) THEN
        rhs0(ix,kz) = -t0(ix,1,kz)/eperao
        phi0(ix,kz) = pot(ix,1,kz)
      ELSEIF ( kz .eq. 0 ) THEN
        rhs0(ix,kz) = 0.0
        phi0(ix,kz) = 0.0
      ELSE
        rhs0(ix,kz) = -rhofair(kz)/eperao
        phi0(ix,kz) = potfair(kz)
      ENDIF

      end do
c      end do
      end do

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
      write (6,101) (iprm(i),i=1,15)
  101 format(/' integer input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2
     +/' nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2,
     +/' method = ',i2, ' work space estimate = ',i7)
      write (6,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option arguments ',
     +/' kcycle = ',i2,
     +/' iprer = ',i2,
     +/' ipost = ',i2,
     +/' intpol = ',i2)
      write(6,103) xa,xb,yc,yd,tolmax
  103 format(/' floating point input parameters ',
     +/' xa = ',f12.3,' xb = ',f12.3,' yc = ',f12.3,' yd = ',f12.3,
     +/' tolerance (error control) =   ',e10.3)
     
      
      ENDIF ! (ndebug)
      tim1 = commas_dtime(timarr) 
c
c     discretize pde
c
      IF ( debug ) write(0,*) 'call mud2sp'


      IF (ndebug .ge. 1) THEN
        write(*,104) intl
      ENDIF
  104 format(/' discretization call to mud3sp', ' intl = ', i2)
      IF ( istretch ) THEN
      call mud2sp(iprm,fprm,work,cfx,cfz,bndc2d,rhs,phi,mgopt,ierror)
      ELSE
      call mud2sp(iprm,fprm,work,cfx,cfzu,bndc2d,rhs,phi,mgopt,ierror)
      ENDIF
      IF (ndebug .ge. 1) THEN
        write (*,105) ierror,iprm(16)
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
      call mud2sp(iprm,fprm,work,cfx,cfz,bndc2d,rhs,phi,mgopt,ierror)
      ELSE
      call mud2sp(iprm,fprm,work,cfx,cfzu,bndc2d,rhs,phi,mgopt,ierror)
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
        write(6,*) k,z1d2(k,1),phi0(nx/2,k),
     :   -0.5/dlz*(phi0(nx/2,k+1)-phi0(nx/2,k-1))*
     :     z1d2(k,3),
     :   t0(nx/2,ny/2,k)
c     :   ,cz
      ENDDO
      k=1
        write(6,*) k,z1d2(k,1),phi0(nxm/2,k),
     :   -0.5/dlz*(phi0(nx/2,k+1)-0.0)*z1d2(k,3),
     :   t0(nx/2,ny/2,k)

      k=1
        write(6,*) k-1,-z1d2(k,1),phi(nxm/2,k)
      
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
         pot(ix,jy,kz) = phi0(ix,kz)
      end do
      end do
      end do


      
      RETURN
      END
      
c ############################################################
c ############################################################

      subroutine bndc2d(kbdy,xory,alfa,gbdy)
c
c     input derivative boundary conditions to mud3sp
c     at lower x,y,z boundaries
c
      implicit none
      real alfa
      integer kbdy
      real xory,yorz,cons,gbdy,x,y,z,pxx,pyy,pzz,px,py,pz,pe
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
!      common/ftmud2sp/xa,xb,yc,yd,ze,zf,tolmax,relmax

      real eztop, pottop
      
      common/bndcz/eztop, pottop
      
      
      IF ( (kbdy .ge. 1 .and. kbdy .le. 3)  ) THEN
        alfa = 0.0
        gbdy = 0.0
      ELSEIF ( kbdy.eq.4 ) THEN 
        alfa = 0.0
        gbdy = 0.0 ! eztop
        gbdy = eztop
      ENDIF
      
      return
      end
