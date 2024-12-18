!-----------------------------------------------------------------------------
!
! ///////////////////// BEGIN \\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\ SUBROUTINE INT_TDLBC ////////////////////
!
! INT_TDLBC interpolates to get time dependent lateral boundary conditions
!-----------------------------------------------------------------------------
! Created by Mike Buban, December 6, 2007
!-----------------------------------------------------------------------------

    SUBROUTINE INT_TDLBC(gd,varname,sintew,sintns,time,nx,ny,nz,kmax)
       
  USE GRID_MODULE
  USE INIT_MODULE, only : athet,aqrat,au,av

  USE COMMASMPI_MODULE
  
  implicit none

#ifdef MPI
   include "mpif.h"
#endif

  TYPE(GRID) :: gd

    integer            :: i,j,k,m,l,n
    integer            :: time,nin(2),ninv(2)
    integer            :: nx,ny,nz
    integer            :: ibeg,iend
    integer            :: kmax

    real, pointer :: dzc(:)
    real, pointer :: dze(:)
    real, pointer :: zcdx(:,:)
    real, pointer :: zedx(:,:)

    real               :: x_sw_loc,y_sw_loc
    real               :: svar(nx,ny,nz),ufl(11),vfl(11)



!
!.... must be saved for three RK steps per one dynamic time step.  Allocated in initsubs.F90.
!

!    real  ::  athet(nxend + 1,nyend + 1,nzend + 1,2)
!    real  ::  aqrat(nxend + 1,nyend + 1,nzend + 1,2)
!    real  ::  au(nxend + 1,nyend + 1,nzend + 1,2)
!    real  ::  av(nxend + 1,nyend + 1,nzend + 1,2)


!    real, save         :: athet(62,62,26,2),aqrat(62,62,26,2)
!    real, save         :: au(62,62,26,2),av(62,62,26,2)

    real               :: tsint(nx,ny,nz,2),tm
    real               :: sint(nx,ny,nz) !,pert(nx,ny,nz)
    real               :: tsintew(ny,nz,2,2),tsintns(nx,nz,2,2)
    real               :: sintew(ny,nz,2),sintns(nx,nz,2)
!    real, allocatable :: sint(:,:,:),pert(:,:,:)
    real, save         :: dxin,dyin,dzin
    real dx,dy,dz,t1,t2
    real grid_xin(62),grid_yin(62),grid_zin(26)
    real f1,f2,f3,f4,f5,f6,f7,f8
    real grid_x(1000),grid_y(1000),grid_z(1000)
    real, save  :: u(100),v(100),hgt(100),p(100),t(100),rh(100)
    real, save  :: reflect(101,101,25,2),uvel(101,101,25,2)
    real, save  :: vvel(101,101,25,2),wvel(101,101,25,2)
    real gxin,gyin,gx,gy,xmx,ymx
    real, save  :: nrint,nlint
    real               :: xmin, ymin ! coordinates (m) of southwest corner of domain
!    real               :: hwlmax,hwlmin,vwlmax,vwlmin
    real,    allocatable :: sd(:,:,:)       ! [nx,ny,nz] standard deviation of smooth perturbations
    real,    pointer     :: xc(:)
    real,    pointer     :: xe(:)
    real,    pointer     :: yc(:)
    real,    pointer     :: ye(:)
    real,    pointer     :: zc(:)
    real,    pointer     :: ze(:)
    real  dt
    integer, pointer     :: member
    integer, save        :: nxin,nyin,nzin
    integer, save        :: nxr,nyr,nzr
    integer ix1,ix2,iy1,iy2,iz1,iz2
    integer, save      :: nrtimes,nltimes,ninold,ninvold
    character(LEN=4)   :: fld
    character(LEN=10)  :: varname    
    character(LEN=70), save  :: infile(100)
    character(LEN=70), save  :: inrad(100)
    character(LEN=50)  :: bc_file

    TYPE(ATTRIBUTE), pointer :: bc_file_name
!    TYPE(ATTRIBUTE), pointer :: member

!   TYPE(VARIABLE)     :: pertu,pertv
!   TYPE(VARIABLE)     :: pertw,pertt,pertq

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  

!-----------------------------------------------------------------------

!   CALL cld_cpu('TDLBC')


    CALL GET_VARIABLE(gd, 'DX', dx)
    CALL GET_VARIABLE(gd, 'DY', dy)
    CALL GET_VARIABLE(gd, 'DZ', dz)


!    CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!    CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
    CALL GET_VARIABLE(gd,'DZC',dzc)
    CALL GET_VARIABLE(gd,'DZE',dze)
    CALL GET_VARIABLE(gd,'XC',xc)
    CALL GET_VARIABLE(gd,'XE',xe)
    CALL GET_VARIABLE(gd,'YC',yc)  
    CALL GET_VARIABLE(gd,'YE',ye)  
    CALL GET_VARIABLE(gd,'ZC',zc)  
    CALL GET_VARIABLE(gd,'ZE',ze)  
    CALL GET_VARIABLE(gd,'DT',dt)

!    print*,varname,sint(1,1,1)


!     print*,'in tdlbc', time, dt, dzc(1),dzc(2),varname

!.....read in list of infile/radar names

      if(time .eq. dt) then

      print*,time,'in tdlbc',' read infile',varname, my_rank

      open(13,file='tdlbc.in',status='old')
      read(13,3101) nrtimes, nrint
 3101 format(1x,i2,1x,f5.0)
      do i=1,nrtimes
      read(13,3012) inrad(i)
      end do
      read(13,3101) nltimes, nlint
 3012 format(1x,a70)
      do i=1,nltimes
      read(13,3012) infile(i)
      end do
      read(13,3013) dxin,dyin,dzin
      read(13,3014) nxin,nyin,nzin
 3013 format(3(1x,f5.0))
 3014 format(3(1x,i3))
      close(13)

!.....read in base state sounding

      open(14,file='pht.snd',status='old')
      read(14,3015)
 3015 format(1x)
      do i=1,25
      read(14,3016) hgt(i),p(i),t(i),rh(i),u(i),v(i)
!      print*,hgt(i),p(i),rh(i),u(i),v(i),i
      end do
 3016 format(f4.0,4x,f7.3,4x,f5.1,3(4x,f4.1))
      close(14)


      ninold=0
      ninvold=0

      end if


      nin(1)= (time/nlint) + 1
      nin(2)= nin(1) +1

      ninv(1)= (time/nrint) + 1
      ninv(2)= ninv(1) +1


      if(nin(2).gt.nltimes) then
      nin(2)=nltimes
      end if

      if(ninv(2).gt.nrtimes) then
      ninv(2)=nrtimes
      end if

!      print*,nin(1),nin(2),ninv(1),ninv(2),nlint,nrint,nltimes,nrtimes,time

      if(ninv(1).gt.ninvold) then
      ninvold=ninv(1)

      print*,time,'in tdlbc', ' read radar', my_rank


      do n=1,2

!.....read in radar fields

      open(50,file=inrad(ninv(n)),status='old')
      read(50,900) fld
  900 format(a12)
      read(50,501) nxr,nyr,nzr
  501 format(3(i3,1x))
      read(50,502) (((reflect(i,j,k,n),i=1,nxr),j=1,nyr),k=1,nzr)
  502 format(9f10.3)
      read(50,503) fld
  503 format(a1)
      read(50,501) nxr,nyr,nzr
      read(50,502) (((uvel(i,j,k,n),i=1,nxr),j=1,nyr),k=1,nzr)
      read(50,503)
      read(50,501) nxr,nyr,nzr
      read(50,502) (((vvel(i,j,k,n),i=1,nxr),j=1,nyr),k=1,nzr)
      read(50,503)
      read(50,501) nxr,nyr,nzr
      read(50,502) (((wvel(i,j,k,n),i=1,nxr),j=1,nyr),k=1,nzr)
      close(50)

      end do


      end if  ! for new set of radar files

      if(nin(1).gt.ninold) then
      ninold=nin(1)

      print*,time,'in tdlbc', ' read lag', my_rank

      do n=1,2

!      print*,infile(n),n,nin(1),nin(2),varname,nx,ny,nz,time
      open(15, file=infile(nin(n)),status='unknown')
      read(15,4023)(((athet(m,l,k,n),m=1,61),l=1,61),k=1,25)
      read(15,4023)(((aqrat(m,l,k,n),m=1,61),l=1,61),k=1,25)
!      read(15,4023)(((au(m,l,k,n),m=1,61),l=1,61),k=1,25)
!      read(15,4023)(((av(m,l,k,n),m=1,61),l=1,61),k=1,25)
      close(15)

      end do

 4023 format (61(f10.4,1x))

      end if


!      print*,nxr,nyr,nzr,time

      do n=1,2  ! bounding times

!......fill in above lag analysis with sounding

      do i=1,nxr
      do j=1,nyr
      do k=12,25
      uvel(i,j,k,n)=u(k)
      vvel(i,j,k,n)=v(k)
!      if(varname.eq.'U') then
!      print*, uvel(i,j,k,n),u(k),i,j,k,n,time
!      end if

      end do
      end do
      end do
      

!.....fill in au,av

      do i=21,81
      do j=1,61
      do k=1,25
      m=i-20
      au(m,j,k,n)=uvel(i,j,k,n) 
      av(m,j,k,n)=vvel(i,j,k,n)
!      if(varname.eq.'V'.and.i.eq.1) then
!      print*,au(m,j,k,n),av(m,j,k,n),i,j,k,n,time
!      end if
      end do 
      end do
      end do


!.....filter bad winds

      do i=1,61
      do j=1,61
      do k=1,25
      if((av(i,j,k,n).lt.0).or.(au(i,j,k,n).lt.-8.).or.(av(i,j,k,n).gt.40.).or.(au(i,j,k,n).gt.30)) then  
      au(i,j,k,n)=u(k)
      av(i,j,k,n)=v(k)
      end if
!      print*,au(i,j,k),av(i,j,k),i,j,k
      end do
      end do
      end do

      end do ! times



!    print*,'before tri-lin', my_rank

! Tri-linear interpolation of input data to model grid
! msb 11/19/07
!
! dxin,dyin,dzin => input grid spacing
! nxin,nyin,nzin => input grid dimensions
!

!  dx=((real(nxin)-1)/(real(nx)-1))*dxin
!  dy=((real(nyin)-1)/(real(ny)-1))*dyin
!  dz=((real(nzin)-1)/(real(nz)-1))*dzin


  do i=1,nxin+1
  grid_xin(i)=(real(i)-1)*dxin
  enddo

  do i=1,nyin+1
  grid_yin(i)=(real(i)-1)*dyin
  enddo

  do i=1,nzin+1
  grid_zin(i)=(real(i)-1)*dzin
  enddo

#ifdef MPI
  do i=nxbeg,nxend
#else  
  do i=1,nx
#endif 
  grid_x(i)=(real(i)-1)*dx
  enddo

#ifdef MPI
  do i=nybeg,nyend
#else  
  do i=1,ny
#endif 
  grid_y(i)=(real(i)-1)*dy
  enddo

!...msb 5/28/08 change for stretched vertical grid


#ifdef MPI
  grid_z(nzbeg)=0.0  
  do i=nzbeg+1,nzend
#else  
  grid_z(1)=0.0
  do i=2,nz
#endif 
  grid_z(i)=grid_z(i-1)+(1/dzc(i-1))
  enddo



   gxin=(real(nxin)-1)*dxin
   gyin=(real(nyin)-1)*dyin

#ifdef MPI   
   gx=(real(nxend)-1)*dx
   gy=(real(nyend)-1)*dy
   xmx=(gxin/gx)*(nxend-1) + 1
   ymx=(gyin/gy)*(nyend-1) + 1
#else   
   gx=(real(nx)-1)*dx
   gy=(real(ny)-1)*dy
   xmx=(gxin/gx)*(nx-1) + 1
   ymx=(gyin/gy)*(ny-1) + 1   
#endif


#ifdef MPI
   do j=jybeg,jyend
    do i=ixbeg,ixend
     do k=kzbeg,kzend
#else     
   do j=1,ny
    do i=1,nx
     do k=1,nz
#endif

     
#ifdef MPI
  ix1= (real(nxin)-1)/(real(nxend)-1)*(gx/gxin)*(i-1) + 1
  
  if(ix1.gt.nxin) then
  ix1=nxin
  end if
  ix2=ix1+1

  iy1= (real(nyin)-1)/(real(nyend)-1)*(gy/gyin)*(j-1) + 1
  if(iy1.gt.nyin) then
  iy1=nyin
  end if
  iy2= iy1+1
  
  iz1=grid_z(k)/dzin + 1
  iz2= iz1 +1
      
    
#else 
    
  ix1= (real(nxin)-1)/(real(nx)-1)*(gx/gxin)*(i-1) + 1
  
  if(ix1.gt.nxin) then
  ix1=nxin
  end if
  ix2=ix1+1
  
  iy1= (real(nyin)-1)/(real(ny)-1)*(gy/gyin)*(j-1) + 1
  if(iy1.gt.nyin) then
  iy1=nyin
  end if
  iy2= iy1+1

  iz1=grid_z(k)/dzin + 1
  iz2= iz1 +1
    
#endif


!     print*,i,j,k,ix1,ix2,iy1,iy2,iz1,iz2,n

     f1=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))
     f2=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))
     f3=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))
     f4=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_zin(iz2)-grid_z(k)))*(1./(dxin*dyin*dzin))

     f5=((grid_xin(ix2)-grid_x(i))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))
     f6=((grid_x(i)-grid_xin(ix1))*(grid_yin(iy2)-grid_y(j))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))
     f7=((grid_xin(ix2)-grid_x(i))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))
     f8=((grid_x(i)-grid_xin(ix1))*(grid_y(j)-grid_yin(iy1))*(grid_z(k)-grid_zin(iz1)))*(1./(dxin*dyin*dzin))

!     print*,f1,f2,f3,f4,f5,f6,f7,f8,'h'


   
     do n=1,2

     if(varname.eq.'TH') then

#ifdef MPI


     if(i.eq.nxbeg) then
     tsintew(j-jybeg+1,k-kzbeg+1,1,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
!     print*,tsintew(j-jybeg+1,k-kzbeg+1,1,n),j-jybeg+1,k-kzbeg+1,n,my_rank,time     
     end if
#else          
     if(i.eq.1) then
     tsintew(j,k,1,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if
#endif


#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if     
#else
     if(i.eq.nx) then   
     tsintew(j,k,2,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if
#endif
     

#ifdef MPI 
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if
#else         
     if(j.eq.1) then
     tsintns(i,k,1,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if
#endif

     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if
#else
     if(j.eq.ny) then
     tsintns(i,k,2,n)=f1*athet(ix1,iy1,iz1,n)+f2*athet(ix2,iy1,iz1,n)+f3*athet(ix1,iy2,iz1,n)+ &
     f4*athet(ix2,iy2,iz1,n)+f5*athet(ix1,iy1,iz2,n)+f6*athet(ix2,iy1,iz2,n)+ &
     f7*athet(ix1,iy2,iz2,n)+f8*athet(ix2,iy2,iz2,n)
     end if
#endif

     elseif(varname.eq.'QV') then

#ifdef MPI
     if(i.eq.nxbeg) then
     tsintew(j-jybeg+1,k-kzbeg+1,1,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#else          
     if(i.eq.1) then
     tsintew(j,k,1,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#endif

     
#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#else          
     if(i.eq.nx) then
     tsintew(j,k,2,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#endif
     

#ifdef MPI
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#else
     if(j.eq.1) then
     tsintns(i,k,1,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#endif
     
     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#else     
     if(j.eq.ny) then
     tsintns(i,k,2,n)=f1*aqrat(ix1,iy1,iz1,n)+f2*aqrat(ix2,iy1,iz1,n)+f3*aqrat(ix1,iy2,iz1,n)+ &
     f4*aqrat(ix2,iy2,iz1,n)+f5*aqrat(ix1,iy1,iz2,n)+f6*aqrat(ix2,iy1,iz2,n)+ &
     f7*aqrat(ix1,iy2,iz2,n)+f8*aqrat(ix2,iy2,iz2,n)
     end if
#endif

     elseif(varname.eq.'U'.or.varname.eq.'ULAG') then


#ifdef MPI
     if(i.eq.nxbeg) then
     tsintew(j-jybeg+1,k-kzbeg+1,1,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)
     end if
#else     
     if(i.eq.1) then
     tsintew(j,k,1,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)
     end if
#endif

     
#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &                          
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)                    
     end if
#else          
     if(i.eq.nx) then
     tsintew(j,k,2,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &                          
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)                    
     end if
#endif


#ifdef MPI 
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &                          
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)                                 
     end if
#else
     if(j.eq.1) then
     tsintns(i,k,1,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &                          
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)                                 
     end if
#endif

     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &                          
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)
     end if
#else     
     if(j.eq.ny) then
     tsintns(i,k,2,n)=f1*au(ix1,iy1,iz1,n)+f2*au(ix2,iy1,iz1,n)+f3*au(ix1,iy2,iz1,n)+ &
     f4*au(ix2,iy2,iz1,n)+f5*au(ix1,iy1,iz2,n)+f6*au(ix2,iy1,iz2,n)+ &                          
     f7*au(ix1,iy2,iz2,n)+f8*au(ix2,iy2,iz2,n)
     end if
#endif

     elseif(varname.eq.'V'.or.varname.eq.'VLAG') then

#ifdef MPI 
     if(i.eq.nxbeg) then 
     tsintew(j-jybeg+1,k-kzbeg+1,1,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ &
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)
     end if
#else
     if(i.eq.1) then
     tsintew(j,k,1,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ &
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)
     end if
#endif


#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ &
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)                                        
     end if
#else
     if(i.eq.nx) then
     tsintew(j,k,2,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ &
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)                                        
     end if
#endif


#ifdef MPI
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ & 
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)                                         
     end if
#else
     if(j.eq.1) then
     tsintns(i,k,1,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ & 
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)                                         
     end if
#endif

     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ &         
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)
     end if
#else
     if(j.eq.ny) then
     tsintns(i,k,2,n)=f1*av(ix1,iy1,iz1,n)+f2*av(ix2,iy1,iz1,n)+f3*av(ix1,iy2,iz1,n)+ &         
     f4*av(ix2,iy2,iz1,n)+f5*av(ix1,iy1,iz2,n)+f6*av(ix2,iy1,iz2,n)+ &
     f7*av(ix1,iy2,iz2,n)+f8*av(ix2,iy2,iz2,n)
     end if
#endif

!     print*,av(ix2,iy1,iz2,n),av(ix2,iy1,iz1,n),av(ix1,iy2,iz1,n),av(ix2,iy2,iz1,n),ix1,iy1,iz1
!     print*,tsint(i,j,k,n),varname,i,j,k,n

     else    !if(varname.eq.'W'.or.varname.eq.'WLAG'.or.varname.eq.'KM') then
#ifdef MPI
     tsintew(j-jybeg+1,k-kzbeg+1,1,n)=0.0
     tsintew(j-jybeg+1,k-kzbeg+1,2,n)=0.0
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n)=0.0
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n)=0.0
#else
     tsintew(j,k,1,n)=0.0
     tsintew(j,k,2,n)=0.0
     tsintns(i,k,1,n)=0.0
     tsintns(i,k,2,n)=0.0
#endif     
     end if

     end do !times

!   CALL cld_cpu('I')

    enddo
      enddo
        enddo


!---------------------------------------------------
!
! End tri-linear interpolation

!--------------------------------------------------
!
! Now do the temporal interpolation between tsint(i,j,k,1)
! and tsint(i,j,k,2)
!
#ifdef MPI
   do j=jybeg,jyend
    do i=ixbeg,ixend
     do k=kzbeg,kzend
#else     
   do j=1,ny
    do i=1,nx
     do k=1,nz
#endif


    if(varname.eq.'U'.or.varname.eq.'ULAG'.or.varname.eq.'V'.or.varname.eq.'VLAG') then
    t1=(ninv(1)-1)*nrint
    t2=(ninv(2)-1)*nrint
    else
    t1=(nin(1)-1)*nlint
    t2=(nin(2)-1)*nlint
    endif
    tm=real(time)    

#ifdef MPI
    if(varname.eq.'U'.or.varname.eq.'ULAG'.or.varname.eq.'V'.or.varname.eq.'VLAG') then

    if(i.eq.nxbeg) then
    sintew(j-jybeg+1,k-kzbeg+1,1)=(1.0/nrint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,1,1)  & 
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,1,2))
         
    end if

    if(i.eq.nxend) then
    sintew(j-jybeg+1,k-kzbeg+1,2)=(1.0/nrint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,2,2))
    end if
    
    if(j.eq.nybeg) then
    sintns(i-ixbeg+1,k-kzbeg+1,1)=(1.0/nrint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,1,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,1,2))
    end if

    if(j.eq.nyend) then
    sintns(i-ixbeg+1,k-kzbeg+1,2)=(1.0/nrint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,2,2))
    end if

    else

    if(i.eq.nxbeg) then
    sintew(j-jybeg+1,k-kzbeg+1,1)=(1.0/nlint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,1,1)  &
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,1,2))
!    print*,tsintew(j-jybeg+1,k-kzbeg+1,1,1),tsintew(j-jybeg+1,k-kzbeg+1,1,2),j-jybeg+1,k-kzbeg+1,n,my_rank,time    
    end if

    if(i.eq.nxend) then
    sintew(j-jybeg+1,k-kzbeg+1,2)=(1.0/nlint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,2,2))
    end if

    if(j.eq.nybeg) then
    sintns(i-ixbeg+1,k-kzbeg+1,1)=(1.0/nlint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,1,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,1,2))
    end if
    
    if(j.eq.nyend) then
    sintns(i-ixbeg+1,k-kzbeg+1,2)=(1.0/nlint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,2,2))
    end if

    end if
#else
    if(varname.eq.'U'.or.varname.eq.'ULAG'.or.varname.eq.'V'.or.varname.eq.'VLAG') then
    sintew(j,k,1)=(1.0/nrint)*((t2-tm)*tsintew(j,k,1,1) + (tm-t1)*tsintew(j,k,1,2))
    sintew(j,k,2)=(1.0/nrint)*((t2-tm)*tsintew(j,k,2,1) + (tm-t1)*tsintew(j,k,2,2))
    sintns(i,k,1)=(1.0/nrint)*((t2-tm)*tsintns(i,k,1,1) + (tm-t1)*tsintns(i,k,1,2))
    sintns(i,k,2)=(1.0/nrint)*((t2-tm)*tsintns(i,k,2,1) + (tm-t1)*tsintns(i,k,2,2))
    else
    sintew(j,k,1)=(1.0/nlint)*((t2-tm)*tsintew(j,k,1,1) + (tm-t1)*tsintew(j,k,1,2))
    sintew(j,k,2)=(1.0/nlint)*((t2-tm)*tsintew(j,k,2,1) + (tm-t1)*tsintew(j,k,2,2))
    sintns(i,k,1)=(1.0/nlint)*((t2-tm)*tsintns(i,k,1,1) + (tm-t1)*tsintns(i,k,1,2))
    sintns(i,k,2)=(1.0/nlint)*((t2-tm)*tsintns(i,k,2,1) + (tm-t1)*tsintns(i,k,2,2))
    end if
#endif



!    print*,'end tri-lin',my_rank
!    if(varname.eq.'TH'.and.j.eq.nyend) then
!    print*,sintns(i-ixbeg+1,k-kzbeg+1,2),tsintns(i-ixbeg+1,k-kzbeg+1,2,1),tsintns(i-ixbeg+1,k-kzbeg+1,2,2),i-ixbeg+1,k-kzbeg+1,n,my_rank,time
!    print*,i-ixbeg+1,k-kzbeg+1
!    end if

    end do
   end do
  end do


!    if(varname.eq.'TH') then
!    write(0,"('--------------ADVECT: after tdlbc --------------')")
!    do j=jybeg,jyend
!    do k=kzbeg,kzend
!    do n=1,2
!    IF( .not. ( sintew(j,k,n) .lt. 350. .and. sintew(j,k,n) .gt. 250. ) ) then
!    print*,sintew(j-jybeg+1,k-kzbeg+1,n),tsintew(j-jybeg+1,k-kzbeg+1,n,1),tsintew(j-jybeg+1,k-kzbeg+1,n,2)
!    print*,j-jybeg+1,k-kzbeg+1,tm,my_rank
!    end if
!    end do
!    end do
!    end do
!
!     end if




!   CALL cld_cpu('J')


!    print*,'end tdlbc',my_rank
!------------------------------------------------------------
!
!  msb 3/3/09
!  Add perturbations to boundary conditions if running EnKF
!
!------------------------------------------------------------

       
!       CALL GET_ATTRIBUTE(gd,'MEMBER', member)

!       print*, member, varname

!     IF(member.gt.0) then

!       CALL GET_ATTRIBUTE(gd,'BC_FILE_NAME', bc_file_name)
!       CALL STRING_LIMITS( bc_file_name%str, ibeg, iend )
!       bc_file=bc_file_name%str(ibeg:iend)
!       tm=real(time)


!     IF( upert .gt. 0.0 .and. varname .eq. 'U') THEN
!       IF( perttype .eq. 1 ) THEN
!         CALL ADDNOISE(sint, upert, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 2 ) THEN
!         CALL ADDBLOBS(sint, upert, 0, rbubh, rbubv, nb, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 3 ) THEN
!        sd(:,:,:) = upert
!        CALL ADD_SMOOTH_PERTS(nx, ny, nz, rbubh, rbubv, sd, xc, yc, zc, sint, xe, yc, zc)
!       ENDIF
!       IF( perttype .eq. 4 ) THEN
!         CALL ADDSINES(sint, upert, 0, rbubh, rbubv, hwlmax, hwlmin, vwlmax, vwlmin, nb, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!         CALL ADDTSINES(tm,sint, rbubh, rbubv, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax, bc_file, 'U', pert)
!       ENDIF
!     ENDIF

!     IF( vpert .gt. 0.0 .and. varname .eq. 'V') THEN
!       IF( perttype .eq. 1 ) THEN
!         CALL ADDNOISE(sint, vpert, nx, ny, nz, xmin, ymin, xc, ye, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 2 ) THEN
!         CALL ADDBLOBS(sint, vpert, 0, rbubh, rbubv, nb, nx, ny, nz, xmin, ymin, xc, ye, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 3 ) THEN
!        sd(:,:,:) = vpert
!        write(*,*) 'ENSINIT:  Perturbing member ', n, ', variable = V '
!        CALL ADD_SMOOTH_PERTS(nx, ny, nz, rbubh, rbubv, sd, xc, yc, zc, sint, xc, ye, zc)
!       ENDIF
!       IF( perttype .eq. 4 ) THEN
!         CALL ADDSINES(sint, vpert, 0, rbubh, rbubv, hwlmax, hwlmin, vwlmax, vwlmin, nb, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!         CALL ADDTSINES(tm,sint, rbubh, rbubv, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax, bc_file, 'V', pert)
!       ENDIF
!     ENDIF

!     IF( wpert .gt. 0.0 .and. varname .eq. 'W') THEN
!       IF( perttype .eq. 1 ) THEN
!         CALL ADDNOISE(sint, wpert, nx, ny, nz, xmin, ymin, xc, yc, ze, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 2 ) THEN
!         CALL ADDBLOBS(sint, wpert, 0, rbubh, rbubv, nb, nx, ny, nz, xmin, ymin, xc, yc, ze, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 3 ) THEN
!        sd(:,:,:) = wpert
!        write(*,*) 'ENSINIT:  Perturbing member ', n, ', variable = W '
!        CALL ADD_SMOOTH_PERTS(nx, ny, nz, rbubh, rbubv, sd, xc, yc, zc, sint, xc, yc, ze)
!       ENDIF
!       w(:,:, 1) = 0.0  ! RESET BOUNDARY
!       w(:,:,nz) = 0.0  ! RESET BOUNDARY
!       IF( perttype .eq. 4 ) THEN
!         CALL ADDSINES(sint, wpert, 0, rbubh, rbubv, hwlmax, hwlmin, vwlmax, vwlmin, nb, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!         CALL ADDTSINES(tm,sint, rbubh, rbubv, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax, bc_file, 'W', pert)
!       ENDIF
!     ENDIF

!     IF( tpert .gt. 0.0 .and. varname .eq. 'TH') THEN
!       IF( perttype .eq. 1 ) THEN
!            print*,'before add noise',sint(1,1,1)
!         CALL ADDNOISE(sint, tpert, nx, ny, nz, xmin, ymin, xc, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!         print*,'after add noise', sint(1,1,1)
!       ENDIF
!       IF( perttype .eq. 2 ) THEN
!         CALL ADDBLOBS(sint, tpert, 1, rbubh, rbubv, nb, nx, ny, nz, xmin, ymin, xc, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 3 ) THEN
!        sd(:,:,:) = tpert
!        write(*,*) 'ENSINIT:  Perturbing member ', n, ', variable = TH '
!        CALL ADD_SMOOTH_PERTS(nx, ny, nz, rbubh, rbubv, sd, xc, yc, zc, sint, xc, yc, zc)
!       ENDIF
!       IF( perttype .eq. 4 ) THEN
!         print*,'before addsines', sint(1,1,1),pert(1,1,1)
!         CALL ADDSINES(sint, tpert, 0, rbubh, rbubv, hwlmax, hwlmin, vwlmax, vwlmin, nb, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!         CALL ADDTSINES(tm,sint, rbubh, rbubv, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax, bc_file, 'T', pert)
!         print*,'after addsines', sint(1,1,1),pert(1,1,1),tm
!       ENDIF
!     ENDIF

!     IF( qpert .gt. 0.0 .and. varname .eq. 'QV') THEN
!       IF( perttype .eq. 1 ) THEN
!         CALL ADDNOISE(sint, qpert, nx, ny, nz, xmin, ymin, xc, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 2 ) THEN
!         CALL ADDBLOBS(sint, qpert, 0, rbubh, rbubv, nb, nx, ny, nz, xmin, ymin, xc, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!       ENDIF
!       IF( perttype .eq. 3 ) THEN
!        sd(:,:,:) = qpert
!        write(*,*) 'ENSINIT:  Perturbing member ', n, ', variable = QV '
!        CALL ADD_SMOOTH_PERTS(nx, ny, nz, rbubh, rbubv, sd, xc, yc, zc, sint, xc, yc, zc)
!       ENDIF
!       IF( perttype .eq. 4 ) THEN
!         CALL ADDSINES(sint, qpert, 0, rbubh, rbubv, hwlmax, hwlmin, vwlmax, vwlmin, nb, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax)
!         CALL ADDTSINES(tm,sint, rbubh, rbubv, nx, ny, nz, xmin, ymin, xe, yc, zc, xbmin, xbmax, ybmin, ybmax, zbmin, zbmax, bc_file, 'Q', pert)
!       ENDIF
!       sint(:,:,:) = max(sint(:,:,:),0.0) ! no negative mixing ratios
!     ENDIF

!    ENDIF ! member


!   CALL cld_cpu('TDLBC')


   END SUBROUTINE INT_TDLBC



!-----------------------------------------------------------------------------
!
! ///////////////////// BEGIN \\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\ SUBROUTINE SET_BASE ////////////////////
!
! INT_TDLBC interpolates to get time dependent lateral boundary conditions
!-----------------------------------------------------------------------------
! Created by Mike Buban, July 5, 2011
!-----------------------------------------------------------------------------

    SUBROUTINE SET_BASE(gd,varname,ndgvarew,ndgvarns,nx,ny,nz,time)
       
  USE GRID_MODULE

  TYPE(GRID) :: gd
  
  integer nx,ny,nz
  integer i,j,k, time
  integer m,l
  real dt
  real :: ndgvarew(ny,nz,2),ndgvarns(nx,nz,2)
  real, pointer :: u1d(:)
  real, pointer :: v1d(:)
  real, pointer :: tz(:)
  real, pointer :: qz(:)

  real, save    :: t3t(150,150,41),q3t(150,150,41)
  real, save    :: u3t(150,150,41),v3t(150,150,41)

  character(LEN=10)  :: varname 

  CALL GET_VARIABLE(gd,'UINIT', u1d)
  CALL GET_VARIABLE(gd,'VINIT', v1d)
  CALL GET_VARIABLE(gd,'THINIT',tz)
  CALL GET_VARIABLE(gd,'QVINIT',qz)
  CALL GET_VARIABLE(gd,'DT',dt)


  if(time.eq.dt) then  

  open(13, file='ideal_def.int',status='unknown')
      
 4022 format (150(f10.4,1x))
       

#ifdef MPI
      read(13,4022)(((t3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      read(13,4022)(((q3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      read(13,4022)(((u3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
      read(13,4022)(((v3t(m,l,k),m=1,nxend),l=1,nyend),k=1,nzend)
#else      
      read(13,4022)(((t3t(m,l,k),m=1,nx),l=1,ny),k=1,nz)
      read(13,4022)(((q3t(m,l,k),m=1,nx),l=1,ny),k=1,nz)
      read(13,4022)(((u3t(m,l,k),m=1,nx),l=1,ny),k=1,nz)
      read(13,4022)(((v3t(m,l,k),m=1,nx),l=1,ny),k=1,nz)
#endif
            
      close(13)
      
!....msb...make format equal to nx

   endif 
     

  do k=1,nz

#ifdef MPI 

  if(varname.eq.'TH') then
  ndgvarew(:,k,1)=tz(k)
  ndgvarew(:,k,2)=t3t(nxend,:,k)
  ndgvarns(:,k,1)=tz(k)
  ndgvarns(:,k,2)=tz(k)
!  print*,ndgvarew(1,k,1),ndgvarew(1,k,2),t3t(nx,4,k),k,time
  elseif(varname.eq.'QV') then
  ndgvarew(:,k,1)=qz(k)
  ndgvarew(:,k,2)=q3t(nxend,:,k)
  ndgvarns(:,k,1)=qz(k)
  ndgvarns(:,k,2)=qz(k)
  elseif(varname.eq.'U') then  
  ndgvarew(:,k,1)=u1d(k)
  ndgvarew(:,k,2)=u3t(nxend,:,k)
  ndgvarns(:,k,1)=u1d(k)
  ndgvarns(:,k,2)=u1d(k)
  elseif(varname.eq.'V') then
  ndgvarew(:,k,1)=v1d(k)
  ndgvarew(:,k,2)=v3t(nxend,:,k)
  ndgvarns(:,k,1)=v1d(k)
  ndgvarns(:,k,2)=v1d(k)
  end if

#else

  if(varname.eq.'TH') then
  ndgvarew(:,k,1)=tz(k)
  ndgvarew(:,k,2)=t3t(nx,:,k)
  ndgvarns(:,k,1)=tz(k)
  ndgvarns(:,k,2)=tz(k)
!  print*,ndgvarew(1,k,1),ndgvarew(1,k,2),t3t(nx,4,k),k,time
  elseif(varname.eq.'QV') then
  ndgvarew(:,k,1)=qz(k)
  ndgvarew(:,k,2)=q3t(nx,:,k)
  ndgvarns(:,k,1)=qz(k)
  ndgvarns(:,k,2)=qz(k)
  elseif(varname.eq.'U') then  
  ndgvarew(:,k,1)=u1d(k)
  ndgvarew(:,k,2)=u3t(nx,:,k)
  ndgvarns(:,k,1)=u1d(k)
  ndgvarns(:,k,2)=u1d(k)
  elseif(varname.eq.'V') then
  ndgvarew(:,k,1)=v1d(k)
  ndgvarew(:,k,2)=v3t(nx,:,k)
  ndgvarns(:,k,1)=v1d(k)
  ndgvarns(:,k,2)=v1d(k)
  end if

#endif

  end do  

   END SUBROUTINE SET_BASE











!-----------------------------------------------------------------------------
!
! ///////////////////// BEGIN \\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\ SUBROUTINE INT_TDLBC_MOV ////////////////////
!
! INT_TDLBC interpolates to get time dependent lateral boundary conditions
!-----------------------------------------------------------------------------
! Created by Mike Buban, August 22, 2011
!-----------------------------------------------------------------------------
!
!  (11/10/11) modified by C. Ziegler and M. Buban
!

!
!.... CLZ(3-5-13): additional heavy debugging
!


    SUBROUTINE INT_TDLBC_MOV(gd,varname,sintew,sintns,    &
    time,nx,ny,nz,kmax,ugrid,vgrid,x_sw_loc,y_sw_loc)
       
  USE GRID_MODULE
  USE INIT_MODULE, only : nmtimes,nmint,nxmeso,nymeso,nzmeso,    &
                          athet,abmeso,au,av,aqrat, &
                          dxmeso,dymeso
  USE COMMASMPI_MODULE
  
  implicit none


  TYPE(GRID) :: gd

    integer            :: i,j,k,m,l,n
    integer            :: time
    integer  :: nin(2)
    integer  :: ninv(2)
    integer            :: nx,ny,nz
    integer            :: ibeg,iend
    integer            :: kmax
    integer            :: n2

    real, pointer :: dzc(:)
    real, pointer :: dze(:)
    real, pointer :: zcdx(:,:)
    real, pointer :: zedx(:,:)

    real               :: x_sw_loc,y_sw_loc
    real               :: ugrid,vgrid

    real ufl(11)
    real vfl(11)

!.... CLZ (3/4/13): nzend+1 should be nxend+1????
!    real               :: athet(nzend+1,nyend+1,nzend+1,2)
!    real               :: aqrat(nzend+1,nyend+1,nzend+1,2)
!    real               :: au(nzend+1,nyend+1,nzend+1,2)
!    real               :: av(nzend+1,nyend+1,nzend+1,2)


!
!.... must be saved for three RK steps per one dynamic time step.  Allocated in initsubs.F90.
!

!    real  ::  athet(nxend + 1,nyend + 1,nzend + 1,2)
!    real  ::  aqrat(nxend + 1,nyend + 1,nzend + 1,2)
!    real  ::  au(nxend + 1,nyend + 1,nzend + 1,2)
!    real  ::  av(nxend + 1,nyend + 1,nzend + 1,2)


!
!.... CLZ (3/5/13): dimension real, allocatable to avoid array stack problem
!

    real, allocatable  :: tsint(:,:,:,:)
    real, allocatable  :: sint(:,:,:)
    real, allocatable  :: tsintew(:,:,:,:)
    real, allocatable  :: tsintns(:,:,:,:)

!
!.... CLZ (3/5/13): sintew and sintns were dimensioned in and passed from the calling routine advectrk
!

    real  :: sintew(ny,nz,2)
    real  :: sintns(nx,nz,2)


    real tm
    real dx
    real dy
    real dz
    real t1
    real t2
    real grid_xin(nxend+1)
    real grid_yin(nyend+1)
    real grid_zin(nzend+1)
    real f1
    real f2
    real f3
    real f4
    real f5
    real f6
    real f7
    real f8
    real grid_x(1000)
    real grid_y(1000)
    real grid_z(1000)
    real, save  :: u(100),v(100),hgt(100),p(100),t(100),rh(100)



!    real, save  :: reflect(101,101,25,2),uvel(101,101,25,2)
!    real, save  :: vvel(101,101,25,2),wvel(101,101,25,2)



    real gxin,gyin,gx,gy,xmx,ymx
    real               :: xmin, ymin ! coordinates (m) of southwest corner of domain
!    real               :: hwlmax,hwlmin,vwlmax,vwlmin
    real,    allocatable :: sd(:,:,:)       ! [nx,ny,nz] standard deviation of smooth perturbations
    real,    pointer     :: xc(:)
    real,    pointer     :: xe(:)
    real,    pointer     :: yc(:)
    real,    pointer     :: ye(:)
    real,    pointer     :: zc(:)
    real,    pointer     :: ze(:)
    real  dt
    integer, pointer     :: member
    integer, save        :: nxr,nyr,nzr
    integer ix1,ix2,iy1,iy2,iz1,iz2
    integer, save      :: ninold,ninvold
    character(LEN=4)   :: fld
    character(LEN=10)  :: varname    
    character(LEN=70), save  :: infile(100)
    character(LEN=70), save  :: inrad(100)
    character(LEN=50)  :: bc_file

    TYPE(ATTRIBUTE), pointer :: bc_file_name
!    TYPE(ATTRIBUTE), pointer :: member

!   TYPE(VARIABLE)     :: pertu,pertv
!   TYPE(VARIABLE)     :: pertw,pertt,pertq

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze

      integer, parameter :: ntot = 50
      double precision  mpitotin(ntot), mpitotout(ntot)

#endif  

!-----------------------------------------------------------------------







!
!.... begin executable code
!








!   CALL cld_cpu('TDLBC')

!
!.... CLZ (3/5/13): dimension real, allocatable to avoid array stack problem
!

!    allocate  ( athet(nxend + 1,nyend + 1,nzend + 1,2) )
!    allocate  ( aqrat(nxend + 1,nyend + 1,nzend + 1,2) )
!    allocate  ( au(nxend + 1,nyend + 1,nzend + 1,2) )
!    allocate  ( av(nxend + 1,nyend + 1,nzend + 1,2) )

    allocate  ( tsint(nx,ny,nz,2) )
    allocate  ( sint(nx,ny,nz) )
    allocate  ( tsintew(ny,nz,2,2) )
    allocate  ( tsintns(nx,nz,2,2) )


    CALL GET_VARIABLE(gd, 'DX', dx)
    CALL GET_VARIABLE(gd, 'DY', dy)
    CALL GET_VARIABLE(gd, 'DZ', dz)


!    CALL GET_VARIABLE(gd, 'ZCDX',   zcdx)
!    CALL GET_VARIABLE(gd, 'ZEDX',   zedx)
    CALL GET_VARIABLE(gd,'DZC',dzc)
    CALL GET_VARIABLE(gd,'DZE',dze)
    CALL GET_VARIABLE(gd,'XC',xc)
    CALL GET_VARIABLE(gd,'XE',xe)
    CALL GET_VARIABLE(gd,'YC',yc)  
    CALL GET_VARIABLE(gd,'YE',ye)  
    CALL GET_VARIABLE(gd,'ZC',zc)  
    CALL GET_VARIABLE(gd,'ZE',ze)  
    CALL GET_VARIABLE(gd,'DT',dt)


!     print*,'in tdlbc', time, dt, dzc(1),dzc(2),varname
!      write(luno,*)'entered int_tdlbc_mov'

      if(time .eq. dt) then

      ninold = 0
      ninvold = 0

      endif


!
!.... CLZ (11/10/11): Do not time-interpolate if nmtimes = 1
!

      if (nmtimes .gt. 1) then
      nin(1)= (time/nmint) + 1
      nin(2)= nin(1) +1

      if(nin(2).gt.nmtimes) then
      nin(2)=nmtimes
      endif

      n2 = 2          ! n2 = 2 here because two time levels are needed


!
!.... CLZ (11/10/11): Do not time-interpolate if nmtimes = 1
!

      elseif (nmtimes .eq. 1) then

      nin(1) = 1
      nin(2) = 1
      n2 = 1          ! n2 = 1 here because only one time level is needed

      endif


!      print*,nin(1),nin(2),ninv(1),ninv(2),nmint,nmint,nmtimes,nmtimes,time



      if(nin(1) .gt. ninold) then

      ninold = nin(1)
!      write(6,*)'in INT_TDLBC_MOV',' time=',time,' my_rank=',my_rank

!
!.... CLZ (11/10/11): Do not time-interpolate if nmtimes = 1 (i.e., only if nmtimes > 1)
!

      do n = 1, n2

      do i = 1, nxmeso
      do j = 1, nymeso
      do k = 1, nzmeso
      
      athet(i,j,k,n) = abmeso(i,j,k,1,nin(n))
      aqrat(i,j,k,n) = abmeso(i,j,k,2,nin(n))
      au(i,j,k,n) = abmeso(i,j,k,3,nin(n)) - ugrid
      av(i,j,k,n) = abmeso(i,j,k,4,nin(n)) - vgrid
      
      enddo
      enddo
      enddo
      
      enddo

!
!.... endif(nin(1) .gt. ninold) then
!

      endif

    

!    write(6,*)'before tri-lin',' my_rank=',my_rank

! Tri-linear interpolation of input data to model grid
! msb 11/19/07
!



#ifdef MPI
  do i = nxbeg, nxend
#else  
  do i = 1, nx
#endif 
  grid_x(i) = x_sw_loc + (real(i) - 1.0) * dx
  enddo

#ifdef MPI
  do i = nybeg, nyend
#else  
  do i = 1, ny
#endif 
  grid_y(i) = y_sw_loc + (real(i) - 1.0) * dy
  enddo


!..... CLZ (3/27/13): mesoscale and model grid vertical levels are now coincident

#ifdef MPI
  do i = nzbeg, nzend - 1
#else  
  do i = 1, nz - 1
#endif
  grid_z(i) = zc(i)
  enddo


! dxmeso,dymeso,dzmeso => input grid spacing
! nxmeso,nymeso,nzmeso => input grid dimensions
!

!  dx=((real(nxmeso)-1)/(real(nx)-1))*dxmeso
!  dy=((real(nymeso)-1)/(real(ny)-1))*dymeso
!  dz=((real(nzmeso)-1)/(real(nz)-1))*dzmeso


  do i = 1, nxmeso + 1
  grid_xin(i) = (real(i) - 1.0) * dxmeso
  enddo

  do i = 1, nymeso + 1
  grid_yin(i) = (real(i) - 1.0) * dymeso
  enddo

!.... CLZ (3/27/13):  soundings, abmeso, and model grid now share same vertical levels

  do i = 1, nzmeso + 1
   grid_zin(i) = grid_z(i)
  enddo

!      write(6,*)'finished grid_zin'

   gxin = (real(nxmeso) - 1.0) * dxmeso
   gyin = (real(nymeso) - 1.0) * dymeso

#ifdef MPI   
   gx=(real(nxend)-1)*dx
   gy=(real(nyend)-1)*dy
   xmx=(gxin/gx)*(nxend-1) + 1
   ymx=(gyin/gy)*(nyend-1) + 1
#else   
   gx=(real(nx)-1)*dx
   gy=(real(ny)-1)*dy
   xmx=(gxin/gx)*(nx-1) + 1
   ymx=(gyin/gy)*(ny-1) + 1   
#endif

!      write(6,*)'finished ymx'

#ifdef MPI
   do j = jybeg, jyend
    do i = ixbeg, ixend

     do k = kzbeg, Min(kzend,nzmeso)
#else     
   do j=1,ny
    do i=1,nx

     do k=1,nz
#endif


     
#ifdef MPI
   ix1 = int( grid_x(i) / dxmeso ) + 1
  
  if(ix1 .gt. nxmeso) then
  ix1 = nxmeso
  endif
  ix2 = ix1 + 1

  iy1 = int( grid_y(j) / dymeso ) + 1

  if(iy1 .gt. nymeso) then
  iy1 = nymeso
  endif
  iy2 = iy1 + 1


   iz1 = k
   iz2 = iz1 + 1      
    
#else 
    
   ix1=int(grid_x(i)/dxmeso) +1
  
  if(ix1.gt.nxmeso) then
  ix1=nxmeso
  end if
  ix2=ix1+1

   iy1 = int(grid_y(j)/dymeso) +1
  
  if(iy1.gt.nymeso) then
  iy1=nymeso
  end if
  iy2= iy1+1


   iz1 = k
   iz2 = iz1 + 1
    
#endif

!     write(6,*) 'spatial interpolation',i,j,k,ix1,ix2,iy1,iy2,iz1,iz2

!.... CLZ (3/1/13): horizontal bilinear interpolation

      f1 = ( (grid_xin(ix2) - grid_x(i)) * (grid_yin(iy2) - grid_y(j)) )/(dxmeso * dymeso)
      f2 = ( (grid_x(i) - grid_xin(ix1)) * (grid_yin(iy2) - grid_y(j)) )/(dxmeso * dymeso)
      f3 = ( (grid_xin(ix2) - grid_x(i)) * (grid_y(j) - grid_yin(iy1)) )/(dxmeso * dymeso)
      f4 = ( (grid_x(i) - grid_xin(ix1)) * (grid_y(j) - grid_yin(iy1)) )/(dxmeso * dymeso)

!.... CLZ (11/10/11): Do not time-interpolate if nmtimes = 1 (i.e., n2 = 1)

     do n = 1, n2

     if(varname .eq. 'TH') then

#ifdef MPI

!      write(luno,*)'before i.eq.nxbeg with n=',n,'varname=',varname

     if(i .eq. nxbeg) then
     tsintew(j-jybeg+1,k-kzbeg+1,1,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + &
                                        f3*athet(ix1,iy2,k,n) + f4*athet(ix2,iy2,k,n)

!      if (j .eq. 30 .and. k .eq. 64) then
!      write(luno,*)'tsw_athet=',tsintew(j,k,1,n)
!      write(luno,*)'f1=',f1,' thet11=',athet(ix1,iy1,k,n)
!      write(luno,*)'f2=',f2,' thet21=',athet(ix2,iy1,k,n)
!      write(luno,*)'f3=',f3,' thet12=',athet(ix1,iy2,k,n)
!      write(luno,*)'f4=',f4,' thet22=',athet(ix2,iy2,k,n)
!      write(luno,*)'Px=',grid_x(i),' Py=',grid_y(j)
!      write(luno,*)'Q11x=',grid_xin(ix1),' Q11y=',grid_yin(iy1)
!      write(luno,*)'Q21x=',grid_xin(ix2),' Q21y=',grid_yin(iy1)
!      write(luno,*)'Q12x=',grid_xin(ix1),' Q12y=',grid_yin(iy2)
!      write(luno,*)'Q22x=',grid_xin(ix2),' Q22y=',grid_yin(iy2)
!     endif


!     print*,tsintew(j-jybeg+1,k-kzbeg+1,1,n),j-jybeg+1,k-kzbeg+1,n,my_rank,time     

     endif
#else          
     if(i .eq. 1) then
     tsintew(j,k,1,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + f3*athet(ix1,iy2,k,n) + &
                        f4*athet(ix2,iy2,k,n)
     endif
#endif


#ifdef MPI

!      write(6,*)'before i.eq.nxend with n=',n,'varname=',varname

     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) +  &
                                        f3*athet(ix1,iy2,k,n) + f4*athet(ix2,iy2,k,n)
     endif     
#else
     if(i.eq.nx) then   
     tsintew(j,k,2,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + f3*athet(ix1,iy2,k,n) +  &
                        f4*athet(ix2,iy2,k,n)
     endif
#endif
     

#ifdef MPI

!      write(6,*)'before j.eq.nybeg with n=',n,'varname=',varname

     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + &
                                        f3*athet(ix1,iy2,k,n) + f4*athet(ix2,iy2,k,n)
     endif
#else         
     if(j.eq.1) then
     tsintns(i,k,1,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + f3*athet(ix1,iy2,k,n) + f4*athet(ix2,iy2,k,n)
     endif
#endif

     
#ifdef MPI

!      write(6,*)'before j.eq.nyend with n=',n,'varname=',varname

     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + &
                                        f3*athet(ix1,iy2,k,n) + f4*athet(ix2,iy2,k,n)
     endif
#else
     if(j.eq.ny) then
     tsintns(i,k,2,n) = f1*athet(ix1,iy1,k,n) + f2*athet(ix2,iy1,k,n) + f3*athet(ix1,iy2,k,n) + f4*athet(ix2,iy2,k,n)
     endif
#endif



     elseif(varname .eq. 'QV') then

#ifdef MPI
     if(i.eq.nxbeg) then
     tsintew(j-jybeg+1,k-kzbeg+1,1,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + &
                                        f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#else          
     if(i.eq.1) then
     tsintew(j,k,1,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#endif

     
#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + &
                                        f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#else          
     if(i.eq.nx) then
     tsintew(j,k,2,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#endif
     

#ifdef MPI
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) +  &
                                        f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#else
     if(j.eq.1) then
     tsintns(i,k,1,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     end if
#endif
     
     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + &
                                        f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#else     
     if(j.eq.ny) then
     tsintns(i,k,2,n) = f1*aqrat(ix1,iy1,k,n) + f2*aqrat(ix2,iy1,k,n) + f3*aqrat(ix1,iy2,k,n) + f4*aqrat(ix2,iy2,k,n)
     endif
#endif



     elseif(varname .eq. 'U') then
!     elseif(varname .eq. 'U' .or. varname .eq. 'ULAG') then


#ifdef MPI
     if(i.eq.nxbeg) then
     tsintew(j-jybeg+1,k-kzbeg+1,1,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)


!      if (j .eq. 30 .and. k .eq. 64) then
!      write(luno,*)'tsw_au=',tsintew(j,k,1,n)
!      write(luno,*)'f1=',f1,' au11=',au(ix1,iy1,k,n)
!      write(luno,*)'f2=',f2,' au21=',au(ix2,iy1,k,n)
!      write(luno,*)'f3=',f3,' au12=',au(ix1,iy2,k,n)
!      write(luno,*)'f4=',f4,' au22=',au(ix2,iy2,k,n)
!      write(luno,*)'Px=',grid_x(i),' Py=',grid_y(j)
!      write(luno,*)'Q11x=',grid_xin(ix1),' Q11y=',grid_yin(iy1)
!      write(luno,*)'Q21x=',grid_xin(ix2),' Q21y=',grid_yin(iy1)
!      write(luno,*)'Q12x=',grid_xin(ix1),' Q12y=',grid_yin(iy2)
!      write(luno,*)'Q22x=',grid_xin(ix2),' Q22y=',grid_yin(iy2)
!     endif


     endif
#else     
     if(i.eq.1) then
     tsintew(j,k,1,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#endif

     
#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#else          
     if(i.eq.nx) then
     tsintew(j,k,2,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#endif


#ifdef MPI 
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#else
     if(j.eq.1) then
     tsintns(i,k,1,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#endif

     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#else     
     if(j.eq.ny) then
     tsintns(i,k,2,n) = f1*au(ix1,iy1,k,n) + f2*au(ix2,iy1,k,n) + f3*au(ix1,iy2,k,n) + f4*au(ix2,iy2,k,n)
     endif
#endif



     elseif(varname .eq. 'V') then
!     elseif(varname .eq. 'V' .or. varname .eq. 'VLAG') then

#ifdef MPI 
     if(i.eq.nxbeg) then 
     tsintew(j-jybeg+1,k-kzbeg+1,1,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#else
     if(i.eq.1) then
     tsintew(j,k,1,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#endif


#ifdef MPI
     if(i.eq.nxend) then
     tsintew(j-jybeg+1,k-kzbeg+1,2,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)                 
     endif
#else
     if(i.eq.nx) then
     tsintew(j,k,2,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#endif


#ifdef MPI
     if(j.eq.nybeg) then
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#else
     if(j.eq.1) then
     tsintns(i,k,1,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#endif

     
#ifdef MPI
     if(j.eq.nyend) then
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#else
     if(j.eq.ny) then
     tsintns(i,k,2,n) = f1*av(ix1,iy1,k,n) + f2*av(ix2,iy1,k,n) + f3*av(ix1,iy2,k,n) + f4*av(ix2,iy2,k,n)
     endif
#endif

!     print*,av(ix2,iy1,iz2,n),av(ix2,iy1,iz1,n),av(ix1,iy2,iz1,n),av(ix2,iy2,iz1,n),ix1,iy1,iz1
!     print*,tsint(i,j,k,n),varname,i,j,k,n



     else

#ifdef MPI
     tsintew(j-jybeg+1,k-kzbeg+1,1,n) = 0.0
     tsintew(j-jybeg+1,k-kzbeg+1,2,n) = 0.0
     tsintns(i-ixbeg+1,k-kzbeg+1,1,n) = 0.0
     tsintns(i-ixbeg+1,k-kzbeg+1,2,n) = 0.0
#else
     tsintew(j,k,1,n) = 0.0
     tsintew(j,k,2,n) = 0.0
     tsintns(i,k,1,n) = 0.0
     tsintns(i,k,2,n) = 0.0
#endif     
     endif


!      write(6,*)'done interpolating field (varname) =',varname



     enddo !times

!
!
!.... CLZ (3-5-13):  commented out this bad line for debug.  Model is now running!!!
!   CALL cld_cpu('I')
!
!

    enddo
      enddo
        enddo


!      write(6,*)'finished tri-linear interpolation to bounding planes'

!---------------------------------------------------
!
! Ended tri-linear interpolation
!
!--------------------------------------------------





!
! Now do the temporal interpolation between tsint(i,j,k,1)
! and tsint(i,j,k,2)
!
#ifdef MPI
   do j=jybeg,jyend
    do i=ixbeg,ixend
     do k=kzbeg,kzend
#else     
   do j=1,ny
    do i=1,nx
     do k=1,nz
#endif



!
!.... CLZ (11/10/11): Do not time-interpolate if nmtimes = 1
!

    if (nmtimes .gt. 1) then


    t1 = (nin(1) - 1) * nmint
    t2 = (nin(2) - 1) * nmint
    tm = real(time)    

#ifdef MPI
    if(varname.eq.'U'.or.varname.eq.'ULAG'.or.varname.eq.'V'.or.varname.eq.'VLAG') then

    if(i.eq.nxbeg) then
    sintew(j-jybeg+1,k-kzbeg+1,1)=(1.0/nmint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,1,1)  & 
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,1,2))
         
    end if

    if(i.eq.nxend) then
    sintew(j-jybeg+1,k-kzbeg+1,2)=(1.0/nmint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,2,2))
    end if
    
    if(j.eq.nybeg) then
    sintns(i-ixbeg+1,k-kzbeg+1,1)=(1.0/nmint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,1,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,1,2))
    end if

    if(j.eq.nyend) then
    sintns(i-ixbeg+1,k-kzbeg+1,2)=(1.0/nmint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,2,2))
    end if

    else

    if(i.eq.nxbeg) then
    sintew(j-jybeg+1,k-kzbeg+1,1)=(1.0/nmint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,1,1)  &
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,1,2))
!    print*,tsintew(j-jybeg+1,k-kzbeg+1,1,1),tsintew(j-jybeg+1,k-kzbeg+1,1,2),j-jybeg+1,k-kzbeg+1,n,my_rank,time    
    end if

    if(i.eq.nxend) then
    sintew(j-jybeg+1,k-kzbeg+1,2)=(1.0/nmint)*((t2-tm)*tsintew(j-jybeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintew(j-jybeg+1,k-kzbeg+1,2,2))
    end if

    if(j.eq.nybeg) then
    sintns(i-ixbeg+1,k-kzbeg+1,1)=(1.0/nmint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,1,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,1,2))
    end if
    
    if(j.eq.nyend) then
    sintns(i-ixbeg+1,k-kzbeg+1,2)=(1.0/nmint)*((t2-tm)*tsintns(i-ixbeg+1,k-kzbeg+1,2,1)  &
    + (tm-t1)*tsintns(i-ixbeg+1,k-kzbeg+1,2,2))
    end if

    end if
#else
    if(varname.eq.'U'.or.varname.eq.'ULAG'.or.varname.eq.'V'.or.varname.eq.'VLAG') then
    sintew(j,k,1)=(1.0/nmint)*((t2-tm)*tsintew(j,k,1,1) + (tm-t1)*tsintew(j,k,1,2))
    sintew(j,k,2)=(1.0/nmint)*((t2-tm)*tsintew(j,k,2,1) + (tm-t1)*tsintew(j,k,2,2))
    sintns(i,k,1)=(1.0/nmint)*((t2-tm)*tsintns(i,k,1,1) + (tm-t1)*tsintns(i,k,1,2))
    sintns(i,k,2)=(1.0/nmint)*((t2-tm)*tsintns(i,k,2,1) + (tm-t1)*tsintns(i,k,2,2))
    else
    sintew(j,k,1)=(1.0/nmint)*((t2-tm)*tsintew(j,k,1,1) + (tm-t1)*tsintew(j,k,1,2))
    sintew(j,k,2)=(1.0/nmint)*((t2-tm)*tsintew(j,k,2,1) + (tm-t1)*tsintew(j,k,2,2))
    sintns(i,k,1)=(1.0/nmint)*((t2-tm)*tsintns(i,k,1,1) + (tm-t1)*tsintns(i,k,1,2))
    sintns(i,k,2)=(1.0/nmint)*((t2-tm)*tsintns(i,k,2,1) + (tm-t1)*tsintns(i,k,2,2))
    end if
#endif


!
!.... CLZ (11/10/11): Do not time-interpolate if nmtimes = 1
!

    elseif (nmtimes .eq. 1) then


#ifdef MPI

    if(i .eq. nxbeg) then
    sintew(j-jybeg+1,k-kzbeg+1,1) = tsintew(j-jybeg+1,k-kzbeg+1,1,1)         
    endif

    if(i .eq. nxend) then
    sintew(j-jybeg+1,k-kzbeg+1,2) = tsintew(j-jybeg+1,k-kzbeg+1,2,1)
    endif
    
    if(j .eq. nybeg) then
    sintns(i-ixbeg+1,k-kzbeg+1,1) = tsintns(i-ixbeg+1,k-kzbeg+1,1,1)
    endif

    if(j .eq. nyend) then
    sintns(i-ixbeg+1,k-kzbeg+1,2) = tsintns(i-ixbeg+1,k-kzbeg+1,2,1)
    endif

#else

    sintew(j,k,1) = tsintew(j,k,1,1)
    sintew(j,k,2) = tsintew(j,k,2,1)
    sintns(i,k,1) = tsintns(i,k,1,1)
    sintns(i,k,2) = tsintns(i,k,2,1)

#endif


    endif  ! endif (nmtimes .gt. 1) then


!      write(6,*)'finished time-interpolation loop of bounding planes'







!    if(varname.eq.'V'.and.j.eq.nyend) then
!    print*,sintns(i-ixbeg+1,k-kzbeg+1,2),tsintns(i-ixbeg+1,k-kzbeg+1,2,1),tsintns(i-ixbeg+1,k-kzbeg+1,2,2),i-ixbeg+1,k-kzbeg+1,n,my_rank,time
!    print*,i-ixbeg+1,k-kzbeg+1,t2,tm,t1
!    end if

    enddo
   enddo
  enddo


!    if(varname.eq.'TH') then
!    write(0,"('--------------ADVECT: after tdlbc --------------')")
!    do j=jybeg,jyend
!    do k=kzbeg,kzend
!    do n=1,2
!    IF( .not. ( sintew(j,k,n) .lt. 350. .and. sintew(j,k,n) .gt. 250. ) ) then
!    print*,sintew(j-jybeg+1,k-kzbeg+1,n),tsintew(j-jybeg+1,k-kzbeg+1,n,1),tsintew(j-jybeg+1,k-kzbeg+1,n,2)
!    print*,j-jybeg+1,k-kzbeg+1,tm,my_rank
!    end if
!    end do
!    end do
!    end do
!
!     end if




!   CALL cld_cpu('TDLBC')


!      write(6,*)'done with INT_TDLBC_MOV'



   END SUBROUTINE INT_TDLBC_MOV





