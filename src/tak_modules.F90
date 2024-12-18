
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module kind_parameters

        implicit none

          public

          integer,parameter :: int4 =kind(0)
          integer,parameter :: real4=kind(0.0)
          integer,parameter :: real8=kind(0.0d0)
!          integer,parameter :: real8=real4 ! kind(0.0d0)
  
      end module kind_parameters

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module takmpi

        use kind_parameters, only: i4=>int4

        implicit none

#if   defined(MPI)
        include 'mpif.h'

#else
          integer(kind=i4),parameter :: MPI_COMM_WORLD       = 0
          integer(kind=i4),parameter :: MPI_STATUS_SIZE      = 1 ! set one to avoid zero-sized array declaration
          integer(kind=i4),parameter :: MPI_DOUBLE_PRECISION = 0
          integer(kind=i4),parameter :: MPI_INTEGER          = 0
          integer(kind=i4),parameter :: MPI_PACKED           = 0
          integer(kind=i4),parameter :: MPI_SUM              = 0
          integer(kind=i4),parameter :: MPI_PROC_NULL        = 0
#endif

          integer(kind=i4) :: rank=0,npcs,root=0,error
          integer(kind=i4) :: status(MPI_STATUS_SIZE)

          integer(kind=i4) :: xnum,ynum
          integer(kind=i4) :: is,ie,js,je,ie1,je1
          integer(kind=i4) :: xs,xe,ys,ye,xdn,xup,ydn,yup

          integer(kind=i4),allocatable :: xpos(:),ypos(:)
          integer(kind=i4),allocatable :: ixs (:),ixe (:),jys(:),jye(:)

      end module takmpi

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module parameters

        use kind_parameters, only: i4=>int4,r8=>real8

        implicit none

          public

          integer(kind=i4)  :: im, jm, km
          integer(kind=i4),parameter  :: kmax = 410
!         integer(kind=i4),parameter :: im= 33,jm= 33,km=61  !12km
!45km     integer(kind=i4),parameter :: im=111,jm=111,km=61  !45km
!90km     integer(kind=i4),parameter :: im=227,jm=227,km=61  !90km
!          integer(kind=i4),parameter :: nx=im+1,ny=jm+1

          real   (kind=r8),parameter :: nocell=0.00d0
          integer(kind=i4)  :: ijk ! =im*jm

          integer(kind=i4),parameter :: stdin=5,stdout=6,stderr=7

      end module parameters

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module comm1

        use kind_parameters, only: r8=>real8, r4=>real4
        use parameters,      only: im,km,kmax
        use takmpi,             only: is,ie,js,je,ie1,je1

        implicit none

          public

          real(kind=r8) :: rhofn(45,2),amas(45,2)
          real(kind=r8) :: rhof (   2)

          real(kind=r8), allocatable :: b  (:,:) ! (im,km)
          real(kind=r8), allocatable :: aww(:,:) ! (km,km)
          real(kind=r8) :: rho0(kmax)

          real(kind=r4),allocatable   ::                                       &
     &                   thpr(:,:,:),qv  (:,:,:),qw   (:,:,:) ,etaft(:,:,:)    &
     &                    ,sigw(:,:,:),sigi(:,:,:),cn   (:,:,:)                &
     &                    ,ci  (:,:,:),cng (:,:,:)                              
!                       ,qww  (:,:,:)                &
!     &                    ,qwf (:,:,:),qwi (:,:,:),pipr (:,:,:)                &
!     &                    ,ak  (:,:,:),fx  (:,:,:),fy   (:,:,:)                &
!     &                    ,fz  (:,:,:),f   (:,:,:)

          contains

            subroutine acomm1

            implicit none
            
            allocate(   b(im,km) )
            allocate( aww(km,km) )

            allocate(thpr (is:ie ,js:je ,km)) ; thpr =0.0d0
            allocate(qv   (is:ie ,js:je ,km)) ; qv   =0.0d0
            allocate(qw   (is:ie ,js:je ,km)) ; qw   =0.0d0
            allocate(sigw (is:ie ,js:je ,km)) ; sigw =0.0d0
                allocate(sigi (is:ie ,js:je ,km)) ; sigi =0.0d0
                allocate(cn   (is:ie ,js:je ,km)) ; cn   =0.0d0
                allocate(ci   (is:ie ,js:je ,km)) ; ci   =0.0d0
                allocate(cng  (is:ie ,js:je ,km)) ; cng  =0.0d0
!                allocate(qww  (is:ie ,js:je ,km)) ; qww  =0.0d0
!                allocate(qwf  (is:ie ,js:je ,km)) ; qwf  =0.0d0
!                allocate(qwi  (is:ie ,js:je ,km)) ; qwi  =0.0d0
!                allocate(pipr (is:ie ,js:je ,km)) ; pipr =0.0d0
!                allocate(ak   (is:ie ,js:je ,km)) ; ak   =0.0d0
!                allocate(fx   (is:ie ,js:je ,km)) ; fx   =0.0d0
!                allocate(fy   (is:ie ,js:je ,km)) ; fy   =0.0d0
!                allocate(fz   (is:ie ,js:je ,km)) ; fz   =0.0d0
!                allocate(f    (is:ie ,js:je ,km)) ; f    =0.0d0
                allocate(etaft(is:ie ,js:je ,km)) ; etaft=0.0d0
            
            return
            
            end subroutine acomm1

      end module comm1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
      
      module comm2

        use kind_parameters, only: r8=>real8
        use parameters,      only: kmax

        use index_module,     only: nr=>ntakrd,np=>ntakpd,ni=>ntakid, nt=>ntakit

        implicit none

          public

          real(kind=r8) :: sr (np) = 0.d0,sx (np),sf (np),srf(np,2),sfr(np),sfl(np)    &
     &                    ,fw (np),vw (np),bew(np),cdw(np),gw(np)              &
     &                    ,szr(np), szf(np,2),szi(ni,nt)                        &
     &                    ,sxinv(np)
          real(kind=r8) :: sri(ni),sfi(ni),di (ni),rhoi1(ni),dichg(ni), dieff(ni),srfchg(np,2)

          real(kind=r8) :: sxf(np,2),ff(np,2),vf(np,2),gf(np,2),sxfmlr(np,2)           &
                          ,sxfinv(np,2),sxfmlrinv(np,2)
          integer       :: isxfmlr(np,2)
          real(kind=r8) :: sxi(ni,nt+1),rhoxi(ni,nt+1),sh(ni,nt+1),fi(ni,nt),vi(ni,nt),gi(ni,nt)
          real(kind=r8) :: sxiinv(ni,nt+1)
          real(kind=r8) :: sh00(ni,nt+1)

          real(kind=r8) :: t0   (kmax),t00  (kmax),p0   (kmax),h0  (kmax)              &
     &                    ,qv0  (kmax),rhoa (kmax),qw0  (kmax),esw0(kmax)              &
     &                    ,esi0 (kmax),qvsw0(kmax),qvsi0(kmax),tv0 (kmax)              &
     &                    ,th0  (kmax),th0pr(kmax),gamw (kmax),gami(kmax)              &
     &                    ,dfuw (kmax),dfui (kmax),geni (kmax),genf(kmax)              &
     &                    ,stick(kmax),tk1  (kmax)

      end module comm2
      
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


      module comm3

        use kind_parameters, only: r8=>real8
        use parameters,      only: km
!        use takmpi,             only: is,ie,js,je
        use index_module,     only: nr=>ntakrd,np=>ntakpd,ni=>ntakid, nt=>ntakit,nmore

        implicit none

          public
          
          real(kind=r8) :: rew   (45),xprw (45),xrimi (45),xfrzf (45)          &
     &                    ,xrimf (45),finx (45),workwp(45),workwm(45),workw(45) &
     &                    ,sum1  (45),sum2 (45),save  (45),savex (45),savemore(45,nmore),savexmore(45,nmore)    &
     &                    ,savesx(45),savec(45),betlp (45),betlm (45)          &
     &                    ,sum1a (45),sum2a(45)
          integer, parameter :: nxtramax = 10
          integer       :: nxtra
          real(kind=r8) :: finxtra(45,nxtramax),savextra(45,nxtramax),scalefacxtra(nxtramax)
          real(kind=r8) :: ftrni (ni)

          real(kind=r8) :: ref   (45, 2),xprf  (45, 2),xprfr (45, 2)           &
     &                    ,termfp(45, 2),termfm(45, 2), sqrtrefnummks(45,2)    &
     &                    ,sqrtreinummks(ni,5),sqrtrewnummks(45),              &
     &                     termf(45, 2)
          real(kind=r8) :: rei   (ni, 5),xpri  (ni, 5),xprir (ni, 5)           &
     &                    ,termip(ni, 5),termim(ni, 5),termkp(ni, 5+1)         &
     &                    ,termkm(ni, 5+1),srim  (ni, 5),frzix (ni, 5)         &
     &                    ,termk(5+1, ni),termi(ni, 5+1)
          real(kind=r8) :: dm    (nr,nr),dk    (nr,nr),xk    (nr,nr)
          real(kind=r8) :: pmxy(nr,nr,nr) ! break-up fragment array
          real(kind=r8) :: pxmy(nr,nr,nr) ! break-up fragment array
          real(kind=r8) :: coaleff(nr,nr),sepeff(nr,nr),dksep(nr,nr),dkall(nr,nr)
          real(kind=r8) :: c_k(nr,nr), ckbott(nr,nr) ! "courant number" and collision kernel in Bott 1998 
          real(kind=r8) :: c_ki(ni,ni,nt), ckibott(ni,ni,nt) ! "courant number" and collision kernel in Bott 1998 (for ice crystals)
          integer       :: ima(nr,nr)  ! Bott index of mass bin k corresponding to sum of masses i and j
          integer       :: imai(ni,ni,nt)  ! Bott index of mass bin k corresponding to sum of masses i and j
          integer       :: mplusy(nr,nr) ! index of mass bin just below the sum of bins m and y
          real(kind=r8) :: dmi   (ni,ni),dki   (ni,ni,nt),xki   (ni,ni,nt)           &
     &                    ,yii   (21,21)
          real(kind=r8) :: yxx   (45,45)
          real(kind=r8) :: yix   (21,45),yixf(21,45)

          real(kind=r8) :: dkwf  (nr,45, 2)
          real(kind=r8) :: dkwf1 (nr,45, 2)
          real(kind=r8) :: dkwfrar(nr,45, 2)
          real(kind=r8) :: dkwi  (nr,ni, 5)
          real(kind=r8) :: dkfi  (45,ni, 5,2)
          real(kind=r8) :: dkfi1 (45,ni, 5,2)
          real(kind=r8) :: dkff  (45,45, 2)
          real(kind=r8) :: dkff1 (45,45, 2)
          integer :: idoflgs(nmore)
          
!          data yii/  &
!      0.90,1.13,1.20,1.15,1.10,1.02,0.97,0.93,0.89,0.88,            &
!      0.86,0.83,0.81,0.79,0.77,0.73,0.71,0.70,0.68,0.67,0.65,   &
!      1.13,1.20,1.35,1.35,1.20,1.13,1.07,1.00,0.97,0.92,            &
!      0.89,0.87,0.85,0.83,0.81,0.79,0.77,0.75,0.73,0.70,0.69,   &
!      1.20,1.35,1.55,1.50,1.40,1.28,1.15,1.08,1.03,0.98,            &
!      0.94,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,   &
!      1.15,1.35,1.50,1.70,1.60,1.45,1.30,1.18,1.10,1.05,            &
!      1.00,0.97,0.92,0.89,0.87,0.85,0.83,0.81,0.80,0.78,0.76,   &
!      1.10,1.20,1.40,1.60,1.73,1.60,1.45,1.30,1.20,1.10,            &
!      1.07,1.01,0.97,0.93,0.90,0.88,0.87,0.85,0.83,0.81,0.80,   &
!      1.02,1.13,1.28,1.45,1.60,1.78,1.60,1.50,1.30,1.20,            &
!      1.10,1.07,1.03,0.98,0.95,0.91,0.89,0.87,0.85,0.84,0.83,   &
!      0.97,1.07,1.15,1.30,1.45,1.60,1.80,1.70,1.50,1.35,            &
!      1.23,1.13,1.08,1.02,1.00,0.97,0.93,0.90,0.88,0.86,0.85,   &
!      0.93,1.00,1.08,1.18,1.30,1.50,1.70,1.80,1.70,1.50,            &
!      1.38,1.25,1.14,1.08,1.05,1.00,0.97,0.93,0.90,0.88,0.87,   &
!      0.89,0.97,1.03,1.10,1.20,1.30,1.50,1.70,1.80,1.70,            &
!      1.52,1.40,1.25,1.14,1.10,1.05,1.02,0.98,0.95,0.92,0.90,   &
!      0.88,0.92,0.98,1.05,1.10,1.20,1.35,1.50,1.70,1.80,            &
!      1.70,1.58,1.45,1.30,1.20,1.12,1.07,1.02,0.99,0.96,0.93,   &
!      0.86,0.89,0.94,1.00,1.07,1.10,1.23,1.38,1.52,1.70,            &
!      1.80,1.75,1.60,1.48,1.37,1.23,1.15,1.08,1.04,1.00,0.97,   &
!      0.83,0.87,0.90,0.97,1.01,1.07,1.13,1.25,1.40,1.58,            &
!      1.75,1.80,1.70,1.60,1.48,1.35,1.23,1.14,1.09,1.03,1.00,   &
!      0.81,0.85,0.88,0.92,0.97,1.03,1.08,1.14,1.25,1.45,            &
!      1.60,1.70,1.80,1.75,1.60,1.47,1.35,1.23,1.14,1.10,1.05,   &
!      0.79,0.83,0.86,0.89,0.93,0.98,1.02,1.08,1.14,1.30,            &
!      1.48,1.60,1.75,1.80,1.75,1.60,1.50,1.38,1.25,1.18,1.10,   &
!      0.77,0.81,0.84,0.87,0.90,0.95,1.00,1.05,1.10,1.20,            &
!      1.37,1.48,1.60,1.75,1.80,1.75,1.62,1.50,1.40,1.28,1.18,   &
!      0.73,0.79,0.82,0.85,0.88,0.91,0.97,1.00,1.05,1.12,            &
!      1.23,1.35,1.47,1.60,1.75,1.80,1.73,1.60,1.50,1.40,1.30,   &
!      0.71,0.77,0.80,0.83,0.87,0.89,0.93,0.97,1.02,1.07,            &
!      1.15,1.23,1.35,1.50,1.62,1.73,1.80,1.75,1.60,1.52,1.40,   &
!      0.70,0.75,0.78,0.81,0.85,0.87,0.90,0.93,0.98,1.02,            &
!      1.08,1.14,1.23,1.38,1.50,1.60,1.75,1.80,1.75,1.62,1.52,   &
!      0.68,0.73,0.76,0.80,0.83,0.85,0.88,0.90,0.95,0.99,            &
!      1.04,1.09,1.14,1.25,1.40,1.50,1.60,1.75,1.80,1.75,1.63,   &
!      0.67,0.70,0.74,0.78,0.81,0.84,0.86,0.88,0.92,0.96,            &
!      1.00,1.03,1.10,1.18,1.28,1.40,1.52,1.62,1.75,1.80,1.77,   &
!      0.65,0.69,0.72,0.76,0.80,0.83,0.85,0.87,0.90,0.93,            &
!      0.97,1.00,1.05,1.10,1.18,1.30,1.40,1.52,1.63,1.77,1.80   /
      

!---

       end module comm3

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module comm4

        use kind_parameters, only: r8=>real8
        use parameters,      only: km,ijk
!        use takmpi,             only: is,ie,js,je
        use index_module,     only: nr=>ntakrd,np=>ntakpd,ni=>ntakid, nt=>ntakit

        implicit none

          public

          real(kind=r8) :: bcon1(45,2),bcon2(45,2)
          real(kind=r8) :: acon1(21,5),acon2(21,5)

          real(kind=r8),allocatable :: cmvf(:,:,:) ! (45,2,km)
          real(kind=r8),allocatable :: cmvi(:,:,:) ! (21,5,km)

          real(kind=r8) :: sum3(nr),sum4(nr),tbk(nr)
          real(kind=r8) :: qbk (nr,nr)

!          real(kind=r8) :: a(ijk)
!          real(kind=r8),allocatable                                            &
!     &                  :: c(:,:)
!---
          contains
            subroutine acomm4

             use index_module,     only: ntakrd,ntakpd,ntakid, ntakit

              implicit none
              
              allocate ( cmvf(ntakpd,2,km) )
              allocate ( cmvi(ntakid,ntakit,km) )
             
             return
            end subroutine acomm4

      end module comm4



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module comm5

        use kind_parameters, only: r8=>real8, r4=>real4
        use parameters,      only: km,kmax
        use takmpi,             only: is,ie,js,je,ie1,je1

        implicit none

          public

!          real(kind=r4),allocatable                                            &
!     &                  :: up    (:,:,:),vp   (:,:,:),wp (:,:,:)               &
!     &                    ,etaft1(:,:,:),qv1  (:,:,:),qvp(:,:,:)               &
!     &                    ,thpr1 (:,:,:),thprp(:,:,:)

!          real(kind=r4),allocatable                                            &
!     &                  :: aheat(:,:,:),thpr1 (:,:,:)
          real(kind=r8) :: dmas(kmax),pai0(kmax),aht0(kmax),aht1(kmax),aht2(kmax)        &
     &                                      ,aht3(kmax),ahtx(kmax),aht (kmax)
!---

          contains

            subroutine acomm5

              implicit none

!                allocate(up    (is:ie1,js:je ,km)) ; up    =0.0d0
!                allocate(vp    (is:ie ,js:je1,km)) ; vp    =0.0d0
!                allocate(wp    (is:ie ,js:je ,km)) ; wp    =0.0d0
!                allocate(etaft1(is:ie ,js:je ,km)) ; etaft1=0.0d0
!                allocate(qv1   (is:ie ,js:je ,km)) ; qv1   =0.0d0
!                allocate(qvp   (is:ie ,js:je ,km)) ; qvp   =0.0d0
!                 allocate(thpr1 (is:ie ,js:je ,km)) ; thpr1 =0.0d0
!                allocate(thprp (is:ie ,js:je ,km)) ; thprp =0.0d0
!                 allocate(aheat (is:ie ,js:je ,km)) ; aheat =0.0d0
 
            return
          end subroutine acomm5

      end module comm5


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module takcommon

        use kind_parameters, only: i4=>int4,r8=>real8
!        use comm1,           only: acomm1
!        use comm3,           only: acomm3
        use comm4,           only: acomm4
        use comm5,           only: acomm5

        implicit none

          public

          real   (kind=r8) :: eps1,eps2,d2t
          integer(kind=i4) :: ipst5,isw5

          real   (kind=r8) :: dt0,dx2,dy2,dz2                                  &
                             ,gamd,accel,abk,bbk,bkw,bka,bmu,bet,tn
          integer(kind=i4) :: jc,isw,isw1,lr41,lr100,ls250

          real   (kind=r8) :: dt,dx,dy,dz,dzsq,d2x,d2y,d2z,hdz                 &
     &                       ,h00,p00,cn0,pi,eps                          &
     &                       ,gg,power,alh,dmu,hlat,ck,ck1                     &
     &                       ,capth,ckapa,cpd,cpe,ceta,czeta,capw,capi         &
     &                       ,rd,rhof0,rhoi,rhoi0,rhoa0,rdry,rlsqw,rlsqi             &
     &                       ,sigs,sr0,sri0,sh0,dj,dji,dji0,sheat,ratio             &
     &                       ,sceta,scu,scrh,scqv,scqvs,scqw,sct,scth          &
     &                       ,scw,scx,sctemp,pastx
          integer(kind=i4) :: itau,munut,msec,ipst,ipst1,ipst2                 &
     &                       ,im1,jm1,km1,imph,jmph                            &
     &                       ,lmax,lfmax,iimax,kimax                     &
     &                       ,lmax1,lfmax1,iimax1,kimax1

          integer(kind=i4) :: kfmax = 2
          integer(kind=i4) :: ims,ime,jms,jme
          real   (kind=r8) :: mtim,mstr,mend
          integer          :: ccntype, shedsmall
          real   (kind=r8) :: cckm,ccne,ccne0,ccnefac,cnexp

!      integer nqsat
!      parameter (nqsat=1000001) ! (nqsat=20001)
!      real fqsat,fqsati
!      parameter (fqsat=0.002,fqsati=1./fqsat)

!      real tabqvs(nqsat),tabqis(nqsat),dtabqvs(nqsat),dtabqis(nqsat)

!---

          contains

            subroutine acommon

              implicit none

!                call acomm1
!                call acomm3
                call acomm4
                call acomm5

            return
          end subroutine acommon

      end module takcommon



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      module clark

        use kind_parameters, only: r8=>real8, r4=>real4
        use parameters,      only: km,kmax
        use takmpi,             only: is,ie,js,je

        implicit none

          public

          real(kind=r8) :: sigw0(kmax),sigi0(kmax)

          real(kind=r4),allocatable                                            &
     &                  :: sigwo (:,:,:),sigwp (:,:,:),sigw1(:,:,:)            &
     &                    ,advs  (:,:,:),difsig(:,:,:),difqv(:,:,:)            &
     &                    ,difth (:,:,:)
!---

          contains

            subroutine aclark

              implicit none

                allocate(sigwo (is:ie,js:je,km)) ; sigwo =0.0d0
                allocate(sigwp (is:ie,js:je,km)) ; sigwp =0.0d0
                allocate(sigw1 (is:ie,js:je,km)) ; sigw1 =0.0d0
                allocate(advs  (is:ie,js:je,km)) ; advs  =0.0d0
                allocate(difsig(is:ie,js:je,km)) ; difsig=0.0d0
                allocate(difqv (is:ie,js:je,km)) ; difqv =0.0d0
                allocate(difth (is:ie,js:je,km)) ; difth =0.0d0

            return
          end subroutine aclark

      end module clark


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------




!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


! MISCELLANEOUS SUBROUTINES


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      subroutine sruletest(sum,itop,save)
! simpson rule integration?
        use kind_parameters, only: i4=>int4,r8=>real8

        implicit none

          integer(kind=i4),intent( in) :: itop
          real   (kind=r8),intent( in) :: save(itop)
          real   (kind=r8),intent(out) :: sum

          integer(kind=i4) :: i,im2
!---
! test version with simple sum
            sum = save(1)
            do i=2,itop
              sum = sum+save(i)
            enddo

!---

        return
      end subroutine sruletest
!--------------------------------------------------------------------------

      subroutine srule(sum,itop,save)
! simpson rule integration?
        use kind_parameters, only: i4=>int4,r8=>real8

        implicit none

          integer(kind=i4),intent( in) :: itop
          real   (kind=r8),intent( in) :: save(itop)
          real   (kind=r8),intent(out) :: sum

          integer(kind=i4) :: i,im2
!---

          sum = save(1)
          if(itop>=4) then
            im2 = itop-2
            do i=2,im2,2
              sum = sum+4.0*save(i)+2.0*save(i+1)
            enddo
          endif

          if(itop==itop/2*2) then
            if(itop==2) then
              sum = (sum+save(itop))*0.5
            else
              sum = (sum-save(itop-1))/3.0+(save(itop-1)+save(itop))*0.5
            endif
          else
              sum = (sum+save(itop)+4.0*save(itop-1))/3.0
          endif
!---

        return
      end subroutine srule


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      subroutine Sumd(buf,n)

        use kind_parameters, only: i4=>int4,r8=>real8
        use takmpi,             only: MPI_COMM_WORLD                              &
     &                            ,MPI_DOUBLE_PRECISION,MPI_SUM                &
     &                            ,npcs,error

        implicit none

          integer(kind=i4),intent( in) :: n
          real   (kind=r8),intent(out) :: buf(n)
!---

#if   defined(MPI)

          real   (kind=r8) :: wk(n)
!---

          if(Npcs==1) return
!---

          call MPI_Allreduce(buf,wk,n,MPI_DOUBLE_PRECISION                     &
     &                               ,MPI_SUM,MPI_COMM_WORLD,Error)
          buf = wk
!---

#else

          if(Npcs==1) return
          buf = buf
!---
#endif

       return
      end subroutine Sumd

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!---

      subroutine Sumi(buf,n)

        use kind_parameters, only: i4=>int4
        use takmpi,             only: MPI_COMM_WORLD,MPI_INTEGER,MPI_SUM          &
     &                            ,npcs,error

        implicit none

          integer(kind=i4),intent( in) :: n
          integer(kind=i4),intent(out) :: buf(n)
!---

#if   defined(MPI)

          integer(kind=i4) :: wk(n)
!---

          if(Npcs==1) return
!---

          call MPI_Allreduce(buf,wk,n,MPI_INTEGER                              &
     &                               ,MPI_SUM,MPI_COMM_WORLD,Error)
          buf = wk
!---

#else

          if(Npcs==1) return
          buf = buf
!---
#endif

       return
      end subroutine Sumi

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

      subroutine const                                                 

        use parameters, only: im,jm,km
        use takcommon,     only: itau,munut,msec,isw,isw1                         &
     &                       ,im1,jm1,km1,lmax,lfmax,kfmax,iimax,kimax         &
     &                       ,lmax1,lfmax1,iimax1,kimax1                       &
     &                       ,dx,dy,dz,dzsq,dx2,dy2,dz2,d2x,d2y,d2z,hdz        &
     &                       ,h00,p00,cn0,pi,eps,gamd,gg,accel            &
     &                       ,power,alh,dmu,hlat,abk,bbk,bkw,bka               &
     &                       ,bmu,bet,tn,ck,ck1,capth,ckapa,cpd,cpe            &
     &                       ,ceta,czeta,capw,capi                             &
     &                       ,rd,rhof0,rhoi,rhoi0,rhoa0,rdry,rlsqw,rlsqi             &
     &                       ,sigs,sr0,sri0,sh0,dj,dji,dji0,sheat,ratio             &
     &                       ,sceta,scu,scw,scrh,scqv,scqvs,scqw,sct           &
     &                       ,scth,scx,sctemp,ccntype
        use takmpi,        only: rank,root
        use index_module,     only: ntakrd,ntakpd,ntakid, ntakit
        use micro_module,   only: takrhoi,taksri0,takkfmax
 
        implicit none
!---

!          if(Rank==Root) print '(1x,a,/)','ENTRY  TO  CONST'
!---
!          ccntype = 1 ! 1=maritime, 2=continental
          itau   = 0
          munut  = 0
          msec   = 0

          isw    = 0
          isw1   = 0

          im1    = im-1
          jm1    = jm-1
          km1    = km-1

          lmax   = ntakrd ! 33
          lfmax  = ntakpd ! 45
          kfmax  = takkfmax
          iimax  = ntakid ! 21
          kimax  = ntakit ! 5
          IF ( rank == root ) write(0,*) 'iimax,kimax,lfmax,lmax = ',iimax,kimax,lfmax,lmax
          lmax1  = lmax-1
          lfmax1 = lfmax-1
          iimax1 = iimax-1
          kimax1 = Max(1,kimax-1)

          dx     = 4.0e4
          dy     = 4.0e4
          dz     = 2.0e4
          dzsq   = dz*dz
          dx2    = dx*dx
          dy2    = dy*dy
          dz2    = dz*dz
          d2x    = 2.0*dx
          d2y    = 2.0*dy
          d2z    = 2.0*dz
          hdz    = dz*0.5

!          t001   = 300.0 ! surface temperature (init.F) (NOT theta!)
          h00    = 0.70  ! init.F: sfc. specific hum.
          p00    = 988.0 ! init.F: sfc. pressure
          cn0    = 1.0e-8 ! N0 in Fletcher function for active ice nuclei (seems to be units of cm**-3)
                          ! normally would be 1.0e-8??? Since factor is 1.e-5 for liter**-1, or 1.e-2 for m**-3
                          ! changed value from 1.e-7 (original) to 1.e-8

          pi     = 3.14159 ! Pi the irrational one
          eps    = 0.622   ! eps, ratio of Rv/Rd
          gamd   = -9.75e-5 ! lapse rate?
          gg     = 980.0    ! accel due to gravity (9.8 m/s**2)
          accel  = gg
          power  = gg/2.87e6 ! for hydrostatic eq. to get base state pressure, p0
                             ! 2.87e6 is Rd (dry gas constant)
          alh    = 595.0     ! not used?
          dmu    = 0.1346    ! kinematic viscosity (ceta/rhoa0)
          hlat   = 80.0     ! freezing latent heat, calories/gram = 334 kJ/kg
!          2260*80/334 = 541 cal/gram latent heat vapor.

          abk    = 62.3    ! "a" breakup constant from Srivastava 1971
          bbk    = 7.0     ! "b" breakup constant from Srivastava 1971
          bkw    = 1.47e-3
          bka    = 6.0e-5
          bmu    = 0.24       ! diffusion constant of water vapor in air
          bet    = 4.8e-7     ! beta 
!          tn     = 300.0
          ck     = 1.26e+6 ! coefficient in droplet fall speed (Tak. 1976b, eq. 10)
          ck1    = 0.4 ! not used in microphysics
          capth  = 298.0
          ckapa  = 0.286
          cpd    = 0.2399    ! calories per degree??  1 Joule = 0.2388 calories
          cpe    = (4.2e7)*cpd
          ceta   = 1.718e-4 ! dynamic viscosity?
          czeta  = 0.6      ! constant "beta" in Fletcher ice nuclei curve
          capw   = 595.0 ! Lv latent heat of vap. (calories/gram)
          capi   = 677.0 ! Ls latent heat of sub. (calories/gram)

          rd     = 0.0685
          rhof0  = 0.3
          rhoi   = takrhoi ! 0.9 ! 0.1
          rhoi0  = 0.1
          rhoa0  = 1.3e-3
          rdry   = 0.0685  ! Dry air gas constant
          rlsqw  = capw*capw*eps/(rdry*cpd)  ! Lv**2 * 0.622/(Rd*Cp) = Lv**2 * Rv/Cp
          rlsqi  = capi*capi*eps/(rdry*cpd)

          sigs   = 5.0
          sr0    = 2.0d-4
          sri0   = taksri0 ! sr0*10.0
          sh0    = 1.0d-3
          dj     = 4.329 ! = 3.0d0/log(2.d0) ! bin constant (hjo) for rain/graupel/hail
          dji    = 2.885 ! hjo for ice crystals
          dji0   = 2.885 ! original hjo for ice crystals (m = d**2)
          sheat  = 0.2399 ! ?
          ratio  = 0.622  ! eps

! print scale factors? (history times)
          sceta  = 1.0e6  ! not used? (seconds)
          scu    = 1.0    ! not used? (seconds)
          scw    = 1.0     ! not used? (seconds)
          scrh   = 1.0e3
          scqv   = 1.0e6
          scqvs  = 1.0e5
          scqw   = 1.0e6
          sct    = 100.0
          scth   = 1.0e2
          scx    = 1.0e4
          sctemp = 1.0e-4
!---
!          if(Rank==Root) print '(1x,a,/)','OUT  FROM  CONST'
!---

        return                                                            
      end subroutine const


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

#include "tak_modules_include.F90"

