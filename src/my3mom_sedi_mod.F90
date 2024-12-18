module my3mom_sedi_mod

!================================================================================!
!  The following subroutines are used by the 'tmom' version of the multimoment   !
!  package.                                                                      !
!                                                                                !
!  Package version:  2.20.0     (internal bookkeeping)                           !
!  Last modified  :  2011-03-03                                                  !
!================================================================================!

   implicit none

  private
  public :: SEDI_ISGH_V33,countColumns_v33,blg5sedi

   contains

!=====================================================================================!
 SUBROUTINE SEDI_ISGH_V33(QX,NX,ZX,cat,Q,T,DE,gamfact,epsQ,epsN,epsZ,afx,bfx,cmx,dmx,    &
                          dtx,cx6,ALFxfix,Noxfix,LXP,npassx,ni,nk,VxMax,DxMax,DZ,SR,     &
                          scheme,ktop_sedi,SS_ON,ixe)

!-------------------------------------------------------------------------------------!
! Sedimentation subroutine for categories whose fall velocity equation is
! V(D) = gamfact * afx * D^bfx
!
!  ***  for my_main_full.ftn90  ***
!-------------------------------------------------------------------------------------!

  use my2mom_fncs_mod
  use my3mom_fncs_mod

  implicit none

! PASSING PARAMETERS:
  real, dimension(ni,nk), intent(inout) :: QX,NX,ZX,Q,T
  real, dimension(ni),    intent(inout) :: SR
  real, dimension(ni,nk), intent(in)    :: DE,DZ
  real,    intent(in)    :: dtx,epsQ,epsN,epsZ,cx6,VxMax,LXP
  real*8,  intent(in)    :: afx,bfx,cmx,dmx,ALFxfix,Noxfix,DxMax
  integer, intent(in)    :: npassx,ni,nk,scheme,cat !,ktop_sedi
  integer, dimension(ni), intent(in) :: ktop_sedi
  logical, intent(in) :: SS_ON 	! Added by DTD.  .false. to suppress size-sorting (i.e. set
  								! all moment fall speeds to mass-weighted fall speed)
  integer, intent(in) :: ixe

! LOCAL PARAMETERS:
  real, dimension(ni,nk) :: VVQ,VVN,VVZ,RHOQX,gamfact
  integer, dimension(nk) :: FLIM
  logical                :: slabHASmass,LOCALLIM,QxPresent
  integer, dimension(ni) :: activeColumn
  real                   :: VqMax,VnMax,Vzmax,cmxSP
  real*8                 :: ALFx,GX2,GX5,ckQx1,ckQx2,ckQx3,iLAMx,iLAMxB0,tmpdp1,tmpdp2,Dx
  integer                :: nnn,a,i,k,counter
  real*8, parameter      :: thrd  = 1.d0/3.d0

!-------------------------------------------------------------------------------------!

  cmxSP= sngl(cmx)
!
  !Determine for which slabs and columns sedimentation should be computes:
   call countColumns_v33(QX,ni,nk,epsQ,counter,activeColumn,slabHASmass,ktop_sedi,ixe)

   IF (slabHASmass) THEN

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    DO nnn= 1,npassx

       RHOQX= DE*QX
       VVQ= 0.;  VVN= 0.;  VVZ= 0.;  VqMax= 0.;  VnMax= 0.;  VzMax= 0.
       do a= 1,counter
         i=activeColumn(a)
!        do k= 1,nk
!         do k= ktop_sedi(i),nk
!		do k = 1,ktop_sedi(i)  ! Changed by DTD: go from bottom to top
         do k=1,nk-1

           if (scheme==1)                    then
              QxPresent = (QX(i,k)>epsQ)
              ALFx      = ALFxfix
              if (cat==2) then  ![ice]
                NX(i,k) = 5.*exp(0.304*(273.15-max(233.,T(i,k))))  !Cooper eqn.
              else if (cat>=3.and.cat<=5) then  ![snow, grpl, or hail]
                tmpdp1  = gammaDP(1.d0+ALFx)
                tmpdp2  = gammaDP(4.d0+ALFx)
                NX(i,k) = (Noxfix*tmpdp1)**(3./(4.+ALFx))*(tmpdp1/tmpdp2*DE(i,k)*     &
                           QX(i,k)/cmx)**((1.+ALFx)/(4.+ALFx))
              endif
           else if (scheme==2.or.scheme==3)  then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN)
              if (QxPresent) then
                 Dx  = (dble(DE(i,k)*QX(i,k)/NX(i,k))/cmx)**thrd
                 ALFx= diagAlpha_v33(Dx,cat)
              endif
           else if (scheme==4)               then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN .and. ZX(i,k)>epsZ)
              if (QxPresent) ALFx= solveAlpha_v33(QX(i,k),NX(i,k),ZX(i,k),cmxSP,DE(i,k))
           endif

           if (QxPresent) then
              GX2      = 1.d0/gammaDP(1.d0+ALFx)
              GX5      = gammaDP(1.d0+ALFx+dmx)
              ckQx1    = afx*gammaDP(1.d0+ALFx+dmx+bfx)/GX5
              ckQx2    = afx*gammaDP(1.d0+ALFx+bfx)*GX2
              ckQx3    = afx*gammaDP(7.d0+ALFx+bfx)/gammaDP(7.d0+ALFx)
              iLAMx    = (dble(QX(i,k)*DE(i,k)/NX(i,k))/(cmx*GX5*GX2))**thrd
              iLAMxB0  = iLAMx**bfx
              tmpdp1   = -gamfact(i,k)*iLAMxB0
              VVQ(i,k) = tmpdp1*ckQx1;   VqMax= max(VxMAX,-VVQ(i,k))
              if(SS_ON) then
                VVN(i,k) = tmpdp1*ckQx2;   VnMax= max(VxMAX,-VVN(i,k))
                VVZ(i,k) = tmpdp1*ckQx3;   VzMax= max(VxMAX,-VVZ(i,k))
              else
                VVN(i,k) = VVQ(i,k);   VnMax= max(VxMAX,-VVN(i,k))
                VVZ(i,k) = VVQ(i,k);   VzMax= max(VxMAX,-VVZ(i,k))
              end if
           endif

         enddo  !k-loop
       enddo    !i-loop
       locallim= (nnn==1)
       call blg5sedi(RHOQX,DZ,VVQ,   nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)
       if (scheme >1)  &
          call blg5sedi(NX,DZ,VVN,   nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)
      if (scheme==4)  &
          call blg5sedi(ZX,DZ,VVZ,   nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)

       DO k = 1,nk-1             ! Changed by WYH.
         DO i = 1,ixe
           QX(i,k)= RHOQX(i,k)/DE(i,k)
         END DO
       END DO

    ! Prevent levels with zero N and nonzero Q and size-limiter:
       IF (scheme>1) THEN
       do a= 1,counter
         i=activeColumn(a)
!        do k= 1,nk
!        do k= ktop_sedi(i),nk
		     !do k = 1,ktop_sedi(i)  ! Changed by DTD: Go from bottom to top
		     do k=1,nk-1
           if (scheme==2.or.scheme==3)  then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN)
           elseif (scheme==4)           then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN .and. ZX(i,k)>epsZ)
           endif
           if (.not. QxPresent) then
              Q(i,k) = Q(i,k)+QX(i,k)
              T(i,k) = T(i,k) - LXP*QX(i,k)   !LCP for rain; LSP for i,s,g,h
              QX(i,k)= 0.;  NX(i,k)= 0.  
              if (scheme == 4) ZX(i,k)= 0.		! Changed by WYH
           else  ! size limiter:
              Dx     = (dble(DE(i,k)*QX(i,k)/NX(i,k))/cmx)**thrd
              tmpdp1 = sngl(max(Dx,DxMAX)/DxMAX)
              NX(i,k)= NX(i,k)*tmpdp1*tmpdp1*tmpdp1
           endif
         enddo
       enddo
       ENDIF !(if scheme>1)

!       SR(:)= SR(:) - cx6*VVQ(:,nk)*DE(:,nk)*QX(:,nk)
       DO i = 1,ixe     ! Changed by WYH
         SR(i)= SR(i) - cx6*VVQ(i,2)*DE(i,2)*QX(i,2)
       END DO

    ENDDO  !nnn-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

   ENDIF  !slabHASmass

 END SUBROUTINE SEDI_ISGH_V33

!=====================================================================================!

 SUBROUTINE countColumns_v33(QX,ni,nk,minQX,counter,activeColumn,slabHASmass,ktop_sedi,ixe)

! Searches the hydrometeor array QX(ni,nk) for non-zero (>minQX) values.
! Returns slabHASmass=TRUE if there is a single non-zero value in the array
! and returns the array if i-indices (activeColumn) for the columns (i)
! which contain at least one non-zero value, as well as the number of such
! columns (counter).

  implicit none

 ! PASSING PARAMETERS:
  integer, intent(in)                   :: ni,nk !,ktop_sedi
  integer, intent(out)                  :: counter
  integer, dimension(ni), intent(in)    :: ktop_sedi
  integer, dimension(ni), intent(out)   :: activeColumn
  real,    dimension(ni,nk), intent(in) :: QX
  real,    intent(in)                   :: minQX
  logical, intent(out)                  :: slabHASmass
  integer, intent(in)                   :: ixe

! LOCAL PARAMETERS:
   integer                              :: i !,k
   integer, dimension(ni)               :: k

   k=0; counter=0; activeColumn=0; slabHASmass=.false.			! DTD: changed back to starting at
   																! k=0 (vertical index order reversed in ARPS)
!   k=ktop_sedi-1; counter=0; activeColumn=0; slabHASmass=.false.
   do i=1,ixe                 ! change to ni-1 for ARPS

      do
         k(i)=k(i)+1
         if (QX(i,k(i))>minQX) then
            counter=counter+1
            activeColumn(counter)=i
            slabHASmass=.true.
            k(i)=0
            exit
         else
            if (k(i)==nk-1) then  ! change to nk-1 for ARPS
			      !if(k(i) == ktop_sedi(i)) then ! Actually go up to the top of sedimentation column now 
										! defined by ktop_sedi
               k(i)=0
               exit
            endif
         endif
      enddo

   enddo  !i-loop


 END SUBROUTINE countColumns_v33

!=====================================================================================!
!This subroutine is modified from S/P BLG.FTN

 SUBROUTINE blg5sedi (RO,DZ,WW,NK,DT,COMPLIM,WMAX,FLIM,counter,activeColumn,ktop_sedi)

 implicit none

!PASSING ARGUMENTS:

 integer, intent(in)                             :: NK,counter
 real, dimension(:,:), intent(inout)             :: ro
 real, dimension(:,:), intent(in)                :: dz,ww
 real, intent(in)                                :: dt,wmax
 integer, dimension(nk)                          :: FLIM
 integer, dimension(size(RO,dim=1)), intent(in)  :: activeColumn
 logical, intent(inout)                          :: COMPLIM
! integer, intent(in)                             :: ktop_sedi
 integer, dimension(size(RO,dim=1))              :: ktop_sedi

! Author
!          C. Girard and A. Plante (2001)
!
! Revisions
!
! 001      A. Glazer and A. Plante (Oct 2001)
!             - introduced complim and ind_limit in computation of precipitation
! 002      A. Plante (June 2003)
!             - IBM conversion, added ro_star. This imporved sedimentaion
!               in mixphase4 by 38%.
!             - change computation limit to insure reentrance in OpenMP
!             - blg2.ftn validates perfectly with blg.ftn of PHY 3.8.
! 003      J. Milbrandt (Nov 2004)
!             - Modified for use in Milbrandt-Yau microphysics scheme
!               + add condition for mass-in-column (activeColumn); for all i-loops, the
!                 line 'do i=1,ni' was replaced with 'do a=1,counter' & 'i=activeColumn(a)'
!               + remove RT calculation (and pass)
!               + hard-wired various options for multimoment.ftn90 (removed unnecessary IF statements)
! 004      J. Milbrandt (Jan 2007)
!             - corrected 'idzmin' initial value (for use in activeColumn i-loops)
!             - added option to exclude upper levels (ktop_sedi)
! 005      J. Milbrandt and R. McTaggart-Cowan (March 2007)
!             - removed i-dependency; cleaned up arrays, etc.  RO, DZ, WW are (ni,nk) arrays but there is
!               only one outer i (a)-loop (which could be moved outside of this subroutine, leaving it
!               as a single column subroutine).
! 006      J. Milbrandt (Nov 2008)
!             - changed 'ktop_sedi' from scalar to i-array
! 007      D. Dawson and Yunheng Wang (CAPS, Jul 2011)
!          Adapted for ARPS system (imported changes from original blg3 subroutine in previous
!          ARPS interface)
!
! Object
!
!
!   Version 1.0
!
!   CALCULATES
!                sedimentation of a quantity of matter ro(kg/m3)
!                falling with negative downward velocity ww(m/s)
!
!   ACCORDING TO
!                the BOX-LAGRANGIAN SCHEME
!                         (Kato,1995:JMSJ,73,241-245)
!             OR
!                the ADJACENT-LAGRANGIAN SCHEME
!                Girard and Plante
!
!   PLUS
!                a conservative two-grid-length filter
!                in order to control noise if needed
!
!Arguments
!
!          - Input/Output -
!
! RO       density in cell before and after sedimentation.
! COMPLIM  logical switch to force recomputing the table FLIM
!
!          - Input -
!
! DZ       cell thickness
! WW       vertical velocity (negative downward)
! NI       number of profiles
! NK       number of grid point in profiles
! DT       timestep
! KF       number of filtering pass
! FLIM     index of hier grid point where to look for mass in sedimentation.
! WMAX     maximum vertical velocity (positive or negative value).
!          This is used to save computation time. For each
!          level, the index of the heighest level from which
!          mass can be received will be computed if COMPLIM is .true.
! IDIR           direction of index:
!                idir=0 if index one is at top of model
!                idir=1 if index one is at surface of model
! ktop_sedi      uppermost level below which sedimentation is computed
! COUNTER        number of columns in (ni,nk) slab with non-zero RO
! ACTIVECOLUMN   array of i-indices with non-zero columns in RO(ni,nk)


! LOCAL VARIABLES AND PARAMETERS:

 integer                :: i,k,l,km1,kp1,ks,kw,a !,ke
 integer, dimension(size(RO,dim=1)) :: ke
 real, dimension(nk)    :: vp,zt,zb,dzi,ro_star
 real, dimension(0:nk)  :: zz
 real                   :: zmax,tempo
 integer                :: idzmin
 real, parameter        :: epsilon = 1.e-2

!---------------------------------------------------------------------------
!     Set parameters related to index direction.

!       ks=1
!       ke=nk
!       kw=1
!       if(idir.eq.0)then
!          ks=nk
!          ke=1
!          kw=-1
!       endif

  !For nk=bottom:  (hard-wired to remove pass of 'idir' parameter)
  !   ks=  nk
  !   ke=  ktop_sedi   !ke=1 (old)
  !   kw= -1
     
  ! For ARPS, k=2 is bottom
  ks=1
  !ke=ktop_sedi
  ke=nk-1
  kw=1

    !------------------------------------------------------------------------
     DO a= 1,counter   !i=1,ni
        i=activeColumn(a)

!     Compute cell height and final position of the top (zt) and bottom (zb)
!     of moving boxes:

        zz(ks)=0.
        do k=ks,ke(i),kw
            zz (k+kw)=zz(k)+dz(i,k)
            dzi(k)=1./dz(i,k)
            vp (k)=0.
        enddo

        do k=ks,ke(i),kw
            zb(k)=zz(k)+ww(i,k)*dt
        enddo

!     Scheme='Girard' (not 'Kato')
        zt(ke(i))=zb(ke(i))+dz(i,ke(i))
        do k=ks,ke(i)-kw,kw
               zb(k)=min(zb(k+kw)-epsilon*dz(i,k),zz(k)+ww(i,k)*dt)
               zt(k)=zb(k+kw)
        enddo

        do k=ks,ke(i),kw    !k=1,nk
            ro_star(k)= ro(i,k)*dz(i,k)/(zt(k)-zb(k))
        enddo

!      Compute limit index where to look for mass:
        if (complim) then
           zmax=abs(wmax*dt)
           do l=ks,ke(i),kw
              flim(l)=l
              do k=l,ke(i),kw
                  if (zmax.ge.zz(k)-zz(l+kw)) flim(l)=k
              enddo
           enddo
        endif

        do l=ks,ke(i),kw
          do k=l,flim(l),kw
               vp(l)= vp(l) + ro_star(k)*max( 0.,min(zz(l+kw),zt(k)) &
                         - max(zz(l),zb(k)) )
          enddo
        enddo

        do k=ks,ke(i),kw
            ro(i,k)=vp(k)*dzi(k)
        enddo


     ENDDO !i-loop
    !------------------------------------------------------------------------
 RETURN

 END SUBROUTINE blg5sedi
!=====================================================================================!

end module my3mom_sedi_mod
