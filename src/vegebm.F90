      subroutine vegebm(tsfc,eflx,fflx,tflx,qflx,tsoil,wsfc,wsoil, &
                       tcanp,wcanp,precip,radsw,radlw,veg,nx,ny,dt)
!-----------------------------------------------------------------------
!     subroutine slab calculates the ground temperature tendency
!     according to the residual of the surface energy budget
!     (blackadar, 1978b).                                  
!-----------------------------------------------------------------------
      
      USE PARAM_MODULE
      USE SPHYS_MODULE
      USE COMMASMPI_MODULE
      
      implicit none
      
      integer nx,ny
      
      real dt

      real tsfc (-ng+1:nx+ng,-ng+1:ny+ng)
      real wsfc (-ng+1:nx+ng,-ng+1:ny+ng)
      real tsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real wsoil(-ng+1:nx+ng,-ng+1:ny+ng)
      real tcanp(-ng+1:nx+ng,-ng+1:ny+ng)
      real wcanp(-ng+1:nx+ng,-ng+1:ny+ng)

      real radsw (-ng+1:nx+ng,-ng+1:ny+ng)
      real radlw (-ng+1:nx+ng,-ng+1:ny+ng)
      real tflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real qflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real eflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real fflx  (-ng+1:nx+ng,-ng+1:ny+ng)
      real precip(-ng+1:nx+ng,-ng+1:ny+ng)
      real veg   (-ng+1:nx+ng,-ng+1:ny+ng)      
      

      logical debug
      real :: d1p = 10.               ! wsfc depth in cm
      real :: d2p = 100.               ! wsoil depth in cm
!......try critical soil moisture value of 0.40
!      data wcrit   / 0.30 /
      real :: wcrit =  0.40 
      real :: cal2j =  4.1868E+04        ! convert cal cm**-2 to J m**-2
      real :: cm2m  = .01               ! convert cm to m
      
      real :: c1, c2, c2p
      real :: thrd, wmax
      integer :: i,j
      real :: rpr, rcg, rc2, alamg, alam2, d1g, d12, rcd
      real :: capg, cap2, rnet, restore, qs
      real :: dtgdt, dt2dt, fac, c1p, prec, dwdt
      
      integer :: ixe,jye
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      debug = .false.
!-----------------------------------------------------------------------
      if(debug) write(*,*) ' VEGEBM: begin '


      ixe = itile
      IF ( myproci == nproci ) ixe = ixend - ixbeg
      jye = jtile
      IF ( myprocj == nprocj ) jye = jyend - jybeg

! --- set constants
!
      c1 = 2. * SQRT(pii)
      c2 = 2. * pii / daysec
      c2p = 0.9 / daysec
      thrd = 1./3.
      wmax = 1.33 * wcrit

!     if(debug) then
!     write(*,*) ' wsfc '
!     do i=1,ixe,(nx-1)/10
!       write(*,'(i3,10f7.4)') i,(wsfc(j,i), j=1,ny-1,(ny-1)/10)
!     enddo

      do j = 1,jye
      do i = 1,ixe
      
!
!-----find slab thermal capacity (capg) (after Deardorff, 1978)
!
        rpr = MIN(1., 0.30 + 0.05 * wsfc(i,j) / (wsoil(i,j)+1.E-6) )
        rcg = 0.27 + wsfc(i,j)
        rc2 = 0.27 + wsoil(i,j)

! --- Sun and Wu (1992)
!       pterm = 0.5 - sign(0.5,0.07 - wsfc(i,j)) ! + for wsfc >= 0.07
!       alamg =      pterm * ( 0.00118 + 0.0062 * wsfc(i,j)) +  &
!              (1.-pterm) * ( 0.0065  + 0.0432 * wsfc(i,j)**thrd )

!       pterm = 0.5 - sign(0.5,0.07 - wsoil(i,j)) ! + for wsoil >= 0.07
!       alam2 =      pterm * ( 0.00118 + 0.0062 * wsoil(i,j)) +  &
!              (1.-pterm) * ( 0.0065  + 0.0432 * wsoil(i,j)**thrd )
! --- Deardorff (1978)
        alamg =  0.001 + 0.004 * wsfc(i,j)**thrd
        alam2 =  0.001 + 0.004 * wsoil(i,j)**thrd

        d1g = SQRT( daysec * alamg / rcg )
        d12 = SQRT( daysec * alam2 / rc2 )

        rcd = rpr * rcg * d1g + (1.-rpr) * rc2 * d12
        capg = c1 / (cal2j * rcd )

        cap2 = 1. / (cal2j * (0.27+wsoil(i,j)) * d2p )
!-----------------------------------------------------------------------
!
!-----compute the surface energy budget
!
!   EQN 36: Ha = -G = Hsg + LEg - Rsw + Rlwup - Rlwdown
!
!-----------------------------------------------------------------------
        rnet   = radlw(i,j) - radsw(i,j)

        restore = c2 * ( tsfc(i,j) - tsoil(i,j) )

        qs     = tflx(i,j) + Lav * qflx(i,j)

!        print*,Lav*qflx(i,j),tflx(i,j),qs,radlw(i,j),radsw(i,j)


!-----------------------------------------------------------------------
!
!-----compute the surface and soil temperature
!
!   EQN 8: Tsfc(n)=Tsfc(n-1)+dt*(-c1*Ha/(rho*cs*d1) 
!                           -c2*(Tsfc(n)-Tsoil(n))
!
!-----------------------------------------------------------------------

        dtgdt  = -1. * (rnet + qs) * capg - restore
!       dtgdt  = 0.
        tsfc(i,j)= tsfc(i,j) + dt * dtgdt

        dt2dt  = -1. * (rnet + qs) * cap2
!       dt2dt  = 0.
        tsoil(i,j)= tsoil(i,j) + dt * dt2dt
!

!        print*,tsfc(i,j),dt,dtgdt,rnet,qs,capg,restore,i,j

!-----surface soil moisture   (Eq. 42)
!
        precip(i,j)= 0.              ! eliminate precip feedback

        fac = MIN( .75, wsfc(i,j) / wmax)
        fac = MAX( .15, fac)
        c1p = 14. - 22.5 * (fac - .15)
        capg = -c1p / (rhow * cm2m * d1p)
        prec = (1. - veg(i,j) ) * precip(i,j)

        dwdt = capg * (qflx(i,j) + 0.1 * eflx(i,j) - prec) -  &
               c2p * (wsfc(i,j) - wsoil(i,j) )
        wsfc(i,j) = wsfc(i,j) + dt * dwdt
        wsfc(i,j) = MAX( 0.0, wsfc(i,j))
!
!-----deep soil moisture  (Eq. 43)
!
        prec = (1. - veg(i,j) ) * precip(i,j) 
        dwdt = - (qflx(i,j) + eflx(i,j) - prec) / (rhow * cm2m * d2p)
        wsoil(i,j) = wsoil(i,j) + dt * dwdt
        wsoil(i,j) = MAX( 0.0, wsoil(i,j))
!
!-----canapy moisture   (Eq. 29)
!
        prec = veg(i,j) * precip(i,j) 
        dwdt = prec - ( fflx(i,j) - eflx(i,j) )
        wcanp(i,j) = wcanp(i,j) + dt * dwdt
        wcanp(i,j) = MAX( 0.0, wcanp(i,j))
        wcanp(i,j) = MIN( .1 * veg(i,j) , wcanp(i,j))
      enddo
      enddo
      if(debug) write(*,*) ' VEGEBM: exit '
!-----------------------------------------------------------------------
      RETURN
      END
