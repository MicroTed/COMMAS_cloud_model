!-----------------------------------------------------------------------------
!
! ///////////////////// BEGIN \\\\\\\\\\\\\\\\\\\\
! \\\\\\\\\\\\\\\\\\\\\ SUBROUTINE CORIOLIS ////////////////////
!
! CORIOLIS computes the Coriolis acceleration for the horizontal winds
!-----------------------------------------------------------------------------
! Created by Louis Wicker, November 11, 1991
! Latest update: 8-27-03
! Last modified by M. Gilmore
!
! Update Notes:
!
! NEW A-GRID VERSION!
!-----------------------------------------------------------------------------
       SUBROUTINE CORIOLIS(u,v,w,fu,fv,fw,ub,vb,ugrid,vgrid,lat,nx,ny,nz)

       USE PARAM_MODULE
       USE COMMASMPI_MODULE

       implicit none

       integer nx,ny,nz
       real ugrid, vgrid
       real :: lat   ! latitude
       real ub(-ng+1:nz+ng),vb(-ng+1:nz+ng)
       real u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
       real v (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
       real w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
       real fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
       real fv(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
       real fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)

! Local variables

       integer i, j, k
       real uavg2v, vavg2u, wavg2u, uavg2w
       real fcor1,fcor2

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
       integer :: ixb, jyb, kzb
       integer :: ixe, jye, kze

!-----------------------------------------------------------------------------

        IF( coriol .eq. 0.0 ) RETURN

!-----------------------------------------------------------------------------
! U-V acceleration
!  Coriolis = -2 (Omega x U) where Omega= Omega cos(phi) j + Omega sin(phi) k
!                            and where U = u i + v j + w k
!
!  Add in grid motion to model winds to get total ground-relative wind.
!  Subtract out reference state since Coriolis acts on the perturbation wind.
!-----------------------------------------------------------------------------

       ! DTD: Added calculation of coriolis parameter based on central latitude of grid

       fcor1 = 2.0*(7.292e-5)*Sin(ATan(1.0)/45.0*lat)
       fcor2 = 2.0*(7.292e-5)*Cos(ATan(1.0)/45.0*lat)       

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg+1

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg+1

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg+1

       do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
       DO k = 1,nz-1
        DO j = 1,ny-1
         DO i = 1,nx-1
#endif
!  First 2 terms are needed in the synoptic scale when U=ui+vj
!  Additional 2 terms needed when w>>0 as in a thunderstorm

!          fu(i,j,k) = fu(i,j,k) + coriol*(v(i,j,k)-w(i,j,k) -vb(k) + vgrid)
!          fv(i,j,k) = fv(i,j,k) - coriol*(u(i,j,k)          -ub(k) + ugrid)
!          fw(i,j,k) = fw(i,j,k) + coriol*(u(i,j,k)          -ub(k) + ugrid)

          !DTD: corrected calculation of coriolis force.  Technically still not 
          !quite correct because v and u are staggered.

          fu(i,j,k) = fu(i,j,k) + fcor1*(v(i,j,k)-vb(k)+vgrid)-fcor2*w(i,j,k)
          fv(i,j,k) = fv(i,j,k) - fcor1*(u(i,j,k)          -ub(k) + ugrid)
          fw(i,j,k) = fw(i,j,k) + fcor2*(u(i,j,k)          -ub(k) + ugrid)

         ENDDO
        ENDDO
       ENDDO

       RETURN
       END SUBROUTINE CORIOLIS
