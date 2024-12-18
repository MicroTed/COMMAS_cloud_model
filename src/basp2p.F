      SUBROUTINE SETCLR
      
      USE PLOT_MODULE, only: dbz_color_option
      
      implicit none
      integer i,j
      real x,dbz,r,g,b
      
        CALL GSCR(1,50,0.75,0.75,0.75)  !  25% gray
        CALL GSCR(1,51,0.80,0.80,0.80)  !  20% gray
C Second foreground color is light blue
c        CALL GSCR(1,52,94.0/255.0,208.0/255.0,1.0)
        CALL GSCR(1,52,65.0/255.0,183.0/255.0,212./255.)

        CALL GSCR(1,61,247./255.,246./255.,25./255.)   ! yellow
        CALL GSCR(1,261,247./255.,246./255.,25./255.)   ! yellow
        CALL GSCR(1,62,1.0,0.5,0.0)   ! orange
        CALL GSCR(1,262,1.0,0.5,0.0)   ! orange
        CALL GSCR(1,53,1.0,1.0,0.5)     ! light yellow
        CALL GSCR(1,253,1.0,1.0,0.5)     ! light yellow
        CALL GSCR(1,54,1.0,0.5,0.5)     ! light pink
C set color to green
        CALL GSCR(1,55,0.0,0.7,0.0)
        CALL GSCR(1,255,0.0,0.7,0.0)
C set color to blue
c        CALL GSCR(1,56,0.0,0.0,1.0)
        CALL GSCR(1,56,25./255.,76./255.,159./255.)
CC set color to red
        CALL GSCR(1,57,0.7,0.0,0.0)
c 29% gray
        call gscr(1,58,0.71,0.71,0.71)
c black
        call gscr(1,59,0.0,0.0,0.0)
c almost black
        call gscr(1,60,0.05,0.05,0.05)

C set color to other red
        CALL GSCR(1,83,249.0/255.,10./255.,2./255.)
C set color to other green
        CALL GSCR(1,85,21.0/255.,98./255.,35./255.)
C set color to darker blue
        CALL GSCR(1,86,29.0/255.,14./255.,130./255.)

c 99% black (for vector default color)
        CALL GSCR(1,99,0.01,0.01,0.01)  !  20% gray

      j=70
      DO i = 0,10
       x = Max(0.0, Min(1.0,1.0 - float(i)/10.))
        call gscr(1,i+j,x,x,x)
      ENDDO
      
      
      DO i = 200,215
       dbz = 5.0*(i - 200+1) - 0.5
   
      IF ( dbz_color_option == 1 ) THEN
       IF  ( dbz .lt. 5.0) THEN
       r = 220
       g = 220
       b = 220
       ELSEIF ( dbz .lt. 10.0) THEN
       r = 156
       g = 156
       b = 156
       ELSEIF ( dbz .lt. 15.0) THEN
       r = 118
       g = 118
       b = 118
       ELSEIF ( dbz .lt. 20.0) THEN
c       r = 170
c       g = 170
c       b = 255
       r = 158
       g = 134
       b = 124
       ELSEIF ( dbz .lt. 25.0) THEN
c       r = 120
c       g = 120
c       b = 240
       r = 10
       g = 80
       b = 160
       ELSEIF ( dbz .lt. 30.0) THEN
c       r = 90
c       g = 90
c       b = 220
       r = 85
       g = 87
       b = 166
       ELSEIF ( dbz .lt. 35.0) THEN
c       r = 0
c       g = 251
c       b = 144
       r = 78
       g = 174
       b = 163
       ELSEIF ( dbz .lt. 40.0) THEN
c       r = 0
c       g = 187
c       b = 0
       r = 17
       g = 142
       b = 52
       ELSEIF ( dbz .lt. 45.0) THEN
       r = 255
       g = 255
       b = 112
       ELSEIF ( dbz .lt. 50.0) THEN
c       r = 208
c       g = 208
c       b = 96
       r = 252
       g = 164
       b = 73
       ELSEIF ( dbz .lt. 55.0) THEN
       r = 255
       g = 96
       b = 96
       ELSEIF ( dbz .lt. 60.0) THEN
       r = 236
       g = 51
       b = 19
       ELSEIF ( dbz .lt. 65.0) THEN
       r = 174
       g = 0
       b = 0
       ELSEIF ( dbz .lt. 70.0) THEN
c       r = 255
c       g = 0
c       b = 255
       r = 223
       g = 90
       b = 159
       ELSEIF ( dbz .lt. 75.0) THEN
       r = 255
       g = 255
       b = 255
       ELSE
       r = 100
       g = 100
       b = 100
       ENDIF

      ELSEIF ( dbz_color_option /= 1 ) THEN
       IF  ( dbz .lt. 5.0) THEN
       r = 220
       g = 220
       b = 220
       ELSEIF ( dbz .lt. 10.0) THEN
       r = 99 
       g = 153 
       b = 194
       ELSEIF ( dbz .lt. 15.0) THEN
       r = 99  
       g = 194
       b = 219
       ELSEIF ( dbz .lt. 20.0) THEN
       r = 51  
       g = 102
       b = 51
       ELSEIF ( dbz .lt. 25.0) THEN
       r = 48  
       g = 153
       b = 102
       ELSEIF ( dbz .lt. 30.0) THEN
       r = 153  
       g = 201
       b = 153
       ELSEIF ( dbz .lt. 35.0) THEN
       r = 245  
       g = 247
       b = 201
       ELSEIF ( dbz .lt. 40.0) THEN
       r = 247  
       g = 245
       b = 51
       ELSEIF ( dbz .lt. 45.0) THEN
       r = 250  
       g = 204
       b = 48
       ELSEIF ( dbz .lt. 50.0) THEN
       r = 222  
       g = 145
       b = 46
       ELSEIF ( dbz .lt. 55.0) THEN
       r = 250  
       g = 99
       b = 18
       ELSEIF ( dbz .lt. 60.0) THEN
       r = 250  
       g = 15
       b = 13
       ELSEIF ( dbz .lt. 65.0) THEN
       r = 204  
       g = 48
       b = 51
       ELSEIF ( dbz .lt. 70.0) THEN
       r = 176  
       g = 10
       b = 51
       ELSEIF ( dbz .lt. 75.0) THEN
       r = 204  
       g = 102
       b = 153
       ELSE
       r = 100
       g = 100
       b = 100
       ENDIF
       
       ENDIF



      
      r = r/255.
      g = g/255.
      b = b/255.
      
c      write(0,*) 'dbz,r,g,b = ',dbz,r,g,b
      
      CALL gscr(1,i,r,g,b)
      ENDDO ! i

      DO i = 1,2
       call gscr(1, 200*(i-1)+20 ,0.142,  0.000,  0.850) ! dark blue
       call gscr(1, 200*(i-1)+21 ,0.097,  0.112,  0.970)
       call gscr(1, 200*(i-1)+22 ,0.160,  0.342,  1.000)
       call gscr(1, 200*(i-1)+23 ,0.240,  0.531,  1.000)
       call gscr(1, 200*(i-1)+24 ,0.340,  0.692,  1.000)
       call gscr(1, 200*(i-1)+25 ,0.460,  0.829,  1.000)
       call gscr(1, 200*(i-1)+26 ,0.600,  0.920,  1.000)
       call gscr(1, 200*(i-1)+27 ,0.740,  0.978,  1.000)     
       call gscr(1, 200*(i-1)+28 ,0.900,  0.900,  0.900) ! last negative
       call gscr(1, 200*(i-1)+29 ,1.000,  1.000,  0.920) ! first positive
       call gscr(1, 200*(i-1)+30 ,1.000,  0.948,  0.740)
       call gscr(1, 200*(i-1)+31 ,1.000,  0.840,  0.600)
       call gscr(1, 200*(i-1)+32 ,1.000,  0.676,  0.460)
       call gscr(1, 200*(i-1)+33 ,1.000,  0.472,  0.340)
       call gscr(1, 200*(i-1)+34 ,1.000,  0.240,  0.240)
       call gscr(1, 200*(i-1)+35 ,0.970,  0.155,  0.210)
       call gscr(1, 200*(i-1)+36 ,0.850,  0.085,  0.187)
       call gscr(1, 200*(i-1)+37 ,0.650,  0.000,  0.130) ! dark red
      ENDDO
      
      RETURN
      END

! #####################################################################
      SUBROUTINE CSTCLR (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
      
      use plot_module, only: maxclr, customcolors
C
        REAL :: XCRA(*),YCRA(*)
        INTEGER :: IAIA(*),IGIA(*)
        
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
c        write(0,*) 'naia = ',naia
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
c          write(*,*) 'IFLL = ',ifll
          IF (IFLL.GT.0.AND.IFLL.LE.maxclr) THEN
c            write(*,*) 'ifll+220 = ',ifll+220,ifll,naia
              CALL GSFACI (customcolors(ifll))
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
c

! #####################################################################
      SUBROUTINE BLURED (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
C
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
        
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
c        write(0,*) 'naia = ',naia
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
c          write(*,*) 'IFLL = ',ifll
          IF (IFLL.GT.0.AND.IFLL.LE.18) THEN
c            write(*,*) 'ifll+220 = ',ifll+220,ifll,naia
              CALL GSFACI (IFLL+220-1)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
c

! #####################################################################

! #####################################################################
      SUBROUTINE REDPOS (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
C
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
c        write(0,*) 'naia = ',naia
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
          IF (IFLL.GT.0.AND.IFLL.LE.9) THEN
c            write(0,*) 'ifll+200 = ',ifll+200
            CALL GSFACI (IFLL+229-1)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
c

! #####################################################################
      SUBROUTINE BLUNEG (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
C
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
c        write(0,*) 'naia = ',naia
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
          IF (IFLL.GT.0.AND.IFLL.LE.9) THEN
c            write(0,*) 'ifll+200 = ',ifll+200
            CALL GSFACI (IFLL+200-1)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
c

! #####################################################################


! #####################################################################
      SUBROUTINE COLDBZ (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
C
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
c        write(0,*) 'naia = ',naia
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
          IF (IFLL.GT.0.AND.IFLL.LE.15) THEN
c            write(0,*) 'ifll+200 = ',ifll+200
            CALL GSFACI (IFLL+200)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
c

! #####################################################################

! #####################################################################

      SUBROUTINE SHDCLR (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
       integer ishdclr
       common /clrstuff/ ishdclr
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C

        IF (ISH.NE.0) THEN
         CALL GSFACI(ishdclr)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END


c
c
c  begin dxlabs
c
c writes distance values at set intervals along the axis
c
      subroutine dxlabs(it,xo,xn,yo,yn,nx,yvalp,ixdiv,dx,dxp,xfac,irev)
      character*7 lab1,dxlab,lab0
      character*24 xlabel
      integer nchar,ilen,ic
      integer irev
      nchar = 6
      xscl = xn-xo

      do 10 ix = 1,nx
      xvalp = xscl*(float(ix-1))/(float(nx-1))
      IF ( irev .eq. 0 ) THEN
        xdist = ((ix-1)*dx + dxp)*xfac
      ELSEIF ( irev .eq. 1 ) THEN  ! reverse direction
        xdist = ((nx-ix)*dx + dxp)*xfac
      ENDIF
      if ( (ix-1)/ixdiv*ixdiv .eq. (ix-1) ) then
      IF ( (xdist - float(int(xdist))) .le. 0.09 ) then
        IF ( xdist .lt. 1e7 ) THEN
          write(lab1,201) Int(xdist)
        ELSE
          write(lab1,201) Int(xdist)
        ENDIF
      ELSE
        write(lab1,200) xdist
      END IF
 200  format(f5.1)
 201  format(i6)
 202  format(1pe7.1)
      axvalp = xo + xvalp
      ayvalp = yo + 2.7*yvalp
      ilen = 0
      DO ic=1,7
       IF ( lab1(ic:ic) .ne. ' ' ) THEN
        ilen = ilen + 1
        write(lab0(ilen:ilen),'(a1)') lab1(ic:ic)
       ENDIF
      ENDDO
c      call pwrit(axvalp,ayvalp,lab1,6,1,0,0)
      call pwrit(axvalp,ayvalp,lab0,ilen,1,0,0)
c      call pwrit(axvalp,ayvalp,' ',1,1,0,0) ! mark the location of the label
      ELSE
!        write(6,*) 'dxlabs: bad interval ix-1,ixdiv = ',ix-1,ixdiv,(ix-1)/ixdiv*ixdiv
      end if
  10  continue
      axvalp = xscl/2. + xo
      xlabel = ' '
      if ( xfac .lt. 1. ) then
      if ( it .eq. 1 ) xlabel = 'X (km)'
      if ( it .eq. 2 ) xlabel = 'X (km)'
      if ( it .eq. 3 ) xlabel = 'Y (km)'
      if ( it .eq. 5 ) THEN 
        xlabel = 'Horizontal Distance (km)'
        nchar = 24
      END IF
      
      IF ( it .eq. 11 ) THEN
       nchar = 8
       IF (xfac .eq. 0.001) xlabel = '(x 1000)'
       IF (xfac .eq. 0.01) xlabel = '(x 100)'
       IF (xfac .eq. 0.1) xlabel = '(x 10)'
      ENDIF
      else
      if ( it .eq. 1 ) xlabel = 'X (M)'
      if ( it .eq. 2 ) xlabel = 'X (M)'
      if ( it .eq. 3 ) xlabel = 'Y (M)'
      if ( it .eq. 5 ) xlabel = 'R (M)'
      IF ( it .eq. 10 ) THEN
       nchar = 10
       xlabel = 'Time (min)'
      ENDIF
      end if
      ayvalp = yo + 2.0*yvalp
      call pwrit(axvalp,ayvalp,xlabel,nchar,1,0,0)
      return
      end
c
c
c





c
c ny is nyp or nzp
c
      subroutine dylabs(it,yo,yn,xo,xn,ny,xvalp,axb
     >                 ,jydiv,dy,dyp,yfac,itlab,templt,zk,nz,kb1)
      implicit none
      character*6 lab1,lab2,ylabel,lab0
      character*15 tlabel
      integer nz,ny
      real templt(nz),zk(nz)
      real yo,yn,xo,xn,axb
      real xvalp,dy,dyp,yfac
      integer it,jydiv,itlab,ic
      real yscl,yvalp,ydist,ayvalp,axvalp,axvalpt
      integer jy,jy2,jy0,kb1
      real yint,dyscl,yvalp2
      real line(2,2)
      integer ierr
      real width,w,yl
      integer ny0
      integer titlen, titlenmx
      
      jy0 = Max(0, kb1-1)
      ny0 = ny + jy0
      w = 2.0
      titlenmx = 0
      
      yscl = yn-yo
      dyscl = yscl/float(ny-1)
      
      DO 10 jy = 1,ny
       yvalp =  yscl*float(jy-1)/float(ny-1)
       yvalp2 = yscl*float(jy)/float(ny-1)
       dyscl = yvalp2 - yvalp
       ydist = ((jy-1)*dy + dyp)*yfac
       
       jy2 = Min(nz-1, jy + Int(dyp/dy))
       
c       write(0,*) 'dylabs: it,itlab,jy,jy2,ydist,dyp,ny,ny0,jy0 ',
c     : it,itlab,jy,jy2,ydist,dyp,templt(jy2),ny,ny0,jy0 
       
       if ( (jy-1)/jydiv*jydiv .eq. (jy-1) ) then
       IF ( (ydist - float(int(ydist))) .le. 0.09 ) then
         write(lab1,199) Int(ydist)
       ELSE
         write(lab1,200) ydist
       END IF
 199  format(i4)
 200  format(f3.0)
      write(lab2,201) templt(jy)
 201  format(f5.1)
 202  format(i3)
      ayvalp = yo + yvalp
      axvalp = xo + 4.5*xvalp
c     CALL PWRIT (PX,PY,CH,NC,IS,IO,IC)
      titlen = 0
      
      DO ic=1,4
       IF ( lab1(ic:ic) .ne. ' ' ) THEN
        titlen = titlen + 1
        write(lab0(titlen:titlen),'(a1)') lab1(ic:ic)
       ENDIF
      ENDDO
      titlenmx = Max( titlen, titlenmx)

c      titlen = len(lab1)
c      CALL strlnth(lab1,titlen)

      call pwrit(axvalp,ayvalp,lab0,titlen,1,0,1)
      if ( (it .eq. 2 .or. it. eq. 3 .or. it .eq. 5 .or. it .ge. 10) 
     :   .and.  itlab .ge. 1 ) then
       IF ( itlab .eq. 1 .or. itlab .eq. 3 ) THEN
         axvalpt = xn + (0.35)*xvalp
             titlen = len(lab2)
             CALL strlnth(lab2,titlen)
c             write(0,*) 'titlen1 = ',titlen
         call pwrit(axvalpt,ayvalp,lab2,titlen,1,0,-1)
       ENDIF
      end if
      end if
       
      IF ( (it .eq. 2 .or. it. eq. 3 .or. it .eq. 5 .or. it .ge. 10) 
     :   .and.  itlab .ge. 1 .and. jy .lt. ny) then
       IF (itlab .eq. 2 .or. itlab .ge. 3 ) THEN
c        write(0,*) '1:xo,xn,yvalp,dyscl = ',xo,xn,yvalp,dyscl
         IF ( templt(jy2) .ge. 0.0 .and. 
     :    templt(Min(ny0,Min(nz,jy2+1))) .le. 0.0 ) THEN
c          write(0,*) 'jy2,templt(jy2),templt(jy2+1) = ',
c     :      jy2,templt(jy2),templt(jy2+1)
         yint = dyscl*(templt(jy2) - 0.0)/(templt(jy2) - templt(jy2+1))
c         write(0,*) 'xo,xn,ayvalp,yint = ',xo,xn,ayvalp,yint
         line(1,1) = axb - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = axb + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call GSPLCI(75)
         CALL GQLWSC (ierr, width)
         CALL GSLWSC(w)
         call CURVE(line(1,1),line(1,2),2)
         line(1,1) = xn - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = xn + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call CURVE(line(1,1),line(1,2),2)
         call GSPLCI(1)
         CALL GSLWSC(width)
           IF ( itlab .eq. 2 ) THEN
             ayvalp = yo + yvalp + yint
             write(lab2,'(i1)') 0
             axvalpt = xn + (1.1)*xvalp
             
             titlen = len(lab2)
             CALL strlnth(lab2,titlen)
c             write(0,*) 'titlen = ',titlen
             call pwrit(axvalpt,ayvalp,lab2,titlen,1,0,-1)
           ENDIF
         ELSEIF (templt(jy2) .ge. -10. .and. 
     :       templt (Min(ny0,Min(nz,jy2+1))) .le. -10.0 )  THEN
         yint = dyscl*(templt(jy2) + 10.0)/(templt(jy2) - templt(jy2+1))
c         write(0,*) 'xo,xn,ayvalp,yint = ',xo,xn,ayvalp,yint
         line(1,1) = axb - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = axb + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         CALL GQLWSC (ierr, width)
         CALL GSLWSC(w)
         call GSPLCI(75)
         call CURVE(line(1,1),line(1,2),2)
         line(1,1) = xn - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = xn + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call CURVE(line(1,1),line(1,2),2)
         call GSPLCI(1)
         CALL GSLWSC(width)
           IF ( itlab .eq. 2 ) THEN
             ayvalp = yo + yvalp + yint
             write(lab2,202) -10
             axvalpt = xn + (0.45)*xvalp
             titlen = len(lab2)
             CALL strlnth(lab2,titlen)
c             write(0,*) 'titlen = ',titlen
             call pwrit(axvalpt,ayvalp,lab2,titlen,1,0,-1)
           ENDIF

         ELSEIF (templt(jy2) .ge. -20. .and. 
     :      templt (Min(ny0,Min(nz,jy2+1))) .le. -20.0 )   THEN
         yint = dyscl*(templt(jy2) + 20.0)/(templt(jy2) - templt(jy2+1))
c         write(0,*) 'xo,xn,ayvalp,yint = ',xo,xn,ayvalp,yint
         line(1,1) = axb - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = axb + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call GSPLCI(75)
         CALL GQLWSC (ierr, width)
         CALL GSLWSC(w)
         call CURVE(line(1,1),line(1,2),2)
         line(1,1) = xn - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = xn + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call CURVE(line(1,1),line(1,2),2)
         call GSPLCI(1)
         CALL GSLWSC(width)
           IF ( itlab .eq. 2 ) THEN
             ayvalp = yo + yvalp + yint
             write(lab2,202) -20
             axvalpt = xn + (0.45)*xvalp
             titlen = len(lab2)
             CALL strlnth(lab2,titlen)
c             write(0,*) 'titlen = ',titlen
             call pwrit(axvalpt,ayvalp,lab2,titlen,1,0,-1)
           ENDIF
         ELSEIF (templt(jy2) .ge. -40. .and.  
     :      templt (Min(ny0,Min(nz,jy2+1))) .le. -40.0 )  THEN
         yint = dyscl*(templt(jy2) + 40.0)/(templt(jy2) - templt(jy2+1))
c         write(0,*) 'xo,xn,ayvalp,yint = ',xo,xn,ayvalp,yint
         line(1,1) = axb - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = axb + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         CALL GQLWSC (ierr, width)
         CALL GSLWSC(w)
         call GSPLCI(75)
         call CURVE(line(1,1),line(1,2),2)
         line(1,1) = xn - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = xn + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call CURVE(line(1,1),line(1,2),2)
         call GSPLCI(1)
         CALL GSLWSC(width)
           IF ( itlab .eq. 2 ) THEN
             ayvalp = yo + yvalp + yint
             write(lab2,202) -40
             axvalpt = xn + (0.45)*xvalp
             titlen = len(lab2)
             CALL strlnth(lab2,titlen)
c             write(0,*) 'titlen = ',titlen
             call pwrit(axvalpt,ayvalp,lab2,titlen,1,0,-1)
           ENDIF

         ENDIF
         IF (templt(jy2) .ge. -30. .and. 
     :      templt (Min(ny0,Min(nz,jy2+1))) .le. -30.0 )   THEN
         yint = dyscl*(templt(jy2) + 30.0)/(templt(jy2) - templt(jy2+1))
c         write(0,*) 'xo,xn,ayvalp,yint = ',xo,xn,ayvalp,yint
         line(1,1) = axb - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = axb + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call GSPLCI(75)
         CALL GQLWSC (ierr, width)
         CALL GSLWSC(w)
         call CURVE(line(1,1),line(1,2),2)
         line(1,1) = xn - 0.35*xvalp
         line(1,2) = yo + yvalp + yint
         line(2,1) = xn + 0.35*xvalp
         line(2,2) = yo + yvalp + yint
         call CURVE(line(1,1),line(1,2),2)
         call GSPLCI(1)
         CALL GSLWSC(width)
           IF ( itlab .eq. 2 ) THEN
             ayvalp = yo + yvalp + yint
             write(lab2,202) -30
             axvalpt = xn + (0.45)*xvalp
             titlen = len(lab2)
             CALL strlnth(lab2,titlen)
c             write(0,*) 'titlen = ',titlen
             call pwrit(axvalpt,ayvalp,lab2,titlen,1,0,-1)
           ENDIF
         ENDIF
       ENDIF
      ENDIF
  10  continue
      ayvalp = yscl/2. + yo
      if ( yfac .lt. 1. ) then
      if ( it .eq. 1 ) ylabel = 'Y (km)'
      if ( it .eq. 2 .or. it .eq. 3 .or. it .ge. 5 ) ylabel = 'Z (km)'
       yl = 3.2 ! 3.5
!       yl = 6.5 - titlenmx
      else
       yl = 2.5
      if ( it .eq. 1 ) ylabel = 'Y (M)'
      if ( it .eq. 2 .or. it .eq. 3 .or. it .eq. 5 ) ylabel = 'Z  (M)'
      end if
      axvalp = xo + yl*xvalp
      call pwrit(axvalp,ayvalp,ylabel,6,1,90,0)
      if ( itlab .eq. 1 .or. itlab .eq. 3 ) then
      if ( it .eq. 2 .or. it .eq. 3 .or. it .ge. 5 ) then
      tlabel = 'TEMPERATURE (C)'
      axvalpt = xn + (3.0)*xvalp 
      call pwrit(axvalpt,ayvalp,tlabel,15,1,90,0)
      end if
      end if
      return
      end
c
c
c


c
c  subroutine dfclr
c
c
      subroutine dfclr(idx,isign,irecmj,irecmn,irectx,ndot)
      integer idx,irecmj,irecmn,irectx,isign,ndot
c
c  idx    is the variable index
c  isign  is the sign of the array
c
c  colors
c
c    1:  black
c    2:  red
c    3:  green
c    4:  blue
c    5:  cyan (aqua-blue)
c    6:  magenta 
c    7:  yellow
c    8:  black
c
c  all color default to: red for positive #'s, blue for positive #'s
c  and black for text
c
      irectx = 1
      if ( isign .le. 0 ) then
      ndot = -682
      irecmj = 4
      end if
      if ( isign .gt. 0 ) then
      ndot = -682
      irecmj = 15
      end if
      if ( isign .gt. 990 ) then
      ndot = -682
      irecmj = 1
      end if
c
c  Special cases
c
c  cloud water
c  
      if ( idx .eq. 11 ) then
      irecmj = 3
      end if
c
c  rain water
c  
      if ( idx .eq. 12 ) then
      irecmj = 3
      end if
c
c  cloud ice water
c  
      if ( idx .eq. 13 ) then
      irecmj = 4
      end if
c
c  snow water
c  
      if ( idx .eq. 14 ) then
      irecmj = 4
      end if
c
c  hail water
c  
      if ( idx .eq. 15 ) then
      irecmj = 2
      end if
c
c  cloudwater and cloudice 
c  
      if ( idx .eq. 16 ) then
      irecmj = 5
      end if
c
c  visible cloud
c  
      if ( idx .eq. 17 ) then
      irecmj = 8
      end if
c
c  total ice
c  
      if ( idx .eq. 18 ) then
      irecmj = 4
      end if
c
c  total liquid
c  
      if ( idx .eq. 19 ) then
      irecmj = 3
      end if
c
c  total precip
c  
      if ( idx .eq. 20 ) then
      irecmj = 3
      end if
c
c  total condenstate
c  
      if ( idx .eq. 21 ) then
      irecmj = 3
      end if
c
c  radar dBZ
c  
      if ( idx .eq. 22 ) then
      irecmj = 4
      end if
c   
c  set irecmn to irecmj
c
      irecmn = irecmj
c
c
c
      return
      end
c
c
c
      subroutine extrms(ia,nor,na,nx,ny,nz
     >  ,imx,imn,jmx,jmn,kmx,kmn,xmax,xmin,v)
      implicit none
      integer ia,nor,na,nx,ny,nz
      integer imx,imn,jmx,jmn,kmx,kmn
      integer ix,jy,kz
      real v(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real xmax,xmin
      xmax = -1.e10
      xmin =  1.e10
      do 10 jy = 1,ny
      do 11 kz = 1,nz
      do 12 ix = 1,nx
      if ( v(ix,jy,kz,ia) .gt. xmax ) then
      imx = ix  
      jmx = jy 
      kmx = kz  
      xmax = v(ix,jy,kz,ia)
      end if
      if ( v(ix,jy,kz,ia) .lt. xmin ) then
      imn = ix  
      jmn = jy 
      kmn = kz  
      xmin = v(ix,jy,kz,ia)
      end if
  12  continue
  11  continue
  10  continue
      return
      end
c
c
c
      subroutine readmp
     >  ( nx, ny, nz, dum )
c
c     implicit none
      character*15 fnmnsv,fnmloc,fnmval
      integer ngs,nch
      parameter(ngs=500,nch=22)
      integer nx,ny,nz
      integer n,mgs,ngscnt,l,ncuse,iotype
      integer igs(ngs),kgs(ngs),jgs
      real val(ngs,nch)
      real dum(nx,ny,nz)
c  
      write(*,*) 'INPUT:  filenames (characters);'
      write(*,*) '    1:  nsvcnt'
      write(*,*) '    2:  fnmloc'
      write(*,*) '    3:  fnmval'
      write(*,*) 'AND 4:  iotype (integer; 1=unform,2=form)'
      write(*,*) 'AND 5:  ncuse  (number of categories)'
      read(22,*)  fnmnsv
      read(22,*)  fnmloc
      read(22,*)  fnmval
      read(22,*)  iotype,ncuse  
c
c
c
c  open files fnmloc and fnmval which contain location and values
c  and number saves:  nsvcnt
c
      open(unit=087,file=fnmnsv,form='formatted',status='old')
c
      if ( iotype .eq. 1 ) then
      open(unit=088,file=fnmloc,form='unformatted',status='old')
      end if
c
      if ( iotype .eq. 2 ) then
      open(unit=088,file=fnmloc,form='formatted',status='old')
      end if
c
      if ( iotype .eq. 1 ) then
      open(unit=089,file=fnmval,form='unformatted',status='old')
      end if
c
      if ( iotype .eq. 2 ) then
      open(unit=089,file=fnmval,form='formatted',status='old')
      end if
c
c  xero dum
c
      do 09991 kz = 1,nz
      do 09992 jy = 1,ny
      do 09993 ix = 1,nx
      dum(ix,jy,kz) = 0.0
09993 continue
09992 continue
09991 continue
c
c  read in data
c
      read(087,*) nsvcnt
c
      do 10000 n = 1,nsvcnt
c
      if ( iotype .eq. 1 ) then
      read(088)  ngscnt
      read(088) (igs(mgs),mgs=1,ngscnt)
      read(088) (kgs(mgs),mgs=1,ngscnt)
      read(088) jgs
      end if
      if ( iotype .eq. 2 ) then
      read(088,1001)  ngscnt
      read(088,1001) (igs(mgs),mgs=1,ngscnt)
      read(088,1001) (kgs(mgs),mgs=1,ngscnt)
      read(088,1001) jgs
1001  format(4i16)
      end if
c
      do 10010 l = 1,ncuse
      if ( iotype .eq. 1 ) then
      read(089) (val(mgs,l),mgs=1,ngscnt)
      end if
      if ( iotype .eq. 2 ) then
      read(089,1000) (val(mgs,l),mgs=1,ngscnt)
1000  format(4e16.7)
      end if
10010 continue
c
      do 10020   l = 1,ncuse
      do 10021 mgs = 1,ngscnt
      dum(igs(mgs),kgs(mgs),jgs) = val(mgs,l)
10021 continue
10020 continue
c
10000 continue
c
      close(unit=087)
      close(unit=088)
      close(unit=089)
c
      return
      end
c
c
c
c
c
c
      subroutine readmc
     >  ( na, nx, ny, nz, dumc )
c
c     implicit none
      character*15 fnmnsv,fnmloc,fnmval
      integer ngs,nch
      parameter(ngs=500,nch=22)
      integer na,nx,ny,nz
      integer n,mgs,ngscnt,l,ncuse,iotype
      integer igs(ngs),kgs(ngs),jgs
      real val(ngs,nch)
      real dumc(nx,ny,nz,na)
c
      write(*,*) 'INPUT:  filenames (characters);'
      write(*,*) '    1:  nsvcnt'
      write(*,*) '    2:  fnmloc'
      write(*,*) '    3:  fnmval'
      write(*,*) 'AND 4:  iotype (integer; 1=unform,2=form)'
      write(*,*) 'AND 5:  ncuse  (number of categories)'
      read(22,*)  fnmnsv
      read(22,*)  fnmloc
      read(22,*)  fnmval
      read(22,*)  iotype,ncuse
c
c
c
c  open files fnmloc and fnmval which contain location and values
c  and number saves:  nsvcnt
c
      open(unit=087,file=fnmnsv,form='formatted',status='old')
c
      if ( iotype .eq. 1 ) then
      open(unit=088,file=fnmloc,form='unformatted',status='old')
      end if
c
      if ( iotype .eq. 2 ) then
      open(unit=088,file=fnmloc,form='formatted',status='old')
      end if
c
      if ( iotype .eq. 1 ) then
      open(unit=089,file=fnmval,form='unformatted',status='old')
      end if
c
      if ( iotype .eq. 2 ) then
      open(unit=089,file=fnmval,form='formatted',status='old')
      end if
c
c  xero dum
c
      do 09990 l  = 1,na
      do 09991 kz = 1,nz
      do 09992 jy = 1,ny
      do 09993 ix = 1,nx
      dumc(ix,jy,kz,l) = 0.0
09993 continue
09992 continue
09991 continue
09990 continue
c
c  read in data
c
      read(087,*) nsvcnt
c
      do 10000 n = 1,nsvcnt
c
      if ( iotype .eq. 1 ) then
      read(088)  ngscnt
      read(088) (igs(mgs),mgs=1,ngscnt)
      read(088) (kgs(mgs),mgs=1,ngscnt)
      read(088) jgs
      end if
      if ( iotype .eq. 2 ) then
      read(088,1001)  ngscnt
      read(088,1001) (igs(mgs),mgs=1,ngscnt)
      read(088,1001) (kgs(mgs),mgs=1,ngscnt)
      read(088,1001) jgs
1001  format(4i16)
      end if
c
      do 10010 l = 1,ncuse
      if ( iotype .eq. 1 ) then
      read(089) (val(mgs,l),mgs=1,ngscnt)
      end if
      if ( iotype .eq. 2 ) then
      read(089,1000) (val(mgs,l),mgs=1,ngscnt)
1000  format(4e16.7)
      end if
10010 continue
c
      do 10020   l = 1,ncuse
      do 10021 mgs = 1,ngscnt
      dumc(igs(mgs),kgs(mgs),l,jgs) = val(mgs,l)
10021 continue
10020 continue
c
10000 continue
c
      close(unit=087)
      close(unit=088)
      close(unit=089)
c
      return
      end
c
c
c



c
c  routine smooth
c
      subroutine smth
     >  (nsmth,ismthx,ismthy,ismthz,ibase,
     >   nx,ny,nz,na,nor,nxt,nyt,nzt,nxl,nyl,nzl,
     >   nxbdo,nybdo,nzbdo,nabdo,
     >   nxedo,nyedo,nzedo,naedo,
     >   t0,at,ab)
c
      integer    ia
      integer    ii
      integer    jj
      integer    kk
      integer    ix
      integer    jy
      integer    kz
      integer    nx
      integer    ny
      integer    nz
      integer    nxl
      integer    nyl
      integer    nzl
      integer    nxt
      integer    nyt
      integer    nzt
      integer    nxbdo
      integer    nybdo
      integer    nzbdo
      integer    nxedo
      integer    nyedo
      integer    nzedo
      integer    na
      integer    nabdo
      integer    naedo
      integer    nsmth
      integer    ismthx
      integer    ismthy
      integer    ismthz
      integer    ibase
      integer    ismth  
c
      real       t0(-nor+1:nxt+nor,-nor+1:nyt+nor,-nor+1:nzt+nor)
      real       ab(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real       at(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor,na)
c
c  smooth data
c
      do 10000 ismth = 1, nsmth
c
      do 10010 ia =    nabdo, naedo
c
      do 11001 ii = 1,ismthx
c
c     print*,'SMOOTHX',ii,ismth
c
c  x-direction
c
      do 10101 kz = nzbdo,nzedo
      do 10102 jy = nybdo,nyedo
      do 10103 ix = nxbdo+1,nxedo-1
      t0(ix,jy,kz) =  0.25*(at(ix+1,jy,kz,ia)+at(ix-1,jy,kz,ia)
     >                + 2.0*at(ix,jy,kz,ia))
10103 continue
10102 continue
10101 continue
c     print*,'SMOOTHXR',ii,ismth
      do 10111 kz = nzbdo,nzedo
      do 10112 jy = nybdo,nyedo
      do 10113 ix = nxbdo+1,nxedo-1
      at(ix,jy,kz,ia) = t0(ix,jy,kz)
10113 continue
10112 continue
10111 continue
c
11001 continue
c
c  y-direction
c
      do 12001 jj = 1,ismthy
c     print*,'SMOOTHY',jj,ismth
      do 10201 kz = nzbdo,nzedo
      do 10202 jy = nybdo+1,nyedo-1
      do 10203 ix = nxbdo,nxedo
      t0(ix,jy,kz) =  0.25*(at(ix,jy+1,kz,ia)+at(ix,jy-1,kz,ia)
     >                + 2.0*at(ix,jy,kz,ia))
10203 continue
10202 continue
10201 continue
      do 10211 kz = nzbdo,nzedo
      do 10212 jy = nybdo+1,nyedo-1
      do 10213 ix = nxbdo,nxedo
      at(ix,jy,kz,ia) = t0(ix,jy,kz)
10213 continue
10212 continue
10211 continue
12001 continue
c
c  z-direction
c
      do 13001 kk = 1,ismthz
c     print*,'SMOOTHZ',kk,ismth
      if ( ibase .eq. 0 ) then
      do 10301 kz = nzbdo+1,nzedo-1
      do 10302 jy = nybdo,nyedo
      do 10303 ix = nxbdo,nxedo
      t0(ix,jy,kz) =  0.25*(at(ix,jy,kz+1,ia)+at(ix,jy,kz-1,ia)
     >                     + 2.0*at(ix,jy,kz,ia))
10303 continue
10302 continue
10301 continue
      do 10311 kz = nzbdo+1,nzedo-1
      do 10312 jy = nybdo,nyedo
      do 10313 ix = nxbdo,nxedo
      at(ix,jy,kz,ia) = t0(ix,jy,kz)
10313 continue
10312 continue
10311 continue
      end if
      if ( ibase .eq. 1 ) then
      do 10391 kz = nzbdo+1,nzedo-1
      do 10392 jy = nybdo,nyedo
      do 10393 ix = nxbdo,nxedo
      t0(ix,jy,kz) =  0.25*(at(ix,jy,kz+1,ia)-ab(ix,jy,kz+1,ia)
     >                     +at(ix,jy,kz-1,ia)-ab(ix,jy,kz-1,ia)
     >               + 2.0*(at(ix,jy,kz,ia)  -ab(ix,jy,kz,ia)))
10393 continue
10392 continue
10391 continue
      do 10381 kz = nzbdo+1,nzedo-1
      do 10382 jy = nybdo,nyedo
      do 10383 ix = nxbdo,nxedo
      at(ix,jy,kz,ia) = t0(ix,jy,kz)+ab(ix,jy,kz,ia)
10383 continue
10382 continue
10381 continue
      end if
c
13001 continue
c
10010 continue
c
10000 continue
c
c 
      return
      end
c
c
c
      SUBROUTINE DRAWCL (XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of DRAWCL draws the polyline defined by the points
C ((XCS(I),YCS(I)),I=1,NCS) if and only if none of the area identifiers
C for the area containing the polyline are negative.  The dash package
C routine CURVED is called to do the drawing.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Turn on drawing.
C
        IDR=1
C
C If any area identifier is negative, turn off drawing.
C
        DO 101 I=1,NAI
          IF (IAI(I).LT.0) IDR=0
  101   CONTINUE
C
C If drawing is turned on, draw the polyline.
C
        IF (IDR.NE.0) CALL CURVED (XCS,YCS,NCS)
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDH1(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
c
c  common
c
c
      call getsh(chsval, chsang)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          IANG = chsang
          CALL SFSETI ('ANGLE',IANG)
          CALL SFSETR ('SPACING',chsval)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END 
c
c
c
      SUBROUTINE SHDH2(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          call getsh(chsval, chsang)
          IANG = chsang
          CALL SFSETI ('ANGLE',IANG)
          CALL SFSETR ('SPACING',chsval)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
          call getsh(chsval, chsang)
          IANG = chsang + 90
          CALL SFSETI ('ANGLE',IANG)
          CALL SFNORM (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDHR(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',45)
          CALL SFSETR ('SPACING',.006)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDHL(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',135)
          CALL SFSETR ('SPACING',.006)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDHH(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',90)
          CALL SFSETR ('SPACING',.006)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE SHDHV(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',0)
          CALL SFSETR ('SPACING',.006)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE SHDHX(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',45)
          CALL SFSETR ('SPACING',.006)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
          CALL SFSETI ('ANGLE',135)
          CALL SFNORM (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDDOT(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          call getsh(chsval, chsang)
          CALL SFSETR ('SP - SPACING OF FILL LINES',chsval)
          CALL SFSETI ('AN - ANGLE OF FILL LINES',0)
          CALL SFSETI ('DO - DOT-FILL FLAG',1)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDD6(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETR ('SP - SPACING OF FILL LINES',.006)
          CALL SFSETI ('AN - ANGLE OF FILL LINES',0)
          CALL SFSETI ('DO - DOT-FILL FLAG',1)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE SHDD4(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETR ('SP - SPACING OF FILL LINES',.004)
          CALL SFSETI ('AN - ANGLE OF FILL LINES',0)
          CALL SFSETI ('DO - DOT-FILL FLAG',1)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE SHDD2(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETR ('SP - SPACING OF FILL LINES',.002)
          CALL SFSETI ('AN - ANGLE OF FILL LINES',0)
          CALL SFSETI ('DO - DOT-FILL FLAG',1)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE SHDD8(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETR ('SP - SPACING OF FILL LINES',.008)
          CALL SFSETI ('AN - ANGLE OF FILL LINES',0)
          CALL SFSETI ('DO - DOT-FILL FLAG',1)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE SHDFIL(XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',45)
          CALL SFSETR ('SPACING',.006)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,1100,IND,1200)
          CALL SFSETI ('ANGLE',135)
          CALL SFNORM (XCS,YCS,NCS-1,DST,1100,IND,1200)
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE COLRAM (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
C
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
          IF (IFLL.GT.0.AND.IFLL.LT.15) THEN
            CALL GSFACI (IFLL+1)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
c
c
c
      SUBROUTINE SHDGRAY1 (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
C Second foreground color is gray
        CALL GSCR(1,21,.8,.8,.8)

        IF (ISH.NE.0) THEN
         CALL GSFACI(21)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END
c

      SUBROUTINE SHDBLUE1 (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
C Second foreground color is blue
        CALL GSCR(1,22,0.85,1.0,1.0)

        IF (ISH.NE.0) THEN
         CALL GSFACI(22)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END
c
      SUBROUTINE SHDBLUE2 (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
C Second foreground color is blue
        CALL GSCR(1,22,0.0,0.0,1.0)

        IF (ISH.NE.0) THEN
         CALL GSFACI(22)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END
c
      SUBROUTINE SHDRED1 (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
C Second foreground color is red
c        CALL GSCR(1,23,1.0,1.0,0.0)   ! yellow
c        CALL GSCR(1,23,1.0,0.5,0.5)     ! light pink
        CALL GSCR(1,23,1.0,1.0,0.5)     ! light yellow

        IF (ISH.NE.0) THEN
         CALL GSFACI(23)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END
c
      SUBROUTINE SHDRED2 (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
C Second foreground color is red
        CALL GSCR(1,23,1.0,0.0,0.0)

        IF (ISH.NE.0) THEN
         CALL GSFACI(23)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END
c
      SUBROUTINE SHDGREEN1 (XCS,YCS,NCS,IAI,IAG,NAI)
C                            (XCS,YCS,NCS,IAI,IAG,NAI)
       DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        DIMENSION DST(1100),IND(1200)
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE
C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
C Second foreground color is green
        CALL GSCR(1,23,0.0,1.0,0.0)

        IF (ISH.NE.0) THEN
         CALL GSFACI(23)
         CALL GFA (NCS,XCS,YCS)
        END IF

        RETURN
C
      END
c
c




c
c  subroutine getsh
c
      subroutine getsh(cv,ca)
      common /hsaw/  chsval, chsang
      cv = chsval
      ca = chsang
      print*,'CA,CV',ca,cv
      return
      end 
c
      subroutine test
      call getsh(ca,cv)
      print*,'TEST: CA,CV',ca,cv
      return
      end                 
c
c
c
      subroutine rastset(idx,ityp,nclr,iform,fac1,fac2,fac3,fac4
     >                  ,rvar,n1x,n2y,ivar,nlab,lab,ifill
     >                  ,rclo,rchi,rcin,idolab,ilablab,n1,n2)
      implicit none

      integer n,n1,n2,i1,i2,ityp,nlab,iform,nclr,idx,ifac
      integer n1x,n2y
      character*6 lab(nlab)
      integer idolab,ilablab
      real fac1,fac2,fac3,fac4
      real rclo,rchi,rcin
      real xvar,rvar(n1x,n2y)
      integer ivar(n1x,n2y),ifill(nlab)
c

         do i1 = 1,nlab
         ifill(i1) = i1+ Max(0, Nint(fac4) - 1)
         end do
      lab(:) = ''
      write(*,*) 'rastset: iform,nlab = ',iform,nlab
      write(*,*) 'rastset: rclo,rchi,rcin,idolab,ilablab = ',
     :  rclo,rchi,rcin,idolab,ilablab
      write(*,*) 'rastset: fac1,fac2,fac3,fac4 = ',fac1,fac2,fac3,fac4
      
      if ( iform .ge. 0 ) then
      do i2 = 1,n2
      do i1 = 1,n1
      ivar(i1,i2) = 0
      if ( rvar(i1,i2) .gt. rclo ) then
      ivar(i1,i2) = 0
      end if
      if ( rvar(i1,i2) .gt. rclo ) then
      xvar = fac1*(rvar(i1,i2)+fac2)/fac3 + fac4
      ivar(i1,i2) = ifix(xvar) ! + 220
      end if
      end do
c     print200,(ivar(i1,i2),i1=1,30)
c     print201,(rvar(i1,i2),i1=1,30)
200   format(1x,30i3)
201   format(1x,30f3.1)
      end do
      end if
c
      ifac = 0 ! fac4+.01
      do n = 1,nlab
      if ( ilablab .eq. 0 ) then
      write(lab(n),'(a)') ' '
      end if
      if ( ilablab .eq. 1 ) then
      write(lab(n+ifac),900) rclo + (n-1)*rcin
      end if
900   format(f6.2)
      end do

      DO n = 1,nlab
        IF ( lab(n) .ne. '' ) write(*,*) 'rastset: n,lab = ',n,lab(n)
      ENDDO
c
      return
      end
c
c
      subroutine dfclrs2(nclr,icrtyp,ishft)
      DIMENSION RGBV18(3,17)
C
C Define the RGB color triples for icrtpy=0 needed below.
C
       DATA RGBV18/ 1.000,  1.000,  1.000 ,
     +              0.142,  0.000,  0.850 ,
     +              0.097,  0.112,  0.970 ,
     +              0.160,  0.342,  1.000 ,
     +              0.240,  0.531,  1.000 ,
     +              0.340,  0.692,  1.000 ,
     +              0.460,  0.829,  1.000 ,
     +              0.600,  0.920,  1.000 ,
     +              0.740,  0.978,  1.000 ,
     +              1.000,  0.948,  0.740 ,
     +              1.000,  0.840,  0.600 ,
     +              1.000,  0.676,  0.460 ,
     +              1.000,  0.472,  0.340 ,
     +              1.000,  0.240,  0.240 ,
     +              0.970,  0.155,  0.210 ,
     +              0.850,  0.085,  0.187 ,
     +              0.650,  0.000,  0.130 /
C
      print*,'IN DFCLRS2',ishft
      if ( icrtyp .eq. 1 ) then
      if ( ishft .eq. 0 ) then
C
C Define NCLRS different color indices, for indices 0 through NCLR.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
C BACKGROUND COLOR
C BLACK
c        CALL GSCR(1,0,0.,0.,0.)
C First foreground color is white
c        CALL GSCR(1,1,1.,1.,1.)
C Second foreground color is gray
        CALL GSCR(1,2,.75,.75,.75)
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./NCLR
        DO  I=1,NCLR
          XHUE=I*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            CALL GSCR(1,NCLR+3-I,RED,GREEN,BLUE)
            ICNT=ICNT+1
          ELSE
            CALL GSCR(1,I-ICNT+2,RED,GREEN,BLUE)
          ENDIF
        END DO
c
      end if
      if ( ishft .eq. 1 ) then
C
C Define NCLRS different color indices, for indices 0 through NCLR.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
C BACKGROUND COLOR
C BLACK
        CALL GSCR(1,0,0.,0.,0.)
C First foreground color is white
        CALL GSCR(1,1,1.,1.,1.)
C Second foreground color is gray
c       CALL GSCR(1,2,.75,.75,.75)
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./NCLR
        DO I=1,NCLR
          XHUE=I*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            CALL GSCR(1,NCLR+2-I,RED,GREEN,BLUE)
            ICNT=ICNT+1
          ELSE
            CALL GSCR(1,I-ICNT+1,RED,GREEN,BLUE)
          ENDIF
        END DO
c
      end if
      if ( ishft .eq. 2 ) then
C
C Define NCLRS different color indices, for indices 0 through NCLR.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
C BACKGROUND COLOR
C BLACK
        CALL GSCR(1,0,0.,0.,0.)
C First foreground color is white
c       CALL GSCR(1,1,1.,1.,1.)
C Second foreground color is gray
c       CALL GSCR(1,2,.75,.75,.75)
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./NCLR
        DO I=1,NCLR
          XHUE=I*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            CALL GSCR(1,NCLR+1-I,RED,GREEN,BLUE)
            ICNT=ICNT+1
          ELSE
            CALL GSCR(1,I-ICNT+0,RED,GREEN,BLUE)
          ENDIF
        END DO
c
      end if
      end if
c
c
      if ( icrtyp .eq. 0 ) then
      if ( nclr .ne. 15 ) then
      write(*,*) 'nclr must be 15 for icrtyp = 0'
      stop
      end if
C
C Define 16 different color indices, for indices 0 through 15.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
c        CALL GSCR (1,0,0.,0.,0.)
C Dont reset colors 0 or 1 because they get redefined depending
c on the output type!  (e.g. ncgm or postscript)
        DO  I=2,17
          CALL GSCR (1,I,RGBV18(1,I),RGBV18(2,I),RGBV18(3,I))
        END DO
      end if
C
C Done.
C
        RETURN
C
      END  
c
c
c
c
c
      subroutine dfclrs(nclr,icrtyp,ishft)
      DIMENSION RGBV15(3,15)
C
C Define the RGB color triples for icrtpy=0 needed below.
C
        DATA RGBV15 / 1.00 , 1.00 , 1.00 ,
     +              0.70 , 0.70 , 0.70 ,
     +              0.75 , 0.50 , 1.00 ,
     +              0.50 , 0.00 , 1.00 ,
     +              0.00 , 0.00 , 1.00 ,
     +              0.00 , 0.50 , 1.00 ,
     +              0.00 , 1.00 , 1.00 ,
     +              0.00 , 1.00 , 0.60 ,
     +              0.00 , 1.00 , 0.00 ,
     +              0.70 , 1.00 , 0.00 ,
     +              1.00 , 1.00 , 0.00 ,
     +              1.00 , 0.75 , 0.00 ,
     +              1.00 , 0.38 , 0.38 ,
     +              1.00 , 0.00 , 0.38 ,
     +              1.00 , 0.00 , 0.00 /
C
      print*,'IN DFCLRS',ishft
      if ( icrtyp .eq. 1 ) then
      if ( ishft .eq. 0 ) then
C
C Define NCLRS different color indices, for indices 0 through NCLR.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
C BACKGROUND COLOR
C BLACK
c        CALL GSCR(1,0,0.,0.,0.)
C First foreground color is white
c        CALL GSCR(1,1,1.,1.,1.)
C Second foreground color is gray
        CALL GSCR(1,2,.75,.75,.75)
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./NCLR
        DO  I=1,NCLR
          XHUE=I*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            CALL GSCR(1,NCLR+3-I,RED,GREEN,BLUE)
            ICNT=ICNT+1
          ELSE
            CALL GSCR(1,I-ICNT+2,RED,GREEN,BLUE)
          ENDIF
        END DO
c
      end if
      if ( ishft .eq. 1 ) then
C
C Define NCLRS different color indices, for indices 0 through NCLR.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
C BACKGROUND COLOR
C BLACK
        CALL GSCR(1,0,0.,0.,0.)
C First foreground color is white
        CALL GSCR(1,1,1.,1.,1.)
C Second foreground color is gray
c       CALL GSCR(1,2,.75,.75,.75)
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./NCLR
        DO I=1,NCLR
          XHUE=I*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            CALL GSCR(1,NCLR+2-I,RED,GREEN,BLUE)
            ICNT=ICNT+1
          ELSE
            CALL GSCR(1,I-ICNT+1,RED,GREEN,BLUE)
          ENDIF
        END DO
c
      end if
      if ( ishft .eq. 2 ) then
C
C Define NCLRS different color indices, for indices 0 through NCLR.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
C BACKGROUND COLOR
C BLACK
        CALL GSCR(1,0,0.,0.,0.)
C First foreground color is white
c       CALL GSCR(1,1,1.,1.,1.)
C Second foreground color is gray
c       CALL GSCR(1,2,.75,.75,.75)
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./NCLR
        DO I=1,NCLR
          XHUE=I*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            CALL GSCR(1,NCLR+1-I,RED,GREEN,BLUE)
            ICNT=ICNT+1
          ELSE
            CALL GSCR(1,I-ICNT+0,RED,GREEN,BLUE)
          ENDIF
        END DO
c
      end if
      end if
c
c
      if ( icrtyp .eq. 0 ) then
      if ( nclr .ne. 15 ) then
      write(*,*) 'nclr must be 15 for icrtyp = 0'
      stop
      end if
C
C Define 16 different color indices, for indices 0 through 15.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
c        CALL GSCR (1,0,0.,0.,0.)
C Dont reset colors 0 or 1 because they get redefined depending
c on the output type!  (e.g. ncgm or postscript)
        DO  I=2,15
          CALL GSCR (1,I,RGBV15(1,I),RGBV15(2,I),RGBV15(3,I))
        END DO
      end if
C
C Done.
C
        RETURN
C
      END  


c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE STRLNTH                    ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE STRLNTH( string, length )
c
c#######################################################################
c      
c     PURPOSE:
c
c     Return the length of the non-blank part of a character string.
c
c#######################################################################
c
c     AUTHOR: Ming Xue
c     11/20/91
c
c     MODIFICATION HISTORY:
c
c     5/05/92 (M. Xue)
c     Added full documentation.
c
c     9/10/94 (Weygandt & Y. Lu)
c     Cleaned up documentation.
c
c#######################################################################
c
c     INPUT:
c
c       string   A character string
c       length   The declared length of the character string 'string'.
c
c     OUTPUT:
c
c       length   The length of the non-blank part of the string.
c
c#######################################################################
c
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      character string*(*)      ! A character string for the name of
                                ! this run.
      integer length            ! The length of the non-blank part
                                ! of a string.

      integer i
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c

      DO 100 i = length,1,-1

        IF(string(i:i) .ne. ' ') GOTO 200

100   CONTINUE

200   CONTINUE

      length = max(1,i) 

      RETURN
      END


C              SUBROUTINE LBFILL (IFTP,XCRA,YCRA,NCRA,INDX)
C                DIMENSION XCRA(*),YCRA(*)
C        integer ifll0
C        common/brparam/ ifll0
C                i = 0
C                IF ( ifll0 .eq. 0 ) THEN
C                  CALL GSFACI (INDX)
C                ELSE
C                  i = ifll0+2 ! (indx - ifll0 + 1)*2 - 1
C!                  IF ( i .le. 0 ) i = 19
C                  CALL GSFACI (i)
C                  ifll0 = ifll0 + 1
C                ENDIF
C                write(0,*) 'lbfill: indx = ',indx,i
C                CALL GFA (NCRA-1,XCRA,YCRA)
C                RETURN
C              END

