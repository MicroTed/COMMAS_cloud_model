c############################################################################
c
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE SNDPARAMS                   ######
c     ######                                                      ######
c     ##################################################################
c
c
c     PURPOSE:
c
c     Calculates a handful of sounding parameters from the ensemble of
c     initial states and outputs them to a file called 'sndparams.out'.
c
c############################################################################
c
c     Author:  Mike Coniglio
c
c     Creation Date:  04 October 2005
c
c############################################################################

      subroutine sndparams(nd)

#ifdef NAG
      USE F90_UNIX_PROC
#endif

c** nd is the number of lines in the sounding files
c** nz is INTERNAL and is the number of vertical levels in the interpolation arrays

      integer nd
      integer,parameter :: nz=25 
      real   ,parameter :: nodata=-9999.9

      integer i,j,k,ptype,klvl
      real, dimension(nd) :: hgt,prs,u,v,tmk,tdk
      real, dimension(nz) :: zlvl,uz,vz

c** Output from parcel subroutine
      real, dimension(nz) :: penergy
      real ca,ci,plpl,zlpl,plcl,zlcl,plfc,zlfc,pel,zel

c** Output from zinterp routine
      real, dimension(nz) :: p_out

c** Output from shearz subroutine
      real, dimension(0:12,0:12) :: bulkshz,totalshz
      real, dimension(nz-1) :: zshz,ushz,vshz
      real effshz,clyrshz,sreh01,sreh03

c** Newly defined variables
      real mlcape,mlcin,mllcl,mllfc,dcape

c** stuff for reading in data files 
      integer n_sndgs
      integer,parameter :: n_max_sndgs=500
      character(LEN=140) command
      character(LEN=80) snd_file(0:n_max_sndgs) 

c** Create a text file that contains a list of the sounding files

      call system('rm profile_file_list.txt')
      write(command,"(A)") 'ls -1 profile* > profile_file_list.txt'
      call system(command)

c** Read in the names of the files

      open(unit=30, file='profile_file_list.txt', status='old')
      n_sndgs = 0
      do
       if(n_sndgs.gt.n_max_sndgs) then
        write(*,*) 'too many soundings, max=', n_max_sndgs
        stop
       endif
       read(30,"(A)",end=999) snd_file(n_sndgs)
       n_sndgs = n_sndgs + 1
      enddo
 999  close(30)

c** Start the main loop for reading in the profiles

      main: do i=0,n_sndgs-1

       open(30,file=snd_file(i))
       read(30,fmt='(65x)')
       read(30,fmt='(66x)')
       do k=1,nd
        read(30,100,end=998) klvl,hgt(k),prs(k),u(k),v(k),tmk(k),tdk(k)
 100    format(I3,3X,F8.1,1X,F8.1,1X,F6.1,5X,F6.1,5X,F6.1,5X,F6.1)
       enddo
 998   close(30)

c**  ptype=1 --> SURFACE PARCEL
c**  ptype=2 --> MOST UNSTABLE SINGLE PARCEL IN LOWEST 400 MB
c**  ptype=3 --> SURFACE LAYER MIXED PARCEL
c**  ptype=4 --> MOST UNSTABLE MIXED PARCEL IN LOWEST 400 MB
c**  ptype=5 --> DOWNDRAFT PARCEL FROM HEIGHT AT MINIMUM THETA-W
c**              OR IN LOWEST 500 MB OR FROM WET-BULB ZERO HEIGHT
c**  ptype=6 --> RUC-2 "BEST" CAPE: IF MAX THETA-E PARCEL
c**              IS IN LOWEST 70 MB, MIX THE LOWEST 70 MB TO FIND THE
c**              PARCEL THETA-E, OTHERWISE USE THE SINGLE MOST
c**              UNSTABLE PARCEL

       ptype=4   ! (MLCAPE)
       call parcel(ptype,nodata,nd,nd,tmk,tdk,prs,hgt
     $  ,penergy
     $  ,ca
     $  ,ci
     $  ,plpl
     $  ,zlpl
     $  ,plcl
     $  ,zlcl
     $  ,plfc
     $  ,zlfc
     $  ,pel
     $  ,zel
     $   )

       mlcape=ca
       mlcin=ci
       mllcl=zlcl
       mllfc=zlfc

       ptype=5   ! (DCAPE)
       call parcel(ptype,nodata,nd,nd,tmk,tdk,prs,hgt
     $  ,penergy
     $  ,ca
     $  ,ci
     $  ,plpl
     $  ,zlpl
     $  ,plcl
     $  ,zlcl
     $  ,plfc
     $  ,zlfc
     $  ,pel
     $  ,zel
     $   )

       dcape=ca

c** set up an array of height levels for interpolation

       do j=1,nz
        zlvl(j)=0.0+(500.0*float(j-1))
       enddo

c** interpolate u and v to 0,500,1000,....,12000 m levels

       call zinterp(nd,nd,nodata,u,v,prs,hgt,nz,zlvl
     $  ,uz,vz,p_out)

       call shearz(nodata,nz,zlvl,uz,vz,zlpl,zlcl,zel
     $  ,bulkshz
     $  ,zshz
     $  ,ushz
     $  ,vshz
     $  ,totalshz
     $  ,effshz
     $  ,clyrshz
     $  ,sreh01
     $  ,sreh03
     $   )

       if(i.eq.0) then
        open(40,file='sndparams.out')
        write(40,*) 'MEM  CAPE    CIN  DCAPE   ZLCL   ZLFC   0-1S   0-6S  0-10S 0-1SRH 0-3SRH'
        write(40,*) '------------------------------------------------------------------------' 
       endif

        write(40,101) i,int(mlcape),int(mlcin),int(dcape),int(mllcl),int(mllfc),bulkshz(0,1)
     $  ,bulkshz(0,6),bulkshz(0,10),sreh01,sreh03
 101    format(I3,5(2X,I5),5(2X,F5.1))

      end do main

      close(40)

      end subroutine sndparams

c**********************************************************************
      subroutine parcel(ptype,nodata,nd,mcnt,tdry,tdew,prst,hgtt
     $  ,penergy
     $  ,ca
     $  ,ci
     $  ,plpl
     $  ,zlpl
     $  ,plcl
     $  ,zlcl
     $  ,plfc
     $  ,zlfc
     $  ,pel
     $  ,zel
     $   )

c**      REQUIRES:  TEMPERATURE IN KELVIN (tdry)
c**                 DEWPOINT IN KELVIN (tdew) - NEED GOOD DATA IN LOWEST 400 MB!!
c**                 PRESSURE IN MB (prst)
c**                 HEIGHT IN M (hgtt)

      implicit none

c**  INPUT PARAMETERS AND ARRAYS
      integer ptype,nd,mcnt
      real nodata
      real, dimension(nd) :: tdry,tdew,prst,hgtt

c**  INTERNAL PARAMETERS
      real dpmpl

c**  INTERNAL VARIABLES AND ARRAYS
      integer i,j,k,nmpl,jthwmax,jthwmin,jthwbz,jdcape
      real pfact,tfact,plvl,tlcl,vap,thwmax,thwmin,aos,tvp,tavg,dz
      real tvp1,tvp2,pwbz,zwbz,zthwmin
      real, dimension(mcnt) :: p,z,tk,dk,q,tv,temp,qmix,pmix,pbot
     $    ,avgcount,qavgcount,pl,tl,tha,the,thw,penergy

c**  INTERNAL FUNCTIONS
      real epres_td,wmix,os,tua,tsa

c**  OUTPUT VARIABLES
      real ca,ci,plpl,zlpl,plcl,zlcl,plfc,zlfc,pel,zel

c**  INITIALIZE THE PRESSURE DIFFERENCE IN MIXED PARCEL LAYERS
      if(ptype.le.5) dpmpl=50.0
      if(ptype.eq.6) dpmpl=70.0  ! Used for RUC-2 BEST CAPE

c**  INITIALIZE ARRAYS
      temp(:)=0.
      qmix(:)=0.
      pmix(:)=0.
      pbot(:)=0.
      avgcount(:)=0.
      qavgcount(:)=0.

c**  MAKE SURE PRESSURES ARE IN MB AND TEMPS ARE IN KELVIN

       if(prst(1).ge.80000.) then
         pfact=0.01
       else
         pfact=1.0
       endif

       if(tdry(1).le.60.) then
         tfact=273.16
       else
         tfact=0.0
       endif

       do k=1,mcnt
        p(k)=pfact*prst(k)
        z(k)=hgtt(k)
        tk(k)=tdry(k)
c        p(k)=pfact*prst(mcnt-k+1)
c        z(k)=hgtt(mcnt-k+1)
c        tk(k)=tdry(mcnt-k+1)
        if(tdew(k).ne.nodata) then
         dk(k)=tdew(k)+tfact
        else
         dk(k)=nodata
        endif
c        if(tdew(mcnt-k+1).ne.nodata) then
c         dk(k)=tdew(mcnt-k+1)+tfact
c        else
c         dk(k)=nodata
c        endif
c        print *,p(k),z(k),tk(k),dk(k)
       enddo   

c**  CALCULATE MIXING RATIO AND VIRTUAL TEMP AT EACH LEVEL

       do k=1,mcnt
        if(dk(k).ne.nodata) then
         vap=epres_td(dk(k))  ! output in mb
         q(k)=wmix(p(k),vap)  ! output in g/kg
         tv(k)=tk(k)*(1.+0.000608*q(k))
c         print *,q(k),tv(k)
        else
         q(k)=nodata
         tv(k)=tk(k)
        endif
       enddo

c#############################################
c Section to compute theta-w for parcels
c#############################################
c**  SET UP THE NUMBER OF PRESSURE LEVELS TO LOOP
c**  THROUGH ACCORDING TO THE TYPE OF PARCEL SPECIFIED
c**  IN THE CALL
c**
c**  ptype=1 --> SURFACE PARCEL
c**  ptype=2 --> MOST UNSTABLE SINGLE PARCEL IN LOWEST 400 MB
c**  ptype=3 --> SURFACE LAYER MIXED PARCEL
c**  ptype=4 --> MOST UNSTABLE MIXED PARCEL IN LOWEST 400 MB
c**  ptype=5 --> DOWNDRAFT PARCEL FROM HEIGHT AT MINIMUM THETA-W
c**              OR IN LOWEST 500 MB OR FROM WET-BULB ZERO HEIGHT
c**  ptype=6 --> CAPE AS CALCULATED IN RUC-2: IF MAX THETA-E PARCEL
c**              IS IN LOWEST 70 MB, MIX THE LOWEST 70 MB TO FIND THE
c**              PARCEL THETA-E, OTHERWISE USE THE SINGLE MOST
c**              UNSTABLE PARCEL

       nmpl=0

       select case (ptype)
        case(1)
            nmpl=1
            pbot(nmpl)=p(1)
        case(2)
            do k=1,mcnt
             if((p(1)-p(k)).le.400.0) then
              nmpl=nmpl+1
              pbot(nmpl)=p(k)
             endif
            enddo
        case(3)
            nmpl=2
            do j=1,nmpl
             pbot(j)=p(1)-(float(j-1)*dpmpl)
            enddo
        case(4)
            do k=1,mcnt
             plvl=p(1)-(float(k-1)*dpmpl)
             if((p(1)-plvl).le.400.0) then
              nmpl=nmpl+1
              pbot(nmpl)=plvl
             endif
            enddo
        case(5)             ! Used for finding wet-bulb zero and theta-w min below
            do k=1,mcnt
             if(p(k).ge.500.0) then
              nmpl=nmpl+1
              pbot(nmpl)=p(k)
             endif
            enddo
        case(6)
            do k=1,mcnt
             if(p(k).ge.500.0) then
              nmpl=nmpl+1
              pbot(nmpl)=p(k)
             endif
            enddo
       end select

c**  IF WE'RE ONLY TAKING A SINGLE PARCEL (PTYPE=1 OR 2), THEN
c**  YOU WON'T NEED TO SUBTRACT 1 FROM THE FOLLOWING 
c**  LOOP COUNTERS (YOU AREN'T AVERAGING ANYTHING) SO SET i=0

       select case (ptype)

        case(1,2)

         i=0
         do j=1,nmpl-i
          do k=1,mcnt
           if( p(k).eq.pbot(j) ) then
            temp(j)=tk(k)
            qmix(j)=q(k)
            pmix(j)=p(k)
            avgcount(j)=1.0
            qavgcount(j)=1.0
           endif
          enddo
         enddo

        case(3,4)

         i=1
         do j=1,nmpl-i
          do k=1,mcnt
           if( p(k).eq.pbot(j) ) then
            temp(j)=temp(j)+tk(k)
            if(q(k).ne.nodata) then
             qmix(j)=qmix(j)+q(k)
             qavgcount(j)=qavgcount(j)+1.0
            endif
            pmix(j)=pmix(j)+p(k)
            avgcount(j)=avgcount(j)+1.0
           elseif( p(k).ge.pbot(j+1).and.p(k).lt.pbot(j) ) then
            temp(j)=temp(j)+tk(k)
            if(q(k).ne.nodata) then
             qmix(j)=qmix(j)+q(k)
             qavgcount(j)=qavgcount(j)+1.0
            endif
            pmix(j)=pmix(j)+p(k)
            avgcount(j)=avgcount(j)+1.0
           endif
          enddo
         enddo

        case(5)

         i=0
         do j=1,nmpl-i
          do k=1,mcnt
           if( p(k).eq.pbot(j) ) then
            temp(j)=tk(k)
            qmix(j)=q(k)
            pmix(j)=p(k)
            avgcount(j)=1.0
            qavgcount(j)=1.0
           endif
          enddo
         enddo

        case(6)

         i=0
         do j=1,nmpl-i
          do k=1,mcnt
           if( p(k).eq.pbot(j) ) then
            temp(j)=tk(k)
            qmix(j)=q(k)
            pmix(j)=p(k)
            avgcount(j)=1.0
            qavgcount(j)=1.0
           endif
          enddo
         enddo

       end select

c**  FIND THE EQUIVALENT WET-BULB POTENTIAL TEMP FOR EACH PARCEL
c**  (THERE MAY NOT ALWAYS BE DATA EVERY 'DPMPL' MB, SO MAKE SURE
c**   YOU SET THETA-W FOR THESE LEVELS EQUAL TO THE BAD VALUE)

       do j=1,nmpl-i
        if( avgcount(j).ne.0.0 ) then
         temp(j)=temp(j)/avgcount(j)
         if(qavgcount(j).ne.0.0) then
          qmix(j)=qmix(j)/qavgcount(j)
         else
          qmix(j)=nodata
         endif
         pmix(j)=pmix(j)/avgcount(j)
              if(qmix(j).ne.nodata) then
                call parcelthw(temp(j),qmix(j),pmix(j),
     $             pl(j),tl(j),tha(j),the(j),thw(j))
              else
                thw(j)=nodata
              endif
        else
         thw(j)=nodata
        endif
c        print *,thw(j)
       enddo

c#######################################################################
c  end of section to compute lifted parcel level, LCL temp., theta,
c  theta-e and theta-w for parcels
c#######################################################################

       cape_or_dcape: if(ptype.le.4.or.ptype.eq.6) then

c**  FIND THE MOST UNSTABLE PARCEL (PARCEL WITH LARGEST THETA-W)

        thwmax=nodata
        jthwmax=0
        do j=1,nmpl-i
         if(thw(j).gt.thwmax) then
           thwmax=thw(j)
           jthwmax=j
         endif
        enddo

        if(jthwmax.eq.0) then
         print *,'No good parcels for this sounding - must remove it'
         stop
        endif

c**  IF YOU WANT THE RUC-2 CAPE, WE NEED TO MIX THE PARCEL IF THE
c**  THETA-W MAX IS IN THE LOWEST 70 MB

        ruc_cape: if(ptype.eq.6) then

         if( (pmix(1)-pmix(jthwmax)).le.70.0 ) then
           ! Reset the mixed parcel arrays
             temp(:)=0.
             qmix(:)=0.
             pmix(:)=0.
             avgcount(:)=0.
             qavgcount(:)=0.
           ! Set the lower and upper boundaries for the mixed layer
             nmpl=2
             do j=1,nmpl
              pbot(j)=p(1)-(float(j-1)*dpmpl)
             enddo
           ! there's only one layer to mix, so set j=1
             j=1
           ! find the mixed parcel in the lowest 70 mb
             do k=1,mcnt
               if( p(k).eq.pbot(j) ) then
                temp(j)=temp(j)+tk(k)
                if(q(k).ne.nodata) then
                 qmix(j)=qmix(j)+q(k)
                 qavgcount(j)=qavgcount(j)+1.0
                endif
                pmix(j)=pmix(j)+p(k)
                avgcount(j)=avgcount(j)+1.0
               elseif( p(k).ge.pbot(j+1).and.p(k).lt.pbot(j) ) then
                temp(j)=temp(j)+tk(k)
                if(q(k).ne.nodata) then
                 qmix(j)=qmix(j)+q(k)
                 qavgcount(j)=qavgcount(j)+1.0
                endif
                pmix(j)=pmix(j)+p(k)
                avgcount(j)=avgcount(j)+1.0
               endif
             enddo
           ! Finally, find the new theta-w of this mixed parcel
             j=1
             if( avgcount(j).ne.0.0 ) then
              temp(j)=temp(j)/avgcount(j)
              if(qavgcount(j).ne.0.0) then
               qmix(j)=qmix(j)/qavgcount(j)
              else
               qmix(j)=nodata
              endif
              pmix(j)=pmix(j)/avgcount(j)
                   if(qmix(j).ne.nodata) then
                     call parcelthw(temp(j),qmix(j),pmix(j),
     $                  pl(j),tl(j),tha(j),the(j),thw(j))
                   else
                     thw(j)=nodata
                   endif
             else
              thw(j)=nodata
             endif
             jthwmax=1
         endif

        end if ruc_cape

c**  DEFINE THE STARTING PARCEL PRESSURE LEVEL, THE LCL PRESSURE, AND THE
c**  POTENTIAL TEMPERATURE OF THE PARCEL DURING UNSATURATED ASCENT

       plpl=pmix(jthwmax)
       plcl=pl(jthwmax)
       tlcl=tl(jthwmax)

c**  COMPUTE THE CAPE AND CIN FOR THE PARCEL

       ca=0.
       ci=0.
       penergy(:)=nodata
       do k=1,mcnt-1
        if(p(k).le.plpl) then
           if(p(k).ge.plcl) then
            tvp=tua(tha(jthwmax),p(k),q(k))  ! temp (K) of parcel for unsaturated ascent
           else
            tvp=tua(tha(jthwmax),p(k),0.0)
           endif
           if(p(k).lt.plcl) tvp=tsa(the(jthwmax),p(k),tlcl)  ! temp (K) of parcel for saturated ascent
           if(q(k).ne.nodata) then
            tvp=tvp*(1.+0.000608*q(k))  ! restore virtual temperature correction to parcel temp
           endif
           tavg=(tv(k+1)+tv(k))/2.
           dz=((-287.04*tavg)/9.8)*alog(p(k+1)/p(k))
           penergy(k)=9.8*((tvp-tv(k))*dz)/tv(k)
c           print *,p(k),penergy(k)
           if(penergy(k).lt.0.0.and.p(k).gt.500.) then
            ci=ci+penergy(k)
           elseif(penergy(k).ge.0.) then
            ca=ca+penergy(k)
           endif
        endif
       enddo

c**  FIND THE PRESSURE AT THE LFC FOR THE PARCEL

       if(ca.le.0.1) then
        plfc=nodata
       else
            plfc=plcl
            do k=1,mcnt-2
             if(p(k).le.plpl.and.p(k).ge.500.0) then
                if(penergy(k).le.0.0.and.penergy(k+1).gt.0.0) then
                 call linterpmike(penergy(k),alog(p(k)),penergy(k+1)
     $             ,alog(p(k+1)),0.0,plfc)
                 plfc=exp(plfc)
                 if(plfc.gt.plcl) plfc=plcl
                endif
             endif
            enddo
       endif

c**  FIND THE PRESSURE AT THE EQUILIBRIUM LEVEL FOR THE PARCEL

       if(plfc.eq.nodata) then
        pel=nodata
       else
            do k=1,mcnt-2
             if(p(k).le.plfc) then
                if(penergy(k).gt.0.0.and.penergy(k+1).le.0.0) then
                 call linterpmike(penergy(k),alog(p(k)),penergy(k+1)
     $             ,alog(p(k+1)),0.0,pel)
                 pel=exp(pel)
                 goto 91
                endif
                pel=p(mcnt-1)
             endif
            enddo
       endif
 91    continue

c**  FIND THE HEIGHT AT THE LIFTED PARCEL LEVEL, THE LCL,LFC, AND EL

       do k=1,mcnt-1
        if(p(k).ge.plpl.and.p(k+1).le.plpl) then
          call linterpmike( alog(p(k)),z(k),alog(p(k+1)),z(k+1),
     $                  alog(plpl),zlpl )
        goto 92
        endif
       enddo
 92    continue

       do k=1,mcnt-1
        if(p(k).ge.plcl.and.p(k+1).le.plcl) then
          call linterpmike( alog(p(k)),z(k),alog(p(k+1)),z(k+1),
     $                  alog(plcl),zlcl )
        goto 93
        endif
       enddo
 93    continue

       zlfc=nodata
       do k=1,mcnt-1
        if(plfc.ne.nodata.and.p(k).ge.plfc.and.p(k+1).le.plfc) then
          call linterpmike( alog(p(k)),z(k),alog(p(k+1)),z(k+1),
     $                  alog(plfc),zlfc )
        goto 94
        endif
       enddo
 94    continue

       zel=nodata
       do k=1,mcnt-1
        if(pel.ne.nodata.and.p(k).ge.pel.and.p(k+1).le.pel) then
          call linterpmike( alog(p(k)),z(k),alog(p(k+1)),z(k+1),
     $                  alog(pel),zel )
        goto 95
        endif
       enddo
 95    continue

      elseif(ptype.eq.5) then     ! cape_or_dcape construct

c** FIND THE WET-BULB ZERO HEIGHT AND USE THE GREATER OF THAT HEIGHT
c** AND THE STARTING HEIGHT OF THE PARCEL WITH THE MINIMUM THETA-W
c** ABOVE 1 KM AS
c** THE STARTING HEIGHT FOR COMPUTING THE DOWNDRAFT CAPE

       jthwbz=0
       do j=1,(nmpl-1)-i
        if(thw(j).ne.nodata) then
         do k=j+1,nmpl-i
          if(thw(k).ne.nodata) then
           tvp1=tsa(the(j),pmix(j),tl(j))
           tvp2=tsa(the(k),pmix(k),tl(k))
            if(tvp1.ge.273.16.and.tvp2.lt.273.16) then
               jthwbz=k
               goto 96
            endif
          endif
         enddo
        endif
       enddo
 96    continue

       thwmin=9999.9
       jthwmin=0
       do j=3,nmpl-i
        if(thw(j).lt.thwmin) then
         thwmin=thw(j)
         jthwmin=j
        endif
       enddo

       if(jthwmin.eq.0) then
        print *,'Cant find a minimum theta-w - check your sounding'
        stop
       endif

       if(jthwbz.eq.0) then
        pwbz=nodata
        zwbz=nodata
        goto 97
c        print *,'Cant find a wet-bulb zero temp - check your sounding'
c        stop
       endif

c** FIND THE WET-BULB ZERO HEIGHT

       call linterpmike(tvp1,alog(pmix(jthwbz-1)),tvp2,alog(pmix(jthwbz)),
     $               273.16,pwbz)
       pwbz=exp(pwbz)

       do k=1,mcnt-1
        if(p(k).ge.pwbz.and.p(k+1).le.pwbz) then
         call linterpmike( alog(p(k)),z(k),alog(p(k+1)),z(k+1),
     $                  alog(pwbz),zwbz )
          goto 97
        endif
       enddo
 97    continue

c** FIND THE STARTING HEIGHT OF THE DOWNDRAFT PARCEL

       do k=1,mcnt-1
        if(p(k).ge.pmix(jthwmin).and.p(k+1).le.pmix(jthwmin)) then
          call linterpmike( alog(p(k)),z(k),alog(p(k+1)),z(k+1),
     $                  alog(pmix(jthwmin)),zthwmin )
          goto 98
        endif
       enddo
 98    continue

c**  DEFINE THE STARTING PARCEL LEVEL

       if(zwbz.gt.zthwmin) then
        jdcape=jthwbz
       else
        jdcape=jthwmin
       endif

       plpl=pmix(jdcape)
       tlcl=tl(jdcape)

c**  COMPUTE THE DOWNDRAFT CAPE

       ci=0.
       ca=0.
       penergy(:)=nodata
       do k=mcnt,2,-1
        if(p(k).ge.plpl) then
           tvp=tsa(the(jdcape),p(k),tlcl)  ! temp (K) of parcel for saturated ascent
           tvp=tvp*(1.+0.000608*q(k))      ! restore virtual temperature correction to parcel temp
           tavg=(tv(k)+tv(k-1))/2.
           dz=((-287.04*tavg)/9.8)*alog(p(k)/p(k-1))
           penergy(k)=9.8*((tv(k)-tvp)*dz)/tv(k)
           if(penergy(k).lt.0.0) then
            ci=ci+penergy(k)
           else
            ca=ca+penergy(k)
           endif
        endif
       enddo

       ci=0.

       zlpl=max(zthwmin,zwbz)
       zlcl=zthwmin
       zlfc=zwbz

       plcl=nodata
       plfc=nodata
       pel=nodata
       zel=nodata

      end if cape_or_dcape

      RETURN
      END SUBROUTINE parcel
c*********************************************************************
      SUBROUTINE parcelthw(t,q,p,pl,tl,tha,the,thw)

      real,parameter :: cp=1004.
      real,parameter :: r=287.04
      real t,q,p                                 ! t in K, q in g/kg
      real e,tl,dzdry,dpdry,tha,pl,the,thw       ! tha,the in K, thw in C
      real epres,tlcl1,theta,theta_e,theta_w     ! Functions

      e=epres(p,q)                               ! vapor pressure in mb
      tl=tlcl1(t,e)                              ! LCL temp. in K
      dzdry=(t-tl)/0.00977                       ! Height diff. between LCL and parcel level
      dpdry=(-9.801*p*(t-tl)) / (0.00977*r*t)    ! Pressure diff. between LCL and parcel level
      tha=theta(p,t,q)                           ! dry air potential temp. in K (accounts for water vapor)
      pl=1000.*((tl/tha)**(cp/r))                ! LCL Pressure in mb
      the=theta_e(tha,tl,q)                      ! Equivalent potential temp. in K
      thw=theta_w(the)                           ! Equivalent wet-bulb potential temp. in C

      RETURN
      END SUBROUTINE parcelthw
c********************************************************************
      SUBROUTINE shearz(nodata,nz,z,u,v,zlpl,zlcl,zel
     $  ,bulkshz
     $  ,zshz
     $  ,ushz
     $  ,vshz
     $  ,totalshz
     $  ,effshz
     $  ,clyrshz
     $  ,sreh01
     $  ,sreh03
     $   )

      implicit none

c** INPUT VARIABLES

      integer nz
      real nodata,zlpl,zlcl,zel
      real, dimension(nz) :: z,u,v

c** INTERNAL VARIABLES

      integer j,j1,j2,l1,l2,k,intshz,j1000,j3000
      real zl1,zl2,shseg,ustorm,vstorm
      real ulpl,vlpl,ul1,vl1,ul2,vl2,sum01,sum03
      real bulkshear  ! FUNCTION

c** OUTPUT VARIABLES

      real, dimension(0:12,0:12) :: bulkshz,totalshz,um,vm
      real, dimension(nz-1) :: dz,zshz,ushz,vshz
      real effshz,clyrshz,sreh01,sreh03,dir0to6,spd0to6

c** CALCULATE AN ARRAY OF BULK SHEAR VALUES IN M/S.
c** THIS WILL CREATE A 13 x 13 SYMMETRIC MARTIX
c** WITH THE ROWS AND COLUMNS DEFINED BY THE HEIGHTS
c** IN KM.  THE DUPLICATE ENTRIES WILL BE DEFINED AS
c** BAD VALUES

      bulkshz(:,:)=nodata
      do l1=0,12
       do l2=0,12
         if(l1.le.l2) then
          zl1=1000.*float(l1)
          zl2=1000.*float(l2)
          do j=1,nz
           if(z(j).eq.zl1) j1=j
           if(z(j).eq.zl2) j2=j
          enddo
          if(u(j1).ne.nodata.and.v(j1).ne.nodata.and.
     $       u(j2).ne.nodata.and.v(j2).ne.nodata)     then
            bulkshz(l1,l2)=bulkshear(u(j1),u(j2),v(j1),v(j2))
          endif
          bulkshz(l2,l1)=bulkshz(l1,l2)  ! Fills in matrix
         endif
       enddo
      enddo

c      do l1=0,12
c       write(*,fmt='(13f8.1)') (bulkshz(l1,l2), l2=0,12)
c      enddo

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c** 10/04/05: These sections below work fine, but I commented
c**           them out because I'm not using any of these variables
c**
c** CREATE AN ARRAY OF HEIGHT DIFFERENCES (IN KM) AND THE HEIGHT
c** LEVELS FOR THE SHEAR CALCULATIONS
c
c      do j=1,nz-1
c       dz(j)=0.001*(z(j+1)-z(j))
c       zshz(j)=(z(j+1)+z(j))/2.0
c      enddo
c
c** CALCULATE U AND V SHEAR PROFILES IN M/(S*KM)
c
c      ushz(:)=nodata
c      vshz(:)=nodata
c
c      do j=1,nz-1
c       if(u(j).ne.nodata.and.u(j+1).ne.nodata.and.
c     $    v(j).ne.nodata.and.v(j+1).ne.nodata)     then
c        ushz(j)=(u(j+1)-u(j))/dz(j)
c        vshz(j)=(v(j+1)-v(j))/dz(j)
c       endif
c      enddo
c  
c** CALCULATE THE TOTAL (OR CUMULATIVE) SHEAR IN 1 KM LAYERS
c** UP TO 12 KM USING EVERY LEVEL
c
c      intshz=1
c      totalshz(:,:)=nodata
c
c      do l1=0,12
c       do l2=0,12
c
c        if(l1.eq.l2) then
c         totalshz(l1,l2)=0.0
c         totalshz(l2,l1)=totalshz(l1,l2)
c        endif
c
c        if(l1.lt.l2) then
c         zl1=1000.*float(l1)
c         zl2=1000.*float(l2)
c         do j=1,nz
c          if(zl1.eq.z(j)) j1=j
c          if(zl2.eq.z(j)) j2=j
c         enddo
c         shseg=0.0
c         do k=j1,j2-intshz,intshz
c          if(u(k).ne.nodata.and.u(k+intshz).ne.nodata.and.
c     $       v(k).ne.nodata.and.v(k+intshz).ne.nodata)     then
c           totalshz(l1,l2)=shseg+bulkshear(u(k),u(k+intshz)
c     $       ,v(k),v(k+intshz))
c           shseg=totalshz(l1,l2)
c          endif
c          totalshz(l2,l1)=totalshz(l1,l2)   !  Fills in matrix
c         enddo
c        endif
c
c       enddo
c      enddo
c
c      do l1=0,12
c       write(*,fmt='(13f5.1)') (totalshz(l1,l2), l2=0,12)
c      enddo
c
c** CALCULATE THE EFFECTIVE SHEAR (SEE THOMPSON ET AL. 2004 22nd SLS)
c** THIS IS THE BULK SHEAR BETWEEN THE LIFTED PARCEL LEVEL OF THE
c** MOST UNSTABLE PARCEL AND THE LEVEL THAT IS 50% OF THE HEIGHT BETWEEN
c** THE LIFTED PARCEL LEVEL AND THE EQUILIBRIUM LEVEL
c
c      effshz=nodata
c
c      if(zel.eq.nodata) then
c       zel=12000.0
c      endif
c
c      zl1=zlpl
c      zl2=(zel/2.)
c
c      do j=1,nz-1
c       if(z(j).le.zl1.and.z(j+1).ge.zl1) then
c          if(u(j).ne.nodata.and.u(j+1).ne.nodata.and.
c     $       v(j).ne.nodata.and.v(j+1).ne.nodata)     then
c           call linterpmike(z(j),u(j),z(j+1),u(j+1),zl1,ul1)
c           call linterpmike(z(j),v(j),z(j+1),v(j+1),zl1,vl1)
c          endif
c       elseif(z(j).le.zl2.and.z(j+1).ge.zl2) then
c          if(u(j).ne.nodata.and.u(j+1).ne.nodata.and.
c     $       v(j).ne.nodata.and.v(j+1).ne.nodata)     then
c           call linterpmike(z(j),u(j),z(j+1),u(j+1),zl2,ul2)
c           call linterpmike(z(j),v(j),z(j+1),v(j+1),zl2,vl2)
c          endif
c       endif
c      enddo
c
c      effshz=bulkshear(ul1,ul2,vl1,vl2)
c      print *,zl1,zl2,effshz
c
c 90   continue
c
c** FIND THE BULK SHEAR BETWEEN THE LCL AND THE EL OR 10 KM
c** (WHICH EVER IS LOWER) AND USE THAT TO DEFINE THE CLOUD 
c** LAYER SHEAR 
c
c      clyrshz=nodata
c
c      if(zel.eq.nodata) then
c       zel=10000.0
c      endif
c
c      zl1=zlcl
c      zl2=min(zel,10000.0)
c
c      do j=1,nz-1
c       if(z(j).le.zl1.and.z(j+1).ge.zl1) then
c          if(u(j).ne.nodata.and.u(j+1).ne.nodata.and.
c     $       v(j).ne.nodata.and.v(j+1).ne.nodata)     then
c           call linterpmike(z(j),u(j),z(j+1),u(j+1),zl1,ul1)
c           call linterpmike(z(j),v(j),z(j+1),v(j+1),zl1,vl1)
c          endif
c       elseif(z(j).le.zl2.and.z(j+1).ge.zl2) then
c          if(u(j).ne.nodata.and.u(j+1).ne.nodata.and.
c     $       v(j).ne.nodata.and.v(j+1).ne.nodata)     then
c           call linterpmike(z(j),u(j),z(j+1),u(j+1),zl2,ul2)
c          call linterpmike(z(j),v(j),z(j+1),v(j+1),zl2,vl2)
c           endif
c       endif
c      enddo
c
c      clyrshz=bulkshear(ul1,ul2,vl1,vl2)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c** CALCULATE THE 0-1 and 0-3 STORM RELATIVE ENVIRONMENTAL HELICITY USING THE
c** 20R85% TECHNIQUE 

      sreh01=nodata
      sreh03=nodata

      ! Find the mean wind speeds and directions

      call meanuvz(nodata,nz,0.,z,u,v,um,vm)
      call uv2dir( um(1,7),vm(1,7),dir0to6,spd0to6 ) 
      spd0to6 = 0.85*spd0to6
      dir0to6 = dir0to6 + 20.0
      if (dir0to6.gt.360.0) dir0to6 = dir0to6 - 360.0
      call dir2uv( dir0to6,spd0to6,270.0,ustorm,vstorm )
       
c      ustorm=10.0
c      vstorm=0.0

      do j=1,nz
       if(z(j).eq.1000.0) j1000=j
       if(z(j).eq.3000.0) j3000=j
      enddo

      if(u(j1000).eq.nodata) goto 93
 
      sum01=0.0
      do j=1,j1000-1
       if(u(j).ne.nodata.and.v(j).ne.nodata) then
        sum01=sum01+(
     $                ((u(j+1)-ustorm)*(v(j)-vstorm))
     $               -((u(j)-ustorm)*(v(j+1)-vstorm))
     $              )
       endif
      enddo
      sreh01=sum01

 93   continue
      if(u(j3000).eq.nodata) goto 94

      sum03=0.0
      do j=1,j3000-1
       if(u(j).ne.nodata.and.v(j).ne.nodata) then
        sum03=sum03+(
     $                ((u(j+1)-ustorm)*(v(j)-vstorm))
     $               -((u(j)-ustorm)*(v(j+1)-vstorm))
     $              )
       endif
      enddo
      sreh03=sum03

94    continue

      RETURN
      END SUBROUTINE shearz
c********************************************************************
      SUBROUTINE meanuvz(nodata,nz,zlow,zzi,uzi,vzi,uzim,vzim)

c** Creates a 13 x 13 array of mean wind speeds (m/s) from 0 to 12 km
c** every 1 km [ the 0-3 km mean wind speed is stored in uzim(1,4) or
c** in uzim(4,1) ] 

      implicit none

      integer nz,j,j12,k,l1,l2,m
      real nodata,zlow,zj,zk
      real, dimension(nz) :: zzi,uzi,vzi
      real, dimension(13,13) :: uzim,vzim

      do j=1,nz
       if(zzi(j).eq.12000.0) j12=j
      enddo

      uzim(:,:)=nodata
      vzim(:,:)=nodata

       do j=1,12
        zj=zlow+(float(j-1)*1000.)
        do l1=1,j12
         if(zzi(l1).eq.zj) then
          do k=j+1,13
           zk=zlow+(float(k-1)*1000.)
            do l2=1,j12
             if(zzi(l2).eq.zk) then
              if(uzi(l1).ne.nodata.and.uzi(l2).ne.nodata) then
                uzim(j,k)= SUM(uzi(l1:l2))/REAL((l2-l1)+1)
                vzim(j,k)= SUM(vzi(l1:l2))/REAL((l2-l1)+1)
                uzim(k,j)= uzim(j,k)
                vzim(k,j)= vzim(j,k)
              endif
             endif
            enddo
          enddo
         endif
        enddo
       enddo

       do j=1,13
        zj=zlow+(float(j-1)*1000.)
        do l1=1,j12
         if(zzi(l1).eq.zj) then
          uzim(j,j)=uzi(l1)
          vzim(j,j)=vzi(l1)
         endif
        enddo
       enddo

       RETURN
       END SUBROUTINE meanuvz
c***************************************************************************
      SUBROUTINE uv2dir(u,v,dirdeg,spd)

      real u,v,pi,rtd,quad,dirrad,dirdeg,spd

      pi=3.1415927
      rtd=180.0/pi
      small=0.1

      spd=SQRT( (u**2.0) + (v**2.0) )

      if(abs(u).le.small) then
        if(abs(v).le.small) then
          dirdeg=0.0
        elseif(v.gt.0.0) then
          dirdeg=180.0
        elseif(v.lt.0.0) then
          dirdeg=360.0
        endif
      else
        dirrad=ATAN(v/u)
        if(u.gt.0.0) then
         dirdeg=270.0-(rtd*dirrad)
        else
         dirdeg=90.0-(rtd*dirrad)
        endif
      endif

      RETURN
      END SUBROUTINE uv2dir
c*************************************************************************
      SUBROUTINE dir2uv(dir,spd,rot,u,v)

      real dir,spd,rot,utmp,vtmp,u,v,dtr
   ! rot defines the rotation angle of the new cartesian coordinate system
   ! rot=270.0 will produce the u,v coordinates relative to the typical 
   ! 270.0 = E-W orientation

      dtr=3.1415927/180.0

      utmp = spd*COS((270.-dir)*dtr)
      vtmp = spd*SIN((270.-dir)*dtr)
      u = vtmp*SIN(dtr*(270.-rot)) + utmp*COS(dtr*(270.-rot))
      v = vtmp*COS(dtr*(270.-rot)) - utmp*SIN(dtr*(270.-rot))

      RETURN
      END SUBROUTINE dir2uv
c*************************************************************************
      SUBROUTINE ZINTERP(ND,NLV,BADVAL,TDRY,TDEW,P,HGT
     $ ,NZ,H_OUT,VAR1_OUT,VAR2_OUT,P_OUT)

      INTEGER ND,NLV,NZ
      REAL BADVAL,ZSFC
      REAL TDRY(ND),TDEW(ND),P(ND),HGT(ND),H_OUT(NZ),VAR1_OUT(NZ),
     $ VAR2_OUT(NZ),P_OUT(NZ)

C**  FIRST CREATE AGL HEIGHTS FROM THE INPUT DATA

      ZSFC=HGT(1)
      DO I=1,NLV
       HGT(I)=HGT(I)-ZSFC
      ENDDO

C**  INTERPOLATE FROM SURFACE ON UP ACCORDING TO LEVELS SPECIFIED
C**  IN THE H_OUT ARRAY

      VAR1_OUT(:)=BADVAL
      VAR2_OUT(:)=BADVAL
      P_OUT(:)=BADVAL

      MLV=0
      DO I=1,NZ
        IF(H_OUT(I).GT.HGT(NLV)) GOTO 99
        DO 30 J=1,NLV
         IF(HGT(J).GE.H_OUT(I)) GOTO 31
 30     CONTINUE
        GOTO 99
 31     IF(H_OUT(I).EQ.HGT(J)) THEN
         VAR1_OUT(I)=TDRY(J)
         VAR2_OUT(I)=TDEW(J)
         P_OUT(I)=P(J)
         MLV=MLV+1
         GOTO 98
        ENDIF
        DO 32 K=J-1,1,-1
         IF(HGT(K).LE.H_OUT(I)) GOTO 33
 32     CONTINUE
        GOTO 99
 33     IF(H_OUT(I).EQ.HGT(K)) THEN
         VAR1_OUT(I)=TDRY(K)
         VAR2_OUT(I)=TDEW(K)
         P_OUT(I)=P(K)
        ELSE
         CALL LINTERPMIKE(HGT(K),TDRY(K),HGT(J),TDRY(J),
     $    H_OUT(I),VAR1_OUT(I))
C** THE DEWPOINTS MIGHT BE MISSING
         IF(TDEW(K).NE.BADVAL.AND.TDEW(J).NE.BADVAL) THEN
          CALL LINTERPMIKE(HGT(K),TDEW(K),HGT(J),TDEW(J),
     $     H_OUT(I),VAR2_OUT(I))
         ENDIF
         CALL LINTERPMIKE(HGT(K),ALOG(P(K)),HGT(J),ALOG(P(J)),
     $    H_OUT(I),P_OUT(I))
         P_OUT(I)=EXP(P_OUT(I))
        ENDIF
        MLV=MLV+1
 98     CONTINUE
      ENDDO

 99   CONTINUE

      IF(MLV.LE.3) THEN
       PRINT *,'WARNING: NUMBER OF INTERPOLATED LEVELS < 3'
       STOP
      ENDIF

      RETURN
      END SUBROUTINE ZINTERP
c********************************************************

      REAL FUNCTION bulkshear(u1,u2,v1,v2)

      real u1,u2,v1,v2,du,dv

      du=u2-u1
      dv=v2-v1
      bulkshear=SQRT( (du**2.0) + (dv**2.0) )

      RETURN
      END FUNCTION bulkshear
c********************************************************************
      REAL FUNCTION relh(t,d)

c t and d should be in Kelvin

      real t,d,vap,satvap

      vap=epres_td(d)
      satvap=esat(t)
      relh=100.*(vap/satvap)

      RETURN
      END FUNCTION relh
c********************************************************************
      real function esat(t)
c**  compute saturation vapor pressure (mb) using temperature (K)

      real t,tc

      tc=t-273.16
      esat=6.112*exp((17.67*tc)/(tc+243.5))
      return
      end
C********************************************************************
      REAL FUNCTION TMR(W,P)
C  TMR(KELVIN),W(GRAMS WATER VAPOR/KILOGRAM DRY AIR),P(MILLIBAR)
C  ALOG10  15   LOG TO THE BASE  TEN.
      X =  ALOG10(   W*P/(622.+ W)  )
      TMR=10.**(.0498646455*X+2.4082965)-7.07475+38.9114*
     X                   ((10.**( .0915*X ) - 1.2035 )**2 )
      RETURN
      END
C********************************************************************
      REAL FUNCTION TUA(OU,P,W)
C  FINDS THE TEMPERATURE (K) ALONG A DRY ADIABAT GIVEN THE THETA
C  ALONG THE ADIABAT (K), THE CURRENT PRESSURE (MB), AND THE MIXING
C  RATIO
      REAL R,CP,A1

      R=287.04
      CP=1004.0
      A1=(1 - (0.00028*W) )
      TUA = OU*( (P/1000.)**((R/CP)*A1) )
      RETURN
      END
c********************************************************************
      REAL FUNCTION TSA(OE,P,Q)
C  SEE BOLTON (1980) MWR FOR A REFERENCE
C  FINDS THE TEMPERATURE (K) ALONG A MOIST ADIABAT GIVEN THE THETA-E
C  ALONG THAT ADIABAT AND THE CURRENT PRESSURE
C  TSA AND OE(KELVIN),P(MILLIBARS)
C  SIGN(A,B) REPLACES THE ALGEBRAIC SIGN OF A WITH THE SIGN OF B

      REAL OE,P,A,TQ,D,THETA1,THETA2,X, Q

      A  =  OE
      TQ  =  253.16
      D =  120
      DO I = 1,20
         D = D/2.
         THETA1=A*EXP(-2.6518986*W(TQ,P)/TQ)
         THETA2=TQ*((1000./P)**.286)
C  IF THE TEMPERATURE DIFFERENCE,X, IS SMALL,EXIT THIS LOOP
         X=THETA1-THETA2
C         X=A*EXP(-2.6518986*W(TQ,P)/TQ)-TQ*((1000./P)**.286)
         IF(ABS(X).LT.0.01) GO TO 10
          TQ = TQ + SIGN(D,X)
      ENDDO
 10   TSA=TQ
      RETURN
      END
C*************************************************************************
      REAL FUNCTION  W(T,P)
C  T IN K
C  W(GRAMS WATER VAPOR/KILOGRAM DRY AIR ), P(MILLIBAR )
      IF(T.GE.999.)GO TO 10
         X  =       ESAT(T)
         W  =  621.97*X/(P-X)
      RETURN
   10 W=0.0
      RETURN
      END
c*************************************************************************
        real function epres(p,w)
c
c       compute water vapor pressure from p (mb) and w (g/kg)
c
        epres=(p*w)/(622+w)
        return
        end
c***********************************************************************
        real function theta(p,t,w)
c
c**     computes the dry air potential temperature from t(K),p(mb) and
c**     w(g/kg).  Theta output in K
c
        ak=0.2859*(1.-0.00028*w)
        theta=t*(1000.0/p)**ak
        return
        end
c************************************************************************
        real function theta_e(theta,tlcl,w)
c
c**     computes the equivalent potential temperature (K) from theta(K),
c**     tlcl(K) and w(g/kg).  Routine to calculate tlcl is below.
c**     See Bolton (1980) MWR 1046-1053
c
        argu=((3.376/tlcl)-0.00254)*w*(1.+0.00081*w)
        theta_e=theta*exp(argu)
        return
        end
c*************************************************************************
        real function tlcl1(t,epres)
c
c**     compute temperature of lcl(k) from t(k) and epres(mb)
c       See Bolton (1980) MWR 1046-1053
c
        tlcl1=(2840./(3.5*alog(t)-alog(epres)-4.805))+55.
        return
        end
c************************************************************************
        real function rl(t)
c
c**     latent heat of vaporization using t(C).
c
        tm2=t
        rl=597.31-0.589533*tm2+0.001005333*tm2**2
        return
        end
c*************************************************************************
        real function epres_td(td)
c
c       compute vapor pressure using dew point temperature (K)
c
        tdew=td-273.16
        epres_td=6.112*exp((17.67*tdew)/(tdew+243.5))
        return
        end
c*************************************************************************
        real function wmix(p,epres)
c
c**     compute mixing ratio wmix using p (mb) and epres (mb).  Output
c**     in g/kg.
c
        rd=287.04
        rv=461.50
        eps=rd/rv
c
        wmix=1000.*eps*epres/(p-epres)
        return
        end
c*************************************************************************
        real function theta_w(thae)
c
c**     compute equivalent potential wet-bulb temperature at 1000 mb
c**     using theta_e (k) and calculating w (g/kg)
c
        rd=287.04
        rv=461.50
        eps=rd/rv
c
        xlnthe=alog(thae)
c
        if(thae.le.283.) then
          tnewk=273.
        else
          tnewk=273.+3.*sqrt(thae-273.-10.)
        end if
        icount=0
c
c       begin newton method iteration
c
        toldk = 200.
        do while((abs(tnewk-toldk).ge.0.01).and.(icount.le.9))
        icount=icount+1
        toldk=tnewk
c
        told=toldk-273.15
        hotph=243.5/(told+243.5)
        e=6.112*exp((17.67*told)/(told+243.5))
c
c       w is in g/kg
c
        w=1000.*eps*e/(1000.-e)
        cotma=(3.376/toldk)-.00254
        onepbr=1.+0.00081*w
        drdt=w*(1000./(1000.-e))*(17.67/(told+243.5))*hotph
        g=alog(toldk)+cotma*onepbr*w
c
        dgdt=(1./toldk)*(1.-(3.376/toldk)*onepbr*w)+
     1    cotma*(1.+0.00162*w)*drdt
c
        tnewk=toldk+(xlnthe-g)/dgdt
        end do
        if(abs(tnewk-toldk).ge.2.0) then
         theta_w=-10.
         print *,' WARNING:Theta-w does not converge',tnewk
        else
         theta_w=tnewk-273.15
        end if
        return
        end
c**********************************************
      subroutine linterpmike(x1,y1,x2,y2,xpt,ypt)

      real x1,y1,x2,y2,xpt,ypt

      slope = (y2-y1)/(x2-x1)
      ypt = y1 + slope*(xpt-x1)
      return
      end
c**********************************************

