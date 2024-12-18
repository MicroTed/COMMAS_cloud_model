      subroutine settitles(ntit,title,units)
       
       implicit none
       integer ntit
      
       character*10 title(2,ntit)
       character*20 units(ntit)
       
       integer k
       
c  titles of plots
c
      title(1,1) = 'THETAP'
      units(1) = "K"
      title(1,2) = 'TEMP'
      units(2) = "C"
      title(1,3) = 'U'
      title(1,4) = 'V'
      title(1,5) = 'W'
      units(3) = "m/s"
      units(4) = "m/s"
      units(5) = "m/s"

      title(1,6) = 'FLASH' ! add only if there was a flash
      title(1,7) = 'FLASH' ! add regardless
      title(1,8) = 'DRHO' ! add only if there was a flash
      title(1,9) = 'DRHO' ! add regardless
      units(6) = "X"
      units(7) = "X"
      units(8) = "nC/m**3"
      units(9) = "nC/m**3"
      title(1,10) = 'CHAN'
      units(10) = 'X'
      title(1,11) = 'WZ'
      title(2,11) = 'VORT'
      units(11) = '1/s(*10^-3)'
      title(1,12) = 'THETA'
      units(12) = 'K'
      title(1,13) = 'THETAV'
      units(13) = 'K'
      title(1,14) = 'TBASE'
      units(14) = 'K'
      title(1,15) = 'TDBASE'
      units(15) = 'K'
      title(1,16) = 'T'
      units(16) = 'K'
      title(1,17) = 'TD'
      units(17) = 'K'
      title(1,18) = 'FLSHI' ! lightning initiation points
      title(1,19) = 'FLSHP' ! lightning summary positive channel count
      title(1,20) = 'FLSHN' ! lightning summary negative channel count
      title(1,21) = 'DRHOM' ! 
      units(18) = "X"
      units(19) = "X"
      units(20) = "X"
      units(21) = "nC/m**3"

      title(1,22) = 'DBZ'
      title(2,22) = 'REF'
      units(22) = 'dBZ'
      
      title(1,23) = 'KM'
      title(2,23) = 'TKE'
      units(23) = 'm**2/s'
      
      title(1,24) = 'SCRN'
      units(24) = "nC/m**3"

      title(1,25) = 'CPION'
      title(2,25) = 'CIONP'
      units(25) = "cm**-3"

      title(1,26) = 'CNION'
      title(2,26) = 'CIONN'
      units(26) = "cm**-3"

      title(1,27) = 'SCION'
      units(27) = "nC/m**3"
      
      title(1,28) = 'X'
      units(28) = 'km'
      
      title(1,29) = 'Y'
      units(29) = 'km'

      title(1,30) = 'GLM_8KM' ! Fake GLM data
      units(30) = "FED"

      title(1,31) = 'FED' ! native grid FED over history interval 'FLSHICGP' ! +CG lightning initiation points
      units(31) = "FED"
      title(1,32) = 'FLSHICGN' ! -CG lightning initiation points
      units(32) = "X"

      title(1,33) = 'FLSH' ! lightning density, positive and negative
      units(33) = 'X'   
      
      title(1,34) = 'VZF' ! z-weighted fall speed
      units(34) = 'm/s'

      title(1,35) = 'UVMAG' ! Horizontal wind speed
      units(35) = 'm/s'

      title(1,45) = 'SSIA'  ! ice supersaturation after microphysics
      units(45) = '%'

      title(1,46) = 'SSWA'  ! water supersaturation after microphysics
      units(46) = '%'

      title(1,50) = 'DBZCMP'  ! composite reflectivity
      units(50) = 'dBZ'

      title(1,51) = 'DBZSFC'  ! surface reflectivity
      units(51) = 'dBZ'
      
      title(1,52) = 'dBZ'  !  Time-Height output
      units(52) = 'dBZ'

      title(1,53) = 'charge'  !  time-height output
      units(53) = 'C/level'
      
      title(1,54) = 'charge'  ! Vernon Cooray charge dump
      units(54) = 'C/m**3'
      
      
      title(1,55) = 'PP'
      title(2,55) = 'PPRES'
      units(55) = 'Pa'
      
      title(1,56) = 'P'
      title(2,56) = 'PRES'
      units(56) = 'mb'

      title(1,57) = 'PI'
      title(2,57) = 'EXNERPI'
      units(57) = 'mb/mb'

      title(1,58) = 'HAILRAT'
      units(58) = 'kg/m^2/s'

      title(1,59) = 'HAILACC'
      units(59) = 'kg/m^2'


      title(1,60) = 'IWP'
      units(60) = 'g/m^2'

      title(1,61) = 'FLSHR'
      units(61) = 'flsh/clm'

      title(1,62) = 'RAINRAT'
      units(62) = 'kg/m^2/s'

      title(1,63) = 'RAINACC'
      units(63) = 'kg/m^2'

      title(1,65) = 'TCANP' ! Canopy potential temperature
      units(65) = 'K'

      title(1,66) = 'WCANP' ! Canopy mixing ratio
      units(66) = 'kg kg-1'
      
      title(1,67) = 'QAV' ! Mixing ratio within foliage
      units(67) = 'kg kg-1' 
      
      title(1,68) = 'VEG' ! Vegetatiohn density (fraction)
      units(68) = 'fraction'        
      
      title(1,69) = 'STYPE' ! Soil type
      units(69) = 'number'        
            
      title(1,70) = 'TSRFC' ! Surface potential temperature
      units(70) = 'K'        
      
      title(1,71) = 'WSFC' ! Surface mixing ratio
      units(71) = 'kg kg-1'
      
      title(1,72) = 'TSOIL' ! Soil potential temperature
      units(72) = 'K'     
    
      title(1,73) = 'WSOIL' ! Soil mixing ratio
      units(73) = 'kg kg-1'        
                   
      title(1,74) = 'VLAI' ! Leaf area index
      units(74) = 'number'
      
      title(1,75) = 'ALBEDO' ! Albedo
      units(75) = 'fraction'
      
      title(1,76) = 'ROUGH' ! Roughness length
      units(76) = 'm'        
                    
      title(1,77) = 'KH' ! Heat eddy mixing coefficient
      units(77) = 'm2 s-1'        

      title(1,78) = 'EFLX' ! Vegetation moisture flux
      units(78) = 'm-2 s-1'

      title(1,79) = 'FFLX' ! Evaporation moisture flux
      units(79) = 'm-2 s-1'
                  
      title(1,80) = 'UFLX' ! U momentum flux
      units(80) = 'm-1 s-2'
                        
      title(1,81) = 'VFLX' ! V momentum flux         
      units(81) = 'm-1 s-2'   

      title(1,82) = 'TFLX' ! SFC sensible heat flux         
      units(82) = 'w m-2'

      title(1,83) = 'QFLX' ! SFC moisture flux         
      units(83) = 'w m-2'        

      title(1,84) = 'RADSW' ! Net SW radiation         
      units(84) = 'w m-2'        
                                    
      title(1,85) = 'RADLW' ! Net LW radiation         
      units(85) = 'w m-2'        

      title(1,86) = 'SFCPP'
      units(86) = 'Pa'
      
      title(1,87) = 'SFCP'
      units(87) = 'mb'
             
!      title(1,86) = 'PERTU' ! EnKF U perturbation         
!      units(86) = 'm s-1'

!      title(1,87) = 'PERTV' ! EnKF V perturbation         
!      units(87) = 'm s-1'      

!      title(1,88) = 'PERTW' ! EnKF W perturbation         
!      units(88) = 'm s-1'      

!      title(1,89) = 'PERTT' ! EnKF T perturbation         
!      units(89) = 'K'      

!      title(1,90) = 'PERTQ' ! EnKF Q perturbation         
!      units(90) = 'm s-1'            

      title(1,91) = 'THETAE'
      units(91) = 'K'

      title(1,92) = 'DIVERG'
      title(2,92) = 'CONVERG'
      units(92) = '1/s(*10^-3)'

      title(1,103) = 'US' ! staggered values
      title(1,104) = 'VS'
      title(1,105) = 'WS'
      units(103) = "m/s"
      units(104) = "m/s"
      units(105) = "m/s"

      title(1,118) = 'FLSHDEN'
      units(118) = 'log_src'
      title(1,119) = 'FLSHPT' ! lightning summary positive channel count, delete above cloud top
      title(1,120) = 'FLSHNT' ! lightning summary negative channel count
      units(119) = "X"
      units(120) = "X"
      title(1,121) = 'FLSHNDN' ! 
      units(121) = 'log_src'
      title(1,122) = 'FLSHPDN' ! 
      units(122) = 'log_src'
      title(1,123) = 'FLSHLMA' ! negative plus 10% of positive
      units(123) = 'log_src'

      title(1,133) = 'FLSHT' ! lightning density, positive and negative
      units(133) = 'X'   

      title(1,241) = 'CFW' ! 'frozen drops'
      title(1,242) = 'CHW' ! 'SMALL HAIL: 0.5<D<2.0 (CM) (N/M**3)'
      title(1,243) = 'CHL' !'LARGE HAIL: 2.0>D (CM) (N/M**3)'
      title(1,244) = 'CHLF' ! number of hail from frozen drops
      title(1,245) = 'CHF' ! number of graupel from frozen drops
      DO k=241,245
        units(k) = "#/m**3"
      END DO

      title(1,250) = 'CPLION'
      title(2,250) = 'CLIONP'
      units(250) = "cm**-3"

      title(1,251) = 'SCPION'
      title(2,251) = 'SCIONP'
      units(251) = "nC/m**3"

      title(1,252) = 'SCPLION'
      title(2,252) = 'SCLIONP'
      units(252) = "nC/m**3"

      title(1,260) = 'CNLION'
      title(2,260) = 'CLIONN'
      units(260) = "cm**-3"

      title(1,261) = 'SCNION'
      title(2,261) = 'SCIONN'
      units(261) = "nC/m**3"

      title(1,262) = 'SCNLION'
      title(2,262) = 'SCLIONN'
      units(262) = "nC/m**3"

      title(1,270) = 'SCLION'
      units(270) = "nC/m**3"
c
      title(1,392) = 'E2ENRGY' ! 'Electrical energy (J) from E^2 only'
      units(392) = "MJ"
      title(1,393) = 'EENERGY' ! 'Electrical energy (J)'
      units(393) = "MJ"
      title(1,394) = 'SCITOT' ! 'TOTAL CLOUD ICE CHARGE (nC/m**3)'
      title(1,395) = 'SCGHTOT' ! 'TOTAL GRAUPEL + HAIL CHARGE (nC/m**3)'
      title(1,396) = 'SCHL' ! 'LARGE HAIL CHARGE (nC/m**3)'
      title(1,397) = 'SCH_TOT' ! total hail+frozen drop charge
      title(1,399) = 'SCCW' ! 'CLOUD WATER CHARGE (nC/m**3)'
      title(1,400) = 'SCCI' ! 'CLOUD ICE-COLUMN CHARGE (nC/m**3)'
      title(1,401) = 'SCIP' ! 'CLOUD ICE-PLATE CHARGE (nC/m**3)'
      title(1,402) = 'SCIR' ! 'CLOUD ICE-RIMED CHARGE (nC/m**3)'
      title(1,403) = 'SCRW' ! 'RAIN WATER CHARGE (nC/m**3)'
      title(1,404) = 'SCS' ! 'SNOW AGGREGATE CHARGE (nC/m**3)'
      title(1,405)  = 'SCGL' ! 'LOW DEN GRAUPEL CHARGE (nC/m**3)' 
      title(1,406)  = 'SCGM' ! 'MED DEN GRAUPEL CHARGE (nC/m**3)'
      title(1,407)  = 'SCGH' ! 'HI DEN GRAUPEL CHARGE (nC/m**3)'
      title(1,408)  = 'SCFD' ! 'FROZEN DROP CHARGE (nC/m**3)'
      title(1,409) = 'SCH' ! 'HAIL CHARGE (nC/m**3)'
      title(1,410) = 'SCNET' ! 'TOTAL CHARGE (nC/m**3) from rho arrary'
      title(2,410) = 'SCTOT' ! 'TOTAL CHARGE (nC/m**3) from rho arrary'
      title(1,411) = 'EX' ! 'ELECTRIC FIELD: X-DIR (kV/m)'
      title(1,412) = 'EY' ! 'ELECTRIC FIELD: Y-DIR (kV/m)'
      title(1,413) = 'EZ' ! 'ELECTRIC FIELD: Z-DIR (kV/m)'
      title(1,414) = 'EMAG' ! 'ELECTRIC FIELD: |E| (kV/m)'
      title(1,415) = 'POT' ! 'ELECTRIC POTENTIAL FIELD'
      title(1,416) = 'SCG_NET' ! 'TOTAL GRAUPEL CHARGE (nC/m**3)'
      title(1,417) = 'SCICE' ! 'TOTAL ICE CHARGE (nC/m**3)'
      title(1,418) = 'SCNET2' ! 'TOTAL CHARGE (nC/m**3) sum up categories'
      title(1,419) = 'EZW' ! 'E-Field: Ez at w-pnts (kV/m)'
      
c      title(1,420) = 'SCLIQ' ! sccw + scrw
      title(1,41) = 'EZSFC' ! 'E-Field: Ez at surface (kV/m)'
      DO k=394,410
        units(k) = "nC/m**3"
      END DO
      DO k=411,414
        units(k) = "kV/m"
      END DO
      units(415) = "MV"
      DO k=416,418
        units(k) = "nC/m**3"
      END DO
      units(419) = "kV/m"
      units(41) = "kV/m"

c
      title(1,420) = 'CCW' ! 'CLOUD WATER CONCENTRATION (#/m**3)'
      title(1,421) = 'CCI' !  'CLOUD ICE-COLUMN CONCENTRATION (#/m**3)'
      title(1,422) = 'CIP' ! 'CLOUD ICE-PLATE CONCENTRATION (#/m**3)'
      title(1,423) = 'CIR' ! 'CLOUD ICE-RIMED CONCENTRATION (#/m**3)'
      title(1,424) = 'CRW' ! 'RAIN WATER CONCENTRATION (#/m**3)'
      title(1,425) = 'CSW'  ! 'SNOW AGGREGATE CONCENTRATION (#/m**3)'
      title(1,426) = 'CGL' ! ''LOW DEN GRAUPEL CONCENTRATION (#/m**3)' 
      title(1,427) = 'CGM' ! ''MED DEN GRAUPEL CONCENTRATION (#/m**3)'
      title(1,428) = 'CGH' ! 'HI DEN GRAUPEL CONCENTRATION (#/m**3)'
      title(1,429) = 'CFW' ! 'FROZEN DROP CONCENTRATION (#/m**3)'
      title(1,430) = 'CHTOT' ! 'TOTAL HAIL CONCENTRATION (#/m**3)'
      title(1,431) = 'CGTOT' ! 'TOTAL GRAUPEL CONCENTRATION (#/m**3)'
      title(1,432) = 'CITOT' ! 'TOTAL ICE CONCENTRATION (#/m**3)'
      
      units(420) = "#/cm**3"
      DO k=421,432
        units(k) = "#/m**3"
      END DO
      
      units(429) = "#/L" ! used now for small frozen drops from the TAK scheme
      
      units(421) = "#/L"
      units(422) = "#/L"
      units(424) = "#/L" ! rain
      units(425) = "#/L" ! snow agg.
      units(432) = "#/L"
      
      title(1,433) = 'CWRAD' ! 'CLOUD WATER Diameter'
      units(433) = "micron"
      title(1,434) = 'SS' ! supersaturation w.r.t. liquid water before microphysics
      units(434) = "%"
      title(1,435) = 'SSI' ! supersaturation w.r.t. ice
      units(435) = "%"
      title(1,436) = 'RWRAD' ! 'Rain drop radius'
      units(436) = "micron"
      title(1,437) = 'SSMX' ! supersaturation
      units(437) = "%"

      title(1,438) = 'CCN' ! number of unactivated CCN
      title(2,438) = 'CCCN' ! number of unactivated CCN
      units(438) = '#/cm**3'

      title(1,439) = 'SSNEW' ! supersaturation w.r.t. liquid _after_ microphys
      units(439) = "%"

      title(1,440) = 'CCHAFF' ! number of chaff fibers
      units(440) = '#/m**3'

      title(1,441) = 'CLNOX' ! lightning NOx concentration
      units(441) = 'ppbv'

      title(1,442) = 'CCNA' ! number of activated CCN
      title(2,442) = 'CCCNA' ! number of activated CCN
      units(442) = '#/cm**3'

      title(1,443) = 'CINA' ! number of activated IN
      title(2,443) = 'CCINA' ! number of activated IN
      units(443) = '#/L'

c
      title(1,464)  = 'CIDIA'
      title(1,465)  = 'SWDIA'
      title(1,466)  = 'GLDIA' 
      title(1,467)  = 'GMDIA'
      title(1,468)  = 'GHDIA'
      title(1,469)  = 'FWDIA'
      title(1,470)  = 'HWDIA'
      title(1,471)  = 'HWVDIA'
      title(1,472)  = 'HLDIA'
      DO k=464,472
        units(k) = "micron"
      END DO

      title(1,474)  = 'HWMAS' ! mean particle mass
      title(1,475)  = 'HLMAS' ! mean particle mass
      title(1,476)  = 'SWMAS' ! mean particle mass
      title(1,477)  = 'RWMAS' ! mean particle mass
      DO k=474,477
        units(k) = "mg"
      END DO
      
      title(1,480)  = 'N0SW'
      title(1,481)  = 'N0GL' 
      title(1,482)  = 'N0GM'
      title(1,483)  = 'N0GH'
      title(1,484)  = 'N0FW'
      title(1,485)  = 'N0HW'
      title(1,486)  = 'N0HL'

      DO k=480,486
        units(k) = "m**-4"
      END DO

      title(1,490)  = 'HMDIA' ! maximum mass diameter
      title(1,491)  = 'HLMDIA'
      title(1,492)  = 'HMWD' !  mass-weighted diameter
      title(1,493)  = 'HLMWD'
      title(1,494)  = 'RWMWD' !  mass-weighted rain diameter
      title(1,495)  = 'RWMDIA' ! maximum mass rain diameter
      title(1,496)  = 'RWD0' ! rain median volume diameter
      title(1,497)  = 'RWDIA' ! mean mass rain diameter
      title(1,498)  = 'HCDIA' ! characteristic diameter
      title(1,499)  = 'HLCDIA' ! characteristic diameter
      DO k=490,499
        units(k) = "micron"
      END DO
      

      title(1,502)  = 'QVP' ! perturbation qv
      title(1,503)  = 'QV'
      title(1,504)  = 'QC'
      title(1,505)  = 'QR'
      title(1,506)  = 'QCLMI' ! column ice
      title(1,507)  = 'QS'
      title(1,508)  = 'QGL'
      title(1,509)  = 'QGM'
      title(1,510)  = 'QGH'
      title(1,511)  = 'QFD' 
      title(1,512)  = 'QH' 
      title(1,513)  = 'QHL' 
      title(1,514)  = 'QIR' 
      title(1,515)  = 'QIP' 
      title(1,516)  = 'QSW' ! water fraction mixing ratio on snow
      title(1,517)  = 'QHW' ! water fraction mixing ratio on graupel
      title(1,518)  = 'QHI' ! ice fractional mixing ratio (graupel)
      title(1,519)  = 'QHLW' ! water fraction mixing ratio on hail
      DO k=502,519
        units(k) = "g/kg"
      END DO
      
      title(1,520)  = 'QHLLOW'
        units(520) = "g/kg"

      title(1,521)  = 'FWS' ! water fraction on snow
      title(1,522)  = 'FWH' ! water fraction on graupel
      title(1,523)  = 'FWHL' ! water fraction on hail

      DO k=521,523
        units(k) = "%"
      END DO


      title(1,524)  = 'QRLOW'
      title(1,525)  = 'QHLOW'
      DO k=524,525
        units(k) = "g/kg"
      END DO
      title(1,526)  = 'CRWLOW'
        units(526) = "#/m**3"

      title(1,527)  = 'CCILOW'
        units(527) = "#/L"

      title(1,528)  = 'CSWLOW'
        units(528) = "#/m**3"

c
c
      title(1,529) = 'QITOT' ! 'TOTAL CLOUD ICE MIXING RATIO (G/KG)'
      title(1,530) = 'QGTOT' ! 'TOTAL GRAUPEL MIXING RATIO (G/KG)'
      title(1,531) = 'QCI' ! 'TOTAL CLOUD MIXING RATIO (G/KG)'
      title(1,532) = 'QCIS' ! 'TOTAL CLOUD + SNOW MIXING RATIO (G/KG)'
      title(1,533) = 'QLWTOT' ! 'TOTAL LIQUID WATER MIXING RATIO (G/KG)'
      title(1,534) = 'QISTOT' ! 'TOTAL CLOUD ICE + SNOW MIXING RATIO (G/KG)'
      title(1,535) = 'QIPRE' ! 'TOTAL PRECIP ICE MIXING RATIO (G/KG)'
      title(1,536) = 'QICE' ! 'TOTAL ICE WATER MIXING RATIO (G/KG)'
      title(1,537) = 'QTOT' ! 'TOTAL CONDENSATE MIXING RATIO (G/KG)'
      title(1,538) = 'QPRE' ! 'TOTAL PRECIPITATION MIXING RATIO (G/KG)'
      title(1,539) = 'QGHTOT' ! 'TOTAL GRAUPEL + HAIL MIXING RATIO (G/KG)'
      title(1,540) = 'QHTOT' ! TOTAL HAIL MIXING RATIO (G/KG)'

      DO k=529,540
        units(k) = "g/kg"
      END DO

c
c integrated particle volumes:
c
      title(1,556)  = 'VCI' 
      title(1,557)  = 'VSW' 
      title(1,558)  = 'VGL' 
      title(1,559)  = 'VGM' 
      title(1,560)  = 'VGH' 
      title(1,561)  = 'VF' 
      title(1,562)  = 'VHW' 
      title(1,563)  = 'VHL' 

      DO k=550,563
        units(k) = "cm**3/m**3"
      END DO

      title(1,564) = 'ALPH'
      title(1,565) = 'ALPHL'
      title(1,566) = 'ALPR'
!      DO k=564,566
!        units(k) = "unit"
!      END DO
c
c reflectivities (6th moments):
c
      title(1,570)  = 'ZRW' 
      title(1,571)  = 'ZSW' 
      title(1,572)  = 'ZGL' 
      title(1,573)  = 'ZGM' 
      title(1,574)  = 'ZGH' 
      title(1,575)  = 'ZF' 
      title(1,576)  = 'ZHW' 
      title(1,577)  = 'ZHL' 
      title(1,578)  = 'ZRW2' 
      title(1,579)  = 'ZCI' 

      DO k=570,579
        units(k) = "dBZ"
      END DO

c
c integrated particle density:
c
      title(1,586)  = 'DNI' 
      title(1,587)  = 'DNS' 
      title(1,588)  = 'DNGL' 
      title(1,589)  = 'DNGM' 
      title(1,590)  = 'DNGH' 
      title(1,591)  = 'DNF' 
      title(1,592)  = 'DNH' 
      title(1,593)  = 'DNHL' 

      DO k=580,593
        units(k) = "kg/m**3"
      END DO

c
c Mass contents
c
      title(1,604)  = 'CWC'
      title(1,605)  = 'RC'
      title(1,606)  = 'CLMIC' ! column ice
      title(1,607)  = 'SC'
      title(1,608)  = 'GLC'
      title(1,609)  = 'GMC'
      title(1,610)  = 'GHC'
      title(1,611)  = 'FDC' 
      title(1,612)  = 'HWC' 
      title(1,613)  = 'HLC' 
      title(1,614)  = 'IRC' 
      title(1,615)  = 'IPC' 
      DO k=602,615
        units(k) = "g/m**3"
      END DO
c
c
      title(1,629) = 'ITOTC' ! 'TOTAL CLOUD ICE CONTENT (G/M**3)'
      title(1,630) = 'GTOTC' ! 'TOTAL GRAUPEL CONTENT (G/M**3)'
      title(1,631) = 'CIC' ! 'TOTAL CLOUD CONTENT (G/M**3)'
      title(1,632) = 'CISC' ! 'TOTAL CLOUD + SNOW CONTENT (G/M**3)'
      title(1,633) = 'LWTOTC' ! 'TOTAL LIQUID WATER CONTENT (G/M**3)'
      title(1,634) = 'ISTOTC' ! 'TOTAL CLOUD ICE + SNOW CONTENT (G/M**3)'
      title(1,635) = 'IPREC' ! 'TOTAL PRECIP ICE CONTENT (G/M**3)'
      title(1,636) = 'ICEC' ! 'TOTAL ICE WATER CONTENT (G/M**3)'
      title(1,637) = 'TOTC' ! 'TOTAL CONDENSATE CONTENT (G/M**3)'
      title(1,638) = 'PREC' ! 'TOTAL PRECIPITATION CONTENT (G/M**3)'
      title(1,639) = 'GHTOTC' ! 'TOTAL GRAUPEL + HAIL CONTENT (G/M**3)'
      title(1,640) = 'HTOTC' ! TOTAL HAIL CONTENT (G/M**3)'

      DO k=629,640
        units(k) = "g/m**3"
      END DO

c miscellaneous rate
      title(1,790) = 'CRFRZ'
      title(1,791) = 'CHCNSI'

      title(1,800) = 'QPROC'
      units(800) = 's^-1'

c noninductive charging rates: snow-ice
!      title(1,800) = 'SCSACI'
      title(1,801) = 'SCSACIP'
      title(1,802) = 'SCSACIR'
      
c noninductive charging rates: low den graup-snow,ice
      title(1,803) = 'SCGLACS'
      title(1,804) = 'SCGLACI'
      title(1,805) = 'SCGLACIR'
      title(1,806) = 'SCGLACIP'
      
c noninductive charging rates: med den graup-snow,ice
      title(1,807) = 'SCGMACS'
      title(1,808) = 'SCGMACI'
      title(1,809) = 'SCGMACIR'
      title(1,810) = 'SCGMACIP'
      
c noninductive charging rates: high den graup-snow,ice
      title(1,811) = 'SCGHACS'
      title(1,812) = 'SCGHACI'
      title(1,813) = 'SCGHACIR'
      title(1,814) = 'SCGHACIP'
      
c noninductive charging rates: hail-snow,ice
      title(1,815) = 'SCHACS'
      title(1,816) = 'SCHACI'
      title(1,817) = 'SCHACIR'
      title(1,818) = 'SCHACIP'
      
c noninductive charging rates: large hail-snow,ice
      title(1,819) = 'SCHLACS'
      title(1,820) = 'SCHLACI'
      title(1,821) = 'SCHLACIR'
      title(1,822) = 'SCHLACIP'

c noninductive charging rates: frozen drops-snow,ice
      title(1,823) = 'SCFACS'
      title(1,824) = 'SCFACI'
      title(1,825) = 'SCFACIR'
      title(1,826) = 'SCFACIP'

c grouped noninductive charging rates: 
      title(1,827) = 'CGI_NET'
      title(1,828) = 'CGS_NET'
      title(1,829) = 'CHI_NET'
      title(1,830) = 'CHS_NET'

c inductive charging rates:
      title(1,831) = 'SCSACW'
      title(1,832) = 'SCGLACW'
      title(1,833) = 'SCGMACW'
      title(1,834) = 'SCGHACW'
      title(1,835) = 'SCHACW'
      title(1,836) = 'SCHLACW'
      title(1,837) = 'SCFACW'

c grouped inductive charging rates:
      title(1,838) = 'CGW_NET'
      title(1,839) = 'CHW_NET'

c rime accretion rates:
      title(1,840) = 'RAR_S'
      title(1,841) = 'RAR_GL'
      title(1,842) = 'RAR_GM'
      title(1,843) = 'RAR_GH'
      title(1,844) = 'RAR_H'
      title(1,845) = 'RAR_HL'
      title(1,846) = 'RAR_FD'

c more group noninductive rates 
      title(1,847) = 'SCGL_I'
      title(1,848) = 'SCGM_I'
      title(1,849) = 'SCGH_I'
      title(1,850) = 'SCH_I'
      title(1,851) = 'SCHL_I'
      title(1,852) = 'SCFD_I'
      title(1,853) = 'SCS_I'
      
      title(1,854) = 'CGHI'
      title(1,855) = 'CGHS'
      title(1,856) = 'CGHIS'
      
      title(1,857) = 'CGHW'
      title(1,858) = 'CGHNET'

      title(1,871) = 'CADDL'
      title(1,872) = 'CNOLIQ'
      units(871) = 'pC/(m**3 s)'
      units(872) = 'pC/(m**3 s)'
      
      title(1,860) = 'QHSHR'
      units(860) = 's^-1'
      

      DO k=801,839
        units(k) = "pC/(m**3 s)"
      END DO
      DO k=840,846
        units(k) = "g/(m**2 s)"
      END DO
      DO k=847,858
        units(k) = "pC/(m**3 s)"
      END DO

      RETURN
      END
