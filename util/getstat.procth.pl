#!/usr/bin/perl
#
# 11/2012: get time-height microphysics process data from .stat file (1-minute integrated from thproc array in micro_driver/icezvd_gs)
#
# a perl program to read in log file(s) and track the max updraft
#   and its location at some interval  (set to 1 minute bins)
# Looks for the lines:
# $infile = 
#open(INPUT, "$run.hlgttyphsll") or die "file not found";
#open(INPUT, "linegrab") or die "file not found";

#open(INPUT, "$run.log") or die "file not found";

#open(OUT, ">$run.flashbin") or die "cannot create file";
 open(OUT, '>-');
$cr = chr(13);  # newline character for Mac
$lf = chr(10);  # newline character for Unix
#$/ = $lf;  # set record separator to Unix newline (default is chr(13) )
# foreach $file (@ARGV) {
#   open INPUT, "<$file";
#   $next = <INPUT>;
#   chomp($next);
#    @parts = split( / +/, $next);
#    print OUT ("$next \n");
#    }
    $binshift = 0;
    $tries = 0;
    $bintry = 0;
    $curbin = 0;
    $ic = 0;
    $cgn = 0;
    $cgp = 0;
    $ictot = 0;
    $cgntot = 0;
    $cgptot = 0;
    $icdistot = 0;
    $icchargetot = 0;
    $cgchargetot = 0;
    $cgpchargetot = 0;
    $cgpchargetot = 0;
    $cgpcharge = 0;
    $cgncharge = 0;
    $triestot = 0;
    $addtime = 0;
    $addstep = 0;
    $nonindtotp = 0;
    $nonindtotn = 0;
    $indtot = 0;
    $nonindp = 0;
    $nonindn = 0;
    $ind = 0;
    $ghsn = 0;
    $ghsp = 0;
    $ghin = 0;
    $ghip = 0;
    $istrz = 0;
    $iwcmx = 0;
    $nrain = 0;
    $layercheck = 0;
# foreach ...
foreach $file (@ARGV) {
    open INPUT, "<$file";
    $start = 0;
#    if ( $file =~ /wk30e4/ ) {
#        $addtime = 900.0;
#        $addstep = 180;
#      } else {
#        $addtime = 0;
#        $addstep = 0;
#      }
    while ( <INPUT> ) {
    
#nx,ny,nz,na =    81   81   81   42
    if ( /^nx,ny,nz,ns =/ ) {
      chomp();
      @parts = split(/ +/);
      $nz = $parts[4];
      
    }
    if ( /^ Layer average, max & min for scalars/ && $layercheck == 0  ) {
       $layercheck = 1; # use to know that we are past the first set of 'amax' prints
       }
    if ( /^ amax/  && $layercheck == 0 ) {
        # only go through the first set of amax prints to check for rain types
        chomp();
        @parts = split(/ +/);
        $len = $#parts;
        $string = $parts[8];
        if ( $string eq 'QRAUTO' ) { $nrain = $nrain + 1 }
        if ( $string eq 'QRSHED' ) { $nrain = $nrain + 1 }
        if ( $string eq 'QRMELT' ) { $nrain = $nrain + 1 }
       # print OUT ( "$string, $len, $nrain \n");
    
    }
      
# istrx,y,z=           0           0           1
#    if ( /^ istrx/ ) {
#      chomp();
#      @parts = split(/ +/);
#      $istrz = $parts[4];
#      }
# STRETCH FAC 
    if ( /^ KZ, ZC/ ) {
#      $next = <INPUT>;
#      chomp($next);
# kz, scal-hgt, w-hgt, mfc, mfe, 0.5*(mfc[k-1]+mfc[k]) 0.5*(mfe[k+1]+mfe[k])
#   1     80.      5.  1.4000  1.4000  1.4000  1.3969
       
     for ($kz = 1; $kz < $nz; $kz++ ) {      
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
       $z[$kz] = $parts[2]*0.001;
  #     print OUT ( "kz,z $kz $z[$kz] \n" );
      }

       }
    
      
    
    if ( /nstop=/ && $start == 0 ) {
      chomp();
      @parts = split( / +/ );
      $nstop = $parts[$#parts];
#      print("nstop=$nstop\n");
      }
    if ( /^nstep=/ && $start == 0 ) {
      $start = 1;
      @parts = split( / +/ );
      $nstep = $parts[1];
      }

#     if ( /^Max dBZ and integrated UDV at time =/ ) {
#       chomp();
#      @parts = split(/ +/);
#      $time[$curbin] = $parts[8];
#      }
     
     if ( /^Layer microphysics processes at time/ ) {
     chomp();
     @parts = split(/ +/);
     $min[1] = $parts[6]/60.0;
#     print OUT ( "minute = $min[1], $parts[6],$parts[5] \n" );
     # skip 2 lines
      $next = <INPUT>;
      $next = <INPUT>;
     
     for ($kz = $nz-1; $kz > 0; $kz-- ) { 
      $next = <INPUT>;
      $_ = $next;
      chomp($next);

      @parts = split(/ +/, $next);
      

#kz,ctghsnz,ctghspz,ctghinz,ctghipz,ctghwnz,ctghwpz
      # $i is the number of extra values (kz, altitude, temperature) before the thproc values start
      if ( /^ / ) { $i = 3 } else { $i = 2 }

       $len = $#parts - $i;
       # if ( $kz == 1 ) {
       #  print OUT ("len = $len\n");
       #  }
#      $crfrzf[$kz] = $parts[1+$i];
#      $ciacrf[$kz] = $parts[2+$i];
#      $chcnsh[$kz] = $parts[3+$i];
#      $chcnih[$kz] = $parts[4+$i];
#      $qhacw[$kz] = $parts[5+$i];
#      $qracw[$kz] = $parts[6+$i];
#      $qrcnw[$kz] = $parts[7+$i];
      
       for ($j = 1; $j < $len+1; $j++ ) {
         $procth[$kz][$j] = $parts[$j+$i];
         } # for j
         
         $qrfrzf[$kz] = $procth[$kz][12];
         $qiacrf[$kz] = $procth[$kz][13];
         $vhacw[$kz] = $procth[$kz][8];
         $qhacw[$kz] = $procth[$kz][5];
         $crfrzf[$kz] = $procth[$kz][1];
         $ciacrf[$kz] = $procth[$kz][2];
     if ( $vhacw[$kz] > 0 ) {
        $rimdens[$kz] = $qhacw[$kz]/$vhacw[$kz];
      } else {
         $rimdens[$kz] = 0;
      }
         $chcnsh[$kz] = $procth[$kz][3];
         $chcnih[$kz] = $procth[$kz][4];
         $chmul1[$kz] = $procth[$kz][10];
         $csplinter[$kz] = $procth[$kz][11];
         $cwctfz[$kz] = $procth[$kz][33];
         $ciint[$kz] = $procth[$kz][35];
       
       if ( $len - $nrain >= 44 ) {
           $qhshr[$kz] = $procth[$kz][42];
           $qfshr[$kz] = $procth[$kz][43];
           $qhlshr[$kz] = $procth[$kz][44];
           $chlcnh[$kz] = $procth[$kz][36];
           $chlcnf[$kz] = $procth[$kz][37];
          } else {
           $qhshr[$kz] = 0;
           $qfshr[$kz] = 0;
           $qhlshr[$kz] = 0;
           $chlcnh[$kz] = 0;
           $chlcnf[$kz] = 0;
          }
       
       if ( $nrain > 0 ) {
          $len = $#parts;
          for ($j = 1; $j <= $nrain; $j++ ) {
          $m = $len-$nrain+$j;
          $raintypefrz[$kz][$j] = $parts[$m];
          }
         } else {
          for ($j = 1; $j <= 5; $j++ ) {
            $raintypefrz[$kz][$j] = 0;
            }
         }
         
      } # for

      if ( $curbin == 0 ) {
#      print OUT ( " Time, altitude, dBZ-max, dBZ-I, udvm5, udv5, udv10, udv20\n");
    print OUT ( "Time (min), altitude,crfrzf,ciacrf,chcnsh,chcnih,qhacw,qracw,qrcnw,vhacw,ptem,....etc....,crfrz+ciarcr,qrfrzf+qiacrf,Rime Density,chcnsh+chcnih,chmul1+csplinter\n" );
      
      }
      $curbin = $curbin + 1;

     if ( $min[1] >= 1 ) {
     for ($kz = 1; $kz < $nz; $kz++ ) {
    printf OUT ( "%7.2f,  %6.3f,  %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n",
      $min[1],$z[$kz],$procth[$kz][1],$procth[$kz][2],$procth[$kz][3],$procth[$kz][4],$procth[$kz][5],$procth[$kz][6],$procth[$kz][7],$procth[$kz][8],$procth[$kz][9],$procth[$kz][10],$procth[$kz][11],
       $procth[$kz][12],$procth[$kz][13],$procth[$kz][14],$procth[$kz][15],$procth[$kz][16],$procth[$kz][17],$procth[$kz][18],$procth[$kz][19],$procth[$kz][20],$procth[$kz][21],$procth[$kz][22],
       $crfrzf[$kz]+$ciacrf[$kz],$qrfrzf[$kz]+$qiacrf[$kz], $rimdens[$kz], $chcnsh[$kz]+$chcnih[$kz], $procth[$kz][10]+$procth[$kz][11],$raintypefrz[$kz][1],$raintypefrz[$kz][2],$raintypefrz[$kz][3],
       $chlcnh[$kz],$chlcnf[$kz],$cwctfz[$kz],$ciint[$kz]);

      }  }
     
     } # if
      
      
    } # while
   } # foreach
   
