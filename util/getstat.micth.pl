#!/usr/bin/perl
# a perl program to read in stat file(s) and get time-height data
#   of max mixing ratio and total mass for each hydrometeor category
# Looks for the lines:
#
# Writes to standard output:
 open(OUT, '>-');
# 
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
    $time[$curbin] = 0.0;
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
    $imicro = 3;
    $nrain = 0;
    $layercheck = 0;
    $qf = 0;
    $qhl = 0;
    $qh = 0;
    $qs = 0;
    $qi = 0;
    $qr = 0;
    $qc = 0;
    $qrauto = 0;
    $qrshed = 0;
    $qrmelt = 0;
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
      $nx = $parts[2];
      $ny = $parts[3];
      $mult = ($nx-1)*($ny-1);
      $nz = $parts[4];
      
      # set zero values
      for ($kz = $nz-1; $kz > 0; $kz-- ) { 
         $qmax[0][$kz] = 0; 
         $mass[0][$kz] = 0;
         }
      
    }
    if ( /^dt,dx,dy,dz =/ ) {
      chomp();
      @parts = split(/ +/);
      $dt = $parts[2];
      $dx = $parts[3];
      $dy = $parts[4];
      $dz = $parts[5];
#      $dv = $dx*$dy*$dz;
      
#   for ($kz = 1; $kz < $nz; $kz++ ) {      
#      $z[kz] = ($kz - 0.5)*$dz*0.001;
#      }
      
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
        if ( $string eq 'QRAUTO' ) {$qrauto = $parts[4]; $nrain = $nrain + 1; } # print OUT ( "qrauto = $qrauto\n"); }
        if ( $string eq 'QRSHED' ) {$qrshed = $parts[4]; $nrain = $nrain + 1; } # print OUT ( "qrshed = $qrshed\n"); }
        if ( $string eq 'QRMELT' ) {$qrmelt = $parts[4]; $nrain = $nrain + 1; } # print OUT ( "qrmelt = $qrmelt\n"); }
        if ( $string eq 'QC' ) { $qc = $parts[4]; }
        if ( $string eq 'QR' ) { $qr = $parts[4]; }
        if ( $string eq 'QI' ) { $qi = $parts[4]; }
        if ( $string eq 'QS' ) { $qs = $parts[4]; }
        if ( $string eq 'QH' ) { 
           $qh = $parts[4]; # print OUT ( "qh = $qh\n"); 
           }
        if ( $string eq 'QF' ) { $qf = $parts[4]; #print OUT ( "qf = $qf\n"); 
            }
        if ( $string eq 'QHL' ) { $qhl = $parts[4]; #print OUT ( "qhl = $qhl\n"); 
          }
       # print OUT ( "$string, $len, $nrain \n");
    
    }

    if ( /^micro =/ ) {
      chomp();
      @parts = split(/ +/);
      $micro = $parts[2];
#      print OUT ("micro = $micro, $dx, $dy, $dz, nx=$nx, ny=$ny, mult=$mult\n");
#      $test =  substr( $micro, 0, 3 );
#      print OUT ("micro = $test\n");
      if ( substr( $micro, 0, 3 ) eq 'ICE' ) {
#        print OUT ("micro = ICE\n");
        $imicro = 2;
      }
      if ( substr( $micro, 0, 1 ) eq 'Z' ) {
#        print OUT ("micro = Z\n");
        $imicro = 3;
      }
#      die;
    }

    if ( /^lqb/ ) {
      chomp();
      @parts = split(/ +/);
      $lqb = $parts[1];
      $lqe = $parts[2];
#      print OUT ("lqb,lqe = $lqb, $lqe\n");
      }

# istrx,y,z=           0           0           1
    if ( /^ istrx/ ) {
      chomp();
      @parts = split(/ +/);
      $istrz = $parts[4];
      }

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
       $dzgt[$kz] = $parts[3];
#       print OUT ( "kz,z $kz $z[$kz], $dzgt[$kz] \n" );
      }
       }
    
# sounding
# kz  z       p        pi     den    theta   tem
    if ( /^ kz  z       p        pi     den    theta/ ) {
       
     for ($kz = 1; $kz < $nz; $kz++ ) {      
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dens[$kz] = $parts[5];
  #    print OUT ( "kz,z $kz $dens[$kz] \n" );
      }
       
       }

      
#     if ( /^ pmax/ ) {
#     if ( /^ Max and Mins for w at level 1/ ) {
     if ( /^ wmax/ ) {
      chomp();
      @parts = split(/ +/);
      $curbin = $curbin + 1;
      if ( $curbin == 1 ) {
       if ( $imicro == 2 ) {
      print OUT ( " Time, Altitude, qcmx, qrmx, qimx, qirmx, qsmx, ",
      "qglmx, qgmmx, qghmx, qfmx, qhmx, qipmx, qhlmx, ",
      "qcmass, qrmass, qimass, qirmass, qsmass, qglmass, qgmmass, ",
      "qghmass, qfmass, qhmass, qipmass, qhlmass, ",
      "qrautomx, qrshedmx,qrmeltmx, qrautomass, qrshedmass, qrmeltmass\n");
       } elsif ( $imicro == 3 ) {
      print OUT ( " Time, Altitude, qcmx, qrmx, qimx, qirmx, qsmx, ",
      "qrautomx, qrshedmx,qrmeltmx, qfmx, qhmx, qipmx, qhlmx, ",
      "qcmass, qrmass, qimass, qirmass, qsmass, qrautomass, qrshedmass, qrmeltmass, ",
      "qfmass, qhmass, qipmass, qhlmass\n");

#      print OUT ( " Time, Altitude, qcmx, qrmx, qimx, qsmx, ",
#      "qgmx, qcmass, qrmass, qimass,  qsmass, qgmass");
#       
       #}
       }
      }
  # $qmax[4][$kz]*$dens[$kz], cwcmx, 
      $time[$curbin] = $parts[2];
      if ( $time[$curbin] > $time[$curbin-1]+$dt) {
#      print OUT ( "time = $time[$curbin] \n" );
     for ($kz = 1; $kz < $nz; $kz++ ) {
# output changed 3/19/07 by erm to correct for index changes
       if ( $imicro == 2 ) {
    printf OUT ( "%5.2f,  %6.3f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n",
      $time[$curbin]/60., $z[$kz],
      $qmax[3][$kz],$qmax[4][$kz],$qmax[5][$kz], $qmax[7][$kz], $qmax[8][$kz],
      $qmax[9][$kz],$qmax[10][$kz],$qmax[11][$kz],$qmax[12][$kz],$qmax[13][$kz],
      $qmax[6][$kz],$qmax[14][$kz],
      $mass[3][$kz],$mass[4][$kz],$mass[5][$kz],$mass[7][$kz],$mass[8][$kz],
      $mass[9][$kz],$mass[10][$kz],$mass[11][$kz],$mass[12][$kz],$mass[13][$kz],
      $mass[6][$kz],$mass[14][$kz]
       );
# old indices:
#       $qmax[4][$kz],$qmax[5][$kz],$qmax[6][$kz], $qmax[7][$kz], $qmax[8][$kz],
#       $qmax[9][$kz],$qmax[10][$kz],$qmax[11][$kz],$qmax[12][$kz],$qmax[13][$kz],
#       $qmax[14][$kz],$qmax[15][$kz],
#       $mass[4][$kz],$mass[5][$kz],$mass[6][$kz],$mass[7][$kz],$mass[8][$kz],
#       $mass[9][$kz],$mass[10][$kz],$mass[11][$kz],$mass[12][$kz],$mass[13][$kz],
#       $mass[14][$kz],$mass[15][$kz]
       } elsif ( $imicro == 3 ) {
    printf OUT ( "%5.2f,  %6.3f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n",
      $time[$curbin]/60., $z[$kz],
      $qmax[$qc][$kz],$qmax[$qr][$kz],$qmax[$qi][$kz],0, $qmax[$qs][$kz], $qmax[$qrauto][$kz], 
      $qmax[$qrshed][$kz], $qmax[$qrmelt][$kz],$qmax[$qf][$kz], $qmax[$qh][$kz],
      0,$qmax[$qhl][$kz],
      $mass[$qc][$kz],$mass[$qr][$kz],$mass[$qi][$kz],0, $mass[$qs][$kz],$mass[$qrauto][$kz], 
      $mass[$qrshed][$kz], $mass[$qrmelt][$kz],$mass[$qf][$kz], $mass[$qh][$kz],
      0,$mass[$qhl][$kz]
      
       );
#      $mass[10][$kz],$mass[11][$kz],$mass[12][$kz],$mass[13][$kz],$mass[14][$kz],
#      $mass[7][$kz],$mass[15][$kz]
      
      }
      
#      $udmf[$kz],$uicmf[$kz],$grvol[$kz],
#      $hlvol[$kz],$grms[$kz],$hlms[$kz],$ticms[$kz],$rnms[$kz],$swms[$kz] );

    }}
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

#
# Get hydrometeor time-height info if 4 <= ia <= 15
#
    if ( /^ Layer average, max & min for scalars/ ) {

#      print OUT ("Read scalars\n");
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $ia = $parts[3];
#      print OUT ("Read scalars, ia = $ia, $lqb, $lqe\n");
      
#      if ( $ia >= $lqb ) {  print OUT ("ia > lqb\n"); }
#      if ( $ia <= $lqe ) {  print OUT ("ia < lqe\n"); }
#      if ( $ia >= $lqb && $ia <= $lqe ) {  print OUT ("ia >< lqb\n"); }
      
      if ( ($ia >= $lqb && $ia <= $lqe) || ($ia >= $qrauto && $ia <= $qrmelt)  ) {
      
 
      for ($kz = $nz-1; $kz > 0; $kz-- ) { 
       $next = <INPUT>;
       chomp($next);

       @parts = split(/ +/, $next);
      
      $qmax[$ia][$kz] = $parts[7];
      # NOTE THAT there is already a factor of density in the 'average' value!
#      $mass[$ia][$kz] = $parts[6]*$dz*$dx*$dy; # *$mult;
      $mass[$ia][$kz] = $parts[6]*$dzgt[$kz]*$dx*$dy*$mult;
#      if ( $ia == $qh ) {
#       print OUT ("kz,qh = $kz, $qmax[$ia][$kz], $mass[$ia][$kz]\n" );
#      }
#      if ($ia >= $qrauto && $ia <= $qrmelt) {
#       print OUT ("kz,$ia = $kz, $qmax[$ia][$kz], $mass[$ia][$kz]\n" );
#      }
      } # for

      }
     } # if
    } # while
   } # foreach
   
