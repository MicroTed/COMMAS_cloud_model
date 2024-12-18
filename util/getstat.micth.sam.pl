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
    if ( /^nx,ny,nz,na =/ ) {
      chomp();
      @parts = split(/ +/);
      $nx = $parts[2];
      $ny = $parts[3];
      $mult = ($nx-1)*($ny-1);
      $nz = $parts[4];
      
      
    }
    if ( /^ dtg/ ) {
      chomp();
      @parts = split(/ +/);
      $dt = $parts[2];
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dx = $parts[2];
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dy = $parts[2];
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dz = $parts[2];
      $dv = $dx*$dy*$dz;
      
#   for ($kz = 1; $kz < $nz; $kz++ ) {      
#      $z[kz] = ($kz - 0.5)*$dz*0.001;
#      }
      
      }
# istrx,y,z=           0           0           1
    if ( /^ istrx/ ) {
      chomp();
      @parts = split(/ +/);
      $istrz = $parts[4];
      }
# STRETCH FAC 
    if ( /^ KZ,ZTRAN,GT/ ) {
#      $next = <INPUT>;
#      chomp($next);
# kz, scal-hgt, w-hgt, mfc, mfe, 0.5*(mfc[k-1]+mfc[k]) 0.5*(mfe[k+1]+mfe[k])
#   1     80.      5.  1.4000  1.4000  1.4000  1.3969
       
     for ($kz = 1; $kz < $nz; $kz++ ) {      
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
       $z[$kz] = $parts[3]*0.001;
       $gt[$kz] = $parts[4];
       $dzgt[$kz] = $dz/$gt[$kz];
   #    print OUT ( "kz,z $kz $dzgt[$kz] \n" );
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
      print OUT ( " Time, Altitude, qcmx, qrmx, qimx, qirmx, qsmx, ",
      "qglmx, qgmmx, qghmx, qfmx, qhmx, qipmx, qhlmx ",
      "qcmass, qrmass, qimass, qirmass, qsmass, qglmass, qgmmass, ",
      "qghmass, qfmass, qhmass, qipmass, qhlmass\n");
      
      }
  # $qmax[4][$kz]*$dens[$kz], cwcmx, 
      $time[$curbin] = $parts[2];
      if ( $time[$curbin] > $time[$curbin-1]+$dt) {
#      print OUT ( "time = $time[$curbin] \n" );
     for ($kz = 1; $kz < $nz; $kz++ ) {
# output changed 3/19/07 by erm to correct for index changes
    printf OUT ( "%5.2f,  %6.3f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n",
      $time[$curbin]/60., $z[$kz],
      $qmax[4][$kz],$qmax[5][$kz],$qmax[6][$kz], $qmax[8][$kz], $qmax[9][$kz],
      $qmax[10][$kz],$qmax[11][$kz],$qmax[12][$kz],$qmax[13][$kz],$qmax[14][$kz],
      $qmax[7][$kz],$qmax[15][$kz],
      $mass[4][$kz],$mass[5][$kz],$mass[6][$kz],$mass[8][$kz],$mass[9][$kz],
      $mass[10][$kz],$mass[11][$kz],$mass[12][$kz],$mass[13][$kz],$mass[14][$kz],
      $mass[7][$kz],$mass[15][$kz]
# old indices:
#       $qmax[4][$kz],$qmax[5][$kz],$qmax[6][$kz], $qmax[7][$kz], $qmax[8][$kz],
#       $qmax[9][$kz],$qmax[10][$kz],$qmax[11][$kz],$qmax[12][$kz],$qmax[13][$kz],
#       $qmax[14][$kz],$qmax[15][$kz],
#       $mass[4][$kz],$mass[5][$kz],$mass[6][$kz],$mass[7][$kz],$mass[8][$kz],
#       $mass[9][$kz],$mass[10][$kz],$mass[11][$kz],$mass[12][$kz],$mass[13][$kz],
#       $mass[14][$kz],$mass[15][$kz]
      
       );
      
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
    if ( /^ Layer average, max & min for scalars on grid/ ) {

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $ia = $parts[3];
      
      if ( $ia >= 4 && $ia <= 15 ) {

      for ($kz = $nz-1; $kz > 0; $kz-- ) { 
       $next = <INPUT>;
       chomp($next);

       @parts = split(/ +/, $next);
      
      $qmax[$ia][$kz] = $parts[7];
      # NOTE THAT there is already a factor of density in the 'average' value!
      $mass[$ia][$kz] = $parts[6]*$dzgt[$kz]*$dx*$dy*$mult;
      } # for

      }
     } # if
    } # while
   } # foreach
   
