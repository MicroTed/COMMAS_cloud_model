#!/usr/bin/perl
#
# 4/2007: Now time-height lightning data is in the stat file!
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
    $bintime = 0;
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
    $offtimeflag = 0;
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
     
     if ( /^ Time Altitude Channels/ ) {
     $title = $_;
     for ($kz = $nz-1; $kz >= 1; $kz-- ) { 
      $next = <INPUT>;
      $_ = $next;
      $dataline[$kz] = $next;
 #     print OUT ( "$dataline[$kz]");
 #     chomp($next);

      chomp($next);
      @parts = split(/ +/, $next);
      $numdat = $#parts;

#kz,ctghsnz,ctghspz,ctghinz,ctghipz,ctghwnz,ctghwpz
      if ( /^ / ) { $i = 1 } else { $i = 0 }
       $min[$kz] = $parts[0+$i];
     #  print OUT ("minute,numdat = $min[$kz], $numdat, $parts[0], @parts[1..$numdat] \n");
#      $z[$kz] = $parts[1+$i];

      for ( $j = $i+1; $j < $numdat+1; $j++ ) {
        $thparts[$kz][$j] = $thparts[$kz][$j] + $parts[$j];
       # print OUT ("thparts $kz, $j, is $thparts[$kz][$j]\n");
        }
      
      if ( $kz == $nz-1 ) {
        $bintime[$curbin] = $min[$kz];
     #   print OUT ("bintime = $bintime[$curbin]\n");
        $offtimeflag = 0;
        if ( $bintime[$curbin] - int($bintime[$curbin]) != 0  ) {
     #     print OUT ("found an off time! $bintime[$curbin], $curbin\n");
          $offtimeflag = 1;
          }
        }
      
      } # for kz

      if ( $curbin == 0 ) {
      print OUT ( "$title" );
     # print OUT ( " Time, altitude, dBZ-max, dBZ-I, udvm5, udv5, udv10, udv20\n");
      
      }
      $curbin = $curbin + 1;

     if ( $min[1] >= 1 && $offtimeflag == 0) {
     for ($kz = 1; $kz < $nz; $kz++ ) {
#      print OUT ( "$dataline[$kz]");
#       print OUT ( "$thparts[$kz][$i..$numdat]\n");
       print OUT (" $min[1] ");
       print OUT "@$_\n" for $thparts[$kz];
#    printf OUT ( "%7.2f,  %6.3f, %9.2f, %9.2f, %e, %e, %e, %e\n",
#      $min[$kz], $z[$kz], $dbzmax[$kz],$dbzsum[$kz],$udvm5[$kz],
#      $udv5[$kz],$udv10[$kz],$udv20[$kz] );

      } 
      @thparts = 0;
      
      }
     
     } # if
      
    } # while
   } # foreach
   
