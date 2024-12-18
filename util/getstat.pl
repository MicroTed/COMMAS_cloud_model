#!/usr/bin/perl
#
# 9/2005: added time-height updraft/downdraft max values
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
    $printflag = 0;
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
    if ( /^STOP HERE/ ) {
      goto "foo";
      }
    
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
#    if ( /^ KZ,ZTRAN,GT/ ) {
#     if ( /^Time, altitude, dBZ-max/ ) {
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
   #    print OUT ( "kz,z $kz $z[$kz] \n" );
      }
       }
     if ( /^Max dBZ and integrated UDV at time =/ ) {
       chomp();
      @parts = split(/ +/);
      $time[$curbin] = $parts[8];
      
      $printflag = 1;

      if ( $curbin == 1 ) {
      print OUT ( " Time, Altitude, udmf, uicmf, grvol, hlvol, grms, hlms, ticms, rnms, snms, iwcmx, wmax, wmin\n");
      
      }
 #     $time[$curbin] = $parts[2];
      if ( $time[$curbin] > $time[$curbin-1]+$dt) {
#      print OUT ( "time = $time[$curbin] \n" );
     for ($kz = 1; $kz < $nz; $kz++ ) {
    printf OUT ( "%5.2f,  %6.3f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n",
      $time[$curbin]/60., $z[$kz], $udmf[$kz],$uicmf[$kz],$grvol[$kz],
      $hlvol[$kz],$grms[$kz],$hlms[$kz],$ticms[$kz],$rnms[$kz],$swms[$kz],$iwcmx[$kz],$wmax[$kz],$wmin[$kz] );

    }}

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
      
     if ( /^ pmax/ && $printflag == 0 ) {
      chomp();
      @parts = split(/ +/);
      if ( $curbin == 1 ) {
      print OUT ( " Time, Altitude, udmf, uicmf, grvol, hlvol, grms, hlms, ticms, rnms, snms, iwcmx, wmax, wmin\n");
      
      }
      $time[$curbin] = $parts[2];
      if ( $time[$curbin] > $time[$curbin-1]+$dt) {
#      print OUT ( "time = $time[$curbin] \n" );
     for ($kz = 1; $kz < $nz; $kz++ ) {
    printf OUT ( "%5.2f,  %6.3f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n",
      $time[$curbin]/60., $z[$kz], $udmf[$kz],$uicmf[$kz],$grvol[$kz],
      $hlvol[$kz],$grms[$kz],$hlms[$kz],$ticms[$kz],$rnms[$kz],$swms[$kz],$iwcmx[$kz],$wmax[$kz],$wmin[$kz] );

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
    if ( /^kz,cgizmx,cgizmn/ ) {
     $curbin = $curbin+1;

     for ($kz = $nz-1; $kz > 0; $kz-- ) { 
#     for ($kz1 = 1; $kz1 < $nz; $kz1++ ) 
      $next = <INPUT>;
      chomp($next);

      @parts = split(/ +/, $next);
      
      $cgizmn[$kz] = $parts[2];
      $cgizmx[$kz] = $parts[3];
      $cgszmn[$kz] = $parts[4];
      $cgszmx[$kz] = $parts[5];
      $chizmn[$kz] = $parts[6];
      $chizmx[$kz] = $parts[7];
      $chszmn[$kz] = $parts[8];
      $chszmx[$kz] = $parts[9];
      $cghwzmn[$kz] = $parts[10];
      $cghwzmx[$kz] = $parts[11];
      } # for


      $next = <INPUT>; # (list of titles)
      $_ = $next;
     
     if ( /^kz,ctghsnz,ctghspz/ ) {
     for ($kz = $nz-1; $kz > 0; $kz-- ) { 
      $next = <INPUT>;
      chomp($next);

      @parts = split(/ +/, $next);

#kz,ctghsnz,ctghspz,ctghinz,ctghipz,ctghwnz,ctghwpz
      
      $ctghsnz[$kz] = $parts[2];
      $ctghspz[$kz] = $parts[3];
      $ctghinz[$kz] = $parts[4];
      $ctghipz[$kz] = $parts[5];
      $ctghwnz[$kz] = $parts[6];
      $ctghwpz[$kz] = $parts[7];
      } # for

      $next = <INPUT>; # Layer volumes/masses/rat
      $next = <INPUT>; # kz,udmf,uicmf,grvol,hlvol,grms,hlms,ticms,rnms
      
      }  # if ( /^kz,ctghsnz/ )
       else { 
      $next = <INPUT>; # kz,udmf,uicmf,grvol,hlvol,grms,hlms,ticms,rnms
       }

     for ($kz = $nz-1; $kz > 0; $kz-- ) {      
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $udmf[$kz] = $parts[2] ; # updraft mass flux
      $uicmf[$kz] = $parts[3]; # upward ice crystal mass flux
      $grvol[$kz] = $parts[4]; # graupel volume
      $hlvol[$kz] = $parts[5]; # hail volume
      $grms[$kz] = $parts[6] ; # graupel mass
      $hlms[$kz] = $parts[7] ; # hail mass
      $ticms[$kz] = $parts[8]; # total ice crystal mass
      $rnms[$kz] = $parts[9] ; # total rain mass
      $swms[$kz] = $parts[10] ; # snow mass
      $iwcmx[$kz] = $parts[11] ; # max ice water content
      } # for

     
     
     } # if
        if ( /^kz,udmf,uicmf,grvol,hlvol/ ) {  # this test will pick up udmf etc. in the non-electrification case
     $curbin = $curbin+1;
     for ($kz = $nz-1; $kz > 0; $kz-- ) {      
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $udmf[$kz] = $parts[2] ; # updraft mass flux
      $uicmf[$kz] = $parts[3]; # upward ice crystal mass flux
      $grvol[$kz] = $parts[4]; # graupel volume
      $hlvol[$kz] = $parts[5]; # hail volume
      $grms[$kz] = $parts[6] ; # graupel mass
      $hlms[$kz] = $parts[7] ; # hail mass
      $ticms[$kz] = $parts[8]; # total ice crystal mass
      $rnms[$kz] = $parts[9] ; # total rain mass
      $swms[$kz] = $parts[10] ; # snow mass
      $iwcmx[$kz] = $parts[11] ; # max ice water content
      } # for
     } # if

# old    if ( /^ Layer average, max & min for velocity/ ) {
    if ( /^ Layer average, min & max for velocity/ ) {
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      if ( $parts[3] == 3 ) {
      for ($kz = $nz-1; $kz > 0; $kz-- ) {      
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $wmin[$kz] = $parts[7] ; # min updraft
      $wmax[$kz] = $parts[8] ; # max updraft
      } # for
      
      
      }
     }
    } # while
   foo:
   } # foreach
   
