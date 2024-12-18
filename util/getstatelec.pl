#!/usr/bin/perl
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
      if ( $curbin == 1 ) {
#      print OUT ( " Time Altitude udmf uicmf grvol hlvol grms hlms ticms rnms\n");
      print OUT ( " Time Altitude cgizmn cgizmx cghwzmn cghwzmx ctghsinz ctghsipz ctghwnz ctghwpz\n");
      
      }
      $time[$curbin] = $parts[8];
      if ( $time[$curbin] > $time[$curbin-1]) {
#      print OUT ( "time = $time[$curbin] \n" );
     for ($kz = 1; $kz < $nz; $kz++ ) {
    printf OUT ( "%5.2f  %6.3f %e %e %e %e %e %e %e %e \n",
      $time[$curbin]/60., $z[$kz], $cgizmn[$kz],$cgizmx[$kz],$cghwzmn[$kz],
      $cghwzmx[$kz],$ctghinz[$kz],$ctghipz[$kz],$ctghwnz[$kz],$ctghwpz[$kz] );
#      $time[$curbin]/60., $z[$kz], $udmf[$kz],$uicmf[$kz],$grvol[$kz],
#      $hlvol[$kz],$grms[$kz],$hlms[$kz],$ticms[$kz],$rnms[$kz] );

    }}
      }
    
    if ( /^kz,cgizmx,cgizmn/ ) {
     $curbin = $curbin+1;

     for ($kz = $nz-1; $kz > 0; $kz-- ) { 
      $next = <INPUT>;
      chomp($next);

      @parts = split(/ +/, $next);
#kz,cgizmx,cgizmn,cgszmx,cgszmn,chizmx,chizmn,chszmx,chszmn,cghwzmx,cghwzmn
      
      $cgizmx[$kz] = $parts[2];
      $cgizmn[$kz] = $parts[3];
      $cgszmx[$kz] = $parts[4];
      $cgszmn[$kz] = $parts[5];
      $chizmx[$kz] = $parts[6];
      $chizmn[$kz] = $parts[7];
      $chszmx[$kz] = $parts[8];
      $chszmn[$kz] = $parts[9];
      $cghwzmx[$kz] = $parts[10];
      $cghwzmn[$kz] = $parts[11];
      } # for

#      $next = <INPUT>; # Layer volumes/masses/rat
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
      $ctghinz[$kz] = $parts[4] + $parts[2];
      $ctghipz[$kz] = $parts[5] + $parts[3];
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
#      $swms[$kz] = $parts[10] ; # snow mass
      } # for

     
     
     } # if
    if ( /^wmax/ ) {
      chomp();
      @parts = split( / +/ );
      $wmax = $parts[2];
      $iwmax = $parts[3];
      $jwmax = $parts[4];
      $kwmax = $parts[5];
      }
    if ( /^ EFIELD--E-MAX/ ) {
      chomp();
      @parts = split( / +/ );
      $efield = $parts[$#parts];
      }
   if ( $nstep < $nstop - 1 )  {
    if ( /ctswin,ctswip/ ) {
      chomp();
      @parts = split(/ +|,/);
      $nonindn = $nonindn + $parts[3];
      $nonindp = $nonindp + $parts[5];

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $nonindn = $nonindn + $parts[3];
      $ghsn = $ghsn + $parts[3];
      $nonindp = $nonindp + $parts[5];
      $ghsp = $ghsp + $parts[5];

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $nonindn = $nonindn + $parts[3];
      $ghin = $ghin + $parts[3];
      $nonindp = $nonindp + $parts[5];
      $ghip = $ghip + $parts[5];

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $indn = $indn + $parts[3];
      $indp = $indp + $parts[5];     
      }
    if ( /^Try/ ) {
     $tries = $tries + 1;
     $bintry = $bintry + 1;
     }
    if ( /IC DISCHARGE/ ) {
      $ic = $ic + 1;
      $type = 1;
     }
    if ( /DISCHARGE IS POSITIVE/ ) {
      $cgp = $cgp + 1;
      $type = 2;
     }
    if ( /DISCHARGE IS NEGATIVE/ ) {
      $type = 3;
      $cgn = $cgn + 1;
     }
    if (/^check totals/ ) {
     @parts = split(/ +/ );
     $pos = $parts[2];
     $neg = $parts[3];
        if ( $type == 1 ) {
           $iccharge = $iccharge + $pos + $neg;
           if ( abs($pos) < abs($neg) ) {
             $icdis = $icdis + abs($pos);
            } else {
             $icdis = $icdis + abs($neg); 
            }
         } elsif ( $type == 2 ) {
           $cgpcharge = $cgpcharge + $pos + $neg;
         } elsif ( $type == 3 ) {
           $cgncharge = $cgncharge + $pos + $neg;
            }
     }    
     } # if ( $nstep < $nstop - 1)  
    } # while
   } # foreach
   
