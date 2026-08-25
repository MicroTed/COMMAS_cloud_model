#!/usr/bin/perl
# a perl program to read in a .out file and total rain and hail rates and accumlations over time
# 
# usage:
# rainfall.pl run.000.000.out > run.rainfall.txt
#
#
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
    $intervalm = 1;   # time discretization (usually 1 minute)
    $intervals = $intervalm*60;  
    $istop = 0;
    $binshift = 0;
    $tries = 0;
    $bintry = 0;
    $newbin = 1;
    $curbin = 0;
    $ic = 0;
    $cgn = 0;
    $cgp = 0;
    $icx = 0;
    $cgnx = 0;
    $cgpx = 0;
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
    $netchg = 0;
    $netchgp = 0;
    $netchgn = 0;
    $irestart = 0;
    $dv = 1.0; 
    $irst2 = 0;
    $readfirst = 0;
    $dt = 4;
      $ionnet   = 0;
      $iontb    = 0;
      $times = 0;
      $timemin = 0;
      $energy = 0;
      $lgtflag = 0;
      
      $raintot = 0;
      $raintotold = 0;
      $hailtot = 0;
      $hailtotold = 0;
      $printflag = 0;

     print ("Time (s), Time (min), Rain Rate, Ice Rate, Hail Rate, Hail-FD Rate, Rain Accum., Ice Accum, Hail Accum., Hail-FD Accum.\n");

# foreach ...
foreach $file (@ARGV) {
    open INPUT, "<$file";
    $start = 0;
    while ( (<INPUT>)  ) {
    if ( /^STOP HERE/ ) {
      goto "foo";
      }
    if ( /^ TIME STEP  / ) {
      chomp();
      @parts = split(/ +/);
      $dt = $parts[4];
#      print("dt = $dt\n");
#      $next = <INPUT>;
#      chomp($next);
#      @parts = split(/ +/, $next);
       }
    if ( /^ NX / ) {
      chomp();
      @parts = split(/ +/);
      $dx = $parts[8];
      $next = <INPUT>;
      $next = <INPUT>;
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dy = $parts[8];
      $next = <INPUT>;
      $next = <INPUT>;
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +/, $next);
      $dz = $parts[8];
 #     print("dx = $dx\n");
 #     print("dy = $dy\n");
 #     print("dz = $dz\n");
      $dv = 1.0; #$dx*$dy*$dz;
      }
    if ( / SIMULATION STOP  TIME/ && $start == 0 ) {
      chomp();
      @parts = split( / +/ );
      $nstop = $parts[$#parts];
#      print("tstop = $nstop\n");
#      print("nstop=$nstop\n");
      }
    if ( /^  COMMAS:   RESTART/ ) {
      $irestart = 1;
      $irst2 = 1;
     }
    if ( /^ COMMAS:  START =/ && $irestart == 1 ) {
      $irestart = 0;
      chomp();
      @parts = split(/ +/);
      $times1 = $parts[4];
      # check if the restart is a continuation or goes back in time
      if ( $times1 < $times ) {
      print ("restart: $times, $times1\n");
        $min = int($times1/$intervals);
        $curbin = $min - $binshift;
        $newbin = $curbin;
        $nel = $#nonindn;
        print ("restart: min,curbin,binshift = $min, $curbin, $binshift\n");
        print ("current number of times = $nel\n");
        @nonindn[$curbin..$nel] = 0;
        @nonindp[$curbin..$nel] = 0;
        @wmax[$curbin..$nel] = 0;
        @ic[$curbin..$nel] = 0;
        @tries[$curbin..$nel] = 0;
        @icdis[$curbin..$nel] = 0;
        $ic = 0;
        $icdis = 0;
        $nonindn = 0;
        $nonindp = 0;
        $tries = 0;
        print ("restart: curbin = $curbin\n");
       } 
    
     }
     
    if ( /^END OF STEP: nstep, time/ ) {
      @parts = split(/ +/);
      $nsteptst = $parts[6];
#      print ("end of step: $nsteptst, $nstop, $min, $curbin, $binshift\n");
      if ( $nsteptst == $nstop && $curbin == 1 ) {
       $min = $min+1;
       $newbin = $min - $binshift;
       }
#       if ( $irst2 == 1 ) {  print ("HERE I AM: $newbin, $curbin\n"); }
     }
    if ( /^NSTEP,/ ) {
     $next = <INPUT>;
     chomp($next);
     @parts = split(/ +/, $next);
     $nstep = $parts[1];
     $nstop = $parts[3];
     $timesec = $parts[4];
    
      }
    
# Total rainfall =   0.0000000E+00
# Total hailfall =   0.0000000E+00
# -----------------------------------------------------------------------
# T =    120.00    MAX                            MIN
    if ( /^ Total rainfall/ ) {

      chomp();
      @parts = split( / +/ );
      $raintot = $parts[4];
      $timemin = $timesec/60.;
      
      $delrain = $raintot-$raintotold;
      $raintotold = $raintot;
      $printflag = 1;

      }
      
    if ( /^ Total hailfall/ ) {
#       $next = <INPUT>;
#       chomp($next);
#       @parts = split(/ +/,$next);
      chomp();
      @parts = split( / +/ );
      $icetot = $parts[4];
      $delice = $icetot-$icetotold;
      $icetotold = $icetot;

      }
    if ( /^ Total hail =/ ) {
#       $next = <INPUT>;
#       chomp($next);
#       @parts = split(/ +/,$next);
      chomp();
      @parts = split( / +/ );
      $hailtot = $parts[4];
      $delhail = $hailtot-$hailtotold;
      $hailtotold = $hailtot;
      }

    if ( /^ Total hail from FD/ ) {
#       $next = <INPUT>;
#       chomp($next);
#       @parts = split(/ +/,$next);
      chomp();
      @parts = split( / +/ );
      $hailfdtot = $parts[6];
      $delhailfd = $hailfdtot-$hailfdtotold;
      $hailfdtotold = $hailfdtot;

      }
      
     if ( /^ T =  /) {
      if ( $printflag > 0 ) {
      print OUT ("$timesec, $timemin, $delrain, $delice, $delhail, $delhailfd, $raintot, $icetot, $hailtot, $hailfdtot\n");
      $printflag = 0;
      }
    
      }
    
    if (  /^ Integration Done/ ) {
      if ( $newbin == $curbin ) { 
         $newbin = $curbin + 1; 
         $min = int(($times)/$intervals); 
         }
     }
    if (  $newbin > $curbin) {

     $bintime[$curbin] = $min*$intervalm;
#     print OUT ("curbin, bintime, min, newbin = $curbin, $bintime[$curbin], $min, $times, $newbin\n");

     $ictot = $ictot + $ic;
     $ic[$curbin] = $ic;
     $icpermin[$curbin] = $ic/$intervalm;
     $ic = 0;
     
     $cgptot = $cgptot + $cgp;
     $cgp[$curbin] = $cgp;
     $cgp = 0;
     
     $cgntot = $cgntot + $cgn;
     $cgn[$curbin] = $cgn;
     $cgn = 0;
     
     $triestot = $triestot + $tries;
     $tries[$curbin] = $tries;
     $tries = 0;
     
     $icdistot = $icdistot + $icdis; # *$dv;
     $icdis[$curbin] = $icdis; # *$dv;
     $icdis = 0;

     $icchargetot = $icchargetot + $iccharge; #*$dv;
     $iccharge[$curbin] = $iccharge; #*$dv;
     $iccharge = 0;
     
     $cgpchargetot = $cgpchargetot + $cgpcharge; # *$dv;
     $cgpcharge[$curbin] = $cgpcharge; # *$dv;
     $cgpcharge = 0;

     $cgnchargetot = $cgnchargetot + $cgncharge; # *$dv;
     $cgncharge[$curbin] = $cgncharge; # *$dv;
     $cgncharge = 0;
     
     $nonindtotn = $nonindtotn + $nonindn*$dt; # *$dv
     $nonindn[$curbin] = $nonindn*$dt; # *$dv
#      print ("nonindn: $nonindn[$curbin], $nonindn ($curbin)\n");
     $nonindn = 0;

     $nonindtotp = $nonindtotp + $nonindp*$dt; # *$dv
     $nonindp[$curbin] = $nonindp*$dt;  # *$dv
     $nonindp = 0;
     
     $ghsn[$curbin] = $ghsn*$dt; # *$dv
     $ghsn = 0;

     $ghsp[$curbin] = $ghsp*$dt;  # *$dv
     $ghsp = 0;

     $ghin[$curbin] = $ghin*$dt;  # *$dv
     $ghin = 0;

     $ghip[$curbin] = $ghip*$dt; # *$dv
     $ghip = 0;

     $indtotn = $indtotn + $indn*$dv*$dt;
     $indn[$curbin] = $indn*$dv*$dt;
     $indn = 0;

     $indtotp = $indtotp + $indp*$dv*$dt;
     $indp[$curbin] = $indp*$dv*$dt;
     $indp = 0;
     
#     $wmax[$curbin] = $wmax;

     if ( $intervalm == 1 ) {
     $wmax[$curbin] = $wmax;
     } else {
     $wmax[$curbin] = $wmaxtot/$intervals;
     $wmaxtot = 0;
     }

     $efield[$curbin] = $efield;

     $efield2[$curbin] = $efield2;
     $efield2 = 0;

     $netchg[$curbin] = $netchg;
     $netchgp[$curbin] = $netchgp;
     $netchgn[$curbin] = $netchgn;
#     print ("$min, $curbin, $times \n");
#       if ( $irst2 == 1 ) {  print ("HERE I AM 3: $newbin, $curbin\n"); }
     $curbin = $newbin;
      }
     
    if ( /^wmax/ ) {
      chomp();
      @parts = split( / +/ );
      $wmax = $parts[2];
      $wmaxtot = $wmaxtot + $wmax*$dt;
      }
     if ( /^net1 pos/ ) {
      chomp();
      @parts = split( / +/ );
      $chgpos = $parts[4];
      $chgneg = $parts[5];
      $chgnet = $parts[6];
      }
     if ( /^Energy info: old,new,difference/ ) {
      chomp();
      @parts = split( / +/ );
      $energytmp = $parts[4];
      $delenergytmp = $parts[5];
      $energylgt = 0;
      $delenergylgt = 0;
       if ( $icyn > 0 || $cgnyn > 0 || $cgpyn > 0 ) {
          $timelgt = $times + 0.1;
          $timelgtmin = $timelgt/60;
      #      $times = $times + 0.1; 
          $lgtflag = 1; 
          $energylgt = $energytmp;
          $delenergylgt = $delenergytmp;
         } else {
          $energy = $energytmp;
          $delenergy = $delenergytmp;
         }
       
      if ( $lgtflag == 1 ) {
#      print ("$times, $energy, $delenergy, $icyn, $cgnyn, $cgpyn, $chgpos, $chgneg, $chgnet, $efield\n");
      }
#        $icyn = 0; $cgnyn = 0; $cgpyn = 0;
      }

    if ( /^Net Ion charge fluxes/ ) {
#     print ;
      $ion = 1;
      chomp();
      @parts = split( / +/ );
      $ionsides = $ionsides + $parts[7] - $parts[8] - $parts[9];
      $ionnet   = $ionnet + $parts[7];
      $iontb    = $iontb + $parts[8] + $parts[9];
      $ionbottom = $ionbottom + $parts[9];
      $iontop = + $iontop + $parts[8];
#      print ("$times, $energy, $delenergy, $icyn, $cgnyn, $cgpyn, $chgpos, $chgneg, $chgnet, $efield\n");
      }

# postlight pos/neg/net charge (C):   1.29770E+03  -1.36479E+03  -6.70915E+01
    if ( /^postlight/ ) {
      chomp();
      @parts = split( / +/ );
      $netchgp = $parts[4];
      $netchgn = $parts[5];
      $netchg = $parts[6];
      
#      $ionchg2 = $ionchg2 + $netchg - $netchg1;
#      $ionchgcum2 = $ionchgcum2 + $netchg - $netchg1;
#      $ionerrtot = $ionerrtot + $netchg - $netchg1 - $scrchg0;
      }
    if ( /^ EFIELD--E-MAX/ ) {
      chomp();
      @parts = split( / +/ );
      $efield = $parts[$#parts];
        if ( $readfirst == 1 ) {
          $readfirst = 0;
          $efield2 = Max($efield2, $efield);
        }
      }

    if ( /^Charge fallout/ ) {
        $readfirst = 1;
      }
    if ( /^ EFIELD2--E-MAX/ ) {
      chomp();
      @parts = split( / +/ );
#      $efield2 = $parts[$#parts];
      }
   if ( $nstep < $nstop  )  {
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
#      print ("ghin/ghip: $ghin, $ghip\n");
#      print ("nonindn/nonindp: $nonindn, $nonindp\n");

      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
      $indn = $indn + $parts[3];
      $indp = $indp + $parts[5];     
      }
    if ( /^Try number/ ) {
     $tries = $tries + 1;
     $bintry = $bintry + 1;
     }
    if ( /IC DISCHARGE/ ) {
      $ic = $ic + 1;
      $type = 1;
     $icyn = 1;
   #   print OUT ("ic = $ic time = $times \n");
     }
    if ( /DISCHARGE IS POSITIVE/ ) {
      $cgp = $cgp + 1;
   #   print OUT ("cgp = $cgp time = $times \n");
      $type = 2;
      $cgpyn = 1;
     }
    if ( /DISCHARGE IS NEGATIVE/ ) {
      $type = 3;
      $cgn = $cgn + 1;
      $cgnyn = 1;
     }
#    if (/^check totals/ ) {
#      @parts = split(/ +/ );
#      $pos = $parts[2];
#      $neg = $parts[3];
#      if ( ($pos < 0.0001) && (($neg) > -0.0001) ) {
#         $pos = $pos*$dv;
#         $neg = $neg*$dv;
# #       print OUT ( "This is odd: pos = $pos, neg = $neg\n");
#         }
#         if ( $type == 1 ) {
#            $iccharge = $iccharge + $pos + $neg;
#            if ( abs($pos) < abs($neg) ) {
#              $icdis = $icdis + abs($pos);
#             } else {
#              $icdis = $icdis + abs($neg); 
#             }
#          } elsif ( $type == 2 ) {
#            $cgpcharge = $cgpcharge + $pos + $neg;
#   #         print OUT ("cgp : $pos $neg\n");
#          } elsif ( $type == 3 ) {
#            $cgncharge = $cgncharge + $pos + $neg;
#             }
#      }
    if (/^ADJ COULOMBS SCNETP/ ) {
     @parts = split(/ +/ );
     $pos = $parts[3];
      $next = <INPUT>;
      chomp($next);
      @parts = split(/ +|,/,$next);
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
  #         print OUT ("cgp : $pos $neg\n");
         } elsif ( $type == 3 ) {
           $cgncharge = $cgncharge + $pos + $neg;
            }
     }
    if (/^ WARNING: Lightning increased the total energy/) {
        if ( $type == 1 ) {
           $iccharge = $iccharge - $pos - $neg;
           $ic = $ic - 1;
           $icx = $icx + 1;
           if ( abs($pos) < abs($neg) ) {
             $icdis = $icdis - abs($pos);
            } else {
             $icdis = $icdis - abs($neg); 
            }
         } elsif ( $type == 2 ) {
           $cgpcharge = $cgpcharge - $pos - $neg;
           $cgp = $cgp - 1;
           $cgpx = $cgpx + 1;
         } elsif ( $type == 3 ) {
           $cgncharge = $cgncharge - $pos - $neg;
           $cgn = $cgn - 1;
           $cgnx = $cgnx + 1;
            }
     }    
     } # if ( $nstep < $nstop - 1)  
    } # while
  foo: 
   } # foreach
   
#   if ( $icx > 0 || $cgpx > 0 || $cgnx > 0 ) {}
#    print OUT ( "Removed $icx ICs, $cgpx +CGs, $cgnx -CGs\n" );
    
    $ictot = 0; 
    $cgptot = 0;  
    $cgntot = 0;
    $triestot = 0;
    $icdistot = 0;
    $cgpchargetot = 0;
    $cgnchargetot = 0;
    
#    print ("curbin = $curbin; newbin = $newbin\n");
    
   for ($i = 1; $i < $curbin; $i++ ) {
     $ictot = $ictot + $ic[$i];
     $cgptot = $cgptot + $cgp[$i];
     $cgntot = $cgntot + $cgn[$i];
     $triestot = $triestot + $tries[$i];
     $icdistot = $icdistot + $icdis[$i];
     $cgpchargetot = $cgpchargetot + $cgpcharge[$i];
     $cgnchargetot = $cgnchargetot + $cgncharge[$i];
    }
    
#    print OUT ( "Totals: $ictot, $cgptot, $cgntot, $triestot, ",
#           "$icdistot ", $icdistot/Max($ictot,1), " $cgpchargetot $cgnchargetot \n" ); 
#    print OUT ( "bin ICs CGPs CGNs Tries ICDIS  chg/ic   CGPchg   CGNchg   Netchg ",
#                "   NONIn    NONIp    INDn     INDp   wmax  Emax\n" );
#   for ($i = 1; $i < $curbin; $i++ ) {
#    printf OUT ( "%03d %3d %3d %3d %3d %8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f %  8.3f %5.2f %e %e\n",
#      $bintime[$i],  $ic[$i],  $cgp[$i],  $cgn[$i],  
#      $tries[$i], $icdis[$i], $icdis[$i]/Max(1,$ic[$i]), $cgpcharge[$i],
#      $cgncharge[$i], $netchg[$i],
#      $nonindn[$i], $nonindp[$i], $indn[$i], 
#      $indp[$i], $wmax[$i], $efield[$i] ,$efield2[$i] );

##    print OUT ( $i+$binshift+1,"  $ic[$i]  $cgp[$i]  $cgn[$i] ", 
##      "$tries[$i] $icdis[$i] $iccharge[$i] $cgcharge[$i]\n" );
#    }

    sub Max {
        my ($max,$tem);
        $max = shift(@_);
        foreach $tem (@_) {
            $max = $tem if $max < $tem;
        }
        return $max;
    }
