#!/usr/bin/perl
# a perl program to read in a flash summary file and report the
#   ICs and CGs in time bins of whatever width. (set to 1 minute bins)
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
    $intervalminv = 1; # I do not recall why this was done -- looks like if set to a value of 2, then it sort of does 30s output, but not really?
    $intervalm = 1/$intervalminv;   # time discretization (usually 1 minute)
    $intervals = $intervalm*60;  
    print ("intervalminv,intervalm,intervals = $intervalminv, $intervalm, $intervals\n");
    $iontop = 0;
    $ionbottom = 0;
    $ionsides = 0;
    $ionchgcum2 = 0;
    $ionchgcum3 = 0;
    $ionerrtot = 0;
    $iondrifterr = 0;
    $first = 0;
    $binshift = 0;
    $tries = 0;
    $bintry = 0;
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
    $nonindgip = 0;
    $nonindgin = 0;
    $nonindhip = 0;
    $nonindhin = 0;
    $chgfall = 0;
    $scrchg = 0;
    $ionchg = 0;
    $ionchgcum = 0;
    $chgfallcum = 0;
    $scrchgcum = 0;
    $ind = 0;
    $ghsn = 0;
    $ghsp = 0;
    $ghin = 0;
    $ghip = 0;
    $gin = 0;
    $gip = 0;
    $gsn = 0;
    $gsp = 0;
    $hsn = 0;
    $hsp = 0;
    $grdn = 0;
    $hldn = 0;
    $postlightchg = 0;
    $netchg = 0;
    $netchgp = 0;
    $netchgn = 0;
    $advchgcum = 0;
    $advchgdiff = 0;
    $sedmichgdiff = 0;
    $sedmichgcum = 0;
    $sedchgdiff = 0;
    $sedchgcum = 0;
    $rnms = 0;
    $liqms = 0;
    $cwms = 0;
    $swms = 0;
    $rnmsa = 0;
    $rnmsb = 0;
    $ion = 1;
    $irestart = 0;
    $irestart2 = 0;
    $dv  = 1.0;
    $com2trmm = 0;
    $negchan = 0;
    $poschan = 0;

    $negchantot = 0;
    $poschantot = 0;
    $times = 0;
    $grmstot = 0;
    $udketottot = 0;
    $grpotetottot = 0;
    $grmstotti = 0; # time integrated
    $grvoltot = 0;
    $hlmstot = 0;
    $hlvoltot = 0;
    $wmaxmax = 0;
    $ssmxmax = 0;
    
# foreach ...
foreach $file (@ARGV) {
    open INPUT, "<$file";
    $start = 0;
    while ( <INPUT> ) {
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
      $irestart2 = 1;
     }
    if ( /^ COMMAS:  START =/ && $irestart == 1 ) {
      $irestart = 0;
      chomp();
      @parts = split(/ +/);
      $times1 = $parts[4];
     @timeparts = split(/\./, $times);
     if ( $timeparts[1] == 1 ) { $times = $times - 0.001 }
      # check if the restart is a continuation or goes back in time
      if ( $times1 < $times ) {
      print ("restart: $times, $times1\n");
#      print ("restart: $times, $times1, $timeparts[0], $timeparts[1]\n");
#        $min = int($times1/60.0);
        $min = int($times1/$intervals);
        $shift = int( ($times - $times1)/60.0 ) + 1;
#        $curbin = $min - $shift;
#        $curbin = $min+1;
        $curbin = $min - $binshift;
        $newbin = $curbin;
        $nel = $#nonindn;
#        print ("restart: min,curbin,shift,newbin = $min, $curbin, $shift, $newbin\n");
        print ("restart: min,curbin,binshift = $min, $curbin, $binshift\n");
        print ("current number of times = $nel\n");

# NOTE: Need to reset other counters here to zero and cumulative values to previous bin's value, as well, to be accurate!!!

        @nonindn[$curbin..$nel] = 0;
        @nonindp[$curbin..$nel] = 0;
        @wmax[$curbin..$nel] = 0;
        @ssmx[$curbin..$nel] = 0;
        @ic[$curbin..$nel] = 0;
        @tries[$curbin..$nel] = 0;
        @icdis[$curbin..$nel] = 0;
        @poschan[$curbin..$nel] = 0;
        @negchan[$curbin..$nel] = 0;
        @bintime[$curbin..$nel] = 0;
        $poschan = 0;
        $negchan = 0;
        $ic = 0;
        $icdis = 0;
        $nonindn = 0;
        $nonindp = 0;
        $tries = 0;
        $times = $times1;
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
    if ( /^END OF STEP: nstep, time/ ) {
      @parts = split(/ +/);
      $nsteptst = $parts[6];
#      print ("end of step: $nsteptst, $nstop, $min, $curbin, $newbin, $binshift\n");
      if ( $nsteptst == $nstop && $curbin == 1 ) {
#       $min = $min+1;
#       $newbin = $min - $binshift;
#       print ("HERE I AM: $newbin, $curbin/n");
       }
     }
    if ( /^NSTEP,/ ) {
     $next = <INPUT>;
     chomp($next);
     @parts = split(/ +/, $next);
     $nstep = $parts[1];
     $nstop = $parts[3];
     $times = $parts[4];
     @timeparts = split(/\./, $times);
     if ( $timeparts[1] == 1 ) { $times = $times - 0.001 }
#     $min = int(($times-$dt)/60.0);
     $min = int(($times-$dt)/$intervals);
     if ($binshift == 0) {
      $binshift = $min - 1;
      $curbin = $min - $binshift;
#      @min = ($min);
     }
     $newbin = $min - $binshift;
 #     if ( $times >= 10790 ) {
 #     print OUT ("$times $min $curbin $newbin $cgp $ic \n");
 #     }
     }
    if (  /^ Integration Done/ ) {
         $min = int(($times)/$intervals); 
         $tmp = (($times)/$intervals);
  # second part of 'if' test added to account for "off" restart times that
  # can result from EnKF timing (oldobsmethod)
      if ( $newbin == $curbin  &&  int(($times)/$intervals) == (($times)/$intervals)) { 
         $newbin = $curbin + 1; 
         $min = int(($times)/$intervals); 
         $tmp = (($times)/$intervals);
         }
     }
    if (  $newbin > $curbin  ) {

#     if ( $irestart2 == 1 ) {
#        print ("newbin,curbin = $newbin, $curbin\n");
        
#       }
       
     $bintime[$curbin] = $min*$intervalm;

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

     $negchantot = $negchantot + $negchan; # *$dv;
     $negchan[$curbin] = $negchan; # *$dv;
     $negchan = 0;

     $poschantot = $poschantot + $poschan; # *$dv;
     $poschan[$curbin] = $poschan; # *$dv;
     $poschan = 0;

     $icchargetot = $icchargetot + $iccharge; #*$dv;
     $icchargecum[$curbin] = $icchargetot;
     $iccharge[$curbin] = $iccharge; #*$dv;
     $iccharge = 0;
     
     $cgpchargetot = $cgpchargetot + $cgpcharge; # *$dv;
     $cgpchargecum[$curbin] = $cgpchargetot;
     $cgpcharge[$curbin] = $cgpcharge; # *$dv;
     $cgpcharge = 0;

     $cgnchargetot = $cgnchargetot + $cgncharge; # *$dv;
     $cgnchargecum[$curbin] = $cgnchargetot;
     $cgncharge[$curbin] = $cgncharge; # *$dv;
     $cgncharge = 0;
     
     $nonindtotn = $nonindtotn + $nonindn*$dv*$dt;
     $nonindn[$curbin] = $nonindn*$dv*$dt;
     $nonindn = 0;

     $nonindtotp = $nonindtotp + $nonindp*$dv*$dt;
     $nonindp[$curbin] = $nonindp*$dv*$dt;
     $nonindp = 0;
     
     $ghsn[$curbin] = $ghsn*$dv*$dt;
     $ghsn = 0;

     $ghsp[$curbin] = $ghsp*$dv*$dt;
     $ghsp = 0;

     $ghin[$curbin] = $ghin*$dv*$dt;
     $ghin = 0;

     $ghip[$curbin] = $ghip*$dv*$dt;
     $ghip = 0;

     $gin[$curbin] = $gin*$dv*$dt;
     $gin = 0;

     $gip[$curbin] = $gip*$dv*$dt;
     $gip = 0;

     $gsn[$curbin] = $gsn*$dv*$dt;
     $gsn = 0;

     $gsp[$curbin] = $gsp*$dv*$dt;
     $gsp = 0;

     $hin[$curbin] = $hin*$dv*$dt;
     $hin = 0;

     $hip[$curbin] = $hip*$dv*$dt;
     $hip = 0;

     $hsn[$curbin] = $hsn*$dv*$dt;
     $hsn = 0;

     $hsp[$curbin] = $hsp*$dv*$dt;
     $hsp = 0;

     $indtotn = $indtotn + $indn*$dv*$dt;
     $indn[$curbin] = $indn*$dv*$dt;
     $indn = 0;

     $indtotp = $indtotp + $indp*$dv*$dt;
     $indp[$curbin] = $indp*$dv*$dt;
     $indp = 0;
     
     $wmax[$curbin] = $wmax;
     $ssmx[$curbin] = $ssmx;
#     if ( $intervalm == 1 ) {
#     $wmax[$curbin] = $wmax;
#     } else {
#     $wmax[$curbin] = $wmaxtot/$intervals;
#     $wmaxtot = 0;
#     }


     $efield[$curbin] = $efield;

     $netchg[$curbin] = $netchg;
     $netchgp[$curbin] = $netchgp;
     $netchgn[$curbin] = $netchgn;
     $grvol[$curbin] = $grvol;
     $hlvol[$curbin] = $hlvol;
      $grms[$curbin] = $grms;
      $hlms[$curbin] = $hlms;
      $ticms[$curbin] = $ticms;
      $rnms[$curbin] = $rnms;
      $rnmsa[$curbin] = $rnmsa;
      $rnmsb[$curbin] = $rnmsb;

      $qcms[$curbin] =  $qcms;
      $qcmsa[$curbin] = $qcmsa;
      $qcmsb[$curbin] = $qcmsb;

      $liqms[$curbin] = $liqms;

      $fdmasstot[$curbin] = $fdmasstot;
      $fdmassa[$curbin] = $fdmassa;
      $fdmassb[$curbin] = $fdmassb;

     $hlnum[$curbin]       = $hlnum      ;
     $hlfdnum[$curbin]     = $hlfdnum    ;
     $hlnumfall[$curbin]   = $hlnumfall  ;
     $hlfdnumfall[$curbin] = $hlfdnumfall;
     $hlnumfall = 0;
     $hlfdnumfall = 0;
     
     $hlfdmasstot[$curbin] = $hlfdmasstot;
     
     $rainfalltot[$curbin] = $rainfalltot;
     $hailfalltot[$curbin] = $hailfalltot;
     $hailonlyfalltot[$curbin] = $hailonlyfalltot;
     $hailfdfalltot[$curbin] = $hailfdfalltot;

     $rainfallautotot[$curbin] = $rainfallautotot;
     $rainfallshedtot[$curbin] = $rainfallshedtot;
     $rainfallmelttot[$curbin] = $rainfallmelttot;

     $grvola[$curbin] = $grvola;
     $hlvola[$curbin] = $hlvola;
      $grmsa[$curbin] = $grmsa;
      $hlmsa[$curbin] = $hlmsa;

     $grvolb[$curbin] = $grvolb;
     $hlvolb[$curbin] = $hlvolb;
      $grmsb[$curbin] = $grmsb;
      $hlmsb[$curbin] = $hlmsb;
      
      $grdn[$curbin] = $grdn;
      $hldn[$curbin] = $hldn;

      $cwms[$curbin] = $cwms;
      $swms[$curbin] = $swms;

      $chgfall[$curbin] = $chgfall;
      $chgfallcum[$curbin] = $chgfallcum;
      $chgfall = 0;
      
      $scrchg[$curbin] = $scrchg;
      $scrchgcum[$curbin] = $scrchgcum;
      $scrchg = 0;
      
      $ionbottom[$curbin] = $ionbottom;
      $ionsides[$curbin] = $ionsides;
      $iontop[$curbin] = $iontop;
      
      $ionchg[$curbin] = $ionchg;
      $ionchgcum[$curbin] = $ionchgcum;
      $ionchg = 0;

      $ionchg2[$curbin] = $ionchg2;
      $ionchgcum2[$curbin] = $ionchgcum2;
      $ionchg2 = 0;

      $ionchg3[$curbin] = $ionchg3;
      $ionchgcum3[$curbin] = $ionchgcum3;
      $ionchg3 = 0;

      $advchgcum[$curbin] = $advchgcum;
      
      $sedmichgcum[$curbin] = $sedmichgcum;
      $sedchgcum[$curbin] = $sedchgcum;

      $udmf00[$curbin] = $udmf00;
      $udmf10[$curbin] = $udmf10;
      $udmf20[$curbin] = $udmf20;
      $udmf30[$curbin] = $udmf30;

      $uicmf10[$curbin] = $uicmf10;
      $uicmf20[$curbin] = $uicmf20;
      $uicmf30[$curbin] = $uicmf30;

      $udv5n[$curbin] = $udv5n;
      $udv05[$curbin] = $udv05;
      $udv10[$curbin] = $udv10;
      $udv20[$curbin] = $udv20;

      $udketot[$curbin] = $udketot;
      $ddketot[$curbin] = $ddketot;

      $grpotctot[$curbin] = $grpotctot;
      $grpotptot[$curbin] = $grpotptot;
      $grpottot[$curbin] = $grpotctot+$grpotptot;

#     print ("$min, $curbin, $times \n");
     $curbin = $newbin;
      }
     
    if ( /^wmax/ ) {
      chomp();
      @parts = split( / +/ );
      $wmax = $parts[2];
      $wmaxtot = $wmaxtot + $wmax*$dt;
      $wmaxmax = Max($wmax, $wmaxmax);
      }

    if ( /^ ssmax/ ) {
      chomp();
      @parts = split( / +/ );
      $ssmx = $parts[3];
      $ssmxtot = $ssmxtot + $ssmx*$dt;
      $ssmxmax = Max($ssmx, $ssmxmax);
      }

#net1 pos/neg/net charge (C):   3.85422E+02  -3.99048E+02  -1.36255E+01     
    if ( /^net1 pos\/neg/ ) {
      chomp();
      @parts = split( / +/ );
      $netchgp1 = $parts[4];
      $netchgn1 = $parts[5];
      $netchg1 = $parts[6];
      if ( $first == 0 ) {
       $first = 1;
       $basechg = $netchg1;
       print OUT ("basechg = $basechg\n");
       }
      }
# postlight pos/neg/net charge (C):   1.29770E+03  -1.36479E+03  -6.70915E+01
    if ( /^postlight/ ) {
      chomp();
      @parts = split( / +/ );
      $netchgp = $parts[4];
      $netchgn = $parts[5];
      $netchg = $parts[6];
      
      $ionchg2 = $ionchg2 + $netchg - $netchg1;
      $ionchgcum2 = $ionchgcum2 + $netchg - $netchg1;
      $ionerrtot = $ionerrtot + $netchg - $netchg1 - $scrchg0;
      }

#premic pos/neg/net charge (C):   1.00463E+03  -1.28107E+03  -2.76433E+02
    if ( /^premic/ ) {
      chomp();
      @parts = split( / +/ );
      $netchg2 = $parts[6];
      if ( $netchg != 0 ) {
        $advchgdiff = $netchg2 - $netchg;
        $advchgcum = $advchgcum + $advchgdiff;
        }
      if ( $com2trmm == 1 ) {
        $netchgp = $parts[4];
        $netchgn = $parts[5];
        $netchg = $parts[6];
       }
      }

#chgnet1,chg2,delta-chg =   -18.2435141292038       -18.3262037913795     
# chgnet1,chg2,delta-chg =   -15.4431444197265       -15.5181581480950     

    if ( /^ chgnet1/ ) {
      chomp();
      @parts = split( / +/ );
      $netchg2a = $parts[3];
      $netchg2b = $parts[4];
      if ( $netchg2a != 0 ) {
        $sedchgdiff = $netchg2a - $netchg2;
        $sedchgcum = $sedchgcum - $sedchgdiff;
#        $advchgcum = $advchgcum + $advchgdiff;
        }
      if ( $com2trmm == 1 ) {
#        $netchgp = $parts[4];
#        $netchgn = $parts[5];
#        $netchg = $parts[6];
       }
      }

    if ( /^postmic/ ) {
      chomp();
      @parts = split( / +/ );
      $netchg3 = $parts[6];
      if ( $netchg3 != 0 ) {
        $sedmichgdiff = $netchg3 - $netchg2a;
        $sedmichgcum = $sedmichgcum - $sedmichgdiff;
#        $advchgcum = $advchgcum + $advchgdiff;
        }
      if ( $com2trmm == 1 ) {
#        $netchgp = $parts[4];
#        $netchgn = $parts[5];
#        $netchg = $parts[6];
       }
      }

#frozen drop masses tot/above/below 0C  2.26704E+09  2.15390E+09  1.13139E+08
    if ( /^frozen drop masses tot/ ) {
      chomp();
      @parts = split( / +/ );
      #$grvola = $parts[4];
      $fdmasstot = $parts[5];
      $fdmassa = $parts[6];
      $fdmassb = $parts[7];
      }
#Hail/HailFD number, fallout  3.57260E+10  2.66921E+10  1.51377E+08  1.07314E+08
    if ( /^Hail\/HailFD number, fallout/ ) {
      chomp();
      @parts = split( / +/ );
      $hlnum = $parts[3];
      $hlfdnum = $parts[4];
      $hlnumfall = $hlnumfall + $parts[5];
      $hlfdnumfall = $hlfdnumfall + $parts[6];
      }
#HailFD mass   9.46950E+07  1.24366E+08
    if ( /^HailFD mass/ ) {
      chomp();
      @parts = split( / +/ );
      $hlfdmasstot = $parts[2];
      }
      
#graupel, hail volumes  7.17800E+03  2.46250E+03
    if ( /^graupel, hail volumes/ ) {
      chomp();
      @parts = split( / +/ );
      $grvol = $parts[3];
      $hlvol = $parts[4];
      }
#graupel, hail, ice masses  1.07262E+10  2.59410E+09  3.41498E+09
    if ( /^graupel, hail, ice masses/ ) {
      chomp();
      @parts = split( / +/ );
      $grms = $parts[4];
      $hlms = $parts[5];
      $ticms = $parts[6];
      }
#graupel/hail volumes above/below 0C  6.91350E+03  2.28350E+03  2.36500E+02  1.71000E+02
#graupel/hail masses above/below 0C  1.03317E+10  2.33682E+09  3.64165E+08
    if ( /^graupel\/hail volumes above/ ) {
      chomp();
      @parts = split( / +/ );
      $grvola = $parts[4];
      $hlvola = $parts[5];
      $grvolb = $parts[6];
      $hlvolb = $parts[7];
      }
    if ( /^graupel\/hail masses above/ ) {
      chomp();
      @parts = split( / +/ );
      $grmsa = $parts[4];
      $hlmsa = $parts[5];
      $grmsb = $parts[6];
      $hlmsb = $parts[7];
      $grmstot = $grmstot + $dt*$parts[4];
      $grmstotti = $grmstotti + $dt*$parts[4]/60;
      $hlmstot = $hlmstot + $dt*$parts[5];
      }
    if ( /^graupel\/hail particle volumes/ ) {
      chomp();
      @parts = split( / +/ );
      $grvoltot = $grvoltot + $dt*$parts[5];
      $hlvoltot = $hlvoltot + $dt*$parts[6];
      $grdn = $parts[7];
      $hldn = $parts[8];
      }
#rain masses tot,above/below 0C  2.78572E+09  4.16217E+08  2.36951E+09
    if ( /^rain masses tot/ ) {
      chomp();
      @parts = split( / +/ );
      $rnms =  $parts[4];
      $rnmsa = $parts[5];
      $rnmsb = $parts[6];
      }

#total liquid on ice   1.78681E+08
    if ( /^total liquid on ice/ ) {
      chomp();
      @parts = split( / +/ );
      $liqms =  $parts[4];
      }

    if ( /^ Total rainfall/ ) { 
      chomp();
      @parts = split( / +/ );
      $rainfalltot = $parts[4];
    }

    if ( /^ Total AUTO rainfall/ ) { 
      chomp();
      @parts = split( / +/ );
      $rainfallautotot = $parts[5];
    }

    if ( /^ Total SHED rainfall/ ) { 
      chomp();
      @parts = split( / +/ );
      $rainfallshedtot = $parts[5];
    }

    if ( /^ Total MELT rainfall/ ) { 
      chomp();
      @parts = split( / +/ );
      $rainfallmelttot = $parts[5];
    }

    if ( /^ Total hailfall/ ) { 
      chomp();
      @parts = split( / +/ );
      $hailfalltot = $parts[4];
    }

    if ( /^ Total hail =/ ) { 
      chomp();
      @parts = split( / +/ );
      $hailonlyfalltot = $parts[4];
    }

    if ( /^ Total hail from FD/ ) { 
      chomp();
      @parts = split( / +/ );
      $hailfdfalltot = $parts[6];
    }

#udmf at 0,-10,-20,and -30:  1.78594E+09  1.96553E+09  1.89928E+09  1.75535E+09
#uicmf at -10, -20, -30:  4.75490E+03  2.62660E+04  4.00441E+05
#udv  -5, 5, 10, 20:  1.72591E+10  1.30640E+11  4.75146E+10  6.53299E+09
    if ( /^udmf at 0/ ) {
      chomp();
      @parts = split( / +/ );
      $udmf00 = $parts[4];
      $udmf10 = $parts[5];
      $udmf20 = $parts[6];
      $udmf30 = $parts[7];
      }
    if ( /^uicmf at -10/ ) {
      chomp();
      @parts = split( / +/ );
      $uicmf10 = $parts[5];
      $uicmf20 = $parts[6];
      $uicmf30 = $parts[7];
      }

#udketot, ddketot :  1.12744E+12  2.46643E+11
#grpotctot, grpotptot :  4.17616E+12  1.16881E+13
    
    if ( /^udketot, ddketot/ ) {
      chomp();
      @parts = split( / +/ );
      $udketot = $parts[3];
      $ddketot = $parts[4];
      $udketottot = $udketottot + $dt*$parts[3];;
     
     }

    if ( /^grpotctot/ ) {
      chomp();
      @parts = split( / +/ );
      $grpotctot = $parts[3];
      $grpotptot = $parts[4];
      $grpotetottot = $grpotetottot + $dt*$parts[4];
     
     }

#qc mass above/below 0C  3.69897E+07  3.71538E+07
    if ( /^qc mass above/ ) {
      chomp();
      @parts = split( / +/ );
      $qcms =  $parts[4]+$parts[5];
      $qcmsa = $parts[4];
      $qcmsb = $parts[5];
      }
    
    if ( /^udv  -5/ ) {
      chomp();
      @parts = split( / +/ );
      $udv5n = $parts[5];
      $udv05 = $parts[6];
      $udv10 = $parts[7];
      $udv20 = $parts[8];
      }

#qctot, qvtot  0.00000E+00  7.19420E+10
    if ( /^qctot, qvtot/ ) {
      chomp();
      @parts = split( / +/ );
      $cwms = $parts[2];
      $vpms = $parts[3];
      }

#snow mass  0.00000E+00
    if ( /^snow mass/ ) {
      chomp();
      @parts = split( / +/ );
      $swms = $parts[2];
      }

#Charge fallout this step (pos,neg,tot):    5.59332E-01  -1.27936E+00  -7.20028E-01
    if ( /^Charge fallout this step/ ) {
      chomp();
      @parts = split( / +/ );
      $chgfall = $chgfall + $parts[7];
      $chgfallcum = $chgfallcum + $parts[7];
      }
#Screening layer charge (pos,neg,tot) 7.67661E+00  -7.62283E+00   5.37712E-02
    if ( /^Screening layer charge/ ) {
#     print ;
      $ion = 0;
      chomp();
      @parts = split( / +/ );
      $scrchg = $scrchg + $parts[6];
      $scrchgcum = $scrchgcum + $parts[6];
#     print ("\n $nstep scrchg = $scrchg \n");      
      }
#difference ion pos/neg/net charge (C):  -9.81570E-04   3.23683E-03   2.24585E-03
    if ( /^difference ion pos/ ) {
      chomp();
      @parts = split( / +/ );
      $ionchg = $ionchg + $parts[7];
      $ionchgcum = $ionchgcum + $parts[7];
      }

#difference2 ion pos/neg/net charge (C):  -9.81570E-04   3.23683E-03   2.24585E-03
    if ( /^difference2 ion pos/ ) {
      chomp();
      @parts = split( / +/ );
      $ionchg3 = $ionchg3 + $parts[7];
      $ionchgcum3 = $ionchgcum3 + $parts[7];
      
      }


#Net Ion charge fluxes (domain,top,bottom):   1.45763E-02 -1.52142E-02  2.97905E-02  2.44317E-02 -1.12446E-02
# 0   1    2       3            4                 5           6             7             8           9
    if ( /^Net Ion charge fluxes/ ) {
#     print ;
      $ion = 1;
      chomp();
      @parts = split( / +/ );
      $ionsides = $ionsides + $parts[7] - $parts[8] - $parts[9];
      $ionbottom = $ionbottom + $parts[9];
      $iontop = $iontop + $parts[8];
      
   #    $scrchg = $scrchg + $ionsides + $ionbottom + $iontop;
   #    $scrchgcum = $scrchgcum + $ionsides + $ionbottom + $iontop;
   #    $scrchg0 = $ionsides + $ionbottom + $iontop;
      $scrchg = $scrchg + $parts[7];
      $scrchgcum = $scrchgcum + $parts[7];
      $scrchg0 = $parts[7];
      }

    if ( /^ EFIELD--E-MAX/ ) {
      chomp();
      @parts = split( / +/ );
      $efield = $parts[$#parts];
      }
   if ( $nstep <= $nstop )  {
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

      $next = <INPUT>;  # graupel-snow alone
      chomp($next);
      @parts = split(/ +|,/,$next);
      $gsn = $gsn + $parts[3];
      $gsp = $gsp + $parts[5];

      $next = <INPUT>;  # graupel-ice alone
      chomp($next);
      @parts = split(/ +|,/,$next);
      $nonindgin = $nonindgin + $parts[3];
      $gin = $gin + $parts[3];
      $nonindgip = $nonindgip + $parts[5];
      $gip = $gip + $parts[5];

      $next = <INPUT>;  # hail-snow alone
      chomp($next);
      @parts = split(/ +|,/,$next);
      $hsn = $hsn + $parts[3];
      $hsp = $hsp + $parts[5];

      $next = <INPUT>;  # hail-ice alone
      chomp($next);
      @parts = split(/ +|,/,$next);
      $nonindhin = $nonindhin + $parts[3];
      $hin = $hin + $parts[3];
      $nonindhip = $nonindhip + $parts[5];
      $hip = $hip + $parts[5];

      }
    if ( /^ COM2TRMM: IC flashes,/ ) {
     $com2trmm = 1;
     $ion = 1;
     @parts = split(/ +/ );
     $ic  = $parts[7];
     $cgn = $parts[8];
     $cgp = $parts[9];
     }
    if ( /^Try number/ ) {
     $tries = $tries + 1;
     $bintry = $bintry + 1;
     }
    if ( /^numchan/ ) {
      chomp();
      @parts = split(/ +/);
      $negchan = $negchan + $parts[1];
      $poschan = $poschan + $parts[2];
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
     if ( $pos < 0.0001 && $neg < 0.0001 ) {
        $pos = $pos*$dv;
        $neg = $neg*$dv;
        }
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
     } # if ( $nstep < $nstop )  
    } # while
  foo: 
   } # foreach
   
    print OUT ( "Removed $icx ICs, $cgpx +CGs, $cgnx -CGs\n" );
#    print OUT ( "Totals: $ictot, $cgptot, $cgntot, $triestot, $nonindtotn, $nonindtotp, ",
#           "$icdistot, $icchargetot, $cgpchargetot, $cgnchargetot \n" ); 
#    print OUT ("Totals: ");

    $ictot = 0; 
    $cgptot = 0;  
    $cgntot = 0;
    $triestot = 0;
    $icdistot = 0;
    $icchargetot = 0;
    $cgpchargetot = 0;
    $cgnchargetot = 0;
    $nonindtotn = 0;
    $nonindtotp = 0;
    $indtotn = 0;
    $indtotp = 0;
    
   for ($i = 1; $i < $curbin; $i++ ) {
     $ictot = $ictot + $ic[$i];
     $cgptot = $cgptot + $cgp[$i];
     $cgntot = $cgntot + $cgn[$i];
     $triestot = $triestot + $tries[$i];
     $icdistot = $icdistot + $icdis[$i];
     $icchargetot = $icchargetot + $iccharge[$i];
     $cgpchargetot = $cgpchargetot + $cgpcharge[$i];
     $cgnchargetot = $cgnchargetot + $cgncharge[$i];
     $nonindtotn = $nonindtotn + $nonindn[$i];
     $nonindtotp = $nonindtotp + $nonindp[$i];
     $indtotn = $indtotn + $indn[$i];
     $indtotp = $indtotp + $indp[$i];
    }


    $grdens = 0;
    if ( $grvoltot > 0 ) {
     $grdens = $grmstot/$grvoltot;
     }
    $hailfalltot[0] = $hailfalltot[1];
    $hailfdfalltot[0] = $hailfdfalltot[1];
    $hailonlyfalltot[0] = $hailonlyfalltot[1];
    
    printf OUT ( "Totals:  %3d, %3d, %3d, %3d,  %10.1f, % 10.1f,%10.1f, % 10.1f, %10.1f, %7.1f, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e,%10.3e, %10.3e, %10.3e, %10.3e, %10.3e,%10.3e,%10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e\n",
                $ictot, $cgptot, $cgntot, $triestot, $nonindtotn, $nonindtotp, $indtotn, $indtotp,
                $icdistot,  $icdistot/Max($ictot,1), $icchargetot, $cgpchargetot, $cgnchargetot, $poschantot, 
                $negchantot, $poschantot + $negchantot,$grmstot,$grvoltot,$grdens,$rainfalltot,$hailfalltot,
                $wmaxmax,$grmstotti, $nonindtotp-$nonindtotn, $udketottot, $grpotetottot,$rainfallautotot,$rainfallshedtot,$rainfallmelttot,$hlmstot,$grmstot+$hlmstot ); 

    if ( $ion == 0 ) {
    print OUT ( "bin ICs CGPs CGNs Tries ICDIS CGPchg CGNchg CGPchgcum CGNchgcum ",
                "NONIn NONIp INDn INDp wmax Emax netcharge netchgp netchgn ",
                "graupelvolume  graupelmass  grvola grvolb ",
                "icemass rainmass chargefallout screenchg advchg udmf00 udmf-10 ",
                "udmf-20 udmf-30 uicmf-10 uicmf-20 uicmf-30 icnetchg ",
                "udv5n udv05 udv10 udv20 \n" );
    } elsif ($ion == 1 ) {
    print OUT ( "Time (min),IC Rate,IC/min,+CG Rate,-CG Rate,Tries,ICDIS,+CG charge,-CG charge,+CG cum. charge,-CG cum. charge,",
                "NONI. neg,NONI. pos,IND. neg,IND. pos.,W-max (m/s),E-max (kV/m),Net Charge,Net Pos. Chg.,Net Neg. Chg.,",
                "Graupel Volume (km**3),Graupel Mass (kg),Graupel Mass T < 0,Graupel Mass T > 0,Graupel Vol. T < 0,Graupel Vol. T > 0,Hail Mass (kg),Hail Vol. T < 0,Hail Vol. T > 0,Graupel+Hail Mass T < 0,",
                "Cloud Ice Mass,Cloud Droplet Mass,Rain Mass Tot,Rain Mass (T < 0),Rain Mass (T > 0),Snow Mass,G-I Neg,G-I Pos,G-S Neg,G-S Pos,H-I Neg,H-I Pos,H-S Neg,H-S Pos,Updraft Mass Flux (T=0) (kg/s),Updraft Mass Flux (T=-10) (kg/s),",
                "Updraft Mass Flux (T=-20) (kg/s),Updraft Mass Flux (T=-30) (kg/s),Crystal Mass Flux (T=-10),Crystal Mass Flux (T=-20),Crystal Mass Flux (T=-30),icnetchg,",
                "Downdraft Volume (< -5m/s),Updraft Volume (> 5m/s),Updraft Volume (> 10m/s),Updraft Volume (> 20m/s),",
                "Charge Fallout (C),Advection Charge Loss (C),Ion Drift Charge (C),Other Ion (C),Unaccounted Charge (C),",
                "Ion Flux: Top,Ion Flux: Bottom,Ion Flux: Sides,Positive Sources,Negative Sources,Total Sources,Mean Graupel Density,",
                "Mean Hail Density,Total Liquid Fraction,Sed/Mic charge (C),Sed Charge2 (C), SSmx (%),UD KE (J),DD KE (J),Grav. Pot (cld) (J),",
                "Grav. Pot. (prec) (J),Grav. Pot (J),Cloud Droplet Mass (T<0),FD Mass Tot.,HailFD Mass Tot.,Hail Num. Tot,HailFD Num. Tot.,",
                "Hail Num Fallout,HailFD Num Fallout,Hail+Graupel Accum.,Hail Accum.,HailFD Accum.,Rain Accum.,Rain-Auto Accum.,Rain-Shed Accum.,Rain-Melt Accum.,",
                "1-min Hail Accum.,1-min HailFD Accum.,Mean Mass Sfc. Hail Diam. (mm)\n" );
     }
   for ($i = 1; $i < $curbin; $i++ ) {
#    printf OUT ( "%03d, %3d, %7.2f, %3d, %3d, %3d, %8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, %5.2f, %5.2f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %9d, %9d, %9d, %e, %e, %e, %e, %e\n",
    if ( $intervalminv == 1 ) {
    
      $meanhailmass = ($hailonlyfalltot[$i]-$hailonlyfalltot[$i-1])/($hlnumfall[$i] + 1.e-8);
      $meanhaildiam = 1000.*($meanhailmass*6.0/(3.14159*900.0) )**(1.0/3.0);
      
    
#    printf OUT ( "%8.1f, %3d, %7.2f, %3d, %3d, %3d, %8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, %5.2f, %5.2f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %9d, %9d, %9d, %e, %e, %e, %e, %e\n",
    printf OUT ( "%8.1f,  %3d,   %7.2f,   %3d,   %3d,   %3d,   %8.3f,   % 8.3f,  % 8.3f,  % 8.3f,  % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, %5.2f,  %5.2f,  %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %9d, %9d, %9d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n",
      $bintime[$i],  $ic[$i], $icpermin[$i], $cgp[$i],  $cgn[$i],  
      $tries[$i], $icdis[$i],  $cgpcharge[$i],
      $cgncharge[$i], $cgpchargecum[$i], $cgnchargecum[$i], 
      $nonindn[$i], $nonindp[$i], $indn[$i], 
      $indp[$i], $wmax[$i], $efield[$i]*0.001,
      $netchg[$i],$netchgp[$i],$netchgn[$i],$grvol[$i], $grms[$i],$grmsa[$i],$grmsb[$i], $grvola[$i], $grvolb[$i], $hlms[$i], $hlvola[$i], $hlvolb[$i], $hlmsa[$i]+$grmsa[$i],
      $ticms[$i],$cwms[$i], $rnms[$i],$rnmsa[$i],$rnmsb[$i], $swms[$i], $gin[$i], $gip[$i], $gsn[$i], $gsp[$i], $hin[$i], $hip[$i], $hsn[$i], $hsp[$i],
      $udmf00[$i], $udmf10[$i] , $udmf20[$i],$udmf30[$i],
      $uicmf10[$i], $uicmf20[$i], $uicmf30[$i], $icchargecum[$i],
      $udv5n[$i], $udv05[$i], $udv10[$i], $udv20[$i],-$chgfallcum[$i], $advchgcum[$i],$ionchgcum3[$i],$ionchgcum2[$i]-$ionchgcum3[$i],
       $cgpchargecum[$i] + $cgnchargecum[$i] - $chgfallcum[$i] + $ionchgcum2[$i] + $advchgcum[$i] - $netchg[$i] + $basechg,
       $iontop[$i],$ionbottom[$i],$ionsides[$i], $poschan[$i], $negchan[$i], $poschan[$i] + $negchan[$i], $grdn[$i], $hldn[$i], $liqms[$i], $sedmichgcum[$i], $sedchgcum[$i],
       $ssmx[$i],$udketot[$i], $ddketot[$i], $grpotctot[$i], $grpotptot[$i], $grpottot[$i],$qcmsa[$i],$fdmasstot[$i],$hlfdmasstot[$i],
       $hlnum[$i],$hlfdnum[$i],$hlnumfall[$i],$hlfdnumfall[$i],$hailfalltot[$i],$hailonlyfalltot[$i],$hailfdfalltot[$i],$rainfalltot[$i],
       $rainfallautotot[$i],$rainfallshedtot[$i],$rainfallmelttot[$i],
       $hailonlyfalltot[$i]-$hailonlyfalltot[$i-1],$hailfdfalltot[$i]-$hailfdfalltot[$i-1]),$meanhaildiam;

#      $scrchgcum[$i], $ionchgcum[$i], $iontop[$i],$ionbottom[$i],$ionsides[$i] );

#      $cgpchargecum[$i] + $cgnchargecum[$i] - $chgfallcum[$i] + $scrchgcum[$i] + $advchgcum[$i] - $netchg[$i] + $basechg, 
#      $cgpchargecum[$i] + $cgnchargecum[$i] - $chgfallcum[$i] + $ionchgcum3[$i] + $advchgcum[$i] - $netchg[$i] + $basechg,
      
      # ,$iontop[$i]+$ionbottom[$i]+$ionsides[$i]
#    print OUT ( $i+$binshift+1,"  $ic[$i]  $cgp[$i]  $cgn[$i] ", 
#      "$tries[$i] $icdis[$i] $iccharge[$i] $cgcharge[$i]\n" );
      } else {

     if ( $bintime[$i] - int($bintime[$i]) == 0  ) {

#    printf OUT ( "%8.1f, %3d, %7.2f, %3d, %3d, %3d, %8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, %5.2f, %5.2f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %9d, %9d, %9d, %e, %e, %e, %e, %e\n",
    printf,OUT ( "%8.1f,  %3d,   %7.2f,   %3d,  %3d,  %3d,  %8.3f,  % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, %5.2f,  %5.2f,  %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e,%e,  %e, %9d, %9d, %9d, %e, %e, %e, %e, %e\n",
      $bintime[$i],  $ic[$i]+$ic[$i-1], ($icpermin[$i]+$icpermin[$i-1])*$intervalm, $cgp[$i]+$cgp[$i-1],  $cgn[$i]+ $cgn[$i-1],  
      $tries[$i]+$tries[$i-1], $icdis[$i]+$icdis[$i-1],  $cgpcharge[$i]+$cgpcharge[$i-1],
      $cgncharge[$i]+$cgncharge[$i-1], $cgpchargecum[$i], $cgnchargecum[$i], 
      $nonindn[$i]+$nonindn[$i-1], $nonindp[$i]+$nonindp[$i-1], $indn[$i]+$indn[$i-1], 
      $indp[$i]+$indp[$i-1], $wmax[$i], $efield[$i]*0.001,
      $netchg[$i],$netchgp[$i],$netchgn[$i],$grvol[$i], $grms[$i],$grmsa[$i],$grmsb[$i], $grvola[$i], $grvolb[$i], $hlms[$i], $hlvola[$i], $hlvolb[$i],
      $ticms[$i],$cwms[$i], $rnms[$i],$rnmsa[$i],$rnmsb[$i], $swms[$i], $gin[$i], $gip[$i], $gsn[$i], $gsp[$i], $hin[$i], $hip[$i], $hsn[$i], $hsp[$i],
      $udmf00[$i], $udmf10[$i] , $udmf20[$i],$udmf30[$i],
      $uicmf10[$i], $uicmf20[$i], $uicmf30[$i], $icchargecum[$i],
      $udv5n[$i], $udv05[$i], $udv10[$i], $udv20[$i],-$chgfallcum[$i], $advchgcum[$i],$ionchgcum3[$i],$ionchgcum2[$i]-$ionchgcum3[$i],
       $cgpchargecum[$i] + $cgnchargecum[$i] - $chgfallcum[$i] + $ionchgcum2[$i] + $advchgcum[$i] - $netchg[$i] + $basechg,
       $iontop[$i],$ionbottom[$i],$ionsides[$i], $poschan[$i]+$poschan[$i-1], $negchan[$i]+$negchan[$i-1], $poschan[$i] + $negchan[$i] +$poschan[$i-1] + $negchan[$i-1], 
       $grdn[$i], $hldn[$i], $liqms[$i], $sedmichgcum[$i], $sedchgcum[$i]);



      }
      
     
     }

    }


    sub Max {
        my ($max,$tem);
        $max = shift(@_);
        foreach $tem (@_) {
            $max = $tem if $max < $tem;
        }
        return $max;
    }
