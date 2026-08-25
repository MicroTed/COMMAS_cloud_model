#!/usr/bin/perl
# a perl program to read in a .out file and get the domain total microphysical
#  processes (numproct=31) in time bins of whatever width. (set to 1 minute bins)
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
    $crfrzf = 0;
    $intervalm = 1;   # time discretization (usually 1 minute)
    $intervals = $intervalm*60;  
    $fdauto = 0;
    $fdshed = 0;
    $fdmelt = 0;
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
    $rnms = 0;
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
    $iqrauto = 0;
    $iqrshed = 0;
    $iqrmelt = 0;
    $nrain = 0;
    $nfdauto = 0;
    $nfdshed = 0;
    $nfdmelt = 0;
    
    
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
      # check if the restart is a continuation or goes back in time
      if ( $times1 < $times ) {
      print ("restart: $times, $times1\n");
        $min = int($times1/60.0);
#        $min = int(($times1-$dt)/$intervals);
        $shift = int( ($times - $times1)/60.0 ) + 1;
#        $curbin = $min - $shift;
        $curbin = $min+1;
        $newbin = $curbin;
        $nel = $#nonindn;
        print ("restart: min,curbin,shift,newbin = $min, $curbin, $shift, $newbin\n");
        print ("current number of times = $nel\n");

# NOTE: Need to reset other counters here to zero and cumulative values to previous bin's value, as well, to be accurate!!!

        @nonindn[$curbin..$nel] = 0;
        @nonindp[$curbin..$nel] = 0;
        @wmax[$curbin..$nel] = 0;
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
        $vhacw = 0;
        $qrcnw = 0;
        $qrcnw = 0;
        $pcond = 0;
        $qcond = 0;
        $pevap = 0;
        $pmlt = 0;
        $pdep = 0;
        $psub = 0;
        $pevapr = 0;
        $qrcnw = 0;
        $qracw = 0;
        $qhlmlr = 0;
        $qhmlr = 0;
        $qhlacr = 0;
        $qhlacw = 0;
        $ciint = 0;
        $ic = 0;
        $cgp = 0;
        $fdauto = 0;
        $fdshed = 0;
        $fdmelt = 0;
        $crfrzf = 0;
        $ciacrf = 0;
        $chcnsh = 0;
        $chcnih = 0;
        $qhacwrsh = 0;
        $qhlacwrsh = 0;
        $qhacw = 0;
        $qhacr = 0;
        $qhshr = 0;
        $qfshr = 0;
        $qhlshr = 0;
        $qhlcnh = 0;
        $qhlcnf = 0;
        $pchld = 0;
        $chlfmlr = 0;
        $nfdauto = 0;
        $nfdshed = 0;
        $nfdmelt = 0;
        $chlcnh = 0;
        $chlcnf = 0;
        $cwctfz = 0;
        $pmlt = 0;
        $pdep = 0;
        $psub = 0;
        $pevapr = 0;
        $pfrz = 0;
        $crfrzs = 0;
        $ciacrs = 0;
        $chmul1 = 0;
        $csplinter = 0;
        $qrfrzf = 0;
        $qiacrf = 0;
        $crcnw = 0;
        $crmltshd = 0;
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
    if ( /^ QRAUTO/ && $iqrauto == 0 ) {
         $iqrauto = 1;
         $nrain = $nrain + 1;
       }
    if ( /^ QRSHED/ && $iqrshed == 0 ) {
         $iqrshed = 1;
         $nrain = $nrain + 1;
       }
    if ( /^ QRMELT/ && $iqrmelt == 0 ) {
         $iqrmelt = 1;
         $nrain = $nrain + 1;
       }
    if ( /^NSTEP,/ ) {
     $next = <INPUT>;
     chomp($next);
     @parts = split(/ +/, $next);
     $nstep = $parts[1];
     $nstop = $parts[3];
     $times = $parts[4];
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
      if ( $newbin == $curbin ) { 
         $newbin = $curbin + 1; 
         $min = int(($times)/$intervals); 
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
     
     $fdautotot = $fdautotot + $fdauto;
     $fdauto[$curbin] = $fdauto;
     $fdauto = 0;
     
     $fdshedtot = $fdshedtot + $fdshed;
     $fdshed[$curbin] = $fdshed;
     $fdshed = 0;
     
     $fdmelttot = $fdmelttot + $fdmelt;
     $fdmelt[$curbin] = $fdmelt;
     $fdmelt = 0;
     
     $crfrzftot = $crfrzftot + $crfrzf;
     $crfrzf[$curbin] = $crfrzf;
     $crfrzf = 0;

     $ciacrftot = $ciacrftot + $ciacrf;
     $ciacrf[$curbin] = $ciacrf;
     $ciacrf = 0;

     $chcnshtot = $chcnshtot + $chcnsh;
     $chcnsh[$curbin] = $chcnsh;
     $chcnsh = 0;

     $chcnihtot = $chcnihtot + $chcnih;
     $chcnih[$curbin] = $chcnih;
     $chcnih = 0;

     $qhacwrshtot = $qhacwrshtot + $qhacwrsh;
     $qhacwrsh[$curbin] = $qhacwrsh;
     $qhacwrsh = 0;

     $qlhacwrshtot = $qhlacwrshtot + $qhlacwrsh;
     $qhlacwrsh[$curbin] = $qhlacwrsh;
     $qhlacwrsh = 0;

     $qhacwtot = $qhacwtot + $qhacw;
     $qhacw[$curbin] = $qhacw;
     $qhacw = 0;

     $qhacrtot = $qhacrtot + $qhacr;
     $qhacr[$curbin] = $qhacr;
     $qhacr = 0;

     $qhshrtot = $qhshrtot + $qhshr;
     $qhshr[$curbin] = $qhshr;
     $qhshr = 0;

     $qfshrtot = $qfshrtot + $qfshr;
     $qfshr[$curbin] = $qfshr;
     $qfshr = 0;

     $qhlshrtot = $qhlshrtot + $qhlshr;
     $qhlshr[$curbin] = $qhlshr;
     $qhlshr = 0;

     $qhlcnhtot = $qhlcnhtot + $qhlcnh;
     $qhlcnh[$curbin] = $qhlcnh;
     $qhlcnh = 0;

     $qhlcnftot = $qhlcnftot + $qhlcnf;
     $qhlcnf[$curbin] = $qhlcnf;
     $qhlcnf = 0;

     $pchldtot = $pchldtot + $pchld;
     $pchld[$curbin] = $pchld;
     $pchld = 0;

     $chlfmlrtot = $chlfmlrtot + $chlfmlr;
     $chlfmlr[$curbin] = $chlfmlr;
     $chlfmlr = 0;

     $nfdautotot = $nfdautotot + $nfdauto;
     $nfdauto[$curbin] = $nfdauto;
     $nfdauto = 0;

     $nfdshedtot = $nfdshedtot + $nfdshed;
     $nfdshed[$curbin] = $nfdshed;
     $nfdshed = 0;

     $nfdmelttot = $nfdmelttot + $nfdmelt;
     $nfdmelt[$curbin] = $nfdmelt;
     $nfdmelt = 0;

     $chlcnhtot = $chlcnhtot + $chlcnh;
     $chlcnh[$curbin] = $chlcnh;
     $chlcnh = 0;

     $chlcnftot = $chlcnftot + $chlcnf;
     $chlcnf[$curbin] = $chlcnf;
     $chlcnf = 0;
     
     $cwctfztot = $cwctfztot + $cwctfz;
     $cwctfz[$curbin] = $cwctfz;
     $cwctfz = 0;

     $ciinttot = $ciinttot + $ciint;
     $ciint[$curbin] = $ciint;
     $ciint = 0;

     $qhlacwtot = $qhlacwtot + $qhlacw;
     $qhlacw[$curbin] = $qhlacw;
     $qhlacw = 0;

     $qhlacrtot = $qhlacrtot + $qhlacr;
     $qhlacr[$curbin] = $qhlacr;
     $qhlacr = 0;


     $qhmlrtot = $qhmlrtot + $qhmlr;
     $qhmlr[$curbin] = $qhmlr;
     $qhmlr = 0;

     $qhlmlrtot = $qhlmlrtot + $qhlmlr;
     $qhlmlr[$curbin] = $qhlmlr;
     $qhlmlr = 0;

     $qracwtot = $qracwtot + $qracw;
     $qracw[$curbin] = $qracw;
     $qracw = 0;

     $qrcnwtot = $qrcnwtot + $qrcnw;
     $qrcnw[$curbin] = $qrcnw;
     $qrcnw = 0;

     $vhacwtot = $vhacwtot + $vhacw;
     $vhacw[$curbin] = $vhacw;
     $vhacw = 0;
     
     if ( $vhacw[$curbin] > 0 ) {
       $rimdens[$curbin] = $qhacwrsh[$curbin]/$vhacw[$curbin];
      } else {
       $rimdens[$curbin] = 0;
      }

     $ptemtot = $ptemtot + $ptem;
     $ptem[$curbin] = $ptem;
     $qrcnw = 0;

     $pcondtot = $pcondtot + $pcond;
     $pcond[$curbin] = $pcond;
     $pcond = 0;

#      $qcond   = $qcond + $parts[17];  # condensation rate
     $qcondtot = $qcondtot + $qcond;
     $qcond[$curbin] = $qcond;
     $qcond = 0;
#      $pevap  = $pevap + $parts[18]; # rain evap
     $pevaptot = $pevaptot + $pevap;
     $pevap[$curbin] = $pevap;
     $pevap = 0;
#      $pmlt   = $pmlt + $parts[19]; # melting rate
     $pmlttot = $pmlttot + $pmlt;
     $pmlt[$curbin] = $pmlt;
     $pmlt = 0;
#      $pdep   = $pdep + $parts[20]; # deposition rate
     $pdeptot = $pdeptot + $pdep;
     $pdep[$curbin] = $pdep;
     $pdep = 0;
#      $psub   = $psub + $parts[21]; # sublimation rate
     $psubtot = $psubtot + $psub;
     $psub[$curbin] = $psub;
     $psub = 0;
     
     $pevaprtot = $pevaprtot + $pevapr;
     $pevapr[$curbin] = $pevapr;
     $pevapr = 0;
     
#      $pfrz   = $pfrz + $parts[22]; # melting rate
     $pfrztot = $pfrztot + $pfrz;
     $pfrz[$curbin] = $pfrz;
     $pfrz = 0;
#      $crfrzs   = $crfrzs + $parts[23]; # rain -> snow number rate (bigg)
     $crfrzstot = $crfrzstot + $crfrzs;
     $crfrzs[$curbin] = $crfrzs;
     $crfrzs = 0;
#      $ciacrs   = $ciacrs + $parts[24]; # rain -> snow number rate (ice capture)
     $ciacrstot = $ciacrstot + $ciacrs;
     $ciacrs[$curbin] = $ciacrs;
     $ciacrs = 0;


     $chmul1tot = $chmul1tot + $chmul1;
     $chmul1[$curbin] = $chmul1;
     $chmul1 = 0;

     $csplintertot = $csplintertot + $csplinter;
     $csplinter[$curbin] = $csplinter;
     $csplinter = 0;

     $qrfrzftot = $qrfrzftot + $qrfrzf;
     $qrfrzf[$curbin] = $qrfrzf;
     $qrfrzf = 0;

     $qiacrftot = $qiacrftot + $qiacrf;
     $qiacrf[$curbin] = $qiacrf;
     $qiacrf = 0;

     $crcnwtot = $crcnwtot + $crcnw;
     $crcnw[$curbin] = $crcnw;
     $crcnw = 0;

     $crmltshdtot = $crmltshdtot + $crmltshd;
     $crmltshd[$curbin] = $crmltshd;
     $crmltshd = 0;
     
     
#     $wmax[$curbin] = $wmax;
     if ( $intervalm == 1 ) {
     $wmax[$curbin] = $wmax;
     } else {
     $wmax[$curbin] = $wmaxtot/$intervals;
     $wmaxtot = 0;
     }


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

#     print ("$min, $curbin, $times \n");
     $curbin = $newbin;
      }
     
    if ( /^wmax/ ) {
      chomp();
      @parts = split( / +/ );
      $wmax = $parts[2];
      $wmaxtot = $wmaxtot + $wmax*$dt;
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
#       print OUT ("basechg = $basechg\n");
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
      }
    if ( /^graupel\/hail particle volumes/ ) {
      chomp();
      @parts = split( / +/ );
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
   if ( $nstep < $nstop )  {
    if ( /^FD source rates/ ) {
      chomp(); # mass rated
      @parts = split(/ +/);
      $fdauto = $fdauto + $parts[3];
      $fdshed = $fdshed + $parts[4];
      $fdmelt = $fdmelt + $parts[5];
       }
    if ( /^processes/ ) {
      chomp();
      @parts = split(/ +/);
      $crfrzf = $crfrzf + $parts[1];
      $ciacrf = $ciacrf + $parts[2];
      $chcnsh = $chcnsh + $parts[3];
      $chcnih = $chcnih + $parts[4];
      $qhacwrsh  = $qhacwrsh + $parts[5];
      $qracw  = $qracw + $parts[6];
      $qrcnw  = $qrcnw + $parts[7];
      $vhacw  = $vhacw + $parts[8]; # actually vhacw + vhacr + vhshdr
      $ptem   = $ptem + $parts[9];
      $chmul1 = $chmul1 + $parts[10];
      $csplinter = $csplinter + $parts[11];
      $qrfrzf = $qrfrzf + $parts[12];
      $qiacrf = $qiacrf + $parts[13];
      $crcnw  = $crcnw + $parts[14];
      $crmltshd = $crmltshd + $parts[15];
      $pcond   = $pcond + $parts[16]; # latent heating condensation + droplet evaporation


      $qcond   = $qcond + $parts[17];  # condensation rate
      $pevap  = $pevap + $parts[18]; # rain evap
      $pmlt   = $pmlt + $parts[19]; # melting rate
      $pdep   = $pdep + $parts[20]; # deposition rate
      $psub   = $psub + $parts[21]; # sublimation rate
      $pfrz   = $pfrz + $parts[22]; # freezing rate
      $crfrzs   = $crfrzs + $parts[23]; # rain -> snow number rate (bigg)
      $ciacrs   = $ciacrs + $parts[24]; # rain -> snow number rate (ice capture)
      
      $qhmlr   = $qhmlr + $parts[25]; # graupel melting rate
      $qhlmlr  = $qhlmlr + $parts[26]; # hail melting rate
      
      $qhlacwrsh  = $qhlacwrsh + $parts[27];
      $qhlacw  = $qhlacw + $parts[28];
      $qhlacr  = $qhlacr + $parts[29];

      $qhacw  = $qhacw + $parts[30];
      $qhacr  = $qhacr + $parts[31];
      $qhlcnh = $qhlcnh + $parts[32];
      $cwctfz = $cwctfz + $parts[33];
      $pevapr  = $pevapr + $parts[34]; # rain evaporation rate
      $ciint = $ciint + $parts[35];
      $chlcnh = $chlcnh + $parts[36];
      $chlcnf = $chlcnf + $parts[37];
      #$qhlcnh = $qhlcnh + $parts[38]; # repeat of 32
      $qhlcnf = $qhlcnf + $parts[39];
      $pchld = $pchld + $parts[40];
      $chlfmlr = $chlfmlr + $parts[41];
      
      $len = $#parts;
    #  print OUT ("len = $len\n");
      if ( $len >= 44 && $nrain > 0 ) {
       # number rates
        $nfdauto = $nfdauto + $parts[$len-2];
        $nfdshed = $nfdshed + $parts[$len-1];
        $nfdmelt = $nfdmelt + $parts[$len-0];
        
        }
      if ( $len - $nrain >= 44 ) {
        # file has shedding terms
          $qhshr = $qhshr + $parts[42];
          $qfshr = $qfshr + $parts[43];
          $qhlshr = $qhlshr + $parts[44];
        }
      
      
#      $qhacw  = $qhacw + $parts[5];
      }
    
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
   
#    print OUT ( "Removed $icx ICs, $cgpx +CGs, $cgnx -CGs\n" );
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



    printf OUT ( "Totals:  %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e\n",
                $crfrzftot, $ciacrftot, $chcnshtot, $chcnihtot, $qhacwrshtot, $qracwtot, $qrcnwtot, $vhacwtot, $ptemtot, $pcondtot, $chmul1tot,  $csplintertot,  $qrfrzftot, 
                $qiacrftot,  $crcnwtot, $crmltshdtot, $qcondtot, $pevaptot, $pmlttot, $pdeptot, $psubtot, $pfrztot, $crfrzstot, $ciacrstot, $qhmlrtot, $qhlmlrtot, $qhlacwrshtot,
                $qhlacwtot, $qhlacrtot, $qhacwtot, $qhacrtot, $qhacwrshtot+$qhlacwrshtot, $qhmlrtot + $qhlmlrtot, $pevaprtot, $qhlcnhtot, $qhlcnftot, $chlcnhtot, $chlcnftot, $pchldtot, $chlfmlrtot,
                $fdautotot, $fdshedtot, $fdmelttot, $qhshrtot, $qfshrtot, $qhlshrtot ); 
    if ( $ion == 0 ) {
    print OUT ( "bin ICs CGPs CGNs Tries ICDIS CGPchg CGNchg CGPchgcum CGNchgcum ",
                "NONIn NONIp INDn INDp wmax Emax netcharge netchgp netchgn ",
                "graupelvolume  graupelmass  grvola grvolb ",
                "icemass rainmass chargefallout screenchg advchg udmf00 udmf-10 ",
                "udmf-20 udmf-30 uicmf-10 uicmf-20 uicmf-30 icnetchg ",
                "udv5n udv05 udv10 udv20 \n" );
    } elsif ($ion == 1 ) {
    print OUT ( "Time (min),crfrzf,ciacrf,chcnsh,chcnih,qhacwrsh,qracw,qrcnw,vhacw,ptem,pcond,chmul1,csplinter,qrfrzf,",
                "qiacrf,crcnw,crmltshd,crfrz+ciarcr,qrfrzf+qiacrf,Rime Density,chcnsh+chcnih,chmul1+csplinter,qcond, pevap,",
                "pmlt, pdep, psub, pfrz, crfrzs, ciacrs, qhmlr, qhlmlr, qhlacwrsh, qhlacw, qhlacr, qhacw, qhacr, qhacwrsh+qhlacwrsh,",
                "qh+hlmlr,cwctfz,ciint,rainevap,qhlcnh,qhlcnf,chlcnh,chlcnf,pchld,chlfmlr,FD-auto,FD-shed,FD-melt,Shed-GR,Shed-FD,Shed-HL\n" );
     }
   for ($i = 1; $i < $curbin; $i++ ) {
    printf OUT ( "%03d, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n",
#    printf OUT ( "%03d, %3d, %7.2f, %3d, %3d, %3d, %8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, % 8.3f, %5.2f, %5.2f, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %9d, %9d, %9d, %e, %e\n",
      $bintime[$i],  $crfrzf[$i], $ciacrf[$i],  $chcnsh[$i], $chcnih[$i],  $qhacwrsh[$i],  $qracw[$i],  $qrcnw[$i], 
       $vhacw[$i],  $ptem[$i], $pcond[$i], $chmul1[$i],  $csplinter[$i],  $qrfrzf[$i],$qiacrf[$i],  $crcnw[$i],  $crmltshd[$i],
       $crfrzf[$i]+$ciacrf[$i],$qrfrzf[$i]+$qiacrf[$i], $rimdens[$i], $chcnsh[$i]+$chcnih[$i], $chmul1[$i]+$csplinter[$i],$qcond[$i], 
       $pevap[$i], $pmlt[$i], $pdep[$i], $psub[$i], $pfrz[$i], $crfrzs[$i], $ciacrs[$i], $qhmlr[$i], $qhlmlr[$i], $qhlacwrsh[$i], $qhlacw[$i],
       $qhlacr[$i], $qhacw[$i], $qhacr[$i], $qhacwrsh[$i]+$qhlacwrsh[$i],$qhmlr[$i]+$qhlmlr[$i],$cwctfz[$i],$ciint[$i],$pevapr[$i],$qhlcnh[$i],
       $qhlcnf[$i],$chlcnh[$i],$chlcnf[$i],$pchld[$i],$chlfmlr[$i],$fdauto[$i],$fdshed[$i],$fdmelt[$i],$qhshr[$i],$qfshr[$i],$qhlshr[$i]);

#      $scrchgcum[$i], $ionchgcum[$i], $iontop[$i],$ionbottom[$i],$ionsides[$i] );

#      $cgpchargecum[$i] + $cgnchargecum[$i] - $chgfallcum[$i] + $scrchgcum[$i] + $advchgcum[$i] - $netchg[$i] + $basechg, 
#      $cgpchargecum[$i] + $cgnchargecum[$i] - $chgfallcum[$i] + $ionchgcum3[$i] + $advchgcum[$i] - $netchg[$i] + $basechg,
      
      # ,$iontop[$i]+$ionbottom[$i]+$ionsides[$i]
#    print OUT ( $i+$binshift+1,"  $ic[$i]  $cgp[$i]  $cgn[$i] ", 
#      "$tries[$i] $icdis[$i] $iccharge[$i] $cgcharge[$i]\n" );
    }


    sub Max {
        my ($max,$tem);
        $max = shift(@_);
        foreach $tem (@_) {
            $max = $tem if $max < $tem;
        }
        return $max;
    }
