#!/usr/bin/perl
# a perl program to read in a set of outputs from bin.ion.pl and print out the totals line
# Assumes a set of binion.txt files with naming convention $prefix + CCN + $postfix
# Arguments are prefix, postfix, and list of ccn vals. For example:
#
# biniontot.pl jun29z blah.binion.txt 50 100 250 500 1000 1500 > output.txt

 open(OUT, '>-');
 $numarg = $#ARGV;
# print("numarg = $numarg\n");
 $prefix = $ARGV[0];
 $postfix = $ARGV[1];
 @ccn = @ARGV[2 .. $#ARGV];
# print ( @ccn );
# print ( "prefix = $prefix\n");
#print OUT ("CCN,crfrzftot, ciacrftot, chcnshtot, chcnihtot, qhacwrshtot, qracwtot, qrcnwtot, vhacwtot, ptemtot, pcondtot, chmul1tot,  csplintertot,  qrfrzftot,  qiacrftot,  crcnwtot, crmltshdtot, qcondtot, pevaptot, pmlttot, pdeptot, psubtot, pfrztot, crfrzstot, ciacrstot, qhmlrtot, qhlmlrtot, qhlacwrshtot, qhlacwtot, qhlacrtot, qhacwtot, qhacrtot, qhacwrshtot+qhlacwrshtot, qhmlrtot+qhlmlrtot\n");
print OUT ("CCN,ICtot,+CG,-CG,Tries,-NIC,+NIC,-IND,+IND,IC discharge (C), Charge per IC,Net IC charge, +CG charge (C), -CG charge (C), Pos. Channel pnts, Neg. Channel pnts, Total Channel pnts,Tot. Graupel Mass (kg s),Tot. Graupel vol. (kg^3 s),grdens,rainfalltot,hailfalltot,wmaxmax,grmstotti, nonindtotp-nonindtotn,Total Updraft KE (J s), Total Prec. Grav. PE (J s),rainfallautotot,rainfallshedtot,rainfallmelttot,Tot. Hail/FD Mass (kg s),Tot. Graupel+Hail Mass (kg s)\n");

#foreach $file (@ARGV) {
foreach $val (@ccn) {
    $file = "$prefix$val$postfix";
    #print OUT ("$file");
    print OUT ("$val,");
    open INPUT, "<$file";
    while ( <INPUT> ) {
    if ( /^Totals:/ ) {
      @parts = split(/Totals:/);
#      print OUT ("$_");
      print OUT ("$parts[1]");
      goto "foo";
      }
      
      
    } # while

    foo: 
    
    } # foreach

