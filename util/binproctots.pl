#!/usr/bin/perl
# a perl program to read in a set of outputs from binproc35.pl and print out the totals line
# Assumes a set of binproc.txt files with naming convention $prefix + CCN + $postfix
# Arguments are prefix, postfix, and list of ccn vals. For example:
#
# binproctots.pl jun29z blah.binproc.txt 50 100 250 500 1000 1500

 open(OUT, '>-');
 $numarg = $#ARGV;
# print("numarg = $numarg\n");
 $prefix = $ARGV[0];
 $postfix = $ARGV[1];
 @ccn = @ARGV[2 .. $#ARGV];
# print ( @ccn );
# print ( "prefix = $prefix\n");
print OUT ("CCN,crfrzftot, ciacrftot, chcnshtot, chcnihtot, qhacwrshtot, qracwtot, qrcnwtot, vhacwtot, ptemtot, pcondtot, chmul1tot,  csplintertot,  qrfrzftot,  qiacrftot,  crcnwtot, crmltshdtot, qcondtot, pevaptot, pmlttot, pdeptot, psubtot, pfrztot, crfrzstot, ciacrstot, qhmlrtot, qhlmlrtot, qhlacwrshtot, qhlacwtot, qhlacrtot, qhacwtot, qhacrtot, qhacwrshtot+qhlacwrshtot, qhmlrtot+qhlmlrtot\n");
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

