This directory contains example run scripts for the hail embryo type paper. It generates a modified Weisman-Klemp sounding with a given surface pressure and temperature. The hodograph is generated as a quarter circle followed by piecewise continuous linear shear. Control parameters are given in Table 2 of Mansell and Straka and in the run scripts in this directory. The scripts here are for the "base" simulations. Sensitivities are 

Options used for the manuscript (Mansell and Straka 2027)

Cloud condensation nuclei concentration is set by the 'ccn' namelist token. 1.3e9 for base (1300CN), and 0.5e9 for the 500CN variant.

Option that controls rain production by shedding
ished2cld :  1: Send shed liquid (from wet growth) to cloud droplets
             2: Reduce collection to offset shedding ("no-shed")
             3: As for 2 but only for T < 0 (turns off wet growth [cold] shedding) ("no-shed cold")
             4: As for 2 but only for T > 0 (turns off warm shedding (graupel/hail collecting cloud droplets)

Frozen drops option:
 frozendrops : 0 = off; 1 = turn on separate FD category (and fraction in hail), used for all but the 'fd2' experiments
               2 = FD are a number fraction within graupel ("fd2") and in hail

Rain sources tracking:
  iraintypes : 0 = off; 1 = on (all simulations used 1)

