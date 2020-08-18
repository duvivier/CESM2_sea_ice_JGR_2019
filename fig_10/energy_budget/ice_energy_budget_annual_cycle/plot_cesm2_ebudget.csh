#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set hemis = ('NH' 'SH')

##############
# start hemisphere loop
##############
set h = 1
while ($h <= 2)  # max: 2
set hemi = $hemis[$h]

# NCL calls that don't require a varcode
   echo 'Plotting ebudget diffs for hemisphere '$hemi
   ncl 'hemi          = "'$hemi'"' \
        ./pi_clim_ebudget.ncl

   # Archive figures nicely
   rm *.ps
@ h ++
end


