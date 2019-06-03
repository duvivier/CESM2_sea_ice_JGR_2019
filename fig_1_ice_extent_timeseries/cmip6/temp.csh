#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('siconc')
set hemis = ('NH' 'SH')

##############
# start loops
set h = 1   
while ($h <= 2)  # max: 2
set v = 1
while ($v <= 1)  # max: 1

set varcode = $vars[$v]
set hemi = $hemis[$h]

##############
# Input into ncl
##############

        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./hist_timeseries_ice_ext.ncl

	    
        # Archive figures nicely
        rm *.ps	

@ v ++
end
@ h ++
end


