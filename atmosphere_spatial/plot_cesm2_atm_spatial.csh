#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('clwvi' 'clivi' 'clt' 'rlds' 'rsds' 'sfcWind' 'psl')
set hemis = ('NH' 'SH')

##############
# start loops
set h = 1   
while ($h <= 2)  # max: 2
set v = 1
while ($v <= 7)  # max: 7

set varcode = $vars[$v]
set hemi = $hemis[$h]

##############
# Input into ncl
##############

        echo 'Plotting cloud term '$varcode 'for hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_spatial_atm_vars.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ v ++
end
@ h ++
end


