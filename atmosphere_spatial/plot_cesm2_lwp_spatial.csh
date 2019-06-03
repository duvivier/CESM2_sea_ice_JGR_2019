#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('cllvi')
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

        echo 'Plotting cloud term '$varcode 'for hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_spatial_lwp.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ v ++
end
@ h ++
end


