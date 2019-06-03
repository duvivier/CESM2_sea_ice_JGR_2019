#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('silwnet' 'siswnet')
set hemis = ('NH' 'SH')

##############
# start loops
set v = 1   
while ($v <= 2)  # max: 2
set h = 1
while ($h <= 1)  # max: 2

set varcode = $vars[$v]
set hemi = $hemis[$h]

##############
# Input into ncl
##############

        echo 'Plotting ebudget term '$varcode 'for hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_clim_ebudget_netrad.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ h ++
end
@ v ++
end

