#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set hemis = ('NH' 'SH')

##############
# start loops
set h = 1
while ($h <= 1)  # max: 2
set hemi = $hemis[$h]

##############
# Input into ncl
##############

        echo 'Plotting net ebudget for hemisphere '$hemi
        ncl 'hemi          = "'$hemi'"' \
            ./pi_clim_ebudget_net.ncl

        echo 'Plotting net ebudget for hemisphere '$hemi
        ncl 'hemi          = "'$hemi'"' \
            ./pi_clim_ebudget_netturb.ncl

        echo 'Plotting albedo for hemisphere '$hemi
        ncl 'hemi          = "'$hemi'"' \
            ./pi_clim_albedo.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ h ++
end


