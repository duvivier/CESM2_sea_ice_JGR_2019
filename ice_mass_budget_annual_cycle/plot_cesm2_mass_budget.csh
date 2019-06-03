#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('sidmassgrowthwat' 'sidmassgrowthbot' 'sidmasssi' 'sidmassmelttop' 'sidmassmeltbot' 'sidmasslat')
#set vars = ('sidmassgrowthwat' 'sidmassgrowthbot' 'sidmasssi' 'sidmassmelttop' 'sidmassmeltbot' 'sidmasslat' 'sidmassevapsubl')
set hemis = ('NH' 'SH')

##############
# start loops
set v = 1   
while ($v <= 7)  # max: 7
set h = 1
while ($h <= 1)  # max: 2

set varcode = $vars[$v]
set hemi = $hemis[$h]

##############
# Input into ncl
##############

        echo 'Plotting  mass budget term '$varcode 'for hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_clim_mass_budget.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ h ++
end
@ v ++
end


