#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars = ('sithick')
set hemis = ('NH' 'SH')

##############
# start loops
set h = 2   
while ($h <= 2)  # max: 2
set v = 1
while ($v <= 1)  # max: 2

set varcode = $vars[$v]
set hemi = $hemis[$h]

##############
# Input into ncl
##############

        echo 'Plotting '$varcode 'for pdf hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_pdf_ice_thick.ncl

        #ncl 'varcode       = "'$varcode'"'\
        #    'hemi          = "'$hemi'"' \
        #    ./hist_pdf_ice_thick.ncl

        #ncl 'varcode       = "'$varcode'"'\
        #    'hemi          = "'$hemi'"' \
        #    ./both_hist_pi_pdf_ice_thick.ncl

        #ncl 'varcode       = "'$varcode'"'\
        #    'hemi          = "'$hemi'"' \
        #    ./hist_pdf_ice_thick_w_icesat.ncl

        # Archive figures nicely
        rm *.ps	

@ v ++
end
@ h ++
end


