#! /bin/tcsh -f

# Script to cycle through regridding dates
# Alice DuVivier- Sept 2016
#################################################
##############
# USER INPUTS
##############
set vars1 = ('clivi' 'clt')
set vars2 = ('swcrf' 'lwcrf')
set hemis = ('NH' 'SH')

##############
# start loops
set h = 1   
while ($h <= 1)  # max: 2
set hemi = $hemis[$h]

##############
# varcode 2 loop - crf
##############
set v2 = 1
while ($v2 <= 2)  # max: 7

set varcode = $vars2[$v2]

        echo 'Plotting crf term '$varcode 'for hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_spatial_atm_crf.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ v2 ++
end

##############
# LWP doesn't require a varcode
##############
   echo 'Plotting LWP for hemisphere '$hemi
   ncl 'hemi          = "'$hemi'"' \
       ./pi_spatial_atm_lwp.ncl

   ncl 'hemi          = "'$hemi'"' \
       ./pi_spatial_atm_lwp_pcnt_diff.ncl

   # Archive figures nicely
    rm *.ps	

##############
# varcode 1 loop - cloud vars
##############
set v1 = 1
while ($v1 <= 2)  # max: 7

set varcode = $vars1[$v1]

        echo 'Plotting cloud term '$varcode 'for hemisphere '$hemi
        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_spatial_atm.ncl

        ncl 'varcode       = "'$varcode'"'\
            'hemi          = "'$hemi'"' \
            ./pi_spatial_atm_pcnt_diff.ncl
	    
        # Archive figures nicely
        rm *.ps	

@ v1 ++
end


@ h ++
end


