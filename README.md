# CESM2_sea_ice_JGR_2020

CREATOR:

Alice DuVivier - duvivier@ucar.edu
August 2020

DESCRIPTION:

Scripts for figures and analysis prepared for manuscript
published by JGR (2020). 

DuVivier, A. K., Holland, M. M., Kay, J. E., Tilmes, S., Gettelman, A., & Bailey, D. A. (2020). Arctic and Antarctic Sea Ice Mean State in the Community Earth System Model Version 2 and the Influence of Atmospheric Chemistry. Journal of Geophysical Research: Oceans, 125(8). https://doi.org/10.1029/2019JC015934


Analysis scripts are in:
(1) NCL with cshell wrapper scripts and cheyenne dav 
submission scripts and (2) python notebooks.

Directories are separated by figure number.
Output is files in png and ps format. Final figures for the
paper were put together in paint to get appropriate paneling.

To run ncl on the command line, do
    >ncl SCRIPTNAME.ncl

To run csh on the command line, do
    >./plot_fig_1.csh
