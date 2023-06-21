# This shell script will take the name of a root file in the worksim directory and run it through
# the hms_foil_rates_carbon_inelastic_sieveholes macro, producing some histograms and getting the rate.
if ( $#argv != 2 ) then
    echo "This script requires arguments of the name of a root file found in the worksim directory and a current value in uA in that order."
else
    set pth = "../worksim/$argv[1].root"
    if (-e "$pth") then
        echo "$pth does exist"
        cd ../pbmodel/
        make
        wait
        cd ..
        clear
        root << EOF
        .L pbmodel/libF1F209.so
        .x examples/hms_foil_rates_carbon_inelastic_sieveholes.C("$argv[1]", $argv[2])
        .q
        EOF
    else
        echo "File name provided ($pth) does not exist or is not found in the worksim directory."
    endif
endif