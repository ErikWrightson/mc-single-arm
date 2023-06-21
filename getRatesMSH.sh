# This shell script will take the name of a root file in the worksim directory and run it through
# the hms_foil_rates_carbon_inelastic_MSH macro, producing some histograms and getting the rate of events originating from each.
if ( $#argv != 4 ) then
    echo "This script requires arguments of the name of a root file found in the worksim directory, a current value in uA,"
    echo "the foil seperation in cm, and whether there are two or three foils in that order."
else
    set pth = "worksim/$argv[1].root"
    if (-e "$pth") then
        echo "$pth does exist"
        cd pbmodel/
        make
        wait
        cd ..
        clear
        root $pth << EOF
        .L pbmodel/libF1F209.so
        .x examples/hms_foil_rates_carbon_inelastic_MSH.C("$argv[1]", $argv[2], $argv[3], $argv[4])
        .q
        EOF
    else
        echo "File name provided ($pth) does not exist or is not found in the worksim directory."
    endif
endif