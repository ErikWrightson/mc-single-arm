# This shell script will take the name of an input file in the infiles directory and run it through
# the mc_single_arm simulator and then convert the compressed file to a root file.
if ( $#argv != 1 ) then
    echo "This script requires arguments of the name of a input file found in the infiles directory."
else
    set pth = "infiles/$argv[1].inp"
    if (-e "$pth") then
        cd src
        echo Here1
        mc_single_arm << endofinput >! ../runout/$argv[1].out
$argv[1]
endofinput
        echo Here2
        cd ..
        source /apps/root/6.10.02/bin/thisroot.csh
        echo Here3
        h2root worksim/$argv[1].rzdat worksim/$argv[1].root
    else
        echo "File name provided ($pth) does not exist or is not found in the worksim directory."
    endif
endif