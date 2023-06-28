# This shell script will take the name of an input file in the infiles directory and run it through the mc_single_arm simulator
# the , producing some histograms and getting the rate.
if ( $#argv != 4 ) then
    echo "This script requires arguments of the name of an input file found in the infiles directory, a current value in uA,"
    echo "the foil seperation in cm, and whether there are two or three foils in that order (true = three)."
else
    set pth = "infiles/$argv[1].inp"
    if (-e "$pth") then
        echo "$pth does exist"
        buildRoot.sh $argv[1]
        wait
        #Change the input file name to all lowercase since the rootfile will always be all lowercase.
        set rootName = `echo $argv[1] | tr '[:upper:]' '[:lower:]'`
        wait
        getRatesMulti.sh $rootName $argv[2] $argv[3] $argv[4]
        wait
        getRatesMSH.sh $rootName $argv[2] $argv[3] $argv[4]
    else
        echo "File name provided ($pth) does not exist or is not found in the worksim directory."
    endif
endif