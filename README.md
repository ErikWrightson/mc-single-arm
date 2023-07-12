mc_single_arm
==============

Hall C single arm Monte Carlo for HMS and SHMS 



To compile type "make" in the src directory.
NOTE: Make sure paths to Cern libararies are defined.
On JLab machines, do "setup cernlib/2005".

Main code is mc_single_arm.f

If include files are changed it is best to do "make Clean" and then "make"
since the Makefile is not smart enough to look for dependency

Running code
------------

Easiest way is to use "buildRoot.sh <filename>" (filename assumed to be in infiles subdirectory - will stop if not)

NOTE: This will generate the compressed .rzdat file in the worksim directory and then automatically uncompress into a .root file leaving both in the worksim directory.

Alternative 1:
use "run_mc_single_arm <filename>" (filename assumed to be in infiles subdirectory)

NOTE: This leaves compressed the .rzdat file in the worksim directory that must be unpacked.

Alternative 2:
cd src
mc_single_arm 
(ask for input file name (assumed in infiles subdirectory))

NOTE: This leaves compressed the .rzdat file in the worksim directory that must be unpacked.

* Input file : infile_name.inp
* Output file is at outfiles/infile_name.out 
* The hbook file is at worksim/infile_name.rzdat
* Hard coded flags in mc_single_arm.f  :
    * hut_ntuple : Write out ntuple (1 for HMS and 1411 for SHMS) described below. (true)
    * spec_ntuple: Write out ntuple 1412 described below (false)
    * decay_flag :	Allow decay of particle(false)
    * use_front_sieve : Pass particles through front sieve (false)   
* Hard coded flags in shms/mc_shms.f :
    * use_sieve : pass particles through the sieve between HB and Q1 (false)
    * use_coll  : pass particles through the collimator between HB and Q1 (true)

Uncompressing a .rzdat file
---------
* source /apps/root/6.10.02/bin/thisroot.csh
* h2root worksim/filename.rzdat worksim/filename.root (assuming the .rzdat file exists)

Code flow
--------- 
* Read in the input file
* Starts event loop
* The seed for the random number generator is picked by the time. each run of the code will be different set of random numbers.
* The lab or beam coordinate system has +z along the target towards the beam dump, +y is vertical pointing down and +x is pointing beam left.
* Randomly selects beam-target interaction point.
* The coordinate system of the SHMS is +z_s along the central ray, +x_s is vertical pointing down and +y_s is horizontal point to large angle. 
* Randomly selects particle delta=(p-pcent)/p, dy_s/dz_s amd dx_s/dz_s.
* Drifts to z_s=0 plane
* Calculated multiple scattering the target using the input radiation length
    * If target length greater than 3cm assume LH2 cryotarget for multiple scattering calculation.
    *  If target length less than 3cm assume simple solid target for multiple scattering calculation.
* calulated mutliple scattering in scattering chamber window, air and spectrometer entrance window.
* Call the mc_shms/mc_hms subroutine to transport the particle through SHMS/HMS apertures and to the focal plane.
* At the focal plane drift the particle through  1st Cerenkov or vacuum pipe (for SHMS only, depending on flag in input file),the drift chamber, scintillators, 2nd cerenkov and calorimeter to check if it hits them . 
* As it passes through the drift chamber randomly change the position according to the presumed DC resolution and fit the positions and calculate the focal plane positions and angles.
* From focal plane calculate the reconstructed target quantities ytar,yptar,xptar and delta using the optics matrix. Note that for now use the known x target position in the optics matrix calculations. To mimic what would be in the analysis, need to calculate xtar from the reconstructed quantities using the raster vertical position and then recalculate the reconstructed quantities with updated xtar.
* If particle makes it to the focal plane and hits all detectors the event variables are put in ntuple.

Multifoil and Optics Sieve
---------------
* Currently the Multifoil and Optics Sieve optics options require some hard coded changes depending on what you would like to do. Hard coded changes can all be done within the mc_single_arm.f file.

NOTE: It is recommended that these be adapted to options allowed to be set in the input file (or having more multifoil options registered throughout the code for target thickness).

* Multifoil : Running with a multifoil target can currently be done for two optics options of target thickness in the input file.
    * Optics 1 - Three foil (foils at z=0, +/- 10 cm)
        * Set target thickness in input file to -3 to access this optics target.
        * This target is currently not in the in use target ladder (check current target).
        * Changing this setting to a different foil seperation just requires changing the 10 initial z calculation to the distance from 0 you want the foils in cm.
        * For reference, the default initial z calculation is "z = (grnd() - 0.5) * foil_tk + foil_nm * 10" and the 10 would be what to change.
    * Optics 2 - Two foil (foils set by default to z+/-3cm)
        * Set target thickness in input file to -2 to acces this optics target.
        * Check your current target to see what multifoil seperation you need.
        * Changing this setting to a different foil seperation change the 3 to whatever seperation from 0 you want and 6 to the total seperation between your two foils. (Works for symmetrically seperated foils).
        * For reference, the default initial z calculation is "z = (grnd() - 0.5) * foil_tk - 3+ foil_nm * 6" and the 3 and 6 would be what to change.

* Optics Sieve : Running with the optics sieve in place or not has two different options.
    * Front Sieve (SHMS) - to use the front sieve in the SHMS the the use_front_sieve logical must be switched from false to true or vice versa run with no sieve.
    * Sieve (HMS) - to use the sieve for the HMS the use_sieve logical must be changed from false to true or vice versa run with no sieve. (currently set to true by default).
    * This could be implemented into being an option in the input file.

Sub Directories
---------------
* src : location of source code for simulator
    * hms: HMS subroutines
    * shms  : SHMS subroutines
* infiles : input files
* outfiles : the output files 
* worksim : rzdat files and root files after uncompressing using scripts
* inelastic_carbon : output files from running HMS carbon foil rates macros
    * multifoil : multifoil HMS inelastic carbon rate pdf outputs
    * multifoilHoles : multifoil with sieve hole seperation HMS inelastic carbon rate pdf outputs
    * sieveData : single foil with sieve hole seperation HMS inelastic carbon rate pdf outputs
    * single : single foil HMS inelastic carbon rate pdf outputs
* runout : collected terminal outputs from running the simulator
* pbmodel : model used to link to the object file libF1F209.so for getting cross-sections for the HMS carbon foil rates macros
* examples : example analysis macros and HMS ideal (calculated using intial conditions, i.e. perfect sieve ID, perfect multifoil ID, etc.) Carbon foil rate macros
* ratesScripts : script files to automatically run any of the HMS inelastic carbon rate macros

Info on infiles
---------------
* Mostly self explanatory
* dp/p down and up should be at least -15.0 and 25.0 
* theta down and up is dy/dz , horizontal angle relative to central ray keep at least -55,55mr
* phi down and up is the dx/dz, vertical angle relative to central ray ,keep at least -50, 50mr
* Remember to keep the Dp/p,theta,phi and ztgt reconstruction cut larger than thrown. 
* The target length just need to set flag if aerogel will be in the detector stack
* Up to the experiment to decide if the 1st Cerenkov detector will be needed for the experiemnt
* If 1st Cerenkov not used then can replace with vacuum pipe. Option of helium bag is availble. this was for study.  

Info on Shell (.sh) scripts
---------------
* buildRoot.sh : 
    * input - 
        * filename (filename.inp file must be in the infiles directory)
    * description -
        * 
    * script flow - 
        * Runs mc-single-arm for for the given .inp file
        * Ensures that the filename given is cast to all lowercase since the .root file will always be all lowercase.
        * Unpacks the filename.rzdat file into a filename.root file in the worksim directory

NOTE: ALL SCRIPTS BELOW FOUND IN ratesScripts DIRECTORY

* buildAndRunSingle :
    * input - 
        * filename - filename.inp file must be in the infiles directory
        * current - current to get the inelastic carbon rates at in uA
    * description -
        * Runs the simulator and gets the pdfs and histograms for a single carbon foil target with and without sieve seperation cuts. 
    * script flow - 
        * Checks that the right amount of arguments was entered and the path to the input file exists
        * Calls buildRoot.sh for filename
        * Ensures that the filename given is cast to all lowercase since the .root file will always be all lowercase.
        * Runs the getRates.sh script (found in ratesScripts)for the filename.root file to get the HMS carbon rates pdfs and histograms for a single carbon foil without sieve seperation cuts.
        * Runs the getRatesSieve.sh script (found in ratesScripts) for the filename.root file to get the HMS carbon rates pdfs and histograms for a single carbon foil WITH sieve seperation cuts.

* buildAndRunMulti.sh :
    * input - 
        * filename - filename.inp file must be in the infiles directory
        * current - current to get the inelastic carbon rates at in uA
        * foilSep - serpation of carbon foils in cm
        * threeFoil - boolean determining if this was a 3 foil run (foil assumed to be at z=0).
    * description -
        * Runs the simulator and gets the pdfs and histograms for a carbon multifoil target with and without sieve seperation cuts.
    * script flow - 
        * Checks that the right amount of arguments was entered and the path to the input file exists
        * Calls buildRoot.sh for filename
        * Ensures that the filename given is cast to all lowercase since the .root file will always be all lowercase.
        * Runs the getRatesMulti.sh script (found in ratesScripts)for the filename.root file to get the HMS carbon rates pdfs and histograms from each foil at the given current and foil seperation for a multifoil target without sieve seperation cuts.
        * Runs the getRatesMSH.sh script (found in ratesScripts) for the filename.root file to get the HMS carbon rates pdfs and histograms from each foil for each sieve hole at the given current and foil seperation for a multifoil target WITH sieve seperation cuts.

* getRates.sh -
    * input -
        * filename - filename.root must be in worksim directory
        * current - current to get the inelastic carbon rates at in uA
    * description - 
        * Gets the total rate and prints select kinematic plots to pdf and root histogram files from a single carbon foil target without any sieve seperation cuts (can be used with or without sieve).
    * script flow - 
        * Ensures the right number of arguments was entered and that the path to filename.root exists.
        * Enters the pbmodel directory and ensures the libF1F209.so is made and linked with the files in the examples directory that uses it.
        * Clears the terminal.
        * Launches root and loads the object file from the pbmodel directory.
        * Executes the hms_foil_rates_carbon_inelastic macro for the input root file and the given current which is found in the example directory to get the rates, make the histograms, and print the pdfs.
        * Quits root.


* getRatesSieve - 
    * input -
        * filename - filename.root must be in worksim directory
        * current - current to get the inelastic carbon rates at in uA
    * description - 
        * Gets the per sieve hole rate (and time to 200 events per hole, etc.) and prints select kinematic plots to pdf and root histogram files for each sieve hole for a single carbon foil target WITH sieve seperation cuts.
    * script flow - 
        * Ensures the right number of arguments was entered and that the path to filename.root exists.
        * Enters the pbmodel directory and ensures the libF1F209.so is made and linked with the files in the examples directory that uses it.
        * Clears the terminal.
        * Launches root and loads the object file from the pbmodel directory.
        * Executes the hms_foil_rates_carbon_inelastic_sieveholes macro for the input root file and the given current which is found in the example directory to get the rates per sieve hole, make the histograms, and print the pdfs.
        * Quits root.

* getRatesMulti - 
    * input -
        * filename - filename.root must be in worksim directory
        * current - current to get the inelastic carbon rates at in uA
        * foilSep - serpation of carbon foils in cm
        * threeFoil - boolean determining if this was a 3 foil run (foil assumed to be at z=0).
    * description - 
        * Gets the total rate for each foil and prints select kinematic plots to pdf and root histogram files from a carbon multifoil target without any sieve seperation cuts (can be used with or without sieve).
    * script flow - 
        * Ensures the right number of arguments was entered and that the path to filename.root exists.
        * Enters the pbmodel directory and ensures the libF1F209.so is made and linked with the files in the examples directory that uses it.
        * Clears the terminal.
        * Launches root and loads the object file from the pbmodel directory.
        * Executes the hms_foil_rates_carbon_inelastic_multifoil macro for the input root file, the given current, the given foilSep, and whether there are three foils or not which is found in the example directory to get the rates per target foil, make the histograms, and print the pdfs.
        * Quits root.

* getRatesMSH.sh
    * input -
        * filename - filename.root must be in worksim directory
        * current - current to get the inelastic carbon rates at in uA
        * foilSep - serpation of carbon foils in cm
        * threeFoil - boolean determining if this was a 3 foil run (foil assumed to be at z=0).
    * description - 
        * Gets the per sieve hole rate from each foil source (and time to 200 events per hole, etc.) and prints select kinematic plots to pdf and root histogram files for each sieve hole for a carbon multifoil target WITH sieve seperation cuts.
    * script flow - 
        * Ensures the right number of arguments was entered and that the path to filename.root exists.
        * Enters the pbmodel directory and ensures the libF1F209.so is made and linked with the files in the examples directory that uses it.
        * Clears the terminal.
        * Launches root and loads the object file from the pbmodel directory.
        * Executes the hms_foil_rates_carbon_inelastic_MSH macro for the input root file, the given current, the given foilSep, and whether there are three foils or not which is found in the example directory to get the rates per sieve hole from each target foil, make the histograms, and print the pdfs.
        * Quits root.

Ntuple variables in SHMS hut ntuple ntuple id = 1411 
---------------------
* psxfp  Focal plane vertical position 
* psyfp  Focal plane horizontal position
* psxpfp Focal plane vertical angle=dx/dz
* psypfp Focal plane horizontal angle=dy/dz
* psztari  Initial random position along the target
* psytari  Initial y (horizontal) position at the z=0 plane perpendicular to central ray of SHMS
* psdeltai Initial random  delta = (p-pcent)/pcent
* psyptari Initial random  horizontal angle=dy/dz
* psxptari Initial random  vertical angle=dx/dz
* psztar  Reconstructed position along the target
* psytar   Reconstructed 
* psdelta   Reconstructed 
* psyptar  Reconstructed 
* psxptar Reconstructed 
* psxtari  Initial x (vertical) position at the z=0 plane perpendicular to central ray of SHMS
* fry   Initial random vertical raster position
* xsnum   sieve slit vertical hole number ( hole number at front sieve if use_front_sieve = true)
* ysnum   sieve slit horizontal hole number ( hole number at front sieve if use_front_sieve = true)
* xsieve   sieve slit vertical position (cm)  ( position at front sieve if use_front_sieve = true)
* ysieve   sieve slit horizontal position (cm) ( position at front sieve if use_front_sieve = true)

Ntuple variable in HMS hut ntuple ntuple id = 1
---------------------
* hsxfp Focal plane vertical position
* hsyfp Focal plane horizontal position
* hsxpfp Focal plane vertical angle=dx/dz
* hsypfp Focal plane horizaontal angle=dy/dz
* hsytari Initial y (horizontal) position at the z=0 plane perpendicular to the central ray of the HMS
* hsdeltai Initial random delta = (p-pcent)/(pcent)
* hsyptari Initial random horizontal angle=dy/dz
* hsxptari Initial random vertical angle=dx/dz
* hsytar Reconstructed horizontal position at the target
* hsdelta Reconstructed
* hsyptar Reconstructed horizontal angle=dy/dz at the target
* hsxptar Reconstructed vertical angle=dx/dz at the target
* fry Initial random vertical raster position
* ztari Initial random position along the target (z)
* stop_id flag marking valid event 0 if valid (can generate with only valid)
* xs_num Sieve hole passed through in the x (vertical) direction (NOT RECONSTRUCTED)
* ys_num Sieve hole passed through in the y (horizontal) direction (NOT RECONSTRUCTED)
* xc_sieve Reconstructed sieve hole passed through in the x (vertical) direction
* yc_sieve Reconstucted sieve hole passed through in the y (horizontal) direction
* wfac weighting factor - solid angle times energy divided by the number of trials (domega * denergy)/n_trials [rad * MeV]
* beam_E Beam Energy [GeV]
* p_spec Spectrometer Momentum  [GeV]
* th_spec Spectometer Angle [deg]