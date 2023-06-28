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
* h2root worksim/<filename>.rzdat worksim/<filename>.root (assuming the .rzdat file exists)

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
** If target length greater than 3cm assume LH2 cryotarget for multiple scattering calculation.
**  If target length less than 3cm assume simple solid target for multiple scattering calculation.
* calulated mutliple scattering in scattering chamber window, air and spectrometer entrance window.
* Call the mc_shms/mc_hms subroutine to transport the particle through SHMS/HMS apertures and to the focal plane.
* At the focal plane drift the particle through  1st Cerenkov or vacuum pipe (for SHMS only, depending on flag in input file),the drift chamber, scintillators, 2nd cerenkov and calorimeter to check if it hits them . 
* As it passes through the drift chamber randomly change the position according to the presumed DC resolution and fit the positions and calculate the focal plane positions and angles.
* From focal plane calculate the reconstructed target quantities ytar,yptar,xptar and delta using the optics matrix. Note that for now use the known x target position in the optics matrix calculations. To mimic what would be in the analysis, need to calculate xtar from the reconstructed quantities using the raster vertical position and then recalculate the reconstructed quantities with updated xtar.
* If particle makes it to the focal plane and hits all detectors the event variables are put in ntuple.



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
        * <filename> (<filename>.inp file must be in the infiles directory)
    * description - 
        * Runs mc-single-arm for for the given .inp file
        * Unpacks the <filename>.rzdat file into a <filename>.root file in the worksim directory
* buildAndRunSingle :
    * input - 
        * <filename> - <filename>.inp file must be in the infiles directory
        * current - current to get the inelastic carbon rates at in uA
    * description -
        * Runs the simulator and gets the pdfs and histograms for a single carbon foil target with and without sieve seperation cuts. 
    * script flow - 
        * Calls buildRoot.sh for <filename>
        * Runs the getRates.sh script (found in ratesScripts)for the <filename>.root file to get the HMS carbon rates pdfs and histograms for a single carbon foil without sieve seperation cuts.
        * Runs the getRatesSieve.sh script (found in ratesScripts) for the <filename>.root file to get the HMS carbon rates pdfs and histograms for a single carbon foil WITH sieve seperation cuts.
* buildAndRunMulti.sh :
    * input - 
        * <filename> - <filename>.inp file must be in the infiles directory
        * current - current to get the inelastic carbon rates at in uA
        * foilSep - serpation of carbon foils in cm
        * threeFoil - boolean determining if this was a 3 foil run (foil assumed to be at z=0).
    * description -
        * Runs the simulator and gets the pdfs and histograms for a carbon multifoil target with and without sieve seperation cuts.
    * script flow - 
        * Calls buildRoot.sh for <filename>
        * Runs the getRatesMulti.sh script (found in ratesScripts)for the <filename>.root file to get the HMS carbon rates pdfs and histograms from each foil at the given current and foil seperation for a multifoil target without sieve seperation cuts.
        * Runs the getRatesMSH.sh script (found in ratesScripts) for the <filename>.root file to get the HMS carbon rates pdfs and histograms from each foil for each sieve hole at the given current and foil seperation for a multifoil target WITH sieve seperation cuts.

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