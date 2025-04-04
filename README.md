# Interactive Spectroscopy Analysis

https://eic.github.io/tutorial-analysis/index.html

This demonstrates how to analyse spectroscopy reactions using the
Reaction AWare Dataframes framework. There are options for both
HepMC generated files and full simulated ePIC data. The code for both
is essentially identical with a different choice of interfacing
high-level dataframe.

# Submodules/Dependencies

## Reaction Aware (R)Dataframes

https://github.com/dglazier/rad

## Setup

      git clone  --recurse-submodules  https://github.com/dglazier/SpectInteractiveSession.git
      cd SpectInteractiveSession
      export RAD=$PWD/rad/
      export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$RAD/include

## Simulated Data

Paths are set to analyse data from the SpectInteractiveSession/simulated_data directory. Please download files to this directory.

You may download some Y simulations files from

https://www.dropbox.com/scl/fo/j0skvooeknlyeqglpki4k/AOUcVUcukKw3FTXFSSEs3hc?rlkey=dsww27zy77oxyu0p75y60hvoj&st=b4oodal6&dl=0

You may download some X simulations files from

https://www.dropbox.com/scl/fo/rzdhayjxg7l1auxxaiduk/AGyHXdLUhdkt_kgRmSZOk1M?rlkey=8sy09il6rcwxy7rf73srfvlbw&st=s3e75x1x&dl=0

The chi_c1(3872) data include 4 possible branches,

     		 J/psi   rho0	     2.8 %
		 J/psi   omega	     4.1 %
		 D0	 D0_bar	pi0  45 %
		 D*0_bar D0	     34 %

The psi(4230) just contains J/psi pi+ pi- events.


## Processing HepMC data

Files are provided in the hepmc_data directory containing generated
event records in the HepMC3 format, converted to .root files.

rad scripts are provided to analyse these and produce histograms
of kinematic distributions. They should run without editing e.g.

      cd analysis/X
      root ProcessHepMCJpsiPipPim.C

And a histogram file will be created in the histos directory.

The rad interface to the hepmc format is the HepMCElectro
high level dataframe, e.g. :

     	   //create reaction dataframe, {treename, filename}
   	   rad::config::HepMCElectro hepmc{"hepmc3_tree", "root_files/jpac_x3872_10_100_halfday.root"};
	   


## Processing ePIC data

Full simulations have been performed for the chi_c1(3872), aka X, and psi(4230) aka Y, production.
You will need to download the files first. See intructions above
Once the files are in the correct place they should run without editing e.g.

      cd analysis/X
      root ePICJpsiPipPim.C

The rad interface to ePIC simulations is the ePICReaction 
high level dataframe, e.g. :

     	   //create reaction dataframe, {treename, filename}
   	   rad::config::ePICReaction epic{"events","../../simulated_data/jpac_x3872_10_100_halfday_*_recon.root"};