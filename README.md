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

      git clone https://github.com/dglazier/SpectInteractiveSession.git
      cd SpectInteractiveSession
      export RAD=$PWD/rad/
      export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$RAD/include


## Processing HepMC data

Files are provided in the hepmc_data directory containing generated
event records in the HepMC3 format, converted to .root files.

rad scripts are provided to analyse these and produce histograms
of kinematic distributions.

The rad interface to the hepmc format is the HepMCElectro
high level dataframe, e.g. :

     	   //create reaction dataframe, {treename, filename}
   	   rad::config::HepMCElectro hepmc{"hepmc3_tree", "root_files/jpac_x3872_10_100_halfday.root"};
	   
The chi_c1(3872) data include 4 possible branches,

     		 J/psi   rho0	     2.8 %
		 J/psi   omega	     4.1 %
		 D0	 D0_bar	pi0  45 %
		 D*0_bar D0	     34 %

The psi(4230) just contains J/psi pi+ pi- events.


## Processing ePIC data

Full simulations have been performed for the chi_c1(3872), aka X, and psi(4230) aka Y, production.

The rad interface to ePIC simulations is the ePICReaction 
high level dataframe, e.g. :

     	   //create reaction dataframe, {treename, filename}
   	   rad::config::ePICReaction epic{"events","/home/dglazier/EIC/data/sim/jpac_z3900_10x100.root"};