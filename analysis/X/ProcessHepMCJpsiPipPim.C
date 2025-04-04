//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.//

#include "HepMCElectro.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

//beam components must be defined even if 0
inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0.,0.,0,0};}
inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0.,0.,0.,0};}

void ProcessHepMCJpsiPipPim(){
  using namespace rad::names::data_type; //for MC()
  

  ///////////////////////////////////////////////////////////
  // Configure the data frame for the reaction
  // ep->e'+ X ( Jpsi(e+e-) + 2Pi (pi+pi-) )
  ///////////////////////////////////////////////////////////

  // we need treename, files
  rad::config::HepMCElectro hepmc{"hepmc3_tree", "../../hepmc_data/jpac_x3872_10_100_twohour.root"};
  hepmc.AliasStableMomentumComponents();
   
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(1);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronIndex(2);

  //give final state particle names,
  hepmc.setParticleIndex("el",rad::indice::useNthOccurance(3,11),{"mc_pid"});
  hepmc.setParticleIndex("po",rad::indice::useNthOccurance(1,-11),{"mc_pid"});
  hepmc.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{"mc_pid"});
  hepmc.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{"mc_pid"});
  hepmc.setParticleIndex("p",rad::indice::useNthOccurance(2,2212),{"mc_pid"});

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  hepmc.setBaryonParticles({"p"});
  
 //Use ParticleCreator to make intermediate states
  hepmc.Particles().Sum("J",{"el","po"});
  hepmc.Particles().Sum("rho",{"pip","pim"});
  hepmc.Particles().Sum("X",{"rho","J"});
  
  hepmc.setMesonParticles({"J","rho"});

  //must call this after all particles are configured
  hepmc.makeParticleMap();
  

  ///////////////////////////////////////////////////////////
  // Do some filtering on particle tracks
  ///////////////////////////////////////////////////////////
  // at least 1 pi+, pi-, e-, e+ and zero photons (remove omega decays)
  hepmc.Filter("(rad::helpers::Count(mc_pid,211)>=1)*(rad::helpers::Count(mc_pid,-211)>=1)*(rad::helpers::Count(mc_pid,11)>=1)*(rad::helpers::Count(mc_pid,-11)>=1)*(rad::helpers::Count(mc_pid,22)==0)","reaction_topo");

  ///////////////////////////////////////////////////////////
  // Declare kinematic calculations you want 
  ///////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(hepmc,"W","{scat_ele}");
  rad::rdf::Mass(hepmc,"Whad","{X,p}");
  rad::rdf::Mass(hepmc,"JMass","{J}");
  rad::rdf::Mass(hepmc,"XMass","{X}");
  rad::rdf::Mass(hepmc,"RhoMass","{rho}");

  //t distribution, column name
  rad::rdf::TBot(hepmc,"tb");
  rad::rdf::TPrimeBot(hepmc,"tbp");
  rad::rdf::TTop(hepmc,"tt");
  rad::rdf::TPrimeTop(hepmc,"ttp");

  //CM production angles
  rad::rdf::CMAngles(hepmc,"CM");
  rad::rdf::Q2(hepmc,"Q2");

  //decay angles
  rad::rdf::gn2s0s0s12::HelicityAngles(hepmc,"Heli");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",hepmc};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in mc_W between 4 and 54 GeV 
  //histo.Splitter().AddRegularDimension(MC()+"W", rad::histo::RegularSplits(10,4,54) );
  //can add as many split dimensions as we like
  //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
  histo.Init({MC()});//will create histograms for mc

  //define histograms
  histo.Create<TH1D,double>({"Q2","Q2",500,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,20.},{"W"});
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+, #pi#pi) [GeV]",100,3,4},{"XMass"});
  histo.Create<TH1D,double>({"JMass","M(e-,e+) [GeV]",100,.3,5.},{"JMass"});
  histo.Create<TH1D,double>({"RhoMass","M(2#pi) [GeV]",100,0,3},{"RhoMass"});
  histo.Create<TH1D,double>({"tb","t(p,p') [GeV^{2}]",100,-2,5},{"tb"});
  histo.Create<TH1D,double>({"tt","t(g,X) [GeV^{2}]",100,-2,5},{"tt"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,double>({"cthHeli","cos(#theta_{hel})",100,-1,1},{"Heli_CosTheta"});
  histo.Create<TH1D,double>({"phCHeli","#phi_{hel})",100,-TMath::Pi(),TMath::Pi()},{"Heli_Phi"});


  gBenchmark->Start("processing");
  //save all histograms to file
  histo.File("histos/HepMCJpsiPipPim_hists.root");

  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
 


  // hepmc.Snapshot("output/HepMCJpsiPipPim_hists.root");

}
