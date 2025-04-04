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

void ProcessHepMCKpKmPipPimPi0Pi0(){
  using namespace rad::names::data_type; //for MC()
  
  ///////////////////////////////////////////////////////////
  // Configure the data frame for the reaction
  // ep->e'+ X ( D0(K-pi+pi0) + D0bar (K+pi-) + pi0 )
  ///////////////////////////////////////////////////////////

  gBenchmark->Start("df");

  //create reaction dataframe
   rad::config::HepMCElectro hepmc{"hepmc3_tree", "../../hepmc_data/jpac_x3872_10_100_twohour.root"};

   hepmc.AliasStableMomentumComponents();
   
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
  hepmc.setBeamIonIndex(1);
  hepmc.setBeamElectronIndex(0);
  hepmc.setScatElectronIndex(2);
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  // hepmc.setParticleIndex("el",4);
  // hepmc.setParticleIndex("po",5);
  // hepmc.setParticleIndex("pim",7);
  // hepmc.setParticleIndex("pip",6);
  hepmc.setParticleIndex("Kp",rad::indice::useNthOccurance(1,321),{"mc_pid"});
  hepmc.setParticleIndex("Km",rad::indice::useNthOccurance(1,-321),{"mc_pid"});
  hepmc.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{"mc_pid"});
  hepmc.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{"mc_pid"});

  //2pi0 -> 4gamma
  hepmc.setParticleIndex("g1",rad::indice::useNthOccurance(1,22),{"mc_pid"});
  hepmc.setParticleIndex("g2",rad::indice::useNthOccurance(2,22),{"mc_pid"});
  hepmc.setParticleIndex("g3",rad::indice::useNthOccurance(3,22),{"mc_pid"});
  hepmc.setParticleIndex("g4",rad::indice::useNthOccurance(4,22),{"mc_pid"});

  hepmc.setParticleIndex("p",rad::indice::useNthOccurance(2,2212),{"mc_pid"});

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  hepmc.setBaryonParticles({"p"});
  
 //Use ParticleCreator to make intermediate states
  hepmc.Particles().Sum("pi0_1",{"g1","g2"});
  hepmc.Particles().Sum("pi0_2",{"g3","g4"});

  //Look for pi0 in other decay
  // hepmc.Particles().Sum("D0bar",{"Kp","pim","pi0_1"});
  //hepmc.Particles().Sum("D0",{"Km","pip"});
  
  hepmc.Particles().Sum("D0bar",{"Kp","pim"});
  hepmc.Particles().Sum("D0",{"Km","pip","pi0_1"});
  
  hepmc.Particles().Sum("X",{"D0","D0bar","pi0_2"});
  
  hepmc.setMesonParticles({"X"});

  //must call this after all particles are configured
  hepmc.makeParticleMap();
  
   ///////////////////////////////////////////////////////////
  // Do some filtering on particle tracks
  ///////////////////////////////////////////////////////////
  // TRUTH : Correct topology
  hepmc.Filter("mc_pid[Kp]==321&&mc_pid[Km]==-321&&mc_pid[pip]==211&&mc_pid[pim]==-211","event_pid"); //could add photons here too


  ///////////////////////////////////////////////////////////
  // Declare kinematic calculations you want 
  ///////////////////////////////////////////////////////////

    //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(hepmc,"W","{scat_ele}");
  rad::rdf::Mass(hepmc,"Whad","{X,p}");
  rad::rdf::Mass(hepmc,"D0Mass","{D0}");
  rad::rdf::Mass(hepmc,"D0barMass","{D0bar}");
  rad::rdf::Mass(hepmc,"Pi0_1Mass","{pi0_2}");
  rad::rdf::Mass(hepmc,"XMass","{X}");

  //t distribution, column name
  rad::rdf::TBot(hepmc,"tb");
  rad::rdf::TPrimeBot(hepmc,"tbp");
  rad::rdf::TTop(hepmc,"tt");
  rad::rdf::TPrimeTop(hepmc,"ttp");

  //CM production angles
  rad::rdf::CMAngles(hepmc,"CM");
  rad::rdf::Q2(hepmc,"Q2");

  rad::rdf::MissMass2(hepmc,"MissMass2","{scat_ele,p,X}");
  rad::rdf::MissMass(hepmc,"MissMassP","{scat_ele,p}");
 
 ///////////////////////////////////////////////////////////
  //Further filtering on calculated variables  
  ///////////////////////////////////////////////////////////
  //TRUTH : invariant mass cuts
  hepmc.Filter("mc_MissMass2<0.1","truthMM2");
  hepmc.Filter("mc_D0Mass>1.8&&mc_D0Mass<1.95","truthD0Mass");
  hepmc.Filter("mc_D0barMass>1.8&&mc_D0barMass<1.95","truthD0barMass");
  hepmc.Filter("mc_XMass>3.8&&mc_XMass<3.9","truthXMass");

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
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+, #pi#pi) [GeV]",500,0,5},{"XMass"});
  histo.Create<TH1D,double>({"D0Mass","M(K-#pi+) [GeV]",100,.3,5.},{"D0Mass"});
  histo.Create<TH1D,double>({"D0barMass","M(K+#pi-#pi0) [GeV]",100,.3,5.},{"D0barMass"});
  histo.Create<TH1D,double>({"Pi0_1Mass","M(2#gamma) [GeV]",100,0,0.5},{"Pi0_1Mass"});
  histo.Create<TH1D,double>({"tb","t(p,p') [GeV^{2}]",100,-2,5},{"tb"});
  histo.Create<TH1D,double>({"tt","t(g,X) [GeV^{2}]",100,-2,5},{"tt"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,double>({"MissMass2","Mmiss [GeV]",1000,-10,10},{"MissMass2"});
  histo.Create<TH1D,double>({"MissMassP","Mmiss [GeV]",1000,-10,10},{"MissMassP"});


  gBenchmark->Start("processing");
  //save all histograms to file
  histo.File("histos/HepMCKpKmPipPimPi0Pi0_hists.root");


  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
  gBenchmark->Stop("df");
  gBenchmark->Print("df");



  // hepmc.Snapshot("output/HepMCKpKmPipPimPi0Pi0_hists.root");

}
