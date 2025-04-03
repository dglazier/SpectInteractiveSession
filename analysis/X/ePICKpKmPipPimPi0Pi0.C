//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.//

#include "ePICReaction.h"
#include "ParticleCreator.h"
#include "ePICParticleCreator.h"
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

void ePICKpKmPipPimPi0Pi0(){
  using namespace rad::names::data_type; //for MC()
  
  gBenchmark->Start("df");

  //create reaction dataframe
  // rad::config::EpicElectro epic{"epic3_tree", "root_files/jpac_x3872_10_100_10day.root"};
  rad::config::ePICReaction epic{"events","/home/dglazier/EIC/data/sim/jpac_x3872_10_100_halfday_*_recon.root"};
   epic.SetBeamsFromMC();

  epic.AliasColumnsAndMatchWithMC();

  //Assign particles names and indices
  //indicing comes from ordering in epic file
  //epic.setBeamIonIndex(1);
  //epic.setBeamElectronIndex(0);
  epic.setScatElectronIndex(0);
 //hack the reconstructed lowQ2 electron to be MCMatched with this 
  rad::epic::ePICParticleCreator epic_particles{epic};
  epic_particles.MCMatchedLowQ2Electron();

  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  // epic.setParticleIndex("el",4);
  // epic.setParticleIndex("po",5);
  // epic.setParticleIndex("pim",7);
  // epic.setParticleIndex("pip",6);
  epic.setParticleIndex("Kp",rad::indice::useNthOccurance(1,321),{"tru_pid"});
  epic.setParticleIndex("Km",rad::indice::useNthOccurance(1,-321),{"tru_pid"});
  epic.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{"tru_pid"});
  epic.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{"tru_pid"});

  //2pi0 -> 4gamma
  epic.setParticleIndex("g1",rad::indice::useNthOccurance(1,22),{"tru_pid"});
  epic.setParticleIndex("g2",rad::indice::useNthOccurance(2,22),{"tru_pid"});
  epic.setParticleIndex("g3",rad::indice::useNthOccurance(3,22),{"tru_pid"});
  epic.setParticleIndex("g4",rad::indice::useNthOccurance(4,22),{"tru_pid"});

  epic.setParticleIndex("p",rad::indice::useNthOccurance(1,2212),{"tru_pid"});

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  epic.setBaryonParticles({"p"});
  
 //Use ParticleCreator to make intermediate states
  epic.Particles().Sum("pi0_1",{"g1","g2"});
  epic.Particles().Sum("pi0_2",{"g3","g4"});
  // epic.Particles().Sum("D0bar",{"Kp","pim","pi0_1"});
  //epic.Particles().Sum("D0",{"Km","pip"});
  epic.Particles().Sum("D0bar",{"Kp","pim"});
  epic.Particles().Sum("D0",{"Km","pip","pi0_1"});
  
  epic.Particles().Sum("X",{"D0","D0bar","pi0_2"});
  
  epic.setMesonParticles({"X"});

  //must call this after all particles are configured
  epic.makeParticleMap();
  
  
  ///////////////////////////////////////////////////////////
  // Do some filtering on particle tracks
  ///////////////////////////////////////////////////////////
  // at least 1 pi+, pi-, e-, e+ and zero photons (remove omega decays)
  epic.Filter("(rad::helpers::Count(tru_pid,211)==1)*(rad::helpers::Count(tru_pid,-211)==1)*(rad::helpers::Count(tru_pid,321)==1)*(rad::helpers::Count(tru_pid,-321)==1)*(rad::helpers::Count(tru_pid,22)==4)","reaction_topo");
  // epic.Filter("rec_scat_ele>-1","lowQ2");
  epic.Filter("tru_pid[Kp]==321&&tru_pid[Km]==-321&&tru_pid[pip]==211&&tru_pid[pim]==-211","event_pid");
  ///////////////////////////////////////////////////////////
  // Declare kinematic calculations you want 
  ///////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Mass(epic,"Whad","{X,p}");
  rad::rdf::Mass(epic,"D0Mass","{D0}");
  rad::rdf::Mass(epic,"D0barMass","{D0bar}");
  rad::rdf::Mass(epic,"Pi0_1Mass","{pi0_2}");
  rad::rdf::Mass(epic,"XMass","{X}");

  //t distribution, column name
  rad::rdf::TBot(epic,"tb");
  rad::rdf::TPrimeBot(epic,"tbp");
  rad::rdf::TTop(epic,"tt");
  rad::rdf::TPrimeTop(epic,"ttp");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");
  rad::rdf::Q2(epic,"Q2");

  rad::rdf::MissMass2(epic,"MissMass2","{scat_ele,p,X}");
  rad::rdf::MissMass(epic,"MissMassP","{scat_ele,p}");
 
  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",epic};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in tru_W between 4 and 54 GeV 
  //histo.Splitter().AddRegularDimension(MC()+"W", rad::histo::RegularSplits(10,4,54) );
  //can add as many split dimensions as we like
  //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
  histo.Init({Truth(),Rec()});//will create histograms for mc

  //define histograms
  histo.Create<TH1D,double>({"Q2","Q2",500,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,50.},{"W"});
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
  histo.File("histos/EpicKpKmPipPimPi0Pi0_hists.root");




  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
  gBenchmark->Stop("df");
  gBenchmark->Print("df");



  // epic.Snapshot("output/EpicKpKmPipPimPi0Pi0_hists.root");

}
