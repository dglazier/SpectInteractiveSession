//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.//

#include "ePICReaction.h"
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

void ePICJpsiPipPim(){
  using namespace rad::names::data_type; //for Rec(), Tru()
  

  ///////////////////////////////////////////////////////////
  // Configure the data frame for the reaction
  // ep->e'+ X ( Jpsi(e+e-) + 2Pi (pi+pi-) )
  ///////////////////////////////////////////////////////////
  rad::config::ePICReaction epic{"events","../../simulated_data/jpac_x3872_10_100_halfday_*_recon.root"};
  epic.SetBeamsFromMC();

  epic.AliasColumnsAndMatchWithMC();

   
  //Assign particles names and indices
  //indicing comes from ordering in epic file
  epic.setScatElectronIndex(0);
 //hack the reconstructed lowQ2 electron to be MCMatched with this 
  rad::epic::ePICParticleCreator epic_particles{epic};
  epic_particles.MCMatchedLowQ2Electron();

  //give final state particle names,
  epic.setParticleIndex("el",rad::indice::useNthOccurance(2,11),{"tru_pid"});
  epic.setParticleIndex("po",rad::indice::useNthOccurance(1,-11),{"tru_pid"});
  epic.setParticleIndex("pim",rad::indice::useNthOccurance(1,-211),{"tru_pid"});
  epic.setParticleIndex("pip",rad::indice::useNthOccurance(1,211),{"tru_pid"});
  epic.setParticleIndex("p",rad::indice::useNthOccurance(1,2212),{"tru_pid"});

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  epic.setBaryonParticles({"p"});
  
 //Use ParticleCreator to make intermediate states
  epic.Particles().Sum("J",{"el","po"});
  epic.Particles().Sum("rho",{"pip","pim"});
  epic.Particles().Sum("X",{"rho","J"});
  
  epic.setMesonParticles({"J","rho"});

  //must call this after all particles are configured
  epic.makeParticleMap();
  
 
  ///////////////////////////////////////////////////////////
  // Do some filtering on particle tracks
  ///////////////////////////////////////////////////////////

  //Veto events with reconstructed photons
  //epic.Filter("rec_Ngamma==0","rec_particles_topo");//rec_pid not reliable, just vero gammas
  
  //Minimum momentum cut on reconstructed particles
  epic.Filter("(rec_pmag[el]>0.1)*(rec_pmag[po]>0.1)*(rec_pmag[pip]>0.1)*(rec_pmag[pim]>0.1)","rec_cut");
  
  // In truth at least 1 pi+, pi-, e-, e+ and zero photons (remove omega decays)
  //epic.Filter("(rad::helpers::Count(tru_pid,211)==1)*(rad::helpers::Count(tru_pid,-211)==1)*(rad::helpers::Count(tru_pid,11)>=1)*(rad::helpers::Count(tru_pid,-11)==1)*(rad::helpers::Count(tru_pid,22)==0)","reaction_topo");

  ///////////////////////////////////////////////////////////
  // Declare kinematic calculations you want 
  ///////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Mass(epic,"Whad","{X,p}");
  rad::rdf::Mass(epic,"JMass","{J}");
  rad::rdf::Mass(epic,"XMass","{X}");
  rad::rdf::Mass(epic,"RhoMass","{rho}");

  //t distribution, column name
  rad::rdf::TBot(epic,"tb");
  rad::rdf::TPrimeBot(epic,"tbp");
  rad::rdf::TTop(epic,"tt");
  rad::rdf::TPrimeTop(epic,"ttp");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");
  rad::rdf::Q2(epic,"Q2");

  //decay angles
  rad::rdf::gn2s0s0s12::HelicityAngles(epic,"Heli");

  ///////////////////////////////////////////////////////////
  //Further filtering on calculated variables  
  ///////////////////////////////////////////////////////////
  epic.Filter("rec_JMass>2.8","recJ");

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
  histo.Create<TH1D,double>({"W","W",100,0,20.},{"W"});
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+, #pi#pi) [GeV]",100,3,4},{"XMass"});
  histo.Create<TH1D,double>({"JMass","M(e-,e+) [GeV]",100,.3,5.},{"JMass"});
  histo.Create<TH1D,double>({"RhoMass","M(2#pi) [GeV]",100,0,3},{"RhoMass"});
  histo.Create<TH1D,double>({"tb","t(p,p') [GeV^{2}]",100,-2,5},{"tb"});
  histo.Create<TH1D,double>({"tt","t(g,X) [GeV^{2}]",100,-2,5},{"tt"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM})",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,float>({"cthHeli","cos(#theta_{hel})",100,-1,1},{"Heli_CosTheta"});
  histo.Create<TH1D,float>({"phCHeli","#phi_{hel})",100,-TMath::Pi(),TMath::Pi()},{"Heli_Phi"});


  rad::histo::Histogrammer histo_res{"resolution",epic};
  histo_res.Init({"res_"});//will create histograms for rec and truth
  histo_res.Create<TH1D,float>({"EleP","#Delta p_{e}",100,-0.2,0.2},{"pmag[scat_ele]"});
  histo_res.Create<TH1D,float>({"pipP","#Delta p_{#pi+}",100,-0.2,0.2},{"pmag[pip]"});
  histo_res.Create<TH1D,float>({"pimP","#Delta p_{#pi-}",100,-0.2,0.2},{"pmag[pim]"});
  histo_res.Create<TH1D,float>({"JelP","#Delta p_{e-}",100,-0.2,0.2},{"pmag[el]"});
  histo_res.Create<TH1D,float>({"JpoP","#Delta p_{e+}",100,-0.2,0.2},{"pmag[po]"});

  epic.Define("res_p_pip","rec_pmag[pip]");//just a hack prefixing with res
  histo_res.Create<TH2D,float,float>({"pipP2D","#Delta p v p #pi+",50,0,10,50,-0.5,0.5},{"p_pip","pmag[pip]"});
  epic.Define("res_p_pim","rec_pmag[pim]");//just a hack prefixing with res
  histo_res.Create<TH2D,float,float>({"pimP2D","#Delta p v p #pi-",50,0,10,50,-0.5,0.5},{"p_pim","pmag[pim]"});

  gBenchmark->Start("processing");
  //save all histograms to file
  histo.File("histos/EpicJpsiPipPim_hists_nogamveto.root");
  histo_res.File("histos/EpicJpsiPipPim_reshists_nogamveto.root");




  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
 


  // epic.Snapshot("output/EpicJpsiPipPim_hists.root");

}
