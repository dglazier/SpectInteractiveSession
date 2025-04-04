//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "ePICReaction.h"
#include "HepMCElectro.h"
#include "Histogrammer.h"
#include "ParticleCreator.h"
#include "ePICParticleCreator.h"
#include "Indicing.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>

//with afterburner need slightly altered energies
inline constexpr std::array<double,4>  rad::beams::InitBotComponents() {return {0,0,0,0};}
inline constexpr std::array<double,4>  rad::beams::InitTopComponents() {return {0,0,0,0};}

void ePICJpsiPipPim(){
  
  using namespace rad::names; //for Scat_Ele()
  using namespace rad::names::data_type; //for Rec(), Truth()

  gBenchmark->Start("df");
  
  ///////////////////////////////////////////////////////////
  // Configure the data frame for the reaction
  // ep->e'+ Y ( Jpsi(e+e-) + 2Pi (pi+pi-) )
  ///////////////////////////////////////////////////////////

  // we need treename, files
  rad::config::ePICReaction epic{"events","../../simulated_data/jpac_y4260_18_275_10day_*_recon.root"};
   
  //as this is likely to be what we have in the experiment
  epic.SetBeamsFromMC();

  //will reorder all rec tracks into truth ordering
  epic.AliasColumnsAndMatchWithMC();
  
  
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file, ignoring beams
  epic.setScatElectronIndex(0);
  
  //hack the reconstructed lowQ2 electron to be MCMatched with this 
  rad::epic::ePICParticleCreator epic_particles{epic};
  epic_particles.MCMatchedLowQ2Electron();
  
   //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  epic.setParticleIndex("Jel",4,11);
  epic.setParticleIndex("Jpo",5,-11);
  epic.setParticleIndex("pip",2,211);
  epic.setParticleIndex("pim",3,-211);
  epic.setParticleIndex("pro",1,2212);
  

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calcualting reaction kinematics
  //e.g. t distributions
  //but must be done after ParticleMap
  //so currently cannot use as baryon...

  rad::config::ParticleCreator particles{epic};
  particles.Sum("Jpsi",{"Jel","Jpo"});
  particles.Sum("Rho",{"pip","pim"});
  particles.Sum("Meson",{"Jpsi","Rho"});
  
  epic.setMesonParticles({"Jpsi","Rho"});
  
  //can also add missing particles
  //And use those in calculations
  //Miss(name,{other final state particles})
  particles.Miss("CalcPro",{ScatEle().data(),"Meson"});
  epic.setBaryonParticles({"CalcPro"});

  //must call this after all particles are configured
  epic.makeParticleMap();

  ///////////////////////////////////////////////////////////
  // Do some filtering on particle tracks
  ///////////////////////////////////////////////////////////

  // tagged and reconstructed Y
  epic.Filter("rec_scat_ele>-1&&rec_pmag[Jel]>0.1&&rec_pmag[Jpo]>0.1&&rec_pmag[pip]>0.1&&rec_pmag[pim]>0.1","tagger");
  
  // make sure e+ e- decay of Jpsi
  epic.Filter("tru_pid[Jpo]==-11&&tru_pid[Jel]==11","epem");
  
  // tagged and untagged and reconstructed Y
  //epic.Filter("rec_pmag[Jel]>0.1&&rec_pmag[Jpo]>0.1&&rec_pmag[pip]>0.1&&rec_pmag[pim]>0.1","tagger");
  
  ///////////////////////////////////////////////////////////
  // Declare kinematic calculations you want 
  ///////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Mass(epic,"Whad","{Meson,pro}");
  rad::rdf::Mass(epic,"JMass","{Jpsi}");
  rad::rdf::Mass(epic,"MesonMass","{Meson}");
  rad::rdf::Mass(epic,"RhoMass","{Rho}");

  //t distribution, column name
  rad::rdf::TTop(epic,"t_meson");
  rad::rdf::TBot(epic,"t_baryon");
  rad::rdf::TPrimeBot(epic,"tp_baryon");
  rad::rdf::TPrimeTop(epic,"tp_meson");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");

 //decay angles
  rad::rdf::gn2s0s0s12::HelicityAngles(epic,"Heli");

  //exlusivity
  rad::rdf::MissMass(epic,"MissMass","{scat_ele,Meson,pro}");
  rad::rdf::MissMass(epic,"MissMassMeson","{scat_ele,Meson}");
  rad::rdf::MissP(epic,"MissP_Meson","{scat_ele,Meson}");
  rad::rdf::MissPt(epic,"MissPt_Meson","{scat_ele,Meson}");
  rad::rdf::MissPz(epic,"MissPz_Meson","{scat_ele,Meson}");
  rad::rdf::MissTheta(epic,"MissTheta_Meson","{scat_ele,Meson}");

  ///////////////////////////////////////////////////////////
  //Further filtering on calculated variables  
  ///////////////////////////////////////////////////////////
  epic.Filter("rec_MissP_Meson>270&&rec_MissP_Meson<276","missP");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",epic};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in tru_W between 4 and 54 GeV 
  // histo.Splitter().AddRegularDimension(Truth()+"W", rad::histo::RegularSplits(10,4,14) );
  //can add as many split dimensions as we like
  //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
  histo.Init({Rec(),Truth()});//will create histograms for rec and truth

  histo.Create<TH1D,double>({"Q2","Q2",500,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"W","W",100,0,150.},{"W"});
  histo.Create<TH1D,double>({"MesonMass","M(e-,e+, #pi+#pi-) [GeV]",100,3,5.},{"MesonMass"});
  histo.Create<TH1D,double>({"JMass","M(e-,e+) [GeV]",100,2.5,3.5},{"JMass"});
  histo.Create<TH1D,double>({"RhoMass","M(#pi+#pi-) [GeV]",100,0,2.},{"RhoMass"});
  histo.Create<TH1D,double>({"tmes","t(g,Meson) [GeV^{2}]",100,-2,5},{"t_meson"});
  histo.Create<TH1D,double>({"tbar","t(p,p') [GeV^{2}]",100,-2,5},{"t_baryon"});
  histo.Create<TH1D,double>({"cthCM","cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"phCM","#phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  histo.Create<TH1D,double>({"MissMass","Mmiss [GeV]",100,-40,40},{"MissMass"});
  histo.Create<TH1D,double>({"MissMassMeson","Mmiss [GeV]",100,-40,40},{"MissMassMeson"});
  histo.Create<TH1D,double>({"missP","p_{miss}(e',Meson)",105,250,305},{"MissP_Meson"});
  histo.Create<TH1D,double>({"missPt","p_{t,miss}(e',Meson)",100,0,3},{"MissPt_Meson"});
  histo.Create<TH1D,double>({"missPz","p_{z,miss}(e',Meson)",105,250,305},{"MissPz_Meson"});
  histo.Create<TH1D,double>({"missTheta","#theta_{miss}(e',Meson)",100,0,1e-2},{"MissTheta_Meson"});
  histo.Create<TH1D,float>({"EleP","p_{e}",100,0,20},{"pmag[scat_ele]"});
  histo.Create<TH1D,float>({"phHeli","#phi_{hel}",100,-TMath::Pi(),TMath::Pi()},{"Heli_Phi"});
  histo.Create<TH1D,float>({"thHeli","cos(#theta_{hel})",100,-1,1},{"Heli_CosTheta"});

  rad::histo::Histogrammer histo_res{"resolution",epic};
  histo_res.Init({"res_"});//will create histograms for rec and truth
  //Make already define resolutions of electron momentum
  histo_res.Create<TH1D,float>({"EleP","#Delta p_{e}",100,-0.2,0.2},{"pmag[scat_ele]"});
  histo_res.Create<TH1D,float>({"EleTh","#Delta #theta_{e}",100,-0.01,0.01},{"theta[scat_ele]"});
  histo_res.Create<TH1D,float>({"ElePhi","#Delta #phi_{e}",100,-TMath::Pi(),TMath::Pi()},{"phi[scat_ele]"});

  //Need to define resolutions for calculated quantities
  epic.Resolution("CM_Phi");
  histo_res.Create<TH1D,double>({"phCM","#Delta #phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  epic.Resolution("MissP_Meson");  
  histo_res.Create<TH1D,double>({"missP","#Delta p_{miss}(e',Meson)",100,-2,2},{"MissP_Meson"});
  epic.Resolution("Heli_Phi");
  histo_res.Create<TH1D,float>({"phHeli","#Delta #phi_{hel}",100,-0.1*TMath::Pi(),0.1*TMath::Pi()},{"Heli_Phi"});
  epic.Resolution("Heli_CosTheta");
  histo_res.Create<TH1D,float>({"thHeli","#Delta cos(#theta_{hel})",100,-0.1,0.1},{"Heli_CosTheta"});
  
  epic.Define("res_p_el","rec_pmag[rec_scat_ele]");
  histo_res.Create<TH2D,float,float>({"EleP2D","#Delta p_{e} v p_{e}",50,6,18,50,-0.05,0.05},{"p_el","pmag[scat_ele]"});

  
  ///////////////////////////////////////////////////////////
  // Save all histograms to file
  ///////////////////////////////////////////////////////////
  histo.File("histos/Jpsi2Pi_hists_missP_epem.root");
  histo_res.File("histos/Jpsi2Pi_reshists_missP_epem.root");


   // epic.Snapshot("Jpsi2Pi.root");
}
