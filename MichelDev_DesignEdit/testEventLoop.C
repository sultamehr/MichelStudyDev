//==============================================================================
// MINERvA Analysis Toolkit "Minimum Adoption" Event Loop Example
//
// "Minimum adoption" means that it only uses the two essential MAT tools:
// Universes and HistWrappers. For an "full adoption" example that additionally
// makes use of Cuts, MacroUtil, and Variable, refer to the example in
// ../MAT_Tutorial/.
//
// This script follows the canonical event-looping structure:
// Setup (I/O, variables, histograms, systematics)
// Loop events
//   loop universes
//     make cuts
//     loop variables
//       fill histograms
// Plot and Save
//==============================================================================


//PlotUtils includes
//No junk from PlotUtils please!  I already
//know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

#include "CVUniverse.h"
#include "Systematics.h"

#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/MacroUtil.h"
#ifndef __CINT__
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/cuts/CCInclusiveCuts.h"
#include "PlotUtils/cuts/CCInclusiveSignal.h"
#include "Categorized.h"
#include "PlotUtils/Cutter.h"
#include "Variable.h"
#include "BestMichelDistance.h"
#include "BestMichelDistance2D.h"
#include "SignalDefinition.h"
#include "q3RecoCut.h"
#include "hasMichel.h"
#include "hasPion.h"
#include "studies/Study.h"
#include "studies/PerMichelVarByGENIELabel.h"
#include "studies/PerMichelEventVarByGENIELabel.h" 
#include "studies/CreateDataHistPerMichelEvent.h"
#include "studies/CreateDataHistPerMichel.h"
#endif //CINT
#include "Michel.h"
#include "Cluster.h"
#include "MichelEvent.h"
#include "MatchedMichel.h"
#include "Pion.h"
#include "PionEvent.h"
//#include "TruthMatchParticles/RecoElectron.h"
//#include "TruthMatchParticles/RecoEPrimaryParent.h"
//#include "Binning.h" //TODO: Fix me
#pragma GCC diagnostic pop
#include <iostream>

class Variable;


// Histogram binning constants
const int nbins = 30;
const double xmin = 0.;
const double xmax = 20.e3;

//==============================================================================
// Cuts
//==============================================================================




//==============================================================================
// Plot
//==============================================================================
void PlotErrorSummary(PlotUtils::MnvH1D* h, std::string label) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  PlotUtils::MnvPlotter mnv_plotter( kCompactStyle);
  TCanvas cE("c1", "c1");

  // Group GENIE bands
    mnv_plotter.error_summary_group_map.clear();
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrAbs_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrAbs_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrCEx_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrCEx_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrElas_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrElas_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrInel_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrPiProd_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back(
        "GENIE_FrPiProd_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_MFP_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_MFP_pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_AGKYxF1pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_AhtBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_BhtBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CCQEPauliSupViaKF");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CV1uBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CV2uBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_EtaNCEL");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MaNCEL");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MaRES");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MvRES");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_NormDISCC");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_NormNCRES");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_RDecBR1gamma");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn1pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn2pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn3pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvp1pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvp2pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Theta_Delta2Npi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_VecFFCCQEshape");

    mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_HighQ2");
    mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_LowQ2");

  const bool do_fractional_uncertainty = true;
  const bool do_cov_area_norm = false;
  const bool include_stat_error = false;
  const std::string do_fractional_uncertainty_str =
      do_fractional_uncertainty ? std::string("Frac") : std::string("Abs");
  const std::string do_cov_area_norm_str =
      do_cov_area_norm ? std::string("CovAreaNorm") : std::string("");

  mnv_plotter.DrawErrorSummary(hist, "TR", include_stat_error, true,
                              0.0, do_cov_area_norm, "",
                              do_fractional_uncertainty);
  mnv_plotter.AddHistoTitle("Event Selection");
  std::string plotname =
      Form("ErrorSummary_%s_%s_%s", do_fractional_uncertainty_str.c_str(),
           do_cov_area_norm_str.c_str(), label.c_str());
  mnv_plotter.MultiPrint(&cE, plotname, "png");
}

void PlotCVAndError(PlotUtils::MnvH1D* h, std::string label) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  if (!gPad)
    throw std::runtime_error("Need a TCanvas. Please make one first.");
  PlotUtils::MnvH1D* datahist = new PlotUtils::MnvH1D(
      "adsf", "", nbins, xmin, xmax);
  bool statPlusSys = true;
  int mcScale = 1.;
  bool useHistTitles = false;
  const PlotUtils::MnvH1D* bkgdHist = NULL;
  const PlotUtils::MnvH1D* dataBkgdHist = NULL;
  mnv_plotter.DrawDataMCWithErrorBand(datahist, hist, mcScale, "TL",
                                     useHistTitles, NULL, NULL, false,
                                     statPlusSys);
  mnv_plotter.AddHistoTitle("Event Selection");
  std::string plotname = Form("CV_w_err_%s", label.c_str());
  mnv_plotter.MultiPrint(&cE, plotname, "png");
  delete datahist;
}

#ifndef __CINT__
//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)
    //util::Categorized<PlotUtils::HistWrapper<CVUniverse>, int> hw_bestpirange) //,PlotUtils::HistWrapper<CVUniverse> hw_emu)
    {
    for (int i=0; i<chain->GetEntries(); ++i) {
    if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;


    int isSignal = 0;
    
    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands) {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes) {
         MichelEvent myevent; // make sure your event is inside the error band loop. 
    
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        
        // THis is where you would Access/create a Michel

          //std::cout << "Best Michel Distance is " << min_dist << std::endl;
          
          if (!michelcuts.isMCSelected(*universe, myevent, 1).all()) continue; //all is another function that will later help me with sidebands
          std::vector<Pion> npiplus = universe->MakeNPrimePions();
          if (npiplus.empty()) continue;
          if (npiplus.size() == 0) continue;
          //For each Variable
          //if(best_michel != -1) //Equivalent to if (return_michels.size() == 0) continue;
        
          for(auto& var: vars)
          {
            for(auto& study: studies) study->SelectedSignal(*universe, myevent, 1); //TODO: Last argument is weight
            //TODO: #include "study.h"

            //#include "MichelEventStudy.h"	
	    (*var->m_bestPionByGENIELabel)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), universe->GetWeight());

            //Fill other per-Variable histograms here
            
          }
          
        //=========================================
        // Cuts in each universe
        //=========================================
        //if(PassesCuts(*universe)) {
        //  hw_emu.univHist(universe)->Fill(universe->GetEmu(), 
                                        //  universe->GetWeight());
          
        //} // End if passed cuts
      } // End band's universe loop
    } // End Band loop
    
    //if (isSignal == 1) std::cout << "THIS EVENT HAS A SIGNAL MICHEL " << std::endl;

  } //End entries loop
}


void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable*> vars,
                                std::vector<Study*> studies,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)

{
      for (auto universe : data_band) {
   		for (int i=0; i<data->GetEntries(); ++i) {
     			universe->SetEntry(i);
     			if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;
     			MichelEvent myevent; 
     			if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
     			for(auto& study: studies) study->Selected(*universe, myevent, 1); 
                        for(auto& var: vars)
     			{
        			(var->dataHist)->FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
     			}
   		}
      }
}



void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    				std::vector<Variable*> vars,
    				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)
{
    for (int i=0; i<truth->GetEntries(); ++i) {
      if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;


      
      //=========================================
      // Systematics loop(s)
      //=========================================
      for (auto band : truth_bands) {
        std::vector<CVUniverse*> truth_band_universes = band.second;
        for (auto universe : truth_band_universes) {
      
          // Tell the Event which entry in the TChain it's looking at
          universe->SetEntry(i);

          if (!michelcuts.isEfficiencyDenom(*universe, 1)) continue; 
          //Fill efficiency denominator now: 
          //TODO ADD PLOTS
	  
          

	}

      }

    }
}




#endif // __CINT__
//==============================================================================
// Main
//==============================================================================
void testEventLoop() {
  TH1::AddDirectory(false);
  // Make a chain of events
  //PlotUtils::ChainWrapper* chain = makeChainWrapperPtr("CCQENu_minervame1A_MC_Inextinguishable_merged.txt", 
  //                                                     "CCQENu");
  //PlotUtils::ChainWrapper* truth = makeChainWrapperPtr("CCQENu_minervame1A_MC_Inextinguishable_merged.txt",
  //						       "Truth");
  //PlotUtils::ChainWrapper* data = makeChainWrapperPtr("CCQENu_minervame1A_DATA_Inextinguishable_merged.txt",
  //						      "CCQENu");
  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr("shorttest.txt",
                                                       "CCQENu");
  PlotUtils::ChainWrapper* truth = makeChainWrapperPtr("shorttest.txt",
                                                       "Truth");
  //PlotUtils::ChainWrapper* data = makeChainWrapperPtr("CCQENu_minervame1A_DATA_Inextinguishable_merged.txt",
  //                                                    "CCQENu");


  
  const std::string mc_file_list("shorttest.txt");
  const std::string data_file_list("shorttest.txt");
  const std::string plist_string("minervame1a");
  const std::string reco_tree_name("CCQENu");
  
  const bool do_truth = false;
  const bool is_grid = false;
  
 // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(false);
  PlotUtils::MinervaUniverse::SetPlaylist("minervame1A");
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(false);
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);

  // Make a map of systematic universes
  std::map< std::string, std::vector<CVUniverse*> > error_bands = 
      GetStandardSystematics(chain);

  std::map< std::string, std::vector<CVUniverse*> > truth_bands =
      GetStandardSystematics(truth);

  // Create a File for output hists:
  //
  //TFile* outDir = TFile::Open("nocuts/RequirePiPlus_MichelFitXY_MCHists.root", "RECREATE");
  TFile* outDir = TFile::Open("ShortTruthSTudies_MCHists.root", "RECREATE");

  TFile* dataDir = TFile::Open("DataMichelStudyHists.root", "RECREATE");
  //ExclusiveVariable1Arg<CVUniverse, Variable> is there to forward a Michel index
  //to CVUniverse::GetPionKE().  It derives from Variable and thus has the
  //histogram(s) you want to fill.
  //std::vector<Variable*> vars = { 
  //new ExclusiveVariable1Arg<CVUniverse, Variable>("tpi", "T#pi [MeV]", 100, 0., 1., &CVUniverse::GetPionKE, &CVUniverse::GetPionKE)
  //};


 
  std::vector<Variable*> vars = {
  new Variable("tpi", "T#pi [GeV]", 500, 0., 1.0, &CVUniverse::GetLowTpi, &CVUniverse::GetLowTpi),
  new Variable("q3", "q3 (GeV)", 100, 0.0, 1.3, &CVUniverse::Getq3, &CVUniverse::Getq3)
  };
  std::vector<Study*> studies;

  #ifndef __CINT__
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> oneVar = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   return evt.m_nmichels[whichMichel]->Best3Ddist;
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> michelE = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   return evt.m_nmichels[whichMichel]->energy;
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
				 {
                                   double micheltime = evt.m_nmichels[whichMichel]->time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
				   return deltat;
				 };

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_dist = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
				{
				  double micheldist = evt.m_nmichels[whichMichel]->up_to_vertex_dist3D;
 				  return micheldist;
				};
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_dist = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_vertex_dist3D;
                                  return micheldist;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clus_range = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_cluster_dist3D;
                                  return micheldist;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clus_range = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_cluster_dist3D;
                                  return micheldist;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clus_mich = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_clus_michel_dist3D;
                                  return micheldist;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clus_mich = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_clus_michel_dist3D;
                                  return micheldist;
                                };

   std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_dist = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
				{
				  double micheldist = evt.m_nmichels[whichMichel]->Best3Ddist;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 3 || match == 4) return micheldist;
				  else return 9999.;
				};
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_dist = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->Best3Ddist;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 1 || match == 2) return micheldist;
                                  else return 9999.;
                                };


 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_mich_dist = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 3 || match == 4){
				    micheldist = evt.m_nmichels[whichMichel]->up_clus_michel_dist3D;
				    return micheldist;
				  }
                                  else return 9999.;
                                };
   std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_XZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->best_XZ;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 1 || match == 2) return micheldist;
				  else return 9999.;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_XZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->best_XZ;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 3 || match == 4) return micheldist;
				  else return 9999.;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_UZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->best_UZ;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 1 || match == 2) return micheldist;
				  else return 9999.;
                                };

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_UZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->best_UZ;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 3 || match == 4) return micheldist;
				  else return 9999.;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_VZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->best_VZ;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 1 || match == 2) return micheldist;
				  else return 9999.;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_VZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->best_VZ;
                                  int match = (evt.m_nmichels[whichMichel]->BestMatch);
                                  if (match == 3 || match == 4) return micheldist;
				  else return 9999.;
                                };


  std::function<double(const CVUniverse&)> low_tpi = [](const CVUniverse& univ)
				{
				   double lowesttpi = 9999.;
                                   int nFSpi = univ.GetNTruePions();
                                   for (int i = 0; i < nFSpi; i++){
         			       int pdg = univ.GetPionPDG(i);
				       int trackid = univ.GetPionParentID(i);
				       if (trackid != 0 ) continue;
  				       if (pdg != 211) continue;
 				       double pionKE = univ.GetPionKE(i);
                                       if (lowesttpi > pionKE) lowesttpi = pionKE;
				   }
                                   return lowesttpi;
				};
  std::function<double(const CVUniverse&, const MichelEvent&)> bestdistvar = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   return evt.m_bestdist;
                                 }; 
  
  std::function<double(const CVUniverse&, const MichelEvent&)> vtxdist = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
  				   int matchtype = evt.m_matchtype;
	                           if (matchtype == 1 || matchtype == 2) return evt.m_bestdist;
	                           else return -9999.;
                                 };
   std::function<double(const CVUniverse&, const MichelEvent&)> clusdist = [](const CVUniverse& univ, const MichelEvent& evt)
                                 { 
 				   int bestidx = evt.m_idx;       
                                   int matchtype = evt.m_matchtype;
                                   if (matchtype == 3 || matchtype == 4) return evt.m_bestdist;
                                   else return -9999.;
                                 };

   std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_XZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_vertex_XZ;
                                  return micheldist;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_XZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_vertex_XZ;
			          return micheldist;
				};
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_UZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_vertex_UZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_UZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_vertex_UZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_VZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_vertex_VZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_VZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_vertex_VZ;
                                  return micheldist;
                                };
  
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_XZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_clus_XZ;
                                  return micheldist;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_XZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_clus_XZ;
                                  return micheldist;
                                };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_UZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_clus_UZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_UZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_clus_UZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_VZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->up_to_clus_VZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_VZ = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel]->down_to_clus_VZ;
                                  return micheldist;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> overlay_frac = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double frac = evt.m_nmichels[whichMichel]->overlay_fraction;
                                  std::cout << "Printing overlay fraction : " << frac << std::endl;
 				  return frac;
                                };
std::function<int(const CVUniverse&, const MichelEvent&, const int)> is_overlay = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int isoverlay = evt.m_nmichels[whichMichel]->is_overlay;
                                  std::cout << "Printing overlay boolean " << isoverlay << std::endl;
                                  return isoverlay;
                                };

 std::function<double(const CVUniverse&, const PionEvent&)> minKE = [](const CVUniverse& univ, const PionEvent& evt)             
                                {
                                  double KE = evt.m_lowKE;;
                                  return KE;
                                };
   
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_x = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->up_vtx_x;
    return x; 
  };
   std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_y = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double y = evt.m_nmichels[whichMichel]->up_vtx_y;
    return y;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_z = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double z = evt.m_nmichels[whichMichel]->up_vtx_z;
    return z;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_x = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->down_vtx_x;
    return x;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_y = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double y = evt.m_nmichels[whichMichel]->down_vtx_y;
    return y;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_z = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double z = evt.m_nmichels[whichMichel]->down_vtx_z;
    return z;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clus_x = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->up_clus_x;
    return x;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clus_y = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double y = evt.m_nmichels[whichMichel]->up_clus_y;
    return y;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clus_z = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double z = evt.m_nmichels[whichMichel]->up_clus_z;
    return z;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clus_x = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->down_clus_x;
    return x;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clus_y = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double y = evt.m_nmichels[whichMichel]->down_clus_y;
    return y;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clus_z = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double z = evt.m_nmichels[whichMichel]->down_clus_z;
    return z;
  };

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> true_x = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->true_initialx;
    return x;
  };
   std::function<double(const CVUniverse&, const MichelEvent&, const int)> true_y = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double y = evt.m_nmichels[whichMichel]->true_initialy;
    return y;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> true_z = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double z = evt.m_nmichels[whichMichel]->true_initialz;
    return z;
  };

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> truerecodiff_x1 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->true_initialx;
    double recx = evt.m_nmichels[whichMichel]->m_x1;
    double diff = x - recx;
    return diff;
  };
   std::function<double(const CVUniverse&, const MichelEvent&, const int)> truerecodiff_y1 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double y = evt.m_nmichels[whichMichel]->true_initialy;
    double recy = evt.m_nmichels[whichMichel]->m_y1;
    double diff = y - recy;
    return diff;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> truerecodiff_z1 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double z = evt.m_nmichels[whichMichel]->true_initialz;
    double recz = evt.m_nmichels[whichMichel]->m_z1;
    double diff = z - recz;
    return diff;
  };

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> truerecodiff_x2 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  {
    double x = evt.m_nmichels[whichMichel]->true_initialx;
    double recx = evt.m_nmichels[whichMichel]->m_x2;
    double diff = x - recx;
    return diff;
  };
   std::function<double(const CVUniverse&, const MichelEvent&, const int)> truerecodiff_y2 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  { 
    double y = evt.m_nmichels[whichMichel]->true_initialy;
    double recy = evt.m_nmichels[whichMichel]->m_y2;
    double diff = y - recy;
    return diff;
  };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> truerecodiff_z2 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
  { 
    double z = evt.m_nmichels[whichMichel]->true_initialz;
    double recz = evt.m_nmichels[whichMichel]->m_z2;
    double diff = z - recz;
    return diff;
  };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff1 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
			   	   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
				   else{
                                   double x = evt.m_nmichels[bestidx]->m_x1;
                                   double tx = evt.m_nmichels[bestidx]->true_initialx;
                                   return tx-x;
				   }
                                 };

 std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff1 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
				   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
		    		   else{
				   double y = evt.m_nmichels[bestidx]->m_y1;
                                   double ty = evt.m_nmichels[bestidx]->true_initialy;
                                   return ty-y;
				   }
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff1 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
				   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
				   else {
				   double z = evt.m_nmichels[bestidx]->m_z1;
                                   double tz = evt.m_nmichels[bestidx]->true_initialz;
                                   return tz-z;
				   }
                                 };


  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff2 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
			           else{
				   double x = evt.m_nmichels[bestidx]->m_x2;
                                   double tx = evt.m_nmichels[bestidx]->true_initialx;
                                   return tx-x;
				   }
                                 };

 std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff2 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
				   else{
				   double y = evt.m_nmichels[bestidx]->m_y2;
                                   double ty = evt.m_nmichels[bestidx]->true_initialy;
                                   return ty-y;
				   }
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff2 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 { 
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
				   else {
				   double z = evt.m_nmichels[bestidx]->m_z2;
                                   double tz = evt.m_nmichels[bestidx]->true_initialz;
                                   return tz-z;
				   }
                                 };

  studies.push_back(new PerMichelVarByGENIELabel(oneVar, "Best Distance", "mm", 100, 0., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(michelE, "Michel energy", "MeV", 20, 0., 100., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(delta_t, "Michel Time Diff", "#mus", 30, 0.0, 9.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_dist, "Best Vertex to Cluster 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_dist, "Best Vertex to Michel 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_dist, "Best Michel to Cluster 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));

  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_dist, "Up Vertex to Michel 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_dist, "Down Vertex to Michel 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clus_range, "Up Cluster to Vertex 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clus_range, "Down Cluster to Vertex 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clus_mich, "Up Cluster to Michel 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clus_mich, "Down Cluster to Michel 3D Distance", "mm", 100, 0.0, 1000.0, error_bands));

  studies.push_back(new PerMichelVarByGENIELabel(vtx_XZ, "Best Vertex to Michel XZ Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_UZ, "Best Vertex to Michel UZ Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_VZ, "Best Vertex to Michel VZ Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_XZ, " Best Cluster to Michel XZ Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_UZ, " Best Cluster to Michel UZ Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_VZ, " Best Cluster to Michel VZ Distance", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_XZ, "Up Vtx-Michel XZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_XZ, "Down Vtx-Michel XZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_UZ, "Up Vtx-Michel UZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_UZ, "Down Vtx-Michel UZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_VZ, "Up Vtx-Michel VZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_VZ, "Down Vtx-Michel VZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_XZ, "Up Clus-Michel XZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_XZ, "Down Clus-Michel XZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_UZ, "Up Clus-Michel UZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_UZ, "Down Clus-Michel UZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_VZ, "Up Clus-Michel VZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_VZ, "Down Clus-Michel VZ", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(overlay_frac, "Overlay Fraction", "no units", 20, -2., 2., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(is_overlay, "Is Overlay", "no units", 5, -2, 2, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_x, "Endpoint 1 X", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_y, "Endpoint 1 Y", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_z, "Endpoint 1 Z", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_x, "Endpoint 2 X", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_y, "Endpoint 2 Y", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_z, "Endpoint 2 Z", "mm", 1000, -10000., 10000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clus_x, "Endpoint 1 X (Clus Matched)", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clus_y, "Endpoint 1 Y (clus Matched)", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clus_z, "Endpoint 1 Z (Clus Matched)", "mm", 100, 4000., 10000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clus_x, "Endpoint 2 X (Clus Matched)", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clus_y, "Endpoint 2 Y (clus Matched)", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clus_z, "Endpoint 2 Z (Clus Matched)", "mm", 100, 4000., 10000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(true_x, "True Initial X", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(true_y, "True Initial Y", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(true_z, "True Initial Z", "mm", 100,0, 10000., error_bands));

  studies.push_back(new PerMichelVarByGENIELabel(truerecodiff_x1, "True Initial - Michel(1) X ", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(truerecodiff_y1, "True Initial - Michel(1) Y ", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(truerecodiff_z1, "True Initial - Michel(1) Z ", "mm", 500, -10000., 10000., error_bands));

  studies.push_back(new PerMichelVarByGENIELabel(truerecodiff_x2, "True Initial - Michel(2) X", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(truerecodiff_y2, "True Initial - Michel(2) Y", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(truerecodiff_z2, "True Initial - Michel(2) Z", "mm", 500, -10000., 10000., error_bands));

  studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff1, "True Initial - Best Michel  X1", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff1, "True Initial - Best Michel  Y1", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff1, "True Initial - Best Michel  Z1", "mm", 500, -5000., 5000., error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff2, "True Initial - Best Michel  X2", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff2, "True Initial - Best Michel  Y2", "mm", 100, -1000., 1000., error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff2, "True Initial - Best Michel  Z2", "mm", 200, -5000., 5000., error_bands));

  //studies.push_back(new PerMichelEventVarByGENIELabel(minKE, "Lowest Energy T#pi", "MeV", 100, 0., 1000., error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestdistvar, "Best Pion Range in Event", "mm", 100, 0., 1000., error_bands));
 // studies.push_back(new PerMichelEventVarByGENIELabel(clusdist, "Vertex - Cluster Distance", "mm", 100, 0., 1000., error_bands));
  //vars.push_back( new Variable("tpi", "T#pi [MeV]", 100, 0., 1., bestdistvar, bestdistvar));
  //std::vector<Variable*> vars = {
  //new ExclusiveVariable1Arg<CVUniverse, Variable>("tpi", "T#pi [MeV]", 100, 0., 1., low_tpi, low_tpi)
  //};
  #endif //__CINT__

  //new ExclusiveVariable1Arg<CVUniverse, Variable>("pionrange", 100, 0., 2400, &Michel::Best3Ddist, &CVUniverse::Best3Ddist)
  
  //using namespace PlotUtils;
  //Creating the single Data universe 
  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list,
                    plist_string, do_truth, is_grid);
  CVUniverse* data_universe = new CVUniverse(util.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv_data"] = data_band;
  
  #ifndef __CINT__
  std::vector<Study*> data_studies;
  // data_studies.push_back(new CreateDataHistPerMichel(oneVar, "Best Distance", "mm", 100, 0, 1000., data_band));
  // data_studies.push_back(new CreateDataHistPerMichel(michelE, "Michel energy", "MeV", 20, 0., 100., data_band));
  // data_studies.push_back(new CreateDataHistPerMichel(delta_t, "Michel Time Diff", "#mus", 30, 0.0, 9.0, data_band));
  // data_studies.push_back(new CreateDataHistPerMichelEvent(bestdistvar, "Best Distance", "mm", 100, 0, 1000., data_band));
   
  data_studies.push_back(new PerMichelVarByGENIELabel(oneVar, "Best Distance DATA", "mm", 100, 0., 1000., data_error_bands));
  data_studies.push_back(new PerMichelVarByGENIELabel(michelE, "Michel energy DATA", "MeV", 20, 0., 100., data_error_bands));
  data_studies.push_back(new PerMichelVarByGENIELabel(delta_t, "Michel Time Diff DATA", "#mus", 30, 0.0, 9.0, data_error_bands));
  data_studies.push_back(new PerMichelEventVarByGENIELabel(bestdistvar, "Best Distance DATA", "mm", 100, 0., 1000., data_error_bands));
  data_studies.push_back(new PerMichelEventVarByGENIELabel(vtxdist, "Vertex - Michel Distance", "mm", 100, 0., 1000., data_error_bands));
  data_studies.push_back(new PerMichelEventVarByGENIELabel(clusdist, "Vertex - Cluster Distance", "mm", 100, 0., 1000., data_error_bands)); 
   #endif //__CINT__


  for(auto& var: vars) var->InitializeMCHists(error_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands;
  auto precuts = reco::GetCCInclusive2DCuts<CVUniverse, MichelEvent>();
  precuts.emplace_back(new Q3RangeReco<CVUniverse, MichelEvent>(0.0, 1.20));
  precuts.emplace_back(new hasMichel<CVUniverse, MichelEvent>());
//  precuts.emplace_back(new BestMichelDistance2D<CVUniverse, MichelEvent>(102.));
  auto signalDefinition = truth::GetCCInclusive2DSignal<CVUniverse>();
  signalDefinition.emplace_back(new Q3Limit<CVUniverse>(1.2));
//  signalDefinition.emplace_back(new hasPion<CVUniverse>);


  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(precuts), std::move(sidebands) , std::move(signalDefinition),std::move(truth::GetCCInclusive2DPhaseSpace<CVUniverse>()));
  // Loop entries and fill
  LoopAndFillEventSelection(chain, error_bands, vars, studies, mycuts);
  //LoopAndFillEffDenom(truth, truth_bands, vars, mycuts);
  std::cout << mycuts << std::endl;
  //mycuts.resetStats();
  //LoopAndFillData(data, data_band,vars, data_studies, mycuts);
  //std::cout << mycuts << std::endl;

  #ifndef __CINT__ //For "auto" c++11 feature because Andrew is lazy.  Added some lambda function action later.
  for(auto& var: vars)
  {
    // You must always sync your HistWrappers before plotting them
    var->SyncCVHistos();

    //Categorized makes sure GetTitle() is the same as the labels you were looping over before
    var->m_bestPionByGENIELabel->visit([](PlotUtils::HistWrapper<CVUniverse>& categ)
                                       {
                                         PlotCVAndError(categ.hist, categ.hist->GetTitle());
                                         PlotErrorSummary(categ.hist, categ.hist->GetTitle());
                                       });
   //var->Write(*outDir); 
  }
  #endif //__CINT__

  //TFile* outDir = TFile::Open("PerMichelHists.root", "RECREATE");
  for(auto& study: studies) study->SaveOrDraw(*outDir);
  for(auto& var: vars) var->Write(*outDir);
  outDir->Write();
  for(auto& study: data_studies) study->SaveOrDraw(*dataDir);  
//  std::cout << mycuts << std::endl; 
  //mycuts.summarizeTruthWithStats(std::cout); 
  std::cout << "Success" << std::endl;
}
/*
void CreateMichels(CVUnivere &univ,std::vector<Michel*> &return_michels){
  unsigned int nmichels = univ->GetNMichels();
  for (unsigned int i = 0; i < nmichels; ++i)
  {
    Michel* current_michel = new Michel(univ, i);
    if (current_michel->is_fitted != 1) continue;
    current_michel->DoesMichelMatchVtx(univ, current_michel);
    current_michel->DoesMichelMatchClus(univ, current_michel);
    return_michels.push_back(current_michel);
  }
  return return_michels;
}

*/
