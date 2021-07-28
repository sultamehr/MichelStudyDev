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
#include "has1PiPlus.h"
#endif //CINT
#include "Michel.h"
#include "Cluster.h"
#include "MichelEvent.h"

//The following are in development and might not be needed
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
    std::vector<Variable*> vars, //Set of Event level Variables to Fill
    std::vector<Study*> studies, // Michel level variables currently built as studies
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
          
          //For each Variable
          //if(best_michel != -1) //Equivalent to if (return_michels.size() == 0) continue;
        
          for(auto& var: vars)
          {
            //TODO: #include "study.h"

            //#include "MichelEventStudy.h"	
            (*var->m_bestPionByGENIELabel)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), universe->GetWeight());

            //Fill other per-Variable histograms here
          }
          for(auto& study: studies) study->SelectedSignal(*universe, myevent, 1); //TODO: Last argument is weight
          
      } // End band's universe loop
    } // End Band loop

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
void pionEventLoop() {

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
  TFile* outDir = TFile::Open("July16_truthcomapre_MCHists.root", "RECREATE");

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
  
    std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   int evttype = evt.eventtype;
                                   double micheltime = evt.m_nmichels[whichMichel].time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
                                   if (evttype == 1) return deltat;
    			           else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   int evttype = evt.eventtype;
                                   double micheltime = evt.m_nmichels[whichMichel].time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
                                   if (evttype == 2) return deltat;
                                   else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   int evttype = evt.eventtype;
                                   double micheltime = evt.m_nmichels[whichMichel].time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
                                   if (evttype == 3) return deltat;
                                   else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   int evttype = evt.eventtype;
                                   double micheltime = evt.m_nmichels[whichMichel].time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
                                   if (evttype == 4) return deltat;
                                   else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   int evttype = evt.eventtype;
                                   double micheltime = evt.m_nmichels[whichMichel].time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
                                   if (evttype == 5) return deltat;
                                   else return -9999.;
                                 };


// Michel - Vertex XZ Distance
 
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_XZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_XZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { 
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_XZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { 
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_XZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { 
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_XZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { 
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_XZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };


// Michel - Vertex UZ Distance


  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_UZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_UZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_UZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_UZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_UZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_UZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

//Michel - Vertex VZ


  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_VZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_VZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_VZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_VZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_up_VZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_VZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };


// Michel Endpoint2 to Vertex XZ
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_XZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_XZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_XZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_XZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_XZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_XZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint2 to Vertex UZ

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_UZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_UZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_UZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_UZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_UZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_UZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 2 to Vertex VZ

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_VZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_VZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_VZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_VZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> vtx_down_VZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_VZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 1 - Vertex 3D distance
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_dist_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_dist3D;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
				};

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_dist_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_dist3D;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_dist_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_dist3D;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_dist_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_dist3D;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_vtx_dist_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_vertex_dist3D;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

//Michel Endpoint 2 - Vertex 3D Distance
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_dist_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_dist3D;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_dist_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_dist3D;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_dist_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_dist3D;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_dist_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { 
				  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_dist3D;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_vtx_dist_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
 				  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_vertex_dist3D;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 1 To Cluster XZ

std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_XZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
				  int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
                                  if (evttype == 1) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_XZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
                                  if (evttype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_XZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
                                  if (evttype == 3) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_XZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
                                  if (evttype == 4) return micheldist;
                                  else return -9999.;
                                };

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_XZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                { int evttype = evt.eventtype;
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_XZ;
                                  if (evttype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 1 to Cluster UZ

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_UZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_UZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_UZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_UZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_UZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_UZ;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 1 to Cluster VZ

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_VZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_VZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_VZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_VZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_up_VZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_clus_VZ;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 2 to Cluster XZ

std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_XZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_XZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_XZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_XZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_XZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_XZ;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };

// Michel Endpoint 2 to Cluster UZ

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_UZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_UZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_UZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_UZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_UZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_UZ;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };



//Michel Endpoint 2 to Cluster VZ

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_VZ_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_VZ_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_VZ_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_VZ_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> clus_down_VZ_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_clus_VZ;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };


// Michel Endpoint 1 to Cluster 3D dist

std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusmich_dist_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_clus_michel_dist3D;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusmich_dist_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_clus_michel_dist3D;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
                                
std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusmich_dist_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_clus_michel_dist3D;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusmich_dist_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_clus_michel_dist3D;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
                                
                                
std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusmich_dist_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_clus_michel_dist3D;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };
                                                                
// Michel Endpoint 2 to Cluster 3D distance

std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusmich_dist_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_clus_michel_dist3D;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusmich_dist_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_clus_michel_dist3D;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusmich_dist_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_clus_michel_dist3D;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusmich_dist_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_clus_michel_dist3D;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusmich_dist_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_clus_michel_dist3D;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };
 
                                
// Clusters matched to Michel Endpoint 1 - Vertex 3D distance

std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusvtx_dist_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_cluster_dist3D;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusvtx_dist_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_cluster_dist3D;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusvtx_dist_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_cluster_dist3D;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusvtx_dist_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_cluster_dist3D;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> up_clusvtx_dist_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].up_to_cluster_dist3D;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };

// Clusters Matched to Michel Endpoint 2 - Vertex 3D distance

std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusvtx_dist_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_cluster_dist3D;
                                  if (evt.eventtype == 1) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusvtx_dist_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_cluster_dist3D;
                                  if (evt.eventtype == 2) return micheldist;
                                  else return -9999.;
                                };
std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusvtx_dist_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_cluster_dist3D;
                                  if (evt.eventtype == 3) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusvtx_dist_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_cluster_dist3D;
                                  if (evt.eventtype == 4) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> down_clusvtx_dist_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                {
                                  double micheldist = evt.m_nmichels[whichMichel].down_to_cluster_dist3D;
                                  if (evt.eventtype == 5) return micheldist;
                                  else return -9999.;
                                };

std::function<double(const CVUniverse&, const MichelEvent&, const int)> bestmichel_range_1pi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
                                   if (evt.eventtype == 1) return micheldist;
                                   else return -9999.;
                                 };

 std::function<double(const CVUniverse&, const MichelEvent&, const int)> bestmichel_range_npi = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
                                   if (evt.eventtype == 2) return micheldist;
                                   else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> bestmichel_range_pi0 = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
                                   if (evt.eventtype == 3) return micheldist;
                                   else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> bestmichel_range_kaons = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
                                   if (evt.eventtype == 4) return micheldist;
                                   else return -9999.;
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&, const int)> bestmichel_range_other = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   double micheldist = evt.m_nmichels[whichMichel].Best3Ddist;
                                   if (evt.eventtype == 5) return micheldist;
                                   else return -9999.;
                                 };                                 
                                                                 
// Comparing True Michel Initial Position to Best Matched position

  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff1_1pi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 { 
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x1;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
				   double diff = tx -x;
                                   if (evt.eventtype == 1) return diff;
                                   else return -9999.;
                                   }
                                 };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff1_npi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x1;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 2) return diff;
                                   else return -9999.;
                                   }
                                 };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff1_pi0 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x1;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 3) return diff;
                                   else return -9999.;
                                   }
                                 };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff1_kaons = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x1;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 4) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff1_other = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x1;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 5) return diff;
                                   else return -9999.;
                                   }
                                 };

// U Diff (true - reco)

 std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff1_1pi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y1;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 1) return diff;
                                   else return -9999.;
                                   }
                                 };
std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff1_npi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y1;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 2) return diff;
                                   else return -9999.;
                                   }
                                 };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff1_pi0 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y1;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 3) return diff;
                                   else return -9999.;
                                   }
                                 };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff1_kaons = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y1;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 4) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff1_other = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y1;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 5) return diff;
                                   else return -9999.;
                                   }
                                 };

// Z diff

std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff1_1pi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z1;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 1) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff1_npi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z1;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 2) return diff;
                                   else return -9999.;
                                   }
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff1_pi0 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z1;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 3) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff1_kaons = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z1;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 4) return diff;
                                   else return -9999.;
                                   }
                                 };
 std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff1_other = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z1;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 5) return diff;
                                   else return -9999.;
                                   }
                                 };

// X diff 2

  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff2_1pi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x2;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 1) return diff;
                                   else return -9999.;
                                   }
                                 };
std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff2_npi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x2;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 2) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff2_pi0 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x2;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 3) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff2_kaons = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x2;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 4) return diff;
                                   else return -9999.;
                                   }
                                 };
std::function<double(const CVUniverse&, const MichelEvent&)> bestxdiff2_other = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double x = evt.m_nmichels[bestidx].m_x2;
                                   double tx = evt.m_nmichels[bestidx].true_initialx;
                                   double diff = tx -x;
                                   if (evt.eventtype == 5) return diff;
                                   else return -9999.;
                                   }
                                 };
// Y diff 2

  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff2_1pi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y2;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 1) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff2_npi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y2;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 2) return diff;
                                   else return -9999.;
                                   }
                                 };

  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff2_pi0 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y2;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 3) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff2_kaons = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y2;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 4) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestydiff2_other = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double y = evt.m_nmichels[bestidx].m_y2;
                                   double ty = evt.m_nmichels[bestidx].true_initialy;
                                   double diff = ty -y;
                                   if (evt.eventtype == 5) return diff;
                                   else return -9999.;
                                   }
                                 };

//Z diff 2
 std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff2_1pi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z2;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 1) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff2_npi = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z2;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 2) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff2_pi0 = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z2;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 3) return diff;
                                   else return -9999.;
                                   }
                                 };  
  std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff2_kaons = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z2;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 4) return diff;
                                   else return -9999.;
                                   }
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&)> bestzdiff2_other = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   int bestidx = evt.m_idx;
                                   if (bestidx < 0) return -9999.;
                                   else{
                                   double z = evt.m_nmichels[bestidx].m_z2;
                                   double tz = evt.m_nmichels[bestidx].true_initialz;
                                   double diff = tz -z;
                                   if (evt.eventtype == 5) return diff;
                                   else return -9999.;
                                   }
                                 };

  // Fill Studies here
  //
  //
  // Michel Time difference with primary muon
 
  studies.push_back(new PerMichelVarByGENIELabel(delta_t_1pi, "micheltime_1pi", "#mus", 30, 0.0, 9.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(delta_t_npi, "micheltime_npi", "#mus", 30, 0.0, 9.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(delta_t_pi0, "micheltime_pi0", "#mus", 30, 0.0, 9.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(delta_t_kaons, "micheltime_kaons", "#mus", 30, 0.0, 9.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(delta_t_other, "micheltime_other", "#mus", 30, 0.0, 9.0, error_bands));

  // Michel XZ distance with Vertec
  
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_XZ_1pi, "michelvtx_up_XZ_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_XZ_npi, "michelvtx_up_XZ_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_XZ_pi0, "michelvtx_up_XZ_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_XZ_kaons, "michelvtx_up_XZ_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_XZ_other, "michelvtx_up_XZ_other", "mm", 100, 0.0, 1000.0, error_bands));

  // Michel UZ distance with Vertex 
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_UZ_1pi, "michelvtx_up_UZ_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_UZ_npi, "michelvtx_up_UZ_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_UZ_pi0, "michelvtx_up_UZ_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_UZ_kaons, "michelvtx_up_UZ_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_UZ_other, "michelvtx_up_UZ_other", "mm", 100, 0.0, 1000.0, error_bands));
  
  // Michel VZ distance with Vertex
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_VZ_1pi, "michelvtx_up_VZ_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_VZ_npi, "michelvtx_up_VZ_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_VZ_pi0, "michelvtx_up_VZ_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_VZ_kaons, "michelvtx_up_VZ_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_up_VZ_other, "michelvtx_up_VZ_other", "mm", 100, 0.0, 1000.0, error_bands));
  
  // Michel Endpoint 2 XZ distance to Vertex
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_XZ_1pi, "michelvtx_down_XZ_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_XZ_npi, "michelvtx_down_XZ_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_XZ_pi0, "michelvtx_down_XZ_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_XZ_kaons, "michelvtx_down_XZ_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_XZ_other, "michelvtx_down_XZ_other", "mm", 100, 0.0, 1000.0, error_bands));

  // Michel Endpoint 2 UZ distance to Vertex
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_UZ_1pi,   "michelvtx_down_UZ_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_UZ_npi,   "michelvtx_down_UZ_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_UZ_pi0,   "michelvtx_down_UZ_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_UZ_kaons, "michelvtx_down_UZ_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_UZ_other, "michelvtx_down_UZ_other", "mm", 100, 0.0, 1000.0, error_bands));

  // Michel Endpoint 2 VZ distance to Vertex
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_VZ_1pi,   "michelvtx_down_VZ_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_VZ_npi,   "michelvtx_down_VZ_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_VZ_pi0,   "michelvtx_down_VZ_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_VZ_kaons, "michelvtx_down_VZ_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(vtx_down_VZ_other, "michelvtx_down_VZ_other", "mm", 100, 0.0, 1000.0, error_bands));

  // Michel Endpoint 1 to Vertex 3D distance
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_dist_1pi, "michelvtx_up_3Ddist_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_dist_npi, "michelvtx_up_3Ddist_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_dist_pi0, "michelvtx_up_3Ddist_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_dist_kaons, "michelvtx_up_3Ddist_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_vtx_dist_other, "michelvtx_up_3Ddist_other", "mm", 100, 0.0, 1000.0, error_bands));

  // Michel Endpoint 2 to Vertex 3D Distance
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_dist_1pi, "michelvtx_down_3Ddist_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_dist_npi, "michelvtx_down_3Ddist_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_dist_pi0, "michelvtx_down_3Ddist_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_dist_kaons, "michelvtx_down_3Ddist_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_vtx_dist_other, "michelvtx_down_3Ddist_other", "mm", 100, 0.0, 1000.0, error_bands));

 // Michel Endpoint 1 to Cluster XZ
 
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_XZ_1pi, "michelclus_XZ_1pi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_XZ_npi, "michelclus_XZ_npi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_XZ_pi0, "michelclus_XZ_pi0", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_XZ_kaons, "michelclus_XZ_kaons", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_XZ_other, "michelclus_XZ_other", "mm", 100, 0.0, 1000., error_bands));

// Michel Endpoint 1 to Cluster UZ

  studies.push_back(new PerMichelVarByGENIELabel(clus_up_UZ_1pi, "michelclus_UZ_1pi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_UZ_npi, "michelclus_UZ_npi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_UZ_pi0, "michelclus_UZ_pi0", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_UZ_kaons, "michelclus_UZ_kaons", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_UZ_other, "michelclus_UZ_other", "mm", 100, 0.0, 1000., error_bands));

// Michel Endpoint1 to Cluster VZ

  studies.push_back(new PerMichelVarByGENIELabel(clus_up_VZ_1pi,   "michelclus_VZ_1pi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_VZ_npi,   "michelclus_VZ_npi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_VZ_pi0,   "michelclus_VZ_pi0", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_VZ_kaons, "michelclus_VZ_kaons", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_up_VZ_other, "michelclus_VZ_other", "mm", 100, 0.0, 1000., error_bands));

//Michel Endpoint 2 to Cluster XZ

  studies.push_back(new PerMichelVarByGENIELabel(clus_down_XZ_1pi,   "michelclus2_XZ_1pi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_XZ_npi,   "michelclus2_XZ_npi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_XZ_pi0,   "michelclus2_XZ_pi0", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_XZ_kaons, "michelclus2_XZ_kaons", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_XZ_other, "michelclus2_XZ_other", "mm", 100, 0.0, 1000., error_bands));
  
// Michel Endpoint 2 to Cluster UZ

  studies.push_back(new PerMichelVarByGENIELabel(clus_down_UZ_1pi,   "michelclus2_UZ_1pi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_UZ_npi,   "michelclus2_UZ_npi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_UZ_pi0,   "michelclus2_UZ_pi0", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_UZ_kaons, "michelclus2_UZ_kaons", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_UZ_other, "michelclus2_UZ_other", "mm", 100, 0.0, 1000., error_bands));


// Michel Endpoint 2 to Cluster VZ

  studies.push_back(new PerMichelVarByGENIELabel(clus_down_VZ_1pi,   "michelclus2_VZ_1pi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_VZ_npi,   "michelclus2_VZ_npi", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_VZ_pi0,   "michelclus2_VZ_pi0", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_VZ_kaons, "michelclus2_VZ_kaons", "mm", 100, 0.0, 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(clus_down_VZ_other, "michelclus2_VZ_other", "mm", 100, 0.0, 1000., error_bands));

// Michel Endpiont 1 to Cluster 3D dist

  studies.push_back(new PerMichelVarByGENIELabel(up_clusmich_dist_1pi,   "michelclus1_3Ddist_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clusmich_dist_npi,   "michelclus1_3Ddist_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clusmich_dist_pi0,   "michelclus1_3Ddist_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clusmich_dist_kaons, "michelclus1_3Ddist_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(up_clusmich_dist_other, "michelclus1_3Ddist_other", "mm", 100, 0.0, 1000.0, error_bands));


// Michel Endpoint 2 to Cluster 3D Dist

studies.push_back(new PerMichelVarByGENIELabel(down_clusmich_dist_1pi,   "michelclus2_3Ddist_1pi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clusmich_dist_npi,   "michelclus2_3Ddist_npi", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clusmich_dist_pi0,   "michelclus2_3Ddist_pi0", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clusmich_dist_kaons, "michelclus2_3Ddist_kaons", "mm", 100, 0.0, 1000.0, error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(down_clusmich_dist_other, "michelclus2_3Ddist_other", "mm", 100, 0.0, 1000.0, error_bands));

// clusters matched to endpoint 1 - vertex 3D distance

 studies.push_back(new PerMichelVarByGENIELabel(up_clusvtx_dist_1pi,   "up_clusvtx_dist_1pi", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(up_clusvtx_dist_npi,   "up_clusvtx_dist_npi", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(up_clusvtx_dist_pi0,   "up_clusvtx_dist_pi0", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(up_clusvtx_dist_kaons,   "up_clusvtx_dist_kaons", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(up_clusvtx_dist_other,   "up_clusvtx_dist_other", "mm", 100, 0.0, 1000.0, error_bands));

// Clusters matched to endpoint 2 - Vertex 3D distance

 studies.push_back(new PerMichelVarByGENIELabel(down_clusvtx_dist_1pi,     "down_clusvtx_dist_1pi", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(down_clusvtx_dist_npi,     "down_clusvtx_dist_npi", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(down_clusvtx_dist_pi0,     "down_clusvtx_dist_pi0", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(down_clusvtx_dist_kaons,   "down_clusvtx_dist_kaons", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(down_clusvtx_dist_other,   "down_clusvtx_dist_other", "mm", 100, 0.0, 1000.0, error_bands));

// Pion Range for Best Michel

 studies.push_back(new PerMichelVarByGENIELabel(bestmichel_range_1pi,     "bestmichel_range_1pi", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(bestmichel_range_npi,     "bestmichel_range_npi", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(bestmichel_range_pi0,     "bestmichel_range_pi0", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(bestmichel_range_kaons,   "bestmichel_range_kaons", "mm", 100, 0.0, 1000.0, error_bands));
 studies.push_back(new PerMichelVarByGENIELabel(bestmichel_range_other,   "bestmichel_range_other", "mm", 100, 0.0, 1000.0, error_bands));
// Get Dfiference between True - reco X for thebest Michel in the event.

studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff1_1pi  , "bestxdiff1_1pi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff1_npi  , "bestxdiff1_npi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff1_pi0  , "bestxdiff1_pi0", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff1_kaons, "bestxdiff1_kaons", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff1_other, "bestxdiff1_other", "mm", 50, -500., 500., error_bands));


 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff1_1pi  , "bestydiff1_1pi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff1_npi  , "bestydiff1_npi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff1_pi0  , "bestydiff1_pi0", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff1_kaons, "bestydiff1_kaons", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff1_other, "bestydiff1_other", "mm", 50, -500., 500., error_bands));

 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff1_1pi  , "bestzdiff1_1pi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff1_npi  , "bestzdiff1_npi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff1_pi0  , "bestzdiff1_pi0", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff1_kaons, "bestzdiff1_kaons", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff1_other, "bestzdiff1_other", "mm", 50, -500., 500., error_bands));
 
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff2_1pi  , "bestxdiff2_1pi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff2_npi  , "bestxdiff2_npi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff2_pi0  , "bestxdiff2_pi0", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff2_kaons, "bestxdiff2_kaons", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestxdiff2_other, "bestxdiff2_other", "mm", 50, -500., 500., error_bands));


 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff2_1pi  , "bestydiff2_1pi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff2_npi  , "bestydiff2_npi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff2_pi0  , "bestydiff2_pi0", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff2_kaons, "bestydiff2_kaons", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestydiff2_other, "bestydiff2_other", "mm", 50, -500., 500., error_bands));

 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff2_1pi  , "bestzdiff2_1pi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff2_npi  , "bestzdiff2_npi", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff2_pi0  , "bestzdiff2_pi0", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff2_kaons, "bestzdiff2_kaons", "mm", 50, -500., 500., error_bands));
 studies.push_back(new PerMichelEventVarByGENIELabel(bestzdiff2_other, "bestzdiff2_other", "mm", 50, -500., 500., error_bands));


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
   
  //data_studies.push_back(new PerMichelVarByGENIELabel(oneVar, "Best Distance DATA", "mm", 100, 0., 1000., data_error_bands));
  //data_studies.push_back(new PerMichelVarByGENIELabel(michelE, "Michel energy DATA", "MeV", 20, 0., 100., data_error_bands));
  //data_studies.push_back(new PerMichelVarByGENIELabel(delta_t, "Michel Time Diff DATA", "#mus", 30, 0.0, 9.0, data_error_bands));
  //data_studies.push_back(new PerMichelEventVarByGENIELabel(bestdistvar, "Best Distance DATA", "mm", 100, 0., 1000., data_error_bands));
  //data_studies.push_back(new PerMichelEventVarByGENIELabel(vtxdist, "Vertex - Michel Distance", "mm", 100, 0., 1000., data_error_bands));
  //data_studies.push_back(new PerMichelEventVarByGENIELabel(clusdist, "Vertex - Cluster Distance", "mm", 100, 0., 1000., data_error_bands)); 
   #endif //__CINT__


  for(auto& var: vars) var->InitializeMCHists(error_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands;
  auto precuts = reco::GetCCInclusive2DCuts<CVUniverse, MichelEvent>();
  precuts.emplace_back(new Q3RangeReco<CVUniverse, MichelEvent>(0.0, 1.20));
  precuts.emplace_back(new hasMichel<CVUniverse, MichelEvent>());
  precuts.emplace_back(new BestMichelDistance2D<CVUniverse, MichelEvent>(102.));
//  precuts.emplace_back(new has1PiPlus<CVUniverse, MichelEvent>(0.0));
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
