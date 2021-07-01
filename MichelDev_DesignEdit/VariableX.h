#ifndef VARIABLEX_H
#define VARIABLEX_H

#ifndef __CINT__ //PlotUtils/VariableBase.h uses std::function which is from c++11
#include "SafeROOTName.cpp" //TODO: This is a very bad idea
#endif //__CINT__
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "Categorized.h"

class VariableX: public PlotUtils::VariableBase<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
    typedef PlotUtils::Hist2DWrapper<CVUniverse> Hist2D;
  public:
    #ifndef __CINT__ //For variadic template parameter ARGS
    template <class ...ARGS>
    VariableX(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...)
    {
    }

    #endif //__CINT__
    //Very dissappointed in Ben for not allowing me to pass universes to the constructor.
    //This is not a solution I like :(
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& error_bands)
    {
      //For example only.  Don't actually use GENIELabels as your backgrounds!
      //You'd end up with a very model-dependent result, and Luke Pickering
      //would frown on your paper ;)
      std::map<int, std::string> GENIELabels = {{1, "QE"},
                                                {8, "2p2h"},
                                                {2, "RES"},
                                                {3, "DIS"}};
 
//      m_bestPionByGENIELabel = new util::Categorized<Hist, int>((GetName() + "_by_GENIE_Label").c_str(),
//                                                                GetName(), GENIELabels,
//                                                                GetNBins(), GetBinVec(), error_bands);
         
      m_eff_num = new Hist((GetName() + "_eff_num").c_str(), (GetName()+"; truth "+GetName()+"; events").c_str(), GetBinVec(), error_bands);
      m_eff_denom = new Hist((GetName() + "_eff_denom").c_str(), (GetName()+"; truth "+GetName()+"; events").c_str(), GetBinVec(), error_bands);
      m_migration = new Hist2D((GetName() + "_migration").c_str(), (GetName()+"; truth "+GetName()+"; events").c_str(), GetBinVec(), GetBinVec(), error_bands);
      //TODO: make m_selected_events_mc later 
 
    }

    #ifndef __CINT__
    //Histots needed for cross sections
    //util::Categorized<Hist, int>* m_bestPionByGENIELabel;
    Hist* m_selected_events; //, m_selection_events_mc TODO: make Mc version need sometimes
    Hist* m_eff_num;
    Hist* m_eff_denom;
    Hist2D* m_migration; 
   // util::Categorized<Hist, int>* m_backgrounds, m_sidebands; 
    // TODO: closure test, m_selectedSignal_events eff numerator with reco variables instead of true 
    #endif //__CINT__

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
	std::vector<double> bins = GetBinVec();
        std::string name = GetName().c_str();
        //dataHist = new Hist(Form("_data_%s", name), name, GetNBins(), bins, universe);
  //	dataHist = new Hist(Form("_data_%s", name), name, GetNBins(), bins, data_error_bands);
        m_selected_events = new Hist(Form("_selected_events_data_reco_%s", name.c_str()), (name + "; reco " + name + "; events").c_str(), GetNBins(), bins, data_error_bands);
    }

    void Write(TFile& file)
    {
      SyncCVHistos();

      file.cd();
      
      m_eff_num->hist->Write();
      m_eff_denom->hist->Write();
      m_migration->hist->Write();
      m_selected_events->hist->Write();
      
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
       m_eff_num->SyncCVHistos();
       m_eff_denom->SyncCVHistos();
       m_migration->SyncCVHistos();
       m_selected_events->SyncCVHistos();   
    }

//Write a destructor
//  ~Variable()
//  {
//    delete dataHist;
//    dataHist = NULL;
//   }
};

#endif //VARIABLEX_H
