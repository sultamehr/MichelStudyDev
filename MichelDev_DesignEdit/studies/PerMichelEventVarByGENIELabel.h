//studies includes
#include "studies/Study.h"

//Mehreen's includes
#include "MichelEvent.h"
#include "Categorized.h"
#include "CVUniverse.h"

//c++ includes
#include <functional> //for std::function

class PerMichelEventVarByGENIELabel: public Study
{
  public:
    //PerMichelVarByGENIELabel fills a histogram with 1 entry per Michel with some variable calculated from that Michel.  Your function will get to see the CVUniverse, the MichelEvent (= reconstructed Michels), and which Michel it's looping over.
    using reco_t = std::function<double(const CVUniverse&, const MichelEvent&)>;
    
    PerMichelEventVarByGENIELabel(reco_t reco, const std::string& varName, const std::string& varUnits, const int nBins, const double minBin, const double maxBin, std::map<std::string, std::vector<CVUniverse*>>& univs): Study(), fReco(reco), dataHist(new HIST(varName.c_str(), (varName + " [" + varUnits + "]").c_str(), nBins, minBin, maxBin, univs))
    {
       // std::map<int, std::string> GENIELabels = {{1, "QE"},
       //                                                // {8, "2p2h"},
       //                                                                                               // {2, "RES"},
       //                                                                                                                                              // {3, "DIS"}}; //other category is built in for free. 
       std::map<int, std::string> GENIELabels = {{0, "All?"}}; //other category is built in for free.  

       m_VarToGENIELabel = new util::Categorized<HIST, int>(varName, varName + " [" + varUnits + "]", GENIELabels, nBins, minBin, maxBin, univs);
    }
    void SaveOrDraw(TDirectory& outDir)
    {
       outDir.cd();
       m_VarToGENIELabel->visit([](HIST& wrapper)
                                {
                                  wrapper.SyncCVHistos();
                                  wrapper.hist->Write();
                                });
       
       dataHist->SyncCVHistos();
       dataHist->hist->Write();
       //TODO: You could do plotting here
       
    }

  private:
    using HIST = PlotUtils::HistWrapper<CVUniverse>;
    reco_t fReco;
    util::Categorized<HIST, int>* m_VarToGENIELabel;
    HIST* dataHist;
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const MichelEvent& evt, const double weight) {
	dataHist->FillUniverse(&univ, fReco(univ, evt), weight);	
    }

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
        
        (*m_VarToGENIELabel)[univ.GetInteractionType()].FillUniverse(&univ, fReco(univ, evt), weight);
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) {}
};
