#include "PlotUtils/Cut.h"
#include "CVUniverse.h" 
template <class UNIVERSE>
class hasPion: public PlotUtils::SignalConstraint<UNIVERSE>
{
 public:
    hasPion(): PlotUtils::SignalConstraint<UNIVERSE>("Event Has True Pion ")
    {
    }

 private:
    bool checkConstraint(const UNIVERSE& univ) const
    {
 
      int npions = univ.GetNTruePions();
      //evt.m_lowKE = 9999.;
      int nprimpions = 0;
      for (int i = 0; i < npions; ++i)
      {

         int pdg = univ.GetPionPDG(i);
         int parentID = univ.GetPionParentID(i);;
         if (pdg == 211 && parentID == 0) 
         {
           nprimpions +=1; 
         }
      }
    return (nprimpions > 0); 
    }
};
