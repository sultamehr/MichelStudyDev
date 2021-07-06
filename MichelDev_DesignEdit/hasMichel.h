//#include "CVUniverse.h"
//#include "Michel.h"
//Cutter written to select Michels 
//PlotUtils includes
#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class hasMichel: public PlotUtils::Cut<UNIVERSE, EVENT>
{
 public:
    hasMichel(): PlotUtils::Cut<UNIVERSE, EVENT>("Event Has Michel ")
    {
    }

  private:
    using Michel = typename std::remove_reference<decltype(*std::declval<EVENT>().m_nmichels.front())>::type;

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {
      int nmichels = univ.GetNMichels();
      evt.m_bestdist = 9999.; // setting some default value for best distance
      for (int i = 0; i < nmichels; ++i)
      {
        Michel current_michel = new Michel(univ, i);
        if (current_michel->is_fitted != 1) continue;
         current_michel->DoesMichelMatchVtx(univ, current_michel);
         current_michel->DoesMichelMatchClus(univ, current_michel);

         current_michel->GetBestMatch(current_michel);
         int is_overlay= current_michel->is_overlay;
         if (is_overlay == 1){
   	    std::cout << "THIS IS AN OVERLAY MICHEL" << std::endl;
           // continue;
         }
         
         double dist = current_michel->Best3Ddist; //getting the minimum pion range 
           if (dist <= evt.m_bestdist) {
           evt.m_bestdist = dist;
           evt.m_idx = i;

           evt.m_best_XZ = current_michel->best_XZ;
           evt.m_best_UZ = current_michel->best_UZ;
           evt.m_best_VZ = current_michel->best_VZ;
           evt.m_matchtype = current_michel->BestMatch;
	       int bmatch = current_michel->BestMatch;
           if(bmatch == 1 || bmatch == 3){
           evt.best_x = current_michel->m_x1;
           evt.best_y = current_michel->m_y1;
           evt.best_z = current_michel->m_z1;
           }
           else if(bmatch == 2 || bmatch == 4){
           evt.best_x = current_michel->m_x2;
           evt.best_y = current_michel->m_y2;
           evt.best_z = current_michel->m_z2;
           }		
           else {continue;}
           evt.b_truex = current_michel->true_initialx;
           evt.b_truey = current_michel->true_initialy;
           evt.b_truez = current_michel->true_initialz;
           }

         evt.m_nmichels.push_back(current_michel);
        
       }
       
       


        return !evt.m_nmichels.empty();
    }

};
