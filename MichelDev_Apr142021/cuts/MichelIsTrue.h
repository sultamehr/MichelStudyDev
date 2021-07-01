#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class MichelIsTrue: public PlotUtils::Cut<UNIVERSE, EVENT>
{
 public:
    MichelIsTrue(): PlotUtils::Cut<UNIVERSE, EVENT>("Event Has Michel ")
    {
    }

  private:
    using Michel = typename std::remove_reference<decltype(*std::declval<EVENT>().m_nmichels.front())>::type;

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {
      int nmichels = univ.GetNMichels();
      for (int i = 0; i < nmichels; ++i)
      {
        Michel* current_michel = new Michel(univ, i);
        if (current_michel->is_fitted != 1) continue;
         current_michel->DoesMichelMatchVtx(univ, current_michel);
         current_michel->DoesMichelMatchClus(univ, current_michel);

         current_michel->GetBestMatch(current_michel);
         int is_overlay= current_michel->is_overlay;
         if (is_overlay == 1){
            std::cout << "THIS IS AN OVERLAY MICHEL" << std::endl;
            continue;
          }

       }

    }
