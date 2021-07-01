#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class hasRecoTrueMichel: public PlotUtils::Cut<UNIVERSE, EVENT>
{
 public:
    hasRecoTrueMichel(): PlotUtils::Cut<UNIVERSE, EVENT>("Event Has Reco True Michel ")
    {
    }
 private:
    using Michel = typename std::remove_reference<decltype(*std::declval<EVENT>().m_nmichels.front())>::type;

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {

       int ntrue = univ.GetNTrueRecoMichels();
       for (int i = 0; i < ntrue; ++i)
       {
         
 
       }


    }

