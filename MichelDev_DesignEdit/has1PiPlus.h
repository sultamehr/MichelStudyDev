//File: has1PiPlus.h

//PlotUtils includes
#include "PlotUtils/Cut.h"
#include "Michel.h"

template <class UNIVERSE, class EVENT>
class has1PiPlus: public PlotUtils::Cut<UNIVERSE, EVENT>
{
  public:
    has1PiPlus(const double maxDistance): PlotUtils::Cut<UNIVERSE, EVENT>("Event has 1 pi+" + std::to_string(maxDistance)), m_maxDistance(maxDistance)
    {
    }

  private:
    double m_maxDistance; //Maximum distance from the vertex that the best Michel can have in mm

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {
       if (evt.eventtype ==4 ) return true;
       else return false;
    }
};
