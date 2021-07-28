//File: BestMichelDistance2D.h
//Brief: Require the the distance to vertex for the best Michel electron is below some value.
//Author: Andrew Olivier aolivier@ur.rochester.edu, Mehreen Sultana msultana@ur.rochester.edu

//PlotUtils includes
#include "PlotUtils/Cut.h"
#include "Michel.h"

template <class UNIVERSE, class EVENT>
class BestMichelDistance2D: public PlotUtils::Cut<UNIVERSE, EVENT>
{
  public:
    BestMichelDistance2D(const double maxDistance): PlotUtils::Cut<UNIVERSE, EVENT>("Per Michel 2D Distance in at least two views is < " + std::to_string(maxDistance) + "mm"), m_maxDistance(maxDistance)
    {
    }

  private:
    double m_maxDistance; //Maximum distance from the vertex that the best Michel can have in mm

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {

      std::vector<Michel> nmichelspass;
      for (unsigned int i = 0; i < evt.m_nmichels.size(); i++)
      {
        double XZ = evt.m_nmichels[i].best_XZ; //this values is the 2D distance for the closest type of match to the michel. (can be either endpoint 1-vertex or endpoint1-cluster or endpoint2 - vertex or endpoint2 - cluster) 
        double UZ = evt.m_nmichels[i].best_UZ;//this values is the 2D distance for the closest type of match to the michel. (can be either endpoint 1-vertex or endpoint1-cluster or endpoint2 - vertex or endpoint2 - cluster) 
        double VZ = evt.m_nmichels[i].best_VZ;//this values is the 2D distance for the closest type of match to the michel. (can be either endpoint 1-vertex or endpoint1-cluster or endpoint2 - vertex or endpoint2 - cluster) 
        bool pass = false;
        if (XZ < m_maxDistance && (UZ < m_maxDistance || VZ < m_maxDistance)) pass = true;
        else if (UZ < m_maxDistance && (XZ < m_maxDistance || VZ < m_maxDistance)) pass = true;
        else if (VZ < m_maxDistance && (XZ < m_maxDistance || UZ < m_maxDistance)) pass = true;
        else pass = false;
        
        if (pass == true) nmichelspass.push_back(evt.m_nmichels[i]);
            
      }
      evt.m_nmichels.clear(); //empty existing vector of Michels
      evt.m_nmichels = nmichelspass; // replace vector of michels with the vector of michels that passed the above cut
      return !evt.m_nmichels.empty();
    }
};
