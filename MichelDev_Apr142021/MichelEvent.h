#ifndef MichelEvent_h
#define MichelEvent_h

#include "Michel.h"
#include "CVUniverse.h"
#include "MatchedMichel.h"

struct MichelEvent {
    int m_idx = -1; // Index for Best Michel in nmichels
    double m_bestdist = 9999.; // in mm 
    std::vector<double> m_best2D; //0: XZ, 1: UZ, 2:VZ   
    double m_best_XZ = 9999.;
    double m_best_UZ = 9999.;
    double m_best_VZ = 9999.;
    int m_matchtype; // 0 NULL 1 UPVTX 2 DOWNVTX 3 UPCLUS 4 DOWNCLUS
    std::vector<Michel*> m_nmichels; //nmatched michels
    double best_x = 9999.;
    double best_y = 9999.;
    double best_z = 9999.;
    double b_truex = 9999.;
    double b_truey = 9999.;
    double b_truez = 9999.;

};
#endif
