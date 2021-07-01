#ifndef PionEvent_h
#define PionEvent_h

#include "Pion.h"
#include "CVUniverse.h"

struct PionEvent {
        int m_idx = -1; // index of the lowest kinetic energy pion
	double m_lowKE = 9999.; // lowestkinetic  energy pion 
        std::vector<Pion*> m_npions;

};

#endif
