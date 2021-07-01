#ifndef PION_H
#define PION_H

class Pion
{
 public: 
     int trackID = -1;
     int parentID = -1;
     int parentPDG = -1;
     int pdg = -1;
     double energy = -9999.;
     double p = -9999.;
     double px = -9999.;
     double py = -9999.;
     double pz = -9999.;
     double fx = -9999.; // final position
     double fy = -9999.;  
     double fz = -9999.; 
     double ke = -9999.; // kinetic energy
     double mass = -9999.;
     
};


#endif // PION_H
