#ifndef RECOELECTRON_H
#define RECOELECTRON_H

class recoElectron;

class recoElectron
{
 public: 
     recoElectron(const CVUniverse &univ, int &i);
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

}

recoElectron::recoElectron(const CVUniverse &univ, int i)
{
  trackID = univ.GetVecElem("FittedMichel_reco_micheltrajectory_trackID", i);
  parentID = univ.GetVecElem("FittedMichel_reco_micheltrajectory_ParentID", i);
  parentPDG = univ.GetVecElem("FittedMichel_reco_micheltrajectory_ParentPDG", i);
  pdg = univ.GetVecElem("FittedMichel_reco_micheltrajectory_pdg", i);
  energy = univ.GetVecElem("FittedMichel_reco_micheltrajectory_energy", i);
  p = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentum", i);
  px = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentumx", i);
  py = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentumy", i);
  pz = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentumz", i);
  fx = univ.GetVecElem("FittedMichel_reco_micheltrajectory_finalx", i);
  fy = univ.GetVecElem("FittedMichel_reco_micheltrajectory_finaly", i);
  fz = univ.GetVecElem("FittedMichel_reco_micheltrajectory_finalz", i);

}

#endif // RECOELECTRON_H
