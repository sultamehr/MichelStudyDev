#ifndef RECOELECTRON_H
#define RECOELECTRON_H

class recoElectron;

class recoElectron
{
 public: 
     recoElectron(const CVUniverse &univ, int &ci);
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
     int reco_idx = -1; // This is the index of the Michel THIS particular recoElectron is association with. The index would reference each electron in the entire physcis event
}

recoElectron::recoElectron(const CVUniverse &univ, int ci)
{
  trackID = univ.GetVecElem("FittedMichel_reco_micheltrajectory_trackID", ci);
  parentID = univ.GetVecElem("FittedMichel_reco_micheltrajectory_ParentID", ci);
  parentPDG = univ.GetVecElem("FittedMichel_reco_micheltrajectory_ParentPDG", ci);
  pdg = univ.GetVecElem("FittedMichel_reco_micheltrajectory_pdg", ci);
  energy = univ.GetVecElem("FittedMichel_reco_micheltrajectory_energy", ci);
  p = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentum", ci);
  px = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentumx", ci);
  py = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentumy", ci);
  pz = univ.GetVecElem("FittedMichel_reco_micheltrajectory_momentumz", ci);
  fx = univ.GetVecElem("FittedMichel_reco_micheltrajectory_finalx", ci);
  fy = univ.GetVecElem("FittedMichel_reco_micheltrajectory_finaly", ci);
  fz = univ.GetVecElem("FittedMichel_reco_micheltrajectory_finalz", ci);
  reco_idx = univ.GetVecElem("FittedMichel_reco_true_michel_indexvec", ci);
}

#endif // RECOELECTRON_H
