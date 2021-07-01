#ifndef OTHER_H
#define OTHER_H

class Other;

class Other
{
 public: 
     Other(const CVUniverse &univ, int &i);
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

Other::Other(const CVUniverse &univ, int i)
{
  trackID = univ.GetVecElem("FittedMichel_all_othertrajectory_trackID", i);
  parentID = univ.GetVecElem("FittedMichel_all_othertrajectory_ParentID", i);
  parentPDG = univ.GetVecElem("FittedMichel_all_othertrajectory_ParentPDG", i);
  pdg = univ.GetVecElem("FittedMichel_all_othertrajectory_pdg", i);
  energy = univ.GetVecElem("FittedMichel_all_othertrajectory_energy", i);
  p = univ.GetVecElem("FittedMichel_all_othertrajectory_momentum", i);
  px = univ.GetVecElem("FittedMichel_all_othertrajectory_momentumx", i);
  py = univ.GetVecElem("FittedMichel_all_othertrajectory_momentumy", i);
  pz = univ.GetVecElem("FittedMichel_all_othertrajectory_momentumz", i);
  fx = univ.GetVecElem("FittedMichel_all_othertrajectory_finalx", i);
  fy = univ.GetVecElem("FittedMichel_all_othertrajectory_finaly", i);
  fz = univ.GetVecElem("FittedMichel_all_othertrajectory_finalz", i);

}

#endif // Other_H
