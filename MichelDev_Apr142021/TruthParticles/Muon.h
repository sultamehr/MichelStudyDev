#ifndef MUON_H
#define MUON_H

class Muon;

class Muon
{
 public: 
     Muon(const CVUniverse &univ, int &i);
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

Muon::Muon(const CVUniverse &univ, int i)
{
  trackID = univ.GetVecElem("FittedMichel_all_muontrajectory_trackID", i);
  parentID = univ.GetVecElem("FittedMichel_all_muontrajectory_ParentID", i);
  parentPDG = univ.GetVecElem("FittedMichel_all_muontrajectory_ParentPDG", i);
  pdg = univ.GetVecElem("FittedMichel_all_muontrajectory_pdg", i);
  energy = univ.GetVecElem("FittedMichel_all_muontrajectory_energy", i);
  p = univ.GetVecElem("FittedMichel_all_muontrajectory_momentum", i);
  px = univ.GetVecElem("FittedMichel_all_muontrajectory_momentumx", i);
  py = univ.GetVecElem("FittedMichel_all_muontrajectory_momentumy", i);
  pz = univ.GetVecElem("FittedMichel_all_muontrajectory_momentumz", i);
  fx = univ.GetVecElem("FittedMichel_all_muontrajectory_finalx", i);
  fy = univ.GetVecElem("FittedMichel_all_muontrajectory_finaly", i);
  fz = univ.GetVecElem("FittedMichel_all_muontrajectory_finalz", i);

}

#endif // MUON_H
