#include "TreeManager.hh"
#include <iostream>
#include <vector>

using std::cout;
using std::cin;
using std::endl;

using namespace std;

TreeManager::TreeManager(InputManager* aInMgr) {

  InMgr = aInMgr;

  string FileName;
  InMgr->GetVariable("Events",FileName);

  file = new TFile(FileName.c_str(),"RECREATE");
  
  tree_Hits = new TTree("Hits", "Data ROOT Tree");
  tree_AllScatEvents = new TTree("All Scat Events", "Data ROOT Tree");
  tree_BeamParticle = new TTree("Beam Particles", "Data ROOT Tree");  ////////////////////////////
  tree_AllEvents = new TTree("All Events", "Data ROOT Tree");

  tree_Hits->Branch("Events" /*branch name*/, &data4 /*pointer to class instance*/, "E_dep:Det_no");
  tree_AllScatEvents->Branch("Events" /*branch name*/, &data /*pointer to class instance*/, "E_cm:E_scat:E_det:Det_no:phideg:thetadeg:x:y:z");
  tree_AllBeamEvents->Branch("Events" /*branch name*/, &data3 /*pointer to class instance*/, "E_cm:E_beam:E_det:Det_no:phideg:thetadeg:x:y:z"); ///////
  tree_AllEvents->Branch("Events" /*branch name*/, &data2 /*pointer to class instance*/, "N_event");

  Hit = false;

}

TreeManager::~TreeManager() {
  file->Write();
}

void TreeManager::DetEvent() {
  tree_Hits->Fill();
}

void TreeManager::EndOfEvent() {

  // if (Hit == true) {
    // tree_Hits->Fill();//fills detector hits tree
    // }

  tree_AllScatEvents->Fill();//fills all events tree

  tree_AllBeamEvents->Fill();/////////////////////////////////////////

  // Hit = false;

}  

//Lets Tree manager know if there is a hit in detector
void TreeManager::SetHit(bool Ht) {
  Hit = Ht;
}

//functions that fill ntuple
void TreeManager::SetE_cm(double E_cm) {
  data.E_cm = E_cm;
}

void TreeManager::SetE_scat(double E_scat) {
  data.E_scat = E_scat;
}

void TreeManager::SetE_det(double E_det) {
  data.E_det = E_det;
}

void TreeManager::SetDet_no_scat(double Det_no_scat) {
  data.Det_no = Det_no_scat;
}

void TreeManager::SetX(double x) {
  data.x = x;
}

void TreeManager::SetY(double y) {
  data.y = y;
}

void TreeManager::SetZ(double z) {
  data.z = z;
}
//recieve angle in rad, set to deg
void TreeManager::SetTheta(double theta) {
  data.thetadeg = theta*180./3.14159;
}

void TreeManager::SetPhi(double phi) {
  data.phideg = phi*180./3.14159;
}
///////////////////////////////////////////////////////////////////////////
void TreeManager::SetE_beam_cm(double E_beam_cm) {
  data3.E_cm = E_beam_cm;
}

void TreeManager::SetE_beam(double E_beam) {
  data3.E_beam = E_beam;
}

void TreeManager::SetE_beam_det(double E_beam_det) {
  data3.E_det = E_beam_det;
}

void TreeManager::SetDet_no_beam(double Det_no_beam) {
  data3.Det_no = Det_no_beam;
}

void TreeManager::SetX_beam(double x_beam) {
  data3.x = x_beam;
}

void TreeManager::SetY_beam(double y_beam) {
  data3.y = y_beam;
}

void TreeManager::SetZ_beam(double z_beam) {
  data3.z = z_beam;
}

void TreeManager::SetThetaBeam(double theta_beam) {
  data3.thetadeg = theta_beam*180./3.14159;
}

void TreeManager::SetPhiBeam(double phi_beam) {
  data3.phideg = phi_beam*180./3.14159;
}
//////////////////////////////////////////////////////
void TreeManager::SetE_dep(double E_dep) {
  data4.E_dep = E_dep;
}

void TreeManager::SetDet_no(double Det_no) {
  data4.Det_no = Det_no;
}

