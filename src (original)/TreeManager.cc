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
  tree_AllEvents = new TTree("All Events", "Data ROOT Tree");

  tree_Hits->Branch("Events" /*branch name*/, &data /*pointer to class instance*/, "E_cm:E_scat:E_det:Det_no:phideg:thetadeg:x:y:z");
  tree_AllScatEvents->Branch("Events" /*branch name*/, &data /*pointer to class instance*/, "E_cm:E_scat:E_det:Det_no:phideg:thetadeg:x:y:z");
  tree_AllEvents->Branch("Events" /*branch name*/, &data2 /*pointer to class instance*/, "N_event");

  Hit = false;

}

TreeManager::~TreeManager() {
  file->Write();
}

void TreeManager::EndOfEvent() {

  if (Hit == true) {
    tree_Hits->Fill();//fills detector hits tree
  }

  tree_AllScatEvents->Fill();//fills all events tree

  Hit = false;

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

void TreeManager::SetDet_no(double Det_no) {
  data.Det_no = Det_no;
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
