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
  file->mkdir("Total_Detector");
  file->mkdir("Recoils_ONLY");
  file->mkdir("Ejectiles_ONLY");

  const char* DetAngles[13] = {"22_5","25","30","35","40","45","55","60","65","75","90","120","135"};

  for(int i=0; i<13; i++) {
    string TotalTree = "angle_";
    TotalTree.append(DetAngles[i]);
    const char * cTotalTree = TotalTree.c_str();
    tree_TotalDetector[i] = new TTree(cTotalTree, "Data ROOT Tree");
    tree_TotalDetector[i]->Branch(DetAngles[i] /*branch name*/, &dataX /*pointer to class instance*/, "E_dep:E_lab:theta_cm:theta_lab:phi_lab:x:y:z:E_beam");

    string RecoilTree = "Recoils_ONLY_angle_";
    RecoilTree.append(DetAngles[i]);
    const char * cRecoilTree = RecoilTree.c_str();
    tree_RecoilDetector[i] = new TTree(cRecoilTree, "Data ROOT Tree");
    tree_RecoilDetector[i]->Branch(DetAngles[i] /*branch name*/, &data0 /*pointer to class instance*/, "E_dep:E_lab:theta_cm:theta_lab:phi_lab:x:y:z:E_beam");

    string EjectileTree = "Ejectiles_ONLY_angle_";
    EjectileTree.append(DetAngles[i]);
    const char * cEjectileTree = EjectileTree.c_str();
    tree_EjectileDetector[i] = new TTree(cEjectileTree, "Data ROOT Tree");
    tree_EjectileDetector[i]->Branch(DetAngles[i] /*branch name*/, &data1 /*pointer to class instance*/, "E_dep:E_lab:theta_cm:theta_lab:phi_lab:x:y:z:E_beam");
  }

}

TreeManager::~TreeManager() {
  for (int i=0; i<13; i++) {
    file->cd("Total_Detector");
    tree_TotalDetector[i]->Write();
    file->cd();
    file->cd("Recoils_ONLY");
    tree_RecoilDetector[i]->Write();
    file->cd();
    file->cd("Ejectiles_ONLY");
    tree_EjectileDetector[i]->Write();
    file->cd();
  }
  file->Close("R");
}

void TreeManager::DetRclEvent(const char* aDetAngle) {
  dataX = data0;
  const char* DetAngles[13] = {"22.5","25","30","35","40","45","55","60","65","75","90","120","135"};
  for (int i=0;i<13;i++) { if (aDetAngle==DetAngles[i]) {
      tree_TotalDetector[i]->Fill();
      tree_RecoilDetector[i]->Fill(); } }
}

void TreeManager::DetEjcEvent(const char* aDetAngle) {
  dataX = data1;
  const char* DetAngles[13] = {"22.5","25","30","35","40","45","55","60","65","75","90","120","135"};
  for (int i=0;i<13;i++) { if (aDetAngle==DetAngles[i]) {
      tree_TotalDetector[i]->Fill();
      tree_EjectileDetector[i]->Fill(); } }
}


//functions that fill ntuple
////////////////////////////////////////////////////////////////////

void TreeManager::Set0E_dep(double E_dep) {
  data0.E_dep = E_dep;
}

void TreeManager::Set0E_lab(double E_lab) {
  data0.E_lab = E_lab;
}

void TreeManager::Set0Tht_cm(double Tht_cm) {
  data0.Tht_cm = Tht_cm*180./3.14159;
}

void TreeManager::Set0Tht_lab(double Tht_lab) {
  data0.Tht_lab = Tht_lab*180./3.14159;
}

void TreeManager::Set0Phi_lab(double Phi_lab) {
  data0.Phi_lab = Phi_lab*180./3.14159;
}

void TreeManager::Set0x(double x) {
  data0.x = x;
}

void TreeManager::Set0y(double y) {
  data0.y = y;
}

void TreeManager::Set0z(double z) {
  data0.z = z;
}

void TreeManager::Set0E_beam(double E_beam) {
  data0.E_beam = E_beam;
}

////////////////////////////////////////////////////////////////////

void TreeManager::Set1E_dep(double E_dep) {
  data1.E_dep = E_dep;
}

void TreeManager::Set1E_lab(double E_lab) {
  data1.E_lab = E_lab;
}

void TreeManager::Set1Tht_cm(double Tht_cm) {
  data1.Tht_cm = Tht_cm*180./3.14159;
}

void TreeManager::Set1Tht_lab(double Tht_lab) {
  data1.Tht_lab = Tht_lab*180./3.14159;
}

void TreeManager::Set1Phi_lab(double Phi_lab) {
  data1.Phi_lab = Phi_lab*180./3.14159;
}

void TreeManager::Set1x(double x) {
  data1.x = x;
}

void TreeManager::Set1y(double y) {
  data1.y = y;
}

void TreeManager::Set1z(double z) {
  data1.z = z;
}

void TreeManager::Set1E_beam(double E_beam) {
  data1.E_beam = E_beam;
}
