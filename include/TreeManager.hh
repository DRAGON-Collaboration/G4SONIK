#ifndef TreeManager_h
#define TreeManager_h 1

#include <TFile.h>
#include <TFolder.h>
#include <TTree.h>

#include "InputManager.hh"

#include <string>

class TreeManager {

public:
  TreeManager(InputManager* aInMgr);
  ~TreeManager();

  void DetRclEvent(const char* aDetAngle);
  void DetEjcEvent(const char* aDetAngle);

  void Set0E_dep(double E_dep);
  void Set0E_lab(double E_lab);
  void Set0Tht_cm(double Tht_cm);
  void Set0Tht_lab(double Tht_lab);
  void Set0Phi_lab(double Phi_lab);
  void Set0x(double x);
  void Set0y(double y);
  void Set0z(double z);
  void Set0E_beam(double E_beam);

  void Set1E_dep(double E_dep);
  void Set1E_lab(double E_lab);

  void Set1Tht_cm(double Tht_cm);
  void Set1Tht_lab(double Tht_lab);
  void Set1Phi_lab(double Phi_lab);
  void Set1x(double x);
  void Set1y(double y);
  void Set1z(double z);
  void Set1E_beam(double E_beam);


private:
  struct DataD {
    Float_t E_dep;
    Float_t E_lab;
    Float_t Tht_cm;
    Float_t Tht_lab;
    Float_t Phi_lab;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t E_beam;
  };

  TFile* file;

  DataD dataX;
  DataD data0;
  DataD data1;

  TTree* tree_TotalDetector[13];
  TTree* tree_RecoilDetector[13];
  TTree* tree_EjectileDetector[13];

  InputManager* InMgr;

};

#endif
