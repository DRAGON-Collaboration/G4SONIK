
#ifndef TreeManager_h
#define TreeManager_h 1

#include <TFile.h>
#include <TFolder.h>
#include <TTree.h>
#include "InputManager.hh"

class TreeManager {

  public:

  TreeManager(InputManager* aInMgr);
 ~TreeManager();

  void EndOfEvent();

  void SetEvent(bool Evnt);
  void SetHit(bool Ht);
  void SetE_cm(double E_cm);
  void SetE_scat(double E_scat);
  void SetE_det(double E_det);
  void SetDet_no(double Det_no);
  void SetX(double x);
  void SetY(double y);
  void SetZ(double z);
  void SetTheta(double theta);
  void SetPhi(double phi);

  private:

  struct Data {
    Float_t E_cm;
    Float_t E_scat;
    Float_t E_det;
    Float_t Det_no;
    Float_t phideg;   // The azimuthal lab angle in degrees
    Float_t thetadeg; // The polar lab angle in degrees
    Float_t x;   // Initial x,y,z coordinates
    Float_t y;
    Float_t z;
  };

  struct Data2 {
    Float_t N_events;
  };

  Int_t N_events;
  TFile* file;
  Data data;
  Data2 data2;
  TTree* tree_Hits;
  TTree* tree_AllScatEvents;
  TTree* tree_AllEvents;
  bool Hit;//detector hit this event
  double N_hit;//number of detected events

  InputManager* InMgr;

};

#endif
