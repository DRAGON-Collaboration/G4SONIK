#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CrossSectionManager.hh"
#include "TreeManager.hh"
#include "InputManager.hh"
#include <fstream>
#include <iomanip>
#include <cstdlib>
using std::rand;
using std::srand;

#include <iostream>
using std::cout;
using std::cin;
using std::endl;

#include <ctime> //prototype for time
using std::time;

#include <string>

#include "TRandom.h"

#include <TFile.h>
#include <TFolder.h>
#include <TTree.h>
#include "TH1.h"

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
  PrimaryGeneratorAction(TreeManager* aTMgr, CrossSectionManager* aCSMgr, InputManager* aInMgr);
  ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event* anEvent);
  void IonBeam(G4Event* anEvent);

private:
  string GeomType;
  int EventCounts;

  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;

  TreeManager* TMgr;
  CrossSectionManager* CSMgr;
  InputManager* InMgr;

  G4double Target_Ion_A, Target_Ion_Z, Beam_Ion_A, Beam_Ion_Z;

};

#endif
