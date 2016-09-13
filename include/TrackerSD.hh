
#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"

#include "TrackerHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "CrossSectionManager.hh"
#include "TreeManager.hh"
#include "InputManager.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include <string>

class G4Step;
class G4HCofThisEvent;

class TrackerSD : public G4VSensitiveDetector {

  public:
  TrackerSD(G4String, TreeManager* aTMgr, InputManager* aInMgr, const char* aDetAngle);
  ~TrackerSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  double GetCopynum();
  
  private:
  G4ParticleTable* particleTable;
  G4double Ion_A, Ion_Z, Ion2_A, Ion2_Z;

  //TrackerHitsCollection* trackerCollection;

  double copynum;
  TreeManager* TMgr;
  InputManager* InMgr;
  const char* DetAngle;

};

#endif

