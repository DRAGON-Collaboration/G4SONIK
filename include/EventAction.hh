#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "TreeManager.hh"
#include "CrossSectionManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4VVisManager.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4AttValue.hh"

class G4Event;

class EventAction : public G4UserEventAction {

public:
  EventAction(TreeManager* aTMgr, CrossSectionManager* aCSMgr);
  ~EventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

private:
  TreeManager* TMgr;
  CrossSectionManager* CSMgr;

};

#endif
