#include "EventAction.hh"

EventAction::EventAction(TreeManager* aTMgr, CrossSectionManager* aCSMgr) {

  TMgr = aTMgr;  
  CSMgr = aCSMgr;

}

EventAction::~EventAction() {

}

void EventAction::BeginOfEventAction(const G4Event*) {

}

void EventAction::EndOfEventAction(const G4Event* evt) {

  G4int event_id = evt->GetEventID();
  
//  G4cout << "end of event" << G4endl;
  
  if(CSMgr->GetEvent() == true) TMgr->EndOfEvent();
  
  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
}
