#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
 
#include "TrackerSD.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "CrossSectionManager.hh"
#include "TreeManager.hh"
#include "InputManager.hh"

#include "ResultsManager.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv) {

  if (argc<2) {// argc should be more than 2 for correct execution
    cerr << "Error in <main>: program " << argv[0] << " requires config filename argument" << endl;
    exit(1);
  }

  std::vector<string> filename;
  std::vector<string>::iterator it;

  for (int i=1; i<argc; i++) {
    filename.push_back(argv[i]);//config input filename(s)
  }

  //Construct custom manager classes
  InputManager* InMgr = new InputManager();//handles config file input

  for (it=filename.begin(); it!=filename.end(); it++) {
    InMgr->ReadFile(it->c_str());
  } 

  TreeManager* TMgr = new TreeManager(InMgr);//handles output root tree
  CrossSectionManager* CSMgr = new CrossSectionManager(InMgr);//handles cross section calc

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  cout << "DetectorConstruction" << endl;
  G4VUserDetectorConstruction* detector = new DetectorConstruction(TMgr,CSMgr,InMgr);
  runManager->SetUserInitialization(detector);
  cout << "PhysicsList" << endl;
  G4VUserPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);

  // set mandatory user action class
  G4VUserPrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction(TMgr,CSMgr,InMgr);
  runManager->SetUserAction(gen_action);

  G4UserRunAction* run_action = new RunAction;
  runManager->SetUserAction(run_action);

  G4UserEventAction* event_action = new EventAction(TMgr,CSMgr);
  runManager->SetUserAction(event_action);
  
  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  #ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
  #endif 

  G4UImanager * UImanager = G4UImanager::GetUIpointer();

//  #ifdef G4UI_USE
  G4UIExecutive * ui = new G4UIExecutive(argc,argv);
//  #ifdef G4VIS_USE

  bool viewer;
  InMgr->GetVariable("viewer",viewer);
  if (viewer == true) UImanager->ApplyCommand("/control/execute vis.mac");     

  char beamOn[30];
  double N_scat = CSMgr->GetN_scat();//gets total number of scattering events
  int N_scatint = int(N_scat);

  cout << "Generating " << N_scatint << " scattering events" << endl << endl;

  sprintf(beamOn,"/run/beamOn %i", N_scatint);

  UImanager->ApplyCommand(beamOn);//generates scattering event N_scat times
  
  delete ui;
  delete visManager;

/*
  #endif
  ui->SessionStart();
  delete ui;
  #endif
     
  #ifdef G4VIS_USE
  delete visManager;
  #endif     
*/
  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //  

  delete runManager;
  delete TMgr;
  delete CSMgr;

  ResultsManager* RMgr = new ResultsManager(InMgr);
  RMgr->CreateTable();
  delete RMgr;  

  delete InMgr;

  return 0;
}


