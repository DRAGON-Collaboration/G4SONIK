#include "TrackerSD.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4AttValue.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "TRandom.h"

namespace {
TRandom* g_random = new TRandom();
}

static int n_hit;//number of detector hits

TrackerSD::TrackerSD(G4String name, TreeManager* aTMgr, InputManager* aInMgr, const char* aDetAngle):G4VSensitiveDetector(name) {

  TMgr = aTMgr;
  InMgr = aInMgr;
  DetAngle = aDetAngle;
  G4String HCname;
  collectionName.insert(HCname=DetAngle);
  n_hit=0;

  InMgr->GetVariable("targetA",Ion_A);
  InMgr->GetVariable("targetZ",Ion_Z);

  InMgr->GetVariable("beamA",Ion2_A); //for beam particle
  InMgr->GetVariable("beamZ",Ion2_Z);

  particleTable = G4ParticleTable::GetParticleTable();

}

TrackerSD::~TrackerSD(){ }

void TrackerSD::Initialize(G4HCofThisEvent* ) {

  //trackerCollection = new TrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 

}

G4bool TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*) {

  G4double res, edep;
  InMgr->GetVariable("DRes",res);
  if (res != 0) {
    edep = g_random->Gaus( aStep->GetTotalEnergyDeposit(), res / 2.355); 
  }
  else {
    edep = aStep->GetTotalEnergyDeposit();
  }

  G4double thres;
  InMgr->GetVariable("Thres",thres);

  if(edep<thres) return false;

  TrackerHit* newHit = new TrackerHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());

  //trackerCollection->insert(newHit);
  
  copynum = aStep->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(0);

  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == particleTable->GetIon(Ion_Z,Ion_A,0)->GetParticleName() ) {
    TMgr->Set0E_dep(edep);
    TMgr->DetRclEvent(DetAngle);
  }

  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == particleTable->GetIon(Ion2_Z,Ion2_A,0)->GetParticleName() ) {
    TMgr->Set1E_dep(edep);
    TMgr->DetEjcEvent(DetAngle);
  }

  n_hit += 1;//number of detector hits
  //if (n_hit>1) cout << '\r';
  cout << "\r" << "Number of detected events: " << n_hit << "      " ;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent* HCE) {
/*
  if (verboseLevel>0) {
    G4int NbHits = trackerCollection->entries();
    G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
            << " hits in the tracker chambers: " << G4endl;
    for (G4int i=0;i<NbHits;i++) 
    {
      (*trackerCollection)[i]->Print();
    }
  } 

  static G4int HCID = -1;
  if(HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double TrackerSD::GetCopynum() {
  
  return copynum;
  
}
