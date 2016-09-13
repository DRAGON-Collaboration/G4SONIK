#include "PrimaryGeneratorAction.hh"

#include "TRandom.h"

namespace {
TRandom* g_random = new TRandom();
}

PrimaryGeneratorAction::PrimaryGeneratorAction(TreeManager* aTMgr, CrossSectionManager* aCSMgr, InputManager* aInMgr) {
  
  TMgr  = aTMgr;
  CSMgr = aCSMgr;
  InMgr = aInMgr;
  srand(time(0)); //seed random number generator

  InMgr->GetVariable("targetA",Ion_A);
  InMgr->GetVariable("targetZ",Ion_Z);

  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  particleTable = G4ParticleTable::GetParticleTable();

}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete particleGun;
}

//------------------------------------------------------------------------
//Generates single scattering event
//------------------------------------------------------------------------

void PrimaryGeneratorAction::IonBeam(G4Event* anEvent) {

  particleGun->SetParticleDefinition(particleTable->GetIon(Ion_Z,Ion_A,0));

  bool Event = CSMgr->GenerateEvent();//genrates postion/energy/direction of event

  if (Event == true) {//if event occurs within config limits

    double E_scat = CSMgr->GetE();

    particleGun->SetParticleEnergy((E_scat)*MeV);

    double theta = CSMgr->GetTheta();
    double phi   = CSMgr->GetPhi();

    theta = g_random->Gaus(theta,0.001);

    //Determine direction vector of event
    double vx = sin(theta)*sin(phi);
    double vy = sin(theta)*cos(phi);
    double vz = cos(theta);

    double TargetOffset = CSMgr->GetTargetOffset();

    double x = CSMgr->GetX();
    double y = CSMgr->GetY();
    double z = CSMgr->GetZ()+TargetOffset;//start of target is offset from axis origin

    G4ThreeVector vect(vx,vy,vz);
    G4ThreeVector pos(x*cm,y*cm,z*cm);

    //Set variables that fill root tree
    TMgr->SetX(x);
    TMgr->SetY(y);
    TMgr->SetZ(z);
    TMgr->SetTheta(theta);//function recieves angle in rad, sets to deg
    TMgr->SetPhi(phi);
    TMgr->SetE_scat(E_scat);
    TMgr->SetE_det(0.);
    TMgr->SetE_cm(CSMgr->GetEcm());
    TMgr->SetDet_no(-1);

    particleGun->SetParticlePosition(pos);
    particleGun->SetParticleMomentumDirection(vect);

    particleGun->GeneratePrimaryVertex(anEvent);

  }

}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  IonBeam(anEvent);

}
