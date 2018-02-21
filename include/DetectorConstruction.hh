#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VSensitiveDetector.hh"
#include "TreeManager.hh"
#include "InputManager.hh"
#include "CrossSectionManager.hh"
#include "TrackerSD.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "iostream"

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {

public:
    DetectorConstruction(TreeManager* aTMgr, CrossSectionManager* aCSMgr, InputManager* aInMgr);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();
    void ConstructSONIK();
    void ConstructRegular();
    void SetTargetLength(double len);
    void SetTargetOffset(double offset);//cm

private:
    // Logical volumes
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* mainCyl_log;
    G4LogicalVolume* innerCyl_log;
    G4LogicalVolume* Aperture_log;

    // Physical volumes
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* mainCyl_phys;
    G4VPhysicalVolume* innerCyl_phys;
    G4VPhysicalVolume* Aperture_phys;
    G4VPhysicalVolume* CylCol01_phys;
    G4VPhysicalVolume* CylCol02_phys;
    G4Box* experimentalHall_box;

    G4Material* Target;
    G4Material* Si;
    G4Material* Al;

    string GeomType;

    InputManager* InMgr;
    TreeManager* TMgr;
    CrossSectionManager* CSMgr;

    double TargetLength;
    double TargetOffset;

};

#endif
