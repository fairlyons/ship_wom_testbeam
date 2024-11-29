//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef OpNoviceDetectorConstruction_h
#define OpNoviceDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Cache.hh"

#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4GenericTrap.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

// #include "LXePMTSD.hh"
// #include "G4GDMLParser.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceDetectorMessenger;

class OpNoviceDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    OpNoviceDetectorConstruction();
    void DefineMaterials();
    void DefineMPTs();
    void DefineSurfaces();
    void DefineSolids();
    void DefineLogicalVolumes();
    void DefineVisAttributes();
    void ConstructVolumes();
    virtual G4VPhysicalVolume* Construct();
    virtual ~OpNoviceDetectorConstruction();
    G4bool intersect_check;

    double ytl;
    double ytr;
    double ybl;
    double ybr;
    std::vector<std::pair<G4double,G4double>>  WOM_coord_vec;

  private:
    G4double fExpHall;

    G4Material* steel;
    G4Material* Al;
    G4Material* Silicon;
    G4Material* ResinSi;
    G4Material* air;
    G4Material* LAB_PPO;
    G4Material* Bis_MSB;
    G4Material* PEMA;
    G4Material* PTP;
    G4Material* WLS_Coat;
    G4Material* PMMA_side;
    G4Material* PMMA_bottom;
    G4Material* BaSO4;

    // geometrical parameters
    G4double SteelZ;
    G4double WallThick;
    G4double Thickness_Reflect;
    G4double ReflectZ;
    G4double SctZ;
    
    G4double Diam_In_In;
    G4double Diam_In_Out;
    G4double Diam_Out_In;
    G4double Diam_Out_Out;
    G4double Diam_WOM_In;
    G4double Diam_WOM_Out;
    G4double Diam_Steel_Add;
    G4double Diam_Hat;
    G4double Length_Out;
    G4double Length_In;
    G4double Thickness_Ring;
    G4double Thickness_Disk;
    G4double Thickness_Hat;
    G4double Thickness_Steel_Add_Top;
    G4double Thickness_Steel_Add_Bot;
    G4double Thickness_Gap;
    G4double Length_WOM;
    G4double Thickness_WLS;
    G4double Length_sipm_box;

    G4double sipmSize;
    G4double sipmSizeSens;
    G4double sipmWindowThickness;
    G4double sipmSensThickness;
    G4double sipmSensThicknessTop;
    G4double sipmBaseThickness;


    // solids
    G4Box* ExpHallBox;
  //-------------------------------------------------------------------
    G4VSolid *sipmBase;
  //-------------------------------------------------------------------
    G4UnionSolid* ScintillatorBox;
    G4Box* sipmBaseBox;
    G4Box* sipmSens;
    G4Box* sipmSensTop;
    G4Box* sipmWindowAll;
    G4Box* sipmHole;
    G4Box* sipmBox;
    G4VSolid *OuterTube;
    G4VSolid *AirGapOut;
    G4VSolid *WLSOut;
    G4VSolid *WOMTube;
    G4VSolid *HoleBox;
    G4VSolid *HoleSct;
    G4VSolid *WLSIn;
    G4VSolid *AirGapIn;
    G4VSolid *InnerTube;
    G4VSolid *PMMARingLower;
    G4VSolid *PMMADisk;
    G4VSolid *PMMAHat;
    G4VSolid *SteelAdd;
    G4VSolid *SctInside;
    G4VSolid *WLSRing;
    G4VSolid *AirRingOut;
    G4VSolid *AirRingIn;
    G4VSolid *PMMARing;
    G4VSolid *EmptySteelBoxWithHole;
    G4VSolid *EmptyReflectBoxWithHole;
    G4SubtractionSolid *ScintillatorBoxWithHole;
    G4SubtractionSolid *sipmWindow;

    // logical volumes
    G4LogicalVolume* ExpHallLV;
    G4LogicalVolume* sipmSensLV;
    G4LogicalVolume* sipmSensTopLV;
    G4LogicalVolume* sipmWindowLV;
    G4LogicalVolume* sipmBaseBoxLV;
    G4LogicalVolume* sipmBoxLV;
    G4LogicalVolume* ScintillatorBoxLV;
    G4LogicalVolume* ReflectBoxLV;
    G4LogicalVolume* SteelBoxLV;
    G4LogicalVolume *OuterTubeLV;
    G4LogicalVolume *WOMTubeLV;
    G4LogicalVolume *InnerTubeLV;
    G4LogicalVolume *PMMARingLowerLV;
    G4LogicalVolume *PMMADiskLV;
    G4LogicalVolume *AirGapOutLV;
    G4LogicalVolume *AirGapInLV;
    G4LogicalVolume *WLSOutLV;
    G4LogicalVolume *WLSInLV;
    G4LogicalVolume *PMMAHatLV;
    G4LogicalVolume *AirRingOutLV;
    G4LogicalVolume *AirRingInLV;
    G4LogicalVolume *PMMARingLV;
    G4LogicalVolume *SteelAddLV;
    G4LogicalVolume *SctInsideLV;

    // physical volumes
    G4VPhysicalVolume* ExpHallPV;
    G4VPhysicalVolume* SteelBoxPV;
    G4VPhysicalVolume* ReflectBoxPV;
    G4VPhysicalVolume* ScintillatorBoxPV;
    G4VPhysicalVolume* sipmBasePV;
    G4VPhysicalVolume* sipmWindowPV;
    G4VPhysicalVolume* sipmSensTopPV;

    std::vector<G4VPhysicalVolume*> OuterTubePV_vect;
    std::vector<G4VPhysicalVolume*> WOMTubePV_vect;
    std::vector<G4VPhysicalVolume*> InnerTubePV_vect;
    std::vector<G4VPhysicalVolume*> PMMARingLowerPV_vect;
    std::vector<G4VPhysicalVolume*> PMMADiskPV_vect;
    std::vector<G4VPhysicalVolume*> AirGapOutPV_vect;
    std::vector<G4VPhysicalVolume*> AirGapInPV_vect;
    std::vector<G4VPhysicalVolume*> WLSOutPV_vect;
    std::vector<G4VPhysicalVolume*> WLSInPV_vect;
    std::vector<G4VPhysicalVolume*> PMMAHatPV_vect;
    std::vector<G4VPhysicalVolume*> AirRingOutPV_vect;
    std::vector<G4VPhysicalVolume*> AirRingInPV_vect;
    std::vector<G4VPhysicalVolume*> PMMARingPV_vect;
    std::vector<G4VPhysicalVolume*> SteelAddPV_vect;
    std::vector<G4VPhysicalVolume*> SctInsidePV_vect;
    std::vector<G4VPhysicalVolume*> sipmSensPV_vect;
    std::vector<G4VPhysicalVolume*> sipmBoxPV_vect;

    // visualisation
    G4Color blue;
    G4Color grey;
    G4Color green;
    G4Color red;
    G4Color white;
    G4Color cyan;
    G4Color magenta;
    G4Color yellow;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
