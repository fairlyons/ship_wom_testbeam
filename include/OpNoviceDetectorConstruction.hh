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

#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
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

    std::vector<std::pair<G4double,G4double>>  WOM_coord_vec;

  private:
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4Material* steel;
    G4Material* air;
    G4Material* LAB_PPO;
    G4Material* Bis_MSB;
    G4Material* PMMA_side;
    G4Material* PMMA_bottom;

//geometrical parameters
    G4double SteelX;
    G4double SteelY;
    G4double SteelZ;
    G4double WallThickness_XY;
    G4double WallThickness_Z_Bottom;
    G4double WallThickness_Z_Cover;
    G4double SctX;
    G4double SctY;
    G4double SctZ;
    
    G4double Additional_Length;
    
    G4double Diam_In_In;
    G4double Diam_In_Out;
    G4double Diam_Out_In;
    G4double Diam_Out_Out;
    G4double Diam_WOM_In;
    G4double Diam_WOM_Out;
    G4double Diam_Hole;
    G4double Diam_Steel_Add;
    G4double Diam_Hat;
    G4double Length_1;
    G4double Length_2;
    G4double Thickness_Ring;
    G4double Thickness_Disk;
    G4double Thickness_Hat;
    G4double Thickness_Steel_Add;
    G4double Length_WOM;
    G4double Thickness_WLS;
    
    G4double delta_X;
    G4double delta_Y;

    G4double OlengthOuter;
    G4double OlengthAG1;
    G4double OlengthInside;
    G4double sipmbasewidth;


//solids
    G4Box* expHall_box;
  //-------------------------------------------------------------------
    G4VSolid *sipm_base;
  //-------------------------------------------------------------------
  G4Box* ScintilatorBox;
  G4Box* SteelBox;
  G4VSolid *Outer_tube;
  G4VSolid *Air_gap1;
  G4VSolid *WLS_tube1;
  G4VSolid *WOM_tube;
  G4VSolid *Hole_box;
  G4VSolid *Hole_sct;
  G4VSolid *WLS_tube2;
  G4VSolid *Air_gap2;
  G4VSolid * Inner_tube;
  G4VSolid * PMMA_Ring;
  G4VSolid * PMMA_disk;
  G4VSolid *PMMA_Hat;
  G4VSolid *SteelAdd;
  G4VSolid *SctInside;
  G4VSolid *WLS_ring;
  G4SubtractionSolid *EmptySteelBoxWithHole;
  G4SubtractionSolid *ScintilatorBoxWithHole;

//logical volumes
    G4LogicalVolume* expHall_log;
    G4LogicalVolume* sipm_base_log;
    G4LogicalVolume* ScintilatorBox_log;
    G4LogicalVolume* SteelBox_log;
    G4LogicalVolume *Outer_tube_log;
    G4LogicalVolume *WOM_tube_log;
    G4LogicalVolume *Inner_tube_log;
    G4LogicalVolume *PMMA_Ring_log;
    G4LogicalVolume *PMMA_disk_log;
    G4LogicalVolume *Air_gap1_log;
    G4LogicalVolume *Air_gap2_log;
    G4LogicalVolume *WLS_tube1_log;
    G4LogicalVolume *WLS_tube2_log;
    G4LogicalVolume *PMMA_Hat_log;
    G4LogicalVolume *Steel_Add_log;
    G4LogicalVolume *Sct_Inside_log;

//physical volumes
    G4VPhysicalVolume* expHall_phys;
    std::vector<G4VPhysicalVolume*> sipm_base_phys_vect;
    G4VPhysicalVolume* SteelBox_phys;
    G4VPhysicalVolume* ScintilatorBox_phys;

  std::vector<G4VPhysicalVolume*> Outer_tube_phys_vect;
  std::vector<G4VPhysicalVolume*> WOM_tube_phys_vect;
  std::vector<G4VPhysicalVolume*> Inner_tube_phys_vect;
  std::vector<G4VPhysicalVolume*> PMMA_Ring_phys_vect;
  std::vector<G4VPhysicalVolume*> PMMA_Disk_phys_vect;
  std::vector<G4VPhysicalVolume*> Air_gap_1_phys_vect;
  std::vector<G4VPhysicalVolume*> Air_gap_2_phys_vect;
  std::vector<G4VPhysicalVolume*> WLS_tube1_phys_vect;
  std::vector<G4VPhysicalVolume*> WLS_tube2_phys_vect;
  std::vector<G4VPhysicalVolume*> PMMA_Hat_phys_vect;
  std::vector<G4VPhysicalVolume*> Steel_Add_phys_vect;
  std::vector<G4VPhysicalVolume*> Sct_Inside_phys_vect;

//visualisation

  G4Color blue;
  G4Color blue_trans;
  G4Color green;
  G4Color red;
  G4Color white;
  G4Color cyan;
  G4Color magenda;
  G4Color DircColor;
  G4Color SensColor;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpNoviceDetectorConstruction_h*/
