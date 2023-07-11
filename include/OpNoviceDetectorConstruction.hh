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
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4GenericTrap.hh"
#include "G4MultiUnion.hh"

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
    G4Material* Al;
    G4Material* Si;
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
    G4double WallThickness_Z_Cover;
    G4double ReflectZ;
    G4double ReflectThick;
    G4double SctZ;
    
    G4double Diam_In_In;
    G4double Diam_In_Out;
    G4double Diam_Out_In;
    G4double Diam_Out_Out;
    G4double Diam_WOM_In;
    G4double Diam_WOM_Out;
    G4double Diam_Hole;
    G4double Diam_Steel_Add;
    G4double Diam_Hat;
    G4double Length_Out;
    G4double Length_In;
    G4double Thickness_Ring;
    G4double Thickness_Disk;
    G4double Thickness_Hat;
    G4double Thickness_Steel_Add;
    G4double Thickness_Gap;
    G4double Thickness_SupportRing;
    G4double Thickness_WLS;
    G4double Length_WOM;

    G4double sipmSize;
    G4double sipmSizeSens;
    G4double sipmWindowThickness;
    G4double sipmSensThickness;
    G4double sipmSensThicknessTop;
    G4double sipmBaseThickness;


    // solids
    G4Box* expHall_box;
  //-------------------------------------------------------------------
    G4VSolid *sipm_base;
  //-------------------------------------------------------------------
    G4GenericTrap* SteelBox;
    G4GenericTrap* ReflectBox;
	G4SubtractionSolid* EmptySteelBox;
    G4GenericTrap* SteelBeamV;
    G4GenericTrap* SteelBeamH;
	G4UnionSolid* SteelBeam;
	G4SubtractionSolid* EmptyReflectBox;
    G4GenericTrap* ScintillatorBoxLT;
    G4GenericTrap* ScintillatorBoxRT;
    G4GenericTrap* ScintillatorBoxLB;
    G4GenericTrap* ScintillatorBoxRB;
    G4MultiUnion* ScintillatorInBeams;
    G4Box* sipmBaseBox;
    G4Box* sipmSens;
    G4Box* sipmSensTop;
    G4Box* sipmWindowAll;
    G4Box* sipmHole;
    G4VSolid* Outer_tube;
    G4VSolid* Air_gap_out;
    G4VSolid* WLS_tube_out;
    G4VSolid* WOM_tube;
    G4VSolid* Hole_box;
    G4VSolid* WLS_tube_in;
    G4VSolid* Air_gap_in;
    G4VSolid* Inner_tube;
    G4VSolid* PMMA_Ring;
    G4VSolid* PMMA_disk;
    G4VSolid* PMMA_Hat;
    G4VSolid* SteelAdd;
    G4VSolid* SctInside;
    G4VSolid* WLS_ring;
    G4VSolid* Air_ring_out;
    G4VSolid* Air_ring_in;
    G4VSolid* PMMA_ring_lower;
    G4SubtractionSolid* EmptySteelBoxWithHole;
    G4SubtractionSolid* SteelBeamVWithHoles;
    G4SubtractionSolid* SteelBeamHWithHoles;
    G4SubtractionSolid* EmptyReflectBoxWithHole;
    G4SubtractionSolid* ScintillatorBoxWithHoleLT;
    G4SubtractionSolid* ScintillatorBoxWithHoleRT;
    G4SubtractionSolid* ScintillatorBoxWithHoleLB;
    G4SubtractionSolid* ScintillatorBoxWithHoleRB;
    G4SubtractionSolid* sipmWindow;

    // logical volumes
    G4LogicalVolume* expHall_log;
    G4LogicalVolume* SteelBox_log;
    G4LogicalVolume* SteelBeam_log;
    G4LogicalVolume* ReflectBox_log;
    G4LogicalVolume* ScintillatorBoxLT_log;
    G4LogicalVolume* ScintillatorBoxRT_log;
    G4LogicalVolume* ScintillatorBoxLB_log;
    G4LogicalVolume* ScintillatorBoxRB_log;
    G4LogicalVolume* ScintillatorInBeams_log;
    G4LogicalVolume* WOM_cell_log;
    G4LogicalVolume* sipmSens_log;
    G4LogicalVolume* sipmSensTop_log;
    G4LogicalVolume* sipmWindow_log;
    G4LogicalVolume* sipmBaseBox_log;
    G4LogicalVolume* Outer_tube_log;
    G4LogicalVolume* WOM_tube_log;
    G4LogicalVolume* Inner_tube_log;
    G4LogicalVolume* PMMA_Ring_log;
    G4LogicalVolume* PMMA_disk_log;
    G4LogicalVolume* Air_gap_out_log;
    G4LogicalVolume* Air_gap_in_log;
    G4LogicalVolume* WLS_tube_out_log;
    G4LogicalVolume* WLS_tube_in_log;
    G4LogicalVolume* PMMA_Hat_log;
    G4LogicalVolume* Air_ring_out_log;
    G4LogicalVolume* Air_ring_in_log;
    G4LogicalVolume* PMMA_ring_lower_log;
    G4LogicalVolume* Steel_Add_log;
    G4LogicalVolume* Sct_Inside_log;

    // physical volumes
    G4VPhysicalVolume* expHall_phys;
    G4VPhysicalVolume* SteelBox_phys;
    G4VPhysicalVolume* SteelBeam_phys;
    G4VPhysicalVolume* ReflectBox_phys;
    G4VPhysicalVolume* ReflectInBeamsV_phys;
    G4VPhysicalVolume* ReflectInBeamsH_phys;
    G4VPhysicalVolume* ScintillatorBoxLT_phys;
    G4VPhysicalVolume* ScintillatorBoxRT_phys;
    G4VPhysicalVolume* ScintillatorBoxLB_phys;
    G4VPhysicalVolume* ScintillatorBoxRB_phys;
    G4VPhysicalVolume* ScintillatorInBeams_phys;
    G4VPhysicalVolume* sipmBase_phys;
    G4VPhysicalVolume* sipmWindow_phys;
    G4VPhysicalVolume* sipmSensTop_phys;

    std::vector<G4VPhysicalVolume*> Outer_tube_phys_vec;
    std::vector<G4VPhysicalVolume*> WOM_tube_phys_vec;
    std::vector<G4VPhysicalVolume*> Inner_tube_phys_vec;
    std::vector<G4VPhysicalVolume*> PMMA_Ring_phys_vec;
    std::vector<G4VPhysicalVolume*> PMMA_Disk_phys_vec;
    std::vector<G4VPhysicalVolume*> Air_gap_out_phys_vec;
    std::vector<G4VPhysicalVolume*> Air_gap_in_phys_vec;
    std::vector<G4VPhysicalVolume*> WLS_tube_out_phys_vec;
    std::vector<G4VPhysicalVolume*> WLS_tube_in_phys_vec;
    std::vector<G4VPhysicalVolume*> PMMA_Hat_phys_vec;
    std::vector<G4VPhysicalVolume*> Air_ring_out_phys_vec;
    std::vector<G4VPhysicalVolume*> Air_ring_in_phys_vec;
    std::vector<G4VPhysicalVolume*> PMMA_ring_lower_phys_vec;
    std::vector<G4VPhysicalVolume*> Steel_Add_phys_vec;
    std::vector<G4VPhysicalVolume*> Sct_Inside_phys_vec;
    std::vector<G4VPhysicalVolume*> sipm_phys_vec;
    std::vector<G4VPhysicalVolume*> WOM_cells_phys_vec;

    // visualisation
    G4Color blue;
    G4Color grey;
    G4Color blue_trans;
    G4Color green;
    G4Color red;
    G4Color white_trans;
    G4Color cyan;
    G4Color magenta;
    G4Color yellow;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
