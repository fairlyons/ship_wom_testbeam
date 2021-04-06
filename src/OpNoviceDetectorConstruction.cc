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

#include "OpNoviceDetectorConstruction.hh"
// #include "OpNoviceDetectorMessenger.hh"s
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
// #include"SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4SDParticleFilter.hh"
#include "G4UnionSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4PhysicalConstants.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction::OpNoviceDetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fExpHall_x = 0.8*m;
  fExpHall_y = 0.6*m;
  fExpHall_z = 0.5*m;

  SteelX = 120*cm;
  SteelY = 80*cm;
  SteelZ = 335*mm;
//  SteelZ = 135*mm;
  WallThickness_XY = 10*mm;
  WallThickness_Z_Bottom = 20*mm;
  WallThickness_Z_Cover = 5*mm;
  SctX = SteelX - 2*WallThickness_XY;
  SctY = SteelY - 2*WallThickness_XY;
  SctZ = SteelZ - WallThickness_Z_Bottom - WallThickness_Z_Cover;
    
  Additional_Length = 0.*mm;
    
  Diam_In_In = 44*mm;
  Diam_In_Out = 50*mm;
  Diam_Out_In = 64*mm;
  Diam_Out_Out = 70*mm;
  Diam_WOM_In = 54*mm;
  Diam_WOM_Out = 60*mm;
  Diam_Hole = 72*mm;
  Diam_Steel_Add = 143*mm;
  Diam_Hat = 145*mm;
  Length_1 = 226*mm + Additional_Length;
  Length_2 = 206*mm + Additional_Length;
  Thickness_Ring = 4*mm;
  Thickness_Disk = 4*mm;
  Thickness_Hat = 6*mm;
  Thickness_Steel_Add = 15*mm;
  Length_WOM = 230*mm + Additional_Length;
//  Length_WOM = 300*mm + Additional_Length;
  Thickness_WLS = 0.05*mm;
    
  delta_X = SteelX/2 - 91.5*mm;
  delta_Y = SteelY/2 - 91.5*mm;

//   WOM_coord_vec = {
//                     {-delta_X, delta_Y},
//                     {delta_X, delta_Y},
//                     {-delta_X, -delta_Y},
//                     {delta_X, -delta_Y}
//                   };
    WOM_coord_vec = {
                    {-400.*mm, 0.},
                    {400.*mm, 0.}
                  };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction::~OpNoviceDetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();

    steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    // G4Material* Al=nist->FindOrBuildMaterial("G4_Al");

  G4double a, z, density;
  G4int nelements, ncomponent, natoms;

// Air
//
    G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

    air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
    air->AddElement(N, 70.*perCent);
    air->AddElement(O, 30.*perCent);


  G4Element* H = new G4Element("Hydrogen", "H", 1 , 1.01*g/mole);

  //Linear alkyl benzene (LAB)
  G4Element* C = new G4Element("Carbon", "C", 6 , 12.01*g/mole);
  G4Material* LAB = new G4Material("LAB",density=0.86*g/cm3,ncomponent=2);
  LAB->AddElement(H,natoms=28);
  LAB->AddElement(C,natoms=17);
  //Diphenyloxazole (PPO)
  G4Material* PPO = new G4Material("PPO",density=1.184*g/cm3,ncomponent=4);
  PPO->AddElement(H,natoms=11);
  PPO->AddElement(C,natoms=15);
  PPO->AddElement(N,natoms=1);
  PPO->AddElement(O,natoms=1);
  //Scintilator (LAB+PPO) 23233 cm^3 23.233 l
      //G4Material* LAB_PPO = new G4Material("LAB_PPO", density=239.5*g/mole, ncomponent=2);
    LAB_PPO = new G4Material("LAB_PPO", density=0.9*g/cm3, ncomponent=2);
    LAB_PPO->AddMaterial(LAB, 87.8*perCent);
    LAB_PPO->AddMaterial(PPO, 12.2*perCent);
  // Bis-MSB WLS
    Bis_MSB = new G4Material("Bis_MSB",density=1.07*g/cm3,ncomponent=2);
    Bis_MSB->AddElement(H,natoms=22);
    Bis_MSB->AddElement(C,natoms=24);
  // PMMA side
    PMMA_side = new G4Material("PMMA",density=1.200*g/cm3,ncomponent=2);
    PMMA_side->AddElement(H,natoms=2);
    PMMA_side->AddElement(C,natoms=4);
  // PMMA bottom
    PMMA_bottom = new G4Material("PMMA",density=1.200*g/cm3,ncomponent=2);
    PMMA_bottom->AddElement(H,natoms=2);
    PMMA_bottom->AddElement(C,natoms=4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineMPTs()
{
  //------------------------------------------------------------------------------
  //----------------------------- LAB_PPO -----------------------------
  //------------------------------------------------------------------------------
    G4double rindex_LAB_PPO[100];
    G4double photon_en_LAB_PPO[100];
    auto Ridndex_LAB_PPO = [=](G4double wl)
    {
        G4double rind=1.;
        G4double B[4], C[4];
        B[0] =0.821384; C[0] =94.7625;
        B[1] =0.311375; C[1] =160.751;
        B[2] =0.0170099;C[2] =219.575;
        B[3] =0.608268; C[3] =9385.54;
        for(int term = 0; term<4; term++)
          rind+=B[term]/( 1.-(C[term]/wl)*(C[term]/wl) );
        return sqrt(rind);
    };

    G4double wl;
    for(int i=0;i<100;i++)
    {
      wl = 250.+5.*i;
      photon_en_LAB_PPO[i]=1240./wl*eV;
      rindex_LAB_PPO[i]=Ridndex_LAB_PPO(wl);
    }
    G4MaterialPropertiesTable *MPT_LAB_PPO = new G4MaterialPropertiesTable();
    MPT_LAB_PPO -> AddConstProperty("SCINTILLATIONYIELD",10800./MeV);
    MPT_LAB_PPO -> AddProperty("RINDEX", photon_en_LAB_PPO, rindex_LAB_PPO, 100)->SetSpline(true);
          
    // emission
    G4double photonWaveLength3[201] = {300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324,
        325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356,
        357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388,
        389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420,
        421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452,
        453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
        485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500};
    G4double photon_en_LAB_PPO_2[201];
    for(int i=0; i<201; i++)
        photon_en_LAB_PPO_2[i] = 1240./photonWaveLength3[i]*eV;
    G4double scintilFast_LAB_PPO[201] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.045,
        0.09, 0.135, 0.18, 0.225, 0.27, 0.315, 0.36, 0.405, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.93, 0.91, 0.89, 0.87, 0.85, 0.83, 0.81,
        0.79, 0.77, 0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.77, 0.74, 0.71, 0.68, 0.65, 0.62, 0.59, 0.56, 0.53, 0.5, 0.4925, 0.485,
        0.4775, 0.47, 0.4625, 0.455, 0.4475, 0.44, 0.4325, 0.425, 0.41, 0.395, 0.38, 0.365, 0.35, 0.335, 0.32, 0.305, 0.29, 0.275, 0.2675, 0.26, 0.2525, 0.245, 0.2375,
        0.23, 0.2225, 0.215, 0.2075, 0.2, 0.1925, 0.185, 0.1775, 0.17, 0.1625, 0.155, 0.1475, 0.14, 0.1325, 0.125, 0.12, 0.115, 0.11, 0.105, 0.1, 0.095, 0.09, 0.085,
        0.08, 0.075, 0.0725, 0.07, 0.0675, 0.065, 0.0625, 0.06, 0.0575, 0.055, 0.0525, 0.05, 0.0475, 0.045, 0.0425, 0.04, 0.0375, 0.035, 0.0325, 0.03, 0.0275, 0.025,
        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.0225, 0.02, 0.0175, 0.015, 0.0125, 0.01, 0.0075, 0.005, 0.0025, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    MPT_LAB_PPO -> AddProperty("FASTCOMPONENT",photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201) -> SetSpline(true);
    MPT_LAB_PPO -> AddProperty("SLOWCOMPONENT",photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201) -> SetSpline(true);
                                                                  
  // transmission
    G4double waveLength2[211] = {320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344,
              345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 360.75, 361.5, 362.25, 363, 363.75, 364.5, 365.25, 366, 366.75, 367.5,
              368, 368.5, 369, 369.5, 370, 370.5, 371, 371.5, 372, 372.5, 373.5, 374.5, 375.5, 376.5, 377.5, 378.5, 379.5, 380.5, 381.5, 382.5, 382.75, 383, 383.25,
              383.5, 383.75, 384, 384.25, 384.5, 384.75, 385, 385.4, 385.8, 386.2, 386.6, 387, 387.4, 387.8, 388.2, 388.6, 389, 389.6, 390.2, 390.8, 391.4, 392, 392.6,
              393.2, 393.8, 394.4, 395, 395.5, 396, 396.5, 397, 397.5, 398, 398.5, 399, 399.5, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413,
              414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443,
              444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473,
              474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500};
    G4double photon_en_LAB_PPO_3[211];
    for(int i=0; i<211; i++)
      photon_en_LAB_PPO_3[i] = 1240./waveLength2[i]*eV;
    G4double transCoef2[211] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.0134, 0.0258, 0.0382, 0.0506, 0.063, 0.0754,
              0.0878, 0.1002, 0.1126, 0.125, 0.1475, 0.17, 0.1925, 0.215, 0.2375, 0.26, 0.2825, 0.305, 0.3275, 0.35, 0.365, 0.38, 0.395, 0.41, 0.425, 0.44, 0.455,
              0.47, 0.485, 0.5, 0.495, 0.49, 0.485, 0.48, 0.475, 0.47, 0.465, 0.46, 0.455, 0.45, 0.4675, 0.485, 0.5025, 0.52, 0.5375, 0.555, 0.5725, 0.59, 0.6075,
              0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.6275, 0.63, 0.6325, 0.635, 0.6375, 0.64, 0.6425, 0.645, 0.6475, 0.65,
              0.6425, 0.635, 0.6275, 0.62, 0.6125, 0.605, 0.5975, 0.59, 0.5825, 0.575, 0.5975, 0.62, 0.6425, 0.665, 0.6875, 0.71, 0.7325, 0.755, 0.7775, 0.8, 0.8, 0.8,
              0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8075, 0.815, 0.8225, 0.83, 0.8375, 0.845, 0.8525, 0.86, 0.8675, 0.875, 0.8825, 0.89, 0.8975, 0.905, 0.9125,
              0.92, 0.9275, 0.935, 0.9425, 0.95, 0.9525, 0.955, 0.9575, 0.96, 0.9625, 0.965, 0.9675, 0.97, 0.9725, 0.975, 0.975, 0.975, 0.975, 0.975, 0.975, 0.975,
              0.975, 0.975, 0.975, 0.975, 0.9725, 0.97, 0.9675, 0.965, 0.9625, 0.96, 0.9575, 0.955, 0.9525, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95,
              0.95,0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.945, 0.94, 0.935, 0.93, 0.925, 0.92, 0.915, 0.91, 0.905, 0.9, 0.9, 0.9, 0.9,
              0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82, 0.81, 0.8};
    G4double absLen2[211];
    for(int i=0; i<211; i++)
      absLen2[i] = -1./std::log(transCoef2[i])*m;
    MPT_LAB_PPO -> AddProperty("ABSLENGTH",photon_en_LAB_PPO_3,absLen2,211) -> SetSpline(true);        
    MPT_LAB_PPO -> AddConstProperty("RESOLUTIONSCALE", 2.0);
    MPT_LAB_PPO -> AddConstProperty("FASTTIMECONSTANT", 5.2*ns);
    MPT_LAB_PPO -> AddConstProperty("SLOWTIMECONSTANT", 18.4*ns);
    MPT_LAB_PPO->AddConstProperty("YIELDRATIO",0.78);
    LAB_PPO -> SetMaterialPropertiesTable(MPT_LAB_PPO);



  //------------------------------------------------------------------------------
  //----------------------------- PMMA -----------------------------
  //------------------------------------------------------------------------------
  const G4int pmma_mpt_entr = 13;
  G4double pmma_wl[pmma_mpt_entr] = {700.,  600.,  550.,  500.,  450.,  400.,  390.,  380.,  370.,  350.,  320., 310.,  300. };
  G4double pmma_rind[pmma_mpt_entr] = {1.489, 1.492, 1.495, 1.498, 1.502, 1.511, 1.512, 1.514, 1.516, 1.522, 1.54, 1.541, 1.542};
  G4double pmma_en[pmma_mpt_entr];
  for(int i=0; i<pmma_mpt_entr; i++ )
  {
    pmma_en[i]=1240./pmma_wl[i]*eV;
  }
  G4double pmma_side_abslen[pmma_mpt_entr] =  { 10.55*mm , 18.23*mm , 24.6*mm, 36.07*mm, 39.7*mm, 42.6*mm, 43.69*mm, 45.77*mm, 52.97*mm, 61.48*mm, 66.44*mm, 70.39*mm, 79.11*mm};
  G4double pmma_bottom_abslen[pmma_mpt_entr] = {  0.01*mm,  0.01*mm,   0.01*mm, 0.01*mm, 1.31*mm, 4.26*mm, 14.98*mm, 24.09*mm, 28.56*mm, 30.35*mm, 32.39*mm, 33.93*mm, 37.32*mm};

  G4MaterialPropertiesTable *MPT_PMMA_side = new G4MaterialPropertiesTable();
  MPT_PMMA_side->AddProperty("RINDEX", pmma_en, pmma_rind, pmma_mpt_entr)->SetSpline(true);
  MPT_PMMA_side->AddProperty("ABSLENGTH", pmma_en, pmma_side_abslen, pmma_mpt_entr)->SetSpline(true);
  PMMA_side->SetMaterialPropertiesTable(MPT_PMMA_side);

  G4MaterialPropertiesTable *MPT_PMMA_bottom = new G4MaterialPropertiesTable();
  MPT_PMMA_bottom->AddProperty("RINDEX", pmma_en, pmma_rind, pmma_mpt_entr)->SetSpline(true);
  MPT_PMMA_bottom->AddProperty("ABSLENGTH", pmma_en, pmma_bottom_abslen, pmma_mpt_entr);
  PMMA_bottom->SetMaterialPropertiesTable(MPT_PMMA_bottom);



  //------------------------------------------------------------------------------
  //----------------------------- Bis_MSB -----------------------------
  //------------------------------------------------------------------------------
  G4MaterialPropertiesTable *MPT_Bis_MSB = new G4MaterialPropertiesTable();
  G4double waveLength3[7] = {300., 340., 380., 400., 420., 460., 500};
  G4double photonEnergy5[7];
  for(int i=0; i<7; i++)
    photonEnergy5[i] = 1240./waveLength3[i]*eV;
  G4double absLen3[7] = {10*nm, 10*nm, 10*nm, 1*mm, 200*m, 200*m, 200*m};
  MPT_Bis_MSB->AddProperty("WLSABSLENGTH", photonEnergy5, absLen3, 7);

  // BIS reemission
  G4double waveLength4[16] = {380., 390., 400., 410., 420., 430., 440.,
                                            450., 460., 470., 480., 490., 500., 510., 520., 530};
  G4double photonEnergy6[16];
  for(int i=0; i<16; i++)
    photonEnergy6[i] = 1240./waveLength4[i]*eV;
  G4double reEmit4[16] = {0., 0., 0.1, 0.8, 1., 0.8, 0.5,
                                        0.45, 0.3, 0.2, 0.15, 0.1, 0.05, 0.05, 0.05, 0.};
  G4double ppckovEmit[8] = { 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.6401*eV , 3.0402*eV , 3.5403*eV , 3.8404*eV};

  G4double rindexWLS[8] = { 1.5, 1.5, 1.5, 1.5, 1.504 , 1.505 , 1.515 , 1.52 };

  MPT_Bis_MSB->AddProperty("WLSCOMPONENT", photonEnergy6, reEmit4, 16);
  MPT_Bis_MSB->AddConstProperty("WLSTIMECONSTANT", 3.*ns);
      //MPT_Bis_MSB-> AddConstProperty("WLSMEANNUMBERPHOTONS",1.0);
  MPT_Bis_MSB->AddProperty("RINDEX", ppckovEmit, rindexWLS, 8)->SetSpline(true);
  Bis_MSB->SetMaterialPropertiesTable(MPT_Bis_MSB);



  //------------------------------------------------------------------------------
  //----------------------------- Air -----------------------------
  //------------------------------------------------------------------------------
  G4double photonEnergy_Air[2] = { 2.*eV,  5.*eV };
  G4double refractiveIndex_Air[2] = { 1.00, 1.00 };
  G4MaterialPropertiesTable* MPT_Air = new G4MaterialPropertiesTable();
  MPT_Air->AddProperty("RINDEX", photonEnergy_Air, refractiveIndex_Air, 2);
  air->SetMaterialPropertiesTable(MPT_Air);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineSurfaces()
{
  //------------------------------------------------------------------------------
  //----------------------------- AirBoxSurface -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* AirBoxSurface = new G4OpticalSurface("AirBoxSurface"); // air-PMMA border
  AirBoxSurface->SetType(dielectric_dielectric);
  AirBoxSurface->SetFinish(polished);
  AirBoxSurface->SetModel(glisur);
         
  const G4int num1 = 2;
  G4double pp1[num1] = {2.*eV, 5.*eV};
  G4double reflectivity1[num1] = {0.999, 0.999}; // by A. K.
  G4MaterialPropertiesTable *MPTsurf_AirBoxSurface = new G4MaterialPropertiesTable();
  MPTsurf_AirBoxSurface->AddProperty("REFLECTIVITY", pp1, reflectivity1, num1);
  AirBoxSurface -> SetMaterialPropertiesTable(MPTsurf_AirBoxSurface);

  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
  {
    new G4LogicalBorderSurface( (std::string("TubeSurface_")+std::to_string(pos)).c_str(), Air_gap_1_phys_vect[pos], Outer_tube_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("TubeSurface_")+std::to_string(pos)).c_str(), Outer_tube_phys_vect[pos], Air_gap_1_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("PMMASurface_")+std::to_string(pos)).c_str(), expHall_phys, PMMA_Disk_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("PMMASurface_")+std::to_string(pos)).c_str(), expHall_phys, Inner_tube_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("TubeSurface_")+std::to_string(pos)).c_str(), PMMA_Ring_phys_vect[pos], Air_gap_1_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("TubeSurface_")+std::to_string(pos)).c_str(), PMMA_Ring_phys_vect[pos], Air_gap_2_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("TubeSurface_")+std::to_string(pos)).c_str(), Inner_tube_phys_vect[pos], Air_gap_2_phys_vect[pos], AirBoxSurface);
    new G4LogicalBorderSurface( (std::string("TubeSurface_")+std::to_string(pos)).c_str(), Air_gap_2_phys_vect[pos], Inner_tube_phys_vect[pos], AirBoxSurface);
  }




  //------------------------------------------------------------------------------
  //----------------------------- Steel -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* SteelBoxSurface = new G4OpticalSurface("SteelBoxSurface");
  SteelBoxSurface -> SetType(dielectric_metal);
  SteelBoxSurface -> SetFinish(ground);
  SteelBoxSurface -> SetModel(glisur);

  G4double waveLength5[18] = {330., 340., 350., 360., 370., 380., 390., 400.,
                           410., 420., 430., 440., 450., 460., 470., 480., 490., 500.};
  G4double photonEnergy7[18];
  for(int i=0; i<18; i++)
    photonEnergy7[i] = 1240./waveLength5[i]*eV;
//           G4double reflectSteel[18] = {0.430, 0.440, 0.449, 0.457, 0.463, 0.469, 0.474, 0.479,
//                          0.483, 0.487, 0.490, 0.493, 0.496, 0.499, 0.501, 0.504, 0.507, 0.508};    // from refractiveindex.info
//            G4double reflectSteel[18] = {0.467, 0.472, 0.477, 0.482, 0.486, 0.491, 0.496, 0.500,
//                           0.506, 0.510, 0.515, 0.519, 0.521, 0.522, 0.530, 0.538, 0.542, 0.546};    // from journal article

  G4double reflectSteel[18] = {0.35, 0.36, 0.365, 0.37, 0.375, 0.38, 0.39, 0.395,
                            0.4, 0.405, 0.41, 0.415, 0.42, 0.42, 0.425, 0.43, 0.43, 0.435};         // from borexino

//  G4double reflectSteel[18] = {0.97752157, 0.9722873 , 0.96680556, 0.97042264, 0.97126233, 0.97271017,
//                               0.97712734, 0.97275713, 0.97658953, 0.97920518, 0.97380716, 0.98033894,
//                               0.9776999 , 0.978477  , 0.97228729, 0.97114843, 0.97375261, 0.97113254};

  G4MaterialPropertiesTable *MPTsurf_Steel = new G4MaterialPropertiesTable();      
  MPTsurf_Steel -> AddProperty("REFLECTIVITY", photonEnergy7, reflectSteel, 18);
  SteelBoxSurface -> SetMaterialPropertiesTable(MPTsurf_Steel);

  G4LogicalSkinSurface* Surface = new G4LogicalSkinSurface("SteelSurface",SteelBox_log,SteelBoxSurface);



  //------------------------------------------------------------------------------
  //----------------------------- PMMA WLS surface -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* FiberSurfaceInside = new G4OpticalSurface("FiberSurfaceInside");
  FiberSurfaceInside->SetType(dielectric_dielectric);
  FiberSurfaceInside->SetFinish(polished);
  FiberSurfaceInside->SetModel(glisur);

  const int nentries_PMMA_WLS_surf = 3;
  G4double photon_en_PMMA_sWLS[nentries_PMMA_WLS_surf] = {2.*eV,3.5*eV,5.*eV};
  G4double reflectivityPMMA[nentries_PMMA_WLS_surf] = {0.999,0.999, 0.999};

  G4MaterialPropertiesTable *MPTsurf_PMMA_WLS = new G4MaterialPropertiesTable();
  MPTsurf_PMMA_WLS->AddProperty("REFLECTIVITY", photon_en_PMMA_sWLS, reflectivityPMMA, nentries_PMMA_WLS_surf);

  FiberSurfaceInside->SetMaterialPropertiesTable(MPTsurf_PMMA_WLS);                                               

  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
  {
    new G4LogicalBorderSurface( (std::string("FiberInnerSurface_")+std::to_string(pos)).c_str(), WOM_tube_phys_vect[pos], WLS_tube1_phys_vect[pos], FiberSurfaceInside);
    new G4LogicalBorderSurface( (std::string("FiberInnerSurface_")+std::to_string(pos)).c_str(), WOM_tube_phys_vect[pos], WLS_tube2_phys_vect[pos], FiberSurfaceInside);
  }



  //------------------------------------------------------------------------------
  //----------------------------- Scintillator_PMMA surface -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* Scintillator_PMMA_Surface = new G4OpticalSurface("Scintillator_PMMA_Surface");
  Scintillator_PMMA_Surface->SetType(dielectric_dielectric);
  Scintillator_PMMA_Surface->SetFinish(polished);
  Scintillator_PMMA_Surface->SetModel(glisur);

  const G4int num  = 3;
  G4double pp[num] = {1.6*eV,3.44*eV, 5.0*eV};
  G4double reflectivitySct[num] = {0.999,0.999, 0.999};

  G4MaterialPropertiesTable *MPTsurf_Scintillator_PMMA = new G4MaterialPropertiesTable();
  MPTsurf_Scintillator_PMMA->AddProperty("REFLECTIVITY", pp, reflectivitySct, num);
  Scintillator_PMMA_Surface->SetMaterialPropertiesTable(MPTsurf_Scintillator_PMMA);

  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
  {
    new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintilatorBox_phys, Outer_tube_phys_vect[pos], Scintillator_PMMA_Surface);
    new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintilatorBox_phys, PMMA_Ring_phys_vect[pos], Scintillator_PMMA_Surface);
    new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintilatorBox_phys, Inner_tube_phys_vect[pos], Scintillator_PMMA_Surface);
    new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), Sct_Inside_phys_vect[pos], Inner_tube_phys_vect[pos], Scintillator_PMMA_Surface);
    new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), Sct_Inside_phys_vect[pos], PMMA_Disk_phys_vect[pos], Scintillator_PMMA_Surface);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineSolids()
{
  expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);
  //-------------------------------------------------------------------
  sipmbasewidth = 1*cm;
  sipm_base  = new G4Tubs("sipm_base", Diam_WOM_In/2 - Thickness_WLS , Diam_WOM_Out/2 + Thickness_WLS, sipmbasewidth/2, 0, 360*deg);
  //-------------------------------------------------------------------
  //Scintillator box
  ScintilatorBox = new G4Box("ScintilatorBox",SctX/2,SctY/2,SctZ/2);
//Steel box
  SteelBox = new G4Box("Steelbox",SteelX/2,SteelY/2,SteelZ/2);
                
  G4double Rin, Rout;
  G4double delta_Z;
    
//Outer tube
  Rin = Diam_Out_In/2;
  Rout = Diam_Out_Out/2;
  OlengthOuter = Length_2 + Thickness_Hat + Thickness_Ring;
  Outer_tube  = new G4Tubs("Outer_tube", Rin, Rout, OlengthOuter/2, 0, 360*deg);
//Air gap1
  Rin = Diam_WOM_Out/2 + Thickness_WLS;
  Rout = Diam_Out_In/2;
  OlengthAG1 = Length_WOM + Thickness_WLS;
  Air_gap1  = new G4Tubs("Air_gap1", Rin, Rout, OlengthAG1/2, 0, 360*deg);
//WLS tube 1
  Rin = Diam_WOM_Out/2;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WLS_tube1  = new G4Tubs("WLS_tube1", Rin, Rout, Length_WOM/2, 0, 360*deg);
//WOM tube
  Rin = Diam_WOM_In/2;
  Rout = Diam_WOM_Out/2;
  WOM_tube  = new G4Tubs("WOM_tube", Rin, Rout, Length_WOM/2, 0, 360*deg);
//Hole in box
  Rin = 0.0*mm;
  Rout = Diam_Out_Out/2;
double OlengthHoleBox = Length_2 + Thickness_Ring - Thickness_Steel_Add;
  Hole_box  = new G4Tubs("Hole_box", Rin, Rout, OlengthHoleBox/2, 0, 360*deg);
//Hole in scintilator
  Rin = 0*mm;
  Rout = Diam_Out_Out/2;
double OlengthHoleSct = Length_2 + Thickness_Ring - Thickness_Steel_Add - WallThickness_Z_Cover;
  Hole_sct = new G4Tubs("Hole_sct", Rin, Rout, OlengthHoleSct/2, 0, 360*deg);
//WLS tube 2
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_In/2;
double OlengthWLS2 = Length_WOM;
  WLS_tube2  = new G4Tubs("WLS_tube2", Rin, Rout, OlengthWLS2/2, 0, 360*deg);
//Air gap2
  Rin = Diam_In_Out/2;
  Rout = Diam_WOM_In/2 - Thickness_WLS;
double OlengthAG2 = Length_WOM + Thickness_WLS;
  Air_gap2  = new G4Tubs("Air_gap2", Rin, Rout, OlengthAG2/2, 0, 360*deg);
//Inner tube
  Rin = Diam_In_In/2;
  Rout = Diam_In_Out/2;
  Inner_tube  = new G4Tubs("Inner_tube", Rin, Rout, Length_1/2, 0, 360*deg);
//PMMA Ring
  Rin = Diam_In_Out/2;
  Rout = Diam_Out_In/2;
  PMMA_Ring  = new G4Tubs("PMMA_Ring", Rin, Rout, Thickness_Ring/2, 0, 360*deg);
//PMMA disk
  Rin = 0*mm;
  Rout = Diam_In_In/2;
  PMMA_disk  = new G4Tubs("PMMA_disk", Rin, Rout, Thickness_Disk/2, 0, 360*deg);
//PMMA "hat"
  Rin = Diam_Out_Out/2;
  Rout = Diam_Hat/2;
  PMMA_Hat = new G4Tubs("PMMA_Hat ", Rin, Rout, Thickness_Hat/2, 0, 360*deg);
//Additional steel
  Rin = Diam_Hole/2;
  Rout = Diam_Steel_Add/2;
double OlengthAdd = Thickness_Steel_Add;
  SteelAdd = new G4Tubs("SteelAdd ", Rin, Rout, OlengthAdd/2, 0, 360*deg);
//LAB&PPO inside tube
  Rin = 0.0*mm;
  Rout = Diam_In_In/2;
  OlengthInside = Length_1 - Thickness_Disk;
  SctInside = new G4Tubs("SctInside", Rin, Rout, OlengthInside/2, 0, 360*deg);
//WLS Ring
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
double OlengthWlsRing = Thickness_WLS;
  WLS_ring = new G4Tubs("WLS_ring ", Rin, Rout, OlengthWlsRing, 0, 360*deg);

    
  delta_Z = (WallThickness_Z_Bottom - WallThickness_Z_Cover)/2;
  G4SubtractionSolid *EmptySteelBox= new G4SubtractionSolid("EmptySteelBox",SteelBox, ScintilatorBox,0,G4ThreeVector(0, 0, delta_Z));
                
  G4double delta_Z_EmptySteelBoxWithHole = SteelZ/2 - OlengthHoleBox/2;
  G4double delta_Z_ScintilatorBoxWithHole = SctZ/2 - OlengthHoleSct/2;


  std::vector<G4SubtractionSolid*> EmptySteelBoxWithHole_tempvec;
  std::vector<G4SubtractionSolid*> ScintilatorBoxWithHole_tempvec;
  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
  {
    if(pos==0)
    {
      EmptySteelBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("EmptySteelBoxWithHole_")+std::to_string(pos)).c_str(),
                                              EmptySteelBox,Hole_box,0,
                                              G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_EmptySteelBoxWithHole)) );
      ScintilatorBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("ScintilatorBoxWithHole_")+std::to_string(pos)).c_str(),
                                                ScintilatorBox, Hole_sct, 0,
                                                G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_ScintilatorBoxWithHole)) );
    }
    else
    {
      EmptySteelBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("EmptySteelBoxWithHole_")+std::to_string(pos)).c_str(),
                                              EmptySteelBoxWithHole_tempvec[pos-1] ,Hole_box,0,
                                              G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_EmptySteelBoxWithHole)) );
      ScintilatorBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("ScintilatorBoxWithHole_")+std::to_string(pos)).c_str(),
                                                ScintilatorBoxWithHole_tempvec[pos-1] , Hole_sct, 0,
                                                G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_ScintilatorBoxWithHole)) );
    }
  }
  EmptySteelBoxWithHole=EmptySteelBoxWithHole_tempvec[WOM_coord_vec.size()-1];
  ScintilatorBoxWithHole = ScintilatorBoxWithHole_tempvec[WOM_coord_vec.size()-1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineLogicalVolumes()
{
  expHall_log = new G4LogicalVolume(expHall_box,air,"World",0,0,0);
  sipm_base_log = new G4LogicalVolume(sipm_base, steel,"sipm_base",0,0,0);
  ScintilatorBox_log = new G4LogicalVolume(ScintilatorBoxWithHole, LAB_PPO,"ScintilatorBoxLV",0,0,0);
  SteelBox_log = new G4LogicalVolume(EmptySteelBoxWithHole,steel,"Steelbox",0,0,0);
  Outer_tube_log = new G4LogicalVolume(Outer_tube, PMMA_side, "Outer_tubeLV");
  WOM_tube_log = new G4LogicalVolume(WOM_tube, PMMA_side, "WOM_tubeLV");
  Inner_tube_log = new G4LogicalVolume(Inner_tube, PMMA_side, "Inner_tubeLV");
  PMMA_Ring_log = new G4LogicalVolume(PMMA_Ring, PMMA_bottom, "PMMA_RingLV");
  PMMA_disk_log = new G4LogicalVolume(PMMA_disk, PMMA_side, "PMMA_diskLV");
  Air_gap1_log = new G4LogicalVolume(Air_gap1, air, "Air_gap1LV");
  Air_gap2_log = new G4LogicalVolume(Air_gap2, air, "Air_gap2LV");
  WLS_tube1_log = new G4LogicalVolume(WLS_tube1, Bis_MSB, "WLS1LV");
  WLS_tube2_log = new G4LogicalVolume(WLS_tube2, Bis_MSB, "WLS2LV");
  PMMA_Hat_log = new G4LogicalVolume(PMMA_Hat, PMMA_side, "PMMA_HatLV");
  Steel_Add_log = new G4LogicalVolume(SteelAdd, steel, "Steel_AddLV");
  Sct_Inside_log = new G4LogicalVolume(SctInside, LAB_PPO, "Sct_InsideLV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::ConstructVolumes()
{

  expHall_phys = new G4PVPlacement(0, G4ThreeVector(),expHall_log,"World",0,false,0);

  G4double delta_Z;
  G4RotationMatrix*RM1=new G4RotationMatrix(0*deg,0*deg,0*deg);
  G4double delta_Z_0 = SteelZ/2 + Thickness_Steel_Add - Length_2;  // upper surface of PMMA ring

  SteelBox_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),SteelBox_log,"SteelBox",expHall_log,false,0, intersect_check);
  
  delta_Z = (WallThickness_Z_Bottom - WallThickness_Z_Cover)/2;
  ScintilatorBox_phys = new G4PVPlacement(0,G4ThreeVector(0, 0, delta_Z),ScintilatorBox_log,"ScintilatorBoxPV",expHall_log,false,0, intersect_check);

//-------------------------------------------------------------------
  G4double delta_Z_sipm = delta_Z_0 + Thickness_WLS + Length_WOM+sipmbasewidth/2;
//-------------------------------------------------------------------

  // PMMA Staff
  // Outer_tube
  G4double delta_Z_Outer_tube = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat - OlengthOuter/2;
  //WOM tube
  G4double delta_Z_WOM = delta_Z_0 + Thickness_WLS + Length_WOM/2;
  //Inner tube
  G4double delta_Z_Inner_tube = delta_Z_0 - Thickness_Ring + Length_1/2;
  //PMMA Ring
  G4double delta_Z_PMMA_Ring = delta_Z_0 - Thickness_Ring + Thickness_Ring/2;
  //PMMA Disk
  G4double delta_Z_PMMA_Disk = delta_Z_0 - Thickness_Ring + Length_1 - Thickness_Disk/2;
  // Air gap 1
  G4double delta_Z_Air_gap_1 = delta_Z_0 + OlengthAG1/2;
  // Air gap 2
  G4double delta_Z_Air_gap_2 = delta_Z_0 + OlengthAG1/2;
  // WLS 1
  G4double delta_Z_WLS_1 = delta_Z_0 + Thickness_WLS + Length_WOM/2;
  // WLS  2
  G4double delta_Z_WLS_2 = delta_Z_0 + Thickness_WLS + Length_WOM/2;
  // PMMA Hat
  G4double delta_Z_PMMA_Hat = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat/2;
  // Additional Steel
  G4double delta_Z_Steel_Add = SteelZ/2 + Thickness_Steel_Add/2;
  // LAB&PPO inside tube
  G4double delta_Z_Sct_Inside = delta_Z_0 - Thickness_Ring + OlengthInside/2;

  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
  {
    sipm_base_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_sipm), sipm_base_log, "sipm_base", expHall_log, false, pos, intersect_check) );
    Outer_tube_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Outer_tube), Outer_tube_log, "Outer_tube", expHall_log, false, 0, intersect_check) );
    WOM_tube_phys_vect.push_back(   new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), WOM_tube_log, "WOM tube", expHall_log, false, 0, intersect_check) );
    Inner_tube_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Inner_tube), Inner_tube_log, "Inner_tube", expHall_log, false, 0, intersect_check) );
    PMMA_Ring_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Ring), PMMA_Ring_log, "PMMA_Ring", expHall_log, false, 0, intersect_check) );
    PMMA_Disk_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Disk), PMMA_disk_log, "PMMA_Disk", expHall_log, false, 0, intersect_check) );
    Air_gap_1_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gap_1), Air_gap1_log, "Air_gap1", expHall_log, false, 0, intersect_check) );
    Air_gap_2_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gap_2), Air_gap2_log, "Air_gap2", expHall_log, false, 0, intersect_check) );
    WLS_tube1_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WLS_1), WLS_tube1_log, "WLS1", expHall_log, false, 0, intersect_check) );
    WLS_tube2_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WLS_2), WLS_tube2_log, "WLS2", expHall_log, false, 0, intersect_check) );
    PMMA_Hat_phys_vect.push_back(   new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Hat), PMMA_Hat_log, "PMMA_Hat", expHall_log, false, 0, intersect_check) );
    Steel_Add_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Steel_Add), Steel_Add_log, "Steel_Add", expHall_log, false, 0, intersect_check) );
    Sct_Inside_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Sct_Inside), Sct_Inside_log, "Sct_Inside", expHall_log, false, 0, intersect_check) );
  }
}

void OpNoviceDetectorConstruction::DefineVisAttributes(){
    G4Color blue        = G4Color(0., 0., 1.);
    G4Color green       = G4Color(0., 1., 0.);
    G4Color red         = G4Color(1., 0., 0.);
    G4Color white       = G4Color(1., 1., 1.);
    G4Color white_trans = G4Color(1., 1., 1., 0.0);
    G4Color cyan        = G4Color(0., 1., 1., 0.3);
    G4Color magenda     = G4Color(1.,0.,1.);

    G4VisAttributes *worldVisAtt = new G4VisAttributes;
    worldVisAtt->SetColor(green);
    worldVisAtt->SetVisibility(true);
    worldVisAtt->SetForceSolid(false);
    expHall_log->SetVisAttributes(worldVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* OpNoviceDetectorConstruction::Construct()
{
  intersect_check = false; // global intersection check
  DefineMaterials();
  DefineMPTs();
  DefineSolids();
  DefineLogicalVolumes();
  ConstructVolumes();
  DefineSurfaces();

  return expHall_phys;
}
