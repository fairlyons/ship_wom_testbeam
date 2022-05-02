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
//#include "OpNoviceDetectorMessenger.hh"
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
//#include"SensitiveDetector.hh"
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

  sipmSize = 3*mm;
  sipmBaseThickness = 1.42*mm;
  sipmWindowThickness = sipmBaseThickness + 0.3*mm;

  SteelX = 80*cm;
  SteelY = 120*cm;
  SteelZ = 35*cm;

  WallThickness_XY = 10*mm;
  WallThickness_Z_Bottom = 20*mm;
  WallThickness_Z_Cover = 5*mm;
  SctX = SteelX - 2*WallThickness_XY;
  SctY = SteelY - 2*WallThickness_XY;
  SctZ = SteelZ - WallThickness_Z_Bottom - WallThickness_Z_Cover;

  //G4double delta_Z_0 = SteelZ/2 + Thickness_Steel_Add - Length_2;  // upper surface of PMMA ring

  //Additional_Length = -100*mm + 21*mm;

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
  Thickness_Gap = 4*mm;
  Thickness_Disk = 4*mm;
  Thickness_Hat = 6*mm;
  Thickness_Steel_Add = 15*mm;
  Length_WOM = 230*mm + Additional_Length;
  Thickness_WLS = 0.02*mm;
  Thickness_Reflect = 0.02*mm;

  delta_X = SteelX/2 - 91.5*mm;
  delta_Y = SteelY/2 - 91.5*mm;

  WOM_coord_vec = {{0., -SctY/6.}, {0., SctY/6.}};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction::~OpNoviceDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();

  steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  Si = nist->FindOrBuildMaterial("G4_Si");

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

  // Linear alkyl benzene (LAB)
  G4Element* C = new G4Element("Carbon", "C", 6 , 12.01*g/mole);
  G4Material* LAB = new G4Material("LAB",density=0.86*g/cm3,ncomponent=2);  //density: http://cpmaindia.com/lab_about.php
  LAB->AddElement(H,natoms=28);
  LAB->AddElement(C,natoms=17);
  // Diphenyloxazole (PPO)
  G4Material* PPO = new G4Material("PPO",density=1.184*g/cm3,ncomponent=4); //??
  PPO->AddElement(H,natoms=11);
  PPO->AddElement(C,natoms=15);
  PPO->AddElement(N,natoms=1);
  PPO->AddElement(O,natoms=1);
  // Scintilator (LAB+PPO) 23233 cm^3 23.233 l
  LAB_PPO = new G4Material("LAB_PPO", density=0.9*g/cm3, ncomponent=2); //??
  LAB_PPO->AddMaterial(LAB, 87.8*perCent); //?? Should be 2 g/L of PPO (First measurement of the surface tension of a liquid scintillator based on Linear Alkylbenzene (HYBLENE 113))
  LAB_PPO->AddMaterial(PPO, 12.2*perCent);
  // Bis-MSB WLS
  Bis_MSB = new G4Material("Bis_MSB",density=1.076*g/cm3,ncomponent=2); // density: http://www.molbase.com/moldata/368101.html
  Bis_MSB->AddElement(H,natoms=22);
  Bis_MSB->AddElement(C,natoms=24);
  // PEMA to build the WLS dye coat  
  PEMA = new G4Material("PEMA",density=1.11 * g / cm3,ncomponent=3); // density: https://polymerdatabase.com/polymers/polyethylmethacrylate.html
  PEMA->AddElement(H,natoms=10);
  PEMA->AddElement(C,natoms=6);
  PEMA->AddElement(O,natoms=242);
  // PTP (para-Terphenyl) to build the scintillator C18H14 and the dye coat
  PTP = new G4Material("PTP",density=1.23 * g / cm3,ncomponent=2); //  density: https://m.molbase.com/moldata/64879.html   
  PTP->AddElement(H,natoms=14);
  PTP->AddElement(C,natoms=18);
  /// WLS Coating  (150g PEMA, 3g PTP 1.5g bis-MSB)
  WLS_Coat = new G4Material("WLSCoat", density=1.1*g/cm3, ncomponent = 3); // density: same as PEMA 
  WLS_Coat->AddMaterial(Bis_MSB, 0.97 * perCent);
  WLS_Coat->AddMaterial(PTP,1.94 * perCent);
  WLS_Coat->AddMaterial(PEMA,  97.09 * perCent);
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
  auto Ridndex_LAB_PPO = [=](G4double wl) {
    G4double rind=1.;
    G4double B[4], C[4];
    B[0] =0.821384; C[0] =94.7625;
    B[1] =0.311375; C[1] =160.751;
    B[2] =0.0170099;C[2] =219.575;
    B[3] =0.608268; C[3] =9385.54;
    for(int term = 0; term<4; term++) rind+=B[term]/( 1.-(C[term]/wl)*(C[term]/wl) ); //formula eand coefficients: https://arxiv.org/pdf/1105.2101.pdf
    return sqrt(rind);
  };

  G4double wl;
  for(int i=0;i<100;i++) {
    wl = 250.+5.*i;
    photon_en_LAB_PPO[i]=1240./wl*eV;
    rindex_LAB_PPO[i]=Ridndex_LAB_PPO(wl);
  }
  G4MaterialPropertiesTable *MPT_LAB_PPO = new G4MaterialPropertiesTable();
  MPT_LAB_PPO -> AddConstProperty("SCINTILLATIONYIELD",10800./MeV); // https://underground.physics.berkeley.edu/WbLS/slides/PennRnD-Grullon.pdf
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
  for(int i=0; i<201; i++) photon_en_LAB_PPO_2[i] = 1240./photonWaveLength3[i]*eV;
  G4double scintilFast_LAB_PPO[201] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.045,
        0.09, 0.135, 0.18, 0.225, 0.27, 0.315, 0.36, 0.405, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.93, 0.91, 0.89, 0.87, 0.85, 0.83, 0.81,
        0.79, 0.77, 0.75, 0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, 0.77, 0.74, 0.71, 0.68, 0.65, 0.62, 0.59, 0.56, 0.53, 0.5, 0.4925, 0.485,
        0.4775, 0.47, 0.4625, 0.455, 0.4475, 0.44, 0.4325, 0.425, 0.41, 0.395, 0.38, 0.365, 0.35, 0.335, 0.32, 0.305, 0.29, 0.275, 0.2675, 0.26, 0.2525, 0.245, 0.2375,
        0.23, 0.2225, 0.215, 0.2075, 0.2, 0.1925, 0.185, 0.1775, 0.17, 0.1625, 0.155, 0.1475, 0.14, 0.1325, 0.125, 0.12, 0.115, 0.11, 0.105, 0.1, 0.095, 0.09, 0.085,
        0.08, 0.075, 0.0725, 0.07, 0.0675, 0.065, 0.0625, 0.06, 0.0575, 0.055, 0.0525, 0.05, 0.0475, 0.045, 0.0425, 0.04, 0.0375, 0.035, 0.0325, 0.03, 0.0275, 0.025,
        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.0225, 0.02, 0.0175, 0.015, 0.0125, 0.01, 0.0075, 0.005, 0.0025, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //https://arxiv.org/abs/1001.3946 
  MPT_LAB_PPO -> AddProperty("FASTCOMPONENT",photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201) -> SetSpline(true);
  MPT_LAB_PPO -> AddProperty("SLOWCOMPONENT",photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201) -> SetSpline(true);

  // transmission
  G4double waveLength[151] = {300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384
  ,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,476,478,480,482,484,486
  ,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,
  590,592,594,596,598,600};
  G4double photon_en_LAB_PPO_3[151];
  for(int i=0; i<211; i++) photon_en_LAB_PPO[i] = 1240./waveLength[i]*eV;
  G4double absLen_unpurified[151] = {0.01279493404851*m,0.0128029524293*m,0.012807751818825*m,0.012816426865377*m,0.01281835788177*m,0.012822230954434*m,0.012830840543912*m,0.012833994111742*m,0.012838615276908*m,0.012849574178154*m,0.01286687734122*m,0.012881105659852*m,0.012890383936614*m,0.012897247581221*m,0.012901565135242*m,0.01290710804677*m,0.012909213986899*m,0.012913052869886*m,0.01292213175909*m,0.012928356570193*m,0.01293558959276*m,0.01294674529108*m,0.012958145564748*m,0.012971261573831*m,0.012986890579723*m,0.013016803788927*m,0.013236885579806*m,0.017250128891369*m,0.028221561196635*m,0.048210452365814*m,0.085384830013475*m,0.153480082825598*m,0.265001590325018*m,0.406926671366136*m,0.541094630094606*m,0.663313263320605*m,0.779908411841287*m,0.875744366988952*m,0.927425529399054*m,0.956557242515996*m,0.9941373021031*m,1.04234271232283*m,1.1046585196898*m,1.12861772574266*m,1.12298134358669*m,1.21010948607421*m,1.40361593174363*m,1.59879984499012*m,1.74682662700547*m,1.85826095244245*m,1.94127546199575*m,2.05636621586457*m,2.23054267956921*m,2.41918635504621*m,2.58972782886322*m,2.74757561840523*m,2.90804104570209*m,3.06634990702638*m,3.219533855462*m,3.36858314280228*m,3.51899905487626*m,3.65200866136557*m,3.79194507519573*m,3.91998342801409*m,4.04589813769501*m,4.16757742933685*m,4.28519387842501*m,4.40871123138021*m,4.55010411037092*m,4.6681139833648*m,4.7917682338105*m,4.90159698959031*m,4.99959477462503*m,5.0911745929865*m,5.19943157180204*m,5.30550194154553*m,5.39615082475944*m,5.49731700891059*m,5.57289275228197*m,5.66214631723293*m,5.73506328468108*m,5.80521857408766*m,5.86604224769804*m,5.94477355531049*m,6.02864030249323*m,6.09293247274225*m,6.15871441142064*m,6.21946196680514*m,6.31098091183809*m,6.42934488959944*m,6.49136793799083*m,6.6118444804822*m,6.63947275377847*m,6.71587945264393*m,6.79526318509296*m,6.83526672561362*m,6.84525178280356*m,6.84830374947254*m,6.90597103165448*m,6.96176069278602*m,7.01469716585332*m,7.11844709300408*m,7.22642533541789*m,7.35800834707499*m,7.43721662692559*m,7.45713654579468*m,7.50711974244003*m,7.55677185990663*m,7.62664064911952*m,7.64365261429804*m,7.66923841145675*m,7.70325134674453*m,7.70586123507162*m,7.7082984257603*m,7.69320756579235*m,7.6061295282424*m,7.5800588175385*m,7.65415458551032*m,7.8726545001728*m,8.13213459072979*m,8.23755943796866*m,8.33114526507463*m,8.33232353493578*m,8.35648742634518*m,8.38645678568186*m,8.43190261832563*m,8.44405991017534*m,8.44097817045948*m,8.48935917676502*m,8.70096539741452*m,8.97482880285293*m,9.45328043541621*m,9.51335775244241*m,9.76515374992161*m,9.48924122323874*m,9.05389571096987*m,8.81555436133421*m,8.85221167762828*m,8.85772575154252*m,8.5724798759478*m,8.36896424710391*m,8.24046357622784*m,8.15949366895364*m,8.11488877012421*m,8.09059945466605*m,8.05778791764816*m,8.08548242289423*m,7.98944770316961*m,7.91287034263266*m,7.68678403498895*m,7.41928480365393*m}; // Patrick measurement Mainz TB_scintillator
  G4double absLen_purified[151] = {0.0128715069987*m,0.012878996045324*m,0.012882419482255*m,0.012889608969039*m,0.012897730665173*m,0.012901203446858*m,0.012901893060039*m,0.012909301261196*m,0.012913555972765*m,0.012919404279236*m,0.012941499851657*m,0.012956981482458*m,0.012962219216222*m,0.012968698731763*m,0.012975348745523*m,0.012978816961298*m,0.012984198598302*m,0.012992079547825*m,0.012996396487048*m,0.012998844829755*m,0.013008275438005*m,0.013019803591978*m,0.013030737069578*m,0.013043343934459*m,0.013063354253155*m,0.013096262643175*m,0.013310009583829*m,0.01737142379557*m,0.02880952553709*m,0.0503991725824*m,0.093267712946957*m,0.181782601946722*m,0.360363924864896*m,0.680589789294733*m,1.14377099881739*m,1.66966271660288*m,2.16375663283064*m,2.5867387338465*m,2.90168745463681*m,3.15134802038665*m,3.36485969966114*m,3.5472216422986*m,3.71820540936771*m,3.87499604590985*m,4.0416225005906*m,4.22693302308812*m,4.41961657539783*m,4.60195158864653*m,4.8012155457967*m,4.98340338198499*m,5.15238124063095*m,5.3394385508064*m,5.53766879317458*m,5.73334467669744*m,5.9379419460511*m,6.13214412943496*m,6.3463075559428*m,6.57082163294415*m,6.78144433425612*m,7.00057529008316*m,7.24933185585757*m,7.46193191606457*m,7.67375314752583*m,7.84805101504057*m,8.02793280531241*m,8.18529383430015*m,8.33199872471845*m,8.50086621234753*m,8.71252334840956*m,8.85479490604375*m,9.00886669573072*m,9.12179277182491*m,9.17977310076775*m,9.22629281746381*m,9.32376570267581*m,9.3663798468712*m,9.42391963661282*m,9.49917662439382*m,9.49239969544417*m,9.53291003815433*m,9.54447464256824*m,9.55647549660189*m,9.53375540026655*m,9.56114765877604*m,9.59797651082656*m,9.56671620921144*m,9.57416795211246*m,9.55705239801755*m,9.59985620139408*m,9.67791627886446*m,9.69593799132894*m,9.80518793240826*m,9.72163465573446*m,9.78343285609755*m,9.80672901676583*m,9.79744658014088*m,9.71997172370098*m,9.62140630852179*m,9.63069786946983*m,9.65338010651176*m,9.66256670125958*m,9.77530373579883*m,9.90904041225772*m,10.0761167690658*m,10.1350472968498*m,10.1133311242014*m,10.086029649153*m,10.0818415254723*m,10.1054177127325*m,10.0802588950859*m,10.0358961669191*m,10.0187720504539*m,9.9643151940587*m,9.91423094765397*m,9.79692373936646*m,9.62746317964499*m,9.54402290594051*m,9.61898800830966*m,9.92788014563769*m,10.231517114707*m,10.3365901053944*m,10.3860868021326*m,10.3153782242223*m,10.2925732890662*m,10.2662283365215*m,10.3319544221576*m,10.3275701553562*m,10.3082593229783*m,10.3800086540752*m,10.670267622928*m,10.9852633716743*m,11.7094782859565*m,11.6636028590624*m,11.9325107778151*m,11.4980811433371*m,10.8002040964469*m,10.4753853636081*m,10.5108807867562*m,10.4907062873286*m,10.1046947101091*m,9.82801139775702*m,9.60814130202656*m,9.48082544118919*m,9.36095897891052*m,9.32945396614867*m,9.31493016410715*m,9.37254639295624*m,9.2987694632999*m,9.27373907607242*m,9.02995755323176*m,8.69251407037729*m}; // Patrick measurement Mainz Column00+ppo
  MPT_LAB_PPO -> AddProperty("ABSLENGTH",photon_en_LAB_PPO,absLen_unpurified,151) -> SetSpline(true);
  MPT_LAB_PPO -> AddConstProperty("RESOLUTIONSCALE", 2.0); //??
  MPT_LAB_PPO -> AddConstProperty("FASTTIMECONSTANT", 5.2*ns); //??
  MPT_LAB_PPO -> AddConstProperty("SLOWTIMECONSTANT", 18.4*ns); //??
  MPT_LAB_PPO->AddConstProperty("YIELDRATIO",0.78);
  LAB_PPO -> SetMaterialPropertiesTable(MPT_LAB_PPO);



  



  //------------------------------------------------------------------------------
  //----------------------------- PMMA -----------------------------
  //------------------------------------------------------------------------------
   const G4int pmma_mpt_entr = 13;
  G4double pmma_side_wl[pmma_mpt_entr] = {700.,  600.,  550.,  500.,  450.,  400.,  390.,  380.,  370.,  350.,  320., 310.,  300. };
  G4double pmma_bottom_wl[75]= {441.465722452963,446.516468959472,451.567215465981,456.61796197249,461.668708478999,466.719454985508,471.770201492017,476.820947998526,481.871694505034,486.922441011543,
  491.973187518052,497.023934024561,500.238045437794,436.414975946454,431.364229439945,426.313482933437,421.262736426928,416.211989920419,411.16124341391,406.110496907401,401.518909174211,397.845638987659,
  396.238583281043,393.805041782452,392.152070198504,390.8926632774,389.39711755859,388.20330474796,387.055407814663,386.366669654684,385.172856844055,384.300455174749,383.382137628111,382.463820081473,
  381.545502534835,380.627184988197,379.662951564227,378.377306998934,377.183494188305,375.609235536925,372.821485841774,368.918636268563,363.867889762054,358.817143255545,353.766396749036,348.715650242527,
  343.664903736018,338.614157229509,333.563410723,328.512664216491,323.461917709983,318.411171203474,313.360424696965,308.309678190456,303.258931683947,298.361238101878,293.157438670929,288.10669216442,
  283.055945657911,278.005199151402,272.954452644894,267.903706138385,262.852959631876,257.802213125367,252.751466618858,247.892036267899,242.64997360584,237.599227099331,232.548480592822,227.497734086313,
  222.446987579804,217.396241073296,212.345494566787,207.294748060278,202.052685398219};
  G4double pmma_rind[pmma_mpt_entr] = {1.489, 1.492, 1.495, 1.498, 1.502, 1.511, 1.512, 1.514, 1.516, 1.522, 1.54, 1.541, 1.542}; // https://refractiveindex.info/?shelf=3d&book=plastics&page=pmma
  G4double pmma_side_en[pmma_mpt_entr];
  G4double pmma_bottom_en[75];
  for(int i=0; i<pmma_mpt_entr; i++ ) pmma_side_en[i]=1240./pmma_side_wl[i]*eV;
  for(int i=0; i<75; i++ ) pmma_bottom_en[i]=1240./pmma_bottom_wl[i]*eV;
  G4double pmma_side_abslen[pmma_mpt_entr] =  { 79.11*mm,70.39*mm,66.44*mm,61.48*mm,52.97*mm,45.77*mm,43.9*mm,42.6*mm,39.7*mm,36.07*mm,24.6*mm,18.23*mm,10.55*mm}; // picture in the folder 
  G4double pmma_bottom_abslen[75] = {3685.07275616407*mm,1975.90212068429*mm,1739.86542909011*mm,1556.52723104106*mm,1409.67117813625*mm,1289.12312649519*mm,1188.17858974627*mm,1102.24118950533*mm,1028.05540850485*mm,963.25152928162*mm,906.063978187016*mm,995.546039820757*mm,825.532796331376*mm,2739.34967914697*mm,3422.40435226653*mm,2615.11915536775*mm,1132.20893793201*mm,562.025169239236*mm,172.116453593441*mm,99.4344973522557*mm,42.1314492053959*mm,23.1587007862628*mm,14.5614893710837*mm,11.0809619119315*mm,8.65636778196572*mm,6.85307504536698*mm,5.55976428458039*mm,4.57033958442137*mm,3.85080777903292*mm,3.37910236494851*mm,2.91311721698663*mm,2.55268573906641*mm,2.27698546816599*mm,2.01668418096676*mm,1.79297764368291*mm,1.59752491160343*mm,1.4038711788544*mm,1.23546305784829*mm,1.06726603411345*mm,0.908835789380408*mm,0.740420832185596*mm,0.566454871207457*mm,0.42886782895817*mm,0.353271102463193*mm,0.307327896210654*mm,0.301885494993155*mm,0.301900827031845*mm,0.301917181987131*mm,0.301934641361397*mm,0.301953291111774*mm,0.301973221802586*mm,0.301994528780971*mm,0.30201731239245*mm,0.302041678261247*mm,0.302067737671548*mm,0.341435928448701*mm,0.302125413985435*mm,0.347533944235681*mm,0.437466744235865*mm,0.556750517609407*mm,0.647685273255554*mm,0.644020575189511*mm,0.535816480007306*mm,0.367918166877008*mm,0.302452281175179*mm,0.320791760065612*mm,0.3025669127656*mm,0.302631272369388*mm,0.302701282266752*mm,0.30277794773123*mm,0.302862615780021*mm,0.302957110157253*mm,0.303063920927376*mm,0.303186471853871*mm,0.321725431864253*mm}; // trasmittance: Transmission curves of plexiglass (PMMA) and optical grease - n: https://refractiveindex.info/?shelf=3d&book=plastics&page=pmma

  G4MaterialPropertiesTable *MPT_PMMA_side = new G4MaterialPropertiesTable();
  MPT_PMMA_side->AddProperty("RINDEX", pmma_side_en, pmma_rind, pmma_mpt_entr)->SetSpline(true);
  MPT_PMMA_side->AddProperty("ABSLENGTH", pmma_side_en, pmma_side_abslen, pmma_mpt_entr)->SetSpline(true);
  PMMA_side->SetMaterialPropertiesTable(MPT_PMMA_side);

  G4MaterialPropertiesTable *MPT_PMMA_bottom = new G4MaterialPropertiesTable();
  MPT_PMMA_bottom->AddProperty("RINDEX", pmma_side_en, pmma_rind, pmma_mpt_entr)->SetSpline(true);
  MPT_PMMA_bottom->AddProperty("ABSLENGTH", pmma_bottom_en, pmma_bottom_abslen,75);
  PMMA_bottom->SetMaterialPropertiesTable(MPT_PMMA_bottom);



  //------------------------------------------------------------------------------
  //----------------------------- WLS_Coat -----------------------------
  //------------------------------------------------------------------------------
  G4MaterialPropertiesTable *MPT_WLSCoat = new G4MaterialPropertiesTable();
  G4double waveLength3[136] = {475,474,473,472,471,470,469,468,467,466,465,464,463,462,461,460,459,458,457,456,455,454,453,452,451,450,449,448,447,446,445,444,443,442,441,440,439,438,437,436,435,434,433,432,431,430,429,428,427,426,425,424,
  423,422,421,420,419,418,417,416,415,414,413,412,411,410,409,408,407,406,405,404,403,402,401,400,399,398,397,396,395,394,393,392,391,390,389,388,387,386,385,384,383,382,381,380,379,378,377,376,375,374,373,372,
  371,370,369,368,367,366,365,364,363,362,361,360,359,358,357,356,355,354,353,352,351,350,349,348,347,346,345,344,343,342,341,340};
  G4double photonEnergy5[136];
  for(int i=0; i<7; i++) photonEnergy5[i] = 1240./waveLength3[i]*eV;
  G4double absLen3[136] = {0.029778641727802*mm,0.025727615980593*mm,0.022665520558419*mm,0.020403678900674*mm,0.018760693300556*mm,0.017577531521066*mm,0.016743287118922*mm,0.016167813685957*mm,0.015763998051761*mm,0.015488194203623*mm,0.015307536362709*mm,0.015191914667119*mm,0.015136095095359*mm,0.015132499525881*mm,0.015160770477501*mm,0.015211299915083*mm,0.01527461409535*mm,0.01533570681267*mm,0.015388969892783*mm,0.015446517621522*mm,0.015514103409318*mm,0.015604008404231*mm,0.015722238215546*mm,0.015842947046885*mm,0.015935673285256*mm,0.015961747565318*mm,0.015876626438453*mm,0.015663034252966*mm,0.0153368750864*mm,0.014921444513454*mm,0.014457211086833*mm,0.014005829236491*mm,0.013596300159085*mm,0.013236174984747*mm,0.012937598346848*mm,0.012704729003243*mm,0.012556254080098*mm,0.012492925033901*mm,0.012503072214779*mm,0.012580455720222*mm,0.012714392487393*mm,0.012884403875863*mm,0.013077390909031*mm,0.013283452954831*mm,0.013498835770703*mm,0.013723246482374*mm,0.01397419169843*mm,0.014241998652728*mm,0.014538624404617*mm,0.014867013478646*mm,0.015207837487861*mm,0.01555093737029*mm,0.015891392059806*mm,0.016191775144739*mm,0.016464263989425*mm,0.016696079973972*mm,0.016868032735573*mm,0.017004398964974*mm,0.017129512051334*mm,0.017253072087858*mm,0.017408464591978*mm,0.017599657075201*mm,0.017827956625514*mm,0.018107546799042*mm,0.018444575338237*mm,0.018811768412869*mm,0.019241846033599*mm,0.019713216602749*mm,0.02016609923033*mm,0.020604550905032*mm,0.021045484562316*mm,0.021456716664703*mm,0.021869584131306*mm,0.022271735128411*mm,0.022584232665911*mm,0.022783558437929*mm,0.022857134705945*mm,0.022745419121801*mm,0.022526778487644*mm,0.022194097391056*mm,0.021754927694089*mm,0.021275832295674*mm,0.020767552670763*mm,0.020200340029712*mm,0.019620327155833*mm,0.019025999593978*mm,0.018418473190498*mm,0.017867045375677*mm,0.017400073171183*mm,0.016956900957442*mm,0.016537832503563*mm,0.01611521038162*mm,0.015685271706403*mm,0.015297671648531*mm,0.014964140123143*mm,0.014666054270214*mm,0.014392594469716*mm,0.014081160494145*mm,0.013793257461952*mm,0.013492336322506*mm,0.013220777232175*mm,0.012963588085775*mm,0.012705359036573*mm,0.012441068320623*mm,0.012181936920422*mm,0.011918775570912*mm,0.011665080611527*mm,0.01142246436451*mm,0.011194560639647*mm,0.010993900870109*mm,0.010807618064254*mm,0.010635604013039*mm,0.010484281619043*mm,0.010359474572177*mm,0.010243379307035*mm,0.010147717929237*mm,0.010063195881395*mm,0.009986378514733*mm,0.009930972141261*mm,0.009887238724068*mm,0.009850512834273*mm,0.009819243355476*mm,0.009791168587455*mm,0.009754555980911*mm,0.009723489799155*mm,0.009690262669493*mm,0.009654553130722*mm,0.009619087536725*mm,0.009589732009137*mm,0.009559515170076*mm,0.009530682070846*mm,0.009502456734169*mm,0.009485526465539*mm,0.009468925475183*mm,0.009469233941075*mm,0.00948309435039*mm}; //?? need to be changed with the Jakobs one 
  MPT_WLSCoat->AddProperty("WLSABSLENGTH", photonEnergy5, absLen3, 136);

  // BIS reemission
  G4double waveLength4[16] = {380., 390., 400., 410., 420., 430., 440.,
                                            450., 460., 470., 480., 490., 500., 510., 520., 530};
  G4double photonEnergy6[16];
  for(int i=0; i<16; i++) photonEnergy6[i] = 1240./waveLength4[i]*eV;
  G4double reEmit4[16] = {0., 0., 0.1, 0.8, 1., 0.8, 0.5,
                                        0.45, 0.3, 0.2, 0.15, 0.1, 0.05, 0.05, 0.05, 0.}; //??
  G4double ppckovEmit[8] = { 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.6401*eV , 3.0402*eV , 3.5403*eV , 3.8404*eV};

  G4double rindexWLS[8] = { 1.5, 1.5, 1.5, 1.5, 1.504 , 1.505 , 1.515 , 1.52 };

  MPT_WLSCoat->AddProperty("WLSCOMPONENT", photonEnergy6, reEmit4, 16);
  MPT_WLSCoat->AddConstProperty("WLSTIMECONSTANT", 3.*ns); //??
  MPT_WLSCoat->AddProperty("RINDEX", pmma_side_en, pmma_rind, pmma_mpt_entr)->SetSpline(true); //??

  WLS_Coat->SetMaterialPropertiesTable(MPT_WLSCoat);


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
//    new G4LogicalBorderSurface( (std::string("WLSSurface_")+std::to_string(pos)).c_str(), WLS_tube1_phys_vect[pos], WOM_tube_phys_vect[pos], AirBoxSurface);
//    new G4LogicalBorderSurface( (std::string("WLSSurface_")+std::to_string(pos)).c_str(), WLS_tube2_phys_vect[pos], WOM_tube_phys_vect[pos], AirBoxSurface);
  }


  G4OpticalSurface* SipmWindowSurface = new G4OpticalSurface("SipmWindowSurface"); // WOM_tube -- sipmWindow border
  SipmWindowSurface->SetType(dielectric_dielectric);
  SipmWindowSurface->SetFinish(polished);
  SipmWindowSurface->SetModel(glisur);

  G4OpticalSurface* SipmBaseSurface = new G4OpticalSurface("SipmBaseSurface"); // sipmWindow -- sipmBase border
  SipmBaseSurface->SetType(dielectric_dielectric);
  SipmBaseSurface->SetFinish(polished);
  SipmBaseSurface->SetModel(glisur);
  G4double reflectivity2[num1] = {0.0999, 0.0999}; // by A. K.

  G4MaterialPropertiesTable *MPTsurf_SipmWindowSurface = new G4MaterialPropertiesTable();
  MPTsurf_SipmWindowSurface->AddProperty("REFLECTIVITY", pp1, reflectivity2, num1);
  SipmWindowSurface -> SetMaterialPropertiesTable(MPTsurf_SipmWindowSurface);

  for(unsigned int sipm_id = 0; sipm_id<80; sipm_id++)
  {
      new G4LogicalBorderSurface( (std::string("SipmWindowSurface_")+std::to_string(sipm_id)).c_str(),
                                  WOM_tube_phys_vect[sipm_id/40], sipm_phys_vect[sipm_id], SipmWindowSurface);
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
  for(int i=0; i<18; i++) photonEnergy7[i] = 1240./waveLength5[i]*eV;

//           G4double reflectSteel[18] = {0.430, 0.440, 0.449, 0.457, 0.463, 0.469, 0.474, 0.479,
//                          0.483, 0.487, 0.490, 0.493, 0.496, 0.499, 0.501, 0.504, 0.507, 0.508};    // from refractiveindex.info

//  G4double reflectSteel[18] = {0.467, 0.472, 0.477, 0.482, 0.486, 0.491, 0.496, 0.500,
//                 0.506, 0.510, 0.515, 0.519, 0.521, 0.522, 0.530, 0.538, 0.542, 0.546};    // from journal article

//  G4double reflectSteel[18] = {0.35, 0.36, 0.365, 0.37, 0.375, 0.38, 0.39, 0.395,
//                            0.4, 0.405, 0.41, 0.415, 0.42, 0.42, 0.425, 0.43, 0.43, 0.435};         // from borexino

  G4double reflectSteel[18] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                               0., 0., 0., 0., 0., 0., 0., 0.,};         // turn off the reflectivity

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

/*   for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
   {
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintilatorBox_phys, Outer_tube_phys_vect[pos], Scintillator_PMMA_Surface);
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintilatorBox_phys, PMMA_Ring_phys_vect[pos], Scintillator_PMMA_Surface);
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintilatorBox_phys, Inner_tube_phys_vect[pos], Scintillator_PMMA_Surface);
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), Sct_Inside_phys_vect[pos], Inner_tube_phys_vect[pos], Scintillator_PMMA_Surface);
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), Sct_Inside_phys_vect[pos], PMMA_Disk_phys_vect[pos], Scintillator_PMMA_Surface);
   } 
*/
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
  
  // Steel box
  std::vector<G4TwoVector> det(4);
  det[0].set(-402.848*mm, 642.936*mm);
  det[1].set(397.515*mm, 589.330*mm);
  det[3].set(-402.484*mm, -603.021*mm);
  det[2].set(397.515*mm, -635.084*mm);
  SteelBox = new G4ExtrudedSolid("Steelbox", det, SteelZ/2.,G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);
  //SteelBox = new G4Box("Steelbox",SteelX/2,SteelY/2,SteelZ/2);

  // Reflectivity box
  std::vector<G4TwoVector> reflect(4);
  reflect[0].set(-402.848*mm + WallThickness_XY, 642.936*mm - WallThickness_XY);
  reflect[1].set(397.515*mm - WallThickness_XY, 589.330*mm - WallThickness_XY);
  reflect[3].set(-402.484*mm + WallThickness_XY, -603.021*mm + WallThickness_XY);
  reflect[2].set(397.515*mm - WallThickness_XY, -635.084*mm + WallThickness_XY);
  ReflectBox = new G4ExtrudedSolid("ReflectBox", reflect, SctZ/2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);

  // Scintillator box
  std::vector<G4TwoVector> scint(4);
  scint[0].set(-402.848*mm + WallThickness_XY + Thickness_Reflect, 642.936*mm - WallThickness_XY - Thickness_Reflect);
  scint[1].set(397.515*mm - WallThickness_XY - Thickness_Reflect, 589.330*mm - WallThickness_XY - Thickness_Reflect);
  scint[3].set(-402.484*mm + WallThickness_XY + Thickness_Reflect, -603.021*mm + WallThickness_XY + Thickness_Reflect);
  scint[2].set(397.515*mm - WallThickness_XY - Thickness_Reflect, -635.084*mm + WallThickness_XY + Thickness_Reflect);
  ScintilatorBox = new G4ExtrudedSolid("ScintilatorBox", scint, (SctZ-Thickness_Reflect)/2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);
  //ScintilatorBox = new G4Box("ScintilatorBox",SctX/2,SctY/2,SctZ/2);

  G4double Rin, Rout;
  G4double delta_Z;

  // Outer tube
  Rin = Diam_Out_In/2;
  Rout = Diam_Out_Out/2;
  OlengthOuter = Length_2 + Thickness_Hat + Thickness_Ring + Thickness_Gap;
  Outer_tube  = new G4Tubs("Outer_tube", Rin, Rout, OlengthOuter/2, 0, 360*deg);
  // Air gap 1
  Rin = Diam_WOM_Out/2 + Thickness_WLS;
  Rout = Diam_Out_In/2;
  OlengthAG = Length_WOM + Thickness_WLS + Thickness_Gap;
  Air_gap1  = new G4Tubs("Air_gap1", Rin, Rout, OlengthAG/2, 0, 360*deg);
  // WLS tube 1
  Rin = Diam_WOM_Out/2;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WLS_tube1  = new G4Tubs("WLS_tube1", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // WOM tube
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WOM_tube  = new G4Tubs("WOM_tube", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // Hole in box
  Rin = 0.0*mm;
  Rout = Diam_Out_Out/2;
  double OlengthHoleBox = Length_2 + Thickness_Ring + Thickness_Gap - Thickness_Steel_Add;
  Hole_box  = new G4Tubs("Hole_box", Rin, Rout, OlengthHoleBox/2, 0, 360*deg);
  // Hole in scintilator
  Rin = 0*mm;
  Rout = Diam_Out_Out/2;
  double OlengthHoleSct = Length_2 + Thickness_Ring + Thickness_Gap - Thickness_Steel_Add - WallThickness_Z_Cover;
  Hole_sct = new G4Tubs("Hole_sct", Rin, Rout, OlengthHoleSct/2, 0, 360*deg);
  // WLS tube 2
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_In/2;
  WLS_tube2  = new G4Tubs("WLS_tube2", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // Air gap 2
  Rin = Diam_In_Out/2;
  Rout = Diam_WOM_In/2 - Thickness_WLS;
  Air_gap2  = new G4Tubs("Air_gap2", Rin, Rout, OlengthAG/2, 0, 360*deg);
  // Inner tube
  Rin = Diam_In_In/2;
  Rout = Diam_In_Out/2;
  Inner_tube  = new G4Tubs("Inner_tube", Rin, Rout, Length_1/2, 0, 360*deg);
  // PMMA ring
  Rin = Diam_In_Out/2;
  Rout = Diam_Out_In/2;
  PMMA_Ring  = new G4Tubs("PMMA_ring", Rin, Rout, Thickness_Ring/2, 0, 360*deg);
  // PMMA disk
  Rin = 0*mm;
  Rout = Diam_In_In/2;
  PMMA_disk  = new G4Tubs("PMMA_disk", Rin, Rout, Thickness_Disk/2, 0, 360*deg);
  // PMMA "hat"
  Rin = Diam_Out_Out/2;
  Rout = Diam_Hat/2;
  PMMA_Hat = new G4Tubs("PMMA_hat", Rin, Rout, Thickness_Hat/2, 0, 360*deg);
  // Additional steel
  Rin = Diam_Hole/2;
  Rout = Diam_Steel_Add/2;
  double OlengthAdd = Thickness_Steel_Add;
  SteelAdd = new G4Tubs("Steel_add", Rin, Rout, OlengthAdd/2, 0, 360*deg);
  // LAB&PPO inside tube
  Rin = 0.0*mm;
  Rout = Diam_In_In/2;
  OlengthInside = Length_1 - Thickness_Disk;
  SctInside = new G4Tubs("Sct_inside", Rin, Rout, OlengthInside/2, 0, 360*deg);
  // WLS ring
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  double OlengthWlsRing = Thickness_WLS;
  WLS_ring = new G4Tubs("WLS_ring ", Rin, Rout, OlengthWlsRing, 0, 360*deg);
  // Air ring outer
  Rin = Diam_WOM_Out/2 - 1*mm;
  Rout = Diam_WOM_Out/2;
  Air_ring1 = new G4Tubs("Air_ring1", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // PMMA ring supporting WOM
  Rin = Diam_WOM_In/2 + 1*mm;
  Rout = Diam_WOM_Out/2 - 1*mm;
  PMMA_ring_lower = new G4Tubs("PMMA_ring_lower", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // Air ring inner
  Rin = Diam_WOM_In/2;
  Rout = Diam_WOM_In/2 + 1*mm;
  Air_ring2 = new G4Tubs("Air_ring2", Rin, Rout, Thickness_Gap/2, 0, 360*deg);

  delta_Z = (WallThickness_Z_Bottom - WallThickness_Z_Cover)/2;
  G4SubtractionSolid *EmptySteelBox = new G4SubtractionSolid("EmptySteelBox",SteelBox,ReflectBox,0,G4ThreeVector(0,0,delta_Z));
  G4SubtractionSolid *EmptyReflectBox = new G4SubtractionSolid("EmptyReflectBox",ReflectBox,ScintillatorBox,0,G4ThreeVector(0,0,delta_Z));

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

  sipmBox = new G4Box("sipmBox", sipmSize/2., sipmSize/2., sipmWindowThickness/2.);
  sipmBaseBox = new G4Box("sipmBase", sipmSize/2., sipmSize/2., sipmBaseThickness/2.);
  WOM_cellBox = new G4Box("wom_cell", SteelX/2,SteelY/2,SteelZ/2 + 15*cm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineLogicalVolumes()
{
  expHall_log = new G4LogicalVolume(expHall_box,air,"World",0,0,0);
  WOM_cell_log = new G4LogicalVolume(WOM_cellBox,air,"wom_cell",0,0,0);
  sipm_base_log = new G4LogicalVolume(sipm_base, steel,"sipm_base",0,0,0);
  ScintilatorBox_log = new G4LogicalVolume(ScintilatorBoxWithHole, LAB_PPO,"ScintilatorBoxLV",0,0,0);
  SteelBox_log = new G4LogicalVolume(EmptySteelBoxWithHole,steel,"Steelbox",0,0,0);
  ReflectBox_log = new G4LogicalVolume(EmptyReflectBox, steel,"Reflectbox", 0,0,0);
  Outer_tube_log = new G4LogicalVolume(Outer_tube, PMMA_side, "Outer_tubeLV");
  WOM_tube_log = new G4LogicalVolume(WOM_tube, PMMA_bottom, "WOM_tubeLV");
  Inner_tube_log = new G4LogicalVolume(Inner_tube, PMMA_side, "Inner_tubeLV");
  PMMA_Ring_log = new G4LogicalVolume(PMMA_Ring, PMMA_bottom, "PMMA_RingLV");
  PMMA_disk_log = new G4LogicalVolume(PMMA_disk, PMMA_side, "PMMA_diskLV");
  Air_gap1_log = new G4LogicalVolume(Air_gap1, air, "Air_gap1LV");
  Air_gap2_log = new G4LogicalVolume(Air_gap2, air, "Air_gap2LV");
  WLS_tube1_log = new G4LogicalVolume(WLS_tube1, Bis_MSB, "WLS1LV");
  WLS_tube2_log = new G4LogicalVolume(WLS_tube2, Bis_MSB, "WLS2LV");
  PMMA_Hat_log = new G4LogicalVolume(PMMA_Hat, PMMA_side, "PMMA_HatLV");
  Air_ring1_log = new G4LogicalVolume(Air_ring1, air, "Air_ring1LV");
  Air_ring2_log = new G4LogicalVolume(Air_ring2, air, "Air_ring2LV");
  PMMA_ring_lower_log = new G4LogicalVolume(PMMA_ring_lower, PMMA_bottom, "PMMA_ring_lowerLV");
  Steel_Add_log = new G4LogicalVolume(SteelAdd, steel, "Steel_AddLV");
  Sct_Inside_log = new G4LogicalVolume(SctInside, LAB_PPO, "Sct_InsideLV");
  sipmBox_log = new G4LogicalVolume(sipmBox, PMMA_side, "sipmBox");
  sipmBaseBox_log = new G4LogicalVolume(sipmBaseBox, Si, "sipmBaseBox");
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

  SteelBox_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),SteelBox_log,"SteelBox",WOM_cell_log,false,0, intersect_check);
  delta_Z = (WallThickness_Z_Bottom - WallThickness_Z_Cover)/2;
  ScintilatorBox_phys = new G4PVPlacement(0,G4ThreeVector(0,0,delta_Z),ScintilatorBox_log,"ScintilatorBoxPV",WOM_cell_log,false,0, intersect_check);

  sipmBase_phys = new G4PVPlacement(0,G4ThreeVector(0, 0, sipmWindowThickness/2. - sipmBaseThickness/2. ),sipmBaseBox_log,"sipmBase",sipmBox_log,false,0, intersect_check);

//-------------------------------------------------------------------
  G4double delta_Z_sipm = delta_Z_0 + Thickness_WLS + Length_WOM + sipmWindowThickness/2;
//-------------------------------------------------------------------

  // PMMA Staff
  // Outer_tube
  G4double delta_Z_Outer_tube = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat - OlengthOuter/2;
  // WOM tube
  G4double delta_Z_WOM = delta_Z_0 + Thickness_WLS + Length_WOM/2;
  G4cout << "########   " << delta_Z_WOM << "########   " << G4endl;
  // Inner tube
  G4double delta_Z_Inner_tube = delta_Z_0 - Thickness_Ring + Length_1/2;
  // PMMA Ring
  G4double delta_Z_PMMA_Ring = delta_Z_0 - Thickness_Ring + Thickness_Ring/2 - Thickness_Gap;
  // PMMA Disk
  G4double delta_Z_PMMA_Disk = delta_Z_0 - Thickness_Ring + Length_1 - Thickness_Disk/2;
  // Air gap
  G4double delta_Z_Air_gap = delta_Z_0 + OlengthAG/2 - Thickness_Gap;
  // WLS 1
  G4double delta_Z_WLS_1 = delta_Z_0 + Thickness_WLS + Length_WOM/2;
  // WLS  2
  G4double delta_Z_WLS_2 = delta_Z_0 + Thickness_WLS + Length_WOM/2;
  // PMMA Hat
  G4double delta_Z_PMMA_Hat = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat/2;
  // Additional Steel
  G4double delta_Z_Steel_Add = SteelZ/2 + Thickness_Steel_Add/2;
  // LAB&PPO inside tube
  G4double delta_Z_Sct_Inside = delta_Z_0 - Thickness_Ring + OlengthInside/2 - Thickness_Gap;
  // Air ring
  G4double delta_Z_upper_ring = delta_Z_0 - Thickness_Ring + Thickness_Ring/2;

  G4int n_sipm = 40;
  G4double radius_sipm = (Diam_WOM_In + Diam_WOM_Out)/4.;
  G4int sipm_id = 0;

  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++) {
    RM1 = new G4RotationMatrix();
    //sipm_base_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_sipm), sipm_base_log, "sipm_base", WOM_cell_log, false, pos, intersect_check) );
    Outer_tube_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Outer_tube), Outer_tube_log, "Outer_tube", WOM_cell_log, false, 0, intersect_check) );
    WOM_tube_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), WOM_tube_log, "WOM tube", WOM_cell_log, false, 0, intersect_check) );
    Inner_tube_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Inner_tube), Inner_tube_log, "Inner_tube", WOM_cell_log, false, 0, intersect_check) );
    PMMA_Ring_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Ring), PMMA_Ring_log, "PMMA_Ring", WOM_cell_log, false, 0, intersect_check) );
    PMMA_Disk_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Disk), PMMA_disk_log, "PMMA_Disk", WOM_cell_log, false, 0, intersect_check) );
    Air_gap_1_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gap), Air_gap1_log, "Air_gap1", WOM_cell_log, false, 0, intersect_check) );
    Air_gap_2_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gap), Air_gap2_log, "Air_gap2", WOM_cell_log, false, 0, intersect_check) );
    WLS_tube1_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(0.,0., 0.), WLS_tube1_log, "WLS1", WOM_tube_log, false, 0, intersect_check) );
    WLS_tube2_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(0.,0., 0.), WLS_tube2_log, "WLS2", WOM_tube_log, false, 0, intersect_check) );
    PMMA_Hat_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Hat), PMMA_Hat_log, "PMMA_Hat", WOM_cell_log, false, 0, intersect_check) );
    Air_ring_1_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_upper_ring), Air_ring1_log, "Air_ring1", WOM_cell_log, false, 0, intersect_check) );
    Air_ring_2_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_upper_ring), Air_ring2_log, "Air_ring2", WOM_cell_log, false, 0, intersect_check) );
    PMMA_ring_lower_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_upper_ring), PMMA_ring_lower_log, "PMMA_ring_lower", WOM_cell_log, false, 0, intersect_check) );
    Steel_Add_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Steel_Add), Steel_Add_log, "Steel_Add", WOM_cell_log, false, 0, intersect_check) );
    Sct_Inside_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Sct_Inside), Sct_Inside_log, "Sct_Inside", WOM_cell_log, false, 0, intersect_check) );

    for(int i = 0; i<n_sipm; i++){
      RM1 = new G4RotationMatrix();
      RM1->rotateZ(- i * 360./n_sipm * deg);
      sipm_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first +  radius_sipm*std::cos(i*2*pi/n_sipm),
                                                                    WOM_coord_vec[pos].second + radius_sipm*std::sin(i*2*pi/n_sipm),
                                                                    delta_Z_sipm),
                                                                    sipmBox_log, "sipm", WOM_cell_log, false, sipm_id++, intersect_check) );
    }
  }

  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 0, intersect_check);
  //RM1 = new G4RotationMatrix();
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(-SteelX/2.,-SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 2, intersect_check) );
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(+SteelX/2.,+SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 0, intersect_check) );    
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(-SteelX/2.,+SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 1, intersect_check) );
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(+SteelX/2.,-SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 3, intersect_check) );     

/*
   for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++)
   {
     RM1 = new G4RotationMatrix();
     sipm_base_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_sipm), sipm_base_log, "sipm_base", expHall_log, false, pos, intersect_check) );
     Outer_tube_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Outer_tube), Outer_tube_log, "Outer_tube", expHall_log, false, 0, intersect_check) );
     WOM_tube_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), WOM_tube_log, "WOM tube", expHall_log, false, 0, intersect_check) );
     Inner_tube_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Inner_tube), Inner_tube_log, "Inner_tube", expHall_log, false, 0, intersect_check) );
     PMMA_Ring_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Ring), PMMA_Ring_log, "PMMA_Ring", expHall_log, false, 0, intersect_check) );
     PMMA_Disk_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Disk), PMMA_disk_log, "PMMA_Disk", expHall_log, false, 0, intersect_check) );
     Air_gap_1_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gap_1), Air_gap1_log, "Air_gap1", expHall_log, false, 0, intersect_check) );
     Air_gap_2_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gap_2), Air_gap2_log, "Air_gap2", expHall_log, false, 0, intersect_check) );
     WLS_tube1_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(0.,0., 0.), WLS_tube1_log, "WLS1", WOM_tube_log, false, 0, intersect_check) );
     WLS_tube2_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(0.,0., 0.), WLS_tube2_log, "WLS2", WOM_tube_log, false, 0, intersect_check) );
     PMMA_Hat_phys_vect.push_back(  new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Hat), PMMA_Hat_log, "PMMA_Hat", expHall_log, false, 0, intersect_check) );
     Steel_Add_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Steel_Add), Steel_Add_log, "Steel_Add", expHall_log, false, 0, intersect_check) );
     Sct_Inside_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Sct_Inside), Sct_Inside_log, "Sct_Inside", expHall_log, false, 0, intersect_check) );

     for(int i = 0; i<n_sipm; i++){
       RM1 = new G4RotationMatrix();
       RM1->rotateZ(- i * 360./n_sipm * deg);
       sipm_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first +  radius_sipm*std::cos(i*2*pi/n_sipm),
                                                                     WOM_coord_vec[pos].second + radius_sipm*std::sin(i*2*pi/n_sipm),
                                                                     delta_Z_sipm),
                                                                     sipmBox_log, "sipm", expHall_log, false, sipm_id++, intersect_check) );
     }
   }
*/
}

void OpNoviceDetectorConstruction::DefineVisAttributes()
{
  blue        = G4Color(0., 0., 1.);
  G4Color grey        = G4Color(0.3, 0.3, 0.3, 0.2);
  blue_trans  = G4Color(0., 0., 1., 0.3);
  green       = G4Color(0., 1., 0., 0.2);
  red         = G4Color(1., 0., 0., 0.2);
  white       = G4Color(1., 1., 1.);
  G4Color white_trans = G4Color(1., 1., 1., 0.2);
  cyan        = G4Color(0., 1., 1., 0.3);
  G4Color magenta     = G4Color(1.,0.,1., 0.3);

  G4VisAttributes *worldVisAtt = new G4VisAttributes;
  worldVisAtt->SetVisibility(false);
  expHall_log->SetVisAttributes(worldVisAtt);
  WOM_cell_log->SetVisAttributes(worldVisAtt);

  G4VisAttributes *steelBoxVisAtt = new G4VisAttributes;
  steelBoxVisAtt->SetVisibility(true);
  steelBoxVisAtt->SetColor(white_trans);
  SteelBox_log->SetVisAttributes(steelBoxVisAtt);
  Steel_Add_log->SetVisAttributes(steelBoxVisAtt);

  G4VisAttributes *sctBoxVisAtt = new G4VisAttributes;
  sctBoxVisAtt->SetColor(blue_trans);
  sctBoxVisAtt->SetVisibility(true);
  ScintilatorBox_log->SetVisAttributes(sctBoxVisAtt);
  Sct_Inside_log->SetVisAttributes(sctBoxVisAtt);

  G4VisAttributes *PMMAVisAtt = new G4VisAttributes;
  PMMAVisAtt->SetVisibility(true);
  PMMAVisAtt->SetColor(grey);
  PMMA_disk_log->SetVisAttributes(PMMAVisAtt);
  PMMA_Ring_log->SetVisAttributes(PMMAVisAtt);
  PMMA_Hat_log->SetVisAttributes(PMMAVisAtt);
  Outer_tube_log->SetVisAttributes(PMMAVisAtt);
  Inner_tube_log->SetVisAttributes(PMMAVisAtt);
  PMMA_ring_lower_log->SetVisAttributes(PMMAVisAtt);

  G4VisAttributes *airVisAtt = new G4VisAttributes;
  airVisAtt->SetColor(green);
  airVisAtt->SetVisibility(true);
  Air_gap1_log->SetVisAttributes(airVisAtt);
  Air_gap2_log->SetVisAttributes(airVisAtt);
  Air_ring1_log->SetVisAttributes(airVisAtt);
  Air_ring2_log->SetVisAttributes(airVisAtt);

  G4VisAttributes *WLSVisAtt = new G4VisAttributes;
  WLSVisAtt->SetColor(red);
  WLSVisAtt->SetVisibility(true);
  WLS_tube1_log->SetVisAttributes(WLSVisAtt);
  WLS_tube2_log->SetVisAttributes(WLSVisAtt);

  G4VisAttributes *WOMVisAtt = new G4VisAttributes;
  WOMVisAtt->SetColor(magenta);
  WOMVisAtt->SetVisibility(true);
  WOM_tube_log->SetVisAttributes(WOMVisAtt);

  G4VisAttributes *sipmVisAtt = new G4VisAttributes;
  sipmVisAtt->SetColor(grey);
  sipmVisAtt->SetVisibility(true);
  sipmBox_log->SetVisAttributes(sipmVisAtt);
  sipmBaseBox_log->SetVisAttributes(sipmVisAtt);
  sipm_base_log->SetVisAttributes(sipmVisAtt);
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
  DefineVisAttributes();
  return expHall_phys;
}
