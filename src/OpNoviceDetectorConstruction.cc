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
  SteelZ = 27*cm;

  WallThickness_XY = 10*mm;
  WallThickness_Z_Cover = 5*mm;
  Thickness_Reflect = 0.05*mm;
  ReflectZ = SteelZ - WallThickness_XY*2;
  SctX = SteelX - 2*WallThickness_XY;
  SctY = SteelY - 2*WallThickness_XY;
  SctZ = SteelZ - WallThickness_XY*2 - Thickness_Reflect*2;

  Diam_In_In = 44*mm;
  Diam_In_Out = 50*mm;
  Diam_Out_In = 64*mm;
  Diam_Out_Out = 70*mm;
  Diam_WOM_In = 53*mm;
  Diam_WOM_Out = 60*mm;
  Diam_Hole = 72*mm;
  Diam_Steel_Add = 120*mm;
  Diam_Hat = 120*mm;
  Thickness_Disk = 5*mm;
  Thickness_Ring = 4*mm;
  Thickness_Hat = 10*mm;
  Thickness_WLS = 0.02*mm;
  Thickness_Gap = 1*mm;
  Thickness_Steel_Add = 15*mm;
  Length_Out = 225*mm;
  Length_WOM = 230*mm;
  Length_In = 225*mm;

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
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);

  // Linear alkyl benzene (LAB)
  G4Element* H = new G4Element("Hydrogen", "H", 1 , 1.01*g/mole);
  G4Element* C = new G4Element("Carbon", "C", 6 , 12.01*g/mole);
  G4Material* LAB = new G4Material("LAB",density=0.856*g/cm3,ncomponent=2);  //density: https://www.knowde.com/stores/sasol/documents/101895
  LAB->AddElement(H,natoms=28);
  LAB->AddElement(C,natoms=17);
  // Diphenyloxazole (PPO)
  G4Material* PPO = new G4Material("PPO",density=1.128*g/cm3,ncomponent=4); // density: https://www.echemi.com/sds/24-diphenyloxazole-pid_Rock24446.html
  PPO->AddElement(H,natoms=11);
  PPO->AddElement(C,natoms=15);
  PPO->AddElement(N,natoms=1);
  PPO->AddElement(O,natoms=1);
  // Scintillator (LAB+PPO) 23233 cm^3 23.233 l
  LAB_PPO = new G4Material("LAB_PPO", density=0.9*g/cm3, ncomponent=2); //??
  LAB_PPO->AddMaterial(LAB, 99.77*perCent); //Calculated considering  be 2 g/L of PPO (First measurement of the surface tension of a liquid scintillator based on Linear Alkylbenzene (HYBLENE 113)). 2g of PPO every liter of LAB
  LAB_PPO->AddMaterial(PPO, 0.23*perCent);
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
  // Barium sulphate (BaSO4) Reflectivity coating
  G4Element* Ba = new G4Element("Barium", "Ba", 56 , 137.327*g/mole);
  G4Element* S = new G4Element("Sulphur", "S", 16 , 32.065*g/mole);
  BaSO4 = new G4Material("BaSO4",density=4.5*g/cm3,ncomponent=3); // https://www.chemeurope.com/en/encyclopedia/Barium_sulfate.html
  BaSO4->AddElement(Ba,natoms=1);
  BaSO4->AddElement(S,natoms=1);
  BaSO4->AddElement(O,natoms=4);
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
  G4double photonWaveLength3[201] = {500,499,498,497,496,495,494,493,492,491,490,489,488,487,486,485,484,483,482,481,480,479,478,477,476,475,474,473,472,471,470,469,468,467,466,465,464,463,462,461,460,459,458,457,456,455,
  454,453,452,451,450,449,448,447,446,445,444,443,442,441,440,439,438,437,436,435,434,433,432,431,430,429,428,427,426,425,424,423,422,421,420,419,418,417,416,415,414,413,412,411,410,409,408,
  407,406,405,404,403,402,401,400,399,398,397,396,395,394,393,392,391,390,389,388,387,386,385,384,383,382,381,380,379,378,377,376,375,374,373,372,371,370,369,368,367,366,365,364,363,362,361,
  360,359,358,357,356,355,354,353,352,351,350,349,348,347,346,345,344,343,342,341,340,339,338,337,336,335,334,333,332,331,330,329,328,327,326,325,324,323,322,321,320,319,318,317,316,315,314,
  313,312,311,310,309,308,307,306,305,304,303,302,301,300};
  G4double photon_en_LAB_PPO_2[201];
  for(int i=0; i<201; i++) photon_en_LAB_PPO_2[i] = 1240./photonWaveLength3[i]*eV;
  G4double scintilFast_LAB_PPO[201] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.0275,
  0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.0725,0.075,0.08,0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12,0.125,
  0.1325,0.14,0.1475,0.155,0.1625,0.17,0.1775,0.185,0.1925,0.2,0.2075,0.215,0.2225,0.23,0.2375,0.245,0.2525,0.26,0.2675,0.275,0.29,0.305,0.32,0.335,0.35,0.365,0.38,0.395,0.41,
  0.425,0.4325,0.44,0.4475,0.455,0.4625,0.47,0.4775,0.485,0.4925,0.5,0.53,0.56,0.59,0.62,0.65,0.68,0.71,0.74,0.77,0.8,0.795,0.79,0.785,0.78,0.775,0.77,0.765,0.76,0.755,0.75,0.77,
  0.79,0.81,0.83,0.85,0.87,0.89,0.91,0.93,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.405,0.36,0.315,0.27,0.225,0.18,0.135,0.09,0.045,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}; //https://arxiv.org/abs/1001.3946 
  MPT_LAB_PPO -> AddProperty("FASTCOMPONENT",photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201) -> SetSpline(true);
  MPT_LAB_PPO -> AddProperty("SLOWCOMPONENT",photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201) -> SetSpline(true);

  // transmission
  G4double waveLength[151] = {600,598,596,594,592,590,588,586,584,582,580,578,576,574,572,570,568,566,564,562,560,558,556,554,552,550,548,546,544,542,540,538,536,534,532,
  530,528,526,524,522,520,518,516,514,512,510,508,506,504,502,500,498,496,494,492,490,488,486,484,482,480,478,476,474,472,470,468,466,464,462,460,458,456,454,452,450,448,446,444,442,
  440,438,436,434,432,430,428,426,424,422,420,418,416,414,412,410,408,406,404,402,400,398,396,394,392,390,388,386,384,382,380,378,376,374,372,370,368,366,364,362,360,358,356,354,352,350,348,
  346,344,342,340,338,336,334,332,330,328,326,324,322,320,318,316,314,312,310,308,306,304,302,300};
  G4double photon_en_LAB_PPO_3[151];
  for(int i=0; i<211; i++) photon_en_LAB_PPO[i] = 1240./waveLength[i]*eV;
  G4double absLen_unpurified[151] = {7.41928480365393*m,7.68678403498895*m,7.91287034263266*m,7.98944770316961*m,8.08548242289423*m,8.05778791764816*m,8.09059945466605*m,8.11488877012421*m,8.15949366895364*m,8.24046357622784*m,8.36896424710391*m,8.5724798759478*m,8.85772575154252*m,8.85221167762828*m,8.81555436133421*m,9.05389571096987*m,9.48924122323874*m,9.76515374992161*m,9.51335775244241*m,9.45328043541621*m,8.97482880285293*m,8.70096539741452*m,8.48935917676502*m,8.44097817045948*m,8.44405991017534*m,8.43190261832563*m,8.38645678568186*m,8.35648742634518*m,8.33232353493578*m,8.33114526507463*m,8.23755943796866*m,8.13213459072979*m,7.8726545001728*m,7.65415458551032*m,7.5800588175385*m,7.6061295282424*m,7.69320756579235*m,7.7082984257603*m,7.70586123507162*m,7.70325134674453*m,7.66923841145675*m,7.64365261429804*m,7.62664064911952*m,7.55677185990663*m,7.50711974244003*m,7.45713654579468*m,7.43721662692559*m,7.35800834707499*m,7.22642533541789*m,7.11844709300408*m,7.01469716585332*m,6.96176069278602*m,6.90597103165448*m,6.84830374947254*m,6.84525178280356*m,6.83526672561362*m,6.79526318509296*m,6.71587945264393*m,6.63947275377847*m,6.6118444804822*m,6.49136793799083*m,6.42934488959944*m,6.31098091183809*m,6.21946196680514*m,6.15871441142064*m,6.09293247274225*m,6.02864030249323*m,5.94477355531049*m,5.86604224769804*m,5.80521857408766*m,5.73506328468108*m,5.66214631723293*m,5.57289275228197*m,5.49731700891059*m,5.39615082475944*m,5.30550194154553*m,5.19943157180204*m,5.0911745929865*m,4.99959477462503*m,4.90159698959031*m,4.7917682338105*m,4.6681139833648*m,4.55010411037092*m,4.40871123138021*m,4.28519387842501*m,4.16757742933685*m,4.04589813769501*m,3.91998342801409*m,3.79194507519573*m,3.65200866136557*m,3.51899905487626*m,3.36858314280228*m,3.219533855462*m,3.06634990702638*m,2.90804104570209*m,2.74757561840523*m,2.58972782886322*m,2.41918635504621*m,2.23054267956921*m,2.05636621586457*m,1.94127546199575*m,1.85826095244245*m,1.74682662700547*m,1.59879984499012*m,1.40361593174363*m,1.21010948607421*m,1.12298134358669*m,1.12861772574266*m,1.1046585196898*m,1.04234271232283*m,0.9941373021031*m,0.956557242515996*m,0.927425529399054*m,0.875744366988952*m,0.779908411841287*m,0.663313263320605*m,0.541094630094606*m,0.406926671366136*m,0.265001590325018*m,0.153480082825598*m,0.085384830013475*m,0.048210452365814*m,0.028221561196635*m,0.017250128891369*m,0.013236885579806*m,0.013016803788927*m,0.012986890579723*m,0.012971261573831*m,0.012958145564748*m,0.01294674529108*m,0.01293558959276*m,0.012928356570193*m,0.01292213175909*m,0.012913052869886*m,0.012909213986899*m,0.01290710804677*m,0.012901565135242*m,0.012897247581221*m,0.012890383936614*m,0.012881105659852*m,0.01286687734122*m,0.012849574178154*m,0.012838615276908*m,0.012833994111742*m,0.012830840543912*m,0.012822230954434*m,0.01281835788177*m,0.012816426865377*m,0.012807751818825*m,0.0128029524293*m,0.01279493404851*m}; // Patrick measurement Mainz TB_scintillator
  G4double absLen_purified[151] = {8.69251407037729*m,9.02995755323176*m,9.27373907607242*m,9.2987694632999*m,9.37254639295624*m,9.31493016410715*m,9.32945396614867*m,9.36095897891052*m,9.48082544118919*m,9.60814130202656*m,9.82801139775702*m,10.1046947101091*m,10.4907062873286*m,10.5108807867562*m,10.4753853636081*m,10.8002040964469*m,11.4980811433371*m,11.9325107778151*m,11.6636028590624*m,11.7094782859565*m,10.9852633716743*m,10.670267622928*m,10.3800086540752*m,10.3082593229783*m,10.3275701553562*m,10.3319544221576*m,10.2662283365215*m,10.2925732890662*m,10.3153782242223*m,10.3860868021326*m,10.3365901053944*m,10.231517114707*m,9.92788014563769*m,9.61898800830966*m,9.54402290594051*m,9.62746317964499*m,9.79692373936646*m,9.91423094765397*m,9.9643151940587*m,10.0187720504539*m,10.0358961669191*m,10.0802588950859*m,10.1054177127325*m,10.0818415254723*m,10.086029649153*m,10.1133311242014*m,10.1350472968498*m,10.0761167690658*m,9.90904041225772*m,9.77530373579883*m,9.66256670125958*m,9.65338010651176*m,9.63069786946983*m,9.62140630852179*m,9.71997172370098*m,9.79744658014088*m,9.80672901676583*m,9.78343285609755*m,9.72163465573446*m,9.80518793240826*m,9.69593799132894*m,9.67791627886446*m,9.59985620139408*m,9.55705239801755*m,9.57416795211246*m,9.56671620921144*m,9.59797651082656*m,9.56114765877604*m,9.53375540026655*m,9.55647549660189*m,9.54447464256824*m,9.53291003815433*m,9.49239969544417*m,9.49917662439382*m,9.42391963661282*m,9.3663798468712*m,9.32376570267581*m,9.22629281746381*m,9.17977310076775*m,9.12179277182491*m,9.00886669573072*m,8.85479490604375*m,8.71252334840956*m,8.50086621234753*m,8.33199872471845*m,8.18529383430015*m,8.02793280531241*m,7.84805101504057*m,7.67375314752583*m,7.46193191606457*m,7.24933185585757*m,7.00057529008316*m,6.78144433425612*m,6.57082163294415*m,6.3463075559428*m,6.13214412943496*m,5.9379419460511*m,5.73334467669744*m,5.53766879317458*m,5.3394385508064*m,5.15238124063095*m,4.98340338198499*m,4.8012155457967*m,4.60195158864653*m,4.41961657539783*m,4.22693302308812*m,4.0416225005906*m,3.87499604590985*m,3.71820540936771*m,3.5472216422986*m,3.36485969966114*m,3.15134802038665*m,2.90168745463681*m,2.5867387338465*m,2.16375663283064*m,1.66966271660288*m,1.14377099881739*m,0.680589789294733*m,0.360363924864896*m,0.181782601946722*m,0.093267712946957*m,0.0503991725824*m,0.02880952553709*m,0.01737142379557*m,0.013310009583829*m,0.013096262643175*m,0.013063354253155*m,0.013043343934459*m,0.013030737069578*m,0.013019803591978*m,0.013008275438005*m,0.012998844829755*m,0.012996396487048*m,0.012992079547825*m,0.012984198598302*m,0.012978816961298*m,0.012975348745523*m,0.012968698731763*m,0.012962219216222*m,0.012956981482458*m,0.012941499851657*m,0.012919404279236*m,0.012913555972765*m,0.012909301261196*m,0.012901893060039*m,0.012901203446858*m,0.012897730665173*m,0.012889608969039*m,0.012882419482255*m,0.012878996045324*m,0.0128715069987*m}; // Patrick measurement Mainz Column00+ppo
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
  G4double pmma_bottom_wl[75]= {500.238045437794,497.023934024561,491.973187518052,486.922441011543,481.871694505034,476.820947998526,471.770201492017,466.719454985508,461.668708478999,456.61796197249,451.567215465981,446.516468959472,
  441.465722452963,436.414975946454,431.364229439945,426.313482933437,421.262736426928,416.211989920419,411.16124341391,406.110496907401,401.518909174211,397.845638987659,396.238583281043,393.805041782452,
  392.152070198504,390.8926632774,389.39711755859,388.20330474796,387.055407814663,386.366669654684,385.172856844055,384.300455174749,383.382137628111,382.463820081473,381.545502534835,380.627184988197,
  379.662951564227,378.377306998934,377.183494188305,375.609235536925,372.821485841774,368.918636268563,363.867889762054,358.817143255545,353.766396749036,348.715650242527,343.664903736018,338.614157229509,
  333.563410723,328.512664216491,323.461917709983,318.411171203474,313.360424696965,308.309678190456,303.258931683947,298.361238101878,293.157438670929,288.10669216442,283.055945657911,278.005199151402,
  272.954452644894,267.903706138385,262.852959631876,257.802213125367,252.751466618858,247.892036267899,242.64997360584,237.599227099331,232.548480592822,227.497734086313,222.446987579804,217.396241073296,
  212.345494566787,207.294748060278,202.052685398219};
  G4double pmma_rind[pmma_mpt_entr] = {1.489, 1.492, 1.495, 1.498, 1.502, 1.511, 1.512, 1.514, 1.516, 1.522, 1.54, 1.541, 1.542}; // https://refractiveindex.info/?shelf=3d&book=plastics&page=pmma
  G4double pmma_side_en[pmma_mpt_entr];
  G4double pmma_bottom_en[75];
  for(int i=0; i<pmma_mpt_entr; i++ ) pmma_side_en[i]=1240./pmma_side_wl[i]*eV;
  for(int i=0; i<75; i++ ) pmma_bottom_en[i]=1240./pmma_bottom_wl[i]*eV;
  G4double pmma_side_abslen[pmma_mpt_entr] =  {79.11*mm,70.39*mm,66.44*mm,61.48*mm,52.97*mm,45.77*mm,43.9*mm,42.6*mm,39.7*mm,36.07*mm,24.6*mm,18.23*mm,10.55*mm}; // picture in the folder 
  G4double pmma_bottom_abslen[75] = {825.532796331375*mm,995.546039820757*mm,906.063978187016*mm,963.25152928162*mm,1028.05540850485*mm,1102.24118950533*mm,1188.17858974627*mm,1289.12312649519*mm,1409.67117813625*mm,1556.52723104106*mm,1739.86542909011*mm,1975.90212068429*mm,3685.07275616407*mm,2739.34967914697*mm,3422.40435226653*mm,2615.11915536775*mm,1132.20893793201*mm,562.025169239236*mm,172.116453593441*mm,99.4344973522557*mm,42.1314492053959*mm,23.1587007862628*mm,14.5614893710837*mm,11.0809619119315*mm,8.65636778196572*mm,6.85307504536698*mm,5.55976428458039*mm,4.57033958442137*mm,3.85080777903292*mm,3.37910236494851*mm,2.91311721698663*mm,2.55268573906641*mm,2.27698546816599*mm,2.01668418096676*mm,1.79297764368291*mm,1.59752491160343*mm,1.4038711788544*mm,1.23546305784829*mm,1.06726603411345*mm,0.908835789380408*mm,0.740420832185596*mm,0.566454871207457*mm,0.42886782895817*mm,0.353271102463193*mm,0.307327896210654*mm,0.301885494993155*mm,0.301900827031845*mm,0.301917181987131*mm,0.301934641361397*mm,0.301953291111774*mm,0.301973221802586*mm,0.301994528780971*mm,0.30201731239245*mm,0.302041678261247*mm,0.302067737671548*mm,0.341435928448701*mm,0.302125413985435*mm,0.347533944235681*mm,0.437466744235865*mm,0.556750517609407*mm,0.647685273255554*mm,0.644020575189511*mm,0.535816480007306*mm,0.367918166877008*mm,0.302452281175179*mm,0.320791760065612*mm,0.3025669127656*mm,0.302631272369388*mm,0.302701282266752*mm,0.30277794773123*mm,0.302862615780021*mm,0.302957110157253*mm,0.303063920927376*mm,0.303186471853871*mm,0.321725431864253*mm}; //measure from CheapCal simulation / Doramas and Andrew measurements 
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
  G4double absLen3[136] = {0.029778641727802*mm,0.025727615980593*mm,0.022665520558419*mm,0.020403678900674*mm,0.018760693300556*mm,0.017577531521066*mm,0.016743287118922*mm,0.016167813685957*mm,0.015763998051761*mm,0.015488194203623*mm,0.015307536362709*mm,0.015191914667119*mm,0.015136095095359*mm,0.015132499525881*mm,0.015160770477501*mm,0.015211299915083*mm,0.01527461409535*mm,0.01533570681267*mm,0.015388969892783*mm,0.015446517621522*mm,0.015514103409318*mm,0.015604008404231*mm,0.015722238215546*mm,0.015842947046885*mm,0.015935673285256*mm,0.015961747565318*mm,0.015876626438453*mm,0.015663034252966*mm,0.0153368750864*mm,0.014921444513454*mm,0.014457211086833*mm,0.014005829236491*mm,0.013596300159085*mm,0.013236174984747*mm,0.012937598346848*mm,0.012704729003243*mm,0.012556254080098*mm,0.012492925033901*mm,0.012503072214779*mm,0.012580455720222*mm,0.012714392487393*mm,0.012884403875863*mm,0.013077390909031*mm,0.013283452954831*mm,0.013498835770703*mm,0.013723246482374*mm,0.01397419169843*mm,0.014241998652728*mm,0.014538624404617*mm,0.014867013478646*mm,0.015207837487861*mm,0.01555093737029*mm,0.015891392059806*mm,0.016191775144739*mm,0.016464263989425*mm,0.016696079973972*mm,0.016868032735573*mm,0.017004398964974*mm,0.017129512051334*mm,0.017253072087858*mm,0.017408464591978*mm,0.017599657075201*mm,0.017827956625514*mm,0.018107546799042*mm,0.018444575338237*mm,0.018811768412869*mm,0.019241846033599*mm,0.019713216602749*mm,0.02016609923033*mm,0.020604550905032*mm,0.021045484562316*mm,0.021456716664703*mm,0.021869584131306*mm,0.022271735128411*mm,0.022584232665911*mm,0.022783558437929*mm,0.022857134705945*mm,0.022745419121801*mm,0.022526778487644*mm,0.022194097391056*mm,0.021754927694089*mm,0.021275832295674*mm,0.020767552670763*mm,0.020200340029712*mm,0.019620327155833*mm,0.019025999593978*mm,0.018418473190498*mm,0.017867045375677*mm,0.017400073171183*mm,0.016956900957442*mm,0.016537832503563*mm,0.01611521038162*mm,0.015685271706403*mm,0.015297671648531*mm,0.014964140123143*mm,0.014666054270214*mm,0.014392594469716*mm,0.014081160494145*mm,0.013793257461952*mm,0.013492336322506*mm,0.013220777232175*mm,0.012963588085775*mm,0.012705359036573*mm,0.012441068320623*mm,0.012181936920422*mm,0.011918775570912*mm,0.011665080611527*mm,0.01142246436451*mm,0.011194560639647*mm,0.010993900870109*mm,0.010807618064254*mm,0.010635604013039*mm,0.010484281619043*mm,0.010359474572177*mm,0.010243379307035*mm,0.010147717929237*mm,0.010063195881395*mm,0.009986378514733*mm,0.009930972141261*mm,0.009887238724068*mm,0.009850512834273*mm,0.009819243355476*mm,0.009791168587455*mm,0.009754555980911*mm,0.009723489799155*mm,0.009690262669493*mm,0.009654553130722*mm,0.009619087536725*mm,0.009589732009137*mm,0.009559515170076*mm,0.009530682070846*mm,0.009502456734169*mm,0.009485526465539*mm,0.009468925475183*mm,0.009469233941075*mm,0.00948309435039*mm}; //Jakobs measure without cloroform  
  MPT_WLSCoat->AddProperty("WLSABSLENGTH", photonEnergy5, absLen3, 136);

  // BIS reemission
  G4double waveLength4[51] = {498.50394264782,495.512169509067,489.3462170527,484.343493743363,479.225925762846,472.348366539066,467.849074868481,464.04990973537,461.338683704906,459.793161903923,457.56412302422,454.795680527534,
  451.768205666593,446.867663596341,441.763591710159,435.933166703427,433.243802562634,431.985904788012,429.963943988042,427.926782673168,426.395623114075,425.126769876294,422.410632741153,420.384395226982,
  417.93569195544,415.411789507185,413.448584543676,412.387719607238,411.249479281543,409.9604817049,406.3089221806,404.020506020059,401.848223327902,399.264457914078,397.452088829242,394.82668707538,
  393.620389802785,392.661003458199,390.948731712797,389.702126722733,388.693734027439,388.169942931827,387.606839417529,386.814987134474,386.022324288337,384.176496975521,382.716212740349,381.014431692937,
  378.83045463701,376.455517150279,373.305908776801};  
  G4double photonEnergy6[51];
  for(int i=0; i<16; i++) photonEnergy6[i] = 1240./waveLength4[i]*eV;
  G4double reEmit4[51] = {0.095367125872878,0.100714589065125,0.13933792816855,0.176288041362039,0.211928140548633,0.245505005810625,0.273370328546945,0.30221148342213,0.333093045311687,0.364223283871252,0.399221677160171,
  0.448451857613,0.514706458856799,0.566173063784637,0.587186163465084,0.619912730325973,0.657430320109198,0.6974131402047,0.747799923939089,0.82411878610154,0.886968238152526,0.927207898674191,0.972182090856316,
  0.995579234732068,0.973067015983444,0.930457495931142,0.880606822664751,0.837907842014882,0.808317266222523,0.778775703403637,0.739211827896116,0.774816618055286,0.818288494562106,0.86029367303487,
  0.836373519169372,0.783204456572165,0.716723640486035,0.653897573842595,0.600839623111776,0.539529302980934,0.482274276482319,0.45143792246339,0.403951196810372,0.355312021195194,0.308055666198011,
  0.257497891356571,0.199434040798557,0.131499952640494,0.074850729151259,0.038245233587249,0.007817590412588}; // (jakob) https://omlc.org/spectra/PhotochemCAD/html/044.html
  G4double ppckovEmit[8] = { 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.6401*eV , 3.0402*eV , 3.5403*eV , 3.8404*eV}; 

  G4double rindexWLS[8] = { 1.5, 1.5, 1.5, 1.5, 1.504 , 1.505 , 1.515 , 1.52 };

  MPT_WLSCoat->AddProperty("WLSCOMPONENT", photonEnergy6, reEmit4, 51);
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

  for(unsigned int sipm_id = 0; sipm_id<80; sipm_id++) new G4LogicalBorderSurface( (std::string("SipmWindowSurface_")+std::to_string(sipm_id)).c_str(),
                                  WOM_tube_phys_vect[sipm_id/40], sipm_phys_vect[sipm_id], SipmWindowSurface);

  //------------------------------------------------------------------------------
  //----------------------------- Steel -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* SteelBoxSurface = new G4OpticalSurface("SteelBoxSurface");
  SteelBoxSurface -> SetType(dielectric_metal);
  SteelBoxSurface -> SetFinish(ground);
  SteelBoxSurface -> SetModel(glisur);

  G4double waveLength5[59] = {838.887127914013,828.817393308936,818.112266950984,808.64910814627,797.918047431058,
  787.812292885341,778.327522116244,768.231853153906,758.120335417687,747.343169178457,737.87136558799,727.743999077889,
  717.623836555917,708.13906578682,698.655735815347,687.914589516756,678.422614759531,668.301011439933,657.604529867735,
  648.073653574619,637.329625680777,627.859262887935,617.112353398842,609.530876293116,598.147134253523,588.662363484425,
  578.540760164827,568.414834452353,558.921418897502,548.170187015532,538.05002449356,527.928421173962,517.811140247241,
  507.062789960522,497.578019191425,487.460738264703,477.336253349855,467.217531625508,457.728438463534,447.612598334438,
  439.389966285336,428.007665043368,418.520012679019,408.401290954673,398.2811284327,388.155202720226,378.599832867476,
  367.563323055644,358.664956920151,347.905080252428,337.790680920958,328.027115811313,318.290925856554,
  308.234158430107,298.119759098637,288.617698758033,278.498977033686,269.063193383857,257.757254416045};
  
  G4double photonEnergy7[59];
  for(int i=0; i<18; i++) photonEnergy7[i] = 1240./waveLength5[i]*eV;


  G4double specular_steel[59] = {0.30216675359479,0.255752492579892,0.211988869838342,0.189438738725759,0.169538407342305,
  0.156267606546984,0.153603551566147,0.131052581909345,0.123084734749157,0.145608032446762,0.131012331786877,
  0.137627607123304,0.137614190415814,0.134950135434977,0.130960342045356,0.120340179523384,0.124304816586462,
  0.125617138287756,0.073899085093508,0.11365865919373,0.105689973489324,0.089768534420655,0.084451325533816,0.0764868325505,
 0.075146000345791,0.072481945364954,0.073794267066248,0.079083803993891,0.084374179465752,0.083034185805262,0.083020769097773,
 0.084333090799066,0.081668197274011,0.077676726795954,0.075012671815117,0.072347778290061,0.076311576808922,0.074972421692649,
  0.076285581938161,0.072294950004322,0.072284048929487,0.069617478315996,0.069604900152724,0.068265745036452,  0.068252328328963,0.071960349606668,0.121685701224445,0.28143546343462,0.268712935986641,0.224746728193772,0.160625740904492,
  0.109660495475257,0.060320585725261,0.031405124642217,0.02550943578365,0.03241234077211,0.039010130683933,0.029993785627795,
  0.019391300509685};  //specular reflectivity Patrick 
  
  G4double other_steel[59] = {0,0,0,0,0,2.6712738969614,6.99912455322394,11.461179103357,13.3303728143136,
  7.76049324144496,10.5470832286767,8.48032079100203,7.89613039576159,9.83548465381472,11.2382332487584,12.9694290049434,
 11.7364909186509,12.943617756497,38.1908081955548,16.6488888058538,19.1187061362638,26.7334655790205,29.7746096076586,        35.7572000618953,36.8950525383479,38.834406796401,40.2088285102441,36.333977297541,31.4552729744343,31.7565672164803,
  31.9252037632262,31.9612664690937,32.7296404495901,35.6382106373499,37.0756802674122,39.1824132558841,36.276610264044,
  36.5778206516682,36.27920975112,39.8568755884458,39.8579656959293,41.7975715172479,42.4680088375628,43.7729884811686,
  45.6145737878838,46.0802460400981,42.111480134302,27.1402731692661,29.2490002940487,35.1512749573081,43.7382070741963,
  50.8422701290832,56.9473252360614,61.6791149803321,63.7743377501613,65.5934703912693,       65.7701657800718,68.177454169658,70.9106514414385}; // neither diffusive or specular reflectivity Patrick 

  G4double photonEnergy9[2] = {1.*eV,5.*eV}; 
  G4double specular_steel2[2] = {0,0}; // this is for not smooth surfaces. In this case probably the specular reflectivity is some of this and some of the other but we cannot know, so I would 						say that we can define just one
  G4MaterialPropertiesTable *MPTsurf_Steel = new G4MaterialPropertiesTable();
  MPTsurf_Steel->AddProperty("SPECULARLOBECONSTANT", photonEnergy7, specular_steel,59); // In order to have diffuse 												reclectivity (Lambertian), it is 													necessary define all the other 												three. The diffuse is 
   MPTsurf_Steel->AddProperty("SPECULARSPIKECONSTANT", photonEnergy9, specular_steel2,2); // 1-other three (in this case 1). 
   MPTsurf_Steel->AddProperty("BACKSCATTERCONSTANT", photonEnergy7, other_steel,59);
  //MPTsurf_Steel -> AddProperty("REFLECTIVITY", photonEnergy7, reflectSteel, 18);
  SteelBoxSurface -> SetMaterialPropertiesTable(MPTsurf_Steel);

  G4LogicalSkinSurface* Surface = new G4LogicalSkinSurface("SteelSurface",SteelBox_log,SteelBoxSurface); //I think this should be SteelBox_log instead of ReflectBox_log (ale)

  //------------------------------------------------------------------------------
  //----------------------------- BaSO4 (reflective coating) -----------------------------
  //------------------------------------------------------------------------------
   G4OpticalSurface* BaSO4_surface = new G4OpticalSurface("BaSO4_surface");
   BaSO4_surface -> SetType(dielectric_metal);
   BaSO4_surface -> SetFinish(ground);
   BaSO4_surface -> SetModel(glisur);
  
   G4double waveLength6[59] = {838.566997868485,826.743136751539,817.005392482906,807.237260496162,798.17412356959,
   788.399914039223,776.524393801489,767.462776260823,757.682489186835,747.874853166547,738.103682407991,
   726.914924599563,718.526395015053,707.343714750247,697.564947062164,687.789218145893,
   677.303936011731,668.216488910671,657.787424055015,647.261118501404,
   636.786472068582,627.715738212483,618.614616638273,607.456246547955,599.772712023618,587.1891579539,578.107788396462,
   567.616428718679,557.825505943352,548.040660711647,537.55385919158,
   528.464892704614,516.58937246688,508.203881654181,497.01968200347,487.941351217843,477.448472154154,466.265791889348,
   459.972495468584,447.398057714299,437.613212482594,427.828367250889,
   417.344604502634,408.257157401573,398.469273398057,388.682908780446,377.440972465324,368.095229760321,358.26480295145,
   347.753691256894,337.962768481567,326.80439839125,317.776207340506,
   308.044540615495,296.855782807067,287.740986759707,277.243549538301,267.499727726046,257.806045648673}; 
   
  G4double photonEnergy8[59];
  for(int i=0; i<18; i++) photonEnergy8[i] = 1240./waveLength6[i]*eV;
  
   G4double specular_coating[59] = { 0.043090222881352,0.031689546391218,0.041461554811332,0.033318214461237,0.02843221025118,0.036575550601275,0.033318214461237,0.038204218671295,
   0.030060878321199,0.038204218671295,0.021917537971104,0.039832886741313,0.033453755370293,0.020371406492993,0.015465525664007,0.033453755370293,0.025277287321981,0.030183168150969,
   0.030183168150969,0.030183168150969,0.025277287321981,0.022006700102656,0,0.031839947655533,0.031839947655533,0.020385162591644,0.031839947655533,0.030203549789264,0.022021560457915,
   0.031839947655533,0.020385162591644,0.025294356190454,0.018748764725372,0.026985273431965,0.022066141523691,0.022066141523691,0.026985273431965,0.030264694704149,0.020426430887599,
   0.020426430887599,0.015507298979326,0.026985273431965,0.01714700961542,0.023721816771648,0.025362631664347,0.02208100187895,0.025362631664347,0.026730849680839,0.03634062423262,
   0.037564094688637,0.035970028431757,0.034756854425251,0.027461227079687,0.020334977012879,0.025054734763179,0.045507018347822,0.045507018347822,0.026627987346615,0.026627987346615}; // specular reflection Patrick data, calculated as 1- diffuse reflection since the specular was very small and the measurements for specular are more difficult to separete from the diffuse ones.

   G4double other_coating[2] = {0,0}; // it is not relevant in our case, it is for not smooth surfaces 
   
   G4double p_coating_refl[19] = { 0.653775342148272*eV,0.687758866578874*eV,0.729108650377501*eV,0.773110023703728*eV,0.827676607343705*eV,0.887209379832684*eV,0.954680932864894*eV,1.03471853868547*eV,1.12204531460793*eV,1.24533694371132*eV,1.3767339264975*eV,1.55649354301667*eV,1.77189963551006*eV,2.09375317772178*eV,2.51208306681883*eV,3.16797881984726*eV,3.50425046717685*eV,3.63215740198218*eV,3.81157924082758*eV};
   
   G4double refl_coating[19] = { 0.861320132013201,0.916105610561056,0.922706270627063,0.920726072607261,0.916765676567657,
   0.931287128712871,0.95042904290429,0.958349834983498,0.962970297029703,0.968250825082508,0.969570957095709,
   0.973531353135314,0.977491749174918,0.978151815181518,0.98013201320132,0.980792079207921,0.970891089108911,
   0.960990099009901,0.953729372937294}; // reflectivity of the coating https://www.optopolymer.de/produktuebersicht/diffuse-reflecting-materials/bariumsulfate-baso4-coating-oprc/
       
   G4MaterialPropertiesTable* MTP_BaSO = new G4MaterialPropertiesTable();
   MTP_BaSO->AddProperty("SPECULARLOBECONSTANT", photonEnergy8, specular_coating,59); // In order to have diffuse reclectivity (Lambertian), it is necessary define all the other three.
   MTP_BaSO->AddProperty("SPECULARSPIKECONSTANT", photonEnergy9, other_coating,2); //  The diffuse is 1-other three (in this case 1-sppecular). 
   MTP_BaSO->AddProperty("BACKSCATTERCONSTANT", photonEnergy9, other_coating,2);
   
   MTP_BaSO->AddProperty("REFLECTIVITY", p_coating_refl, refl_coating,19);
   //MTP_BaSO->AddProperty("EFFICIENCY", <>, <>);

   BaSO4_surface->SetMaterialPropertiesTable(MTP_BaSO); 	
   

   G4LogicalSkinSurface* Surface1 = new G4LogicalSkinSurface("BaSO4_surface", ReflectBox_log, BaSO4_surface);

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

  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++) {
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
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintillatorBox_phys, Outer_tube_phys_vect[pos], Scintillator_PMMA_Surface);
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintillatorBox_phys, PMMA_Ring_phys_vect[pos], Scintillator_PMMA_Surface);
     new G4LogicalBorderSurface( (std::string("FiberBoxSurfaceOut_")+std::to_string(pos)).c_str(), ScintillatorBox_phys, Inner_tube_phys_vect[pos], Scintillator_PMMA_Surface);
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
  ReflectBox = new G4ExtrudedSolid("ReflectBox", reflect, ReflectZ/2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);

  // Scintillator box
  std::vector<G4TwoVector> scint(4);
  scint[0].set(-402.848*mm + WallThickness_XY + Thickness_Reflect, 642.936*mm - WallThickness_XY - Thickness_Reflect);
  scint[1].set(397.515*mm - WallThickness_XY - Thickness_Reflect, 589.330*mm - WallThickness_XY - Thickness_Reflect);
  scint[3].set(-402.484*mm + WallThickness_XY + Thickness_Reflect, -603.021*mm + WallThickness_XY + Thickness_Reflect);
  scint[2].set(397.515*mm - WallThickness_XY - Thickness_Reflect, -635.084*mm + WallThickness_XY + Thickness_Reflect);
  ScintillatorBox = new G4ExtrudedSolid("ScintillatorBox", scint, SctZ/2., G4TwoVector(0., 0.), 1., G4TwoVector(0., 0.), 1.);
  //ScintillatorBox = new G4Box("ScintillatorBox",SctX/2,SctY/2,SctZ/2);

  G4double Rin, Rout;
  G4double delta_Z;

  // Outer tube
  Rin = Diam_Out_In/2;
  Rout = Diam_Out_Out/2;
  Outer_tube  = new G4Tubs("Outer_tube", Rin, Rout, Length_Out/2, 0, 360*deg);
  // Air gap 1
  Rin = Diam_WOM_Out/2 + Thickness_WLS;
  Rout = Diam_Out_In/2;
  Air_gap1  = new G4Tubs("Air_gap1", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // WLS tube 1
  Rin = Diam_WOM_Out/2;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WLS_tube1  = new G4Tubs("WLS_tube1", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // WOM tube
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WOM_tube  = new G4Tubs("WOM_tube", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // WLS tube 2
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_In/2;
  WLS_tube2  = new G4Tubs("WLS_tube2", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // Air gap 2
  Rin = Diam_In_Out/2;
  Rout = Diam_WOM_In/2 - Thickness_WLS;
  Air_gap2  = new G4Tubs("Air_gap2", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // Inner tube
  Rin = Diam_In_In/2;
  Rout = Diam_In_Out/2;
  Inner_tube  = new G4Tubs("Inner_tube", Rin, Rout, Length_In/2, 0, 360*deg);
  // PMMA disk
  Rin = 0.0*mm;
  Rout = Diam_In_Out/2;
  PMMA_disk  = new G4Tubs("PMMA_disk", Rin, Rout, Thickness_Disk/2, 0, 360*deg);
  // PMMA ring
  Rin = Diam_In_In/2;
  Rout = Diam_Out_Out/2;
  PMMA_Ring  = new G4Tubs("PMMA_ring", Rin, Rout, Thickness_Ring/2, 0, 360*deg);
  // WLS ring
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WLS_ring = new G4Tubs("WLS_ring ", Rin, Rout, Thickness_WLS, 0, 360*deg);
  // Air ring outer
  Rin = (Diam_WOM_Out - 1*mm)/2;
  Rout = Diam_WOM_Out/2;
  Air_ring1 = new G4Tubs("Air_ring1", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // PMMA ring supporting WOM
  Rin = (Diam_WOM_In + 1*mm)/2;
  Rout = (Diam_WOM_Out - 1*mm)/2;
  PMMA_ring_lower = new G4Tubs("PMMA_ring_lower", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // Air ring inner
  Rin = Diam_WOM_In/2;
  Rout = (Diam_WOM_In + 1*mm)/2;
  Air_ring2 = new G4Tubs("Air_ring2", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // Hole in box
  Rin = 0.0*mm;
  Rout = Diam_Out_Out/2;
  Hole_box  = new G4Tubs("Hole_box", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // Hole in scintillator
  Rin = 0.0*mm;
  Rout = Diam_Out_Out/2;
  G4double Length_Scint_Hole = Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY;
  Hole_sct = new G4Tubs("Hole_sct", Rin, Rout, Length_Scint_Hole/2, 0, 360*deg);
  // PMMA "hat"
  Rin = Diam_Out_Out/2;
  Rout = Diam_Hat/2;
  PMMA_Hat = new G4Tubs("PMMA_hat", Rin, Rout, Thickness_Hat/2, 0, 360*deg);
  // Additional steel
  Rin = Diam_Hole/2;
  Rout = Diam_Steel_Add/2;
  SteelAdd = new G4Tubs("Steel_add", Rin, Rout, Thickness_Steel_Add/2, 0, 360*deg);
  // LAB&PPO inside tube
  Rin = 0.0*mm;
  Rout = Diam_In_In/2;
  SctInside = new G4Tubs("Sct_inside", Rin, Rout, Length_In/2, 0, 360*deg);

  G4SubtractionSolid *EmptySteelBox = new G4SubtractionSolid("EmptySteelBox",SteelBox,ReflectBox,0,G4ThreeVector(0,0,0));
  G4SubtractionSolid *EmptyReflectBox = new G4SubtractionSolid("EmptyReflectBox",ReflectBox,ScintillatorBox,0,G4ThreeVector(0,0,0));

  G4double delta_Z_EmptySteelBoxWithHole = SteelZ/2 - Length_WOM/2;
  G4double delta_Z_EmptyReflectBoxWithHole = ReflectZ/2 - Length_WOM/2;
  G4double delta_Z_ScintillatorBoxWithHole = SctZ/2 - Length_Scint_Hole/2;

  std::vector<G4SubtractionSolid*> EmptySteelBoxWithHole_tempvec;
  std::vector<G4SubtractionSolid*> EmptyReflectBoxWithHole_tempvec;
  std::vector<G4SubtractionSolid*> ScintillatorBoxWithHole_tempvec;
  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++) {
    if(pos == 0) {
      EmptySteelBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("EmptySteelBoxWithHole_")+std::to_string(pos)).c_str(),
                                              EmptySteelBox,Hole_box,0,
                                              G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_EmptySteelBoxWithHole)) );
      EmptyReflectBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("EmptyReflectBoxWithHole_")+std::to_string(pos)).c_str(),
                                              EmptyReflectBox,Hole_box,0,
                                              G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_EmptySteelBoxWithHole)) );
      ScintillatorBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("ScintillatorBoxWithHole_")+std::to_string(pos)).c_str(),
                                                ScintillatorBox, Hole_sct, 0,
                                                G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_ScintillatorBoxWithHole)) );
    }
    else {
      EmptySteelBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("EmptySteelBoxWithHole_")+std::to_string(pos)).c_str(),
                                              EmptySteelBoxWithHole_tempvec[pos-1] ,Hole_box,0,
                                              G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_EmptySteelBoxWithHole)) );
      EmptyReflectBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("EmptyReflectBoxWithHole_")+std::to_string(pos)).c_str(),
                                              EmptyReflectBoxWithHole_tempvec[pos-1] ,Hole_box,0,
                                              G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_EmptyReflectBoxWithHole)) );
      ScintillatorBoxWithHole_tempvec.push_back( new G4SubtractionSolid( (std::string("ScintillatorBoxWithHole_")+std::to_string(pos)).c_str(),
                                                ScintillatorBoxWithHole_tempvec[pos-1] , Hole_sct, 0,
                                                G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_ScintillatorBoxWithHole)) );
    }
  }
  EmptySteelBoxWithHole = EmptySteelBoxWithHole_tempvec[WOM_coord_vec.size()-1];
  EmptyReflectBoxWithHole = EmptyReflectBoxWithHole_tempvec[WOM_coord_vec.size()-1];
  ScintillatorBoxWithHole = ScintillatorBoxWithHole_tempvec[WOM_coord_vec.size()-1];

  sipmBox = new G4Box("sipmBox", sipmSize/2., sipmSize/2., sipmWindowThickness/2.);
  sipmBaseBox = new G4Box("sipmBase", sipmSize/2., sipmSize/2., sipmBaseThickness/2.);
  WOM_cellBox = new G4Box("wom_cell", SteelX/2,SteelY/2,SteelZ/2 + 15*cm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineLogicalVolumes()
{
  expHall_log = new G4LogicalVolume(expHall_box, air,"World",0,0,0);
  WOM_cell_log = new G4LogicalVolume(WOM_cellBox, air,"wom_cell",0,0,0);
  sipm_base_log = new G4LogicalVolume(sipm_base, steel,"sipm_base",0,0,0);
  ScintillatorBox_log = new G4LogicalVolume(ScintillatorBoxWithHole, LAB_PPO,"ScintillatorBoxLV",0,0,0);
  SteelBox_log = new G4LogicalVolume(EmptySteelBoxWithHole, steel,"Steelbox",0,0,0);
  ReflectBox_log = new G4LogicalVolume(EmptyReflectBoxWithHole, BaSO4,"Reflectbox", 0,0,0);
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
  G4RotationMatrix *RM1 = new G4RotationMatrix(0*deg,0*deg,0*deg);

  SteelBox_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),SteelBox_log,"SteelBox",WOM_cell_log,false,0, intersect_check);
  ReflectBox_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),ReflectBox_log,"ReflectBox",WOM_cell_log,false,0, intersect_check);
  ScintillatorBox_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),ScintillatorBox_log,"ScintillatorBoxPV",WOM_cell_log,false,0, intersect_check);

  sipmBase_phys = new G4PVPlacement(0,G4ThreeVector(0, 0, sipmWindowThickness/2. - sipmBaseThickness/2. ),sipmBaseBox_log,"sipmBase",sipmBox_log,false,0, intersect_check);

//-------------------------------------------------------------------
  G4double delta_Z_sipm = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat + Thickness_WLS + Thickness_Gap + sipmWindowThickness/2;
//-------------------------------------------------------------------

  // PMMA Staff
  // Outer_tube
  G4double delta_Z_Outer_tube = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY) + Length_Out/2 + Thickness_Ring;
  // WOM tube
  G4double delta_Z_WOM = SteelZ/2 - Length_WOM/2 + Thickness_Steel_Add + Thickness_Hat + Thickness_Gap;
  // Inner tube
  G4double delta_Z_Inner_tube = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY) + Length_In/2 + Thickness_Ring;
  // PMMA Ring
  G4double delta_Z_PMMA_Ring = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY) + Thickness_Ring/2;
  // PMMA Disk
  G4double delta_Z_PMMA_Disk = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY) + Length_In + Thickness_Ring + Thickness_Disk/2;
  // PMMA Hat
  G4double delta_Z_PMMA_Hat = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat/2;
  // Additional Steel
  G4double delta_Z_Steel_Add = SteelZ/2 + Thickness_Steel_Add/2;
  // LAB&PPO inside tube
  G4double delta_Z_Sct_Inside = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY) + Length_In/2;
  // Air ring
  G4double delta_Z_upper_ring = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThickness_XY) + Thickness_Ring + Thickness_Gap/2;

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
    Air_gap_1_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), Air_gap1_log, "Air_gap1", WOM_cell_log, false, 0, intersect_check) );
    Air_gap_2_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), Air_gap2_log, "Air_gap2", WOM_cell_log, false, 0, intersect_check) );
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
      sipm_phys_vect.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first +  radius_sipm*std::cos(i*2*pi/n_sipm),WOM_coord_vec[pos].second + radius_sipm*std::sin(i*2*pi/n_sipm),delta_Z_sipm),sipmBox_log, "sipm", WOM_cell_log, false, sipm_id++, intersect_check) );
    }
  }

  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 0, intersect_check);
  //RM1 = new G4RotationMatrix();
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(-SteelX/2.,-SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 2, intersect_check) );
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(+SteelX/2.,+SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 0, intersect_check) );    
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(-SteelX/2.,+SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 1, intersect_check) );
  //WOM_cells_phys_vect.push_back( new G4PVPlacement(RM1, G4ThreeVector(+SteelX/2.,-SteelY/2.,0.), WOM_cell_log, "wom_cell", expHall_log, false, 3, intersect_check) );     
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
  ScintillatorBox_log->SetVisAttributes(sctBoxVisAtt);
  Sct_Inside_log->SetVisAttributes(sctBoxVisAtt);

  G4VisAttributes *PMMAVisAtt = new G4VisAttributes;
  PMMAVisAtt->SetVisibility(true);
  PMMAVisAtt->SetColor(grey);
  PMMA_disk_log->SetVisAttributes(PMMAVisAtt);
  //PMMA_Ring_log->SetVisAttributes(PMMAVisAtt);
  PMMA_Hat_log->SetVisAttributes(PMMAVisAtt);
  Outer_tube_log->SetVisAttributes(PMMAVisAtt);
  //Inner_tube_log->SetVisAttributes(PMMAVisAtt);
  //PMMA_ring_lower_log->SetVisAttributes(PMMAVisAtt);

  G4VisAttributes *airVisAtt = new G4VisAttributes;
  airVisAtt->SetColor(green);
  airVisAtt->SetVisibility(true);
  Air_gap1_log->SetVisAttributes(airVisAtt);
  Air_gap2_log->SetVisAttributes(airVisAtt);
  Air_ring1_log->SetVisAttributes(airVisAtt);
  Air_ring2_log->SetVisAttributes(airVisAtt);

  G4VisAttributes *WLSVisAtt = new G4VisAttributes;
  WLSVisAtt->SetColor(red);
  WLSVisAtt->SetVisibility(false);
  WLS_tube1_log->SetVisAttributes(WLSVisAtt);
  WLS_tube2_log->SetVisAttributes(WLSVisAtt);
  ReflectBox_log->SetVisAttributes(WLSVisAtt);

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
