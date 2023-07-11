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
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4SDParticleFilter.hh"
#include "G4UnionSolid.hh"
#include "G4PhysicalConstants.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction::OpNoviceDetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fExpHall_x = 5*m;
  fExpHall_y = 5*m;
  fExpHall_z = 5*m;

  SteelZ = 22*cm;
  WallThick = 10*mm;
  WallThickness_Z_Cover = 5*mm;
  ReflectZ = SteelZ - WallThick*2;
  ReflectThick = 0.05*mm;
  SctZ = SteelZ - WallThick*2 - ReflectThick*2;

  Diam_In_In = 44*mm;
  Diam_In_Out = 50*mm;
  Diam_Out_In = 64*mm;
  Diam_Out_Out = 70*mm;
  Diam_WOM_In = 54*mm;
  Diam_WOM_Out = 60*mm;
  Diam_Hole = 72*mm;
  Diam_Steel_Add = 120*mm;
  Diam_Hat = 120*mm;
  Length_Out = 185*mm;
  Length_In = 175*mm;
  Thickness_Ring = 4*mm;
  Thickness_Disk = 5*mm;
  Thickness_Hat = 10*mm;
  Thickness_Steel_Add = 15*mm;
  Thickness_Gap = 1*mm;
  Thickness_SupportRing = 1*mm;
  Thickness_WLS = 0.02*mm;
  Length_WOM = 180*mm;
  
  sipmSize = 3.4*mm;
  sipmSizeSens = 3*mm;
  sipmBaseThickness = 1.*mm;
  sipmWindowThickness =  0.15*mm;
  sipmSensThickness = 0.149*mm;
  sipmSensThicknessTop = 0.001*mm;

  WOM_coord_vec = {{-397*mm,321.88*mm}, {-397*mm,953.7975*mm}, {397*mm,285.1*mm}, {397*mm,905.2025*mm}, {-397*mm,-762.54125*mm}, {-397*mm,-240.12*mm}, {397*mm,-778.45875*mm}, {397*mm,-266.88*mm}};
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
  Silicon = nist->FindOrBuildMaterial("G4_Si");
  Al = nist->FindOrBuildMaterial("G4_Al");
  
  G4double a, z, density;
  G4int nelements, ncomponent, natoms;

  // Air
  G4Element* N = new G4Element("Nitrogen","N", z = 7, a = 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen","O", z = 8, a = 16.00*g/mole);

  air = new G4Material("Air", density = 1.29*mg/cm3, nelements = 2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);

  // Linear alkyl benzene (LAB)
  G4Element* H = new G4Element("Hydrogen", "H", 1 , 1.01*g/mole);
  G4Element* C = new G4Element("Carbon", "C", 6 , 12.01*g/mole);
  G4Material* LAB = new G4Material("LAB", density = 0.856*g/cm3, ncomponent = 2);  //density: https://www.knowde.com/stores/sasol/documents/101895
  LAB->AddElement(H, natoms = 28);
  LAB->AddElement(C, natoms = 17);
  // Diphenyloxazole (PPO)
  G4Material* PPO = new G4Material("PPO", density = 1.128*g/cm3, ncomponent = 4); // density: https://www.echemi.com/sds/24-diphenyloxazole-pid_Rock24446.html
  PPO->AddElement(H, natoms = 11);
  PPO->AddElement(C, natoms = 15);
  PPO->AddElement(N, natoms = 1);
  PPO->AddElement(O, natoms = 1);
  // Scintillator (LAB+PPO) 23233 cm^3 23.233 l
  LAB_PPO = new G4Material("LAB_PPO", density = 0.9*g/cm3, ncomponent = 2); //??
  LAB_PPO->AddMaterial(LAB, 99.77*perCent); //Calculated considering  be 2 g/L of PPO (First measurement of the surface tension of a liquid scintillator based on Linear Alkylbenzene (HYBLENE 113)). 2g of PPO every liter of LAB
  LAB_PPO->AddMaterial(PPO, 0.23*perCent);
  // Bis-MSB WLS
  Bis_MSB = new G4Material("Bis_MSB", density = 1.076*g/cm3, ncomponent = 2); // density: http://www.molbase.com/moldata/368101.html
  Bis_MSB->AddElement(H, natoms = 22);
  Bis_MSB->AddElement(C, natoms = 24);
  // PEMA to build the WLS dye coat  
  PEMA = new G4Material("PEMA", density = 1.11*g/cm3, ncomponent = 3); // density: https://polymerdatabase.com/polymers/polyethylmethacrylate.html
  PEMA->AddElement(H, natoms = 10);
  PEMA->AddElement(C, natoms = 6);
  PEMA->AddElement(O, natoms = 242);
  // PTP (para-Terphenyl) to build the scintillator C18H14 and the dye coat
  PTP = new G4Material("PTP", density=1.23*g/cm3, ncomponent = 2); //  density: https://m.molbase.com/moldata/64879.html   
  PTP->AddElement(H, natoms = 14);
  PTP->AddElement(C, natoms = 18);
  /// WLS Coating  (150g PEMA, 3g PTP 1.5g bis-MSB)
  WLS_Coat = new G4Material("WLSCoat", density = 1.1*g/cm3, ncomponent = 3); // density: same as PEMA 
  WLS_Coat->AddMaterial(Bis_MSB, 0.97*perCent);
  WLS_Coat->AddMaterial(PTP, 1.94*perCent);
  WLS_Coat->AddMaterial(PEMA, 97.09*perCent);
  // PMMA side
  PMMA_side = new G4Material("PMMA_side", density = 1.200*g/cm3, ncomponent = 2);
  PMMA_side->AddElement(H, natoms = 2);
  PMMA_side->AddElement(C, natoms = 4);
  // PMMA bottom
  PMMA_bottom = new G4Material("PMMA_bottom", density = 1.200*g/cm3, ncomponent = 2);
  PMMA_bottom->AddElement(H, natoms = 2);
  PMMA_bottom->AddElement(C, natoms = 4);
  // Barium sulphate (BaSO4) Reflectivity coating
  G4Element* Ba = new G4Element("Barium", "Ba", 56 , 137.327*g/mole);
  G4Element* S = new G4Element("Sulphur", "S", 16 , 32.065*g/mole);
  BaSO4 = new G4Material("BaSO4", density = 4.5*g/cm3, ncomponent = 3); // https://www.chemeurope.com/en/encyclopedia/Barium_sulfate.html
  BaSO4->AddElement(Ba, natoms = 1);
  BaSO4->AddElement(S, natoms = 1);
  BaSO4->AddElement(O, natoms = 4);  
  
  // Silicone resin
  G4Element* Si = new G4Element("Silicon", "Si", 14 , 28.0855*g/mole);
  ResinSi = new G4Material("ResinSi",density=3*g/cm3,ncomponent=2); 
  ResinSi->AddElement(Si,natoms=1);
  ResinSi->AddElement(O,natoms=4);
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
    G4double rind = 1.;
    G4double B[4], C[4];
    B[0] = 0.821384;  C[0] = 94.7625;
    B[1] = 0.311375;  C[1] = 160.751;
    B[2] = 0.0170099; C[2] = 219.575;
    B[3] = 0.608268;  C[3] = 9385.54;
    for(int term = 0; term < 4; term++) rind += B[term]/(1. - (C[term]/wl)*(C[term]/wl)); //formula eand coefficients: https://arxiv.org/pdf/1105.2101.pdf
    return sqrt(rind);
  };

  G4double wl;
  for(unsigned int i = 0; i < sizeof(photon_en_LAB_PPO)/sizeof(photon_en_LAB_PPO[0]); i++) {
    wl = 745. - 5.*i;
    photon_en_LAB_PPO[i] = 1240./wl*eV;
    rindex_LAB_PPO[i] = Ridndex_LAB_PPO(wl);
  }
  G4MaterialPropertiesTable *MPT_LAB_PPO = new G4MaterialPropertiesTable();
  MPT_LAB_PPO->AddConstProperty("SCINTILLATIONYIELD", 10800./MeV); // https://underground.physics.berkeley.edu/WbLS/slides/PennRnD-Grullon.pdf
  MPT_LAB_PPO->AddProperty("RINDEX", photon_en_LAB_PPO, rindex_LAB_PPO, 100);//->SetSpline(true);

  // emission
  G4double photonWaveLength3[201] = {500,499,498,497,496,495,494,493,492,491,490,489,488,487,486,485,484,483,482,481,480,479,478,477,476,475,474,473,472,471,470,469,468,467,466,465,464,463,462,461,460,459,458,457,456,455,454,453,452,451,450,449,448,447,446,445,444,443,442,441,440,439,438,437,436,435,434,433,432,431,430,429,428,427,426,425,424,423,422,421,420,419,418,417,416,415,414,413,412,411,410,409,408,407,406,405,404,403,402,401,400,399,398,397,396,395,394,393,392,391,390,389,388,387,386,385,384,383,382,381,380,379,378,377,376,375,374,373,372,371,370,369,368,367,366,365,364,363,362,361,360,359,358,357,356,355,354,353,352,351,350,349,348,347,346,345,344,343,342,341,340,339,338,337,336,335,334,333,332,331,330,329,328,327,326,325,324,323,322,321,320,319,318,317,316,315,314,313,312,311,310,309,308,307,306,305,304,303,302,301,300};
  G4double photon_en_LAB_PPO_2[201];
  for(unsigned int i = 0; i < sizeof(photonWaveLength3)/sizeof(photonWaveLength3[0]); i++) photon_en_LAB_PPO_2[i] = 1240./photonWaveLength3[i]*eV;

  G4double scintilFast_LAB_PPO[201] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06,0.0625,0.065,0.0675,0.07,0.0725,0.075,0.08,0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12,0.125,0.1325,0.14,0.1475,0.155,0.1625,0.17,0.1775,0.185,0.1925,0.2,0.2075,0.215,0.2225,0.23,0.2375,0.245,0.2525,0.26,0.2675,0.275,0.29,0.305,0.32,0.335,0.35,0.365,0.38,0.395,0.41,0.425,0.4325,0.44,0.4475,0.455,0.4625,0.47,0.4775,0.485,0.4925,0.5,0.53,0.56,0.59,0.62,0.65,0.68,0.71,0.74,0.77,0.8,0.795,0.79,0.785,0.78,0.775,0.77,0.765,0.76,0.755,0.75,0.77,0.79,0.81,0.83,0.85,0.87,0.89,0.91,0.93,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45,0.405,0.36,0.315,0.27,0.225,0.18,0.135,0.09,0.045,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}; //https://arxiv.org/abs/1001.3946 
  MPT_LAB_PPO->AddProperty("FASTCOMPONENT", photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201)->SetSpline(true);
  MPT_LAB_PPO->AddProperty("SLOWCOMPONENT", photon_en_LAB_PPO_2, scintilFast_LAB_PPO, 201)->SetSpline(true);

  // transmission
  G4double waveLength[151] = {600,598,596,594,592,590,588,586,584,582,580,578,576,574,572,570,568,566,564,562,560,558,556,554,552,550,548,546,544,542,540,538,536,534,532,530,528,526,524,522,520,518,516,514,512,510,508,506,504,502,500,498,496,494,492,490,488,486,484,482,480,478,476,474,472,470,468,466,464,462,460,458,456,454,452,450,448,446,444,442,440,438,436,434,432,430,428,426,424,422,420,418,416,414,412,410,408,406,404,402,400,398,396,394,392,390,388,386,384,382,380,378,376,374,372,370,368,366,364,362,360,358,356,354,352,350,348,346,344,342,340,338,336,334,332,330,328,326,324,322,320,318,316,314,312,310,308,306,304,302,300};
  G4double photon_en_LAB_PPO_3[151];
  for(unsigned int i = 0; i < sizeof(waveLength)/sizeof(waveLength[0]); i++) photon_en_LAB_PPO_3[i] = 1240./waveLength[i]*eV;

  G4double absLen_unpurified[151] = {7.41928480365393*m,7.68678403498895*m,7.91287034263266*m,7.98944770316961*m,8.08548242289423*m,8.05778791764816*m,8.09059945466605*m,8.11488877012421*m,8.15949366895364*m,8.24046357622784*m,8.36896424710391*m,8.5724798759478*m,8.85772575154252*m,8.85221167762828*m,8.81555436133421*m,9.05389571096987*m,9.48924122323874*m,9.76515374992161*m,9.51335775244241*m,9.45328043541621*m,8.97482880285293*m,8.70096539741452*m,8.48935917676502*m,8.44097817045948*m,8.44405991017534*m,8.43190261832563*m,8.38645678568186*m,8.35648742634518*m,8.33232353493578*m,8.33114526507463*m,8.23755943796866*m,8.13213459072979*m,7.8726545001728*m,7.65415458551032*m,7.5800588175385*m,7.6061295282424*m,7.69320756579235*m,7.7082984257603*m,7.70586123507162*m,7.70325134674453*m,7.66923841145675*m,7.64365261429804*m,7.62664064911952*m,7.55677185990663*m,7.50711974244003*m,7.45713654579468*m,7.43721662692559*m,7.35800834707499*m,7.22642533541789*m,7.11844709300408*m,7.01469716585332*m,6.96176069278602*m,6.90597103165448*m,6.84830374947254*m,6.84525178280356*m,6.83526672561362*m,6.79526318509296*m,6.71587945264393*m,6.63947275377847*m,6.6118444804822*m,6.49136793799083*m,6.42934488959944*m,6.31098091183809*m,6.21946196680514*m,6.15871441142064*m,6.09293247274225*m,6.02864030249323*m,5.94477355531049*m,5.86604224769804*m,5.80521857408766*m,5.73506328468108*m,5.66214631723293*m,5.57289275228197*m,5.49731700891059*m,5.39615082475944*m,5.30550194154553*m,5.19943157180204*m,5.0911745929865*m,4.99959477462503*m,4.90159698959031*m,4.7917682338105*m,4.6681139833648*m,4.55010411037092*m,4.40871123138021*m,4.28519387842501*m,4.16757742933685*m,4.04589813769501*m,3.91998342801409*m,3.79194507519573*m,3.65200866136557*m,3.51899905487626*m,3.36858314280228*m,3.219533855462*m,3.06634990702638*m,2.90804104570209*m,2.74757561840523*m,2.58972782886322*m,2.41918635504621*m,2.23054267956921*m,2.05636621586457*m,1.94127546199575*m,1.85826095244245*m,1.74682662700547*m,1.59879984499012*m,1.40361593174363*m,1.21010948607421*m,1.12298134358669*m,1.12861772574266*m,1.1046585196898*m,1.04234271232283*m,0.9941373021031*m,0.956557242515996*m,0.927425529399054*m,0.875744366988952*m,0.779908411841287*m,0.663313263320605*m,0.541094630094606*m,0.406926671366136*m,0.265001590325018*m,0.153480082825598*m,0.085384830013475*m,0.048210452365814*m,0.028221561196635*m,0.017250128891369*m,0.013236885579806*m,0.013016803788927*m,0.012986890579723*m,0.012971261573831*m,0.012958145564748*m,0.01294674529108*m,0.01293558959276*m,0.012928356570193*m,0.01292213175909*m,0.012913052869886*m,0.012909213986899*m,0.01290710804677*m,0.012901565135242*m,0.012897247581221*m,0.012890383936614*m,0.012881105659852*m,0.01286687734122*m,0.012849574178154*m,0.012838615276908*m,0.012833994111742*m,0.012830840543912*m,0.012822230954434*m,0.01281835788177*m,0.012816426865377*m,0.012807751818825*m,0.0128029524293*m,0.01279493404851*m}; // Patrick measurement Mainz TB_scintillator

  G4double absLen_purified[151] = {8.69251407037729*m,9.02995755323176*m,9.27373907607242*m,9.2987694632999*m,9.37254639295624*m,9.31493016410715*m,9.32945396614867*m,9.36095897891052*m,9.48082544118919*m,9.60814130202656*m,9.82801139775702*m,10.1046947101091*m,10.4907062873286*m,10.5108807867562*m,10.4753853636081*m,10.8002040964469*m,11.4980811433371*m,11.9325107778151*m,11.6636028590624*m,11.7094782859565*m,10.9852633716743*m,10.670267622928*m,10.3800086540752*m,10.3082593229783*m,10.3275701553562*m,10.3319544221576*m,10.2662283365215*m,10.2925732890662*m,10.3153782242223*m,10.3860868021326*m,10.3365901053944*m,10.231517114707*m,9.92788014563769*m,9.61898800830966*m,9.54402290594051*m,9.62746317964499*m,9.79692373936646*m,9.91423094765397*m,9.9643151940587*m,10.0187720504539*m,10.0358961669191*m,10.0802588950859*m,10.1054177127325*m,10.0818415254723*m,10.086029649153*m,10.1133311242014*m,10.1350472968498*m,10.0761167690658*m,9.90904041225772*m,9.77530373579883*m,9.66256670125958*m,9.65338010651176*m,9.63069786946983*m,9.62140630852179*m,9.71997172370098*m,9.79744658014088*m,9.80672901676583*m,9.78343285609755*m,9.72163465573446*m,9.80518793240826*m,9.69593799132894*m,9.67791627886446*m,9.59985620139408*m,9.55705239801755*m,9.57416795211246*m,9.56671620921144*m,9.59797651082656*m,9.56114765877604*m,9.53375540026655*m,9.55647549660189*m,9.54447464256824*m,9.53291003815433*m,9.49239969544417*m,9.49917662439382*m,9.42391963661282*m,9.3663798468712*m,9.32376570267581*m,9.22629281746381*m,9.17977310076775*m,9.12179277182491*m,9.00886669573072*m,8.85479490604375*m,8.71252334840956*m,8.50086621234753*m,8.33199872471845*m,8.18529383430015*m,8.02793280531241*m,7.84805101504057*m,7.67375314752583*m,7.46193191606457*m,7.24933185585757*m,7.00057529008316*m,6.78144433425612*m,6.57082163294415*m,6.3463075559428*m,6.13214412943496*m,5.9379419460511*m,5.73334467669744*m,5.53766879317458*m,5.3394385508064*m,5.15238124063095*m,4.98340338198499*m,4.8012155457967*m,4.60195158864653*m,4.41961657539783*m,4.22693302308812*m,4.0416225005906*m,3.87499604590985*m,3.71820540936771*m,3.5472216422986*m,3.36485969966114*m,3.15134802038665*m,2.90168745463681*m,2.5867387338465*m,2.16375663283064*m,1.66966271660288*m,1.14377099881739*m,0.680589789294733*m,0.360363924864896*m,0.181782601946722*m,0.093267712946957*m,0.0503991725824*m,0.02880952553709*m,0.01737142379557*m,0.013310009583829*m,0.013096262643175*m,0.013063354253155*m,0.013043343934459*m,0.013030737069578*m,0.013019803591978*m,0.013008275438005*m,0.012998844829755*m,0.012996396487048*m,0.012992079547825*m,0.012984198598302*m,0.012978816961298*m,0.012975348745523*m,0.012968698731763*m,0.012962219216222*m,0.012956981482458*m,0.012941499851657*m,0.012919404279236*m,0.012913555972765*m,0.012909301261196*m,0.012901893060039*m,0.012901203446858*m,0.012897730665173*m,0.012889608969039*m,0.012882419482255*m,0.012878996045324*m,0.0128715069987*m}; // Patrick measurement Mainz Column00+ppo
  MPT_LAB_PPO->AddProperty("ABSLENGTH", photon_en_LAB_PPO_3, absLen_purified, 151)->SetSpline(true);
  MPT_LAB_PPO->AddConstProperty("RESOLUTIONSCALE", 1.0); //??
  MPT_LAB_PPO->AddConstProperty("FASTTIMECONSTANT", 2*ns); // Main results 
  MPT_LAB_PPO->AddConstProperty("SLOWTIMECONSTANT", 14*ns); // Main results 
  MPT_LAB_PPO->AddConstProperty("YIELDRATIO", 0.5);
  LAB_PPO->SetMaterialPropertiesTable(MPT_LAB_PPO);

  //------------------------------------------------------------------------------
  //----------------------------- PMMA -----------------------------
  //------------------------------------------------------------------------------
  const G4int pmma_mpt_entr = 13;

  G4double pmma_refl_index_wl[pmma_mpt_entr] = {700.,600.,550.,500.,450.,400.,390.,380.,370.,350.,320.,310.,300.};
  G4double pmma_side_wl[8] = {632,543,440,405,375,320,310,300};
  G4double pmma_bottom_wl[37] = {632,543,440,405,375,358.8171433,353.7663967,348.7156502,343.6649037,338.6141572,333.5634107,328.5126642,323.4619177,318.4111712,313.3604247,308.3096782,303.2589317,298.3612381,293.1574387,288.1066922,283.0559457,278.0051992,272.9544526,267.9037061,262.8529596,257.8022131,252.7514666,247.8920363,242.6499736,237.5992271,232.5484806,227.4977341,222.4469876,217.3962411,212.3454946,207.2947481,202.0526854};
  G4double pmma_rind[pmma_mpt_entr] = {1.489,1.492,1.495,1.498,1.502,1.511,1.512,1.514,1.516,1.522,1.54,1.541,1.542}; // https://refractiveindex.info/?shelf=3d&book=plastics&page=pmma
  G4double pmma_refl_index_en[pmma_mpt_entr];
  G4double pmma_side_en[8];
  G4double pmma_bottom_en[37];
  for(int i = 0; i < pmma_mpt_entr; i++) pmma_refl_index_en[i] = 1240./pmma_refl_index_wl[i]*eV;
  for(unsigned int i = 0; i < sizeof(pmma_side_en)/sizeof(pmma_side_en[0]); i++) pmma_side_en[i] = 1240./pmma_side_wl[i]*eV;
  for(unsigned int i = 0; i < sizeof(pmma_bottom_en)/sizeof(pmma_bottom_en[0]); i++) pmma_bottom_en[i] = 1240./pmma_bottom_wl[i]*eV;
  G4double pmma_side_abslen[8] = {1.5*m,1*m,1*m,1*m,604.44*mm,24.6*mm,18.23*mm,10.55*mm}; // Tecacryl 3-4 / RTP 1-2-3 / Spartech 1-2 /Mc Master 	(From Micheal: https://arxiv.org/pdf/1310.6454.pdf)	for >375nm --- old measurements from the tube for <375nm 

  G4double pmma_bottom_abslen[37] = {17900*mm,125150*mm,2440*mm,262.5*mm,95*mm,0.353271102*mm,0.307327896*mm,0.301885495*mm,0.301900827*mm,0.301917182*mm,0.301934641*mm,0.301953291*mm,0.301973222*mm,0.301994529*mm,0.302017312*mm,0.302041678*mm,0.302067738*mm,0.341435928*mm,0.302125414*mm,0.347533944*mm,0.437466744*mm,0.556750518*mm,0.647685273*mm,0.644020575*mm,0.53581648*mm,0.367918167*mm,0.302452281*mm,0.32079176*mm,0.302566913*mm,0.302631272*mm,0.302701282*mm,0.302777948*mm,0.302862616*mm,0.30295711*mm,0.303063921*mm,0.303186472*mm,0.321725432*mm}; //RPT 4 / Spartech 3-4	 (From Micheal: https://arxiv.org/pdf/1310.6454.pdf) for >375nm --- andrew and doramas measurement for <375nm
  
  G4double pmma_side_wl_new[57] = {550,545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465,460,455,450,445,440,435,430,425,420,415,410,405,400,395,390,385,380,375,370,365,360,355,350,345,340,335,330,325,320,315,310,305,300,295,290,285,280,275,270};
  
  G4double pmma_side_abslen_new[57] = {472.626403500485,472.47431927557,472.442875276582,472.548069065499,472.503941286975,472.717214657718,473.449818351276,473.924385233674,475.069753001741,474.877497105366,475.03999205846,473.890412182362,475.831343905636,473.79358106133,475.253004464295,470.790625727439,470.497047736906,472.259113792486,472.931200356411,472.365745218049,470.881247507505,470.708148061025,469.722655892388,469.20573460488,468.340978652268,467.724990027454,466.424650040469,465.387657347346,464.676776284062,463.446924175645,462.489280782535,461.435833610124,459.974991506271,458.087943295396,455.842768535932,452.224327990847,449.384367527057,444.861880129657,439.189956124801,432.269666813194,423.7477750101,414.172287702685,404.353357130681,390.519697889426,369.419537577731,346.777708977586,327.537561287855,296.568368423133,257.906895924098,214.610379424564,175.424523581044,141.385052228302,112.73196443914,90.6391278206925,73.2747503375892,57.4643292737885,44.7190787748292}; // from John Rack-Helleis (Mainz) measurements and corrected with Christian Scharf (Berlin) thesis formula 
  
  G4double pmma_bottom_wl_new[58] = {545,540,535,530,525,520,515,510,505,500,495,490,485,480,475,470,465,460,455,450,445,440,435,430,425,420,415,410,405,400,395,390,385,380,375,370,365,360,355,350,345,340,335,330,325,320,315,310,305,300,295,290,285,280,275,270,265,260};
  
  G4double pmma_bottom_abslen_new[58] = {351.350063334171,351.999125361121,350.789205816769,350.157544826916,349.964814524976,349.905090650762,349.209562761158,349.2001239779,349.050276048572,348.952244986273,348.460024783965,347.714842442585,347.227536096873,346.878098226731,346.294344025831,346.316453434341,345.45962390967,344.223231991828,343.439096892365,342.746147890309,341.376979566388,341.060656138989,340.221475069691,339.08960889221,337.887636842857,336.098835449247,333.584401979514,328.567896278262,316.064904235738,286.663031783464,229.193328437948,150.621295694047,82.6320843414912,43.3974238845659,23.9716600974233,15.1960657437843,12.4427779773735,12.1974417410274,12.4676517189804,12.9540782532771,13.3653510241641,13.541055915037,13.4195096806019,13.2944731221916,13.1983301819591,13.1537884218875,13.1125467948418,13.1384675907102,13.2390134915731,13.4511160842576,13.8139219043694,14.3241576852257,14.9667556104642,15.7120711674145,16.6545751702108,17.6424558957449,18.2588587216547,18.5796529296663}; // from John Rack-Helleis (Mainz) measurements and corrected with Christian Scharf (Berlin) thesis formula 
  
  G4double pmma_side_en_new[8];
  G4double pmma_bottom_en_new[37];
  for(unsigned int i = 0; i < sizeof(pmma_side_en_new)/sizeof(pmma_side_en_new[0]); i++) pmma_side_en_new[i] = 1240./pmma_side_wl_new[i]*eV;
  for(unsigned int i = 0; i < sizeof(pmma_bottom_en_new)/sizeof(pmma_bottom_en_new[0]); i++) pmma_bottom_en_new[i] = 1240./pmma_bottom_wl_new[i]*eV;

  G4MaterialPropertiesTable *MPT_PMMA_side = new G4MaterialPropertiesTable();
  MPT_PMMA_side->AddProperty("ABSLENGTH", pmma_side_en, pmma_side_abslen, 8)->SetSpline(true);
  MPT_PMMA_side->AddProperty("RINDEX", pmma_refl_index_en, pmma_rind, pmma_mpt_entr)->SetSpline(true);
  PMMA_side->SetMaterialPropertiesTable(MPT_PMMA_side);

  G4MaterialPropertiesTable *MPT_PMMA_bottom = new G4MaterialPropertiesTable();
  MPT_PMMA_bottom->AddProperty("RINDEX", pmma_refl_index_en, pmma_rind, pmma_mpt_entr)->SetSpline(true);
  MPT_PMMA_bottom->AddProperty("ABSLENGTH", pmma_bottom_en, pmma_bottom_abslen, 37);
  PMMA_bottom->SetMaterialPropertiesTable(MPT_PMMA_bottom);

  //------------------------------------------------------------------------------
  //----------------------------- WLS_Coat -----------------------------
  //------------------------------------------------------------------------------
  G4MaterialPropertiesTable *MPT_WLSCoat = new G4MaterialPropertiesTable();
  G4double waveLength3[272] = {560,559,558,557,556,555,554,553,552,551,550,549,548,547,546,545,544,543,542,541,540,539,538,537,536,535,534,533,532,531,530,529,528,527,526,525,524,523,522,521,520,519,518,517,516,515,514,513,512,511,510,509,508,507,506,505,504,503,502,501,500,499,498,497,496,495,494,493,492,491,490,489,488,487,486,485,484,483,482,481,480,479,478,477,476,475,474,473,472,471,470,469,468,467,466,465,464,463,462,461,460,459,458,457,456,455,454,453,452,451,450,449,448,447,446,445,444,443,442,441,440,439,438,437,436,435,434,433,432,431,430,429,428,427,426,425,424,423,422,421,420,419,418,417,416,415,414,413,412,411,410,409,408,407,406,405,404,403,402,401,400,399,398,397,396,395,394,393,392,391,390,389,388,387,386,385,384,383,382,381,380,379,378,377,376,375,374,373,372,371,370,369,368,367,366,365,364,363,362,361,360,359,358,357,356,355,354,353,352,351,350,349,348,347,346,345,344,343,342,341,340,339,338,337,336,335,334,333,332,331,330,329,328,327,326,325,324,323,322,321,320,319,318,317,316,315,314,313,312,311,310,309,308,307,306,305,304,303,302,301,300,299,298,297,296,295,294,293,292,291,290,289};
  G4double photonEnergy5[272];
  for(unsigned int i = 0; i < sizeof(photonEnergy5)/sizeof(photonEnergy5[0]); i++) photonEnergy5[i] = 1240./waveLength3[i]*eV;
  G4double absLen3[272] = {1875.34495*mm,1815.34495*mm,1755.34495*mm,1695.34495*mm,1635.34495*mm,1575.34495*mm,1515.34495*mm,1455.34495*mm,1415.34495*mm,1375.34495*mm,1335.34495*mm,1295.34495*mm,1255.34495*mm,1215.34495*mm,1175.34495*mm,1135.34495*mm,1095.34495*mm,1075.34495*mm,1055.34495*mm,1035.34495*mm,1015.34495*mm,995.34495*mm,975.34495*mm,955.34495*mm,935.34495*mm,915.34495*mm,895.34495*mm,875.34495*mm,855.34495*mm,835.34495*mm,815.34495*mm,795.34495*mm,775.34495*mm,755.34495*mm,735.34495*mm,715.34495*mm,695.34495*mm,675.34495*mm,655.34495*mm,635.34495*mm,625.34495*mm,615.34495*mm,605.34495*mm,595.34495*mm,585.34495*mm,575.34495*mm,565.34495*mm,555.34495*mm,545.34495*mm,535.34495*mm,525.34495*mm,515.34495*mm,505.34495*mm,495.34495*mm,485.34495*mm,475.34495*mm,465.34495*mm,455.34495*mm,445.34495*mm,439.34495*mm,433.34495*mm,427.34495*mm,421.34495*mm,415.34495*mm,409.34495*mm,403.34495*mm,397.34495*mm,391.34495*mm,385.34495*mm,379.34495*mm,373.34495*mm,367.34495*mm,361.34495*mm,355.34495*mm,349.34495*mm,343.34495*mm,337.34495*mm,331.34495*mm,325.34495*mm,319.34495*mm,313.34495*mm,307.34495*mm,301.34495*mm,295.34495*mm,289.34495*mm,283.34495*mm,277.34495*mm,271.34495*mm,265.34495*mm,259.34495*mm,253.34495*mm,247.34495*mm,242.84495*mm,238.34495*mm,233.84495*mm,229.34495*mm,224.84495*mm,220.34495*mm,215.84495*mm,211.34495*mm,206.84495*mm,202.34495*mm,197.84495*mm,193.34495*mm,188.84495*mm,184.34495*mm,179.84495*mm,175.34495*mm,170.84495*mm,166.34495*mm,161.84495*mm,157.34495*mm,152.84495*mm,148.34495*mm,143.84495*mm,139.34495*mm,134.84495*mm,130.34495*mm,125.84495*mm,121.34495*mm,116.84495*mm,112.34495*mm,107.84495*mm,103.34495*mm,98.84495*mm,95.34495*mm,91.84495*mm,88.34495*mm,84.84495*mm,81.34495*mm,77.84495*mm,74.34495*mm,70.84495*mm,67.34495*mm,63.84495*mm,60.34495*mm,56.84495*mm,53.34495*mm,49.84495*mm,46.34495*mm,42.84495*mm,39.34495*mm,35.84495*mm,32.34495*mm,28.84495*mm,26.34495*mm,23.84495*mm,21.34495*mm,18.84495*mm,16.34495*mm,13.84495*mm,11.34495*mm,8.84495*mm,6.34495*mm,3.84495*mm,1.34495*mm,0.4858*mm,0.30961*mm,0.21896*mm,0.15935*mm,0.11815*mm,0.08933*mm,0.06919*mm,0.0549*mm,0.04435*mm,0.03637*mm,0.0303*mm,0.02561*mm,0.02199*mm,0.01915*mm,0.0169*mm,0.01509*mm,0.01362*mm,0.01244*mm,0.01145*mm,0.01062*mm,0.00991*mm,0.00931*mm,0.0088*mm,0.00835*mm,0.00795*mm,0.0076*mm,0.00724*mm,0.00695*mm,0.00668*mm,0.00642*mm,0.00618*mm,0.00596*mm,0.00574*mm,0.00554*mm,0.00536*mm,0.00519*mm,0.00503*mm,0.0049*mm,0.00477*mm,0.00466*mm,0.00457*mm,0.0045*mm,0.00443*mm,0.00437*mm,0.00433*mm,0.0043*mm,0.00427*mm,0.00425*mm,0.00423*mm,0.00422*mm,0.00422*mm,0.00421*mm,0.00422*mm,0.00422*mm,0.00423*mm,0.00425*mm,0.00427*mm,0.00429*mm,0.00433*mm,0.00437*mm,0.00442*mm,0.00448*mm,0.00454*mm,0.00462*mm,0.0047*mm,0.00479*mm,0.00489*mm,0.005*mm,0.00512*mm,0.00524*mm,0.00538*mm,0.00552*mm,0.00568*mm,0.00584*mm,0.00601*mm,0.00619*mm,0.00637*mm,0.00654*mm,0.0067*mm,0.00685*mm,0.00697*mm,0.00707*mm,0.00713*mm,0.00717*mm,0.00717*mm,0.00714*mm,0.00708*mm,0.00698*mm,0.00684*mm,0.00667*mm,0.00647*mm,0.00626*mm,0.00602*mm,0.00579*mm,0.00556*mm,0.00535*mm,0.00515*mm,0.00498*mm,0.00484*mm,0.00474*mm,0.00469*mm,0.00469*mm,0.00473*mm,0.00482*mm,0.00497*mm,0.00517*mm,0.00544*mm,0.00584*mm,0.00628*mm,0.00694*mm,0.00779*mm,0.009*mm,0.01068*mm,0.01369*mm,0.01871*mm,0.03013*mm}; //Jakobs measure without BPEA  
  MPT_WLSCoat->AddProperty("WLSABSLENGTH", photonEnergy5, absLen3, 272);  // jakob thesis -- without BPEA 

  // BIS reemission
  G4double waveLength4[51] = {498.50394264782,495.512169509067,489.3462170527,484.343493743363,479.225925762846,472.348366539066,467.849074868481,464.04990973537,461.338683704906,459.793161903923,457.56412302422,454.795680527534,
  451.768205666593,446.867663596341,441.763591710159,435.933166703427,433.243802562634,431.985904788012,429.963943988042,427.926782673168,426.395623114075,425.126769876294,422.410632741153,420.384395226982,
  417.93569195544,415.411789507185,413.448584543676,412.387719607238,411.249479281543,409.9604817049,406.3089221806,404.020506020059,401.848223327902,399.264457914078,397.452088829242,394.82668707538,
  393.620389802785,392.661003458199,390.948731712797,389.702126722733,388.693734027439,388.169942931827,387.606839417529,386.814987134474,386.022324288337,384.176496975521,382.716212740349,381.014431692937,
  378.83045463701,376.455517150279,373.305908776801};   // ol. the actualy one  
  G4double photonEnergy6[51];
  for(unsigned int i = 0; i < sizeof(photonEnergy6)/sizeof(photonEnergy6[0]); i++) photonEnergy6[i] = 1240./waveLength4[i]*eV;
  G4double reEmit4[51] = {0.095367125872878,0.100714589065125,0.13933792816855,0.176288041362039,0.211928140548633,0.245505005810625,0.273370328546945,0.30221148342213,0.333093045311687,0.364223283871252,0.399221677160171,
  0.448451857613,0.514706458856799,0.566173063784637,0.587186163465084,0.619912730325973,0.657430320109198,0.6974131402047,0.747799923939089,0.82411878610154,0.886968238152526,0.927207898674191,0.972182090856316,
  0.995579234732068,0.973067015983444,0.930457495931142,0.880606822664751,0.837907842014882,0.808317266222523,0.778775703403637,0.739211827896116,0.774816618055286,0.818288494562106,0.86029367303487,
  0.836373519169372,0.783204456572165,0.716723640486035,0.653897573842595,0.600839623111776,0.539529302980934,0.482274276482319,0.45143792246339,0.403951196810372,0.355312021195194,0.308055666198011,
  0.257497891356571,0.199434040798557,0.131499952640494,0.074850729151259,0.038245233587249,0.007817590412588}; // (jakob) https://omlc.org/spectra/PhotochemCAD/html/044.html
  //G4double ppckovEmit[8] = {2.95*eV,2.95*eV,2.95*eV,2.95*eV,2.6401*eV,3.0402*eV,3.5403*eV,3.8404*eV}; 
  //G4double rindexWLS[8] = {1.5,1.5,1.5,1.5,1.504,1.505,1.515,1.52};

  MPT_WLSCoat->AddProperty("WLSCOMPONENT", photonEnergy6, reEmit4, 51);
  MPT_WLSCoat->AddConstProperty("WLSTIMECONSTANT", 2.*ns); // More or less it should be this value
  MPT_WLSCoat->AddProperty("RINDEX", pmma_refl_index_en, pmma_rind, pmma_mpt_entr)->SetSpline(true); //??

  WLS_Coat->SetMaterialPropertiesTable(MPT_WLSCoat);

  //------------------------------------------------------------------------------
  //----------------------------- Air -----------------------------
  //------------------------------------------------------------------------------
  G4double photonEnergy_Air[2] = {2.*eV, 5.*eV};
  G4double refractiveIndex_Air[2] = {1.00, 1.00};
  G4MaterialPropertiesTable* MPT_Air = new G4MaterialPropertiesTable();
  MPT_Air->AddProperty("RINDEX", photonEnergy_Air, refractiveIndex_Air, 2);
  air->SetMaterialPropertiesTable(MPT_Air);
  
  //------------------------------------------------------------------------------
  //----------------------------- Silicon Resin  -----------------------------
  //------------------------------------------------------------------------------
  G4double photonEnergy_ResinSi[3] = {2.*eV, 4.*eV, 5.*eV};
  G4double refractiveIndex_ResinSi[3] = {1.57, 1.57, 1.57}; //Hamamatsu data sheet
  G4MaterialPropertiesTable* MPT_ResinSi = new G4MaterialPropertiesTable();
  MPT_ResinSi->AddProperty("RINDEX", photonEnergy_ResinSi, refractiveIndex_ResinSi, 3);
  ResinSi->SetMaterialPropertiesTable(MPT_ResinSi);
  
  //------------------------------------------------------------------------------
  //----------------------------- Silicon  -----------------------------
  //------------------------------------------------------------------------------
  G4double photonEnergy_Si[3] = {2.76*eV, 2.88*eV, 3.10*eV};
  G4double refractiveIndex_Si[3] = {2.59, 2.75, 2.91}; //from article  2002.04218 
  G4MaterialPropertiesTable* MPT_Si = new G4MaterialPropertiesTable(); 
  MPT_Si->AddProperty("RINDEX", photonEnergy_Si, refractiveIndex_Si, 3);
  Silicon->SetMaterialPropertiesTable(MPT_Si);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineSurfaces()
{
  //------------------------------------------------------------------------------
  //----------------------------- Steel -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* SteelBoxSurface = new G4OpticalSurface("SteelBoxSurface");
  SteelBoxSurface->SetType(dielectric_metal);
  SteelBoxSurface->SetFinish(ground);
  SteelBoxSurface->SetModel(glisur);

  G4double waveLength5[59] = {838.5669979,826.7431368,817.0053925,807.2372605,798.1741236,788.399914,776.5243938,767.4627763,757.6824892,747.8748532,738.1036824,726.9149246,718.526395,707.3437148,697.5649471,687.7892181,677.303936,668.2164889,657.7874241,647.2611185,636.7864721,627.7157382,618.6146166,607.4562465,599.772712,587.189158,578.1077884,567.6164287,557.8255059,548.0406607,537.5538592,528.4648927,516.5893725,508.2038817,497.019682,487.9413512,477.4484722,466.2657919,459.9724955,447.3980577,437.6132125,427.8283673,417.3446045,408.2571574,398.4692734,388.6829088,377.4409725,368.0952298,358.264803,347.7536913,337.9627685,326.8043984,317.7762073,308.0445406,296.8557828,287.7409868,277.2435495,267.4997277,257.8060456};
  
  G4double photonEnergy7[59];
  for(unsigned int i = 0; i < sizeof(photonEnergy7)/sizeof(photonEnergy7[0]); i++) photonEnergy7[i] = 1240./waveLength5[i]*eV;

  G4double specular_steel[59] = {0.30216675359479,0.255752492579892,0.211988869838342,0.189438738725759,0.169538407342305,0.156267606546985,0.153603551566147,0.131052581909346,0.123084734749157,0.145608032446762,0.131012331786877,0.137627607123303,0.137614190415815,0.134950135434977,0.130960342045356,0.120340179523384,0.124304816586462,0.125617138287756,0.073899085093508,0.113658659193729,0.105689973489324,0.089768534420655,0.084451325533816,0.0764868325505,0.075146000345791,0.072481945364954,0.073794267066248,0.079083803993891,0.084374179465752,0.083034185805262,0.083020769097773,0.084333090799066,0.081668197274011,0.077676726795954,0.075012671815117,0.072347778290061,0.076311576808922,0.074972421692649,0.076285581938161,0.072294950004322,0.072284048929487,0.069617478315996,0.069604900152724,0.068265745036452,0.068252328328963,0.071960349606668,0.121685701224445,0.28143546343462,0.268712935986641,0.224746728193772,0.160625740904492,0.109660495475257,0.060320585725261,0.031405124642217,0.02550943578365,0.03241234077211,0.039010130683933,0.029993785627795,0.019391300509685};  //specular reflectivity Patrick 
  G4MaterialPropertiesTable *MPTsurf_Steel = new G4MaterialPropertiesTable();
  MPTsurf_Steel->AddProperty("REFLECTIVITY", photonEnergy7, specular_steel, 59);
  SteelBoxSurface->SetMaterialPropertiesTable(MPTsurf_Steel);

  G4LogicalSkinSurface* Surface_Steel = new G4LogicalSkinSurface("SteelSurface", SteelBox_log, SteelBoxSurface); //surface between the coating layer (that need to be filled with lab_ppo and steel)

  //------------------------------------------------------------------------------
  //----------------------------- BaSO4 (reflective coating) -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* BaSO4_surface = new G4OpticalSurface("BaSO4_surface");
  BaSO4_surface->SetType(dielectric_metal);
  BaSO4_surface->SetFinish(ground);
  BaSO4_surface->SetModel(glisur);
  
  G4double waveLength6[59] = {837.9882615,828.5013155,818.487317,807.9462659,798.45932,787.9182689,778.4313229,767.8902719,758.9303785,748.3893274,739.429434,727.8342778,717.8202793,709.3874384,697.2652297,687.2512312,678.8183903,668.8043918,658.7903933,648.7763948,638.2353437,627.6942927,618.7343992,609.2474533,598.7064022,589.2194563,580.2595628,568.1373541,558.1233556,546.0011469,537.568306,528.6084126,518.0673615,508.053363,499.6205222,487.4983134,478.0113675,467.997369,458.510423,447.9693719,439.5365311,427.4143223,418.4544289,408.4404304,399.480537,387.8853808,378.3984349,368.3844363,357.8433853,347.8293868,338.3424408,327.8013897,318.3144438,308.8274978,299.3405518,287.7453957,278.7855023,269.8256089,258.2304527}; 
   
  G4double photonEnergy8[59];
  for(unsigned int i = 0; i < sizeof(photonEnergy8)/sizeof(photonEnergy8[0]); i++) photonEnergy8[i] = 1240./waveLength6[i]*eV;
  G4double photonEnergy9[2] = {1*eV, 5*eV};
  G4double specular_coating[59] = { 0.01827689,0.006247856,0.016558456,0.00796629,0.002810989,0.011403156,0.00796629,0.01312159,0.004529423,0.01312159,0,0.014840023,0.00796629,0,0,0.011985604,0.003428249,0.008562662,0.008562662,0.008562662,0.003428249,5.30684E-06,0,0.010274133,0.010942009,0,0.010942009,0.009231693,0.000680113,0.010942009,0,0.004100745,0,0.005811061,0.000680113,0.002699077,0.007819659,0.01123338,0.000992216,0.000992216,0,0.007819659,0,0.004405937,0.006781677,0.003370253,0.006781677,0.008487389,0.00857075,0.010279258,0,0,0,0,0,0.001514251,0.001514251,0,0,}; // specular reflection Patrick data, calculated as 1- diffuse reflection since the specular was very small and the measurements for specular are more difficult to separete from the diffuse ones.

  G4double other_coating[2] = {0, 0}; // it is not relevant in our case, it is for not smooth surfaces 
  
  G4double p_coating_refl[19] = {0.653775342148272*eV,0.687758866578874*eV,0.729108650377501*eV,0.773110023703728*eV,0.827676607343705*eV,0.887209379832684*eV,0.954680932864894*eV,1.03471853868547*eV,1.12204531460793*eV,1.24533694371132*eV,1.3767339264975*eV,1.55649354301667*eV,1.77189963551006*eV,2.09375317772178*eV,2.51208306681883*eV,3.16797881984726*eV,3.50425046717685*eV,3.63215740198218*eV,3.81157924082758*eV};
   
  G4double refl_coating[19] = {0.861320132013201,0.916105610561056,0.922706270627063,0.920726072607261,0.916765676567657,
   0.931287128712871,0.95042904290429,0.958349834983498,0.962970297029703,0.968250825082508,0.969570957095709,
   0.973531353135314,0.977491749174918,0.978151815181518,0.98013201320132,0.980792079207921,0.970891089108911,
   0.960990099009901,0.953729372937294}; // reflectivity of the coating https://www.optopolymer.de/produktuebersicht/diffuse-reflecting-materials/bariumsulfate-baso4-coating-oprc/
       
  G4MaterialPropertiesTable* MTP_BaSO = new G4MaterialPropertiesTable();
  MTP_BaSO->AddProperty("SPECULARLOBECONSTANT", photonEnergy8, specular_coating, 59); // In order to have diffuse reclectivity (Lambertian), it is necessary define all the other three.
  MTP_BaSO->AddProperty("SPECULARSPIKECONSTANT", photonEnergy9, other_coating, 2); //  The diffuse is 1-other three (in this case 1-sppecular). 
  MTP_BaSO->AddProperty("BACKSCATTERCONSTANT", photonEnergy9, other_coating, 2);
  
  G4double refl_test[19];
  for(int i = 0; i < 19; ++i) refl_test[i] = 0.65*refl_coating[i];
 
  MTP_BaSO->AddProperty("REFLECTIVITY", p_coating_refl, refl_coating, 19);

  BaSO4_surface->SetMaterialPropertiesTable(MTP_BaSO); 	
   
  G4LogicalSkinSurface* Surface_Refl = new G4LogicalSkinSurface("BaSO4_surface", ReflectBox_log, BaSO4_surface);

  //------------------------------------------------------------------------------
  //----------------------------- Absorbent surface SiPMs -----------------------------
  //------------------------------------------------------------------------------
  G4OpticalSurface* sipms_surface = new G4OpticalSurface("sipms_surface");
  sipms_surface->SetType(dielectric_metal);
  sipms_surface->SetFinish(ground);
  sipms_surface->SetModel(glisur);
  
  G4double waveLength7[2] = {1000, 100}; 
   
  G4double photonEnergy10[2];
  for(unsigned int i = 0; i < sizeof(photonEnergy10)/sizeof(photonEnergy10[0]); i++) photonEnergy10[i] = 1240./waveLength7[i]*eV;
  G4double abs_sipm[2] = {0, 0}; // we want that the surface of the sipms sides is completely absorbent
       
  G4MaterialPropertiesTable* MTP_sipms = new G4MaterialPropertiesTable();
  MTP_sipms->AddProperty("REFLECTIVITY", photonEnergy10, abs_sipm, 2);

  sipms_surface->SetMaterialPropertiesTable(MTP_sipms); 	
   
  G4LogicalSkinSurface* Surface_abs_sipms = new G4LogicalSkinSurface("sipms_surface", sipmSens_log, sipms_surface);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineSolids()
{
  expHall_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------

  double x = 800*mm;
  double ybr = 1035*mm; // bottom right
  double ytr = 1179*mm; // top right
  double ytl = 1287*mm; // top left
  double ybl = 1013*mm; // bottom left
  double ytm = 1233*mm; // top middle
  double ybm = 1024*mm; // bottom middle
  double slopet = 0.0675; // slope of top edge
  double slopeb = 0.01375; // slope of bottom edge
  double slopem = 0.04; // slope of middle beam
  double holeR = 15*mm; // radius of holes in beams
  
  // Steel box
  std::vector<G4TwoVector> det = {G4TwoVector(-x, -ybl), G4TwoVector(-x, ytl), G4TwoVector(x, ytr), G4TwoVector(x, -ybr), G4TwoVector(-x, -ybl), G4TwoVector(-x, ytl), G4TwoVector(x, ytr), G4TwoVector(x, -ybr)};
  SteelBox = new G4GenericTrap("SteelBox", SteelZ/2., det);
  std::vector<G4TwoVector> refl = 
  {G4TwoVector(-x+WallThick, -ybl+WallThick), G4TwoVector(-x+WallThick, ytl-WallThick), G4TwoVector(x-WallThick, ytr-WallThick), G4TwoVector(x-WallThick, -ybr+WallThick),
  G4TwoVector(-x+WallThick, -ybl+WallThick), G4TwoVector(-x+WallThick, ytl-WallThick), G4TwoVector(x-WallThick, ytr-WallThick), G4TwoVector(x-WallThick, -ybr+WallThick)};
  ReflectBox = new G4GenericTrap("ReflectBox", ReflectZ/2., refl);
  EmptySteelBox = new G4SubtractionSolid("EmptySteelBox", SteelBox, ReflectBox, 0, G4ThreeVector(0,0,0));
  
  // Steel vertical beam
  std::vector<G4TwoVector> beamv = 
  {G4TwoVector(-WallThick/2., -ybm+WallThick+slopeb*(WallThick/2.)), G4TwoVector(-WallThick/2., ytm-WallThick+slopet*(WallThick/2.)), 
  G4TwoVector(WallThick/2., ytm-WallThick-slopet*(WallThick/2.)),    G4TwoVector(WallThick/2., -ybm+WallThick-slopeb*(WallThick/2.)),
  G4TwoVector(-WallThick/2., -ybm+WallThick+slopeb*(WallThick/2.)),  G4TwoVector(-WallThick/2. ,ytm-WallThick+slopet*(WallThick/2.)), 
  G4TwoVector(WallThick/2., ytm-WallThick-slopet*(WallThick/2.)),    G4TwoVector(WallThick/2., -ybm+WallThick-slopeb*(WallThick/2.))};
  SteelBeamV = new G4GenericTrap("SteelBeamV", ReflectZ/2., beamv);
  // Holes in vertical beam
  G4Tubs* BeamHoleV = new G4Tubs("BeamHoleV", 0*mm, holeR, WallThick/2. + ReflectThick + 10*mm, 0*deg, 360*deg);
  std::vector<G4ThreeVector> SteelBeamV_Holes = 
  {G4ThreeVector(0, -ybm+WallThick, ReflectZ/2.), G4ThreeVector(0, -ybm+WallThick, -ReflectZ/2.),
  G4ThreeVector(0, ytm-WallThick, ReflectZ/2.),   G4ThreeVector(0, ytm-WallThick, -ReflectZ/2.)};
  G4RotationMatrix* RotV = new G4RotationMatrix(); RotV->rotateY(pi/2);
  std::vector<G4SubtractionSolid*> SteelBeamV_tempvec;
  for(unsigned int pos = 0; pos < SteelBeamV_Holes.size(); pos++) {
    if(pos == 0) SteelBeamV_tempvec.push_back(new G4SubtractionSolid((std::string("SteelBeamV_")+std::to_string(pos)).c_str(), SteelBeamV, BeamHoleV, RotV, SteelBeamV_Holes[0]));
    else SteelBeamV_tempvec.push_back(new G4SubtractionSolid((std::string("SteelBeamV_")+std::to_string(pos)).c_str(), SteelBeamV_tempvec[pos-1], BeamHoleV, RotV, SteelBeamV_Holes[pos]));
  }
  SteelBeamVWithHoles = SteelBeamV_tempvec[SteelBeamV_Holes.size()-1];
  
  // Steel horizontal Beam
  std::vector<G4TwoVector> beamh = 
  {G4TwoVector(-x+WallThick, slopem*(x-WallThick)-WallThick/2.), G4TwoVector(-x+WallThick, slopem*(x-WallThick)+WallThick/2.),
  G4TwoVector(x-WallThick, -slopem*(x-WallThick)+WallThick/2.),  G4TwoVector(x-WallThick, -slopem*(x-WallThick)-WallThick/2.),
  G4TwoVector(-x+WallThick, slopem*(x-WallThick)-WallThick/2.),  G4TwoVector(-x+WallThick, slopem*(x-WallThick)+WallThick/2.),
  G4TwoVector(x-WallThick, -slopem*(x-WallThick)+WallThick/2.),  G4TwoVector(x-WallThick, -slopem*(x-WallThick)-WallThick/2.)};
  SteelBeamH = new G4GenericTrap("SteelBeamH",ReflectZ/2.,beamh);
  // Holes in horizontal beam
  G4CutTubs* BeamHoleH = new G4CutTubs("BeamHoleH", 0*mm, holeR, WallThick/2.+ReflectThick+10*mm, 0*deg, 360*deg, G4ThreeVector(-0.03996, 0, -0.999), G4ThreeVector(0.03996, 0, 0.999));
  std::vector<G4ThreeVector> SteelBeamH_Holes = 
  {G4ThreeVector(x-WallThick, -slopem*(x-WallThick), ReflectZ/2.),   G4ThreeVector(-x+WallThick, slopem*(x-WallThick), ReflectZ/2.),
  G4ThreeVector(x-WallThick, -slopem*(x-WallThick), -ReflectZ/2.),   G4ThreeVector(-x+WallThick, slopem*(x-WallThick), -ReflectZ/2.),
  G4ThreeVector(WallThick/2., -slopem*(WallThick/2.), ReflectZ/2.),  G4ThreeVector(-WallThick/2., slopem*(WallThick/2.), ReflectZ/2.),
  G4ThreeVector(WallThick/2., -slopem*(WallThick/2.), -ReflectZ/2.), G4ThreeVector(-WallThick/2., slopem*(WallThick/2.), -ReflectZ/2.)};
  G4RotationMatrix* RotH = new G4RotationMatrix(0, -pi/2, 0);
  std::vector<G4SubtractionSolid*> SteelBeamH_tempvec;
  for(unsigned int pos = 0; pos < SteelBeamH_Holes.size(); pos++) {
    if(pos == 0) SteelBeamH_tempvec.push_back(new G4SubtractionSolid((std::string("SteelBeamH_")+std::to_string(pos)).c_str(), SteelBeamH, BeamHoleH, RotH, SteelBeamH_Holes[0]));
    else SteelBeamH_tempvec.push_back(new G4SubtractionSolid((std::string("SteelBeamH_")+std::to_string(pos)).c_str(), SteelBeamH_tempvec[pos-1], BeamHoleH, RotH, SteelBeamH_Holes[pos]));
  }
  SteelBeamHWithHoles = SteelBeamH_tempvec[SteelBeamH_Holes.size()-1];
  // Unified beams
  SteelBeam = new G4UnionSolid("SteelBeam", SteelBeamVWithHoles, SteelBeamHWithHoles, 0, G4ThreeVector(0, 0, 0));
  
  // Reflectivity box
  EmptyReflectBox = new G4SubtractionSolid("EmptyReflectBox", ReflectBox, SteelBeam, 0, G4ThreeVector(0, 0, 0));

  // Scintillator box left top
  std::vector<G4TwoVector> scintLT = 
  {G4TwoVector(-x+WallThick+ReflectThick, slopem*(x-WallThick-ReflectThick)+WallThick/2.+ReflectThick),   G4TwoVector(-x+WallThick+ReflectThick, ytl-WallThick-ReflectThick),
  G4TwoVector(-WallThick/2.-ReflectThick, ytm-WallThick-ReflectThick+slopet*(WallThick/2.+ReflectThick)), G4TwoVector(-WallThick/2.-ReflectThick, WallThick/2.+ReflectThick+slopem*(WallThick/2.+ReflectThick)),
  G4TwoVector(-x+WallThick+ReflectThick, slopem*(x-WallThick-ReflectThick)+WallThick/2.+ReflectThick),    G4TwoVector(-x+WallThick+ReflectThick, ytl-WallThick-ReflectThick),
  G4TwoVector(-WallThick/2.-ReflectThick, ytm-WallThick-ReflectThick+slopet*(WallThick/2.+ReflectThick)), G4TwoVector(-WallThick/2.-ReflectThick, WallThick/2.+ReflectThick+slopem*(WallThick/2.+ReflectThick))};
  ScintillatorBoxLT = new G4GenericTrap("ScintillatorBoxLT", SctZ/2., scintLT);

  // Scintillator box right top
  std::vector<G4TwoVector> scintRT = 
  {G4TwoVector(WallThick/2.+ReflectThick, WallThick/2.+ReflectThick-slopem*(WallThick/2.+ReflectThick)), G4TwoVector(WallThick/2.+ReflectThick, ytm-WallThick-ReflectThick-slopet*(WallThick/2.+ReflectThick)),
  G4TwoVector(x-WallThick-ReflectThick, ytr-WallThick-ReflectThick),                                     G4TwoVector(x-WallThick-ReflectThick, -slopem*(x-WallThick-ReflectThick)+WallThick/2.+ReflectThick),
  G4TwoVector(WallThick/2.+ReflectThick, WallThick/2.+ReflectThick-slopem*(WallThick/2.+ReflectThick)),  G4TwoVector(WallThick/2.+ReflectThick, ytm-WallThick-ReflectThick-slopet*(WallThick/2.+ReflectThick)),
  G4TwoVector(x-WallThick-ReflectThick, ytr-WallThick-ReflectThick),                                     G4TwoVector(x-WallThick-ReflectThick, -slopem*(x-WallThick-ReflectThick)+WallThick/2.+ReflectThick)};
  ScintillatorBoxRT = new G4GenericTrap("ScintillatorBoxRT", SctZ/2., scintRT);

  // Scintillator box left bottom
  std::vector<G4TwoVector> scintLB = 
  {G4TwoVector(-x+WallThick+ReflectThick, -ybl+WallThick+ReflectThick),                                   G4TwoVector(-x+WallThick+ReflectThick, slopem*(x-WallThick-ReflectThick)-WallThick/2.-ReflectThick),
  G4TwoVector(-WallThick/2.-ReflectThick, -WallThick/2.-ReflectThick+slopem*(WallThick/2.+ReflectThick)), G4TwoVector(-WallThick/2.-ReflectThick, -ybm+WallThick+ReflectThick+slopeb*(WallThick/2.+ReflectThick)),
  G4TwoVector(-x+WallThick+ReflectThick, -ybl+WallThick+ReflectThick),                                    G4TwoVector(-x+WallThick+ReflectThick, slopem*(x-WallThick-ReflectThick)-WallThick/2.-ReflectThick),
  G4TwoVector(-WallThick/2.-ReflectThick, -WallThick/2.-ReflectThick+slopem*(WallThick/2.+ReflectThick)), G4TwoVector(-WallThick/2.-ReflectThick, -ybm+WallThick+ReflectThick+slopeb*(WallThick/2.+ReflectThick))};
  ScintillatorBoxLB = new G4GenericTrap("ScintillatorBoxLB", SctZ/2., scintLB);
  
  // Scintillator box right bottom
  std::vector<G4TwoVector> scintRB =
  {G4TwoVector(WallThick/2.+ReflectThick, -ybm+WallThick+ReflectThick-slopeb*(WallThick/2.+ReflectThick)), G4TwoVector(WallThick/2.+ReflectThick, -WallThick/2.-ReflectThick-slopem*(WallThick/2.+ReflectThick)),
  G4TwoVector(x-WallThick-ReflectThick, -slopem*(x-WallThick-ReflectThick)-WallThick/2.-ReflectThick),     G4TwoVector(x-WallThick-ReflectThick, -ybr+WallThick+ReflectThick),
  G4TwoVector(WallThick/2.+ReflectThick, -ybm+WallThick+ReflectThick-slopeb*(WallThick/2.+ReflectThick)),  G4TwoVector(WallThick/2.+ReflectThick, -WallThick/2.-ReflectThick-slopem*(WallThick/2.+ReflectThick)),
  G4TwoVector(x-WallThick-ReflectThick, -slopem*(x-WallThick-ReflectThick)-WallThick/2.-ReflectThick),     G4TwoVector(x-WallThick-ReflectThick, -ybr+WallThick+ReflectThick)};
  ScintillatorBoxRB = new G4GenericTrap("ScintillatorBoxRB", SctZ/2., scintRB);

  // Scintillator in beam holes
  std::vector<G4ThreeVector> ScintillatorV_Holes = 
  {G4ThreeVector(0, -ybm+WallThick+ReflectThick, SctZ/2.), G4ThreeVector(0, -ybm+WallThick+ReflectThick, -SctZ/2.),
  G4ThreeVector(0, ytm-WallThick-ReflectThick, SctZ/2.),   G4ThreeVector(0, ytm-WallThick-ReflectThick, -SctZ/2.)};
  std::vector<G4Transform3D> TrV =
  {G4Transform3D(G4RotationMatrix(0, pi/2, pi/2),  ScintillatorV_Holes[2]), G4Transform3D(G4RotationMatrix(0, -pi/2, pi/2),  ScintillatorV_Holes[3]),
  G4Transform3D(G4RotationMatrix(0, pi/2, -pi/2),  ScintillatorV_Holes[0]), G4Transform3D(G4RotationMatrix(0, -pi/2, -pi/2), ScintillatorV_Holes[1])};
  G4double gap = 1*mm;
  std::vector<G4ThreeVector> ScintillatorH_Holes = 
  {G4ThreeVector(x-WallThick-ReflectThick-gap, -slopem*(x-WallThick-ReflectThick-gap), SctZ/2.-gap),   G4ThreeVector(-x+WallThick+ReflectThick+gap, slopem*(x-WallThick-ReflectThick-gap), SctZ/2.-gap),
  G4ThreeVector(x-WallThick-ReflectThick-gap, -slopem*(x-WallThick-ReflectThick-gap), -SctZ/2.+gap),   G4ThreeVector(-x+WallThick+ReflectThick+gap, slopem*(x-WallThick-ReflectThick-gap), -SctZ/2.+gap),
  G4ThreeVector(WallThick/2.+ReflectThick+gap, -slopem*(WallThick/2.+ReflectThick+gap), SctZ/2.-gap),  G4ThreeVector(-WallThick/2.-ReflectThick-gap, slopem*(WallThick/2.+ReflectThick+gap), SctZ/2.-gap),
  G4ThreeVector(WallThick/2.+ReflectThick+gap, -slopem*(WallThick/2.+ReflectThick+gap), -SctZ/2.+gap), G4ThreeVector(-WallThick/2.-ReflectThick-gap, slopem*(WallThick/2.+ReflectThick+gap), -SctZ/2.+gap)};
  std::vector<G4Transform3D> TrH1 =
  {G4Transform3D(G4RotationMatrix(0,pi/2,pi),ScintillatorH_Holes[0]), G4Transform3D(G4RotationMatrix(0,pi/2,0),ScintillatorH_Holes[1]),
  G4Transform3D(G4RotationMatrix(0,pi/2,0),ScintillatorH_Holes[4]),   G4Transform3D(G4RotationMatrix(0,pi/2,pi),ScintillatorH_Holes[5])}; 
  std::vector<G4Transform3D> TrH2 =
  {G4Transform3D(G4RotationMatrix(0, -pi/2, pi),ScintillatorH_Holes[2]), G4Transform3D(G4RotationMatrix(0, -pi/2, 0),ScintillatorH_Holes[3]),
  G4Transform3D(G4RotationMatrix(0, -pi/2, 0),ScintillatorH_Holes[6]),   G4Transform3D(G4RotationMatrix(0, -pi/2, pi),ScintillatorH_Holes[7])};
  G4Transform3D TrV1 = G4Transform3D(G4RotationMatrix(0, 0, 0),G4ThreeVector(0, ReflectThick, 0));
  G4Transform3D TrV2 = G4Transform3D(G4RotationMatrix(0, 0, 0),G4ThreeVector(0, -ReflectThick, 0));
  std::vector<G4TwoVector> beamvrefl =
  {G4TwoVector(-WallThick/2.-ReflectThick, -ybm+WallThick+slopeb*(WallThick/2.+ReflectThick)), G4TwoVector(-WallThick/2.-ReflectThick, ytm-WallThick+slopet*(WallThick/2.+ReflectThick)),
  G4TwoVector(WallThick/2.+ReflectThick, ytm-WallThick-slopet*(WallThick/2.+ReflectThick)),    G4TwoVector(WallThick/2.+ReflectThick, -ybm+WallThick-slopeb*(WallThick/2.+ReflectThick)),
  G4TwoVector(-WallThick/2.-ReflectThick, -ybm+WallThick+slopeb*(WallThick/2.+ReflectThick)), G4TwoVector(-WallThick/2.-ReflectThick, ytm-WallThick+slopet*(WallThick/2.+ReflectThick)),
  G4TwoVector(WallThick/2.+ReflectThick, ytm-WallThick-slopet*(WallThick/2.+ReflectThick)),    G4TwoVector(WallThick/2.+ReflectThick, -ybm+WallThick-slopeb*(WallThick/2.+ReflectThick))};
  G4GenericTrap* SteelBeamVWithRefl = new G4GenericTrap("SteelBeamVWithRefl", ReflectZ/2., beamvrefl);
  G4Tubs* ScintHoleV = new G4Tubs("ScintHoleV", 0*mm, holeR-(ReflectThick*2), WallThick/2.+ReflectThick, 0*deg, 360*deg);
  G4CutTubs* ScintillatorHoleH1 = new G4CutTubs("BeamHoleH1", 0*mm, holeR-(ReflectThick*2)-gap*2, WallThick/2.+ReflectThick, 0*deg, 90*deg, G4ThreeVector(-0.03996, 0, -0.999), G4ThreeVector(0.03996, 0, 0.999));
  G4CutTubs* ScintillatorHoleH2 = new G4CutTubs("BeamHoleH2", 0*mm, holeR-(ReflectThick*2)-gap*2, WallThick/2.+ReflectThick, 0*deg, 90*deg, G4ThreeVector(0.03996, 0, -0.999), G4ThreeVector(-0.03996, 0, 0.999));
  std::vector<G4IntersectionSolid*> ScintillatorHoleV;
  for(unsigned int pos = 0; pos < ScintillatorV_Holes.size(); pos++) 
  ScintillatorHoleV.push_back(new G4IntersectionSolid((std::string("ScintillatorHoleV_")+std::to_string(pos)).c_str(), SteelBeamVWithRefl, ScintHoleV, RotV, ScintillatorV_Holes[pos]));
  ScintillatorInBeams = new G4MultiUnion("ScintillatorInBeams");
  for(unsigned int pos = 0; pos < ScintillatorV_Holes.size()/2; pos++) ScintillatorInBeams->AddNode(*ScintillatorHoleV[pos], TrV1);
  for(unsigned int pos = ScintillatorV_Holes.size()/2; pos < ScintillatorV_Holes.size(); pos++) ScintillatorInBeams->AddNode(*ScintillatorHoleV[pos], TrV2);
  for(unsigned int pos = 0; pos < TrH1.size(); pos++) { ScintillatorInBeams->AddNode(*ScintillatorHoleH1, TrH1[pos]); ScintillatorInBeams->AddNode(*ScintillatorHoleH2, TrH2[pos]); }
  ScintillatorInBeams->Voxelize();
  
  G4double Rin, Rout;

  // Outer tube
  Rin = Diam_Out_In/2;
  Rout = Diam_Out_Out/2;
  Outer_tube  = new G4Tubs("Outer_tube", Rin, Rout, Length_Out/2, 0, 360*deg);
  // Air gap out
  Rin = Diam_WOM_Out/2 + Thickness_WLS;
  Rout = Diam_Out_In/2;
  Air_gap_out  = new G4Tubs("Air_gap_out", Rin, Rout, (Length_WOM+Thickness_Gap)/2, 0, 360*deg);
  // WLS tube out
  Rin = Diam_WOM_Out/2;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  WLS_tube_out  = new G4Tubs("WLS_tube_out", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // WOM tube
  Rin = Diam_WOM_In/2;
  Rout = Diam_WOM_Out/2;
  WOM_tube  = new G4Tubs("WOM_tube", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // WLS tube in
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = Diam_WOM_In/2;
  WLS_tube_in  = new G4Tubs("WLS_tube_in", Rin, Rout, Length_WOM/2, 0, 360*deg);
  // Air gap in
  Rin = Diam_In_Out/2;
  Rout = Diam_WOM_In/2 - Thickness_WLS;
  Air_gap_in  = new G4Tubs("Air_gap_in", Rin, Rout, (Length_WOM+Thickness_Gap)/2, 0, 360*deg);
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
  // Air ring outer
  Rin = (Diam_WOM_Out - Thickness_SupportRing)/2;
  Rout = Diam_WOM_Out/2 + Thickness_WLS;
  Air_ring_out = new G4Tubs("Air_ring_out", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // PMMA ring supporting WOM
  Rin = (Diam_WOM_In + Thickness_SupportRing)/2;
  Rout = (Diam_WOM_Out - Thickness_SupportRing)/2;
  PMMA_ring_lower = new G4Tubs("PMMA_ring_lower", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // Air ring inner
  Rin = Diam_WOM_In/2 - Thickness_WLS;
  Rout = (Diam_WOM_In + Thickness_SupportRing)/2; 
  Air_ring_in = new G4Tubs("Air_ring_in", Rin, Rout, Thickness_Gap/2, 0, 360*deg);
  // Hole in box
  Rin = 0.0*mm;
  Rout = Diam_Out_Out/2; 
  Hole_box  = new G4Tubs("Hole_box", Rin, Rout, (Length_Out + Thickness_Ring)/2, 0, 360*deg);
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
  SctInside = new G4Tubs("Sct_inside", Rin, Rout, (Length_In + Thickness_Ring)/2, 0, 360*deg);

  G4double delta_Z_Hole = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_Out/2 + Thickness_Ring/2;

  std::vector<G4SubtractionSolid*> EmptySteelBoxWithHole_tempvec;
  std::vector<G4SubtractionSolid*> EmptyReflectBoxWithHole_tempvec;
  std::vector<G4SubtractionSolid*> ScintillatorBoxWithHoleLT_tempvec;
  std::vector<G4SubtractionSolid*> ScintillatorBoxWithHoleRT_tempvec;
  std::vector<G4SubtractionSolid*> ScintillatorBoxWithHoleLB_tempvec;
  std::vector<G4SubtractionSolid*> ScintillatorBoxWithHoleRB_tempvec;
  
  for(unsigned int pos = 0; pos<WOM_coord_vec.size(); pos++) {
    if(pos == 0) {
      EmptySteelBoxWithHole_tempvec.push_back(new G4SubtractionSolid((std::string("EmptySteelBoxWithHole_")+std::to_string(pos)).c_str(), EmptySteelBox, Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      EmptyReflectBoxWithHole_tempvec.push_back(new G4SubtractionSolid((std::string("EmptyReflectBoxWithHole_")+std::to_string(pos)).c_str(), EmptyReflectBox, Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      ScintillatorBoxWithHoleLT_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleLT_")+std::to_string(pos)).c_str(), ScintillatorBoxLT, Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
    }
    if(pos != 0) {
      EmptySteelBoxWithHole_tempvec.push_back(new G4SubtractionSolid((std::string("EmptySteelBoxWithHole_")+std::to_string(pos)).c_str(), EmptySteelBoxWithHole_tempvec[pos-1], Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      EmptyReflectBoxWithHole_tempvec.push_back(new G4SubtractionSolid((std::string("EmptyReflectBoxWithHole_")+std::to_string(pos)).c_str(), EmptyReflectBoxWithHole_tempvec[pos-1], Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 1) ScintillatorBoxWithHoleLT_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleLT_")+std::to_string(pos)).c_str(), ScintillatorBoxWithHoleLT_tempvec[0], Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 2) ScintillatorBoxWithHoleRT_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleRT_")+std::to_string(pos)).c_str(), ScintillatorBoxRT, Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 3) ScintillatorBoxWithHoleRT_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleRT_")+std::to_string(pos)).c_str(), ScintillatorBoxWithHoleRT_tempvec[0], Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 4) ScintillatorBoxWithHoleLB_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleLB_")+std::to_string(pos)).c_str(), ScintillatorBoxLB, Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 5) ScintillatorBoxWithHoleLB_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleLB_")+std::to_string(pos)).c_str(), ScintillatorBoxWithHoleLB_tempvec[0], Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 6) ScintillatorBoxWithHoleRB_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleRB_")+std::to_string(pos)).c_str(), ScintillatorBoxRB, Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
      if(pos == 7) ScintillatorBoxWithHoleRB_tempvec.push_back(new G4SubtractionSolid((std::string("ScintillatorBoxWithHoleRB_")+std::to_string(pos)).c_str(), ScintillatorBoxWithHoleRB_tempvec[0], Hole_box, 0, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Hole)));
    }
  }
  EmptySteelBoxWithHole = EmptySteelBoxWithHole_tempvec[WOM_coord_vec.size()-1];
  EmptyReflectBoxWithHole = EmptyReflectBoxWithHole_tempvec[WOM_coord_vec.size()-1];
  ScintillatorBoxWithHoleLT = ScintillatorBoxWithHoleLT_tempvec[1];
  ScintillatorBoxWithHoleRT = ScintillatorBoxWithHoleRT_tempvec[1];
  ScintillatorBoxWithHoleLB = ScintillatorBoxWithHoleLB_tempvec[1];
  ScintillatorBoxWithHoleRB = ScintillatorBoxWithHoleRB_tempvec[1];

  // SIPMs 14161 - 3050HS
  sipmHole = new G4Box("sipmHole", sipmSizeSens/2., sipmSizeSens/2., (sipmSensThickness + sipmSensThicknessTop)/2);
  sipmSensTop = new G4Box("sipmSensTop", sipmSizeSens/2., sipmSizeSens/2., sipmSensThicknessTop/2.);
  sipmSens = new G4Box("sipmSens", sipmSizeSens/2., sipmSizeSens/2., sipmSensThickness/2.);
  sipmWindowAll = new G4Box("sipmAll", sipmSize/2., sipmSize/2., sipmWindowThickness);
  sipmWindow = new G4SubtractionSolid("sipmWindow", sipmWindowAll, sipmHole, 0, G4ThreeVector(0, 0, (sipmSensThickness + sipmSensThicknessTop)/2)); 
  sipmBaseBox = new G4Box("sipmBase", sipmSize/2., sipmSize/2., sipmBaseThickness/2.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::DefineLogicalVolumes()
{ 
  expHall_log = new G4LogicalVolume(expHall_box, air, "WorldLV", 0, 0, 0);
  SteelBox_log = new G4LogicalVolume(EmptySteelBoxWithHole, steel, "SteelBoxLV");
  SteelBeam_log = new G4LogicalVolume(SteelBeam, steel, "SteelBeamLV");
  ReflectBox_log = new G4LogicalVolume(EmptyReflectBoxWithHole, steel, "SteelBeamLV");
  ScintillatorBoxLT_log = new G4LogicalVolume(ScintillatorBoxWithHoleLT, LAB_PPO, "ScintillatorBoxLTLV");
  ScintillatorBoxRT_log = new G4LogicalVolume(ScintillatorBoxWithHoleRT, LAB_PPO, "ScintillatorBoxRTLV");
  ScintillatorBoxLB_log = new G4LogicalVolume(ScintillatorBoxWithHoleLB, LAB_PPO, "ScintillatorBoxLBLV");
  ScintillatorBoxRB_log = new G4LogicalVolume(ScintillatorBoxWithHoleRB, LAB_PPO, "ScintillatorBoxRBLV");
  ScintillatorInBeams_log = new G4LogicalVolume(ScintillatorInBeams, LAB_PPO, "ScintillatorInBeamsLV");
  Outer_tube_log = new G4LogicalVolume(Outer_tube, PMMA_side, "Outer_tubeLV");
  WOM_tube_log = new G4LogicalVolume(WOM_tube, PMMA_bottom, "WOM_tubeLV");
  Inner_tube_log = new G4LogicalVolume(Inner_tube, PMMA_side, "Inner_tubeLV");
  PMMA_Ring_log = new G4LogicalVolume(PMMA_Ring, PMMA_bottom, "PMMA_RingLV");
  PMMA_disk_log = new G4LogicalVolume(PMMA_disk, PMMA_side, "PMMA_diskLV");
  Air_gap_out_log = new G4LogicalVolume(Air_gap_out, air, "Air_gap_outLV");
  Air_gap_in_log = new G4LogicalVolume(Air_gap_in, air, "Air_gap_inLV");
  WLS_tube_out_log = new G4LogicalVolume(WLS_tube_out, WLS_Coat, "WLS_outLV");
  WLS_tube_in_log = new G4LogicalVolume(WLS_tube_in, air, "WLS_inLV"); //we should keep it without
  PMMA_Hat_log = new G4LogicalVolume(PMMA_Hat, PMMA_side, "PMMA_HatLV");
  Air_ring_out_log = new G4LogicalVolume(Air_ring_out, air, "Air_ring_outLV");
  Air_ring_in_log = new G4LogicalVolume(Air_ring_in, air, "Air_ring_inLV");
  PMMA_ring_lower_log = new G4LogicalVolume(PMMA_ring_lower, PMMA_bottom, "PMMA_ring_lowerLV");
  Steel_Add_log = new G4LogicalVolume(SteelAdd, steel, "Steel_AddLV");
  Sct_Inside_log = new G4LogicalVolume(SctInside, LAB_PPO, "Sct_InsideLV");
  sipmSens_log = new G4LogicalVolume(sipmSens, Silicon, "sipmSensLV");  ////?????
  sipmSensTop_log = new G4LogicalVolume(sipmSensTop, Silicon, "sipmSensTopLV");  ////?????
  sipmWindow_log = new G4LogicalVolume(sipmWindow, ResinSi, "sipmWindowLV");
  sipmBaseBox_log = new G4LogicalVolume(sipmBaseBox, Al, "sipmBaseBoxLV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstruction::ConstructVolumes()
{
  expHall_phys = new G4PVPlacement(0, G4ThreeVector(), expHall_log, "WorldPV", 0, false, 0);

  G4RotationMatrix *RM1 = new G4RotationMatrix(0*deg, 0*deg, 0*deg);

  SteelBox_phys = new G4PVPlacement(0, G4ThreeVector(), SteelBox_log, "SteelBoxPV", expHall_log, false, 101, intersect_check);
  SteelBeam_phys = new G4PVPlacement(0, G4ThreeVector(), SteelBeam_log, "SteelBeamPV", expHall_log, false, 102, intersect_check);
  ReflectBox_phys = new G4PVPlacement(0, G4ThreeVector(), ReflectBox_log, "ReflectBoxPV", expHall_log, false, 201, intersect_check);
  ScintillatorBoxLT_phys = new G4PVPlacement(0, G4ThreeVector(), ScintillatorBoxLT_log, "ScintillatorBoxLTPV", ReflectBox_log, false, 301, intersect_check);
  ScintillatorBoxRT_phys = new G4PVPlacement(0, G4ThreeVector(), ScintillatorBoxRT_log, "ScintillatorBoxRTPV", ReflectBox_log, false, 302, intersect_check);
  ScintillatorBoxLB_phys = new G4PVPlacement(0, G4ThreeVector(), ScintillatorBoxLB_log, "ScintillatorBoxLBPV", ReflectBox_log, false, 303, intersect_check);
  ScintillatorBoxRB_phys = new G4PVPlacement(0, G4ThreeVector(), ScintillatorBoxRB_log, "ScintillatorBoxRBPV", ReflectBox_log, false, 304, intersect_check);
  ScintillatorInBeams_phys = new G4PVPlacement(0, G4ThreeVector(), ScintillatorInBeams_log, "ScintillatorInBeamsPV", expHall_log, false, 305, intersect_check);

//-------------------------------------------------------------------
  G4double delta_Z_sipm = SctZ/2 + Thickness_Steel_Add + Thickness_Hat + WallThick + Thickness_Ring + Thickness_Gap + sipmWindowThickness; // no small gap of air 
//-------------------------------------------------------------------

  // PMMA Staff
  // Outer tube
  G4double delta_Z_Outer_tube = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_Out/2 + Thickness_Ring;
  // Air gaps
  G4double delta_Z_Air_gaps = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_WOM/2 + Thickness_Gap/2 + Thickness_Ring;
  // WOM tube
  G4double delta_Z_WOM = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_WOM/2 + Thickness_Ring + Thickness_Gap;
  // Inner tube
  G4double delta_Z_Inner_tube = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_In/2 + Thickness_Ring;
  // PMMA Ring
  G4double delta_Z_PMMA_Ring = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Thickness_Ring/2;
  // PMMA Disk
  G4double delta_Z_PMMA_Disk = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_In + Thickness_Ring + Thickness_Disk/2;
  // PMMA Hat
  G4double delta_Z_PMMA_Hat = SteelZ/2 + Thickness_Steel_Add + Thickness_Hat/2;
  // Additional Steel
  G4double delta_Z_Steel_Add = SteelZ/2 + Thickness_Steel_Add/2;
  // LAB&PPO inside tube
  G4double delta_Z_Sct_Inside = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Length_In/2 + Thickness_Ring/2;
  // Air ring
  G4double delta_Z_upper_ring = SctZ/2 - (Length_WOM - Thickness_Steel_Add - Thickness_Hat - WallThick) + Thickness_Ring + Thickness_Gap/2;

  G4int n_sipm = 40;
  G4double radius_sipm = (Diam_WOM_In + Diam_WOM_Out)/4.;
  G4int sipm_id = 0;

  for(unsigned int pos = 0; pos < WOM_coord_vec.size(); pos++) {
    RM1 = new G4RotationMatrix();
    Outer_tube_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Outer_tube), Outer_tube_log, "Outer_tubePV", expHall_log, false, 401, intersect_check));
    Air_gap_out_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gaps), Air_gap_out_log, "Air_gap_outPV", expHall_log, false, 501, intersect_check));
    WLS_tube_out_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), WLS_tube_out_log, "WLS_outPV", expHall_log, false, 601, intersect_check));
    WOM_tube_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), WOM_tube_log, "WOM_tubePV", expHall_log, false, 701, intersect_check));
    WLS_tube_in_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_WOM), WLS_tube_in_log, "WLS_inPV", expHall_log, false, 602, intersect_check));
    Air_gap_in_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Air_gaps), Air_gap_in_log, "Air_gap_inPV", expHall_log, false, 502, intersect_check));
    Inner_tube_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Inner_tube), Inner_tube_log, "Inner_tubePV", expHall_log, false, 402, intersect_check));
    PMMA_Ring_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Ring), PMMA_Ring_log, "PMMA_RingPV", expHall_log, false, 403, intersect_check));
    Air_ring_out_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_upper_ring), Air_ring_out_log, "Air_ring_outPV", expHall_log, false, 503, intersect_check));
    PMMA_ring_lower_phys_vec.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_upper_ring), PMMA_ring_lower_log, "PMMA_ring_lowerPV", expHall_log, false, 404, intersect_check));
    Air_ring_in_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_upper_ring), Air_ring_in_log, "Air_ring_inPV", expHall_log, false, 504, intersect_check));
    PMMA_Disk_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Disk), PMMA_disk_log, "PMMA_DiskPV", expHall_log, false, 405, intersect_check));
    PMMA_Hat_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_PMMA_Hat), PMMA_Hat_log, "PMMA_HatPV", expHall_log, false, 406, intersect_check));
    Steel_Add_phys_vec.push_back( new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Steel_Add), Steel_Add_log, "Steel_AddPV", expHall_log, false, 103, intersect_check));
    Sct_Inside_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(WOM_coord_vec[pos].first, WOM_coord_vec[pos].second, delta_Z_Sct_Inside), Sct_Inside_log, "Sct_InsidePV", expHall_log, false, 205, intersect_check));
    if(pos == 0) sipm_id = 1000;
    if(pos == 2) sipm_id = 2000;
    if(pos == 4) sipm_id = 3000;
    if(pos == 6) sipm_id = 4000;
  for(int i = 0; i < n_sipm; i++) {
    RM1 = new G4RotationMatrix();
    RM1->rotateZ(-(i+0.5)*360./n_sipm*deg);
    G4double Xrotation = WOM_coord_vec[pos].first + radius_sipm*std::cos((i + 0.5)*2*pi/n_sipm);
    G4double Yrotation = WOM_coord_vec[pos].second + radius_sipm*std::sin((i + 0.5)*2*pi/n_sipm);
    sipm_phys_vec.push_back(new G4PVPlacement(RM1, G4ThreeVector(Xrotation, Yrotation, sipmSensThicknessTop + sipmSensThickness/2 + delta_Z_sipm), sipmSens_log, "sipmSensPV", expHall_log, false, sipm_id++, intersect_check));
    sipmSensTop_phys = new G4PVPlacement(RM1, G4ThreeVector(Xrotation, Yrotation, sipmSensThicknessTop/2 + delta_Z_sipm), sipmSensTop_log, "sipmSensTopPV", expHall_log, false, 802, intersect_check);
    sipmWindow_phys = new G4PVPlacement(RM1, G4ThreeVector(Xrotation, Yrotation, delta_Z_sipm), sipmWindow_log, "sipmWindowPV", expHall_log, false, 803, intersect_check);
    sipmBase_phys = new G4PVPlacement(RM1, G4ThreeVector(Xrotation, Yrotation, sipmWindowThickness + sipmBaseThickness/2. + delta_Z_sipm), sipmBaseBox_log, "sipmBasePV", expHall_log, false, 801, intersect_check);
    } 
  }

}

void OpNoviceDetectorConstruction::DefineVisAttributes()
{
  blue        = G4Color(0., 0., 1.);
  grey        = G4Color(0.3, 0.3, 0.3, 0.2);
  blue_trans  = G4Color(0., 0., 1., 0.1);
  green       = G4Color(0., 1., 0., 0.2);
  red         = G4Color(1., 0., 0., 0.2);
  white_trans = G4Color(1., 1., 1., 0.25);
  cyan        = G4Color(0., 1., 1., 0.3);
  magenta     = G4Color(1.,0.,1., 0.3);
  yellow     = G4Color(1.,1.,.8, 0.3);

  G4VisAttributes *worldVisAtt = new G4VisAttributes;
  worldVisAtt->SetVisibility(false);
  expHall_log->SetVisAttributes(worldVisAtt);
  
  G4VisAttributes *steelBoxVisAtt = new G4VisAttributes;
  steelBoxVisAtt->SetVisibility(true);
  steelBoxVisAtt->SetColor(white_trans);
  SteelBox_log->SetVisAttributes(steelBoxVisAtt);
  SteelBeam_log->SetVisAttributes(steelBoxVisAtt);
  G4VisAttributes *steelBoxVisAtt1 = new G4VisAttributes;
  steelBoxVisAtt1->SetVisibility(true);
  steelBoxVisAtt1->SetColor(grey);
  Steel_Add_log->SetVisAttributes(steelBoxVisAtt1);

  G4VisAttributes *sctBoxVisAtt = new G4VisAttributes;
  sctBoxVisAtt->SetColor(blue_trans);
  sctBoxVisAtt->SetVisibility(true);
  ScintillatorBoxLT_log->SetVisAttributes(sctBoxVisAtt);
  ScintillatorBoxRT_log->SetVisAttributes(sctBoxVisAtt);
  ScintillatorBoxLB_log->SetVisAttributes(sctBoxVisAtt);
  ScintillatorBoxRB_log->SetVisAttributes(sctBoxVisAtt);
  ScintillatorInBeams_log->SetVisAttributes(sctBoxVisAtt);
  G4VisAttributes *sctBoxVisAtt1 = new G4VisAttributes;
  sctBoxVisAtt1->SetColor(blue_trans);
  sctBoxVisAtt1->SetVisibility(true);
  Sct_Inside_log->SetVisAttributes(sctBoxVisAtt1);
  
  G4VisAttributes *reflectBoxVisAtt = new G4VisAttributes;
  reflectBoxVisAtt->SetColor(cyan);
  reflectBoxVisAtt->SetVisibility(false);
  ReflectBox_log->SetVisAttributes(reflectBoxVisAtt);

  G4VisAttributes *PMMAVisAtt = new G4VisAttributes;
  PMMAVisAtt->SetVisibility(true);
  PMMAVisAtt->SetColor(yellow);
  PMMA_disk_log->SetVisAttributes(PMMAVisAtt);
  PMMA_Ring_log->SetVisAttributes(PMMAVisAtt);
  PMMA_Hat_log->SetVisAttributes(PMMAVisAtt);
  Outer_tube_log->SetVisAttributes(PMMAVisAtt);
  Inner_tube_log->SetVisAttributes(PMMAVisAtt);
  PMMA_ring_lower_log->SetVisAttributes(PMMAVisAtt);

  G4VisAttributes *airVisAtt = new G4VisAttributes;
  airVisAtt->SetColor(green);
  airVisAtt->SetVisibility(false);
  Air_gap_out_log->SetVisAttributes(airVisAtt);
  Air_gap_in_log->SetVisAttributes(airVisAtt);
  Air_ring_out_log->SetVisAttributes(airVisAtt);
  Air_ring_in_log->SetVisAttributes(airVisAtt);

  G4VisAttributes *WLSVisAtt = new G4VisAttributes;
  WLSVisAtt->SetColor(red);
  WLSVisAtt->SetVisibility(false);
  WLS_tube_out_log->SetVisAttributes(WLSVisAtt);
  WLS_tube_in_log->SetVisAttributes(WLSVisAtt);

  G4VisAttributes *WOMVisAtt = new G4VisAttributes;
  WOMVisAtt->SetColor(magenta);
  WOMVisAtt->SetVisibility(true);
  WOM_tube_log->SetVisAttributes(WOMVisAtt);

  G4VisAttributes *sipmVisAtt1 = new G4VisAttributes;
  sipmVisAtt1->SetColor(blue);
  sipmVisAtt1->SetVisibility(true);
  sipmWindow_log->SetVisAttributes(sipmVisAtt1);
  G4VisAttributes *sipmVisAtt2 = new G4VisAttributes;
  sipmVisAtt2->SetColor(blue);
  sipmVisAtt2->SetVisibility(true);
  sipmBaseBox_log->SetVisAttributes(sipmVisAtt2);
  G4VisAttributes *sipmVisAtt3 = new G4VisAttributes;
  sipmVisAtt3->SetColor(blue);
  sipmVisAtt3->SetVisibility(true);
  sipmSens_log->SetVisAttributes(sipmVisAtt3);
  G4VisAttributes *sipmVisAtt4 = new G4VisAttributes;
  sipmVisAtt4->SetColor(blue);
  sipmVisAtt4->SetVisibility(true);
  sipmSensTop_log->SetVisAttributes(sipmVisAtt4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* OpNoviceDetectorConstruction::Construct()
{
  intersect_check = true; // global intersection check
  DefineMaterials();
  DefineMPTs();
  DefineSolids();
  DefineLogicalVolumes();
  ConstructVolumes();
  DefineSurfaces();
  DefineVisAttributes();
  return expHall_phys;
}
