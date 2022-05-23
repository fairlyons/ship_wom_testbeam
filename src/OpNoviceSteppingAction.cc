//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
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
// $Id: OpNoviceSteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file OpNoviceSteppingAction.cc
/// \brief Implementation of the OpNoviceSteppingAction class

#include "OpNoviceSteppingAction.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Geantino.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include <cmath>
#include "g4root.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4LogicalVolume.hh"
#include "OpNoviceDetectorConstruction.hh"
#include "G4VSolid.hh"
#include <vector>
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "OpNoviceEventAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceSteppingAction::OpNoviceSteppingAction(OpNoviceEventAction* EvAct):G4UserSteppingAction(),fEventAction(EvAct)
{
  fEventNumber = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceSteppingAction::~OpNoviceSteppingAction()
{}


G4int sipm_detection(G4double wl)
{
  std::vector<double> wl_vec = {283.32062824769,289.582356019053,291.851293147889,291.715757914742,291.542247040993,291.147787961894,293.34486564184,293.282467883949,295.507643349451,294.293612273865,296.340897745862
  ,297.015056654157,307.013465028043,298.568341976658,300.113443634526,300.967315603349,302.9596835306,311.78302720843,315.340827078522,321.744866844688,329.572026558891,338.110746247113,
  347.598212567359,361.592225389723,372.867478761135,377.958104792148,384.362144558314,389.705046550231,399.304904012702,408.397059236272,420.289721059469,431.325102843534,446.97942227194,
  461.479007289261,476.864941180716,485.737129093121,493.527720069697,502.481100245381,508.986791436407,515.101986504042,532.408950163217,519.558539621825,526.674139362009,539.482218894342,
  544.463138712471,554.42497834873,562.963698036952,571.502417725173,580.041137413395,588.579857101617,597.183264060204,612.363050704067,602.018304344688,620.600055932448,630.561895568707,
  640.523735204965,650.485574841224,660.447414477483,671.645180788106,680.37109375,690.332933386259,704.822881948089,720.218452295035,735.872771723441,751.527091151847,767.181410580254,
  782.83573000866,798.60864276607,810.133758102824,822.683088553695,838.337407982101,853.991727410508,869.646046838914,885.300366267321};
  std::vector<double> eff_vec = {0.026526786824411,0.042202926083421,0.099461599466405,0.077280017273974,0.060746384220011,0.052571693483039,0.134999961095994,0.118459329017847,0.170466786513831,0.153436957527306,0.198733559050505,
  0.185564020054936,0.300076779492663,0.221644583298286,0.244494495093093,0.262674077462957,0.27963118570716,0.318640406447815,0.3324116751661,0.347682245368915,0.361821294310075,0.375023432468283,
  0.388101669952467,0.39946901483621,0.411514245710859,0.42670217815887,0.442491341724927,0.453393011986188,0.470229259321717,0.481924633423035,0.494084426931976,0.501906240893377,0.509188412034585, 0.507662157652688,0.497041799436511,0.485277195273939,0.472187986094495,0.456045925166636,0.4445350205978,0.432226649587411,0.401911114819206,0.424042352885693,0.411459999903534,0.391695977817467,
  0.380843018476912,0.365951864733511,0.352614224890481,0.339378270020637,0.326549055043531,0.3136036286685,0.299313719619317,0.281472503125024,0.293666519075303,0.267462403056718,0.253497308890536,
  0.239706531821242,0.226090071848837,0.213332746139666,0.201155690203397,0.191486486991037,0.180237081246375,0.166175861979421,0.153000281895251,0.139768724905833,0.127424600409666,0.116023372937575,
  0.105970047833846,0.096379388185778,0.088385388686368,0.079750594524847,0.070789806409192,0.061995411886021,0.053589269078647,0.045737771579554}; // from the datasheet S14160-3050HS
  G4int veclen = wl_vec.size();
  G4int ind = 0;
  while(wl>wl_vec[ind]) ind++;
  if( (ind==veclen)||(ind==0) ) return 0;
  G4double rnd = double(rand())/RAND_MAX;

  if( rnd < eff_vec[ind] ) return 1;
  else return 0;
}

void OpNoviceSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) fEventNumber = eventNumber;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4String ParticleName = aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
      

  if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    G4Track* track = aStep->GetTrack();

    G4int stepnum = track->GetCurrentStepNumber();
    G4int parentid = track->GetParentID();

    G4int trackid = track->GetTrackID();

    G4int process = 0;
    std::string processname = "";

    if(track->GetCreatorProcess()) std::string processname = track->GetCreatorProcess()->GetProcessName();
    if(track->GetCreatorProcess()) processname = track->GetCreatorProcess()->GetProcessName();
    if(processname == "Scintillation") process = 1;
    else if( processname == "OpWLS") process = 3;
    else if( processname == "Cerenkov") process = 2;

    G4VPhysicalVolume* prevolumephys = NULL;
    prevolumephys = aStep->GetPreStepPoint()->GetPhysicalVolume();
    G4String prephysvolname = prevolumephys->GetName();
    G4int pre_copynum = prevolumephys->GetCopyNo();
    G4VPhysicalVolume* postvolumephys = NULL;
    postvolumephys = aStep->GetPostStepPoint()->GetPhysicalVolume();
    if(track->GetTrackStatus() == fStopAndKill) {
      analysisManager->FillNtupleIColumn(6,0, pre_copynum);
      analysisManager->FillNtupleDColumn(6,1, 1.24e-3 / track->GetKineticEnergy());
      analysisManager->AddNtupleRow(6);
    }
    if(!postvolumephys) return;
    G4String postphysvolname = postvolumephys->GetName();
    G4int post_copynum = postvolumephys->GetCopyNo();
    const G4VTouchable* posttouchable = NULL;
    posttouchable = aStep->GetPostStepPoint()->GetTouchable();
    if(!posttouchable) return;

    if(stepnum == 1) {
      if(process == 1) fEventAction->scintillation_photons++;
      if(process == 2) fEventAction->cherenkov_photons++;
    }

    //---------------- 2. Born in WLS
    if(process == 3) {
      PhotonInfo info = {parentid, process, pre_copynum, post_copynum, 1.24e-3 / track -> GetKineticEnergy(), stepnum};
      fEventAction->map_bornWLS[trackid] = info;
    }
    //---------------- 2. Born in WLS END

    //---------------- 4. photons that reached SiPM's
    if(postphysvolname == "sipmBase") {
      G4cout << "sipm hit !!!!!!!!! " << posttouchable->GetCopyNumber(1) << G4endl;
      analysisManager->FillNtupleDColumn(0,0, aStep-> GetPostStepPoint()-> GetPosition().getX());
      analysisManager->FillNtupleDColumn(0,1, aStep-> GetPostStepPoint()-> GetPosition().getY());
      analysisManager->FillNtupleIColumn(0,2, process);
      if(posttouchable->GetCopyNumber(1) < 41) analysisManager->FillNtupleIColumn(0,3,1); // WOM number
      if(posttouchable->GetCopyNumber(1) > 40) analysisManager->FillNtupleIColumn(0,3,2); // WOM number
      analysisManager->FillNtupleDColumn(0,4, 1.24e-3 / track->GetKineticEnergy());
      analysisManager->FillNtupleDColumn(0,5, track -> GetGlobalTime());
      analysisManager->FillNtupleIColumn(0,6, sipm_detection(1.24e-3 / track->GetKineticEnergy()));
      analysisManager->FillNtupleIColumn(0,7, eventNumber);
      analysisManager->FillNtupleIColumn(0,8, posttouchable->GetCopyNumber(1)); //sipm number
      analysisManager->AddNtupleRow(0);
      track->SetTrackStatus(fStopAndKill);
    }
    //---------------- 4. photons that reached SiPM's END

    //---------------- 3. Absorbed in WLS
    if((postphysvolname == "WLS1") || (postphysvolname == "WLS2") || (postphysvolname == "WLSring")) {
      PhotonInfo info = {parentid, process, pre_copynum, post_copynum, 1.24e-3 / track -> GetKineticEnergy(), stepnum};
      fEventAction->map_absorbedWLS[trackid] = true;
      fEventAction->map_absorbedWLS_info[trackid] = info;
    }
    else if(fEventAction->map_absorbedWLS.find(trackid) != fEventAction->map_absorbedWLS.end()) fEventAction->map_absorbedWLS[trackid] = false;
    //---------------- 3. Absorbed in WLS END 

    //---------------- 1. Fallen on WOM
    if((postphysvolname == "WOM_tube")) {
      PhotonInfo info = {parentid, process, pre_copynum, post_copynum, 1.24e-3 / track -> GetKineticEnergy(), stepnum};
      fEventAction->map_entersWOM[trackid] = info;
    }
    //---------------- 1. Fallen on WOM END 

    //---------------- 5. Fallen on PMMA vessel from scintillator
    if((postphysvolname == "Outer_tube") || (postphysvolname == "Inner_tube") || (postphysvolname == "PMMA_Ring")  || (postphysvolname == "PMMA_Disk") || (postphysvolname == "PMMA_ring_lower")) {
      if((prephysvolname == "ScintillatorBoxPV") || (prephysvolname == "SteelBox") || (prephysvolname == "Sct_Inside") || (prephysvolname == "ReflectBox")) {
        PhotonInfo info = {parentid, process, pre_copynum, post_copynum, 1.24e-3 / track -> GetKineticEnergy(), stepnum};
        fEventAction->map_entersPMMAvessel[trackid] = info;
      }
    }
    //---------------- 5. Fallen on PMMA vessel from scintillator END
  }

  else {
    // non-optical particles
    G4int volume_index;
    G4VPhysicalVolume* prevolumephys = NULL;
    prevolumephys = aStep->GetPreStepPoint()->GetPhysicalVolume();
    G4String prephysvolname = prevolumephys->GetName();
    if(prephysvolname == "ScintillatorBoxPV") volume_index = 1;
    else if(prephysvolname == "SteelBox") volume_index = 2;
    else volume_index = 0;
    fEventAction->addEdep(volume_index, aStep -> GetTotalEnergyDeposit());
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
