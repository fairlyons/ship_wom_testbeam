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
// $Id: OpNoviceSteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file OpNoviceSteppingAction.cc
/// \brief Implementation of the OpNoviceSteppingAction class

#include "OpNoviceSteppingAction.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include <cmath>
#include "g4root.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
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
{ ; }


G4int sipm_detection(G4double wl)
{
  G4double filling_factor = 0.67;
  std::vector<double> wl_vec = {320., 330., 340., 350., 360., 370., 380., 390., 400., 410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550.};
  std::vector<double> eff_vec = {0.03, 0.09, 0.19, 0.27, 0.31, 0.34, 0.39, 0.43, 0.45, 0.46, 0.48, 0.49, 0.49, 0.50, 0.50, 0.50, 0.50, 0.49, 0.48, 0.47, 0.46, 0.44, 0.42, 0.41};
  G4int veclen = wl_vec.size();
  G4int ind = 0;
  while(wl>wl_vec[ind])
    ind++;
  if( (ind==veclen)||(ind==0) )
    return 0;
  G4double rnd = double(rand())/RAND_MAX;

  if( rnd < eff_vec[ind]*filling_factor )
  {
    return 1;
  }
  else return 0;
}

void OpNoviceSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4int eventNumber = G4RunManager::GetRunManager()->
                                              GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber)
  {
     fEventNumber = eventNumber;
  }
  

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4String ParticleName = aStep->GetTrack()->GetDynamicParticle()->
                                  GetParticleDefinition()->GetParticleName();
      

  if(aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    G4Track* track = aStep->GetTrack();

    G4int stepnum=track->GetCurrentStepNumber();
    G4int parentid=track -> GetParentID();


    G4int trackid = track -> GetTrackID();

    G4int process = 0;
    std::string processname = track->GetCreatorProcess()->GetProcessName();

    if(processname == "Scintillation")
      process = 1;
    else if( processname == "OpWLS") 
      process = 3;
    else if( processname == "Cerenkov")
      process = 2;
    
    G4VPhysicalVolume* postvolumephys = NULL;
    postvolumephys = aStep->GetPostStepPoint()->GetPhysicalVolume();
    if(!postvolumephys) return;
    G4String postphysvolname = postvolumephys->GetName();

    G4int post_copynum = postvolumephys->GetCopyNo();


    if(stepnum ==1)
    {
      if(process==1) fEventAction->scintillation_photons++;
      if(process==2) fEventAction->cherenkov_photons++;


      //---------------- 2. Born in WLS
      G4VPhysicalVolume* prevolumephys = NULL;
      prevolumephys = aStep->GetPreStepPoint()->GetPhysicalVolume();
      G4String prephysvolname = prevolumephys->GetName();
      
      
      
      if( (prephysvolname == "WLS1") || (prephysvolname == "WLS2") || (prephysvolname == "WLSring") )
      {
        PhotonInfo info = {parentid, process, post_copynum, 1.24e-3 / track -> GetKineticEnergy()};
        fEventAction->map_bornWLS[trackid] = info;
        // if(process == 0)
        // {
        //   std::cout<<"process: "<<processname<<"\n";
        // }
      }
      //---------------- 2. Born in WLS END
    }




    //---------------- 4. photons that reached SiPM's

    if(postphysvolname == "sipm_base")
    {
      analysisManager->FillNtupleDColumn(0,0, aStep -> GetPostStepPoint() -> GetPosition().getX() );
      analysisManager->FillNtupleDColumn(0,1, aStep -> GetPostStepPoint() -> GetPosition().getY() );
      analysisManager->FillNtupleIColumn(0,2, process );
      analysisManager->FillNtupleIColumn(0,3, post_copynum );
      // analysisManager->FillNtupleIColumn(0,2, parentid );
      // analysisManager->FillNtupleIColumn(0,3, trackid);
      analysisManager->FillNtupleDColumn(0,4, 1.24e-3 / track -> GetKineticEnergy());
      analysisManager->FillNtupleDColumn(0,5, track -> GetGlobalTime() );
      analysisManager->FillNtupleIColumn(0,6, sipm_detection(1.24e-3 / track -> GetKineticEnergy()));
      analysisManager->FillNtupleIColumn(0,7, eventNumber);
      analysisManager->AddNtupleRow(0);
      track->SetTrackStatus(fStopAndKill);
    }
    //---------------- 4. photons that reached SiPM's END

    //---------------- 3. Fallen on WOM + absorbed in WLS


    if( (postphysvolname == "WLS1") || (postphysvolname == "WLS2") || (postphysvolname == "WLSring") )
    {
      PhotonInfo info = {parentid, process, post_copynum, 1.24e-3 / track -> GetKineticEnergy()};
      if( fEventAction->map_entersWOM.find(trackid) == fEventAction->map_entersWOM.end() )
      {
        fEventAction->map_entersWOM[trackid] = info;
      }
      fEventAction->map_absorbedWLS[trackid] = true;
      fEventAction->map_absorbedWLS_info[trackid] = info;
    }
    else
    {
      if( fEventAction->map_absorbedWLS.find(trackid) != fEventAction->map_absorbedWLS.end() )
      {
        fEventAction->map_absorbedWLS[trackid] = false;
      }
    }


    //---------------- 3. Fallen on WOM + absorbed in WLS END 

    //---------------- 5. Fallen on PMMA vessel from scintillator
    if( (postphysvolname == "Outer_tube") || (postphysvolname == "Inner_tube") || (postphysvolname == "PMMA_Ring")  || (postphysvolname == "PMMA_Disk") )
    {
      G4VPhysicalVolume* prevolumephys = NULL;
      prevolumephys = aStep->GetPreStepPoint()->GetPhysicalVolume();
      G4String prephysvolname = prevolumephys->GetName();
      if( (prephysvolname == "ScintilatorBoxPV") || (prephysvolname == "SteelBox") || (prephysvolname == "Sct_Outside") )
      {
        PhotonInfo info = {parentid, process, post_copynum, 1.24e-3 / track -> GetKineticEnergy()};
        fEventAction->map_entersPMMAvessel[trackid] = info;
      }
    }
    //---------------- 5. Fallen on PMMA vessel from scintillator END
  }
  else
  {
    // non-optical particles
      G4int volume_index;
      G4VPhysicalVolume* prevolumephys = NULL;
      prevolumephys = aStep->GetPreStepPoint()->GetPhysicalVolume();
      G4String prephysvolname = prevolumephys->GetName();
      if(prephysvolname == "ScintilatorBoxPV")
      {
        volume_index =1;
      }
      else if(prephysvolname == "SteelBox")
      {
        volume_index =2;
      }
      else
      {
        volume_index =0;
      }

    fEventAction->addEdep(volume_index, aStep -> GetTotalEnergyDeposit());
  }

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
