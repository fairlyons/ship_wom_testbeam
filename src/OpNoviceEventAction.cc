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
/// \file electromagnetic/TestEm4/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 75839 2013-11-06 17:27:26Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "OpNoviceEventAction.hh"
#include "G4RunManager.hh"
#include "g4root.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceEventAction::OpNoviceEventAction()
:G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceEventAction::BeginOfEventAction(const G4Event*)
{ 
  G4cout << "\n---> Begin of Event: " << G4RunManager::GetRunManager()->
  GetCurrentEvent()->GetEventID() << G4endl;

  scintillation_photons = 0;
  cherenkov_photons = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceEventAction::EndOfEventAction(const G4Event*)
{                           
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  /*
  std::cout<<"fallen_on_steel= "<<fallen_on_steel<<"\n";
  std::cout<<"absorbed_by_steel= "<<absorbed_by_steel<<"\n";
  std::cout<<"absorbtion coef= "<<double(absorbed_by_steel)/fallen_on_steel<<"\n";

  map<G4int, PhotonInfo>::iterator it;
  for(it=map_entersWOM.begin(); it!=map_entersWOM.end(); it++) {
    analysisManager->FillNtupleIColumn(1,0, it->first);
    analysisManager->FillNtupleIColumn(1,1, it->second.parentID );
    analysisManager->FillNtupleIColumn(1,2, it->second.process);
    analysisManager->FillNtupleIColumn(1,3, it->second.WOWnumber );
    analysisManager->FillNtupleDColumn(1,4, it->second.waveLen );
    analysisManager->FillNtupleIColumn(1,5, eventNumber);
    analysisManager->AddNtupleRow(1);
  }
  for(it=map_entersPMMAvessel.begin(); it!=map_entersPMMAvessel.end(); it++) {
    analysisManager->FillNtupleIColumn(5,0, it->first);
    analysisManager->FillNtupleIColumn(5,1, it->second.parentID );
    analysisManager->FillNtupleIColumn(5,2, it->second.process);
    analysisManager->FillNtupleIColumn(5,3, it->second.WOWnumber );
    analysisManager->FillNtupleDColumn(5,4, it->second.waveLen );
    analysisManager->FillNtupleIColumn(5,5, eventNumber);
    analysisManager->AddNtupleRow(5);
  }
  for(it=map_bornWLS.begin(); it!=map_bornWLS.end(); it++) {
    analysisManager->FillNtupleIColumn(2,0, it->first);
    analysisManager->FillNtupleIColumn(2,1, it->second.parentID );
    analysisManager->FillNtupleIColumn(2,2, it->second.process);
    analysisManager->FillNtupleIColumn(2,3, it->second.WOWnumber );
    analysisManager->FillNtupleDColumn(2,4, it->second.waveLen );
    analysisManager->FillNtupleIColumn(2,5, eventNumber);
    analysisManager->AddNtupleRow(2);
  }
  map<G4int, G4bool>::iterator it2;

  for(it2=map_absorbedWLS.begin(); it2!=map_absorbedWLS.end(); it2++) {
    if(it2->second) {
      PhotonInfo info = map_absorbedWLS_info[it2->first];
      analysisManager->FillNtupleIColumn(3,0, it2->first);
      analysisManager->FillNtupleIColumn(3,1, info.parentID );
      analysisManager->FillNtupleIColumn(3,2, info.process);
      analysisManager->FillNtupleIColumn(3,3, info.WOWnumber );
      analysisManager->FillNtupleDColumn(3,4, info.waveLen );
      analysisManager->FillNtupleIColumn(3,5, eventNumber);
      analysisManager->AddNtupleRow(3);
    }
  }
  */

  map_bornWLS.clear();
  map_absorbedWLS.clear();
  map_absorbedWLS_info.clear();
  map_entersWOM.clear();
  map_entersPMMAvessel.clear();

  analysisManager->FillNtupleIColumn(1,0, scintillation_photons );
  analysisManager->FillNtupleIColumn(1,1, cherenkov_photons );
  analysisManager->FillNtupleDColumn(1,2, fEnergyDeposit[0] );
  analysisManager->FillNtupleDColumn(1,3, fEnergyDeposit[1] );
  analysisManager->FillNtupleDColumn(1,4, fEnergyDeposit[2] );
  analysisManager->FillNtupleIColumn(1,5, eventNumber);
  analysisManager->AddNtupleRow(1);

  fEnergyDeposit.clear();
}

OpNoviceEventAction::~OpNoviceEventAction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


