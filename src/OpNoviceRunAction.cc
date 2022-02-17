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

// Make this appear first!
#include "G4Timer.hh"

#include "OpNoviceRunAction.hh"
#include "OpNoviceDetectorConstruction.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"
#include "Randomize.hh"
#include "G4HCofThisEvent.hh"

#include "OpNoviceEventAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceRunAction::OpNoviceRunAction()
 : G4UserRunAction(),
   fTimer(0)
{
  fTimer = new G4Timer;
  // automatic (time-based) random seeds for each run
     G4cout << "*******************" << G4endl;
     G4cout << "*** AUTOSEED ON ***" << G4endl;
     G4cout << "*******************" << G4endl;
     long seeds[1];
     time_t systime = time(NULL);
     //seeds[0] = (long) systime;
     G4Random::showEngineStatus();

     // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName("LScin");
     


//--------------------------------------------------------------------
  analysisManager -> CreateNtuple("Detected","Results");
  analysisManager -> CreateNtupleDColumn("x");
  analysisManager -> CreateNtupleDColumn("y");
  analysisManager -> CreateNtupleIColumn("process");
  analysisManager -> CreateNtupleIColumn("WOWnumber");
  analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleDColumn("time");
  analysisManager -> CreateNtupleIColumn("detection");
  // analysisManager -> CreateNtupleIColumn("Steel_reflections");
  // analysisManager -> CreateNtupleDColumn("Scint_totalpath");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> CreateNtupleIColumn("sipmNumber");

  analysisManager -> FinishNtuple(0);
//--------------------------------------------------------------------

//--------------------------------------------------------------------
  analysisManager -> CreateNtuple("EventStat","Results");
  analysisManager -> CreateNtupleIColumn("scintillation_photons");
  analysisManager -> CreateNtupleIColumn("Cerenkov_photons");
  analysisManager -> CreateNtupleDColumn("Edep_other");
  analysisManager -> CreateNtupleDColumn("Edep_scintillator");
  analysisManager -> CreateNtupleDColumn("Edep_walls");
  analysisManager -> CreateNtupleIColumn("eventNumber");

  analysisManager -> FinishNtuple(1);
//--------------------------------------------------------------------


//--------------------------------------------------------------------
  // analysisManager -> CreateNtuple("PotentiallyDetected","Results");
  // analysisManager -> CreateNtupleIColumn("photonID");
  // analysisManager -> CreateNtupleIColumn("parentID");
  // analysisManager -> CreateNtupleIColumn("process");
  // analysisManager -> CreateNtupleDColumn("waveLen");
  // analysisManager -> CreateNtupleIColumn("WOM_1");
  // analysisManager -> CreateNtupleIColumn("WOM_2");
  // analysisManager -> CreateNtupleIColumn("WOM_3");
  // analysisManager -> CreateNtupleIColumn("WOM_4");
  // // "WOM_*" variable content:
  // // +10 if a photon fell on corresponding PMMA vessel
  // // +100 if a photon fell on corresponding WOM
  // // +1000 if a photon was absorbed in WLS
  // // +10000 if a photon was born in WLS
  // analysisManager -> CreateNtupleIColumn("eventNumber");
  // analysisManager -> FinishNtuple(1);
//--------------------------------------------------------------------

  // analysisManager -> CreateNtuple("MassiveParticles","Results");
  // analysisManager -> CreateNtupleIColumn("trackID");
  // analysisManager -> CreateNtupleIColumn("parentID");
  // analysisManager -> CreateNtupleIColumn("charge");
  // analysisManager -> CreateNtupleDColumn("postX");
  // analysisManager -> CreateNtupleDColumn("postY");
  // analysisManager -> CreateNtupleDColumn("postZ");
  // analysisManager -> CreateNtupleDColumn("totEn");
  // analysisManager -> CreateNtupleDColumn("Edep");
  // analysisManager -> CreateNtupleDColumn("time");
  // analysisManager -> CreateNtupleIColumn("stepnum");
  // analysisManager -> CreateNtupleIColumn("volume_index");
  // analysisManager -> CreateNtupleIColumn("eventNumber");
  // analysisManager -> FinishNtuple(4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceRunAction::~OpNoviceRunAction()
{
  delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


    // Open an output file
    // The default file name is set in RunAction::RunAction(),
    // it can be overwritten in a macro
    analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
