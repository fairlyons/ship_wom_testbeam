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
#include "G4AnalysisManager.hh"
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
  //analysisManager->SetFileName("LScin");  

//--------------------------------------------------------------------
  analysisManager -> CreateNtuple("Detected","Results");
  analysisManager -> CreateNtupleDColumn("x");
  analysisManager -> CreateNtupleDColumn("y");
  //analysisManager -> CreateNtupleIColumn("process");
  //analysisManager -> CreateNtupleIColumn("WomNo");
  //analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleDColumn("time");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> CreateNtupleIColumn("sipmNumber");
  analysisManager -> CreateNtupleIColumn("WOMNumber");
  analysisManager -> CreateH3("evt_quadrant_time","",10,0,10,64,0,64,1024,0,300);
  analysisManager -> FinishNtuple(0);
//--------------------------------------------------------------------

//--------------------------------------------------------------------
/*  analysisManager -> CreateNtuple("EventStat","Results");
  analysisManager -> CreateNtupleIColumn("scintillation_photons");
  analysisManager -> CreateNtupleIColumn("Cerenkov_photons");
  analysisManager -> CreateNtupleDColumn("Edep_other");
  analysisManager -> CreateNtupleDColumn("Edep_scintillator");
  analysisManager -> CreateNtupleDColumn("Edep_walls");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> FinishNtuple(1);*/
//--------------------------------------------------------------------
/*
  analysisManager -> CreateNtuple("entersWOM","Results");
  analysisManager -> CreateNtupleIColumn("first");
  analysisManager -> CreateNtupleIColumn("parentID");
  analysisManager -> CreateNtupleIColumn("process");
  analysisManager -> CreateNtupleIColumn("PreCopyNo");
  analysisManager -> CreateNtupleIColumn("PostCopyNo");
  analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleIColumn("stepNumber");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> FinishNtuple(2);

  analysisManager -> CreateNtuple("entersPMMA","Results");
  analysisManager -> CreateNtupleIColumn("first");
  analysisManager -> CreateNtupleIColumn("parentID");
  analysisManager -> CreateNtupleIColumn("process");
  analysisManager -> CreateNtupleIColumn("PreCopyNo");
  analysisManager -> CreateNtupleIColumn("PostCopyNo");
  analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleIColumn("stepNumber");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> FinishNtuple(3);

  analysisManager -> CreateNtuple("bornWLS","Results");
  analysisManager -> CreateNtupleIColumn("first");
  analysisManager -> CreateNtupleIColumn("parentID");
  analysisManager -> CreateNtupleIColumn("process");
  analysisManager -> CreateNtupleIColumn("PreCopyNo");
  analysisManager -> CreateNtupleIColumn("PostCopyNo");
  analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleIColumn("stepNumber");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> FinishNtuple(4);

  analysisManager -> CreateNtuple("absorbedWLS","Results");
  analysisManager -> CreateNtupleIColumn("first");
  analysisManager -> CreateNtupleIColumn("parentID");
  analysisManager -> CreateNtupleIColumn("process");
  analysisManager -> CreateNtupleIColumn("PreCopyNo");
  analysisManager -> CreateNtupleIColumn("PostCopyNo");
  analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleIColumn("eventNumber");
  analysisManager -> CreateNtupleIColumn("stepNumber");
  analysisManager -> FinishNtuple(5);

  analysisManager -> CreateNtuple("Dies","Results");
  analysisManager -> CreateNtupleIColumn("CopyNo");
  analysisManager -> CreateNtupleDColumn("waveLen");
  analysisManager -> CreateNtupleDColumn("x");
  analysisManager -> CreateNtupleDColumn("y");
  analysisManager -> FinishNtuple(6); 
*/
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
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << " " << *fTimer << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
