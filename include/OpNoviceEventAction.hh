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
/// \file electromagnetic/TestEm4/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 75839 2013-11-06 17:27:26Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef OpNoviceEventAction_h
#define OpNoviceEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include <string>
#include <map>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;

typedef struct
{
    G4int parentID;
    G4int process;
    G4int WOWnumber;
    G4double waveLen;
} PhotonInfo;

class OpNoviceEventAction : public G4UserEventAction
{
  public:
    OpNoviceEventAction();
   ~OpNoviceEventAction();
  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    void addEdep(G4int volume_index, G4double Edep) {fEnergyDeposit[volume_index] += Edep;};
    //G4double GetEnergyDeposit()     {return fTotalEnergyDeposit;};

    map<G4int, PhotonInfo> map_bornWLS;
    map<G4int, G4bool> map_absorbedWLS;
    map<G4int, PhotonInfo> map_absorbedWLS_info;
    map<G4int, PhotonInfo> map_entersWOM;

    map<G4int, PhotonInfo> map_entersPMMAvessel;

    //long long int fallen_on_steel;
    G4int scintillation_photons;
    G4int cherenkov_photons;
    
  private:
    map<G4int, G4double>  fEnergyDeposit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
