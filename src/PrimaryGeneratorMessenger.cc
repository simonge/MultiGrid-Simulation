// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify primary event parameters
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "Randomize.hh"

//--------------------------------------------------------------------------
PrimaryGeneratorMessenger::
  PrimaryGeneratorMessenger(PrimaryGeneratorAction* pga)
  : G4UImessenger(),
    fPGAction(pga)
{
  fGunDir = new G4UIdirectory("/MG/gun/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fModeCmd = new G4UIcmdWithAnInteger("/MG/gun/Mode",this);
  fModeCmd->SetGuidance("Set generator mode");
  fModeCmd->SetGuidance("0:GPS, 1:Alpha scan");
  fModeCmd->SetParameterName("mode",true);
  fModeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fFNameCmd = new G4UIcmdWithAString("/MG/gun/File-name",this);
  fFNameCmd->SetGuidance("Set root file name");
  fFNameCmd->SetParameterName("fname",true);
  fFNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fHNameCmd = new G4UIcmdWithAString("/MG/gun/Hist-name",this);
  fHNameCmd->SetGuidance("Set root histogram name");
  fHNameCmd->SetParameterName("hname",true);
  fHNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPolarCmd = new G4UIcmdWithoutParameter("/MG/gun/PhotonPolar",this);
  fPolarCmd->SetGuidance("Set linear polarization");
  fPolarCmd->AvailableForStates(G4State_Idle);
}

//--------------------------------------------------------------------------
PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fModeCmd;
  delete fGunDir;
}

//--------------------------------------------------------------------------
void PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{
  if( command == fModeCmd ) {
    G4int mode = fModeCmd->GetNewIntValue(newValue);
    fPGAction->SetMode(mode);
  }
  else if( command == fFNameCmd ) {
    fPGAction->SetFName(newValue);
  }
  else if( command == fHNameCmd ) {
    fPGAction->SetHName(newValue);
  }
  if( command == fPolarCmd ) {
    fPGAction->SetOptPhotonPolar();
  }
}
