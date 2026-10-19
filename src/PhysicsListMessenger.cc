// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Set variable parameters for Physics List
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

//--------------------------------------------------------------------------
PhysicsListMessenger::
  PhysicsListMessenger(PhysicsList* pPhys) 
  : G4UImessenger(),
    fPhysicsList(pPhys)
{
  fMGDir = new G4UIdirectory("/MG/");
  fMGDir->SetGuidance("UI commands of this example");

  fPhysDir = new G4UIdirectory("/MG/phys/");
  fPhysDir->SetGuidance("Physics List control");
 
  fVerboseCmd = new G4UIcmdWithAnInteger("/MG/phys/verbose",this);
  fVerboseCmd->SetGuidance("set verbose for physics processes");
  fVerboseCmd->SetParameterName("verbose",true);
  fVerboseCmd->SetDefaultValue(1);
  fVerboseCmd->SetRange("verbose>=0");
  fVerboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fCerenkovCmd =
           new G4UIcmdWithAnInteger("/MG/phys/cerenkovMaxPhotons",this);
  fCerenkovCmd->SetGuidance("set max # of photons per step");
  fCerenkovCmd->SetParameterName("MaxNumber",false);
  fCerenkovCmd->SetRange("MaxNumber>=0");
  fCerenkovCmd->AvailableForStates(G4State_PreInit);
}

//--------------------------------------------------------------------------
PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fVerboseCmd;
  delete fCerenkovCmd;
  delete fPhysDir;
  delete fMGDir;
}

//--------------------------------------------------------------------------
void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{
  if( command == fVerboseCmd ){
    fPhysicsList->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if( command == fCerenkovCmd ){
    fPhysicsList->
              SetNbOfPhotonsCerenkov(fCerenkovCmd->GetNewIntValue(newValue));
  }
}
