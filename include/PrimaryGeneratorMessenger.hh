// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify primary event parameters
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction* );
    virtual ~PrimaryGeneratorMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:
  PrimaryGeneratorAction* fPGAction;
  G4UIdirectory*          fGunDir;
  G4UIcmdWithAnInteger*   fModeCmd;
  G4UIcmdWithAString*     fFNameCmd;
  G4UIcmdWithAString*     fHNameCmd;
  G4UIcmdWithoutParameter* fPolarCmd;
};

#endif
