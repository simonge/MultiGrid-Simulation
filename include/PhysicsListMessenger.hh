// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Set variable parameters for Physics List
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
//--------------------------------------------------------------------------
class PhysicsListMessenger: public G4UImessenger
{
  public:
    PhysicsListMessenger(PhysicsList* );
    virtual ~PhysicsListMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:
  PhysicsList*  fPhysicsList;
  
  G4UIdirectory*        fMGDir;
  G4UIdirectory*        fPhysDir;
  G4UIcmdWithAnInteger* fVerboseCmd;
  G4UIcmdWithAnInteger* fCerenkovCmd;
};
//--------------------------------------------------------------------------
#endif
