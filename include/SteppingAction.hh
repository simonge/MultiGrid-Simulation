// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify stepping
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef SteppingAction_H
#define SteppingACtion_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

//---------------------------------------------------------------------------

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  ~SteppingAction();
  virtual void UserSteppingAction(const G4Step*);
};

#endif

