// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify run action...introduce timer
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef RunAction_h
#define RunAction_h 1
#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    G4Timer* fTimer;
};

#endif /*RunAction_h*/
