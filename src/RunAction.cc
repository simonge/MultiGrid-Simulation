// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify run action...introduce timer
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "G4Timer.hh"
#include "RunAction.hh"
#include "G4Run.hh"

//--------------------------------------------------------------------------
RunAction::RunAction()
 : G4UserRunAction(),
   fTimer(0)
{
  fTimer = new G4Timer;
}

//--------------------------------------------------------------------------
RunAction::~RunAction()
{
  delete fTimer;
}

//--------------------------------------------------------------------------
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
}

//--------------------------------------------------------------------------
void RunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
}
