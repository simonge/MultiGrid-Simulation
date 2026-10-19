// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify stepping
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "SteppingAction.hh"
#include "DetectorSD.hh"
#include "DetectorConstruction.hh"//added this one

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//---------------------------------------------------------------------------

SteppingAction::SteppingAction()
{
}

//---------------------------------------------------------------------------

SteppingAction::~SteppingAction()
{
}

//---------------------------------------------------------------------------

void SteppingAction::UserSteppingAction(const G4Step * aStep)
{
  G4Track* track = aStep->GetTrack();
  // abort track if over max time
  if(track->GetGlobalTime() > 0.1e9*CLHEP::ns)
    track->SetTrackStatus(fStopAndKill);
}













