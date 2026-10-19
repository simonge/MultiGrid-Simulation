// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Write hits to ROOT file at end of event
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

//---------------------------------------------------------------------------

class G4Event;
class AnalysisManager;
class PrimaryGeneratorAction;
class DetectorConstruction;
class DetectorHit;

//---------------------------------------------------------------------------

class EventAction : public G4UserEventAction
{
  public:
    EventAction(AnalysisManager*);
   ~EventAction();

  public:
    void     BeginOfEventAction(const G4Event*);
    void     EndOfEventAction(const G4Event*);

  private:
  AnalysisManager* fAM;
  PrimaryGeneratorAction* fPGA;
  DetectorConstruction* fDC;
  DetectorHit* fDetHit;
  G4int* fHits;
};
#endif

//---------------------------------------------------------------------------

    
