// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Record hits in sensitive detector
// Discard G4 hits collection stuff
// JRM Annand 20/08/2018      
//--------------------------------------------------------------------------
#ifndef DetectorSD_h
#define DetectorSD_h 1

#include "DetectorHit.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class G4HCofThisEvent;
class G4Step;

//---------------------------------------------------------------------------

class DetectorSD : public G4VSensitiveDetector
{
public:
  
  DetectorSD(G4String,G4int,G4int,G4int,G4int=0);
  ~DetectorSD();
  
  void   Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  void   EndOfEvent(G4HCofThisEvent*);
  G4int* GetHits(){ return fHits; }
  //G4int* GetHitID(){ return fHitID; }
  G4int  GetNhits(){ return fNhits; }
  DetectorHit* GetDetHits(){ return fDetectorHits; }
  void ZeroHits(){
    for (G4int i=0;i<fNhits;i++) {
      //fHitID[fHits[i]] = -1;
      fHits[i]         = 0;
    }
    fNhits = 0;
  }
private:
  DetectorHit*             fDetectorHits;
  G4int                    fNmax;
  G4int                    fNcell;
  G4int                    fNvess;
  G4int                    fNbox;
  G4int                    fNhits;
  //G4int*                   fHitID;
  G4int*                   fHits;
  //G4double*                fPixd;
};
#endif

//---------------------------------------------------------------------------


