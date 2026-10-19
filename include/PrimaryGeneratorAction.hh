// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Modify primary event parameters
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh" 
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class AnalysisManager;
class DetectorConstruction;
class G4GeneralParticleSource;
class TFile;
class TTree;

enum { EPGA_GPS, EPGA_Chop, EPGA_Alpha};
const G4double hp0 = 6.6260755e-34;   // h in J/s
const G4double mn0 = 1.50535074e-10;  // m neutron in J
const G4double mn1 = 939.56563;       // in MeV
const G4double c0 = 299792458;        // light vel m/s
const G4double csq = c0*c0;

//--------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
  void GeneratePrimaries (G4Event*);
  void GenChoppedBeam(G4Event*);
  void GenAlpha(G4Event*);
  void SetFName(G4String);
  void SetHName(G4String);
  void SetMode(G4int);
  void GetDCinfo();
  G4int GetNEvents(); 
  G4int GetMode()                     { return fGenMode; }
  G4ThreeVector GetVertex()           { return G4ThreeVector(fVx,  fVy,  fVz);  }
  G4ThreeVector GetDirection()        { return G4ThreeVector(fPxp, fPyp, fPzp); }
  G4double      GetEnergy()           { return fEp; }
  G4double      GetTime()             { return fTp; }
  G4ParticleDefinition* GetPrimPDef() { return fPDef; }
  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);
  void SetAnalysisManager(AnalysisManager* am){ fAM = am; }
  void SetDetectorConstruction(DetectorConstruction* dc){ fDC = dc; GetDCinfo();}
  G4double GetTn(){ return fTn; }
  G4double GetEn(){ return fEn; }
private:
  G4ParticleGun* fParticleGun;
  G4GeneralParticleSource* fGPS;
  PrimaryGeneratorMessenger* fGunMessenger;
  AnalysisManager*           fAM;
  DetectorConstruction*      fDC;
  G4ParticleTable*           fParticleTable;
  G4ParticleDefinition*      fPDef;
  G4int                      fGenMode;   
  G4int                      fNGenBranches;   
  G4int                      fNevent;
  G4float                    fVx;
  G4float                    fVy;
  G4float                    fVz;
  G4float                    fPxp;
  G4float                    fPyp;
  G4float                    fPzp;
  G4float                    fEp;
  G4float                    fTp;
  char*                      fFName;
  char*                      fHName;
  G4ThreeVector              fPos;
  G4ThreeVector              fPos0;
  G4ThreeVector              fDPos;
  G4double                   fLamda;
  G4double                   fFP;
  G4double                   fdTOF;
  G4double                   fSkew;
  G4double                   fT0n;
  G4double                   fEn;
  G4double                   fTn;
  G4int                      fRunMode;
  G4double*                  fAlphaLines;
  G4int                      fNAlphaLines;
  G4double                   fScPx,fScPy,fScPz;
  G4double                   fScX,fScY,fScZ;
};

#endif /*PrimaryGeneratorAction_h*/
