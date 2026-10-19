// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Write hits info to root file
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "globals.hh"
#include "AnalysisMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include <string>

//---------------------------------------------------------------------------

class TTree;
class TFile;
class TH1D;
class TF1;

//---------------------------------------------------------------------------

class AnalysisManager {

public:

  AnalysisManager();
  ~AnalysisManager();

  void InitOutput();
  void ZeroArray();
  void FillPDArray( G4int);

  void FillTree( );
  
  void SetOutFileName(const std::string&, G4bool lockForCommandLine = false);
  void SetInput1D (const char*, const char* ) ;

  void SetPDPDef  (G4ParticleDefinition* pdp) { fPDpdef = pdp;    }
  void SetPDP3    (G4ThreeVector  pdmom)      { fPDp3   = pdmom;  }
  void SetPDPos   (G4ThreeVector  pdpos)      { fPDpos  = pdpos;  }
  void SetPDVertex(G4ThreeVector  pdvert)     { fPDvertex  = pdvert; }
  void SetPDTime  (G4double       pdtime)     { fPDtime = pdtime; }
  void SetPDke    (G4double       pdke)       { fPDke = pdke; }
  void SetPDEdep  (G4double       pdedep)     { fPDedep = pdedep; }
  void SetPDVID   (G4int          pdvid)      { fPDvid_single  = pdvid;  }
  void SetEscale  (G4double       escale)     { fEscale  = escale;  }
  void SetPDT0n   (G4double       t0n)        { fPDT0n  = (G4float)t0n;  }
  void SetPDE0n   (G4double       e0n)        { fPDE0n  = (G4float)e0n;  }
  TH1D* GetH1Input(){ return fH1Input; }
  void SetSkGaus(G4double, G4double, G4double );
  G4double SampleH1();
  G4double SampleTOF();
  // G4double SkPk(G4double*,G4double*);

private:
  
  AnalysisMessenger*    fAnaMessenger;
  std::string           fOutFileName;
  G4bool                fOutFileNameLocked = false;
  TFile*                fROOTfile = nullptr;
  TTree*                fROOTtree = nullptr;
  TH1D*                 fH1Input = nullptr;
  TF1*                  fSkGaus = nullptr;

  // Single hit variables
  G4ParticleDefinition* fPDpdef;
  G4ThreeVector         fPDp3;
  G4ThreeVector         fPDpos;
  G4ThreeVector         fPDvertex;
  G4double              fPDtime;
  G4double              fPDke;
  G4double              fPDedep;
  G4int                 fPDvid_single;
  G4double              fEscale;
  G4float               fPDT0n;
  G4float               fPDE0n;
  G4double              fGausSig = 0;

  // Event array variables
  //static const G4int    kMaxHits = 4096;

  G4int                 fPDNhits = 0;
  G4int                 fPDpdg[EMaxHits];
  G4int                 fPDvid[EMaxHits];
  G4float               fPDm[EMaxHits];
  G4float               fPDt[EMaxHits];
  G4float               fPDx[EMaxHits];
  G4float               fPDy[EMaxHits];
  G4float               fPDz[EMaxHits];
  G4float               fPDvx[EMaxHits];
  G4float               fPDvy[EMaxHits];
  G4float               fPDvz[EMaxHits];
  G4float               fPDpx[EMaxHits];
  G4float               fPDpy[EMaxHits];
  G4float               fPDpz[EMaxHits];
  G4float               fPDth[EMaxHits];
  G4float               fPDph[EMaxHits];
  G4float               fPDp[EMaxHits];
  G4float               fPDE[EMaxHits];
  G4float               fPDEd[EMaxHits];

};

#endif

//---------------------------------------------------------------------------
