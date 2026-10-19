// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Generate primary event
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"
#include "G4String.hh"
/*
#include "TObjArray.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
*/
//--------------------------------------------------------------------------
PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  fGPS = new G4GeneralParticleSource();  //isotropic source

  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);

  //default kinematic
  //
  fParticleTable = G4ParticleTable::GetParticleTable();
  fPDef = fParticleTable->FindParticle("neutron");

  fParticleGun->SetParticleDefinition(fPDef);
  fParticleGun->SetParticleTime(0.0*ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.)); //x, y, z directions
  fParticleGun->SetParticleEnergy(500.0*keV);
}

//--------------------------------------------------------------------------
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
  delete fGPS;
}

//--------------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //G4double th,phi,ux,uy,uz;
  switch(fGenMode) {

  case EPGA_GPS:
    fPDef        = fGPS->GetParticleDefinition();
    fGPS->GeneratePrimaryVertex(anEvent); //fParticleGun
    fVx          = fGPS->GetParticlePosition().getX();
    fVy          = fGPS->GetParticlePosition().getY();
    fVz          = fGPS->GetParticlePosition().getZ();
    fPxp         = fGPS->GetParticleMomentumDirection().getX();
    fPyp         = fGPS->GetParticleMomentumDirection().getY();
    fPzp         = fGPS->GetParticleMomentumDirection().getZ();
    fEp          = fGPS->GetParticleEnergy()/GeV;
    fTp          = fGPS->GetParticleTime();
    break;
    
  case EPGA_Chop:
    GenChoppedBeam(anEvent);
    break;
  case EPGA_Alpha:
    GenAlpha(anEvent);
    break;
  default:
    G4cerr << "Unknown mode given to PrimiaryGeneratorAction" << G4endl; 
    exit(1);
 }	
}

//--------------------------------------------------------------------------
void PrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

//--------------------------------------------------------------------------
void PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 fParticleGun->SetParticlePolarization(polar);
}

//--------------------------------------------------------------------------
G4int PrimaryGeneratorAction::GetNEvents() {
  //return static_cast<G4int>(fGenTree->GetEntries());
  return 0;
}

//--------------------------------------------------------------------------
void PrimaryGeneratorAction::SetMode(G4int mode){
  fGenMode = mode;
}
//--------------------------------------------------------------------------
void PrimaryGeneratorAction::SetFName(G4String name){ 
  char* fn = (char*)name.data();
  fFName = new char(strlen(fn)+1);
  strcpy(fFName,fn);
}
//--------------------------------------------------------------------------
void PrimaryGeneratorAction::SetHName(G4String name){ 
  char* fn = (char*)name.data();
  fHName = new char(strlen(fn)+1);
  strcpy(fHName,fn);
  fAM->SetInput1D(fFName,fHName);
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GenChoppedBeam(G4Event* anEvent)
{
  // Chopped-beam time
  //
  fVx = fPos.getX();
  fVy = fPos.getY();
  fVz = fPos.getZ();
  fVx = fVx + fDPos.getX()*(-1 + 2*G4UniformRand());
  fVy = fVy + fDPos.getY()*(-1 + 2*G4UniformRand());
  fVz = fVz + fDPos.getZ()*(-1 + 2*G4UniformRand());
  G4double v0,t0,lam;
  // fixed lambda
  if(fLamda > 0){
    if(fSkew == 0){
      t0 = G4RandGauss::shoot(fT0n,fdTOF); 
    }
    else
      t0 = fAM->SampleTOF();
    v0 = 1000000*fFP/t0;
  }
  else if(fLamda == -8888){
    t0 = fAM->SampleH1();
    v0 = 1000*fFP/t0;
  }
  // lambda distribution
  else if(fLamda == -9999){
    lam = fAM->SampleH1();
    v0 = csq*hp0/(mn0*lam*1e-10);
    t0 = 1000000*fFP/v0;
  }
  // flat lambda distribution
  else{
    lam = G4UniformRand()*-fLamda;
    v0 = csq*hp0/(mn0*lam*1e-10);
    t0 = 1000000*fFP/v0;
  }
  G4double e0 = 1000000*0.5*mn1*v0*v0/csq;
  fTn = t0;
  fEn = e0*eV;
  fParticleGun->SetParticlePosition(G4ThreeVector(fVx,fVy,fVz));
  fParticleGun->SetParticleEnergy(fEn);
  G4ThreeVector pdir = G4ThreeVector(0,0,1);
  fParticleGun->SetParticleMomentumDirection( pdir );
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GenAlpha(G4Event* anEvent)
{
  fVx = fScPx;
  fVy = fScPy;
  fVz = fScPz;
  fVx = fVx + fScX*(-1 + 2*G4UniformRand());
  fVy = fVy + fScY*(-1 + 2*G4UniformRand());
  fVz = fVz + fScZ*(-1 + 2*G4UniformRand());
  G4double cosTh = 2*G4UniformRand() - 1.0;
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4double sinTh = sqrt(1. - cosTh*cosTh);
  G4double ux = sinTh*cos(phi);
  G4double uy = sinTh*sin(phi);
  G4double uz = cosTh;
  fParticleGun->SetParticlePosition(G4ThreeVector(fVx,fVy,fVz));
  G4int ie = (int)(fNAlphaLines*G4UniformRand());
  fEn = fAlphaLines[ie]*MeV;
  fParticleGun->SetParticleEnergy(fEn);
  G4ThreeVector pdir = G4ThreeVector(ux,uy,uz);
  fParticleGun->SetParticleMomentumDirection( pdir );
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GetDCinfo()
{
  // Get parameters from DetectorConstruction read from file SoNDe.dat
  //
  fRunMode = fDC->GetRunMode();
  //if((fRunMode == 5) || (fRunMode == 6)){
  //fPDef = fParticleTable->FindParticle("alpha");
  fPDef = fParticleGun->GetParticleDefinition();
   // }
  fPos = fPos0 = fDC->GetPos0();
  fDPos = fDC->GetDPos();
  fLamda = fDC->GetLamda();
  if(fLamda > 0)
    fGenMode = EPGA_Chop;
  else if((fLamda == -9999)||(fLamda == -8888)){
    fAM->SetInput1D(fDC->GetFName(), fDC->GetHName());
    fAM->SetEscale(1.0);
    fGenMode = EPGA_Chop;
  }
  else
    fGenMode = EPGA_Chop;
  fFP = fDC->GetFP()/1000;              // in m
  fSkew = fDC->GetSkew();               // skew in non-normal distribution
  fdTOF = fDC->GetdTOF();               // sigma in us
  G4double vn = csq*hp0/(mn0*fLamda*1e-10);
  fT0n = 1.0e6*fFP/vn;
  if(fSkew == 0.0)
    fdTOF = fdTOF/2.355;                // was FWHM of Gaussian
  else
    fAM->SetSkGaus(fT0n,fdTOF,fSkew);
  fNAlphaLines = fDC->GetNAlphaLines();
  if(fNAlphaLines > 0.0){
    fGenMode = EPGA_Alpha;
    fAlphaLines = fDC->GetAlphaLines();
    fDC->GetScatter(fScPx,fScPy,fScPz,fScX,fScY,fScZ);
    fPDef = fParticleTable->FindParticle("alpha");
    fParticleGun->SetParticleDefinition(fPDef);
  }
}
