// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Record hits in sensitive detector
// Discard G4 hits collection stuff
// JRM Annand 20/08/2018
//--------------------------------------------------------------------------
#include "G4RunManager.hh"
#include "DetectorSD.hh"
#include "DetectorHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleTypes.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scintillation.hh" //additions from here
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVReplica.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4MaterialTable.hh"
#include "G4String.hh"
#include "globals.hh"

#include "G4ParticleDefinition.hh"

#include <cmath>
//---------------------------------------------------------------------------

DetectorSD::DetectorSD(G4String name,G4int Nmax,G4int Ncell,G4int Ncol,G4int Nbox):
  G4VSensitiveDetector(name)
{
  fNmax  = Nmax;
  fNcell = Ncell;                // # voxels in single MG
  fNvess = Ncell * Ncol;         // # voxels in column of MG
  fNbox = fNvess * Nbox;         // # voxels in box (6 TREX, 0 otherwise)
  fNhits      = 0;
  fHits       = new G4int[fNmax];
  //fHitID       = new G4int[fNmax];
  fDetectorHits = new DetectorHit[fNmax];
  //for(G4int i=0; i<fNmax; i++) fHitID[i] = -1;
  for(G4int i=0; i<fNmax; i++) fHits[i]  = 0;
}

//---------------------------------------------------------------------------

DetectorSD::~DetectorSD()
{;}

//---------------------------------------------------------------------------

void DetectorSD::Initialize(G4HCofThisEvent*)
{
  return;
}

//---------------------------------------------------------------------------

G4bool DetectorSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  // Multigrid related stuff
  if(fNhits >= fNmax)
    return false;
  G4Track* aTrack = aStep->GetTrack();
  G4int pdg = aTrack->GetDynamicParticle()->GetPDGcode();
  G4ThreeVector pos;
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* gas = theTouchable->GetVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();     // energy
  G4int vgas = gas->GetCopyNo();
  G4int vid,vcell,vgrid,vvess,vbox,st1,st2;
  switch(vgas){
  case 1000:
    // copyN0 1000 is ArCO2 inside multigrid cell
    // copyNo 1001 - 1096 cells in single grid
    // copyNo 2000 - N for N grids per column
    // copyNo 100 - 999 no of columns
    // copyNo 10 - 99 no of boxes (TREX only)
    if ((edep/keV == 0.)) return false;
    vcell = theTouchable->GetVolume(1)->GetCopyNo() - 1001;  // voxel #
    vgrid = theTouchable->GetVolume(2)->GetCopyNo() - 2000;  // grid #
    vvess = theTouchable->GetVolume(3)->GetCopyNo() - 100;   // column #
    vbox = theTouchable->GetVolume(4)->GetCopyNo() - 10;     // box #
    vid = vbox*fNbox + vvess*fNvess + vgrid*fNcell + vcell;  // voxel ID
    break;
  case 1001:
    // 3He Tube
    if ((edep/keV == 0.)) return false;
    vid = theTouchable->GetVolume(1)->GetCopyNo();
    break;
  case 10:
    // Scattering sample
    if(pdg != 2112)
      return false;
    st1 = aStep->GetPreStepPoint()->GetStepStatus();
    st2 = aStep->GetPostStepPoint()->GetStepStatus();
    if((st1==1)&&(st2==1))
       return false;
    //pos = aStep->GetPreStepPoint()->GetPosition();
    //pos = aStep->GetPostStepPoint()->GetPosition();
    vid = 1000000;
    break;
  case 20:
  case 21:
  case 22:
  case 23:
    // Test scattering sample
    if(pdg == 2112){
      st1 = aStep->GetPreStepPoint()->GetStepStatus();
      st2 = aStep->GetPostStepPoint()->GetStepStatus();
      if((st1==1)&&(st2==1))
	return false;
      vid = 3;
    }
    else{
      if ((edep/keV == 0.)) return false;
      vid = 1;
    }
    break;
  case 30:
    // Test detection sphere
    vid = 2;
    if ((edep/keV == 0.)) return false;
    break;
  case 40:
    // Beam Monitor
    vid = 2000000;
    if(edep/keV == 0.0)
       return false;
    break;
  case 5:
    // gas outside of multigrid
    if(pdg ==  1000020040){
      vid = -1;
      if (edep/keV == 0.)
	return false;
    }
    else
      vid = -2;
    break;
  default:
    vid = vgas;
    break;
  }
  pos = aStep->GetPostStepPoint()->GetPosition();
  DetectorHit* hit = fDetectorHits + fNhits;
  hit->SetPDG(pdg);
  hit->SetVid(vid);
  hit->AddEnergy(edep);
  hit->SetMomentum(aTrack->GetMomentum()/MeV);
  hit->SetKE(aTrack->GetKineticEnergy()/MeV + edep);
  // hit->SetKE(aTrack->GetVertexKineticEnergy()/MeV);
  hit->SetPDef(aTrack->GetDefinition());
  hit->SetPos(pos);
  hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
  hit->SetParentID(aTrack->GetParentID());
  fHits[fNhits]=vid;
  fNhits++;
  return true;
}

//---------------------------------------------------------------------------

void DetectorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  // Hits zeroing is done from Event Action 
  return;
}

//---------------------------------------------------------------------------
