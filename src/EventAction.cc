// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Write hits to ROOT file at end of event
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "EventAction.hh"
#include "DetectorHit.hh"
#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

//---------------------------------------------------------------------------

EventAction::EventAction(AnalysisManager* ana)
{
  fAM = ana;
  fPGA = (PrimaryGeneratorAction*)
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  fDC = (DetectorConstruction*)
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
}

//---------------------------------------------------------------------------

EventAction::~EventAction()
{;}

//---------------------------------------------------------------------------

void EventAction::BeginOfEventAction(const G4Event* evt)
{ 
  // Init some pointers if 1st event
  if( evt->GetEventID() == 0){
    fAM->InitOutput();
    fDetHit = fDC->GetDetSD()->GetDetHits();
    fHits = fDC->GetDetSD()->GetHits();
  }
}

//---------------------------------------------------------------------------

void EventAction::EndOfEventAction(const G4Event* evt)
{
  // Transfer SD hits info over to Analysis Manager
  // where its written to root file
  G4int event_id = evt->GetEventID();  
  if (event_id%10000 == 0)
    printf("MG event: %d\n",event_id);
  
  fAM->ZeroArray();
  DetectorHit* hit;
  G4int nhits = fDC->GetDetSD()->GetNhits();
  for(G4int j = 0; j < nhits; j++) {
    hit = fDetHit + j;
    fAM->SetPDT0n(fPGA->GetTn());
    fAM->SetPDE0n(fPGA->GetEn());
    fAM->SetPDPDef ( (G4ParticleDefinition*)    hit->GetPDef() );
    fAM->SetPDP3   ( (G4ThreeVector)            hit->GetMom()  );
    fAM->SetPDPos  ( (G4ThreeVector)            hit->GetPos()  );
    fAM->SetPDVertex( (G4ThreeVector)           fPGA->GetVertex() );
    fAM->SetPDTime ( (G4double)                 hit->GetTime() );
    fAM->SetPDke  ( (G4double)                  hit->GetKE() );
    fAM->SetPDEdep ( (G4double)                 hit->GetEdep() );
    fAM->SetPDVID  ( (G4int)                    hit->GetVid()  );
    fAM->FillPDArray( j ); 
  }
  // record data in TTree if there any hits
  if(nhits > 0)
    fAM->FillTree();
  fDC->GetDetSD()->ZeroHits();
}
  

//---------------------------------------------------------------------------
