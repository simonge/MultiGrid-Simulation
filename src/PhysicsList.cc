// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Setup Physics List
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
//
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC_AllHP.hh"
#include "QGSP_BIC_AllHP_EmLow.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPThermalScattering.hh"
//
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4ProcessVector.hh"
//
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4Decay.hh"
#include <iomanip>
#include <assert.h>
#include "DetectorConstruction.hh"

PhysicsList::PhysicsList(DetectorConstruction* dc) 
 : G4VUserPhysicsList(),
   fVerboseLevel(0), fMaxNumPhotonStep(25)
{
  fDC = dc;
  fRunMode = dc->GetRunMode();
  fHadronMode = dc->GetHadronMode();
  fMessenger = new PhysicsListMessenger(this);
  fParticleIterator = NULL;
  fHadronList = NULL;
  fEMPhysicsList = fEMXPhysicsList = NULL; 
  fParticleList = new G4DecayPhysics("decays");
  switch(fHadronMode){
  case EFTFP_BERT_HP:
    fHadronList = new FTFP_BERT_HP();
    break;
  case EQGSP_BERT_HP:
    fHadronList = new QGSP_BIC_AllHP_EmLow();
    break;
  case EQGSP_BIC_HP:
    fHadronList = new QGSP_BIC_AllHP();
    break;
  default:
    break;
  }
}

//--------------------------------------------------------------------------
PhysicsList::~PhysicsList() {
}

//--------------------------------------------------------------------------
void PhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
}

//--------------------------------------------------------------------------
void PhysicsList::ConstructProcess()
{
  if(fHadronMode == ENoHadron){
      AddTransportation();
      //fEMPhysicsList = new G4EmStandardPhysics();
      fEMPhysicsList = new G4EmLowEPPhysics();
      fEMXPhysicsList = new G4EmExtraPhysics("extra EM");
      fEMPhysicsList->ConstructProcess();
      fEMXPhysicsList->ConstructProcess();
      fParticleList->ConstructProcess();
  }
  else{
    fHadronList->ConstructProcess();
  }
  if(fDC->IsThermal())
    AddThermalNeutron();
  if(fDC->IsOptical())
    ConstructOptical();
}

//--------------------------------------------------------------------------
void PhysicsList::ConstructOptical()
{
  fParticleIterator = GetParticleIterator();
  //
  G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");
  cerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
  cerenkovProcess->SetMaxBetaChangePerStep(10.0);
  cerenkovProcess->SetTrackSecondariesFirst(true);
  G4Scintillation* scintillationProcess = new G4Scintillation("Scintillation");
  //scintillationProcess->SetScintillationYieldFactor(1.);
  scintillationProcess->SetTrackSecondariesFirst(true);
  // scintillationProcess->SetScintillationTime(45.*ns);//check
  G4OpAbsorption* absorptionProcess = new G4OpAbsorption();
  G4OpRayleigh* rayleighScatteringProcess = new G4OpRayleigh();
  G4OpMieHG* mieHGScatteringProcess = new G4OpMieHG();
  G4OpBoundaryProcess* boundaryProcess = new G4OpBoundaryProcess();

  cerenkovProcess->SetVerboseLevel(fVerboseLevel);
  scintillationProcess->SetVerboseLevel(fVerboseLevel);
  absorptionProcess->SetVerboseLevel(fVerboseLevel);
  rayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  mieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
  boundaryProcess->SetVerboseLevel(fVerboseLevel);
  
  // Use Birks Correction in the Scintillation process
  if(!G4Threading::IsWorkerThread())
  {
    G4EmSaturation* emSaturation =
              G4LossTableManager::Instance()->EmSaturation();
      scintillationProcess->AddSaturation(emSaturation);
  }

  fParticleIterator->reset();
  while( (*fParticleIterator)() ){
    G4ParticleDefinition* particle = fParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (cerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(cerenkovProcess);
      pmanager->SetProcessOrdering(cerenkovProcess,idxPostStep);
    }
    if (scintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(scintillationProcess);
      pmanager->SetProcessOrderingToLast(scintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(scintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(absorptionProcess);
      pmanager->AddDiscreteProcess(rayleighScatteringProcess);
      pmanager->AddDiscreteProcess(mieHGScatteringProcess);
      pmanager->AddDiscreteProcess(boundaryProcess);
    }
  }
}

//--------------------------------------------------------------------------
void PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  //
  //SetCutsWithDefault();
  SetCutValue(0.01*mm, "gamma");
  SetCutValue(0.01*mm, "e-");
  SetCutValue(0.01*mm, "e+");
  SetCutValue(0.0*mm,"proton");
  SetCutValue(0.0*mm,"neutron");

  if (verboseLevel>0) DumpCutValuesTable();
}
 
//----------------------------------------------------------------------------
void PhysicsList::AddThermalNeutron()
{
  // Find the neutron elastic process
  // and the high-precision low-energy scattering model (both required)
  // Add thermal neutron scattering scattering to the process
  // Set min HP energy and max thermal energy to 4 eV
  //
  G4ParticleDefinition* nPD = G4Neutron::Definition();
  // Find the neutron elastic process
  G4ProcessVector* pvec = nPD->GetProcessManager()->GetProcessList();
  G4HadronicProcess* nEl = 0;
  for (G4int i=0; i<pvec->size(); i++) {
    if ((*pvec)[i]->GetProcessSubType() != fHadronElastic) continue;
    nEl = dynamic_cast<G4HadronicProcess*>((*pvec)[i]);
    break;
  }
  assert(nEl != 0); //ensure its found
  G4HadronicInteraction* nElHP =
    G4HadronicInteractionRegistry::Instance()->FindModel("NeutronHPElastic");
  assert(nElHP != 0);
  nElHP->SetMinEnergy(4.0*eV);
  G4NeutronHPThermalScattering* nTherm = new G4NeutronHPThermalScattering;
  nTherm->SetMaxEnergy(4.0*eV);
  nEl->RegisterMe(nTherm);
  nEl->AddDataSet(new G4NeutronHPThermalScatteringData);
  return;
}
