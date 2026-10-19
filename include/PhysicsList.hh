// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Setup Physics List
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class PhysicsListMessenger;
class G4VPhysicsConstructor;
class G4VModularPhysicsList;
class DetectorConstruction;

//--------------------------------------------------------------------------
class PhysicsList : public G4VUserPhysicsList
{
  public:
  PhysicsList(DetectorConstruction*);
  virtual ~PhysicsList();

  public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

  //these methods Construct physics processes and register them
  void ConstructOptical();
  void AddThermalNeutron();

  //for the Messenger 
  void SetVerbose(G4int verb){ fVerboseLevel = verb; }
  void SetNbOfPhotonsCerenkov(G4int nc){ fMaxNumPhotonStep = nc; }
 
  private:
  G4int                fVerboseLevel;
  PhysicsListMessenger* fMessenger;
  G4VModularPhysicsList* fHadronList;
  G4VPhysicsConstructor* fParticleList;
  G4VPhysicsConstructor* fEMPhysicsList;
  G4VPhysicsConstructor* fEMXPhysicsList;
  DetectorConstruction* fDC;
  G4int fMaxNumPhotonStep;
  G4int fRunMode;
  G4int fHadronMode;
  G4ParticleTable::G4PTblDicIterator* fParticleIterator;
};


#endif /* PhysicsList_h */
