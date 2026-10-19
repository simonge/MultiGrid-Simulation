// Geant-4 simulation of the Multigrid 
// for the European Spallation Source in Lund Sweden
// Simulation includes calculation of
// Gamma-ray response
// Fast-neutron response
// Thermal-neutron response
// Transport of electrons through MWPC gas
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------

//#include "NXSG4.hh"
#include "G4NCrystal/G4NCrystal.hh"
#include "G4GDMLParser.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
//#include "StackingAction.hh"
#include "AnalysisManager.hh"
#include "EventAction.hh"
#include "CommandLineParser.h"

#define G4VIS_USE 1
#define G4UI_USE 1

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include <exception>
#include <fstream>
#include <iostream>
#include <string>

//--------------------------------------------------------------------------
int main(int argc,char** argv)
{
  MG::CommandLineOptions options;
  try {
    options = MG::ParseCommandLine(argc, argv);
  } catch (const std::exception& error) {
    G4cerr << "Error: " << error.what() << G4endl;
    MG::PrintUsage(std::cerr, argv[0]);
    return 2;
  }
  if (options.help) {
    MG::PrintUsage(std::cout, argv[0]);
    return 0;
  }

  const G4String macro = options.macro;
  const G4String session = options.session;
  const G4String geometryFile = options.geometryFile;
  G4long rseed = static_cast<G4long>(options.seed);
  // NCrystal::libClashDetect();//Detect broken installation

  // Choose the Random engine
  // and set the seed
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  if(rseed == 0){
    time_t systime = time(NULL);
    rseed = (long)systime;
  }
  G4Random::setTheSeed(rseed);
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;


  // Set mandatory initialization classes
  //
  // Detector construction
  DetectorConstruction* dc;
  if(!geometryFile.empty())
    dc = new DetectorConstruction(geometryFile.c_str());
  else
    dc = new DetectorConstruction();

  runManager-> SetUserInitialization(dc);

  if (!options.gdmlOutputFile.empty()) {
    runManager->InitializeGeometry();
    G4GDMLParser gdmlParser;
    gdmlParser.Write(options.gdmlOutputFile, dc->GetWorldVolume());
    std::ifstream gdmlFile(options.gdmlOutputFile, std::ios::binary | std::ios::ate);
    if (!gdmlFile || gdmlFile.tellg() == std::streampos(0)) {
      G4cerr << "Failed to write GDML geometry to "
             << options.gdmlOutputFile << G4endl;
      delete runManager;
      return 1;
    }
    delete runManager;
    return 0;
  }

  // Physics list
  runManager-> SetUserInitialization(new PhysicsList(dc));
  // User action initialization
  PrimaryGeneratorAction* pga = new PrimaryGeneratorAction();
  runManager->SetUserAction(pga);
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new SteppingAction());
  //runManager->SetUserAction(new StackingAction());
  AnalysisManager* am = new AnalysisManager();
  if (options.hasOutputFile)
    am->SetOutFileName(options.outputFile, true);
  runManager->SetUserAction(new EventAction(am));
  pga->SetAnalysisManager(am);
  pga->SetDetectorConstruction(dc);

  // Initialize G4 kernel moved to run macros
  runManager->Initialize();
  //Install NXSG4:
  //NXSG4::installOnDemand();
  //Install G4NCrystal:
  if(dc->IsCrysMat())
    G4NCrystal::installOnDemand();

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   
  if (!macro.empty()) {
     // Batch mode
     G4String command = "/control/execute ";
     UImanager->ApplyCommand(command+macro);
  }
  else if (options.events >= 0) {
     UImanager->ApplyCommand("/run/beamOn " + std::to_string(options.events));
  }
  else // Define UI session for interactive mode
  {
#ifdef G4UI_USE
      char* uiArgv[] = {argv[0], nullptr};
      G4UIExecutive * ui = new G4UIExecutive(1,uiArgv,session);
#ifdef G4VIS_USE
       UImanager->ApplyCommand("/control/execute macros/visMG.mac");
#endif
     ui->SessionStart();
     delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete am;
  //G4NCrystal::Manager::cleanup();//delete manager singleton, unref cached ncrystal objects (for valgrind).

  return 0;
}
