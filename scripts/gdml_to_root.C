#include <TGeoManager.h>
#include <TSystem.h>
#include <TError.h>

void gdml_to_root()
{
  const char* gdmlFile = gSystem->Getenv("MG_GDML_INPUT");
  const char* rootFile = gSystem->Getenv("MG_ROOT_OUTPUT");
  if (!gdmlFile || !*gdmlFile || !rootFile || !*rootFile) {
    Error("gdml_to_root", "MG_GDML_INPUT and MG_ROOT_OUTPUT must be set");
    gSystem->Exit(1);
    return;
  }

  TGeoManager* geometry = TGeoManager::Import(gdmlFile);
  if (!geometry || !geometry->GetTopVolume()) {
    Error("gdml_to_root", "Could not import GDML geometry from %s", gdmlFile);
    gSystem->Exit(1);
    return;
  }

  if (geometry->Export(rootFile) <= 0) {
    Error("gdml_to_root", "Could not export ROOT geometry to %s", rootFile);
    gSystem->Exit(1);
  }
}
