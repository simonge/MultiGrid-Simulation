// Geant-4 simulation of Multigrid type detectors
// for the European Spallation Source in Lund Sweden
// Detector geometry, materials including optical properties
// JRM Annand 26/11/23
//--------------------------------------------------------------------------
#include "DetectorConstruction.hh"
#include "DetectorSD.hh"

//#include "NXSG4.hh"
#include "G4NCrystal/G4NCrystal.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scintillation.hh" //additions from here
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVReplica.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4MaterialTable.hh"
#include "G4String.hh"
#include "G4MultiUnion.hh"
#include "G4NistManager.hh"
#include "globals.hh"
// Colours
G4VisAttributes* black   = new G4VisAttributes( G4Colour(0.0,0.0,0.0) );
G4VisAttributes* white   = new G4VisAttributes( G4Colour(1.0,1.0,1.0) );
G4VisAttributes* red     = new G4VisAttributes( G4Colour(1.0,0.0,0.0) );
G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0) );
G4VisAttributes* purple  = new G4VisAttributes( G4Colour(1.0,0.0,1.0) );
G4VisAttributes* green   = new G4VisAttributes( G4Colour(0.0,1.0,0.0) );
G4VisAttributes* grey    = new G4VisAttributes( G4Colour(0.5,0.5,0.5) );
G4VisAttributes* yellow  = new G4VisAttributes( G4Colour(1.0,1.0,0.0) );
G4VisAttributes* cyan    = new G4VisAttributes( G4Colour(0.0,1.0,1.0) );

G4double pixd[] = {6.25, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.25};
FILE* CoordFile;

DetectorConstruction::DetectorConstruction(const char* readfile)
 : G4VUserDetectorConstruction()
{
  // Default initialisation of some dimensions
  fExpX = fExpY = 500*cm; fExpZ = 500*cm;
  fRunMode = ERunThermal;   // default optical transport
  fHadronMode = ENoHadron;
  fIsCrysMat = false;
  fLamda = fFP = fdTOF = 0.0;
  // Containment vessel default off
  fVx = fVy = fVz =  fVth = 0;
  fNcol = 1;
  fIsMB = fIsTrex = false;
  fMBth = fMBth1 = 0.0;
  //
  CoordFile = fopen("CoordFile.txt","w+");
  ReadParameters(readfile);    // dimensions from file
  fIsThermal = fIsOptical = false;
  if(fIsTarget && (fNy==0))
    fIsScatter = false;
  switch(fRunMode){
  case ERunOptical:     // optical only
  case ERunFastOpt:
    fIsOptical = true;
    break;
  case ERunFast:       // fast neutrons, gamma rays
    break;
  case ERunThermal:
    fIsThermal = true;
    break;
  case ERunThermalOpt:
    fIsThermal = true;
    fIsOptical = true;
    break;
  }
}

DetectorConstruction::~DetectorConstruction(){;}

//---------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Construct the geometry
  // Define materials 1st
  Materials();
  // Outer concrete on exp volume
  // The experimental volume, vacuum (0) or air@STP (1) filled
  G4Box* expBox = new G4Box("World",fExpX,fExpY,fExpZ);
  if(fExpMat == 0)
    fExpLog = new G4LogicalVolume(expBox,fVac,"World",0,0,0);
  else if(fExpMat == 1)
    fExpLog = new G4LogicalVolume(expBox,fAir,"World",0,0,0);
  fExpPhys = new G4PVPlacement(0,G4ThreeVector(),fExpLog,"World",0,false,0);
  //
  if(fIsEMMAtest){
    EMMAtest();
    for(G4int n=0; n<fNvess; n++)
      MultiGrid(n);
  }
  else if(fIsTrex)
    TRbox();
  else if(fNy){
    if(fVx)
      Vessel();
    // Components
    MultiGrid();
  }
  else if(fNHe3){
    // He3 tubes
    He3Tube();
  }
  else if(fIsScatter)
    // Scattering test
    ScatterTest();
  //
  if(fColSz > 0)
    Collimator();           // insert collimator
  if(fBMsz > 0)
    BeamMonitor();          // insert beam monitor
  if(fTabSx > 0)            // extra shielding
    EMMAshield();
  // Sensitive detector
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4int nbox = 0;
  if(fIsTrex)
    nbox = 6;
  fDetSD = new DetectorSD("MGDetectorSD",EMaxHits,fNx*fNz,fNy,nbox);
  SDman->AddNewDetector(fDetSD );
  if(fNy){
    for(G4int i=0; i<nVxSz; i++){
      if(flogArCO2[i])
	flogArCO2[i]->SetSensitiveDetector( fDetSD );
    }
  }
  else if(fNHe3)
    flogHe3->SetSensitiveDetector( fDetSD );
  if(fIsTarget)
    flogTarget->SetSensitiveDetector( fDetSD );
  else if(fIsScatter && (fNy==0)){
    for(G4int im=0; im<fNmat; im++)
      flogScat[im+1]->SetSensitiveDetector( fDetSD );
    if(flogScDet)
      flogScDet->SetSensitiveDetector( fDetSD );
  }
  for(G4int n=0; n<fNmg; n++){
    if(flogVGas[n])
      flogVGas[n]->SetSensitiveDetector( fDetSD );
  }
  if(flogBM)
    flogBM->SetSensitiveDetector( fDetSD );
  if(fnT3)
    flogHe3->SetSensitiveDetector( fDetSD );
  //
  //---------------------------------------------------------------------------
  // Set Visualisation Attributes (wireframe colour)
  //---------------------------------------------------------------------------
  //World
  //fExpLog->SetVisAttributes(G4VisAttributes(false));
  fExpLog->SetVisAttributes(black);
  // Grid
  for(G4int n=0; n<fNmg; n++){
    //if(flogMG[n])
    //flogMG[n]->SetVisAttributes(black);
    if(flogVess[n])
      flogVess[n]->SetVisAttributes(black);
    if(flogVGas[n])
      flogVGas[n]->SetVisAttributes(black);
  }
  if(flogTarget)
    flogTarget->SetVisAttributes(black);
  if(fNHe3){
    flogHe3->SetVisAttributes(cyan);
    flogTube->SetVisAttributes(black);
    flogWire->SetVisAttributes(black);
  }

//always return the physical World
  fclose(CoordFile);
  return fExpPhys;
}

//----------------------------------------------
void DetectorConstruction::Vessel(G4int nmg)
{
  G4double thfr = fVth;
  if(fVfr > 0)
    thfr = fVfr;
  G4double dzV = fVz+0.5*(fVbk+thfr);
  if(fIsTrex){
    G4double dt = fdTh;
    G4double dx = dzV*tan(dt);
    G4double vx = fRad*sin(0.5*dt); // 1/2 chord length
    G4Trd* vBox =
      new G4Trd("vbox",vx,vx+dx,fVy+fVth,fVy+fVth,dzV);
    flogVess[nmg] = new G4LogicalVolume(vBox,fAl,"vLog",0,0,0);
  }
  else{
    G4Box* vBox = new G4Box("vbox",fVx+fVth,fVy+fVth,dzV);
    flogVess[nmg] = new G4LogicalVolume(vBox,fAl,"vLog",0,0,0);
  }
  G4Box* vgBox = new G4Box("vgbox",fVx,fVy,fVz);
  G4double dz = 0.5*(thfr - fVbk);
  flogVGas[nmg] = new G4LogicalVolume(vgBox,fArCO2,"vgLog",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,dz),flogVGas[nmg],"vgExp",flogVess[nmg],
		    false,5);
  G4double th = fThSt;
  G4double dth = fdTh;
  G4int nv = 0;
  for(G4int n=0; n<fNvess; n++){
    G4RotationMatrix* yR = new G4RotationMatrix;
    G4double x = (fRad+dzV) * sin(th);
    G4double y = (fRad+dzV) * cos(th);
    yR->rotateY(-th);
    G4PVPlacement* vess = 
      new G4PVPlacement(yR,G4ThreeVector(x,0,y),flogVess[nmg],"vExp",fExpLog,
			false,n+100);
    G4ThreeVector xyz = -vess->GetFrameTranslation();
    fprintf(CoordFile,"Vess: %d %g %g %g %g\n",
	    n+100,xyz.x(),xyz.y(),xyz.z(),-th);

    th += dth;
    nv++;
    if(n == fNvGap)
      th += fThGap*deg;
    if(nv == fNvGap1){
      if(n != fNvGap)
	th += fThGap1*deg;
      nv = 0;
    }
  }
  fprintf(CoordFile,"Gas:  0 0 0 %g\n",dz);
  Target();
}
//----------------------------------------------
void DetectorConstruction::EMMAtest()
{
  G4double thfr = fVth;
  if(fVfr > 0)
    thfr = fVfr;
  G4double th = 0*deg;
  G4double dzV = fVz+0.5*(fVbk+thfr);
  auto yR = new G4RotationMatrix;
  auto vBox = new G4Box("vbox",fVx+fVth,fVy+fVth,dzV);
  auto vMB0 = new G4Box("vMB0",fVx+fVth,fVy+fVth+thMB,dzV);
  auto vMB1 = new G4Box("vMB1",fVx+fVth+thMB,fVy+fVth+thMB,dzV);
  auto vMBs = new G4SubtractionSolid("vMBs",vMB1,vMB0,0,G4ThreeVector(0,-thMB,0));
  auto vgBox = new G4Box("vgbox",fVx,fVy,fVz);
  G4double dz = 0.5*(thfr - fVbk);
  // yR->rotateY(-90*deg);
  yR->rotateZ(fETrz*deg);
  yR->rotateY(fETry*deg);
  yR->rotateX(fETrx*deg);
  auto lMBs = new G4LogicalVolume(vMBs,fMB,"lMBs",0,0,0);
  lMBs->SetVisAttributes(red);
  for(G4int n=0; n<fNvess; n++){
    flogVess[n] = new G4LogicalVolume(vBox,fAl,"vLog",0,0,0);
    flogVGas[n] = new G4LogicalVolume(vgBox,fArCO2,"vgLog",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,dz),flogVGas[n],"vgExp",flogVess[n],
		      false,5);
    G4PVPlacement* vess = 
      new G4PVPlacement(yR,G4ThreeVector(fETx[n],fETy[n],fETz[n]),flogVess[n],
			"vExp",fExpLog,false,n+100);
    new G4PVPlacement(yR,G4ThreeVector(fETx[n],fETy[n],fETz[n]-thMB*0.0),lMBs,
			"pMBs",fExpLog,false,0);
    auto xyz = -vess->GetFrameTranslation();
    fprintf(CoordFile,"Vess: %d %g %g %g %g\n",
	    n+100,xyz.x(),xyz.y(),xyz.z(),th);

  }
  fprintf(CoordFile,"Gas:  0 0 0 %g\n",dz);
  Target();
  if(fnT3)
    EMMAHe3();
}
//------------------------------------------------------------
void DetectorConstruction::Target(){
  if(fIsTarget){
    G4Material* tgmat;
    //   if(!strcmp("Vanadium",fTgMat))
    //  tgmat = fV;
    // else if(!strcmp("Graphite",fTgMat))
    //  tgmat = fGraf;
    // else if(!strcmp("Polyethylene",fTgMat))
    //   tgmat = fCH2;
    tgmat = fScMatList[fTgMat-1];
    G4RotationMatrix* rot = new G4RotationMatrix;
    if(fnmTg > 0){
      rot->rotateY(fTgAng*deg);
      G4double sz = fmTgSz[0]+fmTgSz[1]+fmTgSz[2];
      auto targ = new G4Box("Target",fmTgSx[0],fmTgSy[0],sz);
      flogTarget = new G4LogicalVolume(targ,fAir,"logTarg",0,0,0);
      G4double z0 = -sz + fmTgSz[0];
      for(G4int n=0; n<fnmTg; n++){
	G4Material* tmat = fScMatList[fmTgMat[n]-1];
	auto tb = new G4Box("tb",fmTgSx[n],fmTgSy[n],fmTgSz[n]);
	auto ltb = new G4LogicalVolume(tb,tmat,"ltb",0,0,0);
	ltb->SetVisAttributes(black);
	new G4PVPlacement(0,G4ThreeVector(0,0,z0),ltb,"ptb",flogTarget,
			  false,10);
	z0 += fmTgSz[n] + fmTgSz[n+1];
      }
    }
    else if(fTgSx == 0){
      auto targ = new G4Tubs("Target",fTgRadI,fTgRad,fTgLen,0.0,360.0*deg);
      flogTarget = new G4LogicalVolume(targ,tgmat,"logTarg",0,0,0);
      rot->rotateX(90*deg);
      //rot->rotateY(90*deg);
    }
    else{
      auto targ = new G4Box("Target",fTgSx,fTgSy,fTgSz);
      flogTarget = new G4LogicalVolume(targ,tgmat,"logTarg",0,0,0);
      rot->rotateY(fTgAng*deg);
    }
    new G4PVPlacement(rot,G4ThreeVector(0,0,0),flogTarget,"pTarg",fExpLog,false,10);
  }
}

//------------------------------------------------------------
void DetectorConstruction::CryoStat(){
  // cryostat, dimensions mm
  if(!fIsCryoStat)
    return;
  G4RotationMatrix* rot = new G4RotationMatrix;
  rot->rotateX(90*deg);
  G4double ir[3]={25,42.5,50};
  G4double irw0=0.5,irw1=1.5,irl=60,irl1=240,irlb,irl2=5,irt=245,endl=2.0;
  char nam[64];
  for(G4int i=0; i<3; i++){
    sprintf(nam,"cr%d",i);
    auto cr = new G4Tubs(nam,ir[i],ir[i]+irw0,ir[i],0.0,360*deg);
    sprintf(nam,"lcr%d",i);
    auto lcr = new G4LogicalVolume(cr,fAl,nam,0,0,0);
    lcr->SetVisAttributes(red);
    sprintf(nam,"pcr%d",i);
    new G4PVPlacement(rot,G4ThreeVector(0,0,0),lcr,nam,fExpLog,false,0);
    //
    sprintf(nam,"cra%d",i);
    auto cra = new G4Tubs(nam,ir[i],ir[i]+irw1,irl,0.0,360*deg);
    sprintf(nam,"lcra%d",i);
    auto lcra = new G4LogicalVolume(cra,fAl,nam,0,0,0);
    lcra->SetVisAttributes(blue);
    sprintf(nam,"pcra%d",i);
    G4double yoff = ir[i] + irl;
    new G4PVPlacement(rot,G4ThreeVector(0,-yoff,0),lcra,nam,fExpLog,false,0);
    sprintf(nam,"end%d",i);
    auto enda = new G4Tubs(nam,0,ir[i]+irw1,endl,0.0,360*deg);
    sprintf(nam,"lend%d",i);
    auto lenda = new G4LogicalVolume(enda,fAl,nam,0,0,0);
    lenda->SetVisAttributes(blue);
    yoff += irl + endl;
    sprintf(nam,"pend%d",i);
    new G4PVPlacement(rot,G4ThreeVector(0,-yoff,0),lenda,nam,fExpLog,false,0);
    //
    sprintf(nam,"crb%d",i);
    irlb = 0.5*(irl1 - ir[i]);
    auto crb = new G4Tubs(nam,ir[i],ir[i]+irw1,irlb,0.0,360*deg);
    sprintf(nam,"lcrb%d",i);
    auto lcrb = new G4LogicalVolume(crb,fAl,nam,0,0,0);
    lcrb->SetVisAttributes(blue);
    sprintf(nam,"pcrb%d",i);
    yoff = ir[i]+irlb;
    new G4PVPlacement(rot,G4ThreeVector(0,yoff,0),lcrb,nam,fExpLog,false,0);
  }
  auto endt = new G4Tubs("endt",0,irt,irl2,0.0,360*deg);
  auto lendt = new G4LogicalVolume(endt,fAl,"lendt",0,0,0);
  lendt->SetVisAttributes(blue);
  new G4PVPlacement(rot,G4ThreeVector(0,irl1+irl2,0),lendt,"pendt",fExpLog,false,0);
}

//----------------------------------------------
void DetectorConstruction::He3Tube()
{
  fWrR = 0.025*mm;
  G4Tubs* Tube = new G4Tubs("He3Tube",0.0,fHe3r+fHe3th,fHe3l,0.0,360.0*deg);
  G4Tubs* He3Gas = new G4Tubs("He3Gas",0.0,fHe3r,fHe3l,0.0,360.0*deg);
  G4Tubs* Wire = new G4Tubs("Wire",0.0,fWrR,fHe3l,0.0,360.0*deg);  
  flogHe3 = new G4LogicalVolume(He3Gas,fHe3CO2,"logHe3",0,0,0);
  flogTube = new G4LogicalVolume(Tube,fFe,"logTube",0,0,0);
  flogWire = new G4LogicalVolume(Wire,fW,"logWire",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),flogWire,"WireExp",flogHe3,false,999);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),flogHe3,"GasExp",flogTube,false,1001);
  G4double th = fThSt;
  G4double dth = fdTh;
  for(G4int n=0; n<fNHe3; n++){
    G4RotationMatrix* xR = new G4RotationMatrix;
    G4double x = fRad * sin(th);
    G4double y = fRad * cos(th);
    xR->rotateX(90*deg);
    G4PVPlacement* phe3 = new G4PVPlacement(xR,G4ThreeVector(x,0,y),flogTube,"vExp",fExpLog,false,n+200);
    G4ThreeVector xyz = -phe3->GetFrameTranslation();
    fprintf(CoordFile,"He3: %d %g %g %g\n",n,xyz.x(),xyz.y(),xyz.z());
    th += dth;
  }
  if(fIsScatter)
    ScatterTest();
}

//----------------------------------------------
void DetectorConstruction::MultiGrid(G4int nmg, G4LogicalVolume* lgas,
				     G4double x0, G4double z0, G4double th0)
{
  G4double gx,gz,hx,hy,hz,hx1=0,mgthr;
  G4LogicalVolume* mgCont;
  char phynm[32];
  if(lgas){
    mgCont = lgas;
    fTRsp0 = z0;
  }
  else if(flogVess[nmg])
    mgCont = flogVGas[nmg];
  else
    mgCont = fExpLog;
  //
  if(flogMG[nmg]==NULL){                      // if column log vol not int
    if(nmg==0)
      mgthr = fMGthr;
    else
      mgthr = fMGthr1;
    if(fIsMB)
      Shield(nmg);
    if(flogMG[nmg]==NULL){
      Grid(mgthr,fIsScat);                    // 1st grid
    }
    gx = fNx*(fMGx + mgthr);
    gz = fNz*(fMGz + fMGth);
    if(flogMB)                                 // extend for rear shield
      gz += fMBth2;
    if(flogMB1)                                // extend for side shield
      gx += 2*fMBth1;
    //
    if(fHspmat >= 0){                          // if H spacers included
      hx = gx;
      hy = fMGy+0.5*fMGsp;
      hz = fTRsp0 + fHsz1;
      HSpacer();
      hx1 = 2*(fHsx - fHsx1);
    }
    fGrx = gx+mgthr+hx1;                       // set container for column
    fGrz = gz+fMGth;
    fGry = fNy*(fMGy + fMGsp);
    auto bxCol = new G4Box("bxCol",fGrx,fGry,fGrz);
    flogMG[nmg] = new G4LogicalVolume(bxCol,fArCO2,"lbxCol",0,0,0);
    flogMG[nmg]->SetVisAttributes(red); 
    
    fphysMG = new G4VPhysicalVolume*[fNy];
    G4int nn = fNy/2;
    G4int ir = fNy%2;  // odd or even
    G4double yoff,xoff;
    if(fNy == 1)
      yoff = 0;
    else{
      yoff = (nn-0.5)*(2*fMGy + fMGsp);
      if( (ir) && (nn) ){
	yoff += fMGy + 0.5*fMGsp;
      }
    }
    G4double yoff0 = yoff;
    nn = fNcol/2;
    ir = fNcol%2;  // odd or even
    G4double mgx = fNx*(fMGx + fMGsp);
    if(fNcol == 1)
      xoff = x0;
    else{
      xoff = x0+ (nn-0.5)*(2*mgx + fMGhsp);
      if( (ir) && (nn) ){
	xoff += mgx + 0.5*fMGhsp;
      }
    }
    nn = 2000;
    //if(fIsTrex)
    //  nn = fNTRcol+2000;
    //
    for(G4int j=0; j<fNcol; j++){
      yoff = yoff0;
      for(G4int i=0; i<fNy; i++){
	sprintf(phynm,"MGphys%d",i);
	new G4PVPlacement(0,G4ThreeVector(0,yoff,-fMBth2),
			  flogGrid,phynm,flogMG[nmg],false,nn);
	//G4ThreeVector xyz = -fphysMG[i]->GetFrameTranslation();
	//fprintf(CoordFile,"MG: %d %g %g %g\n",nn,xyz.x(),xyz.y(),xyz.z());
	//}
	//xoff -= (2*mgx + fMGhsp);
	if(flogMB)
	  //new G4PVPlacement(0,G4ThreeVector(0,yoff,fMposz),flogMB,
	  new G4PVPlacement(0,G4ThreeVector(0,yoff,fMposz-fMBth2),flogMB,
			    "pMB",flogMG[nmg],false,0);
	if(flogMB1){
	  new G4PVPlacement(0,G4ThreeVector(fMposx,yoff,0),flogMB1,
			    "pMB1",flogMG[nmg],false,0);
	  new G4PVPlacement(0,G4ThreeVector(-fMposx,yoff,0),flogMB1,
			    "pMB1",flogMG[nmg],false,0);
	}
	yoff -= (2*fMGy + fMGsp);
	nn++;
      }
    }
  }
  G4RotationMatrix* mgR = new G4RotationMatrix;
  mgR->rotateY(th0);
  sprintf(phynm,"pMGcol%d",fNTRcol);
  new G4PVPlacement(mgR,G4ThreeVector(x0,0,fTRsp0),
		    flogMG[nmg],phynm,mgCont,false,fNTRcol);
  fNTRcol++;
}
/*------------------------
  new G4PVPlacement(0,G4ThreeVector(0,0,-fMBth2),logmg,
			  "RadphysA",logmg1,false,9999);
  if(fHspmat >= 0){
    G4RotationMatrix* hR = new G4RotationMatrix;
    hR->rotateZ(180*deg);
    new G4PVPlacement(0,G4ThreeVector(hx,hy,hz),
		      flogHsp,"pHsp",logmg1,false,0);
    new G4PVPlacement(hR,G4ThreeVector(-hx,hy,hz),
		      flogHsp,"pHsp",logmg1,false,0);
  }
  if(flogMG[0]==NULL){
    flogMG[0] = logmg1;
  }
  else{
    flogMG[1] = logmg1;
  }
*/
//----------------------------------------------
void DetectorConstruction::Grid(G4double mgthr, G4int isScat)
{
  G4double gx = fNx*(fMGx + mgthr);
  G4double gz = fNz*(fMGz + fMGth);
  G4Box* bx1 = new G4Box("bx1",gx+mgthr,fMGy,gz+fMGth);
  flogGrid = new G4LogicalVolume(bx1,fAl,"lGrid",0,0,0);
  flogGrid->SetVisAttributes(black);
  G4Box* bx2 = new G4Box("bx2",fMGx,fMGy,fMGz);
  G4Box* bx3[nVxSz];
  G4Box* bx4[4] = {NULL};
  G4LogicalVolume* logbx4[4] = {NULL};
  // If radial blade is different material(s) to Al Normal blades
  if(isScat){
    G4double sx = mgthr;   // radial blade thickness
    ScatterTest();
    for(G4int im=0; im<fNmat; im++){
      sx -= fScX[im];
      bx4[im] = new G4Box("bx2",sx,fMGy,gz+fMGth);
      logbx4[im] = new G4LogicalVolume(bx4[im],fScatMat[im],"radiallog",0,0,0);
      if(im==0)
	logbx4[im]->SetVisAttributes(black);
      else
	logbx4[im]->SetVisAttributes(black);
      if(im > 0)
	new G4PVPlacement(0,G4ThreeVector(0,0,0),logbx4[im],"RadPhys1",
			  logbx4[im-1],false,0);
    }
  }
  if(flogbx2[0] == NULL){
    for(G4int i=0; i<nVxSz; i++){
      flogbx2[i] = new G4LogicalVolume(bx2,fB4Ce,"b4clog",0,0,0);
      flogbx2[i]->SetVisAttributes(black);
      G4double z = fMGz - (fBt1[i] + fBt2[i])*0.5;
      G4double dz = -0.5*(fBt2[i] - fBt1[i]);
      bx3[i] = new G4Box("bx3",fMGx-fBthx,fMGy,z);
      flogArCO2[i] = new G4LogicalVolume(bx3[i],fArCO2,"ArCO2log",0,0,0);
      flogArCO2[i]->SetVisAttributes(black);
      new G4PVPlacement(0,G4ThreeVector(0,0,dz),flogArCO2[i],"ArCO2phys",
			flogbx2[i],false,1000);
    }
  }
  //
  G4double xoff = -gx + mgthr + fMGx;
  //+ 2*(fMGx + fMGth);
  G4int l,nn=1001;
  for(G4int i=0; i<fNx; i++){
    G4double zoff = -gz + fMGth + fMGz;
    //+ 2*(fMGz + fMGth);
    for(G4int j=0; j<fNz; j++){
      if(j < fn1) l=0;
      else if(j==fn1) l=1;
      else if(j<fn2) l=2;
      else if(j==fn2) l=3;
      else if(j<fn3) l=4;
      else if(j==fn3) l=5;
      else l=6;
      G4PVPlacement* bxp = 
	new G4PVPlacement(0,G4ThreeVector(xoff,0,zoff),flogbx2[l],"MGphys",
			  flogGrid,false,nn);
      G4ThreeVector xyz = -bxp->GetFrameTranslation();
      fprintf(CoordFile,"Cell: %d %g %g %g\n",nn,xyz.x(),xyz.y(),xyz.z());
      zoff += 2*(fMGz + fMGth);
      nn++;
    }
    if(isScat){
      new G4PVPlacement(0,G4ThreeVector(xoff-(fMGx+mgthr),0,0),logbx4[0],
			"Radphys",flogGrid,false,0);
      if(i == fNx-1)
	new G4PVPlacement(0,G4ThreeVector(xoff+(fMGx+mgthr),0,0),logbx4[0],
			  "Radphys",flogGrid,false,0);
    }
    xoff += 2*(fMGx + mgthr);
  }
  return;
}

//---------------------------------------------------------------------------
void DetectorConstruction::TRbox()
{
  // Hardwired model of single TREX box, containing 6 multigrids
  // Box(es) placed at end of procedure
  //
  G4double r0=3000*mm,r1,r3,r4,r5,dr0=18*mm;       // radii
  G4double wx1=1.5*mm,wx2,wy=1170*mm,wz=173*mm,wzi;// wedge struts wy was 1192
  G4double xi,zi,x1i,th,thmin,thmax,thr,dzt;
  G4double bx,by=1170*mm,bz=27.0*mm;              // back plate
  G4double bbx=52*mm,bby=1090*mm,bbz=10*mm;       // back bone
  G4double cax=52*mm,cay=75*mm,caz=16*mm,cz=10*mm;// carriage
  G4double pcay[] = {-900,-300,300,900};          // carriage y pos
  G4double rlx=30*mm,rly=1090*mm,rlz,rlz1;        // rail..slide MG
  G4double gz = fNz*(fMGz + fMGth)+ fMBth2;       // length of MG
  G4double tp1=40*mm,tp2=30*mm,tp3=20*mm;         // for top/bot plate
  fGrz = gz+fMGth;                                // length of MG
  r1 = r0+2*wx1;                                  //
  r4 = r0+2*wz+bz;
  //  thmax = 3*fdTh;                                 // angle limits of box
  thmax = 8.85*deg;
  thmin = 90*deg-thmax;
  thr = 2*thmax;                                   // max angle
  dzt = 0.5*(r4-r0*cos(thmax));                    // max extent in z
  bx = xi = r4*tan(thmax);
  x1i = r0*sin(thmax);
  // trapezoid has outline of box and is filled with Ar/CO2
  // box will be placed in this mother volume
  auto mggas = new G4Trd("mggas",x1i,xi,by+2*tp2,by+2*tp2,dzt);
  auto lmggas = new G4LogicalVolume(mggas,fArCO2,"lmggas",0,0,0);
  lmggas->SetVisAttributes(cyan);
  // top/bottom plate
  auto tbp0 = new G4Trd("tbp0",x1i,xi,tp2,tp2,dzt);
  auto tbp1 = new G4Trd("tbp1",x1i-tp1,xi-tp1,tp2,tp2,dzt-tp1);
  auto tbp=new G4SubtractionSolid("tbp",tbp0,tbp1,0,G4ThreeVector(0,tp3,0));
  auto ltbp = new G4LogicalVolume(tbp,fAl,"ltbp",0,0,0);
  ltbp->SetVisAttributes(blue);
  // Start building box as a multi-union shape
  auto FrontP = new G4Tubs("FrontP",r0,r1,wy,thmin,thr); // front window
  auto BackP = new G4Box("BackP",bx,by,bz);              // back plate
  wx2 = wx1 + wz*tan(fdTh);
  auto strut = new G4Trd("strut",wx1,wx2,wy,wy,wz);      // central strut
  auto trbox = new G4MultiUnion("trbox");
  auto pos = G4ThreeVector(0.,0.,0.);
  auto rotm = G4RotationMatrix();
  auto tr3d = G4Transform3D(rotm,pos);
  trbox->AddNode(*BackP,tr3d);            // back plate
  pos.set(0,0,-wz);
  tr3d = G4Transform3D(rotm,pos);
  trbox->AddNode(*strut,tr3d);            // central strut
  pos.set(0,0,-r0-2*wz);
  rotm.rotateX(90*deg);
  tr3d = G4Transform3D(rotm,pos);
  trbox->AddNode(*FrontP,tr3d);           // front window
  rotm.rotateX(-90*deg);
  // Add the other 6 wedge struts
  G4Trd* struti[3];
  auto rotm1 = G4RotationMatrix();
  for(G4int i=1; i<=3; i++){
    th = i*fdTh;  // bx = (r4-bz)*sin(thmax);
    rotm.rotateY(fdTh);
    rotm1.rotateY(-fdTh);
    r5 = (r0+2*wz)/cos(th);
    wzi = 0.5*(r5-r0);
    r3 = r0+wzi;
    xi = r3*sin(th);
    zi = -wzi*cos(th);
    auto ppi = G4ThreeVector(xi,0,zi);
    auto tri = G4Transform3D(rotm,ppi);
    //if(i==3)
    //  dt *= 0.5;
    struti[i] = new G4Trd("struti",wx1,wx1+wzi*tan(fdTh),wy,wy,wzi);
    trbox->AddNode(*struti[i],tri);
    ppi.set(-xi,0,zi);
    auto tri1 = G4Transform3D(rotm1,ppi);
    trbox->AddNode(*struti[i],tri1);
  }
  //
  trbox->Voxelize();
  rotm = G4RotationMatrix();
  pos.set(0.,0.,-dzt+bz);
  tr3d = G4Transform3D(rotm,pos);
  auto trboxa = new G4IntersectionSolid("trboxa",trbox,mggas,tr3d);
  //auto trboxb = trboxa->CreatePolyhedron();
  auto ltrbox = new G4LogicalVolume(trboxa,fAl,"ltrbox",0,0,0);
  //auto ltrbox = new G4LogicalVolume(trbox,fAl,"ltrbox",0,0,0);
  ltrbox->SetVisAttributes(blue);
  new G4PVPlacement(0,G4ThreeVector(0,0,+dzt-bz),ltrbox,
		      "ptrbox",lmggas,false,0);
  // Add MG mounting assy on box
  auto car1 = new G4Box("car1",cax,cay,caz);
  auto car2 = new G4Box("car2",cax-15,cay,caz);
  auto cars = new G4SubtractionSolid("cars",car1,car2,0,G4ThreeVector(0,0,cz));
  auto lcar1 = new G4LogicalVolume(cars,fAl,"lcars",0,0,0);
  //auto rail1 = new G4Box("rail1",rlx*2,rly+5,rlz1);
  lcar1->SetVisAttributes(red);
  auto bb = new G4Box("backbone",bbx,bby,bbz);
  auto lbb = new G4LogicalVolume(bb,fAl,"lbackbone",0,0,0);
  lbb->SetVisAttributes(blue);
  G4RotationMatrix rotm3 = G4RotationMatrix();
  G4RotationMatrix rotm4 = G4RotationMatrix();
  th = -2.5*fdTh;
  rotm4.rotateY(th);
  G4LogicalVolume* lrail[6];
  char nme[64];
  for(int i=0; i<6; i++){
    // Backbone
    r5 = r0 + dr0 + 2*fGrz + bbz;
    xi = r5*sin(th);
    zi = -r4+r5*cos(th)+dzt;
    auto roti = new G4RotationMatrix();
    roti->rotateY(-th);
    sprintf(nme,"pBackBone%d",i);
    new G4PVPlacement(roti,G4ThreeVector(xi,0,zi),lbb,
		      nme,lmggas,false,0);
    // Slider
    r5 += caz + bbz;
    xi = r5*sin(th);
    zi = -r4+r5*cos(th)+dzt;
    for(int j=0; j<4; j++){
      sprintf(nme,"pCar%d%d",i,j);
      new G4PVPlacement(roti,G4ThreeVector(xi,pcay[j],zi),lcar1,
		      nme,lmggas,false,0);
    }
    // Rail
    rlz = 0.5*(dzt-2*bz-zi);
    rlz1 = rlz - rlx*tan(th);
    auto rail = new G4Trap("rail",2*rly,2*rlx,2*rlz,2*rlz1);
    lrail[i] = new G4LogicalVolume(rail,fAl,"lrail",0,0,0);
    lrail[i]->SetVisAttributes(blue);
    auto rotr = new G4RotationMatrix();
    rotr->rotateY(-90*deg);
    rotr->rotateX(90*deg);
    rlz1 = 0.5*(rlz+rlz1);
    r5 = r0 + dr0 + 2*fGrz + 2*bbz + rlz;
    xi = r5*sin(th);
    zi = dzt - 2*bz - rlz1;
    sprintf(nme,"pRail%d",i);
    new G4PVPlacement(rotr,G4ThreeVector(xi,0,zi),lrail[i],
		      nme,lmggas,false,0);
    // Add MultiGrid
    r5 = r0 + dr0 + fGrz;
    xi = r5*sin(th);
    zi = r4-r5*cos(th)-dzt;
    MultiGrid(0,lmggas,xi,-zi,-th);
    rotm4.rotateY(fdTh);
    th += fdTh;
  }
  //Place top/bottom plate
  new G4PVPlacement(0,G4ThreeVector(0,by+tp2,0),ltbp,
		    "ptop",lmggas,false,0);
  auto rotb = new G4RotationMatrix();
  rotb->rotateZ(180*deg);
  new G4PVPlacement(rotb,G4ThreeVector(0,-(by+tp2),0),ltbp,
		    "ptop",lmggas,false,0);
  // Place TREX boxes
  r5 = r0*cos(thmax)+dzt;
  th = fThSt + thmax;
  for(int n=0; n<fNvess; n++){
    th = fThBox[n];
    auto rottr = new G4RotationMatrix();
    rottr->rotateY(-th);
    xi = r5*sin(th);
    zi = r5*cos(th);
    new G4PVPlacement(rottr,G4ThreeVector(xi,0,zi),lmggas,
		    "ptrgas",fExpLog,false,n+10);
    //th += 2*thmax;
  }
  // Optional shielding vanes between boxes
  if(fShTh > 0){
    G4double vnz = 215;
    G4int nsw;
    auto shbx = new G4Box("shbx",fShTh,by,vnz);
    auto lshbx = new G4LogicalVolume(shbx,fScMatList[fShMat-1],"lshbx",0,0,0);
    lshbx->SetVisAttributes(green);
    //r5 += dr0+10;
    r5 -= vnz + 120;
    G4double dth = 0.5*(fThBox[1]-fThBox[0]);
    for(int n=0; n<fNvess+2; n++){
      if(n<fNvess-1){
	if((fThBox[n+1]<0)&&(fThBox[n]>0)){  // switch sides?
	  th = fThBox[n] + dth;
	  nsw = n+1;
	}
	else
	  th = 0.5*(fThBox[n+1]+fThBox[n]);
      }
      else if(n == fNvess-1)
	th = fThBox[fNvess-1] - dth;   //end
      else if(n == fNvess)
	th = fThBox[nsw] + dth;
      else if(n==fNvess+1)
	th = fThBox[0] - dth;
      auto rottr = new G4RotationMatrix();
      rottr->rotateY(-th);
      xi = r5*sin(th);
      zi = r5*cos(th);
      new G4PVPlacement(rottr,G4ThreeVector(xi,0,zi),lshbx,
			"pshbx",fExpLog,false,0);
    } 
  }
  // Make and place optional target
  Target();
  // Make and place optional cryostat around the target
  CryoStat();
}

//---------------------------------------------------------------------------
void DetectorConstruction::Trex()
{
  G4double sx = fNx * (fMGx + fMGthr);
  G4double sy = fNy * (fMGy + fMGsp);
  G4Box* trb1 = new G4Box("trb1",sx,sy,fTRz1);
  G4Box* trb1a = new G4Box("trb1a",sx-fMGthr,sy,fTRz1-fMGth);
  G4RotationMatrix* rot = new G4RotationMatrix();
  G4SubtractionSolid* trb1b =
    new G4SubtractionSolid("trb1b",trb1,trb1a,rot,G4ThreeVector(0,0,-fMGth));
  flogTR1 = new G4LogicalVolume(trb1b,fAl,"ltrb1",0,0,0);
  G4Box* trb2 = new G4Box("trb2",fTRx2,sy,fTRz2);
  flogTR2 = new G4LogicalVolume(trb2,fAl,"ltrb2",0,0,0);
  G4Box* trb3 = new G4Box("trb3",fTRx3,sy,fTRz3);
  flogTR3 = new G4LogicalVolume(trb3,fAl,"ltrb3",0,0,0);
  flogTR1->SetVisAttributes(blue);
  flogTR2->SetVisAttributes(blue);
  flogTR3->SetVisAttributes(blue);
}
  
//---------------------------------------------------------------------------
void DetectorConstruction::Shield(G4int nmg)
{
  G4Material* mb = fScMatList[fMBmat[nmg]-1];
  G4Material* mb1 = fScMatList[fMBmat1[nmg]-1];
  G4Material* mb2 = fScMatList[fMBmat2[nmg]-1];
  G4double mgthr;
  if(nmg==0)
    mgthr = fMGthr;
  else
    mgthr = fMGthr1;
  // Back shield
  if(fMBth){
    G4double sx = fNx * (fMGx + mgthr) - mgthr;
    //G4double sy = fNy * (fMGy + fMGsp);
    //fMposz = fNz*(fMGz + fMGth) + fMBth + fMGth;
    fMposz = fNz*(fMGz + fMGth) + fMBth2 + fMGth;
    G4Box* box = new G4Box("Shield",sx,fMGy,fMBth);
    G4LogicalVolume* logMB = new G4LogicalVolume(box,mb,"lShield",0,0,0);
    logMB->SetVisAttributes(purple);
    sx += 2*mgthr;
    G4Box* box2 = new G4Box("RBlock",sx,fMGy,fMBth2);
    flogMB = new G4LogicalVolume(box2,mb2,"lRBlock",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,fMBth-fMBth2),logMB,
		      "pShield",flogMB,false,10000);
    flogMB->SetVisAttributes(black);
  }
  // side shield
  if(fMBth1){
    G4double sz = fNz * (fMGz + fMGth);
    fMposx = fNx*(fMGx + mgthr) + fMBth1 + mgthr;
    G4Box* box1 = new G4Box("Shield1",fMBth1,fMGy,sz);
    flogMB1 = new G4LogicalVolume(box1,mb1,"lShield1",0,0,0);
    flogMB1->SetVisAttributes(purple);
  }
}
//---------------------------------------------------------------------------
void DetectorConstruction::ScatterTest()
{
  // not incorporated in multigrid or He3 tube array?
  if((fNy == 0)&&(fNHe3 == 0)){
    // 4-pi sphere detector?
    if(fScDetSx == 0){
      G4Sphere* sph = new G4Sphere("ScSphere",fScDetSz*mm,(fScDetSz+25)*mm,
				   0.0*deg,360*deg,2.5*deg,175*deg);
      flogScDet = new G4LogicalVolume(sph,fScMatList[fScDetMat-1],
				      "lScDet",0,0,0);
      flogScDet->SetVisAttributes(black);
    }
    else{
      G4double szAl = 0.5;
      G4double szAr = 2.0;
      auto scboxA = new G4Box("ScBoxAl",fScDetSx+1,fScDetSy+1,fScDetSz+1);
      auto scboxB = new G4Box("ScBoxB",fScDetSx,fScDetSy,fScDetSz);
      //auto scboxD = new G4Box("ScBoxAr",fScDetSx,fScDetSy,szAr);
      auto logScDetA = new G4LogicalVolume(scboxA,fAl,"lScDetAl",0,0,0);
      auto logScDetB = new G4LogicalVolume(scboxB,fScMatList[fScDetMat-1],
					   "lScDetB",0,0,0);
      flogScDet = new G4LogicalVolume(scboxA,fArCO2,"lScDetAr",0,0,0);
      G4double dz = fScDetZ - (szAl + szAr + 2*fScDetSz);
      new G4PVPlacement(0,G4ThreeVector(fScDetX,fScDetY,dz),logScDetA,
		      "pScDetAl",fExpLog,false,0);
      G4double dz1 = fScDetZ - (szAr + fScDetSz);
      new G4PVPlacement(0,G4ThreeVector(fScDetX,fScDetY,dz1),logScDetB,
		      "pScDetB",fExpLog,false,0);
      logScDetA->SetVisAttributes(black); 
      logScDetB->SetVisAttributes(blue); 
      flogScDet->SetVisAttributes(cyan);
    }
    new G4PVPlacement(0,G4ThreeVector(fScDetX,fScDetY,fScDetZ),flogScDet,
		      "pScDet",fExpLog,false,30);
  }
  G4RotationMatrix* yR = new G4RotationMatrix;
  yR->rotateY(fScTh[0]*deg);
  // ArCO2?  put scatterer in detector
  if(fScDetSx > 0)
    flogScat[0] = flogScDet;
  else
    flogScat[0] = fExpLog;
  char name[16];
  for(G4int im=0; im<fNmat; im++){
    fScatMat[im] = fScMatList[fScMat[im]-1];
    sprintf(name,"lScat%d",im);
    if(fNy > 0)
      continue;                    // Trex simulation, skip scatter sample
    G4Box* box = new G4Box("Scat",fScX[im],fScY[im],fScZ[im]);
    flogScat[im+1] = new G4LogicalVolume(box,fScatMat[im],name,0,0,0);
    flogScat[im+1]->SetVisAttributes(black);
    sprintf(name,"pScat%d",im);
    new G4PVPlacement(yR,G4ThreeVector(fScPx,fScPy,fScPz),flogScat[im+1],name,
				flogScat[im],false,20+im);
    yR = NULL;
  }
}
//----------------------------------------------------------------
void DetectorConstruction::Collimator()
{
  G4LogicalVolume* lcbx0=NULL;
  G4LogicalVolume* lcbx0p=NULL;
  G4double lconc;
  if(fWallTh){
    auto cbx0 = new G4Box("ColliBox0",fExpX,fExpY,fWallTh);
    auto cbx0a = new G4Tubs("ColliBox0a",0,fColSx,fWallTh,0,360*deg);
    auto cbx0s = new G4SubtractionSolid("Wall",cbx0,cbx0a);
    auto cbx0p = new G4Tubs("Pipe",fColSy,fColSx,2*fWallTh,0,360*deg);
    lcbx0 = new G4LogicalVolume(cbx0s,fConcrete,"lcbx0",0,0,0);
    lcbx0->SetVisAttributes(black);
    lcbx0p = new G4LogicalVolume(cbx0p,fAl,"lcbx0p",0,0,0);
    lcbx0p->SetVisAttributes(black);
  }
  auto cbx1 = new G4Box("ColliBox1",fColSx,fColSx,fColSz);
  auto cbx1a = new G4Box("ColliBox1a",fColHx,fColHy,fColSz);
  auto cbx1s = new G4SubtractionSolid("ColliBox",cbx1,cbx1a);
  auto lcbx1 = new G4LogicalVolume(cbx1s,fScMatList[fColMat-1],
					      "lcbx1",0,0,0);
  lcbx1->SetVisAttributes(black);
  //G4double zoff = fColZ - fWallTh + fColSz;
  G4double zoff = fColZ - 2*fWallTh - fColSz;
  if(lcbx0){
    new G4PVPlacement(0,G4ThreeVector(0,0,zoff),lcbx0,
		      "pcbx0",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,zoff),lcbx0p,
		      "pcbx0p",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,fWallZ1),lcbx0,
			"pcbx0a",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,fWallZ1),lcbx0p,
		      "pcbx0pa",fExpLog,false,0);
    //lconc = 0.5*(fWallZ1 - fColZ-fWallTh-fColSz);
    lconc = 0.5*(fWallZ1 - zoff - 2*fWallTh);
  }
  //new G4PVPlacement(0,G4ThreeVector(0,0,fColZ+fWallTh+fColSz*2),lcbx1,
  //		      "pcbx1",fExpLog,false,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,fColZ),lcbx1,
		      "pcbx1",fExpLog,false,0);
  if(fWallTh){
    G4double offs = 500;
    G4double xcub = fExpX - offs;
    //G4double ycub = fExpY - offs;
    G4Box* cbx3 = new G4Box("ConcBox0",xcub,fExpY,lconc);
    //G4Box* cbx3 = new G4Box("ConcBox0",fExpX,ycub,lconc);
    G4Box* cbx3a = new G4Box("ConcBox0a",xcub-2.0*fWallTh,
			     fExpY-2.0*fWallTh,lconc);
    //G4Box* cbx3a = new G4Box("ConcBox0a",fExpX-2.0*fWallTh,
    //ycub-2.0*fWallTh,lconc);
    G4SubtractionSolid* cbx3s = new G4SubtractionSolid("Wall2",cbx3,cbx3a);
    G4LogicalVolume* lcbx3 = new G4LogicalVolume(cbx3s,fConcrete,"lcbx3",0,0,0);
    lcbx3->SetVisAttributes(black);
    //G4double zz = fColZ + fColSz + lconc;
    G4double zz = fColZ + lconc -fWallTh - fColSz;
    new G4PVPlacement(0,G4ThreeVector(offs,0,zz),lcbx3,
    		      "pcbx3",fExpLog,false,0);
    //new G4PVPlacement(0,G4ThreeVector(0,-offs,zz),lcbx3,
    //		      "pcbx3",fExpLog,false,0);
  }
}
//----------------------------------------------------------------
void DetectorConstruction::BeamMonitor()
{
  // In beam monitor
  auto xR = new G4RotationMatrix;
  xR->rotateY(90.0*deg);
  //xR->rotateX(90.0*deg);
  //auto bmcase = new G4Tubs("BMcase",fBMsx+1,fBMsx+2,fBMsy+5,0.0,360.0*deg);
  auto bmcase = new G4Tubs("BMcase",fBMsx+1,fBMsx+1.5,fBMsy+5,0.0,360.0*deg);
  auto bmbx = new G4Box("BMBox1",fBMsx,fBMsy,fBMsz);
  flogBM = new G4LogicalVolume(bmbx,fScMatList[fBMmat-1],
			       "lbmbx",0,0,0);
  auto logBMcase = new G4LogicalVolume(bmcase,fAl,"lbmcase",0,0,0);
  flogBM->SetVisAttributes(black);
  logBMcase->SetVisAttributes(black);
  new G4PVPlacement(0,G4ThreeVector(fBMx,fBMy,fBMz),flogBM,
		      "pbmbx",fExpLog,false,40);
  new G4PVPlacement(xR,G4ThreeVector(fBMx,fBMy,fBMz),logBMcase,
		      "pbmcase",fExpLog,false,0);
}
//----------------------------------------------------------------
void DetectorConstruction::EMMAshield()
{
  // extra shielding installed at EMMA
  //fTabSx = 20;
  //fTabSy = 500;
  //fTabSz = 500;
  //fTabX = 1400;
  //fTabY = 0;
  //fTabZ = 250;
  G4double ha=670,hb=630,hc=730,hd=790,wab=382.5,thab=5.0;
  G4double hp=800,wpo=35.0,wpi=20.0;
  G4double xa = fTabX - ha - fTabSx;
  G4double xb = fTabX - hb - fTabSx;
  G4double xc = fTabX - hc - fTabSx;
  G4double xd = fTabX - hd - fTabSx;
  G4double xp = fTabX - hp - fTabSx;
  G4double xe = fTabX - fTabSx - 100;
  G4double d0 = wab + thab;
  G4double d1 = wab + wpo;
  // Table top
  auto tab = new G4Box("Table",fTabSx,fTabSy,fTabSz);
  auto ltab = new G4LogicalVolume(tab,fAl,"ltab",0,0,0);
  ltab->SetVisAttributes(black);
  new G4PVPlacement(0,G4ThreeVector(fTabX,fTabY,fTabZ),ltab,
		      "pTab",fExpLog,false,0);
  if(fShL1 > 0){
    // Al Pillar
    auto bxO = new G4Box("PilO",hp,wpo,wpo);
    auto bxI = new G4Box("PilI",hp,wpi,wpi);
    auto lbxO = new G4LogicalVolume(bxO,fAl,"lbxO",0,0,0);
    auto lbxI = new G4LogicalVolume(bxI,fAir,"lbxI",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),lbxI,"pbxI",lbxO,false,0);
    lbxO->SetVisAttributes(black);
    lbxI->SetVisAttributes(black);
    //
    auto bxA = new G4Box("ShA",ha,thab,wab);
    auto bxB = new G4Box("ShB",hb,wab,thab);
    auto bxC = new G4Box("ShC",hc,thab,wab);
    auto bxD = new G4Box("ShD",hd,wab,thab);
    auto bxE = new G4Box("ShE",thab,wab,wab);
    auto lbxA = new G4LogicalVolume(bxA,fMB,"lbxA",0,0,0);
    auto lbxB = new G4LogicalVolume(bxB,fMB,"lbxB",0,0,0);
    auto lbxC = new G4LogicalVolume(bxC,fMB,"lbxC",0,0,0);
    auto lbxD = new G4LogicalVolume(bxD,fMB,"lbxD",0,0,0);
    auto lbxE = new G4LogicalVolume(bxE,fMB,"lbxE",0,0,0);
    lbxA->SetVisAttributes(black);
    lbxB->SetVisAttributes(black);
    lbxC->SetVisAttributes(black);
    lbxD->SetVisAttributes(black);
    lbxE->SetVisAttributes(black);
    new G4PVPlacement(0,G4ThreeVector(xa,-d0+fTabY,fTabZ),lbxA,"pbxA",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xb,fTabY,fTabZ-d0),lbxB,"pbxB",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xc,d0+fTabY,fTabZ),lbxC,"pbxC",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xd,fTabY,fTabZ+d0),lbxD,"pbxD",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xe,fTabY,fTabZ),lbxE,"pbxE",fExpLog,false,0);
    //
    new G4PVPlacement(0,G4ThreeVector(xp,-d1+fTabY,fTabZ+d1),lbxO,"pbxO",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xp,-d1+fTabY,fTabZ-d1),lbxO,"pbx1",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xp,d1+fTabY,fTabZ+d1),lbxO,"pbx2",fExpLog,false,0);
    new G4PVPlacement(0,G4ThreeVector(xp,d1+fTabY,fTabZ-d1),lbxO,"pbx3",fExpLog,false,0);
  }
  G4double shsx,shsy,shsz;
  if(fShL0 > 0){
    shsz = fShL0*fabs(fBMz);
    //shsx = fColHx+5;
    //shsy = fColHy+5;
    shsx = 25.0;
    shsy = 26.0;
    //G4Box* shbm0 = new G4Box("ShBox0",shsx+5,shsy+5,shsz);
    //G4Box* shbm0a = new G4Box("ShBox0a",shsx,shsy,shsz);
    //G4SubtractionSolid* shbm0s = new G4SubtractionSolid("EMSh0",shbm0,shbm0a);
    //G4LogicalVolume* lshbm0 = new G4LogicalVolume(shbm0s,fB4C,"lshbm0",0,0,0);
    auto shbm0 = new G4Tubs("EMSh0",shsx,shsy,shsz,0.0,360.0*deg);
    auto lshbm0 = new G4LogicalVolume(shbm0,fCd,"lshbm0",0,0,0);
    lshbm0->SetVisAttributes(black);
    new G4PVPlacement(0,G4ThreeVector(0,0,0.5*fBMz),lshbm0,
		      "BmShield",fExpLog,false,0);
  }
  if(fShL2 > 0){
    auto yR = new G4RotationMatrix;
    yR->rotateX(90*deg);
    yR->rotateY(180*deg);
    yR->rotateZ(180*deg);
    //shsz = fShL2*fabs(fBMz);
    //shsx = 100.0;
    //shsy = 2.5;
    shsy = 2*fShL2*fabs(fBMz);
    //shsx = 5.0;
    //shsz = 200.0;
    //auto shbm2 = new G4Box("ShBox2",shsx,shsy,shsz);
    auto shbm2 = new G4Trap("ShBox2",fSh2sx,shsy,fSh2sz,fSh2sz1);
    auto lshbm2 = new G4LogicalVolume(shbm2,fMB,"lshbm2",0,0,0);
    lshbm2->SetVisAttributes(black);
    new G4PVPlacement(yR,G4ThreeVector(fSh2z,-fSh2x,0.5*fBMz),lshbm2,
		      "BmShield2a",fExpLog,false,0);
    new G4PVPlacement(yR,G4ThreeVector(fSh2z,fSh2x,0.5*fBMz),lshbm2,
		      "BmShield2b",fExpLog,false,0);
  }
  /*
  if(fShL1 > 0){
    if(fShZ > 0)
      shsz = fShZ;
    else
      shsz = fVx+fVth+20;
    shsx = fTabX*fShL1;
    if(fShY > 0)
      shsy = fShY;
    else
      shsy = fVy+fVth+20;
    G4Subdouble shpx = fTabX -fTabSx - shsx;
    G4Box* shbm1 = new G4Box("ShBox1",shsx,shsy+50,shsz+50);
    G4Box* shbm1a = new G4Box("ShBox1a",shsx,shsy,shsz);
    G4SubtractionSolid* shbm1s = new G4SubtractionSolid("EMSh1",shbm1,shbm1a);
    G4LogicalVolume* lshbm1 = new G4LogicalVolume(shbm1s,fB4C,"lshbm1",0,0,0);
    lshbm1->SetVisAttributes(black);
    new G4PVPlacement(0,G4ThreeVector(shpx,0,0),lshbm1,
		      "BmShield1",fExpLog,false,0);
  }
  */
}
//----------------------------------------------
void DetectorConstruction::EMMAHe3()
{
  G4double wR = 0.025*mm;
  auto tube = new G4Tubs("he3Tube",0.0,fT3ro,fT3l,0.0,360.0*deg);
  auto he3Gas = new G4Tubs("he3Gas",0.0,fT3ri,fT3l,0.0,360.0*deg);
  auto wire = new G4Tubs("wire",0.0,wR,fT3l,0.0,360.0*deg);  
  flogHe3 = new G4LogicalVolume(he3Gas,fHe3Ar,"logHe3",0,0,0);
  flogT3 = new G4LogicalVolume(tube,fFe,"logT3",0,0,0);
  auto logWire = new G4LogicalVolume(wire,fW,"logWire",0,0,0);
  flogT3->SetVisAttributes(black);
  flogHe3->SetVisAttributes(blue);
  logWire->SetVisAttributes(black);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),logWire,"wireExp",flogHe3,false,999);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),flogHe3,"gasExp",flogT3,false,1001);
  auto xR = new G4RotationMatrix;
  G4double y = fT3y - (fnT3-1)*fT3ro;
  for(G4int n=0; n<fnT3; n++){
    xR->rotateX(90*deg);
    auto phe3 = new G4PVPlacement(xR,G4ThreeVector(fT3x,y,fT3z),flogT3,"vT3",fExpLog,false,n+1440);
    y += fT3ro*2;
    G4ThreeVector xyz = -phe3->GetFrameTranslation();
    fprintf(CoordFile,"He3: %d %g %g %g\n",n,xyz.x(),xyz.y(),xyz.z());
  }
}
//----------------------------------------------
void DetectorConstruction::HSpacer()
{
  // H-Spacers between grids
  //auto yR = new G4RotationMatrix;
  //yR->rotateX(90*deg);
  auto bx1 = new G4Box("bx1",fHsx,fHsy,fHsz);
  auto bx2 = new G4Box("bx2",fHsx1,fHsy1,fHsz);
  G4double dy = fHsy - fHsy1;
  G4double dx = fHsx - fHsx1;
  auto hbx = new G4SubtractionSolid("hbx",bx1,bx2,0,G4ThreeVector(-dx,dy,0));
  auto hbx1 = new G4SubtractionSolid("hbx1",hbx,bx2,0,G4ThreeVector(-dx,-dy,0));
  flogHsp = new G4LogicalVolume(hbx1,fAl2O3,"logHsp",0,0,0);
  flogHsp->SetVisAttributes(red);
}

//----------------------------------------------------------------
void DetectorConstruction::Materials()
{
  // Define materials of Sonde detector and surroundings
  //
  char nccfg[256];
  G4NistManager* nistM = G4NistManager::Instance();
  // Vaccum
  fVac = nistM->FindOrBuildMaterial("G4_Galactic");
  // Air
  fAir = nistM->FindOrBuildMaterial("G4_AIR");
  // Water
  fConcrete = nistM->FindOrBuildMaterial("G4_CONCRETE");
  if(!fIsCrysMat){
    fLucite = nistM->FindOrBuildMaterial("G4_LUCITE");
    fWater = nistM->FindOrBuildMaterial("G4_WATER");
  }
  // ArCO2 80/20%
  G4Material* CO2  = nistM->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  //G4Material* gAr = new G4Material("gAr",18.,39.95*g/mole,1.784*mg/cm3);
  G4Material* gAr = nistM->FindOrBuildMaterial("G4_Ar");
  fArCO2 = new G4Material("Ar-CO2",fArCO2dens,2);
  fArCO2->AddMaterial(CO2,fCO2Frac); 
  fArCO2->AddMaterial(gAr,fArFrac);
  G4cout << fArCO2 << G4endl;
  //
  G4Element* H  = nistM->FindOrBuildElement("H");
  G4Element* O  = nistM->FindOrBuildElement("O");
  G4Element* C  = nistM->FindOrBuildElement("C");
  G4Element* enHe = new G4Element("enHe","enHe",2);
  enHe->AddIsotope(new G4Isotope("He3", 2, 3, 3.0160*g/mole), 0.995);
  enHe->AddIsotope(new G4Isotope("He4", 2, 4, 4.0026*g/mole), 0.005);
  //fHe3CO2 = new G4Material("He3CO2", fHe3pr*0.495*mg/cm3,3);
  //fHe3CO2 = new G4Material("He3CO2", fHe3pr*0.219*mg/cm3,3);
  fHe3CO2 = new G4Material("He3CO2", fHe3pr*0.18015*mg/cm3,3);
  //fHe3CO2->AddElement(enHe,4);
  // fHe3CO2->AddElement(C,1);
  // fHe3CO2->AddElement(O,2);
  //fHe3CO2->AddElement(enHe,57);
  //fHe3CO2->AddElement(C,1);
  //fHe3CO2->AddElement(O,2);
  fHe3CO2->AddElement(enHe,117);
  fHe3CO2->AddElement(C,1);
  fHe3CO2->AddElement(O,2);
  G4cout << fHe3CO2 << G4endl;
  // 10 bar He3 + 3 bar Ar for EMMA test He3 array
  G4double rHe3Ar = (10.0*rHe3 + 3.0*rAr)*g/cm3;
  fHe3Ar = new G4Material("He3Ar",rHe3Ar,2);
  fHe3Ar->AddElement(enHe,275);
  G4Element* eAr  = nistM->FindOrBuildElement("Ar");
  fHe3Ar->AddElement(eAr,67);
  G4cout << fHe3Ar << G4endl;
  //
  // Boron Carbide converter and shielding
  //G4Element* enB = new G4Element("enB", "enB", 2);
  //enB->AddIsotope(new G4Isotope("B10", 5, 10, 10.012*g/mole), fB10frac);
  G4double b11frac = 1.0 - fB10frac;
  //enB->AddIsotope(new G4Isotope("B11", 5, 11, 11.009*g/mole), b11frac);
  //scale density according to b10/b11 ratio
  G4double rb4c = rB4C*(0.2 + 0.8*(fB10frac*10.01 + b11frac*11.01)/10.811);
  //fB4Ce = new G4Material("B4Ce", rb4c*g/cm3 ,2);
  //fB4Ce->AddElement(enB,4);
  //fB4Ce->AddElement(C,1);
  // NCrystal B4C
  sprintf(nccfg,"B4C_sg166_BoronCarbide.ncmat;atomdb=B:is:%g:B10:%g:B11",fB10frac,1-fB10frac);
  fB4Ce =  G4NCrystal::createMaterial(nccfg);
  G4cout << fB4Ce << std::endl;
  // pure boron-10
  G4Element* b10 = new G4Element("b10", "b10", 1);
  b10->AddIsotope(new G4Isotope("B10", 5, 10, 10.012*g/mole), 1.0);
  fB10 = new G4Material("B10pure", 2.5*g/cm3 ,1);
  fB10->AddElement(b10,1);
  G4cout << fB10 << std::endl;
  // natural B4C
  fB4C = nistM->FindOrBuildMaterial("G4_BORON_CARBIDE");
  G4cout << fB4C << G4endl;
  // B4C in epoxy matrix
  G4double vmix = W2V(fDripMix,rb4c,rEpoxy);    // volume ratios
  //G4double rmb = vmix*rb4c + (1-vmix)*rEpoxy;
  G4Material* glue = new G4Material("glue",rEpoxy*g/cm3,3);
  glue->AddElement(C,11);
  glue->AddElement(O,3);
  glue->AddElement(H,12);
  fMB = new G4Material("Mirrobor",rMb*g/cm3,2);
  fMB->AddMaterial(fB4C,fDripMix);
  fMB->AddMaterial(glue,(1.0-fDripMix));
  G4cout << fMB << std::endl;
  // NCrystal dripped mixture
  sprintf(nccfg,"phases<%g*B4C_sg166_BoronCarbide.ncmat&%g*Epoxy_Araldite506_C18H20O3.ncmat>;atomdb=B:is:%g:B10:%g:B11",vmix,1-vmix,fB10frac,1-fB10frac);
  fDrip =  G4NCrystal::createMaterial(nccfg);
  G4cout << fDrip << std::endl;
  // Al,Fe etc. poly-crystaline
  if(fIsCrysMat){
    fAl = G4NCrystal::createMaterial("Al_sg225.ncmat");
    fFe = G4NCrystal::createMaterial("Fe_sg229_Iron-alpha.ncmat");
    fV = G4NCrystal::createMaterial("V_sg229.ncmat");
    fGraf = G4NCrystal::createMaterial("C_sg194_pyrolytic_graphite.ncmat");
    fCH2 = G4NCrystal::createMaterial("Polyethylene_CH2.ncmat");
    fCr  = G4NCrystal::createMaterial("Cr_sg229.ncmat");
    fNi  = G4NCrystal::createMaterial("Ni_sg225.ncmat");
    fB4Cnc = G4NCrystal::createMaterial("B4C_sg166_BoronCarbide.ncmat");
    fLucite = G4NCrystal::createMaterial("AcrylicGlass_C5O2H8.ncmat");
    fWater = G4NCrystal::createMaterial("LiquidWaterH2O_T293.6K.ncmat");
    fPb = G4NCrystal::createMaterial("Pb_sg225.ncmat");
    fW = G4NCrystal::createMaterial("W_sg229.ncmat");
    fC8H8 = G4NCrystal::createMaterial("Polystyrene_C8H8.ncmat");
  }
  // Al,Fe etc. amorphous
  else{
    fAl = nistM->FindOrBuildMaterial("G4_Al");
    fFe = nistM->FindOrBuildMaterial("G4_Fe");
    fV = nistM->FindOrBuildMaterial("G4_V");
    fCr = nistM->FindOrBuildMaterial("G4_Cr");
    fNi = nistM->FindOrBuildMaterial("G4_Ni");
  }
  // Nat-Boron-Al mixture...solid sheets
  vmix = W2V(fAlMix,rAl,rB4C);
  G4double rbal = vmix*rAl + (1-vmix)*rB4C;
  fBorAl = new G4Material("BorAl",rbal*g/cm3,2);
  fBorAl->AddMaterial(fB4C,(1-fAlMix));
  fBorAl->AddMaterial(fAl,fAlMix);
  G4cout << fBorAl << std::endl;
  // NCrystal BorAl
  sprintf(nccfg,"phases<%g*Al_sg225.ncmat&%g*B4C_sg166_BoronCarbide.ncmat>;density=%ggcm3",vmix,1-vmix,rbal);
  fBorAlnc =  G4NCrystal::createMaterial(nccfg);
  G4cout << fBorAlnc << std::endl;
  // Nat-Boron-Al mixture...spray mixture
  vmix = W2V(fAlMixS,rAl,rB4C);
  rbal = vmix*rAl + (1-vmix)*rB4C;
  fBorAlS = new G4Material("BorAlS",rbal*g/cm3,2);
  fBorAlS->AddMaterial(fB4C,(1-fAlMixS));
  fBorAlS->AddMaterial(fAl,fAlMixS);
  G4cout << fBorAlS << std::endl;
  // Nat-Boron-Ni spray mixture
  vmix = W2V(fNiMix,rNi,rB4C);
  G4double rbni = vmix*rNi + (1-vmix)*rB4C;
  fBorNi = new G4Material("BorNi",rbni*g/cm3,2);
  fBorNi->AddMaterial(fB4C,(1-fNiMix));
  fBorNi->AddMaterial(fNi,fNiMix);
  G4cout << fBorNi << std::endl;
  // Ni plating with P additive
  //G4double rNiP = rNi*(1-fPhosMix) + rPhos*fPhosMix;
  //fNiP = new G4Material("NiP",rNiP*g/cm3,2);
  //fPhos = nistM->FindOrBuildMaterial("G4_P");
  //fNiP->AddMaterial(fPhos, fPhosMix);
  //fNiP->AddMaterial(fNi, (1-fPhosMix));
  vmix = W2V(fPhosMix,rPhos,rNi);
  sprintf(nccfg,"phases<%g*Ni_sg225.ncmat&%g*solid::P/P_is_1.0_P31/%ggcm3>",
	  1-vmix,vmix,rPhos);
  fNiP = G4NCrystal::createMaterial(nccfg);
  G4cout << fNiP << std::endl;
  // CH2
  // Cadmium
  fCd = nistM->FindOrBuildMaterial("G4_Cd");
  // Li glasses
  G4Element* enLi = new G4Element("enLi", "enLi", 2);
  enLi->AddIsotope(new G4Isotope("Li6",3,6,6.015*g/mole), 0.95);
  enLi->AddIsotope(new G4Isotope("Li7",3,7,7.016*g/mole), 0.05);
  G4Element* Ce  = nistM->FindOrBuildElement("Ce");
  G4Material* Ce2O3 = new G4Material("Ce2O3", 6.2*g/cm3 , 2);
  Ce2O3->AddElement(Ce,2);
  Ce2O3->AddElement(O,3);
  G4Material* Al2O3 = nistM->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  G4Material* SiO2 = nistM->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Li2O = nistM->FindOrBuildMaterial("G4_LITHIUM_OXIDE");
  G4Material* MgO = nistM->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
  fGS1 = new G4Material("GS1",2.64*g/cm3,5);
  fGS1->AddMaterial(SiO2,55.0*perCent);
  fGS1->AddMaterial(Al2O3,11.0*perCent);
  fGS1->AddMaterial(MgO,24.0*perCent);
  fGS1->AddMaterial(Ce2O3,4.0*perCent);
  fGS1->AddMaterial(Li2O,6.0*perCent);
  G4cout << fGS1 << std::endl;
  G4Material* enLi2O = new G4Material("enLi2O",2.013*g/cm3,2);
  enLi2O->AddElement(enLi,2);
  enLi2O->AddElement(O,1);
  fGS20 = new G4Material("GS20",2.50*g/cm3,5);
  fGS20->AddMaterial(SiO2,56.2*perCent);
  fGS20->AddMaterial(Al2O3,18.2*perCent);
  fGS20->AddMaterial(MgO,4.2*perCent);
  fGS20->AddMaterial(Ce2O3,3.2*perCent);
  fGS20->AddMaterial(enLi2O,18.2*perCent);
  G4cout << fGS20 << std::endl;
  // Corundum
  fAl2O3 = G4NCrystal::createMaterial("Al2O3_sg167_Corundum.ncmat");
  G4cout << fAl2O3 << std::endl;
  //
  fScMatList[0] = fBorAl;
  fScMatList[1] = fAl;
  fScMatList[2] = fB4Ce;
  fScMatList[3] = fCr;
  fScMatList[4] = fB4C;
  fScMatList[5] = fNi;
  fScMatList[6] = fMB;
  fScMatList[7] = fBorNi;
  fScMatList[8] = fB10;
  fScMatList[9] = fBorAlS;
  fScMatList[10] = fHe3CO2;
  fScMatList[11] = fArCO2;
  fScMatList[12] = fGS1;
  fScMatList[13] = fGS20;
  fScMatList[14] = fB4Cnc;
  fScMatList[15] = fDrip;
  fScMatList[16] = fNiP;
  fScMatList[17] = fBorAlnc;
  fScMatList[18] = fLucite;
  fScMatList[19] = fV;
  fScMatList[20] = fWater;
  fScMatList[21] = fAl2O3;
  fScMatList[22] = fCd;
  fScMatList[23] = fPb;
  fScMatList[24] = fW;
  fScMatList[25] = fCH2;
  fScMatList[26] = fC8H8;
  fScMatList[27] = fGraf;
  return;

}

//---------------------------------------------------------------------------
void DetectorConstruction::ReadParameters(const char* file)
{
  //
  // Initialise detector parameters
  //
  const char* keylist[] = {"MG-Dim:","Run-Mode:","Energy-Hist:",
    "B4C-Thick:", "Crys-Material:", "Vessel:", "Column:", "He3-Dim:",
    "Shield:", "T-Rex:", "Chopper:", "Target:", "Scatter:",
    "Scatter-Detector:","Gas:","B10frac:","Alpha-Spectrum:", "World:",
    "Collimator:", "Beam-Monitor:", "EMMA-Shield:", "EMMA-Test:",
    "EMMA-He3:","M-Target:","H-spacer:","Cryostat:",NULL };
  enum { EMG_Dim, ERun_Mode, EEnergy_Hist, EB4C_Thick, ECrysMat,
    EVessel, EColumn, EHe3_Dim, EShield, ETrex, EChop, ETarget, EScatter,
    EScatDet,EGas, EB10frac, EAlphaSpec, EWorld, ECollimator, EBeamMonitor,
    EEMMAshield,EEMMAtest,EEMMAHe3,EMTarget,EHspacer,ECryoStat,ENULL };
  G4int ikey, iread;
  G4int ierr = 0;
  char line[256];
  char delim[64];
  char hname[32];
  char fname[32];
  G4double t1,t2,t3,t4,x,y,z,dx,dy,dz;
  FILE* pdata;
  if( (pdata = fopen(file,"r")) == NULL ){
    printf("Error opening detector parameter file: %s\n",file);
    return;
  }
  while( fgets(line,256,pdata) ){
    if( line[0] == '#' ) continue;       // comment #
    printf("%s\n",line);                 // record datum to log
    sscanf(line,"%s",delim);
    for(ikey=0; ikey<ENULL; ikey++)
      if(!strcmp(keylist[ikey],delim)) break;
    switch( ikey ){
    default:
      printf("Unrecognised delimiter: %s\n",delim);
      ierr++;
      break;
    case EWorld:
      iread = sscanf(line,"%*s%lf%lf%lf%d",&fExpX,&fExpY,&fExpZ,&fExpMat);
      if( iread < 3 ) ierr++;
      break;
    case EMG_Dim:
      iread = sscanf(line,"%*s%d%d%d%lf%lf%lf%lf%lf%lf%d%lf%d",
		     &fNx,&fNy,&fNz,&fMGx,&fMGy,&fMGz,&fMGth,&fMGsp,&fMGthr,&fIsScat,&fMGthr1,&fIsScat1);
      if( iread < 8 ) ierr++;
      else if(iread == 8)
	fMGthr = fMGth;
      break;
    case EHe3_Dim:
      iread = sscanf(line,"%*s%d%lf%lf%lf%lf%lf%lf%lf",
		     &fNHe3,&fHe3r,&fHe3l,&fHe3th,&fHe3pr,&fRad,&fdTh,&fThSt);
      if( iread != 8 ) ierr++;
      fdTh *= deg;   // convert to radians
      fThSt *= deg;
      break;
    case ERun_Mode:
      iread = sscanf(line,"%*s%d%s",&fRunMode,hname);
      if( iread != 2 ) ierr++;
      if(!strcmp(hname,"FTFP_BERT_HP")) {
	fHadronMode = EFTFP_BERT_HP;
      }
      else if (!strcmp(hname,"QGSP_BERT_HP")) {
	fHadronMode = EQGSP_BERT_HP;
      }
      else if (!strcmp(hname,"QGSP_BIC_HP")) {
	fHadronMode = EQGSP_BIC_HP;
      }
      else if (!strcmp(hname,"NoHadronic")) {
	fHadronMode = ENoHadron;
      }
      else {
	printf("Hadronic physics list %s is not recognised\n",hname);
	ierr++;
      }
      break;
    case EEnergy_Hist:
      iread = sscanf(line,"%*s%s%s",fname,hname);
      if(iread < 2)
	ierr++;
      fFName = new char[strlen(fname)+1];
      strcpy(fFName,fname);
      fHName = new char[strlen(hname)+1];
      strcpy(fHName,hname);
      break;
    case EB4C_Thick:
      iread = sscanf(line,"%*s%lf%d%lf%d%lf%d%lf%lf",
		     &fBthx,&fn1,&t1,&fn2,&t2,&fn3,&t3,&t4);
      if( iread != 8 ) ierr++;
      fBt1[0] = fBt2[0] = t1;
      fBt1[1] = t1; fBt2[1] = t2;
      fBt1[2] = fBt2[2] = t2;
      fBt1[3] = t2; fBt2[3] = t3;
      fBt1[4] = fBt2[4] = t3;
      fBt1[5] = t3; fBt2[5] = t4;
      fBt1[6] = fBt2[6] = t4;
      break;
    case ECrysMat:
      fIsCrysMat = true;
      break;
    case EShield:
      fIsMB = true;
      iread = sscanf(line,"%*s%lf%lf%lf%d%d%d%d%d%d", &fMBth, &fMBth1, &fMBth2,
		     fMBmat, fMBmat1, fMBmat2, fMBmat+1, fMBmat1+1, fMBmat2+1);
      if( iread < 6 ) ierr++;
      break;
    case ETrex:
      fIsTrex = true;
      iread = sscanf(line,"%*s%d%lf%lf%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     &fNvess,&fRad,&fdTh,&fShTh,&fShMat,fThBox,fThBox+1,fThBox+2,fThBox+3,
		     fThBox+4,fThBox+5,fThBox+6,fThBox+7,fThBox+8,fThBox+9);
      if( iread < fNvess+4 ) ierr++;
      fdTh = fdTh*deg;
      fThSt = fThSt*deg;
      fNcol = 1;
      fNTRcol = 100;
      for(int n=0; n<fNvess; n++)
	fThBox[n] = fThBox[n]*deg;
      break;
    case EVessel:
      if(fNHe3 > 0){
	printf("Vessel not compatible with He3 counters...ignore line\n");
	break;
      }
      if(fIsTrex){
	printf("Vessel not compatible with TREX\n");
	break;
      } 
      iread = sscanf(line,"%*s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%lf%d%lf",
		     &fNvess,&fRad,&fdTh,&fThSt,&fVx,&fVy,&fVz,&fVth,&fVbk,&fThGap,&fNvGap,&fThGap1,&fNvGap1,&fVfr);
      if( iread < 9 ) ierr++;
      fdTh *= deg;
      fThSt *= deg;
      break;
    case EColumn:
      iread = sscanf(line,"%*s%d%lf", &fNcol, &fMGhsp);
      if( iread != 2 ) ierr++;
      break;
    case EChop:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     &x,&y,&z,&dx,&dy,&dz,&fLamda,&fFP,&fdTOF,&fSkew);
      if( iread < 9 ) ierr++;
      fprintf(CoordFile,"Beam: 0 %g %g %g %g\n",x,y,z,fFP);
      fPos0.set(x,y,z);
      fDPos.set(dx,dy,dz);
      break;
    case ETarget:
      iread = sscanf(line,"%*s%lf%lf%lf%d%lf%lf%lf%lf",
		     &fTgLen, &fTgRad, &fTgRadI,&fTgMat, &fTgSx, &fTgSy, &fTgSz, &fTgAng);
      if( iread < 4 ) ierr++;
      fIsTarget = true;
      break;
    case EMTarget:
      iread = sscanf(line,"%*s%d%lf%d%lf%lf%lf%d%lf%lf%lf%d%lf%lf%lf",
		     &fnmTg,&fTgAng,fmTgMat,fmTgSx,fmTgSy,fmTgSz,
		     fmTgMat+1,fmTgSx+1,fmTgSy+1,fmTgSz+1,
		     fmTgMat+2,fmTgSx+2,fmTgSy+2,fmTgSz+2);
      if( iread < 10 ) ierr++;
      fIsTarget = true;
      break;
    case ECryoStat:
      fIsCryoStat = true;
      break;
    case EScatter:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%d",&fScPx,&fScPy,&fScPz,
		     fScX+fNmat,fScY+fNmat,fScZ+fNmat,fScTh+fNmat,fScMat+fNmat);
      if( iread != 8 ) ierr++;
      fNmat++;
      fIsScatter = true;
      break;
    case EScatDet:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%d",&fScDetX,&fScDetY,&fScDetZ,
		     &fScDetSx,&fScDetSy,&fScDetSz,&fScDetMat);
      if( iread != 7 ) ierr++;
      break;
    case EGas:
      iread = sscanf(line,"%*s%lf%lf%lf", &fArFrac, &fCO2Frac, &fArCO2dens);
      if( iread != 3 ) ierr++;
      fArFrac *= perCent;
      fCO2Frac *= perCent;
      fArCO2dens *= mg/cm3;
      break;
    case EB10frac:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf",
		     &fB10frac,&fAlMix,&fDripMix,&fNiMix,&fAlMixS);
      if( iread < 4 ) ierr++;
      break;
    case ECollimator:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf%d%lf%lf",
		     &fColX,&fColY,&fColZ,&fColSx,&fColSy,&fColSz,
		     &fColHx,&fColHy,&fColMat,&fWallTh,&fWallZ1);
      if( iread < 11 ) ierr++;
      break;
    case EAlphaSpec:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     fAlphaLines,fAlphaLines+1,fAlphaLines+2,fAlphaLines+3,
		     fAlphaLines+4,fAlphaLines+5,fAlphaLines+6,fAlphaLines+7,
		     fAlphaLines+8,fAlphaLines+9,fAlphaLines+10,fAlphaLines+11,
		     fAlphaLines+12,fAlphaLines+13,fAlphaLines+14,fAlphaLines+15);
      if( iread < 1 ) ierr++;
      fNAlphaLines = iread;
      break;
    case EBeamMonitor:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%d",
		     &fBMx,&fBMy,&fBMz,&fBMsx,&fBMsy,&fBMsz,&fBMmat);
      if( iread < 7 ) ierr++;
      break;
    case EEMMAshield:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     &fTabX,&fTabY,&fTabZ,&fTabSx,&fTabSy,&fTabSz,
		     &fShL0,&fShL1,&fShY,&fShZ,
		     &fShL2,&fSh2x,&fSh2z,&fSh2sx,&fSh2sz,&fSh2sz1);
      if( iread < 8 ) ierr++;
      break;
    case EEMMAtest:
      fIsEMMAtest = true;
      fNvess = fNmg = 2;
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		     fETx,fETy,fETz,fETx+1,fETy+1,fETz+1,
		     &fETrx,&fETry,&fETrz);
      if( iread < 9 ) ierr++;
      if((fETx[1]==-1.0)&&(fETy[1]==-1.0)&&(fETz[1]==-1))
	fNvess = 1;
      break;
    case EEMMAHe3:
      iread = sscanf(line,"%*s%lf%lf%lf%d%lf%lf%lf",
		     &fT3x,&fT3y,&fT3z,&fnT3,&fT3ri,&fT3ro,&fT3l);
      if( iread < 7 ) ierr++;
      break;
    case EHspacer:
      iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%d",
		     &fHsx,&fHsy,&fHsz,&fHsx1,&fHsy1,&fHsz1,&fHspmat);
      if( iread < 7 ) ierr++;
      break;
    }
    if( ierr ){
      printf("Fatal Error: invalid read of parameter line %s\n %s\n",
	     keylist[ikey],line);
      exit(-1);
    }
  }
}
