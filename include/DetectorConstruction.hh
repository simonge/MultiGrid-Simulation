// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Detector geometry, materials including optical properties
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GDMLParser.hh"
#include "DetectorSD.hh"

// Mode of simulation run
enum {ERunOptical, ERunFast, ERunFastOpt, ERunThermal, ERunThermalOpt};
// Hadronic Physics List
enum {EFTFP_BERT_HP, EQGSP_BERT_HP, EQGSP_BIC_HP, ENoHadron};
enum {EMaxHits = 16384};
//enum {EBoral=1,EAl,EB4Ce,ECr,EB4C,ENi,EMB,EBorNi,EB10,EBorAlS,EHe3CO2};
//
const G4double rdn = 57.29578;
const G4int nMat = 4;
const G4int nVxSz = 7;
const G4int nMGtype = 2;
class G4LogicalVolume;
class G4VPhysicalVolume;
class DetectorSD;
class G4SubtractionSolid;
class G4Material;

const G4double rAl = 2.699;    // dens Al
const G4double rNi = 8.907;    // dens Ni
const G4double rPhos = 2.20;   // dens P (additive to Ni)
//const G4double rNiP = 8.37;    // dens Ni-P with 8% wgt P
const G4double rB4C = 2.52;    // dens nat B4C
const G4double rEpoxy = 1.20;  // dens epoxy matrix
const G4double rMb = 1.55;     // dens Mirrobor
const G4double rDrip = 1.55;    // dens dripped mixture
const G4double rHe3 = 0.000166; // He3 gas @ STP
const G4double rCO2 = 0.001977; // CO2 gas @ STP
const G4double rAr = 0.001784;  // Ar gas @ STP
const G4double thMB = 5.0;      // thickness Mirrobor sheet
//
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(const char* ="macros/MG.dat");
    virtual ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
  void ReadParameters(const char*);
  void Grid(G4double,G4int);
  void MultiGrid(G4int=0,G4LogicalVolume* =NULL,
		 G4double=0, G4double=0, G4double=0);
  void He3Tube();
  void ScatterTest();
  void Vessel(G4int=0);
  void Target();
  void EMMAtest();
  void EMMAHe3();
  void Materials();
  void ShieldCH2();
  void ShieldPb();
  void EMMAshield();
  void Shield(G4int=0); // internal MG shielding
  void Trex();
  void Collimator();    // neutron collimation
  void BeamMonitor();   // beam monitor
  void HSpacer();       // spacers between grids
  void TRbox();         // TREX box for 6 MG columns
  void CryoStat();      // cryostat around target
  inline DetectorSD*        GetDetSD()  { return fDetSD; }
  inline G4VPhysicalVolume* GetWorldVolume() const { return fExpPhys; }
  inline G4int GetRunMode(){ return fRunMode; }
  inline G4int GetHadronMode(){ return fHadronMode; }
  inline G4bool IsOptical(){ return fIsOptical; }
  inline G4bool IsThermal(){ return fIsThermal; }
  inline G4double W2V(G4double w1, G4double r1, G4double r2){
    // weight fraction to volume fraction
    G4double v1 = w1/r1;
    G4double v2 = (1-w1)/r2;
    v1 = v1/(v1+v2);
    return v1;
  }
  inline G4double V2W(G4double v1, G4double r1, G4double r2){
    // volume fraction to weight fraction
    G4double w1 = v1*r1;
    G4double w2 = (1-v1)*r2;
    w1 = w1/(w1+w2);
    return w1;
  }
  G4LogicalVolume* fExpLog;
  G4LogicalVolume* fShellLog=NULL;
  G4VPhysicalVolume* fExpPhys;
  G4VPhysicalVolume* fShellPhys=NULL;
  //
  private:
  G4double fExpX;    // experiment dimensions
  G4double fExpY;
  G4double fExpZ;
  G4double fWallTh=0.0; // any concrete wall thickness
  G4double fWallZ1=0.0; // 2nd wall?
  G4int    fExpMat=0; // experiment volume material...air,vacuum(default)....
  G4double fMGx;    // cell dimensions
  G4double fMGy;
  G4double fMGz;
  G4double fMGth;   // Al normal blade thickness
  G4double fMGthr;  // radial blade thickness
  G4double fMGthr1=0;  // alternative radial blade thickness
  G4double fMGsp;   // MG vertical spacing
  G4double fMGhsp;  // MG column horizontal spacing
  G4double fGrx,fGry,fGrz; // Grid volume dimensions
  G4double fHe3r;   // He3 radius
  G4double fHe3l;   // He3 length
  G4double fHe3th;  // He3 tube wall thickness
  G4double fHe3pr = 10.0;  // He3 gas pressure 10 bar default
  G4double fWrR;    // Wire radius
  G4int    fNx;     // # cells horizontally
  G4int    fNy;     // # grids vertically
  G4int    fNz;     // # cells depth-wise
  G4int    fNcol;   // # multigrid columns
  G4int    fNTRcol=100;   // multigrid column index
  G4int    fNvess=0;// # multigrid vessels
  G4int    fNHe3=0; // # He3 tubes
  G4double fRad;    // vessel placement radius
  G4double fdTh;    // vessel placement angle spacing
  G4double fShTh=0.0;   // thickness shield plate between boxes
  G4int    fShMat=18;  // material of shielding plate
  G4double fThBox[10]; // TREX box angle, max 10 TREX boxes
  G4double fThSt;   // vessel placement start angle
  G4double fThGap = 0;  // direct beam gap angle
  G4double fThGap1= 0;  // box-box gap angle
  G4int    fNvGap = 0;  // vessel # where beam gap occurs
  G4int    fNvGap1= 0;  // vessel # where box gap occurs
  G4double fVfr= 0;   // front vessel thickness
  G4double fBthx;   // B4C coating
  G4double fBt1[nVxSz]; // B4C coating
  G4double fBt2[nVxSz]; // B4C coating
  G4int    fn1,fn2,fn3; // for setting up coating
  G4double fVx;     // containment vessel internal x,y,z dimensions
  G4double fVy;
  G4double fVz;
  G4double fVth;    // containment vessel thickness
  G4double fVbk;    // backplate thickness
  G4double fMBth=0;   // back shield 1/2 thickness
  G4double fMBth1=0;  // side shield 1/2 thickness
  G4double fMBth2=0;  // back block 1/2 thickness
  G4int    fMBmat[nMGtype] = {12,0};
  G4int    fMBmat1[nMGtype] = {2,0};
  G4int    fMBmat2[nMGtype] = {2,0};
  G4double fTRsp0,fTRsp1;
  G4double fTRx1,fTRx2,fTRx3,fTRz1,fTRz2,fTRz3;
  G4double fLamda,fFP,fdTOF,fSkew=0;            // chopped beam parameters
  G4double fScPx=0,fScPy=0,fScPz=0;             // centre of scatterer
  G4double fScX[nMat],fScY[nMat],fScZ[nMat],fScTh[nMat]; // scattering sample 
  G4int    fScMat[nMat];
  G4double fScDetX=0,fScDetY=0,fScDetZ=0;
  G4double fScDetSx=0,fScDetSy=0,fScDetSz=500.0;
  G4double fAlphaLines[16] = {0.0};
  G4double fETx[nMGtype] = {0.0};
  G4double fETy[nMGtype] = {0.0};
  G4double fETz[nMGtype] = {0.0};
  G4double fETrx=0.0,fETry=0.0,fETrz=0.0;
  G4int    fNAlphaLines = 0;
  G4int    fScDetMat;
  G4int    fNmat = 0;
  G4int    fNmg = 1;        // # different configurations of grids, default 1
  G4int    fRunMode;
  G4int    fHadronMode;
  G4bool   fIsOptical;
  G4bool   fIsThermal;
  G4bool   fIsCrysMat;
  G4bool   fIsMB;
  G4bool   fIsTrex;
  G4bool   fIsTarget = false;
  G4bool   fIsCryoStat = false;
  G4bool   fIsScatter = false;
  G4int    fIsScat=1,fIsScat1=0;
  G4bool   fIsEMMAtest = false;
  G4double fBMx=0,fBMy=0,fBMz=0;
  G4double fBMsx=0,fBMsy=0,fBMsz=0;
  G4int    fBMmat=-1;
  G4int    fHspmat=-1;
  G4double fHsx,fHsy,fHsz,fHsx1,fHsy1,fHsz1;
  G4double fWr;     // beamline window radius, thickness dist
  G4double fWz;
  G4double fWd;
  G4double fRspot;  // beam spot radius
  G4double fMposz;  // rear mirrorbor position
  G4double fMposx;  // side mirrorbor position
  G4double fEscale; // energy scale
  G4double fColX,fColY,fColZ;           // neutron colli position
  G4double fColSx=0,fColSy=0,fColSz=0;  // neutron colli 1/2 size
  G4double fColHx=0,fColHy=0;           // colli aperture 1/2 size
  G4int    fColMat=5;                   // B4C material
  G4double fTabX,fTabY,fTabZ;           // Al table at EMMA
  G4double fTabSx=0,fTabSy,fTabSz;        // Al table at EMMA
  G4double fShL0,fShL1,fShL2=0,fShY=0,fShZ=0;
  G4double fSh2x=0,fSh2z=0,fSh2sx,fSh2sz,fSh2sz1;
  G4double fT3x,fT3y,fT3z,fT3ri,fT3ro,fT3l;
  G4int    fnT3=0;
  G4int    fnmTg=0;
  G4double fmTgSx[4] = {0.0};
  G4double fmTgSy[4] = {0.0};
  G4double fmTgSz[4] = {0.0};
  G4int    fmTgMat[4] = {-1};
  G4VPhysicalVolume** fphysMG;
  G4LogicalVolume* flogMG[nMGtype]={NULL};
  G4LogicalVolume* flogArCO2[nVxSz]={NULL};
  G4LogicalVolume* flogbx2[nVxSz]={NULL};
  G4LogicalVolume* flogGrid=NULL;
  G4LogicalVolume* flogVess[nMGtype]={NULL};
  G4LogicalVolume* flogVGas[nMGtype]={NULL};
  G4LogicalVolume* flogHe3=NULL;
  G4LogicalVolume* flogTube=NULL;
  G4LogicalVolume* flogWire=NULL;
  G4LogicalVolume* flogMB=NULL;
  G4LogicalVolume* flogMB1=NULL;
  G4LogicalVolume* flogTR1=NULL;
  G4LogicalVolume* flogTR2=NULL;
  G4LogicalVolume* flogTR3=NULL;
  G4LogicalVolume* flogTarget=NULL;
  G4LogicalVolume* flogScat[5]={NULL};
  G4LogicalVolume* flogScDet=NULL;
  G4LogicalVolume* flogBM=NULL;
  G4LogicalVolume* flogT3=NULL;
  G4LogicalVolume* flogHsp=NULL;
  DetectorSD* fDetSD;
  G4Material* fAir;    // air at STP
  G4Material* fVac;    // high vacuum
  G4Material* fWater;  // water
  G4Material* fAl;
  G4Material* fFe;
  G4Material* fAlpc;
  G4Material* fPb;
  G4Material* fW;
  G4Material* fV;
  G4Material* fCd;
  G4Material* fCH2;    // polythene ncrystal
  G4Material* fB4C;    // natural boron carbide
  G4Material* fB4Ce;   // boron carbide enriched 10B
  G4Material* fB4Cnc;   // boron carbide in NCrystal
  G4Material* fArCO2;  // MG gas
  G4Material* fHe3CO2; // He3 tube gas
  G4Material* fHe3Ar;  // He3 tube gas
  G4Material* fMB;     // MirroBor B based shielding
  G4Material* fDrip;   // Dripped B based shielding
  G4Material* fGraf;   // ncrystal graphite
  G4Material* fC8H8;   // ncrystal polystyrene
  G4Material* fBorAl;  // Al - B4C mixture solid sheets
  G4Material* fBorAlnc;// Al - B4C mixture solid sheets (B4C in NCrystal)
  G4Material* fBorAlS; // Al - B4C mixture spray
  G4Material* fBorNi;  // Ni - B4C mixture spray
  G4Material* fCr;     // Cromium plating
  G4Material* fNi;     // Ni plating
  G4Material* fPhos;   // Phosphorous additive
  G4Material* fNiP;    // Ni plating with 8% wgt of P
  G4Material* fB10;    // hypothetical boron-10
  G4Material* fGS1;    // Nat-Li scint glass
  G4Material* fGS20;   // 95% Li6 scint glass
  G4Material* fConcrete;// Concrete
  G4Material* fLucite; // Perspex
  G4Material* fAl2O3;  // Corundum
  G4Material* fScMatList[48]={NULL};
  // 83% B4C, 0.85% H,
  // glue (20%) 70% C, 25% O, 5% H
  // dens 1.55 gcm-3
  G4Material* fScatMat[nMat] = {NULL};
  G4double fArFrac = 78.4*CLHEP::perCent;
  G4double fCO2Frac = 21.6*CLHEP::perCent;
  G4double fArCO2dens = 1.796*CLHEP::mg/CLHEP::cm3;
  G4double fTgLen,fTgRad,fTgRadI;
  G4double fTgSx=0, fTgSy=0, fTgSz=0, fTgAng=0;
  G4int    fTgMat;
  G4double fB10frac = 0.199; // 10B abundance in natural boron
  G4double fAlMix = 0.69;    // fraction (weight) Al in BorAl sheet mixture
  G4double fAlMixS = 0.52;   // fraction (weight) Al in BorAl spray mixture
  G4double fNiMix = 0.779;   // fraction (weight) Ni in BorNi spray mixture
  G4double fDripMix = 0.83;  // fraction (weight) BorAl in epoxy mixture
  G4double fPhosMix = 0.08;  // fraction phosphorous (wgt) in Ni plating
 //
  // Stuff for stepped scans
  G4ThreeVector fPos0;
  G4ThreeVector fDPos;
  G4double      fDTh;
  G4int         fNscan;
  // Stuff for sampling from root histogram
  char*         fFName;
  char*         fHName;
public:
  G4ThreeVector GetPos0(){ return fPos0; }
  G4ThreeVector GetDPos(){ return fDPos; }
  G4double      GetDTh(){ return fDTh; }
  G4double      GetRspot(){ return fRspot;}
  G4int         GetNscan(){ return fNscan; }
  char*         GetFName(){ return fFName; }
  char*         GetHName(){ return fHName; }
  G4int         GetNCell(){ return fNx*fNz; }
  G4bool        IsCrysMat(){ return fIsCrysMat; }
  G4double      GetLamda(){ return fLamda; }
  G4double      GetFP(){ return fFP; }
  G4double      GetdTOF(){ return fdTOF; }
  G4double      GetSkew(){ return fSkew; }
  G4int         GetNAlphaLines(){ return fNAlphaLines; }
  G4double*     GetAlphaLines(){ return fAlphaLines; }
  void          GetScatter(G4double& x, G4double& y, G4double& z,
			   G4double& dx, G4double& dy, G4double& dz){
    x = fScPx; y = fScPy; z = fScPz;
    dx = fScX[fNmat-1]; dy = fScY[fNmat-1]; dz = fScZ[fNmat-1];
  }
};

#endif /*DetectorConstruction_h*/
