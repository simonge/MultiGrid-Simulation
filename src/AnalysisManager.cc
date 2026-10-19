// Geant-4 simulation of the Solid state Neutron Detector (SoNDe)
// for the European Spallation Source in Lund Sweden
// Write hits info to root file
// JRM Annand 13/03/2018
//--------------------------------------------------------------------------
#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Exception.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include <boost/math/distributions/skew_normal.hpp>
using namespace boost::math;

//---------------------------------------------------------------------------

AnalysisManager::AnalysisManager()
{
  fAnaMessenger = new AnalysisMessenger(this);
  fOutFileName = "output/out_default.root";
  fEscale = 1.0;
}

//---------------------------------------------------------------------------

AnalysisManager::~AnalysisManager()
{
   if (fROOTtree)
     fROOTtree->Write();
   if (fROOTfile) {
     fROOTfile->Close();
     delete fROOTfile;
   }
   delete fAnaMessenger;
}

//---------------------------------------------------------------------------

void AnalysisManager::InitOutput()
{
  // Setup the ROOT TTree file
  fROOTfile = new TFile(fOutFileName.c_str(),"RECREATE","fROOTfile",1);
  if (!fROOTfile || fROOTfile->IsZombie()) {
    G4cerr << "Could not create ROOT output file: " << fOutFileName << G4endl;
    G4Exception("AnalysisManager::InitOutput", "MGOutput001", FatalException,
                "Failed to create the requested ROOT output file.");
    return;
  }
  fROOTtree = new TTree("PD","Optical Photons");
  fROOTtree->SetAutoSave();

  // Set Photon Detector Branches
  fROOTtree->Branch("PD_Nhits", &fPDNhits, "PD_Nhits/I");
  fROOTtree->Branch("PD_T0n", &fPDT0n,   "PD_T0n/F");
  fROOTtree->Branch("PD_E0n", &fPDE0n,   "PD_E0n/F");
  fROOTtree->Branch("PD_pdg", &fPDpdg,   "PD_pdg[PD_Nhits]/I" );
  fROOTtree->Branch("PD_vid", &fPDvid,   "PD_vid[PD_Nhits]/I" );
  fROOTtree->Branch("PD_m",   &fPDm,     "PD_m[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_t",   &fPDt,     "PD_t[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_x",   &fPDx,     "PD_x[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_y",   &fPDy,     "PD_y[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_z",   &fPDz,     "PD_z[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_vx",   &fPDvx,    "PD_vx[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_vy",   &fPDvy,    "PD_vy[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_vz",   &fPDvz,    "PD_vz[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_px",  &fPDpx,    "PD_px[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_py",  &fPDpy,    "PD_py[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_pz",  &fPDpz,    "PD_pz[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_th",  &fPDth,    "PD_th[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_ph",  &fPDph,    "PD_ph[PD_Nhits]/F"  );
  fROOTtree->Branch("PD_p",   &fPDp,     "PD_p[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_E",   &fPDE,     "PD_E[PD_Nhits]/F"   );
  fROOTtree->Branch("PD_Ed",  &fPDEd,    "PD_Ed[PD_Nhits]/F"  );
}

//---------------------------------------------------------------------------

void AnalysisManager::ZeroArray()
{
  // Zero the data arrays for tree branches
  G4ThreeVector zero(0.,0.,0.);
  fPDpdef  = NULL;
  fPDpos   = (zero);
  fPDp3    = (zero);
  fPDtime  = 0;
  fPDedep  = 0;
  for ( Int_t i = 0; i < fPDNhits; i++ ) {
    fPDpdg[i] = 0;
    fPDvid[i] = 0;
    fPDm[i] = 0.;
    fPDt[i] = 0.;
    fPDx[i] = 0.;
    fPDy[i] = 0.;
    fPDz[i] = 0.;
    fPDvx[i] = 0.;
    fPDvy[i] = 0.;
    fPDvz[i] = 0.;
    fPDpx[i] = 0.;
    fPDpy[i] = 0.;
    fPDpz[i] = 0.;
    fPDth[i] = 0.;
    fPDph[i] = 0.;
    fPDp[i] = 0.;
    fPDE[i] = 0.;
    fPDEd[i] = 0.;
  }
  fPDNhits = 0;
}

//---------------------------------------------------------------------------

void AnalysisManager::FillPDArray( Int_t hitn )
{
  if(hitn >= EMaxHits)
    return;
  // Increment hit number
  fPDNhits++;

  // Fill branch arrays from info transferred from SD hit
  fPDpdg[hitn]   = (Int_t)  fPDpdef->GetPDGEncoding();
  fPDvid[hitn]   = (Int_t)  fPDvid_single;
  fPDm[hitn]     = (Float_t)fPDpdef->GetPDGMass();                     
  fPDt[hitn]     = (Float_t)fPDtime;                                   
  fPDx[hitn]     = (Float_t)fPDpos.getX();                             
  fPDy[hitn]     = (Float_t)fPDpos.getY();                             
  fPDz[hitn]     = (Float_t)fPDpos.getZ();                             
  fPDvx[hitn]    = (Float_t)fPDvertex.getX();                             
  fPDvy[hitn]    = (Float_t)fPDvertex.getY();                             
  fPDvz[hitn]    = (Float_t)fPDvertex.getZ();                             
  fPDpx[hitn]    = (Float_t)fPDp3.getX();                             
  fPDpy[hitn]    = (Float_t)fPDp3.getY();                             
  fPDpz[hitn]    = (Float_t)fPDp3.getZ();                             
  fPDth[hitn]    = (Float_t)fPDp3.getTheta();                         
  fPDph[hitn]    = (Float_t)fPDp3.getPhi();                          
  fPDp[hitn]     = (Float_t)fPDp3.mag();                             
  //  fPDE[hitn]     = TMath::Sqrt( fPDp[hitn]*fPDp[hitn] + fPDm[hitn]*fPDm[hitn]);
  fPDE[hitn]     = (Float_t)fPDke;
  fPDEd[hitn]    = (Float_t)fPDedep;
}

//---------------------------------------------------------------------------

void AnalysisManager::FillTree()
{              
  fROOTtree->Fill();
}

//---------------------------------------------------------------------------
void AnalysisManager::SetOutFileName(const std::string& fname,
                                     G4bool lockForCommandLine) {
  if (fOutFileNameLocked && !lockForCommandLine) {
    G4cout << "Ignoring macro output-file setting; command-line output takes precedence."
           << G4endl;
    return;
  }
  fOutFileName = fname;
  fOutFileNameLocked = lockForCommandLine;
  G4cout << "Output file name " << fOutFileName << G4endl;
}

//---------------------------------------------------------------------------
void AnalysisManager::SetInput1D(const char* fname, const char* hname) {
  // Setup data sampling from ROOT 1D hist held in TFile fname
  TFile* file = new TFile(fname);
  if(!file){
    printf("Error opening file %s\n",fname);
    exit(1);
  }
  file->ls();
  //Get 1D hist
  fH1Input = (TH1D*)file->Get(hname);
  if(!fH1Input){
    printf("Setup 1D histogram input from file %s:no 1D histogram %s\n",fname,hname);
    file->ls();
    exit(1);
  }
}

//---------------------------------------------------------------------------
G4double AnalysisManager::SampleH1() {
  // Setup data sampling from ROOT 1D hist held in TFile fname
  return fEscale * fH1Input->GetRandom();
}

//--------------------------------------------------------------------------
G4double SkPk(Double_t* x, Double_t* par)
{
  // Inverted (about mean mu) Landau
  Double_t mu = par[0];
  Double_t sig = par[1];
  Double_t xx = 2*mu - x[0];
  Double_t func = TMath::Landau(xx,mu,sig);
  return func;
}
//---------------------------------------------------------------------------
void AnalysisManager::SetSkGaus(Double_t mean, Double_t sig, Double_t skew )
{
  // skewed Gaussian
  // peak offset correction
  //  Double_t sk = fabs(skew);
  //  Double_t corr = 0.692 - 0.093*sk + 0.0070*sk*sk -0.00027*sk*sk*sk +
  //    4.03e-6*sk*sk*sk*sk;
  
  //Double_t lim1 = mean - 35*sig;
  Double_t lim1 = mean - 90*sig;
  Double_t lim2 = mean + 10*sig;
  //fSkGaus = new TF1("skGaus","pdf(skew_normal([0],[1],[2]), x)",lim1,lim2);
  //fSkGaus->SetParameters(mean,sig,skew);
  fSkGaus = new TF1("skk",SkPk,lim1,lim2,2);
  fSkGaus->SetParameters(mean,sig);
  return;
}
//
//---------------------------------------------------------------------------
G4double AnalysisManager::SampleTOF()
{
  return fSkGaus->GetRandom();
}
