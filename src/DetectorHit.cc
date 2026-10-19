#include "DetectorHit.hh"

#include "G4Color.hh"
#include "G4VisAttributes.hh"

//---------------------------------------------------------------------------

G4Allocator<DetectorHit> DetectorHitAllocator;

DetectorHit::DetectorHit()
{
  fIDX  = 0;
  fIDY  = 0;
  fPID  = 0;
  fEdep = 0;
  fPDef = 0;
  fTime = 0;
  fPos.set(0,0,0);
  fMom.set(0,0,0);
}

//---------------------------------------------------------------------------

DetectorHit::~DetectorHit()
{;}

//---------------------------------------------------------------------------

DetectorHit::DetectorHit(const DetectorHit& right)
  :G4VHit()
{
  fVid  = right.fVid;
  fIDX  = right.fIDX;
  fIDY  = right.fIDY;
  fPDG  = right.fPDG;
  fPID  = right.fPID;
  fEdep = right.fEdep;
  fPDef = right.fPDef;
  fTime = right.fTime;
  fPos  + right.fPos;
  fMom  + right.fMom;
  fKE   = right.fKE;
}

//---------------------------------------------------------------------------

const DetectorHit& DetectorHit::operator=(const DetectorHit& right)
{
  fVid  =right.fVid;
  fIDX  =right.fIDX;
  fIDY  =right.fIDY;
  fPDG  =right.fPDG;
  fPID  =right.fPID;
  fEdep =right.fEdep;
  fPDef =right.fPDef;
  fTime =right.fTime;
  fPos  +right.fPos;
  fMom  +right.fMom;
  fKE   =right.fKE;
  return *this;
}

//---------------------------------------------------------------------------
void DetectorHit::Draw()
{;}

//---------------------------------------------------------------------------
void DetectorHit::Print()
{;}

//---------------------------------------------------------------------------







