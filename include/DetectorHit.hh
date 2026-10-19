#ifndef DetectorHit_h
#define DetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"

//---------------------------------------------------------------------------

class DetectorHit : public G4VHit
{
public:
  
  DetectorHit();
  ~DetectorHit();
  DetectorHit(const DetectorHit&);
  const DetectorHit& operator=(const DetectorHit&);
  inline int operator==(const DetectorHit&) const{ return 0; }
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  virtual void Draw();
  void Print();

protected:
  G4int                  fVid;
  G4int                  fPDG;
  G4int                  fIDX; 
  G4int                  fIDY; 
  G4int                  fPID; 
  G4double               fEdep;  
  G4ParticleDefinition*  fPDef; 
  G4double               fTime;
  G4ThreeVector          fPos; 
  G4ThreeVector          fMom;  
  G4double               fKE;
public:
  
  inline void SetVid      (G4int iv)                 { fVid = iv;   };
  inline void SetIDX      (G4int ix)                 { fIDX = ix;   };
  inline void SetIDY      (G4int iy)                 { fIDY = iy;   };
  inline void SetPDG      (G4int pdg)                { fPDG = pdg;  };
  inline void SetParentID (G4int j)                  { fPID = j;    };
  inline void AddEnergy   (G4double de)              { fEdep = de; };
  inline void SetPDef     (G4ParticleDefinition* pd) { fPDef = pd;  };
  inline void SetTime     (G4double t)               { fTime = t;   };
  inline void SetPos      (G4ThreeVector pos)        { fPos=pos;    };
  inline void SetMomentum (G4ThreeVector p3)         { fMom = p3;   };
  inline void SetKE       (G4double ke)              { fKE = ke;   };
  
  inline G4int                 GetVid()              { return fVid; };
  inline G4int                 GetIDX()              { return fIDX; };
  inline G4int                 GetIDY()              { return fIDY; };
  inline G4int                 GetPDG()              { return fPDG; };
  inline G4int                 GetParentID()         { return fPID; };
  inline G4ParticleDefinition* GetPDef()             { return fPDef;};
  inline G4double              GetTime()             { return fTime;};
  inline G4double              GetEdep()             { return fEdep;};
  inline G4ThreeVector         GetPos()              { return fPos; };
  inline G4ThreeVector         GetMom()              { return fMom; };
  inline G4double              GetKE()               { return fKE; }
};

//---------------------------------------------------------------------------

typedef G4THitsCollection<DetectorHit> DetectorHitsCollection;

extern G4Allocator<DetectorHit> DetectorHitAllocator;


inline void* DetectorHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) DetectorHitAllocator.MallocSingle();
  return aHit;
}


inline void DetectorHit::operator delete(void* aHit)
{
  DetectorHitAllocator.FreeSingle((DetectorHit*) aHit);
}

#endif

//---------------------------------------------------------------------------




