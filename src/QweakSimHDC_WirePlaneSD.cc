//=============================================================================
// 
//   ---------------------------
//  | Doxygen File Information |
//  ---------------------------
// 
/**
 
   \file QweakSimHDC_WirePlaneSD.cc

   $Revision: 1.3 $	
   $Date: 2006/05/05 21:39:12 $

   \author Klaus Hans Grimm   

*/
//=============================================================================

//=============================================================================
//   -----------------------
//  | CVS File Information |
//  -----------------------
// 
//  Last Update:      $Author: grimm $
//  Update Date:      $Date: 2006/05/05 21:39:12 $
//  CVS/RCS Revision: $Revision: 1.3 $
//  Status:           $State: Exp $
// 
// ===================================
//  CVS Revision Log at end of file !!
// ===================================
//
//============================================================================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "QweakSimHDC_WirePlaneSD.hh"

// user includes
#include "QweakSimSensitiveDetectorDefinition.hh"
#include "QweakSimHDC_WirePlaneHit.hh"
#include "QweakSimTrackInformation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QweakSimHDC_WirePlaneSD::QweakSimHDC_WirePlaneSD(G4String name)
:G4VSensitiveDetector(name)
{
//  G4cout << G4endl << "###### Calling QweakSimHDC_WirePlaneSD::QweakSimHDC_WirePlaneSD() " << G4endl << G4endl; 

  collectionName.insert("HDCWirePlaneCollection"); //collectionName[0]
   HCID = -1;

//  G4cout << G4endl << "###### Leaving QweakSimHDC_WirePlaneSD::QweakSimHDC_WirePlaneSD() " << G4endl << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
QweakSimHDC_WirePlaneSD::~QweakSimHDC_WirePlaneSD()
{
  //G4cout << G4endl << "###### Calling QweakSimHDC_WirePlaneSD::~QweakSimHDC_WirePlaneSD() " << G4endl << G4endl; 
  //delete hitsCollection;
  //G4cout << G4endl << "###### Leaving QweakSimHDC_WirePlaneSD::~QweakSimHDC_WirePlaneSD() " << G4endl << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QweakSimHDC_WirePlaneSD::Initialize(G4HCofThisEvent* HCE)
{
  //G4cout << G4endl << "###### Calling QweakSimHDC_WirePlaneSD::Initialize() " << G4endl << G4endl; 

    hitsCollection     = new QweakSimHDC_WirePlane_HitsCollection(SensitiveDetectorName,collectionName[0]);
    
     HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
  // HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); 

     HCE->AddHitsCollection( HCID, hitsCollection );

  //G4cout << G4endl << "###### Leaving QweakSimHDC_WirePlaneSD::Initialize() " << G4endl << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool QweakSimHDC_WirePlaneSD::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
  //G4cout << G4endl << "###### Calling QweakSimHDC_WirePlaneSD::ProcessHits() " << G4endl << G4endl; 


  G4double  charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  G4String  particlename = aStep->GetTrack()->GetDefinition()->GetParticleName();
  G4int     particletype = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  
//   G4cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Charge =" << charge << G4endl;

  // reject non-charged particle hits
  if (fabs(charge)<0.1) 
      return false;

  G4StepPoint*        preStepPoint  = aStep->GetPreStepPoint();
  G4StepPoint*        postStepPoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4VPhysicalVolume*  physVol      = theTouchable->GetVolume();

  // we are only interested if the particle is crossing a boundary (glass->vaccuum)
  if (postStepPoint->GetStepStatus() != fGeomBoundary) { 
//       G4cout << "=======================================================================" << G4endl; 
//       G4cout << ">>>> QweakSimHDC_WirePlaneSD: We are NOT crossing a boundary  <<<<<<" << G4endl; 
//       G4cout << ">>>>       Aborting QweakSimHDC_WirePlaneSD::ProcessHits      <<<<<<" << G4endl; 
//       G4cout << "=======================================================================" << G4endl;   

    return false;
  }

//   G4cout << "================== We are in volume :" << physVol->GetName() << G4endl << G4endl;


  if(strcmp(physVol->GetName(),"HDC_WirePlane_Physical")==0) {

//     G4cout << "=============================================================" << G4endl;
//     G4cout << ">>>>>>> Particle crossing  : HDC_WirePlane_Physical <<<<<<<<<" << G4endl;
//     G4cout << "=============================================================" << G4endl;

//     G4int trackID = aStep->GetTrack()->GetTrackID();
//     G4cout << "====> Track ID     : " <<  trackID << G4endl;
    
//     G4int parentID = aStep->GetTrack()->GetParentID();
//     G4cout << "====> Parent ID    : " <<  parentID << G4endl;
  
  

//     G4int    MotherCopyNo = theTouchable->GetVolume(1)->GetCopyNo();  // Several MotherVolumes
    G4int WirePlaneCopyNo = theTouchable->GetVolume()->GetCopyNo();   // but only one WirePlane per MV
//     G4int MotherReplicaNo = theTouchable->GetReplicaNumber(1);        // Several MotherVolumes

    G4int    MotherCopyNo2 = theTouchable->GetVolume(2)->GetCopyNo();  // Several MotherVolumes
//     G4int MotherReplicaNo2 = theTouchable->GetReplicaNumber(2);        // Several MotherVolumes


//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane MV    CopyNumber  :" << MotherCopyNo    << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane       CopyNumber  :" << WirePlaneCopyNo << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane    ReplicaNumber  :" << MotherReplicaNo << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane MV    CopyNumber2 :" << MotherCopyNo2    << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane    ReplicaNumber2 :" << MotherReplicaNo2 << G4endl;
    

    G4ThreeVector worldPos = preStepPoint->GetPosition();
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane WorldPos.X [cm]:" << worldPos.x()/cm << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane WorldPos.Y [cm]:" << worldPos.y()/cm << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane WorldPos.Z [cm]:" << worldPos.z()/cm << G4endl;
    
    G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane LocalPos.X [cm]:" << localPos.x()/cm << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane LocalPos.Y [cm]:" << localPos.y()/cm << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane LocalPos.Z [cm]:" << localPos.z()/cm << G4endl;
  
    G4ThreeVector worldMomentum = preStepPoint->GetMomentum();
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane WorldMomentum.X [GeV]:" << worldMomentum.x()/GeV << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane WorldMomentum.Y [GeV]:" << worldMomentum.y()/GeV << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane WorldMomentum.Z [GeV]:" << worldMomentum.z()/GeV << G4endl;
    
    G4ThreeVector localMomentum = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldMomentum);
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane LocalMomentum.X [GeV]:" << localMomentum.x()/GeV << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane LocalMomentum.Y [GeV]:" << localMomentum.y()/GeV << G4endl;
//     G4cout << "%%%%%%%%%%%%%%%%%%%  WirePlane LocalMomentum.Z [GeV]:" << localMomentum.z()/GeV << G4endl;



 // Vertex: where this track was generated
  G4ThreeVector originVertexPosition          = aStep->GetTrack()->GetVertexPosition();
  G4ThreeVector originVertexMomentumDirection = aStep->GetTrack()->GetVertexMomentumDirection();
  G4double      originVertexKineticEnergy     = aStep->GetTrack()->GetVertexKineticEnergy();

//   G4cout << "====> originVertexPosition_X [cm]     : " <<  originVertexPosition.x()/cm << G4endl;
//   G4cout << "====> originVertexPosition_Y [cm]     : " <<  originVertexPosition.y()/cm << G4endl;
//   G4cout << "====> originVertexPosition_Z [cm]     : " <<  originVertexPosition.z()/cm << G4endl;

//   G4cout << "====> originVertexMomentumDirection_X : " <<  originVertexMomentumDirection.x() << G4endl;
//   G4cout << "====> originVertexMomentumDirection_Y : " <<  originVertexMomentumDirection.y() << G4endl;
//   G4cout << "====> originVertexMomentumDirection_Z : " <<  originVertexMomentumDirection.z() << G4endl;

//   G4cout << "====> originVertexKineticEnergy [Mev] : " <<  originVertexKineticEnergy/MeV << G4endl;


  G4double      currentKineticEnergy     = preStepPoint->GetKineticEnergy();
  G4double      currentTotalEnergy       = preStepPoint->GetTotalEnergy();

  // get User Track Info
  //QweakSimTrackInformation* info = (QweakSimTrackInformation*) (aStep->GetTrack()->GetUserInformation());
  //info->Print();

//   G4double primaryQ2                      = info->GetPrimaryQ2();
//   G4double crossSection                   = info->GetCrossSection();
//   G4double crossSectionWeight             = info->GetCrossSectionWeight();
//   G4int    primaryEventNumber             = info->GetPrimaryEventNumber();

//   G4cout << " @@@@@@@@@@@@@@@@@@@@@@@ Original Track Q2                 : " << primaryQ2                       << G4endl;
//   G4cout << " @@@@@@@@@@@@@@@@@@@@@@@ Original Track Event Number       : " << primaryEventNumber              << G4endl;
//   G4cout << " @@@@@@@@@@@@@@@@@@@@@@@ Original CrossSection Weight (ub) : " << crossSectionWeight              << G4endl << G4endl;



  //G4cout << " =====> Storing HDC hit information into aHit" << G4endl;

    //   QweakSimHDCHit* aHit = new QweakSimHDCHit(MotherCopyNo); // there is only one plane per motherVolume
    QweakSimHDC_WirePlaneHit* aHit = new QweakSimHDC_WirePlaneHit(); // there is only one plane per motherVolume

    aHit->StoreHDCID(MotherCopyNo2);         // 0: Front HDC       , 1: Back HDC
    aHit->StoreWirePlaneID(WirePlaneCopyNo);    // [0,5] : WirePlane #0 to #5 (6 wire plane in total per HDC)

    aHit->StoreWorldPosition(worldPos);
    aHit->StoreLocalPosition(localPos);
    aHit->StoreGlobalTime(preStepPoint->GetGlobalTime());
    
    aHit->StoreOriginVertexPosition(originVertexPosition);
    aHit->StoreOriginVertexMomentumDirection(originVertexMomentumDirection);
    aHit->StoreOriginVertexKineticEnergy(originVertexKineticEnergy);

//     aHit->StorePrimaryQ2(primaryQ2);
//     aHit->StoreCrossSection(crossSection);
//     aHit->StoreCrossSectionWeight(crossSectionWeight);
//     aHit->StorePrimaryEventNumber(primaryEventNumber);

    aHit->StoreWorldMomentum(worldMomentum);
    aHit->StoreLocalMomentum(localMomentum);

    aHit->StoreParticleName(particlename);
    aHit->StoreParticleType(particletype);
    
    aHit->StoreKineticEnergy(currentKineticEnergy);
    aHit->StoreTotalEnergy(currentTotalEnergy);

    // now store all the hit information into the hit collection
    hitsCollection->insert(aHit);
  
  } // end of  if(strcmp(physVol->GetName(),"HDC_WirePlane_Physical")==0)


  //G4cout << G4endl << "###### Leaving QweakSimHDC_WirePlaneSD::ProcessHits() " << G4endl << G4endl; 

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void QweakSimHDC_WirePlaneSD::EndOfEvent(G4HCofThisEvent*  HCE)
void QweakSimHDC_WirePlaneSD::EndOfEvent(G4HCofThisEvent* )
{
  //G4cout << G4endl << "###### Calling QweakSimHDC_WirePlaneSD::EndOfEvent() " << G4endl << G4endl; 

//      G4int NbHits = hitsCollection->entries();
//       G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
//              << " hits in the Wire Planes: " << G4endl;
//       for (G4int i=0;i<NbHits;i++) (*hitsCollection)[i]->Print();

   //G4cout << G4endl << "###### Leaving QweakSimHDC_WirePlaneSD::EndOfEvent() " << G4endl << G4endl; 

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//=======================================================
//   -----------------------
//  | CVS File Information |
//  -----------------------
// 
//      $Revisions$  
//      $Log: QweakSimHDC_WirePlaneSD.cc,v $
//      Revision 1.3  2006/05/05 21:39:12  grimm
//      Records now the kinetic and total energy.
//
//      Revision 1.2  2005/12/27 19:12:16  grimm
//      - Redesign of Doxygen header containing CVS info like revision and date
//      - Added CVS revision log at the end of file
//
// 
