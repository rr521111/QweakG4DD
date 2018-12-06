//=============================================================================
//
//   ---------------------------
//  | Doxygen File Information |
//  ---------------------------
//
/**

   \file QweakSimEventAction.cc

   $Date: Fri Jul  3 11:38:14 CDT 2009 $

   \author Klaus Hans Grimm
   \author Jie Pan

*/
//=============================================================================
/*
Buddhiin 01-27-2015 - Strip down the root output for PMT DD analysis. We dont need;
Target, Region 1, Region 2, Region 3, Trigger Scintillator and Lumi

*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "QweakSimEventAction.hh"
#include "QweakSimEventActionMessenger.hh"

// geant4 includes
#include "G4Event.hh"
#include "G4UserEventAction.hh"
#include "G4SDManager.hh"

// user includes
#include "QweakSimAnalysis.hh"
#include "QweakSimUserInformation.hh"
//#include "QweakSimTarget_DetectorHit.hh"
//#include "QweakSimGEM_WirePlaneHit.hh"
//#include "QweakSimHDC_WirePlaneHit.hh"
//#include "QweakSimVDC_WirePlaneHit.hh"
//#include "QweakSimVDC_DriftCellHit.hh"
//#include "QweakSimTriggerScintillator_DetectorHit.hh"
//#include "QweakSimTriggerScintillator_PMTHit.hh"
// #include "QweakSimLeadGlass_DetectorHit.hh"
// #include "QweakSimLeadGlass_PMTHit.hh"
// #include "QweakSimPMTOnly_DetectorHit.hh"
// #include "QweakSimPMTOnly_PMTHit.hh"
#include "QweakSimCerenkov_DetectorHit.hh"
#include "QweakSimCerenkovDetector_PMTHit.hh"
//#include "QweakSimLumi_DetectorHit.hh"
#include "QweakSimTrajectory.hh"
#include "QweakSimUserMainEvent.hh"

#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
QweakSimEventAction::QweakSimEventAction(std::vector<double>*asInfo,QweakSimAnalysis* AN, QweakSimUserInformation* UI)
: analysis(AN),myUserInfo(UI),asymInfo(asInfo)
{
//---------------------------------------------------------------------------------------------
//! Constructor of QweakSimEventAction
    /*!

    \param  AN     class containing the Geamt4 hit data structure
    \param  myUI   class containing user information like Q2 for this event or QE of some PMTs
                                         which is needed for processing/saving hit information


    */
//---------------------------------------------------------------------------------------------


  //TargetDetector_CollID              = -1;
    //GEM_WirePlane_CollID               = -1;
    //HDC_WirePlane_CollID               = -1;
    //VDC_WirePlane_CollID               = -1;
    //VDC_DriftCellFront_CollID          = -1;
    //VDC_DriftCellBack_CollID           = -1;
    //TriggerScintillatorDetector_CollID = -1;
    //TriggerScintillatorPMT_CollID      = -1;
  //  LeadGlassDetector_CollID           = -1;
    //LeadGlassPMT_CollID                = -1;
    //PMTOnlyDetector_CollID	       = -1;
    //PMTOnlyDetectorPMT_CollID	       = -1;
    CerenkovDetector_CollID            = -1;
    CerenkovDetectorPMT_CollID         = -1;
    //LumiDetector_CollID                = -1;

    // Event action messenger
    fEventActionMessenger = new QweakSimEventActionMessenger(this);

    // Initialize map from string to trigger mode
    fTrigger.resize(kNumTriggers, false);
    fTriggerName.resize(kNumTriggers);
    if (kMapTriggerMode.size() == 0) {
        kMapTriggerMode["all"]   = kTriggerAll;
        fTriggerName[kTriggerAll] = "all";
        kMapTriggerMode["4fold"] = kTrigger4Fold;
        fTriggerName[kTrigger4Fold] = "4fold";
        kMapTriggerMode["3fold"] = kTrigger3Fold;
        fTriggerName[kTrigger3Fold] = "3fold";
        //kMapTriggerMode["scint"] = kTriggerScint;
	// fTriggerName[kTriggerScint] = "scint";
        //kMapTriggerMode["leadglass"] = kTriggerLeadGlass;  // trigger for the lead glass
        //fTriggerName[kTriggerLumi] = "lumi";
        //kMapTriggerMode["lumi"] = kTriggerLumi;  // trigger for the lead glass
	// fTriggerName[kTriggerLeadGlass] = "leadglass";
        // kMapTriggerMode["gem"]   = kTriggerGEM;
        // fTriggerName[kTriggerGEM] = "gem";
        kMapTriggerMode["cer"]   = kTriggerCer;
        fTriggerName[kTriggerCer] = "cer";
	//kMapTriggerMode["hdc"] = kTriggerHDC;
	//fTriggerName[kTriggerHDC] = "hdc";
	//	kMapTriggerMode["pmtonly"] = kTriggerPMTOnly;
	//fTriggerName[kTriggerHDC] = "pmtonly";
    }
    if (kMapTriggerMode.size() != kNumTriggers)
        G4cout << "Number of software triggers is not defined correctly!" << G4endl;

    // Initialize software trigger to false
    for (size_t iTrigger = 0; iTrigger < fTrigger.size(); iTrigger++)
        fTrigger[iTrigger] = false;

    // By default enable only cerenkov trigger
    fTrigger[kTriggerCer] = true;

    // By default, enable print out of hit information
    //printhits = true;
    printhits = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
QweakSimEventAction::~QweakSimEventAction(){
    // Delete the event action messenger
    if (fEventActionMessenger) delete fEventActionMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QweakSimEventAction::ShowTrigger()
{
    for (size_t iTrigger = 0; iTrigger < fTrigger.size(); iTrigger++)
        G4cout << (fTrigger[iTrigger]? "Enabled":"Disabled")
        << " software trigger " << fTriggerName[iTrigger] << "." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QweakSimEventAction::SetTrigger(const G4String value, const G4bool status)
{
    // No error checking here...
    //std::transform(value.begin(), value.end(), value.begin(), std::tolower);
    fTrigger[kMapTriggerMode[value]] = status;
    G4cout << (status? "Enabled":"Disabled") << " software trigger " << value << "." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QweakSimEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();

    // // check for existing PMTOnly Collection ID (if it's -1 it will be assigned)
    // if (TargetDetector_CollID==-1) {
    //     TargetDetector_CollID = SDman->GetCollectionID("TargetSD/TargetCollection");
    // }

    // check for existing GEM_WirePlane Collection ID (if it's -1 it will be assigned)
//     if (GEM_WirePlane_CollID==-1) {
//         GEM_WirePlane_CollID = SDman->GetCollectionID("GEMWirePlaneSD/GEMWirePlaneCollection");
//     }

    // // check for existing HDC_WirePlane Collection ID (if it's -1 it will be assigned)
    // if (HDC_WirePlane_CollID==-1) {
    //     HDC_WirePlane_CollID = SDman->GetCollectionID("HDCWirePlaneSD/HDCWirePlaneCollection");
    // }

    // // check for existing VDC_WirePlane Collection ID (if it's -1 it will be assigned)
    // if (VDC_WirePlane_CollID==-1) {
    //     VDC_WirePlane_CollID = SDman->GetCollectionID("VDCWirePlaneSD/VDCWirePlaneCollection");
    // }

    // // check for existing VDC_DriftCellFront Collection ID (if it's -1 it will be assigned)
    // if (VDC_DriftCellFront_CollID==-1) {
    //     VDC_DriftCellFront_CollID = SDman->GetCollectionID("VDCDriftCellFrontSD/DriftCellFrontCollection");
    // }

    // // check for existing VDC_DriftCellBack Collection ID (if it's -1 it will be assigned)
    // if (VDC_DriftCellBack_CollID==-1) {
    //     VDC_DriftCellBack_CollID = SDman->GetCollectionID("VDCDriftCellBackSD/DriftCellBackCollection");
    // }


    // // check for existing TriggerScintillator Collection ID (if it's -1 it will be assigned)
    // if (TriggerScintillatorDetector_CollID==-1) {
    //     TriggerScintillatorDetector_CollID = SDman->GetCollectionID("TriggerScintillatorSD/TriggerScintillatorCollection");
    // }

    // // check for existing LeadGlass Collection ID (if it's -1 it will be assigned)
    // if (LeadGlassDetector_CollID==-1) {
    //     LeadGlassDetector_CollID = SDman->GetCollectionID("LeadGlassSD/LeadGlassCollection");
    // }
    
    // // check for existing PMTOnly Collection ID (if it's -1 it will be assigned)
    // if (PMTOnlyDetector_CollID==-1) {
    //     PMTOnlyDetector_CollID = SDman->GetCollectionID("PMTOnlySD/PMTOnlyCollection");
    // }
    
    //     // check for existing PMTOnly Collection ID (if it's -1 it will be assigned)
    // if (PMTOnlyDetectorPMT_CollID==-1) {
    //     PMTOnlyDetectorPMT_CollID = SDman->GetCollectionID("PMTOnly_PMTSD/PMTHitCollection");
    // }
	
    // check for existing CerenkovDetector Collection ID (if it's -1 it will be assigned)
    if (CerenkovDetector_CollID==-1) {
        CerenkovDetector_CollID = SDman->GetCollectionID("CerenkovDetectorSD/CerenkovDetectorCollection");
    }

    // check for existing CerenkovDetectorPMT Collection ID (if it's -1 it will be assigned)
    if (CerenkovDetectorPMT_CollID==-1) {
        CerenkovDetectorPMT_CollID = SDman->GetCollectionID("CerenkovPMTSD/PMTHitCollection");
    }

    // // check for existing LumiDetector Collection ID (if it's -1 it will be assigned)
    // if (LumiDetector_CollID==-1) {
    //     // Do we want to change this so that both US and DS lumis have the same SD?
    //     LumiDetector_CollID = SDman->GetCollectionID("USLumiSD/LumiCollection");
    // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QweakSimEventAction::EndOfEventAction(const G4Event* evt) {

  const G4bool debugPrint=false;
  if(debugPrint)
    G4cout<<G4endl<<"Start EndOfEvenAction"<<G4endl;
  
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;  
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  if(debugPrint)
    G4cout<<" n_trajectories: "<<n_trajectories<<G4endl;

  //FIXME cg -- not sure why but this block crashes the code when running without optical photons
  // // extract the trajectories and draw them
  // if (G4VVisManager::GetConcreteInstance()) {
  //   G4cout<<" apparently there is a vis manager"<<G4endl;
  //   for (G4int i=0; i<n_trajectories; i++) {
  //     G4cout<<"   "<<i<<" "<<G4endl;
  //     if((*(evt->GetTrajectoryContainer()))[i]){
  // 	QweakSimTrajectory* trj = (QweakSimTrajectory*) ((*(evt->GetTrajectoryContainer()))[i]);
  // 	trj->DrawTrajectory(50);
  //     }
  //   }
  // }

  if(debugPrint)
    G4cout<<" done drawing trajectories"<<G4endl;

  //-----------------------------------------------------------------------------

  // preset variables for hit collection
  Initialize();    
  
  // Get event number
  G4int myEventCounter = myUserInfo->GetPrimaryEventNumber();


  // Get hit collection
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  
  // initialize HitsCollection pointers
  //QweakSimTarget_DetectorHitsCollection*              TargetDetector_HC              = 0;
  //QweakSimGEM_WirePlane_HitsCollection*               GEM_WirePlane_HC               = 0;
  // QweakSimHDC_WirePlane_HitsCollection*               HDC_WirePlane_HC               = 0;
  // QweakSimVDC_WirePlane_HitsCollection*               VDC_WirePlane_HC               = 0;
  // QweakSimVDC_DriftCellHitsCollection*                VDC_DriftCellFront_HC          = 0;
  // QweakSimVDC_DriftCellHitsCollection*                VDC_DriftCellBack_HC           = 0;
  // QweakSimTriggerScintillator_DetectorHitsCollection* TriggerScintillatorDetector_HC = 0;
  //QweakSimTriggerScintillator_PMTHitsCollection*      TriggerScintillatorPMT_HC      = 0;
  //QweakSimLeadGlass_DetectorHitsCollection*           LeadGlassDetector_HC           = 0;
  //QweakSimLeadGlass_PMTHitsCollection*                LeadGlassPMT_HC                = 0;
  //QweakSimPMTOnly_DetectorHitsCollection*             PMTOnlyDetector_HC             = 0;
  //QweakSimPMTOnly_PMTHitsCollection*                  PMTOnlyPMT_HC		       = 0;
  QweakSimCerenkovDetectorHitsCollection*             CerenkovDetector_HC            = 0;
  QweakSimCerenkovDetector_PMTHitsCollection*         CerenkovDetectorPMT_HC         = 0;
  //QweakSimLumi_DetectorHitsCollection*                LumiDetector_HC                = 0;
  
  //G4int n_hitTarget = 0;
  //G4int n_GEMhitWirePlane = 0;
  //G4int n_HDChitWirePlane = 0;
  // G4int n_VDChitWirePlane = 0;
  // G4int n_VDChitDCFront = 0;
  // G4int n_VDChitDCBack = 0;
  //G4int n_hitTriggerScintillator = 0;
  //G4int n_hitTriggerScintillatorPMT = 0;
  //G4int n_hitLeadGlass = 0;
  //G4int n_hitLeadGlassPMT = 0;
  //G4int n_hitPMTOnly = 0;
  //G4int n_hitPMTOnlyPMT = 0;
  G4int n_hitCerenkov = 0;
  G4int n_hitCerenkovPMT = 0;
  //G4int n_hitLumi = 0;
  //G4int n_hitLumiPMT = 0;
  
  if (HCE) {
    
    // // get  PMTOnly Hit Collector pointer
    // if (TargetDetector_CollID > -1) {
    //     TargetDetector_HC  = (QweakSimTarget_DetectorHitsCollection*)(HCE->GetHC(TargetDetector_CollID));
    //     n_hitTarget        = TargetDetector_HC -> entries();
    // }
    
    // get  GEM_WirePlane Hit Collector pointer
    //GEM_WirePlane_HC       = (QweakSimGEM_WirePlane_HitsCollection*)(HCE->GetHC(GEM_WirePlane_CollID));
    //n_GEMhitWirePlane      = GEM_WirePlane_HC -> entries();
    
    // // get  HDC_WirePlane Hit Collector pointer
    // if (HDC_WirePlane_CollID > -1) {
    //     HDC_WirePlane_HC       = (QweakSimHDC_WirePlane_HitsCollection*)(HCE->GetHC(HDC_WirePlane_CollID));
    //     n_HDChitWirePlane      = HDC_WirePlane_HC  -> entries();
    // }
    
    // // get  VDC_WirePlane Hit Collector pointer
    // if (VDC_WirePlane_CollID > -1) {
    //     VDC_WirePlane_HC       = (QweakSimVDC_WirePlane_HitsCollection*)(HCE->GetHC(VDC_WirePlane_CollID));
    //     n_VDChitWirePlane      = VDC_WirePlane_HC -> entries();
    // }
    
    // // get  VDC_DriftCellFront Hit Collector pointer
    // if (VDC_DriftCellFront_CollID > -1) {
    //     VDC_DriftCellFront_HC  = (QweakSimVDC_DriftCellHitsCollection*)(HCE->GetHC(VDC_DriftCellFront_CollID));
    //     n_VDChitDCFront        = VDC_DriftCellFront_HC -> entries();
    // }
    
    // // get  VDC_DriftCellFront Hit Collector pointer
    // if (VDC_DriftCellBack_CollID > -1) {
    //     VDC_DriftCellBack_HC   = (QweakSimVDC_DriftCellHitsCollection*)(HCE->GetHC(VDC_DriftCellBack_CollID));
    //     n_VDChitDCBack         = VDC_DriftCellBack_HC -> entries();
    // }
    
    // // get  TriggerScintillator Hit Collector pointer
    // if (TriggerScintillatorDetector_CollID > -1) {
    //     TriggerScintillatorDetector_HC = (QweakSimTriggerScintillator_DetectorHitsCollection*)(HCE->GetHC(TriggerScintillatorDetector_CollID));
    //     n_hitTriggerScintillator       = TriggerScintillatorDetector_HC -> entries();
    // }
    
    // get  TriggerScintillatorPMT Hit Collector pointer
    //if (TriggerScintillatorPMT_CollID > -1) {
    //    TriggerScintillatorPMT_HC   = (QweakSimTriggerScintillator_PMTHitsCollection*)(HCE->GetHC(TriggerScintillatorPMT_CollID));
    //    n_hitTriggerScintillatorPMT = TriggerScintillatorPMT_HC -> entries();
    //}
    
    // // get  LeadGlass Hit Collector pointer
    // if (LeadGlassDetector_CollID > -1) {
    //     LeadGlassDetector_HC  = (QweakSimLeadGlass_DetectorHitsCollection*)(HCE->GetHC(LeadGlassDetector_CollID));
    //     n_hitLeadGlass        = LeadGlassDetector_HC -> entries();
    // }
    
    // get  LeadGlassPMT Hit Collector pointer
    //if (LeadGlassPMT_CollID > -1) {
    //    LeadGlassPMT_HC    = (QweakSimLeadGlass_PMTHitsCollection*)(HCE->GetHC(LeadGlassPMT_CollID));
    //    n_hitLeadGlassPMT  = LeadGlassPMT_HC -> entries();
    //}
    
    // // get  PMTOnly Hit Collector pointer
    // if (PMTOnlyDetector_CollID > -1) {
    //     PMTOnlyDetector_HC  = (QweakSimPMTOnly_DetectorHitsCollection*)(HCE->GetHC(PMTOnlyDetector_CollID));
    //     n_hitPMTOnly        = PMTOnlyDetector_HC -> entries();
    // }
    
    // // get PMTOnly_PMT Hit Collection
    // if (PMTOnlyDetectorPMT_CollID > -1) {
    //     PMTOnlyPMT_HC  = (QweakSimPMTOnly_PMTHitsCollection*)(HCE->GetHC(PMTOnlyDetectorPMT_CollID));
    //     n_hitPMTOnlyPMT        = PMTOnlyPMT_HC -> entries();
    // }
    
    // get  CerenkovDetector Hit Collector pointer
    if (CerenkovDetector_CollID > -1) {
      CerenkovDetector_HC    = (QweakSimCerenkovDetectorHitsCollection*)(HCE->GetHC(CerenkovDetector_CollID));
      n_hitCerenkov          = CerenkovDetector_HC -> entries();
    }
    
    // get  CerenkovDetectorPMT Hit Collector pointer
    if (CerenkovDetectorPMT_CollID > -1) {
      CerenkovDetectorPMT_HC = (QweakSimCerenkovDetector_PMTHitsCollection*)(HCE->GetHC(CerenkovDetectorPMT_CollID));
      n_hitCerenkovPMT       = CerenkovDetectorPMT_HC -> entries();
    }
    
    // get  LumiDetector Hit Collector pointer
    // if (LumiDetector_CollID > -1) {
    //     LumiDetector_HC    = (QweakSimLumi_DetectorHitsCollection*)(HCE->GetHC(LumiDetector_CollID));
    //     n_hitLumi          = LumiDetector_HC -> entries();
    // }
  }
  
  if (printhits) {
    //G4cout << "Target " << n_hitTarget
    // << ",\tHDC " << n_HDChitWirePlane
    //        << ",\tVDC_Front " << n_VDChitDCFront
    //        << ",\tVDC_Back " << n_VDChitDCBack
    //        << ",\tTS " << n_hitTriggerScintillator
    G4cout //<< ",\tLeadGlass " << n_hitLeadGlass
      //<< ",\tPMTOnly " << n_hitPMTOnly
      //<< ",\tPMTOnlyPMT " << n_hitPMTOnlyPMT
      << ",\tCerenkov " << n_hitCerenkov
      //<< ",\tLumi " << n_hitLumi
      << "\tCerenkovPMT " << n_hitCerenkovPMT << G4endl;
  }
  
  
  // Initialize/Clear Event variables in target scattering window
  // analysis->fRootEvent->Target.Detector.Initialize();
  
  // Initialize/Clear Event variables, initialize Cerenkov Detector with NoHit Flag
  analysis->fRootEvent->Cerenkov.Detector.Initialize();
  analysis->fRootEvent->Cerenkov.PMT.Initialize();
  
  // // Initialize/Clear Event variables in Region 1
  // analysis->fRootEvent->Region1.ChamberFront.WirePlane.Initialize();
  // analysis->fRootEvent->Region1.ChamberBack.WirePlane.Initialize();
  
  // // Initialize/Clear Event variables in Region 2
  // analysis->fRootEvent->Region2.ChamberFront.WirePlane1.Initialize();
  // analysis->fRootEvent->Region2.ChamberFront.WirePlane2.Initialize();
  // analysis->fRootEvent->Region2.ChamberFront.WirePlane3.Initialize();
  // analysis->fRootEvent->Region2.ChamberFront.WirePlane4.Initialize();
  // analysis->fRootEvent->Region2.ChamberFront.WirePlane5.Initialize();
  // analysis->fRootEvent->Region2.ChamberFront.WirePlane6.Initialize();
  // //
  // analysis->fRootEvent->Region2.ChamberBack.WirePlane1.Initialize();
  // analysis->fRootEvent->Region2.ChamberBack.WirePlane2.Initialize();
  // analysis->fRootEvent->Region2.ChamberBack.WirePlane3.Initialize();
  // analysis->fRootEvent->Region2.ChamberBack.WirePlane4.Initialize();
  // analysis->fRootEvent->Region2.ChamberBack.WirePlane5.Initialize();
  // analysis->fRootEvent->Region2.ChamberBack.WirePlane6.Initialize();
  
  // // Initialize Region 3 wire planes (2: u,v ) with NoHit Flag
  // analysis->fRootEvent->Region3.ChamberFront.WirePlaneU.Initialize();
  // analysis->fRootEvent->Region3.ChamberFront.WirePlaneV.Initialize();
  // //
  // analysis->fRootEvent->Region3.ChamberBack.WirePlaneU.Initialize();
  // analysis->fRootEvent->Region3.ChamberBack.WirePlaneV.Initialize();
  
  // // Initialize DriftCells with NoHit Flag
  // analysis->fRootEvent->Region3.ChamberFront.DriftCell.Initialize();
  // //
  // analysis->fRootEvent->Region3.ChamberBack.DriftCell.Initialize();
  
  // Initialize TriggerScintillator with NoHit Flag
  // analysis->fRootEvent->TriggerScintillator.Detector.Initialize();
  
  // Initialize LeadGlass and PMTOnly //--- with NoHit Flag
  // analysis->fRootEvent->LeadGlass.Detector.Initialize();
  // analysis->fRootEvent->PMTOnly.Detector.Initialize();
  // analysis->fRootEvent->PMTOnly.PMT.Initialize();
  //analysis->fRootEvent->Lumi.Detector.Initialize();
  
  //#########################################################################################################################
  //#########################################################################################################################
  //
  //                                 ================================================================
  //                                  The Main "Software DAQ Trigger": setting the coincidence level
  //
  //                                  or: what is required for filling the Root ntuple for this event
  //                                 =================================================================
  //
  //##########################################################################################################################
  //##########################################################################################################################
  //
  if(debugPrint){
    G4cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ana: "<<" "<<myUserInfo->GetPrimaryEventNumber()<<" "
	  <<fTrigger[kTriggerAll]<<" "<<fTrigger[kTriggerCer]<<" "<<n_hitCerenkov<<G4endl;
  }
  //       <<" "<<fTrigger[kTriggerCer]<<" "<<n_hitCerenkov<<G4endl;//FIXME
  if ( fTrigger[kTriggerAll] /* Trigger on every event */
       //|| (fTrigger[kTrigger4Fold] && (n_VDChitWirePlane == 4) && (n_VDChitDCFront > 0) && (n_VDChitDCBack > 0) && (n_hitCerenkov > 0) ) /* 4-fold coincidence */
       // || (fTrigger[kTrigger3Fold] && (n_VDChitWirePlane >= 2) && (n_VDChitDCFront > 0) && (n_VDChitDCBack > 0) ) /* 3-fold coincidence */
       //|| (fTrigger[kTriggerScint] && (n_hitTriggerScintillator > 0) ) /* Qweak trigger on a hit in the trigger scintillator */
       // || (fTrigger[kTriggerHDC]   && (n_HDChitWirePlane >= 6))         /* HDC Trigger */
       //|| (fTrigger[kTriggerLeadGlass] && (n_hitLeadGlass >0))         /* a hit in the LeadGlass */
       // || (fTrigger[kTriggerPMTOnly] && (n_hitPMTOnly >0))         /* a hit in the PMTOnly */
       || (fTrigger[kTriggerCer]   && (n_hitCerenkov > 0) )          /* Triggering on Main Detector */
       //|| (fTrigger[kTriggerCer])           /* Triggering on Main Detector */	
       //|| (fTrigger[kTriggerLumi] && (n_hitLumi > 0))
       ) {
    
    //========================================
    // Store Primary Information into /Primary
    //========================================
    
    //-------------------------------------------------------------------------------------------
        G4int    PrimaryEventNumber = myUserInfo->GetPrimaryEventNumber();
        G4int    ReactionRegion = myUserInfo->GetReactionRegion();
        G4int    ReactionType = myUserInfo->GetReactionType();
        G4int    PDGcode = myUserInfo->GetPDGcode();
        G4int    TrackID = myUserInfo->GetTrackID();
        G4double GlobalTime = myUserInfo->GetGlobalTime();
        G4double PrimaryQ2 = myUserInfo->GetPrimaryQ2();
        G4double CrossSection = myUserInfo->GetCrossSection();
        G4double CrossSectionWeight = myUserInfo->GetCrossSectionWeight();
        G4double CrossSectionBornTotal     = myUserInfo->GetCrossSectionBornTotal();
        G4double CrossSectionBornInelastic = myUserInfo->GetCrossSectionBornInelastic();
        G4double CrossSectionBornQE        = myUserInfo->GetCrossSectionBornQE();
        G4double CrossSectionRadTotal      = myUserInfo->GetCrossSectionRadTotal();
        G4double CrossSectionRadElastic    = myUserInfo->GetCrossSectionRadElastic();
        G4double CrossSectionRadQE         = myUserInfo->GetCrossSectionRadQE();
        G4double CrossSectionRadDIS        = myUserInfo->GetCrossSectionRadDIS();
        G4double CrossSectionRadTotalIntOnly   = myUserInfo->GetCrossSectionRadTotalIntOnly();
        G4double CrossSectionRadElasticIntOnly = myUserInfo->GetCrossSectionRadElasticIntOnly();
        G4double CrossSectionRadQEIntOnly      = myUserInfo->GetCrossSectionRadQEIntOnly();
        G4double CrossSectionRadDISIntOnly     = myUserInfo->GetCrossSectionRadDISIntOnly();
        G4double CrossSectionRadElasticPeak    = myUserInfo->GetCrossSectionRadElasticPeak();
	G4double Asymmetry = myUserInfo->GetAsymmetry();
        G4double OriginVertexPositionX = myUserInfo->GetOriginVertexPositionX();
        G4double OriginVertexPositionY = myUserInfo->GetOriginVertexPositionY();
        G4double OriginVertexPositionZ = myUserInfo->GetOriginVertexPositionZ();
        G4double OriginVertexThetaAngle = myUserInfo->GetOriginVertexThetaAngle();
        G4double OriginVertexPhiAngle = myUserInfo->GetOriginVertexPhiAngle();
        G4double OriginVertexMomentumDirectionX = myUserInfo->GetOriginVertexMomentumDirectionX();
        G4double OriginVertexMomentumDirectionY = myUserInfo->GetOriginVertexMomentumDirectionY();
        G4double OriginVertexMomentumDirectionZ = myUserInfo->GetOriginVertexMomentumDirectionZ();
	G4double PreScatteringKineticEnergy = myUserInfo->GetPreScatteringKineticEnergy();
        G4double OriginVertexKineticEnergy = myUserInfo->GetOriginVertexKineticEnergy();
        G4double OriginVertexTotalEnergy = myUserInfo->GetOriginVertexTotalEnergy();

        G4int NumberOfEventToBeProcessed = myUserInfo->GetNumberOfEventToBeProcessed();
        G4double PhiAngle_Min            = myUserInfo->GetPhiAngle_Min();
        G4double PhiAngle_Max            = myUserInfo->GetPhiAngle_Max();
        G4double ThetaAngle_Min          = myUserInfo->GetThetaAngle_Min();
        G4double ThetaAngle_Max          = myUserInfo->GetThetaAngle_Max();
        G4double EPrime_Min              = myUserInfo->GetEPrime_Min();
        G4double EPrime_Max              = myUserInfo->GetEPrime_Max();
        G4double BeamEnergy              = myUserInfo->GetBeamEnergy();
        G4double Luminosity              = myUserInfo->GetLuminosity();
        G4double PhaseSpace              = myUserInfo->GetPhaseSpace();
        CalculateKinematicVariables();
        G4double OriginVertexKinematicNu = myUserInfo->GetOriginVertexKinematicNu();
        G4double OriginVertexKinematicQ2 = myUserInfo->GetOriginVertexKinematicQ2();
        G4double OriginVertexKinematicX  = myUserInfo->GetOriginVertexKinematicX();
        G4double OriginVertexKinematicW  = myUserInfo->GetOriginVertexKinematicW();
        G4double EffectiveKinematicNu    = myUserInfo->GetEffectiveKinematicNu();
        G4double EffectiveKinematicQ2    = myUserInfo->GetEffectiveKinematicQ2();
        G4double EffectiveKinematicX     = myUserInfo->GetEffectiveKinematicX();
        G4double EffectiveKinematicW     = myUserInfo->GetEffectiveKinematicW();

        analysis->fRootEvent->Primary.StoreTrackID((Int_t) TrackID);
        analysis->fRootEvent->Primary.StoreGlobalTime((Float_t) GlobalTime);
        analysis->fRootEvent->Primary.StoreOriginVertexPositionX((Float_t) OriginVertexPositionX);
        analysis->fRootEvent->Primary.StoreOriginVertexPositionY((Float_t) OriginVertexPositionY);
        analysis->fRootEvent->Primary.StoreOriginVertexPositionZ((Float_t) OriginVertexPositionZ);
        analysis->fRootEvent->Primary.StoreOriginVertexMomentumDirectionX((Float_t) OriginVertexMomentumDirectionX);
        analysis->fRootEvent->Primary.StoreOriginVertexMomentumDirectionY((Float_t) OriginVertexMomentumDirectionY);
        analysis->fRootEvent->Primary.StoreOriginVertexMomentumDirectionZ((Float_t) OriginVertexMomentumDirectionZ);
        analysis->fRootEvent->Primary.StoreOriginVertexThetaAngle((Float_t) OriginVertexThetaAngle);
        analysis->fRootEvent->Primary.StoreOriginVertexPhiAngle((Float_t) OriginVertexPhiAngle);
	analysis->fRootEvent->Primary.StorePreScatteringKineticEnergy((Float_t) PreScatteringKineticEnergy);
        analysis->fRootEvent->Primary.StoreOriginVertexKineticEnergy((Float_t) OriginVertexKineticEnergy);
        analysis->fRootEvent->Primary.StoreOriginVertexTotalEnergy((Float_t) OriginVertexTotalEnergy);
        analysis->fRootEvent->Primary.StorePrimaryQ2((Float_t) PrimaryQ2);
        analysis->fRootEvent->Primary.StoreCrossSection((Float_t) CrossSection);
        analysis->fRootEvent->Primary.StoreCrossSectionWeight((Float_t) CrossSectionWeight);
        analysis->fRootEvent->Primary.StoreCrossSectionBornTotal    ((Float_t) CrossSectionBornTotal);
        analysis->fRootEvent->Primary.StoreCrossSectionBornInelastic((Float_t) CrossSectionBornInelastic);
        analysis->fRootEvent->Primary.StoreCrossSectionBornQE       ((Float_t) CrossSectionBornQE);
        analysis->fRootEvent->Primary.StoreCrossSectionRadTotal     ((Float_t) CrossSectionRadTotal);
        analysis->fRootEvent->Primary.StoreCrossSectionRadElastic   ((Float_t) CrossSectionRadElastic);
        analysis->fRootEvent->Primary.StoreCrossSectionRadQE        ((Float_t) CrossSectionRadQE);
        analysis->fRootEvent->Primary.StoreCrossSectionRadDIS       ((Float_t) CrossSectionRadDIS);
        analysis->fRootEvent->Primary.StoreCrossSectionRadTotalIntOnly  ((Float_t) CrossSectionRadTotalIntOnly);
        analysis->fRootEvent->Primary.StoreCrossSectionRadElasticIntOnly((Float_t) CrossSectionRadElasticIntOnly);
        analysis->fRootEvent->Primary.StoreCrossSectionRadQEIntOnly     ((Float_t) CrossSectionRadQEIntOnly);
        analysis->fRootEvent->Primary.StoreCrossSectionRadDISIntOnly    ((Float_t) CrossSectionRadDISIntOnly);
        analysis->fRootEvent->Primary.StoreCrossSectionRadElasticPeak   ((Float_t) CrossSectionRadElasticPeak);
	analysis->fRootEvent->Primary.StoreAsymmetry((Float_t) Asymmetry);
        analysis->fRootEvent->Primary.StorePrimaryEventNumber((Int_t) PrimaryEventNumber);
        analysis->fRootEvent->Primary.StoreReactionRegion((Int_t) ReactionRegion);
        analysis->fRootEvent->Primary.StoreReactionType((Int_t) ReactionType);
        analysis->fRootEvent->Primary.StorePDGcode((Int_t) PDGcode);
        analysis->fRootEvent->Primary.StoreNumberOfEventToBeProcessed((Int_t) NumberOfEventToBeProcessed);
        analysis->fRootEvent->Primary.StorePhiAngle_Min((Float_t) PhiAngle_Min/degree);
        analysis->fRootEvent->Primary.StorePhiAngle_Max((Float_t) PhiAngle_Max/degree);
        analysis->fRootEvent->Primary.StoreThetaAngle_Min((Float_t) ThetaAngle_Min/degree);
        analysis->fRootEvent->Primary.StoreThetaAngle_Max((Float_t) ThetaAngle_Max/degree);
        analysis->fRootEvent->Primary.StoreEPrime_Min((Float_t) EPrime_Min);
        analysis->fRootEvent->Primary.StoreEPrime_Max((Float_t) EPrime_Max);
        analysis->fRootEvent->Primary.StoreBeamEnergy((Float_t) BeamEnergy);
        analysis->fRootEvent->Primary.StoreLuminosity((Float_t) Luminosity);
        analysis->fRootEvent->Primary.StorePhaseSpace((Float_t) PhaseSpace);
        analysis->fRootEvent->Primary.StoreOriginVertexKinematicNu((Float_t) OriginVertexKinematicNu);
        analysis->fRootEvent->Primary.StoreOriginVertexKinematicQ2((Float_t) OriginVertexKinematicQ2);
        analysis->fRootEvent->Primary.StoreOriginVertexKinematicX((Float_t) OriginVertexKinematicX);
        analysis->fRootEvent->Primary.StoreOriginVertexKinematicW((Float_t) OriginVertexKinematicW);
        analysis->fRootEvent->Primary.StoreEffectiveKinematicNu((Float_t) EffectiveKinematicNu);
        analysis->fRootEvent->Primary.StoreEffectiveKinematicQ2((Float_t) EffectiveKinematicQ2);
        analysis->fRootEvent->Primary.StoreEffectiveKinematicX((Float_t) EffectiveKinematicX);
        analysis->fRootEvent->Primary.StoreEffectiveKinematicW((Float_t) EffectiveKinematicW);

	////////
	// store energy loss variables into rootfile
	//  all Elosses are in MeV (see QweakSimSteppingAction.cc)
	analysis->fRootEvent->Primary.StoredEIonIn((Float_t) myUserInfo->GetdEIonIn());
	analysis->fRootEvent->Primary.StoredEIonOut((Float_t) myUserInfo->GetdEIonOut());
	analysis->fRootEvent->Primary.StoredEIonTot((Float_t) myUserInfo->GetdEIonTot());
	analysis->fRootEvent->Primary.StoredEBremIn((Float_t) myUserInfo->GetdEBremIn());
	analysis->fRootEvent->Primary.StoredEBremOut((Float_t) myUserInfo->GetdEBremOut());
	analysis->fRootEvent->Primary.StoredEBremTot((Float_t) myUserInfo->GetdEBremTot());
	analysis->fRootEvent->Primary.StoredEMscIn((Float_t) myUserInfo->GetdEMscIn());
	analysis->fRootEvent->Primary.StoredEMscOut((Float_t) myUserInfo->GetdEMscOut());
	analysis->fRootEvent->Primary.StoredEMscTot((Float_t) myUserInfo->GetdEMscTot());
	analysis->fRootEvent->Primary.StoredETotIn((Float_t) myUserInfo->GetdETotIn());
	analysis->fRootEvent->Primary.StoredETotOut((Float_t) myUserInfo->GetdETotOut());
	analysis->fRootEvent->Primary.StoredETot((Float_t) myUserInfo->GetdETot());
	///////

	///////
	// store asymmetry information
	Double_t pp=asymInfo->at(0);
	Double_t pm=asymInfo->at(1);
	analysis->fRootEvent->Primary.StoreAsymPrim((pp-pm)/(pp+pm));
	analysis->fRootEvent->Primary.StoreAsymPlus(pp);
	analysis->fRootEvent->Primary.StoreAsymNomi(pp-pm);
	analysis->fRootEvent->Primary.StoreAsymDeno(pp+pm);
	///////
	
        //==========================================================================================

        //===========================================
        // Store Number Of Hits of each Detector
        //===========================================

        // Store Number of Hits for: Cerenkov Detector
        analysis->fRootEvent->Cerenkov.Detector.StoreDetectorNbOfHits(n_hitCerenkov);
		
        // Store Number of Hits for: Cerenkov PMT Detector
        analysis->fRootEvent->Cerenkov.PMT.StoreDetectorNbOfHits(n_hitCerenkovPMT);



        //===============================================================================================================

        //=========================================================
        // Store Cerenkov Detector hits into /Cerenkov
        //=========================================================

        if (n_hitCerenkov > 0) {

            // loop over hits
            for (int i1 = 0; i1 < n_hitCerenkov; i1++) {

                QweakSimCerenkov_DetectorHit* aHit = (*CerenkovDetector_HC)[i1];

                G4int octantID = (Int_t) aHit->GetDetectorID() + 1;

                if (print_Cerenkov_DetectorHit) aHit->Print();

                // get local position of hit
                G4ThreeVector localPosition  = aHit->GetLocalPosition();
                Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
                Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
                Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;

                // get world position of hit
                G4ThreeVector globalPosition  = aHit->GetWorldPosition();
                Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
                Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
                Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;

                // get local Momentum of hit
                G4ThreeVector localMomentum  = aHit->GetLocalMomentum();
                Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
                Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
                Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;
                Float_t rLocalPhiAngle   = (Float_t) localMomentum.phi() / degree - 90.0;
                Float_t rLocalThetaAngle = (Float_t) localMomentum.theta() / degree;

                // get world Momentum of hit
                G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
                Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
                Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
                Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;
                Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;
                Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;

                G4ThreeVector localExitPosition = myUserInfo->GetLocalCerenkovExitPosition();
                Float_t rLocalExitPositionX = (Float_t) localExitPosition.x() / cm;
                Float_t rLocalExitPositionY = (Float_t) localExitPosition.y() / cm;
                Float_t rLocalExitPositionZ = (Float_t) localExitPosition.z() / cm;

                G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
                Float_t rOriginVertexPositionX = (Float_t) originVertexPosition.x() / cm;
                Float_t rOriginVertexPositionY = (Float_t) originVertexPosition.y() / cm;
                Float_t rOriginVertexPositionZ = (Float_t) originVertexPosition.z() / cm;

                G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
                Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
                Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
                Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
                Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
                Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;

                Float_t rOriginVertexKineticEnergy = (Float_t) aHit->GetOriginVertexKineticEnergy() / MeV;
                Float_t rOriginVertexTotalEnergy   = (Float_t) aHit->GetOriginVertexTotalEnergy() / MeV;

                Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;

                TString rParticleName = TString(aHit->GetParticleName());
                G4int rParticleType = (Int_t) aHit->GetParticleType();
		G4int rTrackId = (Int_t) aHit->GetTrackID();
		G4int rParentId = (Int_t) aHit->GetParentID();
		TString rCreatorName = TString(aHit->GetCreatorProcessName());
		
                // get total Energy of hit
                Float_t rTotalEnergy     = (Float_t) aHit->GetTotalEnergy() / MeV;

                // get kinetic Energy of hit
                Float_t rKineticEnergy     = (Float_t) aHit->GetKineticEnergy() / MeV;

                // get polarization of hit
                G4ThreeVector rGlobalPolarizationVector = aHit->GetPolarization();
                Float_t rGlobalPolarizationX = rGlobalPolarizationVector.x();
                Float_t rGlobalPolarizationY = rGlobalPolarizationVector.y();
                Float_t rGlobalPolarizationZ = rGlobalPolarizationVector.z();
                // longitudinal polarization is along the momentum direction
                Float_t rLongitudinalPolarization (0);
 		if(globalMomentum.mag()>0 && rGlobalPolarizationVector.mag()>0)
 		  rLongitudinalPolarization = rGlobalPolarizationVector.dot(globalMomentum)/
 		    globalMomentum.mag()/rGlobalPolarizationVector.mag();
                // transverse polarization is all except longitudinal
		G4ThreeVector rTransversePolarizationVector(0);
		if(globalMomentum.mag()>0)
		  rTransversePolarizationVector = rGlobalPolarizationVector -
		    rLongitudinalPolarization*globalMomentum/globalMomentum.mag();
                Float_t rTransversePolarization = rTransversePolarizationVector.mag();
                Float_t rTransversePolarizationX = rTransversePolarizationVector.x();
                Float_t rTransversePolarizationY = rTransversePolarizationVector.y();
                Float_t rTransversePolarizationZ = rTransversePolarizationVector.z();
                Float_t rTransversePolarizationPhiAngle = rTransversePolarizationVector.phi() / degree;

                //  edgeEvent = myUserInfo->GetEdgeEventDetected();

                //==========================================================
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorID(octantID);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorHasBeenHit(5);
                analysis->fRootEvent->Cerenkov.Detector.StoreParticleName(rParticleName);
                analysis->fRootEvent->Cerenkov.Detector.StoreParticleType(rParticleType);
		analysis->fRootEvent->Cerenkov.Detector.StoreParticleID(rTrackId);
                analysis->fRootEvent->Cerenkov.Detector.StoreParentID(rParentId);
		analysis->fRootEvent->Cerenkov.Detector.StoreCreatorProcessName(rCreatorName);
		analysis->fRootEvent->Cerenkov.Detector.StoreGlobalTimeOfHit(rGlobalTime);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexPositionX(rOriginVertexPositionX);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexPositionY(rOriginVertexPositionY);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexPhiAngle(rOriginVertexPhiAngle);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexThetaAngle(rOriginVertexThetaAngle);

                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
                analysis->fRootEvent->Cerenkov.Detector.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);

                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorLocalPositionX(rLocalPositionX);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorLocalPositionY(rLocalPositionY);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorLocalPositionZ(rLocalPositionZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorLocalExitPositionX(rLocalExitPositionX);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorLocalExitPositionY(rLocalExitPositionY);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorLocalExitPositionZ(rLocalExitPositionZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorGlobalPositionX(rGlobalPositionX);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorGlobalPositionY(rGlobalPositionY);
                analysis->fRootEvent->Cerenkov.Detector.StoreDetectorGlobalPositionZ(rGlobalPositionZ);

                // store polarization of hit
                analysis->fRootEvent->Cerenkov.Detector.StorePolarizationX(rGlobalPolarizationX);
                analysis->fRootEvent->Cerenkov.Detector.StorePolarizationY(rGlobalPolarizationY);
                analysis->fRootEvent->Cerenkov.Detector.StorePolarizationZ(rGlobalPolarizationZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreLongitudinalPolarization(rLongitudinalPolarization);
                analysis->fRootEvent->Cerenkov.Detector.StoreTransversePolarization(rTransversePolarization);
                analysis->fRootEvent->Cerenkov.Detector.StoreTransversePolarizationX(rTransversePolarizationX);
                analysis->fRootEvent->Cerenkov.Detector.StoreTransversePolarizationY(rTransversePolarizationY);
                analysis->fRootEvent->Cerenkov.Detector.StoreTransversePolarizationZ(rTransversePolarizationZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreTransversePolarizationPhiAngle(rTransversePolarizationPhiAngle);

                // store local momentum and track angles Phi and Theta
                analysis->fRootEvent->Cerenkov.Detector.StoreLocalMomentumX(rLocalMomentumX);
                analysis->fRootEvent->Cerenkov.Detector.StoreLocalMomentumY(rLocalMomentumY);
                analysis->fRootEvent->Cerenkov.Detector.StoreLocalMomentumZ(rLocalMomentumZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreLocalPhiAngle(rLocalPhiAngle);
                analysis->fRootEvent->Cerenkov.Detector.StoreLocalThetaAngle(rLocalThetaAngle);

                // store global momentum and track angles Phi and Theta
                analysis->fRootEvent->Cerenkov.Detector.StoreGlobalMomentumX(rGlobalMomentumX);
                analysis->fRootEvent->Cerenkov.Detector.StoreGlobalMomentumY(rGlobalMomentumY);
                analysis->fRootEvent->Cerenkov.Detector.StoreGlobalMomentumZ(rGlobalMomentumZ);
                analysis->fRootEvent->Cerenkov.Detector.StoreGlobalPhiAngle(rGlobalPhiAngle);
                analysis->fRootEvent->Cerenkov.Detector.StoreGlobalThetaAngle(rGlobalThetaAngle);

                // store total+kinetic energy of a hit
                analysis->fRootEvent->Cerenkov.Detector.StoreTotalEnergy(rTotalEnergy);
                analysis->fRootEvent->Cerenkov.Detector.StoreKineticEnergy(rKineticEnergy);

                Int_t edgeEvent;
                if (GetDistance(localPosition,localExitPosition)/cm < 1.15)
                    edgeEvent = 1;
                else
                    edgeEvent = 0;

                analysis->fRootEvent->Cerenkov.Detector.StoreEdgeEventFlag(edgeEvent);
            } // end  for(int i1=0;i1<n_hitCerenkov;i1++)
        } // end    if (n_hitCerenkov >0)

        //=========================================================
        // Store Number of Photoelectrons of Cerenkov Detector hits
        //=========================================================

        if (n_hitCerenkovPMT > 0) {

            std::vector<G4int>   PmtHasBeenHit(PmtMaxSize);

            std::vector<G4float> PmtTime(PmtMaxSize);
            std::vector<G4float> PmtEnergy(PmtMaxSize);
            std::vector<G4int>   PmtOctant(PmtMaxSize);

            std::vector<G4int>   PmtHitsLeft(PmtMaxSize);
            std::vector<G4int>   PmtHitsRight(PmtMaxSize);
            std::vector<G4int>   PmtHitsTotal(PmtMaxSize);
            std::vector<G4float> PmtNPELeft(PmtMaxSize);
            std::vector<G4float> PmtNPERight(PmtMaxSize);
            std::vector<G4float> PmtNPETotal(PmtMaxSize);

            std::vector<G4float> PmtRateLeft(PmtMaxSize);
            std::vector<G4float> PmtRateRight(PmtMaxSize);
            std::vector<G4float> PmtRateTotal(PmtMaxSize);
            std::vector<G4float> PmtRateLeftEL(PmtMaxSize);
            std::vector<G4float> PmtRateRightEL(PmtMaxSize);
            std::vector<G4float> PmtRateTotalEL(PmtMaxSize);
            std::vector<G4float> PmtRateLeftDIS(PmtMaxSize);
            std::vector<G4float> PmtRateRightDIS(PmtMaxSize);
            std::vector<G4float> PmtRateTotalDIS(PmtMaxSize);
            std::vector<G4float> PmtRateLeftQE(PmtMaxSize);
            std::vector<G4float> PmtRateRightQE(PmtMaxSize);
            std::vector<G4float> PmtRateTotalQE(PmtMaxSize);
            std::vector<G4float> PmtRateLeftELPeak(PmtMaxSize);
            std::vector<G4float> PmtRateRightELPeak(PmtMaxSize);
            std::vector<G4float> PmtRateTotalELPeak(PmtMaxSize);

            std::vector<G4float> PmtYieldLeft(PmtMaxSize);
            std::vector<G4float> PmtYieldRight(PmtMaxSize);
            std::vector<G4float> PmtYieldTotal(PmtMaxSize);
            std::vector<G4float> PmtYieldLeftEL(PmtMaxSize);
            std::vector<G4float> PmtYieldRightEL(PmtMaxSize);
            std::vector<G4float> PmtYieldTotalEL(PmtMaxSize);
            std::vector<G4float> PmtYieldLeftDIS(PmtMaxSize);
            std::vector<G4float> PmtYieldRightDIS(PmtMaxSize);
            std::vector<G4float> PmtYieldTotalDIS(PmtMaxSize);
            std::vector<G4float> PmtYieldLeftQE(PmtMaxSize);
            std::vector<G4float> PmtYieldRightQE(PmtMaxSize);
            std::vector<G4float> PmtYieldTotalQE(PmtMaxSize);
            std::vector<G4float> PmtYieldLeftELPeak(PmtMaxSize);
            std::vector<G4float> PmtYieldRightELPeak(PmtMaxSize);
            std::vector<G4float> PmtYieldTotalELPeak(PmtMaxSize);

            // loop over hits
            for (int i1 = 0; i1 < n_hitCerenkovPMT; i1++) {

                QweakSimCerenkovDetector_PMTHit* aHit = (*CerenkovDetectorPMT_HC)[i1];
                G4int octantID = (Int_t) aHit->GetDetectorID() + 1;

                if (octantID > (Int_t) PmtHitsLeft.size()) {
                  G4cerr << "octantID is larger than size of vectorPmtHitsLeft in QweakSimEventAction!" << G4endl;
                  break;
                }

                //------------------------------------------------------------------------
                PmtTime.push_back(aHit->GetHitTime());
                PmtEnergy.push_back(aHit->GetPhotonEnergy());
                PmtOctant.push_back(0);

                //------------------------------------------------------------------------
                if (aHit->GetPMTID() == 0) { // left PMT
                    //if(aHit->IsHitValid())
                    {
                        PmtOctant.back() = -octantID;
                        PmtHitsLeft[octantID]++;
                        PmtNPELeft[octantID] += myUserInfo->GetNumberOfPhotoelectronsS20(aHit->GetPhotonEnergy()*1.0e6);
			//G4cout<<"   pmtNPELeft: "<<i1<<" "<<aHit->GetPhotonEnergy()<<" "<<myUserInfo->GetNumberOfPhotoelectronsS20(aHit->GetPhotonEnergy()*1.0e6)<<" "<<PmtNPELeft[octantID]<<" "<<PmtHitsLeft[octantID]<<G4endl;
                    }
                }

                if (aHit->GetPMTID() == 1) { // right PMT
                    //if(aHit->IsHitValid())
                    {
                        PmtOctant.back() = +octantID;
                        PmtHitsRight[octantID]++;
                        PmtNPERight[octantID] += myUserInfo->GetNumberOfPhotoelectronsS20(aHit->GetPhotonEnergy()*1.0e6);
                        // G4cout<<"pmtNPERight: "<<pmtNPERight[octantID]<<G4endl;
                    }
                }

                PmtHasBeenHit[octantID] = 5;
                PmtHitsTotal[octantID]++;
                PmtNPETotal[octantID] += myUserInfo->GetNumberOfPhotoelectronsS20(aHit->GetPhotonEnergy()*1.0e6);
                // G4cout<<"pmtNPETTotal: "<<pmtNPETotal[octantID]<<G4endl;

                //------------------------------------------------------------------------

            } // end for (int i1 = 0; i1 < n_hitCerenkovPMT; i1++)

            // Loop over all octants
            for (int octantID = 0; octantID < PmtMaxSize; octantID++) {
                if (PmtHasBeenHit[octantID] != 5) continue; // skip octants without hits

                PmtRateTotal[octantID] = CalculateRate( myUserInfo->GetCrossSection(), (Bool_t)PmtNPETotal[octantID] );
                PmtRateLeft[octantID] = CalculateRate( myUserInfo->GetCrossSection(), (Bool_t)PmtNPELeft[octantID] );
                PmtRateRight[octantID] = CalculateRate( myUserInfo->GetCrossSection(), (Bool_t)PmtNPERight[octantID] );

                PmtRateTotalEL[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElastic(), (Bool_t)PmtNPETotal[octantID] );
                PmtRateLeftEL[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElastic(), (Bool_t)PmtNPELeft[octantID] );
                PmtRateRightEL[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElastic(), (Bool_t)PmtNPERight[octantID] );

                PmtRateTotalDIS[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadDIS(), (Bool_t)PmtNPETotal[octantID] );
                PmtRateLeftDIS[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadDIS(), (Bool_t)PmtNPELeft[octantID] );
                PmtRateRightDIS[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadDIS(), (Bool_t)PmtNPERight[octantID] );

                PmtRateTotalQE[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadQE(), (Bool_t)PmtNPETotal[octantID] );
                PmtRateLeftQE[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadQE(), (Bool_t)PmtNPELeft[octantID] );
                PmtRateRightQE[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadQE(), (Bool_t)PmtNPERight[octantID] );

                PmtRateTotalELPeak[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElasticPeak(), (Bool_t)PmtNPETotal[octantID] );
                PmtRateLeftELPeak[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElasticPeak(), (Bool_t)PmtNPELeft[octantID] );
                PmtRateRightELPeak[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElasticPeak(), (Bool_t)PmtNPERight[octantID] );

                PmtYieldTotal[octantID] = CalculateRate( myUserInfo->GetCrossSection(), PmtNPETotal[octantID] );
                PmtYieldLeft[octantID] = CalculateRate( myUserInfo->GetCrossSection(), PmtNPELeft[octantID] );
                PmtYieldRight[octantID] = CalculateRate( myUserInfo->GetCrossSection(), PmtNPERight[octantID] );

                PmtYieldTotalEL[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElastic(), PmtNPETotal[octantID] );
                PmtYieldLeftEL[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElastic(), PmtNPELeft[octantID] );
                PmtYieldRightEL[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElastic(), PmtNPERight[octantID] );

                PmtYieldTotalDIS[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadDIS(), PmtNPETotal[octantID] );
                PmtYieldLeftDIS[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadDIS(), PmtNPELeft[octantID] );
                PmtYieldRightDIS[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadDIS(), PmtNPERight[octantID] );

                PmtYieldTotalQE[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadQE(), PmtNPETotal[octantID] );
                PmtYieldLeftQE[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadQE(), PmtNPELeft[octantID] );
                PmtYieldRightQE[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadQE(), PmtNPERight[octantID] );

                PmtYieldTotalELPeak[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElasticPeak(), PmtNPETotal[octantID] );
                PmtYieldLeftELPeak[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElasticPeak(), PmtNPELeft[octantID] );
                PmtYieldRightELPeak[octantID] = CalculateRate( myUserInfo->GetCrossSectionRadElasticPeak(), PmtNPERight[octantID] );

            } // end for (int octantID = 0; octantID < PmtMaxSize; octantID++)

            // has the detector been hit
            analysis->fRootEvent->Cerenkov.PMT.StorePMTHasBeenHit(PmtHasBeenHit);

            //---------------------------------------------
            // store number of hits for left and right PMT
            //---------------------------------------------
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTimeOfHits(PmtTime);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTEnergyOfHits(PmtEnergy);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTOctantOfHits(PmtOctant);

            //---------------------------------------------
            // store number of hits for left and right PMT
            //---------------------------------------------
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftNbOfHits(PmtHitsLeft);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightNbOfHits(PmtHitsRight);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalNbOfHits(PmtHitsTotal);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftNbOfPEs(PmtNPELeft);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightNbOfPEs(PmtNPERight);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalNbOfPEs(PmtNPETotal);

            //---------------------------------------------
            // store rates for left and right PMT
            //---------------------------------------------
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftRate(PmtRateLeft);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightRate(PmtRateRight);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalRate(PmtRateTotal);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftRateEL(PmtRateLeftEL);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightRateEL(PmtRateRightEL);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalRateEL(PmtRateTotalEL);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftRateDIS(PmtRateLeftDIS);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightRateDIS(PmtRateRightDIS);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalRateDIS(PmtRateTotalDIS);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftRateQE(PmtRateLeftQE);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightRateQE(PmtRateRightQE);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalRateQE(PmtRateTotalQE);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftRateELPeak(PmtRateLeftELPeak);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightRateELPeak(PmtRateRightELPeak);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalRateELPeak(PmtRateTotalELPeak);

            //---------------------------------------------
            // store yields for left and right PMT
            //---------------------------------------------
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftYield(PmtYieldLeft);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightYield(PmtYieldRight);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalYield(PmtYieldTotal);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftYieldEL(PmtYieldLeftEL);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightYieldEL(PmtYieldRightEL);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalYieldEL(PmtYieldTotalEL);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftYieldDIS(PmtYieldLeftDIS);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightYieldDIS(PmtYieldRightDIS);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalYieldDIS(PmtYieldTotalDIS);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftYieldQE(PmtYieldLeftQE);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightYieldQE(PmtYieldRightQE);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalYieldQE(PmtYieldTotalQE);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTLeftYieldELPeak(PmtYieldLeftELPeak);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTRightYieldELPeak(PmtYieldRightELPeak);
            analysis->fRootEvent->Cerenkov.PMT.StorePMTTotalYieldELPeak(PmtYieldTotalELPeak);

            //==============================================================================

        } // end if (n_hitCerenkovPMT > 0)

        // //==============================
        // // Store HDC hits into /Region2
        // //==============================

        // if (n_HDChitWirePlane > 0) {

        //     // loop over wire plane hits
        //     for (int i1=0;i1<n_HDChitWirePlane;i1++) {

        //         int HDC_Chamber_Plane_NbOfHits[2][6];
        //         for (int chamber = 0; chamber < 2; chamber++)
        //           for (int plane = 0; plane < 6; plane++)
        //             HDC_Chamber_Plane_NbOfHits[chamber][plane] = 0;

        //         // get hit pointer for each hit
        //         QweakSimHDC_WirePlaneHit* aHit = (*HDC_WirePlane_HC)[i1];

        //         // get local position of hit
        //         G4ThreeVector localPosition  = aHit->GetLocalPosition();
        //         Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
        //         Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
        //         Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;

        //         // get world position of hit
        //         G4ThreeVector globalPosition  = aHit->GetWorldPosition();
        //         Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
        //         Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
        //         Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;

        //         // get local Momentum of hit
        //         G4ThreeVector localMomentum = aHit->GetLocalMomentum();
        //         Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
        //         Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
        //         Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;

        //         // get world Momentum of hit
        //         G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
        //         Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
        //         Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
        //         Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;
        //         Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;
        //         Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;

        //         G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
        //         Float_t rOriginVertexPositionX      = (Float_t) originVertexPosition.x() / cm;
        //         Float_t rOriginVertexPositionY      = (Float_t) originVertexPosition.y() / cm;
        //         Float_t rOriginVertexPositionZ      = (Float_t) originVertexPosition.z() / cm;


        //         G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
        //         Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
        //         Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
        //         Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
        //         Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
        //         Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;

        //         Float_t rOriginVertexKineticEnergy = (Float_t) aHit->GetOriginVertexKineticEnergy() / MeV;
        //         Float_t rOriginVertexTotalEnergy   = (Float_t) aHit->GetOriginVertexTotalEnergy() / MeV;

        //         Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;
        //         TString rParticleName = TString(aHit->GetParticleName());
        //         Int_t rParticleType = (Int_t) aHit->GetParticleType();

        //         // get total Energy of hit
        //         Float_t rTotalEnergy     = (Float_t) aHit->GetTotalEnergy() / MeV;

        //         // get kinetic Energy of hit
        //         Float_t rKineticEnergy     = (Float_t) aHit->GetKineticEnergy() / MeV;

        //         //-----------------------------------
        //         int iHDCID = aHit->GetHDCID();
        //         int iHDC_Chamber = -1; //0 Corresponds to Front, 1 Corresponds to Back
        //         QweakSimUserHDC_SingleHDCEvent* single_hdc_event = 0;
        //         if (iHDCID == 0 || iHDCID == 2){
        //           iHDC_Chamber = 0;
        //           single_hdc_event = &(analysis->fRootEvent->Region2.ChamberFront);
        //         }
        //         if (iHDCID == 1 || iHDCID == 3){
        //           iHDC_Chamber = 1;
        //           single_hdc_event = &(analysis->fRootEvent->Region2.ChamberBack);
        //         }

        //         //-----------------------------------
        //         if (single_hdc_event == 0) {
        //             G4cerr << "HDC hit with incorrect chamber ID." << G4endl;
        //             break;
        //         }
        //         int iHDCpackage = -1; // 0 correponds to package 1, 1 corresponds to package 2
        //                    if(iHDCID == 0 || iHDCID == 1)
        //                    	iHDCpackage = 1;
        //                    if(iHDCID == 2 || iHDCID == 3)
        //                    	iHDCpackage = 2;
        //         //-----------------------------------
        //         int iWirePlaneID = aHit->GetWirePlaneID();
        //         QweakSimUserHDC_WirePlaneEvent* wire_plane_event = 0;
        //         if (iWirePlaneID == 0)
        //           wire_plane_event = &(single_hdc_event->WirePlane1);
        //         if (iWirePlaneID == 1)
        //           wire_plane_event = &(single_hdc_event->WirePlane2);
        //         if (iWirePlaneID == 2)
        //           wire_plane_event = &(single_hdc_event->WirePlane3);
        //         if (iWirePlaneID == 3)
        //           wire_plane_event = &(single_hdc_event->WirePlane4);
        //         if (iWirePlaneID == 4)
        //           wire_plane_event = &(single_hdc_event->WirePlane5);
        //         if (iWirePlaneID == 5)
        //           wire_plane_event = &(single_hdc_event->WirePlane6);

        //         //-----------------------------------
        //         if (wire_plane_event == 0) {
        //             G4cerr << "HDC hit with incorrect plane ID." << G4endl;
        //             break;
        //         }

        //         //-----------------------------------
        //         // store number of hits
        //         HDC_Chamber_Plane_NbOfHits[iHDC_Chamber][iWirePlaneID]++;
        //         wire_plane_event->StoreNbOfHits(HDC_Chamber_Plane_NbOfHits[iHDC_Chamber][iWirePlaneID]);

        //         // mark wire plane as been hit
        //         wire_plane_event->StorePlaneHasBeenHit(5);

        //         // store package hit occurred in
        //         wire_plane_event->StorePackageID(iHDCpackage);

        //         // store origin vertex info
        //         wire_plane_event->StoreOriginVertexPositionX(rOriginVertexPositionX);
        //         wire_plane_event->StoreOriginVertexPositionY(rOriginVertexPositionY);
        //         wire_plane_event->StoreOriginVertexPositionZ(rOriginVertexPositionZ);
        //         wire_plane_event->StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
        //         wire_plane_event->StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
        //         wire_plane_event->StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
        //         wire_plane_event->StoreOriginVertexPhiAngle( rOriginVertexPhiAngle );
        //         wire_plane_event->StoreOriginVertexThetaAngle( rOriginVertexThetaAngle );

        //         wire_plane_event->StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
        //         wire_plane_event->StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);

        //         wire_plane_event->StoreGlobalTimeOfHit(rGlobalTime);
        //         wire_plane_event->StoreParticleName(rParticleName);
        //         wire_plane_event->StoreParticleType(rParticleType);

        //         wire_plane_event->StorePlaneLocalPositionX(rLocalPositionX);
        //         wire_plane_event->StorePlaneLocalPositionY(rLocalPositionY);
        //         wire_plane_event->StorePlaneLocalPositionZ(rLocalPositionZ);

        //         wire_plane_event->StorePlaneGlobalPositionX(rGlobalPositionX);
        //         wire_plane_event->StorePlaneGlobalPositionY(rGlobalPositionY);
        //         wire_plane_event->StorePlaneGlobalPositionZ(rGlobalPositionZ);

        //         // store wire plane hit momentum
        //         wire_plane_event->StorePlaneLocalMomentumX(rLocalMomentumX);
        //         wire_plane_event->StorePlaneLocalMomentumY(rLocalMomentumY);
        //         wire_plane_event->StorePlaneLocalMomentumZ(rLocalMomentumZ);

        //         wire_plane_event->StorePlaneGlobalMomentumX(rGlobalMomentumX);
        //         wire_plane_event->StorePlaneGlobalMomentumY(rGlobalMomentumY);
        //         wire_plane_event->StorePlaneGlobalMomentumZ(rGlobalMomentumZ);

        //         // store global track angles Phi and Theta
        //         wire_plane_event->StoreGlobalPhiAngle(rGlobalPhiAngle);
        //         wire_plane_event->StoreGlobalThetaAngle(rGlobalThetaAngle);

        //         // store total+kinetic energy of hit
        //         wire_plane_event->StoreTotalEnergy(rTotalEnergy);
        //         wire_plane_event->StoreKineticEnergy(rKineticEnergy);

        //         //-----------------------------------

        //     } // end of  for(int i1=0;i1<n_HDChitWirePlane;i1++){

        // } // end of if (n_HDChitWirePlane > 0)


        //===============================================================================================================

        //==============================
        // Store GEM hits into /Region1
        //==============================

/* start GEM comments

        if (n_GEMhitWirePlane > 0) {
            //========================================
            // Store GEM Hit Information into /Region1
            //========================================

            int GEM_ChamberFront_WirePlane_NbOfHits = 0;
            int GEM_ChamberBack_WirePlane_NbOfHits = 0;

            // loop over wire plane hits
            // up to now there should be only one GEM per octant
            for (int i1=0;i1<n_GEMhitWirePlane;i1++) {

                // get hit pointer for each hit
                QweakSimGEM_WirePlaneHit* aHit = (*GEM_WirePlane_HC)[i1];

                // G4cout << G4endl << "###### Printing GEM hit info within QweakSimEventAction::EndOfEventAction() " << G4endl << G4endl;
                if (print_GEM_WirePlaneHit) aHit->Print();

                // get local position of hit
                G4ThreeVector localPosition  = aHit->GetLocalPosition();
                Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
                Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
                Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;

                // get world position of hit
                G4ThreeVector globalPosition  = aHit->GetWorldPosition();
                Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
                Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
                Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;

                // get local Momentum of hit
                G4ThreeVector localMomentum = aHit->GetLocalMomentum();
                Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
                Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
                Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;

                // get world Momentum of hit
                G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
                Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
                Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
                Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;
                Float_t rGlobalThetaAngle = globalMomentum.theta() / degree;
                Float_t rGlobalPhiAngle   = globalMomentum.phi() / degree - 90.0;

                G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
                Float_t rOriginVertexPositionX      = (Float_t) originVertexPosition.x() / cm;
                Float_t rOriginVertexPositionY      = (Float_t) originVertexPosition.y() / cm;
                Float_t rOriginVertexPositionZ      = (Float_t) originVertexPosition.z() / cm;

                G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
                Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
                Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
                Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
                Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
                Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;

                Float_t rOriginVertexKineticEnergy = (Float_t) aHit->GetOriginVertexKineticEnergy() / MeV;
                Float_t rOriginVertexTotalEnergy = (Float_t) aHit->GetOriginVertexTotalEnergy() / MeV;

                Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;

                // get total Energy of hit
                Float_t rtotalEnergy     = (Float_t) aHit->GetTotalEnergy() / MeV;

                // get kinetic Energy of hit
                Float_t rkineticEnergy     = (Float_t) aHit->GetKineticEnergy() / MeV;

                //-----------------------------------

                if ((aHit->GetGEMID()==0) && (aHit->GetWirePlaneID()==0)) {

                    GEM_ChamberFront_WirePlane_NbOfHits++;
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreNbOfHits(GEM_ChamberFront_WirePlane_NbOfHits++);

                    // mark wire plane as been hit
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneHasBeenHit(5);


                    // store origin vertex info
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexPositionX(rOriginVertexPositionX);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexPositionY(rOriginVertexPositionY);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexPhiAngle( rOriginVertexPhiAngle );
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexThetaAngle( rOriginVertexThetaAngle );
                    //------------------------------------------------------------------------------------------------------------------------------------------

                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);

                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreGlobalTimeOfHit(rGlobalTime);


                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneLocalPositionX(rLocalPositionX);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneLocalPositionY(rLocalPositionY);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneLocalPositionZ(rLocalPositionZ);

                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneGlobalPositionX(rGlobalPositionX);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneGlobalPositionY(rGlobalPositionY);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneGlobalPositionZ(rGlobalPositionZ);

                    // store wire plane hit momentum
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneLocalMomentumX(rLocalMomentumX);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneLocalMomentumY(rLocalMomentumY);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneLocalMomentumZ(rLocalMomentumZ);

                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneGlobalMomentumX(rGlobalMomentumX);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneGlobalMomentumY(rGlobalMomentumY);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StorePlaneGlobalMomentumZ(rGlobalMomentumZ);

                    // store global track angles Phi and Theta
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreGlobalPhiAngle(rGlobalPhiAngle);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreGlobalThetaAngle(rGlobalThetaAngle);

                    // store total+kinetic energy of hit
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreTotalEnergy(rtotalEnergy);
                    analysis->fRootEvent->Region1.ChamberFront.WirePlane.StoreKineticEnergy(rkineticEnergy);

                } //end of if((aHit->GetGEMID()==0) && (aHit->GetWirePlaneID()==0))


                //-----------------------------------

                if ((aHit->GetGEMID()==1) && (aHit->GetWirePlaneID()==0)) {

                    GEM_ChamberBack_WirePlane_NbOfHits++;
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreNbOfHits(GEM_ChamberBack_WirePlane_NbOfHits++);

                    // mark wire plane as been hit
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneHasBeenHit(5);

                    // store origin vertex info
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexPositionX(rOriginVertexPositionX);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexPositionY(rOriginVertexPositionY);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexPhiAngle( rOriginVertexPhiAngle );
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexThetaAngle( rOriginVertexThetaAngle );

                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);

                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreGlobalTimeOfHit(rGlobalTime);



                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneLocalPositionX(rLocalPositionX);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneLocalPositionY(rLocalPositionY);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneLocalPositionZ(rLocalPositionZ);

                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneGlobalPositionX(rGlobalPositionX);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneGlobalPositionY(rGlobalPositionY);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneGlobalPositionZ(rGlobalPositionZ);

                    // store wire plane hit momentum
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneLocalMomentumX(rLocalMomentumX);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneLocalMomentumY(rLocalMomentumY);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneLocalMomentumZ(rLocalMomentumZ);

                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneGlobalMomentumX(rGlobalMomentumX);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneGlobalMomentumY(rGlobalMomentumY);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StorePlaneGlobalMomentumZ(rGlobalMomentumZ);

                    // store global track angles Phi and Theta
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreGlobalPhiAngle(rGlobalPhiAngle);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreGlobalThetaAngle(rGlobalThetaAngle);

                    // store total+kinetic energy of hit
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreTotalEnergy(rtotalEnergy);
                    analysis->fRootEvent->Region1.ChamberBack.WirePlane.StoreKineticEnergy(rkineticEnergy);

                } // end of if((aHit->GetGEMID()==1) && (aHit->GetWirePlaneID()==0)) {

                //-----------------------------------


            } // end of  for(int i1=0;i1<n_GEMhitWirePlane;i1++){

        } // end of if ( (n_GEMhitWirePlane == 1)

*/  // end of GEM comments


        // //===========================================================
        // // Store Target hits into /Target
        // //===========================================================

        // if (n_hitTarget >0) {
        //     // initialize deposited energy
        //     Float_t rTotalDepositedEnergy = 0.0;

        //     //--- loop over hits
        //     for (int i1=0;i1<n_hitTarget;i1++) {

        //         QweakSimTarget_DetectorHit* aHit = (*TargetDetector_HC)[i1];

        //         //--- track ID
        //         Float_t rTrackID = (Float_t) aHit->GetTrackID();

        //         //--- particle name & type
        //         TString rParticleName = TString(aHit->GetParticleName());
        //         Int_t rParticleType = (Int_t) aHit->GetParticleType();

        //         //--- get global time of hit
        //         Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;

        //         //--- get world position of hit
        //         G4ThreeVector globalPosition  = aHit->GetWorldPosition();
        //         Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
        //         Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
        //         Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;

        //         //--- get local position of hit
        //         G4ThreeVector localPosition  = aHit->GetLocalPosition();
        //         Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
        //         Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
        //         Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;

        //         //--- get origin vertex position
        //         G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
        //         Float_t rOriginVertexPositionX = (Float_t) originVertexPosition.x() / cm;
        //         Float_t rOriginVertexPositionY = (Float_t) originVertexPosition.y() / cm;
        //         Float_t rOriginVertexPositionZ = (Float_t) originVertexPosition.z() / cm;

        //         //--- get world momentum of hit
        //         G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
        //         //Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
        //         //Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
        //         //Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;

        //         //--- get global theta & phi angle
        //         Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;
        //         Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;

        //         //--- get local momentum of hit
        //         //G4ThreeVector localMomentum = aHit->GetLocalMomentum();
        //         //Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
        //         //Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
        //         //Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;

        //         //--- get local vertex momentum direction of hit
        //         G4ThreeVector localVertexMomentumDirection = aHit->GetMomentumDirection();
        //         Float_t rLocalVertexMomentumDirectionX = (Float_t) localVertexMomentumDirection.x();
        //         Float_t rLocalVertexMomentumDirectionY = (Float_t) localVertexMomentumDirection.y();
        //         Float_t rLocalVertexMomentumDirectionZ = (Float_t) localVertexMomentumDirection.z();

        //         //--- get origin vertex momentum direction of hit
        //         G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
        //         Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
        //         Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
        //         Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();

        //         //--- get origin vertex theta & phi angle
        //         Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;
        //         Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;

        //         //--- get origin vertex kinetic energy & total energy
        //         Float_t rOriginVertexKineticEnergy = (Float_t ) aHit->GetOriginVertexKineticEnergy() / MeV;
        //         Float_t rOriginVertexTotalEnergy   = (Float_t ) aHit->GetOriginVertexTotalEnergy() / MeV;

        //         //--- get total energy & total energy of hit
        //         Float_t rKineticEnergy = (Float_t) aHit->GetKineticEnergy() / MeV;
        //         Float_t rTotalEnergy = (Float_t) aHit->GetTotalEnergy() / MeV;

        //         //--- get PMTOnly deposited energy
        //         Float_t rDepositedEnergy = (Float_t) aHit->GetHitDepositedEnergy() / MeV;
        //         rTotalDepositedEnergy += rDepositedEnergy;

        //         // beam energy (before energy loss in target), explicitly in MeV
        //         G4double beamEnergy = myUserInfo->GetBeamEnergy() / MeV;

        //         // incoming and outgoing momentum (taking into account rastering)
        //         G4ThreeVector p_in = myUserInfo->GetNormMomentum();
        //         G4ThreeVector p_out = globalMomentum;
        //         // scattering angle between incoming and outgoing momentum direction
        //         G4double scatteringAngle = acos(p_out.dot(p_in)/p_in.mag()/p_out.mag());

        //         //--- determine elastic cross section
        //         G4double rElasticCrossSection = 0;
        //         G4double rElasticMomentumTransfer = 0;
        //         G4double rElasticScatteredEnergy = 0;
        //         myUserInfo->GetEPEvent()->Elastic_Cross_Section_Proton(
        //             /* input  */ beamEnergy,
        //                          scatteringAngle,
        //             /* output */ rElasticCrossSection,
        //                          rElasticMomentumTransfer,
        //                          rElasticScatteredEnergy);

        //         //--- store Primary Event Number
        //         analysis->fRootEvent->Target.Detector.StorePrimaryEventNumber((Int_t) PrimaryEventNumber);

        //         //--- store track ID
        //         analysis->fRootEvent->Target.Detector.StoreTrackID(rTrackID);

        //         //--- store particle name & type
        //         analysis->fRootEvent->Target.Detector.StoreParticleName(rParticleName);
        //         analysis->fRootEvent->Target.Detector.StoreParticleType(rParticleType);

        //         //--- store global time of hit
        //         analysis->fRootEvent->Target.Detector.StoreGlobalTimeOfHit(rGlobalTime);

        //         //--- mark target detector as been hit
        //         analysis->fRootEvent->Target.Detector.StoreDetectorHasBeenHit(5);

        //         //--- store global position
        //         analysis->fRootEvent->Target.Detector.StoreDetectorGlobalPositionX(rGlobalPositionX);
        //         analysis->fRootEvent->Target.Detector.StoreDetectorGlobalPositionY(rGlobalPositionY);
        //         analysis->fRootEvent->Target.Detector.StoreDetectorGlobalPositionZ(rGlobalPositionZ);

        //         //--- store local position
        //         analysis->fRootEvent->Target.Detector.StoreDetectorLocalPositionX(rLocalPositionX);
        //         analysis->fRootEvent->Target.Detector.StoreDetectorLocalPositionY(rLocalPositionY);
        //         analysis->fRootEvent->Target.Detector.StoreDetectorLocalPositionZ(rLocalPositionZ);

        //         //--- store origin vertex position
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexPositionX(rOriginVertexPositionX);
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexPositionY(rOriginVertexPositionY);
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexPositionZ(rOriginVertexPositionZ);

        //         //--- store local vertex momentum direction
        //         analysis->fRootEvent->Target.Detector.StoreLocalVertexMomentumDirectionX(rLocalVertexMomentumDirectionX);
        //         analysis->fRootEvent->Target.Detector.StoreLocalVertexMomentumDirectionY(rLocalVertexMomentumDirectionY);
        //         analysis->fRootEvent->Target.Detector.StoreLocalVertexMomentumDirectionZ(rLocalVertexMomentumDirectionZ);

        //         //--- store origin vertex momentum direction
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);

        //         //--- store origin theta & phi angle
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexThetaAngle(rOriginVertexThetaAngle);
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexPhiAngle(rOriginVertexPhiAngle);

        //         //--- store origin kinetic energy & total energy
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
        //         analysis->fRootEvent->Target.Detector.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);

        //         //--- store local vertex kinetic & total energy
        //         analysis->fRootEvent->Target.Detector.StoreDetectorLocalVertexKineticEnergy(rKineticEnergy);
        //         analysis->fRootEvent->Target.Detector.StoreDetectorLocalVertexTotalEnergy(rTotalEnergy);

        //         //--- store global track theta & phi angle
        //         analysis->fRootEvent->Target.Detector.StoreGlobalThetaAngle(rGlobalThetaAngle);
        //         analysis->fRootEvent->Target.Detector.StoreGlobalPhiAngle(rGlobalPhiAngle);

        //         //--- store target deposited energy
        //         analysis->fRootEvent->Target.Detector.StoreDepositedEnergy(rDepositedEnergy);

        //         //--- store elastic cross section
        //         analysis->fRootEvent->Target.Detector.StoreElasticCrossSection(rElasticCrossSection);
        //         analysis->fRootEvent->Target.Detector.StoreElasticScatteredEnergy(rElasticScatteredEnergy);
        //         analysis->fRootEvent->Target.Detector.StoreElasticMomentumTransfer(rElasticMomentumTransfer);

        //     } // end  for(int i1=0;i1<n_hitTarget;i1++)

        //     // Put the total energy into the thing.
        //     analysis->fRootEvent->Target.Detector.StoreTotalEnergyDeposit(rTotalDepositedEnergy);

        // } // end    if (n_hitTarget >0)

        //===========================================================


        // //===========================================================
        // // Store Trigger Scintillator hits into /TriggerScintillator
        // //===========================================================

        // if (n_hitTriggerScintillator >0) {
        //     // initialize deposited energy
        //     Float_t rTotalDepositedEnergy = 0.0;

        //     // loop over hits
        //     for (int i1=0;i1<n_hitTriggerScintillator;i1++) {

        //         QweakSimTriggerScintillator_DetectorHit* aHit = (*TriggerScintillatorDetector_HC)[i1];

        //         if (print_TriggerScintillator_DetectorHit) aHit->Print();

        //         // get local position of hit
        //         G4ThreeVector localPosition  = aHit->GetLocalPosition();
        //         Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
        //         Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
        //         Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;

        //         // get world position of hit
        //         G4ThreeVector globalPosition  = aHit->GetWorldPosition();
        //         Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
        //         Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
        //         Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;

        //         // get local Momentum of hit
        //         //G4ThreeVector localMomentum = aHit->GetLocalMomentum();
        //         //Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
        //         //Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
        //         //Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;

        //         // get world Momentum of hit
        //         G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
        //         //Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
        //         //Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
        //         //Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;

        //         // get global theta & phi angle
        //         Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;
        //         Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;

        //         G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
        //         Float_t rOriginVertexPositionX      = (Float_t) originVertexPosition.x() / cm;
        //         Float_t rOriginVertexPositionY      = (Float_t) originVertexPosition.y() / cm;
        //         Float_t rOriginVertexPositionZ      = (Float_t) originVertexPosition.z() / cm;

        //         G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
        //         Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
        //         Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
        //         Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
        //         Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
        //         Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;

        //         Float_t rOriginVertexKineticEnergy =  (Float_t) aHit->GetOriginVertexKineticEnergy() / MeV;
        //         Float_t rOriginVertexTotalEnergy = (Float_t) aHit->GetOriginVertexTotalEnergy() / MeV;

        //         Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;
        //         TString rParticleName = TString(aHit->GetParticleName());
        //         Int_t rParticleType = (Int_t) aHit->GetParticleType();
                
        //         //--- get TriggerScintillator deposited energy
        //         Float_t rDepositedEnergy = (Float_t) aHit->GetHitDepositedEnergy() / MeV;
        //         rTotalDepositedEnergy += rDepositedEnergy;

        //         // mark TriggerScintillator detector as been hit
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorHasBeenHit(5);

        //         // store global time of hit
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreGlobalTimeOfHit(rGlobalTime);

        //         //--- store particle name & type
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreParticleName(rParticleName);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreParticleType(rParticleType);

        //         // store origin vertex info
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexPositionX(rOriginVertexPositionX);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexPositionY(rOriginVertexPositionY);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexPhiAngle( rOriginVertexPhiAngle );
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexThetaAngle( rOriginVertexThetaAngle );
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);

        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorLocalPositionX(rLocalPositionX);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorLocalPositionY(rLocalPositionY);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorLocalPositionZ(rLocalPositionZ);

        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorGlobalPositionX(rGlobalPositionX);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorGlobalPositionY(rGlobalPositionY);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorGlobalPositionZ(rGlobalPositionZ);

        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDetectorLocalVertexTotalEnergy((Float_t ) aHit->GetTotalEnergy()/MeV);

        //         // store global track angles Phi and Theta
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreGlobalPhiAngle(rGlobalPhiAngle);
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreGlobalThetaAngle(rGlobalThetaAngle);
                
        //         // Store Energy deposited per hit
        //         analysis->fRootEvent->TriggerScintillator.Detector.StoreDepositedEnergy(rDepositedEnergy);

        //     } // end  for(int i1=0;i1<n_hitTriggerScintillator;i1++)
            
        //     // Store total deposited energy
        //     analysis->fRootEvent->TriggerScintillator.Detector.StoreTotalEnergyDeposit(rTotalDepositedEnergy);

        // } // end    if (n_hitTriggerScintillator >0)


////////////////////////////////////////////
   //      //===========================================================
		
//         //===========================================================
//         // Store Lumi hits into /Lumi
//         //===========================================================
		
//         if (n_hitLumi > 0) {
			
//             //--- loop over hits
//             for (int i1=0;i1<n_hitLumi;i1++) {
				
//                 QweakSimLumi_DetectorHit* aHit = (*LumiDetector_HC)[i1];

//                 //--- track ID
//                 Float_t rTrackID = (Float_t) aHit->GetTrackID();
				
//                 //--- particle name & type
//                 TString rParticleName = TString(aHit->GetParticleName());
//                 Int_t rParticleType = (Int_t) aHit->GetParticleType();
				
//                 //--- get global time of hit
//                 Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;
				
//                 //--- GetHasBeenHit();
//                 //GetEdgeEventFlag();
//                 //n_hitLumi <---> GetNbOfHits();
				
//                 //--- get world position of hit
//                 G4ThreeVector globalPosition  = aHit->GetWorldPosition();
//                 Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
//                 Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
//                 Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;
				
//                 //--- get local position of hit
//                 G4ThreeVector localPosition  = aHit->GetLocalPosition();
//                 Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
//                 Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
//                 Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;
				
//                 //--- get local exit position of hit
//                 //G4ThreeVector localExitPosition = aHit->GetLocalExitPosition();
//                 //Float_t rLocalExitPositionX = (Float_t) localExitPosition.x() / cm;
//                 //Float_t rLocalExitPositionY = (Float_t) localExitPosition.y() / cm;
//                 //Float_t rLocalExitPositionZ = (Float_t) localExitPosition.z() / cm;
				
//                 //--- get origin vertex position
//                 G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
//                 Float_t rOriginVertexPositionX = (Float_t) originVertexPosition.x() / cm;
//                 Float_t rOriginVertexPositionY = (Float_t) originVertexPosition.y() / cm;
//                 Float_t rOriginVertexPositionZ = (Float_t) originVertexPosition.z() / cm;
				
//                 //--- get world momentum of hit
//                 G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
//                 //Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
//                 //Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
//                 //Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;

//                 //--- get global theta & phi angle
//                 Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;
//                 Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;
				
//                 //--- get local momentum of hit  
//                 //G4ThreeVector localMomentum = aHit->GetLocalMomentum();
//                 //Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
//                 //Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
//                 //Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;
				
//                 //--- get local vertex momentum direction of hit
//                 G4ThreeVector localVertexMomentumDirection = aHit->GetMomentumDirection();
//                 Float_t rLocalVertexMomentumDirectionX = (Float_t) localVertexMomentumDirection.x();
//                 Float_t rLocalVertexMomentumDirectionY = (Float_t) localVertexMomentumDirection.y();
//                 Float_t rLocalVertexMomentumDirectionZ = (Float_t) localVertexMomentumDirection.z();
				
//                 //--- get origin vertex momentum direction of hit
//                 G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
//                 Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
//                 Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
//                 Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
				
//                 //--- get origin vertex theta & phi angle
//                 Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;
//                 Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
				
//                 //--- get origin vertex kinetic energy & total energy
//                 Float_t rOriginVertexKineticEnergy = (Float_t ) aHit->GetOriginVertexKineticEnergy() / MeV;
//                 Float_t rOriginVertexTotalEnergy   = (Float_t ) aHit->GetOriginVertexTotalEnergy() / MeV;
				
//                 //--- get total energy & total energy of hit
//                 Float_t rKineticEnergy = (Float_t) aHit->GetKineticEnergy() / MeV;
//                 Float_t rTotalEnergy = (Float_t) aHit->GetTotalEnergy() / MeV;
				
//                 //--- get Lumi deposited energy
//                 Float_t rDepositedEnergy = (Float_t) aHit->GetHitDepositedEnergy() / MeV;
//                 //--------------------------------------------------------------------------------------------

//                 //--- store Primary Event Number
//                 analysis->fRootEvent->Lumi.Detector.StorePrimaryEventNumber((Int_t) PrimaryEventNumber);
				
//                 //--- store track ID
//                 analysis->fRootEvent->Lumi.Detector.StoreTrackID(rTrackID);
				
//                 //--- store particle name & type
//                 analysis->fRootEvent->Lumi.Detector.StoreParticleName(rParticleName);
//                 analysis->fRootEvent->Lumi.Detector.StoreParticleType(rParticleType);
				
//                 //--- store global time of hit
//                 analysis->fRootEvent->Lumi.Detector.StoreGlobalTimeOfHit(rGlobalTime);
				
//                 //--- mark Lumi detector as been hit
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorHasBeenHit(5);
//                 //--- edge event flag
//                 //--- Store Nb of hits  --- Done in previous code
//                 //analysis->fRootEvent->Lumi.Detector.StoreDetectorNbOfHits(n_hitLumi);
				
//                 //--- store global position
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorGlobalPositionX(rGlobalPositionX);
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorGlobalPositionY(rGlobalPositionY);
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorGlobalPositionZ(rGlobalPositionZ);
				
//                 //--- store local position
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalPositionX(rLocalPositionX);
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalPositionY(rLocalPositionY);
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalPositionZ(rLocalPositionZ);
				
//                 //analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalExitPositionX(rLocalExitPositionX);
//                 //analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalExitPositionY(rLocalExitPositionY);
//                 //analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalExitPositionZ(rLocalExitPositionZ);
				
//                 //--- store origin vertex position
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexPositionX(rOriginVertexPositionX);
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexPositionY(rOriginVertexPositionY);
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
				
//                 //--- store local vertex momentum direction
//                 analysis->fRootEvent->Lumi.Detector.StoreLocalVertexMomentumDirectionX(rLocalVertexMomentumDirectionX);
//                 analysis->fRootEvent->Lumi.Detector.StoreLocalVertexMomentumDirectionY(rLocalVertexMomentumDirectionY);
//                 analysis->fRootEvent->Lumi.Detector.StoreLocalVertexMomentumDirectionZ(rLocalVertexMomentumDirectionZ);
				
//                 //--- store origin vertex momentum direction
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
				
//                 //--- store origin theta & phi angle
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexThetaAngle(rOriginVertexThetaAngle);
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexPhiAngle(rOriginVertexPhiAngle);
				
//                 //--- store origin kinetic energy & total energy
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
//                 analysis->fRootEvent->Lumi.Detector.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);
				
//                 //--- store local vertex kinetic & total energy
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalVertexKineticEnergy(rKineticEnergy);
//                 analysis->fRootEvent->Lumi.Detector.StoreDetectorLocalVertexTotalEnergy(rTotalEnergy);

//                 //--- store global track theta & phi angle
//                 analysis->fRootEvent->Lumi.Detector.StoreGlobalThetaAngle(rGlobalThetaAngle);
//                 analysis->fRootEvent->Lumi.Detector.StoreGlobalPhiAngle(rGlobalPhiAngle);
				
//                 //--- store Lumi deposited energy
//                 analysis->fRootEvent->Lumi.Detector.StoreDepositedEnergy(rDepositedEnergy);
// 		// -- TotalDepositedEnergy is evaluated in UserLumi class

// 		//--------------------------------------------------------------------------------------------
				
//             } // end  for(int i1=0;i1<n_hitLumi;i1++)	   
//         } // end    if (n_hitLumi >0)
// ////////////////////////////////////////////
		
		
        // //===========================================================
		
        // //===========================================================
        // // Store LeadGlass hits into /LeadGlass
        // //===========================================================
		
        // if (n_hitLeadGlass >0) {
			
        //     //--- loop over hits
        //     for (int i1=0;i1<n_hitLeadGlass;i1++) {
				
        //         QweakSimLeadGlass_DetectorHit* aHit = (*LeadGlassDetector_HC)[i1];

        //         //--- track ID
        //         Float_t rTrackID = (Float_t) aHit->GetTrackID();
				
        //         //--- particle name & type
        //         TString rParticleName = TString(aHit->GetParticleName());
        //         Int_t rParticleType = (Int_t) aHit->GetParticleType();
				
        //         //--- get global time of hit
        //         Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;
				
        //         //--- GetHasBeenHit();
        //         //GetEdgeEventFlag();
        //         //n_hitLeadGlass <---> GetNbOfHits();
				
        //         //--- get world position of hit
        //         G4ThreeVector globalPosition  = aHit->GetWorldPosition();
        //         Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
        //         Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
        //         Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;
				
        //         //--- get local position of hit
        //         G4ThreeVector localPosition  = aHit->GetLocalPosition();
        //         Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
        //         Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
        //         Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;
				
        //         //--- get local exit position of hit
        //         //G4ThreeVector localExitPosition = aHit->GetLocalExitPosition();
        //         //Float_t rLocalExitPositionX = (Float_t) localExitPosition.x() / cm;
        //         //Float_t rLocalExitPositionY = (Float_t) localExitPosition.y() / cm;
        //         //Float_t rLocalExitPositionZ = (Float_t) localExitPosition.z() / cm;
				
        //         //--- get origin vertex position
        //         G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
        //         Float_t rOriginVertexPositionX = (Float_t) originVertexPosition.x() / cm;
        //         Float_t rOriginVertexPositionY = (Float_t) originVertexPosition.y() / cm;
        //         Float_t rOriginVertexPositionZ = (Float_t) originVertexPosition.z() / cm;
				
        //         //--- get world momentum of hit
        //         G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
        //         //Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
        //         //Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
        //         //Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;
				
        //         //--- get global theta & phi angle
        //         Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;
        //         Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;

        //         //--- get local momentum of hit  
        //         //G4ThreeVector localMomentum = aHit->GetLocalMomentum();
        //         //Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
        //         //Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
        //         //Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;
				
        //         //--- get local vertex momentum direction of hit
        //         G4ThreeVector localVertexMomentumDirection = aHit->GetMomentumDirection();
        //         Float_t rLocalVertexMomentumDirectionX = (Float_t) localVertexMomentumDirection.x();
        //         Float_t rLocalVertexMomentumDirectionY = (Float_t) localVertexMomentumDirection.y();
        //         Float_t rLocalVertexMomentumDirectionZ = (Float_t) localVertexMomentumDirection.z();
				
        //         //--- get origin vertex momentum direction of hit
        //         G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
        //         Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
        //         Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
        //         Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
				
        //         //--- get origin vertex theta & phi angle
        //         Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;
        //         Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
				
        //         //--- get origin vertex kinetic energy & total energy
        //         Float_t rOriginVertexKineticEnergy = (Float_t ) aHit->GetOriginVertexKineticEnergy() / MeV;
        //         Float_t rOriginVertexTotalEnergy   = (Float_t ) aHit->GetOriginVertexTotalEnergy() / MeV;
				
        //         //--- get total energy & total energy of hit
        //         Float_t rKineticEnergy = (Float_t) aHit->GetKineticEnergy() / MeV;
        //         Float_t rTotalEnergy = (Float_t) aHit->GetTotalEnergy() / MeV;

        //         //--- get LeadGlass deposited energy
        //         Float_t rDepositedEnergy = (Float_t) aHit->GetHitDepositedEnergy() / MeV;
        //         //--------------------------------------------------------------------------------------------

        //         //--- store Primary Event Number
        //         analysis->fRootEvent->LeadGlass.Detector.StorePrimaryEventNumber((Int_t) PrimaryEventNumber);
				
        //         //--- store track ID
        //         analysis->fRootEvent->LeadGlass.Detector.StoreTrackID(rTrackID);
				
        //         //--- store particle name & type
        //         analysis->fRootEvent->LeadGlass.Detector.StoreParticleName(rParticleName);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreParticleType(rParticleType);
				
        //         //--- store global time of hit
        //         analysis->fRootEvent->LeadGlass.Detector.StoreGlobalTimeOfHit(rGlobalTime);
				
        //         //--- mark LeadGlass detector as been hit
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorHasBeenHit(5);
        //         //--- edge event flag
        //         //--- Store Nb of hits  --- Done in previous code
        //         //analysis->fRootEvent->LeadGlass.Detector.StoreDetectorNbOfHits(n_hitLeadGlass);
				
        //         //--- store global position
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorGlobalPositionX(rGlobalPositionX);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorGlobalPositionY(rGlobalPositionY);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorGlobalPositionZ(rGlobalPositionZ);
				
        //         //--- store local position
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalPositionX(rLocalPositionX);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalPositionY(rLocalPositionY);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalPositionZ(rLocalPositionZ);
				
        //         //analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalExitPositionX(rLocalExitPositionX);
        //         //analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalExitPositionY(rLocalExitPositionY);
        //         //analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalExitPositionZ(rLocalExitPositionZ);
				
        //         //--- store origin vertex position
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexPositionX(rOriginVertexPositionX);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexPositionY(rOriginVertexPositionY);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
				
        //         //--- store local vertex momentum direction
        //         analysis->fRootEvent->LeadGlass.Detector.StoreLocalVertexMomentumDirectionX(rLocalVertexMomentumDirectionX);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreLocalVertexMomentumDirectionY(rLocalVertexMomentumDirectionY);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreLocalVertexMomentumDirectionZ(rLocalVertexMomentumDirectionZ);
				
        //         //--- store origin vertex momentum direction
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
				
        //         //--- store origin theta & phi angle
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexThetaAngle(rOriginVertexThetaAngle);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexPhiAngle(rOriginVertexPhiAngle);
				
        //         //--- store origin kinetic energy & total energy
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);
				
        //         //--- store local vertex kinetic & total energy
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalVertexKineticEnergy(rKineticEnergy);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDetectorLocalVertexTotalEnergy(rTotalEnergy);
				
        //         //--- store global track theta & phi angle
        //         analysis->fRootEvent->LeadGlass.Detector.StoreGlobalThetaAngle(rGlobalThetaAngle);
        //         analysis->fRootEvent->LeadGlass.Detector.StoreGlobalPhiAngle(rGlobalPhiAngle);
				
        //         //--- sotre LeadGlass deposited energy
        //         analysis->fRootEvent->LeadGlass.Detector.StoreDepositedEnergy(rDepositedEnergy);
	// 	// -- TotalDepositedEnergy is evaluated in UserLeadGlass class

	// 	//--------------------------------------------------------------------------------------------
				
        //     } // end  for(int i1=0;i1<n_hitLeadGlass;i1++)	   
        // } // end    if (n_hitLeadGlass >0)
	
	
	//===========================================================
        // Store PMTOnly hits into /PMTOnly
        //===========================================================
		
        // if (n_hitPMTOnly >0) {
        //     // initialize deposited energy
        //     Float_t rTotalDepositedEnergy = 0.0;

        //     //--- loop over hits
        //     for (int i1=0;i1<n_hitPMTOnly;i1++) {
				
        //         QweakSimPMTOnly_DetectorHit* aHit = (*PMTOnlyDetector_HC)[i1];
				
        //         //--- track ID
        //         Float_t rTrackID = (Float_t) aHit->GetTrackID();
				
        //         //--- particle name & type
        //         TString rParticleName = TString(aHit->GetParticleName());
        //         Int_t rParticleType = (Int_t) aHit->GetParticleType();
				
        //         //--- get global time of hit
        //         Float_t rGlobalTime = (Float_t) aHit->GetGlobalTime() / ns;
				
        //         //--- get world position of hit
        //         G4ThreeVector globalPosition  = aHit->GetWorldPosition();
        //         Float_t rGlobalPositionX = (Float_t) globalPosition.x() / cm;
        //         Float_t rGlobalPositionY = (Float_t) globalPosition.y() / cm;
        //         Float_t rGlobalPositionZ = (Float_t) globalPosition.z() / cm;
				
        //         //--- get local position of hit
        //         G4ThreeVector localPosition  = aHit->GetLocalPosition();
        //         Float_t rLocalPositionX = (Float_t) localPosition.x() / cm;
        //         Float_t rLocalPositionY = (Float_t) localPosition.y() / cm;
        //         Float_t rLocalPositionZ = (Float_t) localPosition.z() / cm;
				
        //         //--- get origin vertex position
        //         G4ThreeVector originVertexPosition  = aHit->GetOriginVertexPosition();
        //         Float_t rOriginVertexPositionX = (Float_t) originVertexPosition.x() / cm;
        //         Float_t rOriginVertexPositionY = (Float_t) originVertexPosition.y() / cm;
        //         Float_t rOriginVertexPositionZ = (Float_t) originVertexPosition.z() / cm;
				
        //         //--- get world momentum of hit
        //         G4ThreeVector globalMomentum  = aHit->GetWorldMomentum();
        //         //Float_t rGlobalMomentumX = (Float_t) globalMomentum.x() / MeV;
        //         //Float_t rGlobalMomentumY = (Float_t) globalMomentum.y() / MeV;
        //         //Float_t rGlobalMomentumZ = (Float_t) globalMomentum.z() / MeV;

        //         //--- get global theta & phi angle
        //         Float_t rGlobalThetaAngle = (Float_t) globalMomentum.theta() / degree;
        //         Float_t rGlobalPhiAngle   = (Float_t) globalMomentum.phi() / degree - 90.0;
				
        //         //--- get local momentum of hit  
        //         //G4ThreeVector localMomentum = aHit->GetLocalMomentum();
        //         //Float_t rLocalMomentumX = (Float_t) localMomentum.x() / MeV;
        //         //Float_t rLocalMomentumY = (Float_t) localMomentum.y() / MeV;
        //         //Float_t rLocalMomentumZ = (Float_t) localMomentum.z() / MeV;
				
        //         //--- get local vertex momentum direction of hit
        //         G4ThreeVector localVertexMomentumDirection = aHit->GetMomentumDirection();
        //         Float_t rLocalVertexMomentumDirectionX = (Float_t) localVertexMomentumDirection.x();
        //         Float_t rLocalVertexMomentumDirectionY = (Float_t) localVertexMomentumDirection.y();
        //         Float_t rLocalVertexMomentumDirectionZ = (Float_t) localVertexMomentumDirection.z();
				
        //         //--- get origin vertex momentum direction of hit
        //         G4ThreeVector originVertexMomentumDirection = aHit->GetOriginVertexMomentumDirection();
        //         Float_t rOriginVertexMomentumDirectionX = (Float_t) originVertexMomentumDirection.x();
        //         Float_t rOriginVertexMomentumDirectionY = (Float_t) originVertexMomentumDirection.y();
        //         Float_t rOriginVertexMomentumDirectionZ = (Float_t) originVertexMomentumDirection.z();
				
        //         //--- get origin vertex theta & phi angle
        //         Float_t rOriginVertexThetaAngle = (Float_t) originVertexMomentumDirection.theta() / degree;
        //         Float_t rOriginVertexPhiAngle   = (Float_t) originVertexMomentumDirection.phi() / degree;
				
        //         //--- get origin vertex kinetic energy & total energy
        //         Float_t rOriginVertexKineticEnergy = (Float_t ) aHit->GetOriginVertexKineticEnergy() / MeV;
        //         Float_t rOriginVertexTotalEnergy   = (Float_t ) aHit->GetOriginVertexTotalEnergy() / MeV;
				
        //         //--- get total energy & total energy of hit
        //         Float_t rKineticEnergy = (Float_t) aHit->GetKineticEnergy() / MeV;
        //         Float_t rTotalEnergy = (Float_t) aHit->GetTotalEnergy() / MeV;
				
        //         //--- get PMTOnly deposited energy
        //         Float_t rDepositedEnergy = (Float_t) aHit->GetHitDepositedEnergy() / MeV;
        //         rTotalDepositedEnergy += rDepositedEnergy;

        //         //--- store Primary Event Number
        //         analysis->fRootEvent->PMTOnly.Detector.StorePrimaryEventNumber((Int_t) PrimaryEventNumber);
				
        //         //--- store track ID
        //         analysis->fRootEvent->PMTOnly.Detector.StoreTrackID(rTrackID);
				
        //         //--- store particle name & type
        //         analysis->fRootEvent->PMTOnly.Detector.StoreParticleName(rParticleName);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreParticleType(rParticleType);
				
        //         //--- store global time of hit
        //         analysis->fRootEvent->PMTOnly.Detector.StoreGlobalTimeOfHit(rGlobalTime);
				
        //         //--- mark PMTOnly detector as been hit
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorHasBeenHit(5);
				
        //         //--- store global position
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorGlobalPositionX(rGlobalPositionX);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorGlobalPositionY(rGlobalPositionY);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorGlobalPositionZ(rGlobalPositionZ);
				
        //         //--- store local position
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorLocalPositionX(rLocalPositionX);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorLocalPositionY(rLocalPositionY);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorLocalPositionZ(rLocalPositionZ);
				
        //         //--- store origin vertex position
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexPositionX(rOriginVertexPositionX);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexPositionY(rOriginVertexPositionY);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexPositionZ(rOriginVertexPositionZ);
				
        //         //--- store local vertex momentum direction
        //         analysis->fRootEvent->PMTOnly.Detector.StoreLocalVertexMomentumDirectionX(rLocalVertexMomentumDirectionX);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreLocalVertexMomentumDirectionY(rLocalVertexMomentumDirectionY);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreLocalVertexMomentumDirectionZ(rLocalVertexMomentumDirectionZ);
				
        //         //--- store origin vertex momentum direction
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexMomentumDirectionX(rOriginVertexMomentumDirectionX);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexMomentumDirectionY(rOriginVertexMomentumDirectionY);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexMomentumDirectionZ(rOriginVertexMomentumDirectionZ);
				
        //         //--- store origin theta & phi angle
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexThetaAngle(rOriginVertexThetaAngle);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexPhiAngle(rOriginVertexPhiAngle);
				
        //         //--- store origin kinetic energy & total energy
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexKineticEnergy(rOriginVertexKineticEnergy);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreOriginVertexTotalEnergy(rOriginVertexTotalEnergy);
				
        //         //--- store local vertex kinetic & total energy
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorLocalVertexKineticEnergy(rKineticEnergy);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDetectorLocalVertexTotalEnergy(rTotalEnergy);
				
        //         //--- store global track theta & phi angle
        //         analysis->fRootEvent->PMTOnly.Detector.StoreGlobalThetaAngle(rGlobalThetaAngle);
        //         analysis->fRootEvent->PMTOnly.Detector.StoreGlobalPhiAngle(rGlobalPhiAngle);
				
        //         //--- store PMTOnly deposited energy
        //         analysis->fRootEvent->PMTOnly.Detector.StoreDepositedEnergy(rDepositedEnergy);
				
        //     } // end  for(int i1=0;i1<n_hitPMTOnly;i1++)	   
            
        //     // Put the total energy into the thing. 
        //     analysis->fRootEvent->PMTOnly.Detector.StoreTotalEnergyDeposit(rTotalDepositedEnergy);
            
        // } // end    if (n_hitPMTOnly >0)

        // if (n_hitPMTOnlyPMT >0) {	//--- loop over hits  
            
        //     G4int PMTOnlyPMTHasBeenHit = 0;
        //     G4int PMTOnlyPMTHits = 0;
        //     G4float PMTOnlyNPE = 0;

        //     for (int i1=0;i1<n_hitPMTOnlyPMT;i1++) {
                
        //         QweakSimPMTOnly_PMTHit* aHit = (*PMTOnlyPMT_HC)[i1]; // This line causes a seg fault it seems.
                
        //         PMTOnlyPMTHasBeenHit = 5;
        //         PMTOnlyPMTHits++;
        //         PMTOnlyNPE += myUserInfo->GetNumberOfPhotoelectronsS20(aHit->GetPhotonEnergy()*1.0e6);

        //     } // end	for(n_hitPMTOnlyPMT)
                        
        //     // Write the PMT results to the ROOTfile
        //     analysis->fRootEvent->PMTOnly.PMT.StorePMTHasBeenHit(PMTOnlyPMTHasBeenHit);
        //     analysis->fRootEvent->PMTOnly.PMT.StorePMTTotalNbOfHits(PMTOnlyPMTHits);
        //     analysis->fRootEvent->PMTOnly.PMT.StorePMTTotalNbOfPEs(PMTOnlyNPE);
          
        // } // end	if (n_hitPMTOnlyPMT >0)
        

	//G4cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << G4endl;

        // Finally fill our event ntuple
        analysis->FillRootNtuple();

        //jpan@nuclear.uwinnipeg.ca
        //clear vector contents
        // analysis->fRootEvent->Region3.Clear();
        // analysis->fRootEvent->Region2.Clear();
	//        analysis->fRootEvent->Region1.Clear();
	
    } //end of if( (n_hitWirePlane == 2)&&(n_hitFront >0)&&(n_hitBack >0)&&(n_hitCerenkov >0) )

    myUserInfo->ResetCerenkovSecondaryParticleInfo();

    // print the Eloss for diagnostics
    if(ELOSS_DEBUG)    myUserInfo->PrintELoss();
    // clear the Eloss variables for even PrimaryEventNumber events only, 
    // after they have been stored in the rootfile
    //
    // odd PrimaryEventNumber events are used to generate physics 
    // up to the scattering vertex Z
    //
    // NOTE:: in QweakSimPrimaryGenerator.cc, the PrimaryEventNumber is increased 
    //   at the end of GeneratePrimaries. GeneratePrimaries is called before EnfOfEventAction.
    //   Therefore, even though, even events are used to generate physics events up to the 
    //   scattering vertex Z in GeneratePrimaries, the permutation is flipped inside this
    //   routine.
    myUserInfo->ClearELoss();



    //=======================================================================
    // Save the Ntuple periodically so we have some data in case of a crash

    G4int eventNumber = evt->GetEventID();

    if (eventNumber%1000 == 1)
        analysis->AutoSaveRootNtuple();

    //reset asymInfo	
    asymInfo->at(0)=-2;
    asymInfo->at(2)=-2;
//=======================================================================

} // end of  QweakSimEventAction::EndOfEventAction()

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QweakSimEventAction::Initialize()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double QweakSimEventAction::GetDistance(G4ThreeVector p1,G4ThreeVector p2) {
    return sqrt((p1.x()-p2.x())*(p1.x()-p2.x())+
                (p1.y()-p2.y())*(p1.y()-p2.y())+
                (p1.z()-p2.z())*(p1.z()-p2.z()));
}


////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  Jim Dowd
// ---------------------------------------------------------
//      Calculates and stores the kinematic variables Nu, Q2,
//      x, and W.
// ---------------------------------------------------------
void QweakSimEventAction::CalculateKinematicVariables() 
{
  // Vertex Kinematic variables
  // These variables are only valid for events generated with the standard
  // generator (not the lookup table).
  G4double mp    = 938.2796*MeV;
  G4double E_in  = myUserInfo->GetPreScatteringKineticEnergy() + 0.511*MeV;
  G4double E_out = myUserInfo->GetOriginVertexTotalEnergy();
  G4double theta = myUserInfo->GetOriginVertexThetaAngle();

  G4double nu    = E_in - E_out;
  G4double q2    = 4.0*E_in*E_out*sin(theta*degree/2.0)*sin(theta*degree/2.0);
  G4double x     = q2/(2.0*mp*nu);
  G4double w     = sqrt(mp*mp+2.0*mp*nu-q2);

  myUserInfo->StoreOriginVertexKinematicNu(nu);
  myUserInfo->StoreOriginVertexKinematicQ2(q2*0.000001);
  myUserInfo->StoreOriginVertexKinematicX(x);
  myUserInfo->StoreOriginVertexKinematicW(w);

  //G4cout << "==== Vertex Kinematics ====" << G4endl;
  //G4cout << "E_in:         " << E_in  << G4endl;
  //G4cout << "E_out:        " << E_out << G4endl;
  //G4cout << "Theta:        " << theta << G4endl;
  //G4cout << "Nu:           " << nu << G4endl;
  //G4cout << "Q2:           " << q2 << G4endl;
  //G4cout << "X:            " << x << G4endl;
  //G4cout << "W:            " << w << G4endl;

  
  // "Effective" kinematic variables
  //  These kinematic variables are not calculated correctly yet.
  //  Need values for E_out and theta grabbed from just downstream 
  //  of the target, instead of the vertex quantities used below.
  //  They are still valid for any events generated using the 
  //  lookup table.
  E_in  = myUserInfo->GetBeamEnergy();
  E_out = myUserInfo->GetOriginVertexTotalEnergy(); 
  theta = myUserInfo->GetOriginVertexThetaAngle(); 

  nu    = E_in - E_out;
  q2    = 4.0*E_in*E_out*sin(theta*degree/2.0)*sin(theta*degree/2.0);
  x     = q2/(2.0*mp*nu);
  w     = sqrt(mp*mp+2.0*mp*nu-q2);

  myUserInfo->StoreEffectiveKinematicNu(nu);
  myUserInfo->StoreEffectiveKinematicQ2(q2*0.000001);
  myUserInfo->StoreEffectiveKinematicX(x);
  myUserInfo->StoreEffectiveKinematicW(w);

  //G4cout << "==== Effective Kinematics ====" << G4endl;
  //G4cout << "E_in:         " << E_in  << G4endl;
  //G4cout << "E_out:        " << E_out << G4endl;
  //G4cout << "Theta:        " << theta << G4endl;
  //G4cout << "Nu:           " << nu << G4endl;
  //G4cout << "Q2:           " << q2 << G4endl;
  //G4cout << "X:            " << x << G4endl;
  //G4cout << "W:            " << w << G4endl;
  //G4cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  Jim Dowd
// ---------------------------------------------------------
//      Calculates a rate based on the generation limits and
//      target type.  
// ---------------------------------------------------------
G4double QweakSimEventAction::CalculateRate(G4double xsec, G4int PEs) 
{
  G4double rate = 2.0*myUserInfo->GetPhaseSpace()*myUserInfo->GetLuminosity()*xsec*PEs/myUserInfo->GetNumberOfEventToBeProcessed();
  return rate;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

