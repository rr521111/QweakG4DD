#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "QweakSimUserMainEvent.hh"
#include "QweakSimSystemOfUnits.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "G4ThreeVector.hh"

using namespace std;

int evNr;
int primary;//0 secondary, 1 primary
float x,y,z,E,phi,th;
float angX,angY,angXi,angYi;
float xi,yi,zi;
float polT,polTi;
double asymPM,asymPP;
float polXraw,polYraw,polZraw,PolMag,polX,polY,polZ;
float polXiraw,polYiraw,polZiraw,PoliMag,polXi,polYi,polZi;
float momXi,momYi,momZi,momX,momY,momZ,MomMag;
float ScatterXraw,ScatterYraw,ScatterZraw,ScatterMag,ScatX,ScatY,ScatZ;
float SpinXraw,SpinYraw,SpinZraw,SpinMag,SpinX,SpinY,SpinZ;
float theta2,phi2sgn,phi2;
int fixedPos(0);
std::vector<double> xI,yI,zI,aXi,aYi,polXI,polYI,polZI,polI;

void findInt(std::vector<int> &inter,std::vector<int> &val, int trackID,int parent, int &hasPar, int &nInt);
void processOne(TTree *QweakSimG4_Tree,TTree *tout, long &nrEvts);
void readInitial(string fnm);


int main(int argc, char** argv){

  if( argc < 3 ) {
    cout<<" usage: build/QweakAna [path to infile with list of output QweakSimG4 trees] [suffix for outfile] [optional: 1 for fixed position]"<<endl;
    return 1;
  }

  string files(argv[1]);
  string suffix(argv[2]);
  if(argc==4)
    fixedPos=atoi(argv[3]);

  TFile *fout=new TFile(Form("o_hits_%s.root",suffix.c_str()),"RECREATE");
  TH1D *hEntries=new TH1D("hEntries","number of processed events",2,0,2);
  TTree *tout=new TTree("t","Stripped QweakSimG4 tree for hits");
  tout->Branch("evNr",&evNr,"evNr/I");
  tout->Branch("primary",&primary,"primary/I");
  tout->Branch("x",&x,"x/F");
  tout->Branch("y",&y,"y/F");
  tout->Branch("z",&z,"z/F");
  tout->Branch("E",&E,"E/F");
  tout->Branch("phi",&phi,"phi/F");
  tout->Branch("th",&th,"th/F");
  tout->Branch("angX",&angX,"angX/F");
  tout->Branch("angY",&angY,"angY/F");
  tout->Branch("polT",&polT,"polT/F");
  //tout->Branch("polXraw",&polXraw,"polXraw/F");	// comment us!
  //tout->Branch("polYraw",&polYraw,"polYraw/F");	//
  //tout->Branch("polZraw",&polZraw,"polZraw/F");	//
  //tout->Branch("MomMag",&MomMag,"MomMag/F");		//
  tout->Branch("xi",&xi,"xi/F");
  tout->Branch("yi",&yi,"yi/F");
  tout->Branch("zi",&zi,"zi/F");
  tout->Branch("angXi",&angXi,"angXi/F");
  tout->Branch("angYi",&angYi,"angYi/F");
  //tout->Branch("polXiraw",&polXiraw,"polXiraw/F");	//
  //tout->Branch("polYiraw",&polYiraw,"polYiraw/F");	//
  //tout->Branch("polZiraw",&polZiraw,"polZiraw/F");	//
  tout->Branch("polTi",&polTi,"polTi/F");
  tout->Branch("asymPM",&asymPM,"asymPM/D");
  tout->Branch("asymPP",&asymPP,"asymPP/D");
  tout->Branch("theta2",&theta2,"theta2/F");
  tout->Branch("phi2",&phi2,"phi2/F");
  
  long totEv=0;
  long simEvts(0);
  if ( files.find(".root") < files.size() ){
    if(!fixedPos) readInitial(files);
    TFile *fin=new TFile(files.c_str(),"READ");
    if(!fin->IsOpen()){
      cout<<"Problem: can't find file: "<<files<<endl;
      fout->Close();
      return -3;
    }

    TTree *QweakSimG4_Tree=(TTree*)fin->Get("QweakSimG4_Tree");
    cout<<"processing only one file: "<<files.c_str()<<endl;    
    
    long evts=QweakSimG4_Tree->GetEntries();
    totEv+=evts;
    cout<<" total nr ev: "<<evts<<" "<<totEv<<endl;
    
    processOne(QweakSimG4_Tree,tout,simEvts);
    simEvts=std::ceil(simEvts/1000.)*1000;
    
    cout<<"Done looping!"<<endl;
    fin->Close();
    delete fin;
  }else{
    ifstream ifile(files.c_str());
    string data;
    while(ifile>>data){
      if(!fixedPos) readInitial(data);
      TFile *fin=new TFile(data.c_str(),"READ");
      if(!fin->IsOpen()){
	cout<<"Problem: can't find file: "<<data<<endl;
	continue;
      }
      TTree *QweakSimG4_Tree=(TTree*)fin->Get("QweakSimG4_Tree");
      cout<<"processing : "<<data.c_str()<<endl;          
      long evts=QweakSimG4_Tree->GetEntries();
      totEv+=evts;
      cout<<" total nr ev: "<<evts<<" "<<totEv<<endl;      
      processOne(QweakSimG4_Tree,tout,simEvts);      
      simEvts=std::ceil(simEvts/1000.)*1000;
      fin->Close();
      delete fin;
    }
    ifile.close();
  }
  
  cout<<"Processed "<<totEv<<" events"<<endl;
  cout<<"\tfrom "<<simEvts<<" simulated events"<<endl;
  
  hEntries->SetBinContent(1,totEv);
  hEntries->SetBinContent(2,simEvts);
  fout->cd();
  hEntries->Write();
  tout->Write();
  fout->Close();
  return 0;
}

void processOne(TTree *QweakSimG4_Tree, TTree *tout, long &nrEvts){

  std::vector<int> interaction;
  std::vector<int> trackID;

  //set addresses of leafs
  QweakSimUserMainEvent* event = 0;
  QweakSimG4_Tree->SetBranchAddress("QweakSimUserMainEvent",&event);    

  for (int i = 0; i < QweakSimG4_Tree->GetEntries(); i++) {
    QweakSimG4_Tree->GetEntry(i);
    if(i%10000==1) cout<<"   at event: "<<i<<endl;

    int nThrown = event->Primary.GetPrimaryEventNumber()-1;
    evNr = nrEvts + nThrown + 1;
    double asymPpM = event->Primary.GetAsymDeno();
    double asymPmM = event->Primary.GetAsymNomi();
    asymPP = (asymPpM + asymPmM)/2;
    asymPM = (asymPpM - asymPmM)/2;

    if(!fixedPos){
      xi=xI[nThrown];
      yi=yI[nThrown];
      zi=zI[nThrown];
      angXi=aXi[nThrown];
      angYi=aYi[nThrown];
      polXiraw=polXI[nThrown];
      polYiraw=polYI[nThrown];
      polZiraw=polZI[nThrown];
      PoliMag=sqrt(polXI[nThrown]*polXI[nThrown]+polYI[nThrown]*polYI[nThrown]+polZI[nThrown]*polZI[nThrown]);
      polXi=polXiraw/PoliMag;
      polYi=polYiraw/PoliMag;
      polZi=polZiraw/PoliMag;
      momXi=tan(angXi*pi/180.0)/sqrt(1 + tan(angXi*pi/180.0)*tan(angXi*pi/180.0) + tan(angYi*pi/180.0)*tan(angYi*pi/180.0));
      momYi=tan(angYi*pi/180.0)/sqrt(1 + tan(angXi*pi/180.0)*tan(angXi*pi/180.0) + tan(angYi*pi/180.0)*tan(angYi*pi/180.0));
      momZi=1/sqrt(1 + tan(angXi*pi/180.0)*tan(angXi*pi/180.0) + tan(angYi*pi/180.0)*tan(angYi*pi/180.0));
      polTi=polI[nThrown];
    }else{
      xi=0.;
      yi=335.0;
      zi=560.;
      angXi=0.;
      angYi=0.;
      polTi=1.;
    }
    
    interaction.clear();
    trackID.clear();
    
    for (int hit = 0; hit < event->Cerenkov.Detector.GetDetectorNbOfHits(); hit++) {	
      if(event->Cerenkov.Detector.GetDetectorID()[hit]!=3) continue;
      int pTypeHit=event->Cerenkov.Detector.GetParticleType()[hit];
      if(abs(pTypeHit)!=11) continue;
      E=event->Cerenkov.Detector.GetTotalEnergy()[hit];
      
      int parentID=event->Cerenkov.Detector.GetParentID()[hit];
      int tID     =event->Cerenkov.Detector.GetParticleID()[hit];
      
      int hasParent(-1),nInt(-1);
      findInt(interaction,trackID,tID,parentID,hasParent,nInt);
      if(nInt!=1 || hasParent==1) continue;

      primary=0;
      if(tID==1 && parentID==0) primary=1;
      G4ThreeVector pol(event->Cerenkov.Detector.GetPolarizationX()[hit],
			event->Cerenkov.Detector.GetPolarizationY()[hit],
			event->Cerenkov.Detector.GetPolarizationZ()[hit]);

      G4ThreeVector mom(event->Cerenkov.Detector.GetGlobalMomentumX()[hit],
			event->Cerenkov.Detector.GetGlobalMomentumY()[hit],
			event->Cerenkov.Detector.GetGlobalMomentumZ()[hit]);
      //pol.rotateUz(mom.unit());
      polT = pol.unit() * mom.unit();
      MomMag = mom.mag();
      PolMag = pol.mag();

      polXraw = (event->Cerenkov.Detector.GetPolarizationX()[hit]);
      polYraw = (event->Cerenkov.Detector.GetPolarizationY()[hit]);
      polZraw = (event->Cerenkov.Detector.GetPolarizationZ()[hit]);

      polX = polXraw/PolMag;
      polY = polYraw/PolMag;
      polZ = polZraw/PolMag;

      x=event->Cerenkov.Detector.GetDetectorGlobalPositionX()[hit];
      y=event->Cerenkov.Detector.GetDetectorGlobalPositionY()[hit];
      z=event->Cerenkov.Detector.GetDetectorGlobalPositionZ()[hit];
      
      double Gphi   = event->Cerenkov.Detector.GetGlobalPhiAngle()[hit];
      double Gtheta = event->Cerenkov.Detector.GetGlobalThetaAngle()[hit];	
      double _Gtheta=Gtheta/180.*pi;
      double _Gphi=(Gphi+90)/180.*pi; //+90 to account for the offset in the output
      phi = Gphi;
      th  = Gtheta;
      angX = atan2(sin(_Gtheta)*cos(_Gphi),cos(_Gtheta)) * 180.0 / pi;
      angY = atan2(sin(_Gtheta)*sin(_Gphi),cos(_Gtheta)) * 180.0 / pi;

      momX=tan(angX*pi/180.0)/sqrt(1 + tan(angX*pi/180.0)*tan(angX*pi/180.0) + tan(angY*pi/180.0)*tan(angY*pi/180.0));
      momY=tan(angY*pi/180.0)/sqrt(1 + tan(angX*pi/180.0)*tan(angX*pi/180.0) + tan(angY*pi/180.0)*tan(angY*pi/180.0));
      momZ=1/sqrt(1 + tan(angX*pi/180.0)*tan(angX*pi/180.0) + tan(angY*pi/180.0)*tan(angY*pi/180.0));

      theta2=acos(momX*momXi + momY*momYi + momZ*momZi)*180.0/pi;

      ScatterXraw=momYi*momZ - momZi*momY;
      ScatterYraw=-1.0*(momXi*momZ - momZi*momX);
      ScatterZraw=momXi*momY - momYi*momX;

      ScatterMag=sqrt(ScatterXraw*ScatterXraw + ScatterYraw*ScatterYraw + ScatterZraw*ScatterZraw);

      ScatX=ScatterXraw/ScatterMag;
      ScatY=ScatterYraw/ScatterMag;
      ScatZ=ScatterZraw/ScatterMag;

      SpinXraw=momYi*polZi - momZi*polYi;
      SpinYraw=-1.0*(momXi*polZi - momZi*polXi);
      SpinZraw=momXi*polYi - momYi*polXi;

      SpinMag=sqrt(SpinXraw*SpinXraw + SpinYraw*SpinYraw + SpinZraw*SpinZraw);

      SpinX=SpinXraw/SpinMag;
      SpinY=SpinYraw/SpinMag;
      SpinZ=SpinZraw/SpinMag;

      phi2sgn=(momXi*(polYi*ScatZ - polZi*ScatY) - momYi*(polXi*ScatZ - polZi*ScatX) + momZi*(polXi*ScatY - polYi*ScatX))/abs(momXi*(polYi*ScatZ - polZi*ScatY) - momYi*(polXi*ScatZ - polZi*ScatX) + momZi*(polXi*ScatY - polYi*ScatX));

      phi2=phi2sgn*acos((polXi*ScatX + polYi*ScatY + polZi*ScatZ)/SpinMag)*180.0/pi;

      tout->Fill();
    }//nhit
  }//tree entries
  nrEvts = evNr;
}

void readInitial(string fnm){
  xI.clear();yI.clear();zI.clear();aXi.clear();aYi.clear();polXI.clear();polYI.clear();polZI.clear();polI.clear();

  string posfnm=fnm.substr(0,fnm.rfind("/")+1)+"positionMomentum.in";
  string polfnm=fnm.substr(0,fnm.rfind("/")+1)+"polarization.in";
  ifstream finpos(posfnm.c_str());
  ifstream finpol(polfnm.c_str());
  float tmpx,tmpy,tmpz,tmpax,tmpay,tmppx,tmppy,tmppz;
  while(finpos>>tmpx>>tmpy>>tmpz>>tmpax>>tmpay){
    finpol>>tmppx>>tmppy>>tmppz;
    xI.push_back(tmpx);
    yI.push_back(tmpy);
    zI.push_back(tmpz);
    aXi.push_back(tmpax);
    aYi.push_back(tmpay);
    polXI.push_back(tmppx);
    polYI.push_back(tmppy);
    polZI.push_back(tmppz);
    polI.push_back(sqrt(pow(tmppx,2)+pow(tmppy,2)));
  }
  finpos.close();
  finpol.close();
}

void findInt(std::vector<int> &inter,std::vector<int> &val, int trackID,int parent, int &hasPar, int &nInt){
  int found=0;
  nInt=-2;
  int findParent=0;
  for(unsigned int i=0;i<val.size();i++){
    if(trackID==val[i]){
      inter[i]++;
      nInt=inter[i];
      found++;
    }
    if(parent==val[i])
      findParent++;
  }
  
  if(findParent) hasPar=1;
  else hasPar=0;
  
  if(!found){
    val.push_back(trackID);
    inter.push_back(1);
    nInt=1;
  }
  if(found>1){
    cout<<"multiple entries for track "<<trackID<<endl;
  }
}
