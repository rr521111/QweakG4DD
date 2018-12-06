#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//                                                                                                        //
#//      New sampling script for polarization maps.							   //
#//                                                                                                        //
#//      Robert Radloff, Ohio University, 2018                                                             //
#//                                                                                                        //
#////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "/u/home/rradloff/dd/QweakG4DD/src/interpolatePOLsNEW.cc"

void readDist();
void sampleDist(int nevents, int vPol, int nDist, string Pol, int nOct, string outpath);
void drawDist();
double getPbPos(double pos,double ang);
double getAngY(double posY);

TH3D *hIn;

// which distribution to sample:
//     1Yx from data
//     20x from data
//   x = oct number;
//   Y either 0 or 1 for data refers to either the first or second tracking run (see readDist)
int nDist;
double getAngYa, getAngYb;
int draw=0;


void samplePrimaryDist3(int nevents, int vPol, int chooseDist=0, string Pol="ini", string _outpath = ""){
  
  if(chooseDist!=0)
    nDist=chooseDist;
  
  if(gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  if(nDist<200)
    //use for files with Hit_Map_Tracks ... DA tracking files
    hIn=new TH3D("hIn","input h",15,320,350,20,-100,100,14,0.34,0.48);
  else    
    //use this for MC files from qelog Det 117 (JP)
    hIn=new TH3D("hIn","input h",24,323,347,200,-100,100,14,0.34,0.48);
  int nOct=nDist%100 - 1;
  readDist();

  if(draw)
    drawDist();
  else
    sampleDist(nevents, nDist, vPol, Pol, nOct, _outpath.c_str());
}

void readDist(){
  int nOct=nDist%100 - 1;
  int nData= ( (nDist-nDist%100)/10 ) %10;
  string fnm;
  int trackingRuns[8][2]={{15121,13681},{13671,0},{13679,13676},{13674,0},{13681,15121},{13671,0},{13676,13679},{13674,0}};

  //for equation a*X+b =>JP det elog 117
  double a[8]={1.37e-3, 1.37e-3, 1.37e-3, 1.37e-3, 1.38e-3, 1.37e-3, 1.37e-3, 1.37e-3};
  double b[8]={-1.8e-4, -1.9e-4,  1.8e-4,  3.1e-4, -2.0e-5,  1.6e-4,  2.2e-4,  1.1e-4};
  getAngYa=a[nOct];
  getAngYb=b[nOct];
  
  if(nDist<200){
    if(trackingRuns[nOct][nData]==0){
      cout<<"Can't find tracking data run: "<<nOct<<" "<<nData<<" "<<endl;
      gSystem->Exit(1);
    }else
      //these are at the MD => projection required
      fnm=Form("/u/home/rradloff/dd/QweakG4DD/input/Hit_Map_Tracks_%d_MD%d_.txt",trackingRuns[nOct][nData],nOct+1);
  }else
    //these are the face of the Pb -- no projection needed
    fnm=Form("/u/home/rradloff/dd/QweakG4DD/input/MC_HitMap_Oct%d.txt",nOct+1);

  ifstream fin(fnm.c_str());

  double x,y,xs,val;
  while(fin>>x>>y>>xs>>val){
    if(val==0)continue;    
    if(nDist<200)
      hIn->SetBinContent(x+1,y+1,xs+1,val);//if using tracking data
    else
      hIn->SetBinContent(x,y,xs,val);//if using sim data
  }
}

void sampleDist(int nevents, int nDist, int vPol, string Pol, int nOct, string outpath){
  ofstream fout(Form("%s/positionMomentum.in",outpath.c_str()),std::ofstream::out);
  ofstream fpol(Form("%s/polarization.in",outpath.c_str()),std::ofstream::out);
  int nv=0;
  
  do{
    double x(-1);//x position (across bar) 
    double y(-1);//y position (along bar)  
    double z(-1);//x angle (across bar)    
    hIn->GetRandom3(x,y,z);

    double pbZpos=571.9;//cm
    double deg=180./3.14159265359;
    double pbXang=getAngY(y);
    double pbYang=z;
    
    double pbXpos,pbYpos;
    if(nDist<200){
      pbYpos=getPbPos(x,pbYang);//this projects from the face of the quartz
      pbXpos=getPbPos(y,pbXang);//this projects from the face of the quartz
    }else{
      pbXpos=y;
      pbYpos=x;
    }

    if( (pow(sin(pbYang),2)+pow(sin(pbXang),2)) > 1 ) continue;
    if( pbYpos<326 || pbYpos>344 ) continue;
    if( fabs(pbXpos)>=100 ) continue;

    fout<<pbXpos<<" "<<pbYpos<<" "<<pbZpos<<" "<<pbXang*deg<<" "<<pbYang*deg<<endl;

    double pX = 0;
    double pY = 0;
    double pZ = 0;

    if(Pol == "H" || Pol == "mH"){
      interpolatePOLsNEW(y,x,z,pY,pX,pZ,Form("polmap%d_H",nOct+1));
    }
    
    if(Pol == "V" || Pol == "mV"){
      interpolatePOLsNEW(y,x,z,pY,pX,pZ,Form("polmap%d_V",nOct+1));
    }
    
    if(Pol == "L" || Pol == "mL"){
      interpolatePOLsNEW(y,x,z,pY,pX,pZ,Form("polmap%d_L",nOct+1));
    }
    
    if(pY>1) pY=1;
    fpol<<vPol*pX<<" "<<vPol*pY<<" "<<vPol*pZ<<endl;

    nv++;
  }while(nv<=nevents);

  fout.close();
  fpol.close();
}

void drawDist(){
  TCanvas *c1=new TCanvas("c1","c1",1600,1400);
  c1->Divide(2);
  c1->cd(1);
  TH1 *h1=hIn->Project3D("xy");
  //TH2D *h1=hIn->Project3D("xy");
  h1->DrawCopy("colz");
  c1->cd(2);
  TH1 *h2=hIn->Project3D("y");
  //TH1D *h2=hIn->Project3D("y");
  h2->DrawCopy("colz");
  gPad->SetGridx(1);
  gPad->SetGridy(1);
}

double getPbPos(double pos,double ang){
  return 5.375*tan(ang)+pos;//[cm]
}

double getAngY(double posY){//[posY]=cm

  if(nDist<200)
    return (1.375e-3 * posY + 0.01); //from DA - data [rad]
  else
    return getAngYa * posY + getAngYb;//from JP - det elog 117
}

