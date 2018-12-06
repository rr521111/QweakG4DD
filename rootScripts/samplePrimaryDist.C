#include <iostream>
#include <fstream>

void readDist();
void sampleDist(int nevents, int vPol, int nDist, string Pol, int nOct);
void drawDist();
void getPol(double pos, double &polX, double &polY, double &polZ, string Pol, int nOct);
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


void samplePrimaryDist(int nevents, int vPol, int chooseDist=0, string Pol="ini"){
  
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
    sampleDist(nevents, nDist, vPol, Pol, nOct);
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
      fnm=Form("../input/Hit_Map_Tracks_%d_MD%d_.txt",trackingRuns[nOct][nData],nOct+1);
  }else
    //these are the face of the Pb -- no projection needed
    fnm=Form("../input/MC_HitMap_Oct%d.txt",nOct+1);

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

void sampleDist(int nevents, int nDist, int vPol, string Pol, int nOct){
  ofstream fout(Form("positionMomentum_%d.in",nDist),std::ofstream::out);
  ofstream fpol(Form("polarization_%d.in",nDist),std::ofstream::out);
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
    getPol(y, pX, pY, pZ, Pol, nOct);
    
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

void getPol(double pos, double &polX, double &polY, double &polZ, string Pol, int nOct){
  //coefficients for along bar(X) and across bar(Y), and forwar beam direction (Z) fits as a function of X position.
  double polXV[3][8] = {{0.9889890676,0.4093218090,0.0002040015,-0.4144668259,-0.9890941219,-0.4089285142,-0.0001610524,0.4110589110}, {-0.0000016391,0.0008429048,-0.0015633406,0.0007949866,0.0000025756,-0.0008331212,0.0015667701,-0.0008248167}, {-0.0000060990,0.0000069103,0.0000000621,-0.0000063397,0.0000061075,-0.0000070491,-0.0000000595,0.0000066280}};
  double polYV[3][8] = {{0.0005595607,-0.6890668688,-0.5973680930,-0.6897161791,-0.0003519301,0.6890262365,0.5971136353,0.6893593540}, {0.0013505283,-0.0010833076,-0.0000090140,0.0010776093,-0.0013515278,0.0010818298,0.0000087267,-0.0010818899}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
  double polZV[3][8] = {{0.0003836084,-0.5519215775,-0.7848460373,-0.5495935204,-0.0004660876,0.5521707751,0.7847000508,0.5513049323}, {0.0031107007,0.0021428279,0.0000053672,-0.0021118644,-0.0031091568,-0.0021378100,-0.0000043785,0.0021365252}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
  double polXH[3][8] = {{0.0000220349,-0.4085240817,-0.9889386626,-0.4124544184,0.0002140031,0.4095578669,0.9889345370,0.4103398378}, {-0.0015693989,0.0010592527,-0.0000022981,-0.0011045102,0.0015706064,-0.0010671500,0.0000029354,0.0010843153}, {-0.0000000532,-0.0000068347,0.0000060471,-0.0000063123,-0.0000000350,0.0000069056,-0.0000060575,0.0000065934}};
  double polYH[3][8] = {{-0.5979350144,-0.6872327977,0.0006418237,0.6869351798,0.5974431329,0.6874576007,-0.0008940756,-0.6871486917}, {0.0000080322,0.0010976883,-0.0013427450,0.0010958530,0.0000049411,-0.0010942259,0.0013481661,-0.0011014535}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
  double polZH[3][8] = {{0.7840770649,0.5524548490,-0.0006835132,-0.5509005813,-0.7844660853,-0.5515524500,0.0007839123,0.5517334626}, {0.0000046206,0.0022611748,0.0031021917,0.0022921491,0.0000057415,-0.0022651743,-0.0031067622,-0.0022771484}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};

  //[pos]=cm=position along bar
  if(Pol=="V"||Pol=="mV"){
    polX = polXV[0][nOct] + pos*polXV[1][nOct] + pos*pos*polXV[2][nOct];
    polY = polYV[0][nOct] + pos*polYV[1][nOct] + pos*pos*polYV[2][nOct];
    polZ = polZV[0][nOct] + pos*polZV[1][nOct] + pos*pos*polZV[2][nOct];
  }
  else if(Pol=="H"||Pol=="mH"){
    polX = polXH[0][nOct] + pos*polXH[1][nOct] + pos*pos*polXH[2][nOct];
    polY = polYH[0][nOct] + pos*polYH[1][nOct] + pos*pos*polYH[2][nOct];
    polZ = polZH[0][nOct] + pos*polZH[1][nOct] + pos*pos*polZH[2][nOct];
  }
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

