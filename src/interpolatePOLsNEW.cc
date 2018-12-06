//////////////////////////////////////////////////////////////////////////////////////////////////////////
//													//
//	Attempt at a more efficient and transparent interpolatePOLs macro.				//
//													//
//	Robert Radloff, Ohio University, 2018								//
//													//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;
	
bool debug = false;

int nbinsx = 200;
int nbinsy = 24;
int nbinsz = 14;
	
float xmin = -100.0;
float xmax = 100.0;	

float ymin = 323.0;
float ymax = 347.0;

float zmin = 0.34;
float zmax = 0.48;
	
float xstep = (xmax-xmin)/nbinsx;
float ystep = (ymax-ymin)/nbinsy;
float zstep = (zmax-zmin)/nbinsz;

float polmapin[201][25][15][3];

float x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;

void interpolatePOLsNEW(double xin, double yin, double zin, double &pac, double &pal, double &pz, string polModel = "ini"){
	
	//CHOOSING INPUT MAP
	
	if(debug)
		cout << "Starting interpolation with:" << endl << "x y z    " << xin << " " << yin << " " << zin << endl << endl;
	
	string path;
  	if("polmap1_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD1_H.txt";
  	}else if("polmap2_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD2_H.txt";
  	}else if("polmap3_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD3_H.txt";
  	}else if("polmap4_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD4_H.txt";
  	}else if("polmap5_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD5_H.txt";
  	}else if("polmap6_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD6_H.txt";
  	}else if("polmap7_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD7_H.txt";
  	}else if("polmap8_H" == polModel) {
    		path = "/u/home/rradloff/dd/QweakG4DD/input/polmap_MD8_H.txt";
  	}else{
    		cout << "Invalid choice of polarization map!" << endl << polModel << endl;	
  	}
	
	ifstream fin(path.c_str());
	if(!fin.is_open()) {
		cout << "Cannot read polmap file!" << endl;
	}
	
	//READING MAP

	int i = 0;
	int xbin = 0;
	int ybin = 0;
	int zbin = 0;
	while(fin>>x0>>x1>>x2>>x3>>x4>>x5>>x6>>x7>>x8>>x9>>x10>>x11){
		
		//cout<<xbin<<" "<<ybin<<" "<<zbin<<" "<<x6<<" "<<i<<endl;
		
		polmapin[xbin][ybin][zbin][0] = x6;
		polmapin[xbin][ybin][zbin][1] = x8;
		polmapin[xbin][ybin][zbin][2] = x10;

		zbin++;

		if(zbin%(nbinsz+1) == 0){
			zbin = 0;
			ybin++;
		}
		if(ybin%(nbinsy+1) == 0 && zbin%(nbinsz+1) == 0){
			ybin = 0;
			xbin++;								
		}
		i++;
	}
	
	if(debug)
		cout << "Found " << i << " corners." << endl << endl;
	
	//INDEXING BINS FOR FASTER SEARCH
	
	std::vector<float> xcornerlist;
	for(int i = 0; i < nbinsx; i++){
		xcornerlist.push_back(xmin+(xstep*i));
	}

	std::vector<float> ycornerlist;
	for(int i = 0; i < nbinsy; i++){
		ycornerlist.push_back(ymin+(ystep*i));
	}

	std::vector<float> zcornerlist;
	for(int i = 0; i < nbinsz; i++){
		zcornerlist.push_back(zmin+(zstep*i));
	}
	
	//FINDING NEAREST CORNERS. STARTS BY BISECTING EACH LIST.
	
	int xhitindex = 0;
	if(xin>=(xmax-xmin)/2){
		for(int i = 100; i < nbinsx; i++){
			if(xin<xcornerlist[i]){
				xhitindex = i-1;
				break;
			}
		}
	}else{
		for(int i = 0; i < nbinsx; i++){
			if(xin<xcornerlist[i]){
				xhitindex = i-1;
				break;
			}
		}
	}
	
	int yhitindex = 0;
	if(yin>=(ymax-ymin)/2){
		for(int i = 12; i < nbinsy; i++){
			if(yin<ycornerlist[i]){
				yhitindex = i-1;
				break;
			}
		}
	}else{
		for(int i = 0; i < nbinsy; i++){
			if(yin<ycornerlist[i]){
				yhitindex = i-1;
				break;
			}
		}
	}
	
	int zhitindex = 0;
	if(zin>=(zmax-zmin)/2){
		for(int i = 7; i < nbinsz; i++){
			if(zin<zcornerlist[i]){
				zhitindex = i-1;
				break;
			}
		}
	}else{
		for(int i = 0; i < nbinsz; i++){
			if(zin<zcornerlist[i]){
				zhitindex = i-1;
				break;
			}
		}
	}
	
	if(debug)
		cout << "Chose corners:" << endl << xcornerlist[xhitindex] << " " << xcornerlist[xhitindex+1] << endl << ycornerlist[yhitindex] << " " << ycornerlist[yhitindex+1] << endl << zcornerlist[zhitindex] << " " << zcornerlist[zhitindex+1] << endl << endl;
	
	//FINDING DISTANCE TO EACH CORNER FOR WEIGHTED AVERAGE
	
	float acpol = 0;
	float alpol = 0;
	float zpol = 0;
	
	float distance[8];
	
	float pols[8][3];
	
	distance[0] = sqrt(pow(xcornerlist[xhitindex]-xin,2)+pow(ycornerlist[yhitindex]-yin,2)+pow(zcornerlist[zhitindex]-zin,2));
	distance[1] = sqrt(pow(xcornerlist[xhitindex+1]-xin,2)+pow(ycornerlist[yhitindex]-yin,2)+pow(zcornerlist[zhitindex]-zin,2));
	distance[2] = sqrt(pow(xcornerlist[xhitindex]-xin,2)+pow(ycornerlist[yhitindex+1]-yin,2)+pow(zcornerlist[zhitindex]-zin,2));
	distance[3] = sqrt(pow(xcornerlist[xhitindex]-xin,2)+pow(ycornerlist[yhitindex]-yin,2)+pow(zcornerlist[zhitindex+1]-zin,2));
	distance[4] = sqrt(pow(xcornerlist[xhitindex+1]-xin,2)+pow(ycornerlist[yhitindex+1]-yin,2)+pow(zcornerlist[zhitindex]-zin,2));
	distance[5] = sqrt(pow(xcornerlist[xhitindex]-xin,2)+pow(ycornerlist[yhitindex+1]-yin,2)+pow(zcornerlist[zhitindex+1]-zin,2));
	distance[6] = sqrt(pow(xcornerlist[xhitindex+1]-xin,2)+pow(ycornerlist[yhitindex]-yin,2)+pow(zcornerlist[zhitindex+1]-zin,2));
	distance[7] = sqrt(pow(xcornerlist[xhitindex+1]-xin,2)+pow(ycornerlist[yhitindex+1]-yin,2)+pow(zcornerlist[zhitindex+1]-zin,2));
	
	for(int i = 0; i < 8; i++){
		if(distance[i] == 0){
			distance[i] = 0.000001;
		}
	}
	
	//LOADING POLARIZATIONS AT CORNERS
	
	for(int i = 0; i < 3; i++){
		pols[0][i] = polmapin[xhitindex][yhitindex][zhitindex][i];
		pols[1][i] = polmapin[xhitindex+1][yhitindex][zhitindex][i];
		pols[2][i] = polmapin[xhitindex][yhitindex+1][zhitindex][i];
		pols[3][i] = polmapin[xhitindex][yhitindex][zhitindex+1][i];
		pols[4][i] = polmapin[xhitindex+1][yhitindex+1][zhitindex][i];
		pols[5][i] = polmapin[xhitindex][yhitindex+1][zhitindex+1][i];
		pols[6][i] = polmapin[xhitindex+1][yhitindex][zhitindex+1][i];
		pols[7][i] = polmapin[xhitindex+1][yhitindex+1][zhitindex+1][i];
	}
	
	//AVERAGING POLARIZATION FOR ALL CORNERS
	
	float num = 0;
	float denom =0;
	for(int i = 0; i < 8; i++){
		num+=pols[i][0]*(1/distance[i]);
		denom+=(1/distance[i]);
	}
	acpol = num/denom;
	
	num = 0;
	denom =0;
	for(int i = 0; i < 8; i++){
		num+=pols[i][1]*(1/distance[i]);
		denom+=(1/distance[i]);
	}
	alpol = num/denom;
	
	num = 0;
	denom =0;
	for(int i = 0; i < 8; i++){
		num+=pols[i][2]*(1/distance[i]);
		denom+=(1/distance[i]);
	}
	zpol = num/denom;
	
	if(debug)
		cout << acpol << " " << alpol << " " << zpol << endl << endl;
	
	pac = acpol;
	pal = alpol;
	pz = zpol;
}
