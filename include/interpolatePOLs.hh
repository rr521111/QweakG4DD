#ifndef interpolatePOLs_hh
#define interpolatePOLs_hh 1

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

class interpolatePOLs
{
private:
  double xPosLowLimit, xPosHighLimit;
  double yPosLowLimit, yPosHighLimit;//xAngLowLimit, xAngHighLimit;
  double yAngLowLimit, yAngHighLimit;//energyLowLimit, energyHighLimit;

  static const int dimension = 6;//3 DoF + 3 POL values
  std::vector<double> scanPoints[dimension];
  std::string polModel;
  
  void getCorners(int lowerIndex, int upperIndex, int depth, std::vector<double> point,
		  std::vector<double> points[dimension]);

  void getPOLs(std::vector<double> in[dimension], std::vector<double> pt,
	      double &outAcross, double &outAlong, double &outZ);

public:
  interpolatePOLs(std::string polmap = "default",int POLuncert = 0);
  ~interpolatePOLs(){};

  void setPolMap(std::string polmap);
  int verbosity;
  int polUncert;
  int getPOLs(double x, double y, double AngY, double &outac, double &outal, double &outz);
  void readScan();

};

#endif
