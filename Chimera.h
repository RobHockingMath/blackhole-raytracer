#ifndef CHIMERA_H
#define CHIMERA_H

#include <vector>
#include "Camera.h"

using namespace std;

class ChimeraCameraArray {
public:
	ChimeraCameraArray(double R, double H_m, double H_M, double UP[3],double arclength, double holo_W, double holo_H, double hogel_size,double zoom);
	virtual void get_rays(int i,int j, int I, int J, const int, const int,vector< vector<double> >*,vector< vector<double> >*){return;};
	vector<vector<Pinhole_Camera*>> cameras;
	double p[3],up[3],l[3];

	double holoW, holoH, hogelSize; //holoW and holoH are in cm.  hogelSize is in micrometers.
	int w,h; // the width and height in pixels of the images created by the array.

};



#endif

