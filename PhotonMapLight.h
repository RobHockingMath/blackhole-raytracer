#ifndef PHOTON_MAP_LIGHT_H
#define PHOTON_MAP_LIGHT_H


#include "Shape.h"
#include "invisibility.h"
#include <vector>


class PhotonMapLight{
public:
	PhotonMapLight(double,int[]);
	virtual void get_Photons(int n_photons, vector< vector<double> >* P, vector< vector<double> >* V){return;};
protected:
	int color[3];
	double power;
};




class spherePhotonMapLight : public PhotonMapLight{
public:
	spherePhotonMapLight(double[],double,double,int[]);
	void get_Photons(int n_photons, vector< vector<double> >* P, vector< vector<double> >* V);
private:
	double pos[3];
	double R;
};

#endif