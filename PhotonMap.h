#ifndef PHOTON_MAP_H
#define PHOTON_MAP_H


#include "Shape.h"
#include "invisibility.h"
#include <vector>
#include "photon.hpp"

class PhotonMap{
public:
	PhotonMap(int num_photons, int num_objects);
	void storePhoton(double x[], double v[], int object_id, double power);
protected:
	int num_photons;
	int num_objects;
};



#endif