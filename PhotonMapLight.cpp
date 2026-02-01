
#include <math.h>
#include <iostream>
using namespace std;
#include "PhotonMapLight.h"
#include <random>

#define PI 3.14159265359

PhotonMapLight::PhotonMapLight(double b,int C[]){
	color[0]=C[0];
	color[1]=C[1];
	color[2]=C[2];
	power=b;
}



spherePhotonMapLight::spherePhotonMapLight(double P[], double b, double R_,int C[])
	: PhotonMapLight(b,C)
{
	pos[0]=P[0];pos[1]=P[1];pos[2]=P[2];
	R=R_;
}

void spherePhotonMapLight::get_Photons(int n_photons, vector< vector<double> >* P, vector< vector<double> >* V){
	if(R<=0){
		// Define the range [min, max]
		double min = -1.0;
		double max = 1.0;

		P->clear();
		V->clear();
		
		vector<double> position;
		position.push_back(pos[0]);
		position.push_back(pos[1]);
		position.push_back(pos[2]);

		// Create a random number generator
		std::random_device rd;  // Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> distrib(min, max);
		vector<double> P1;
		P1.clear();
		int count=0;
		while(count<n_photons){
			
			double x = distrib(gen);
			double y = distrib(gen);
			double z = distrib(gen);
			
			double r2 = x*x+y*y+z*z;
			
			if(0<r2 && r2<=1){
				double r1 = sqrt(r2);
				x = x/r1;
				y = y/r1;
				z = z/r1;
				P1.push_back(x);
				P1.push_back(y);
				P1.push_back(z);
				P->push_back(position);
				V->push_back(P1);
				P1.clear();
				count++;
			}
			
			
		}
	}
	else{
		// Define the range [min, max]
		double min = -1.0;
		double max = 1.0;

		P->clear();
		V->clear();

		// Create a random number generator
		std::random_device rd;  // Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> distrib(min, max);
		vector<double> P1,V1;
		P1.clear();
		V1.clear();
		int count=0;
		while(count<n_photons){
			
			double x = distrib(gen);
			double y = distrib(gen);
			double z = distrib(gen);
			
			double r2 = x*x+y*y+z*z;
			
			if(0<r2 && r2<=1){
				double r1 = sqrt(r2);
				x = x/r1;
				y = y/r1;
				z = z/r1;
				P1.push_back(pos[0]+R*x);
				P1.push_back(pos[1]+R*y);
				P1.push_back(pos[2]+R*z);
				P->push_back(P1);
				P1.clear();
				count++;
			}
			
			double e1 = (1+distrib(gen))/2.0; // uniform [0,1]
			double e2 = (1+distrib(gen))/2.0; // uniform [0,1]
			double theta = acos(sqrt(e1));
			double phi = 2*PI*e2;
			
			double x1,y1,z1;
			if(fabs(z)<0.5){
				x1=-y;y1=x;z1=0;
			}
			else{
				x1=0;y1=-z;z1=y;
			}
			double length = sqrt(x1*x1+y1*y1+z1*z1);
			x1=x1/length;
			y1=y1/length;
			z1=z1/length;
			double x2 = y*z1-y1*z;
			double y2 = z*x1-z1*x;
			double z2 = x*y1-x1*y;
			
			double vx = x*cos(theta)+x1*sin(theta)*sin(phi)+x2*sin(theta)*cos(phi);
			double vy = y*cos(theta)+y1*sin(theta)*sin(phi)+y2*sin(theta)*cos(phi);
			double vz = z*cos(theta)+z1*sin(theta)*sin(phi)+z2*sin(theta)*cos(phi);
			V1.push_back(vx);
			V1.push_back(vy);
			V1.push_back(vz);
			V->push_back(V1);
			V1.clear();
			
			
		}
	}
}