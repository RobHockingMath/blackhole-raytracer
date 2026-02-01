#include "Adaptive_Runge_Kutta.h"
#include <math.h>
#include <iostream>

namespace Adaptive_Runge_Kutta
{
	bool Adaptive_Step(double*& y,bool& y_changed,int n,double& h,double tol,bool F(double*& yp,double *y,void * dummy),double norm(double *x,double *y,void * dummy),void * dummy){
		double* k1=new double[n];
		double* k2=new double[n];
		double* k3=new double[n];
		double* k4=new double[n];
		double* k5=new double[n];
		double* k6=new double[n];
		
		double y1[n],Delta4[n],Delta5[n],diff54[n];
		
		if( !F(k1,y,dummy) ){return false;} // attempts to evaluate f and returns false if y is outside Domain(f)
		for(int i=0;i<n;i++){k1[i]*=h;  
							 y1[i]=y[i]+0.25*k1[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !F(k2,y1,dummy) ){return false;}
		for(int i=0;i<n;i++){k2[i]*=h;  
							 y1[i]=y[i]+(3.0/32)*k1[i]+(9.0/32)*k2[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;					 
		if( !F(k3,y1,dummy) ){return false;}
		for(int i=0;i<n;i++){k3[i]*=h;  
							 y1[i]=y[i]+(1932.0/2197)*k1[i]-(7200.0/2197)*k2[i]+(7296.0/2197)*k3[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !F(k4,y1,dummy) ){return false;}
		for(int i=0;i<n;i++){k4[i]*=h;  
							 y1[i]=y[i]+(439.0/216)*k1[i]-8*k2[i]+(3680.0/513)*k3[i]-(845.0/4104)*k4[i];}		
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !F(k5,y1,dummy) ){return false;}
		for(int i=0;i<n;i++){k5[i]*=h;  
							 y1[i]=y[i]-(8.0/27)*k1[i]+2*k2[i]-(3544.0/2565)*k3[i]+(1859.0/4104)*k4[i]-(11.0/40)*k5[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;				 
		if( !F(k6,y1,dummy) ){return false;}	
		for(int i=0;i<n;i++){k6[i]*=h;
							 Delta4[i]=(25.0/216)*k1[i]+(1408.0/2565)*k3[i]+(2197.0/4104)*k4[i]-(1.0/5)*k5[i];
							 Delta5[i]=(16.0/135)*k1[i]+(6656.0/12825)*k3[i]+(28561.0/56430)*k4[i]-(9.0/50)*k5[i]+(2.0/55)*k6[i];
							 diff54[i]=Delta5[i]-Delta4[i];
			}
		double q,eh;

		eh=norm(y,diff54,dummy)/h;

		if(eh==0){q=2;}
		else{q=pow(tol/eh,0.25);}

	
		h=h*q;
		// accept the step if the estimated truncation error doesn't exceed the tolerance by too much
		
		if(eh<1.1*tol){
			for(int i=0;i<n;i++){y[i]+=Delta5[i];}
			y_changed=true;
		}
		else{ y_changed=false;}
		
		delete[] k1;delete[] k2;delete[] k3;delete[] k4;delete[] k5;delete[] k6;

		return true;		
	}
}