#include "multi_worm_hole.h"

Multi_worm_hole::Multi_worm_hole(vector<array<double, 3>> positions_,double tol_, vector<double> as_,vector<int> stitching_)
	: Spacetime(tol_), positions(positions_), as(as_), stitching(stitching_), num_wormholes(positions_.size())
{
}

int Multi_worm_hole::inside_universe(double *p){
	
	// it's important to know that *p is actually a vector with six components, [p,v]
	
	//iterate over all black holes
	
	for(int i=0;i<num_wormholes;i++){
	
		double r_i=sqrt((p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]));
		
		double a_i = as[i];
		
		if(r_i<=a_i){
			int j = stitching[i];
			//teleport(i,j,p); // means this is a wormhole.
			return 2+encode(i,j); // flag for teleporting
		}
		
	}
	return 1; // we've neither fallen into a worm hole nor teleported.
}

void Multi_worm_hole::teleport(int i,int j,double *p){
	// move from worm hole i to j.
	double dp_i[3] = {p[0]-positions[i][0],p[1]-positions[i][1],p[2]-positions[i][2]};
	
	double delta = sqrt(dp_i[0]*dp_i[0]+dp_i[1]*dp_i[1]+dp_i[2]*dp_i[2]);
	dp_i[0]=dp_i[0]/delta;
	dp_i[1]=dp_i[1]/delta;
	dp_i[2]=dp_i[2]/delta;
	
	double a_j = as[j];
	
	p[0] = positions[j][0]-scale_a*a_j*dp_i[0];
	p[1] = positions[j][1]-scale_a*a_j*dp_i[1];
	p[2] = positions[j][2]-scale_a*a_j*dp_i[2];
	
//	p[3]=-dp_i[0];
//	p[4]=-dp_i[1];
//	p[5]=-dp_i[2];

	
}

double Multi_worm_hole::get_h0(double *p){
	
	//iterate over all black holes
	
	double h0 = 0.1;
	
	for(int i=0;i<num_wormholes;i++){
	
		double r_i=sqrt((p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]));
		
		double a_i = as[i];
		
		if(fabs(1-r_i/a_i)<h0){
			h0=fabs(1-r_i/a_i);
		}
	}
	return h0;

}

// We think of p as being a point and x a vector.  
// This format is useful as it can be used directly in adaptive Runge Kutta routines.

double Multi_worm_hole::error_norm(double *p,double *x){
	//return sqrt((x[0]-pos[0])*(x[0]-pos[0])+(x[1]-pos[1])*(x[1]-pos[1])+(x[2]-pos[2])*(x[2]-pos[2]));
	return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]));
	//double r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	//return (fabs(x[0])+fabs(x[1])+fabs(x[2]))/fabs(1-r/M);
	
}

double Multi_worm_hole::min_step_size(double* y){
	
	//
	double buffer1 = 6.0;
	double buffer2 = 3.0;
	double buffer3 = 1.5;
	double buffer4 = 0.5;
	double buffer5 = 0.1;
	
	double h_min = 1.0;
	
	for(int i=0;i<num_wormholes;i++){
		
	
		double r_i=sqrt( (y[0]-positions[i][0])*(y[0]-positions[i][0])+(y[1]-positions[i][1])*(y[1]-positions[i][1])+(y[2]-positions[i][2])*(y[2]-positions[i][2]) );
		double a_i = as[i];
		h_min=fmin(h_min,fabs(1-a_i/r_i));
	}
	return h_min;
	/**for(int i=0;i<num_black_holes;i++){
	
		double r_i=sqrt((y[0]-positions[i][0])*(y[0]-positions[i][0])+(y[1]-positions[i][1])*(y[1]-positions[i][1])+(y[2]-positions[i][2])*(y[2]-positions[i][2]));
		double m_i = masses[i];
		double h0=fabs(1-0.9*r_i/m_i);
		if(h0<h_min){
			h_min = h0;
		}
		/**double m_i = masses[i];
		if(r_i<m_i+buffer1){
			h_min = min(h_min,1.0);
		}
		if(r_i<m_i+buffer2){
			h_min = min(h_min,1.0);
		}
		if(r_i<m_i+buffer3){
			h_min = min(h_min,0.1);
		}
		if(r_i<m_i+buffer4){
			h_min = min(h_min,0.01);
		}
		if(r_i<m_i+buffer5){
			h_min = min(h_min,0.001);
			return h_min;
		}**/
	//}
	//return h_min;
}

bool Multi_worm_hole::accel_1st_order(double*& yp,double* y, double E){
	
	
	double dPhi_over_Phi[3] = {0.0,0.0,0.0};
	
	
	for(int i=0;i<num_wormholes;i++){
		
	
			double r_i2=(y[0]-positions[i][0])*(y[0]-positions[i][0])+(y[1]-positions[i][1])*(y[1]-positions[i][1])+(y[2]-positions[i][2])*(y[2]-positions[i][2]);
			
			double r_i = sqrt(r_i2);
			
			double a_i = as[i];

			//double Phi_i = (r_i2-a_i*a_i)/(r_i2);
			double Phi_i = (r_i-a_i)*(r_i+a_i)/(r_i2);
			
			double dPhidx = 2*a_i*a_i*(y[0]-positions[i][0])/(r_i2*r_i2);
			double dPhidy = 2*a_i*a_i*(y[1]-positions[i][1])/(r_i2*r_i2);
			double dPhidz = 2*a_i*a_i*(y[2]-positions[i][2])/(r_i2*r_i2);
			dPhi_over_Phi[0]+=dPhidx/Phi_i;
			dPhi_over_Phi[1]+=dPhidy/Phi_i;
			dPhi_over_Phi[2]+=dPhidz/Phi_i;
		
	}
	
	// y and yp are giant vectors of the form y=[p,v] where p is position and v is velocity.
	
	// position derivative is velocity.
	yp[0]=y[3]; 
	yp[1]=y[4]; 
	yp[2]=y[5]; 
	
	yp[3]=dPhi_over_Phi[0]*y[3];
	yp[4]=dPhi_over_Phi[1]*y[4];
	yp[5]=dPhi_over_Phi[2]*y[5];


	
	return inside_universe(yp);
	
}

double Multi_worm_hole::compute_energy(double *p){

	return 1.0;
}

void Multi_worm_hole::scale_velocity(double v[],double *p){
	
	double Phi=1;
	
	
	for(int i=0;i<num_wormholes;i++){
		
	
			double r_i2=(p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]);
			
			double a_i = as[i];

			double Phi_i = (r_i2-a_i*a_i)/(r_i2);
			
			Phi*=Phi_i;
		
	}
	Phi = sqrt(Phi);
	v[0]=v[0]/Phi;
	v[1]=v[1]/Phi;
	v[2]=v[2]/Phi;
}




