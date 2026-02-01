#include "multi_black_hole.h"

#define PI 3.14159

Multi_black_hole::Multi_black_hole(vector<array<double, 3>> positions_,double tol_, vector<double> masses_,vector<int> stitching_)
	: Spacetime(tol_), positions(positions_), masses(masses_), stitching(stitching_), num_black_holes(positions_.size())
{
	is_wormhole = true;
}

int Multi_black_hole::inside_universe(double *p){
	
	// it's important to know that *p is actually a vector with six components, [p,v]
	
	//iterate over all black holes
	
	for(int i=0;i<num_black_holes;i++){
	
		double r_i=sqrt((p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]));
		
		double m_i = masses[i];
		
		if(r_i<m_i && p[6]>=0){
			if(stitching[i]<0){
				return 0; // means this is a black hole.
			}
			else{
				int j = stitching[i];
				//teleport(i,j,p); // means this is a wormhole.
				return 2+encode(i,j); // flag for teleporting
			}
		}
		
	}
	return 1; // we've neither fallen into a black hole nor teleported.
}

void Multi_black_hole::rotate(double x[], double v[], double theta){
	double c = cos(theta/2);
	double s = sin(theta/2);
	
	double X = v[0]; double Y = v[1]; double Z = v[2];
	double P0 = (2*(X*X-1)*s*s+1)*(x[0])+(2*X*Y*s*s-2*Z*c*s)*(x[1])+(2*X*Z*s*s+2*Y*c*s)*(x[2]);
	double P1 = (2*X*Y*s*s+2*Z*c*s)*(x[0])+(2*(Y*Y-1)*s*s+1)*(x[1])+(2*Y*Z*s*s-2*X*c*s)*(x[2]);
	double P2 = (2*X*Z*s*s-2*Y*c*s)*(x[0])+(2*Y*Z*s*s+2*X*c*s)*(x[1])+(2*(Z*Z-1)*s*s+1)*(x[2]);
	x[0]=P0;
	x[1]=P1;
	x[2]=P2;
	
}

void Multi_black_hole::teleport(int i,int j,double *p){
	// move from black hole i to j.
	double dp_i[3] = {p[0]-positions[i][0],p[1]-positions[i][1],p[2]-positions[i][2]};
	
	double delta = sqrt(dp_i[0]*dp_i[0]+dp_i[1]*dp_i[1]+dp_i[2]*dp_i[2]);
	dp_i[0]=dp_i[0]/delta;
	dp_i[1]=dp_i[1]/delta;
	dp_i[2]=dp_i[2]/delta;
	
	double rotation_vector[3]={0,1,0};
	//double rotation_angle = -(PI/4);
	
	double rotation_angle = -(PI/2)*(tanh(-5*(dp_i[2]-0.5))+1)/2;
	this->rotate(dp_i,rotation_vector,rotation_angle);
	
	double m_j = masses[j];
	
	double scale = 1.1;
	
	p[0] = positions[j][0]+m_j*dp_i[0];
	p[1] = positions[j][1]+m_j*dp_i[1];
	p[2] = positions[j][2]+m_j*dp_i[2];
	
	//p[3]=-dp_i[0];
	//p[4]=-dp_i[1];
	//p[5]=-dp_i[2];
	
	/**double rotation_angle = (PI/2)*(tanh(-5*(dp_i[2]-0.5))+1)/2;
	
	double rotation_vector[3]={-p[4],p[3],0};
	
	double velocity[3]={p[3],p[4],p[5]};
	
	this->rotate(velocity,rotation_vector,rotation_angle);
	
	p[3]=velocity[0];
	p[4]=velocity[1];
	p[5]=velocity[2];**/
	
	//double polar_angle = acos(dp_i[2]);
	
	//double Vx = p[3]*cos(polar_angle)+p[5]*sin(polar_angle);
	//double Vz = p[3]*sin(polar_angle)-p[5]*cos(polar_angle);
	
	//p[5]=p[5];
	
	//p[3]=Vx;
	//p[5]=Vz;
	
	//p[3]=-dp_i[0];
	//p[4]=-dp_i[1];
	//p[5]=-dp_i[2];
	
	//double speed = sqrt(p[3]*p[3]+p[4]*p[4]+p[5]*p[5]);
	//p[3]=p[3]/speed;
	//p[4]=p[4]/speed;
	//p[5]=p[5]/speed;
	
	/*p[0] = positions[j][0]-scale*m_j*dp_i[0];
	p[1] = positions[j][1]-scale*m_j*dp_i[1];
	p[2] = positions[j][2]-scale*m_j*dp_i[2];
	double speed = sqrt(p[3]*p[3]+p[4]*p[4]+p[5]*p[5]);
	p[3]=p[3]/speed;
	p[4]=p[4]/speed;
	p[5]=p[5]/speed;*/
	
}

double Multi_black_hole::get_h0(double *p){
	
	//iterate over all black holes
	
	double h0 = 0.1;
	
	for(int i=0;i<num_black_holes;i++){
	
		double r_i=sqrt((p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]));
		
		double m_i = masses[i];
		
		if(fabs(1-r_i/m_i)<h0){
			h0=fabs(1-r_i/m_i);
		}
	}
	return h0;

}

// We think of p as being a point and x a vector.  
// This format is useful as it can be used directly in adaptive Runge Kutta routines.

double Multi_black_hole::error_norm(double *p,double *x){
	//return sqrt((x[0]-pos[0])*(x[0]-pos[0])+(x[1]-pos[1])*(x[1]-pos[1])+(x[2]-pos[2])*(x[2]-pos[2]));
	return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]));
	//double r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	//return (fabs(x[0])+fabs(x[1])+fabs(x[2]))/fabs(1-r/M);
	
}

double Multi_black_hole::min_step_size(double* y){
	
	//
	double buffer1 = 6.0;
	double buffer2 = 3.0;
	double buffer3 = 1.5;
	double buffer4 = 0.5;
	double buffer5 = 0.1;
	
	double h_min = 1.0;
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

bool Multi_black_hole::accel_1st_order(double*& yp,double* y, double E){
	
	
	double U = 1.0;
	
	double r0 = 2.0;
	double dUdx = 0.0;
	double dUdy = 0.0;
	double dUdz = 0.0;
	
	double zero_active = true;
	double one_active = true;
	
	if(E<-10000){
		E = 10000+E;
		zero_active = false;
	}
	else if (E<0){
		one_active = false;
	}
	
	for(int i=0;i<num_black_holes;i++){
		
		if((i==0 && zero_active) || (i==1 && one_active)){
	
			double r_i=sqrt((y[0]-positions[i][0])*(y[0]-positions[i][0])+(y[1]-positions[i][1])*(y[1]-positions[i][1])+(y[2]-positions[i][2])*(y[2]-positions[i][2]));
			
			double m_i = masses[i];
			m_i = m_i*tanh(r0-r_i);
			if(m_i<=0 && y[0]>=2){
				m_i=0;
			}
			U=U+m_i/r_i;
			dUdx = dUdx - m_i*(y[0]-positions[i][0])/(r_i*r_i*r_i);
			dUdy = dUdy - m_i*(y[1]-positions[i][1])/(r_i*r_i*r_i);
			dUdz = dUdz - m_i*(y[2]-positions[i][2])/(r_i*r_i*r_i);
		
		}
		
	}
	
	// y and yp are giant vectors of the form y=[p,v] where p is position and v is velocity.
	
	double speed2 = y[3]*y[3]+y[4]*y[4]+y[5]*y[5];
	
	// position derivative is velocity.
	yp[0]=y[3]; 
	yp[1]=y[4]; 
	yp[2]=y[5]; 
	

	
	// acceleration is zero.
	double dU_dot_v = dUdx*y[3]+dUdy*y[4]+dUdz*y[5];
	yp[3]=((E*E+speed2)*dUdx-2*dU_dot_v*y[3])/U;
	yp[4]=((E*E+speed2)*dUdy-2*dU_dot_v*y[4])/U;
	yp[5]=((E*E+speed2)*dUdz-2*dU_dot_v*y[5])/U;


	
	return inside_universe(yp);
	
}

double Multi_black_hole::compute_energy(double *p){

	double U = 1.0;	
	double r0 = 2.0;

	for(int i=0;i<num_black_holes;i++){
	
		double r_i=sqrt((p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]));
		
		double m_i = masses[i];
		U=U+m_i/r_i;		
	}

	return 1.0/U;
}

void Multi_black_hole::scale_velocity(double v[],double *p){
	
	double U = 1.0;	
	double r0 = 2.0;

	for(int i=0;i<num_black_holes;i++){
	
		double r_i=sqrt((p[0]-positions[i][0])*(p[0]-positions[i][0])+(p[1]-positions[i][1])*(p[1]-positions[i][1])+(p[2]-positions[i][2])*(p[2]-positions[i][2]));
		
		double m_i = masses[i];
		U=U+m_i/r_i;		
	}
	v[0]=v[0]/U;
	v[1]=v[1]/U;
	v[2]=v[2]/U;
}




