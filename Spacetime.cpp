#include "Spacetime.h"

#define PI 4*atan(1.0)

Spacetime::Spacetime(double tol_){
	tolerance = tol_;
	phi = 0; // used for winding number.
	has_plane_coordinates = false;
	radius = 0;
	is_wormhole=false;
	is_collapse=false;
	has_dustcloud=false;
}

void Spacetime::recompute_orthonormal_frame(double o[],double p[],double u[],double e1_new[],double e2_new[],double e3_new[]){

	e1_new[0]=o[0]-p[0];
	e1_new[1]=o[1]-p[1];
	e1_new[2]=o[2]-p[2];
	double L = sqrt(e1_new[0]*e1_new[0]+e1_new[1]*e1_new[1]+e1_new[2]*e1_new[2]);
	e1_new[0]/=L;
	e1_new[1]/=L;
	e1_new[2]/=L;
	
	e2_new[0]=e1_new[1]*u[2]-e1_new[2]*u[1];
	e2_new[1]=e1_new[2]*u[0]-e1_new[0]*u[2];
	e2_new[2]=e1_new[0]*u[1]-e1_new[1]*u[0];

	double len=sqrt(e2_new[0]*e2_new[0]+e2_new[1]*e2_new[1]+e2_new[2]*e2_new[2]);
	e2_new[0]=e2_new[0]/len;e2_new[1]=e2_new[1]/len;e2_new[2]=e2_new[2]/len;

	e3_new[0]=e1_new[1]*e2_new[2]-e1_new[2]*e2_new[1];
	e3_new[1]=e1_new[2]*e2_new[0]-e1_new[0]*e2_new[2];
	e3_new[2]=e1_new[0]*e2_new[1]-e1_new[1]*e2_new[0];	

}

void Spacetime::recompute_orthonormal_frame2(double o[],double p[],double e1_new[],double e2_new[],double e3_new[]){

	e1_new[0]=o[0]-p[0];
	e1_new[1]=o[1]-p[1];
	e1_new[2]=o[2]-p[2];
	double L = sqrt(e1_new[0]*e1_new[0]+e1_new[1]*e1_new[1]+e1_new[2]*e1_new[2]);
	e1_new[0]/=L;
	e1_new[1]/=L;
	e1_new[2]/=L;
	
	orthonormalBasis(e1_new,e2_new,e3_new);
	
	//double e3_new[3] = {e1_new[1] * v[2] - e1_new[2] * v[1],e1_new[2] * v[0] - e1_new[0] * v[2],e1_new[0] * v[1] - e1_new[1] * v[0]};
	
	//double L = sqrt(e3_new[0]*e3_new[0]+e3_new[1]*e3_new[1]+e3_new[2]*e3_new[2]);
	//e1_new[0]/=L;
	//e1_new[1]/=L;
	//e1_new[2]/=L;

}

void Spacetime::orthonormalBasis(double e1[3], double e2[3], double e3[3]) {
    // Compute the norm of e1 and check for degeneracy.
    double norm_e1 = std::sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    if (norm_e1 < 1e-10) {
        std::cerr << "Error: e1 is the zero vector.\n";
        return;
    }
    
    // Normalize e1 (store in n).
    double n[3] = { e1[0] / norm_e1, e1[1] / norm_e1, e1[2] / norm_e1 };

    // Choose an auxiliary vector 'a' that is not parallel to n.
    // We pick the coordinate axis corresponding to the smallest absolute component of n.
    double a[3];
    if (std::fabs(n[0]) <= std::fabs(n[1]) && std::fabs(n[0]) <= std::fabs(n[2])) {
        // n[0] is the smallest component in magnitude: choose a = (1,0,0)
        a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;
    } else if (std::fabs(n[1]) <= std::fabs(n[0]) && std::fabs(n[1]) <= std::fabs(n[2])) {
        // n[1] is smallest: choose a = (0,1,0)
        a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;
    } else {
        // Otherwise, n[2] is smallest: choose a = (0,0,1)
        a[0] = 0.0; a[1] = 0.0; a[2] = 1.0;
    }

    // Compute e2 = n x a (cross product).
    e2[0] = n[1]*a[2] - n[2]*a[1];
    e2[1] = n[2]*a[0] - n[0]*a[2];
    e2[2] = n[0]*a[1] - n[1]*a[0];

    // Normalize e2.
    double norm_e2 = std::sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
    if (norm_e2 < 1e-10) {
        std::cerr << "Error: Failed to compute a valid e2.\n";
        return;
    }
    e2[0] /= norm_e2; e2[1] /= norm_e2; e2[2] /= norm_e2;

    // Compute e3 = n x e2.
    e3[0] = n[1]*e2[2] - n[2]*e2[1];
    e3[1] = n[2]*e2[0] - n[0]*e2[2];
    e3[2] = n[0]*e2[1] - n[1]*e2[0];

    // Normalize e3 (should already be unit length, but we normalize to be safe).
    double norm_e3 = std::sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);
    if (norm_e3 < 1e-10) {
        std::cerr << "Error: Failed to compute a valid e3.\n";
        return;
    }
    e3[0] /= norm_e3; e3[1] /= norm_e3; e3[2] /= norm_e3;
}

bool Spacetime::take_step(double P[],double V[],double P2[],double V2[],double &E,double &H,bool &step_taken, bool force_step){
		
			// initial position and velocity as a single vector
			double* y=new double[7];
			y[0]=P[0];y[1]=P[1];y[2]=P[2];y[3]=V[0];y[4]=V[1];y[5]=V[2];y[6]=P[3];//y[6]=E;

			int n=6; // length of y
		
			// Adaptive_Runge_Kutta needs an H_max
			
			// we're going to assume P has a 4th component which is positive if we are outside of the dust cloud and negative otherwise.
		
		
			bool in_universe=Adaptive_Step(y,step_taken,n,E,H,this->tolerance,force_step);
			
			//if(H>0.01*(Universe->u_max-Universe->u_min)){H=0.01*(Universe->u_max-Universe->u_min);}
			
			// ughh.  Ok, in_universe will be false if anywhere in Adaptive_Runge_Kutta we tried to call accel_1st_order somewhere outside the universe,
			// but it may be that y on exit is outside the universe, and the format of the adaptive runge kutta function does not allow us to check this
			// inside it
			
			int jump;
			
			if(in_universe){jump=this->inside_universe(y);}

			if(jump>=2 && step_taken && is_wormhole){
				int N = jump-2; // jump index
				int i,j;
				decode(N, i, j);
				//cout<<"i="<<i<<" "<<"j="<<j<<endl;
				//cout<<"y pre = "<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl;
				teleport(i,j,y);
				if(i==0){
					E = -fabs(E);
				}
				if(i==1){
					E = -10000-fabs(E);
				}

				//cout<<"y post = "<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl;				
			}
			//if(jump==2 && is_collapse){
			//	E=-fabs(E);
			//}

			// update particle position and velocity
			P2[0]=y[0];P2[1]=y[1];P2[2]=y[2];P2[3]=P[3];
			V2[0]=y[3];V2[1]=y[4];V2[2]=y[5];//V2[3]=V[3];	
			
			//if(P[3]<0 && P2[3]<0 && V[0]==-1 && !(V2[0]==-1)){
			//	cout<<"when we took a step, we somehow got V[0] V2[0] = "<<V[0]<<" "<<V2[0]<<endl;
			//}

			if(has_dustcloud && step_taken){
				if(P[3]>0){ // means we began outside of the dust cloud.
				
					double tol = 1e-4;
					int maxIters = 1000;
					double crossing_time = findDustBoundaryCrossingBisection(P,V,P2,V2,tol,maxIters); // will return -1 if there is no crossing, otherwise is in [0,1]
					bool inside = inside_dust_cloud(P2[1],P2[0]);
					double dust_t = dust_time(P2[1]);
					//cout<<"inside = "<<inside<<endl;
					//cout<<"r = "<<P2[1]<<"dust_time(r) = "<<dust_t<<endl;
					
				
					if(crossing_time>=0 && inside){ // means a crossing was detected
						//cout<<"entering gas cloud! crossing_time = "<<crossing_time<<endl;
						//cout<<" P = "<<P[0]<<" "<<P[1]<<" "<<P[2]<<" "<<P[3]<<endl;
						//cout<<" P2 = "<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<" "<<P[3]<<endl;
						//cout<<" P = "<<P[0]<<" "<<P[1]<<" "<<P[2]<<endl;
						
						double P_crossing[] = {P[0]*(1-crossing_time)+P2[0]*crossing_time,P[1]*(1-crossing_time)+P2[1]*crossing_time,P[2]*(1-crossing_time)+P2[2]*crossing_time,fabs(P[3])*(1-crossing_time)+fabs(P2[3])*crossing_time};
						double V_crossing[] = {V[0]*(1-crossing_time)+V2[0]*crossing_time,V[1]*(1-crossing_time)+V2[1]*crossing_time,V[2]*(1-crossing_time)+V2[2]*crossing_time};
						//cout<<" P_crossing = "<<P_crossing[0]<<" "<<P_crossing[1]<<" "<<P_crossing[2]<<" "<<P_crossing[3]<<endl;
						//cout<<"entering gas cloud at r = "<<P_crossing[1]<<endl;
						//cout<<"2M = "<<2*M<<endl;
						//cout<<"redshift factor of (1-2M/r) = "<<(1-2*M/P_crossing[1])<<endl;
						//cout<<"total red_shift = "<<get_redshift(P_crossing[1]);
						
						double P2_interior[4]={0,0,0,0};
						double V2_interior[3]={0,0,0};
						to_interior_coordinates(P_crossing,V_crossing,P2_interior,V2_interior);
						P2[0]=P2_interior[0];
						P2[1]=P2_interior[1];
						P2[2]=P2_interior[2];
						//P2[3]=-P_crossing[1]; // encodes that we are now inside.
						//P2[3]=-1;
						P2[3]=P2_interior[3];
						if(P2[3]<0){
							P2[3]=-get_redshift(P_crossing[1]);
						}
						//cout<<" P2_interior = "<<P2_interior[0]<<" "<<P2_interior[1]<<" "<<P2_interior[2]<<" "<<P2_interior[3]<<endl;
						// let's print a(tau)
						//cout<<"interior coordinates = "<<P2[0]<<" "<<P2[1]<<endl;
						//P2[3]=P2_interior[3];
						V2[0]=V2_interior[0];
						V2[1]=V2_interior[1];
						V2[2]=V2_interior[2];
						H*=crossing_time;
					}
				}
				else if(P[3]<0){ // means we began inside the dust cloud.
				
					/**double eps = 1e-5; // your chosen threshold
					// If X = P2[1] is small or negative, do reflection
					if(P2[1] < eps) {
						// reflection bounce
						//std::cout<<"Bouncing at the center!\n";
						P2[1] = eps; // set X to a small positive
						V2[1] = -V2[1]; // flip radial velocity

						P2[2] += M_PI;  // shift phi by pi
						V2[2] = -V2[2]; // flip angular velocity
					}**/
						
				
					if(P[0]>=0 && P2[0]<0){
						double tol = 1e-4;
						int maxIters = 1000;
						double crossing_time = findTauZeroCrossing(P,P2,tol,maxIters);
						//double crossing_time = findEtaZeroCrossing(P,P2,tol,maxIters);
						if(crossing_time>=0){
							//cout<<"exiting gas cloud through tau=0!"<<endl;
							double P_crossing[] = {P[0]*(1-crossing_time)+P2[0]*crossing_time,P[1]*(1-crossing_time)+P2[1]*crossing_time,P[2]*(1-crossing_time)+P2[2]*crossing_time,-fabs(P[3])*(1-crossing_time)-fabs(P2[3])*crossing_time};
							double V_crossing[] = {V[0]*(1-crossing_time)+V2[0]*crossing_time,V[1]*(1-crossing_time)+V2[1]*crossing_time,V[2]*(1-crossing_time)+V2[2]*crossing_time};
							double P2_exterior[4]={0,0,0,0};
							double V2_exterior[3]={0,0,0};
							//to_exterior_coordinates_through_tau_equals_zero(P_crossing,V_crossing,P2_exterior,V2_exterior);
							to_exterior_coordinates_through_tau_equals_zero(P_crossing,V_crossing,P2_exterior,V2_exterior);
							P2[0]=P2_exterior[0];
							P2[1]=P2_exterior[1];
							P2[2]=P2_exterior[2];
							P2[3]=P2_exterior[3]; // encodes that we are now outside.
							V2[0]=V2_exterior[0];
							V2[1]=V2_exterior[1];
							V2[2]=V2_exterior[2];
							H*=crossing_time;
						}
					}
					else{
							
						double tol = 1e-4;
						int maxIters = 1000;
						double crossing_time = findDustExitCrossingBisection(P,V,P2,V2,tol,maxIters); // will return -1 if there is no crossing, otherwise is in [0,1]					

						if(crossing_time>=0){ // means a crossing was detected
							//cout<<"exiting gas cloud! crossing_time = "<<crossing_time<<endl;
							//cout<<" P = "<<P[0]<<" "<<P[1]<<" "<<P[2]<<" "<<P[3]<<endl;
							double P_crossing[] = {P[0]*(1-crossing_time)+P2[0]*crossing_time,P[1]*(1-crossing_time)+P2[1]*crossing_time,P[2]*(1-crossing_time)+P2[2]*crossing_time,-fabs(P[3])*(1-crossing_time)-fabs(P2[3])*crossing_time};
							double V_crossing[] = {V[0]*(1-crossing_time)+V2[0]*crossing_time,V[1]*(1-crossing_time)+V2[1]*crossing_time,V[2]*(1-crossing_time)+V2[2]*crossing_time};
							double P2_exterior[4]={0,0,0,0};
							double V2_exterior[3]={0,0,0};
							to_exterior_coordinates(P_crossing,V_crossing,P2_exterior,V2_exterior);
							P2[0]=P2_exterior[0];
							P2[1]=P2_exterior[1];
							P2[2]=P2_exterior[2];
							P2[3]=P2_exterior[3]; // encodes that we are now outside.
							V2[0]=V2_exterior[0];
							V2[1]=V2_exterior[1];
							V2[2]=V2_exterior[2];
							//cout<<" P_crossing = "<<P_crossing[0]<<" "<<P_crossing[1]<<" "<<P_crossing[2]<<" "<<P_crossing[3]<<endl;
							//cout<<" P2_exterior = "<<P2_exterior[0]<<" "<<P2_exterior[1]<<" "<<P2_exterior[2]<<" "<<P2_exterior[3]<<endl;
							H*=crossing_time;
						}
					}
				}
			}
			
			
			delete[] y;		

			/**if(has_plane_coordinates){
				double dx = (P2[0]-P[0])*e1[0]+(P2[1]-P[1])*e1[1]+(P2[2]-P[2])*e1[2];
				double dy = (P2[0]-P[0])*e2[0]+(P2[1]-P[1])*e2[1]+(P2[2]-P[2])*e2[2];
				double dphi = (cos(phi)*dy-sin(phi)*dx)/radius;
				double dr = cos(phi)*dx+sin(phi)*dy;
				phi = phi+dphi;
				radius = radius+dr;
			}**/
			
			return in_universe;	
	}
	
bool Spacetime::Adaptive_Step(double*& y,bool& y_changed,int n,double E,double& h,double tol,bool force_step){
		double* k1=new double[n+1];
		double* k2=new double[n+1];
		double* k3=new double[n+1];
		double* k4=new double[n+1];
		double* k5=new double[n+1];
		double* k6=new double[n+1];
		double y_n_backup=y[n];
		k1[n]=y[n];
		k2[n]=y[n];
		k3[n]=y[n];
		k4[n]=y[n];
		k5[n]=y[n];
		k6[n]=y[n];
		
		//double h_min = 1.0;
		//if( danger_zone(y) ){
		//	h_min=0.1;
		//}
		double h_min = min_step_size(y);
		
		double y1[n+1],Delta4[n],Delta5[n],diff54[n];
		
		bool inside_universe = true;
		
		if( !accel_1st_order(k1,y,E) ){inside_universe = false;} // attempts to evaluate f and returns false if y is outside Domain(f)
		for(int i=0;i<n;i++){k1[i]*=h;  
							 y1[i]=y[i]+0.25*k1[i];}
							 y1[n]=y_n_backup;
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !accel_1st_order(k2,y1,E) ){inside_universe = false;}
		for(int i=0;i<n;i++){k2[i]*=h;  
							 y1[i]=y[i]+(3.0/32)*k1[i]+(9.0/32)*k2[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;					 
		if( !accel_1st_order(k3,y1,E) ){inside_universe = false;}
		for(int i=0;i<n;i++){k3[i]*=h;  
							 y1[i]=y[i]+(1932.0/2197)*k1[i]-(7200.0/2197)*k2[i]+(7296.0/2197)*k3[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !accel_1st_order(k4,y1,E) ){inside_universe = false;}
		for(int i=0;i<n;i++){k4[i]*=h;  
							 y1[i]=y[i]+(439.0/216)*k1[i]-8*k2[i]+(3680.0/513)*k3[i]-(845.0/4104)*k4[i];}		
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !accel_1st_order(k5,y1,E) ){inside_universe = false;}
		for(int i=0;i<n;i++){k5[i]*=h;  
							 y1[i]=y[i]-(8.0/27)*k1[i]+2*k2[i]-(3544.0/2565)*k3[i]+(1859.0/4104)*k4[i]-(11.0/40)*k5[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;				 
		if( !accel_1st_order(k6,y1,E) ){inside_universe = false;}	
		for(int i=0;i<n;i++){k6[i]*=h;
							 Delta4[i]=(25.0/216)*k1[i]+(1408.0/2565)*k3[i]+(2197.0/4104)*k4[i]-(1.0/5)*k5[i];
							 Delta5[i]=(16.0/135)*k1[i]+(6656.0/12825)*k3[i]+(28561.0/56430)*k4[i]-(9.0/50)*k5[i]+(2.0/55)*k6[i];
							 diff54[i]=Delta5[i]-Delta4[i];
			}
		
		/**if(y_n_backup<0){
			Delta4[0]=-h;
			Delta5[0]=-h;

			Delta4[3]=0;
			Delta5[3]=0;

			diff54[0]=0;
			diff54[3]=0;
		}**/

		double q,eh;

		eh=error_norm(y,diff54)/h;

		if(eh==0){q=2;}
		else{q=pow(tol/eh,0.25);}
	
		h=h*q;
		if(true){
			
			h = min(h,h_min);
		}
		// accept the step if the estimated truncation error doesn't exceed the tolerance by too much
		
		if(eh<1.1*tol || force_step){
			double L_Delta5 = sqrt(Delta5[0]*Delta5[0]+Delta5[1]*Delta5[1]+Delta5[2]*Delta5[2]);
			//if(L_Delta5<=h){
				for(int i=0;i<n;i++){y[i]+=Delta5[i];}
				y[n]=y_n_backup;
				//if(y_n_backup<0){y[3]=-1;}
			//}
			//else{
			//	for(int i=0;i<n;i++){y[i]+=Delta5[i]*h/L_Delta5;}
			//}
			y_changed=true;
		}
		else{ y_changed=false;}
		
		delete[] k1;delete[] k2;delete[] k3;delete[] k4;delete[] k5;delete[] k6;

		return inside_universe;		
	}


