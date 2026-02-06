#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
#include "World.h"
#include "bmp.h"
#define PI 4*atan(1.0)
#include <omp.h>

using namespace nanoflann;

// Define the KDTree type
typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<float, PhotonCloud>,
    PhotonCloud,
    3 /* dimension */
> PhotonKDTree;

World::World(double R) {
	//invis_R=R;
	invis_R=0;
	env_map_ID=0;
	//this should really be uncommented out later...
	//I should do this properly at some point or I'll be
	//super confused later...
	/**image = new int**[w];
	for(int i=0;i<w;i++){
		image[i]=new int*[h];
		for(int j=0;j<h;j++){
			image[i][j]=new int[3];
		}
	}**/
	
	// Create a random number generator
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> distrib(0.0, 1.0);
}

World::~World() {
	//delete[] image;
}

void World::init(){
	//int i,j;
	/**for(i=0;i<w;i++){
		for(j=0;j<h;j++){
			
			image[i][j][0]=255;
			image[i][j][1]=255;
			image[i][j][2]=255;
			//if(0.25*w<i && i<0.75*w && 0.25*h<j && j<0.75*h){
			//	image[i][j][0]=0;image[i][j][1]=0;image[i][j][2]=0;
			//}
		}
	}**/
}

//void World::writeData(char* filename){
//  save_bmp(w,h,filename,image);
	/**ofstream myfile;
	myfile.open("World.txt");
	myfile << w <<" "<< h << "\n";
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			myfile << i <<" "<<j<<" "<<image[i][j][0]<<" "<<image[i][j][1]<<" "<<image[i][j][2]<<"\n";
		}
	}
	myfile.close();**/
//}

void World::addShape(Shape* S){
	shapes.push_back(S);
}

void World::addLight(Light* L){
	lights.push_back(L);
}

void World::addPhotonMapLight(PhotonMapLight* L){
	photon_lights.push_back(L);
}

void World::addTexture(Texture* T){
	textures.push_back(T);
}

void World::addPhotonMapTexture(Texture* T){
	photon_maps.push_back(T);
}


void World::read_env_map(int i,int j,int C[]){

	int W=textures[env_map_ID]->get_width()-1;
	int H=textures[env_map_ID]->get_height()-1;

	int ii=int(i*float(W)/w);
	int jj=int(j*float(H)/h);
	
	textures[env_map_ID]->get_at(ii,jj,C);
}

void World::read_sphere_env_map(double V[],int C[]){
	double r = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
	double u[2]={0,0};
	u[0]=(atan2(V[1],V[0])+PI)/(2*PI);
	u[1]=1-acos(V[2]/r)/PI;
	if(u[0]<0){u[0]=0;}
	if(u[1]<0){u[1]=0;}
	if(u[0]>1){u[0]=1;}
	if(u[1]>1){u[1]=1;}
	
	int W,H;
	
	//# pragma omp critical
	{
		W=textures[env_map_ID]->get_width();
		H=textures[env_map_ID]->get_height();
	}
	
	int ii=floor(u[0]*float(W-1));
	int ip=min(ii+1,W-1);
	int jj=floor(u[1]*float(H-1));
	int jp=min(jj+1,H-1);
	
	int C00[3]={C[0],C[1],C[2]};
	int C01[3]={C[0],C[1],C[2]};
	int C10[3]={C[0],C[1],C[2]};
	int C11[3]={C[0],C[1],C[2]};
	
	//# pragma omp critical
	{
		textures[env_map_ID]->get_at(ii,jj,C00);
		textures[env_map_ID]->get_at(ip,jj,C10);
		textures[env_map_ID]->get_at(ii,jp,C01);
		textures[env_map_ID]->get_at(ip,jp,C11);
	}
	
	double s = u[0]*float(W-1)-ii;
	double t = u[1]*float(H-1)-jj;
	
	C[0]=int((1-s)*(1-t)*C00[0]+s*(1-t)*C10[0]+(1-s)*t*C01[0]+s*t*C11[0]);
	C[1]=int((1-s)*(1-t)*C00[1]+s*(1-t)*C10[1]+(1-s)*t*C01[1]+s*t*C11[1]);
	C[2]=int((1-s)*(1-t)*C00[2]+s*(1-t)*C10[2]+(1-s)*t*C01[2]+s*t*C11[2]);
	
	return;
	
}

void World::rotate(double x[], double v[], double theta, double o[]){
	double c = cos(theta/2);
	double s = sin(theta/2);
	
	double X = v[0]; double Y = v[1]; double Z = v[2];
	double P0 = (2*(X*X-1)*s*s+1)*(x[0]-o[0])+(2*X*Y*s*s-2*Z*c*s)*(x[1]-o[1])+(2*X*Z*s*s+2*Y*c*s)*(x[2]-o[2]);
	double P1 = (2*X*Y*s*s+2*Z*c*s)*(x[0]-o[0])+(2*(Y*Y-1)*s*s+1)*(x[1]-o[1])+(2*Y*Z*s*s-2*X*c*s)*(x[2]-o[2]);
	double P2 = (2*X*Z*s*s-2*Y*c*s)*(x[0]-o[0])+(2*Y*Z*s*s+2*X*c*s)*(x[1]-o[1])+(2*(Z*Z-1)*s*s+1)*(x[2]-o[2]);
	x[0]=P0+o[0];
	x[1]=P1+o[1];
	x[2]=P2+o[2];
	
}


void World::non_linear_photon_emission(double P[],double V[],vector< vector<Photon>>& photons,double photon_energy,int depth,double H,double E,int last_index,bool caustic){
	
	double dist = 0;
	bool intersection_found=false;
	int steps=0;
	int max_steps = 5000;
	//double P_prev[3]={P[0],P[1],P[2]};
	while(dist<MAX_DIST && intersection_found==false && steps<max_steps){
		//if(ii==288 && jj>=876){
		//	cout<<"steps "<<steps<<endl;
		//}
		double P2[3]={0,0,0};
		double V2[3]={0,0,0};
		bool step_taken = false;
		bool in_universe = true;
		int max_tries = 20;
		int try_count = 0;
		bool force_step = false;
		while(step_taken==false){
			//if(ii==288 && jj>=876){
			//	cout<<"try_count "<<try_count<<endl;
			//	cout<<"step_taken "<<step_taken<<endl;
			//	cout<<"force_step "<<force_step<<endl;
			//}
			in_universe = this->space_time->take_step(P,V,P2,V2,E,H,step_taken,force_step);
			try_count++;
			if(try_count>max_tries){
				force_step = true;
			}
		}
		steps++;
		if(steps>max_steps ){//} || this->space_time->WindingNumberTooHigh()){
			return;
		}
		double rat = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
		V[0]=V[0]/rat;
		V[1]=V[1]/rat;
		V[2]=V[2]/rat;
		
		if(!in_universe){
			return;
		}
		double step_size = sqrt((P2[0]-P[0])*(P2[0]-P[0])+(P2[1]-P[1])*(P2[1]-P[1])+(P2[2]-P[2])*(P2[2]-P[2]));
		P2[0]=P[0]+step_size*V[0];
		P2[1]=P[1]+step_size*V[1];
		P2[2]=P[2]+step_size*V[2];
		//double step_size2 = sqrt((P2[0]-P[0])*(P2[0]-P_prev[0])+(P2[1]-P_prev[1])*(P2[1]-P_prev[1])+(P2[2]-P_prev[2])*(P2[2]-P_prev[2]));
	
		double t=-1;
		double t_i;
		double x[4],y[4],n[3]; //forth coord of x & y
		//records which sub-body of a compound shape is hit
		
		//cout<<" inside raytracer "<<endl;

		n[0]=0;n[1]=0;n[2]=0;

		double Pn[3],Vn[3];

		Shape* candidate = shapes[0];
		bool winner = false;
		int winner_i = -1;
		double winner_t = -1;
	  
		double epsilon = 1e-3;
	  
		int candidate_i = -1;
		//double y_prev[4];
	  
		for(int i=0; i<shapes.size();i++){
			if(i!=last_index){ // we ignore the object we just hit.
				if(!shapes[i]->reflection_only || depth>0){
					if(invis_R==0){t_i=shapes[i]->get_Intersection(P,V,x);}
					else{t_i=shapes[i]->get_Intersection2(P,V,invis_R,x);}

					if((t_i<t && t_i>-epsilon) || (t_i>-epsilon && t==-1)){
						t=t_i;
						candidate=shapes[i];
						candidate_i = i;
						y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=x[3];
					}
				}
			}
		}
		
		//we hit something but it might be reflective or something
		if( ( t>-epsilon && t<=step_size*(1+epsilon) ) ){
			intersection_found=true;
			if(candidate->reflectivity==0 && (!caustic || depth>0)){
				Photon photon(y[0],y[1],y[2],V[0],V[1],V[2],photon_energy);
				photons[candidate_i].push_back(photon);
				return;
			}
			if(candidate->reflectivity!=0 && depth<MAX_DEPTH){
				
				// we chose a random number to decide whether to reflect or store.
				double rand_num = distrib(gen);
				//if(rand_num<=candidate->reflectivity){
				if(!caustic){
					// absorb.
					Photon photon(y[0],y[1],y[2],V[0],V[1],V[2],photon_energy);
					photons[candidate_i].push_back(photon);
					return;
				}
				else{
					double len;
					//candidate->get_reflection(y,V);
					
					//use bump map on reflection.
					candidate->get_Normal(y,n);
					//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
					//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
					if(candidate->bump_ID>-1){
						double dn[3]={0,0,0};
						//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
						//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
						//#pragma omp critical (water_reflection)
						{
						candidate->get_Bump_Perturbation(y,this->textures,dn);
						}
						n[0]=n[0]-1.0*dn[0];
						n[1]=n[1]-1.0*dn[1];
						n[2]=n[2]+1.0*dn[2];
						len = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
						n[0]/=len;
						n[1]/=len;
						n[2]/=len;
						//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
					}
					if(n[2]<0.98){
						//len = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
						//V[0]/=len;V[1]/=len;V[2]/=len;
						// rescale V
						//V[0]*=rat;
						//V[1]*=rat;
						//V[2]*=rat;
						double vdotn = V[0]*n[0]+V[1]*n[1]+V[2]*n[2];
						//cout<<"V before = "<<V[0]<<" "<<V[1]<<" "<<V[2]<<endl;
						V[0]=V[0]-2*vdotn*n[0];
						V[1]=V[1]-2*vdotn*n[1];
						V[2]=V[2]-2*vdotn*n[2];
						//cout<<"V after = "<<V[0]<<" "<<V[1]<<" "<<V[2]<<endl;
						
						double HH = this->space_time->get_h0(y);
						this->space_time->scale_velocity(V,y);
						double EE = this->space_time->compute_energy(y);
						
						this->non_linear_photon_emission(y,V,photons,photon_energy,depth+1,HH,EE,candidate_i,caustic);
						return;
					}
					return;
				}
				
			}
		
		}
		else{
			// this code is reached if a ray has not intersected
			// anything on this step.
			dist = dist+step_size;
			// reset last_index to -1
			last_index=-1;
			
			if(dist<MAX_DIST){
				P[0]=P2[0];
				P[1]=P2[1];
				P[2]=P2[2];
				V[0]=V2[0];
				V[1]=V2[1];
				V[2]=V2[2];
			}
			else{
				return;
			}
		}
	}
	return;
	
	
}


double World::get_far_time(double P[], double V[],double H,double E,double &phi_far){

	double dist = 0;
	bool intersection_found=false;
	int steps=0;
	int max_steps = 50000;
	//double P_prev[3]={P[0],P[1],P[2]};
	while(dist<MAX_DIST && intersection_found==false && steps<max_steps){
		//if(ii==288 && jj>=876){
		//	cout<<"steps "<<steps<<endl;
		//}
		double P2[4]={1,1,1,1};
		double V2[3]={0,0,0};
		bool step_taken = false;
		bool in_universe = true;
		int max_tries = 200;
		int try_count = 0;
		bool force_step = false;
		double E_backup=E;
		while(step_taken==false){
			//if(ii==288 && jj>=876){
			//	cout<<"try_count "<<try_count<<endl;
			//cout<<"step_taken "<<step_taken<<endl;
			//	cout<<"force_step "<<force_step<<endl;
			//}
			in_universe = this->space_time->take_step(P,V,P2,V2,E,H,step_taken,force_step);
			try_count++;
			if(try_count>max_tries){
				//cout<<"forcing step"<<endl;
				force_step = true;
			}
			//cout<<"P[0]= "<<P[0]<<" P[1]= "<<P[1]<<" P[2]= "<<P[2]<<endl;
			//if(P2[3]<0){
			//	cout<<"P2[0]= "<<P2[0]<<" P2[1]= "<<P2[1]<<" P2[2]= "<<P2[2]<<" P2[3]= "<<P2[3]<<endl;
			//}
			//cout<<"V[0]= "<<V[0]<<" V[1]= "<<V[1]<<" V[2]= "<<V[2]<<endl;
			//cout<<"step taken : "<<step_taken<<endl;
			//cout<<"try_count : "<<try_count<<"max_tries : "<<max_tries<<endl;
			//cout<<"H= "<<H<<endl;
		}
		steps++;
		if(steps>max_steps ){//} || this->space_time->WindingNumberTooHigh()){
			cout<<"max_steps_hit"<<endl;
			//file.close();
			return 0;
		}
		
		if(P2[3]<0){
			double r_crossing=-P2[3]; // just return entree time into the dust cloud.
			cout<<"r_crossing = "<<r_crossing<<endl;
			return this->space_time->get_redshift(r_crossing);
		}



		double R_big = 100;
		
		double exp1 = P[1]*P[1]-18*P[1]*cos(P[2]);
		double exp2 = P2[1]*P2[1]-18*P2[1]*cos(P2[2]);

		if(P[3]>0 && P2[3]>0 && exp1<=R_big*R_big && exp2>R_big*R_big){
		
			//double s = (R_big*R_big-exp1)/(exp2-exp1);
			//double t_far = P[0]*(1-s)+P2[0]*s;
			//return t_far;
			return 1;
		}
		else{
			// this code is reached if a ray has not yet exited u=1e-3
			dist = dist+0;
			// reset last_index to -1
			
			if(dist<MAX_DIST){
				//P_prev[0]=P[0];
				//P_prev[1]=P[1];
				//P_prev[2]=P[2];
				P[0]=P2[0];
				P[1]=P2[1];
				P[2]=P2[2];
				P[3]=P2[3];
				V[0]=V2[0];
				V[1]=V2[1];
				V[2]=V2[2];
			}
			else{
				return 0;
			}
		}
	}
	return 0;
}

//P and V are as in L(t)=P+tV, C is the color array returned
//ii and jj are just the image coords of the pixel - needed sometimes.
void World::non_linear_ray_trace(int ii,int jj,double P[], double V[], int C[],int depth,double H,double E,double x_ij[],int last_index,double o[], double e1[],double e2[],double e3[],double psi){

	//std::ostringstream filename;
	//filename << "Trajectories/test_light_trajectory_"<<ii<<"_"<<jj<<".txt";
	
	// Open the file and write to it
    //std::ofstream file(filename.str());
	
	// last_id is the index in the object array of the object that was last hit, or -1 if there is no such object.  
	// last_id resets to -1 after 1 step.  The point is to prevent roundoff errors where you hit an object, bounce off of it of
	// get refracted or otherwise transmitted through, only to immediately hit it again due to rounding errors.
	
	//std::ofstream output_file("test_trajectory.txt");
	
	double dist = 0;
	bool intersection_found=false;
	int steps=0;
	int max_steps = 25000;
	double max_step_size=0;
	//double P_prev[3]={P[0],P[1],P[2]};
	while(dist<MAX_DIST && intersection_found==false && steps<max_steps){
		//if(ii==288 && jj>=876){
		//	cout<<"steps "<<steps<<endl;
		//}
		double P2[4]={1,1,1,1};
		double V2[3]={0,0,0};
		bool step_taken = false;
		bool in_universe = true;
		int max_tries = 200;
		int try_count = 0;
		bool force_step = false;
		double E_backup=E;
		while(step_taken==false){
			//if(ii==288 && jj>=876){
			//	cout<<"try_count "<<try_count<<endl;
			//cout<<"step_taken "<<step_taken<<endl;
			//	cout<<"force_step "<<force_step<<endl;
			//}
			in_universe = this->space_time->take_step(P,V,P2,V2,E,H,step_taken,force_step);
			try_count++;
			if(try_count>max_tries){
				//cout<<"forcing step"<<endl;
				force_step = true;
			}
			/**cout<<"P[0]= "<<P[0]<<" P[1]= "<<P[1]<<" P[2]= "<<P[2]<<endl;
			cout<<"P2[0]= "<<P[0]<<" P2[1]= "<<P[1]<<" P2[2]= "<<P[2]<<endl;
			cout<<"V[0]= "<<V[0]<<" V[1]= "<<V[1]<<" V[2]= "<<V[2]<<endl;
			cout<<"step taken : "<<step_taken<<endl;
			cout<<"try_count : "<<try_count<<"max_tries : "<<max_tries<<endl;
			cout<<"H= "<<H<<endl;**/
		}
		steps++;
		if(steps>max_steps ){//} || this->space_time->WindingNumberTooHigh()){
			cout<<"max_steps_hit"<<endl;
			C[0]=255;
			C[1]=0;
			C[2]=0;
			//file.close();
			return;
		}
		
		//check if we have redshifted to black
		double redshift = fabs(P2[3]);

		//cout<<"P[3] P2[3] = "<<P[3]<<" "<<P2[3]<<endl;

		//cout<<"redshift = "<<redshift<<endl;

		if(redshift<=1.0/255){
			C[0]=0;
			C[1]=0;
			C[2]=0;
			return;
		}
		
		
		//if(P[3]<0 && V[0]!=-1){
		//	cout<<"P[3], V[0], depth = "<<P[3]<<" "<<V[0]<<" "<<depth<<endl;
		//}
		//cout<<"P = "<<P[0]<<" "<<P[1]<<" "<<P[2]<<" "<<P[3]<<endl;
		//cout<<"P2 = "<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<" "<<P2[3]<<endl;
		//output_file << P[0] << " " << P[1] << " " << P[2] << " " << P[3] << "\n"; // Save the current point (x, y) to the file
		
		
		double P1_cartesian[3], P2_cartesian[3], V_cartesian[3];
		this->space_time->local_to_cartesian(o,e1,e2,e3,P,P1_cartesian,psi);
		this->space_time->local_to_cartesian(o,e1,e2,e3,P2,P2_cartesian,psi);
		
		//if(ii==0 && jj==0){
		//file << P[0] << "\t" << P[1] << "\t" << P[2] << "\t"  << P2[0] << "\t" << P2[1] << "\t" << P2[2] << "\t" << P1_cartesian[0] << "\t" << P1_cartesian[1] << "\t" << P1_cartesian[2] << "\n";
		//}

		V_cartesian[0]=P2_cartesian[0]-P1_cartesian[0];
		V_cartesian[1]=P2_cartesian[1]-P1_cartesian[1];
		V_cartesian[2]=P2_cartesian[2]-P1_cartesian[2];
		//double rat = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
		//V[0]=V[0]/rat;
		//V[1]=V[1]/rat;
		//V[2]=V[2]/rat;
		double rat = sqrt(V_cartesian[0]*V_cartesian[0]+V_cartesian[1]*V_cartesian[1]+V_cartesian[2]*V_cartesian[2]);
		V_cartesian[0]=V_cartesian[0]/rat;
		V_cartesian[1]=V_cartesian[1]/rat;
		V_cartesian[2]=V_cartesian[2]/rat;
		
		bool jump = false;
		bool hit_flat_spacetime = false;
		//if(E<0 || E_backup<0){
		//	cout<<"E = "<<E<<"E_backup = "<<E_backup<<endl;
		//}
		if(E_backup>0 && E<0 && this->space_time->is_wormhole){
			cout<<"we are jumping"<<endl;
			jump=true;
		}
		if((P[0]<=1e-3 || P2[0]<=1e-3) && this->space_time->is_collapse){
			hit_flat_spacetime = true;
			steps = max_steps; // exit loop immediately after a flat spacetime raytrace.
		}
		//if(jump){
		//	cout<<"jump = "<<jump;
		//}
		/**double VV[3]={P2[0]-P_prev[0],P2[1]-P_prev[1],P2[1]-P_prev[1]};
		double rat2 = sqrt(VV[0]*VV[0]+VV[1]*VV[1]+VV[2]*VV[2]);
		VV[0]=VV[0]/rat2;
		VV[1]=VV[1]/rat2;
		VV[2]=VV[2]/rat2;**/
		if(!in_universe){
			//we've fallen in the black hole
			cout<<"exited the universe"<<endl;
			C[0]=0;
			C[1]=255;
			C[2]=0;
			//file.close();
			return;
		}
		double step_size=0;
		//if(jump){C[0]=max(int(-V[0]*255),0);C[1]=max(int(-V[1]*255),0);C[2]=max(int(-V[2]*255),0);return;}
		if(!jump){
			//cout<<"we are not jumping"<<endl;
			//if(!hit_flat_spacetime){
				//step_size = sqrt((P2[0]-P[0])*(P2[0]-P[0])+(P2[1]-P[1])*(P2[1]-P[1])+(P2[2]-P[2])*(P2[2]-P[2]));
				step_size = sqrt((P2_cartesian[0]-P1_cartesian[0])*(P2_cartesian[0]-P1_cartesian[0])+(P2_cartesian[1]-P1_cartesian[1])*(P2_cartesian[1]-P1_cartesian[1])+(P2_cartesian[2]-P1_cartesian[2])*(P2_cartesian[2]-P1_cartesian[2]));
				//cout<<"step_size = "<<step_size<<endl;
				//P2[0]=P[0]+step_size*V[0];
				//P2[1]=P[1]+step_size*V[1];
				//P2[2]=P[2]+step_size*V[2];
				P2_cartesian[0]=P1_cartesian[0]+step_size*V_cartesian[0];
				P2_cartesian[1]=P1_cartesian[1]+step_size*V_cartesian[1];
				P2_cartesian[2]=P1_cartesian[2]+step_size*V_cartesian[2];
				if(step_size>max_step_size){max_step_size=step_size;}
			//}
			//cout<<"P1_cartesian = "<<P1_cartesian[0]<<" "<<P1_cartesian[1]<<" "<<P1_cartesian[2]<<endl;
			//cout<<"P1_cartesian = "<<P1_cartesian[0]<<" "<<P2_cartesian[1]<<" "<<P2_cartesian[2]<<endl;
			//cout<<"V_cartesian = "<<V_cartesian[0]<<" "<<V_cartesian[1]<<" "<<V_cartesian[2]<<endl;
			//double step_size2 = sqrt((P2[0]-P[0])*(P2[0]-P_prev[0])+(P2[1]-P_prev[1])*(P2[1]-P_prev[1])+(P2[2]-P_prev[2])*(P2[2]-P_prev[2]));
		
			double t=-1;
			double t_i;
			double x[4],y[4],n[3]; //forth coord of x & y
			//records which sub-body of a compound shape is hit
			
			//cout<<" inside raytracer "<<endl;
			
			//double v_backup[3] = {V[0],V[1],V[2]};
			double v_backup[3] = {V_cartesian[0],V_cartesian[1],V_cartesian[2]};

			n[0]=0;n[1]=0;n[2]=0;

			double Pn[3],Vn[3];

			Shape* candidate = shapes[0];
			bool winner = false;
			int winner_i = -1;
			double winner_t = -1;
		  
			double epsilon = 1e-3;
		  
			int candidate_i = -1;
			//double y_prev[4];
		  
			for(int i=0; i<shapes.size();i++){
				if(i!=last_index){ // we ignore the object we just hit.
					if(!shapes[i]->reflection_only || depth>0){
						//if(invis_R==0){t_i=shapes[i]->get_Intersection(P,V,x);}
						if(invis_R==0){t_i=shapes[i]->get_Intersection(P1_cartesian,V_cartesian,x);}
						else{t_i=shapes[i]->get_Intersection2(P,V,invis_R,x);}
						//cout<<"rat i="<<i<<" true ="<<true<<" is wireframe = "<<shapes[i]->is_wireframe<<" wireframe condition = "<<shapes[i]->wireframe_condition(x)<<endl;
						if(shapes[i]->is_wireframe && !(shapes[i]->wireframe_condition(x))){t_i=-1;}
						if(shapes[i]->is_wireframe && shapes[i]->wireframe_condition(x) && t_i>0.01){
							winner = true;
							if(winner_t==-1){winner_t = t_i; winner_i = i;}
							else if(winner_t!=-1 && t_i<winner_t){winner_i=i;winner_t = t_i;}
						}

						if((t_i<t && t_i>-epsilon) || (t_i>-epsilon && t==-1)){
							t=t_i;
							candidate=shapes[i];
							candidate_i = i;
							y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=x[3];
						}
					}
				}
			}
			//cout<<"t = "<<t<<"shapes.size() = "<<shapes.size()<<endl;
			
			
			/**Shape* candidate_prev = shapes[0];
			int candidate_prev_i = -1;
			double t_prev = -1;
			
			for(int i=0; i<shapes.size();i++){
				if(i!=last_index){ // we ignore the object we just hit.
					if(!shapes[i]->reflection_only || depth>-epsilon){
						if(invis_R==0){t_i=shapes[i]->get_Intersection(P_prev,VV,x);}
						else{t_i=shapes[i]->get_Intersection2(P_prev,VV,invis_R,x);}
						//cout<<"rat i="<<i<<" true ="<<true<<" is wireframe = "<<shapes[i]->is_wireframe<<" wireframe condition = "<<shapes[i]->wireframe_condition(x)<<endl;
						if(shapes[i]->is_wireframe && !(shapes[i]->wireframe_condition(x))){t_i=-1;}
						if(shapes[i]->is_wireframe && shapes[i]->wireframe_condition(x) && t_i>0.01){
							winner = true;
							if(winner_t==-1){winner_t = t_i; winner_i = i;}
							else if(winner_t!=-1 && t_i<winner_t){winner_i=i;winner_t = t_i;}
						}

						if((t_i<t_prev && t_i>-epsilon) || (t_i>-epsilon && t_prev==-1)){
							t_prev=t_i;
							candidate_prev=shapes[i];
							candidate_prev_i = i;
							y_prev[0]=x[0];y_prev[1]=x[1];y_prev[2]=x[2];y_prev[3]=x[3];
						}
					}
				}
			}**/
			if(t!=-1 && (t<-epsilon || t>step_size*(1+epsilon)) ){
				//cout<<" t = "<<t<<" step_size = "<<step_size<<"epsilon = "<<epsilon<<endl;
			}				
			// we hit something within the step, and reflections etc are impossible
			if(winner && abs(t-winner_t)<0.001 && t<=step_size*(1+epsilon)){
				x_ij[0]=y[0];x_ij[1]=y[1];x_ij[2]=y[2];
				shapes[winner_i]->get_Color(y,this->textures,C);
				intersection_found=true;
				//file.close();
				return;
			}
			//we hit something but it might be reflective or something
			if( ( t>-epsilon && t<=step_size*(1+epsilon) ) || ( t>-epsilon && hit_flat_spacetime ) ){// || ( t_prev>-epsilon && t_prev<=step_size2*(1+epsilon) ) ){
				/**if(candidate_i!=candidate_prev_i){
					candidate_i = candidate_prev_i;
					candidate = candidate_prev;
					y[0]=y_prev[0];
					y[1]=y_prev[1];
					y[2]=y_prev[2];
					y[3]=y_prev[3];
				}**/
				//double t_far=(P[1]+t*V[1]+P[0]+t*V[0])/2.0;
				double t_far=(P[1]+P[0])/2.0;
				//cout<<"intersection found, candidate_i = "<<candidate_i<<endl;
				intersection_found=true;
				x_ij[0]=y[0];x_ij[1]=y[1];x_ij[2]=y[2];
				//cout<<"intersection found!"<<endl;
				candidate->get_Color(y,this->textures,C);
				//cout<<"C = "<<C[0]<<" "<<C[1]<<" "<<C[2]<<endl;
				//if(P[3]<0 or P2[3]<0){
				//	C[2]=255;
				//}
				int C2[3]={0,0,0};
				//cout<<"y = "<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl;
				candidate->get_Shading(y,this->photon_maps,candidate_i,C2);

				//cout<<"C2 = "<<C2[0]<<" "<<C2[1]<<" "<<C2[2]<<endl;
				//int C2[3]={C[0],C[1],C[2]};
				//C[0]=min(C[0]/10+int(floor(2*C[0]*C2[0]/255.0)),255);
				//C[1]=min(C[1]/10+int(floor(2*C[1]*C2[1]/255.0)),255);
				//C[2]=min(C[2]/10+int(floor(2*C[2]*C2[2]/255.0)),255);
				//cout<<" C = "<<C[0]<<" "<<C[1]<<" "<<C[2]<<endl;
				
				
				C[0]=int(floor(redshift*C[0]*C2[0]/255.0));
				C[1]=int(floor(redshift*C[1]*C2[1]/255.0));
				C[2]=int(floor(redshift*C[2]*C2[2]/255.0));
				
				//return;
				
				//cout<<" C = "<<C[0]<<" "<<C[1]<<" "<<C[2]<<endl;
				//cout<<" candidate_id = "<<candidate_i<<endl;
				//if (hit_flat_spacetime){
				//	C[2]=0;
				//}
				//C[3]=int(floor(t_far*10000));
				//cout<<"saving t_far="<<t_far<<" "<<"C[3]="<<C[3]<<endl;
				
				//if(hit_flat_spacetime){
				//	C[0]=0;
				//}
				
				
				//cout<<"C[0] = "<<C[0]<<" C[1] = "<<C[1]<<" C[2] = "<<C[2]<<endl;
				//C[0]=C2[0];
				//C[1]=C2[1];
				//C[2]=C2[2];
				candidate->get_Normal(y,n);
				if(candidate->bump_ID>-1 && !(candidate->self_illuminating)){
					//#pragma omp critical (bump_map)
					{
					//cout<<"hit main fractal"<<endl;
					double dn[3]={0,0,0};
					//if(jj==255){
					//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
					candidate->get_Bump_Perturbation(y,this->textures,dn);
					//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
					n[0]=n[0]+5*dn[0];
					n[1]=n[1]+5*dn[1];
					n[2]=n[2]+5*dn[2];
					double len = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
					n[0]/=len;
					n[1]/=len;
					n[2]/=len;
					//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
					//}
					}
				}
				//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
				/**if(ii>4 && jj>4){
					//double x1[3]={X_data[ii-4][jj][0],X_data[ii-4][jj][1],X_data[ii-4][jj][2]};
					//double x2[3]={X_data[ii][jj-4][0],X_data[ii][jj-4][1],X_data[ii][jj-4][2]};
					double x1[3]={x_im[0]-x_ij[0],x_im[1]-x_ij[1],x_im[2]-x_ij[2]};
					double x2[3]={x_jm[0]-x_ij[0],x_jm[1]-x_ij[1],x_jm[2]-x_ij[2]};
					n[0]=x1[1]*x2[2]-x1[2]*x2[1];
					n[1]=x1[2]*x2[0]-x1[0]*x2[2];
					n[2]=x1[0]*x2[1]-x1[1]*x2[0];
					double len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
					n[0]=n[0]/len;
					n[1]=n[1]/len;
					n[2]=n[2]/len;
					double n_dot_v = n[0]*v_backup[0]+n[1]*v_backup[1]+n[2]*v_backup[2];
					if(n_dot_v>=0){
						n[0]=-n[0];
						n[1]=-n[1];
						n[2]=-n[2];
					}
				}**/
				// temporarily disabling this.
				/**if(!(candidate->is_wireframe)){
					double reflection[3];
					if(!candidate->self_illuminating && false){
						//cout<<"shading..."<<endl;
						//cout<<"n[0] = "<<n[0]<<"n[1] = "<<n[1]<<"n[2]= "<<n[2]<<endl;
						double dot_prod = n[0]*V[0]+n[1]*V[1]+n[2]*V[2];
						reflection[0]=V[0]-2*dot_prod*n[0];
						reflection[1]=V[1]-2*dot_prod*n[1];
						reflection[2]=V[2]-2*dot_prod*n[2];
						//cout<<"reflection computed"<<endl;
						this->shade(y,n,reflection,C);
						//cout<<"shading succeeded"<<endl;
						return;
					}
					else if(candidate->reflectivity==0 && candidate->transparency==0){
						return;
					}
				}**/

				//transport_tangent(P,t,invis_R,V);


				// refraction also disabled temporarily.
				/**if(candidate->transparency!=0 || candidate->alpha_ID>=0 && depth<MAX_DEPTH){
					int Ctrans[3];
					double xx[3];
					if( candidate->get_refracted_ray(y,V) ){depth=depth+1;}
					y[0]=y[0]+0.01*V[0];
					y[1]=y[1]+0.01*V[1];
					y[2]=y[2]+0.01*V[2];
					///double HH = this->space_time->get_h0(y);
					//this->space_time->scale_velocity(V,y);
					//double EE = this->space_time->compute_energy(y);
					V[0]=V[0]*rat;
					V[1]=V[1]*rat;
					V[2]=V[2]*rat;
					double y_copy[3]={y[0],y[1],y[2]};
					double V_copy[3]={V[0],V[1],V[2]};
					this->non_linear_ray_trace(ii,jj,y_copy,V_copy,Ctrans,depth+1,H,E,xx,candidate_i);
					double tr;
					if(candidate->alpha_ID>=0){
						tr = candidate->get_Alpha(y, textures);
					}
					else{
						tr=candidate->transparency;
					}
					C[0]=int((1-tr)*C[0]+tr*Ctrans[0]);
					C[1]=int((1-tr)*C[1]+tr*Ctrans[1]);
					C[2]=int((1-tr)*C[2]+tr*Ctrans[2]);
					//C[0]=int(tr*255);
					//C[1]=int(tr*255);
					//C[2]=int(tr*255);
					return;
				}**/
				// reflections we want.  Not disabling.
				if(candidate->reflectivity!=0 && depth<MAX_DEPTH){
					int Cnew[3];
					double xx[3];
					double len;
					//candidate->get_reflection(y,V);
					
					//use bump map on reflection.
					candidate->get_Normal(y,n);
					//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
					//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
					if(candidate->bump_ID>-1){
						double dn[3]={0,0,0};
						//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
						//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
						//#pragma omp critical (water_reflection)
						{
						candidate->get_Bump_Perturbation(y,this->textures,dn);
						}
						n[0]=n[0]-0.1*dn[0];
						n[1]=n[1]-0.1*dn[1];
						n[2]=n[2]+0.1*dn[2];
						len = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
						n[0]/=len;
						n[1]/=len;
						n[2]/=len;
						//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
					}
					//len = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
					//V[0]/=len;V[1]/=len;V[2]/=len;
					// rescale V
					//V[0]*=rat;
					//V[1]*=rat;
					//V[2]*=rat;
					
					//V_cartesian[0]*=rat;
					//V_cartesian[1]*=rat;
					//V_cartesian[2]*=rat;
					
					//double vdotn = V[0]*n[0]+V[1]*n[1]+V[2]*n[2];
					//V[0]=V[0]-2*vdotn*n[0];
					//V[1]=V[1]-2*vdotn*n[1];
					//V[2]=V[2]-2*vdotn*n[2];
					
					// disable reflections for now.
					
					// need to edit a few things.
					
					//double P_intersection[]={P1_cartesian[0]+t*V_cartesian[0],P1_cartesian[1]+t*V_cartesian[1],P1_cartesian[2]+t*V_cartesian[2]};
					double P_intersection[]={y[0],y[1],y[2]};


					//cout<<" V_cartesian before = "<<V_cartesian[0]<<" "<<V_cartesian[1]<<" "<<V_cartesian[2]<<endl;
					//cout<<" n = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;

					double vdotn = V_cartesian[0]*n[0]+V_cartesian[1]*n[1]+V_cartesian[2]*n[2];
					V_cartesian[0]=V_cartesian[0]-2*vdotn*n[0];
					V_cartesian[1]=V_cartesian[1]-2*vdotn*n[1];
					V_cartesian[2]=V_cartesian[2]-2*vdotn*n[2];
					
					//P_intersection[0]=P_intersection[0]+0.001*n[0];
					//P_intersection[1]=P_intersection[1]+0.001*n[1];
					//P_intersection[2]=P_intersection[2]+0.001*n[2];
					
					//cout<<" P_intersection = "<<P_intersection[0]<<" "<<P_intersection[1]<<" "<<P_intersection[2]<<endl;
					//cout<<" V_cartesian = "<<V_cartesian[0]<<" "<<V_cartesian[1]<<" "<<V_cartesian[2]<<endl;

					// we need to build a new 2D plane containing the intersection point and the origin of our system.
					double u[]={0,0,1};
					double e1_new[]={0,0,0};
					double e2_new[]={0,0,0};
					double e3_new[]={0,0,0};
					
					//return; // exit before a reflection.

					this->space_time->recompute_orthonormal_frame2(o,P_intersection,e1_new,e2_new,e3_new);
					double r_new = sqrt((P_intersection[0]-o[0])*(P_intersection[0]-o[0])
					                    +(P_intersection[1]-o[1])*(P_intersection[1]-o[1])
										+(P_intersection[2]-o[2])*(P_intersection[2]-o[2]));
					
					
					/**if(P[3]>0){
						r_new = P[1]+t*V[1];
					}
					else{
						double tau = P[0]+t*V[0]+0;
						double X_new = asin(r_new/this->space_time->get_a(tau));
						//X_new = P[1]*(1-t)+P2[1]*t;
						double r_new = X_new;
					}**/
					
					double s = t/step_size;
					
					r_new = P[1]*(1-s)+P2[1]*s;
					//r_new = P[1]+(t/step_size)*V[1];
					//double t_new = P[0]+(t/step_size)*V[0]+0;
					double t_new = P[0]*(1-s)+P2[0]*s+0.0;
					//cout<<"P[3], P2[3] = "<<P[3]<<" "<<P2[3]<<endl;
					//cout<<"P[0], t, V[0], t_new = "<<P[0]<<" "<<t<<" "<<V[0]<<" "<<t_new<<endl;
					//if(P[3]*P2[3]<0){return;}
					/**double t_new;
					if(P[3]>0){
						t_new = P[0]+t*V[0];
					}
					else{
						tau_new = this->spacetime->tau_of_eta(P[0])+t*
					}**/
					double p_local_new[4]={t_new,r_new,0,P[3]};
					double v_obs_local[]={0,0,0}; // doesn't do anything.
					double v_local_new[]={0,0,0}; // to be computed.
					double HH = this->space_time->get_h0(p_local_new);
					double EE = this->space_time->compute_energy(p_local_new);
					double psi_new;
					this->space_time->get_null_vector_in_direction(e1_new,e2_new,e3_new,p_local_new,v_obs_local,V_cartesian,v_local_new,psi_new);
					//this->non_linear_ray_trace(ii,jj,y,V,Cnew,depth+1,HH,EE,xx,candidate_i);
					//cout<<"doing reflection"<<endl;
					//cout<<"candidate->texture_ID = "<<candidate->texture_ID<<endl;
					if(depth>1){
						cout<<"depth = "<<depth<<endl;
					}
					this->non_linear_ray_trace(ii,jj,p_local_new,v_local_new,Cnew,depth+1,HH,EE,xx,candidate_i,o,e1_new,e2_new,e3_new,psi_new);
					
					//cout<<"reflection done, C, Cnew = "<<C[0]<<" "<<C[1]<<" "<<C[2]<<" "<<Cnew[0]<<" "<<Cnew[1]<<" "<<Cnew[2]<<endl;
					//double y_local[3];
					//double V_local[3];
					//this->space_time->cartesian_to_local(y,y_local);
					//this->space_time->cartesian_to_local_vector(V_local,V_cartesian);
					
					//double HH = this->space_time->get_h0(y);
					//this->space_time->scale_velocity(V,y);
					//double EE = this->space_time->compute_energy(y);
					
					//this->non_linear_ray_trace(ii,jj,y,V,Cnew,depth+1,HH,EE,xx,candidate_i);
					//double HH = this->space_time->get_h0(y_local);
					//this->space_time->scale_velocity(V_local,y_local);
					//double EE = this->space_time->compute_energy(y_local);
					
					//this->non_linear_ray_trace(ii,jj,y_local,V_local,Cnew,depth+1,HH,EE,xx,candidate_i);
					double rf=candidate->reflectivity;
					C[0]=int((1-rf)*C[0]+rf*redshift*Cnew[0]);
					C[1]=int((1-rf)*C[1]+rf*redshift*Cnew[1]);
					C[2]=int((1-rf)*C[2]+rf*redshift*Cnew[2]);
					//cout<<"C_combined = "<<C[0]<<" "<<C[1]<<" "<<C[2]<<endl;
					//file.close();
					return;
				}
				return;
			}
			else{
				// this code is reached if a ray has not intersected
				// anything on this step.
				dist = dist+step_size;
				// reset last_index to -1
				last_index=-1;
				
				if(dist<MAX_DIST){
					//P_prev[0]=P[0];
					//P_prev[1]=P[1];
					//P_prev[2]=P[2];
					P[0]=P2[0];
					P[1]=P2[1];
					P[2]=P2[2];
					P[3]=P2[3];
					V[0]=V2[0];
					V[1]=V2[1];
					V[2]=V2[2];
				}
				else{
					//cout<<"nothing hit, dist ="<<dist<<" steps = "<<steps<<" max_step_size = "<<max_step_size<<"max_step_size*steps = "<<max_step_size*steps<<endl;
					if(env_map_ID==-1){
						C[0]=0;C[1]=0;C[2]=0;
						//file.close();
						return;
						//if(depth==0){C[0]=170;C[1]=170;C[2]=170;}
					}
					else{
						//cout<<"hit star texture"<<endl;
						//# pragma omp critical
						{
							//cout<<"reading sphere map"<<endl;
							this->read_sphere_env_map(V,C);
							//file.close();
							return;
							//cout<<"sphere map read"<<endl;
						}
					}
				}
			}
		}
		else{
			// this code is reached if a ray has not intersected
			// anything on this step.
			dist = dist+step_size;
			// reset last_index to -1
			last_index=-1;
			
			if(dist<MAX_DIST){
				//P_prev[0]=P[0];
				//P_prev[1]=P[1];
				//P_prev[2]=P[2];
				P[0]=P2[0];
				P[1]=P2[1];
				P[2]=P2[2];
				P[3]=P2[3];
				V[0]=V2[0];
				V[1]=V2[1];
				V[2]=V2[2];
			}
			else{
				//cout<<"nothing hit"<<endl;
				if(env_map_ID==-1){
					cout<<"reached maximum distance.  Distance = "<<dist<<endl;
					C[0]=0;C[1]=0;C[2]=0;
					//file.close();
					return;
					//if(depth==0){C[0]=170;C[1]=170;C[2]=170;}
				}
				else{
					//cout<<"hit star texture"<<endl;
					//# pragma omp critical
					{
						//cout<<"reading sphere map"<<endl;
						//file.close();
						this->read_sphere_env_map(V,C);
						return;
						//cout<<"sphere map read"<<endl;
					}
				}
			}
		}
	}
	//if(steps>=max_steps){cout<<"hit maximum steps without hitting anything"<<endl;}
	// default to reading the sphere map if we exit for an unknown reason.
	//cout<<"final outer"<<endl;
	//cout<<"dist ="<<dist<<" steps = "<<steps<<" max_step_size = "<<max_step_size<<"max_step_size*steps = "<<max_step_size*steps<<endl;
	C[0]=255;
	C[1]=255;
	C[2]=0;
	//file.close();
	//this->read_sphere_env_map(V,C);
	return;
}

bool World::linear_ray_hits_something(int ii,int jj,double P[], double V[]){
	double t=-1;
	double t_i;
	double x[4],y[4],n[3]; //forth coord of x & y
	//records which sub-body of a compound shape is hit
	
	//cout<<" inside raytracer "<<endl;
	
	double v_backup[3] = {V[0],V[1],V[2]};

	n[0]=0;n[1]=0;n[2]=0;

	double Pn[3],Vn[3];

	Shape* candidate = shapes[0];
	bool winner = false;
	int winner_i = -1;
	double winner_t = -1;
  
	for(int i=0; i<shapes.size();i++){
		if(!shapes[i]->reflection_only){
			if(invis_R==0){t_i=shapes[i]->get_Intersection(P,V,x);}
			else{t_i=shapes[i]->get_Intersection2(P,V,invis_R,x);}
			//cout<<"rat i="<<i<<" true ="<<true<<" is wireframe = "<<shapes[i]->is_wireframe<<" wireframe condition = "<<shapes[i]->wireframe_condition(x)<<endl;

			if((t_i<t && t_i>0.01) || (t_i>0.01 && t<0)){
				t=t_i;
				candidate=shapes[i];
				y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=x[3];
			}
		}
	}

	if(t>0.01){
		return true;
	}
	return false;
}


//P and V are as in L(t)=P+tV, C is the color array returned
//ii and jj are just the image coords of the pixel - needed sometimes.
void World::ray_trace(int ii,int jj,double P[], double V[], int C[],int depth,double x_ij[]){
	double t=-1;
	double t_i;
	double x[4],y[4],n[3]; //forth coord of x & y
	//records which sub-body of a compound shape is hit
	
	//cout<<" inside raytracer "<<endl;
	
	double v_backup[3] = {V[0],V[1],V[2]};

	n[0]=0;n[1]=0;n[2]=0;

	double Pn[3],Vn[3];

	Shape* candidate = shapes[0];
	bool winner = false;
	int winner_i = -1;
	double winner_t = -1;
  
	for(int i=0; i<shapes.size();i++){
		if(!shapes[i]->reflection_only || depth>0){
			if(invis_R==0){t_i=shapes[i]->get_Intersection(P,V,x);}
			else{t_i=shapes[i]->get_Intersection2(P,V,invis_R,x);}
			//cout<<"rat i="<<i<<" true ="<<true<<" is wireframe = "<<shapes[i]->is_wireframe<<" wireframe condition = "<<shapes[i]->wireframe_condition(x)<<endl;
			if(shapes[i]->is_wireframe && !(shapes[i]->wireframe_condition(x))){t_i=-1;}
			if(shapes[i]->is_wireframe && shapes[i]->wireframe_condition(x) && t_i>0.01){
				winner = true;
				if(winner_t==-1){winner_t = t_i; winner_i = i;}
				else if(winner_t!=-1 && t_i<winner_t){winner_i=i;winner_t = t_i;}
			}

			if((t_i<t && t_i>0.01) || (t_i>0.01 && t<0)){
				t=t_i;
				candidate=shapes[i];
				y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=x[3];
			}
		}
	}
	if(winner && abs(t-winner_t)<0.001){
		x_ij[0]=y[0];x_ij[1]=y[1];x_ij[2]=y[2];
		shapes[winner_i]->get_Color(y,this->textures,C);
		return;
	}
	//if(t>0.01 && !winner){
	if(t>0.01){
		x_ij[0]=y[0];x_ij[1]=y[1];x_ij[2]=y[2];
		//cout<<"intersection found!"<<endl;
		candidate->get_Color(y,this->textures,C);
		candidate->get_Normal(y,n);
		if(candidate->bump_ID>-1 && !(candidate->self_illuminating)){
			//#pragma omp critical (bump_map)
			{
			//cout<<"hit main fractal"<<endl;
			double dn[3]={0,0,0};
			//if(jj==255){
			//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
			candidate->get_Bump_Perturbation(y,this->textures,dn);
			//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
			n[0]=n[0]+5*dn[0];
			n[1]=n[1]+5*dn[1];
			n[2]=n[2]+5*dn[2];
			double len = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
			n[0]/=len;
			n[1]/=len;
			n[2]/=len;
			//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
			//}
			}
		}
		//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
		/**if(ii>4 && jj>4){
			//double x1[3]={X_data[ii-4][jj][0],X_data[ii-4][jj][1],X_data[ii-4][jj][2]};
			//double x2[3]={X_data[ii][jj-4][0],X_data[ii][jj-4][1],X_data[ii][jj-4][2]};
			double x1[3]={x_im[0]-x_ij[0],x_im[1]-x_ij[1],x_im[2]-x_ij[2]};
			double x2[3]={x_jm[0]-x_ij[0],x_jm[1]-x_ij[1],x_jm[2]-x_ij[2]};
			n[0]=x1[1]*x2[2]-x1[2]*x2[1];
			n[1]=x1[2]*x2[0]-x1[0]*x2[2];
			n[2]=x1[0]*x2[1]-x1[1]*x2[0];
			double len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
			n[0]=n[0]/len;
			n[1]=n[1]/len;
			n[2]=n[2]/len;
			double n_dot_v = n[0]*v_backup[0]+n[1]*v_backup[1]+n[2]*v_backup[2];
			if(n_dot_v>=0){
				n[0]=-n[0];
				n[1]=-n[1];
				n[2]=-n[2];
			}
		}**/
		if(!(candidate->is_wireframe)){
			double reflection[3];
			if(!candidate->self_illuminating){
				double dot_prod = n[0]*V[0]+n[1]*V[1]+n[2]*V[2];
				reflection[0]=V[0]-2*dot_prod*n[0];
				reflection[1]=V[1]-2*dot_prod*n[1];
				reflection[2]=V[2]-2*dot_prod*n[2];
				this->shade(y,n,reflection,C);
			}
		}

		transport_tangent(P,t,invis_R,V);

		if(candidate->transparency!=0 || candidate->alpha_ID>=0){
			int Ctrans[3];
			double xx[3];
			if( candidate->get_refracted_ray(y,V) ){depth=depth+1;}
			this->ray_trace(ii,jj,y,V,Ctrans,depth,xx);
			double tr;
			if(candidate->alpha_ID>=0){
				tr = candidate->get_Alpha(y, textures);
			}
			else{
				tr=candidate->transparency;
			}
			C[0]=int((1-tr)*C[0]+tr*Ctrans[0]);
			C[1]=int((1-tr)*C[1]+tr*Ctrans[1]);
			C[2]=int((1-tr)*C[2]+tr*Ctrans[2]);
		}
		if(candidate->reflectivity!=0 && depth<MAX_DEPTH){
			int Cnew[3];
			double xx[3];
			double len;
			//candidate->get_reflection(y,V);
			
			//use bump map on reflection.
			candidate->get_Normal(y,n);
			//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
			//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
			if(candidate->bump_ID>-1){
				double dn[3]={0,0,0};
				//cout<<"dn = "<<dn[0]<<" "<<dn[1]<<" "<<dn[2]<<endl;
				//cout<<"n before = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
				//#pragma omp critical (water_reflection)
				{
				candidate->get_Bump_Perturbation(y,this->textures,dn);
				}
				n[0]=n[0]-0.1*dn[0];
				n[1]=n[1]-0.1*dn[1];
				n[2]=n[2]+0.1*dn[2];
				len = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
				n[0]/=len;
				n[1]/=len;
				n[2]/=len;
				//cout<<"n after = "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
			}
			//len = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			//V[0]/=len;V[1]/=len;V[2]/=len;
			double vdotn = V[0]*n[0]+V[1]*n[1]+V[2]*n[2];
			//cout<<"V before = "<<V[0]<<" "<<V[1]<<" "<<V[2]<<endl;
			V[0]=V[0]-2*vdotn*n[0];
			V[1]=V[1]-2*vdotn*n[1];
			V[2]=V[2]-2*vdotn*n[2];
			//cout<<"V after = "<<V[0]<<" "<<V[1]<<" "<<V[2]<<endl;
			
			
			
			this->ray_trace(ii,jj,y,V,Cnew,depth+1,xx);
			double rf=candidate->reflectivity;
			C[0]=int((1-rf)*C[0]+rf*Cnew[0]);
			C[1]=int((1-rf)*C[1]+rf*Cnew[1]);
			C[2]=int((1-rf)*C[2]+rf*Cnew[2]);
			
		}
		
	}
	else{
		// this code is reached if a ray has not intersected
		//anything.  
		//cout<<"nothing hit"<<endl;
		if(env_map_ID==-1){
			C[0]=0;C[1]=0;C[2]=0;
			//if(depth==0){C[0]=170;C[1]=170;C[2]=170;}
		}
		else{
			//cout<<"hit star texture"<<endl;
			//# pragma omp critical
			{
				//cout<<"reading sphere map"<<endl;
				this->read_sphere_env_map(V,C);
				//cout<<"sphere map read"<<endl;
			}
		}
		//else{this->read_env_map(ii,jj,C);}
	}
	return;
}


/**void World::render_hogel(Camera* cam, char* filename, int a, int N_a, int N_b, int N_I, double h_c[], double Ly, double Lz, double r, double theta_SLM){

	cout<<"rendering hogel a = "<<a<<endl;

	//initialize image
	hogel_data = new int**[N_I];
	for(int i=0;i<N_I;i++){
		hogel_data[i]=new int*[N_b];
		for(int b=0;b<N_b;b++){
			hogel_data[i][b]=new int[3];
			hogel_data[i][b][0]=0;
			hogel_data[i][b][1]=0;
			hogel_data[i][b][2]=0;
		}
	}
	
	cout<<"initialized hogel image"<<endl;
	
	double p[3];
	
	double h_x = 0.3*(1/10.0);
	double h_y = 0.266*(1/10.0);
	
	int N_rays = 4;
	
	double i_r,j_r,k_r; //these will hold small random perturbations in position index
	
	srand( time(NULL) );
	

	for(int b=0;b<N_b;b++){
		cout<<"rendering hogel b = "<<b<<endl;
		// first let's get our position on the hologram plane.
		//double dy = -Ly/2+(a-1)*Ly/(N_a-1);
		//double dz = -Lz/2+(b-1)*Lz/(N_b-1);
		for(int i=0;i<N_I;i++){
			//cout<<"rendering i = "<<i<<endl;
			
			for(int ii=0;ii<N_rays;ii++){
				for(int jj=0;jj<N_rays;jj++){
					for(int kk=0;kk<N_rays;kk++){
					
						//i_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
						//j_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
						//k_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
						i_r=((rand() % 100))/100.0; //i.e. rand(0,1.0)
						j_r=((rand() % 100))/100.0; //i.e. rand(0,1.0)
						k_r=((rand() % 100))/100.0; //i.e. rand(0,1.0)
						
						//double dy = -Ly/2+h_x/2+(a+float(kk+k_r)/N_rays-1)*h_x;
						//double dz = -Lz/2+h_y/2+(b+float(jj+j_r)/N_rays-1)*h_y;
						double dy = -Ly/2+h_x/2+(a+float(kk+k_r)/N_rays-0.5)*h_x;
						double dz = -Lz/2+h_y/2+(b+float(jj+j_r)/N_rays-0.5)*h_y;
						
						// staggering
						if(b%2==1){
							dy=dy+h_x/2;
						}
						
						p[0]=h_c[0];
						p[1]=h_c[1]+dy;
						p[2]=h_c[2]+dz;
						
						
						
						int C[3]={0,0,0};
						double x[3]={0,0,0};
						
						double dy2 = -tan(theta_SLM/2)+(i+float(ii+i_r)/N_rays-1)*2*tan(theta_SLM/2)/(N_I-1);
						double v[2]={1,dy2};
						//double A = h_c[0]*h_c[0]+h_c[1]*h_c[1];
						double A = v[0]*v[0]+v[1]*v[1];
						double B = 2*v[0]*(p[0]-h_c[0])+2*v[1]*(p[1]-h_c[1]);
						double CC = (p[0]-h_c[0])*(p[0]-h_c[0])+(p[1]-h_c[1])*(p[1]-h_c[1])-r*r;
						double desc = B*B-4*A*CC;
						if(desc<0){
							cout<<"desc = "<<desc<<endl;
						}
						double t1 = (-B+sqrt(desc))/(2*A);
						double t2 = (-B-sqrt(desc))/(2*A);
						double t = fmax(t1,t2);
						double p_c[3]={p[0]+t*v[0],p[1]+t*v[1],h_c[2]};
						double V[3]={p[0]-p_c[0],p[1]-p_c[1],p[2]-p_c[2]}; // this is the ray we should ray trace.
						
						double length = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
						V[0]=V[0]/length;
						V[1]=V[1]/length;
						V[2]=V[2]/length;
						
						
						double H = this->space_time->get_h0(p_c);
						this->space_time->scale_velocity(V,p_c);
						double E = this->space_time->compute_energy(p_c);
						int last_index = -1;
						if(linear_ray_hits_something(i,b,p_c,V)){
							this->non_linear_ray_trace(i,b,p_c,V,C,0,H,E,x,last_index);
						}
						else{
							C[0]=0;
							C[1]=0;
							C[2]=0;
						}
						hogel_data[i][b][0]+=C[0];
						hogel_data[i][b][1]+=C[1];
						hogel_data[i][b][2]+=C[2];
					}
				}
			}
			hogel_data[i][b][0]/=(N_rays*N_rays*N_rays);
			hogel_data[i][b][1]/=(N_rays*N_rays*N_rays);
			hogel_data[i][b][2]/=(N_rays*N_rays*N_rays);
		}
	}
	cout<<"done making hogel image"<<endl;
	
	save_bmp(N_I,N_b,filename,hogel_data);	
	for(int i=0;i<N_I;i++){
		for(int b=0;b<N_b;b++){
			delete[] hogel_data[i][b];
		}
		delete[] hogel_data[i];
	}
	delete[] hogel_data;


	
}


void World::render_full_parallax_hogel(Camera* cam, char* filename, int a, int b, int N_a, int N_b, int N_I, int N_J, double h_c[], double Ly, double Lz, double r, double theta_SLM_h,double theta_SLM_v){

	cout<<"rendering hogel a = "<<a<<endl;

	//initialize image
	hogel_data = new int**[N_I];
	for(int i=0;i<N_I;i++){
		hogel_data[i]=new int*[N_J];
		for(int j=0;j<N_J;j++){
			hogel_data[i][j]=new int[3];
			hogel_data[i][j][0]=0;
			hogel_data[i][j][1]=0;
			hogel_data[i][j][2]=0;
		}
	}
	
	cout<<"initialized hogel image"<<endl;
	
	double p[3];
	
	for(int j=0;j<N_J;j++){
		cout<<"rendering j = "<<j<<endl;
		// first let's get our position on the hologram plane.
		double dy = -Ly/2+(a-1)*Ly/(N_a-1);
		double dz = -Lz/2+(b-1)*Lz/(N_b-1);
		p[0]=h_c[0];
		p[1]=h_c[1]+dy;
		p[2]=h_c[2]+dz;
		int C[3]={0,0,0};
		double x[3]={0,0,0};
		double dz2 = -tan(theta_SLM_v/2)+(j-1)*2*tan(theta_SLM_v/2)/(N_J-1);
		for(int i=0;i<N_I;i++){
			//cout<<"rendering i = "<<i<<endl;
			double dy2 = -tan(theta_SLM_h/2)+(i-1)*2*tan(theta_SLM_h/2)/(N_I-1);
			double v[3]={1,dy2,dz2};
			v[0]=v[0]/sqrt(1+dy2*dy2+dz2*dz2);
			v[1]=v[1]/sqrt(1+dy2*dy2+dz2*dz2);
			v[2]=v[2]/sqrt(1+dy2*dy2+dz2*dz2);
			//double A = h_c[0]*h_c[0]+h_c[1]*h_c[1];
			double A = v[0]*v[0]+v[1]*v[1];
			double B = 2*v[0]*(p[0]-h_c[0])+2*v[1]*(p[1]-h_c[1]);
			double CC = (p[0]-h_c[0])*(p[0]-h_c[0])+(p[1]-h_c[1])*(p[1]-h_c[1])-r*r;
			double desc = B*B-4*A*CC;
			if(desc<0){
				cout<<"desc = "<<desc<<endl;
			}
			double t1 = (-B+sqrt(desc))/(2*A);
			double t2 = (-B-sqrt(desc))/(2*A);
			double t = fmax(t1,t2);
			double p_c[3]={p[0]+t*v[0],p[1]+t*v[1],p[2]+t*v[2]};
			double V[3]={p[0]-p_c[0],p[1]-p_c[1],p[2]-p_c[2]}; // this is the ray we should ray trace.
			
			double length = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0]=V[0]/length;
			V[1]=V[1]/length;
			V[2]=V[2]/length;
			
			
			double H = this->space_time->get_h0(p_c);
			this->space_time->scale_velocity(V,p_c);
			double E = this->space_time->compute_energy(p_c);
			int last_index = -1;
			if(linear_ray_hits_something(i,j,p_c,V)){
				this->non_linear_ray_trace(i,j,p_c,V,C,0,H,E,x,last_index);
			}
			else{
				C[0]=0;
				C[1]=0;
				C[2]=0;
			}
			hogel_data[i][j][0]=C[0];
			hogel_data[i][j][1]=C[1];
			hogel_data[i][j][2]=C[2];
		}
	}
	cout<<"done making hogel image"<<endl;
	
	save_bmp(N_I,N_J,filename,hogel_data);	
	for(int i=0;i<N_I;i++){
		for(int j=0;j<N_J;j++){
			delete[] hogel_data[i][j];
		}
		delete[] hogel_data[i];
	}
	delete[] hogel_data;


	
}**/


void World::inpaint_hogel_image(int*** hogel_data, int*** hogel_mask, int N_I, int N_J, const char* original_filename, const char* mask_filename, const char* inpainted_filename) {
    // Create OpenCV matrices for the image and mask
    cv::Mat image(N_J, N_I, CV_8UC3);  // Swap N_I and N_J here to match OpenCV's row/column convention
    cv::Mat mask(N_J, N_I, CV_8UC3);   // Same for the mask

    // Fill the cv::Mat with the data from the int*** arrays
    for (int i = 0; i < N_I; i++) {
        for (int j = 0; j < N_J; j++) {
            // Copy the RGB values from the 3D array into the Mat (OpenCV uses BGR by default)
            image.at<cv::Vec3b>(j, i) = cv::Vec3b(  // Swap indices i and j
                static_cast<uchar>(hogel_data[i][j][2]), // B
                static_cast<uchar>(hogel_data[i][j][1]), // G
                static_cast<uchar>(hogel_data[i][j][0])  // R
            );

            // Copy the mask values (if the mask is just 0 or 255, it will work the same way)
            mask.at<cv::Vec3b>(j, i) = cv::Vec3b(  // Swap indices i and j
                static_cast<uchar>(hogel_mask[i][j][2]),
                static_cast<uchar>(hogel_mask[i][j][1]),
                static_cast<uchar>(hogel_mask[i][j][0])
            );
        }
    }

    // Save the original image to disk
    cv::imwrite(original_filename, image);
    std::cout << "Original image saved to " << original_filename << std::endl;

    // Save the mask to disk
    cv::imwrite(mask_filename, mask);
    std::cout << "Mask saved to " << mask_filename << std::endl;

    // Convert mask to grayscale since OpenCV expects a single-channel mask for inpainting
    cv::Mat grayscale_mask;
    cv::cvtColor(mask, grayscale_mask, cv::COLOR_BGR2GRAY);

    // Threshold the mask to make sure it is binary (255 for mask areas, 0 for others)
    cv::threshold(grayscale_mask, grayscale_mask, 1, 255, cv::THRESH_BINARY);

    // Inpainting using the Telea method (you can also use INPAINT_NS for Navier-Stokes method)
    cv::Mat inpainted_image;
    cv::inpaint(image, grayscale_mask, inpainted_image, 3, cv::INPAINT_NS);

    // After inpainting, convert the result back to the int*** format (if needed)
    for (int i = 0; i < N_I; i++) {
        for (int j = 0; j < N_J; j++) {
            cv::Vec3b color = inpainted_image.at<cv::Vec3b>(j, i);  // Swap indices i and j
            hogel_data[i][j][0] = color[2];  // R
            hogel_data[i][j][1] = color[1];  // G
            hogel_data[i][j][2] = color[0];  // B
        }
    }

    // Save the inpainted image to disk
    cv::imwrite(inpainted_filename, inpainted_image);
    std::cout << "Inpainted image saved to " << inpainted_filename << std::endl;
}


void World::render_full_parallax_hogel_tilted(char* filename,char* filename2,char* filename3, int a, int b, int N_a, int N_b, int N_I, int N_J, double distance2, double h_c[],double o[],double e0[],double e1[],double e2[], double Ly, double Lz, double r, double theta_SLM_h,double theta_SLM_v){

	//  A few clarifications.
	// h_c is the point that the camera is looking at.
	// o is the black hole position.  We need it to move back and forth between black hole coordinates and cartesian coordinates.

	int N_rays=2;

	double o_copy[3]={o[0],o[1],o[2]};

	double o1[3],o2[3],o3[3]; //camera basis
	//this->space_time->recompute_orthonormal_frame(o,cam->p,cam->up,e1,e2,e3);
	
	
	double time_interval = 2*3.7;
	
	cout<<"rendering hogel a = "<<a<<endl;

	//initialize image
	hogel_data = new int**[N_I];
	for(int i=0;i<N_I;i++){
		hogel_data[i]=new int*[N_J];
		for(int j=0;j<N_J;j++){
			hogel_data[i][j]=new int[3];
			hogel_data[i][j][0]=0;
			hogel_data[i][j][1]=0;
			hogel_data[i][j][2]=0;
		}
	}
	
	//initialize mask
	hogel_mask = new int**[N_I];
	for(int i=0;i<N_I;i++){
		hogel_mask[i]=new int*[N_J];
		for(int j=0;j<N_J;j++){
			hogel_mask[i][j]=new int[3];
			hogel_mask[i][j][0]=0;
			hogel_mask[i][j][1]=0;
			hogel_mask[i][j][2]=0;
		}
	}

	double V_obs_loc[3]={0,0,0};
	double V_loc[4];
	
	cout<<"initialized hogel image"<<endl;
	
	double p[3];
	
	bool time_upwards=false;
	
	double hy = Ly/N_a;
	double hz = Lz/N_b;
	
	for(int j=0;j<N_J;j++){
		cout<<"rendering j = "<<j<<endl;
		
		double time2;
		if(time_upwards){		
			time2 = 36+2.0+time_interval/2-j*time_interval/(N_J-1);
		}
		else{
			time2 = 36+2.0-time_interval/2+j*time_interval/(N_J-1);
		}
		
		// first let's get our position on the hologram plane.
		double x[3]={0,0,0};
		
		for(int i=0;i<N_I;i++){
			
			//double dy = -Ly/2+h_x/2+(a+float(kk+k_r)/N_rays-0.5)*h_x;
			
			double dy = -Ly/2+(a)*Ly/(N_a);
			double dz = -Lz/2+(b)*Lz/(N_b);
			// staggering
			if(b%2==1){
				dy=dy+hy/2;
			}
			//p[0]=h_c[0];
			//p[1]=h_c[1]+dy;
			//p[2]=h_c[2]+dz;
			p[0]=h_c[0]+dy*e1[0]+dz*e2[0];
			p[1]=h_c[1]+dy*e1[1]+dz*e2[1];
			p[2]=h_c[2]+dy*e1[2]+dz*e2[2];
			for(int ii=0;ii<N_rays;ii++){
				for(int jj=0;jj<N_rays;jj++){
					
					int i_r=((rand() % 100))/100.0; //i.e. rand(0,1.0)
					int j_r=((rand() % 100))/100.0; //i.e. rand(0,1.0)
					
					double d_dy2=2*tan(theta_SLM_h/2)/(N_I-1);
					double d_dz2=2*tan(theta_SLM_v/2)/(N_J-1);
					double delta1 = (float(ii+i_r)/N_rays-0.5)*d_dy2;
					double delta2 = (float(jj+j_r)/N_rays-0.5)*d_dz2;
					//cout<<"rendering i = "<<i<<endl;
					double dy2 = -tan(theta_SLM_h/2)+i*2*tan(theta_SLM_h/2)/(N_I-1)+delta1;
					double dz2 = -tan(theta_SLM_v/2)+j*2*tan(theta_SLM_v/2)/(N_J-1)+delta2;
					//double v[3]={1,dy2,dz2};
					double v[3];
					//v[0]=-e0[0]+dy2*e1[0]+dz2*e2[0];
					//v[1]=-e0[1]+dy2*e1[1]+dz2*e2[1];
					//v[2]=-e0[2]+dy2*e1[2]+dz2*e2[2];
					v[0]=e0[0]+dy2*e1[0]+dz2*e2[0];
					v[1]=e0[1]+dy2*e1[1]+dz2*e2[1];
					v[2]=e0[2]+dy2*e1[2]+dz2*e2[2];
					v[0]=v[0]/sqrt(1+dy2*dy2+dz2*dz2);
					v[1]=v[1]/sqrt(1+dy2*dy2+dz2*dz2);
					v[2]=v[2]/sqrt(1+dy2*dy2+dz2*dz2);
					
					double pmhc_dot_e0=(p[0]-h_c[0])*e0[0]+(p[1]-h_c[1])*e0[1]+(p[2]-h_c[2])*e0[2];
					double pmhc_dot_e1=(p[0]-h_c[0])*e1[0]+(p[1]-h_c[1])*e1[1]+(p[2]-h_c[2])*e1[2];
					double v_dot_e0=v[0]*e0[0]+v[1]*e0[1]+v[2]*e0[2];
					double v_dot_e1=v[0]*e1[0]+v[1]*e1[1]+v[2]*e1[2];
					
					double A = v_dot_e0*v_dot_e0+v_dot_e1*v_dot_e1;
					double B = 2*(pmhc_dot_e0*v_dot_e0+pmhc_dot_e1*v_dot_e1);
					double CC = pmhc_dot_e0*pmhc_dot_e0+pmhc_dot_e1*pmhc_dot_e1-r*r;
					
					//double A = v[0]*v[0]+v[1]*v[1];
					//double B = 2*v[0]*(p[0]-h_c[0])+2*v[1]*(p[1]-h_c[1]);
					//double CC = (p[0]-h_c[0])*(p[0]-h_c[0])+(p[1]-h_c[1])*(p[1]-h_c[1])-r*r;
					double desc = B*B-4*A*CC;
					if(desc<0){
						cout<<"desc = "<<desc<<endl;
					}
					double t1 = (-B+sqrt(desc))/(2*A);
					double t2 = (-B-sqrt(desc))/(2*A);
					double t;
					double l1_dot_e0=(p[0]+t1*v[0])*e0[0]+(p[1]+t1*v[1])*e0[1]+(p[2]+t1*v[2])*e0[2];
					double l2_dot_e0=(p[0]+t2*v[0])*e0[0]+(p[1]+t2*v[1])*e0[1]+(p[2]+t2*v[2])*e0[2];
					if(l1_dot_e0>=l2_dot_e0){
						t=t1;
					}
					else{
						t=t2;
					}
					//double t = fmax(t1,t2);
					double p_c[3]={p[0]+t*v[0],p[1]+t*v[1],p[2]+t*v[2]};
					
					// distance from the black hole is not r.
					double r_cylinder = sqrt((p_c[0]-o[0])*(p_c[0]-o[0])+(p_c[1]-o[1])*(p_c[1]-o[1])+(p_c[2]-o[2])*(p_c[2]-o[2]));
					//P_obs_loc[1] = r_cylinder;
					//double time_delta = fabs(distance2-r);
					double time_delta = fabs(distance2-r);
					double P_obs_loc[4] = {time2+time_delta,distance2,0,1};
					
					double V[3]={p[0]-p_c[0],p[1]-p_c[1],p[2]-p_c[2]}; // this is the ray we should ray trace.
					o[0]=o_copy[0];
					o[1]=o_copy[1];
					o[2]=o_copy[2];
					this->space_time->recompute_orthonormal_frame2(o,p_c,o1,o2,o3);
					double length = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
					V[0]=V[0]/length;
					V[1]=V[1]/length;
					V[2]=V[2]/length;
					// perturb V, for anti-aliasing.
			
			
					int C[3]={0,0,0};
					
					
					//double V_perturbed[3]={V[0]+delta1*e1[0]+delta2*e2[0],V[1]+delta1*e1[1]+delta2*e2[1],V[2]+delta1*e1[2]+delta2*e2[2]};
					//double V_perturbed[3]={V[0]+delta1*o2[0]+delta2*o3[0],V[1]+delta1*o2[1]+delta2*o3[1],V[2]+delta1*o2[2]+delta2*o3[2]};
					//double length = sqrt(V_perturbed[0]*V_perturbed[0]+V_perturbed[1]*V_perturbed[1]+V_perturbed[2]*V_perturbed[2]);
					//V_perturbed[0]=V_perturbed[0]/length;
					//V_perturbed[1]=V_perturbed[1]/length;
					//V_perturbed[2]=V_perturbed[2]/length;

					double psi;
					this->space_time->get_null_vector_in_direction(o1,o2,o3,P_obs_loc,V_obs_loc,V,V_loc,psi);
					
					double H = this->space_time->get_h0(P_obs_loc);
					//this->space_time->scale_velocity(V,p_c);
					double E = this->space_time->compute_energy(P_obs_loc);
					int last_index = -1;
					try{
						if(linear_ray_hits_something(i,j,p_c,V)){
							this->non_linear_ray_trace(i,j,P_obs_loc,V_loc,C,0,H,E,x,last_index,o,o1,o2,o3,psi);
						}
						else{
							C[0]=0;
							C[1]=0;
							C[2]=0;
						}
					}
					catch(const std::exception& e){
						hogel_mask[i][j][0]=255;
						hogel_mask[i][j][1]=255;
						hogel_mask[i][j][2]=255;
					}
					hogel_data[i][j][0]+=C[0];
					hogel_data[i][j][1]+=C[1];
					hogel_data[i][j][2]+=C[2];
				}
			}
			hogel_data[i][j][0]/=(N_rays*N_rays);
			hogel_data[i][j][1]/=(N_rays*N_rays);
			hogel_data[i][j][2]/=(N_rays*N_rays);
		}
	}
	cout<<"done making hogel image"<<endl;
	
	inpaint_hogel_image(hogel_data,hogel_mask,N_I,N_J,filename,filename2,filename3);
	
	//save_bmp(N_I,N_J,filename,hogel_data);	
	//save_bmp(N_I,N_J,filename2,hogel_mask);	
	for(int i=0;i<N_I;i++){
		for(int j=0;j<N_J;j++){
			delete[] hogel_data[i][j];
		}
		delete[] hogel_data[i];
	}
	delete[] hogel_data;
	for(int i=0;i<N_I;i++){
		for(int j=0;j<N_J;j++){
			delete[] hogel_mask[i][j];
		}
		delete[] hogel_mask[i];
	}
	delete[] hogel_mask;


	
}


//renders the world as seen from camera cam as a w x h image
//saved under filename

void World::render(Camera* cam, char* filename,int w, int h){

//make image
	im_data = new int**[w];
	for(int i=0;i<w;i++){
		im_data[i]=new int*[h];
		for(int j=0;j<h;j++){
			im_data[i][j]=new int[3];
			im_data[i][j][0]=0;
			im_data[i][j][1]=0;
			im_data[i][j][2]=0;
		}
	}

		

//render

vector< vector<double> > Ps,Vs;
double P[4],V[3];
double P_loc[4],V_loc[3];
double x[3]={0,0,0};
int n,C[3];
C[0]=0;
C[1]=0;
C[2]=0;

double o[3],e1[3],e2[3],e3[3]; //camera basis
cam->get_basis(o,e1,e2,e3);
double psi;

	int thread_id = omp_get_thread_num();
	//iterate over pixels
    for(int i=0;i<w;i++){	
        std::cout << "thread <<"<<thread_id<<" rendering . . . i="<<i<<"/"<<w<<std::endl;
        for(int j=0;j<h;j++){
			//if(i>=288){
			//std::cout << "thread <<"<<thread_id<<" rendering . . . j="<<j<<"/"<<h<<std::endl;
			//}
			im_data[i][j][0]=0;im_data[i][j][1]=0;im_data[i][j][2]=0;
			//each pixel will in general correspond to several rays
			//cout<<" getting rays"<<endl;

            //cam->get_rays(i,j+int(0.5*(w-h)),w,w,&Ps,&Vs);
			//cam->get_rays_geola(i,j,w,h,&Ps,&Vs);
			//#pragma omp critical
			{
			cam->get_rays_chimera(i,j,w,h,&Ps,&Vs);
			}
			//cout<<" got rays"<<endl;
			
			n=Ps.size();
			//cout<<" num rays = "<<n<<endl;
			//iterate over the rays
			//#pragma omp critical
			{
			for(int k=0;k<n;k++){
				
				P[0]=Ps[k][0];P[1]=Ps[k][1];P[2]=Ps[k][2];
				V[0]=Vs[k][0];V[1]=Vs[k][1];V[2]=Vs[k][2];
				//std::cout << "thread <<"<<thread_id<<" rendering . . . k="<<k<<"/"<<n<<std::endl;
				try{
					this->space_time->cartesian_to_local(e1,e2,e3,P,V,P_loc,V_loc,psi);
					//double H = this->space_time->get_h0(P);
					double H = this->space_time->get_h0(P_loc);
					//this->space_time->scale_velocity(V,P);
					this->space_time->scale_velocity(V_loc,P_loc);
					//double E = this->space_time->compute_energy(P);
					double E = this->space_time->compute_energy(P_loc);
					int last_index = -1;
					//if(linear_ray_hits_something(i,j,P,V)){
					this->non_linear_ray_trace(i,j,P_loc,V_loc,C,0,H,E,x,last_index,o,e1,e2,e3,psi);
					
					//}
					//else{
					//	C[0]=0;
					//	C[1]=0;
					//	C[2]=0;
					//}
				}
				catch(...){
					cout<<"ray tracing failed!"<<endl;
				}
				
				for(int l=0;l<3;l++){
					im_data[i][j][l]+=C[l];
				}
			}
			//average the rays
			im_data[i][j][0]/=n;im_data[i][j][1]/=n;im_data[i][j][2]/=n;
			}
		}
    }
			
	/**for(int i=0;i<w;i++){

		cout << "rendering . . . "<<i<<"/"<<w<<endl;
		for(int j=0;j<h;j++){
			cam->get_ray(i,j,w,h,P,V);
			this->ray_trace(P,V,image[i][j],0);

			
		}
	}**/
	
//save

save_bmp(w,h,filename,im_data);	
for(int i=0;i<w;i++){
	for(int j=0;j<h;j++){
		delete[] im_data[i][j];
	}
	delete[] im_data[i];
}
delete[] im_data;

}




void World::encode(double x, int &R, int &G, int &B) {
    // Define the range maximum
    const double L = 100.0; // Adjust this according to your actual L
    
    // Clamp x to [0, L] just in case
    if (x < 0.0) x = 0.0;
    if (x > L)   x = L;
    
    // Normalize x to [0,1]
    double normalized = x / L;
    
    // Scale up to 24-bit integer range: 0 to 16777215 (2^24 - 1)
    unsigned int maxVal = 16777215; // 0xFFFFFF
    unsigned int value = static_cast<unsigned int>(std::floor(normalized * maxVal));

    // Extract R, G, B from the 24-bit value
    R = (value >> 16) & 0xFF;
    G = (value >> 8)  & 0xFF;
    B = value & 0xFF;
}


//renders the world as seen from camera cam as a w x h image
//saved under filename


void World::spherically_symmetric_fartime_table(double P_obs_loc[], char* filename){
	
	std::ofstream outFile(filename);
	if (!outFile) {
		std::cerr << "Error opening file!" << std::endl;
		return;
	}
	
	double V_loc[3]={0,0,0};
	double P_loc_bu[4]={P_obs_loc[0],P_obs_loc[1],P_obs_loc[2],P_obs_loc[3]};
	double near_angle;
	int N_angles = 500;
	//double a = 1;
	//double R = (P_loc[1]-P_loc[0])/2.0;
	double phi_far;
	
	for(int i=0;i<N_angles;i++){
	
		near_angle = (PI/2)*(i)/(N_angles-1);
	
		P_obs_loc[0]=P_loc_bu[0];
		P_obs_loc[1]=P_loc_bu[1];
		P_obs_loc[2]=P_loc_bu[2];
		P_obs_loc[3]=P_loc_bu[3];
	
		double r = P_obs_loc[1];
		
		int time_dir = -1;
		
		V_loc[0]=time_dir/sqrt(1-2*this->space_time->M/r);
		V_loc[1]=-cos(near_angle)*sqrt(1-2*this->space_time->M/r);
		V_loc[2]=sin(near_angle)/r;
	
		double H = this->space_time->get_h0(P_obs_loc);
		double E = this->space_time->compute_energy(P_obs_loc);

		double far_time = get_far_time(P_obs_loc,V_loc,H,E,phi_far);

		double h = 0.01;

		P_obs_loc[0]=P_loc_bu[0]+h;
		P_obs_loc[1]=P_loc_bu[1];
		P_obs_loc[2]=P_loc_bu[2];
		P_obs_loc[3]=P_loc_bu[3];
		
		V_loc[0]=time_dir/sqrt(1-2*this->space_time->M/r);
		V_loc[1]=-cos(near_angle)*sqrt(1-2*this->space_time->M/r);
		V_loc[2]=sin(near_angle)/r;
		
		H = this->space_time->get_h0(P_obs_loc);
		E = this->space_time->compute_energy(P_obs_loc);
		
		double far_time2 = get_far_time(P_obs_loc,V_loc,H,E,phi_far);
		
		outFile << near_angle << "\t" << far_time <<"\t"<<far_time2<<"\n";
	}
	outFile.close();
	return;
}


void World::spherically_symmetric_dimmness_image(Camera* cam, char* filename,char *time_name,int w, int h,double P_obs_loc[],double V_obs_loc[]){
	
	double V_loc[3]={0,0,0};
	double P_loc_bu[4]={P_obs_loc[0],P_obs_loc[1],P_obs_loc[2],P_obs_loc[3]};
	double near_angle;
	int N_angles = 500;
	//double a = 1;
	//double R = (P_loc[1]-P_loc[0])/2.0;
	double phi_far;
	
	double dimness[500];
	
	for(int i=0;i<N_angles;i++){
	
		near_angle = (PI/2)*(i)/(N_angles-1);
	
		P_obs_loc[0]=P_loc_bu[0];
		P_obs_loc[1]=P_loc_bu[1];
		P_obs_loc[2]=P_loc_bu[2];
		P_obs_loc[3]=P_loc_bu[3];
	
		double r = P_obs_loc[1];
		
		int time_dir = -1;
		
		V_loc[0]=time_dir/sqrt(1-2*this->space_time->M/r);
		V_loc[1]=-cos(near_angle)*sqrt(1-2*this->space_time->M/r);
		V_loc[2]=sin(near_angle)/r;
	
		double H = this->space_time->get_h0(P_obs_loc);
		double E = this->space_time->compute_energy(P_obs_loc);

		dimness[i] = get_far_time(P_obs_loc,V_loc,H,E,phi_far);
		
	}

	cout<<" dimness[1] = "<<dimness[1]<<endl;

//make image
	im_data = new int**[w];
	for(int i=0;i<w;i++){
		im_data[i]=new int*[h];
		for(int j=0;j<h;j++){
			im_data[i][j]=new int[3];
			im_data[i][j][0]=0;
			im_data[i][j][1]=0;
			im_data[i][j][2]=0;
		}
	}
	
//make image
	/**time_data = new int**[w];
	for(int i=0;i<w;i++){
		time_data[i]=new int*[h];
		for(int j=0;j<h;j++){
			time_data[i][j]=new int[3];
			time_data[i][j][0]=0;
			time_data[i][j][1]=0;
			time_data[i][j][2]=0;
		}
	}**/

		

//render

	vector< vector<double> > Ps,Vs;
	double P[4],V[3];
	double x[3]={0,0,0};
	int n,C[3];
	C[0]=0;
	C[1]=0;
	C[2]=0;

	double o[3],e1[3],e2[3],e3[3]; //camera basis
	cam->get_basis(o,e1,e2,e3);
	double psi;
	o[0]=2.0;
	o[1]=0.0;
	o[2]=-3.0;

	int thread_id = omp_get_thread_num();
	
	//iterate over pixels
    for(int i=0;i<w;i++){	
        std::cout << "thread <<"<<thread_id<<" rendering . . . i="<<i<<"/"<<w<<std::endl;
        for(int j=0;j<h;j++){
			//if(i>=288){
			//std::cout << "thread <<"<<thread_id<<" rendering . . . j="<<j<<"/"<<h<<std::endl;
			//}
			im_data[i][j][0]=0;im_data[i][j][1]=0;im_data[i][j][2]=0;
			//each pixel will in general correspond to several rays
			//cout<<" getting rays"<<endl;

            //cam->get_rays(i,j+int(0.5*(w-h)),w,w,&Ps,&Vs);
			//cam->get_rays_geola(i,j,w,h,&Ps,&Vs);
			//#pragma omp critical
			{
			cam->get_rays_chimera(i,j,w,h,&Ps,&Vs);
			}
			//cout<<" got rays"<<endl;
			
			n=Ps.size();
			//cout<<" num rays = "<<n<<endl;
			//iterate over the rays
			//#pragma omp critical
			{
			for(int k=0;k<n;k++){
				
				P[0]=Ps[k][0];P[1]=Ps[k][1];P[2]=Ps[k][2];
				V[0]=Vs[k][0];V[1]=Vs[k][1];V[2]=Vs[k][2];
				
				
				
				
				//std::cout << "thread <<"<<thread_id<<" rendering . . . k="<<k<<"/"<<n<<std::endl;
				try{

					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					double v_dot_e1 = V[0]*e1[0]+V[1]*e1[1]+V[2]*e1[2];
					double v_dot_e2 = V[0]*e2[0]+V[1]*e2[1]+V[2]*e2[2];
					double v_dot_e3 = V[0]*e3[0]+V[1]*e3[1]+V[2]*e3[2];
		
					psi =atan2(v_dot_e3,v_dot_e2);
					double epsi[3] = {cos(psi)*e2[0]+sin(psi)*e3[0],cos(psi)*e2[1]+sin(psi)*e3[1],cos(psi)*e2[2]+sin(psi)*e3[2]};
					double v_dot_epsi = V[0]*epsi[0]+V[1]*epsi[1]+V[2]*epsi[2];
					double phi = atan2(v_dot_epsi,v_dot_e1);
					
					double I = (phi/(PI/2))*500;
					int Im = int(floor(I));
					int Ip = min(Im+1,500);
					double s = I-Im;
					double brightness = dimness[Im]*(1-s)+dimness[Ip]*s;
					C[0]=brightness*255;
					C[1]=brightness*255;
					C[2]=brightness*255;
					
					
				}
				catch(...){
					cout<<"ray tracing failed!"<<endl;
				}
				
				for(int l=0;l<3;l++){
					im_data[i][j][l]+=C[l];
				}
				//time_data[i][j][0]=int(floor(255*dt/(far_time-far_time_prev)));
				//time_data[i][j][1]=int(floor(255*dt/(far_time-far_time_prev)));
				//time_data[i][j][2]=int(floor(255*dt/(far_time-far_time_prev)));
				//int R,G,B;
				//this->encode(far_time, R, G, B);
				//time_data[i][j][0]=R;
				//time_data[i][j][1]=G;
				//time_data[i][j][2]=B;
			}
			//average the rays
			im_data[i][j][0]/=n;im_data[i][j][1]/=n;im_data[i][j][2]/=n;
			}
		}
    }
	//for(int i=0;i<w;i++){
	//	for(int j=0;j<h;j++){
	//		time_data[i][j][0]/=n;time_data[i][j][1]/=n;time_data[i][j][2]/=n;
	//	}
	//}
	/**for(int i=0;i<w;i++){

		cout << "rendering . . . "<<i<<"/"<<w<<endl;
		for(int j=0;j<h;j++){
			cam->get_ray(i,j,w,h,P,V);
			this->ray_trace(P,V,image[i][j],0);

			
		}
	}**/
	
//save

save_bmp(w,h,filename,im_data);	
for(int i=0;i<w;i++){
	for(int j=0;j<h;j++){
		delete[] im_data[i][j];
	}
	delete[] im_data[i];
}
delete[] im_data;


}

void World::render_spherically_symmetric_spacetime(Camera* cam,double o[], char* filename,char *time_name,int w, int h,double P_obs_loc[],double V_obs_loc[]){

//make image
	im_data = new int**[w];
	for(int i=0;i<w;i++){
		im_data[i]=new int*[h];
		for(int j=0;j<h;j++){
			im_data[i][j]=new int[3];
			im_data[i][j][0]=0;
			im_data[i][j][1]=0;
			im_data[i][j][2]=0;
		}
	}
	
//make image
	/**time_data = new int**[w];
	for(int i=0;i<w;i++){
		time_data[i]=new int*[h];
		for(int j=0;j<h;j++){
			time_data[i][j]=new int[3];
			time_data[i][j][0]=0;
			time_data[i][j][1]=0;
			time_data[i][j][2]=0;
		}
	}**/

		

//render

vector< vector<double> > Ps,Vs;
double P[4],V[3];
double V_loc[4];
double x[3]={0,0,0};
double P_obs_loc_bu[4] = {P_obs_loc[0],P_obs_loc[1],P_obs_loc[2],1};
double V_obs_loc_bu[3] = {V_obs_loc[0],V_obs_loc[1],V_obs_loc[2]};
int n,C[3];
C[0]=0;
C[1]=0;
C[2]=0;
double e1[3],e2[3],e3[3]; //camera basis
//this->space_time->recompute_orthonormal_frame(o,cam->p,cam->up,e1,e2,e3);
this->space_time->recompute_orthonormal_frame2(o,cam->p,e1,e2,e3);
double psi;
//o[2]=-2.0;
//o[2]=-1.5;


double dt = 0.01;
double far_time,far_time_prev;
double phi = fabs(atan2(e1[1],e1[0]));

	int thread_id = omp_get_thread_num();
	
	int w_l = int(floor(w/4.0));
	int h_l = int(floor(h/4.0));
	int w_u = int(floor(3*w/4.0));
	int h_u = int(floor(3*h/4.0));
	int w_m = int(floor(w/2.0));
	int h_m = int(floor(h/2.0));
	
	//iterate over pixels
    for(int i=0;i<w;i++){	
        std::cout << "thread <<"<<thread_id<<" rendering . . . i="<<i<<"/"<<w<<std::endl;
        for(int j=0;j<h;j++){
			//if(i>=288){
			//std::cout << "thread <<"<<thread_id<<" rendering . . . j="<<j<<"/"<<h<<std::endl;
			//}
			im_data[i][j][0]=0;im_data[i][j][1]=0;im_data[i][j][2]=0;
			//each pixel will in general correspond to several rays
			//cout<<" getting rays"<<endl;

            //cam->get_rays(i,j+int(0.5*(w-h)),w,w,&Ps,&Vs);
			//cam->get_rays_geola(i,j,w,h,&Ps,&Vs);
			//#pragma omp critical
			{
			cam->get_rays_chimera(i,j,w,h,&Ps,&Vs);
			}
			//cout<<" got rays"<<endl;
			
			n=Ps.size();
			//cout<<" num rays = "<<n<<endl;
			//iterate over the rays
			//#pragma omp critical
			{
			for(int k=0;k<n;k++){
				
				P[0]=Ps[k][0];P[1]=Ps[k][1];P[2]=Ps[k][2];
				V[0]=Vs[k][0];V[1]=Vs[k][1];V[2]=Vs[k][2];
				
				
				
				
				//std::cout << "thread <<"<<thread_id<<" rendering . . . k="<<k<<"/"<<n<<std::endl;
				try{
					P_obs_loc[0]=P_obs_loc_bu[0];
					P_obs_loc[1]=P_obs_loc_bu[1];
					P_obs_loc[2]=P_obs_loc_bu[2];
					P_obs_loc[3]=P_obs_loc_bu[3];
					V_obs_loc[0]=V_obs_loc_bu[0];
					V_obs_loc[1]=V_obs_loc_bu[1];
					V_obs_loc[2]=V_obs_loc_bu[2];
					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					this->space_time->get_null_vector_in_direction(e1,e2,e3,P_obs_loc,V_obs_loc,V,V_loc,psi);
					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					//double H = this->space_time->get_h0(P);
					double H = this->space_time->get_h0(P_obs_loc);
					double P_copy[3]={P[0],P[1],P[2]};
					double V_copy[3]={V[0],V[1],V[2]};
					
					//double E = phi;
					//this->space_time->scale_velocity(V,P);
					//this->space_time->scale_velocity(V_loc,P_loc);
					double E = this->space_time->compute_energy(P_obs_loc);
					//double E = this->space_time->compute_energy(P_loc);
					
					int last_index = -1;
					if(linear_ray_hits_something(i,j,P_copy,V_copy)){
					
						this->non_linear_ray_trace(i,j,P_obs_loc,V_loc,C,0,H,E,x,last_index,o,e1,e2,e3,psi);
					}
					
					//P_loc[0]=P_loc_bu[0];
					//P_loc[1]=P_loc_bu[1];
					//P_loc[2]=P_loc_bu[2];
					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					//this->space_time->get_null_vector_in_direction(e1,e2,e3,P_loc,V,V_loc,psi);
					
					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					//double H = this->space_time->get_h0(P);
					//H = this->space_time->get_h0(P_loc);
					
					
					//E = phi;
					//this->space_time->scale_velocity(V,P);
					//this->space_time->scale_velocity(V_loc,P_loc);
					//double E = this->space_time->compute_energy(P);
					//double E = this->space_time->compute_energy(P_loc);
					
					//last_index = -1;
					//if(linear_ray_hits_something(i,j,P,V)){
					
					//far_time = get_far_time(i,j,P_loc,V_loc,0,H,E,x,last_index,o,e1,e2,e3,psi);
					
					/**P_loc[0]=P_loc_bu[0]-dt;
					P_loc[1]=P_loc_bu[1]-dt;
					P_loc[2]=P_loc_bu[2];
					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					this->space_time->get_null_vector_in_direction(e1,e2,e3,P_loc,V,V_loc,psi);
					
					//cout<<"V_loc[0] = "<<V_loc[0]<<"V_loc[1] = "<<V_loc[1]<<"V_loc[2] = "<<V_loc[2]<<endl;
					//double H = this->space_time->get_h0(P);
					H = this->space_time->get_h0(P_loc);
					
					
					E = phi;
					//this->space_time->scale_velocity(V,P);
					//this->space_time->scale_velocity(V_loc,P_loc);
					//double E = this->space_time->compute_energy(P);
					//double E = this->space_time->compute_energy(P_loc);
					
					last_index = -1;
					//if(linear_ray_hits_something(i,j,P,V)){
					
					far_time_prev = get_far_time(i,j,P_loc,V_loc,0,H,E,x,last_index,o,e1,e2,e3,psi);**/
					
					//cout<<"far_time, C[3], dt, ratio = "<<far_time<<" "<<double(C[3])/10000<<" "<<dt<<" "<<dt/(far_time-double(C[3])/10000)<<endl;
					
					//}
					//else{
					//	C[0]=0;
					//	C[1]=0;
					//	C[2]=0;
					//}
				}
				catch(...){
					cout<<"ray tracing failed!"<<endl;
				}
				
				for(int l=0;l<3;l++){
					im_data[i][j][l]+=C[l];
				}
				//time_data[i][j][0]=int(floor(255*dt/(far_time-far_time_prev)));
				//time_data[i][j][1]=int(floor(255*dt/(far_time-far_time_prev)));
				//time_data[i][j][2]=int(floor(255*dt/(far_time-far_time_prev)));
				//int R,G,B;
				//this->encode(far_time, R, G, B);
				//time_data[i][j][0]=R;
				//time_data[i][j][1]=G;
				//time_data[i][j][2]=B;
			}
			//average the rays
			im_data[i][j][0]/=n;im_data[i][j][1]/=n;im_data[i][j][2]/=n;
			}
		}
    }
	//for(int i=0;i<w;i++){
	//	for(int j=0;j<h;j++){
	//		time_data[i][j][0]/=n;time_data[i][j][1]/=n;time_data[i][j][2]/=n;
	//	}
	//}
	/**for(int i=0;i<w;i++){

		cout << "rendering . . . "<<i<<"/"<<w<<endl;
		for(int j=0;j<h;j++){
			cam->get_ray(i,j,w,h,P,V);
			this->ray_trace(P,V,image[i][j],0);

			
		}
	}**/
	
//save

save_bmp(w,h,filename,im_data);	
for(int i=0;i<w;i++){
	for(int j=0;j<h;j++){
		delete[] im_data[i][j];
	}
	delete[] im_data[i];
}
delete[] im_data;

//save_bmp(w,h,time_name,time_data);	
//for(int i=0;i<w;i++){
//	for(int j=0;j<h;j++){
//		delete[] time_data[i][j];
//	}
//	delete[] time_data[i];
//}
//delete[] time_data;

}


//makes and renders the photon map as w x h textures for each object in the scene.

void World::make_photon_map(int w, int h){

	int n_photons = 2500000;
	int n_objects = shapes.size();
	double photon_energy = 15000/double(n_photons);
	int last_index = -1;
	bool caustic = false;
	
	vector< vector<Photon>> photons;
	photons.clear();
	for(int i=0;i<n_objects;i++){
		vector<Photon> empty_vector;
		empty_vector.clear();
		photons.push_back(empty_vector);
	}

	for(int i=0;i<photon_lights.size();i++){
		
		cout<<"making photon map....i = "<<i<<" of "<<photon_lights.size()<<endl;
		
		vector< vector<double> > P;
		vector< vector<double> > V;
		P.clear();
		V.clear();
		photon_lights[i]->get_Photons(n_photons,&P,&V);
		for(int j=0;j<P.size();j++){
			double P1[3] = {P[j][0],P[j][1],P[j][2]};
			double V1[3] = {V[j][0],V[j][1],V[j][2]};
			double H = this->space_time->get_h0(P1);
			this->space_time->scale_velocity(V1,P1);
			double E = this->space_time->compute_energy(P1);
			this->non_linear_photon_emission(P1,V1,photons,photon_energy,0,H,E,last_index,caustic);
		}
	}
	
	//make image
	im_data = new int**[w];
	for(int i=0;i<w;i++){
		im_data[i]=new int*[h];
		for(int j=0;j<h;j++){
			im_data[i][j]=new int[3];
			im_data[i][j][0]=0;
			im_data[i][j][1]=0;
			im_data[i][j][2]=0;
		}
	}	
	
	// now let's make the textures
	for(int k=0;k<n_objects;k++){
		cout<<"length of photons[k] = "<<photons[k].size()<<endl;
		string texture_name = "./PhotonMap/Object_"+to_string(k)+"_Geola.bmp";
		char *filename = new char[texture_name.length() + 1];
		strcpy(filename, texture_name.c_str());
		string texture_name2 = "./PhotonMap/ObjectDots_"+to_string(k)+"_Geola.bmp";
		char *filename2 = new char[texture_name2.length() + 1];
		strcpy(filename2, texture_name2.c_str());
		
		// Initialize the PhotonCloud with the vector of photons
		const PhotonCloud cloud(photons[k]);

		// Create a KD-tree index
		//PhotonKDTree index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
		PhotonKDTree index(3 /*dim*/, cloud, {10 /* max leaf */});
		index.buildIndex();
		
		double alpha = 0.45;
		
		// Query point
		float xx[3] = {0, 0, 0};
		
		//size_t num_results = size_t(floor(pow(photons[k].size(),alpha)));
		size_t num_results = 400;
		cout<<"num_results = "<<num_results<<endl;
		
		for(int i=0;i<w;i++){
			
			//cout<<"Baking photon map....object = "<<k<<" of "<<n_objects<<" i = "<<i<<endl;
			
			for(int j=0;j<h;j++){
				shapes[k]->inverse_texture(i,j,w,h,xx);
				const float query_pt[3] = {xx[0], xx[1], xx[2]};
				// Number of nearest neighbors to find
				num_results = size_t(floor(pow(photons[k].size(),alpha)));
				std::vector<uint32_t> ret_index(num_results);
				std::vector<float> out_dist_sqr(num_results);


				//KNNResultSet<float> resultSet(num_results);
				//resultSet.init(&ret_index[0], &out_dist_sqr[0]);
				//index.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(5000));
				num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);
				
				// Find the maximum squared distance
				float max_dist_sqr = 0;
				for (size_t ii = 0; ii < num_results; ++ii) {
					if (out_dist_sqr[ii] > max_dist_sqr) {
						max_dist_sqr = out_dist_sqr[ii];
					}
				}
				// Find the maximum value in out_dist_sqr
				/**float max_dist_sqr = 0;
				auto max_it = std::max_element(out_dist_sqr.begin(), out_dist_sqr.end());
				if (max_it != out_dist_sqr.end()) {
					float max_dist_sqr = *max_it;
					//std::cout << "The maximum squared distance is: " << max_value << std::endl;
				} else {
					//std::cout << "The vector is empty." << std::endl;
				}**/
				
				// Calculate the total energy of the nearest neighbors
				//float total_energy = 0.0f;
				float total_energy = num_results*photon_energy;
				for (size_t ii = 0; ii < num_results; ++ii) {
					float dp2 = (photons[k][ret_index[ii]].x-query_pt[0])*(photons[k][ret_index[ii]].x-query_pt[0])+(photons[k][ret_index[ii]].y-query_pt[1])*(photons[k][ret_index[ii]].y-query_pt[1])+(photons[k][ret_index[ii]].z-query_pt[2])*(photons[k][ret_index[ii]].z-query_pt[2]);
					float wg = 1.818*(1-(1-exp(-1.953*dp2/(2*max_dist_sqr)))/(1-exp(-1.953)));
					total_energy += wg*photons[k][ret_index[ii]].energy;
				}
				int brightness = min(int(floor(total_energy/(PI*max_dist_sqr))),255);
				im_data[i][j][0]=brightness;
				im_data[i][j][1]=brightness;
				im_data[i][j][2]=brightness;
			}
		}
		save_bmp(w,h,filename,im_data);	
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++){
				im_data[i][j][0]=0;
				im_data[i][j][1]=0;
				im_data[i][j][2]=0;
			}
		}
		for(int ii=0;ii<photons[k].size();ii++){
			double xx[3]={photons[k][ii].x,photons[k][ii].y,photons[k][ii].z};
			shapes[k]->place_dot(xx,im_data,w,h);
		}
		save_bmp(w,h,filename2,im_data);
	}
	
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			delete[] im_data[i][j];
		}
		delete[] im_data[i];
	}
	delete[] im_data;
	

}













void World::shade(double x[],double n[],double r[], int col[]){

double term[3],sum[3];
double specular_term[3], specular_sum[3];
sum[0]=0;sum[1]=0;sum[2]=0;
specular_sum[0]=0;specular_sum[1]=0;specular_sum[2]=0;
double dp[3]={0,0,0};
for(int i=0;i<lights.size();i++){
	if(!soft_shadows_on){
		if(lights[i]->is_Visible(x,invis_R,dp,shapes) || !shadows_on){
			// dp is now x-light_pos
			lights[i]->get_Contribution(n,x,dp,textures,invis_R,term);
			lights[i]->get_Specular_Contribution(r,x,specular_term);
			sum[0]+=term[0];
			sum[1]+=term[1];
			sum[2]+=term[2];
			specular_sum[0]=+specular_term[0];
			specular_sum[1]=+specular_term[1];
			specular_sum[2]=+specular_term[2];
		}
	}
	if(soft_shadows_on){
		double visibility = lights[i]->percent_Visible(x,invis_R,shapes);
		lights[i]->get_Contribution(n,x,dp,textures,invis_R,term);
		lights[i]->get_Specular_Contribution(r,x,specular_term);
		sum[0]+=visibility*term[0];
		sum[1]+=visibility*term[1];
		sum[2]+=visibility*term[2];
		specular_sum[0]=+visibility*specular_term[0];
		specular_sum[1]=+visibility*specular_term[1];
		specular_sum[2]=+visibility*specular_term[2];
	}
}

for(int i=0;i<3;i++){
	col[i]=int(min(int(col[i]*sum[i]),min(int(2*col[i]),255)));
	col[i]=int(col[i]+specular_sum[i]);
	//col[i]=int(col[i]*sum[i]);
	if(col[i]>255){col[i]=255;}
	if(col[i]<0){col[i]=0;}
	}

}
