#include <iostream>
#include <cstdlib>
#include <math.h>
#include "World.h"
#include "Shape.h"
#include "Sphere.h"
#include "Camera.h"
#include "Texture.h"
#include "Quaternion_Julia_Set.h"
#include "Menger_Sponge_4DSlice.h"
#include "Menger_Sponge_4DSlice_Type2.h"
#include "Menger_Sponge.h"
#include "Menger_Sponge_Type2.h"
#include "Mandelbulb.h"
#include "Minkowski.h"
#include "Schwarzschild.h"
#include <string>
#include <cstring>
#include <omp.h>

#define pi 3.1415
#define SQRT3 1.7321

using namespace std;


int main(int argc,char* argv[])
{
//default values for output image filename 
//and width and height, if not entered by command line
char* filename="world.bmp";
int w=1920;
int h=1080;
int eye=1;
double offset = 0;
//attempt to read in these things
if(argc>1){filename=argv[1];}
if(argc>2){w=atoi(argv[2]);}
if(argc>3){h=atoi(argv[3]);}
//if(argc>4){eye=atoi(argv[4]);}
if(argc>4){offset=stod(argv[4]);}

int rat = 30;
string base = "rat";
string extension = to_string(rat);
cout <<(base+extension).c_str()<<endl;




// Set OMP_STACKSIZE to 250MB using _putenv_s
    if (_putenv_s("OMP_STACKSIZE", "250MB") != 0) {
        std::cerr << "Failed to set OMP_STACKSIZE" << std::endl;
        return 1;
    }

    // Verify that OMP_STACKSIZE is set
    char* stacksize = std::getenv("OMP_STACKSIZE");
    if (stacksize) {
        std::cout << "OMP_STACKSIZE is set to: " << stacksize << std::endl;
    } else {
        std::cout << "OMP_STACKSIZE is not set." << std::endl;
    }
	
// Set OMP_DISPLAY_ENV to true using _putenv_s
    if (_putenv_s("OMP_DISPLAY_ENV", "true") != 0) {
        std::cerr << "Failed to set OMP_DISPLAY_ENV" << std::endl;
        return 1;
    }
	
    // Verify that OMP_DISPLAY_ENV is set
    char* display_env = std::getenv("OMP_DISPLAY_ENV");
    if (display_env) {
        std::cout << "OMP_DISPLAY_ENV is set to: " << display_env << std::endl;
    } else {
        std::cout << "OMP_DISPLAY_ENV is not set." << std::endl;
    }

//omp_set_stacksize(2L * 1024 * 1024 * 1024);


//char* stacksize = std::getenv("OMP_STACKSIZE");
//if (stacksize) {
//    std::cout << "OMP_STACKSIZE is set to: " << stacksize << std::endl;
//} else {
//    std::cout << "OMP_STACKSIZE is not set." << std::endl;
//}
//initialize world

//define a camera
	//double p[]={-2,4,7};
	//double p[] = {-2,2,3};
	//double p[]={3,3,3};

	//double l[]={0.5,0.5,0.5};

	//double u[]={0,1,0};
	
	
	//double p[]={2-0.5,-4+0.5,-7+0.5};
	
	//Lense_Camera my_Camera(p,l,u,50,0.1,0.1,6.0622,3); //can also use lense cameras for depth of field
	
	//Lense_Camera my_Camera(p,l,u,50,0.1,0.1,5.39,3); //can also use lense cameras for depth of field
	
	double R=0;double I=0.75;double angle=57.5*pi/180;
	
	//Quaternion_Julia_Set *Q1=new Quaternion_Julia_Set(-2,0,-1,R,I,angle,0.001*1200.0/w,255,0,0);

	//Quaternion_Julia_Set *Q2=new Quaternion_Julia_Set(0,0.25,2,-1,0.25,0,0.001*1200.0/w,0,255,0);

	//Quaternion_Julia_Set *Q3=new Quaternion_Julia_Set(3,0,-3,-1,-0.5,0,0.001*1200.0/w,0,0,255);
	
	//Mandelbulb *M=new Mandelbulb(1,0.5,0,0.001,8,0,255,0);
	double c1=0,c2=0*0.5,c3=0,scale = 1,E=0.001;
	int depth = 4;
	E=0.0001;
	//double offset = 1.0;
	

//define a parallelogram
//this will be the floor
	//double x2[]={-0.5774,0.7887,-0.2113};double x1[]={-0.5774,-0.2113,0.7887};
	//for(int i=0;i<3;i++){
	//	x1[i]*=2;
	//	x2[i]*=2;
	//}
	
	//P1->texture_ID=0;
	//P1->reflectivity=0.5;

//define a bunch of lights

	/**double pos[3]={5,5,5};
	double pos2[3]={0,5,3};
	int color[3]={255,255,255};
	int color2[3]={255,255,255};
	//Light *L=new pointLight(pos,4.5,color);
	Light *L=new pointLight(pos,15,color);
	Light *L2=new ambientLight(0.15,color);
	Light *L3=new pointLight(pos2,10,color);
	
	double pos3[3]={4,0,-1};
	double pos3a[3]={-4,4,-2};
	double pos4[3]={-3,1,1.5};	//double p[] = {-2,2,3};
	double pos4a[3]={-3,1,-1.5};	//double p[] = {-2,2,3};
	
	Light* L4=new pointLight(pos3,6.5,color);
	//Light* L5=new pointLight(pos4,0.5,color);
	Light* L5=new pointLight(pos4,3.0,color);
	
	double pos5[3]={0,0,0};
	double pos6[3]={sqrt(2.0/9)*5,sqrt(2/3.0)*5,(1.0/3)*5};
	double pos6a[3]={sqrt(2.0/9)*5,-sqrt(2/3.0)*5,(1.0/3)*5};
	double pos7[3]={0,0,-3};
	
	Light* L6=new pointLight(pos5,0.5,color);
	Light* L7=new pointLight(pos6,5.0,color);
	Light* L8=new pointLight(pos4a,5.0,color);
	Light* L9=new pointLight(pos7,15.0,color);
	Light* L10=new pointLight(pos3a,5.0,color);
	Light* L11=new pointLight(pos6a,7.0,color);**/
	
	
	double p[3]={3,3,-3};
	double v[3]={1/SQRT3,-1/SQRT3,-1/SQRT3};
	
	const int numThreads = 1;  // Explicitly set the number of threads
	//const int numThreads = 1;  // Explicitly set the number of threads
	std::vector<Texture*> skies(numThreads, nullptr);
	std::vector<Texture*> metals(numThreads, nullptr);
	//std::vector<Texture*> metal_bumps(numThreads, nullptr);
	std::vector<Texture*> gold_spirals(numThreads, nullptr);
	std::vector<Texture*> map_0(numThreads, nullptr);
	std::vector<Texture*> map_1(numThreads, nullptr);
	std::vector<Texture*> map_2(numThreads, nullptr);
	std::vector<Texture*> map_3(numThreads, nullptr);
	std::vector<Texture*> map_4(numThreads, nullptr);
	std::vector<Texture*> map_5(numThreads, nullptr);
	std::vector<Texture*> map_6(numThreads, nullptr);
	std::vector<Texture*> map_7(numThreads, nullptr);
	std::vector<Texture*> map_8(numThreads, nullptr);
	std::vector<Texture*> map_9(numThreads, nullptr);
	
	#pragma omp parallel num_threads(numThreads)
    {
		cout<<"pre-computing space texture..."<<endl;
		//char* space_panorama="BlueSkyAurora.txt";
		char* space_panorama="Tiles_crop.txt";
		char* metal_texture="Tile_wrap.txt";
		char* bump_texture="metal_bump3.txt";
		char* gold_spiral_texture="RoofTile3.txt";
		
		char* photon_map_0 = "Object_0_Geola.txt";
		char* photon_map_1 = "Object_1_Geola.txt";
		char* photon_map_2 = "Object_2_Geola.txt";
		char* photon_map_3 = "Object_3_Geola.txt";
		char* photon_map_4 = "Object_4_Geola.txt";
		char* photon_map_5 = "Object_5_Geola.txt";
		char* photon_map_6 = "Object_6_Geola.txt";
		char* photon_map_7 = "Object_7_Geola.txt";
		char* photon_map_8 = "Object_8_Geola.txt";
		char* photon_map_9 = "Object_9_Geola.txt";
		
		#pragma omp critical
			{
			int threadId = omp_get_thread_num();
			cout<<"thread number "<<threadId<<" of "<<numThreads<<endl;
			skies[threadId] = new Texture(space_panorama);  // Allocate and copy base texture
			metals[threadId] = new Texture(metal_texture);  // Allocate and copy base texture
			//metal_bumps[threadId] = new Texture(bump_texture);  // Allocate and copy base texture
			gold_spirals[threadId] = new Texture(gold_spiral_texture);  // Allocate and copy base texture
			if (skies[threadId] == nullptr) {
				throw std::runtime_error("Failed to allocate memory for thread texture");
			}
			map_0[threadId] = new Texture(photon_map_0);
			map_1[threadId] = new Texture(photon_map_1);
			map_2[threadId] = new Texture(photon_map_2);
			map_3[threadId] = new Texture(photon_map_3);
			map_4[threadId] = new Texture(photon_map_4);
			map_5[threadId] = new Texture(photon_map_5);
			map_6[threadId] = new Texture(photon_map_6);
			map_7[threadId] = new Texture(photon_map_7);
			map_8[threadId] = new Texture(photon_map_8);
			map_9[threadId] = new Texture(photon_map_9);
			//if (metals[threadId] == nullptr) {
			//	throw std::runtime_error("Failed to allocate memory for thread texture");
			//}
			//if (metal_bumps[threadId] == nullptr) {
			//	throw std::runtime_error("Failed to allocate memory for thread texture");
			//}
		}
    }
	
	
	
	
	/**double pos5[3]={2,2,-4};
	double pos6[3]={0,1,-4};	
	
	Light* L6=new pointLight(pos5,4.5,color);
	Light* L7=new pointLight(pos6,0.5,color);**/
	
	
//define some textures
//I do textures in kind of a weird way, by translating images into text files
//and then loading the text files into arrays.

	//char* marble_texture="marble.txt";
	//char* space_texture="Space.txt";
	//char* space_panorama="sky.txt";
	//char* space_panorama="sky2_4096_2048.txt";
	/** move textures */
	//char* space_panorama="BlueSkySmall.txt";
	//char* metal_texture="3metal_texture_big_100617.txt";
	//char* bump_texture="metal_bump3.txt";
	//char* water_texture="Water.txt";
	//char* water_bump="Water_bump5smooth.txt";

	//Texture *T1=new Texture(marble_texture);
	//Texture *T2=new Texture(space_texture);
	//Texture *T1=new Texture(space_panorama);
	/**Texture *T2=new Texture(metal_texture);
	Texture *T3=new Texture(bump_texture);
	Texture *T4=new Texture(water_texture);
	Texture *T5=new Texture(water_bump);**/
	
	/**end texture move**/ 
	
	/** move lights
	int planeColor[3]={0,255,255};
	Light* L12=new planeLight(2.5, 5, planeColor,T4);
	
	
	int planeColor2[3]={255,255,255};
	Light* L13=new planeLight(-5, 5, planeColor2,T4);
	
	end lights **/

//add shapes and light and textures to the world

	//earth.addShape(Q1);
	//earth.addShape(Q2);
	//earth.addShape(Q3);
	/**double o_w[3]={0,0,0};
	double o_x[3]={1,0,0};
	double o_y[3]={0,1,0};
	double o_z[3]={0,0,1};
	double e11[3]={1,0,0};
	double e22[3]={0,1,0};
	double e33[3]={0,0,1};
	Parallelogram *z0=new Parallelogram(o_w,e11,e22,255,255,0);
	z0->is_wireframe=true;
	Parallelogram *z00=new Parallelogram(o_w,e11,e22,127,130,187);
	z00->transparency=0.25;

	Parallelogram *z1=new Parallelogram(o_z,e11,e22,255,255,0);
	z1->is_wireframe=true;
	z1->rat=1;
	Parallelogram *x0=new Parallelogram(o_w,e22,e33,255,255,0);
	x0->is_wireframe=true;
	Parallelogram *x00=new Parallelogram(o_w,e22,e33,127,130,187);
	x00->transparency=0.25;
	Parallelogram *x11=new Parallelogram(o_x,e22,e33,255,255,0);
	x11->is_wireframe=true;
	x11->rat=1;
	Parallelogram *y0=new Parallelogram(o_w,e11,e33,255,255,0);
	y0->is_wireframe=true;
	Parallelogram *y00=new Parallelogram(o_w,e11,e33,127,130,187);
	y00->transparency=0.25;
	Parallelogram *y1=new Parallelogram(o_y,e11,e33,255,255,0);
	y1->is_wireframe=true;
	y1->rat=1;
	
	double p1[3]={-0.25,-0.25,-0.25};
	double p2[3]={1.25,1.25,1.25};
	Cylinder *C = new Cylinder(p1,p2,0.01,100,100,100);**/
	
	

			//double l[]={0.5,0.5,0.5};
	double l[3]={0,0,0};

			//double u[]={0,1,0};
	double u[3]={1/SQRT3,-1/SQRT3,-1/SQRT3};
	
	

//render world and write result to file

	//string filename_base = "./calcul/";
	string filename_base = "./Loic_HP_Hogel_Staggered/";
	string filename_middle = "hogelOut";
	string filename_extension = ".bmp";
	string underscore = "_";
	
	double p_back[3] = {-9.41805, -18.2564, -3.16031};
	//Light* L6 = new pointLight(p_back,16*6.0,color);


	//for(int i=620;i<621;i++){ 
	int N = 768;
	//int M = 240;
	//double dM = double(floor(M/47.0));
	//#pragma omp parallel for private(p,l,u)

	std::vector<Texture*> waters(numThreads, nullptr);
	std::vector<Texture*> water_bumps(numThreads, nullptr);
	//std::vector<Texture*> fires(numThreads, nullptr);
	//for(int ii=0;ii<48;ii++){
	//int N_a = 100;
	int N_a = 150;
	int N_I = 1920;
try { 
	//#pragma omp parallel for private(p,l,u,w,h) num_threads(numThreads) schedule(dynamic)
	for(int i=0;i<15;i++){ 
	//for(int i=126;i<168;i++){ 
		
		string angle = to_string(i);
		if(i<10){
			angle = "00000"+angle;
		}
		else if(i<100){
			angle = "0000"+angle;
		}
		else{
			angle = "000"+angle;
		}
		
		cout<<i<<" of "<<768<<endl;
		
		#pragma omp parallel num_threads(numThreads)
    {
		int i_t = i%179;
		int i_t2 = i%539;
		cout<<"pre-computing water textures..."<<endl;
		string water_str = "./Water2/"+to_string(i_t)+".txt";
		string water_bump_str = "./Water2/bump"+to_string(i_t)+".txt";
		char *water_texture2 = new char[water_str.length() + 1];
		strcpy(water_texture2, water_str.c_str());
		cout<<"water_texture2 = "<<water_texture2<<endl;
			
		char* water_texture="Water.txt";
		char *water_bump2 = new char[water_bump_str.length() + 1];
		strcpy(water_bump2, water_bump_str.c_str());
		cout<<"water_bump2 = "<<water_bump2<<endl;
		char* water_bump="Water_bump5smooth.txt";
		
		//string fire_str = "./Fire/"+to_string(i_t2)+".txt";
		//char *fire_texture = new char[fire_str.length() + 1];
		//strcpy(fire_texture, fire_str.c_str());
		
		double s_r = 0.52;
		double s_g = 0.8168;
		double s_b = 0.7293;
		
		#pragma omp critical
			{
			int threadId = omp_get_thread_num();
			cout<<"thread number "<<threadId<<" of "<<numThreads<<endl;
			waters[threadId] = new Texture(water_texture2);  // Allocate and copy base texture
			waters[threadId]->scale_data(s_r,s_g,s_b);
			water_bumps[threadId] = new Texture(water_bump2);  // Allocate and copy base texture
			//fires[threadId] = new Texture(fire_texture);
			
		}
    }
		
		
		//#pragma omp parallel for private(p,l,u,w,h) num_threads(numThreads) schedule(dynamic)
		for(int j=120;j<=120;j++){
		string slash = "/";
		string folder = to_string(j)+slash;
		//folder = std::string(3 - std::min(3, folder.length()), '0') + folder+'/';
		if(j<10){
			folder = "00"+folder;
		}
		else if(j<100){
			folder = "0"+folder;
		}
		else{
			folder = folder;
		}
		//int i = 384+int(floor(-125+(ii-1)*dM));
		
	
		//if(!(250<=i && i<=386) && !(500<=i && i<=614) && !(750<=i && i<=838) && !(1000<=i && i<=1072) && !(1250<=i && i<=1550) && !(1750<=i && i<=1799) && !(2000<=i && i<=2057) 
	//&&	!(2250<=i && i<=32315) && !(2500<=i && i<=2573) && !(2750<=i && i<=2839) && !(3000<=i && i<=3100) && !(3250<=i && i<=3369) && !(3500<=i && i<=3646) && !(3750<=i && i<=3899)){
		//if((i==2639)){

		
		
			//string filename_str = filename_base+folder+filename_middle+angle+filename_extension;
			//string filename_str = filename_base+folder+angle+filename_extension;
			string filename_str = filename_base+to_string(i)+filename_extension;
			
			p[0]=3;p[1]=3;p[2]=-3;
			
			
			World earth(0);
			
			/** let's not have a Chimera array and instead imbed it's logic here.**/
			
				double R0 = 5.1962;
				//double H_m = -2.1523;
				//double H_m = -2.6;
				double H_m = -3.5225;
				//double H_m = -4.29;
				//double H_M = 2.1523;
				//double H_M = 1.73;
				double H_M = 1.0;
				//double H_M = 0.5;
				//double UP[3]={1/SQRT3,-1/SQRT3,-1/SQRT3};
				double UP[3]={0,0,1};
				double arclength = 120*pi/180;
				double holo_W=20;
				double holo_H=27;
				double hogel_size = 250;
				double zoom = 7;
			
			
				// all cameras in the array look at the same point and have the same up vector.
				l[0]=0;l[1]=0;l[2]=0;
				double up[3];
				up[0]=UP[0];up[1]=UP[1];up[2]=UP[2];
				//we align our cylinder relative to the up vector.  This will break if the up vector is parallel to e1.
				double len = sqrt(up[0]*up[0]+up[1]*up[1]+up[2]*up[2]);
				up[0]/=len;up[1]/=len;up[2]/=len;
				
				double e1[3]={1-(1*up[0])*1,0,0};
				len = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
				e1[0]/=len;e1[1]/=len;e1[2]/=len;
				double e2[3]={up[1]*e1[2]-e1[1]*up[2],e1[0]*up[2]-up[0]*e1[2],up[0]*e1[1]-e1[0]*up[1]}; 
				
				//cout<<e1<<endl;
				//cout<<e2<<endl;
				//cout<<up<<endl;
				
				// zoom is a percentage, so that a zoom of 100 does nothing.
				
				double F=2*tan(holo_W/(2*R0))/(zoom/100);
				// let's build the array of positions and location frames for each camera in the area.
				//vector<Pinhole_Camera*> column;
				//for(int i = 0 ; i<768;i++){
				//	cout<<"building Chimera array...column "<<i<<" of "<<768<<endl;
				double theta = arclength/2-i*arclength/767;
				//	column.clear();
				//	for(int j=45;j<90;j++){
				double z = H_m + (j-1)*(H_M-H_m)/179;
				p[0] = R0*cos(theta)*e1[0]+R0*sin(theta)*e2[0]+z*up[0];
				p[1] = R0*cos(theta)*e1[1]+R0*sin(theta)*e2[1]+z*up[1];
				p[2] = R0*cos(theta)*e1[2]+R0*sin(theta)*e2[2]+z*up[2];
						
				//		Pinhole_Camera* C = new Pinhole_Camera(p,l,up,F);
				//		column.push_back(C);
				//		
				//	}
				//	cameras.push_back(column);
				//}
				double holoW=holo_W;
				double holoH=holo_H;
				double hogelSize = hogel_size;

				// now let's compute the width and height of the images in pixels.
				// notice the coversion factor of 10000 because holoW is in cm but
				// hogelSize is in micrometers.
				
				// also notice the scale factor of 1.1 as 10% of each image is cut away.
				
				w = int(floor(1.1*holoW*10000/(hogel_size)));
				h = int(floor(1.1*holoH*10000/(hogel_size)));
				
				//w = 692;
				//h = 934;
				w = 2502;
				h = 750;
				//w = 1834;
				//h = 375;
			
			cout<<"w = "<<w<<" h = "<<h<<endl;
			
			cout<<"F = "<<F<<endl;
			
			/** end logic **/
			
			/** move textures here **/
			//Texture *T1, *T2,*T3,*T4,*T5,*T6;
			Texture *T1,*T2, *T4,*T5,*T6;
			Texture *TT0,*TT1,*TT2,*TT3,*TT4,*TT5,*TT6,*TT7,*TT8,*TT9;
			cout<<"making texture names"<<endl;
			#pragma omp critical
			{
			int thread_id = omp_get_thread_num();
			/**
			int i_t = i%179;
			//char* space_panorama="BlueSky.txt";
			//char* space_panorama="BlueSkyAurora.txt";
			//char* metal_texture="3metal_texture_big_100617.txt";
			//char* metal_texture="Gold.txt";
			char* bump_texture="metal_bump3.txt";
			//char* bump_texture="Gold_bump.txt";
			string water_str = "./Water2/"+to_string(i_t)+".txt";
			string water_bump_str = "./Water2/bump"+to_string(i_t)+".txt";
			char *water_texture2 = new char[water_str.length() + 1];
			strcpy(water_texture2, water_str.c_str());
			cout<<"water_texture2 = "<<water_texture2<<endl;
			
			char* water_texture="Water.txt";
			char *water_bump2 = new char[water_bump_str.length() + 1];
			strcpy(water_bump2, water_bump_str.c_str());
			cout<<"water_bump2 = "<<water_bump2<<endl;
			char* water_bump="Water_bump5smooth.txt";**/

			cout<<"making textures"<<endl;
			//Texture *T1=new Texture(marble_texture);
			//Texture *T2=new Texture(space_texture);
			//T1=new Texture(space_panorama);
			T1=skies[thread_id];
			T2=metals[thread_id];
			//T3=new Texture(bump_texture);
			//T3=metal_bumps[thread_id];
			T4=waters[thread_id];
			T5=water_bumps[thread_id];
			//T6=fires[thread_id];
			T6=gold_spirals[thread_id];
			
			TT0=map_0[thread_id];
			TT1=map_1[thread_id];
			TT2=map_2[thread_id];
			TT3=map_3[thread_id];
			TT4=map_4[thread_id];
			TT5=map_5[thread_id];
			TT6=map_6[thread_id];
			TT7=map_7[thread_id];
			TT8=map_8[thread_id];
			TT9=map_9[thread_id];
			
			
			cout<<"current thread id = "<<thread_id<<endl;
			
			//T1->save(thread_id);
			
			}
			cout<<"textures made"<<endl;
			
			/**end texture move **/
			
			
			/** move lights to here **/
			
			double pos[3]={0,0,40};
			int color[3]={255,255,255};
			Light *L=new pointLight(pos,10,color);
			L->illumination_map_ID = 2;
			L->pp=0;
			
			/**double pos[3]={5,5,5};
			double pos2[3]={0,5,3};
			int color[3]={255,255,255};
			int color2[3]={100,255,255};
			int color3[3]={255,247,158};
			//Light *L=new pointLight(pos,4.5,color);
			Light *L=new pointLight(pos,15,color);
			Light *L2=new ambientLight(0.15,color);
			//Light *L3=new pointLight(pos2,10,color);
			
			double pos3[3]={4,0,-1};
			double pos3a[3]={-4,4,-2};
			double pos4[3]={-3,1,1.5};	//double p[] = {-2,2,3};
			double pos4a[3]={-3,1,-1.5};	//double p[] = {-2,2,3};
			
			//Light* L4=new pointLight(pos3,6.5,color);
			//Light* L5=new pointLight(pos4,0.5,color);
			//Light* L5=new pointLight(pos4,3.0,color);
			
			double r0=0.65;
			double r1=1;
			
			double pos5[3]={0,0,0};
			double pos6[3]={sqrt(2.0/9)*3,sqrt(2/3.0)*3,(1.0/3)*3};  // this light should be made soft
			double e61[3]={-0.8165,0.5469,-0.1850};
			double e62[3]={-0.3333,-0.1850,0.9245};
			double pos60p[3],pos60m[3],pos6p0[3],pos6m0[3];
			for(int kk=0;kk<3;kk++){
				pos60p[kk]=pos6[kk]+r1*e61[kk];
				pos60m[kk]=pos6[kk]-r1*e61[kk];
				pos6p0[kk]=pos6[kk]+r1*e62[kk];
				pos6m0[kk]=pos6[kk]-r1*e62[kk];
			}
			double pos6a[3]={sqrt(2.0/9)*3,-sqrt(2/3.0)*3,(1.0/3)*3};  // this light should be made soft
			double e6a1[3]={0.8165,0.5469,0.1850};
			double e6a2[3]={-0.3333,0.1850,0.9245};
			double pos6a0p[3],pos6a0m[3],pos6ap0[3],pos6am0[3];
			for(int kk=0;kk<3;kk++){
				pos6a0p[kk]=pos6a[kk]+r1*e6a1[kk];
				pos6a0m[kk]=pos6a[kk]-r1*e6a1[kk];
				pos6ap0[kk]=pos6a[kk]+r1*e6a2[kk];
				pos6am0[kk]=pos6a[kk]-r1*e6a2[kk];
			}
			
			double delta = 2.5;
			
			double pos7[3]={0,0,-delta}; // this light should be made soft
			
			double pos70p[3]={0,r0,-delta}; 
			double pos70m[3]={0,-r0,-delta};
			double pos7p0[3]={r0,0,-delta}; 
			double pos7m0[3]={-r0,0,-delta};
			double pos7pp[3]={r0,r0,-delta}; 
			double pos7pm[3]={r0,-r0,-delta};
			double pos7mp[3]={-r0,r0,-delta}; 
			double pos7mm[3]={-r0,-r0,-delta};
			
			double pos8[3]={0,0,1.9}; // this light should be made soft
			double pos80p[3]={0,r0,1.9}; 
			double pos80m[3]={0,-r0,1.9};
			double pos8p0[3]={r0,0,1.9}; 
			double pos8m0[3]={-r0,0,1.9};
			double pos8pp[3]={r0,r0,1.9}; 
			double pos8pm[3]={r0,-r0,1.9};
			double pos8mp[3]={-r0,r0,1.9}; 
			double pos8mm[3]={-r0,-r0,1.9};
			
			Light* L6=new pointLight(pos5,2.5,color3);
			//L6->illumination_map_ID = 5;
			//L6->pp=0.0;
			double divizzor=4;
			Light* L7=new pointLight(pos6,5.0/divizzor,color);
			
			Light* L70p=new pointLight(pos60p,2*5.0/divizzor,color);
			Light* L70m=new pointLight(pos60m,2*5.0/divizzor,color);
			Light* L7p0=new pointLight(pos6p0,2*5.0/divizzor,color);
			Light* L7m0=new pointLight(pos6m0,2*5.0/divizzor,color);
			
			double divizor=5;
			//Light* L8=new pointLight(pos4a,5.0,color);
			Light* L9=new pointLight(pos7,15.0/divizor,color);
			Light* L90p=new pointLight(pos70p,15.0/divizor,color);
			Light* L90m=new pointLight(pos70m,15.0/divizor,color);
			Light* L9p0=new pointLight(pos7p0,15.0/divizor,color);
			Light* L9m0=new pointLight(pos7m0,15.0/divizor,color);
			Light* L9pp=new pointLight(pos7pp,15.0/divizor,color);
			Light* L9pm=new pointLight(pos7pm,15.0/divizor,color);
			Light* L9mp=new pointLight(pos7mp,15.0/divizor,color);
			Light* L9mm=new pointLight(pos7mm,15.0/divizor,color);
			Light* L14=new pointLight(pos8,15.0/6,color2);
			//L14->illumination_map_ID = 3;
			//L14->is_flat=true;
			Light* L140p=new pointLight(pos80p,15.0/6,color2);
			//L140p->illumination_map_ID = 3;
			//L140p->is_flat=true;
			Light* L140m=new pointLight(pos80m,15.0/6,color2);
			//L140m->illumination_map_ID = 3;
			//L140m->is_flat=true;
			Light* L14p0=new pointLight(pos8p0,15.0/6,color2);
			//L14p0->illumination_map_ID = 3;
			//L14p0->is_flat=true;
			Light* L14m0=new pointLight(pos8m0,15.0/6,color2);
			//L14m0->illumination_map_ID = 3;
			//L14m0->is_flat=true;
			Light* L14pp=new pointLight(pos8pp,15.0/6,color2);
			//L14pp->illumination_map_ID = 3;
			//L14pp->is_flat=true;
			Light* L14pm=new pointLight(pos8pm,15.0/6,color2);
			//L14pm->illumination_map_ID = 3;
			//L14pm->is_flat=true;
			Light* L14mp=new pointLight(pos8mp,15.0/6,color2);
			//L14mp->illumination_map_ID = 3;
			//L14mp->is_flat=true;
			Light* L14mm=new pointLight(pos8mm,15.0/6,color2);
			//L14mm->illumination_map_ID = 3;
			//L14mm->is_flat=true;
			
			
			
			//Light* L10=new pointLight(pos3a,5.0,color);
			Light* L11=new pointLight(pos6a,7.0/divizzor,color);
			
			Light* L110p=new pointLight(pos6a0p,2*7.0/divizzor,color);
			Light* L110m=new pointLight(pos6a0m,2*7.0/divizzor,color);
			Light* L11m0=new pointLight(pos6am0,2*7.0/divizzor,color);
			Light* L11p0=new pointLight(pos6ap0,2*7.0/divizzor,color);
			
			//int planeColor[3]={0,255,255};
			//Light* L12=new planeLight(2.5, 5, planeColor,T4);
	
	
			//int planeColor2[3]={255,255,255};
			//Light* L13=new planeLight(-5, 5, planeColor2,T4);**/
			
			/** end lights **/
			
			
			//ChimeraCameraArray Chimera(R0,H_m,H_M,UP,arclength,holo_W,holo_H,hogel_size,zoom);
			
			char *filename_i = new char[filename_str.length() + 1];
			strcpy(filename_i, filename_str.c_str());
			cout<<filename_i<<endl;
			//offset=1.1-i*1.2/(N-1);
			//offset=(offset-0)/2+2;
			double offset4D=1-0.25+i*0.5/(768-1);
			double offset3D=(offset4D)/2+2-0;  // subtract 1 for top, 0 for bottom.

			//double o[]={offset3D/3-0.5*x1[0]-0.5*x2[0],offset3D/3-0.5*x1[1]-0.5*x2[1],offset3D/3-0.5*x1[2]-0.5*x2[2]};
			//double o[3]={-6,-6,5};
			//double x1[3]={60,0,0};
			//double x2[3]={0,60,0};
			//Parallelogram *P1=new Parallelogram(o,x1,x2,216,222,47);
			//P1->texture_ID=3;
			
			//InfinitePlane* P0 = new InfinitePlane(-6,0,255,255);
			InfinitePlane* P0 = new InfinitePlane(-6,0,255,255);
			P0->normal_dir = -1;
			P0->texture_ID=4;
			P0->self_illuminating=true;
			//P0->tile_size=12;
			P0->tile_size=8;
			P0->upper_bound_coordinate=0;
			//P0->upper_bound=3.0;
			P0->upper_bound=2.0;
			P0->texture_center[0]=0;
			//P0->bump_ID=2;
			
			//water
			InfinitePlane* P1 = new InfinitePlane(0,0,255,255);
			P1->texture_ID=2;
			P1->bump_ID=3;
			P1->self_illuminating=true;
			P1->reflectivity=0.5;
			P1->upper_bound_coordinate=0;
			//P1->upper_bound=3.0;
			P1->upper_bound=2.0;
			//Parallelogram *P1=new Parall
			
			InfinitePlane* P2 = new InfinitePlane(0,255,0,0);
			P2->texture_ID=4;
			P2->self_illuminating=true;
			P2->alpha_ID=4;
			//P2->tiling=false;
			//P1->transparency=0.5;
			
			InfinitePlane* P3 = new InfinitePlane(4,0,255,255);
			P3->coordinate = 0;
			P3->texture_ID=0;
			P3->self_illuminating=true;
			P3->upper_bound_coordinate=0;
			P3->upper_bound=2.0;
			//P3->tile_size=6;
			//InfinitePlane* P4 = new InfinitePlane(-6,0,255,255);
			InfinitePlane* P4 = new InfinitePlane(-2,0,255,255);
			P4->normal_dir = -1;
			P4->coordinate = 0;
			P4->texture_ID=0;
			P4->tile_size=4;
			P4->self_illuminating=true;
			InfinitePlane* P5 = new InfinitePlane(4,0,255,255);
			P5->coordinate = 1;
			P5->texture_ID=0;
			P5->tile_size=4;
			P5->self_illuminating=true;
			P5->upper_bound_coordinate=0;
			//P5->upper_bound=3;
			P5->upper_bound=2;
			InfinitePlane* P6 = new InfinitePlane(-4,0,255,255);
			P6->normal_dir=-1;
			P6->coordinate = 1;
			P6->texture_ID=0;
			P6->tile_size=4;
			P6->upper_bound_coordinate=0;
			//P6->upper_bound=3;
			P6->upper_bound=2;
			P6->self_illuminating=true;
			
			double c_x,c_y,c_z;
			c_x=2;
			c_y=3;
			c_z=-3;
			double Lxx,Lyy,Lzz;
			Lxx = 0.1;
			Lyy = 2;
			Lzz = 6;
			
			Cube* Cu = new Cube(c_x, c_y, c_z, Lxx, Lyy, Lzz, 255, 0, 0);
			Cu->self_illuminating=true;
			Cube* Cu2 = new Cube(c_x, -c_y, c_z, Lxx, Lyy, Lzz, 255, 0, 0);
			Cu2->self_illuminating=true;
			
			double z_t=2;
			double z_b=-6;
			double x00 = 1.0;
			double y00 = 1.5;
			double y00_1 = 0.95;
			double rr = 0.35;
			Upright_Cylinder_Infinite* C1 = new Upright_Cylinder_Infinite(z_t,z_b,x00,y00,rr,E,200,200,200);
			C1->self_illuminating=true;
			C1->texture_ID=1;
			Upright_Cylinder_Infinite* C2 = new Upright_Cylinder_Infinite(z_t,z_b,x00,-y00,rr,E,200,200,200);
			C2->self_illuminating=true;
			C2->texture_ID=1;
			
			Upright_Cylinder_Infinite* C3 = new Upright_Cylinder_Infinite(z_t,z_b,-x00,y00_1,rr,E,200,200,200);
			C3->self_illuminating=true;
			C3->texture_ID=1;
			Upright_Cylinder_Infinite* C4 = new Upright_Cylinder_Infinite(z_t,z_b,-x00,-y00_1,rr,E,200,200,200);
			C4->self_illuminating=true;
			C4->texture_ID=1;
			//Upright_Cylinder* C2 = new Upright_Cylinder(zb,zt,-x00,y00,rr,255,255,0);
			//Upright_Cylinder* C3 = new Upright_Cylinder(zb,zt,-x00,-y00,rr,0,255,0);
			//Upright_Cylinder* C4 = new Upright_Cylinder(zb,zt,x00,-y00,rr,0,0,255);
			//C2->self_illuminating=true;

			//double l[]={0.5,0.5,0.5};
			//l[0]=0;l[1]=0;l[2]=0;

			//double u[]={0,1,0};
			//u[0]=1/SQRT3;
			//u[1]=-1/SQRT3;
			//u[2]=-1/SQRT3;
			//double u[3]={1/SQRT3,-1/SQRT3,-1/SQRT3};
			//for(int k=0;k<3;k++){
			//	p[k]-=2.5*u[k];
			//}
			
			
			
			
			/**double s=0.1*eye;
			double e1[3]={l[0]-p[0],l[1]-p[1],l[2]-p[2]};
			double e2[3]={e1[1]*u[2]-u[1]*e1[2],e1[2]*u[0]-e1[0]*u[2],e1[0]*u[1]-u[0]*e1[1]};
			double Ll=sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]);
			e2[0]/=Ll;e2[1]/=Ll;e2[2]/=Ll;
			p[0]+=e2[0]*s;p[1]+=e2[1]*s;p[2]+=e2[2]*s;
			l[0]+=e2[0]*s;l[1]+=e2[1]*s;l[2]+=e2[2]*s;**/
			
			
			//double theta = 2*i*2*pi/(360.0);
			
			/**double v[3];
			v[0]=p[1]*u[2]-p[2]*u[1];
			v[1]=p[2]*u[0]-p[0]*u[2];
			v[2]=p[0]*u[1]-p[0]*u[1];
			double LL = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
			v[0]=v[0]/LL;v[1]=v[1]/LL;v[2]=v[2]/LL;**/
			
			
			
			//earth.rotate(p,v,theta);
			//earth.rotate(l,v,theta);
			//earth.rotate(u,v,theta);
			
			double Lc = 16.0169/2;
			//double Lc = 2*10.6780/2;
			
			//double p_y = -Lc+2*Lc*i/1002.0;
			double theta_chimera = i*pi/1200;
			double r_geola = 9;//2*8.3155;
			//double p_y = -Lc+2*Lc*i/668.0;
			//Pinhole_Camera my_Camera(p,l,u,50);
			//Pinhole_Camera my_Camera(p,l,up,F);
			//double p_goala[3]={12.4732+2,p_y,-3};
			double p_chimera[3]={2+r_geola*sin(theta_chimera),r_geola*cos(theta_chimera),-3};
			//double p_goala[3]={2*8.3155+2,p_y,-3};//831.55
			//double l_goala[3]={2,p_y,-3};
			double l_chimera[3]={2,0,-3};
			//double up_goala[3] = {0,0,1};
			double up_chimera[3] = {0,0,1};
			double F_goala = 77.49;
			double F_chimera = 77.49;
			
			int w_chimera = 1200+120;
			int h_chimera = 1600+160;
			//double F_goala = 82.86;
			Pinhole_Camera my_Camera(p_chimera,l_chimera,up_chimera,F_chimera);
			//my_Camera.h_m=-3;
			//my_Camera.h_M=3;
			//my_Camera.h_c = 0;
			//Pinhole_Camera* my_Camera = new Pinhole_Camera(p,l,u,F);

			//double black_hole_pos[3]={-1.5,0,-2};
			double black_hole_pos[3]={2.0,0,-3};

			Minkowski* FlatSpace = new Minkowski(0.01);
			double error_tol=1e-6;
			//double black_hole_radius=0.4; // was 0.25
			double black_hole_radius=0.25; // was 0.25
			int dimension = 4;
			Schwarzschild* BlackHole = new Schwarzschild(black_hole_pos,error_tol,black_hole_radius,dimension);
			earth.space_time=BlackHole;
			//earth.space_time=FlatSpace;


			/**double C = offset4D;
			int type = 1;
			depth = 0;
			int depth2 = 4;
			double scale = 1;
			Menger_Sponge_4DSlice *M3= new Menger_Sponge_4DSlice(0,0,0,scale,depth2,type,E,offset4D,200,200,200);
			M3->is_wireframe=false;
			M3->texture_ID=1;
			M3->bump_ID=2;
			//M3->reflectivity=0.2;
			scale = 1.01;
			//Menger_Sponge_4DSlice *M2 = new Menger_Sponge_4DSlice(0,0,0,scale,depth,type,E,offset4D,255,255,0);
			//M2->is_wireframe=true;
			//Menger_Sponge_4DSlice *M= new Menger_Sponge_4DSlice(c1,c2,c3,scale,depth,type,E,offset,0,100,200);
			//M->transparency=0.95;
			//double shift1 = -4;
			//shift = offset;
			scale = 0.5;
			Menger_Sponge_4DSlice *M2 = new Menger_Sponge_4DSlice(0,0,0,scale,depth,type,E,offset4D,255,255,0);
			M2->reflection_only = true;
			M2->texture_ID=5;
			M2->alpha_ID=5;**/
			
			//Menger_Sponge *M2p5 = new Menger_Sponge(c1,c2,c3,scale,depth,shift,E,255,255,0);
			//M2p5->is_wireframe=true;
			//M2p5->rat = true;
			
			
			//depth = 3;
			//shift = 3*offset;
			
			//double p1[3]={0.3,0.2,0.1};
			//double p2[3]={0.1,0.8,0.1};
			//double p1[3]={shift/5,3*shift/5,shift/5};
			//double p2[3]={0,3*shift/4,shift/4};
			//Cylinder *C = new Cylinder(p1,p2,0.01,100,100,100);
			//Sphere *S1= new Sphere(p1[0],p1[1],p1[2],0.025,100,100,100);
			//Sphere *S2= new Sphere(p2[0],p2[1],p2[2],0.025,100,100,100);
			
			//Menger_Sponge *M3 = new Menger_Sponge(c1,c2,c3,scale,depth,offset3D,E,255,168,0);
			//M3->is_wireframe=true;
			//depth = 3;
			//Menger_Sponge_4DSlice *M2= new Menger_Sponge_4DSlice(c1,c2,c3,scale,depth,type,E,offset,200,200,200);
			///Menger_Sponge *M2= new Menger_Sponge(c1,c2,c3,scale,depth,type,E,200,200,200);
			//Menger_Sponge_4DSlice_Type2 *M2= new Menger_Sponge_4DSlice_Type2(c1,c2,c3,scale,depth,type,E,offset,200,200,200);
			
			 //forget about the 0, it's related to some invisibility cloak stuff
			earth.env_map_ID=-1; //the texture ID of the background image - explained later
			//earth.env_map_ID=3; //the texture ID of the background image - explained later
		
			earth.w=w; //this is sloppy, but I have the world know how big the image is going to be
				  //so it can resize the background image accordingly.    
			earth.h=h;

			//earth.addShape(M);
			//earth.addShape(M2);
			//earth.addShape(M2p5);

			//earth.addShape(M3);
			earth.addShape(P1);
			earth.addShape(P0);
			earth.addShape(P3);
			earth.addShape(P4);
			earth.addShape(P5);
			earth.addShape(P6);
			//earth.addShape(P2);
			earth.addShape(C1);
			earth.addShape(C2);
			earth.addShape(C3);
			earth.addShape(C4);
			//earth.addShape(Cu);
			//earth.addShape(Cu2);
			earth.addLight(L);
			//earth.addShape(S1);
			//earth.addShape(S2);
			//earth.addLight(L);
			//earth.addLight(L2);
			
			//earth.addLight(L3);
			//earth.addLight(L4);
			//earth.addLight(L5);
			//earth.addLight(L6);
			//earth.addLight(L6);
			//earth.addLight(L7);
			//earth.addLight(L70p);
			//earth.addLight(L70m);
			//earth.addLight(L7p0);
			//earth.addLight(L7m0);
			
			//earth.addLight(L8);
			/**earth.addLight(L9);
			earth.addLight(L90p);
			earth.addLight(L90m);
			earth.addLight(L9p0);
			earth.addLight(L9m0);
			earth.addLight(L9pp);
			earth.addLight(L9pm);
			earth.addLight(L9mp);
			earth.addLight(L9mm);
			
			earth.addLight(L14);
			earth.addLight(L140p);
			earth.addLight(L140m);
			earth.addLight(L14p0);
			earth.addLight(L14m0);
			earth.addLight(L14mp);
			earth.addLight(L14pm);
			earth.addLight(L14pp);
			earth.addLight(L14mm);
			//earth.addLight(L10);
			earth.addLight(L11);
			earth.addLight(L110p);
			earth.addLight(L110m);
			earth.addLight(L11p0);
			earth.addLight(L11m0);**/
			//earth.addLight(L12);
			//earth.addLight(L13);
			earth.addTexture(T1); //this will have texture ID 0, since it is added first
			earth.addTexture(T2); //this will have texture ID 1, since it is added second
			//earth.addTexture(T3); //this will have texture ID 2, since it is added third
			earth.addTexture(T4); //this will have texture ID 3, since it is added fourth
			earth.addTexture(T5); //this will have texture ID 4, since it is added fifth
			earth.addTexture(T6); //this will have texture ID 5, since it is added sixth
			
			earth.addPhotonMapTexture(TT0);
			earth.addPhotonMapTexture(TT1);
			earth.addPhotonMapTexture(TT2);
			earth.addPhotonMapTexture(TT3);
			earth.addPhotonMapTexture(TT4);
			earth.addPhotonMapTexture(TT5);
			earth.addPhotonMapTexture(TT6);
			earth.addPhotonMapTexture(TT7);
			earth.addPhotonMapTexture(TT8);
			earth.addPhotonMapTexture(TT9);
			
			bool full_column = false;
			double theta_SLM = 120*pi/180;
			//earth.renderChimeraColumn(&Chimera,filename_base, i,full_column);
			//cout<<"frame "<<i<<" about to render!"<<endl;
			//earth.render(&my_Camera,filename_i,w_chimera,h_chimera);
			double Lz = 6.0;
			double Ly = 4.5;
			//int N_b = 134;
			int N_b = 225;
			int a = i;
			earth.render_hogel(&my_Camera, filename_i, a, N_a, N_b, N_I, black_hole_pos, Ly, Lz, r_geola, theta_SLM);
			//cout<<"frame "<<i<<" rendered!"<<endl;
			//free(M2);
			//free(M3);
			//free(L6);
			//free(P1);
			//delete T1;
			//delete T2;
			//delete T3;
			//delete T4;
			//delete T5;
			
			//delete M2;
			//delete M3;
			delete P1;
			delete L;
			/**delete L2;
			delete L6;
			delete L7;
			delete L70p;
			delete L70m;
			delete L7p0;
			delete L7m0;
			
			delete L9;
			delete L90p;
			delete L90m;
			delete L9p0;
			delete L9m0;
			delete L9pp;
			delete L9pm;
			delete L9mp;
			delete L9mm;
			
			delete L14;
			delete L140p;
			delete L140m;
			delete L14p0;
			delete L14m0;
			delete L14mp;
			delete L14pm;
			delete L14pp;
			delete L14mm;
			
			delete L11;
			delete L110p;
			delete L110m;
			delete L11p0;
			delete L11m0;**/
		//}
		
	}
	//clear the texture vectors
	for(int ii=0;ii<numThreads;ii++){
		delete waters[ii];
		delete water_bumps[ii];
		//delete fires[ii];
	}
	waters.clear();
	water_bumps.clear();
	//fires.clear();
	}
	//}

}
catch (std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception caught" << std::endl;
    }

    std::cout << "Program completed successfully." << std::endl;
}