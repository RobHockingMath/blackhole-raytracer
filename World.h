#include <vector>
#include "Sphere.h"
#include "Camera.h"
#include "Chimera.h"
#include "Light.h"
#include "PhotonMapLight.h"
#include "Texture.h"
#include <string>
#include <cstring>
#include "Spacetime.h"
#include "photon.hpp"
#include <random>

#include <opencv2/opencv.hpp>

class World {
public:
	World(double);
	~World();
	//void writeData(char*);
	void init();
	void addShape(Shape*);
	void addLight(Light*);
	void addPhotonMapLight(PhotonMapLight*);
	void addTexture(Texture*);
	void addPhotonMapTexture(Texture*);
	//void render_hogel(Camera* cam, char* filename, int a, int N_a, int N_b, int N_I, double h_c[], double Ly, double Lz, double r, double theta_SLM);
	//void render_full_parallax_hogel(Camera* cam, char* filename, int a, int b, int N_a, int N_b, int N_I, int N_J, double h_c[], double Ly, double Lz, double r, double theta_SLM_h,double theta_SLM_v);
	void encode(double x, int &R, int &G, int &B);
	void render_spherically_symmetric_spacetime(Camera* cam,double o[], char* filename,char *time_name,int w, int h,double P_obs_loc[],double V_obs_loc[]);
	void spherically_symmetric_dimmness_image(Camera* cam, char* filename,char *time_name,int w, int h,double P_obs_loc[],double V_obs_loc[]);
	void render(Camera*,char*,int,int);
	void make_photon_map(int w, int h);
	bool linear_ray_hits_something(int ii,int jj,double P[], double V[]);
	void ray_trace(int ii,int jj,double P[], double V[], int C[],int depth,double x[]);
	void non_linear_photon_emission(double P[],double V[],vector< vector<Photon>>& photons,double photon_energy,int depth,double H,double E,int last_index,bool caustic);
	double get_far_time(double P[], double V[],double H,double E,double &phi_far);
	void spherically_symmetric_fartime_table(double P_loc[], char* filename);
	void non_linear_ray_trace(int ii,int jj,double P[], double V[], int C[],int depth,double H,double E,double x_ij[],int last_index,double o[], double e1[],double e2[],double e3[],double psi);
	void read_env_map(int,int,int[]);
	void read_sphere_env_map(double V[],int C[]);
	void shade(double[],double[],double[],int[]);
	void rotate(double x[], double v[], double theta, double o[]);
	void render_full_parallax_hogel_tilted(char* filename,char* filename2,char* filename3, int a, int b, int N_a, int N_b, int N_I, int N_J, double distance2, double h_c[],double o[],double e0[],double e1[],double e2[], double Ly, double Lz, double r, double theta_SLM_h,double theta_SLM_v);
	void inpaint_hogel_image(int*** hogel_data, int*** hogel_mask, int N_I, int N_J, const char* original_filename, const char* mask_filename, const char* inpainted_filename);
	std::vector<Shape*> shapes;
	Spacetime *space_time;
	int env_map_ID;
	int w;
	int h;
private:
	//0;
	
	static const bool shadows_on=false;
	static const bool soft_shadows_on=false;
	static const int MAX_DEPTH=5;//number of reflections per ray
	double MAX_DIST=250;

	double invis_R;

	std::uniform_real_distribution<> distrib;
	std::mt19937 gen;

	//int image[w][h][3];

	//int ***image;
	int ***im_data;
	int ***time_data;
	int ***hogel_data;
	int ***hogel_mask;

	
	std::vector<Light*> lights;
	std::vector<PhotonMapLight*> photon_lights;
	std::vector<Texture*> textures;
	std::vector<Texture*> photon_maps;
};

