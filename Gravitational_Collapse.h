#ifndef GRAVITATIONAL_COLLAPSE_H
#define GRAVITATIONAL_COLLAPSE_H

#include "Spacetime.h"

#include <vector>
#include <string>
#include <array>
//#include <map>

using namespace std;

/**

	This class models, in uv double null coordinates, the spacetime of a spherically symmetric massless scalar field.
	lines of constant u are called columns and indexed by j, whereas lines of constant v are called rows and indexed by i.
	The spacetime is constructed by numerical solution of the Einstein Equations, as described in the paper
	"The spherically symmetric collapse of a massless scalar field" by Hamade and Stewart, 1995.
	unlike them, we use a limited form of adaptive mesh refinement where cell sizes are constant in the v direction, but
	are allowed to get smaller as u increases.  This is done in such a way that the cells are always square, in order to ensure
	stability of the underlining discretization.
	
	The numerical solution to the Einstein Equations is the spacetime metric (a,r), it's derivatives (c,d,f,g), the scalar field s, 
	and it's derivatives (p,q) on a grid within the set U_min<=u<=v<=U_max.	 This solution is extended via linear interpolation to the whole set 
	U_min<=u<=v<=U_max.  Such an extension is necessary because we want to be able to have observer and light ray objects interact with our spacetime,
	and this requires the above functions to be defined everywhere.  The functions "get_a(u,v)", etc return these functions extended by linear interpolation.

**/

class Gravitational_Collapse : public Spacetime{

public:

	// the imploding scalar field is a gaussian with amplitude Psi0, center r0, width delta.
	// Depending on what these parameters are, a black hole may or may not form.
	// n controls the number of grid cells in the first column.
	// evolution proceeds until either the full set U_min<=u<=v<=U_max has been filled (in the case where no black hole forms)
	// or else if a black hole forms and grid cells are getting smaller without bound, then evolution terminates once the ratio of the original meshsize
	// to the current meshsize exceeds R_max.

	Gravitational_Collapse(double Psi0,double r0,double delta,double U_min,double U_max,int n,double R_max,double tol_);
	
	bool readSpacetimeFromDisk(string filename);
	void load_path(string filename);
	
	/**2024 new routines **/
	void get_null_vector_in_direction(double e1[],double e2[],double e3[],double P_obs_loc[],double V_obs_loc[],double V[],double V_loc[],double &psi);
	void cartesian_to_local(double e1[],double e2[], double e3[], double p_cart[],double v_cart[],double p_local[],double v_local[],double &psi);
	void local_to_cartesian(double o[],double e1[],double e2[], double e3[], double p_local[],double p_cart[],double psi);
	
	
	/** These routines are related to outputting data for visualization **/
	
	// since the solution is not on a uniform mesh, visualizing it is difficult in standard applications like Matlab.
	// Instead I output the solution as a .obj file, so that I can look at it in 3dsmax.
	void write_obj(string id,string rescale);
	void write_objs();
	
	// these are helper routines, for normalizing the output data if desired.
	
	double max_abs(vector<double*> X);
	double get_max_abs(string id);
	
	// write the solution on axis to a txt file, for visualization in matlab.
	void axis_dump();
	void dump_row_or_column(string type,int k);
	void write_stats_to_disk();
	
	/** These routines are related to the numerical solution of the Einstein Equations **/
	
	void Solve_Einstein_Equations();
	bool stop_evolution();
	
	void first_column();
	void add_column(); //adds a new column to the solution of the Einstein equations, i.e. a slice with u constant and v from u to V_max.
	
	//we decide the value of h in the new row based on information in the current row ... there are more routines related to this in there own section
	double get_new_h(); 
	void Boundary_Conditions(int i);
	void get_qhat_and_dhat(double &qhat,double &dhat,int i,int j);
	void Populate(int i,int j);
	double get_phat(int i,int j,double ghat,double qhat,double fhat,double rhat);
	double get_fhat(int i,int j,double ghat,double ahat,double rhat);	
	double get_shat(int i,int j,double qhat);	
	double get_rhat(int i,int j,double ghat);	
	double get_ghat(int i,int j,double qhat,double dhat);
	double get_ahat(int i,int j,double dhat);
	double get_chat(int i,int j);

	/** These routines are related to interpolation and extrapolation of the numerical solution **/
	
	// generally I have overloaded the functions to accept either a pt (u,v) in spacetime or else
	// a column index i in place of u.
	double extrapolate(vector<double*> X,int i,int j);
	double interpolate(vector<double*> X,double u,double v);
	double interpolate(vector<double*> X,int i,double v);
	double get(double u,double v,string id);
	double get(int i,double v,string id); //id is one of "q","d","c","a","g","r","s","f","p";
	int find_i(double u);
	
	// get_q(u,v) and get(u,v,"q") do the same thing...I just found it convenient to also be able to call the function by either name
	double get_q(double u,double v);
	double get_d(double u,double v);
	double get_a(double u,double v);
	double get_g(double u,double v);
	double get_r(double u,double v);
	double get_s(double u,double v);
	double get_f(double u,double v);
	double get_p(double u,double v);
	double get_c(double u,double v);
	double get_a_u(double u,double v);
	double get_q(int i,double v);
	double get_d(int i,double v);
	double get_a(int i,double v);
	double get_g(int i,double v);
	double get_r(int i,double v);
	double get_s(int i,double v);
	double get_f(int i,double v);
	double get_p(int i,double v);
	double get_c(int i,double v);
	
	void get_acdgfr(double u, double v, double &A, double &C, double &D, double &G, double &F, double &R);
	
	/** These routines are related to determining what size the mesh spacing should be in the next column, and extracting
		the local meshsize later **/
	
	double get_isocontour_projection(double u,double v);
	double get_isocontour_projection(int i,double v);
	double minimum_projection(int i);
	double get_local_mesh_size(double u,double v);
	double min_step_size(double* y);
	/** These routines are for determining the trajectories of light rays and observers through the spacetime after it has
	been constructed.  **/
	double get_hawking_mass(double u,double v);
	bool inside_universe(double u,double v);
	void accel(double &udd,double &vdd,double &phidd,double u,double v,double phi,double udot,double vdot,double phidot,double phi0); // the local acceleration of a freely falling particle
	bool accel_1st_order(double*& yp,double* y, double E); // acceleration as a 1st order system - returns false if we attempt to compute at a point outside of the universe
	bool constant_r_flow(double*& yp,double* y); // same as constant_r_flow_vector but formatted in a way that adaptive RK will like.
	void constant_r_flow_vector(double u,double v,double &ud,double &vd); //
	void constant_r_trajectory(double u0,double v0,double tau_max,vector<double>& u_path,vector<double>& v_path,vector<double>& r_path, vector<double>& tau_path, vector<double>& m_path);
	bool Adaptive_Step_Constant_R(double* y,bool& y_changed,int n,double& h,double tol,bool force_step,double r_desired);
	double error_norm(double *p,double *x); 
	
	//void buildHashTable();
	
	/** Now everything else is member variables **/
	
	
	/** 
		column_count is literally the number of columns...NOT the index of the last column, which is
		column_count-1.
	**/
	int column_count;

	/** The solution to the Einstein equations, a system of 9 PDEs in the variables q,d,c,a,g,r,s,f,p
	are stored in 9 2D arrays of the same name
	the first index i corresponds to u, while j corresponds to v.  As i increases, the number of grid cells in a 
	'column' of constant i increases, due to adaptive mesh refinement.  Thus the total number of columns is unknown at the
	outset, which is the reason for the use of a vector data structure
	**/
	
	vector<double*> q,d,c,a,g,r,s,f,p; 
	
	vector<double> h; // h[i] is the step size in column i
					  // *** and also the space between column i-1 and i ***
	vector<int> N;    // N[i] is the number of data points in column i
	vector<double> U; //  U[i] is the value of u in the column i 
	//map<double,int> lookup_table;

	double u_min,u_max,v_min,v_max;
	double Psi0,r0,delta; // parameters for the massless scalar field

	bool contains_blackhole;
	double u_horizon; // u coordinate of the event horizon, if present
	double m_bh;
	
	double R_max; // maximum allowable ratio of smallest to largest grid size - evolution will terminate before cells smaller than this are made
	vector<array<double, 5>> reference_data; //[tau, u_smooth, v_smooth, u_dot_smooth,v_dot_smooth];
	int reference_data_length;
};

#endif