#include "Gravitational_Collapse.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#define pi 3.14159265359

using namespace std;

	Gravitational_Collapse::Gravitational_Collapse(double Psi0_,double r0_,double delta_,double U_min,double U_max,int n,double R_MAX,double tol_)
		: Spacetime(tol_), Psi0(Psi0_), r0(r0_), delta(delta_), R_max(R_MAX), u_min(U_min), u_max(U_max)
		{

		// strictly speaking we always have v_min=u_min, v_max=u_max, but I found it makes formulas easier to read if we have these
		// redundant variables
		
		v_min=U_min;
		v_max=U_max;
		
		double h0=(v_max-v_min)/(n-1);  // the initial step size
		h.push_back(h0);
		N.push_back(n);
		U.push_back(u_min);
		column_count=0;
		
		contains_blackhole=false;
		u_horizon=u_max; // initial value - will be updated if a black hole forms
		is_collapse=true;
	}
	
	bool Gravitational_Collapse::readSpacetimeFromDisk(string filename)
	{
			std::ifstream in(filename, std::ios::binary);
		if(!in) {
			// handle error if you want
			return false;
		}

		// 1) Read metadata in same order
		in.read(reinterpret_cast<char*>(&u_min), sizeof(u_min));
		in.read(reinterpret_cast<char*>(&u_max), sizeof(u_max));
		in.read(reinterpret_cast<char*>(&v_min), sizeof(v_min));
		in.read(reinterpret_cast<char*>(&v_max), sizeof(v_max));

		in.read(reinterpret_cast<char*>(&Psi0),  sizeof(Psi0));
		in.read(reinterpret_cast<char*>(&r0),    sizeof(r0));
		in.read(reinterpret_cast<char*>(&delta), sizeof(delta));

		cout<<"u_min, u_max, v_min, v_max= "<<u_min<<" "<<u_max<<" "<<v_min<<" "<<v_max<<endl;
		cout<<"Psi0, r0, delta"<<Psi0<<" "<<r0<<" "<<delta<<endl;

		{
			char b;
			in.read(&b, 1);
			contains_blackhole = (b != 0);
		}

		in.read(reinterpret_cast<char*>(&u_horizon), sizeof(u_horizon));
		in.read(reinterpret_cast<char*>(&R_max),     sizeof(R_max));

		// 2) Read column_count
		in.read(reinterpret_cast<char*>(&column_count), sizeof(column_count));

		// 3) Read N
		N.resize(column_count);
		in.read(reinterpret_cast<char*>(N.data()), column_count * sizeof(int));

		// 4) Read U
		U.resize(column_count);
		in.read(reinterpret_cast<char*>(U.data()), column_count * sizeof(double));

		// 5) Read h
		h.resize(column_count);
		in.read(reinterpret_cast<char*>(h.data()), column_count * sizeof(double));

		// 6) Prepare pointer columns
		q.resize(column_count);
		d.resize(column_count);
		c.resize(column_count);
		a.resize(column_count);
		g.resize(column_count);
		r.resize(column_count);
		s.resize(column_count);
		f.resize(column_count);
		p.resize(column_count);
		// etc. for all your columns

		std::vector<std::vector<double*>*> columns = { &q,&d,&c,&a,&g,&r,&s,&f,&p};

		// 7) For i in [0..column_count-1], read each pointer array
		for(int i = 0; i < column_count; ++i) {
			int length = N[i];
			for(auto colPtr : columns) {
				double* arr = new double[length];
				in.read(reinterpret_cast<char*>(arr), length * sizeof(double));
				(*colPtr)[i] = arr;
			}
		}

		in.close();
		return true;
	}
	
	
	void Gravitational_Collapse::load_path(string filename){
		std::ifstream infile(filename);
		
		if (!infile.is_open()) {
			std::cerr << "Error: Could not open file " << filename << std::endl;
			return;
		}
		
		std::string line;
		
		while (std::getline(infile, line)) {
			// Skip empty lines if they exist
			if (line.empty()) {
				continue;
			}
			
			std::stringstream ss(line);
			
			// We only need the first three columns
			double col1, col2, col3, col4, col5;
			
			// Attempt to extract the first three columns
			if (ss >> col1 >> col2 >> col3 >> col4 >> col5) {
				// Push those three values into the vector
				reference_data.push_back({col1, col2, col3, col4, col5});
				
				// If there are more columns, we just ignore them
				// because we don't read them into any variables.
			}
			// Otherwise, if the line doesn't have at least three columns, we skip it.
		}
		
		infile.close();
		
		// Check how many rows we successfully read
		std::cout << "Read " << reference_data.size() << " rows of 3 columns each.\n";
		
		// (Optional) Print out the first few rows
		for (size_t i = 0; i < reference_data.size() && i < 5; ++i) {
			std::cout << reference_data[i][0] << " " << reference_data[i][1] << " " << reference_data[i][2] << " " << reference_data[i][3] << " " << reference_data[i][4] << "\n";
		}
		reference_data_length=reference_data.size();
		return;
	}
	
	double Gravitational_Collapse::min_step_size(double* y){
	
		if(y[0]>1e-3)
			return 1;//0.5*get_local_mesh_size(y[0],y[1]);
		else{
			return 10;
		}
	
	}

	// new function 2024 
	// the cartesian vector V specifies the direction.  V_loc will give us the vector.
	void Gravitational_Collapse::get_null_vector_in_direction(double e1[],double e2[],double e3[],double P_obs_loc[],double V_obs_loc[],double V[],double V_loc[],double &psi){
	
		//cout<<"getting null vector"<<endl;
	
		double v_dot_e1 = V[0]*e1[0]+V[1]*e1[1]+V[2]*e1[2];
		double v_dot_e2 = V[0]*e2[0]+V[1]*e2[1]+V[2]*e2[2];
		double v_dot_e3 = V[0]*e3[0]+V[1]*e3[1]+V[2]*e3[2];
		
		psi =atan2(v_dot_e3,v_dot_e2);
		double epsi[3] = {cos(psi)*e2[0]+sin(psi)*e3[0],cos(psi)*e2[1]+sin(psi)*e3[1],cos(psi)*e2[2]+sin(psi)*e3[2]};
		double v_dot_epsi = V[0]*epsi[0]+V[1]*epsi[1]+V[2]*epsi[2];
		double phi = atan2(v_dot_epsi,v_dot_e1);
		
		this->phi=phi;
		
		
		/**double a = get_a(P_loc[0],P_loc[1]);
		a = 1;
		double R = get_r(P_loc[0],P_loc[1]);
		R = (P_loc[1]-P_loc[0])/2.0;
	
		//cout<<"psi, phi, a, r = "<<psi<<" "<<phi<<" "<<a<<" "<<R<<endl;
	
		V_loc[0]=(-1.0+cos(phi))/a;
		V_loc[1]=(-1.0-cos(phi))/a;
		V_loc[2]=sin(phi)/R;**/
		
		double t_dir = -1;
		
		double a = get_a(P_obs_loc[0],P_obs_loc[1]);
		double R = get_r(P_obs_loc[0],P_obs_loc[1]);
		
		V_loc[0]=t_dir*V_obs_loc[0]+sqrt(V_obs_loc[0]/V_obs_loc[1])*cos(phi)/a;
		V_loc[1]=t_dir*V_obs_loc[1]-sqrt(V_obs_loc[1]/V_obs_loc[0])*cos(phi)/a;
		V_loc[2]=sin(phi)/R;
		
		
		//V_loc[3]=phi; // useful to store this for later.
	
	}

	// new function 2024
	void Gravitational_Collapse::cartesian_to_local(double e1[],double e2[], double e3[], double p_cart[],double v_cart[],double p_local[],double v_local[],double &psi){
		// we've simulatenously converting point (x,y,z) and vector cartesian coordinates into the manifold coordinate system (u,v,phi)
		// we exploit spherical symmetry by orienting our coordinate system so that this is always an equatorial orbit in the local manifold coordinates.
		// Notice that we assume we have temporal knowledge of our position, but not of our velocity (as the latter can be deduced from the null condition).
		
		// we assume the components of v_cart are expressed in global cartesian coordinates, 
		// rather than with respect to the local basis e1, e2, e3.  We assume e1 points at the black hole.
		// that means e1=(1/a)*e_u-(1/a)*e_v
		
		//double p_dot_e1 = p_cart[0]*e1[0]+p_cart[1]*e1[1]+p_cart[2]*e1[2];
		//double p_dot_e2 = p_cart[0]*e2[0]+p_cart[1]*e2[1]+p_cart[2]*e2[2];
		//double p_dot_e3 = p_cart[0]*e3[0]+p_cart[1]*e3[1]+p_cart[2]*e3[2];
		
		double v_cart_dot_e1 = v_cart[0]*e1[0]+v_cart[1]*e1[1]+v_cart[2]*e1[2];
		double v_cart_dot_e2 = v_cart[0]*e2[0]+v_cart[1]*e2[1]+v_cart[2]*e2[2];
		double v_cart_dot_e3 = v_cart[0]*e3[0]+v_cart[1]*e3[1]+v_cart[2]*e3[2];
		
		
		// first we solve for u & v from our radius and time.
		double t = 0;
		double R = sqrt(p_cart[0]*p_cart[0]+p_cart[1]*p_cart[1]+p_cart[2]*p_cart[2]);
		double u = t-R;
		double v = t+R;
		
		// we find the angle necessary to make this all happen in the equatorial plane in local coordinates.
		
		psi = atan2(v_cart_dot_e3,v_cart_dot_e2);		
		double epsi[3] = {cos(psi)*e2[0]+sin(psi)*e3[0],cos(psi)*e2[1]+sin(psi)*e3[1],cos(psi)*e2[2]+sin(psi)*e3[2]};
		
		double phi = 0; 
		p_local[0]=u;
		p_local[1]=v;
		p_local[2]=phi;
		
		// now we move onto the vector bit.
		
		// we assume that the three components of v_cart are v.e1+v.e2+v.e3.
		// we then solve for v.e0 using v.e1+v.e2+v.e3+v.e0=0;
		
		double v_cart_dot_e0 = -v_cart_dot_e1-v_cart_dot_e2-v_cart_dot_e3;
		
		// we find the angle necessary to make this all happen in the equatorial plane in local coordinates.
		
		
		double v_cart_dot_epsi = v_cart[0]*epsi[0]+v_cart[1]*epsi[1]+v_cart[2]*epsi[2];
		
		double a = this->get(u,v,"a");
		double r = this->get(u,v,"r");
		
		double v_dot_eu = a*(v_cart_dot_e1+v_cart_dot_e0);
		double v_dot_ev = a*(v_cart_dot_e1-v_cart_dot_e0);
		double v_dot_ephi = r*(v_cart_dot_epsi);
		
		v_local[0]=v_dot_eu;
		v_local[1]=v_dot_ev;
		v_local[2]=v_dot_ephi;
	
	}
	
	// new function 2024
	void Gravitational_Collapse::local_to_cartesian(double o[],double e1[],double e2[], double e3[], double p_local[],double p_cart[],double psi){
		// t = (v+u)/2
		// r = (v-u)/2
		double R;
		if(p_local[0]>0){
			R = get_r(p_local[0],p_local[1]);
		}
		else{
			R = (p_local[1]-p_local[0])/2.0; // assume flat space.
		}
	
		double weight = (1+tanh(r0-R))/2;
		weight = 0;
		R = weight*R+(1-weight)*(p_local[1]-p_local[0])/2.0;
		double phi = p_local[2];
		double x = R*cos(phi);
		double y = R*sin(phi);
		
		
		double epsi[3] = {cos(psi)*e2[0]+sin(psi)*e3[0],cos(psi)*e2[1]+sin(psi)*e3[1],cos(psi)*e2[2]+sin(psi)*e3[2]};
		
		p_cart[0]=o[0]-x*e1[0]+y*epsi[0];
		p_cart[1]=o[1]-x*e1[1]+y*epsi[1];
		p_cart[2]=o[2]-x*e1[2]+y*epsi[2];
		
	}

	
	
	
	// if rescale=="rescale", normalize maximum height to 1 before saving to .obj
	// id is one of "a","r","c","d","f","g","s"
	void Gravitational_Collapse::write_obj(string id,string rescale){

		string filename;
		FILE *fp;
		
		double S=1; //normalization factor
		if(rescale.compare("rescale")==0){
			S=get_max_abs(id);
			filename="obj_files/"+id+"_rescaled.obj";
		}
		else{
			filename="obj_files/"+id+".obj";
		}
		fp = fopen(filename.c_str(),"w");
	
		for(int i=0;i<column_count-2;i++){
		
			//cout<<"we are on column i="<<i<<".  U[i]="<<U[i]<<" "<<"U[i+1]="<<U[i+1]<<" h[i]="<<h[i]<<endl;
		
			// the first element is a triangle.
			fprintf(fp,"v %f %f %f \n",U[i],get(U[i],U[i],id.c_str())/S,U[i]);
			fprintf(fp,"v %f %f %f \n",U[i],get(U[i],U[i+1],id.c_str())/S,U[i+1]);
			fprintf(fp,"v %f %f %f \n",U[i+1],get(U[i+1],U[i+1],id.c_str())/S,U[i+1]);
			// the remaining elements are squares.
			for(int j=1;j<N[i]-1;j++){
				fprintf(fp,"v %f %f %f \n",U[i],get(U[i],U[i]+j*h[i],id.c_str())/S,U[i]+j*h[i]);
				fprintf(fp,"v %f %f %f \n",U[i],get(U[i],U[i]+(j+1)*h[i],id.c_str())/S,U[i]+(j+1)*h[i]);
				fprintf(fp,"v %f %f %f \n",U[i]+h[i+1],get(U[i]+h[i+1],U[i]+(j+1)*h[i],id.c_str())/S,U[i]+(j+1)*h[i]);
				fprintf(fp,"v %f %f %f \n",U[i]+h[i+1],get(U[i]+h[i+1],U[i]+j*h[i],id.c_str())/S,U[i]+j*h[i]);
			}
		
		}
		int count=1;
		for(int i=0;i<column_count-2;i++){
		
			// the first element is a triangle.
			fprintf(fp,"f %i %i %i \n",count,count+1,count+2);
			count+=3;
			// the remaining elements are squares.
			for(int j=1;j<N[i]-1;j++){
				fprintf(fp,"f %i %i %i %i \n",count,count+1,count+2,count+3);
				count+=4;
			}
		
		}			
		fclose(fp);		
	
	}	
	
	void Gravitational_Collapse::write_objs(){	
		write_obj("q","no");
		write_obj("d","no");
		write_obj("a","no");
		write_obj("g","no");
		write_obj("r","no");
		write_obj("s","no");
		write_obj("f","no");
		write_obj("p","no");	
	}

	// this is just a helper routine to get the maximum
	// absolute value of a 2D array X
	double Gravitational_Collapse::max_abs(vector<double*> X){
	
		double max_so_far=-1;
		bool max_so_far_set=false;
	
		for(int i=0;i<column_count;i++){
			for(int j=0;j<N[i];j++){
				if(!max_so_far_set || max_so_far<fabs(X[i][j])){max_so_far=fabs(X[i][j]); max_so_far_set=true;}
			}
		}
		return max_so_far;
	}		
	
	double Gravitational_Collapse::get_max_abs(string id){
		if(id.compare("q")==0){
			return max_abs(q);
		}
		else if(id.compare("d")==0){
			return max_abs(d);
		}
		else if(id.compare("a")==0){
			return max_abs(a);
		}
		else if(id.compare("g")==0){
			return max_abs(g);
		}
		else if(id.compare("r")==0){
			return max_abs(r);
		}
		else if(id.compare("s")==0){
			return max_abs(s);
		}
		else if(id.compare("f")==0){
			return max_abs(f);
		}
		else if(id.compare("p")==0){
			return max_abs(p);
		}
		else if(id.compare("c")==0){
			return max_abs(c);
		}
		else{
			cout<<"variable name not found"<<endl;
			return 0;
		}
	}
	
	void Gravitational_Collapse::axis_dump(){

		string filename="Matlab/axis_dump.txt";

		FILE *fp;
		fp=fopen(filename.c_str(),"w");
	
		for(int i=0;i<column_count-3;i++){
			//cout<<i<<" "<<U[i]<<" "<<d[i][0]<<" "<<q[i][0]<<" "<<a[i][0]<<" "<<g[i][0]<<" "<<r[i][0]<<" "<<s[i][0]<<" "<<f[i][0]<<" "<<p[i][0]<<endl;
			fprintf(fp,"%f %f %f %f %f %f %f %f %f \n",U[i],d[i][0],q[i][0],a[i][0],g[i][0],r[i][0],s[i][0],f[i][0],p[i][0]);	
		}
		fclose(fp);
	}
	
	void Gravitational_Collapse::dump_row_or_column(string type,int k){

		// first make the filename...involves converting k to a string...which is annoying in C++
	
		stringstream k_as_string;
		k_as_string<<k;
	
		string filename="Matlab/"+type+"_"+k_as_string.str()+".txt";

		FILE *fp;
		fp=fopen(filename.c_str(),"w");
	
		if(type.compare("row")==0){ // dump row j=k
	
			for(int i=0;i<column_count-3;i++){
				//cout<<i<<" "<<U[i]<<" "<<d[i][N[i]-1]<<" "<<q[i][N[i]-1]<<" "<<a[i][N[i]-1]<<" "<<g[i][N[i]-1]<<" "<<r[i][N[i]-1]<<" "<<s[i][N[i]-1]<<" "<<f[i][N[i]-1]<<" "<<p[i][N[i]-1]<<endl;
				fprintf(fp,"%f %f %f %f %f %f %f %f %f \n",U[i],d[i][k],q[i][k],a[i][k],g[i][k],r[i][k],s[i][k],f[i][k],p[i][k]);	
			}
		}
		else{ // dump column j=k

			for(int j=0;j<N[k];j++){
				fprintf(fp,"%f %f %f %f %f %f %f %f %f \n",U[k]+j*h[k],d[k][j],q[k][j],a[k][j],g[k][j],r[k][j],s[k][j],f[k][j],p[k][j]);
			}				
		}
		
		fclose(fp);
	}

	//currently this just writes u_max and u_horizon to disk
	void Gravitational_Collapse::write_stats_to_disk(){
		FILE *fp;
		fp = fopen("Matlab/Manifold_Stats.txt","w");
		fprintf(fp,"%f %f \n",u_max,u_horizon);
		fclose(fp);
	}	

	void Gravitational_Collapse::Solve_Einstein_Equations(){
	
		cout<<"Solving Einstein Equations..."<<endl;
	
		first_column(); // The Initial Conditions
	
		while(!stop_evolution()){
			cout<<"i=="<<column_count<<endl;
			add_column();
		}
		
		cout<<"Done Solving Einstein Equations"<<endl;
		
		//cout<<"Building Hash table for U"<<endl;
		//buildHashTable();
		//cout<<"Hash table built"<<endl;
	
	}
	
	// stopping conditions for the evolution.  I also check here to see if a blackhole has formed.
	
	bool Gravitational_Collapse::stop_evolution(){

		/** for technical reasons, we cannot evolve all the way up to u=u_max.
		 instead we stop when u_i>=u_max-2*h_i (this is due to the way BC are handled 
		 
		 We might also want to stop if mesh spacing has gotten really small, i.e. ratio of original to current meshsize is greater than R_max
		 
		**/
	
		if( U[column_count-1]>=u_max-2*h[column_count-1] || h[0]/h[column_count-1]>R_max){
		
			// if the mesh spacing is really small, it's because a black hole has formed.  The current u value is a good approximation
			// to that of the event horizon
			if(h[0]/h[column_count-1]>R_max){
				contains_blackhole=true;
				u_horizon=U[column_count-1];
				m_bh = get_hawking_mass(U[column_count-2],0.5*(v_min+v_max));
			}
		
			return true;
		}
		return false;
	
	}

	void Gravitational_Collapse::first_column(){
		//allocate space for the first column
		q.push_back(new double[N[0]]);
		d.push_back(new double[N[0]]);
		a.push_back(new double[N[0]]);
		g.push_back(new double[N[0]]);
		r.push_back(new double[N[0]]);
		s.push_back(new double[N[0]]);
		f.push_back(new double[N[0]]);
		p.push_back(new double[N[0]]);
		c.push_back(new double[N[0]]);
		// fill in j=0 using BC
		Boundary_Conditions(0);
		for(int j=1;j<N[0];j++){  // fill in j>0
			Populate(0,j);
		}
		column_count++;
	
	}
	
	void Gravitational_Collapse::add_column(){
	
		// we are adding a row which will have INDEX i equal to the current column_count.
		// after we have added it the column_count will increment by one.
		int i=column_count;
	
		/** First calcuate the new step size, the u value for the new column, and the number of cells.
		    if v_max-U[i] is not divisible by h[i], we always overshoot and use one more cell,
		    rather than undershoot.	**/	
		h.push_back(get_new_h());
		U.push_back(U[i-1]+h[i]);
		N.push_back(ceil((v_max-U[i])/h[i]+1));
		
		//allocate space for the next column
		q.push_back(new double[N[i]]);
		d.push_back(new double[N[i]]);
		a.push_back(new double[N[i]]);
		g.push_back(new double[N[i]]);
		r.push_back(new double[N[i]]);
		s.push_back(new double[N[i]]);
		f.push_back(new double[N[i]]);
		p.push_back(new double[N[i]]);
		c.push_back(new double[N[i]]);
		
		Boundary_Conditions(i); // fill in j=0 using BC
		
		for(int j=1;j<N[i];j++){  // fill in j>0
			Populate(i,j);
		}
		
		column_count++;
	
	}
	
	double Gravitational_Collapse::get_new_h(){
		
		cout<<" current mesh ratio "<<minimum_projection(column_count-1)<<" ";
		return min(h[column_count-1],h[0]*minimum_projection(column_count-1));
		
	}	
	
	// fill in cell j=0 of column i, using boundary conditions.
	
	void Gravitational_Collapse::Boundary_Conditions(int i){
	
		double u=U[i];
		double v=u;
	
		// obviously r=0 at u=v which is really the origin.
	
		r[i][0]=0;
		
		// we are trying to set da/dr=ds/dr=0 at the boundary, using
		// a second order one sided difference stencil - but due
		// to nonuniform mesh size, some interpolation is needed
		// (calls to get_a etc) to generate virtual cells.
		
		// another complication is that we need to be on at least the third column (which has index 2, 
		// to further complicate things), for this to work.  In columns 0 and 1 we instead assume a=1
		// (flat space) and s=0 (scalar field negligible)

		if(i==0 || i==1){a[i][0]=1;s[i][0]=0;}
		else{
			a[i][0]=(4*get_a(i-1,v+h[i])-get_a(u-2*h[i],v+2*h[i]))/3;
			s[i][0]=(4*get_s(i-1,v+h[i])-get_s(u-2*h[i],v+2*h[i]))/3;
		}

		// next f=-g=a/2
		
		g[i][0]=0.5*a[i][0];
		f[i][0]=-g[i][0];
		
		
		//here we are making the assumption about things going to 0 fast
		
		double qhat,dhat;
		get_qhat_and_dhat(qhat,dhat,i,0);
		
		if(i==0){  // IC
			q[i][0]=qhat;
			p[i][0]=q[i][0];
			d[i][0]=dhat;
		}
		else{
			//here we are making the assumption about things going to 0 fast
			//cout<<"u-h[i] "<<u-h[i]<<endl;
			q[i][0]=0.5*(qhat+get_q(i-1,v));
			p[i][0]=q[i][0];
			d[i][0]=0.5*(dhat+get_d(i-1,v)-h[i]*qhat*p[i][0]);
		}
		c[i][0]=d[i][0];
		
	}	

	/**
		These are computed together for efficiency reasons.
		q & d are the only variable that depend on the previous column i-1.
		Due to the non-uniform step size, cells in row i-1 due not line
		up with cells in the current row i.  For this reason virtual cells
		are created by interpolating the previous row - hence the 
		calls to "get_q" etc.
	**/
	
	void Gravitational_Collapse::get_qhat_and_dhat(double &qhat,double &dhat,int i,int j){
	
		double v=U[i]+h[i]*j;
		
		if(i==0){ // use initial conditions
			qhat=2*Psi0*v*v*(3-2*(v-r0)*v/(delta*delta))*exp(-(v-r0)*(v-r0)/(delta*delta));
			dhat=0;
		}
		else{
			double qtemp=get_q(i-1,v); // q(u-du,v)
			double ftemp=get_f(i-1,v);
			double gtemp=get_g(i-1,v);
			double ptemp=get_p(i-1,v);
			double rtemp=get_r(i-1,v);
			double atemp=get_a(i-1,v);
			double dtemp=get_d(i-1,v);
			qhat=qtemp-h[i]*(ftemp*qtemp+gtemp*ptemp)/rtemp;
			dhat=dtemp+h[i]*((ftemp*gtemp+0.25*atemp*atemp)/(rtemp*rtemp)-ptemp*qtemp);
		}
	}

	void Gravitational_Collapse::Populate(int i,int j){
		
		double v=U[i]+j*h[i];
		   
		// first obtain qhat and dhat
		// (done together for efficiency reasons
			double qhat,dhat;
			get_qhat_and_dhat(qhat,dhat,i,j);
		
		// now get the other hatted variables
			double ahat=get_ahat(i,j,dhat);
			double ghat=get_ghat(i,j,qhat,dhat);
			double rhat=get_rhat(i,j,ghat);
			double shat=get_shat(i,j,qhat);
			double fhat=get_fhat(i,j,ghat,ahat,rhat);
			double phat=get_phat(i,j,ghat,qhat,fhat,rhat);
			double chat=get_chat(i,j);
			
			if(i==0){ // IC
				q[i][j]=qhat;
				d[i][j]=dhat;
			}
			else{
				q[i][j]=0.5*(qhat+get_q(i-1,v)-h[i]*(fhat*qhat+ghat*phat)/rhat);
				d[i][j]=0.5*(dhat+get_d(i-1,v)+h[i]*((fhat*ghat+0.25*ahat*ahat)/(rhat*rhat)-phat*qhat));
			}
			
			a[i][j]=0.5*(ahat+a[i][j-1]+h[i]*ahat*dhat);
			g[i][j]=0.5*(ghat+g[i][j-1]+h[i]*(2*dhat*ghat-rhat*qhat*qhat));
			r[i][j]=0.5*(rhat+r[i][j-1]+h[i]*ghat);
			s[i][j]=0.5*(shat+s[i][j-1]+h[i]*qhat);
			f[i][j]=0.5*(fhat+f[i][j-1]-h[i]*(fhat*ghat+0.25*ahat*ahat)/rhat);
			p[i][j]=0.5*(phat+p[i][j-1]-h[i]*(fhat*qhat+ghat*phat)/rhat);
			c[i][j]=0.5*(chat+c[i][j-1]+h[i]*((fhat*ghat+0.25*ahat*ahat)/(rhat*rhat)-phat*qhat));
		
		
	}

	double Gravitational_Collapse::get_phat(int i,int j,double ghat,double qhat,double fhat,double rhat){	
		//here we are making the goes to 0 assumption again...
		double chunk;
		if(j==1){
			chunk=0;
		}
		else{
			chunk=(f[i][j-1]*q[i][j-1]+g[i][j-1]*p[i][j-1])/r[i][j-1];
		}
		return (p[i][j-1]-0.5*h[i]*(chunk+fhat*qhat/rhat))/(1+0.5*h[i]*ghat/rhat);
	}

	double Gravitational_Collapse::get_fhat(int i,int j,double ghat,double ahat,double rhat){	
		//here we are making the goes to 0 assumption again...
		double chunk;
		if(j==1){
			chunk=0;
		}
		else{
			chunk=(f[i][j-1]*g[i][j-1]+0.25*a[i][j-1]*a[i][j-1])/r[i][j-1];
		}
		return (f[i][j-1]-0.5*h[i]*(chunk+0.25*ahat*ahat/rhat))/(1+0.5*h[i]*ghat/rhat);
	}		
	double Gravitational_Collapse::get_shat(int i,int j,double qhat){	
		return s[i][j-1]+0.5*h[i]*(q[i][j-1]+qhat);
	}		
	double Gravitational_Collapse::get_rhat(int i,int j,double ghat){	
		return r[i][j-1]+0.5*h[i]*(g[i][j-1]+ghat);
	}
	double Gravitational_Collapse::get_ghat(int i,int j,double qhat,double dhat){	
		return (g[i][j-1]+0.5*h[i]*(2*d[i][j-1]*g[i][j-1]-r[i][j-1]*q[i][j-1]*q[i][j-1]-qhat*qhat*(r[i][j-1]+0.5*h[i]*g[i][j-1])))/(1-h[i]*dhat+0.25*h[i]*h[i]*qhat*qhat);
	}	
	double Gravitational_Collapse::get_ahat(int i,int j,double dhat){	
		return a[i][j-1]*(1+0.5*h[i]*d[i][j-1])/(1-0.5*h[i]*dhat);
	}
	double Gravitational_Collapse::get_chat(int i,int j){
		
			if(j>1){
				return c[i][j-1]+h[i]*((f[i][j-1]*g[i][j-1]+0.25*a[i][j-1]*a[i][j-1])/(r[i][j-1]*r[i][j-1])-p[i][j-1]*q[i][j-1]);
			}
			else{
				return c[i][j-1]+h[i]*(-p[i][j-1]*q[i][j-1]);
			}
	}

	

	

	void Gravitational_Collapse::get_acdgfr(double u, double v, double &A, double &C, double &D, double &G, double &F, double &R){
		
		// check that the range of u & v is valid.  we don't simple call "inside universe" because we allow some extrapolation in the v direction.
		if(u<0 || u>u_horizon || v < u){
			cout<<"invalid interpolation value (u,v)=("<<u<<" "<<v<<")"<<endl;
			return;
		}
		
		//first get the i for u.
		int i=-1;
		int left = 0;
		int right = column_count - 1;

		while (left <= right) {
			int mid = left + (right - left) / 2;

			// Check the range condition
			if (U[mid] <= u && u < U[mid + 1]) {
				i=mid;
				break;
			}
			// If u is less than U[mid], search left half
			else if (u < U[mid]) {
				right = mid - 1;
			}
			// Otherwise, search right half
			else {
				left = mid + 1;
			}
		}

		if(i==-1){
			std::cout << "i not found" << std::endl;
			return;
		}
		
		
		
	
		
		// ok, so we have U[i]<=u<U[i+1].  So long as v>U[i+1], i.e. we are not in the first (triangle shaped) cell, things are easy.
		
		if(v>U[i+1]){	
			double t=(u-U[i])/(U[i+1]-U[i]);
			
			double interpolateAi, interpolateAip1;
			double interpolateCi, interpolateCip1;
			double interpolateDi, interpolateDip1;
			double interpolateGi, interpolateGip1;
			double interpolateFi, interpolateFip1;
			double interpolateRi, interpolateRip1;
			
			// get interpolateAi, etc
			// v is between j*h[i] and (j+1)*h[i], where
			int j=floor((v-U[i])/h[i]);
		
			double tj=(v-U[i]-j*h[i])/h[i];
			
			int j2=floor((v-U[i+1])/h[i+1]);
		
			double tj2=(v-U[i+1]-j2*h[i+1])/h[i+1];
	
			// for technical reasons, it may be that we need values of X 
			// up to two cells beyond the end of the array in the j direction...
			// we call extrapolate(X,i,j) instead of X[i][j] to take care of this.
		
			// on the other hand, we do NOT allow j<0 or i<0 or i>=column_count
			if(j<0 || i<0 || i>=column_count){
				cout<<"invalid interpolation value (i,j,U[i],v)="<<i<<" "<<j<<" "<<U[i]<<" "<<v<<endl;
				return;
			}
			double extrapolateAij,extrapolateAijp1,extrapolateAip1j,extrapolateAip1jp1;
			double extrapolateCij,extrapolateCijp1,extrapolateCip1j,extrapolateCip1jp1;
			double extrapolateDij,extrapolateDijp1,extrapolateDip1j,extrapolateDip1jp1;
			double extrapolateGij,extrapolateGijp1,extrapolateGip1j,extrapolateGip1jp1;
			double extrapolateFij,extrapolateFijp1,extrapolateFip1j,extrapolateFip1jp1;
			double extrapolateRij,extrapolateRijp1,extrapolateRip1j,extrapolateRip1jp1;
			
			// get extrapolateAij etc
			if(j<N[i]){extrapolateAij=a[i][j];
					   extrapolateCij=c[i][j];
					   extrapolateDij=d[i][j];
					   extrapolateGij=g[i][j];
					   extrapolateFij=f[i][j];
					   extrapolateRij=r[i][j];
					  }
			else{
				double dA=(3*a[i][N[i]-1]-4*a[i][N[i]-2]+a[i][N[i]-3])/2;
				extrapolateAij=a[i][N[i]-1]+(j-N[i]+1)*dA;
				double dC=(3*c[i][N[i]-1]-4*c[i][N[i]-2]+c[i][N[i]-3])/2;
				extrapolateCij=c[i][N[i]-1]+(j-N[i]+1)*dC;
				double dD=(3*d[i][N[i]-1]-4*d[i][N[i]-2]+d[i][N[i]-3])/2;
				extrapolateDij=d[i][N[i]-1]+(j-N[i]+1)*dD;
				double dG=(3*g[i][N[i]-1]-4*g[i][N[i]-2]+g[i][N[i]-3])/2;
				extrapolateGij=g[i][N[i]-1]+(j-N[i]+1)*dG;
				double dF=(3*f[i][N[i]-1]-4*f[i][N[i]-2]+f[i][N[i]-3])/2;
				extrapolateFij=f[i][N[i]-1]+(j-N[i]+1)*dF;
				double dR=(3*r[i][N[i]-1]-4*r[i][N[i]-2]+r[i][N[i]-3])/2;
				extrapolateRij=r[i][N[i]-1]+(j-N[i]+1)*dR;
			}
			// get extrapolateAijp1 etc
			if(j+1<N[i]){extrapolateAijp1=a[i][j+1];
					   extrapolateCijp1=c[i][j+1];
					   extrapolateDijp1=d[i][j+1];
					   extrapolateGijp1=g[i][j+1];
					   extrapolateFijp1=f[i][j+1];
					   extrapolateRijp1=r[i][j+1];
					  }
			else{
				double dA=(3*a[i][N[i]-1]-4*a[i][N[i]-2]+a[i][N[i]-3])/2;
				extrapolateAijp1=a[i][N[i]-1]+(j+1-N[i]+1)*dA;
				double dC=(3*c[i][N[i]-1]-4*c[i][N[i]-2]+c[i][N[i]-3])/2;
				extrapolateCijp1=c[i][N[i]-1]+(j+1-N[i]+1)*dC;
				double dD=(3*d[i][N[i]-1]-4*d[i][N[i]-2]+d[i][N[i]-3])/2;
				extrapolateDijp1=d[i][N[i]-1]+(j+1-N[i]+1)*dD;
				double dG=(3*g[i][N[i]-1]-4*g[i][N[i]-2]+g[i][N[i]-3])/2;
				extrapolateGijp1=g[i][N[i]-1]+(j+1-N[i]+1)*dG;
				double dF=(3*f[i][N[i]-1]-4*f[i][N[i]-2]+f[i][N[i]-3])/2;
				extrapolateFijp1=f[i][N[i]-1]+(j+1-N[i]+1)*dF;
				double dR=(3*r[i][N[i]-1]-4*r[i][N[i]-2]+r[i][N[i]-3])/2;
				extrapolateRijp1=r[i][N[i]-1]+(j+1-N[i]+1)*dR;
			}
			// get extrapolateAip1j etc
			if(j2<N[i+1]){extrapolateAip1j=a[i+1][j2];
					   extrapolateCip1j=c[i+1][j2];
					   extrapolateDip1j=d[i+1][j2];
					   extrapolateGip1j=g[i+1][j2];
					   extrapolateFip1j=f[i+1][j2];
					   extrapolateRip1j=r[i+1][j2];
					  }
			else{
				double dA=(3*a[i+1][N[i+1]-1]-4*a[i+1][N[i+1]-2]+a[i+1][N[i+1]-3])/2;
				extrapolateAip1j=a[i+1][N[i+1]-1]+(j2-N[i+1]+1)*dA;
				double dC=(3*c[i+1][N[i+1]-1]-4*c[i+1][N[i+1]-2]+c[i+1][N[i+1]-3])/2;
				extrapolateCip1j=c[i+1][N[i+1]-1]+(j2-N[i+1]+1)*dC;
				double dD=(3*d[i+1][N[i+1]-1]-4*d[i+1][N[i+1]-2]+d[i+1][N[i+1]-3])/2;
				extrapolateDip1j=d[i+1][N[i+1]-1]+(j2-N[i+1]+1)*dD;
				double dG=(3*g[i+1][N[i+1]-1]-4*g[i+1][N[i+1]-2]+g[i+1][N[i+1]-3])/2;
				extrapolateGip1j=g[i+1][N[i+1]-1]+(j2-N[i+1]+1)*dG;
				double dF=(3*f[i+1][N[i]-1]-4*f[i+1][N[i+1]-2]+f[i+1][N[i+1]-3])/2;
				extrapolateFip1j=f[i+1][N[i+1]-1]+(j2-N[i+1]+1)*dF;
				double dR=(3*r[i+1][N[i+1]-1]-4*r[i+1][N[i+1]-2]+r[i+1][N[i+1]-3])/2;
				extrapolateRip1j=r[i+1][N[i+1]-1]+(j2-N[i+1]+1)*dR;
			}
			// get extrapolateAip1jp1 etc
			if(j2+1<N[i+1]){extrapolateAip1jp1=a[i+1][j2+1];
					   extrapolateCip1jp1=c[i+1][j2+1];
					   extrapolateDip1jp1=d[i+1][j2+1];
					   extrapolateGip1jp1=g[i+1][j2+1];
					   extrapolateFip1jp1=f[i+1][j2+1];
					   extrapolateRip1jp1=r[i+1][j2+1];
					  }
			else{
				double dA=(3*a[i+1][N[i+1]-1]-4*a[i+1][N[i+1]-2]+a[i+1][N[i+1]-3])/2;
				extrapolateAip1jp1=a[i+1][N[i+1]-1]+(j2+1-N[i+1]+1)*dA;
				double dC=(3*c[i+1][N[i+1]-1]-4*c[i+1][N[i+1]-2]+c[i+1][N[i+1]-3])/2;
				extrapolateCip1jp1=c[i+1][N[i+1]-1]+(j2+1-N[i+1]+1)*dC;
				double dD=(3*d[i+1][N[i+1]-1]-4*d[i+1][N[i+1]-2]+d[i+1][N[i+1]-3])/2;
				extrapolateDip1jp1=d[i+1][N[i+1]-1]+(j2+1-N[i+1]+1)*dD;
				double dG=(3*g[i+1][N[i+1]-1]-4*g[i+1][N[i+1]-2]+g[i+1][N[i+1]-3])/2;
				extrapolateGip1jp1=g[i+1][N[i+1]-1]+(j2+1-N[i+1]+1)*dG;
				double dF=(3*f[i+1][N[i]-1]-4*f[i+1][N[i+1]-2]+f[i+1][N[i+1]-3])/2;
				extrapolateFip1jp1=f[i+1][N[i+1]-1]+(j2+1-N[i+1]+1)*dF;
				double dR=(3*r[i+1][N[i+1]-1]-4*r[i+1][N[i+1]-2]+r[i+1][N[i+1]-3])/2;
				extrapolateRip1jp1=r[i+1][N[i+1]-1]+(j2+1-N[i+1]+1)*dR;
			}
			
			interpolateAi=(1-tj)*extrapolateAij+tj*extrapolateAijp1;
			interpolateCi=(1-tj)*extrapolateCij+tj*extrapolateCijp1;
			interpolateDi=(1-tj)*extrapolateDij+tj*extrapolateDijp1;
			interpolateGi=(1-tj)*extrapolateGij+tj*extrapolateGijp1;
			interpolateFi=(1-tj)*extrapolateFij+tj*extrapolateFijp1;
			interpolateRi=(1-tj)*extrapolateRij+tj*extrapolateRijp1;
			
			interpolateAip1=(1-tj2)*extrapolateAip1j+tj2*extrapolateAip1jp1;
			interpolateCip1=(1-tj2)*extrapolateCip1j+tj2*extrapolateCip1jp1;
			interpolateDip1=(1-tj2)*extrapolateDip1j+tj2*extrapolateDip1jp1;
			interpolateGip1=(1-tj2)*extrapolateGip1j+tj2*extrapolateGip1jp1;
			interpolateFip1=(1-tj2)*extrapolateFip1j+tj2*extrapolateFip1jp1;
			interpolateRip1=(1-tj2)*extrapolateRip1j+tj2*extrapolateRip1jp1;
			//interpolateAip1=(1-t2)*extrapolateAip1j+t2*extrapolateAip1jp1;
			
			
			A=(1-t)*interpolateAi+t*interpolateAip1;
			C=(1-t)*interpolateCi+t*interpolateCip1;
			D=(1-t)*interpolateDi+t*interpolateDip1;
			G=(1-t)*interpolateGi+t*interpolateGip1;
			F=(1-t)*interpolateFi+t*interpolateFip1;
			R=(1-t)*interpolateRi+t*interpolateRip1;
			
			
			
		}
		// if not things are a little more complicated - essentially we are interpolating on a triangle instead of a square.
		else{
		
		
			// we use barycentric coordinates...straight from wikipedia
			double x1=U[i],y1=U[i];
			double x2=U[i+1],y2=U[i+1];
			double x3=U[i],y3=U[i]+h[i];
			
			double x=u;double y=v;
			
			double det=(y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
			
			double l1=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/det;
			double l2=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/det;
			double l3=1-l1-l2;
			A=l1*a[i][0]+l2*a[i+1][0]+l3*a[i][1];
			C=l1*c[i][0]+l2*c[i+1][0]+l3*c[i][1];
			D=l1*d[i][0]+l2*d[i+1][0]+l3*d[i][1];
			G=l1*g[i][0]+l2*g[i+1][0]+l3*g[i][1];
			F=l1*f[i][0]+l2*f[i+1][0]+l3*f[i][1];
			R=l1*r[i][0]+l2*r[i+1][0]+l3*r[i][1];
			
		}
		
	}
	
	

	

	
	/** X[i] is an array of length N[i].
	 Thus X[i][j] is only defined for j<N[i].
	 But sometimes we might need to add one or two ghost cell's X[i][N[i]],X[i][N[i]+1]
	 This is done to 2nd order using a second order one sided stencil for the derivative 
	 
	 Calling this function is like calling X[i][j], but with ghost cells used if j overshoots
	 
	 **/
	
	double Gravitational_Collapse::extrapolate(vector<double*> X,int i,int j){
		if(j<N[i]){return X[i][j];}
		else{
			double dX=(3*X[i][N[i]-1]-4*X[i][N[i]-2]+X[i][N[i]-3])/2;
			return X[i][N[i]-1]+(j-N[i]+1)*dX;
		}
	}
	
	double Gravitational_Collapse::interpolate(vector<double*> X,int i,double v){
	
		// v is between j*h[i] and (j+1)*h[i], where
		int j=floor((v-U[i])/h[i]);
		
		double t=(v-U[i]-j*h[i])/h[i];
	
		// for technical reasons, it may be that we need values of X 
		// up to two cells beyond the end of the array in the j direction...
		// we call extrapolate(X,i,j) instead of X[i][j] to take care of this.
		
		// on the other hand, we do NOT allow j<0 or i<0 or i>=column_count
		if(j<0 || i<0 || i>=column_count){
					cout<<"invalid interpolation value (i,j,U[i],v)="<<i<<" "<<j<<" "<<U[i]<<" "<<v<<endl;
					return 0;
		}
		
		return (1-t)*extrapolate(X,i,j)+t*extrapolate(X,i,j+1);
	}
	
	double Gravitational_Collapse::interpolate(vector<double*> X,double u,double v){
	
		// check that the range of u & v is valid.  we don't simple call "inside universe" because we allow some extrapolation in the v direction.
		if(u<0 || u>u_horizon || v < u){
			cout<<"invalid interpolation value (u,v)=("<<u<<" "<<v<<")"<<endl;
			return 0;
		}
	
		// u is between U[i] and U[i+1], for some i.
		int i=find_i(u); // just a linear search
		
		// ok, so we have U[i]<=u<U[i+1].  So long as v>U[i+1], i.e. we are not in the first (triangle shaped) cell, things are easy.
		
		if(v>U[i+1]){	
			double t=(u-U[i])/(U[i+1]-U[i]);
			return (1-t)*interpolate(X,i,v)+t*interpolate(X,i+1,v);
		}
		// if not things are a little more complicated - essentially we are interpolating on a triangle instead of a square.
		else{
		
		
			// we use barycentric coordinates...straight from wikipedia
			double x1=U[i],y1=U[i];
			double x2=U[i+1],y2=U[i+1];
			double x3=U[i],y3=U[i]+h[i];
			
			double x=u;double y=v;
			
			double det=(y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
			
			double l1=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/det;
			double l2=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/det;
			double l3=1-l1-l2;
			return l1*X[i][0]+l2*X[i+1][0]+l3*X[i][1];
			
		}
	}
	
	// find the i such that U[i]<=u<U[i+1]
	
	int Gravitational_Collapse::find_i(double u){
	
		//for(int i=0;i<column_count-1;i++){
		//	if(U[i+1]>u){return i;}
		//}
		//cout<<"i not found"<<endl;
		//return 0; // not found
		int left = 0;
		int right = column_count - 1;

		while (left <= right) {
			int mid = left + (right - left) / 2;

			// Check the range condition
			if (U[mid] <= u && u < U[mid + 1]) {
				return mid;
			}
			// If u is less than U[mid], search left half
			else if (u < U[mid]) {
				right = mid - 1;
			}
			// Otherwise, search right half
			else {
				left = mid + 1;
			}
		}

		// If no valid index is found
		std::cout << "i not found" << std::endl;
		return 0;
		
		//auto it = lookup_table.lower_bound(u);

        // Handle edge cases when `u` is smaller than the first element
        //if (it == lookup_table.begin() && u < it->first) {
        //    cout << "i not found (below range)" << endl;
        //    return 0;
        //}

        // If lower_bound points to an element > u, step back
        //if (it == lookup_table.end() || it->first > u) {
        //    --it;
        //}

        // Check if `u` fits the correct interval
        //int i = it->second;
        //if (U[i] <= u && u < U[i+1]) {
        //    return i;
        //}

        //cout << "i not found" << endl;
        //return 0;  // Not found
	
	}
	

	
	
	double Gravitational_Collapse::get(double u,double v,string id){
		if(id.compare("q")==0){
			return get_q(u,v);
		}
		else if(id.compare("d")==0){
			return get_d(u,v);
		}
		else if(id.compare("a")==0){
			return get_a(u,v);
		}
		else if(id.compare("g")==0){
			return get_g(u,v);
		}
		else if(id.compare("r")==0){
			return get_r(u,v);
		}
		else if(id.compare("s")==0){
			return get_s(u,v);
		}
		else if(id.compare("f")==0){
			return get_f(u,v);
		}
		else if(id.compare("p")==0){
			return get_p(u,v);
		}
		else if(id.compare("c")==0){
			return get_c(u,v);
		}
		else{
			cout<<"variable name not found"<<endl;
			return -1;
		}
	}
	
	double Gravitational_Collapse::get(int i,double v,string id){
		if(id.compare("q")==0){
			return get_q(i,v);
		}
		else if(id.compare("d")==0){
			return get_d(i,v);
		}
		else if(id.compare("a")==0){
			return get_a(i,v);
		}
		else if(id.compare("g")==0){
			return get_g(i,v);
		}
		else if(id.compare("r")==0){
			return get_r(i,v);
		}
		else if(id.compare("s")==0){
			return get_s(i,v);
		}
		else if(id.compare("f")==0){
			return get_f(i,v);
		}
		else if(id.compare("p")==0){
			return get_p(i,v);
		}	
		else if(id.compare("c")==0){
			return get_c(i,v);
		}		
		else{
			cout<<"variable name not found"<<endl;
			return -1;
		}
	}	
	
	double Gravitational_Collapse::get_q(double u,double v){
		return interpolate(q,u,v);
	}
	double Gravitational_Collapse::get_d(double u,double v){
		return interpolate(d,u,v);
	}
	double Gravitational_Collapse::get_a(double u,double v){
		return interpolate(a,u,v);
	}	
	double Gravitational_Collapse::get_g(double u,double v){
		return interpolate(g,u,v);
	}	
	double Gravitational_Collapse::get_r(double u,double v){
		return interpolate(r,u,v);
	}
	double Gravitational_Collapse::get_s(double u,double v){
		return interpolate(s,u,v);
	}
	double Gravitational_Collapse::get_f(double u,double v){
		return interpolate(f,u,v);
	}	
	double Gravitational_Collapse::get_p(double u,double v){
		return interpolate(p,u,v);
	}
	double Gravitational_Collapse::get_c(double u,double v){
		return interpolate(c,u,v);
	}	
	double Gravitational_Collapse::get_q(int i,double v){
		return interpolate(q,i,v);
	}
	double Gravitational_Collapse::get_d(int i,double v){
		return interpolate(d,i,v);
	}
	double Gravitational_Collapse::get_a(int i,double v){
		return interpolate(a,i,v);
	}	
	double Gravitational_Collapse::get_g(int i,double v){
		return interpolate(g,i,v);
	}	
	double Gravitational_Collapse::get_r(int i,double v){
		return interpolate(r,i,v);
	}
	double Gravitational_Collapse::get_s(int i,double v){
		return interpolate(s,i,v);
	}
	double Gravitational_Collapse::get_f(int i,double v){
		return interpolate(f,i,v);
	}	
	double Gravitational_Collapse::get_p(int i,double v){
		return interpolate(p,i,v);
	}
	double Gravitational_Collapse::get_c(int i,double v){
		return interpolate(c,i,v);
	}		
	
	double Gravitational_Collapse::get_a_u(double u,double v){
	
		int i=find_i(u);
	
		if(u<=2*h[i]){ //flat space on first band of width h_i
			return 0;
		}
		else if(v>=u+h[i]){ // centered differences if we are far enough from u=v
			return (interpolate(a,u+h[i],v)-interpolate(a,u-h[i],v))/(2*h[i]);
		}
		else if(v-u>=0){ // we're right near u=v, so we have to use a different O(h^2) stencil
			return (3*interpolate(a,u,v)-4*interpolate(a,u-h[i],v)+interpolate(a,u-2*h[i],v))/(2*h[i]);
		}
		else{   // I think we aren't even in a valid domain in this case
			return 0;
		}
	
	}	
	
	// returns the cosine of the angle that a constant areal radius isocontour passing through
	// (u,v) makes with the u axis.
	
	double Gravitational_Collapse::get_isocontour_projection(double u,double v){
	
		return 1./sqrt(1+pow(get_f(u,v)/get_g(u,v),2));
	
	}
	
	// returns the cosine of the angle that a constant areal radius isocontour passing through
	// column i at v makes with the u axis.
	
	double Gravitational_Collapse::get_isocontour_projection(int i,double v){
	
		return 1./sqrt(1+pow(get_f(i,v)/get_g(i,v),2));
	
	}	
	
	// computes the minimum of the above function over a column of constant i
	
	double Gravitational_Collapse::minimum_projection(int i){
	
		double minimum=1; // the minimum will always be less than 1
		
		double u=U[i];
		double v,candidate;
		
		for(int j=0;j<N[i];j++){
			v=U[i]+j*h[i];
			candidate=get_isocontour_projection(i,v);
			if(candidate<minimum){minimum=candidate;}
		}
	
		return minimum;
	
	}

	double Gravitational_Collapse::get_local_mesh_size(double u,double v){
	
		int i=find_i(u);
		return h[i];
	
	}		
	
	
	
	
	
	
	// check if we are in the portion of spacetime on which we have solved the Einstein equations.
	// We do not include boundaries, as this could cause interpolation errors
	
	bool Gravitational_Collapse::inside_universe(double u,double v){
			//if(u<U[0]){return 0;} 
			//else if(u>U[column_count-1] || v>v_max || v<u){return 0;}
			//else{return 1;}
            //if(v>=u && v<= v_max && u>=U[0] && u<=U[column_count-1]){return true;}
			if(v>=u && v<= v_max && u<=U[column_count-1]){return true;}
            else{
				/**cout<<"(u,v)="<<u<<" "<<v<<endl;
				cout<<"U[0]="<<U[0]<<endl;
				cout<<"U[column_count-1]="<<U[column_count-1]<<endl;
				cout<<"U[column_count-1]="<<U[column_count-1]<<endl;
				cout<<"v_max="<<v_max<<endl;
				cout<<"column_count = "<<column_count<<endl;**/
				return false;}
	}
	
	double Gravitational_Collapse::get_hawking_mass(double u,double v){
		double a = get_a(u,v);
		return 0.5*get_r(u,v)*(1.0+4*get_f(u,v)*get_g(u,v)/(a*a));
	}
	
	void Gravitational_Collapse::accel(double &udd,double &vdd,double &phidd,double u,double v,double phi,double udot,double vdot,double phidot,double phi0){
		double A,C,D,G,F,R;
		double A_flat,C_flat,D_flat,G_flat,F_flat,R_flat;
		//double A2,C2,D2,G2,F2,R2;
		
		R_flat=(v-u)/2;
		A_flat=1;
		C_flat=0;
		D_flat=0;
		F_flat=-0.5;
		G_flat=0.5;
		
		
		if(u>1e-3){
			get_acdgfr(u,v,A,C,D,G,F,R);
		}
		else{
			R=R_flat;
			A=A_flat;
			C=C_flat;
			D=D_flat;
			F=F_flat;
			G=G_flat;
		}
		
		double r0 = 4.0;
	
		//double weight = (1+tanh(r0-R))/2;
		
		double y = fabs(R*sin(phi));
		double weight = (1+tanh(r0-y))/2;
		
		//weight = 1;
		double R_blend = R*weight+R_flat*(1-weight);
		double A_blend = A*weight+A_flat*(1-weight);
		double C_blend = C*weight+C_flat*(1-weight);
		double D_blend = D*weight+D_flat*(1-weight);
		double F_blend = F*weight+F_flat*(1-weight);
		double G_blend = G*weight+G_flat*(1-weight);
		
		//double x = R*cos(phi);
		//double y = R*sin(phi);
		
		//double phi0p = atan2(y,x);
		//double phi1 = this->phi;
		//cout<<"phi1"<<phi1<<endl;
		//if(-pi+phi0<=phi0p && phi0p<=phi0 && R > 2.0){
		/**if(abs(phi1)>1e-1){
			//assume we are in flat space
			//udd=-0.5*(v-u)*phidot*phidot;
			//vdd=0.5*(v-u)*phidot*phidot;
			//phidd=2*((udot-vdot)/(v-u))*phidot;
			udd = -1*(-1+cos(phi1))*(C*udot+D*vdot)/A;
			vdd = -1*(-1-cos(phi1))*(C*udot+D*vdot)/A;
			phidd = -sin(phi1)*(F*udot+G*vdot)/(R*R);
			return;
		}**/
		
		
		//A2 = get_a(u,v);
		//C2 = get_c(u,v);
		//D2 = get_d(u,v);
		///G2 = get_g(u,v);
		//F2 = get_f(u,v);
		//R2 = get_r(u,v);
		//if(!(A2==A) || !(C2==C) || !(D2==D) || !(G2==G) || !(F2==F) || !(R2==R)){
		//	cout<<"A A2 C C2 D D2 G G2 F F2 R R2 = "<<A<<" "<<A2<<" "<<C<<" "<<C2<<" "<<D<<" "<<D2<<" "<<G<<" "<<G2<<" "<<F<<" "<<F2<<" "<<R<<" "<<R2<<endl;
		//}
		
		udd=-2*C_blend*udot*udot-2*R_blend*G_blend*phidot*phidot/(A_blend*A_blend);
		vdd=-2*D_blend*vdot*vdot-2*R_blend*F_blend*phidot*phidot/(A_blend*A_blend);
		phidd=-2*G_blend*vdot*phidot/R_blend-2*F_blend*udot*phidot/R_blend;
		
		//udd=-2*C*udot*udot-2*R*G*phidot*phidot/(A*A);
		//vdd=-2*D*vdot*vdot-2*R*F*phidot*phidot/(A*A);
		//phidd=-2*G*vdot*phidot/R-2*F*udot*phidot/R;
		//if(u>0){
		//	udd=-2*get_c(u,v)*udot*udot-2*get_r(u,v)*get_g(u,v)*phidot*phidot/pow(get_a(u,v),2);
		//	vdd=-2*get_d(u,v)*vdot*vdot-2*get_r(u,v)*get_f(u,v)*phidot*phidot/pow(get_a(u,v),2);
		//	phidd=-2*get_g(u,v)*vdot*phidot/get_r(u,v)-2*get_f(u,v)*udot*phidot/get_r(u,v);
		//}
		/**else{
			//assume we are in flat space
			udd=-0.5*(v-u)*phidot*phidot;
			vdd=0.5*(v-u)*phidot*phidot;
			phidd=2*((udot-vdot)/(v-u))*phidot;
		}**/
	}

	void Gravitational_Collapse::constant_r_flow_vector(double u,double v,double &ud,double &vd){
	
		ud=sqrt(-get_g(u,v)/get_f(u,v))/get_a(u,v);
		//vd=sqrt(-get_f(u,v)/get_g(u,v))/get_a(u,v);
		vd=1.0/(pow(get_a(u,v),2)*ud);
	
	}	
	
	// it is useful to convert the acceleration into a first order system, so that methods like 
	// adaptive runge kutta can be used directly
	
	// here y=[ u,v,phi,udot,vdot,phidot], yp=[udot,vdot,phidot,udotdot,vdotdot,phidotdot]
	
	// returns true if u,v is part of the universe and false otherwise
	
	// because C++ is a bit of a pain the function needs to take a void * parameter that is really just "self",
	// because non-static member function pointers are a pain in the ass to deal with, so we have to "trick" C++
	// into thinking the function is static in a convoluted way.
	
	bool Gravitational_Collapse::accel_1st_order(double*& yp,double* y, double E){ 
	

	
		if(!this->inside_universe(y[0],y[1])){return false;}
	
		yp[0]=y[3]; // u'=udot
		yp[1]=y[4]; // v'=vdot
		yp[2]=y[5]; // phi'=phidot
		double udd,vdd,phidd;
		// I have rigged things so that E is actually phi0.
		this->accel(udd,vdd,phidd,y[0],y[1],y[2],y[3],y[4],y[5],E);
		yp[3]=udd;
		yp[4]=vdd;
		yp[5]=phidd;
		return true;
	}
	
	// same deal, really.
	
	bool Gravitational_Collapse::constant_r_flow(double*& yp,double* y){
	
		// y[0], y[1], y[2]    = u, v, phi
		// yp[0], yp[1], yp[2] = udot, vdot, phidot
	
	
		if(!this->inside_universe(y[0],y[1])){return false;}
	
		this->constant_r_flow_vector(y[0],y[1],yp[0],yp[1]);
		yp[2]=0;
		return true;
	
	}
	
	bool Gravitational_Collapse::Adaptive_Step_Constant_R(double* y,bool& y_changed,int n,double& h,double tol,bool force_step,double r_desired){
		double* k1=new double[n+1];
		double* k2=new double[n+1];
		double* k3=new double[n+1];
		double* k4=new double[n+1];
		double* k5=new double[n+1];
		double* k6=new double[n+1];

		
		//double h_min = 1.0;
		//if( danger_zone(y) ){
		//	h_min=0.1;
		//}
		double alpha = 0.3;
		double h_min = alpha*get_local_mesh_size(y[0],y[1]);
		
		double y1[n],Delta4[n],Delta5[n],diff54[n];
		
		bool inside_universe = true;
		
		if( !constant_r_flow(k1,y) ){inside_universe = false;} // attempts to evaluate f and returns false if y is outside Domain(f)
		for(int i=0;i<n;i++){k1[i]*=h;  
							 y1[i]=y[i]+0.25*k1[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !constant_r_flow(k2,y1) ){inside_universe = false;}
		for(int i=0;i<n;i++){k2[i]*=h;  
							 y1[i]=y[i]+(3.0/32)*k1[i]+(9.0/32)*k2[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;					 
		if( !constant_r_flow(k3,y1) ){inside_universe = false;}
		for(int i=0;i<n;i++){k3[i]*=h;  
							 y1[i]=y[i]+(1932.0/2197)*k1[i]-(7200.0/2197)*k2[i]+(7296.0/2197)*k3[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !constant_r_flow(k4,y1) ){inside_universe = false;}
		for(int i=0;i<n;i++){k4[i]*=h;  
							 y1[i]=y[i]+(439.0/216)*k1[i]-8*k2[i]+(3680.0/513)*k3[i]-(845.0/4104)*k4[i];}		
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;
		if( !constant_r_flow(k5,y1) ){inside_universe = false;}
		for(int i=0;i<n;i++){k5[i]*=h;  
							 y1[i]=y[i]-(8.0/27)*k1[i]+2*k2[i]-(3544.0/2565)*k3[i]+(1859.0/4104)*k4[i]-(11.0/40)*k5[i];}
		//std::cout<<"y1 =="<<y1[0]<<" "<<y1[1]<<" "<<y1[2]<<" "<<y1[3]<<" "<<y1[4]<<" "<<y1[5]<<std::endl;				 
		if( !constant_r_flow(k6,y1) ){inside_universe = false;}	
		for(int i=0;i<n;i++){k6[i]*=h;
							 Delta4[i]=(25.0/216)*k1[i]+(1408.0/2565)*k3[i]+(2197.0/4104)*k4[i]-(1.0/5)*k5[i];
							 Delta5[i]=(16.0/135)*k1[i]+(6656.0/12825)*k3[i]+(28561.0/56430)*k4[i]-(9.0/50)*k5[i]+(2.0/55)*k6[i];
							 diff54[i]=Delta5[i]-Delta4[i];
			}
		double q,eh;

		double r = get_r(y[0]+Delta5[0],y[1]+Delta5[1]);

		//eh=error_norm(y,diff54)/h;
		eh = abs(r-r_desired)/h;


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
	
	
	void Gravitational_Collapse::constant_r_trajectory(double u0,double v0,double tau_max,vector<double>& u_path,vector<double>& v_path, vector<double>& r_path, vector<double>& tau_path, vector<double>& m_path){
	
		u_path.clear();
		v_path.clear();
		r_path.clear();
		tau_path.clear();
		m_path.clear();
		u_path.push_back(u0);
		v_path.push_back(v0);
		r_path.push_back(get_r(u0,v0));
		tau_path.push_back(0.0);
		m_path.push_back(get_hawking_mass(u0,v0));

		double tau=0;
		
		bool inside_universe=true;
		bool step_taken;
		
		double u0p,v0p;
		
		double r_desired = get_r(u0,v0);
		
		constant_r_flow_vector(u0,v0,u0p,v0p);
		
		//double H=get_stepsize(u,v,ud,vd);
		double alpha = 0.3;
		//alpha*Universe->get_local_mesh_size(U,V)/max(fabs(Udot),fabs(Vdot));
		double H = alpha*get_local_mesh_size(u0,v0)/max(fabs(u0p),fabs(v0p));
		double H_old;
		
		double y[2]={u0,v0};
		
		double tol = 0.01;
		bool force_step = false;
		
		while(tau<tau_max && inside_universe){
			
			H_old=H;
			//inside_universe=take_constant_r_adaptive_step(H,step_taken); // tries to take a step but will not if the error is too high, in this case H is reduced and we try again.
			inside_universe=Adaptive_Step_Constant_R(y,step_taken,2,H,tol,force_step,r_desired);
			//if(H<get_stepsize(u,v,ud,vd)){H=get_stepsize(u,v,ud,vd);}
			if(step_taken){
				u_path.push_back(y[0]);
				v_path.push_back(y[1]);
				r_path.push_back(get_r(y[0],y[1]));
				m_path.push_back(get_hawking_mass(y[0],y[1]));
				tau+=H_old;
				tau_path.push_back(tau);
				//cout<<"y[0]="<<y[0]<<" "<<"y[1]="<<y[1]<<" "<<"tau = "<<tau<<endl;
			}
			
		
		}

	
	}
	
	
	/** this is a bit of a subtle point
	
	Sometimes it is useful to have a notion of when two points in spacetime are nearby.  A concrete example is 
	when we are numerically solving for trajectories through spacetime.  Simply computing
	the Euclidean norm of the coordinates (u1-u2)^2+(v1-v2)^2 is not good enough, since we know that larger and larger regions of the
	physical universe get squeezed into smaller and smaller patches of uv coordinate space as the event horizon is approached.  The spacetime
	distance is also not good, because (for example) null separated points are not 'the same point' even though the distance between them is zero.
	So instead what I do is the following:
	
	The reason the Euclidean norm is failing is because a(u,v) is going to infinity as the event horizon is approached.  We need a norm that reflects this.
	
	dS=a(u,v)*{ |du|+|dv| }
	
	should be fine.
	
	finally, we assume that p=(u,v) and x=(du,dv).  We think of p as being a point and x a vector  This format is useful as it can be used
	directly in adaptive Runge Kutta routines.
	
	**/
	
	// see "accel_1st_order" for why we need to do the whole dummy thing
	
	double Gravitational_Collapse::error_norm(double *p,double *x){
	
		double u = p[0];
		if(u>1e-3){		
			return this->get_a(p[0],p[1])*( fabs(x[0])+fabs(x[1]) );
		}
		else{
			return 1.0*( fabs(x[0])+fabs(x[1]) );
		}
	
	}
	
	//void Gravitational_Collapse::buildHashTable() {
    //    for (int i = 0; i < column_count - 1; i++) {
    //        lookup_table[U[i]] = i;
    //    }
    //}
	
	
	
	
	
	
	