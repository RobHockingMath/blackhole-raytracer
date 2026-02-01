#include "Oppenheimer_Snyder_Collapse.h"
#include <cmath>
#include <iostream>
#include <iomanip>   // <-- Add this line

using namespace std;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
Oppenheimer_Snyder_Collapse::Oppenheimer_Snyder_Collapse(double pos_[3], double tol_, double mass, 
                                     double R0_)
: Spacetime(tol_), 
  M(mass), R0(R0_), 
  X0(std::asin(std::sqrt(2.0*M/R0))),
  am(std::sqrt((R0*R0*R0)/(2*M))),
  rho0(M / ((4.0/3.0)*M_PI*R0*R0*R0)),
  C((8*M_PI/3.0)*rho0*am*am*am)
{
    // Copy initial position array, if needed
    for(int i=0; i<3; i++){
        pos[i] = pos_[i];
    }
	
	has_dustcloud=true;
	
	//double eps = 1e-9;
	//h = R0/1000000;
	//integrate_a(eps,h);
	
	double tau_collapse = 0.5*am*M_PI;
	
	cout<<"dust cloud collapses in "<<tau_collapse<<" units of proper time."<<endl;
	//cout<<"dust cloud collapses in "<<m_tau[m_tau.size()-1]<<" units of proper time."<<endl;
	
	// I am assuming we always start outside the dust cloud.
}

void Oppenheimer_Snyder_Collapse::integrate_a(double eps,double h){

	double aInit = am - eps;
    if(aInit <= 0.0) {
        std::cerr << "Warning: aInit <= 0. Increase eps?\n";
        aInit = am*0.999; // fallback
    }

    // 4) do a simple ODE integration 
    double tau=0.0;
    double a=aInit;

    m_tau.clear();
    m_a.clear();
    m_tau.push_back(tau);
    m_a.push_back(a);

    // integrate until a <= 0
    while(a>0){
        // compute dadtau = - sqrt(C/a -1)
        double deriv = C/a - 1.0;
        if(deriv < 0.0) {
            // means a> C in principle, or star is not collapsing?
            // we can break or something
            std::cerr << "Encountered deriv<0 => no real collapse\n";
            break;
        }
        double dadtau = -std::sqrt(deriv);

		double aMid   = a + 0.5*h*dadtau;

		double derivMid   = C/aMid - 1.0;
		double dadtauMid  = -std::sqrt(derivMid);

		//double aNext = a + h*dadtauMid;
		
		double aNext = a + h*dadtauMid;


        //double aNext = a + h*dadtau;
        double tauNext = tau + h;

        // store
        m_tau.push_back(tauNext);
        m_a.push_back(aNext);

        // update
        tau = tauNext;
        a   = aNext;

        if(aNext<=0.0) {
            // done
            break;
        }
    }

}


double Oppenheimer_Snyder_Collapse::tau_of_eta(double eta)
{
    // tau(\eta) = (a_m/2) * [ eta + sin(eta) ]
    return 0.5*am*(eta + std::sin(eta));
}

double Oppenheimer_Snyder_Collapse::a_of_eta(double eta)
{
    // a(\eta) = (a_m/2) * [1 + cos(eta)]
    return 0.5*am*(1.0 + std::cos(eta));
}

double Oppenheimer_Snyder_Collapse::get_a_old2(double tau)
{
    // 1) The total collapse time (when eta = pi):
    //    tau_collapse = (a_m/2)*(pi + sin(pi)) = (a_m/2)*pi.
    double tau_collapse = 0.5*am*M_PI;

    // If tau is beyond collapse, the star is fully collapsed => a=0
    if(tau >= tau_collapse) {
        return 0.0;
    }
    // If tau <= 0, at earliest times => a(\tau)=a_m
    // (We can clamp, or handle negative tau if you wish, but typically OS is 0..tau_coll)
    if(tau <= 0.0) {
        return am;
    }

    // 2) Solve tau = tau_of_eta(eta, a_m) for eta in [0, pi]
    //    We'll do a simple bisection on f(eta)= tau_of_eta(eta)- tau.

    double etaLo = 0.0;         // f(etaLo) = tau_of_eta(0)=0
    double etaHi = M_PI;        // f(etaHi) = (a_m/2)*(pi + sin(pi))= a_m*pi/2 => tau_collapse

    double fLo = tau_of_eta(etaLo) - tau; // ~ -tau, negative if tau>0
    double fHi = tau_of_eta(etaHi) - tau; // ~ tau_collapse - tau, should be >0 if tau< tau_collapse

    // Sanity check
    if(fLo>0.0 || fHi<0.0) {
        // This would mean something is off. Possibly user gave tau < 0 or > tau_coll
        // We'll just clamp, but you could also throw an exception.
        if(fLo>0.0) return am;      // tau<=0 => a= a_m
        if(fHi<0.0) return 0.0;      // tau> tau_coll => a=0
    }

    const int maxIters=100;
    const double tol=1e-6;

    double etaMid=0.0;
    for(int iter=0; iter<maxIters; iter++){
        etaMid = 0.5*(etaLo+etaHi);
        double fMid = tau_of_eta(etaMid) - tau;

        if(std::fabs(fMid)< tol){
            break;
        }
        // sign-based bisection
        if(fLo * fMid>0.0){
            // means fLo and fMid have same sign => shift lo
            etaLo= etaMid;
            fLo= fMid;
        } else {
            // shift hi
            etaHi= etaMid;
            fHi= fMid;
        }
    }

    // 3) having found etaMid => we compute a(\etaMid)
    return a_of_eta(etaMid);
}

double Oppenheimer_Snyder_Collapse::get_a(double tau){

    double eta = eta_from_tau_newton(tau);
    return 0.5*am*(1+std::cos(eta));

}

double Oppenheimer_Snyder_Collapse::eta_from_tau_newton(double tau) {

    double tol=1e-6;

    // This function inverts: tau = 0.5 * am * (eta + sin(eta))
    // on the domain eta in [0, pi].
    // We solve for eta using Newton’s method until |f(eta)| < tol,
    // or until a fixed max number of iterations.

    // 1) Check that tau is in the valid range.
    //    - Minimum is tau=0 at eta=0
    //    - Maximum is tau=0.5*am*(pi + sin(pi)) = 0.5 * am * pi, since sin(pi)=0
    if (tau < 0.0 || tau > 0.5 * am * M_PI) {
        throw std::runtime_error("invertTau: tau out of range [0, 0.5*am*pi].");
    }

    // 2) Initial guess:
    //    If we ignore sin(eta), then  tau ~ 0.5*am*(eta),
    //    => eta ~ 2*tau/am.  That generally puts us near the correct value.
    //    This will be between 0 and pi if tau is in [0, 0.5*am*pi].
    double eta  = 2.0 * tau / am;

    // 3) Iteration
    const int maxIter = 50;  // typically ~5-10 is plenty
    for (int i = 0; i < maxIter; i++)
    {
        // f(eta)   = 0.5*am*(eta + sin(eta)) - tau
        // f'(eta)  = 0.5*am*(1 + cos(eta))
        double f   = 0.5 * am * (eta + std::sin(eta)) - tau;
        double fp  = 0.5 * am * (1.0 + std::cos(eta));

        // Check convergence by residual
        if (std::fabs(f) < tol) {
            break;
        }

        // If derivative is extremely small near eta=pi => Newton can stall.
        // But if tau ~ 0.5*am*pi, the solution is indeed close to pi.
        // We'll do a basic safeguard to avoid dividing by zero:
        if (std::fabs(fp) < 1e-14) {
            // fallback (one step of bisection or just clamp):
            eta = M_PI; // the only solution if tau is at max
            break;
        }

        double dEta = f / fp;
        eta -= dEta;

        // Keep eta in [0, pi]
        if (eta < 0.0)    eta = 0.0;
        if (eta > M_PI)   eta = M_PI;
    }

    return eta;

}


double Oppenheimer_Snyder_Collapse::eta_from_tau(double tau) {
    // tau_of_eta(eta) = 0.5 * am * (eta + sin(eta))
    // Let tau_collapse be the proper time when eta = π.
    double tau_collapse = 0.5 * am * M_PI; // since sin(pi) = 0

    // Clamp tau to the physically allowed range.
    if(tau <= 0.0) {
        return 0.0;
    }
    if(tau >= tau_collapse) {
        return M_PI;
    }

    // Bisection bounds: for a collapsing dust cloud, eta is in [0, π].
    double etaLo = 0.0;
    double etaHi = M_PI;
    const int maxIters = 100;
    const double tol = 1e-6;
    double etaMid = 0.0;

    // We assume that tau_of_eta is defined as:
    // double Oppenheimer_Snyder_Collapse::tau_of_eta(double eta) {
    //     return 0.5 * am * (eta + sin(eta));
    // }

    // Begin bisection
    for (int iter = 0; iter < maxIters; iter++) {
        etaMid = 0.5 * (etaLo + etaHi);
        double fMid = 0.5 * am * (etaMid + sin(etaMid)) - tau;
        if (fabs(fMid) < tol) {
            break;
        }
        // Compute f at the lower bound:
        double fLo = 0.5 * am * (etaLo + sin(etaLo)) - tau;
        if (fLo * fMid > 0) {
            // Same sign: move the lower bound upward.
            etaLo = etaMid;
        } else {
            // Otherwise, move the upper bound downward.
            etaHi = etaMid;
        }
    }

    return etaMid;
}


double Oppenheimer_Snyder_Collapse::get_a_old(double tau){
	if(tau>=0){
		if(tau>m_tau[m_tau.size()-1]){
			return 0;
		}
		double i = tau/h;
		int i_m = max(int(floor(i)),0);
		int i_p = min(i_m+1,int(m_tau.size()-1));
		double s = i-i_m;
		return m_a[i_m]*(1-s)+m_a[i_p]*s;
	}
	else{
		return am;
	}
}	

double Oppenheimer_Snyder_Collapse::get_redshift(double r){

	double static_component = sqrt(fabs(1-2*M/r));
	double E0 = sqrt(1-2*M/R0);
	double dtdr = E0*(sqrt(2/M)*(M/sqrt(r)+0.5*sqrt(r))+2*M*sqrt(2*M)/((r-2*M)*sqrt(r)));
	double myBeta = fabs(1.0/dtdr);
	double myGamma = 1/sqrt(1-myBeta*myBeta);
	//double z_factor = myGamma*(1-myBeta)/static_component;
	double z_factor = 1.0/static_component;
	return fmax(1.0/(z_factor),1.0/256); // clamp to prevent round-down to zero.

}


double Oppenheimer_Snyder_Collapse::compute_energy(double p_loc[]){double r=p_loc[1];return t_dir*sqrt(1-2*M/r);}

double Oppenheimer_Snyder_Collapse::error_norm(double *p,double *x){
	//double r = p[1];
	
	return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]));
	
	if(p[6]>0){
		double r=p[1];
		return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]))/fabs(r-2*M);
	}
	else{
		return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]));
	}
	/**if(p[4]>0){
		return sqrt((x[0])*(x[0])*(1-2*M/r)+(x[1])*(x[1])/(1-2*M/r)+(x[2])*(x[2])*r*r);
	}
	else{
		return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]));
	}**/
	
}



double Oppenheimer_Snyder_Collapse::min_step_size(double* y){
	return 1.0;
//double r = y[1];
//double r_dot = y[4];
//double phi_dot = y[5];
//return 1.0/sqrt(r_dot*r_dot+r*r*phi_dot*phi_dot);

}

void Oppenheimer_Snyder_Collapse::get_null_vector_in_direction(double e1[],double e2[],double e3[],double P_obs_loc[],double V_obs_loc[],double V[],double V_loc[],double &psi){

		//cout<<"getting null vector"<<endl;
	
		double v_dot_e1 = V[0]*e1[0]+V[1]*e1[1]+V[2]*e1[2];
		double v_dot_e2 = V[0]*e2[0]+V[1]*e2[1]+V[2]*e2[2];
		double v_dot_e3 = V[0]*e3[0]+V[1]*e3[1]+V[2]*e3[2];
		
		psi =atan2(v_dot_e3,v_dot_e2);
		double epsi[3] = {cos(psi)*e2[0]+sin(psi)*e3[0],cos(psi)*e2[1]+sin(psi)*e3[1],cos(psi)*e2[2]+sin(psi)*e3[2]};
		double v_dot_epsi = V[0]*epsi[0]+V[1]*epsi[1]+V[2]*epsi[2];
		double phi = atan2(v_dot_epsi,v_dot_e1);
		
		
		//(t,r,phi)
		if(P_obs_loc[3]>0){ // external
			double r = P_obs_loc[1];
			
			V_loc[0]=t_dir/sqrt(1-2*M/r);
			V_loc[1]=-cos(phi)*sqrt(1-2*M/r);
			V_loc[2]=sin(phi)/r;
		}
		else{
            double a_of_tau;
            if(use_eta){
				cout<<"we shouldn't go in here"<<endl;
                double eta = P_obs_loc[0];
                a_of_tau=a_of_eta(eta);
            }
            else{
                double tau = P_obs_loc[0];
                double a_of_tau=get_a(tau);
            }
			double X = P_obs_loc[1];
			V_loc[0]=t_dir;
			V_loc[1]=-cos(phi)/(a_of_tau);
			V_loc[2]=sin(phi)/(a_of_tau*sin(X));
		}

}

void Oppenheimer_Snyder_Collapse::local_to_cartesian(double o[],double e1[],double e2[], double e3[], double p_local[],double p_cart[],double psi){
		// t = (v+u)/2
		// r = (v-u)/2
		//double R;
		//if(p_local[0]>0){
		//	R = get_r(p_local[0],p_local[1]);
		//}
		//else{
		//	R = (p_local[1]-p_local[0])/2.0; // assume flat space.
		//}
	
		//double weight = (1+tanh(r0-R))/2;
		//weight = 0;
		//R = weight*R+(1-weight)*(p_local[1]-p_local[0])/2.0;
		//double phi = p_local[2];
		//double x = R*cos(phi);
		//double y = R*sin(phi);
		
		double r,phi,x,y;
		
		if(p_local[3]>0){ //external		
			r = p_local[1];
			phi = p_local[2];	
		}
		else{ //internal
            double a_of_tau;
            if(use_eta){
				cout<<"we shouldn't go in here"<<endl;
                double eta = p_local[0];
                a_of_tau = a_of_eta(eta);
            }
            else{
                double tau=p_local[0];
                a_of_tau=get_a(tau);
            }
			double X = p_local[1];
			phi = p_local[2];
			r = a_of_tau*sin(X);
		}
		x = r*cos(phi);
		y = r*sin(phi);	
		
		
		
		double epsi[3] = {cos(psi)*e2[0]+sin(psi)*e3[0],cos(psi)*e2[1]+sin(psi)*e3[1],cos(psi)*e2[2]+sin(psi)*e3[2]};
		
		p_cart[0]=o[0]-x*e1[0]+y*epsi[0];
		p_cart[1]=o[1]-x*e1[1]+y*epsi[1];
		p_cart[2]=o[2]-x*e1[2]+y*epsi[2];
		
	}

bool Oppenheimer_Snyder_Collapse::accel_1st_order(double*& yp,double* y, double E){
	
	// y and yp are giant vectors of the form y=[p,v] where p is position and v is velocity.
	
	double MM = M;
	if(y[0]<0){
		MM=0; // black hole is turned off for times less than zero.
	}
	
	if(y[6]>0){ // encodes outside.
	
	// position derivative is velocity.
	yp[0]=y[3]; 
	yp[1]=y[4]; 
	yp[2]=y[5]; 

	// Debug check: if \dot{t} is about to become positive, print and stop
	if(yp[0] > 0.0)
	{
		// Print with high precision
		std::cerr << std::setprecision(16) 
				  << "ERROR: dot{t} is positive!\n"
				  << "r = " << y[1] << "\n"
				  << "t = " << y[0] << "\n"
				  << "r_dot = " << y[4] << "\n"
				  << "phi_dot = " << y[5] << "\n"
				  << "dot{t} = " << y[3] << std::endl;

		// Abort the program
		std::abort();  
		// or std::exit(EXIT_FAILURE);
		// or throw std::runtime_error("dot{t} is positive!")
	}

	double r = y[1];

	yp[3]=-2*MM*E*y[4]/((y[1]-2*MM)*(y[1]-2*MM));
	yp[4]=-MM*E*E/(r*(r-2*MM))+MM*y[4]*y[4]/(r*(r-2*MM))+(r-2*MM)*y[5]*y[5];
	yp[5]=-2*y[4]*y[5]/y[1];
	
	}
    else if (y[6] < 0 && !use_eta) {
        // Our state in τ–coordinates:
        // y[0] = τ, y[1] = X, y[2] = φ,
        // y[3] = τ̇, y[4] = Ẋ, y[5] = φ̇.
        ///cout<<"we shouldn't go in here"<<endl;
        double tau = y[0];
        double X   = y[1];
        double phi = y[2];
        double tau_dot = y[3];
        double X_dot   = y[4];
        double phi_dot = y[5];
        
        // Step 1: Convert τ to conformal time η.
        double eta = eta_from_tau(tau);  // You must implement this inversion.
        
        // Step 2: Compute the conformal scale factor and its derivative.
        double a_eta = 0.5 * am * (1.0 + cos(eta));     // a(η)
        double a_eta_prime = -0.5 * am * sin(eta);         // derivative w.r.t. η.
        
        // Step 3: Compute η̇ from τ̇:
        double eta_dot = tau_dot / a_eta;
		//double eta_dot = -1;
        
        // Step 4: Compute the acceleration in the conformal system.
        // For η (conformal time) we assume:
        double eta_ddot = 0;  // Enforce zero conformal time acceleration.
        // For spatial coordinates, use the conformal geodesic equations:
        double X_ddot = sin(X) * cos(X) * (phi_dot * phi_dot);
        double phi_ddot = -2.0 * (1.0 / tan(X)) * X_dot * phi_dot;
        
        // Step 5: Convert the conformal time acceleration back to the τ system.
        double tau_ddot = a_eta * eta_ddot + a_eta_prime * (eta_dot * eta_dot);
        
        // Step 6: Fill in the derivative vector.
        // Copy the velocities (positions' derivatives):
        yp[0] = tau_dot;  // (or equivalently, y[3])
        yp[1] = X_dot;    // y[4]
        yp[2] = phi_dot;  // y[5]
        
        // Set the accelerations:
        yp[3] = tau_ddot;
        yp[4] = X_ddot;
        yp[5] = phi_ddot;
    }
	else if(y[6]<0 && false){ // old form
		
		
		
		// position derivative is velocity.
		yp[0]=y[3]; 
		yp[1]=y[4]; 
		yp[2]=y[5]; 
		
		double tau = y[0];
		
		double a = get_a(tau);
		double a_dot = -sqrt(C/a-1);
		
		// we are inside the dust cloud.
		//tau dot dot = - a(tau)*a'(tau)*(X_dot^2+sin(X)^2*phi_dot^2)
		yp[3]=-a*a_dot*(y[4]*y[4]+sin(y[1])*sin(y[1])*y[5]*y[5]);
		//X dot dot = sin(X)*cos(X)*phi_dot^2-2*(a_dot/a)*tau_dot*X_dot
		yp[4]=sin(y[1])*cos(y[1])*y[5]*y[5]-2*(a_dot/a)*y[3]*y[4];
		//phi dot dot = -2cot(X)*X_dot*phi_dot-2(a/a_dot)*tau_dot*phi_dot
		yp[5]=-2*(1.0/tan(y[1]))*y[4]*y[5]-2*(a_dot/a)*y[3]*y[5];
	}
    else if(y[6]<0 && use_eta){
		
		cout<<"we shouldn't go in here"<<endl;
		
		// position derivative is velocity.
		yp[0]=-1; 
		yp[1]=y[4]; 
		yp[2]=y[5]; 

		
		// we are inside the dust cloud.
		yp[3]=0;
		//X dot dot = sin(X)*cos(X)*phi_dot^2
		yp[4]=sin(y[1])*cos(y[1])*y[5]*y[5];
		//phi dot dot = -2cot(X)*X_dot*phi_dot
		yp[5]=-2*(1.0/tan(y[1]))*y[4]*y[5];
	}
	else{
		cout<<"invalid y[6]!, y = "<<y[0]<<" "<<y[1]<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<" "<<y[5]<<" "<<y[6]<<endl;
	}


	return inside_universe(y);
	
}

int Oppenheimer_Snyder_Collapse::inside_universe(double *p){
	
	return 1;
	
	/**if(p[3]>0){ //external
	
		if( p[1] > 2*M){
			return 1;
		}
		else{
			return 0;
		}
	}
	else{ //internal
		if(get_a(p[0])>0){return 1;}
		else{return 0;}
	}**/
	
}

// this is just the radial infall trajectory of a particle starting
// at r=R0.  See https://www.reed.edu/physics/courses/Physics411/html/411/page2/files/Lecture.31.pdfhttps://www.reed.edu/physics/courses/Physics411/html/411/page2/files/Lecture.31.pdf
double Oppenheimer_Snyder_Collapse::dust_time(double r){

	double E0 = sqrt(1-2*M/R0);
	
	// Compute terms in the t(r) formula
    double sqrt_r = sqrt(r);
    double sqrt_R0 = sqrt(R0);

    double term1 = (E0/3)*sqrt(2/M)*( sqrt_r * (6 * M + r) - sqrt_R0 * (6 * M + R0) );

    double log_term_r = log((1 - sqrt(2 * M / r)) / (1 + sqrt(2 * M / r)));
    double log_term_R0 = log((1 - sqrt(2 * M / R0)) / (1 + sqrt(2 * M / R0)));

    double term2 = 2 * M * E0 * (log_term_r - log_term_R0);

    // Combine terms
    double t_r = term1 + term2;
	
	t_r = -t_r;

    return t_r;

}

bool Oppenheimer_Snyder_Collapse::inside_dust_cloud(double r,double t){

	double t_dust = dust_time(r);
    return (t_dust > t);

}

bool Oppenheimer_Snyder_Collapse::exiting_dust_cloud(double X,double Xdot){

	if(X>=X0 and Xdot>0){
		return true;
	}
	else{
		return false;
	}

}

double Oppenheimer_Snyder_Collapse::findDustExitCrossingBisection(
    double y[4],    // old step: y[0]=tau,y[1]=X
    double yp[3],   // old step velocities: yp[0]=tau_dot, yp[1]=X_dot
    double y_new[4],// new step
    double yp_new[3],
    double tol,
    int maxIters)
{
	
    // Evaluate exit condition at start:
    bool exitOld = exiting_dust_cloud(y[1], yp[1]);
    // Evaluate exit condition at end:
    bool exitNew = exiting_dust_cloud(y_new[1], yp_new[1]);

    // If no bracket, return -1.0
    if(exitOld == exitNew){
        return -1.0;
    }

    // define f(s)= +1 if exiting_dust_cloud(...)=true, else -1
    auto f = [&](double s){
        // linear interpolation for X
        double X_s = (1.0 - s)*y[1] + s*y_new[1];
        // linear interpolation for X_dot
        double Xdot_s = (1.0 - s)*yp[1] + s*yp_new[1];

        bool e = exiting_dust_cloud(X_s, Xdot_s);
        return e? +1.0 : -1.0;
    };

    double sLo = 0.0;
    double sHi = 1.0;
    double fLo = f(sLo); // should be +1 or -1
    double fHi = f(sHi);

    // we require fLo*fHi<0
    // but we do know exitOld != exitNew => so fLo * fHi = -1
    // let's do a simple bisection for s in [0,1]
    for(int iter=0; iter<maxIters; iter++){
        double mid = 0.5*(sLo + sHi);
        double fMid = f(mid);

        if(fabs(sHi - sLo) < tol){
            return mid;
        }

        // sign-based approach
        if(fLo * fMid > 0.0){
            sLo = mid;
            fLo = fMid;
        } else {
            sHi = mid;
            fHi = fMid;
        }
    }

    // if not converged, return midpoint
    return 0.5*(sLo + sHi);
}

double Oppenheimer_Snyder_Collapse::findDustBoundaryCrossingBisection(
    double y[4],      // old step: y[0]=t, y[1]=r
	double yp[4],
    double y_new[4],  // new step
	double yp_new[4],
    double tol,
    int maxIters
)
{
	double rdot = yp[1];
	double rdot_new = yp_new[1];
//if(rdot>0 || rdot_new>0 || y[0]<0 || y_new[0]<0){return -1;}
    // Step 1) compute function g(0) = dust_time(r(0)) - t(0)
    double r0 = y[1];
    double t0 = y[0];
    double f0 = dust_time(r0) - t0;
	//cout<<"r0, dust_time(r0), t0 = "<<r0<<" "<<dust_time(r0)<<" "<<t0<<endl;

    // Step 2) compute g(1) = dust_time(r(1)) - t(1)
    double r1 = y_new[1];
    double t1 = y_new[0];
    double f1 = dust_time(r1) - t1;
	//cout<<"r1, dust_time(r1), t1 = "<<r1<<" "<<dust_time(r1)<<" "<<t1<<endl;

    // If no sign change => no crossing => return -1
    if(f0*f1>0.0) {
        return -1.0;  
    }

    // standard bisection on [0,1]
    double sLo=0.0, sHi=1.0;
    double fLo=f0, fHi=f1;

    for(int iter=0; iter<maxIters; iter++){
        double mid = 0.5*(sLo + sHi);

        // Interpolate to get r(mid), t(mid)
        double t_mid = (1.0-mid)*t0   + mid*t1;
        double r_mid = (1.0-mid)*r0   + mid*r1;

        double fMid = dust_time(r_mid) - t_mid;

        // check if done
        // (we do a double approach: either sign-based or we see if interval < tol)
        if(std::fabs(sHi - sLo)< tol || std::fabs(fMid)<1e-14) {
			double r_mid = r0*(1-mid)+r1*mid;
			double t_mid =t0*(1-mid)+t1*mid;
			//cout<<"r_mid, t_mid, dust_time"<<r_mid<<" "<<t_mid<<" "<<dust_time(r_mid)<<endl;
			
            return mid;
        }

        // Bisection sign check
        if(fLo * fMid > 0.0) {
            sLo=mid;
            fLo=fMid;
        } else {
            sHi=mid;
            fHi=fMid;
        }
    }
	//cout<<"not converged"<<endl;
	double mid = 0.5*(sLo + sHi);
	double r_mid = r0*(1-mid)+r1*mid;
	double t_mid =t0*(1-mid)+t1*mid;
	//cout<<"r_mid, t_mid, dust_time"<<r_mid<<" "<<t_mid<<" "<<dust_time(r_mid)<<endl;

    // not converged => return midpoint
    return 0.5*(sLo + sHi);
}


double Oppenheimer_Snyder_Collapse::findDustBoundaryCrossing(double y[4],
                                                             double y_new[4],
                                                             double tol,
                                                             int maxIters)
{
    // 
    // 1) Check old step inside/outside status:
    //
    bool insideOld = inside_dust_cloud(y[1], y[0]);       // y[1]=r, y[0]=t
    bool insideNew = inside_dust_cloud(y_new[1],y_new[0]);

    // If there's no crossing, we can either throw or return -1.0 or something
    if(insideOld == insideNew){
        // no sign change => no crossing
        // you can handle it as you see fit:
        return -1.0; 
    }

    // We do a bisection in "s" from 0 to 1:
    // define f(s)= +1 if inside_dust_cloud(...) =true, -1 otherwise
    // we want f(s)=0 => crossing. We'll do sign approach.

    // Helper lambda for f(s)
    auto f = [&](double s){
        // linearly interpolate in [0,1]
        double t_s = (1.0 - s)*y[0]     + s*y_new[0];
        double r_s = (1.0 - s)*y[1]     + s*y_new[1];
        
        bool inside = inside_dust_cloud(r_s, t_s);
        return inside ? +1.0 : -1.0;
    };

    double sLo=0.0, sHi=1.0;
    double fLo = f(sLo);
    double fHi = f(sHi);

    // we know fLo * fHi <0.0 because inside vs outside are different
    for(int iter=0; iter< maxIters; iter++){
        double mid = 0.5*(sLo + sHi);
        double fMid = f(mid);
        
        // If we are close enough, or no sign
        // However "fMid" is discrete +1/-1. So let's see if
        // it maybe never is "0." We want sign to flip from +1 to -1. 
        // We'll do a tolerance approach in "s."

        if(std::fabs(sHi - sLo) < tol){
            return mid;
        }

        // Bisection step:
        if(fLo * fMid > 0.0){
            sLo = mid;
            fLo = fMid;
        } else {
            sHi = mid;
            fHi = fMid;
        }
    }
    // if we didn't converge, return midpoint
    return 0.5*(sLo + sHi);
}

double Oppenheimer_Snyder_Collapse::findTauZeroCrossing(double y[4],
                                                        double y_new[4],
                                                        double tol,
                                                        int maxIters)
{
    // We assume that, at the start of the step, y[0] >= 0,
    // and at the end of the step, y_new[0] < 0.
    // We want to find s in [0,1] such that interpolated tau(s) = 0.

    double tau0 = y[0];       // >= 0
    double tau1 = y_new[0];   // < 0

    // Quick check if there's no sign change in tau:
    if(tau0 < 0.0 || tau1 >= 0.0)
    {
        // No bracket => no valid crossing => return sentinel
        return -1.0;
    }

    // Define f(s) = tau(s) = (1-s)*tau0 + s*tau1
    // We want f(s) = 0 => crossing time fraction
    auto f = [&](double s){
        return (1.0 - s)*tau0 + s*tau1; 
    };

    // Evaluate at endpoints
    double fLo = f(0.0); // = tau0 >= 0
    double fHi = f(1.0); // = tau1 < 0

    // Ensure there's indeed a bracket (fLo >= 0, fHi < 0)
    if(fLo*fHi > 0.0) 
    {
        // No sign change => no crossing
        return -1.0;
    }

    // Bisection in [0,1]
    double sLo = 0.0;
    double sHi = 1.0;
    for(int iter = 0; iter < maxIters; ++iter)
    {
        double sMid = 0.5*(sLo + sHi);
        double fMid = f(sMid);

        // If the interval is small enough or f(sMid) is very close to 0
        if(std::fabs(sHi - sLo) < tol || std::fabs(fMid) < 1e-15)
        {
            return sMid;
        }

        // Sign-based bisection
        if(fLo * fMid > 0.0)
        {
            sLo = sMid;
            fLo = fMid;
        }
        else
        {
            sHi = sMid;
            fHi = fMid;
        }
    }

    // Not fully converged => return midpoint
    return 0.5*(sLo + sHi);
}

double Oppenheimer_Snyder_Collapse::findEtaZeroCrossing(double y[4],
                                                        double y_new[4],
                                                        double tol,
                                                        int maxIters)
{
    // We assume that, at the start of the step, y[0] = eta >= 0,
    // and at the end of the step, y_new[0] = eta < 0.
    // We want to find s in [0,1] such that the interpolated eta(s) = 0.
    
    double eta0 = y[0];      // starting eta, should be >= 0
    double eta1 = y_new[0];  // ending eta, should be < 0

    // Quick check: if there's no sign change in eta, return -1 (no valid crossing).
    if(eta0 < 0.0 || eta1 >= 0.0)
    {
        return -1.0;
    }

    // Define f(s) as the linear interpolation between eta0 and eta1.
    auto f = [&](double s) {
        return (1.0 - s) * eta0 + s * eta1;
    };

    double fLo = f(0.0);  // f(0) = eta0 (>= 0)
    double fHi = f(1.0);  // f(1) = eta1 (< 0)

    // Ensure there is a bracket (fLo >= 0 and fHi < 0)
    if (fLo * fHi > 0.0)
    {
        return -1.0;
    }

    // Bisection in the interval [0, 1] to find s such that f(s)=0.
    double sLo = 0.0;
    double sHi = 1.0;
    for (int iter = 0; iter < maxIters; ++iter)
    {
        double sMid = 0.5 * (sLo + sHi);
        double fMid = f(sMid);
        
        // If the interval is small enough or f(sMid) is extremely close to 0, return sMid.
        if (std::fabs(sHi - sLo) < tol || std::fabs(fMid) < 1e-15)
        {
            return sMid;
        }
        
        // Use the sign of f(s) to update the bracket.
        if (fLo * fMid > 0.0)
        {
            sLo = sMid;
            fLo = fMid;
        }
        else
        {
            sHi = sMid;
            fHi = fMid;
        }
    }
    
    // If not fully converged, return the midpoint.
    return 0.5 * (sLo + sHi);
}


void Oppenheimer_Snyder_Collapse::to_interior_coordinates(double y[],double yp[],double y_new[],double yp_new[])
{
	double r = y[1];
	double r_dot = yp[1];
	double phi_dot = yp[2];
    double X = X0;

	//cout<<"X0 =  "<<X0<<endl;
	//cout<<"r = "<<r<<endl;
	//double arg = fmax(fmin((2.0 / am) * (r / sin(X0)) - 1.0,1.0),-1.0);
	double arg = (2.0 / am) * (r / sin(X0))-1.0;
	//cout<<"arg = "<<arg;

    // 1) Invert boundary condition to find eta
    //double eta = acos((2.0 / am) * (r / sin(X0)) - 1.0);
	double eta = acos(arg);

	//cout<<"(2.0 / am) * (r / sin(X0)) - 1.0 = "<<(2.0 / am) * (r / sin(X0)) - 1.0<<endl;

	//cout<<"eta =  "<<eta<<endl;

	double tau = 0.5*am*(eta+sin(eta));

    // 2) a(tau), da/dtau
    double a_of_tau = 0.5 * am * (1.0 + cos(eta));
	
	//cout<<"a_of_tau =  "<<a_of_tau<<endl;
	
	r = a_of_tau*sin(X0);

	//cout<<"C/a_of_tau-1 =  "<<C/a_of_tau-1<<endl;
	

	double dadtau = -sqrt(C/a_of_tau-1);

    // 3) A, B
    double A = (dadtau * sin(X0)) / (a_of_tau * cos(X0));
    double B = r_dot / (a_of_tau * cos(X0));

    // 4) Quadratic coefficients alpha2, alpha1, alpha0
    //    Corrected alpha0 => plus sign, include a_of_tau^2 for phi_dot^2 term
    double alpha2 = -1.0 + (a_of_tau * a_of_tau) * (A * A);
    double alpha1 = -2.0 * (a_of_tau * a_of_tau) * A * B;
    double alpha0 = (B * B)
                    + (a_of_tau * a_of_tau) * (sin(X0) * sin(X0))
                      * (phi_dot * phi_dot);

    // 5) Solve the quadratic
    double rad = alpha1 * alpha1 - 4.0 * alpha2 * alpha0;
    // Check rad >= 0 to avoid sqrt of negative
    if (rad < 0.0) {
		cout<<"rad <0, rad =  "<<rad<<endl;
        // No real solutions => possibly no crossing
        return;
    }
	//cout<<"rad =  "<<rad<<endl;
    double sqrt_rad = sqrt(rad);

    double tau_dot_p = (-alpha1 + sqrt_rad) / (2.0 * alpha2);
    double tau_dot_m = (-alpha1 - sqrt_rad) / (2.0 * alpha2);

    // 6) Compute the corresponding X_dot for each root
    double X_dot_p = (r_dot - dadtau * sin(X0) * tau_dot_p) 
                     / (a_of_tau * cos(X0));
    double X_dot_m = (r_dot - dadtau * sin(X0) * tau_dot_m) 
                     / (a_of_tau * cos(X0));

    // 7) Evaluate validity of each root
    bool p_root_valid = false;
    bool m_root_valid = false;

    // Suppose we want inward crossing => X_dot<0
    // If time_dir=+1 => we want tau_dot>0 ; if time_dir=-1 => tau_dot<0
    // So "tau_dot_p * time_dir" should be >0, not ==1
    if (X_dot_p < 0.0 && (tau_dot_p * t_dir > 0.0)) {
        p_root_valid = true;
    }
    if (X_dot_m < 0.0 && (tau_dot_m * t_dir > 0.0)) {
        m_root_valid = true;
    }

	double tau_dot,X_dot;
    // 8) Decide
    if (p_root_valid && m_root_valid) {
		cout<<"both roots valid"<<endl;
		cout<<"X_dot_p, X_dot_m = "<<X_dot_p<<" "<<X_dot_m<<endl;
		cout<<"tau_dot_p, tau_dot_m = "<<tau_dot_p<<" "<<tau_dot_m<<endl;
        // Could be 2 possible solutions => handle how you wish
        // Possibly pick the one with largest magnitude? or something else
        // ...
    } 
    else if (p_root_valid) {
        tau_dot = tau_dot_p;
        X_dot   = X_dot_p;
    }
    else if (m_root_valid) {
        tau_dot = tau_dot_m;
        X_dot   = X_dot_m;
    }
    else {
		//cout<<"neither roots valid"<<endl;
		//cout<<"X_dot_p, X_dot_m r_dot = "<<X_dot_p<<" "<<X_dot_m<<" "<<r_dot<<endl;
		//cout<<"tau_dot_p, tau_dot_m = "<<tau_dot_p<<" "<<tau_dot_m<<endl;
        // Neither root valid => maybe grazing or no crossing 
		// we'll assume it's a mistake and cancel going in.
		y_new[0]=y[0];
		y_new[1]=y[1];
		y_new[2]=y[2];
		y_new[3]=fabs(y[3]); // inside
		
		yp_new[0]=yp[0];
		yp_new[1]=yp[1];
		yp_new[2]=yp[2];
        return;
    }

	//cout<<"tau = "<<tau<<" X = "<<X<<" y[2]= "<<y[2]<<endl;

    y_new[0]=tau;
	y_new[1]=X;
	y_new[2]=y[2];
	y_new[3]=-fabs(y[3]); // inside
	
	yp_new[0]=tau_dot;
	yp_new[1]=X_dot;
	yp_new[2]=yp[2];
	
}

void Oppenheimer_Snyder_Collapse::to_interior_coordinates_eta(double y[], double yp[],
                                                               double y_new[], double yp_new[])
{
    // Extract exterior data.
    double r         = y[1];
    double r_dot     = yp[1];
    double phi_dot   = yp[2];
    double X         = X0;  // Boundary value (assumed constant)

    // Compute conformal time η from the boundary condition.
    double arg = (2.0 / am) * (r / sin(X0)) - 1.0;
    double eta = acos(arg);  // Interior time coordinate in conformal form

    // Compute proper time τ from η.
    double tau = 0.5 * am * (eta + sin(eta));
    
    // Compute the scale factor a(τ) = 0.5*am*(1 + cos(η)).
    double a_of_tau = 0.5 * am * (1.0 + cos(eta));
    //cout<<" (going in) eta = "<<eta<<endl;
    //cout << "a_of_tau = " << a_of_tau << endl;

    // Reconstruct interior r.
    r = a_of_tau * sin(X0);

    // Compute derivative of a with respect to τ.
    double dadtau = -sqrt(C / a_of_tau - 1);

    // Coefficients from the relation: r_dot = (dadtau * sin(X0))*τ̇ + a(τ)*cos(X0)*Ẋ.
    double commonDenom = a_of_tau * cos(X0);
    double A = (dadtau * sin(X0)) / commonDenom;
    double B = r_dot / commonDenom;

    // Set up quadratic for τ̇: alpha2 * (τ̇)² + alpha1 * τ̇ + alpha0 = 0.
    double alpha2 = -1.0 + (a_of_tau * a_of_tau) * (A * A);
    double alpha1 = -2.0 * (a_of_tau * a_of_tau) * A * B;
    double alpha0 = (B * B) + (a_of_tau * a_of_tau) * (sin(X0) * sin(X0)) * (phi_dot * phi_dot);

    double rad = alpha1 * alpha1 - 4.0 * alpha2 * alpha0;
    if (rad < 0.0) {
        cout << "rad < 0, rad = " << rad << endl;
        // No real solutions: possibly no crossing.
        return;
    }
    double sqrt_rad = sqrt(rad);

    double tau_dot_p = (-alpha1 + sqrt_rad) / (2.0 * alpha2);
    double tau_dot_m = (-alpha1 - sqrt_rad) / (2.0 * alpha2);

    // Compute corresponding Ẋ for each τ̇ root.
    double X_dot_p = (r_dot - dadtau * sin(X0) * tau_dot_p) / commonDenom;
    double X_dot_m = (r_dot - dadtau * sin(X0) * tau_dot_m) / commonDenom;

    // Convert τ̇ to conformal time derivative: η̇ = τ̇ / a(τ).
    double eta_dot_p = tau_dot_p / a_of_tau;
    double eta_dot_m = tau_dot_m / a_of_tau;

    // Validate the roots: we want inward crossing (Ẋ < 0) and proper sign for η̇.
    bool p_root_valid = (X_dot_p < 0.0 && (eta_dot_p * t_dir > 0.0));
    bool m_root_valid = (X_dot_m < 0.0 && (eta_dot_m * t_dir > 0.0));

    double tau_dot, X_dot, eta_dot;
    if (p_root_valid && m_root_valid) {
        cout << "both roots valid" << endl;
        cout << "X_dot_p, X_dot_m = " << X_dot_p << " " << X_dot_m << endl;
        cout << "tau_dot_p, tau_dot_m = " << tau_dot_p << " " << tau_dot_m << endl;
        // Choose the one with larger |η̇|.
        if (fabs(eta_dot_p) > fabs(eta_dot_m)) {
            tau_dot = tau_dot_p;
            X_dot   = X_dot_p;
            eta_dot = eta_dot_p;
        } else {
            tau_dot = tau_dot_m;
            X_dot   = X_dot_m;
            eta_dot = eta_dot_m;
        }
    }
    else if (p_root_valid) {
        tau_dot = tau_dot_p;
        X_dot   = X_dot_p;
        eta_dot = eta_dot_p;
    }
    else if (m_root_valid) {
        tau_dot = tau_dot_m;
        X_dot   = X_dot_m;
        eta_dot = eta_dot_m;
    }
    else {
        // No valid crossing: return the unchanged state.
        y_new[0] = y[0];
        y_new[1] = y[1];
        y_new[2] = y[2];
        y_new[3] = 1; // inside flag
        yp_new[0] = yp[0];
        yp_new[1] = yp[1];
        yp_new[2] = yp[2];
        return;
    }

    // Now store the interior coordinates as (η, X, φ).
    y_new[0] = eta;   // Interior time is now conformal time.
    y_new[1] = X;
    y_new[2] = y[2];
    y_new[3] = -1;    // inside flag.

    // For numerical stability we define the conformal time derivative to be -1.
    eta_dot = -1;
    // Convert the derivative: note that we keep the computed X_dot and φ_dot.
    yp_new[0] = eta_dot;   // dη/dλ = -1.
    yp_new[1] = X_dot;
    yp_new[2] = yp[2];
}

void Oppenheimer_Snyder_Collapse::to_exterior_coordinates(double y[],double yp[],double y_new[],double yp_new[])
{
	// 1) Invert boundary condition to find eta
	
	double tau = y[0];
	double tau_dot = yp[0];
	double X_dot = yp[1];
	double phi_dot = yp[2];
	double a_of_tau = get_a(tau);
	double dadtau=-sqrt(C/a_of_tau-1);
	
	double r = a_of_tau*sin(X0);
	double t = dust_time(r);
    
	
	double r_dot = dadtau*sin(X0)*tau_dot+a_of_tau*cos(X0)*X_dot;
	double t_dot = t_dir*sqrt(r_dot*r_dot/(1-2*M/r)+r*r*phi_dot*phi_dot)/sqrt(1-2*M/r);
	
	y_new[0]=t;
	y_new[1]=r;
	y_new[2]=y[2];
	y_new[3]=fabs(y[3]); //outside
	
	yp_new[0]=t_dot;
	yp_new[1]=r_dot;
	yp_new[2]=yp[2];
	
}

void Oppenheimer_Snyder_Collapse::to_exterior_coordinates_eta(double y[], double yp[],
                                                            double y_new[], double yp_new[])
{
    // --- Read the interior state in conformal coordinates ---
    // y[0] = η, y[1] = X, y[2] = φ,
    // y[3] = dη/dλ (should be -1 by our convention),
    // y[4] = dX/dλ, y[5] = dφ/dλ.
    double eta     = y[0];      // interior conformal time
    double X       = y[1];      // comoving radial coordinate (should equal X0 at the boundary)
    double phi     = y[2];
    
    double eta_dot = yp[0];      // by our convention, eta_dot should be -1.
    double X_dot   = yp[1];
    double phi_dot = yp[2];
    
    // --- Convert conformal time to proper time ---
    // Using: τ = 0.5 * am * (η + sin(η))
    double tau = 0.5 * am * (eta + sin(eta));
    
    // --- Compute the scale factor ---
    // Either use your function get_a(tau) or compute it directly.
    // Here we compute it directly as: a(τ) = 0.5 * am * (1 + cos(η))
    //cout<<" (going out) eta = "<<eta<<endl;
    double a_of_tau = 0.5 * am * (1.0 + cos(eta));
    //cout << "a_of_tau = " << a_of_tau << endl;
    
    // --- Compute the interior r coordinate ---
    // The interior coordinate relation is r = a(τ) * sin(X₀)
    double r = a_of_tau * sin(X0);
    
    // --- Compute the exterior time coordinate ---
    // For OS collapse, you use a dust time function, e.g.:
    double t = dust_time(r);
    
    // --- Compute τ̇ from the conversion ---
    // Since τ = τ(η) with dτ/dη = a(τ), and we force dη/dλ = eta_dot = -1,
    // we have: τ̇ = a(τ) * eta_dot. Thus:
    double tau_dot = a_of_tau * eta_dot;  // should be -a_of_tau.
    
    // --- Compute dadtau from the Friedmann equation ---
    double dadtau = -sqrt(C / a_of_tau - 1);
    
    // --- Compute ṙ using the chain rule ---
    // The relation used is:
    //   r_dot = (dadtau * sin(X0)) * τ̇ + a(τ) * cos(X0) * X_dot.
    double r_dot = dadtau * sin(X0) * tau_dot + a_of_tau * cos(X0) * X_dot;
    
    // --- Compute ṫ in the exterior coordinates ---
    // As before, we use the formula:
    //   t_dot = t_dir * sqrt( r_dot^2/(1-2*M/r) + r^2 * φ̇^2 ) / sqrt(1-2*M/r)
    double t_dot = t_dir * sqrt( (r_dot * r_dot)/(1 - 2 * M / r) + r * r * (phi_dot * phi_dot) )
                   / sqrt(1 - 2 * M / r);
    
    // --- Store the exterior coordinates and their derivatives ---
    y_new[0] = t;      // exterior time t
    y_new[1] = r;      // exterior radial coordinate r
    y_new[2] = phi;    // φ remains the same
    y_new[3] = 1;      // flag: outside
    
    yp_new[0] = t_dot;
    yp_new[1] = r_dot;
    yp_new[2] = phi_dot;
}

void Oppenheimer_Snyder_Collapse::to_exterior_coordinates_through_tau_equals_zero(double y[],double yp[],double y_new[],double yp_new[])
{
	// 1) Invert boundary condition to find eta
	
	double tau = y[0];
	double X = y[1];
	double tau_dot = yp[0];
	double X_dot = yp[1];
	double phi_dot = yp[2];
	double a_of_tau = am;
	double dadtau=0;
	
	double r = a_of_tau*sin(X);
	double t = -1e-6; // slightly below zero.
    
	
	double r_dot = a_of_tau*cos(X)*X_dot;
	double t_dot = t_dir*sqrt(r_dot*r_dot/(1-2*M/r)+r*r*phi_dot*phi_dot)/sqrt(1-2*M/r);
	
	y_new[0]=t;
	y_new[1]=r;
	y_new[2]=y[2];
	y_new[3]=fabs(y[3]); //outside
	
	yp_new[0]=t_dot;
	yp_new[1]=r_dot;
	yp_new[2]=yp[2];
	
}

void Oppenheimer_Snyder_Collapse::to_exterior_coordinates_through_tau_equals_zero_eta(double y[],
                                                                                   double yp[],
                                                                                   double y_new[],
                                                                                   double yp_new[])
{
    // --- Read the interior state in conformal coordinates ---
    // y[0] = eta, y[1] = X, y[2] = phi.
    // yp[0] = d(eta)/d(lambda) (typically -1), yp[1] = dX/d(lambda), yp[2] = d(phi)/d(lambda).
    double eta     = y[0];  // near 0 (since tau ~ 0)
    double X       = y[1];  // the comoving coordinate of the geodesic (may be interior, i.e. X < X0)
    double phi     = y[2];
    
    double eta_dot = yp[0]; // typically -1
    double X_dot   = yp[1];
    double phi_dot = yp[2];
    
    // --- Convert conformal time to proper time ---
    // tau = 0.5 * am * (eta + sin(eta))
    double tau = 0.5 * am * (eta + sin(eta));  // at eta=0, tau = 0.
    
    // --- Compute the scale factor ---
    // a(tau) = 0.5 * am * (1 + cos(eta))
    // For eta near 0, cos(eta) ~ 1, so a ~ am.
    double a_of_tau = 0.5 * am * (1.0 + cos(eta));
    
    // At the transition (tau = 0), we assume the dust cloud ceases, so set da/dtau = 0.
    double dadtau = 0;
    
    // --- Compute the exterior radial coordinate ---
    // Instead of forcing the boundary (r = a * sin(X0)), use the actual interior comoving coordinate:
    double r = a_of_tau * sin(X);   // r = am * sin(X) when a_of_tau = am.
    
    // --- Choose an exterior time coordinate ---
    // Set t to a small negative value to mark the exit.
    double t = -1e-6;  // slightly below zero.
    
    // --- Compute r_dot ---
    // Since r = a * sin(X) and a is constant (with da/dtau = 0),
    // we have: r_dot = a * cos(X) * X_dot.
    double r_dot = a_of_tau * cos(X) * X_dot;
    
    // --- Compute t_dot in the exterior coordinates ---
    // Use the Schwarzschild metric relation:
    // t_dot = t_dir * sqrt( r_dot^2/(1-2*M/r) + r^2*(phi_dot)^2 ) / sqrt(1-2*M/r)
    double t_dot = t_dir * sqrt( (r_dot * r_dot) / (1 - 2 * M / r) + r * r * (phi_dot * phi_dot) )
                   / sqrt(1 - 2 * M / r);
    
    // --- Store the exterior coordinates and derivatives ---
    y_new[0] = t;     // exterior time
    y_new[1] = r;     // exterior radial coordinate
    y_new[2] = phi;   // phi is unchanged
    y_new[3] = 1;     // flag: outside
    
    yp_new[0] = t_dot;
    yp_new[1] = r_dot;
    yp_new[2] = phi_dot;
}