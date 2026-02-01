#ifndef ADAPTIVE_RUNGE_KUTTA_H
#define ADAPTIVE_RUNGE_KUTTA_H


// routines for solving the system of ODEs y'=f(y).

// in general, f may not be defined everywhere.  The "adaptive step" routine returns a bool telling you whether the given step went outside of Domain(f).
// true means everything is ok, false means at some point during the step we attempted to go outside of Domain(f).

namespace Adaptive_Runge_Kutta
{
	// f is encoded in the function F, as is a function "norm" which measures the length of vector v at point y (the standard norm may not be appropriate in
	// GR applications).
	
	// it is assumed that F is of the form z=F(y',y) where if y is in Domain(f) then on exit z=true and y'=f(y), whereas if y is not in the domain(f), then
	// on exit y' is unchanged and z=false - we call the latter case an 'illegal' call to F
	
	// when Adaptive_Step exits, if no illegal calls were made then y is changed to y_new after one step, h is changed to the recommended new h,
	// and the function returns true.  Note that y_new=y if the estimated truncation error is too high (in this case h_new is smaller than h in compensation)
	// the bool y_changed tells us on exit whether or not y was modified.
	// note that the single step is made with the ORIGINAL value of h, and not the value it is changed too internally, which is a recommendation for the next step
	// otherwise, no parameters are changed and the function returns false.
	
	// the final parameter void * dummy in Adaptive_Step, F, and norm is basically a hack to allow us to use functions that are non-static member
	// functions of classes
	
	bool Adaptive_Step(double*& y,bool& y_changed,int n,double& h,double tol,bool F(double*& yp,double *y,void * dummy),double norm(double *y,double *v,void *dummy),void * dummy);
}

#endif