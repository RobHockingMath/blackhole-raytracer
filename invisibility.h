#ifndef INVISIBILITY_H
#define INVISIBILITY_H


#include <math.h>



void expand(double[],double,double[]);
void contract(double[],double,double[]);
void expand_derivative(double[],double[],double,double[]);
void contract_derivative(double[],double[],double,double[]);
void initial_value(double[],double[],double,double[],double[]);
void boundary_value(double[],double[],double,double[],double[]);
void get_direction(double[],double[],double,double[]);
void transport_tangent(double[],double,double,double[]);
void flow(double[],double[],double,double);

#endif
