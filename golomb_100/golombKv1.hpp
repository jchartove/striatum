// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// generic conductance template
// you can either fill this out yourself
// or use the conductance class within MATLAB
// to auto-generate C++ code 
#ifndef GOLOMBKV1
#define GOLOMBKV1
#include "conductance.hpp"
//inherit conductance class spec
class golombKv1: public conductance {
public:
// specify parameters + initial conditions
golombKv1(double g_, double E_, double m_, double h_)
{
gbar = g_;
E = E_;
m = m_;
h = h_;
// defaults 
if (isnan(gbar)) { gbar = 3; }
if (isnan (m)) { m = 0; }
if (isnan (h)) { h = 1; }
if (isnan (E)) { E = -90; }
p = 3;
q = 1;
}
double m_inf(double V, double Ca);
double h_inf(double V, double Ca);
double tau_m(double V, double Ca);
double tau_h(double V, double Ca);
string getClass(void);
};
string golombKv1::getClass(){return "golombKv1";}
double golombKv1::m_inf(double V, double Ca) {return  1/(1+exp(-(V+50)/20));}
double golombKv1::h_inf(double V, double Ca) {return  1/(1+exp(-(V+65)/-6));}
double golombKv1::tau_m(double V, double Ca) {return  2;}
double golombKv1::tau_h(double V, double Ca) {return  150;}
#endif
