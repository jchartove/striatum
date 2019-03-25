// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// generic conductance template
// you can either fill this out yourself
// or use the conductance class within MATLAB
// to auto-generate C++ code 
#ifndef GOLOMBNA
#define GOLOMBNA
#include "conductance.hpp"
//inherit conductance class spec
class golombNa: public conductance {
public:
// specify parameters + initial conditions
golombNa(double g_, double E_, double m_, double h_)
{
gbar = g_;
E = E_;
m = m_;
h = h_;
// defaults 
if (isnan(gbar)) { gbar = 112.5; }
if (isnan (m)) { m = 0; }
if (isnan (h)) { h = 1; }
if (isnan (E)) { E = 50; }
p = 3;
q = 1;
}
double m_inf(double V, double Ca);
double h_inf(double V, double Ca);
//double tau_m(double V, double Ca);
double tau_h(double V, double Ca);
string getClass(void);
};
string golombNa::getClass(){return "golombNa";}
double golombNa::m_inf(double V, double Ca) {return  1/(1+exp(-(V+24)/11.5));}
double golombNa::h_inf(double V, double Ca) {return  1/(1+exp(-(V+58.3)/-6.7));}
//double golombNa::tau_m(double V, double Ca) {return  0;}
double golombNa::tau_h(double V, double Ca) {return  0.5+14/(1+exp(-(V+60)/-12));}
#endif
