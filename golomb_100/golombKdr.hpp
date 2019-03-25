// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// generic conductance template
// you can either fill this out yourself
// or use the conductance class within MATLAB
// to auto-generate C++ code 
#ifndef GOLOMBKDR
#define GOLOMBKDR
#include "conductance.hpp"
//inherit conductance class spec
class golombKdr: public conductance {
public:
// specify parameters + initial conditions
golombKdr(double g_, double E_, double m_)
{
gbar = g_;
E = E_;
m = m_;
//h = h_;
// defaults 
if (isnan(gbar)) { gbar = 225; }
if (isnan (m)) { m = 0; }
//if (isnan (h)) { h = 1; }
if (isnan (E)) { E = -90; }
p = 2;
//q = 0;
}
double m_inf(double V, double Ca);
// double h_inf(double V, double Ca);
double tau_m(double V, double Ca);
// double tau_h(double V, double Ca);
string getClass(void);
};
string golombKdr::getClass(){return "golombKdr";}
double golombKdr::m_inf(double V, double Ca) {return  1/(1+exp(-(V+12.4)/6.8));}
// double golombKdr::h_inf(double V, double Ca) {return  0;}
double golombKdr::tau_m(double V, double Ca) {return  (0.087+11.4/(1+exp((V+14.6)/8.6)))*(0.087+11.4/(1+exp(-(V-1.3)/18.7)));}
// double golombKdr::tau_h(double V, double Ca) {return  0;}
#endif
