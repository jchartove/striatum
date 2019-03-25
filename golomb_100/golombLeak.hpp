// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// Slow Calcium conductance
#ifndef GOLOMBLEAK
#define GOLOMBLEAK
#include "conductance.hpp"

//inherit conductance class spec
class golombLeak: public conductance {

public:

    // specify parameters + initial conditions
    golombLeak(double g_, double E_)
    {
        gbar = g_;
        g = gbar; // this is important as integrate doesn't do anything in the leak channels
        E = E_;

        if (isnan (E)) { E = -70; }

    }


    void integrate(double, double);
    void integrateMS(int, double, double);

    string getClass(void);

};

string golombLeak::getClass(){return "golombLeak";}

void golombLeak::integrate(double V, double Ca) {
    // do nothing
    // because there is nothing to integrate
}

// Runge-Kutta 4 integrator
void golombLeak::integrateMS(int k, double V, double Ca)
{

    // do nothing
    // because there is nothing to integrate
}



#endif
