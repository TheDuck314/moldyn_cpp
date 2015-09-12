#ifndef VELOCITY_VERLET_INTEGRATOR_H
#define VELOCITY_VERLET_INTEGRATOR_H

#include "Sim.h"

class VelocityVerletIntegrator
{
public:
    VelocityVerletIntegrator(Sim &sim, double dt);
    void Step();
    void Run(int num_steps);

private:
    Sim &sim;
    double dt;
    int iter;
    double time;
};

#endif