#include "VelocityVerletIntegrator.h"

VelocityVerletIntegrator::VelocityVerletIntegrator(Sim &sim, double dt)
    : sim(sim), dt(dt), iter(0), time(0)
{
}

void VelocityVerletIntegrator::Step()
{
    sim.ComputeHalfStepVelocities(dt);
    sim.ComputeNewPositions(dt);
    sim.ComputeAccelerations();
    sim.ComputeNewVelocities(iter, dt);

    ++iter;
    time += dt;
    sim.EndStep(iter, time);
}

void VelocityVerletIntegrator::Run(int num_steps)
{
    for (int i = 0; i < num_steps; ++i) {
	Step();
    }
}
