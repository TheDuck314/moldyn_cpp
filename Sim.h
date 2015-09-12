#ifndef SIM_H
#define SIM_H

class Sim
{
public:
    virtual void ComputeHalfStepVelocities(double dt) = 0;
    virtual void ComputeNewPositions(double dt) = 0;
    virtual void ComputeAccelerations() = 0;
    virtual void ComputeNewVelocities(int iter, double dt) = 0;
    virtual void EndStep(int iter, double time) = 0;
};

#endif