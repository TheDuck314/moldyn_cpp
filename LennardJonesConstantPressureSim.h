#ifndef LENNARD_JONES_CONSTANT_PRESSURE_SIM
#define LENNARD_JONES_CONSTANT_PRESSURE_SIM

#include "Sim.h"
#include "Vector2.h"
#include <vector>
#include <cstdio>

class LennardJonesConstantPressureSim : public Sim
{
public:
    LennardJonesConstantPressureSim(int N, double target_temp, double pressure,
	double initial_radius, int thermostat_interval, 
	const char* data_filename);
    ~LennardJonesConstantPressureSim();

    void ComputeHalfStepVelocities(double dt);
    void ComputeNewPositions(double dt);
    void ComputeAccelerations();
    void ComputeNewVelocities(int iter, double dt);
    void EndStep(int iter, double time);    

private:
    void SampleVelocitiesFromMaxwellBoltzmann(double T);
    double MeasureTemperature();

    // state
    std::vector<Vector2> pos;
    std::vector<Vector2> vel;
    std::vector<Vector2> half_vel;
    std::vector<Vector2> acc;
    double wall_radius;
    double wall_vel;
    double wall_half_vel;
    double wall_acc;

    // parameters
    int N;
    double pressure;
    double target_temp;
    double wall_mass;
    double wall_stiffness;
    double r2_cutoff;
    double r_cutoff;
    int thermostat_interval;

    // statistics
    double volume;
    double total_interatom_PE;
    double total_atom_KE;
    double total_wall_PE;
    double wall_KE;
    double pressure_PE;

    const char *data_filename;
};

#endif