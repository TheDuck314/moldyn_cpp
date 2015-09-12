#include "VelocityVerletIntegrator.h"
#include "LennardJonesConstantPressureSim.h"

#include <omp.h>

std::vector<double> GetRange(double min, double max, int count)
{
    std::vector<double> ret(count);
    if (count == 1) {
	ret[0] = count;
    } else {
	for (int i = 0; i < count; ++i) {
	    ret[i] = min + (max - min) * i / (count - 1);
	}
    }
    return ret;
}

int main(int argc, char** argv)
{
    //int N = 500;
    double initial_radius = 60;
    int thermostat_interval = 3000;
    double sim_time = 2000000;

    std::vector<double> temperatures = GetRange(0.47, 0.55, 16);
    std::vector<double> pressures{ 0.06 };
    std::vector<double> Ns{ 5000 };

    int num_sims = Ns.size() * temperatures.size();

    //omp_set_num_threads(9);

    int sim_counter = 0;

#pragma omp parallel for
    for (int i = 0; i < num_sims; ++i) {
	int sim_index;
#pragma omp critical
	{
	    sim_index = sim_counter;
	    sim_counter++;
	}
	int indT = sim_index % temperatures.size();
	sim_index /= temperatures.size();
	int indP = sim_index % pressures.size();
	sim_index /= pressures.size();
	int indN = sim_index % Ns.size();

	double T = temperatures[indT];
	double P = pressures[indP];
	int N = Ns[indN];

	char filename[512];
	sprintf(filename, "transition_P0.06/run2/sim_N%d_T%0.4f_P%0.6f.dat", N, T, P);
	LennardJonesConstantPressureSim sim(N, T, P, initial_radius, thermostat_interval, filename);

	double dt = 0.005;
	int num_steps = (int)(sim_time / dt);
	VelocityVerletIntegrator integrator(sim, dt);
	printf("starting %d-step sime of %s\n", num_steps, filename);
	integrator.Run(num_steps);
    }
}