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

std::vector<double> GetLogRange(double min, double max, int count)
{
    std::vector<double> ret(count);
    double ratio = std::pow(max/min, 1.0/(count-1));
    if (count == 1) {
	ret[0] = count;
    } else {
	for (int i = 0; i < count; ++i) {
	    ret[i] = min * std::pow(ratio, i);
	}
    }
    return ret;
}

void RunLiquidToSolidSim(int N, double P, double T_start, double T_end, double total_time, int num_temps)
{
  int thermostat_interval = 3000;
  char filename[512];
  sprintf(filename, "liquid_to_solid/run7/sim_N%d_Tstart%0.4f_Tend%0.4f_time%d_P%0.6f.dat", N, T_start, T_end, (int)total_time, P);
  LennardJonesConstantPressureSim sim(N, T_start, P, -1, thermostat_interval, filename);
  double dt = 0.005;
  double time_per_temp = total_time / num_temps;
  int steps_per_temp = (int)(time_per_temp / dt);
  double T_decrement = (T_start - T_end) / (num_temps - 1);

  double T = T_start;

  printf("starting run N=%d, P=%f, %d temperatures starting at %f and decreasing by %f\n", N, P, num_temps, T_start, T_decrement);

  VelocityVerletIntegrator integrator(sim, dt);
  for (int i = 0; i < num_temps; ++i) {
    sim.SetParams(T, P);
    integrator.Run(steps_per_temp);
    T -= T_decrement;
  }
}

void RunPressureRangeSim(int N, double T, double P_start, double P_end, double total_time, int num_P)
{
  int thermostat_interval = 3000;
  char filename[512];
  sprintf(filename, "vaporization/run2/sim_N%d_T_%0.4f_Pstart%0.6f_Tend%0.6f_time%d.dat", N, T, P_start, P_end, (int)total_time);
  LennardJonesConstantPressureSim sim(N, T, P_start, -1, thermostat_interval, filename);
  double dt = 0.01; //0.005;
  std::vector<double> Ps = GetLogRange(P_start, P_end, num_P);
  double time_per_P = total_time / num_P;
  int steps_per_P = (int)(time_per_P / dt);

  printf("starting run N=%d, T=%f, %d pressures starting at %f and ending at %f\n", N, T, num_P, P_start, P_end);

  VelocityVerletIntegrator integrator(sim, dt);
  for (int i = 0; i < num_P; ++i) {
    sim.SetParams(T, Ps[i]);
    integrator.Run(steps_per_P);
  }
}

int main(int argc, char** argv)
{
  /*
    //int N = 500;
    double initial_radius = 50;
    int thermostat_interval = 3000;
    double sim_time = 50000; 

    std::vector<double> temperatures = GetRange(0.45, 0.65, 9);
    std::vector<double> pressures{ 0.06 };
    std::vector<double> Ns{ 2000 };

    int num_sims = Ns.size() * temperatures.size();

    omp_set_num_threads(9);

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
	sprintf(filename, "transition_P0.06/run3/sim_N%d_T%0.4f_P%0.6f.dat", N, T, P);
	LennardJonesConstantPressureSim sim(N, T, P, initial_radius, thermostat_interval, filename);

	double dt = 0.005;
	int num_steps = (int)(sim_time / dt);
	VelocityVerletIntegrator integrator(sim, dt);
	printf("starting %d-step sime of %s\n", num_steps, filename);
	integrator.Run(num_steps);
    }
    */

  
  /*
    std::vector<int> Ns { 100, 400, 1600 };
    std::vector<double> total_times { 10000, 40000 };
    //std::vector<double> Ps { 0.001, 0.003, 0.07, 0.01, 0.02, 0.03, 0.04 };
    std::vector<double> Ps { 0.02 };
    double T_start = 0.33;
    double T_end = 0.43;
    double num_temps = 50;

    int sim_counter = 0;
    int num_sims = 2 * Ns.size() * total_times.size() * Ps.size();
#pragma omp parallel for
    for (int i = 0; i < num_sims; ++i) {
	int sim_index;
#pragma omp critical
	{
	    sim_index = sim_counter;
	    sim_counter++;
	}
        int indN = sim_index % Ns.size();
        sim_index /= Ns.size();
        int indTT = sim_index % total_times.size();
        sim_index /= total_times.size();
        int indP = sim_index / 2;
        bool reverse = (sim_index % 2) == 1;
        int N = Ns[indN];
        double total_time = total_times[indTT];
        double P = Ps[indP];

        if (!reverse) RunLiquidToSolidSim(N, P, T_start, T_end, total_time, num_temps);
        else          RunLiquidToSolidSim(N, P, T_end, T_start, total_time, num_temps);
    }*/

    std::vector<int> Ns { 1000 };
    std::vector<double> total_times { 500000 };
    std::vector<double> Ts { 0.45, 0.5, 0.55, 0.6, 0.65, 0.70 };
    double num_P = 50;

    int sim_counter = 0;
    int num_sims = 2 * Ns.size() * total_times.size() * Ts.size();
#pragma omp parallel for
    for (int i = 0; i < num_sims; ++i) {
	int sim_index;
#pragma omp critical
	{
	    sim_index = sim_counter;
	    sim_counter++;
	}
        int indN = sim_index % Ns.size();
        sim_index /= Ns.size();
        int indTT = sim_index % total_times.size();
        sim_index /= total_times.size();
        int indT = sim_index / 2;
        bool reverse = (sim_index % 2) == 1;
        int N = Ns[indN];
        double total_time = total_times[indTT];
        double T = Ts[indT];

        double estPvap = 8*std::exp(-2.7/T);
        printf("T = %f, estPvap = %f\n", T, estPvap);
        double Plow = estPvap / 3;
        double Phigh = estPvap * 3;

        double P_start = (reverse ? Phigh : Plow);
        double P_end = (reverse ? Plow : Phigh);

        RunPressureRangeSim(N, T, P_start, P_end, total_time, num_P);
    }
}

