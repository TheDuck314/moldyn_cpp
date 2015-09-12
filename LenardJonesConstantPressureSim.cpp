#include "LennardJonesConstantPressureSim.h"
#include "RNG.h"
#include <cassert>
#include <sys/time.h>
#include <time.h>

LennardJonesConstantPressureSim::LennardJonesConstantPressureSim(int N, double target_temp, double pressure,
    double initial_radius, int thermostat_interval,
    const char* data_filename)
    : pos(N), vel(N), half_vel(N), acc(N),
    N(N), pressure(pressure), target_temp(target_temp),
    thermostat_interval(thermostat_interval),
    data_filename(data_filename)
{
    wall_radius = initial_radius;
    wall_mass = 50;
    wall_stiffness = 5;

    // place atoms
    for (int i = 0; i < N; ++i) {
	bool collision;
	int tries = 0;
	do {
	    pos[i].x = 2 * (-0.5 + RNG::UnitInterval()) * wall_radius;
	    pos[i].y = 2 * (-0.5 + RNG::UnitInterval()) * wall_radius;
	    collision = false;
	    for (int j = 0; j < i; ++j) {
		if (dist(pos[i], pos[j]) < 1) {
		    collision = true;
		    break;
		}
	    }
	    tries += 1;
	    assert(tries < 10000);
	} while (collision);
    }

    SampleVelocitiesFromMaxwellBoltzmann(target_temp);
    
    r2_cutoff = 9;
    r_cutoff = std::sqrt(r2_cutoff);

    FILE *f = fopen(data_filename, "w");
    if (!f) {
	printf("couldn't open data file %s\n", data_filename);
	exit(-1);
    }
    fclose(f);
}

LennardJonesConstantPressureSim::~LennardJonesConstantPressureSim()
{
}

void LennardJonesConstantPressureSim::ComputeHalfStepVelocities(double dt)
{
    for (int i = 0; i < N; ++i) {
	half_vel[i] = vel[i] + (0.5 * dt) * acc[i];
    }
    wall_half_vel = wall_vel + 0.5 * dt * wall_acc;
}

void LennardJonesConstantPressureSim::ComputeNewPositions(double dt)
{
    for (int i = 0; i < N; ++i) {
	pos[i] += dt * half_vel[i];
    }
    wall_radius += dt * wall_half_vel;
}

void LennardJonesConstantPressureSim::ComputeAccelerations()
{
    // interatomic forces
    for (int i = 0; i < N; ++i) {
	acc[i].Zero();
    }

    /*total_interatom_PE = 0;
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < i; ++j) {
	    Vector2 r = pos[i] - pos[j];
	    double r2 = norm2(r);
	    if (r2 <= r2_cutoff) {
		//if (i == 0 || j == 0) printf("old: %d <-> %d\n", i, j);
		double rm6 = 1 / (r2*r2*r2);
		double FoverR = 12 * rm6 * (rm6 - 1) / r2;
		total_interatom_PE += rm6 * (rm6 - 2);
		Vector2 F = FoverR * r;
		acc[i] += F;
		acc[j] -= F;
	    }
	}
    }
    std::vector<Vector2> oldacc(N);
    for (int i = 0; i < N; ++i) {
	oldacc[i] = acc[i];
	acc[i].Zero();
    }*/


//    struct timeval start_time;
//    gettimeofday(&start_time, NULL);    
//    static int best_Ncell = 1;
//    static long best_elapsed_us = 999999999;
//    static int test_counter = 0;
//    static int test_interval = 1000;

    // decide on cell size
    int max_Ncell = (int)(2 * wall_radius / r_cutoff) - 1;
    if (max_Ncell > std::sqrt(N)) max_Ncell = (int)std::sqrt(N);
    //printf("max_Ncell = %d, best_Ncell = %d\n", max_Ncell, best_Ncell);
    int Ncell = max_Ncell;
//    if (test_counter  < max_Ncell) {
//	Ncell = 1 + test_counter;
//    }
//    test_counter = (test_counter + 1) % test_interval;
    double cell_size = 2 * wall_radius / Ncell;
    //printf("Ncell = %d, wall_radius = %f, cell_size = %f, r_cutoff = %f\n", Ncell, wall_radius, cell_size, r_cutoff);
    assert(cell_size > r_cutoff);

    std::vector<std::vector<std::vector<int> > > cell_members(Ncell,
	std::vector<std::vector<int> >(Ncell));

    // put atoms in cells
    for (int i = 0; i < N; ++i) {
	int cellX = (int)(wall_radius + pos[i].x) / cell_size;
	int cellY = (int)(wall_radius + pos[i].y) / cell_size;
	if (cellX < 0) cellX = 0;
	if (cellY < 0) cellY = 0;
	if (cellX >= Ncell) cellX = Ncell - 1;
	if (cellY >= Ncell) cellY = Ncell - 1;
	cell_members[cellX][cellY].push_back(i);
    }
    // compute forces between atoms
    int neighborOffsetsX[] = { 0, 1, 1, 0, -1 };
    int neighborOffsetsY[] = { 0, 0, 1, 1, 1 };
    total_interatom_PE = 0;
    for (int cellX = 0; cellX < Ncell; ++cellX) {
	for (int cellY = 0; cellY < Ncell; ++cellY) {
	    std::vector<int> &this_cell_members = cell_members[cellX][cellY];
	    for (int off = 0; off < 5; ++off) {
		int neighbor_cellX = (cellX + neighborOffsetsX[off]);
		int neighbor_cellY = (cellY + neighborOffsetsY[off]);
		if (neighbor_cellX < 0 || neighbor_cellY < 0 || neighbor_cellX >= Ncell || neighbor_cellY >= Ncell) continue;
		std::vector<int> &neighbor_cell_members = cell_members[neighbor_cellX][neighbor_cellY];
		for (unsigned this_atom = 0; this_atom < this_cell_members.size(); ++this_atom) {
		    for (unsigned neighbor_atom = 0; neighbor_atom < neighbor_cell_members.size(); ++neighbor_atom) {
			int i = this_cell_members[this_atom];
			int j = neighbor_cell_members[neighbor_atom];
			if (off == 0 && neighbor_atom >= this_atom) break; // avoid self-force and double-counting forces within a cell

			Vector2 r = pos[i] - pos[j];
			double r2 = norm2(r);
			if (r2 <= r2_cutoff) {
			    //if (i == 0 || j == 0) printf("new: %d <-> %d\n", i, j);
			    double rm6 = 1 / (r2*r2*r2);
			    double FoverR = 12 * rm6 * (rm6 - 1) / r2;
			    total_interatom_PE += rm6 * (rm6 - 2);
			    Vector2 F = FoverR * r;
			    acc[i] += F;
			    acc[j] -= F;
			}
		    }
		}
	    }
	}
    }
    
//    struct timeval end_time;
//    gettimeofday(&end_time, NULL);
//    long elapsed_us = (end_time.tv_sec - start_time.tv_sec) * 1000000 + (end_time.tv_usec - start_time.tv_usec);
    //printf("elapsed_us = %ld\n", elapsed_us);

//    if (Ncell == best_Ncell) {
//	best_elapsed_us = elapsed_us;
//    } else if (best_Ncell > max_Ncell || elapsed_us < best_elapsed_us) {
//	best_Ncell = Ncell;
//	best_elapsed_us = elapsed_us;
//    }
//    if (test_counter == 100) printf("test_counter = %d, Ncell = %d, max_Ncell = %d, elapsed_us = %ld, volume = %f\n", test_counter, Ncell, max_Ncell, elapsed_us, wall_radius * wall_radius * 4);

    /*for (int i = 0; i < N; ++i) {
	//printf("acc[%d] = (%f,%f), oldacc[%d] = (%f,%f)\n", i, acc[i].x, acc[i].y, i, oldacc[i].x, oldacc[i].y);
	assert(norm2(acc[i] - oldacc[i]) < 1e-12);
    }*/

    // force of pressure on wall
    double wall_area = 8 * wall_radius;
    double wall_force = -pressure * wall_area;
    wall_acc = wall_force / wall_mass;
    volume = 4 * wall_radius * wall_radius;
    pressure_PE = pressure * volume;

    // force between wall and atoms
    total_wall_PE = 0;
    for (int i = 0; i < N; ++i) {
	if (pos[i].x < -wall_radius) {
	    double d = (pos[i].x + wall_radius);
	    double F = -wall_stiffness * d;
	    acc[i].x += F;
	    wall_acc += F / wall_mass;
	    total_wall_PE += 0.5 * wall_stiffness * d * d;
	}
	if (pos[i].x > wall_radius) {
	    double d = (pos[i].x - wall_radius);
	    double F = wall_stiffness * d;
	    acc[i].x -= F;
	    wall_acc += F / wall_mass;
	    total_wall_PE += 0.5 * wall_stiffness * d * d;
	}
	if (pos[i].y < -wall_radius) {
	    double d = (pos[i].y + wall_radius);
	    double F = -wall_stiffness * d;
	    acc[i].y += F;
	    wall_acc += F / wall_mass;
	    total_wall_PE += 0.5 * wall_stiffness * d * d;
	}
	if (pos[i].y > wall_radius) {
	    double d = (pos[i].y - wall_radius);
	    double F = wall_stiffness * d;
	    acc[i].y -= F;
	    wall_acc += F / wall_mass;
	    total_wall_PE += 0.5 * wall_stiffness * d * d;
	}
    }
}

void LennardJonesConstantPressureSim::ComputeNewVelocities(int iter, double dt)
{
    if ((iter % thermostat_interval) == 0) {
	SampleVelocitiesFromMaxwellBoltzmann(target_temp);
    } else {
	for (int i = 0; i < N; ++i) {
	    vel[i] = half_vel[i] + (0.5 * dt) * acc[i];
	}
	wall_vel = wall_half_vel + 0.5 * dt * wall_acc;
    }

    total_atom_KE = 0;
    for (int i = 0; i < N; ++i) {
	total_atom_KE += 0.5 * norm2(vel[i]);
    }

    wall_KE = 0.5 * wall_mass * wall_vel * wall_vel;
}

void LennardJonesConstantPressureSim::EndStep(int iter, double time)
{
    // write data
    if (iter % 2000 == 0) {
	double total_E = total_interatom_PE + total_atom_KE + total_wall_PE + wall_KE + pressure_PE;
	double temp = MeasureTemperature();
	
	
	FILE *data_file = fopen(data_filename, "a");
	fprintf(data_file, "%d\t%d\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\n",
	    N, iter, time, volume, temp, total_E, total_interatom_PE, total_atom_KE, total_wall_PE, wall_KE, pressure_PE);
	fclose(data_file);
    }
}


void LennardJonesConstantPressureSim::SampleVelocitiesFromMaxwellBoltzmann(double T)
{
    double std = std::sqrt(T);
    for (int i = 0; i < N; ++i) {
	vel[i].x = RNG::Gaussian(std);
	vel[i].y = RNG::Gaussian(std);
    }

    double wall_std = std::sqrt(T / wall_mass);
    wall_vel = RNG::Gaussian(wall_std);
}

double LennardJonesConstantPressureSim::MeasureTemperature()
{
    // average kinetic energy at temperature T is equal to T (in two dimensions)
    double mean_KE = 0;
    for (int i = 0; i < N; ++i) {
	mean_KE += 0.5 * norm2(vel[i]);
    }
    mean_KE /= N;
    return mean_KE;
}