#pragma once
#include <mpi.h>
#include <array>

/**
 * @brief Data structure that holds input parameters for 
 * transfer between MPI processes
 *
 */
struct Parameters {
    double xlength; /* length of the domain x-dir. */
    double ylength; /* length of the domain y-dir. */
    double nu;      /* viscosity */
    double dt;      /* time step */
    double omg;     /* relaxation factor */
    double eps;     /* accuracy bound for pressure */
    double tau;     /* safety factor for time step */
    double gamma;   /* uppwind differencing factor */
    double UI;      /* velocity x-direction */
    double VI;      /* velocity y-direction */
    double PI;      /* pressure */
    double TI;      /* initial temperature */
    double GX;      /* gravitation x-direction */
    double GY;      /* gravitation y-direction */
    double UIN;     /* inlet velocity x-direction */
    double VIN;     /* inlet velocity y-direction */
    double TIN;     /* inlet temperature */
    double beta;    /* thermal expansion coefficient */
    double alpha;   /* thermal diffusivity */
    double t_end;   /* Simulation time */
    double output_freq; /* Solution file outputting frequency */
    int itermax;    /* max. number of iterations for pressure per time step */
    int imax;       /* number of cells x-direction */
    int jmax;       /* number of cells y-direction */
    int iproc;      /* number of subdomains in x direction */
    int jproc;      /* number of subdomains in y direction */
    int num_walls;  /* number of walls */
    bool energy_eq; /* if energy equation is turned on */

    // Currently support 4 walls (Change if needed) change blocklengths of Case.cpp
    std::array<int, 4> wall_idx = {-1, -1, -1, -1};         /* Wall indices (Supports upto 4 walls) */
    std::array<double, 4> wall_vel = {0.0, 0.0, 0.0, 0.0};  /* Wall velocity against the wall index */
    std::array<double, 4> wall_temp = {0.0, 0.0, 0.0, 0.0}; /* Wall temperature against the wall index */
};


