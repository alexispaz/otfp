#pragma once
#ifndef OTFP_ENGINE_H
#define OTFP_ENGINE_H

typedef struct otfp_engine otfp_engine;

otfp_engine * otfp_engine_create ( int ncenters, int ncvs, int nrestr, long seed );
void otfp_engine_free ( otfp_engine * engine );

int otfp_engine_add_center ( otfp_engine * engine, int natoms, const int * indices, const double * masses );
int otfp_engine_set_center_position ( otfp_engine * engine, int center, const double xyz[3] );
const double * otfp_engine_center_force ( otfp_engine * engine, int center );

int otfp_engine_add_cv ( otfp_engine * engine, const char * type, int ncenters, const int * centers,
    double zmin, double zmax, const char * boundary, double boundary_k );
int otfp_engine_add_restraint ( otfp_engine * engine, double k, double target,
    int ncvs, const double * coefficients, const char * function,
    double zmin, double zmax, const char * boundary, double boundary_k );
int otfp_engine_add_tamd ( otfp_engine * engine, int restraint,
    double gamma, double kt, double dt, int chapeau_id, int chapeau_dim );
int otfp_engine_add_smd ( otfp_engine * engine, int restraint,
    double target, int first_step, int last_step );
int otfp_engine_add_adam ( otfp_engine * engine, int restraint,
    double step_size, double beta1, double beta2, double epsilon,
    double decay, double kt, double min_step );

int otfp_engine_compute ( otfp_engine * engine, int first_step, int step, double external_bias );
double otfp_engine_bias_energy ( otfp_engine * engine );

#endif
