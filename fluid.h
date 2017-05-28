#pragma once

void add_source(const int N, float* x, float* s, const float dt);

void set_bnd(const int N, const int b, float* x);

void vel_step(const int N,
	      float* u, float* v, float* u0, float* v0,
	      const float visc, const float dt);

void dens_step(const int N,
	       float* x, float* x0, float* u, float* v,
	       const float diff,
	       const float dt);

