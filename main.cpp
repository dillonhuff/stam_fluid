#include "utils.h"
#include "visualize.h"

void add_source(const int N, float* x, float* s, const float dt) {
  int size = (N+2)*(N+2);

  for (int i = 0; i < size; i++) {
    x[i] += dt*s[i];
  }

}

// Conditions
// 1 = Setting x bounds, symmetric
// 2 = Setting y bounds, symmetric
// 0 = Continuous
void set_bnd(const int N, const int b, float* x) {

  // Set edges
  for (int i = 1; i <= N; i++) {
    x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];

    x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
  }

  // Set corners
  x[IX(0, 0)] = 0.5*( x[IX(1, 0)] + x[IX(0, 1)] );
  x[IX(0, N + 1)] = 0.5*( x[IX(0, N)] + x[IX(1, N + 1)] );
  x[IX(N + 1, 0)] = 0.5*( x[IX(N, 0)] + x[IX(N + 1, 1)] );
  x[IX(N + 1, N + 1)] = 0.5*( x[IX(N + 1, N)] + x[IX(N, N + 1)] );

}

void diffuse(const int N, const int b, float* x, float* x0, float diff, float dt) {
  float a = dt*diff*N*N;

  int num_iterations = 20;
  for (int k = 0; k < num_iterations; k++) {

    for (int i = 1; i <= N; i++) {
      for (int j = 1; j <= N; j++) {
	x[IX(i, j)] = (x[IX(i, j)] +
		       a*(x[IX(i - 1, j)] +
			  x[IX(i + 1, j)] +
			  x[IX(i, j - 1)] +
			  x[IX(i, j + 1)])) / (1 + 4*a);
      }
    }

    set_bnd(N, b, x);
  }
}

void advect(const int N, const int b,
	    float* d,
	    float* d0,
	    float* u, // X velocities?
	    float* v, // Y velocities?
	    const float dt) {
  float dt0 = dt*N;

  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      float x = i - dt0*u[IX(i, j)];
      float y = j - dt0*v[IX(i, j)];

      if (x < 0.5) {
	x = 0.5;
      }
      if (x > N + 0.5) {
	x = N + 0.5;
      }

      int i0 = (int) x;
      int i1 = i0 + 1;

      if (y < 0.5) {
	y = 0.5;
      }
      if (y > N + 0.5) {
	y = N + 0.5;
      }

      int j0 = (int) y;
      int j1 = i0 + 1;

      int s1 = x - i0;
      int s0 = 1 - s1;
      int t1 = y - j0;
      int t0 = 1 - t1;

      d[IX(i, j)] = s0*(t0*d0[IX(i0, j0)] + t1*d0[IX(i0, j1)]) +
	s1*(t0*d0[IX(i1, j0)] + t1*d0[IX(i1, j1)]);
    }
  }

  set_bnd(N, b, d);
}

void dens_step(const int N,
	       float* x, float* x0, float* u, float* v,
	       const float diff,
	       const float dt) {
  add_source(N, x, x0, dt);

  SWAP(x0, x);
  diffuse(N, 0, x, x0, diff, dt);
  SWAP(x0, x);
  advect(N, 0, x, x0, u, v, dt);
}

void project(const int N,
	     float* u,
	     float* v,
	     float* p,
	     float* div) {

  float h = 1.0 / N;

  // Set the divergence
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      div[IX(i, j)] = -0.5*h*(u[IX(i + 1, j)] - u[IX(i - 1, j)] +
			      v[IX(i, j+1)] - v[IX(i, j - 1)]);
      p[IX(i, j)] = 0;
    }
  }

  set_bnd(N, 0, div);
  set_bnd(N, 0, p);

  
  // Solve for pressure
  for (int k = 0; k < 20; k++) {

    for (int i = 1; i <= N; i++) {
      for (int j = 1; j <= N; j++) {
	p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
		       p[IX(i, j + 1)] + p[IX(i, j - 1)]) / 4.0;
      }
    }

    set_bnd(N, 0, p);
  }

  // Update original quantities
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      u[IX(i, j)] -= 0.5*(p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
      v[IX(i, j)] -= 0.5*(p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
    }
  }

  set_bnd(N, 1, u);
  set_bnd(N, 2, v);

}

void vel_step(const int N,
	      float* u, float* v, float* u0, float* v0,
	      const float visc, const float dt) {
  add_source(N, u, u0, dt);
  add_source(N, v, v0, dt);

  SWAP(u0, u);
  diffuse(N, 1, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(N, 2, v, v0, visc, dt);

  project(N, u, v, u0, v0);

  SWAP(u0, u);
  SWAP(v0, v);

  // Advection of velocity?
  advect(N, 1, u, u0, u0, v0, dt);
  advect(N, 2, v, v0, v0, v0, dt);

  project(N, u, v, u0, v0);
  
}

int main() {
  int N = 20;
  int size = (N+2)*(N+2);
  //int u[size], v[size], u_prev[size], v_prev[size];
  float dens[size]; //, dens_prev[size];

  visualize_density(N, dens);

  return 0;
}
