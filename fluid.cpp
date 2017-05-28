#include "fluid.h"

#include "utils.h"
#include "visualize.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

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

  check_nans(N, x, "in bound b = " + std::to_string(b));

  // Set corners
  x[IX(0, 0)] = 0.5*( x[IX(1, 0)] + x[IX(0, 1)] );
  x[IX(0, N + 1)] = 0.5*( x[IX(0, N)] + x[IX(1, N + 1)] );
  x[IX(N + 1, 0)] = 0.5*( x[IX(N, 0)] + x[IX(N + 1, 1)] );
  x[IX(N + 1, N + 1)] = 0.5*( x[IX(N + 1, N)] + x[IX(N, N + 1)] );

  check_nans(N, x, "in bound b after corner = " + std::to_string(b));

}

void diffuse(const int N, const int b, float* x, float* x0, float diff, float dt) {
  float a = dt*diff*N*N;

  int num_iterations = 20;
  for (int k = 0; k < num_iterations; k++) {

    for (int i = 1; i <= N; i++) {
      for (int j = 1; j <= N; j++) {
	x[IX(i, j)] = (x0[IX(i, j)] +
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

  check_nans(N, d0, "d0");

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
      int j1 = j0 + 1;

      int s1 = x - i0;
      int s0 = 1 - s1;
      int t1 = y - j0;
      int t0 = 1 - t1;

      if (i0 < 0) {
	cout << "X = " << x << endl;
	cout << "dt0 = " << dt0 << endl;
	cout << "u[IX(i, j)] = " << u[IX(i, j)] << endl;
	cout << "dt0*u[IX(i, j)] = " << dt0*u[IX(i, j)] << endl;
	cout << "i = " << i << endl;
	cout << "j = " << j << endl;
      }

      assert(i0 >= 0);
      assert(j0 >= 0);
      

      // cout << "i0 = " << i0 << endl;
      // cout << "j0 = " << j0 << endl;

      double s0t0 = s0*t0*d0[IX(i0, j0)];
      double s0t1 = s0*t1*d0[IX(i0, j1)];
      double s1t0 = s1*t0*d0[IX(i1, j0)];
      double s1t1 = s1*t1*d0[IX(i1, j1)];
      d[IX(i, j)] = s0t0 + s0t1 + s1t0 + s1t1;

      if (isnan(d[IX(i, j)])) {
	cout << "s0 = " << s0 << endl;
	cout << "t0 = " << t0 << endl;
	cout << "s1 = " << s1 << endl;
	cout << "t1 = " << t1 << endl;
	cout << "d0[i0, j0] = " << d0[IX(i0, j0)] << endl;
	cout << "( " << i << ", " << j << " )" << endl;
	cout << "s0t0 = " << s0t0 << endl;
	cout << "s0t1 = " << s0t1 << endl;
	cout << "s1t0 = " << s1t0 << endl;
	cout << "s1t1 = " << s1t1 << endl;
	cout << "d[i, j] = " << d[IX(i, j)] << endl;
	assert(!isnan(d[IX(i, j)]));
      }
      
      // d[IX(i, j)] = s0*(t0*d0[IX(i0, j0)] + t1*d0[IX(i0, j1)]) +
      // 	s1*(t0*d0[IX(i1, j0)] + t1*d0[IX(i1, j1)]);
    }
  }

  cout << "Set entries" << endl;
  set_bnd(N, b, d);

  cout << "Done setting bounds" << endl;
}

void dens_step(const int N,
	       float* x, float* x0, float* u, float* v,
	       const float diff,
	       const float dt) {
  add_source(N, x, x0, dt);

  SWAP(x0, x);

  cout << "starting diffuse" << endl;
  diffuse(N, 0, x, x0, diff, dt);
  SWAP(x0, x);

  cout << "starting advect" << endl;
  advect(N, 0, x, x0, u, v, dt);

  cout << "Done with advect" << endl;
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


  cout << "max div = " << array_max(N, div) << endl;
  cout << "min div = " << array_min(N, div) << endl;

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

  check_nans(N, p, "p in project before u, v set");
  check_nans(N, div, "div in project before u, v, set");

  cout << "h = " << h << endl;
  
  // Update original quantities
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      double old_v = v[IX(i, j)];
      u[IX(i, j)] -= 0.5*(p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
      v[IX(i, j)] -= 0.5*(p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;

      if (isinf(v[IX(i, j)])) {
	cout << "old v = " << old_v << endl;
	cout << "subtracted = " << 0.5*(p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h << endl;
	assert(!isinf(v[IX(i, j)]));
	assert(!isinf(u[IX(i, j)]));
      }
    }
  }

  check_nans(N, v, "v in project before set_bnd");
  check_nans(N, u, "u in project before set_bnd");
  
  set_bnd(N, 1, u);
  set_bnd(N, 2, v);

  check_nans(N, v, "v in project");
  check_nans(N, u, "u in project");

}

void vel_step(const int N,
	      float* u, float* v, float* u0, float* v0,
	      const float visc, const float dt) {

  cout << "Adding source" << endl;

  add_source(N, u, u0, dt);
  add_source(N, v, v0, dt);

  check_nans(N, u, "u");
  check_nans(N, v, "v");

  check_nans(N, u0, "u0");
  check_nans(N, v0, "v0");

  SWAP(u0, u);
  diffuse(N, 1, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(N, 2, v, v0, visc, dt);

  check_nans(N, u, "u after diffuse");
  check_nans(N, v, "v after diffuse");

  check_nans(N, u0, "u0 after diffuse");
  check_nans(N, v0, "v0 after diffuse");
  
  project(N, u, v, u0, v0);

  cout << "Done with diffuse and project" << endl;

  SWAP(u0, u);
  SWAP(v0, v);

  check_nans(N, u, "u after diffuse project");
  check_nans(N, v, "v");

  check_nans(N, u0, "u0");
  check_nans(N, v0, "v0");
  
  // Advection of velocity?
  cout << "Advect 1" << endl;
  advect(N, 1, u, u0, u0, v0, dt);

  cout << "Done with advect 1" << endl;

  check_nans(N, u, "u");
  //check_nans(N, v);

  //check_nans(N, u0);
  //check_nans(N, v0);

  cout << "Advect 2" << endl;
  advect(N, 2, v, v0, u0, v0, dt);

  project(N, u, v, u0, v0);
  
}

void read_values(const int N, float* dens, float* u, float* v) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {

      dens[IX(i, j)] = 0;

      u[IX(i, j)] = 0;
      v[IX(i, j)] = 0;
    }

  }
  
}
