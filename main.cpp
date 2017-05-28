#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "fluid.h"
#include "utils.h"

#include <iostream>

using namespace std;

TEST_CASE("One step on 20 by 20 grid, constant velocity") {
  int N = 200;
  int size = (N+2)*(N+2);

  float u[size], v[size], u_prev[size], v_prev[size];
  double dt = 0.1;
  double visc = 0.00;
  double diff = 100.00;

  int cube_bl_x = 50;
  int cube_bl_y = 150;
  int cube_length = 30;

  float cube_vel_x = 10.0;
  float cube_vel_y = 0.0;

  float dens[size], dens_prev[size];

  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {

      u_prev[IX(i, j)] = cube_vel_x;
      v_prev[IX(i, j)] = 0.0;

      u[IX(i, j)] = cube_vel_x;
      v[IX(i, j)] = 0.0;
      
      if ((cube_bl_x <= i && i <= cube_bl_x + cube_length) &&
	  (cube_bl_y <= j && j <= cube_bl_y + cube_length)) {

	dens_prev[IX(i, j)] = 255.0;
	dens[IX(i, j)] = 255.0;

      } else {

	dens_prev[IX(i, j)] = 0;
	dens[IX(i, j)] = 0;

      }
    }

  }

  set_bnd(N, 0, dens_prev);
  set_bnd(N, 1, u_prev);
  set_bnd(N, 2, v_prev);
  
  set_bnd(N, 0, dens);
  set_bnd(N, 1, u);
  set_bnd(N, 2, v);

  cout << "Initialized density" << endl;
  //visualize_density(N, dens);

  cout << "Initialized u" << endl;
  //visualize_density(N, u);

  cout << "Initialized v" << endl;
  //visualize_density(N, v);
  
  vel_step ( N, u, v, u_prev, v_prev, visc, dt );

  // visualize_density(N, u);
  // visualize_density(N, v);

  cout << "Did vel step" << endl;

  dens_step ( N, dens, dens_prev, u, v, diff, dt );

  //read_values( N, dens_prev, u_prev, v_prev );

  cout << "Did dens step" << endl;

  //visualize_density(N, dens);

  //visualize_density(N, dens);

}

TEST_CASE("Full simulation 200 by 200 grid") {
  int N = 200;
  int size = (N+2)*(N+2);

  float u[size], v[size], u_prev[size], v_prev[size];
  double dt = 0.1;
  double visc = 0.00;
  double diff = 100.00;

  int cube_bl_x = 50;
  int cube_bl_y = 150;
  int cube_length = 30;

  float cube_vel_x = 0.1;
  float cube_vel_y = 0.1;

  float dens[size], dens_prev[size];

  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {

      // u_prev[IX(i, j)] = cube_vel_x; //cube_vel_x;//random_float(0, 1);
      // v_prev[IX(i, j)] = -1*cube_vel_y; //cube_vel_y;//-1*u_prev[i];//random_float(0, 
      //      1);

      u_prev[IX(i, j)] = cube_vel_x + random_float(0, cube_vel_x / 3.0);//random_float(0, 1);
      v_prev[IX(i, j)] = -1*cube_vel_y + random_float(0, cube_vel_y / 3.0);//random_float(0      
      u[IX(i, j)] = cube_vel_x + random_float(0, cube_vel_x / 3.0);//random_float(0, 1);
      v[IX(i, j)] = -1*cube_vel_y + random_float(0, cube_vel_y / 3.0);//random_float(0, 1);
      
      if ((cube_bl_x <= i && i <= cube_bl_x + cube_length) &&
	  (cube_bl_y <= j && j <= cube_bl_y + cube_length)) {

	dens_prev[IX(i, j)] = 255.0; //255.0; //random_float(0, 255);
	dens[IX(i, j)] = 255.0; //random_float(0, 255);

	// u_prev[IX(i, j)] = cube_vel_x / 2.0; //cube_vel_x;//random_float(0, 1);
	// v_prev[IX(i, j)] = -1*cube_vel_y / 2.0; //cube_vel_y;//-1*u_prev[i];//random_float(0, 1);

	// u[IX(i, j)] = cube_vel_x;//random_float(0, 1);
	// v[IX(i, j)] = -1*cube_vel_y;//random_float(0, 1);

      } else {

	dens_prev[IX(i, j)] = 0;
	dens[IX(i, j)] = 0;

	// u_prev[IX(i, j)] = 0;
	// v_prev[IX(i, j)] = 0;

	// u[IX(i, j)] = 0;
	// v[IX(i, j)] = 0;
      }
    }

  }

  set_bnd(N, 0, dens_prev);
  set_bnd(N, 1, u_prev);
  set_bnd(N, 2, v_prev);
  
  set_bnd(N, 0, dens);
  set_bnd(N, 1, u);
  set_bnd(N, 2, v);

  cout << "Initialized density" << endl;
  //visualize_density(N, dens);

  cout << "Initialized u" << endl;
  //visualize_density(N, u);

  cout << "Initialized v" << endl;
  //visualize_density(N, v);
  
  int i = 0;
  int max = 10;
  while (i < max) {
    vel_step ( N, u, v, u_prev, v_prev, visc, dt );

    // visualize_density(N, u);
    // visualize_density(N, v);

    cout << "Did vel step" << endl;

    dens_step ( N, dens, dens_prev, u, v, diff, dt );

    //read_values( N, dens_prev, u_prev, v_prev );

    cout << "Did dens step" << endl;

    i++;

    cout << "i = " << i << endl;

    //visualize_density(N, dens);

  }

  //visualize_density(N, dens);

  
}
