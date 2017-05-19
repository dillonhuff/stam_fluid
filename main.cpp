
#define IX(i, j) ((i) + (N+2)*(j))

void add_source(const int N, float* x, float* s, const float dt) {
  int size = (N+2)*(N+2);

  for (int i = 0; i < size; i++) {
    x[i] += dt*s[i];
  }

  
}

void set_bnd(const int N, const double d, float* x) {
  
}

void diffuse(const int N, const int b, float* x, float* x0, float diff, float dt) {
  float a = dt*diff*N*N;

  for (int k = 0; k < 20; k++) {

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

int main() {
  // int N = 5;
  // int size = (N+2)*(N+2);
  // int u[size], v[size], u_prev[size], v_prev[size];
  // int dens[size], dens_prev[size];

  

  return 0;
}
