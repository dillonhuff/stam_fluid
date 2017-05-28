#include "utils.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

void check_nans(const int N, float* x, const std::string& s) {
  for (int i = 0; i < N + 2; i++) {
    for (int j = 0; j < N + 2; j++) {
      if (isnan(x[IX(i, j)])) {
	cout << "NaN in " << s << " ";
	cout << "( " << i << " , " << j << " )" << endl;
	assert(false);
      }

      if (isinf(x[IX(i, j)])) {
	cout << "Infinity in " << s << " ";
	cout << "( " << i << " , " << j << " )" << endl;
	assert(false);
      }
    }
  }
}

float array_min(const int N, float* densities) {
  auto min_elem = min_element(densities, densities + (N+2)*(N+2));
  return *min_elem;
}

float array_max(const int N, float* densities) {
  auto max_elem = max_element(densities, densities + (N+2)*(N+2));
  return *max_elem;
}

float random_float(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

