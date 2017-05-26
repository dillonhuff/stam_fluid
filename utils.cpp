#include "utils.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

void check_nans(const int N, float* x) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      if (isnan(x[i])) {
	cout << "( " << i << " , " << j << " )" << endl;
	assert(false);
      }
    }
  }
}
