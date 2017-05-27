#include "utils.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

void check_nans(const int N, float* x, const std::string& s) {
  for (int i = 0; i < N + 2; i++) {
    for (int j = 0; j < N + 2; j++) {
      if (isnan(x[IX(i, j)])) {
	cout << "In " << s << " ";
	cout << "( " << i << " , " << j << " )" << endl;
	assert(false);
      }
    }
  }
}
