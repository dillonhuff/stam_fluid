#pragma once

#define IX(i, j) ((i) + (N+2)*(j))
// Replace with std::swap?
#define SWAP(x, x0) {float* tmp = x; x0 = x; x = tmp;}

#include <string>

void check_nans(const int N, float* x, const std::string& s);

float array_min(const int N, float* densities);

float array_max(const int N, float* densities);
