#pragma once

#define IX(i, j) ((i) + (N+2)*(j))
// Replace with std::swap?
#define SWAP(x, x0) {float* tmp = x; x0 = x; x = tmp;}

