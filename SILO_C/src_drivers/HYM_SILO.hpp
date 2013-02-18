//============================================================================//
/*

Clayton Myers
HYM_SILO.hpp
Created:  18 August 2009
Modified: 17 September 2009

Header file to centralize HYM_SILO parameters and function declarations.

*/
//============================================================================//
//============================================================================//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <unistd.h>  // For chdir and getcwd in Basic_Functions.cpp
using namespace std;

//============================================================================//

#define DB_USE_LEGACY_DTPTR  // Required to work with SILO 4.7.1 and later
#include <silo.h>
#include <Basic_Functions.hpp>

//============================================================================//
// Basic constants:
const int ndims = 3;  // Number of dimensions
const int nvars = 5;  // Number of variables (p,n,B,v,J,...)
const int Nghost = 0; // Number of SILO ghost zones added on each end (in phi)
const int exclude_origin=0; //Whether to interpolate vectors around r=0 or not. If r=0 not part of domain, then =0
const int half_cyl=1; //whether to relabel the final s co-ord as pi in a half cylindrical geometry

const float pi = 4.0*atan(1.0); // Natural constant

//============================================================================//
//============================================================================//
