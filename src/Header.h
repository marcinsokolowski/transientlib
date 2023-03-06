// Copyright 2003 by Bogumil Pilecki and Dorota Szczygiel

#ifndef HEAD_H
#define HEAD_H

#include <iostream>
#include <fstream>
#include <sstream>	// -> Astroccd.cpp
#include <iomanip>	// -> Astroccd.cpp
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <ctime>	// -> Astroccd.cpp

// optionaly :
#include "mathdefs.h"
// #define PI_VALUE 3.14159265358;

using namespace std;

typedef signed short fshort;


const int DEFAULT_STAR_LIST_SIZE = 200;


// const double PI = 3.14159265358;

typedef char map_type;

const map_type BOUNDARY_STAR = 2;	// this pixel is near the border of an image
const map_type HOT_PIXEL = 1;		// this pixel always has a high value


// number of stars that that are used to find tranformation between starlists

const int MIN_STARS_TO_TRANSFORM      =  1;
const int PREFERED_STARS_TO_TRANSFORM = 20;
const int MAX_STARS_TO_TRANSFORM      = 50;

// x, y, mag tolerance before transformation

const float X_TOLERANCE   = 2.5;
const float Y_TOLERANCE   = 2.5;
const float MAG_TOLERANCE = 0.2;


// x, y, mag tolerance after transformation

const float MAX_X_TOLERANCE   = 1.5;
const float MAX_Y_TOLERANCE   = 1.5;
const float MIN_MAG_TOLERANCE = 0.15;  //0.15;
const float MAX_MAG_TOLERANCE = 1.5; //0.35;


extern const double UP_NORTH;
extern const double UP_SOUTH;
extern const double UP_WEST;
extern const double UP_EAST;
         

#endif
