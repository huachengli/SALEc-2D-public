//
// Created by huacheng on 2020/12/22.
//

#ifndef SALEC_VARIABLE_H
#define SALEC_VARIABLE_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>


#define NIpB 4
#define NIpV 8
#define NGI2d 4
#define NGI3d 8
#define BETA 0.57735026918962576451

// points of Gaussian integral for 2*2 in 2d
extern const double GIPS2d[NGI2d][2];
// points of Gaussian integral for 2*2*2 in 3d
extern const double GIPS3d[NGI3d][3];
// weights of Gaussian integral for 2*2 in 2d
extern const double GIWS2d[NGI2d];
// weights of Gaussian integral for 2*2*2 in 3d
extern const double GIWS3d[NGI3d];

typedef double (*Interpolation2d9_t)(const double *_x);

// interpolate function for a Plane with 4 points
extern Interpolation2d9_t IpB[NIpB];
extern Interpolation2d9_t IpB_X[NIpB][2];
extern Interpolation2d9_t IpB_XT[2][NIpB];
// interpolation function for Volume with in 8 points
extern Interpolation2d9_t IpV[NIpV];
extern Interpolation2d9_t IpV_X[NIpV][3];

#endif //SALEC_VARIABLE_H
