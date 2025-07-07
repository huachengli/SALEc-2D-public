//
// Created by huach on 2/7/2023.
//
// functions used for calculate parameters of Nine-point stencil for Laplacian operator
// at non-uniform grids
// suffix 2d9 means nine points/shapes function in 2d
// shape functions are generated in Lagrange Polynomial
// the interpolate_fn can also been rewritten in macros (future)
//

#ifndef SALE_REBUILD_INTERPOLATE_2D9_H
#define SALE_REBUILD_INTERPOLATE_2D9_H

#include <stdio.h>
#include <math.h>

typedef double (*Interpolation2d9_t)(const double *_x);

// interpolation function for 2d9
extern Interpolation2d9_t Ip2d9[9];
extern Interpolation2d9_t Ip2d9_X[9][2];
extern Interpolation2d9_t Ip2d9_XT[2][9];
extern Interpolation2d9_t Ip2d9_XX[9][3];
// 2d9 Gaussian quadrature
extern const double GIPS2d9[9][2];
extern const double GIWS2d9[9];

// following functions from intepolate_2d9 have been tested.
double Int2d9(const double *Xi,const double *Yi);
double Laplace2d9(const double * Xi, const double *Yi, double gamma);
double Laplace2d9cyl(const double * Xi, const double *Yi, double gamma);

double Interpolate2d9_Fdx(const double *Xi,const double *Yi, const double *xl,int dim);

// gaussian quadrature points
extern const double GIPS2d16[16][2];
extern const double GIWS2d16[16];

#endif //SALE_REBUILD_INTERPOLATE_2D9_H
