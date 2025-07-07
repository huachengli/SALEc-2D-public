//
// Created by huacheng on 11/10/22.
//

#ifndef SALE_REBUILD_INTERPOLATE_H
#define SALE_REBUILD_INTERPOLATE_H
#include "interpolate_fn.h"
#include "linear.h"

double Integrate2d(const double *Xi,const double *Yi);
double Interpolate2d(const double * Xi, const double * xl);
double Interpolate2d_k(double * Xi[], const double * xl, int k);
double Integrate2dcyl(const double *Xi,const double *Yi);
void CrossLine2dcyl(const double *Xi1, const double *Xi2, double *_s);
void Polygon2dcyl(const double *Xi, int n, double * _p);
void Polygon2dRotate(double *Xi,int n);
void Polygon2dRotate2(double *Xi,int n,int m);
double Integrate2d_Fdx(const double *Xi, const double *Fi, int dim);
double Integrate2dcyl_Fdx(const double *Xi, const double *Fi, int dim);
void GradVec2dcyl(const double *Xi, const double *Fi, double * _grad);
void GradVec2d(const double *Xi, const double *Fi, double * _grad);
double DivVec2dcyl(const double *Xi, const double *Fi);
double DivVec2d(const double *Xi, const double *Fi);
void Grad2d(const double *Xi, const double *Fi, double * _grad);
void Grad2dcyl(const double *Xi, const double *Fi, double * _grad);

void Interpolate2dvec2(const double * Xi, const double * xl, double * r, int n);
void Interpolate2dvec(const double * Xi, const double * xl, double * r);
#endif //SALE_REBUILD_INTERPOLATE_H
