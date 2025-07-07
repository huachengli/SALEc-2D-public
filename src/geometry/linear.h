//
// Created by huacheng on 11/10/22.
//

#ifndef SALE_REBUILD_LINEAR_H
#define SALE_REBUILD_LINEAR_H
#include <math.h>
void liner_op(double *x,int n,const double *a,const double *b,double pa,double pb);
void scaler_add(double *x,int n,const double *a,double pa);
void scaler_move(double *x,int n,const double *a,double pa);
void vec_zero(double *x, int n);
void vec_number(double *x, int n, double number);
int vec_isnan(double *x, int n);
int vec_isinf(double *x,int n);
void vec_copy(const double *x, double *y, int n);
void vec_scaler(double *x, int n, double px);
double vec_sum(const double * x,int n);
double vec_dot(const double *x,const double *y, int n);
double det2d(const double *x);
void vec_norm(double *x, int n);
void vec_norm_2(double *x,int n, double tol);
double vec_len(double *x, int n);
double vec_dis(double *x, double *y, int n);
int vec_cmp(const double *x, double y, int n);
double second_inv(double * t);

#define X 0
#define Y 1
#define Z 2

int bisearch(double x, double * x_list, int n);
int abilinear(const double *x,const double *y,double *xl);
#endif //SALE_REBUILD_LINEAR_H
