//
// Created by huacheng on 11/10/22.
//

#include "interpolate.h"
void XgIpB(double *Xg, int n,const double *Xi, const double *xl)
{
    vec_zero(Xg,n);
    for(int k=0;k<NIpB;++k)
    {
        scaler_add(Xg,n,Xi + k*n,IpB[k](xl));
    }
}

void XgIpB_X(double *Xg, int n,const double *Xi, const double *xl, int dim)
{
    vec_zero(Xg,n);
    vec_zero(Xg,n);
    for(int k=0;k<NIpB;++k)
    {
        scaler_add(Xg,n,Xi + k*n,IpB_X[k][dim](xl));
    }
}

double Interpolate_fn_2d(const double * Xi, Interpolation2d9_t * fn, const double * xl)
{
    double res_value = 0;
    if(NULL == fn) fn = IpB;
    for(int k=0;k<NIpB;++k)
    {
        res_value+= fn[k](xl)*Xi[k];
    }
    return res_value;
}

double Interpolate2d(const double * Xi, const double * xl)
{
    double res_value = 0;
    for(int k=0;k<NIpB;++k)
    {
        res_value+= IpB[k](xl)*Xi[k];
    }
    return res_value;
}

double Interpolate2d_k(double * Xi[], const double * xl, int k)
{
    return Interpolate2d((double []){Xi[0][k],Xi[1][k],Xi[2][k],Xi[3][k]},xl);
}


double jacobi2d(const double *Xi, const double *xl)
{
    double jacobi_matrix[4];
    jacobi_matrix[0] = Interpolate_fn_2d((double []){Xi[0],Xi[2],Xi[4],Xi[6]},IpB_XT[0],xl);
    jacobi_matrix[1] = Interpolate_fn_2d((double []){Xi[0],Xi[2],Xi[4],Xi[6]},IpB_XT[1],xl);
    jacobi_matrix[2] = Interpolate_fn_2d((double []){Xi[1],Xi[3],Xi[5],Xi[7]},IpB_XT[0],xl);
    jacobi_matrix[3] = Interpolate_fn_2d((double []){Xi[1],Xi[3],Xi[5],Xi[7]},IpB_XT[1],xl);
    return det2d(jacobi_matrix);
}

double Integrate2d(const double *Xi,const double *Yi)
{
    double res_value = 0;
    for(int j=0;j<NGI2d;++j)
    {
        res_value += GIWS2d[j]* Interpolate2d(Yi,GIPS2d[j]) * jacobi2d(Xi,GIPS2d[j]);
    }
    return res_value;
}

double Integrate2dcyl(const double *Xi,const double *Yi)
{
    /*
     * the difference between Integrate2dcyl and Integrate2d is
     * Integrate2dcyl calculate integral of y*x*dx*dy
     * Integrate calculate integral of y*dx*dy
     */
    double res_value = 0;
    for(int j=0;j<NGI2d;++j)
    {
        double local_r = Interpolate2d((double []){Xi[0],Xi[2],Xi[4],Xi[6]},GIPS2d[j]);
        local_r = fabs(local_r);
        res_value += GIWS2d[j]*Interpolate2d(Yi,GIPS2d[j])*jacobi2d(Xi,GIPS2d[j])*local_r;
    }
    return res_value;
}

void CrossLine2dcyl(const double *Xi1, const double *Xi2, double *_s)
{
    double r0 = (Xi2[0] + Xi1[0])/2.0;
    /*
     * r0 should be always positive
     * it is the distance between origin and midpoint of this line in x direction
     * the normal vector is written into _s
     */
    r0 = fabs(r0);

    _s[0] = -1.0*(Xi2[1] - Xi1[1])*r0;
    _s[1] = (Xi2[0] - Xi1[0])*r0;
    /*
     * _s is the surface represented by line (Xi2 - Xi1) in x-z plane
     */
}

void Polygon2dcyl(const double *Xi, int n, double * _p)
{
    vec_zero(_p,2);
    double tmp_line[2];
    for(int k=1;k<n;++k)
    {
        CrossLine2dcyl(Xi+2*(k-1),Xi+2*k,tmp_line);
        scaler_add(_p,2,tmp_line,1.0);
    }
}

void Polygon2dRotate(double *Xi,int n)
{
    double * Yi = (double *) malloc(sizeof(double)*n*2);
    for(int k=0;k<n;k++)
    {
        int ia = k;
        int ib = (k+1)%n;
        Yi[ia*2] = Xi[ib*2];
        Yi[ia*2 + 1] = Xi[ib*2 + 1];
    }
    for(int k=0;k<n;k++)
    {
        Xi[2*k] = Yi[2*k];
        Xi[2*k+1] = Yi[2*k+1];
    }
    free(Yi);
}

void Polygon2dRotate2(double *Xi,int n,int m)
{
    double * Yi = (double *) malloc(sizeof(double)*n*2);
    for(int k=0;k<n;k++)
    {
        int ia = k;
        int ib = (k+m)%n;
        Yi[ia*2] = Xi[ib*2];
        Yi[ia*2 + 1] = Xi[ib*2 + 1];
    }
    for(int k=0;k<n;k++)
    {
        Xi[2*k] = Yi[2*k];
        Xi[2*k+1] = Yi[2*k+1];
    }
    free(Yi);
}

double Integrate2d_Fdx(const double *Xi, const double *Fi, int dim)
{
    double res_value = 0;
    double local_field[8];
    vec_copy(Xi,local_field,8);
    for(int k=0;k<4;++k)
    {
        local_field[k*2 + dim] = Fi[k];
    }

    for(int j=0;j<NGI2d;++j)
    {
        res_value += GIWS2d[j] * jacobi2d(local_field,GIPS2d[j]);
    }
    return res_value;
}

double Integrate2dcyl_Fdx(const double *Xi, const double *Fi, int dim)
{
    double res_value = 0;
    double local_field[8];
    vec_copy(Xi,local_field,8);
    for(int k=0;k<4;++k)
    {
        local_field[k*2 + dim] = Fi[k];
    }

    for(int j=0;j<NGI2d;++j)
    {
        double local_r = Interpolate2d((double []){Xi[0],Xi[2],Xi[4],Xi[6]},GIPS2d[j]);
        res_value += GIWS2d[j] * jacobi2d(local_field,GIPS2d[j]) * fabs(local_r);
    }
    return res_value;
}

void GradVec2dcyl(const double *Xi, const double *Fi, double * _grad)
{

    _grad[0] = Integrate2dcyl_Fdx(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]},0);
    _grad[1] = Integrate2dcyl_Fdx(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]},1);
    _grad[2] = Integrate2dcyl_Fdx(Xi,(double []){Fi[1],Fi[3],Fi[5],Fi[7]},0);
    _grad[3] = Integrate2dcyl_Fdx(Xi,(double []){Fi[1],Fi[3],Fi[5],Fi[7]},1);
    _grad[4] = Integrate2d(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]});
}

void GradVec2d(const double *Xi, const double *Fi, double * _grad)
{
    _grad[0] = Integrate2d_Fdx(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]},0);
    _grad[1] = Integrate2d_Fdx(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]},1);
    _grad[2] = Integrate2d_Fdx(Xi,(double []){Fi[1],Fi[3],Fi[5],Fi[7]},0);
    _grad[3] = Integrate2d_Fdx(Xi,(double []){Fi[1],Fi[3],Fi[5],Fi[7]},1);
}

double DivVec2dcyl(const double *Xi, const double *Fi)
{
    double res_value = 0;
    res_value += Integrate2dcyl_Fdx(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]},0);
    res_value += Integrate2dcyl_Fdx(Xi,(double []){Fi[1],Fi[3],Fi[5],Fi[7]},1);
    res_value += Integrate2d(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]});
    return res_value;
}

double DivVec2d(const double *Xi, const double *Fi)
{
    double res_value = 0;
    res_value += Integrate2d_Fdx(Xi,(double []){Fi[0],Fi[2],Fi[4],Fi[6]},0);
    res_value += Integrate2d_Fdx(Xi,(double []){Fi[1],Fi[3],Fi[5],Fi[7]},1);
    return res_value;
}

void Grad2d(const double *Xi, const double *Fi, double * _grad)
{
    _grad[0] = Integrate2d_Fdx(Xi,Fi,0);
    _grad[1] = Integrate2d_Fdx(Xi,Fi,1);
}

void Grad2dcyl(const double *Xi, const double *Fi, double * _grad)
{
    _grad[0] = Integrate2dcyl_Fdx(Xi,Fi,0);
    _grad[1] = Integrate2dcyl_Fdx(Xi,Fi,1);
}

void Interpolate2dvec2(const double * Xi, const double * xl, double * r, int n)
{
    vec_zero(r,n);
    for(int k=0;k<n;++k)
    {
        r[k] = Interpolate2d((double []){Xi[k+ 0*n], Xi[k+ 1*n], Xi[k+ 2*n], Xi[k+ 3*n]},xl);
    }
}

void Interpolate2dvec(const double * Xi, const double * xl, double * r)
{
    r[0] = Interpolate2d((double []){Xi[0], Xi[2], Xi[4], Xi[6]}, xl);
    r[1] = Interpolate2d((double []){Xi[1], Xi[3], Xi[5], Xi[7]}, xl);
}