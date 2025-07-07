//
// Created by huacheng on 11/10/22.
//

#include "linear.h"

void liner_op(double *x,int n,const double *a,const double *b,double pa,double pb)
{
    // x = pa*a + pb*b
    for(int k=0;k<n;++k)
    {
        x[k] = pa*a[k] + pb*b[k];
    }
}

void scaler_add(double *x,int n,const double *a,double pa)
{
    // x += pa*a
    for(int k=0;k<n;++k)
    {
        x[k] += pa*a[k];
    }
}


void scaler_move(double *x,int n,const double *a,double pa)
{
    // x = pa*a
    for(int k=0;k<n;++k)
    {
        x[k] = pa*a[k];
    }
}

void vec_scaler(double *x, int n, double px)
{
    // x = px*x
    for(int k=0;k<n;++k) x[k] *= px;
}

void vec_zero(double *x, int n)
{
    for(int k=0;k<n;++k)
    {
        x[k] = 0.0;
    }
}

void vec_number(double *x, int n, double number)
{
    for(int k=0;k<n;++k)
    {
        x[k] = number;
    }
}

int vec_isnan(double *x, int n)
{
    int nan_detected = 0;
    for(int k=0;k<n;++k)
    {
        if(isnan(x[k]))
        {
            nan_detected = 1;
            break;
        }
    }
    return nan_detected;
}

int vec_isinf(double *x,int n)
{
    int inf_detected = 0;
    for(int k=0;k<n;++k)
    {
        if(isinf(x[k]))
        {
            inf_detected = 1;
            break;
        }
    }
    return inf_detected;
}

double vec_sum(const double * x,int n)
{
    double z_res = 0.0;
    for(int k=0;k<n;++k)
    {
        z_res += x[k];
    }
    return z_res;
}


double vec_dot(const double *x,const double *y, int n)
{
    double z_res = 0;
    for(int k=0;k<n;++k)
    {
        z_res += x[k] * y[k];
    }
    return z_res;
}

void vec_copy(const double *x, double *y, int n)
{
    for(int k=0;k<n;++k) y[k] = x[k];
}

double det2d(const double *x)
{
    return x[0]*x[3] - x[1]*x[2];
}

void vec_norm(double *x, int n)
{
    double length = 0.;
    for(int k=0;k<n;++k)
    {
        length += x[k]*x[k];
    }

    length = sqrt(length);

    for(int k=0;k<n;++k)
    {
        x[k] /= length;
        if(isnan(x[k])) x[k] = 1.0;
    }
}

void vec_norm_2(double *x,int n, double tol)
{
    double y[3];
    scaler_move(y,n,x,tol);
    double length = 0.;
    for(int k=0;k<n;++k)
    {
        length += y[k]*y[k];
    }

    length = sqrt(length);
    for(int k=0;k<n;++k)
    {
        y[k] /= length;
        if(isnan(y[k])) y[k] = 1.0;
    }

    vec_copy(y,x,n);
}

double vec_len(double *x, int n)
{
    double res = 0.;
    for(int k=0;k<n;++k)
    {
        res += x[k] * x[k];
    }

    return sqrt(res);
}


double vec_dis(double *x, double *y, int n)
{
    double res = 0.;
    for(int k=0;k<n;++k)
    {
        res += (x[k] - y[k])*(x[k] - y[k]);
    }
    return sqrt(res);
}

int vec_cmp(const double *x, double y, int n)
{
    int ret = 0;
    for(int k=0;k<n;++k)
    {
        if(x[k] > y)
        {
            ret ++;
        }
    }
    return ret;
}

double second_inv(double * t)
{
    // ((t(XX)-t(YY))**2+(t(YY)-t(TH))**2+(t(TH)-t(XX))**2)/6.d0 + t(XY)**2
    // return second invariant of tensor
    double tXX = t[0], tYY = t[2], tXY = t[1], tTH = t[3];
    return ((tXX - tYY)*(tXX - tYY) + (tYY - tTH)*(tYY - tTH) + (tTH - tXX)*(tTH - tXX))/6.0 + tXY*tXY;
}

int bisearch(double x, double * x_list, int n)
{
    // if the x searched is out of the scope of x_list
    // bisearch will return the leftmost/rightmost interval
    int left = 0;
    int right = n - 1;

    if(x < x_list[left]) return left;
    if(x >= x_list[right]) return right - 1;

    while(left + 1 < right)
    {
        int mid = (left + right)/2;
        if(x_list[mid] > x)
        {
            right = mid;
        } else
        {
            left = mid;
        }
    }
    return left;
}

int abilinear(const double *x,const double *y,double *xl)
{
    const double bfactor[4][5] = {
            { 1.0,-1.0, 1.0,-1.0, 0.0},
            {-1.0, 1.0, 1.0,-1.0, 0.0},
            {-1.0,-1.0, 1.0, 1.0, 0.0},
            { 1.0, 1.0, 1.0, 1.0,-4.0}
    };

    double a[4] = {0.}, b[4] = {0.};
    for(int k=0;k<4;++k)
    {
        for(int j=0;j<5;++j)
        {
            a[k] += bfactor[k][j] * x[j];
            b[k] += bfactor[k][j] * y[j];
        }
    }

    double tA = -a[0]*b[2] + a[2]*b[0];
    double tB = -a[0]*b[3] + a[3]*b[0] - a[1]*b[2] + a[2]*b[1];
    double tC = -a[1]*b[3] + a[3]*b[1];

    double tDelta = tB*tB - 4.0*tA*tC;

    double xlt[4] = {0.0};
    int ret_var = 0;

    if(tDelta < 0.)
    {
        ret_var = 0;
    } else if(fabs(tDelta) < 1e-24)
    {
        xlt[1] = -tB/(2.*tA);
        xlt[0] = -1.0 * (a[2]*xlt[1] + a[3])/(a[0]*xlt[1] + a[1]);
        ret_var = 1;
        if(fabs(xlt[0]) > 1.0 || fabs(xlt[1]) > 1.0)
        {
            xlt[0] = 0.0;
            xlt[1] = 0.0;
            ret_var --;
        }

    } else
    {
        xlt[1] = (-tB + sqrt(tDelta))/(2.*tA);
        xlt[0] = -1.0 * (a[2]*xlt[1] + a[3])/(a[0]*xlt[1] + a[1]);

        xlt[3] = (-tB - sqrt(tDelta))/(2.*tA);
        xlt[2] = -1.0 * (a[2]*xlt[3] + a[3])/(a[0]*xlt[3] + a[1]);

        ret_var = 2;
        if(fabs(xlt[0]) > 1.0 || fabs(xlt[1]) > 1.0)
        {
            xlt[0] = 0.0;
            xlt[1] = 0.0;
            ret_var --;
        }

        if(fabs(xlt[2]) > 1.0 || fabs(xlt[3]) > 1.0)
        {
            xlt[2] = 0.;
            xlt[3] = 0.;
            ret_var --;
        }
    }

    if((fabs(xlt[2]) + fabs(xlt[3])) > (fabs(xlt[1]) + fabs(xlt[0])))
    {
        xl[0] = xlt[2];
        xl[1] = xlt[3];
    } else
    {
        xl[0] = xlt[0];
        xl[1] = xlt[1];
    }

    return ret_var;
}