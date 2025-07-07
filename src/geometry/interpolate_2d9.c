//
// Created by huach on 2/7/2023.
//

#include "interpolate_2d9.h"

// macros result example
//double Ip2d9i0(const double *_x)
//{
//    double xi  = _x[0];
//    double eta = _x[1];
//    int i0= xli2d9[0][0], i1 = xli2d9[0][1], i2= xli2d9[0][2];
//    int j0= xlj2d9[0][0], j1 = xlj2d9[0][1], j2= xlj2d9[0][2];
//    double xf = Lagrange(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange(xi ,xl2d9[i0][0],xl2d9[i2][0]);
//    double yf = Lagrange(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange(eta,xl2d9[j0][1],xl2d9[j2][1]);
//    return  xf*yf;
//}
//double Ip2d9f0_X1(const double *_x)
//{
//    double xi = _x[0], eta = _x[1];
//    int i0 = xli2d9[0][0], i1 = xli2d9[0][1], i2 = xli2d9[0][2];
//    int j0 = xlj2d9[0][0], j1 = xlj2d9[0][1], j2 = xlj2d9[0][2];
//    double xf_x = 1.0 / ((xl2d9[i0][0]) - (xl2d9[i1][0]))*(((xi) - (xl2d9[i2][0])) / ((xl2d9[i0][0]) - (xl2d9[i2][0])))
//    + (((xi) - (xl2d9[i1][0])) / ((xl2d9[i0][0]) - (xl2d9[i1][0]))) * 1.0 / ((xl2d9[i0][0]) - (xl2d9[i2][0])));
//    double yf = (((eta) - (xl2d9[j1][1])) / ((xl2d9[j0][1]) - (xl2d9[j1][1]))) *
//                (((eta) - (xl2d9[j2][1])) / ((xl2d9[j0][1]) - (xl2d9[j2][1])));
//    return xf_x * yf;
//}

#define Lagrange(x,x0,x1) (((x) - (x1))/((x0) - (x1)))
#define Lagrange_x(x,x0,x1) (1.0/((x0) - (x1)))
const double xl2d9[9][2] =
        {
                {-1.,-1.},{0.,-1.},{1.,-1.},
                {-1., 0.},{0., 0.},{1., 0.},
                {-1., 1.},{0., 1.},{1., 1.},
        };
const int xli2d9[9][3] = {
        {0,1,2},
        {1,2,0},
        {2,0,1},
        {3,4,5},
        {4,5,3},
        {5,3,4},
        {6,7,8},
        {7,8,6},
        {8,6,7}
};

const int xlj2d9[9][3] = {
        {0,3,6},
        {1,4,7},
        {2,5,8},
        {3,0,6},
        {4,7,1},
        {5,8,2},
        {6,0,3},
        {7,1,4},
        {8,2,5}
};

#define Ip2d9_list(func_d) \
    func_d(0)              \
    func_d(1)              \
    func_d(2)              \
    func_d(3)              \
    func_d(4)              \
    func_d(5)              \
    func_d(6)              \
    func_d(7)              \
    func_d(8)

#define Ip2d9_array(func_d) \
    func_d(0),              \
    func_d(1),              \
    func_d(2),              \
    func_d(3),              \
    func_d(4),              \
    func_d(5),              \
    func_d(6),              \
    func_d(7),              \
    func_d(8)



#define _Ip2d9fi_setij(k) double xi  = _x[0],eta = _x[1]; \
    int i0= xli2d9[k][0], i1 = xli2d9[k][1], i2= xli2d9[k][2]; \
    int j0= xlj2d9[k][0], j1 = xlj2d9[k][1], j2= xlj2d9[k][2];

#define _Ip2d9fi(k) \
    double xf = Lagrange(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange(xi ,xl2d9[i0][0],xl2d9[i2][0]); \
    double yf = Lagrange(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange(eta,xl2d9[j0][1],xl2d9[j2][1]); \
    return xf*yf;
#define Ip2d9fi(k) double Ip2d9i##k(const double * _x){_Ip2d9fi_setij(k) _Ip2d9fi(k) }

#define _Ip2d9fi_X1(k) \
    double xf_x = Lagrange_x(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange(xi ,xl2d9[i0][0],xl2d9[i2][0]) \
    + Lagrange(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange_x(xi ,xl2d9[i0][0],xl2d9[i2][0]); \
    double yf = Lagrange(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange(eta,xl2d9[j0][1],xl2d9[j2][1]); \
    return xf_x*yf;

#define _Ip2d9fi_X2(k) \
    double xf = Lagrange(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange(xi ,xl2d9[i0][0],xl2d9[i2][0]); \
    double yf_x = Lagrange_x(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange(eta,xl2d9[j0][1],xl2d9[j2][1])\
    + Lagrange(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange_x(eta,xl2d9[j0][1],xl2d9[j2][1]);            \
    return xf*yf_x;

#define Ip2d9fi_X1(k) double Ip2d9i##k##_X1(const double *_x){ \
    _Ip2d9fi_setij(k)                                        \
    _Ip2d9fi_X1(k)                                           \
}

#define Ip2d9fi_X2(k) double Ip2d9i##k##_X2(const double *_x){ \
    _Ip2d9fi_setij(k)                                        \
    _Ip2d9fi_X2(k)                                           \
}

#define _Ip2d9fi_XX(k) \
    double x2f_xx = 2.0*Lagrange_x(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange_x(eta,xl2d9[j0][1],xl2d9[j2][1]); \
    double x1f_xx = 2.0*Lagrange_x(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange_x(xi ,xl2d9[i0][0],xl2d9[i2][0]); \
    double x1f_x = Lagrange_x(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange(xi ,xl2d9[i0][0],xl2d9[i2][0]) \
    + Lagrange(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange_x(xi ,xl2d9[i0][0],xl2d9[i2][0]); \
    double x2f_x = Lagrange_x(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange(eta,xl2d9[j0][1],xl2d9[j2][1])\
    + Lagrange(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange_x(eta,xl2d9[j0][1],xl2d9[j2][1]);             \
    double x2f = Lagrange(eta,xl2d9[j0][1],xl2d9[j1][1])* Lagrange(eta,xl2d9[j0][1],xl2d9[j2][1]);            \
    double x1f = Lagrange(xi ,xl2d9[i0][0],xl2d9[i1][0])* Lagrange(xi ,xl2d9[i0][0],xl2d9[i2][0]);

#define Ip2d9fi_X1X1(k) double Ip2d9i##k##_X1X1(const double *_x){ \
    _Ip2d9fi_setij(k) \
    _Ip2d9fi_XX(k)                                                 \
    return x1f_xx*x2f;                                              \
    }

#define Ip2d9fi_X2X2(k) double Ip2d9i##k##_X2X2(const double *_x){ \
    _Ip2d9fi_setij(k) \
    _Ip2d9fi_XX(k)                                                 \
    return x1f*x2f_xx;                                              \
    }

#define Ip2d9fi_X1X2(k) double Ip2d9i##k##_X1X2(const double *_x){ \
    _Ip2d9fi_setij(k) \
    _Ip2d9fi_XX(k)                                                 \
    return x1f_x*x2f_x;                                            \
    }

Ip2d9_list(Ip2d9fi)
Ip2d9_list(Ip2d9fi_X1)
Ip2d9_list(Ip2d9fi_X2)
Ip2d9_list(Ip2d9fi_X1X1)
Ip2d9_list(Ip2d9fi_X2X2)
Ip2d9_list(Ip2d9fi_X1X2)

#define _Interpolation_Ip2d9(k) Ip2d9i##k

Interpolation2d9_t Ip2d9[9] = {
    Ip2d9_array(_Interpolation_Ip2d9)
};

Interpolation2d9_t Ip2d9_XT[2][9] = {
    {Ip2d9i0_X1,
    Ip2d9i1_X1,
    Ip2d9i2_X1,
    Ip2d9i3_X1,
    Ip2d9i4_X1,
    Ip2d9i5_X1,
    Ip2d9i6_X1,
    Ip2d9i7_X1,
    Ip2d9i8_X1},
    {Ip2d9i0_X2,
    Ip2d9i1_X2,
    Ip2d9i2_X2,
    Ip2d9i3_X2,
    Ip2d9i4_X2,
    Ip2d9i5_X2,
    Ip2d9i6_X2,
    Ip2d9i7_X2,
    Ip2d9i8_X2}
};

#define _Interpolation_Ip2d9_X(k) {Ip2d9i##k##_X1, Ip2d9i##k##_X2}

Interpolation2d9_t Ip2d9_X[9][2]= {
    Ip2d9_array(_Interpolation_Ip2d9_X)
};

#define _Interpolation_Ip2d9_XX(k) {Ip2d9i##k##_X1X1, Ip2d9i##k##_X1X2, Ip2d9i##k##_X2X2}

Interpolation2d9_t Ip2d9_XX[9][3]={
        Ip2d9_array(_Interpolation_Ip2d9_XX)
};

#define SQRT_3O5 0.77459666924148337703585307995648
const double GIPS2d9[9][2] = {
        {-SQRT_3O5, -SQRT_3O5},
        {-SQRT_3O5, 0.},
        {-SQRT_3O5,SQRT_3O5},
        {0., -SQRT_3O5},
        {0., 0.},
        {0.,SQRT_3O5},
        {SQRT_3O5, -SQRT_3O5},
        {SQRT_3O5, 0.},
        {SQRT_3O5,SQRT_3O5},
};

const double GIWS2d9[9] = {
        (5.0/9.0) * (5.0/9.0),
        (5.0/9.0) * (8.0/9.0),
        (5.0/9.0) * (5.0/9.0),
        (8.0/9.0) * (5.0/9.0),
        (8.0/9.0) * (8.0/9.0),
        (8.0/9.0) * (5.0/9.0),
        (5.0/9.0) * (5.0/9.0),
        (5.0/9.0) * (8.0/9.0),
        (5.0/9.0) * (5.0/9.0)
};

double Interpolate2d9(const double * Xi, const double * xl){
    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += Ip2d9[k](xl)*Xi[k];
    }
    return rtn_value;
}

double Interpolate2d9_X(const double * Xi, const double * xl,int dim){
    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += Ip2d9_X[k][dim](xl)*Xi[k];
    }
    return rtn_value;
}

double Interpolate2d9_XX(const double * Xi, const double * xl,int dim){
    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += Ip2d9_XX[k][dim](xl)*Xi[k];
    }
    return rtn_value;
}

double vInterpolate2d9(const double * Xi, int s,int f,const double * xl){
    /*
     * just interpolate Xi at shift s*k + f, k in [0,9)
     */
    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += Ip2d9[k](xl)*Xi[k*s + f];
    }
    return rtn_value;
}

double vInterpolate2d9_X(const double * Xi, int s,int f,const double * xl,int dim){
    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += Ip2d9_X[k][dim](xl)*Xi[k*s + f];
    }
    return rtn_value;
}

double vInterpolate2d9_XX(const double * Xi, int s,int f,const double * xl,int dim){
    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += Ip2d9_XX[k][dim](xl)*Xi[k*s + f];
    }
    return rtn_value;
}


double Interpolate2d9_Fdx(const double *Xi,const double *Yi, const double *xl,int dim)
{
    /*
     * jacobi_maxtrix is partial(x,y)/partial(xi,eta)
     * inverse_jacobi_matrix is partial(xi,eta)/partial(x,y)
     */

    if(NULL == Xi) Xi = xl2d9[0];

    double jacobi_matrix[4] = {0.};
    jacobi_matrix[0] = vInterpolate2d9_X(Xi,2,0,xl,0);
    jacobi_matrix[1] = vInterpolate2d9_X(Xi,2,0,xl,1);
    jacobi_matrix[2] = vInterpolate2d9_X(Xi,2,1,xl,0);
    jacobi_matrix[3] = vInterpolate2d9_X(Xi,2,1,xl,1);

    double detJ = jacobi_matrix[0]*jacobi_matrix[3] - jacobi_matrix[1]*jacobi_matrix[2];
    double inverse_jacobi_matrix[4] = {0.};
    inverse_jacobi_matrix[0] = jacobi_matrix[3]/detJ;
    inverse_jacobi_matrix[3] = jacobi_matrix[0]/detJ;
    inverse_jacobi_matrix[1] = -jacobi_matrix[2]/detJ;
    inverse_jacobi_matrix[2] = -jacobi_matrix[1]/detJ;

    /*
     * Ydx1 = partial Y/ partial xi
     * Ydx2 = partial Y/ partial eta
     *
     * partial Y/partial x = (Y_xi)(xi_x) + (Y_eta)(eta_x)
     */
    double Ydx1 = Interpolate2d9_X(Yi,xl,0);
    double Ydx2 = Interpolate2d9_X(Yi,xl,1);

    return Ydx1*inverse_jacobi_matrix[0+dim] + Ydx2*inverse_jacobi_matrix[2+dim];
}

double Interpolate2d9_Fd2x(const double *Xi,const double *Yi, const double *xl,int dim)
{
    double jacobi_matrix[4] = {0.};
    jacobi_matrix[0] = vInterpolate2d9_X(Xi,2,0,xl,0);
    jacobi_matrix[1] = vInterpolate2d9_X(Xi,2,0,xl,1);
    jacobi_matrix[2] = vInterpolate2d9_X(Xi,2,1,xl,0);
    jacobi_matrix[3] = vInterpolate2d9_X(Xi,2,1,xl,1);

    double detJ = jacobi_matrix[0]*jacobi_matrix[3] - jacobi_matrix[1]*jacobi_matrix[2];
    double inverse_jacobi_matrix[4] = {0.};
    inverse_jacobi_matrix[0] = jacobi_matrix[3]/detJ;
    inverse_jacobi_matrix[3] = jacobi_matrix[0]/detJ;
    inverse_jacobi_matrix[1] = -jacobi_matrix[2]/detJ;
    inverse_jacobi_matrix[2] = -jacobi_matrix[1]/detJ;

    double Ydx1x1 = Interpolate2d9_XX(Yi,xl,0);
    double Ydx1x2 = Interpolate2d9_XX(Yi,xl,1);
    double Ydx2x2 = Interpolate2d9_XX(Yi,xl,2);

    if(0==dim) //dxdx
    {
        return Ydx1x1*inverse_jacobi_matrix[0]*inverse_jacobi_matrix[0]
               + 2.0*Ydx1x2*inverse_jacobi_matrix[0]*inverse_jacobi_matrix[2]
               + Ydx2x2*inverse_jacobi_matrix[2]*inverse_jacobi_matrix[2];
    }
    else if(1==dim) //dxdy
    {
        return Ydx1x1*inverse_jacobi_matrix[0]*inverse_jacobi_matrix[1]
        + Ydx1x2*inverse_jacobi_matrix[2]*inverse_jacobi_matrix[1]
        + Ydx1x2*inverse_jacobi_matrix[0]*inverse_jacobi_matrix[3]
        + Ydx2x2*inverse_jacobi_matrix[2]*inverse_jacobi_matrix[3];
    }
    else if(2==dim) //dydy
    {
        return Ydx1x1*inverse_jacobi_matrix[1]*inverse_jacobi_matrix[1]
        + 2.0*Ydx1x2*inverse_jacobi_matrix[1]*inverse_jacobi_matrix[3]
        + Ydx2x2*inverse_jacobi_matrix[3]*inverse_jacobi_matrix[3];
    }
    else
    {
        fprintf(stderr,"invalid dimension in %s\n",__func__ );
        return 0.;
    }
}



double jacobi2d9(const double *Xi, const double *xl)
{
    double jacobi_matrix[4] = {0.};
    jacobi_matrix[0] = vInterpolate2d9_X(Xi,2,0,xl,0);
    jacobi_matrix[1] = vInterpolate2d9_X(Xi,2,0,xl,1);
    jacobi_matrix[2] = vInterpolate2d9_X(Xi,2,1,xl,0);
    jacobi_matrix[3] = vInterpolate2d9_X(Xi,2,1,xl,1);

    return jacobi_matrix[0]*jacobi_matrix[3] - jacobi_matrix[1]*jacobi_matrix[2];
}

double Int2d9(const double *Xi,const double *Yi)
{
    double rtn_value = 0.;
    if(NULL == Xi) Xi = xl2d9[0];

    for(int k=0;k<9;++k)
    {
        rtn_value += GIWS2d9[k]*Interpolate2d9(Yi,GIPS2d9[k])*jacobi2d9(Xi,GIPS2d9[k]);
    }
    return rtn_value;
}

double Int2d9Fdx(const double *Xi,const double *Yi, int dim)
{
    double local_field[18];
    for(int k=0;k<18;++k)
    {
        if(k%9 == dim)
            local_field[k] = Yi[k];
        else
            local_field[k] = Xi[k];
    }

    double rtn_value = 0.;
    for(int k=0;k<9;++k)
    {
        rtn_value += GIWS2d9[k]* jacobi2d9(local_field,GIPS2d9[k]);
    }
    return rtn_value;
}

double Laplace2d9(const double * Xi, const double *Yi, double gamma)
{
    /*
     * the intergate interval is [-gamma,gamma]x[-gamma,gamma]
     * when gamma = sqrt(3/4) = Oono-Puri nine point
     * gamma = sqrt(1.0/2) = Patra-Karttunen (recommended)
     * check at wiki/Nine-point stencil
     */

    if(gamma < 0.) gamma*= -1.;
    if(gamma > 1.0) fprintf(stderr,"[%s]:warning: gamma = %f (this value should <= 1)\n",__func__ ,gamma);

    const double * local_Xi = xl2d9[0];
    if(NULL != Xi) local_Xi = Xi;

    double rtn_value = 0.;
    double rtn_area = 0.;
    for(int k=0;k<9;++k)
    {
        double local_xl[2] = {GIPS2d9[k][0]*gamma,GIPS2d9[k][1]*gamma};
        double Dxx = Interpolate2d9_Fd2x(local_Xi,Yi,local_xl,0);
        double Dyy = Interpolate2d9_Fd2x(local_Xi,Yi,local_xl,2);
        rtn_value += (gamma*gamma)*GIWS2d9[k]*(Dxx + Dyy)* jacobi2d9(local_Xi,local_xl);
        rtn_area += (gamma*gamma)*GIWS2d9[k]*1.0*jacobi2d9(local_Xi,local_xl);
    }
    return rtn_value/rtn_area;
}

double Laplace2d9cyl(const double * Xi, const double *Yi, double gamma)
{
    /*
     * calculate/integrate Laplacian in cylindrical coordinate.
     * Laplacian = D_rr + D_zz + D_r/r
     * the terms related with theta is ignored.
     */

    if(gamma < 0.) gamma*= -1.;
    if(gamma > 1.0) fprintf(stderr,"[%s]:warning: gamma = %f (this value should <= 1)\n",__func__ ,gamma);

    if(NULL == Xi) Xi = xl2d9[0];

    double rtn_value_right = 0.;
    double rtn_value_left = 0.;
    double rtn_area_right = 0.;
    double rtn_area_left = 0.;

    // if (0,0) is in integrate section
    // Laplacian is discontinued
    // split the domain into [-gamma,0]x[-gamma,gamma[ and [0,gamma]x[-gamma,gamma]


    const int num_pts = 16;
    for(int k=0;k<num_pts;++k)
    {

        // the [0,gamma]x[-gamma,gamma] section
        double local_xl[2] = {(GIPS2d16[k][0] + 1.0)*gamma/2.,GIPS2d16[k][1]*gamma};
        {
            double Dxx = Interpolate2d9_Fd2x(Xi,Yi,local_xl,0);
            double Dyy = Interpolate2d9_Fd2x(Xi,Yi,local_xl,2);
            double Dx  = Interpolate2d9_Fdx(Xi,Yi,local_xl,0);
            double R  = vInterpolate2d9(Xi,2,0,local_xl);

            if(R < 0.)
            {
                R  *= -1.;
                Dx *= -1.;
            }

            rtn_value_right  += (gamma/2.*gamma)*GIWS2d16[k]*((Dxx + Dyy)*R + Dx)* jacobi2d9(Xi,local_xl);
            rtn_area_right   += (gamma/2.*gamma)*GIWS2d16[k]*1.0*R*jacobi2d9(Xi,local_xl);
        }

        // the [-gamma,0]x[-gamma,gamma] section
        local_xl[0] = (GIPS2d16[k][0] - 1.0)*gamma/2.0;
        {
            double Dxx = Interpolate2d9_Fd2x(Xi,Yi,local_xl,0);
            double Dyy = Interpolate2d9_Fd2x(Xi,Yi,local_xl,2);
            double Dx  = Interpolate2d9_Fdx(Xi,Yi,local_xl,0);
            double R  = vInterpolate2d9(Xi,2,0,local_xl);

            if(R < 0.)
            {
                R *= -1.;
                Dx *= -1.;
            }

            rtn_value_left  += (gamma/2.*gamma)*GIWS2d16[k]*((Dxx + Dyy)*R + Dx)* jacobi2d9(Xi,local_xl);
            rtn_area_left   += (gamma/2.*gamma)*GIWS2d16[k]*1.0*R*jacobi2d9(Xi,local_xl);
        }

    }

    return (rtn_value_left + rtn_value_right)/(rtn_area_left + rtn_area_right);
}

#define GIP2d16_a __builtin_sqrt(3./7. - 2./7.*__builtin_sqrt(6./5.))
#define GIP2d16_b __builtin_sqrt(3./7. + 2./7.*__builtin_sqrt(6./5.))
#define GIP2d16_c (18. + __builtin_sqrt(30.))/36.
#define GIP2d16_d (18. - __builtin_sqrt(30.))/36.

#define _GIP2d16_array_r(func_d,x,y1,y2) \
    func_d(x,-y2), func_d(x,-y1),        \
    func_d(x, y1), func_d(x, y2)

#define _GIP2d16_array_l(func_d1,func_d2,x1,x2,y1,y2) \
    func_d1(func_d2,-x2,y1,y2),                   \
    func_d1(func_d2,-x1,y1,y2),                    \
    func_d1(func_d2, x1,y1,y2),                   \
    func_d1(func_d2, x2,y1,y2)

#define _GIP2d16_pts(x,y) {x,y}
#define _GIP2d16_wts(x,y) __builtin_fabs((x)*(y))

const double GIPS2d16[16][2] =
        {
                _GIP2d16_array_l(_GIP2d16_array_r,_GIP2d16_pts,GIP2d16_a,GIP2d16_b,GIP2d16_a,GIP2d16_b)
        };

const double GIWS2d16[16] =
        {
                _GIP2d16_array_l(_GIP2d16_array_r,_GIP2d16_wts,GIP2d16_c,GIP2d16_d,GIP2d16_c,GIP2d16_d)
        };

#undef _GIP2d16_array_r
#undef _GIP2d16_array_l
#undef _GIP2d16_pts
#undef _GIP2d16_wts