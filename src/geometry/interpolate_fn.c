//const double BETA = 0.57735026918962576451;
#include "interpolate_fn.h"

const double GIPS2d[NGI2d][2] = {
        {-BETA, -BETA},
        {-BETA, BETA},
        {BETA,  -BETA},
        {BETA,  BETA}};
const double GIPS3d[NGI3d][3] = {
        {-BETA, -BETA, -BETA},
        {-BETA, -BETA, BETA},
        {-BETA, BETA,  -BETA},
        {-BETA, BETA,  BETA},
        {BETA,  -BETA, -BETA},
        {BETA,  -BETA, BETA},
        {BETA,  BETA,  -BETA},
        {BETA,  BETA,  BETA}};
const double GIWS2d[NGI2d] = {1.0, 1.0, 1.0, 1.0};
const double GIWS3d[NGI3d] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// functions used in 4 points interpolate
double Ip2d0(const double *_x)
{
    return 0.25 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip2d0_X1(const double *_x)
{
    return 0.25 * (1 + _x[1]);
}

double Ip2d0_X2(const double *_x)
{
    return 0.25 * (1 + _x[0]);
}

double Ip2d1(const double *_x)
{
    return 0.25 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip2d1_X1(const double *_x)
{
    return -0.25 * (1 + _x[1]);
}

double Ip2d1_X2(const double *_x)
{
    // correct by huacheng
    return 0.25 * (1 - _x[0]);
}

double Ip2d2(const double *_x)
{
    return 0.25 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip2d2_X1(const double *_x)
{
    return -0.25 * (1 - _x[1]);
}

double Ip2d2_X2(const double *_x)
{
    return -0.25 * (1 - _x[0]);
}


// functions used in 4 points interpolate
double Ip2d3(const double *_x)
{
    return 0.25 * (1 + _x[0]) * (1 - _x[1]);
}

double Ip2d3_X1(const double *_x)
{
    return 0.25 * (1 - _x[1]);
}

double Ip2d3_X2(const double *_x)
{
    return -0.25 * (1 + _x[0]);
}
// the functions used in 4 points for a boundary
// see the array IpB and IpB_X

// functions used in 8 points interpolate
// see the function array IpV and IpV_B
double Ip3d0(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d1(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d4(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d5(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d6(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d7(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d0_X1(const double *_x)
{
    return 0.125 * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d1_X1(const double *_x)
{
    return -0.125 * (1 + _x[1]) * (1 - _x[2]);
}

double Ip3d2_X1(const double *_x)
{
    return -0.125 * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d3_X1(const double *_x)
{
    return 0.125 * (1 - _x[1]) * (1 - _x[2]);
}

double Ip3d4_X1(const double *_x)
{
    return 0.125 * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d5_X1(const double *_x)
{
    return -0.125 * (1 + _x[1]) * (1 + _x[2]);
}

double Ip3d6_X1(const double *_x)
{
    return -0.125 * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d7_X1(const double *_x)
{
    return 0.125 * (1 - _x[1]) * (1 + _x[2]);
}

double Ip3d0_X2(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[2]);
}

double Ip3d1_X2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[2]);
}

double Ip3d2_X2(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 - _x[2]);
}

double Ip3d3_X2(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 - _x[2]);
}

double Ip3d4_X2(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[2]);
}

double Ip3d5_X2(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[2]);
}

double Ip3d6_X2(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 + _x[2]);
}

double Ip3d7_X2(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 + _x[2]);
}

double Ip3d0_X3(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip3d1_X3(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip3d2_X3(const double *_x)
{
    return -0.125 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip3d3_X3(const double *_x)
{
    return -0.125 * (1 + _x[0]) * (1 - _x[1]);
}

double Ip3d4_X3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 + _x[1]);
}

double Ip3d5_X3(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 + _x[1]);
}

double Ip3d6_X3(const double *_x)
{
    return 0.125 * (1 - _x[0]) * (1 - _x[1]);
}

double Ip3d7_X3(const double *_x)
{
    return 0.125 * (1 + _x[0]) * (1 - _x[1]);
}

Interpolation2d9_t IpB[NIpB] = {Ip2d2, Ip2d3, Ip2d0, Ip2d1};
Interpolation2d9_t IpB_X[NIpB][2] = {
        {Ip2d2_X1, Ip2d2_X2},
        {Ip2d3_X1, Ip2d3_X2},
        {Ip2d0_X1, Ip2d0_X2},
        {Ip2d1_X1, Ip2d1_X2},
};

Interpolation2d9_t IpB_XT[2][NIpB] = {
        {Ip2d2_X1,Ip2d3_X1,Ip2d0_X1,Ip2d1_X1},
        {Ip2d2_X2,Ip2d3_X2,Ip2d0_X2,Ip2d1_X2}
};

Interpolation2d9_t IpV[NIpV] = {Ip3d0, Ip3d1, Ip3d2, Ip3d3, Ip3d4, Ip3d5, Ip3d6, Ip3d7};
Interpolation2d9_t IpV_X[NIpV][3] = {
        {Ip3d0_X1, Ip3d0_X2, Ip3d0_X3},
        {Ip3d1_X1, Ip3d1_X2, Ip3d1_X3},
        {Ip3d2_X1, Ip3d2_X2, Ip3d2_X3},
        {Ip3d3_X1, Ip3d3_X2, Ip3d3_X3},
        {Ip3d4_X1, Ip3d4_X2, Ip3d4_X3},
        {Ip3d5_X1, Ip3d5_X2, Ip3d5_X3},
        {Ip3d6_X1, Ip3d6_X2, Ip3d6_X3},
        {Ip3d7_X1, Ip3d7_X2, Ip3d7_X3},
};


