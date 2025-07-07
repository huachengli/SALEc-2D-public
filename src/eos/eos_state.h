//
// Created by huacheng on 1/11/22.
//

#ifndef EOSTOOL_EOS_STATE_H
#define EOSTOOL_EOS_STATE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "InputParser.h"
#include "linear.h"

struct ANEOSTable
{
    /*
     * structured table, interpolate the value linearly
     * between data points
     */
    int nTem;
    int nDen;
    double Pnorm;
    double Dnorm;
    double Tnorm;
    double *yDen;
    double *xTem;
    double *** Data;
    // data is (density*energy,pressure,csound)
    // not same with iSALE2d
};

#define EOSTEM -2
#define EOSDEN -1
#define EOSENG 0
#define EOSPRE 1
#define EOSCSD 2
#define CAP_MINPRESSURE 1

#define WITHOUTCORRECT 0
#define CONSTSHIFT 1

struct StateReference
{
    /*
     * The parameters shared in one materials
     * which can be used in EOS or strength model.
     */
    int MaterialId;
    double VaporDen;
    double MeltDen;
    double GaCoef; // Constant to convert bulk modulus to shear modulus

    // parameters for simple ROCK strength model
    double Yint0;
    double Yfricint;
    double Ylimint;
    double Ydam0;
    double Yfricdam;
    double Ylimdam;

    char yfunc[MaxStrLen];
    char tfunc[MaxStrLen];
    double Yten0;

    double IvanA;
    double IvanB;
    double IvanC;
    double Pbd;
    double Pbp;


    char dSDf[MaxStrLen];

    double minPint;
    double minPdam;

    // parameters for Johnson-Cook strength model
    double JCA;
    double JCB;
    double JCN;
    double JCC;
    double JCM;
    double JCTREF;
    double JCMaxPlastic;
    double JCMinP;

    // soft parameters
    // simon (a,c,T) is the solidus
    double Asimon;
    double Csimon;
    double Tmelt0;
    double Tfrac;
    double Tdelta;
    double Viscosity;

    // liquids
    int PolyLiq; // flag for whether use poly liquids
    double Tlidc1;
    double Tlidc2;
    double Tlidc3;

    // solids
    int PolySol;
    double Tsolc1;
    double Tsolc2;
    double Tsolc3;

    // Acoustic Fluidization parameters
    double Toff;
    double Cvib;
    double VibMax;
    double Tdamp;
    double Pvlim;
    double Acvis;
    double GammaEta;
    double GammaBeta;

    // porosity model, Wünnemann et al. (2006)
    // input parameters
    char porosityfunc[MaxStrLen];
    double alpha0;
    double alphax;
    double epse0;
    double chi;
    double kappa;
    double avel; // activate speed
    // derived parameters
    double epsec;
    double epsex;
    double alphamax;
    double alphamax2;
    double alphae;

    // other constants
    double gamfrac; // gamma/(csd0**2)
    double csd0;

    // additional information for temperature/gravity profile
    double thermexp;
    double c_heat;
    double den0;

    // args for melting attach in thermal soft
    double melt_viscosity;

    uintptr_t shear_strength_func;
    uintptr_t tensile_strength_func;
    uintptr_t damage_increment_func;
    uintptr_t porosity_increment_func;
};

void AllocateANEOS(struct ANEOSTable *_t);
void UnAllocateANEOS(struct ANEOSTable *_t);
void LoadANEOS(struct ANEOSTable *_t, FILE *fp);

void ANEOSInitStateRef(struct ANEOSTable *_t, struct StateReference *_s);
double ANEOSInterpolateTD(struct ANEOSTable *_t,double tTem, double tDen, int DataId);
double ANEOSInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId);
void ANEOSWrite(struct ANEOSTable * _t,const char fname[], const char comment[]);

struct AirEOSTable
{
    int nEng;
    int nDen;
    double * xEng;
    double * yDen;
    double *** Data;
    // (All the variables is log10. form)
    // Data = log.10 (Pressure,Temperature,Speed of sound)
};

int BilinearInverse(const double xi[],const double yi[], double xl[]);
void AllocateAirEOS(struct AirEOSTable * _air);
void UnAllocateAirEOS(struct AirEOSTable * _air);
void LoadAirEOS(struct AirEOSTable * _t, FILE *fp);

struct TillotsonTable
{
    double TLRho0;
    double TLCv;
    double TLA;
    double TLB;
    double TLE0;
    double TLa;
    double TLb;
    double TLAlpha;
    double TLBeta;
    double TLEiv;
    double TLEcv;
    double TLTref;

    double * xDen;
    double * yEng;
    double StepDen;
    int nDen;
};

double TillPres(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs);
double TillPres_split(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs);

double TillTemp(const struct TillotsonTable * _t,double EngIn, double RhoIn);
double TillCold(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs);
double TillHot(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs);

void TillColdEnergy(struct TillotsonTable * _t);
double LoadTillEOS(struct TillotsonTable * _t, FILE * fp);
void TillInitStateRef(struct TillotsonTable *_t, struct StateReference *_s);
double TillEOSInterpolateTP(struct TillotsonTable *_t, double tTem, double tPre, int DataId);
double TillEOSInterpolateTD(struct TillotsonTable *_t, double tTem, double tDen, int DataId);
void UnAllocateTillEOS(struct TillotsonTable *_t);


// some addition function
double Wind(double x, double limitL, double limitR);
void Over(FILE * fp,int n);
double Max(double a, double b);
double Min(double a, double b);
double Sqrt_trunc(double x);
void ANEOSLowDenCorrect(struct ANEOSTable * _t, double plimit);

typedef struct StateReference state_reference;
typedef struct ANEOSTable aneos_table;
typedef struct TillotsonTable tillotson_table;
typedef enum EoSType
{
    ANEOS,
    Tillotson,
    ANEOSSieDen,
    Unknown
} eos_type;

typedef struct EoSTable
{
    state_reference * ref;
    void * table;
    eos_type type;
} eos_table;

double calculate_state(eos_table * _e,double den,double eng,double * tem,double * pre,double * csd);
double calculate_state_correct(eos_table * _e,double den,double eng,double plim,double * tem,double * pre,double * csd, double *pty);
void aneos_calculate_state(aneos_table * _t,double den,double eng,double * tem,double * pre,double * csd);
void tillotson_calculate_state(tillotson_table * _t,double den,double eng,double * tem,double * pre,double * csd);
void load_eos_table(eos_table * _e, FILE * fp);
double interpolate_eos_tp(eos_table * _e,double tem,double pre, int DataId);
double interpolate_eos_td(eos_table * _e,double tem,double pre, int DataId);

// some function add for aneos_sieden
//void AllocateANEOSSieDen(struct ANEOSTable *_t);
//void UnAllocateANEOSSieDen(struct ANEOSTable *_t);
void LoadANEOSSieDen(struct ANEOSTable *_t, FILE *fp);

void ANEOSSieDenInitStateRef(struct ANEOSTable *_t, struct StateReference *_s);
double ANEOSSieDenInterpolateTD(struct ANEOSTable *_t,double tTem, double tDen, int DataId);
double ANEOSSieDenInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId);
double aneossieden_calculate_state(aneos_table * _t,double den,double eng,double * tem,double * pre,double * csd);
#endif //EOSTOOL_EOS_STATE_H
