//
// Created by huacheng on 3/15/23.
//

#ifndef SALE_REBUILD_STRENGTH_H
#define SALE_REBUILD_STRENGTH_H
#include "sale2d.h"
#include "mpi_vec.h"

typedef double (*strength_func)(double*,const state_reference *);
void update_e_stress(sale2d_var * _sale, node_type _nt, double dt);
void select_strength_model(sale2d_var * _sale);
void update_strength_condition(sale2d_var * _sale);

// strenth function
double JohnsonCook2(double *_paras, const state_reference *_s);
double JohnsonCook1(double * _paras, const state_reference *_s);
double SimpleRock(double * _paras, const state_reference * _s);
double IvanovShearDam(double *_paras, const state_reference *_s);
double CollinsShearDam(double * _paras, const state_reference * _s);

// some expression from SALEc
double Ylundborg(double p,double y0,double fric,double ylim);
double Ydrucker(double p,double y0,double fric,double ylim);
double LowdensitySoft(double Density, double MeltDensity);
double JNCKSoft(double tem,double pre, const state_reference *_s);
double SimonMelt0(double pre, const state_reference *_s);
double SimonMelt(double pre, const state_reference *_s);
double SolidusTem(double pre, const state_reference * _s);
double PolyLiquids(double pre, const state_reference * _s);
double OhnakaSoft(double tem, const state_reference *_s, double Tmelt);
double LowdensitySoft(double Density, double MeltDensity);
double VonMises1(double * _paras, const state_reference * _s);
double VonMises2(double * _paras, const state_reference * _s);
void block_vibration(sale2d_var * _sale, node_type _nt,double dt,double time);
double PorosityCsdRatio(double alpha, double alpha0, double chi);

typedef struct StrFuncPair
{
    strength_func f;
    char name[256];
} str_func_pair;

extern str_func_pair str_func_list[64];
#endif //SALE_REBUILD_STRENGTH_H
