//
// Created by huacheng on 4/25/23.
//
/*
 * target_profile.c/h is adopt from SALEc/Planet.c/h
 */



#ifndef SALE_REBUILD_TARGET_PROFILE_H
#define SALE_REBUILD_TARGET_PROFILE_H

#include "sale2d.h"

typedef struct TargetProfile
{
    double * pre;
    double * rad;
    double * den;
    double * tem;
    double * gra;
    double * dam;
    double * eng;
    int * mat;
    int n;
    double dz;
    target_type type;
    target_info * info;
} target_profile;

typedef target_info target_definition;


int AllocateTargetProfile(target_profile * _prof);
void UnallocateTargetProfile(target_profile * _prof);
void WriteTargetProfile(FILE * fp,target_profile * _prof);
void WriteTargetProfile_All(FILE * fp,target_profile * _prof, state_reference * _s, target_definition * _pdef);
void ConsistentLith(target_definition * _tdef);
void SetTargetProfile_Tem(target_profile * _prof, state_reference * _s, target_definition * _tdef);
void SetTargetProfile_Mat(target_profile * _prof, double * depth, int * IdTarget, int NumTarget);
void SetTargetProfile_Gra(target_profile * _prof, double gz);
void CorrectDensity(target_profile * _prof, eos_table * _t, state_reference * _s);
void CorrectionDensity1(target_profile * _prof, eos_table * _t, state_reference * _s);
void CorrectionDensity2(target_profile * _prof, eos_table * _t, state_reference * _s);
void IntegrateGravity(target_profile * _prof);
void LoadTargetProfile(InputFile * ifp, target_definition * _tdef);
// function set the initial profile
void init_material_distribution(mesh2d_info * _minfo, sale2d_var * _sale);
void init_material_plane_distribution(mesh2d_info * _minfo, sale2d_var * _sale);
void init_material_sphere_distribution(mesh2d_info * _minfo, sale2d_var * _sale);
void build_target_profile(target_profile * _tprof, mesh2d_info * _minfo, sale2d_var * _sale);

#endif //SALE_REBUILD_TARGET_PROFILE_H
