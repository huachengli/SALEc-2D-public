//
// Created by huacheng on 4/25/23.
//

#include "target_profile.h"

#include <strength.h>

#define is_null(_ptr) (_ptr == NULL)

int AllocateTargetProfile(target_profile * _prof)
{
    if(_prof->n > 0)
    {
        _prof->pre = (double *) malloc(sizeof(double)*_prof->n);
        _prof->rad = (double *) malloc(sizeof(double)*_prof->n);
        _prof->den = (double *) malloc(sizeof(double)*_prof->n);
        _prof->tem = (double *) malloc(sizeof(double)*_prof->n);
        _prof->gra = (double *) malloc(sizeof(double)*_prof->n);
        _prof->dam = (double *) malloc(sizeof(double)*_prof->n);
        _prof->eng = (double *) malloc(sizeof(double)*_prof->n);
        _prof->mat = (int *) malloc(sizeof(int)*_prof->n);

        return is_null(_prof->pre) ? 0 : _prof->n;

    } else
    {
        return -1;
    }

}
void UnallocateTargetProfile(target_profile * _prof)
{
    free(_prof->pre);
    free(_prof->rad);
    free(_prof->den);
    free(_prof->tem);
    free(_prof->mat);
    free(_prof->gra);
    free(_prof->dam);
    free(_prof->eng);
}

void WriteTargetProfile(FILE * fp,target_profile * _prof)
{
    fputs("      Depth, Matid,    Density,    Gravity, Pressure,Temperature\n",fp);
    for(int k= _prof->n-1; k>= 0; --k)
    {
        fprintf(fp,"%10.5f,   %3d, %10.5f, %10.5f,%10.5e,%10.5e\n",_prof->rad[k], _prof->mat[k],_prof->den[k],
                _prof->gra[k] ,_prof->pre[k],_prof->tem[k]);
    }
}

void WriteTargetProfile_All(FILE * fp,target_profile * _prof, state_reference * _s, target_definition * _pdef)
{
    fputs("      Depth, Matid,    Density,    Gravity, Pressure,Temperature,solidus\n",fp);
    for(int k= _prof->n-1; k>= 0; --k)
    {
        state_reference * local_sref = _s + _prof->mat[k];
        double melt_tem = SolidusTem(_prof->pre[k],local_sref);
        fprintf(fp,"%10.5f,   %3d, %10.5f, %10.5f,%10.5e,%10.5e, %10.5e\n",_prof->rad[k], _prof->mat[k],_prof->den[k],
                _prof->gra[k] ,_prof->pre[k],_prof->tem[k],
                melt_tem);
    }
}

void ConsistentLith(target_definition * _pdef)
{
    if(_pdef->t_lith >= 0 && _pdef->d_lith < 0)
    {
        // calculated d_lith automatically
        double Delta = 1.0 - 2.0*(_pdef->t_lith - _pdef->surf_tem)/(_pdef->r_planet*_pdef->surf_dtdz);
        Delta = Max(0.0,Delta);
        _pdef->d_lith = _pdef->r_planet * (1.0 - sqrt(Delta));
    } else if(_pdef->t_lith < 0 && _pdef->d_lith >= 0)
    {
        _pdef->t_lith = _pdef->surf_tem + _pdef->surf_dtdz*_pdef->d_lith*(1.0 - 0.5*_pdef->d_lith/_pdef->r_planet);
    } else
    {
        fprintf(stdout,"%s:cannot set d_lith and t_lith at same time!\n",__func__);
        exit(0);
    }
}


void SetTargetProfile_Tem(target_profile * _prof, state_reference * _s, target_definition * _pdef)
{
    // first set the surf_pressure as initial pressure
    for(int k=0;k<_prof->n;++k)
    {
        _prof->pre[k] = _pdef->surf_pre;
    }

    // check the target type, set the density
    if(_pdef->layer_type == sphere_target)
    {
        for(int k=0;k<_prof->n;++k)
        {
            _prof->den[k] = (_s + _prof->mat[k])->den0;
            _prof->rad[k] = _pdef->r_planet*k/(_prof->n - 1);
        }
    } else if(_pdef->layer_type == plane_target)
    {
        for(int k=0;k<_prof->n;++k)
        {
            _prof->den[k] = (_s + _prof->mat[k])->den0;
            _prof->rad[k] = _pdef->max_depth*k/(_prof->n - 1);
        }
    } else
    {
        fprintf(stderr,"Invalid target layer_type in TargetProfile!\n");
        exit(0);
    }

    // set the initial temperature
    switch (_pdef->tem_type)
    {
        case const_profile:
            for(int k=0;k<_prof->n;++k)
            {
                _prof->tem[k] = _pdef->surf_tem;
            }
            break;
        case conductive_profile:
            for(int k=_prof->n-1;k>=0;--k)
            {
                double k_depth = _prof->rad[_prof->n-1] - _prof->rad[k];
                if(_prof->n== (k+1))
                {
                    _prof->tem[k] = _pdef->surf_tem;
                    continue;
                }

                double dtdz_cond = _pdef->surf_dtdz * (1.0 + k_depth/_pdef->r_planet);
                _prof->tem[k] = _prof->tem[k+1] + (dtdz_cond)*(_prof->rad[k+1] - _prof->rad[k]);
            }
            break;
        case convective_profile:
        case condconv_profile:
        case lunar_profile:
            for(int k=_prof->n-1;k>=0;--k)
            {
                double k_depth = _prof->rad[_prof->n-1] - _prof->rad[k];
                if(_prof->n== (k+1))
                {
                    _prof->tem[k] = _pdef->surf_tem;
                    continue;
                }

                double dtdz_cond = _pdef->surf_dtdz * (1.0 + k_depth/_pdef->r_planet);
                double conv_coef = _s[_prof->mat[k]].thermexp * _pdef->surf_gravity / _s[_prof->mat[k]].c_heat;
                double dtdz_conv = _pdef->t_lith * conv_coef * exp(conv_coef * (k_depth - _pdef->d_lith));
                // double dtdz_conv = _prof->tem[k+1] * conv_coef;
                double conduction = 1.0 - (_prof->tem[k+1] - 0.9*_pdef->t_lith)/(0.11*_pdef->t_lith);
                conduction = Wind(conduction,0.0,1.0);
                _prof->tem[k] = _prof->tem[k+1] + (dtdz_cond*conduction + dtdz_conv*(1.0 - conduction))*(_prof->rad[k+1] - _prof->rad[k]);

                state_reference * local_sref = _s + _prof->mat[k];
                // double melt_tem = local_sref->Tmelt0 * pow(local_sref->Asimon*_prof->pre[k] + 1.0,local_sref->Csimon);
                // replace melt_tem as solidus; if polysol is defined
                double melt_tem = SolidusTem(_prof->pre[k],local_sref);
                // if(_prof->tem[k] > melt_tem) _prof->tem[k] = melt_tem;
            }
            break;
        default:
            fprintf(stderr,"Invalid temperature prof type in SetTem\n");
    }
}

void SetTargetProfile_Mat(target_profile * _prof, double * depth, int * IdTarget, int NumTarget)
{
    int cur_cell = _prof->n - 1;
    for(int iTar=0;iTar<NumTarget;iTar++)
    {
        int TarIntervals = (int)floor(depth[iTar]/_prof->dz);
        for(int j=0;j<TarIntervals && cur_cell > 0;++j)
        {
            _prof->mat[cur_cell--] = IdTarget[iTar];
        }
    }

    for(;cur_cell>=0;--cur_cell)
    {
        _prof->mat[cur_cell] = IdTarget[NumTarget-1];
    }
}

void SetTargetProfile_Gra(target_profile * _prof, double gz)
{
    for(int k=0;k<_prof->n;++k)
    {
        _prof->gra[k] = gz;
        _prof->dam[k] = 0.;
    }
}


/*
 * the difference between CorrectDensity and Cap*Profile is the variables set
 * CorrectDensity recalculate the density according to (temperature,pressure)
 * Cap*Profile initialize the energy profile at after loop (IntergateGravity,CorrectDensity)
 */

void CorrectDensity(target_profile * _prof, eos_table * _t, state_reference * _s)
{
    if(_prof->info->tem_type != convective_profile)
    {
        CorrectionDensity1(_prof,_t,_s);
    }
    else
    {
        CorrectionDensity2(_prof,_t,_s);
    }
}


void CorrectionDensity1(target_profile * _prof, eos_table * _t, state_reference * _s)
{
    // not update the temperature, only cap
    for(int k=0;k<_prof->n;++k)
    {
        // cap the temperature
        state_reference * local_sref = _s + _prof->mat[k];
        // double melt_tem = local_sref->Tmelt0 * pow(local_sref->Asimon*_prof->pre[k] + 1.0,local_sref->Csimon);
        // replace melt_tem as solidus; if polysol is defined
        double melt_tem = SolidusTem(_prof->pre[k],local_sref);
        double local_tem = Min(_prof->tem[k],melt_tem);
        double local_pty = local_sref->alpha0;
        _prof->tem[k] = local_tem;
        eos_table * local_eos = _t + _prof->mat[k];
        _prof->eng[k] = interpolate_eos_tp(local_eos,local_tem,_prof->pre[k]*local_pty,EOSENG)/local_pty;
        _prof->den[k] = interpolate_eos_tp(local_eos,local_tem,_prof->pre[k]*local_pty,EOSDEN)/local_pty;
    }
}

void CorrectionDensity2(target_profile * _prof, eos_table * _t, state_reference * _s)
{
    // conduction lithosphere + convection mantle
    for(int k=_prof->n -1 ;k >= 0; k--)
    {
        // fix the surface temperature; just skip temperature updating of [n-1] index
        double k_depth = _prof->rad[_prof->n-1] - _prof->rad[k];
        target_info * _pdef = _prof->info;
        // recalculate temperature according to the gravity/thermexp, using local gravity not the surf_gravity
        double dtdz_cond = _pdef->surf_dtdz * (1.0 + k_depth/_pdef->r_planet);
        _prof->tem[k] = _pdef->surf_tem;
        if( (_prof->n - 1) != k)
        {
            double conv_coef = _s[_prof->mat[k]].thermexp * fabs(_prof->gra[k]) / _s[_prof->mat[k]].c_heat;
            double dtdz_conv = _prof->tem[k+1] * conv_coef;
            // double dtdz_conv = _pdef->t_lith * conv_coef * exp(conv_coef * (k_depth - _pdef->d_lith));
            // double conduction = 1.0 - (_prof->tem[k+1] - 0.9*_pdef->t_lith)/(0.10*_pdef->t_lith);
            double conduction = 1.0 - (k_depth - 0.9*_pdef->d_lith)/(0.1*_pdef->d_lith);
            // double conduction = 1.0 - (_prof->tem[k+1] - 0.9*_pdef->t_lith)/(0.11*_pdef->t_lith);
            conduction = Wind(conduction,0.0,1.0);

            if(_pdef->d_mantle >= _pdef->d_lith && k_depth < _pdef->d_mantle)
            {
                conduction = 1.0;
            }
            _prof->tem[k] = _prof->tem[k+1] + (dtdz_cond*conduction + dtdz_conv*(1.0 - conduction))*(_prof->rad[k+1] - _prof->rad[k]);

        }

        // cap the temperature greater than melt temperature
        state_reference * local_sref = _s + _prof->mat[k];
        // double melt_tem = local_sref->Tmelt0 * pow(local_sref->Asimon*_prof->pre[k] + 1.0,local_sref->Csimon);
        // replace melt_tem as solidus; if polysol is defined
        double melt_tem = SolidusTem(_prof->pre[k],local_sref);
        double local_tem = _prof->tem[k];

        if(_pdef->d_mantle >= _pdef->d_lith && k_depth > _pdef->d_mantle)
        {
            local_tem = _prof->tem[k];
        }
        else
        {
            local_tem = Min(_prof->tem[k],melt_tem);
        }


         _prof->tem[k] = local_tem;

        double local_pty = local_sref->alpha0;
        eos_table * local_eos = _t + _prof->mat[k];
        // refresh the energy/density
        _prof->eng[k] = interpolate_eos_tp(local_eos,local_tem,_prof->pre[k]*local_pty,EOSENG)/local_pty;
        _prof->den[k] = interpolate_eos_tp(local_eos,local_tem,_prof->pre[k]*local_pty,EOSDEN)/local_pty;
    }
}



void IntegrateGravity(target_profile * _prof)
{
    // just used for planet
    if(_prof->type == sphere_target)
    {
        _prof->gra[0] = 0.0;
        for(int k=1;k<_prof->n;++k)
        {
            _prof->gra[k] = _prof->gra[k-1] - 4.0*M_PI*6.674e-11*(_prof->den[k]+_prof->den[k-1])
                                              *(_prof->rad[k] + _prof->rad[k-1])*(_prof->rad[k] + _prof->rad[k-1])/8.0
                                              *(_prof->rad[k] - _prof->rad[k-1]);
        }
        for(int k=1;k<_prof->n;++k)
        {
            _prof->gra[k] = _prof->gra[k]/(_prof->rad[k]*_prof->rad[k]);
        }

        // _prof->pre[_prof->n-1] = 0.0;
        /*
         * _prof->pre[_prof->n-1] (surface pressure has been set by surf_pre)
         * do not overwrite it
         */
        for(int k=_prof->n-2;k>=0;--k)
        {
            _prof->pre[k] = _prof->pre[k+1] + (_prof->den[k] + _prof->den[k+1]) * (_prof->gra[k] + _prof->gra[k+1])
                                              * (_prof->rad[k] - _prof->rad[k+1])/4.0;
        }
    } else if(_prof->type == plane_target)
    {
        // just integrate from surface to the bottom
        // the density related to the pressure/energy is not considered
        // the correct density function will take this into calculations
        // _prof->pre[_prof->n-1] = 0.0;
        /*
         * _prof->pre[_prof->n-1] (surface pressure has been set by surf_pre)
         * do not overwrite it
         */

        for(int k=_prof->n-2;k>=0;--k)
        {
            // double drhogra = 0.5*(_prof->den[k+1]*_prof->gra[k] + _prof->den[k]*_prof->gra[k+1]) + (_prof->den[k] - _prof->den[k+1])*(_prof->gra[k] - _prof->gra[k+1])/3.0;
            double drhogra = 0.25*(_prof->den[k+1] + _prof->den[k])*(_prof->gra[k] + _prof->gra[k+1]);
            _prof->pre[k] = _prof->pre[k+1] + drhogra * (_prof->rad[k] - _prof->rad[k+1]);
        }
    } else
    {
        fprintf(stderr,"Invalid target type in %s!\n",__func__);
        exit(0);
    }
}

void load_target_info(InputFile * ifp, target_info * _tinfo)
{
    const int n_profile = 3001;
    _tinfo->n_profile = n_profile;
    _tinfo->n_layer = GetValueI(ifp,"target.number","1");

    char TargetTypeOpt[MaxStrLen];
    GetValueS(ifp,"target.type",TargetTypeOpt,"plane");
    if(strcasecmp("plane",TargetTypeOpt) == 0)
    {
        _tinfo->layer_type = plane_target;
    } else if(strcasecmp("sphere",TargetTypeOpt) == 0)
    {
        _tinfo->layer_type = sphere_target;
        // load r_planet info
        _tinfo->r_planet = GetValueD(ifp,"target.r_planet","1738.0e3");
        _tinfo->center[Y] = -_tinfo->r_planet;
        _tinfo->center[X] =  0.0;
    } else
    {
        _tinfo->layer_type = unknown_target;
    }

    // load name of materials,
    char material_name[MAXMAT][80];
    int nmat = GetValueI(ifp,"material.nm","3");
    for(int matid=0;matid<nmat;++matid)
    {
        GetValueSk(ifp,"material.name",material_name[matid],matid,"vacuum_");
    }

    double layer_depth[MAXMAT] = {0.};

    for(int k_layer=0;k_layer<_tinfo->n_layer;++k_layer)
    {
        _tinfo->depth[k_layer] = GetValueDk(ifp,"target.depth",k_layer,"1.");
        layer_depth[k_layer] = _tinfo->depth[k_layer];
        if(k_layer >= 1)
        {
            layer_depth[k_layer] += layer_depth[k_layer - 1];
        }

        // check the material index of layer
        char tmp_name[80];
        GetValueSk(ifp,"target.material",tmp_name,k_layer,"unknown");
        _tinfo->material[k_layer] = -1;
        for(int matid=0;matid<nmat;++matid)
        {
            if(0==strcasecmp(tmp_name,material_name[matid]))
            {
                _tinfo->material[k_layer] = matid;
                break;
            }
        }

        if(material_name[k_layer] < 0)
        {
            fprintf(stderr,"cannot get material id of %s in target layer\n",tmp_name);
            exit(0);
        }
    }

    _tinfo->max_depth = GetValueD(ifp,"target.max_depth","1.");
    if( _tinfo->layer_type == sphere_target)
    {
        _tinfo->max_depth = _tinfo->r_planet;
    }

    // set TEM_PROF_TYPE
    // PLANET/CONSTANT

    char TemProfOpt[MaxStrLen];
    GetValueSk(ifp,"target.temperature",TemProfOpt,0,"const");

    _tinfo->surf_pre = GetValueD(ifp,"target.surf_pressure","1.0e-5");

    if(strcasecmp("const",TemProfOpt)==0)
    {
        _tinfo->tem_type = const_profile;
        _tinfo->surf_tem = GetValueD(ifp,"target.surf_temperature","-1.0");
        if(_tinfo->surf_tem < 0.)
        {
            _tinfo->surf_tem = GetValueDk(ifp,"target.temperature",1,"293.0");
        }
    } else if(strcasecmp("surface",TemProfOpt)==0)
    {
        _tinfo->tem_type = const_profile;
        _tinfo->surf_tem = GetValueDk(ifp,"target.temperature",1,"293.0");
    }
    else if(strcasecmp("conductive",TemProfOpt)==0)
    {
        _tinfo->tem_type = conductive_profile;
        _tinfo->surf_tem = GetValueD(ifp,"target.surf_temperature","-1.0");
        if(_tinfo->surf_tem < 0.)
        {
            _tinfo->surf_tem = GetValueDk(ifp,"target.temperature",1,"293.0");
        }
        _tinfo->surf_dtdz = GetValueD(ifp,"target.surf_dtdz","0.0");
    } else if(strcasecmp("convective",TemProfOpt)==0)
    {
        _tinfo->tem_type = convective_profile;
        _tinfo->surf_gravity = GetValueD(ifp,"target.surf_gravity","-1.63");
        _tinfo->surf_dtdz = GetValueD(ifp,"target.surf_dtdz","0.0");
        _tinfo->t_lith = GetValueD(ifp,"target.t_lith","-1.0");
        _tinfo->d_lith = GetValueD(ifp,"target.d_lith","-1.0");
        _tinfo->r_planet = GetValueD(ifp,"target.r_planet","1738.0e3");
        _tinfo->d_mantle = GetValueD(ifp,"target.d_mantle","-1.0");
        _tinfo->t_mantle = GetValueD(ifp,"target.t_mantle","-1.0");
        _tinfo->surf_tem = GetValueD(ifp,"target.surf_temperature","-1.0");
        if(_tinfo->surf_tem < 0.)
        {
            _tinfo->surf_tem = GetValueDk(ifp,"target.temperature",1,"293.0");
        }
        _tinfo->center[Y] = -_tinfo->r_planet;
        _tinfo->center[X] =  0.0;
        ConsistentLith(_tinfo);
    } else if(strcasecmp("lunar",TemProfOpt)==0)
    {
        _tinfo->tem_type = lunar_profile;
        _tinfo->surf_gravity = GetValueD(ifp,"target.surf_gravity","-1.63");
        _tinfo->surf_dtdz = GetValueD(ifp,"target.surf_dtdz","0.0");
        _tinfo->t_lith = GetValueD(ifp,"target.t_lith","-1.0");
        _tinfo->d_lith = GetValueD(ifp,"target.d_lith","-1.0");
        _tinfo->r_planet = GetValueD(ifp,"target.r_planet","1738.0e3");
        _tinfo->d_mantle = GetValueD(ifp,"target.d_mantle","560.0e3");
        _tinfo->t_mantle = GetValueD(ifp,"target.t_mantle","-1.0");
        _tinfo->surf_tem = GetValueD(ifp,"target.surf_temperature","-1.0");
        if(_tinfo->surf_tem < 0.)
        {
            _tinfo->surf_tem = GetValueDk(ifp,"target.temperature",1,"293.0");
        }
        _tinfo->center[Y] = -_tinfo->r_planet;
        _tinfo->center[X] =  0.0;
        ConsistentLith(_tinfo);
    }
    else
    {
        fprintf(stderr,"unknown target temperature profile: %s \n",TemProfOpt);
    }
}

void load_target_addition_info(InputFile * ifp, target_addition_info * _tainfo)
{
    _tainfo->nobject = GetValueI(ifp,"addition.nobject","0");
    // load name of materials,
    char material_name[MAXMAT][80];
    int nmat = GetValueI(ifp,"material.nm","3");
    for(int matid=0;matid<nmat;++matid)
    {
        GetValueSk(ifp,"material.name",material_name[matid],matid,"vacuum_");
    }

    for(int k_object=0; k_object<_tainfo->nobject;++k_object)
    {
        // check the material index of addition
        char tmp_name[80];
        GetValueSk(ifp,"addition.material",tmp_name,k_object,"unknown");
        _tainfo->material[k_object] = -1;
        for(int matid=0;matid<nmat;++matid)
        {
            if(0==strcasecmp(tmp_name,material_name[matid]))
            {
                _tainfo->material[k_object] = matid;
                break;
            }
        }

        if(material_name[k_object] < 0)
        {
            fprintf(stderr,"cannot get material id of %s in target addition\n",tmp_name);
            exit(0);
        }

        char tmp_type[80], parameter_name[80];
        GetValueSk(ifp,"addition.type",tmp_type,k_object,"unknown");
        GetValueSk(ifp,"addition.parameter",parameter_name,k_object,"unknown");
        char parameter_handle[100];
        snprintf(parameter_handle,100,"addition.%s",parameter_name);
        int r = SearchInput(ifp, parameter_handle);
        if(r<0)
        {
            fprintf(stdout,"can not find parameter:%s\n", parameter_handle);
            exit(0);
        }

        if(strcasecmp("cube",tmp_type) == 0)
        {
            // left-bottom position (lbx,lby)
            // height vector (hx,hy)
            // width vector (wx,wy)
            // velocity (vx,vy)
            // temperature (t)
            // pressure (p)
            // damage (d)
            // porosity
            _tainfo->addition_type[k_object] = cube_addition;
            for(int j=0;j<12;++j)
            {
                _tainfo->parmeters[k_object][j] = GetValueDk(ifp,parameter_handle,j,"-1");
            }
        }
        else if(strcasecmp("sphere",tmp_type)==0)
        {
            // center (cx,cy)
            // radius (r)
            // velcoity (vx,vy)
            // temperature (t)
            // pressure (p)
            // damage (d)
            // porosity
            _tainfo->addition_type[k_object] = sphere_addition;
            for(int j=0;j<9;++j)
            {
                _tainfo->parmeters[k_object][j] = GetValueDk(ifp,parameter_handle,j,"-1");
            }
        }
        else
        {
            _tainfo->addition_type[k_object] = unknown_addition;
        }
    }
}

void init_material_distribution(mesh2d_info * _minfo, sale2d_var * _sale)
{
    if(_minfo->tinfo.layer_type == plane_target)
    {
        init_material_plane_distribution(_minfo,_sale);
    }
    else if(_minfo->tinfo.layer_type == sphere_target)
    {
        init_material_sphere_distribution(_minfo,_sale);
    }
    else
    {
        fprintf(stdout,"undefined layer type in %s\n",__func__);
    }
    init_material_target_addition(_minfo,_sale);
}
void init_material_plane_distribution(mesh2d_info * _minfo, sale2d_var * _sale)
{
    /*
     * the position of every vertex (v_pos) and element (e_pos) has been set in init_mesh_pos().
     * parameters (e_suf, e_vof, e_grad) related to position is initialized at the same time.
     *
     * this subroutine set the initialize values for
     * STEP 1/ e_vof, v_vel, e_pre, e_tem, e_dam FROM mesh info
     * STEP 2/ m_vof = e_vof, m_den, m_eng FROM interpolate_eos_*
     * STEP 3/ e_csd, e_den, e_eng, e_vel FROM *
     */

    // build the target profile
    target_profile targetProfile;
    build_target_profile(&targetProfile,_minfo,_sale);

    // region, STEP 1
    proc_info_2d * e_proc = _sale->e_vof->proc_self;
    int e_nx = e_proc->nx, e_ny = e_proc->ny;
    int nmat = _sale->nmat;


    for(int k=0;k<e_nx;++k)
    {
        for(int j=0;j<e_ny;++j)
        {
            fv_id_2d elid = {.x=k,.y=j};
            // in the initialization stage, all the elements in this progress,
            // none is skipped recording to the node_type
            char * epos_cptr, *evof_cptr, *epre_cptr, *etem_cptr,*edam_cptr,*ewpt_cptr;
            _sale->foe->get_data(&elid,&epos_cptr,_sale->e_pos);
            _sale->foe->get_data(&elid,&evof_cptr,_sale->e_vof);
            _sale->foe->get_data(&elid,&epre_cptr,_sale->e_pre);
            _sale->foe->get_data(&elid,&etem_cptr,_sale->e_tem);
            _sale->foe->get_data(&elid,&edam_cptr,_sale->e_dam);
            _sale->foe->get_data(&elid,&ewpt_cptr,_sale->e_wpt);
            double * epos = (double *)epos_cptr;
            double * evof = (double *)evof_cptr;
            double * epre = (double *)epre_cptr;
            double * etem = (double *)etem_cptr;
            double * edam = (double *)edam_cptr;
            double * ewpt = (double *)ewpt_cptr;


            double * evel = NULL, *mden = NULL, *meng=NULL;
            _sale->foe->get_double(&elid,&mden,_sale->m_den);
            _sale->foe->get_double(&elid,&meng,_sale->m_eng);
            _sale->foe->get_double(&elid,&evel,_sale->e_vel);
            vec_zero(evel,2);
            double *etps,*evib;
            _sale->foe->get_double(&elid,&etps,_sale->e_tps);
            _sale->foe->get_double(&elid,&evib,_sale->e_vib);
            etps[0] = 0.;
            evib[0] = 0.;

            double *mpty_ptr,*ecvs_ptr;
            _sale->foe->get_double(&elid,&mpty_ptr,_sale->m_pty);
            _sale->foe->get_double(&elid,&ecvs_ptr,_sale->e_cvs);

            vec_number(mpty_ptr,nmat,1.0);
            *ecvs_ptr = 0.;

            double projectile_dis = vec_dis(epos,_minfo->pinfo.center,2);
            if(projectile_dis <= _minfo->pinfo.radiu)
            {
                // this element is in projectile
                // set the temperature, pressure, damage FROM projectile_info
                vec_zero(evof,nmat);
                int projectile_id = _minfo->pinfo.material;
                evof[projectile_id] = 1.0;
                *epre = _minfo->pinfo.pressure;
                *etem = _minfo->pinfo.temperature;
                *edam = _minfo->pinfo.damage;

                vec_zero(mden,nmat);
                vec_zero(meng,nmat);

                double projectile_pty = _sale->etb[projectile_id].ref->alpha0;
                mden[projectile_id] = interpolate_eos_tp(_sale->etb+projectile_id,
                                                         _minfo->pinfo.temperature,
                                                         _minfo->pinfo.pressure*projectile_pty,
                                                         EOSDEN)/projectile_pty;
                meng[projectile_id] = interpolate_eos_tp(_sale->etb+projectile_id,
                                                         _minfo->pinfo.temperature,
                                                         _minfo->pinfo.pressure*projectile_pty,
                                                         EOSENG)/projectile_pty;
                // velocity
                vec_copy(_minfo->pinfo.velocity, evel, 2);
                // set wpt
                ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = -0.25;
                mpty_ptr[projectile_id] = projectile_pty;

            } else if(epos[Y] <= 0.)
            {
                // this element is in the target layers

                int layer_index = (int) floor(fabs(epos[Y])/_minfo->tinfo.max_depth * (targetProfile.n - 1)) + 1;
                layer_index = targetProfile.n - 1 - layer_index;

                layer_index = MIN(layer_index,targetProfile.n - 1);
                layer_index = MAX(layer_index,0);

                int layer_id = targetProfile.mat[layer_index];

                vec_zero(evof,nmat);
                vec_zero(mden,nmat);
                vec_zero(meng,nmat);

                *epre = targetProfile.pre[layer_index];
                *etem = targetProfile.tem[layer_index];
                *edam = targetProfile.dam[layer_index];

                evof[layer_id] = 1.0;
                mden[layer_id] = targetProfile.den[layer_index];
                meng[layer_id] = targetProfile.eng[layer_index];
                mpty_ptr[layer_id] =  _sale->etb[layer_id].ref->alpha0;
                ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = -0.25;
            } else
            {
                // this element is in vaccum, just set e_vof
                vec_zero(evof,nmat);
                evof[0] = 1.0;
                ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = 1.0;
            }
        }
    }

    // loop in vertexs , set gravity
    // estimate the transient crater using gravity regime
    // Holsapple, Keith A., and Kevin R. Housen. “A Crater and Its Ejecta: An Interpretation of Deep Impact."
    double t_a = _minfo->pinfo.radiu;
    double t_g = vec_len(_minfo->gravity,2);
    double t_U = fabs(_minfo->pinfo.velocity[Y]);
    double t_delta = _sale->etb[_minfo->pinfo.material].ref->den0;
    int n_transient = (int) floor(20.0*t_a/targetProfile.dz);
    double t_rho = vec_sum(targetProfile.den-n_transient+targetProfile.n,n_transient-1)/(n_transient-1);
    double transient_R = 1.03*t_a * pow(t_g*t_a/(t_U*t_U),-0.170)*pow(t_delta/t_rho,0.332);

    for(int k=0;k<_minfo->v_nx;++k)
    {
        for(int j=0;j<_minfo->v_ny;++j)
        {
            double * vbf = NULL;
            proc_info_2d * proc_domain =  (proc_info_2d *)(_sale->v_bf->domain);
            fv_id_2d vlid = {.x=k,.y=j}, vgid = {0,0};
            _sale->fov->lid2gid(_sale->v_bf->rank,&vlid,&vgid,_sale->v_bf);
            _sale->fov->get_double(&vlid,&vbf,_sale->v_bf);
            vec_copy(_minfo->gravity,vbf,2);

            if(vgid.y <= _minfo->ypad) vbf[Y] = 0.;

            double * vlaplace_ptr;
            double * vpos_ptr;
            _sale->fov->get_double(&vlid,&vlaplace_ptr,_sale->v_laplace);
            _sale->fov->get_double(&vlid,&vpos_ptr,_sale->v_pos);
            double * vdebug_ptr;
            _sale->fov->get_double(&vlid,&vdebug_ptr,_sale->v_debug);

            // reduce the anc effect when cell is far (0,0)
            double anc_reduce = 0.;
            double ancd = _minfo->anc_depth;
            double ancw = _minfo->anc_width;
            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};
            double *esigma_ptr[4];
            for(int i=0;i<4;++i)
            {
                _sale->foe->get_double(elid[i],&esigma_ptr[i],_sale->e_sigma);

            }
            // calculate damp factor on cells
            double vsigma[2] = {0.};
            for(int i=0;i<4;++i)
            {
                if(NULL == esigma_ptr[i]) continue;
                scaler_add(vsigma,2,esigma_ptr[i],0.25);
            }

            if(strcasecmp(_minfo->anc_pattern,"pml") == 0)
            {
                if(_minfo->damp_y >= 1  && vgid.y <= 2*_minfo->damp_y)
                {
                    double local_x = 0.5*vgid.y/_minfo->damp_y;
                    anc_reduce = 2.0*Max(anc_reduce,2.0*hanning(0.5*(1.0 - local_x)));
                }

                if(_minfo->damp_x>= 1 && vgid.x >= proc_domain->nx - 2*_minfo->damp_x)
                {
                    double local_x = 0.5*(proc_domain->nx - 2*_minfo->damp_x - vgid.x)/_minfo->damp_x;
                    anc_reduce = 2.0*Max(anc_reduce,2.0*hanning(0.5*local_x));
                }
            }
            else if(strcasecmp(_minfo->anc_pattern,"uniform") == 0)
            {
                anc_reduce = 1.0;
            }
            else if(strcasecmp(_minfo->anc_pattern,"scaling") == 0)
            {
                if(ancd > 0. && ancw > 0.)
                {
                    double ancX = vpos_ptr[X]/transient_R;
                    double ancY = vpos_ptr[Y]/transient_R;
                    if( ancX <= (ancw + 1) &&  ancY <= -1.0*ancd)
                    {
                        anc_reduce = 1.0;
                        if(ancX >= ancw)
                        {
                            anc_reduce *= 1.0 - hanning(0.5*(ancX - ancw));
                        }

                        if(ancY >= -1.0*(ancd + 1))
                        {
                            anc_reduce *= hanning(0.5*(-ancY-ancd));
                        }
                    }
                }
                else
                {
                    anc_reduce = 1.0;
                }

                // part of pml, ..., unstable if commented
                if(_minfo->damp_y >= 1  && vgid.y <= _minfo->damp_y)
                {
                    anc_reduce = Max(anc_reduce,2.0*hanning(0.5 - 0.5*vgid.y/_minfo->damp_y));
                }

                if(_minfo->damp_x>= 1 && vgid.x >= proc_domain->nx - _minfo->damp_x)
                {
                    anc_reduce = Max(anc_reduce,2.0*hanning(0.5*(proc_domain->nx - _minfo->damp_x - vgid.x)/_minfo->damp_x));
                }
            }
            else if(strcasecmp(_minfo->anc_pattern,"none") == 0)
            {
                anc_reduce = 0.;
            }
            else
            {
                fprintf(stdout,"undefined anc pattern:%s\n",_minfo->anc_pattern);
                exit(0);
            }
            vec_scaler(vlaplace_ptr,9,anc_reduce);
        }
    }

    UnallocateTargetProfile(&targetProfile);
}
void init_material_sphere_distribution(mesh2d_info * _minfo, sale2d_var * _sale)
{
    /*
     * init the materials distribution as sphere target
     * the position of every vertex (v_pos) and element (e_pos) has been set in init_mesh_pos().
     * parameters (e_suf, e_vof, e_grad) related to position is initialized at the same time.
     *
     * this subroutine set the initialize values for
     * STEP 1/ e_vof, v_vel, e_pre, e_tem, e_dam FROM mesh info
     * STEP 2/ m_vof = e_vof, m_den, m_eng FROM interpolate_eos_*
     * STEP 3/ e_csd, e_den, e_eng, e_vel FROM *
     */

    // build the target profile
    target_profile targetProfile;
    build_target_profile(&targetProfile,_minfo,_sale);

    // region, STEP 1
    proc_info_2d * e_proc = _sale->e_vof->proc_self;
    int e_nx = e_proc->nx, e_ny = e_proc->ny;
    int nmat = _sale->nmat;


    for(int k=0;k<e_nx;++k)
    {
        for(int j=0;j<e_ny;++j)
        {
            fv_id_2d elid = {.x=k,.y=j};
            // in the initialization stage, all the elements in this progress,
            // none is skipped recording to the node_type
            char * epos_cptr, *evof_cptr, *epre_cptr, *etem_cptr,*edam_cptr,*ewpt_cptr;
            _sale->foe->get_data(&elid,&epos_cptr,_sale->e_pos);
            _sale->foe->get_data(&elid,&evof_cptr,_sale->e_vof);
            _sale->foe->get_data(&elid,&epre_cptr,_sale->e_pre);
            _sale->foe->get_data(&elid,&etem_cptr,_sale->e_tem);
            _sale->foe->get_data(&elid,&edam_cptr,_sale->e_dam);
            _sale->foe->get_data(&elid,&ewpt_cptr,_sale->e_wpt);
            double * epos = (double *)epos_cptr;
            double * evof = (double *)evof_cptr;
            double * epre = (double *)epre_cptr;
            double * etem = (double *)etem_cptr;
            double * edam = (double *)edam_cptr;
            double * ewpt = (double *)ewpt_cptr;


            double * evel = NULL, *mden = NULL, *meng=NULL;
            _sale->foe->get_double(&elid,&mden,_sale->m_den);
            _sale->foe->get_double(&elid,&meng,_sale->m_eng);
            _sale->foe->get_double(&elid,&evel,_sale->e_vel);
            vec_zero(evel,2);
            double *etps,*evib;
            _sale->foe->get_double(&elid,&etps,_sale->e_tps);
            _sale->foe->get_double(&elid,&evib,_sale->e_vib);
            etps[0] = 0.;
            evib[0] = 0.;

            double *mpty_ptr,*ecvs_ptr;
            _sale->foe->get_double(&elid,&mpty_ptr,_sale->m_pty);
            _sale->foe->get_double(&elid,&ecvs_ptr,_sale->e_cvs);

            vec_number(mpty_ptr,nmat,1.0);
            *ecvs_ptr = 0.;

            double projectile_dis = vec_dis(epos,_minfo->pinfo.center,2);
            double planet_dis = vec_dis(epos,_minfo->tinfo.center,2);

            if(projectile_dis <= _minfo->pinfo.radiu)
            {
                // this element is in projectile
                // set the temperature, pressure, damage FROM projectile_info
                vec_zero(evof,nmat);
                int projectile_id = _minfo->pinfo.material;
                evof[projectile_id] = 1.0;
                *epre = _minfo->pinfo.pressure;
                *etem = _minfo->pinfo.temperature;
                *edam = _minfo->pinfo.damage;

                vec_zero(mden,nmat);
                vec_zero(meng,nmat);

                double projectile_pty = _sale->etb[projectile_id].ref->alpha0;
                mden[projectile_id] = interpolate_eos_tp(_sale->etb+projectile_id,
                                                         _minfo->pinfo.temperature,
                                                         _minfo->pinfo.pressure*projectile_pty,
                                                         EOSDEN)/projectile_pty;
                meng[projectile_id] = interpolate_eos_tp(_sale->etb+projectile_id,
                                                         _minfo->pinfo.temperature,
                                                         _minfo->pinfo.pressure*projectile_pty,
                                                         EOSENG)/projectile_pty;
                // velocity
                vec_copy(_minfo->pinfo.velocity, evel, 2);
                // set wpt
                ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = -0.25;
                mpty_ptr[projectile_id] = projectile_pty;

            }
            else if(planet_dis <= _minfo->tinfo.r_planet)
            {
                // this element is in the target layers
                int layer_index = (int) floor(planet_dis/targetProfile.dz);
                layer_index = MIN(layer_index,targetProfile.n - 1);
                layer_index = MAX(layer_index,0);

                int layer_id = targetProfile.mat[layer_index];

                vec_zero(evof,nmat);
                vec_zero(mden,nmat);
                vec_zero(meng,nmat);

                *epre = targetProfile.pre[layer_index];
                *etem = targetProfile.tem[layer_index];
                *edam = targetProfile.dam[layer_index];

                evof[layer_id] = 1.0;
                mden[layer_id] = targetProfile.den[layer_index];
                meng[layer_id] = targetProfile.eng[layer_index];
                mpty_ptr[layer_id] =  _sale->etb[layer_id].ref->alpha0;
                ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = -0.25;
            } else
            {
                // this element is in vaccum, just set e_vof
                vec_zero(evof,nmat);
                evof[0] = 1.0;
                ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = 1.0;
            }
        }
    }

    // loop in vertexs ,
    for(int k=0;k<_minfo->v_nx;++k)
    {
        for(int j=0;j<_minfo->v_ny;++j)
        {
            double * vbf = NULL;
            proc_info_2d * proc_domain =  (proc_info_2d *)(_sale->v_bf->domain);
            fv_id_2d vlid = {.x=k,.y=j}, vgid = {0,0};
            _sale->fov->lid2gid(_sale->v_bf->rank,&vlid,&vgid,_sale->v_bf);
            _sale->fov->get_double(&vlid,&vbf,_sale->v_bf);

            double * vlaplace_ptr;
            double * vpos_ptr;
            _sale->fov->get_double(&vlid,&vlaplace_ptr,_sale->v_laplace);
            _sale->fov->get_double(&vlid,&vpos_ptr,_sale->v_pos);
            // reduce the anc effect when cell is far from (0,0), !!! anc is closed in sphere model
            double anc_reduce = 0.;
            vec_scaler(vlaplace_ptr,9,anc_reduce);

            // set gravity
            double planet_dis = vec_dis(vpos_ptr,_minfo->tinfo.center,2);
            double local_gravity[2] = {0., 0.};
            liner_op(local_gravity,2,vpos_ptr,_minfo->tinfo.center,1.0,-1.0);
            vec_norm(local_gravity,2);
            int layer_index = (int) floor(planet_dis/targetProfile.dz);
            layer_index = MIN(layer_index,targetProfile.n - 1);
            layer_index = MAX(layer_index,0);
            scaler_move(vbf,2,local_gravity,targetProfile.gra[layer_index]);
            if(planet_dis >= _minfo->tinfo.r_planet)
            {
                vec_scaler(vbf,2,(_minfo->tinfo.r_planet/planet_dis)*(_minfo->tinfo.r_planet/planet_dis));
            }
        }
    }

    UnallocateTargetProfile(&targetProfile);
}

void init_material_target_addition(mesh2d_info * _minfo, sale2d_var * _sale)
{
    /*
     * init the materials distribution of target addition
     * the position of every vertex (v_pos) and element (e_pos) has been set in init_mesh_pos().
     * parameters (e_suf, e_vof, e_grad) related to position is initialized at the same time.
     *
     * this subroutine set the initialize values for
     * STEP 1/ e_vof, v_vel, e_pre, e_tem, e_dam FROM mesh info
     * STEP 2/ m_vof = e_vof, m_den, m_eng FROM interpolate_eos_*
     * STEP 3/ e_csd, e_den, e_eng, e_vel FROM *
     */

    // target addition info in _minfo->tainfo
    // region, STEP 1
    proc_info_2d * e_proc = _sale->e_vof->proc_self;
    int e_nx = e_proc->nx, e_ny = e_proc->ny;
    int nmat = _sale->nmat;

    for(int k=0;k<e_nx;++k)
    {
        for(int j=0;j<e_ny;++j)
        {
            fv_id_2d elid = {.x=k,.y=j};
            // in the initialization stage, all the elements in this progress,
            // none is skipped recording to the node_type
            char * epos_cptr, *evof_cptr, *epre_cptr, *etem_cptr,*edam_cptr,*ewpt_cptr;
            _sale->foe->get_data(&elid,&epos_cptr,_sale->e_pos);
            _sale->foe->get_data(&elid,&evof_cptr,_sale->e_vof);
            _sale->foe->get_data(&elid,&epre_cptr,_sale->e_pre);
            _sale->foe->get_data(&elid,&etem_cptr,_sale->e_tem);
            _sale->foe->get_data(&elid,&edam_cptr,_sale->e_dam);
            _sale->foe->get_data(&elid,&ewpt_cptr,_sale->e_wpt);
            double * epos = (double *)epos_cptr;
            double * evof = (double *)evof_cptr;
            double * epre = (double *)epre_cptr;
            double * etem = (double *)etem_cptr;
            double * edam = (double *)edam_cptr;
            double * ewpt = (double *)ewpt_cptr;


            double * evel = NULL, *mden = NULL, *meng=NULL;
            _sale->foe->get_double(&elid,&mden,_sale->m_den);
            _sale->foe->get_double(&elid,&meng,_sale->m_eng);
            _sale->foe->get_double(&elid,&evel,_sale->e_vel);
            vec_zero(evel,2);
            double *etps,*evib;
            _sale->foe->get_double(&elid,&etps,_sale->e_tps);
            _sale->foe->get_double(&elid,&evib,_sale->e_vib);
            etps[0] = 0.;
            evib[0] = 0.;

            double *mpty_ptr,*ecvs_ptr;
            _sale->foe->get_double(&elid,&mpty_ptr,_sale->m_pty);
            _sale->foe->get_double(&elid,&ecvs_ptr,_sale->e_cvs);

            vec_number(mpty_ptr,nmat,1.0);
            *ecvs_ptr = 0.;


            for(int k_object=0;k_object<_minfo->tainfo.nobject;++k_object)
            {
                int oid = _minfo->tainfo.material[k_object];
                double tem,pre,dam,pty,vel[2];
                if(_minfo->tainfo.addition_type[k_object] == cube_addition)
                {
                    double lbp[2] = {_minfo->tainfo.parmeters[k_object][0],_minfo->tainfo.parmeters[k_object][1]};
                    double xa[2] = {_minfo->tainfo.parmeters[k_object][2],_minfo->tainfo.parmeters[k_object][3]};
                    double ya[2] = {_minfo->tainfo.parmeters[k_object][4],_minfo->tainfo.parmeters[k_object][5]};
                    vel[0] = _minfo->tainfo.parmeters[k_object][6];
                    vel[1] = _minfo->tainfo.parmeters[k_object][7];
                    tem = _minfo->tainfo.parmeters[k_object][8];
                    pre = _minfo->tainfo.parmeters[k_object][9];
                    dam = _minfo->tainfo.parmeters[k_object][10];
                    pty = _minfo->tainfo.parmeters[k_object][11];

                    double x0[2] = {0,0};
                    liner_op(x0,2,epos,lbp,1,-1);
                    double lxa = vec_len(xa,2);
                    double lya = vec_len(ya,2);
                    double xl = vec_dot(x0,xa,2)/lxa/lxa;
                    double yl = vec_dot(x0,ya,2)/lya/lya;


                    if(xl*(1-xl) <0 || yl*(1-yl)<0)
                        continue;
                }
                else if(_minfo->tainfo.addition_type[k_object] == sphere_addition)
                {
                    double center[2] = {_minfo->tainfo.parmeters[k_object][0],_minfo->tainfo.parmeters[k_object][1]};
                    double radius = _minfo->tainfo.parmeters[k_object][2];
                    vel[0] = _minfo->tainfo.parmeters[k_object][3];
                    vel[1] = _minfo->tainfo.parmeters[k_object][4];
                    tem = _minfo->tainfo.parmeters[k_object][5];
                    pre = _minfo->tainfo.parmeters[k_object][6];
                    dam = _minfo->tainfo.parmeters[k_object][7];
                    pty = _minfo->tainfo.parmeters[k_object][8];

                    double object_distance = vec_dis(center, epos, 2);
                    if(object_distance > radius)
                        continue;
                }
                else if(_minfo->tainfo.addition_type[k_object] == unknown_addition)
                {
                    fprintf(stdout,"no method to initialize for unknown_addition type %d\n",k_object);
                    exit(0);
                }
                else
                {
                    fprintf(stdout,"no method to initialize for addition %d\n",k_object);
                    exit(0);
                }

                tem = tem>=0? tem: etem[0];
                pre = pre>=0? pre: epre[0];
                pty = pty>=0? pty: mpty_ptr[oid];
                dam = dam>=0? dam: edam[0];

                if(oid == VACUUM)
                {
                    vec_zero(evof,nmat);
                    vec_zero(mden,nmat);
                    vec_zero(meng,nmat);
                    evof[VACUUM] = 1.0;
                    ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = 1.0;
                }
                else
                {
                    vec_zero(evof,nmat);
                    vec_zero(mden,nmat);
                    vec_zero(meng,nmat);

                    *epre = pre;
                    *etem = tem;
                    *edam = dam;

                    evof[oid] = 1.0;
                    mden[oid] = interpolate_eos_tp(_sale->etb+oid,
                                                             tem,
                                                             pre,
                                                             EOSDEN)/pty;
                    meng[oid] = interpolate_eos_tp(_sale->etb+oid,
                                                             tem,
                                                             pre,
                                                             EOSENG)/pty;
                    mpty_ptr[oid] =  pty;
                    ewpt[0] = ewpt[1] = ewpt[2] = ewpt[3] = -0.25;
                    vec_copy(vel, evel, 2);
                }
            }
        }
    }

}



void build_target_profile(target_profile * _tprof, mesh2d_info * _minfo, sale2d_var * _sale)
{
    target_info * _tinfo = &(_minfo->tinfo);
    state_reference * _ref = _minfo->ref;
    eos_table * _etb = _sale->etb;
    // build the profile used for initialize

    _tprof->dz   = Min(_minfo->dy * 0.5, _tinfo->max_depth/_tinfo->n_profile);
    // _tprof->dz = _minfo->dy*0.5;
    _tprof->n = (int) ceil(_tinfo->max_depth/_tprof->dz);
    _tprof->type = _tinfo->layer_type;
    _tprof->info = _tinfo; // the info struct have directly put in target_profile
    AllocateTargetProfile(_tprof);
    SetTargetProfile_Mat(_tprof,_tinfo->depth,_tinfo->material,_tinfo->n_layer);
    SetTargetProfile_Tem(_tprof,_ref,_tinfo);
    SetTargetProfile_Gra(_tprof,_minfo->gravity[Y]);

    int MaxItr = 2000;
    // if(_sale->porosity_model) MaxItr*=10;

    for(int cycle=0;cycle<MaxItr;++cycle)
    {
        IntegrateGravity(_tprof);
        CorrectDensity(_tprof,_etb,_ref);
    }

    IntegrateGravity(_tprof);

    if(_sale->e_vof->rank == 1)
    {
        FILE * fp = fopen("profile.txt","w");
        //WriteTargetProfile(fp,_tprof);
        WriteTargetProfile_All(fp,_tprof,_ref,_tinfo);
        fclose(fp);
    }
}