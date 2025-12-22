//
// Created by huacheng on 3/15/23.
//

#include "strength.h"

void update_e_stress(sale2d_var * _sale, node_type _nt, double dt)
{
    /*
     * calculate the increment of stress (only elastic stress)
     * VARIABLES required pre-calculated
     *      e_grad_vel
     * VARIABLES overwritten:
     *      e_stress
     *
     */

    field_var * _var = _sale->e_ste;
    field_opt * _opt = _sale->foe;

    proc_info_2d * _proc = (proc_info_2d *) _var->proc_self;
    int nx = _proc->nx, ny = _proc->ny;


    for(int index=0;index<nx*ny;++index)
    {
        int k = INDEX2LIDX(index,nx,ny);
        int j = INDEX2LIDY(index,nx,ny);
        fv_id_2d lid = {k,j};
        node_type lnt = useless_section;
        _opt->get_nty(&lid,&lnt,_var);
        if(!(_nt&lnt)) continue;

        double *evof_ptr, *este_ptr, *eden_ptr;
        _opt->get_double(&lid,&evof_ptr,_sale->e_vof);
        _opt->get_double(&lid,&este_ptr,_sale->e_ste);
        _opt->get_double(&lid,&eden_ptr,_sale->e_den);

        double * esta2s_ptr;
        _opt->get_double(&lid,&esta2s_ptr,_sale->e_sta2s);
        double * este2s_ptr;
        _opt->get_double(&lid,&este2s_ptr,_sale->e_ste2s);
        double * etps_ptr;
        _opt->get_double(&lid,&etps_ptr,_sale->e_tps);
        double *estrength_ptr;
        _opt->get_double(&lid,&estrength_ptr,_sale->e_strength);
        double * edam_ptr=NULL;
        _opt->get_double(&lid,&edam_ptr,_sale->e_dam);

        if(evof_ptr[VACUUM] > TOLVOF)
        {
            vec_zero(este_ptr,4);
            este2s_ptr[0] = 0.;
            esta2s_ptr[0] = 0.;
            estrength_ptr[0] = 0.;
            continue;
        }


        // check the state of local cell
        double vapor_den_ref = 0.;
        double melt_den_ref  = 0.;
        for(int matid=1;matid<_sale->nmat;++matid)
        {
            if(evof_ptr[matid] < TOLVOF) continue;
            vapor_den_ref += evof_ptr[matid]*_sale->etb[matid].ref->VaporDen;
            melt_den_ref  += evof_ptr[matid]*_sale->etb[matid].ref->MeltDen;
        }

        // check the porosity
        double material_porosity = 0.0;
        if(_sale->porosity_model)
        {
            double material_vof = 0.0;
            double *mpty_ptr=NULL;
            _opt->get_double(&lid,&mpty_ptr,_sale->m_pty);
            for (int matid=1;matid<_sale->nmat;++matid)
            {
                if(evof_ptr[matid] < TOLVOF) continue;
                material_porosity += evof_ptr[matid]*mpty_ptr[matid];
                material_vof += evof_ptr[matid];
            }
            material_porosity /= material_vof;
        }
        else
        {
            material_porosity = 1.0;
        }
        material_porosity = Max(material_porosity, 1.0);

        // apply porosity to vapor/melt ref
        vapor_den_ref /= material_porosity;
        melt_den_ref /= material_porosity;

        int solid_cell = 0;
        int melt_cell  = 0;
        int vapor_cell = 0;
        double este2s_old = este2s_ptr[0]; // before update the stress record the j2 of old stress
        if(eden_ptr[0] <= vapor_den_ref)
        {
            // this cell is full of vapor
            // clear tps and damage
            vec_zero(este_ptr,4);
            este2s_ptr[0] = 0.;
            etps_ptr[0] = 0.;
            vapor_cell = 1;
            edam_ptr[0] = 0.;
            estrength_ptr[0] = 0.;
        } else if(eden_ptr[0] <= melt_den_ref)
        {
            // melt cell
            vec_zero(este_ptr,4);
            este2s_ptr[0] = 0.;
            melt_cell = 1;
            double *egradvel_ptr;
            _opt->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);

            double strain_increment[4] = {egradvel_ptr[0],0.5*egradvel_ptr[1]+0.5*egradvel_ptr[2],
                                          egradvel_ptr[3],egradvel_ptr[4]};

            double strain_tr = strain_increment[0] + strain_increment[2] + strain_increment[3];
            strain_increment[0] -= strain_tr/3.0;
            strain_increment[2] -= strain_tr/3.0;
            strain_increment[3] -= strain_tr/3.0;

            /*
             * set strain rate and total plastic rate
             */
            esta2s_ptr[0] = sqrt(second_inv(strain_increment));
            etps_ptr[0] = Min(etps_ptr[0] + esta2s_ptr[0] * dt,10.);
            edam_ptr[0] = 1.;
            estrength_ptr[0] = 0.;
        } else
        {
            solid_cell = 1;
            double *egradvel_ptr, *emu_ptr;
            _opt->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);
            _opt->get_double(&lid,&emu_ptr,_sale->e_mu);

            double strain_increment[4] = {egradvel_ptr[0],0.5*egradvel_ptr[1]+0.5*egradvel_ptr[2],
                                          egradvel_ptr[3],egradvel_ptr[4]};

            double strain_tr = strain_increment[0] + strain_increment[2] + strain_increment[3];
            strain_increment[0] -= strain_tr/3.0;
            strain_increment[2] -= strain_tr/3.0;
            strain_increment[3] -= strain_tr/3.0;

            /*
             * wilkins rotation
             */
#ifdef WIKINS_RORARTION
            double sin2Omega = (egradvel_ptr[2] - egradvel_ptr[1])*dt;
            sin2Omega = Wind(sin2Omega,-1.,1.);
            double cos2Omega1 = sqrt(1.0 - sin2Omega*sin2Omega) - 1.0;
            double ss2 = (este_ptr[0] - este_ptr[2])/2.0;
            double stress_rot[4] = {ss2*cos2Omega1 + este_ptr[1]*sin2Omega,
                                    este_ptr[1]*cos2Omega1 - ss2*sin2Omega,
                                    -ss2*cos2Omega1 - este_ptr[1]*sin2Omega,
                                    0.0};
            scaler_add(este_ptr,4,stress_rot,1.0);
#endif
            scaler_add(este_ptr,4,strain_increment,2.0*(*emu_ptr)*dt);


            /*
             * set strain rate and total plastic rate
             */
            este2s_ptr[0] = sqrt(second_inv(este_ptr));
            esta2s_ptr[0] = sqrt(second_inv(strain_increment));
        }

        /*
         * update damage or strength
         * merged from update_e_strength
         */

        double * epre_ptr=NULL;
        _opt->get_double(&lid,&epre_ptr,_sale->e_pre);
        double *etem_ptr=NULL;
        _opt->get_double(&lid,&etem_ptr,_sale->e_tem);
        double * evib_ptr;
        _opt->get_double(&lid,&evib_ptr,_sale->e_vib);
        double * ecsd_ptr;
        _opt->get_double(&lid,&ecsd_ptr,_sale->e_csd);
        double * epos_ptr;
        _opt->get_double(&lid,&epos_ptr,_sale->e_pos);

        if(vapor_cell)
        {
            edam_ptr[0] = 0.;
        } else if(melt_cell)
        {
            edam_ptr[0] = 1.0;
        } else
        {
            // calculate yield strength
            double local_shear_strength = -1.0;
            if(_sale->strength_mod == avg_strength)
            {
                local_shear_strength = 0.0;
            }
            
            double local_tensile_strength = 0.;
            for(int matid=1;matid<_sale->nmat;++matid)
            {
                if(evof_ptr[matid] < TOLVOF) continue;
                state_reference * local_ref = _sale->etb[matid].ref;
                strength_func local_strength_func = (strength_func)(local_ref->shear_strength_func);
                strength_func local_tensile_func  = (strength_func)(local_ref->tensile_strength_func);
                /*
                 * components of _paras in strength_func
                 * _paras[0] : pressure
                 * _paras[1] : damage
                 * _paras[2] : total plastic strain/ tps
                 * _paras[3] : strain rate second invariant/ sta2s
                 * _paras[4] : reserved
                 * _paras[5] : reserved
                 *
                 * thermal soft or low density model is considered in strength_func.
                 */
                double mat_y = local_strength_func((double []){*epre_ptr, *edam_ptr, *etps_ptr, *esta2s_ptr, *etem_ptr}, local_ref);
                if(_sale->acoustic_fluid && epre_ptr[0] >= 0.)
                {
                    // calculate vib pressure
                    double vibpressure = Min(eden_ptr[0]*ecsd_ptr[0]* sqrt(evib_ptr[0]),epre_ptr[0]);
                    double melt_temp = SimonMelt(epre_ptr[0],local_ref);
                    double soft_factor = OhnakaSoft(etem_ptr[0], local_ref, melt_temp);
                    soft_factor = Wind(soft_factor,0.0,1.0);
                    double local_pre = epre_ptr[0] - vibpressure*edam_ptr[0]*edam_ptr[0]* sqrt(soft_factor);
                    // substitute the pressure in strength_func + the viscosity of acoustic fluid
                    double mat_yvib = local_strength_func((double []){local_pre, *edam_ptr, *etps_ptr, *esta2s_ptr, *etem_ptr}, local_ref)
                                      + eden_ptr[0]*local_ref->Acvis*esta2s_ptr[0];
                    mat_y = Min(mat_yvib,mat_y);
                }


                // use average strength of materials in element (default)
                switch(_sale->strength_mod)
                {
                    case avg_strength:
                        local_shear_strength += evof_ptr[matid]* Max(mat_y,0.0);
                        break;
                    case min_strength:
                        if(local_shear_strength < 0.)
                            local_shear_strength = mat_y;
                        else
                            local_shear_strength = Min(local_shear_strength,mat_y);

                        break;
                    case max_strength:
                        if(local_shear_strength < 0.)
                            local_shear_strength = mat_y;
                        else
                            local_shear_strength = Max(local_shear_strength,mat_y);
                        break;
                    default:
                        local_shear_strength += evof_ptr[matid]* Max(mat_y,0.0);
                        break;
                }

                double mat_t = local_tensile_func((double []){*epre_ptr, *edam_ptr, *etps_ptr, *esta2s_ptr, *etem_ptr}, local_ref);
                local_tensile_strength += evof_ptr[matid] * mat_t;
            }

            // check for strength on tracer
            if(_sale->tracer_lubrication >= 1)
            {
                double * ecnd_ptr;
                _opt->get_double(&lid,&ecnd_ptr,_sale->e_cnd);

                if(*ecnd_ptr >= 1.0)
                {
                    double lub_strength = Min((*epre_ptr)*(_sale->tracer_lub_frac) + _sale->tracer_lub_y0,0.0);
                    local_shear_strength = Min(local_shear_strength,lub_strength);
                }
            }

            // save strength, not necessary
            *estrength_ptr = local_shear_strength;

            if(local_shear_strength < este2s_ptr[0])
            {

                /*
                 * when the materials break, change este2s to plastic strain rate
                 * if the previous step didn't break, the plastic rate < the strain rate
                 * .. break, este2d_old is the previous strength ~
                 * the plastic rate should be >= 0
                 */

                esta2s_ptr[0] = esta2s_ptr[0]/ Max(0.5,(este2s_ptr[0] - este2s_old)/(este2s_ptr[0] - local_shear_strength));


                /*
                 * shear failure condition, cut the stress tensor
                 */
                double sfac = local_shear_strength/este2s_ptr[0];
                double *este_ptr=NULL;
                _opt->get_double(&lid,&este_ptr,_sale->e_ste);
                vec_scaler(este_ptr,4,sfac);
                este2s_ptr[0] = local_shear_strength;


                /*
                 * update the damage
                * The general damage increase model
                * LINEAR　is default
                * The exponential model will be added in future code
                 * adopt from SALEc/manual D = Max(e_p/e_f,1.0)
                */
                double dSDamage = 0.;
                for(int matid=1;matid<_sale->nmat;++matid)
                {
                    state_reference * local_ref = _sale->etb[matid].ref;
                    strength_func dam_func = (strength_func) local_ref->damage_increment_func;
                    dSDamage += evof_ptr[matid] * dam_func((double []){*epre_ptr,*edam_ptr},local_ref);
                }
                edam_ptr[0] = Wind(edam_ptr[0]+ esta2s_ptr[0] * dSDamage*dt,0.0,1.0);
                etps_ptr[0] = Min(etps_ptr[0] + esta2s_ptr[0] * dt,10.);
            }

            /*
             * tensile strength if enable
             * add tensile failure and damage
             */
            if(_sale->tensile_failure)
            {
                double Sxx = este_ptr[0], Szz = este_ptr[2], Sxz = este_ptr[1], Sth = este_ptr[3];
                double Ts  = sqrt(Sxz*Sxz + 0.25*(Sxx - Szz)*(Sxx - Szz));
                double Tmax = Max(0.5*(Sxx + Szz) - Ts,0.5*(Sxx + Szz) + Ts);
                Tmax = Max(Tmax,Sth);

                if((Tmax - epre_ptr[0]) > local_tensile_strength)
                {
                    double dam_inc = 0.4 * ecsd_ptr[0] / Min(epos_ptr[1],epos_ptr[2]);
                    double dam_old = edam_ptr[0];
                    double tdamthird = pow(dam_old,1.0/3.);
                    tdamthird += dam_inc;
                    edam_ptr[0] = Wind(pow(tdamthird,3.0),0.,1.);
                    double sfac = dam_old < 0.98 ? (1.0 - edam_ptr[0])/(1.0 - dam_old) : 0.0;
                    vec_scaler(este_ptr,4,sfac);
                    este2s_ptr[0] *= sfac;
                }
            }
        }

    }

}

void update_strength_condition(sale2d_var * _sale)
{
    if(_sale->tracer_lubrication <= 0) return;
    // detect some tracers that reduce the strength of cells,
    // write results in e_cnd, <conditions of cells>
    clean_fv(_sale->e_cnd); // set all condition to default(0)
    field_list_iterator * fit = fl_iterator_new(_sale->tracer_list,field_list_head);
    field_list_node * node = NULL;

    // loop on tracers, just read,
    // any writing of tracers would change the node_state to dirty/delete,
    // dirty/delete state should be processed by fl_check
    while(NULL != (node = fl_iterator_next(fit)))
    {
        if(node->state != field_list_used) continue;
        int k = node->id.x, j = node->id.y;
        fv_id_2d eid = {k, j};
        double *ecnd_ptr = NULL;
        _sale->foe->get_double(&eid,&ecnd_ptr,_sale->e_cnd);
        if(node->data.gy < _sale->tracer_lub_gy1 && node->data.gy >= _sale->tracer_lub_gy0)
        {
            *ecnd_ptr += 1.0;
        }
    }
    fl_iterator_del(fit);
}


double JohnsonCook2(double *_paras, const state_reference *_s)
{
    double pres = _paras[0];
    double PlasticStrain    = Wind(_paras[2],0.0,_s->JCMaxPlastic);
    double DotPlasticStrain = _paras[3];
    double temp = _paras[4];

    double JcY1 = _s->JCA + _s->JCB * pow(PlasticStrain,_s->JCN);
    double JcY2 = 1.0;
    if(DotPlasticStrain > 1.0 && _s->JCC > 0.0)
    {
        JcY2 = 1.0 + _s->JCC * log(DotPlasticStrain);
    }

    return JcY1*JcY2*JNCKSoft(temp,pres,_s)*BETA;
}

double JohnsonCook1(double * _paras, const state_reference *_s)
{
    double pres = _paras[0];
    double PlasticStrain    = Wind(_paras[2],0.0,_s->JCMaxPlastic);
    double DotPlasticStrain = _paras[3];

    double JcY1 = _s->JCA + _s->JCB * pow(PlasticStrain,_s->JCN);
    double JcY2 = 1.0;
    if(DotPlasticStrain > 1.0 && _s->JCC > 0.0)
    {
        JcY2 = 1.0 + _s->JCC * log(DotPlasticStrain);
    }
    return JcY1*JcY2*BETA;
}


double SimpleRock(double * _paras, const state_reference * _s)
{
    double pre = _paras[0];
    double dam = _paras[1];

    /*
     * this is a example for the strength function
     * This code encourage user to write the Strength function
     * accord to different materials.
     * The parameters required by this function for different materials can be added
     * in the struct StateReference.
     */

    // Intact strength uses the lundborg approximation
    double Yint = Ylundborg(pre,_s->Yint0,_s->Yfricint,_s->Ylimint);
    double Ydam = Ydrucker(pre,_s->Ydam0,_s->Yfricdam,_s->Ylimdam);
    Ydam = Min(Yint,Ydam);
    return Ydam * dam + Yint * (1.0 - dam);
}

double VonMises1(double * _paras, const state_reference * _s)
{
    double pre = _paras[0];
    double dam = _paras[1];

    /*
     * the von mise strength model, (constant)
     */

    // Intact strength uses the lundborg approximation
    double Yint = _s->Yint0;
    double Ydam = _s->Ydam0;
    Ydam = Min(Yint,Ydam);
    return Ydam * dam + Yint * (1.0 - dam);
}

double VonMises2(double * _paras, const state_reference * _s)
{
    // von mises strength with onnakasoft
    double pre = _paras[0];
    double dam = _paras[1];
    double PlasticStrain    = Wind(_paras[2],0.0,_s->JCMaxPlastic);
    double DotPlasticStrain = _paras[3];
    double temp = _paras[4];

    // Intact/Damaged strength as const
    double Yint = _s->Yint0;
    double Ydam = _s->Ydam0;
    Ydam = Min(Yint,Ydam);

    double melt_temp = SimonMelt(pre,_s);
    double soft_factor = OhnakaSoft(temp, _s, melt_temp);
    soft_factor = Wind(soft_factor,0.0,1.0);
    return (Ydam * dam + Yint * (1.0 - dam))*soft_factor;
}



double IvanovShearDam(double *_paras, const state_reference *_s)
{
    /*
     * this function return the increment of damage caused by shear failure
     * This function return a reverse of the counterpart of iSALE2d
     * Example for IVANOV damage model
     */
    double pre = _paras[0];
    return 1.0/Max(_s->IvanA, _s->IvanB * (pre - _s->IvanC));
}

double CollinsShearDam(double * _paras, const state_reference * _s)
{
    /*
     * calculate the increment of damage
     * adopt from Collins 2004
     */
    double pre = _paras[0];
    double dam = _paras[1];
    double dam_inc = 0.;

    if(pre < _s->Pbd)
    {
        dam_inc = Max(0.01 , 0.01 + 0.04 * pre/_s->Pbd);
    } else if(pre < _s->Pbp)
    {
        dam_inc = 0.05 + 0.05 * (pre - _s->Pbd)/(_s->Pbp - _s->Pbd);
    } else
    {
        dam_inc = 0.1 + 0.5 *(pre - _s->Pbp)/_s->Pbp;
    }

    return 1.0/dam_inc;
}

double Ylundborg(double p,double y0,double fric,double ylim)
{
    /*
     * Smooth lundborg function from y0 to ylim
     */
    return y0 + p*fric/(1 + p*fric/(ylim-y0));
}

double Ydrucker(double p,double y0,double fric,double ylim)
{
    /*
     * Linear with pressure from y0 up to ylim
     */
    return Min(y0 + fric*p,ylim);
}

double JNCKSoft(double tem,double pre, const state_reference *_s)
{
    double Tmelt = SimonMelt(pre,_s);
    double MeltFrac = (tem - _s->JCTREF)/(Tmelt - _s->JCTREF);
    MeltFrac = Wind(MeltFrac,0.0,1.0);
    return Wind(1.0 - pow(MeltFrac,_s->JCM),0.0,1.0);
}

double SimonMelt0(double pre, const state_reference *_s)
{
    return _s->Tmelt0 * pow(_s->Asimon*pre + 1.0,_s->Csimon);
}

double PolyLiquids(double pre, const state_reference * _s)
{
    double p_gpga = pre * 1.0e-9;
    double T_liquids = 0.0;
    T_liquids = _s->Tlidc1 + _s->Tlidc2*p_gpga + _s->Tlidc3*p_gpga*p_gpga;
    return T_liquids;
}

double PolySolidus(double pre, const state_reference * _s)
{
    double p_gpga = pre * 1.0e-9;
    double T_solidus = _s->Tsolc1 + _s->Tsolc2*p_gpga + _s->Tsolc3*p_gpga*p_gpga;
    return T_solidus;
}

double SimonMelt(double pre, const state_reference *_s)
{
    // compare the old version SimonMelt;
    // this version will return the liquid temperature if necessary
    double solidus_temp = SimonMelt0(pre,_s);
    if(_s->PolyLiq <= 0 && _s->PolySol <= 0) return solidus_temp;
    if( _s->PolySol > 0) solidus_temp = PolySolidus(pre,_s);
    double liquids_temp  = PolyLiquids(pre,_s);
    return Max(solidus_temp,liquids_temp);
}

double SolidusTem(double pre, const state_reference * _s) {
    // try use simon approximatin
    double solidus_temp = SimonMelt0(pre,_s);
    if( _s->PolySol > 0) solidus_temp = PolySolidus(pre,_s);
    return solidus_temp;
}

double OhnakaSoft(double tem, const state_reference *_s, double Tmelt)
{
    return tanh(((Tmelt+_s->Tdelta)/tem - 1.0)*_s->Tfrac);
}

double LowdensitySoft(double Density, double MeltDensity)
{
    return pow(Density/MeltDensity,4.0);
}


strength_func str2func(const char s[])
{
    int cur_func = -1;
    for(int k=0;str_func_list[k].f != NULL; ++k)
    {
        if(strcasecmp(s,str_func_list[k].name) == 0)
        {
            cur_func = k;
            break;
        }
    }

    if(cur_func >0)
    {
        return str_func_list[cur_func].f;
    } else
    {
        fprintf(stderr,"undefined function: %s \n", s);
        exit(0);
        return NULL;
    }
}


void select_strength_model(sale2d_var * _sale)
{
    /*
     * set strength/damge model function in state reference
     * according name loaded in yfunc/dSDf
     */

    for(int matid=1;matid<_sale->nmat;++matid)
    {
        // select yield strength
        _sale->etb[matid].ref->shear_strength_func =(uintptr_t) str2func(_sale->etb[matid].ref->yfunc);

        // select damage increment
        _sale->etb[matid].ref->damage_increment_func = (uintptr_t) str2func(_sale->etb[matid].ref->dSDf);

        // select tensile strength
        _sale->etb[matid].ref->tensile_strength_func = (uintptr_t) str2func(_sale->etb[matid].ref->tfunc);

        // porosity setting
        _sale->etb[matid].ref->porosity_increment_func = (uintptr_t) str2func(_sale->etb[matid].ref->porosityfunc);
    }
}

void block_vibration(sale2d_var * _sale, node_type _nt,double dt,double time)
{
    /*
     * e_vib is the energy of element vibration,
     * this variable require flux between cells, so the send/internal section should be updated separately
     */

    if(_sale->acoustic_fluid == 0) return;

    field_opt * _opt = _sale->foe;
    field_var * _var = _sale->e_vib;

    proc_info_2d * _proc = (proc_info_2d *) _var->proc_self;
    int nx = _proc->nx, ny = _proc->ny;
    int nmat = _sale->nmat;

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {k,j};
            node_type local_node_type = useless_section;
            _opt->get_nty(&lid,&local_node_type,_var);
            if(!(_nt&local_node_type)) continue;

            double *evof_ptr, *edam_ptr, *evib_ptr;
            _opt->get_double(&lid,&evof_ptr,_sale->e_vof);
            _opt->get_double(&lid,&edam_ptr,_sale->e_dam);
            _opt->get_double(&lid,&evib_ptr,_sale->e_vib);

            if(evof_ptr[VACUUM] > TOLVOF || edam_ptr[0] < 0.5)
            {
                evib_ptr[0] = 0.;
                continue;
            }

            double Tdeacy = -1.0, Cvib=0.0, Toff=-1.0, VibMax=0.0, Pvlim=0.0;
            for(int matid=1;matid<nmat;++matid)
            {
                if(evof_ptr[matid] < TOLVOF) continue;
                state_reference * _sref = _sale->etb[matid].ref;
                Tdeacy = Max(_sref->Tdamp, Tdeacy);
                Cvib   = Max(_sref->Cvib ,   Cvib);
                Toff   = Max(_sref->Toff ,  Toff );
                VibMax = Max(_sref->VibMax,VibMax);
                Pvlim  = Max(_sref->Pvlim,  Pvlim);
            }

            if(Tdeacy <=0. || Toff <= 0.)
            {
                evib_ptr[0] = 0.;
                continue;
            }

            double * evel_ptr, *eq_ptr;
            _opt->get_double(&lid,&evel_ptr,_sale->e_vel);
            _opt->get_double(&lid,&eq_ptr,_sale->e_q);

            double new_vib = Cvib * vec_len(evel_ptr,2);
            double new_vib2 = new_vib*new_vib;

            if( (time <= Toff) && (new_vib2 >= evib_ptr[0]) && (new_vib >= 1.0) &&  (eq_ptr[0] > 0.))
            {
                evib_ptr[0] = Min(new_vib2,VibMax*VibMax);
            } else
            {
                evib_ptr[0] = Max(evib_ptr[0]*(1.0 - 2.0*dt/Tdeacy + 2.0*(dt/Tdeacy)*(dt/Tdeacy)),0.0);
            }
        }
    }
}

double RockStrengthT(double * _paras, const state_reference * _s)
{
    // similar like SimpleRock
    double pre = _paras[0];
    double dam = _paras[1];
    double PlasticStrain    = Wind(_paras[2],0.0,_s->JCMaxPlastic);
    double DotPlasticStrain = _paras[3];
    double temp = _paras[4];
    /*
     * adopt from Simple Rock, add thermal soft model
     */

    // Intact strength uses the lundborg approximation
    double Yint = Ylundborg(pre,_s->Yint0,_s->Yfricint,_s->Ylimint);
    double Ydam = Ydrucker(pre,_s->Ydam0,_s->Yfricdam,_s->Ylimdam);
    Ydam = Min(Yint,Ydam);

    double melt_temp = SimonMelt(pre,_s);
    double soft_factor = OhnakaSoft(temp, _s, melt_temp);

    soft_factor = Wind(soft_factor,0.0,1.0);

    double ymed = (Ydam * dam + Yint * (1.0 - dam))*soft_factor;
    double yvis = DotPlasticStrain*_s->melt_viscosity*(1.0 - soft_factor);
    return  Max(ymed,yvis);
}

double IceMelt(double pressure) {
    /*
     * melt temperature for ice
     * adopt from isale2d
     * Datchi, F., Louberyre, P., and LeToullec, R. (2000)
     * "Extended and accurate determination of the melting curves of
     * argon, helium, ice (H2O) and hydrogen (H2)". Phys. Rev. B
     * Vol. 61, #10, pp. 6535-6546
     */
    const double pice0 =    0.0, pice1 = 0.2e9, pice2 = 0.34e9;
    const double tice0 =  273.0, tice1 = 252.0, tice2 = 257.0;
    const double pice3 = 0.61e9, pice4 = 2.17e9;
    const double tice3 =  274.0, tice4 = 354.8;

    // Calculate melt temp as piecewise series of functions of pressure...
    if (pressure <= pice1) {
        return tice0 + (tice1 - tice0) * (pressure - pice0) / (pice1 - pice0);
    } else if (pressure > pice1 && pressure <= pice2) {
        return tice1 + (tice2 - tice1) * (pressure - pice1) / (pice2 - pice1);
    } else if (pressure > pice2 && pressure <= pice3) {
        return tice2 + (tice3 - tice2) * (pressure - pice2) / (pice3 - pice2);
    } else if (pressure > pice3 && pressure <= pice4) {
        return tice3 + (tice4 - tice3) * (pressure - pice3) / (pice4 - pice3);
    } else {
        // High-pressure trend is a Simon approximation (as for rock)
        return tice4 * pow(((pressure - pice4) / 1.253e9) + 1, 1.0 / 3.0);
    }
}



double IceStrength(double * _paras, const state_reference * _s)
{
    double pre = _paras[0];
    double dam = _paras[1];
    double PlasticStrain    = Wind(_paras[2],0.0,_s->JCMaxPlastic);
    double DotPlasticStrain = _paras[3];
    double temp = _paras[4];

    // Intact/Damaged strength uses the lundborg approximation
    double Yint = Ylundborg(pre,_s->Yint0,_s->Yfricint,_s->Ylimint);
    double Ydam = Ylundborg(pre,_s->Ydam0,_s->Yfricdam,_s->Ylimdam);
    Ydam = Min(Yint,Ydam);

    double melt_temp = IceMelt(pre);
    double soft_factor = OhnakaSoft(temp, _s, melt_temp);
    soft_factor = Wind(soft_factor,0.0,1.0);
    return (Ydam * dam + Yint * (1.0 - dam))*soft_factor;
}



double SimpleTensile(double * _paras, const state_reference * _s)
{
    double pre = _paras[0];
    double dam = _paras[1];
    double melt_temp = SimonMelt(pre,_s);
    double temp = _paras[4];
    double soft_factor = OhnakaSoft(temp, _s, melt_temp);
    soft_factor = Wind(soft_factor,0.0,1.0);
    return _s->Yten0 * (1.0 - dam) * soft_factor;
}


double ZeroFunc(double *_paras, const state_reference *_s)
{
    return 0.0;
}


double PorosityCsdRatio(double alpha, double alpha0, double chi)
{
    // return the sound speed ratio,
    // the porous material at zero pressure to the soild matrix
    return 1.0 + (chi - 1.0) * fmax(0.0, fmin(1.0, (alpha - 1.0) / (alpha0 - 1.0)));
}

double WunnemannPorosity(double * _paras, const state_reference * _s)
{
    // return the \frac{d\alpha}{d t} term
    double pty = _paras[0]; // distension, alpha
    double cvs = _paras[1]; // cold volume strain
    double pre = _paras[2];
    double dcvs = _paras[3]; // time derivation of cold volume strain
    double evx  = _paras[4];
    double evy  = _paras[5];

    // skip some compression, or over tension cells ?
    if(pty <= 1.0 || pty > _s->alphamax || pre <= 0.) return 0.;
    double evel2 = evx*evx + evy*evy;
    if( evel2 <= _s->avel*_s->avel) return 0.;

    // the critical volume strain, porous materials are divied to
    // compaction at low strain and compression regime
    double cvs_c = 0.;
    if(pty > _s->alphae)
    {
        cvs_c = _s->epse0;
    }
    else if(pty > _s->alphax)
    {
        cvs_c = _s->epse0 + log(pty/_s->alphae)/_s->kappa;
    }
    else
    {
        cvs_c = _s->epsec - sqrt((pty - 1.0)/(_s->alphax - 1.0))*(_s->epsec - _s->epsex);
    }

    // the \frac{d\alpha}{d\epsilon} term
    double CsdRatio = PorosityCsdRatio(pty,_s->alphamax,_s->chi);
    double dade = pty*(1.0 - CsdRatio*CsdRatio);
    if(dcvs < 0. && cvs < cvs_c)
    {
        // compression
        if(pty > _s->alphax)
        {
            dade = _s->kappa * pty;
        }
        else
        {
            dade = sqrt((pty - 1.)/(_s->alphax - 1.0)) * _s->kappa * _s->alphax;
        }
    }
    else
    {
        // compaction within low strain
        dade = pty*(1.0 - CsdRatio*CsdRatio);
    }

    return dade*dcvs;
}


/*
 * register all the strength related function
 */
str_func_pair str_func_list[64] = {
        {.name = "SimpleRock",.f = SimpleRock},
        {.name = "SimpleShear",.f = IvanovShearDam},
        {.name = "JohnsonCook1", .f = JohnsonCook1},
        {.name = "JohnsonCook2", .f = JohnsonCook2},
        {.name = "RockStrengthT", .f = RockStrengthT},
        {.name = "IceStrength", .f = IceStrength},
        {.name = "Collins", .f = CollinsShearDam},
        {.name = "Ivanov", .f= IvanovShearDam},
        {.name = "SimpleTensile",.f=SimpleTensile},
        {.name = "None", .f= ZeroFunc},
        {.name = "Wunnemann", .f= WunnemannPorosity},
        {.name = "VonMises1", .f= VonMises1},
        {.name = "VonMises2", .f= VonMises2},
        {.name = "xyz", .f=NULL}
};