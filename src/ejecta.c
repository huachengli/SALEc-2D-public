//
// Created by huach on 12/7/2023.
//

#include "ejecta.h"


void trace_profile(sale2d_var * _sale)
{
    /*
     * trace the surface of ejecta,
     * the materials at (r,z) where z > threshold is marked as ejecta materials
     * the position is recorded in y_profile
     */

    field_var * _var = _sale->e_pos;
    field_opt * _opt = _sale->foe;

    proc_info_2d * _proc = (proc_info_2d *) _var->proc_self;
    int nx = _proc->nx, ny = _proc->ny;

    double threshold_z = _sale->threshold_z;

    // before find the profile position, reset the profile of previous step
    for(int k=0;k<_sale->y_profile->data_size;++k)
    {
        double * yprof_ptr = (double *)(_sale->y_profile->data + _sale->y_profile->bytes_per_unit * k);
        *yprof_ptr = threshold_z*2.0;
    }

    for(int index=0;index<nx*ny;++index)
    {
        int k = INDEX2LIDX(index, nx, ny);
        int j = INDEX2LIDY(index, nx, ny);
        fv_id_2d lid = {k, j};
        node_type lnt = useless_section;
        _opt->get_nty(&lid, &lnt, _var);

        if(useless_section == lnt) continue;

        double * evof_ptr;
        _opt->get_double(&lid,&evof_ptr,_sale->e_vof);
        double * eejt_ptr;
        double * epos_ptr;
        _opt->get_double(&lid,&eejt_ptr,_sale->e_ejt);
        _opt->get_double(&lid,&epos_ptr,_sale->e_pos);

        double * yprof_ptr;
        _sale->foi->get_double(&lid,&yprof_ptr,_sale->y_profile);

        // find element have vacuum, check whether this is below profile
        if(evof_ptr[VACUUM] > 0.1)
        {
            *yprof_ptr = MIN(epos_ptr[Y],*yprof_ptr);
        }

        // for element z> threshold_z, set the ejecta mass fraction is 1
        if(epos_ptr[Y] > threshold_z && evof_ptr[VACUUM] < TOLVOF)
        {
            *eejt_ptr = 1.0;
        }
    }
}

void ejecta_flow(sale2d_var * _sale, node_type _nt,double dt)
{

    // negative friction is invalid turn off resistance flow of ejecta
    if(_sale->flow_resistance <= 0. && _sale->flow_friction <= 0.) return;

    field_var * _var = _sale->e_pos;
    field_opt * _opt = _sale->foe;

    proc_info_2d * _proc = (proc_info_2d *) _var->proc_self;
    int nx = _proc->nx, ny = _proc->ny;

    for(int index=0;index<nx*ny;++index)
    {
        int k = INDEX2LIDX(index, nx, ny);
        int j = INDEX2LIDY(index, nx, ny);
        fv_id_2d lid = {k, j};
        fv_id_2d vlid[2] = {{k - 1, j - 1},
                            {k,     j - 1}};
        node_type lnt = useless_section;
        _opt->get_nty(&lid, &lnt, _var);

        if(useless_section & lnt) continue;
        if(!(lnt&_nt)) continue;

        double *eejt_ptr;
        double *yprof_ptr;
        double *epos_ptr;
        double *dvel_ptr;
        double *evel_ptr;
        double *este_ptr;
        _sale->foe->get_double(&lid, &eejt_ptr, _sale->e_ejt);
        _sale->foe->get_double(&lid, &epos_ptr, _sale->e_pos);
        _sale->foi->get_double(&lid, &yprof_ptr, _sale->y_profile);

        double egravity_ptr[2] = {0., 0.};
        double *vbf_ptr[2];
        _sale->fov->get_double(vlid + 0, &(vbf_ptr[0]), _sale->v_bf);
        _sale->fov->get_double(vlid + 1, &(vbf_ptr[1]), _sale->v_bf);
        if(NULL == vbf_ptr[0] || NULL == vbf_ptr[1])
        {
            fprintf(stderr,"%d,%d, => null vbf\n",k,j);
        }

        scaler_add(egravity_ptr,2,vbf_ptr[0],0.5);
        scaler_add(egravity_ptr,2,vbf_ptr[1],0.5);


        if(epos_ptr[Y] <= yprof_ptr[0] && eejt_ptr[0] > 0.1)
        {
            _sale->foe->get_double(&lid,&dvel_ptr,_sale->e_dvel);
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);
            _sale->foe->get_double(&lid,&este_ptr,_sale->e_ste);

            if(_sale->flow_resistance > 0.)
                scaler_add(dvel_ptr,2,evel_ptr, -_sale->flow_resistance * dt);

            if(_sale->flow_friction > 0.)
            {
                vec_zero(este_ptr,4);
                if(evel_ptr[X] >= 0.)
                    dvel_ptr[X] += -1.0* _sale->flow_friction* fabs(egravity_ptr[Y])*dt;
                else
                    dvel_ptr[X] +=  1.0* _sale->flow_friction* fabs(egravity_ptr[Y])*dt;
            }

        }
    }

}