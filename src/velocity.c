/*
 * functions for calculate acceleration/displacement
 * damping on vertex is also included
 */

#include "sale2d.h"


void calculate_acc_default(sale2d_var * _sale, node_type _nt, double dt)
{
    /*
     * calculate the acceleration of vertex
     * VARIABLES need to be synchronized in advance
     *      e_den, e_pre
     *      e_stress
     * VARIABLES overwritten AFTER
     *      v_acc
     */
    proc_info_2d* proc_acc = (proc_info_2d*) _sale->v_acc->proc_self;
    int nx = proc_acc->nx,ny=proc_acc->ny,myrank=_sale->v_acc->rank;

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j,0,0};
            node_type lnt = useless_section;
            _sale->fov->get_nty(&lid,&lnt,_sale->v_acc);

            if(!(lnt&_nt)) continue;

            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};
            double * acc_ptr=NULL, *bf_ptr=NULL;
            _sale->fov->get_double(&lid,&acc_ptr,_sale->v_acc);
            _sale->fov->get_double(&lid,&bf_ptr,_sale->v_bf);
            vec_zero(acc_ptr,2);

            double *epre_ptr[4], *eden_ptr[4];
            int v_offset[4] = {2, 3, 0, 1};
            double v_mass = 0.0, v_vol = 0.0;
            double * estress_ptr[4], *eq_ptr[4];

            for(int i=0;i<4;++i)
            {
                // physical variables
                _sale->foe->get_double(elid[i],&epre_ptr[i],_sale->e_pre);
                _sale->foe->get_double(elid[i],&eden_ptr[i],_sale->e_den);
                _sale->foe->get_double(elid[i],&estress_ptr[i],_sale->e_ste);
                _sale->foe->get_double(elid[i],&eq_ptr[i],_sale->e_q);
            }

            // add subvol of vertex, geometrical parameters
            double *vsubvol_ptr,*vvol_ptr;
            _sale->fov->get_double(&lid,&vsubvol_ptr,_sale->v_subvol);
            _sale->fov->get_double(&lid,&vvol_ptr,_sale->v_vol);
            double *vsuf_ptr;
            _sale->fov->get_double(&lid,&vsuf_ptr,_sale->v_suf);

            for(int i=0;i<4;i++)
            {
                /*
                 * note esuf_ptr[] is the surface at neighbour element
                 * direction is inward
                 */

                // acceleration for pressure
                double f_pre[2] = {0.};
                double local_pressure = epre_ptr[i][0] + eq_ptr[i][0]; // pressure + artificial pressure
                scaler_add(f_pre,2,vsuf_ptr + 2*i,local_pressure);
                f_pre[0] += vsubvol_ptr[i] * local_pressure;

                // acceleration for stress
                double Sxx = estress_ptr[i][0], Szz = estress_ptr[i][2], Sxz = estress_ptr[i][1], Sth = estress_ptr[i][3];
                double f_stress[2] = {0.};

                // the integrate term on (x,z)
                f_stress[0] = vec_dot(vsuf_ptr + 2*i, (double []){Sxx, Sxz}, 2);
                f_stress[1] = vec_dot(vsuf_ptr + 2*i, (double []){Sxz, Szz}, 2);

                // the term related to the cylindrical coordinate
                f_stress[0] += vsubvol_ptr[i] * Sth;

                // notice, the direction of surface vector is point into the vertex
                double local_vertex_vol = vvol_ptr[1+i];
                scaler_add(acc_ptr,2,f_stress,-1.0);
                scaler_add(acc_ptr,2,f_pre,1.0);

                v_mass += local_vertex_vol * eden_ptr[i][0];
                v_vol  += local_vertex_vol;
            }

            // calculate the dissipation of velocity (ANC)
            // v_acc = anc*(sound speed)*sqrt(vol)*Laplacian(vec)
            fv_id_2d vlid[9] = {
                    {k-1,j-1},
                    {k  ,j-1},
                    {k+1,j-1},
                    {k-1,j},
                    {k  ,j},
                    {k+1,j},
                    {k-1,j+1},
                    {k  ,j+1},
                    {k+1,j+1}
            };

            double tmp_vx[9] = {0.};
            double tmp_vy[9] = {0.};
            for(int kv=0;kv<9;++kv)
            {
                double * tmp_vel = NULL;
                _sale->fov->get_double(vlid+kv,&tmp_vel,_sale->v_vel);
                tmp_vx[kv] = tmp_vel[0];
                tmp_vy[kv] = tmp_vel[1];
            }
            double *vlaplace_ptr;
            _sale->fov->get_double(&lid,&vlaplace_ptr,_sale->v_laplace);

            double acc_anc[2] = {0.};

            if(lnt & pad_section){
                vec_zero(acc_anc, 2);
            }
            else
            {
                double csd_approx = 5.0e3;
                double Lvx = vec_dot(vlaplace_ptr,tmp_vx,9);
                double Lvy = vec_dot(vlaplace_ptr,tmp_vy,9);
                acc_anc[X] = csd_approx*sqrt(fabs(vsubvol_ptr[0]) + fabs(vsubvol_ptr[2]))*Lvx;
                acc_anc[Y] = csd_approx*sqrt(fabs(vsubvol_ptr[0]) + fabs(vsubvol_ptr[2]))*Lvy;
            }

            double *vmass_ptr;
            _sale->fov->get_double(&lid,&vmass_ptr,_sale->v_mass);
            // vmass_ptr[1] is the void corners around this vertex
            if(vmass_ptr[1] >= 1.) vec_zero(acc_anc, 2);

            if(v_mass > TOLRHO * v_vol)
            {
                vec_scaler(acc_ptr,2,1.0/v_mass);
                scaler_add(acc_ptr,2,bf_ptr,1.0);
                scaler_add(acc_ptr,2,acc_anc,_sale->anc);
            } else
            {
                vec_zero(acc_ptr,2);
            }

            // at last remove small accelerations
#ifdef ACCMIN_CUTOFF
            if(fabs(acc_ptr[X]) < _sale->acceleration_min*ACCMIN_CUTOFF) acc_ptr[X] = 0.;
            if(fabs(acc_ptr[Y]) < _sale->acceleration_min*ACCMIN_CUTOFF) acc_ptr[Y] = 0.;
#endif
        }
    }
}


void calculate_acc_pml(sale2d_var * _sale, node_type _nt, double dt)
{
    /*
     * calculate the acceleration of vertex
     * VARIABLES need to be synchronized in advance
     *      e_den, e_pre, e_q
     *      e_stress
     * VARIABLES overwritten AFTER
     *      v_acc, v_alpha
     *  overwrite args has been deleted.
     */

    proc_info_2d* proc_acc = (proc_info_2d*) _sale->v_acc->proc_self;
    int nx = proc_acc->nx,ny=proc_acc->ny,myrank=_sale->v_acc->rank;

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j,0,0};
            node_type lnt = useless_section;
            _sale->fov->get_nty(&lid,&lnt,_sale->v_acc);

            if(!(lnt&_nt)) continue;

            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};
            double * acc_ptr=NULL, *bf_ptr=NULL, *vel_ptr;
            _sale->fov->get_double(&lid,&acc_ptr,_sale->v_acc);
            _sale->fov->get_double(&lid,&bf_ptr,_sale->v_bf);
            _sale->fov->get_double(&lid,&vel_ptr,_sale->v_vel);
            vec_zero(acc_ptr,2); //before used of acc, clean its values


            double *epre_ptr[4], *eden_ptr[4];
            int v_offset[4] = {2, 3, 0, 1};
            double v_mass = 0.0, v_vol = 0.0;
            double * estress_ptr[4], *eq_ptr[4], *esigma_ptr[4], *edampref_ptr[4];

            for(int i=0;i<4;++i)
            {
                // physical variables
                _sale->foe->get_double(elid[i],&epre_ptr[i],_sale->e_pre);
                _sale->foe->get_double(elid[i],&eden_ptr[i],_sale->e_den);
                _sale->foe->get_double(elid[i],&estress_ptr[i],_sale->e_ste);
                _sale->foe->get_double(elid[i],&eq_ptr[i],_sale->e_q);
                // damp factor
                _sale->foe->get_double(elid[i],&esigma_ptr[i],_sale->e_sigma);
                _sale->foe->get_double(elid[i],&edampref_ptr[i],_sale->e_dampref);
            }

            // add subvol of vertex, geometrical parameters
            double *vsubvol_ptr,*vvol_ptr;
            _sale->fov->get_double(&lid,&vsubvol_ptr,_sale->v_subvol);
            _sale->fov->get_double(&lid,&vvol_ptr,_sale->v_vol);
            double *vsuf_ptr;
            _sale->fov->get_double(&lid,&vsuf_ptr,_sale->v_suf);
            double * valpha_ptr;
            _sale->fov->get_double(&lid,&valpha_ptr,_sale->v_alpha);

            // calculate damp factor on cells
            double vsigma[2] = {0.};
            for(int i=0;i<4;++i)
            {
                scaler_add(vsigma,2,esigma_ptr[i],0.25);
            }


            for(int i=0;i<4;i++)
            {
                /*
                 * note esuf_ptr[] is the surface at neighbour element
                 * direction is inward
                 */

                // acceleration for pressure
                double lp = epre_ptr[i][0] + eq_ptr[i][0]; // pressure + artificial pressure
                double lpr = edampref_ptr[i][2];
                // acceleration for stress
                double Sxx = estress_ptr[i][0] - lp, Szz = estress_ptr[i][2] - lp, Sxz = estress_ptr[i][1], Sth = estress_ptr[i][3] - lp;
                double f_stress[2] = {0.};

                // the integrate term on (x,z)
                f_stress[0] = vec_dot(vsuf_ptr + 2*i, (double []){Sxx, Sxz}, 2);
                f_stress[1] = vec_dot(vsuf_ptr + 2*i, (double []){Sxz, Szz}, 2);

                // the term related to the cylindrical coordinate
                f_stress[0] += vsubvol_ptr[i] * Sth;
                scaler_add(acc_ptr,2,f_stress,-1.0);


                if(vsigma[X] > 0. || vsigma[Y] > 0.)
                {
                    Sxx += lpr;
                    Szz += lpr;
                    Sth += lpr;

                    valpha_ptr[X] += vec_dot(vsuf_ptr + 2*i, (double []){vsigma[Y]*Sxx, vsigma[X]*Sxz}, 2)*dt;
                    valpha_ptr[Y] += vec_dot(vsuf_ptr + 2*i, (double []){vsigma[Y]*Sxz, vsigma[X]*Szz}, 2)*dt;
                    valpha_ptr[X] += (vsigma[X] + vsigma[Y])* vsubvol_ptr[i] * Sth * dt;
                }

                // notice, the direction of surface vector is point into the vertex
                double local_vertex_vol = vvol_ptr[1+i];
                v_mass += local_vertex_vol * eden_ptr[i][0];
                v_vol  += local_vertex_vol;
            }

            if(vsigma[X] > 0. || vsigma[Y] > 0.)
            {
                // remove dependence of damping rate on sigma
                valpha_ptr[X] -= _sale->alpha_damp_x * valpha_ptr[X] * dt;
                valpha_ptr[Y] -= _sale->alpha_damp_y * valpha_ptr[Y] * dt;
                scaler_add(acc_ptr,2,valpha_ptr,-1.0);
                // the term: (dx+dy)*V*dt
                // scaler_add(acc_ptr,2,vel_ptr,-1.0*(vsigma[X]+vsigma[Y])*v_mass);
            }


            // calculate the dissipation of velocity (ANC)
            // v_acc = anc*(sound speed)*sqrt(vol)*Laplacian(vec)
            fv_id_2d vlid[9] = {
                    {k-1,j-1},
                    {k  ,j-1},
                    {k+1,j-1},
                    {k-1,j},
                    {k  ,j},
                    {k+1,j},
                    {k-1,j+1},
                    {k  ,j+1},
                    {k+1,j+1}
            };

            double tmp_vx[9] = {0.};
            double tmp_vy[9] = {0.};
            for(int kv=0;kv<9;++kv)
            {
                double * tmp_vel = NULL;
                _sale->fov->get_double(vlid+kv,&tmp_vel,_sale->v_vel);
                tmp_vx[kv] = tmp_vel[0];
                tmp_vy[kv] = tmp_vel[1];
            }
            double *vlaplace_ptr;
            _sale->fov->get_double(&lid,&vlaplace_ptr,_sale->v_laplace);

            double acc_anc[2] = {0.};


            if(lnt & pad_section){
                vec_zero(acc_anc, 2);
            }
            else
            {
                double csd_approx = 5.0e3;
                double Lvx = vec_dot(vlaplace_ptr,tmp_vx,9);
                double Lvy = vec_dot(vlaplace_ptr,tmp_vy,9);
                acc_anc[X] = csd_approx*sqrt(fabs(vsubvol_ptr[0]) + fabs(vsubvol_ptr[2]))*Lvx;
                acc_anc[Y] = csd_approx*sqrt(fabs(vsubvol_ptr[0]) + fabs(vsubvol_ptr[2]))*Lvy;
            }

            double *vmass_ptr;
            _sale->fov->get_double(&lid,&vmass_ptr,_sale->v_mass);
            // vmass_ptr[1] is the void corners around this vertex
            if(vmass_ptr[1] >= 1.) vec_zero(acc_anc, 2);

            if(v_mass > TOLRHO * v_vol)
            {
                vec_scaler(acc_ptr,2,1.0/v_mass);
                scaler_add(acc_ptr,2,bf_ptr,1.0);
                scaler_add(acc_ptr,2,acc_anc,_sale->anc);
            } else
            {
                vec_copy(bf_ptr,acc_ptr,2);
            }


#ifdef ACCMIN_CUTOFF
            if(fabs(acc_ptr[X]) < _sale->acceleration_min*ACCMIN_CUTOFF) acc_ptr[X] = 0.;
            if(fabs(acc_ptr[Y]) < _sale->acceleration_min*ACCMIN_CUTOFF) acc_ptr[Y] = 0.;
#endif
        }
    }
}

void calculate_acc(damp_type_t dampType,sale2d_var * _sale, node_type _nt, double  dt)
{
    switch(dampType)
    {
        case pml3:
        case pml2:
            calculate_acc_pml(_sale,_nt,dt);
            break;
        case pml:
        default:
            calculate_acc_default(_sale,_nt,dt);
    }
}

void damp_on_cell(damp_type_t dampType,sale2d_var * _sale,double dt,double t)
{
    switch(dampType)
    {
        case pml:
            damp_on_cell_pml(_sale,dt,t);
            break;
        case pml2:
            damp_on_cell_pml2(_sale,dt,t);
            break;
        case sbc:
            damp_on_cell_sbc(_sale,dt,t);
            break;
        case pml3:
            damp_on_cell_pml3(_sale,dt,t);
            break;
        case ctlg:
        default:
            break;
    }
}

void damp_on_cell_pml(sale2d_var * _sale,double dt,double t)
{
    /*
     * use perfect match layer algorithm to remove numerical reflection of impact wave
     * the reference stable solution is recorded before cycles
     * at section far from interested, the equation can be simplified to a wave function
     * on density or pressure when the stress/strain term is omitted.
     *
     * ! damp_on_cell_pml only add damping terms related to pressure
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    if(_sale->damp_time <= t)
    {
        return;
    }

    for(int j = 0; j < ny; j++) {
        for (int k = 0; k < nx; ++k) {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if(lnt == useless_section) continue;
            double * evof, *mden, *meng;
            _sale->foe->get_double(&lid,&evof,_sale->e_vof);
            _sale->foe->get_double(&lid,&mden,_sale->m_den);
            _sale->foe->get_double(&lid,&meng,_sale->m_eng);

            /*
             * skip the cell have vacuum
             * absorbing boundary should not be applied to ejector
             */
            if(evof[VACUUM] > TOLVOF) continue;

            double * eden_ptr, *egradvel_ptr,*evel_ptr, *epre_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);

            double * edvel_ptr;
            _sale->foe->get_double(&lid,&edvel_ptr,_sale->e_dvel);

            double *ealpha_ptr,*esigma_ptr;
            _sale->foe->get_double(&lid,&ealpha_ptr,_sale->e_alpha);
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);

            double *edampref_ptr;
            _sale->foe->get_double(&lid,&edampref_ptr,_sale->e_dampref);

            double *mpty_ptr;
            _sale->foe->get_double(&lid,&mpty_ptr,_sale->m_pty);

            double csd_approx = 5000.0;
            // damp on y direction
            if(esigma_ptr[Y] > 0.0)
            {
                // update auxiliary variables
                ealpha_ptr[Y] += (esigma_ptr[Y]*egradvel_ptr[0] + esigma_ptr[Y]*egradvel_ptr[4] - _sale->alpha_damp_y*ealpha_ptr[Y])*dt;
                edvel_ptr[Y]  += -1.0*esigma_ptr[Y]*evel_ptr[Y]*dt;

                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] < 1.0 - TOLVOF) continue;
                    double mden_damp = (-esigma_ptr[Y]*(mden[matid] - edampref_ptr[0]) - edampref_ptr[0]*ealpha_ptr[Y])*dt;
                    double meng_damp = (-esigma_ptr[Y]*(meng[matid] - edampref_ptr[1]) - epre_ptr[0]*ealpha_ptr[Y])*dt;
                    mden[matid] += mden_damp/mpty_ptr[matid];
                    meng[matid] += meng_damp/mpty_ptr[matid];
                    break;
                }
            }

            //damp on X direction
            if(esigma_ptr[X] > 0.0)
            {
                // update auxiliary variables
                if(evof[VACUUM] > TOLVOF)
                {
                    ealpha_ptr[X] += -1.0*_sale->alpha_damp_x*ealpha_ptr[X]*dt;
                    edvel_ptr[X]  += -1.0*esigma_ptr[X]*evel_ptr[X]*dt;
                }
                else
                {
                    ealpha_ptr[X] += (esigma_ptr[X]*egradvel_ptr[3] + esigma_ptr[X]*egradvel_ptr[4] - _sale->alpha_damp_x*ealpha_ptr[X])*dt;
                    edvel_ptr[X] += -1.0*esigma_ptr[X]*evel_ptr[X]*dt;
                }

                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] > 1.0 - TOLVOF)
                    {
                        double mden_damp = (-esigma_ptr[X]*(mden[matid] - edampref_ptr[0]) - edampref_ptr[0]*(ealpha_ptr[X] + ealpha_ptr[Y]))*dt;
                        double meng_damp = (-esigma_ptr[X]*(meng[matid] - edampref_ptr[1]) - epre_ptr[0]*(ealpha_ptr[X]+ealpha_ptr[Y]))*dt;
                        mden[matid] += mden_damp/mpty_ptr[matid];
                        meng[matid] += meng_damp/mpty_ptr[matid];
                        break;
                    }
                }
            }
        }
    }
}



void damp_on_cell_pml2(sale2d_var * _sale,double dt,double t)
{
    /*
     * use perfect match layer algorithm to remove numerical reflection of impact wave
     * the reference stable solution is recorded before cycles
     * at section far from interested, the equation can be simplified to a wave function
     * on density or pressure when the stress/strain term is omitted.
     *
     * add damping term both of pressure and stress
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    if(_sale->damp_time <= t)
    {
        return;
    }

    for(int j = 0; j < ny; j++) {
        for (int k = 0; k < nx; ++k) {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if(lnt & useless_section) continue;
            double * evof, *mden, *meng;
            _sale->foe->get_double(&lid,&evof,_sale->e_vof);
            _sale->foe->get_double(&lid,&mden,_sale->m_den);
            _sale->foe->get_double(&lid,&meng,_sale->m_eng);

            /*
             * skip the cell have vacuum
             * absorbing boundary should not be applied to ejector
             */
            if(evof[VACUUM] > TOLVOF) continue;

            double *eden_ptr, *egradvel_ptr,*evel_ptr, *epre_ptr, *este_ptr, *eq_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);
            _sale->foe->get_double(&lid,&este_ptr,_sale->e_ste);
            _sale->foe->get_double(&lid,&eq_ptr,_sale->e_q);

            double * edvel_ptr;
            _sale->foe->get_double(&lid,&edvel_ptr,_sale->e_dvel);

            double *ealpha_ptr,*esigma_ptr;
            _sale->foe->get_double(&lid,&ealpha_ptr,_sale->e_alpha);
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);

            double *edampref_ptr;
            _sale->foe->get_double(&lid,&edampref_ptr,_sale->e_dampref);
            double *mpty_ptr;
            _sale->foe->get_double(&lid,&mpty_ptr,_sale->m_pty);

            double csd_approx = 5000.0;

            if(esigma_ptr[Y] > 0.0 || esigma_ptr[X] > 0.0)
            {
                // update increment of velocity
                edvel_ptr[X]  += -1.0*(esigma_ptr[X] + esigma_ptr[Y])*evel_ptr[X]*dt;
                edvel_ptr[Y]  += -1.0*(esigma_ptr[X] + esigma_ptr[Y])*evel_ptr[Y]*dt;

                // update auxiliary variables
                if(evof[VACUUM] > TOLVOF)
                {
                    ealpha_ptr[X] += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ealpha_ptr[X]*dt;
                    ealpha_ptr[Y] += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ealpha_ptr[Y]*dt;
                }
                else
                {
                    ealpha_ptr[X] += (esigma_ptr[Y]*egradvel_ptr[0] + esigma_ptr[Y]*egradvel_ptr[4] - _sale->alpha_damp_y*ealpha_ptr[X])*dt
                                     + (esigma_ptr[X]*egradvel_ptr[3] + esigma_ptr[X]*egradvel_ptr[4] - _sale->alpha_damp_x*ealpha_ptr[X])*dt;
                    double local_pre = epre_ptr[0] + eq_ptr[0] - edampref_ptr[2];
                    double ste_damp = (este_ptr[0]-local_pre) * esigma_ptr[Y] *egradvel_ptr[0]
                                    + este_ptr[1]*(esigma_ptr[Y]*egradvel_ptr[1] + esigma_ptr[X]*egradvel_ptr[2])
                                    + (este_ptr[2]-local_pre)*esigma_ptr[X]*egradvel_ptr[3]
                                    + (esigma_ptr[Y] + esigma_ptr[X])*este_ptr[3]*egradvel_ptr[4];
                    ealpha_ptr[Y] += -1.0*ste_damp*dt -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ealpha_ptr[Y]*dt;
                }

                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] < 1.0 - TOLVOF) continue;
                    double mpty0 = _sale->etb[matid].ref->alpha0;
                    double mden_damp = (-(esigma_ptr[X]+esigma_ptr[Y])*(mden[matid] - edampref_ptr[0]) - edampref_ptr[0]*ealpha_ptr[X])*dt;
                    double meng_damp = (-(esigma_ptr[X]+esigma_ptr[Y])*(meng[matid] - edampref_ptr[1]) - ealpha_ptr[Y])*dt;
                    mden[matid] += mden_damp/mpty_ptr[matid];
                    meng[matid] += meng_damp/mpty_ptr[matid];
                    // mden[matid] += mden_damp;
                    // meng[matid] += meng_damp;
                    // mpty_ptr[matid] -= mpty_ptr[matid]*(mden_damp/mden[matid]);
                    break;
                }
            }
        }
    }

}


void damp_on_cell_pml3(sale2d_var * _sale,double dt,double t)
{
    /*
     * use perfect match layer algorithm to remove numerical reflection of impact wave
     * the reference stable solution is recorded before cycles
     * at section far from interested, the equation can be simplified to a wave function
     * on density or pressure when the stress/strain term is omitted.
     *
     * add damping term both of pressure and stress
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    if(_sale->damp_time <= t)
    {
        return;
    }

    for(int j = 0; j < ny; j++) {
        for (int k = 0; k < nx; ++k) {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if(lnt & useless_section) continue;
            double * evof, *mden, *meng;
            _sale->foe->get_double(&lid,&evof,_sale->e_vof);
            _sale->foe->get_double(&lid,&mden,_sale->m_den);
            _sale->foe->get_double(&lid,&meng,_sale->m_eng);

            /*
             * skip the cell have vacuum
             * absorbing boundary should not be applied to ejector
             */
            if(evof[VACUUM] > TOLVOF) continue;

            double *eden_ptr, *egradvel_ptr,*evel_ptr, *epre_ptr, *este_ptr, *eq_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);
            _sale->foe->get_double(&lid,&este_ptr,_sale->e_ste);
            _sale->foe->get_double(&lid,&eq_ptr,_sale->e_q);

            double * edvel_ptr;
            _sale->foe->get_double(&lid,&edvel_ptr,_sale->e_dvel);

            double *ealpha_ptr,*esigma_ptr;
            _sale->foe->get_double(&lid,&ealpha_ptr,_sale->e_alpha);
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);

            double *edampref_ptr;
            _sale->foe->get_double(&lid,&edampref_ptr,_sale->e_dampref);

            double csd_approx = 5000.0;

            if(esigma_ptr[Y] > 0.0 || esigma_ptr[X] > 0.0)
            {
                // update increment of velocity
                edvel_ptr[X]  += -1.0*(esigma_ptr[X] + esigma_ptr[Y])*evel_ptr[X]*dt;
                edvel_ptr[Y]  += -1.0*(esigma_ptr[X] + esigma_ptr[Y])*evel_ptr[Y]*dt;

                // update auxiliary variables
                if(evof[VACUUM] > TOLVOF)
                {
                    ealpha_ptr[X] += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ealpha_ptr[X]*dt;
                    ealpha_ptr[Y] += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ealpha_ptr[Y]*dt;
                }
                else
                {
                    ealpha_ptr[X] += (esigma_ptr[Y]*egradvel_ptr[0] + esigma_ptr[Y]*egradvel_ptr[4] - _sale->alpha_damp_y*ealpha_ptr[X])*dt
                                     + (esigma_ptr[X]*egradvel_ptr[3] + esigma_ptr[X]*egradvel_ptr[4] - _sale->alpha_damp_x*ealpha_ptr[X])*dt;
                    double local_pre = epre_ptr[0] + eq_ptr[0] - edampref_ptr[2];
                    double ste_damp = (este_ptr[0]-local_pre) * esigma_ptr[Y] *egradvel_ptr[0]
                                      + este_ptr[1]*(esigma_ptr[Y]*egradvel_ptr[1] + esigma_ptr[X]*egradvel_ptr[2])
                                      + (este_ptr[2]-local_pre)*esigma_ptr[X]*egradvel_ptr[3]
                                      + (esigma_ptr[Y] + esigma_ptr[X])*este_ptr[3]*egradvel_ptr[4];
                    ealpha_ptr[Y] += -1.0*ste_damp*dt -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ealpha_ptr[Y]*dt;
                }

                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] < 1.0 - TOLVOF) continue;
                    double mden_damp = (-(esigma_ptr[X]+esigma_ptr[Y])*(mden[matid] - edampref_ptr[0]))*dt; // different from pml2
                    double meng_damp = (-(esigma_ptr[X]+esigma_ptr[Y])*(meng[matid] - edampref_ptr[1]) - ealpha_ptr[Y])*dt;
                    mden[matid] += mden_damp;
                    meng[matid] += meng_damp;
                    break;
                }
            }
        }
    }
}


void damp_on_cell_sbc(sale2d_var * _sale,double dt,double t)
{
    /*
     * this implement sponge boundary condition
     * the variables on cells need to be damped is
     * m_den, m_eng
     * the velocity is on the vertex is not damped ...
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    if(_sale->damp_time <= t)
    {
        return;
    }

    for(int j = 0; j < ny; j++) {
        for (int k = 0; k < nx; ++k) {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if(lnt & useless_section) continue;
            double * evof, *mden, *meng;
            _sale->foe->get_double(&lid,&evof,_sale->e_vof);
            _sale->foe->get_double(&lid,&mden,_sale->m_den);
            _sale->foe->get_double(&lid,&meng,_sale->m_eng);

            /*
             * skip the cell have vacuum
             * absorbing boundary should not be applied to ejector
             */
            if(evof[VACUUM] > TOLVOF) continue;

            double * eden_ptr, *egradvel_ptr,*evel_ptr, *epre_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);

            double *ealpha_ptr,*esigma_ptr;
            _sale->foe->get_double(&lid,&ealpha_ptr,_sale->e_alpha);
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);

            double *edampref_ptr;
            _sale->foe->get_double(&lid,&edampref_ptr,_sale->e_dampref);

            // damp on y direction
            if(esigma_ptr[Y] < 1.0)
            {
                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] < TOLVOF) continue;
                    // calculate delta variables
                    double delta_damp_den = mden[matid] - edampref_ptr[0];
                    double delta_damp_eng = meng[matid] - edampref_ptr[1];
                    mden[matid] += delta_damp_den*esigma_ptr[Y] + edampref_ptr[0];
                    meng[matid] += delta_damp_eng*esigma_ptr[Y] + edampref_ptr[1];
                }
            }

            //damp on X direction
            if(esigma_ptr[X] < 1.0)
            {
                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] < TOLVOF) continue;
                    // calculate delta variables
                    double delta_damp_den = mden[matid] - edampref_ptr[0];
                    double delta_damp_eng = meng[matid] - edampref_ptr[1];
                    mden[matid] = delta_damp_den*esigma_ptr[X] + edampref_ptr[0];
                    meng[matid] = delta_damp_eng*esigma_ptr[X]*esigma_ptr[X] + edampref_ptr[1];
                }
            }
        }
    }

}
