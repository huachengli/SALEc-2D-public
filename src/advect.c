//
// Created by huacheng on 3/2/23.
//
/*
 * functions on advect/calculate flux between cells
 * advect_flux
 * calculate_flux
 * and some functions which were piece/block of old version calculate/advect_flux
 */


#include "sale2d.h"

void advect_flux(damp_type_t dampType,sale2d_var * _sale, node_type _nt, double dt)
{
    switch(dampType)
    {
        case pml3:
            advect_flux_pml3(_sale,_nt,dt);
            break;
        default:
            advect_flux_default(_sale,_nt);
    }
}

void advect_flux_default(sale2d_var * _sale, node_type _nt)
{
    /*
     * update the energy and density through advect&compress
     * rewritten with get_double
     *
     * advect_flux is coupled with _calculate_flux_set_bphi
     * the former update element state according the values stored in dphi
     */

    field_var * _var = _sale->e_vof;
    field_opt * _opt = _sale->foe;

    proc_info_2d * proc_self = (proc_info_2d *) _var->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int nmat = _sale->nmat;

    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k, .y=j};
            node_type lnt = useless_section;
            _opt->get_nty(&lid,&lnt,_var);
            if(!(_nt&lnt)) continue;

            fv_id_2d bxl = {k,j}, bxr={k+1,j}, byb={k,j}, byu={k,j+1};

            double *bphi_dptr[4];
            _sale->fobx->get_double(&bxl,&bphi_dptr[0],_sale->bx_phi);
            _sale->fobx->get_double(&bxr,&bphi_dptr[2],_sale->bx_phi);
            _sale->foby->get_double(&byb,&bphi_dptr[1],_sale->by_phi);
            _sale->foby->get_double(&byu,&bphi_dptr[3],_sale->by_phi);

            double direction[4] = {1.0,1.0,-1.0,-1.0};

            double dphi[MAXPHI];
            /*
             * variables in dphi,
             * 3*nmat variables flux with volume + (phi_length - 3*nmat) variables flux with mass
             * vol(matid)   , 0 <= matid < nmat
             * mass(matid)  ,
             * energy(matid),
             * momentum(X),  momentum(Y),
             * stress(XX), stress(XY), stress(YY), stress(TH)
             * damage,
             * tps,
             * vibration,
             * ejecta
             * cold volume strain
             * porosity
             */

            int phi_len = _sale->phi_length;

            if(_sale->phi_length > MAXPHI)
            {
                fprintf(stderr,"Error: %s[%d]: too many materials, phi_length = %d\n",__func__ ,__LINE__,phi_len);
            }

            vec_zero(dphi,phi_len);

            // calculate sum of flux from bx,by
            double dphi_mass = 0.;

            for(int i=0;i<4;++i) // iteration on different direction
            {
                for(int matid=0;matid<nmat;++matid)
                {
                    double mat_vol_i = direction[i]*bphi_dptr[i][matid + 0*nmat];
                    dphi[matid + 0*nmat] += mat_vol_i; // vol
                    dphi[matid + 1*nmat] += mat_vol_i*bphi_dptr[i][matid + 1*nmat]; // mass
                    dphi[matid + 2*nmat] += mat_vol_i*bphi_dptr[i][matid + 2*nmat]; // energy
                    dphi_mass += mat_vol_i*bphi_dptr[i][matid + 1*nmat]; // total mass
                    scaler_add(dphi + 3*nmat, 10 ,bphi_dptr[i] + 3*nmat, mat_vol_i*bphi_dptr[i][matid + 1*nmat]);
                    if(_sale->porosity_model)
                    {
                        dphi[        3*nmat + 10] += (bphi_dptr[i][        3*nmat + 10]) * (mat_vol_i*bphi_dptr[i][matid + 1*nmat]);
                        dphi[matid + 3*nmat + 11] += (bphi_dptr[i][matid + 3*nmat + 11]) * (mat_vol_i*bphi_dptr[i][matid + 1*nmat]);
                    }
                }
            }

            // construct phi in element
            double * mden, *meng, *evof;
            _opt->get_double(&lid,&mden,_sale->m_den);
            _opt->get_double(&lid,&meng,_sale->m_eng);
            _opt->get_double(&lid,&evof,_sale->m_vof);

            double *edvel_ptr;
            _opt->get_double(&lid,&edvel_ptr,_sale->e_dvel);

            double *evol_ptr;
            _opt->get_double(&lid,&evol_ptr,_sale->e_vol);

            double *este;
            _opt->get_double(&lid,&este,_sale->e_ste);

            double *edam_ptr;
            _opt->get_double(&lid,&edam_ptr,_sale->e_dam);

            double *etps_ptr;
            _opt->get_double(&lid,&etps_ptr,_sale->e_tps);

            double * eden_ptr = NULL;
            _opt->get_double(&lid,&eden_ptr,_sale->e_den);

            double * evib_ptr;
            _opt->get_double(&lid,&evib_ptr,_sale->e_vib);

            double * eejt_ptr;
            _opt->get_double(&lid,&eejt_ptr,_sale->e_ejt);

            double * ecvs_ptr, *mpty_ptr;
            _opt->get_double(&lid,&ecvs_ptr,_sale->e_cvs);
            _opt->get_double(&lid,&mpty_ptr,_sale->m_pty);

            double evol = *evol_ptr;
            double emass = (*evol_ptr) * (*eden_ptr);
            double sphi[MAXPHI];
            vec_zero(sphi,phi_len);

            for(int matid=0;matid<nmat;++matid)
            {
                sphi[matid + 0*nmat] += evof[matid] * evol; // material volume
                sphi[matid + 1*nmat] += evof[matid] * evol * mden[matid]; // material mass
                sphi[matid + 2*nmat] += evof[matid] * evol * meng[matid]; // material energy
            }

            /*
             * [deleted] evel is changed to the velocity increment after advect
             * velocity increment has been moved into e_dvel
             * should not be add to sphi
             * scaler_add(sphi + 3*nmat,2,evel,emass);
             *
             * sphi =
             * 3*nmat : {material volume, material mass, material energy} * nmat
             * 2 : velocity/ momentum increment
             * 4 : stress tensor {rr,rz,zz,th}
             * 1 : damage
             * 1 : tps/total plastic strain
             * 1 : vibration
             * 1 : ejt
             * 1 : cvs
             * nmat: pty
             */

            scaler_add(sphi+3*nmat+2,4,este,emass);
            sphi[3*nmat+6] = emass * (*edam_ptr);
            sphi[3*nmat+7] = emass * (*etps_ptr);
            sphi[3*nmat+8] = emass * (*evib_ptr);
            sphi[3*nmat+9] = emass * (*eejt_ptr);

            if(_sale->porosity_model)
            {
                sphi[3*nmat + 10] = emass * (*ecvs_ptr);
                for(int matid=0;matid<nmat;++matid)
                {
                    sphi[3*nmat + 11 + matid] = evof[matid] * evol * mden[matid] * mpty_ptr[matid];
                }
            }
            // add dphi to sphi
            scaler_add(sphi,phi_len,dphi,1.0);

            double Zeta = 0.0; // sum of vof except for vaccum
            emass = 0.;
            for(int matid=0;matid<nmat;++matid)
            {
                if(sphi[matid]<0.)
                {
                    sphi[matid] = 0.;
                    sphi[matid+nmat] = 0.;
                }

                if(VACUUM == matid) continue;
                Zeta += sphi[matid];
                emass += sphi[matid+nmat];
            }

            if(Zeta > evol*(1.0 - TOLVOF) || (Zeta > 0. && sphi[0] < evol * TOLVOF))
            {
                // compress/expand materials in full elements
                scaler_move(evof,nmat,sphi,1.0/Zeta);
                evof[0] = 0.0;
            } else if(Zeta < evol * TOLVOF)
            {
                vec_zero(evof,nmat);
                evof[0] = 1.0;
            } else
            {
                // compress vacuum and other materials consistently
                if(_sale->partpres)
                {
                    scaler_move(evof,nmat,sphi,1.0/(Zeta+sphi[0]));
                } else
                {
                    sphi[0] = evol - Zeta;
                    scaler_move(evof,nmat,sphi,1.0/evol);
                }
            }


            for(int matid=1;matid<nmat;++matid)
            {
                if(evof[matid] < TOLVOF || sphi[matid + 1*nmat] < TOLRHO*TOLVOF*evol)
                {
                    meng[matid] = 0.0;
                    mden[matid] = 0.0;
                    evof[matid] = 0.0;
                    if(_sale->porosity_model) mpty_ptr[matid] = 1.0;
                } else
                {
                    meng[matid] = sphi[matid + 2*nmat] / (evof[matid] * evol);
                    mden[matid] = sphi[matid + 1*nmat] / (evof[matid] * evol);
                    if(_sale->porosity_model) mpty_ptr[matid] = sphi[3*nmat + 11 + matid]/sphi[matid + 1*nmat];
                }
            }

            // update the velocity/stree/damage in element
            if(emass > TOLRHO * evol)
            {
                scaler_move(edvel_ptr,2,sphi+3*nmat,1.0/emass);
                scaler_move(este,4,sphi+3*nmat+2,1.0/emass);
                *edam_ptr = sphi[3*nmat+6]/emass;
                *etps_ptr = sphi[3*nmat+7]/emass;
                *evib_ptr = sphi[3*nmat+8]/emass;
                *eejt_ptr = sphi[3*nmat+9]/emass;
                if(_sale->porosity_model) *ecvs_ptr = sphi[3*nmat + 10]/emass;
            } else
            {
                vec_zero(edvel_ptr,2);
                vec_zero(este,4);
                *edam_ptr = 0.;
                *etps_ptr = 0.;
                *evib_ptr = 0.;
                *eejt_ptr = 0.;
                if(_sale->porosity_model) *ecvs_ptr = 0.;
            }

            // update the density of element
            *eden_ptr = emass/(*evol_ptr);
            // cutoff tps
            *etps_ptr = Wind(*etps_ptr,0.0,10.0);
            *edam_ptr = Wind(*edam_ptr,0.0,1.0);
            *eejt_ptr = Wind(*eejt_ptr,0.0,1.0);

            if(_sale->porosity_model){
                for(int matid=1;matid<nmat;++matid)
                {
                    mpty_ptr[matid] = Wind(mpty_ptr[matid],1.0,_sale->etb[matid].ref->alphamax2);
                }
            }

            // remove element with low density
            if(*eden_ptr < _sale->density_cutoff || *eden_ptr < 10.0*_sale->density_cutoff* (1.0-evof[VACUUM]))
            {
                vec_zero(evof,nmat); evof[VACUUM] = 1.0;
                vec_zero(mden,nmat);
                vec_zero(meng,nmat);
                vec_zero(edvel_ptr,2);
                vec_zero(este,4);
                *edam_ptr = 0.;
                *etps_ptr = 0.;
                *evib_ptr = 0.;
                *eejt_ptr = 0.;
                if(_sale->porosity_model)
                {
                    *ecvs_ptr = 0.;
                    for(int matid=1;matid<nmat;++matid)
                    {
                        mpty_ptr[matid] = 1.0;
                    }
                }
            }
        }
    }
}

void advect_flux_pml3(sale2d_var * _sale, node_type _nt, double dt)
{
    /*
     * update the energy and density through advect&compress
     * rewritten with get_double
     *
     * advect_flux is coupled with _calculate_flux_set_bphi
     * the former update element state according the values stored in dphi
     *
     * add PML auxiliary variables e_alpha[0:1]
     */

    field_var * _var = _sale->e_vof;
    field_opt * _opt = _sale->foe;

    proc_info_2d * proc_self = (proc_info_2d *) _var->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int nmat = _sale->nmat;

    for(int j=0;j<ny;j++)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k, .y=j};
            node_type lnt = useless_section;
            _opt->get_nty(&lid,&lnt,_var);
            if(!(_nt&lnt)) continue;

            fv_id_2d bxl = {k,j}, bxr={k+1,j}, byb={k,j}, byu={k,j+1};

            double *bphi_dptr[4];
            _sale->fobx->get_double(&bxl,&bphi_dptr[0],_sale->bx_phi);
            _sale->fobx->get_double(&bxr,&bphi_dptr[2],_sale->bx_phi);
            _sale->foby->get_double(&byb,&bphi_dptr[1],_sale->by_phi);
            _sale->foby->get_double(&byu,&bphi_dptr[3],_sale->by_phi);

            double direction[4] = {1.0,1.0,-1.0,-1.0};

            double *ealpha_ptr,*ebeta_ptr,*esigma_ptr;
            _sale->foe->get_double(&lid,&ealpha_ptr,_sale->e_alpha);
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);
            _sale->foe->get_double(&lid,&ebeta_ptr,_sale->e_beta);
            double pml_para[4] = {-esigma_ptr[Y],-esigma_ptr[X],-esigma_ptr[Y],-esigma_ptr[X]};

            double dphi[MAXPHI];
            /*
             * variables in dphi,
             * 3*nmat variables flux with volume + (phi_length - 3*nmat) variables flux with mass
             * vol(matid)   , 0 <= matid < nmat
             * mass(matid)  ,
             * energy(matid),
             * momentum(X),  momentum(Y),
             * stress(XX), stress(XY), stress(YY), stress(TH)
             * damage,
             * tps,
             * vibration,
             * ejecta
             */


            int phi_len = _sale->phi_length;

            if(_sale->phi_length > MAXPHI)
            {
                fprintf(stderr,"Error: %s[%d]: too many materials, phi_length = %d\n",__func__ ,__LINE__,phi_len);
            }

            vec_zero(dphi,phi_len);

            // calculate sum of flux from bx,by
            double dphi_mass = 0.;

            for(int i=0;i<4;++i) // iteration on different direction
            {
                for(int matid=0;matid<nmat;++matid)
                {
                    double mat_vol_i = direction[i]*bphi_dptr[i][matid + 0*nmat];
                    dphi[matid + 0*nmat] += mat_vol_i; // vol
                    dphi[matid + 1*nmat] += mat_vol_i*bphi_dptr[i][matid + 1*nmat]; // mass
                    dphi[matid + 2*nmat] += mat_vol_i*bphi_dptr[i][matid + 2*nmat]; // energy
                    dphi_mass += mat_vol_i*bphi_dptr[i][matid + 1*nmat]; // total mass
                    scaler_add(dphi + 3*nmat, phi_len - 3*nmat,bphi_dptr[i] + 3*nmat, mat_vol_i*bphi_dptr[i][matid + 1*nmat]);

                    // insert auxiliary of flux equation
                    if(esigma_ptr[X] > 0. || esigma_ptr[Y] > 0.)
                    {
                        //
                        // flux mass =  mat_vol_i*bphi_dptr[i][matid + 1*nmat]
                        // pml attenuation parameter: pml_para[i]
                        // velocity  = [X] = dphi + 3*nmat+0,[Y]= dphi + 3*nmat+1
                        // scaler_add(ealpha_ptr,2,dphi + 3*nmat,mat_vol_i*bphi_dptr[i][matid + 1*nmat]*pml_para[i]);
                        ebeta_ptr[0] += mat_vol_i*bphi_dptr[i][matid + 1*nmat] * pml_para[i];
                        ebeta_ptr[1] += mat_vol_i*bphi_dptr[i][matid + 2*nmat] * pml_para[i];
                        ebeta_ptr[2] += mat_vol_i*bphi_dptr[i][matid + 1*nmat] * pml_para[i] * bphi_dptr[i][3*nmat + X];
                        ebeta_ptr[3] += mat_vol_i*bphi_dptr[i][matid + 1*nmat] * pml_para[i] * bphi_dptr[i][3*nmat + Y];
                    }
                }
            }

            // construct phi in element
            double * mden, *meng, *evof;
            _opt->get_double(&lid,&mden,_sale->m_den);
            _opt->get_double(&lid,&meng,_sale->m_eng);
            _opt->get_double(&lid,&evof,_sale->m_vof);

            double *edvel_ptr;
            _opt->get_double(&lid,&edvel_ptr,_sale->e_dvel);

            double *evol_ptr;
            _opt->get_double(&lid,&evol_ptr,_sale->e_vol);

            double *este;
            _opt->get_double(&lid,&este,_sale->e_ste);

            double *edam_ptr;
            _opt->get_double(&lid,&edam_ptr,_sale->e_dam);

            double *etps_ptr;
            _opt->get_double(&lid,&etps_ptr,_sale->e_tps);

            double * eden_ptr = NULL;
            _opt->get_double(&lid,&eden_ptr,_sale->e_den);

            double * evib_ptr;
            _opt->get_double(&lid,&evib_ptr,_sale->e_vib);

            double * eejt_ptr;
            _opt->get_double(&lid,&eejt_ptr,_sale->e_ejt);

            double evol = *evol_ptr;
            double emass = (*evol_ptr) * (*eden_ptr);
            double sphi[MAXPHI];
            vec_zero(sphi,phi_len);

            for(int matid=0;matid<nmat;++matid)
            {
                sphi[matid + 0*nmat] += evof[matid] * evol; // material volume
                sphi[matid + 1*nmat] += evof[matid] * evol * mden[matid]; // material mass
                sphi[matid + 2*nmat] += evof[matid] * evol * meng[matid]; // material energy
            }

            /*
             * [deleted] evel is changed to the velocity increment after advect
             * velocity increment has been moved into e_dvel
             * should not be add to sphi
             * scaler_add(sphi + 3*nmat,2,evel,emass);
             *
             * sphi =
             * 3*nmat : {material volume, material mass, material energy} * nmat
             * 2 : velocity/momentum increment : index = 3*nmat+0 and 3*nmat + 1
             * 4 : stress tensor {rr,rz,zz,th}
             * 1 : damage
             * 1 : tps/total plastic strain
             * 1 : vibration
             * 1 : ejt
             */


            scaler_add(sphi+3*nmat+2,4,este,emass);
            sphi[3*nmat+6] = emass * (*edam_ptr);
            sphi[3*nmat+7] = emass * (*etps_ptr);
            sphi[3*nmat+8] = emass * (*evib_ptr);
            sphi[3*nmat+9] = emass * (*eejt_ptr);

            // add dphi to sphi
            scaler_add(sphi,phi_len,dphi,1.0);
            // add pml damping term
            if(esigma_ptr[X]  + esigma_ptr[Y] > 0.)
            {
                // damp on auxiliary vars
                ebeta_ptr[0]  += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ebeta_ptr[0]*dt;
                ebeta_ptr[1]  += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ebeta_ptr[1]*dt;
                ebeta_ptr[2]  += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ebeta_ptr[2]*dt;
                ebeta_ptr[3]  += -1.0*(_sale->alpha_damp_x + _sale->alpha_damp_y)*ebeta_ptr[3]*dt;

                if(sphi[VACUUM] < TOLVOF)
                {
                    // first term: -(dx+dy)*Rho*Vel*dt, shared with damp_on_cell;
                    // second term: related with auxiliary terms
                    sphi[3*nmat+X] += ebeta_ptr[2]*dt;
                    sphi[3*nmat+Y] += ebeta_ptr[3]*dt;
                }

                // first term has been included in damp on cell
                // damp on mass/energy
                for(int matid=1;matid<nmat;++matid)
                {
                    if(sphi[matid]>TOLVOF)
                    {
                        sphi[matid+1*nmat] += ebeta_ptr[0]*dt;
                        sphi[matid+2*nmat] += ebeta_ptr[1]*dt;
                        break;
                    }
                }
            }

            double Zeta = 0.0; // sum of vof except for vaccum
            emass = 0.;
            for(int matid=0;matid<nmat;++matid)
            {
                if(sphi[matid]<0. || sphi[matid+nmat]<0.)
                {
                    sphi[matid] = 0.;
                    sphi[matid+nmat] = 0.;
                }

                if(VACUUM == matid) continue;
                Zeta += sphi[matid];
                emass += sphi[matid+nmat];
            }

            if(Zeta > evol*(1.0 - TOLVOF) || (Zeta > 0. && sphi[0] < evol * TOLVOF))
            {
                // compress/expand materials in full elements
                scaler_move(evof,nmat,sphi,1.0/Zeta);
                evof[0] = 0.0;
            } else if(Zeta < evol * TOLVOF)
            {
                vec_zero(evof,nmat);
                evof[0] = 1.0;
            } else
            {
                // compress vacuum and other materials consistently
                if(_sale->partpres)
                {
                    scaler_move(evof,nmat,sphi,1.0/(Zeta+sphi[0]));
                } else
                {
                    sphi[0] = evol - Zeta;
                    scaler_move(evof,nmat,sphi,1.0/evol);
                }
            }


            for(int matid=1;matid<nmat;++matid)
            {
                if(evof[matid] < TOLVOF)
                {
                    meng[matid] = 0.0;
                    mden[matid] = 0.0;
                    evof[matid] = 0.0;
                } else
                {
                    meng[matid] = sphi[matid + 2*nmat] / (evof[matid] * evol);
                    mden[matid] = sphi[matid + 1*nmat] / (evof[matid] * evol);
                }
            }

            // update the velocity/stree/damage in element
            if(emass > TOLRHO * evol)
            {
                scaler_move(edvel_ptr,2,sphi+3*nmat,1.0/emass);
                scaler_move(este,4,sphi+3*nmat+2,1.0/emass);
                *edam_ptr = sphi[3*nmat+6]/emass;
                *etps_ptr = sphi[3*nmat+7]/emass;
                *evib_ptr = sphi[3*nmat+8]/emass;
                *eejt_ptr = sphi[3*nmat+9]/emass;
            } else
            {
                vec_zero(edvel_ptr,2);
                vec_zero(este,4);
                *edam_ptr = 0.;
                *etps_ptr = 0.;
                *evib_ptr = 0.;
                *eejt_ptr = 0.;
            }

            // update the density of element
            *eden_ptr = emass/(*evol_ptr);
            // cutoff tps
            *etps_ptr = Wind(*etps_ptr,0.0,10.0);
            *edam_ptr = Wind(*edam_ptr,0.0,1.0);
            *eejt_ptr = Wind(*eejt_ptr,0.0,1.0);

            // remove element with low density
            if(*eden_ptr < _sale->density_cutoff)
            {
                vec_zero(evof,nmat); evof[VACUUM] = 1.0;
                vec_zero(mden,nmat);
                vec_zero(meng,nmat);
                vec_zero(edvel_ptr,2);
                vec_zero(este,4);
                *edam_ptr = 0.;
                *etps_ptr = 0.;
                *evib_ptr = 0.;
                *eejt_ptr = 0.;
            }
        }
    }
}



void clean_flux(sale2d_var * _sale)
{
    proc_info_2d * proc_bx = _sale->bx_phi->proc_self;
    int bx_nx = proc_bx->nx, bx_ny = proc_bx->ny;

    for(int j=0;j<bx_ny;++j)
    for(int k=0;k<bx_nx;k++)
    {
        {
            double * bxphi;
            fv_id_2d lid = {.x=k,.y=j};
            _sale->fobx->get_double(&lid,&bxphi,_sale->bx_phi);
            vec_zero(bxphi,_sale->bx_phi->bytes_per_unit/ sizeof(double));
        }
    }

    proc_info_2d * proc_by = _sale->by_phi->proc_self;
    int by_nx = proc_by->nx, by_ny = proc_by->ny;

    for(int j=0;j<by_ny;++j)
    for(int k=0;k<by_nx;++k)
    {
        {
            double *byphi;
            fv_id_2d lid = {.x=k,.y=j};
            _sale->foby->get_double(&lid,&byphi,_sale->by_phi);
            vec_zero(byphi,_sale->by_phi->bytes_per_unit/ sizeof(double));
        }
    }
}

void calculate_flux_primary(sale2d_var * _sale, node_type _nt)
{
    /*
     * select a material as the primary
     * except for this material, all other material is mixed uniformly
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    for (int j = 0; j < ny; ++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d elid = {k, j};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&elid, &lnt, _sale->e_vof);
            if (!(lnt & _nt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};

            double *vpos_dptr[4], *vdis_dptr[4], *vvof_dptr[4];
            for (int i = 0; i < 4; ++i) {
                _sale->fov->get_double(vlid[i], &vpos_dptr[i], _sale->v_pos);
                _sale->fov->get_double(vlid[i], &vdis_dptr[i], _sale->v_dis);
                _sale->fov->get_double(vlid[i], &vvof_dptr[i], _sale->v_vof);
            }

            double l_vol = get_flux_volcyl(vpos_dptr[0],vpos_dptr[3],vdis_dptr[0],vdis_dptr[3]);
            double r_vol = get_flux_volcyl(vpos_dptr[2],vpos_dptr[1],vdis_dptr[2],vdis_dptr[1]);
            double b_vol = get_flux_volcyl(vpos_dptr[1],vpos_dptr[0],vdis_dptr[1],vdis_dptr[0]);
            double u_vol = get_flux_volcyl(vpos_dptr[3],vpos_dptr[2],vdis_dptr[3],vdis_dptr[2]);

            double *egrad_ptr = NULL;
            _sale->foe->get_double(&elid,&egrad_ptr,_sale->e_grad);

            char *epos_cptr;
            double *escale_dptr;
            _sale->foe->get_data2(&elid, &epos_cptr, _sale->e_pos);
            escale_dptr = (double *) (epos_cptr) + 2;

            double *evof_ptr;
            _sale->foe->get_double(&elid, &evof_ptr, _sale->e_vof);

            double *ewpt_ptr;
            _sale->foe->get_double(&elid,&ewpt_ptr,_sale->e_wpt);

            double *evol_ptr;
            _sale->foe->get_double(&elid,&evol_ptr,_sale->e_vol);

            double flux_vol[MAXMAT][4] = {0};

            int primary_mat = 0; // the default value is vacuum

            if(evof_ptr[0] < TOLVOF)
            {
                /*
                 * this element is full of materials except vacuum
                 * change the primary materials to the max vof material
                 */
                for(int cutid = 1;cutid<nmat;++cutid)
                {
                    if(evof_ptr[cutid] > evof_ptr[primary_mat])
                    {
                        primary_mat = cutid;
                    }
                }
            }

            // if the cell is full of primary materials
            // no special treated needed
            if(evof_ptr[primary_mat] > 1.0 - TOLVOF)
            {
                for(int cutid=0;cutid<nmat;++cutid)
                {
                    if(cutid==primary_mat)
                    {
                        flux_vol[cutid][0] =      MIN(l_vol,0);
                        flux_vol[cutid][2] = -1.0*MIN(r_vol,0);
                        flux_vol[cutid][1] =      MIN(b_vol,0);
                        flux_vol[cutid][3] = -1.0*MIN(u_vol,0);
                    }
                    else
                    {
                        flux_vol[cutid][0] = 0.;
                        flux_vol[cutid][2] = 0.;
                        flux_vol[cutid][1] = 0.;
                        flux_vol[cutid][3] = 0.;
                    }
                }

                if(primary_mat == 0)
                {
                    ewpt_ptr[0] = ewpt_ptr[1] = ewpt_ptr[2] = ewpt_ptr[3] = 1.0;
                } else
                {
                    ewpt_ptr[0] = ewpt_ptr[1] = ewpt_ptr[2] = ewpt_ptr[3] = -0.25;
                }
            }
            else
            {
                // branch for mixed cells
                double local_vof_grad[2];
                double local_vvof[4] = {vvof_dptr[0][primary_mat], vvof_dptr[1][primary_mat],
                                        vvof_dptr[2][primary_mat], vvof_dptr[3][primary_mat]};
                local_vof_grad[0] = vec_dot(local_vvof,egrad_ptr + 0, 4);
                local_vof_grad[1] = vec_dot(local_vvof,egrad_ptr + 4, 4);
                vec_norm(local_vof_grad,2);

                double cut_a = escale_dptr[0] * fabs(local_vof_grad[0]);
                double cut_b = escale_dptr[1] * fabs(local_vof_grad[1]);

                if (cut_a < cut_b) {
                    double tmp = cut_a;
                    cut_a = cut_b;
                    cut_b = tmp;
                }

                double dstar = cut_length(evof_ptr[primary_mat], (double[]) {cut_a, cut_b});
                double local_wpt[4], sum_wpt = 0., max_wpt;
                double dis_wpt[4];
                for (int i = 0; i < 4; ++i) {
                    local_wpt[i] = vec_dot(local_vof_grad, vpos_dptr[i], 2);
                    sum_wpt += local_wpt[i];
                    dis_wpt[i] = vec_dot(local_vof_grad,vdis_dptr[i],2);
                }

                max_wpt = MAX(MAX(local_wpt[0], local_wpt[1]), MAX(local_wpt[2], local_wpt[3]));
                for (int i = 0; i < 4; ++i) {
                    local_wpt[i] = local_wpt[i] - max_wpt + dstar;
                }

                if(primary_mat == 0)
                {
                    vec_copy(local_wpt,ewpt_ptr,4);
                } else
                {
                    ewpt_ptr[0] = ewpt_ptr[1] = ewpt_ptr[2] = ewpt_ptr[3] = -0.25;
                }


                double flux_primary[4] = {0.};
                flux_primary[0] = geometric_flux_ratio(local_wpt[3], local_wpt[0]);
                flux_primary[2] = geometric_flux_ratio(local_wpt[2], local_wpt[1]);
                flux_primary[1] = geometric_flux_ratio(local_wpt[0], local_wpt[1]);
                flux_primary[3] = geometric_flux_ratio(local_wpt[3], local_wpt[2]);

                // before this point we just calculate the ratio of primary materials should be fluxed
                for(int cutid=0;cutid<nmat;++cutid)
                {
                    if(cutid == primary_mat)
                    {
                        flux_vol[cutid][0] =      MIN(l_vol,0.0)*flux_primary[0];
                        flux_vol[cutid][2] = -1.0*MIN(r_vol,0.0)*flux_primary[2];
                        flux_vol[cutid][1] =      MIN(b_vol,0.0)*flux_primary[1];
                        flux_vol[cutid][3] = -1.0*MIN(u_vol,0.0)*flux_primary[3];
                    } else
                    {
                        flux_vol[cutid][0] =      MIN(l_vol,0.0)*(1.0-flux_primary[0])*evof_ptr[cutid]/(1.0-evof_ptr[primary_mat]);
                        flux_vol[cutid][2] = -1.0*MIN(r_vol,0.0)*(1.0-flux_primary[2])*evof_ptr[cutid]/(1.0-evof_ptr[primary_mat]);
                        flux_vol[cutid][1] =      MIN(b_vol,0.0)*(1.0-flux_primary[1])*evof_ptr[cutid]/(1.0-evof_ptr[primary_mat]);
                        flux_vol[cutid][3] = -1.0*MIN(u_vol,0.0)*(1.0-flux_primary[3])*evof_ptr[cutid]/(1.0-evof_ptr[primary_mat]);
                    }
                }
            }

            // cutoff flux
            for(int cutid=1;cutid<nmat;++cutid)
            {
                double cutid_vol = -1.0*vec_dot(flux_vol[cutid],(double []){1.0,1.0,-1.0,-1.0},4);
                if(cutid_vol > evof_ptr[cutid]*(*evol_ptr))
                {
                    vec_scaler(flux_vol[cutid],4,evof_ptr[cutid]*(*evol_ptr)/cutid_vol);
                }
            }

            double lbru_vol[4] = {l_vol,b_vol,r_vol,u_vol};
            _calculate_flux_set_bphi(_sale,&elid,lbru_vol,flux_vol);

        }
    }
}

mpi_vec_state sync_all(sale2d_var * _sale)
{
    // this function is just used for test
    // BLOCKED SYNC (using carefully)

    sync_start(_sale->m_vof);
    sync_complete(_sale->m_vof);


// region, 1
    sync_start(_sale->e_vel);
    sync_start(_sale->e_den);
    sync_start(_sale->e_wpt);
    sync_complete(_sale->e_vel);
    sync_complete(_sale->e_den);
    sync_complete(_sale->e_wpt);

    cycle_control rand_ccl;
    rand_ccl.dt = 1.0*rand()/RAND_MAX;
    MPI_Allreduce(&(rand_ccl.local_dt),&(rand_ccl.dt),1,MPI_DOUBLE,MPI_MIN,_sale->e_pre->comm);
    sync_start(_sale->e_pre);
    sync_start(_sale->e_ste);
    sync_start(_sale->e_q);

    sync_complete(_sale->e_pre);
    sync_complete(_sale->e_ste);
    sync_complete(_sale->e_q);



    sync_start(_sale->m_den);
    sync_start(_sale->m_eng);
    sync_start(_sale->e_vel);
    sync_start(_sale->e_tps);
    sync_start(_sale->e_dam);
    sync_start(_sale->e_vib);


    sync_complete(_sale->m_den);
    sync_complete(_sale->m_eng);
    sync_complete(_sale->e_vel);
    sync_complete(_sale->e_tps);
    sync_complete(_sale->e_dam);
    sync_complete(_sale->e_vib);

}

void calculate_flux_mixed(sale2d_var * _sale, node_type _nt)
{
    /*
     * mixed implemented of algebraic and geometric internal boundary reconstruction
     * to capture the boundary between materials
     * select a material as the primary, geometrical will be used for primary
     * non-primary materials will use algebraic method
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    for (int j = 0; j < ny; ++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d elid = {k, j};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&elid, &lnt, _sale->e_vof);
            if (!(lnt & _nt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};
            fv_id_2d l_elid ={k-1,j}, r_elid={k+1,j}, b_elid={k,j-1}, u_elid={k,j+1};
            fv_id_2d *nei_elid[4] = {&l_elid,&b_elid,&r_elid,&u_elid};

            double *vpos_dptr[4], *vdis_dptr[4], *vvof_dptr[4];
            for (int i = 0; i < 4; ++i) {
                _sale->fov->get_double(vlid[i], &vpos_dptr[i], _sale->v_pos);
                _sale->fov->get_double(vlid[i], &vdis_dptr[i], _sale->v_dis);
                _sale->fov->get_double(vlid[i], &vvof_dptr[i], _sale->v_vof);
            }

            double l_vol = get_flux_volcyl(vpos_dptr[0],vpos_dptr[3],vdis_dptr[0],vdis_dptr[3]);
            double r_vol = get_flux_volcyl(vpos_dptr[2],vpos_dptr[1],vdis_dptr[2],vdis_dptr[1]);
            double b_vol = get_flux_volcyl(vpos_dptr[1],vpos_dptr[0],vdis_dptr[1],vdis_dptr[0]);
            double u_vol = get_flux_volcyl(vpos_dptr[3],vpos_dptr[2],vdis_dptr[3],vdis_dptr[2]);
            double lbru_vol[4] = {l_vol,b_vol,r_vol,u_vol};

            double *egrad_ptr = NULL;
            _sale->foe->get_double(&elid,&egrad_ptr,_sale->e_grad);

            double *epos_dptr;
            double *escale_dptr;
            _sale->foe->get_double(&elid, &epos_dptr, _sale->e_pos);
            escale_dptr = (double *) (epos_dptr) + 2;

            double *evof_ptr;
            _sale->foe->get_double(&elid, &evof_ptr, _sale->e_vof);

            double *ewpt_ptr;
            _sale->foe->get_double(&elid,&ewpt_ptr,_sale->e_wpt);

            double *evol_ptr;
            _sale->foe->get_double(&elid,&evol_ptr,_sale->e_vol);

            double flux_vol[MAXMAT][4] = {0};

            int primary_mat = 0; // the default value is vacuum

            if(evof_ptr[0] < TOLVOF)
            {
                /*
                 * this element is full of materials except vacuum
                 * change the primary materials to the max vof material
                 */
                for(int cutid = 1;cutid<nmat;++cutid)
                {
                    if(evof_ptr[cutid] > evof_ptr[primary_mat])
                    {
                        primary_mat = cutid;
                    }
                }
            }

            // if the cell is full of primary materials
            // no special treated needed
            // the in-flux will be capped, only process the out-flux
            if(evof_ptr[primary_mat] > 1.0 - TOLVOF)
            {
                for(int cutid=0;cutid<nmat;++cutid)
                {
                    if(cutid==primary_mat)
                    {
                        flux_vol[cutid][0] =      MIN(l_vol,0);
                        flux_vol[cutid][2] = -1.0*MIN(r_vol,0);
                        flux_vol[cutid][1] =      MIN(b_vol,0);
                        flux_vol[cutid][3] = -1.0*MIN(u_vol,0);
                    }
                    else
                    {
                        flux_vol[cutid][0] = 0.;
                        flux_vol[cutid][2] = 0.;
                        flux_vol[cutid][1] = 0.;
                        flux_vol[cutid][3] = 0.;
                    }
                }

                if(primary_mat == 0)
                {
                    ewpt_ptr[0] = ewpt_ptr[1] = ewpt_ptr[2] = ewpt_ptr[3] = 1.0;
                } else
                {
                    ewpt_ptr[0] = ewpt_ptr[1] = ewpt_ptr[2] = ewpt_ptr[3] = -0.25;
                }
            }
            else
            {
                // branch for mixed cells
                double local_vof_grad[2];
                double local_vvof[4] = {vvof_dptr[0][primary_mat], vvof_dptr[1][primary_mat],
                                        vvof_dptr[2][primary_mat], vvof_dptr[3][primary_mat]};
                local_vof_grad[0] = vec_dot(local_vvof,egrad_ptr + 0, 4);
                local_vof_grad[1] = vec_dot(local_vvof,egrad_ptr + 4, 4);
                vec_norm(local_vof_grad,2);

                double cut_a = escale_dptr[0] * fabs(local_vof_grad[0]);
                double cut_b = escale_dptr[1] * fabs(local_vof_grad[1]);

                if (cut_a < cut_b) {
                    double tmp = cut_a;
                    cut_a = cut_b;
                    cut_b = tmp;
                }

                double dstar = cut_length(evof_ptr[primary_mat], (double[]) {cut_a, cut_b});
                double local_wpt[4], sum_wpt = 0., max_wpt;
                double dis_wpt[4];
                for (int i = 0; i < 4; ++i) {
                    local_wpt[i] = vec_dot(local_vof_grad, vpos_dptr[i], 2);
                    sum_wpt += local_wpt[i];
                    dis_wpt[i] = vec_dot(local_vof_grad,vdis_dptr[i],2);
                }

                max_wpt = MAX(MAX(local_wpt[0], local_wpt[1]), MAX(local_wpt[2], local_wpt[3]));
                for (int i = 0; i < 4; ++i) {
                    local_wpt[i] = local_wpt[i] - max_wpt + dstar;
                }

                if(primary_mat == 0)
                {
                    vec_copy(local_wpt,ewpt_ptr,4);
                } else
                {
                    ewpt_ptr[0] = ewpt_ptr[1] = ewpt_ptr[2] = ewpt_ptr[3] = -0.25;
                }


                double flux_primary[4] = {0.};
                flux_primary[0] = geometric_flux_ratio(local_wpt[3], local_wpt[0]);
                flux_primary[2] = geometric_flux_ratio(local_wpt[2], local_wpt[1]);
                flux_primary[1] = geometric_flux_ratio(local_wpt[0], local_wpt[1]);
                flux_primary[3] = geometric_flux_ratio(local_wpt[3], local_wpt[2]);
                // before this point we just calculate the ratio of primary materials should be fluxed

                // calculate courant number for alg method
                double cell_courant = MIN(l_vol,0.0) + MIN(r_vol,0) + MIN(b_vol,0) + MIN(u_vol,0);
                cell_courant = fabs(cell_courant)/evol_ptr[0];
                // fetch the vof value of number cells
                double * nei_evof_ptr[4];
                for(int i=0;i<4;++i)
                {
                    _sale->foe->get_double(nei_elid[i],&(nei_evof_ptr[i]),_sale->e_vof);
                    if(NULL == nei_evof_ptr[i]) nei_evof_ptr[i] = evof_ptr;
                }
                // load normal vector of boundary
                double *bound_normal[4];
                fv_id_2d bxl = {k,j}, bxr={k+1,j}, byb={k,j}, byu={k,j+1};
                _sale->fobx->get_double(&bxl,&bound_normal[0],_sale->bx_normal);
                _sale->fobx->get_double(&bxr,&bound_normal[2],_sale->bx_normal);
                _sale->foby->get_double(&byb,&bound_normal[1],_sale->by_normal);
                _sale->foby->get_double(&byu,&bound_normal[3],_sale->by_normal);

                double flux_nonprimary[MAXMAT][4] = {0.};
                double sum_nonprimary[4] = {0.};
                for(int matid=1;matid<nmat;++matid)
                {
                    if(primary_mat == matid) continue;
                    double nonpri_vof_grad[2] = {0.};
                    double nonpri_vvof[4] = {vvof_dptr[0][matid], vvof_dptr[1][matid],
                                            vvof_dptr[2][matid], vvof_dptr[3][matid]};
                    nonpri_vof_grad[0] = vec_dot(nonpri_vvof,egrad_ptr + 0, 4);
                    nonpri_vof_grad[1] = vec_dot(nonpri_vvof,egrad_ptr + 4, 4);
                    vec_norm(nonpri_vof_grad,2);
                    double nonpri_vof_upwind[4]   = { nei_evof_ptr[0][matid], nei_evof_ptr[1][matid],
                                                      nei_evof_ptr[2][matid], nei_evof_ptr[3][matid]};
                    double nonpri_vof_acceptor[4] = { nei_evof_ptr[2][matid], nei_evof_ptr[3][matid],
                                                      nei_evof_ptr[0][matid], nei_evof_ptr[1][matid]};
                    // the donor is fixed on {k,j}

                    /*
                     * the bound_normal may be reversed in some case
                     * but the CosTf*CosTf term makes no difference at reverse direction
                     */
                    for(int i=0;i<4;i++)
                    {
                        flux_nonprimary[matid][i] = algebraic_flux_ratio(nonpri_vof_upwind[i],evof_ptr[matid],nonpri_vof_acceptor[i],
                                                                         cell_courant, bound_normal[i],nonpri_vof_grad, AVOFALG);
                        sum_nonprimary[i] += flux_nonprimary[matid][i];
                    }
                }

                // check the sum_nonprimary[i] if approx zero,
                // use evof_ptr in flux other than flux_nonprimary
                // in this branch (1.0 - evof_ptr[primary_mat]) >= TOLVOF

                for(int i=0;i<4;++i)
                {
                    if(sum_nonprimary[i] > TOLVOF)
                    {
                        for(int matid=1;matid<nmat;++matid)
                        {
                            flux_nonprimary[matid][i] = flux_nonprimary[matid][i]/sum_nonprimary[i];
                        }
                    } else
                    {
                        for(int matid=1;matid<nmat;++matid)
                        {
                            flux_nonprimary[matid][i] = evof_ptr[matid]/(1.0 - evof_ptr[primary_mat]);
                        }
                    }
                }


                // insert flux_nonprimary into flux_ratio

                for(int cutid=0;cutid<nmat;++cutid)
                {
                    if(cutid == primary_mat)
                    {
                        flux_vol[cutid][0] =      MIN(l_vol,0.0)*flux_primary[0];
                        flux_vol[cutid][2] = -1.0*MIN(r_vol,0.0)*flux_primary[2];
                        flux_vol[cutid][1] =      MIN(b_vol,0.0)*flux_primary[1];
                        flux_vol[cutid][3] = -1.0*MIN(u_vol,0.0)*flux_primary[3];
                    } else
                    {
                        flux_vol[cutid][0] =      MIN(l_vol,0.0)*(1.0-flux_primary[0])*flux_nonprimary[cutid][0];
                        flux_vol[cutid][2] = -1.0*MIN(r_vol,0.0)*(1.0-flux_primary[2])*flux_nonprimary[cutid][2];
                        flux_vol[cutid][1] =      MIN(b_vol,0.0)*(1.0-flux_primary[1])*flux_nonprimary[cutid][1];
                        flux_vol[cutid][3] = -1.0*MIN(u_vol,0.0)*(1.0-flux_primary[3])*flux_nonprimary[cutid][3];
                    }
                }
            }

            // cutoff flux, avoid negative vof/density values
            for(int cutid=1;cutid<nmat;++cutid)
            {
                double cutid_vol = -1.0*vec_dot(flux_vol[cutid],(double []){1.0,1.0,-1.0,-1.0},4);
                if(cutid_vol > evof_ptr[cutid]*(*evol_ptr))
                {
                    vec_scaler(flux_vol[cutid],4,evof_ptr[cutid]*(*evol_ptr)/cutid_vol);
                }
            }
#ifdef CHECK_NAN
            if(vec_isnan(flux_vol[0],4*nmat) || vec_isinf(flux_vol[0],4*nmat))
            {
                fprintf(stderr,"Warning %s, Invalid flux vol\n",__func__ );
                for(int d=0;d<4;++d)
                {
                    for(int m=0;m<nmat;++m)
                    {
                        fprintf(stderr,"flux_vol[%d,%d]=%f,vvof[]=%f,",m,d,flux_vol[m][d],vvof_dptr[d][m]);
                    }
                    fprintf(stderr,"wpt(%d)=%f,lbru=%10.5e,pos(%f,%f) dis(%e,%e);",d,ewpt_ptr[d],lbru_vol[d],vpos_dptr[d][0],vpos_dptr[d][1],
                            vdis_dptr[d][0],vdis_dptr[d][1]);
                    fprintf(stderr,"\n");
                }
            }
#endif
            _calculate_flux_set_bphi(_sale,&elid,lbru_vol,flux_vol);

        }
    }
}

void _calculate_flux_set_bphi(sale2d_var * _sale, fv_id_2d * _elid, double *lbru_vol,double flux_vol[][4])
{
    /*
     * just the ending part of calculate_flux
     * used for fill bphi_ptr according to information
     * bphi[] record the volume of differnt materials across this boundary and the state of its upwind element.
     * the state of upwind element may be changed ! here we need copy the values not pointer
     */
    int k = _elid->x;
    int j = _elid->y;
    int nmat = _sale->nmat;

    fv_id_2d bxl = {k,j}, bxr={k+1,j}, byb={k,j}, byu={k,j+1};
    double * bphi_ptr[4];
    _sale->fobx->get_double(&bxl,&bphi_ptr[0],_sale->bx_phi);
    _sale->fobx->get_double(&bxr,&bphi_ptr[2],_sale->bx_phi);
    _sale->foby->get_double(&byb,&bphi_ptr[1],_sale->by_phi);
    _sale->foby->get_double(&byu,&bphi_ptr[3],_sale->by_phi);

    double * mden = NULL, *meng = NULL, *evel = NULL;
    _sale->foe->get_double(_elid,&mden,_sale->m_den);
    _sale->foe->get_double(_elid,&meng,_sale->m_eng);
    _sale->foe->get_double(_elid,&evel,_sale->e_vel);
    double * mpty=NULL, *ecvs = NULL;
    _sale->foe->get_double(_elid,&mpty,_sale->m_pty);
    _sale->foe->get_double(_elid,&ecvs,_sale->e_cvs);

    double *este=NULL, *edam=NULL, *etps;
    _sale->foe->get_double(_elid,&este,_sale->e_ste);
    _sale->foe->get_double(_elid,&edam,_sale->e_dam);
    _sale->foe->get_double(_elid,&etps,_sale->e_tps);

    double *evib_ptr;
    _sale->foe->get_double(_elid,&evib_ptr,_sale->e_vib);

    double *eejt_ptr;
    _sale->foe->get_double(_elid,&eejt_ptr,_sale->e_ejt);

    for(int i=0;i<4;++i)
    {
        if(lbru_vol[i] >=0.) continue;
        double flux_mass = 0.0;

        /*
         * variables advect by volume
         * [0,3*nmat-1]{volume, density, energy} for every materials
         */
        for(int cutid=0;cutid<nmat;++cutid)
        {
            bphi_ptr[i][cutid + 0*nmat] = flux_vol[cutid][i];
            bphi_ptr[i][cutid + 1*nmat] = mden[cutid];
            bphi_ptr[i][cutid + 2*nmat] = meng[cutid];
            flux_mass += mden[cutid]*flux_vol[cutid][i];
        }
        /*
         * variables advect by mass
         * length, offset             description
         * -----------------------------------------------
         * 2     ,[3*nmat, 3*nmat+1]  velocity_{r,z}
         * 4     ,[3*nmat+2,3*nmat+5] stress_{rr,rz,zz,theta}
         * 1     ,[3*nmat+6]          damage
         * 1     ,[3*nmat+7]          total plastic strain, tps
         * 1     ,[3*nmat+8]          vibration energy
         * 1     ,[3*nmat+9]          ejecta mass fraction
         * 1  ,[3*namt+10, *]        cold volume strain
         * nmat  ,[3*namt+11, *]      porosity of materials
         * -----------------------------------------------
         * 11+nmat summary : phi_length - 3*nmat
        */
        flux_mass = 1.0; // flux_mass is removed from bphi
        vec_copy(evel,bphi_ptr[i]+3*nmat,2);
        scaler_move(bphi_ptr[i]+3*nmat+2,4,este,flux_mass);
        bphi_ptr[i][3*nmat + 6] = flux_mass * (*edam);
        bphi_ptr[i][3*nmat + 7] = flux_mass * (*etps);
        bphi_ptr[i][3*nmat + 8] = flux_mass * (*evib_ptr);
        bphi_ptr[i][3*nmat + 9] = flux_mass * (*eejt_ptr);

        if(_sale->porosity_model)
        {
            bphi_ptr[i][3*nmat + 10] = flux_mass * (*ecvs);
            for(int cutid=0;cutid<nmat;++cutid)
            {
                bphi_ptr[i][3*nmat + 11 + cutid] = mpty[cutid];
            }
        }
    }
}

double algebraic_flux_ratio(double VOFU, double VOFD, double VOFA, double CourantD,
                            double bnormal[],double grad_vof[],
                            double (*vofAlgorithm)(double,double,double))
{
    double tmp = VOFA - VOFU;
    double nVOFD = (VOFD - VOFU) / tmp;

    double boundaryVOF = 0.0;

    if(fabs(tmp) < TOLVOF)
    {
        boundaryVOF = VOFD;
    } else if(nVOFD < TOLVOF || nVOFD > 1.0 - TOLVOF)
    {
        boundaryVOF = VOFD;
    } else
    {
        double CosTf;
        CosTf = vec_dot(grad_vof, bnormal,2) / vec_len(grad_vof,2);
        CosTf = CosTf * CosTf;
        double Gf, nVOFf;
        Gf = Min(CISSAMK1 * CosTf * CosTf, 1.0);
        nVOFf = vofAlgorithm(nVOFD, Gf, CourantD);
        boundaryVOF = VOFU + nVOFf * (VOFA - VOFU);
    }
    return boundaryVOF;
}


/*
 * MSTACS, CICSAM, STACS, MHRIC used in algebraic reconstruction
 * just copy from SALEc
 */
double vofMSTACS(double nVOFD, double Gf, double CourantD)
{
    double nVOFCDSMSTACS,nVOFfSTOIC;
    if(CourantD < 1.0/3.0)
        nVOFCDSMSTACS = Min(nVOFD/CourantD,1.0);
    else
        nVOFCDSMSTACS = Min(3.0*nVOFD,1.0);

    if(nVOFD < 1.0/5.0)
        nVOFfSTOIC = 3.0 * nVOFD;
    else if(nVOFD < 1.0/2.0)
        nVOFfSTOIC = (1.0 + nVOFD)/2.0;
    else
        nVOFfSTOIC = Min((6.0*nVOFD + 3.0)/8.0,1.0);
    return Gf*nVOFCDSMSTACS + (1.0 - Gf)*nVOFfSTOIC;
}

double vofCICSAM(double nVOFD, double Gf, double CourantD)
{
    double nVOFCBC, nVOFUQ;
    nVOFCBC = Min(nVOFD / CourantD, 1.0);
    nVOFUQ  = Min(CourantD * nVOFD + (1.0 - CourantD) * (6.0 * nVOFD + 3.0) / 8.0, nVOFCBC);
    return Gf * nVOFCBC + (1.0 - Gf) * nVOFUQ;
}

double vofSTACS(double nVOFD, double Gf)
{
    double nVOFfSUP,nVOFfSTOIC;
    if(nVOFD < 1.0/3.0)
        nVOFfSUP = 2.0*nVOFD;
    else if(nVOFD < 1.0/2.0)
        nVOFfSUP = (1.0 + nVOFD)/2.0;
    else
        nVOFfSUP = Min(3.0*nVOFD/2.0,1.0);
    if(nVOFD < 1.0/5.0)
        nVOFfSTOIC = 3.0 * nVOFD;
    else if(nVOFD < 1.0/2.0)
        nVOFfSTOIC = (1.0 + nVOFD)/2.0;
    else
        nVOFfSTOIC = Min((6.0*nVOFD + 3.0)/8.0,1.0);

    return Gf*nVOFfSUP + (1.0 - Gf)*nVOFfSTOIC;
}

double vofMHRIC(double nVOFD, double Gf)
{
    // MHRIC
    double nVOFs1,nVOFs2;
    nVOFs1 = Min(2.0*nVOFD,1.0);
    nVOFs2 = Min((6.0*nVOFD + 3.0)/8.0,nVOFs1);
    return Gf*nVOFs1 + (1.0 - Gf)*nVOFs2;
}

double geometric_flux_ratio(double  _p1, double _p2)
{
    double a = _p1, b=_p2;
    if(b > a)
    {
        a = _p2; b=_p1;
    }

    if(b >= 0.0)
    {
        return 1.0;
    } else if(a <= 0.0)
    {
        return 0.0;
    } else
    {
        return a / (a  + fabs(b));
    }
}

void calculate_flux(sale2d_var * _sale, node_type _nt)
{
    switch (_sale->flux_opt)
    {
        case material_priority:
            calculate_flux_primary(_sale,_nt);
            break;
        case mixed_strategy:
            calculate_flux_mixed(_sale,_nt);
            break;
        default:
            calculate_flux_upwind(_sale,_nt);
    }
}

double get_flux_volcyl(const double * base_a,const double * base_b, const double * dis_a, const double * dis_b)
{
    double xa = base_a[0], ya=base_a[1], xb=base_b[0], yb=base_b[1];
    double x1 = -1.0*dis_a[0], y1= -1.0*dis_a[1], x2= -1.0*dis_b[0], y2= -1.0*dis_b[1];

    double vol = x1*(ya - yb - y2) + x2*(ya - yb + y1) + (y1 + y2)*(xb - xa);
    return vol*fabs(2.0*xa + 2.0*xb + x1 + x2)/8.0;
}

double cut_length_l(double _vf,const double * _p)
{
    /*
     * b << a, b/a approx 0.
     */
    double a = _p[0],b=_p[1];
    double dstar = 0.0;
    if(_vf <= 0.0)
        dstar = 0.0;
    else if(_vf <= 1.0)
        dstar = _vf*a;
    else
        dstar = a;
    return dstar + b;
}

double cut_length_s(double _vf,const double * _p)
{
    /*
     * b < a,
     */
    double a = _p[0],b=_p[1];
    double dstar;
    double crf1 = 0.5*b/a;
    if(_vf <= crf1)
        dstar = sqrt(2.0*a*b*_vf);
    else if(_vf <= 1.0 - crf1)
        dstar = 0.5*b + a*_vf;
    else if(_vf <= 1.0)
        dstar = a + b - sqrt(2.0*a*b*(1.0-_vf));
    else
        dstar = a + b;
    return dstar;
}

double cut_length(double _vf, double * _p)
{
    if(_p[0] < _p[1])
    {
        double tmp = _p[0];
        _p[0] = _p[1];
        _p[1] = tmp;
    }


    if(1.0e-5 * _p[0] > _p[1])
    {
        return cut_length_l(_vf,_p);
    } else
    {
        return cut_length_s(_vf,_p);
    }
}