//
// Created by huach on 14/11/2022.
//

#include "sale2d.h"
void init_cycle_control(InputFile * ifp, cycle_control * _ccl)
{
    _ccl->start_time = GetValueD(ifp,"cycle.time0","0.0");
    _ccl->end_time   = GetValueD(ifp,"cycle.time1","100.0");
    _ccl->max_dt     = GetValueD(ifp,"cycle.maxdt","1.0");
    _ccl->dt = GetValueD(ifp,"cycle.dt0","0.0");
    _ccl->max_step  = GetValueI(ifp,"cycle.maxstep","200");
    _ccl->local_dt = _ccl->dt;

    _ccl->i_step = 0;
    _ccl->k_out_step = 0;
    _ccl->i_out_step = 0;
    _ccl->step_time = _ccl->start_time;

    _ccl->init_step = GetValueI(ifp,"output.init","0");
    _ccl->out_time  = GetValueD(ifp,"output.time","1.0");
    _ccl->out_cycle = GetValueI(ifp,"output.out_cycle","-1");

    GetValueS(ifp, "output.prefix", _ccl->prefix, "ParaTest");

    _ccl->out_tracer = 0;
    _ccl->cycle_continue = 1;
    _ccl->out_log = 0;
    _ccl->out_vts = 0;
    _ccl->out_init = 0;
    _ccl->out_tracer_txt = 0;

    MPI_Comm_rank(MPI_COMM_WORLD,&(_ccl->rank));
    int npgy = GetValueI(ifp,"processor.npgy","2");
    _ccl->x_rank = _ccl->rank / npgy;
    _ccl->y_rank = _ccl->rank % npgy;

#ifdef DEBUG_MODE
    int name_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &name_len);
    fprintf(stderr,"Rank %d running on host %s has PID %d\n", _ccl->rank, hostname, getpid());
#endif

    // prepare directory for vts/vtp output
    // create a FILE to write progress info
    if(_ccl->rank == 1)
    {
        make_empty_dir("./vts");
        make_empty_dir("./vtp");
        make_empty_dir("./txt");
        make_empty_dir("./post");

        _ccl->log = fopen("timestep.txt","w");
        fputs("  STEP,           DT,            T, output,     id\n",_ccl->log);
        FILE * node_info_fp = fopen("node_info.txt","w");
        put_cpuinfo(node_info_fp);
        fclose(node_info_fp);
    }
    else
    {
        _ccl->log = NULL;
    }
}

void clean_cycle_control(cycle_control * _ccl)
{
    if(NULL != _ccl->log) fclose(_ccl->log);
}

void update_cycle_control(cycle_control * _ccl, FILE * fp)
{
    /*
     * check/update cycle control
     */

    // update cycle time
    _ccl->step_time += _ccl->dt;
    _ccl->i_step ++;
    // check the stoping conditions

    _ccl->cycle_continue = 1;

    if(_ccl->i_step > _ccl->max_step)
    {
        _ccl->cycle_continue = 0;
    }

    if(_ccl->step_time > _ccl->end_time + 2*_ccl->max_dt)
    {
        _ccl->cycle_continue = 0;
    }

    // check output state
    // if(fabs(_ccl->step_time - _ccl->out_time*_ccl->k_out_step + _ccl->out_time*_ccl->i_out_step) < 0.5*_ccl->out_time)
    if(_ccl->step_time >= _ccl->out_time*(_ccl->k_out_step - _ccl->i_out_step))
    {
        _ccl->k_out_step ++;
        _ccl->out_vts = 1;
    } else if(_ccl->i_step <= _ccl->init_step)
    {
        _ccl->k_out_step++;
        _ccl->i_out_step++;
        _ccl->out_vts = 1;
    } else if(_ccl->out_cycle > 0 && _ccl->i_step % _ccl->out_cycle == 0)
    {
        _ccl->k_out_step++;
        _ccl->i_out_step++;
        _ccl->out_vts = 1;
    } else
    {
        _ccl->out_vts = 0;
    }

    // check tracer output state
    _ccl->out_tracer = 0;

    // check log output state
    if((_ccl->i_step % 100 == 0) || (_ccl->i_step <= _ccl->init_step))
    {
        _ccl->out_log = 1;
    } else
    {
        _ccl->out_log = 0;
    }


    if(_ccl->i_step % 100 == 1)
    {
        _ccl->sync_vel = 1;
    } else
    {
        _ccl->sync_vel = 0;
    }

}

void sale2d_write_log(FILE * fp, cycle_control * _ccl)
{
    //write information into log file
    if(_ccl->rank == 1)
    {
        if(_ccl->out_vts)
        {
            fprintf(fp,"[*%5d|%4d] dt = %10.5es, time = %10.5es\n",
                    _ccl->i_step,_ccl->k_out_step-1,_ccl->local_dt, _ccl->step_time);
        } else if(_ccl->out_log)
        {
            fprintf(fp,"[ %5d|    ] dt = %10.5es, time = %10.5es\n",
                    _ccl->i_step,_ccl->local_dt, _ccl->step_time);
        }
    }

    if(NULL != _ccl->log)
    {
        if(_ccl->out_vts)
        {
            fprintf(_ccl->log,"%6d, %10.5es, %10.5es, %6d, %6d\n",_ccl->i_step, _ccl->local_dt, _ccl->step_time,
                    _ccl->out_vts,_ccl->k_out_step-1);
        }
        else
        {
            fprintf(_ccl->log,"%6d, %10.5es, %10.5es, %6d,       \n",_ccl->i_step, _ccl->local_dt, _ccl->step_time,
                    _ccl->out_vts);
        }

    }
}


void init_sale2d_var(sale2d_var * _sale, mesh2d_info * _minfo)
{
    // only responsible for allocate memory for parallel environment
    // other codes have been moved to init_sale_mat(eos)
    // all the flags should be set before allocate memory for sale2d_var
    int nmat = _sale->nmat;
    // allocate space for operation collection
    field_opt * local_opt = (field_opt *) malloc(sizeof(field_opt)*5);
    _sale->foe  = local_opt + 0;
    _sale->fov  = local_opt + 1;
    _sale->fobx = local_opt + 2;
    _sale->foby = local_opt + 3;
    _sale->foi  = local_opt + 4;
    // set operation for element, vertex, and bx by
    // init_field_opt_2d_default is designed for element
    init_field_opt_2d_default(_sale->foe);
    init_field_opt_2d_default(_sale->fov);
    _sale->fov->lid2nty = node_type_2d_default_vertex;
    init_field_opt_2d_default(_sale->fobx);
    init_field_opt_2d_default(_sale->foby);
    _sale->fobx->lid2nty = node_type_2d_default_vertex;
    _sale->foby->lid2nty = node_type_2d_default_vertex;

    // allocate space for variables
    int n_fv_num = 64;
    int field_ptr_offset = 0;

    field_var * local_fv = (field_var *) malloc(sizeof(field_var)*n_fv_num);
    if(NULL == local_fv)
    {
        fprintf(stderr,"fail to allocate local field \n");
    }
    // register the variable for parallel send/recv
    // region, (7 variables), variables on elements related to velocity
    //    field_ptr_offset ++;
    _sale->v_vel = local_fv + field_ptr_offset++;
    _sale->v_bf  = local_fv + field_ptr_offset++;
    _sale->v_acc = local_fv + field_ptr_offset++;
    _sale->e_vel = local_fv + field_ptr_offset++;
    _sale->v_dis = local_fv + field_ptr_offset++;
    _sale->e_div_vel = local_fv + field_ptr_offset++;
    _sale->e_grad_vel = local_fv + field_ptr_offset++;
    _sale->e_ste = local_fv + field_ptr_offset++;
    _sale->e_mu = local_fv + field_ptr_offset++;
    _sale->e_dvel = local_fv + field_ptr_offset++;

    field_var * dup_v = _sale->v_vel;
    field_var * dup_e = _sale->e_div_vel;


    init_field_2d(&(field_2d){.npgx=_minfo->npgx,.npgy=_minfo->npgy,.xpad=_minfo->xpad,.ypad=_minfo->ypad,
                              .nx=_minfo->nx+1,.ny=_minfo->ny+1,.bytes_per_unit=2*sizeof(double)/sizeof(char),
                              .padding=vertex_pad},
                  _sale->v_vel,_sale->fov);

    duplicate_field_2d(dup_v,_sale->v_bf ,_sale->fov);
    duplicate_field_2d(dup_v,_sale->v_acc,_sale->fov);
    duplicate_field_2d(dup_v,_sale->v_dis,_sale->fov);

    init_field_2d(&(field_2d){.npgx=_minfo->npgx,.npgy=_minfo->npgy,.xpad=_minfo->xpad,.ypad=_minfo->ypad,
                          .nx=_minfo->nx,.ny=_minfo->ny,.bytes_per_unit=sizeof(double)/sizeof(char),
                          .padding=element_pad},
                  _sale->e_div_vel,_sale->foe);

    duplicate_field_2d_n(dup_e,_sale->e_grad_vel,_sale->foe,5);
    duplicate_field_2d_n(dup_e,_sale->e_vel,_sale->foe,2);
    duplicate_field_2d_n(dup_e,_sale->e_ste,_sale->foe,4);
    duplicate_field_2d(dup_e,_sale->e_mu,_sale->foe);
    duplicate_field_2d_n(dup_e,_sale->e_dvel,_sale->foe,2);


    // region, (7 variables), state variables (one double)
    _sale->e_pre = local_fv + field_ptr_offset++;
    _sale->e_tem = local_fv + field_ptr_offset++;
    _sale->e_den = local_fv + field_ptr_offset++;
    _sale->e_eng = local_fv + field_ptr_offset++;
    _sale->e_csd = local_fv + field_ptr_offset++;
    _sale->e_dam = local_fv + field_ptr_offset++;
    _sale->e_cvs = local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_pre,_sale->foe);
    duplicate_field_2d(dup_e,_sale->e_tem,_sale->foe);
    duplicate_field_2d(dup_e,_sale->e_den,_sale->foe);
    duplicate_field_2d(dup_e,_sale->e_eng,_sale->foe);
    duplicate_field_2d(dup_e,_sale->e_csd,_sale->foe);
    duplicate_field_2d(dup_e,_sale->e_dam,_sale->foe);
    duplicate_field_2d(dup_e,_sale->e_cvs,_sale->foe);

    // region, (5 variables), state variables for different materials (nmat double)
    _sale->m_vof = local_fv + field_ptr_offset++;
    _sale->m_den = local_fv + field_ptr_offset++;
    _sale->m_eng = local_fv + field_ptr_offset++;
    _sale->e_wpt = local_fv + field_ptr_offset++;
    _sale->m_pty = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->m_vof,_sale->foe,nmat);
    duplicate_field_2d_n(dup_e,_sale->m_den,_sale->foe,nmat);
    duplicate_field_2d_n(dup_e,_sale->m_eng,_sale->foe,nmat);
    duplicate_field_2d_n(dup_e,_sale->e_wpt,_sale->foe,4);
    duplicate_field_2d_n(dup_e,_sale->m_pty,_sale->foe,nmat);

    // region, (5 variable), variables defined before loop
    _sale->v_pos  = local_fv + field_ptr_offset++;
    _sale->e_pos  = local_fv + field_ptr_offset++;
    _sale->e_vol  = local_fv + field_ptr_offset++;
    _sale->e_suf  = local_fv + field_ptr_offset++;
    _sale->e_grad = local_fv + field_ptr_offset++;
    _sale->e_subvol = local_fv + field_ptr_offset++;

    duplicate_field_2d(dup_v,_sale->v_pos,_sale->fov);
    duplicate_field_2d(dup_e,_sale->e_vol,_sale->foe);
    duplicate_field_2d_n(dup_e,_sale->e_pos,_sale->foe,4);
    duplicate_field_2d_n(dup_e,_sale->e_suf,_sale->foe,8);
    duplicate_field_2d_n(dup_e,_sale->e_grad,_sale->foe,12);
    duplicate_field_2d_n(dup_e,_sale->e_subvol,_sale->foe,8);
    // region, (2 variables), VARIABLES ON BOUNDARY

    _sale->bx_phi = local_fv + field_ptr_offset++;
    _sale->by_phi = local_fv + field_ptr_offset++;

    // consist of phi_length
    // 3*nmat = (vol,density,energy) of materials
    // velocity
    //
    // ..
    // damage
    // ejecta
    // porosity(nmat)
    // cold volumetric strain
    int phi_length = 3*nmat + 8 + 1 + 1;
    if(_sale->porosity_model){
        // porosity is opened
        phi_length = 3*nmat + 8 + 1 + 1 + (nmat + 1);
    }
    _sale->phi_length = phi_length;

    init_field_2d(&(field_2d){.npgx=_minfo->npgx,.npgy=_minfo->npgy,.xpad=_minfo->xpad,.ypad=_minfo->ypad,
                          .nx=_minfo->nx+1,.ny=_minfo->ny,.bytes_per_unit=phi_length*sizeof(double)/sizeof(char),
                          .padding=bx_pad},
                  _sale->bx_phi, _sale->fobx);
    init_field_2d(&(field_2d){.npgx=_minfo->npgx,.npgy=_minfo->npgy,.xpad=_minfo->xpad,.ypad=_minfo->ypad,
                          .nx=_minfo->nx,.ny=_minfo->ny+1,.bytes_per_unit=phi_length*sizeof(double)/sizeof(char),
                          .padding=by_pad},
                  _sale->by_phi, _sale->foby);
    // normal vector of boundary
    _sale->bx_normal = local_fv + field_ptr_offset++;
    _sale->by_normal = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(_sale->bx_phi,_sale->bx_normal,_sale->fobx,-2);
    duplicate_field_2d_n(_sale->by_phi,_sale->by_normal,_sale->foby,-2);

    // region, (2 variabls), variables for flux
    _sale->v_vof = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_v,_sale->v_vof,_sale->fov,-nmat);
    _sale->e_vof = _sale->m_vof;
    _sale->v_mass = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_v,_sale->v_mass,_sale->fov,1);

    //region, variables for debug
    _sale->e_debug = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_debug,_sale->foe,8);
    _sale->v_debug = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_v,_sale->v_debug,_sale->fov,4);

    //region, variables for stability
    _sale->e_q = local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_q,_sale->foe);

    _sale->e_strength = local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_strength,_sale->foe);

    _sale->e_cnd = local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_cnd,_sale->foe);

    _sale->e_tps =  local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_tps,_sale->foe);

    _sale->e_ste2s = local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_ste2s,_sale->foe);

    _sale->e_sta2s = local_fv + field_ptr_offset++;
    duplicate_field_2d(dup_e,_sale->e_sta2s,_sale->foe);

    _sale->v_vol = local_fv + field_ptr_offset++;
    _sale->v_subvol = local_fv + field_ptr_offset++;
    _sale->v_laplace = local_fv + field_ptr_offset++;

    duplicate_field_2d_n(dup_v,_sale->v_vol,_sale->fov,-5); // when the last args is negative, it represents abs(n)*double
    duplicate_field_2d_n(dup_v,_sale->v_subvol,_sale->fov,2);
    duplicate_field_2d_n(dup_v,_sale->v_laplace,_sale->fov,-9);

    _sale->v_suf = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_v,_sale->v_suf,_sale->fov,4);

    _sale->e_alpha = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_alpha,_sale->foe,2);

    _sale->e_beta = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_beta,_sale->foe,4);

    _sale->e_sigma = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_sigma,_sale->foe,2);

    _sale->v_alpha = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_v,_sale->v_alpha,_sale->fov,1);

    _sale->v_beta = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_v,_sale->v_beta,_sale->fov,1);

    _sale->e_dampref = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_dampref,_sale->foe,3);

    _sale->e_vib = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_vib,_sale->foe,1);

    _sale->e_ejt = local_fv + field_ptr_offset++;
    duplicate_field_2d_n(dup_e,_sale->e_ejt,_sale->foe,1);

    _sale->y_profile = local_fv + field_ptr_offset++;
    init_index_2d(&(field_2d){.npgx=_minfo->npgx,.npgy=_minfo->npgy,.xpad=_minfo->xpad,.ypad=_minfo->ypad,
            .nx=_minfo->nx,.ny=-_minfo->ny,.bytes_per_unit=sizeof(double)/sizeof(char),
            .padding=element_pad},_sale->y_profile);
    init_index_opt_2d_default(_sale->foi);
    // record the fv_pointer allocated in sale2d
    _sale->fv_counter = field_ptr_offset;

    if(1 == dup_e->rank)
    {
        fprintf(stderr,"%d field_var have been initialize.\n",field_ptr_offset);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // check whether rank is consistent
    int rank_consistent = 0;
    for(int fv_off=0;fv_off<field_ptr_offset;++fv_off)
    {
        if(local_fv[fv_off].rank != local_fv[0].rank) rank_consistent += 1;
    }
    if(0 == rank_consistent)
    {
        fprintf(stderr,"rank(%d) consistent check pass.\n",local_fv[0].rank);
    }
    else
    {
        fprintf(stderr,"rank(%d) consistent check fail.\n",local_fv[0].rank);
    }

}


void clean_sale_var(sale2d_var * _sale)
{
    for(int k=0;k<_sale->fv_counter;++k)
    {
        clean_field(_sale->v_vel + k);
    }

    free(_sale->v_vel);
    free(_sale->foe);
}

void init_sale_mat(sale2d_var * _sale, mesh2d_info * _minfo)
{
    /*
     * old version is init_sale_eos
     * some code of init_sale2d_var have been inserted here.
     */

    // some variables copy from minfo
    _sale->nmat = _minfo->n_materials;
    _sale->avis2 = _minfo->avis2;
    _sale->avis1 = _minfo->avis1;
    _sale->anc   = _minfo->anc;
    _sale->density_cutoff = _minfo->density_cutoff;

    // set args related to ejecta
    _sale->threshold_z = _minfo->ejecta_threshold * _minfo->pinfo.radiu;
    _sale->flow_resistance = _minfo->ejecta_resistance;
    _sale->flow_friction = _minfo->ejecta_friction;

    if (_minfo->velocity_max > 0.)
    {
        _sale->velocity_max = _minfo->velocity_max;
    } else
    {
        _sale->velocity_max = fabs(_minfo->velocity_max) * vec_len(_minfo->pinfo.velocity, 2);
    }

    if(_minfo->velocity_min > 0.)
    {
        _sale->velocity_min = _minfo->velocity_min;
    }
    else
    {
        _sale->velocity_min = fabs(_minfo->velocity_min) * vec_len(_minfo->pinfo.velocity,2);
    }

    if(_minfo->acceleration_min > 0.)
    {
        _sale->acceleration_min = _minfo->acceleration_min;
    }
    else
    {
        _sale->acceleration_min = fabs(_minfo->acceleration_min) * vec_len(_minfo->gravity,2);
    }

    _sale->partpres = _minfo->partpres;
    _sale->damp_time = _minfo->damp_t;

    // set the mixed strategy as default
    _sale->flux_opt = FLUXSTRATEGY;
    // set strength mod
    _sale->strength_mod = STRENGTHMOD;

    int nmat = _sale->nmat;
    eos_table * etb = (eos_table *) malloc(sizeof(eos_table)*nmat);
    _sale->etb = etb;

    // check whether acfl is opened
    _sale->acoustic_fluid = 0;
    // close the tensile failure default
    _sale->tensile_failure = 0;
    // close the porosity model default
    // check whether the porosity model is opened
    // default = 0 (turn off)
    _sale->porosity_model = 0;
    for(int matid=1;matid<nmat;++matid) // 0 is vacuum, skip it
    {
        FILE * fp = fopen(_minfo->ename[matid],"r");
        if(NULL == fp)
        {
            fprintf(stdout,"cannot open eos at %s\n", _minfo->ename[matid]);
            exit(0);
        }
        etb[matid].type = _minfo->etype[matid];
        etb[matid].table= NULL;
        etb[matid].ref  = _minfo->ref + matid;

        if(etb[matid].ref->Toff > 0.0)
        {
            _sale->acoustic_fluid = 1;
        }
        load_eos_table(etb+matid,fp);

        if(strcasecmp(_sale->etb[matid].ref->porosityfunc,"None") != 0) _sale->porosity_model++;
    }
}

void init_sale2d_other(sale2d_var * _sale, InputFile * ifp)
{
    // most parameters are set through mesh2d_info
    // but this struct is redundancy, this part directly set some parameters from InputFile
    // this function complete some work of load_mesh_info + init_sale_mat

    // set tracers lubrication, default 0 represent no lubrication tracers on the target surface
    _sale->tracer_lubrication = GetValueI(ifp,"tracer.lubrication","0");
    // value of tracer_lub_y0/y1 is set when the position of vertexes are set (init_sale_pos)
    _sale->tracer_lub_y0   = GetValueD(ifp,"tracer.lub_y0","0.");
    _sale->tracer_lub_frac = GetValueD(ifp,"tracer.lub_frac","0.");

    // flags for tensile failure
    _sale->tensile_failure = GetValueI(ifp,"numerical.tensile_failure","0");
}


void update_velocity(sale2d_var * _sale, node_type _nt,cycle_control * _ccl)
{
    // update the velocity accord to the accerlation /pressure/strain
    proc_info_2d* proc_acc = (proc_info_2d*) _sale->v_acc->proc_self;
    int nx = proc_acc->nx,ny=proc_acc->ny,myrank=_sale->v_acc->rank;

    // update the displacement and velocity

    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;++j)
        {
            fv_id_2d lid = {.x=k,.y=j,0,0};
            node_type lnt = useless_section;
            _sale->fov->lid2nty(myrank,&lid,&lnt,_sale->v_acc);

            if(!(_nt & lnt)) continue;
            if(lnt == receive_section) continue;
            char * acc_cptr, *vel_cptr, *dis_cptr;
            _sale->foe->get_data2(&lid,&acc_cptr,_sale->v_acc);
            _sale->foe->get_data2(&lid,&vel_cptr,_sale->v_vel);
            _sale->foe->get_data2(&lid,&dis_cptr,_sale->v_dis);
            vec_zero((double *)dis_cptr,2);
            scaler_add((double *)dis_cptr,2,(double *)vel_cptr,0.5*_ccl->dt);
            scaler_add((double *)vel_cptr,2,(double *)acc_cptr,_ccl->dt);
            scaler_add((double *)dis_cptr,2,(double *)vel_cptr,0.5*_ccl->dt);
        }
    }
}


void update_energy(sale2d_var * _sale, node_type _nt, double dt)
{
    /*
     * update the energy in elements (internal energy I and kinetic energy)
     *
     */

    proc_info_2d * proc_evel = (proc_info_2d *)_sale->e_vel->proc_self;
    int nx = proc_evel->nx, ny = proc_evel->ny;
    int myrank = _sale->e_vel->rank;
    int nmat = _sale->nmat;

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid,&lnt,_sale->e_vel);
            if(!(lnt&_nt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};

            double * vvel_ptr[4];
            double vvelx[4] = {0.},vvely[4]={0.};

            for(int i=0;i<4;++i)
            {
                _sale->fov->get_double(vlid[i],&(vvel_ptr[i]),_sale->v_vel);
                vvelx[i] = vvel_ptr[i][0];
                vvely[i] = vvel_ptr[i][1];
            }

            double * egrad_ptr;
            _sale->foe->get_double(&lid,&egrad_ptr,_sale->e_grad);

            /*
             * update energy
             * VARIABLES need to be synchronized in advance:
             *      e_div_vel, e_grad_vel, e_stress, e_pre
             * VARIABLES overwritten:
             *      m_eng
             */


            // skip the empty/non-full cells
            double *eden_ptr,*evof_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&evof_ptr,_sale->e_vof);
            double *mpty_ptr,*ecvs_ptr;
            _sale->foe->get_double(&lid,&mpty_ptr,_sale->m_pty);
            _sale->foe->get_double(&lid,&ecvs_ptr,_sale->e_cvs);
            if(evof_ptr[VACUUM] > TOLVOF || eden_ptr[0] < TOLRHO)
            {
                continue;
            }

            // pressure and divergence of velocity
            double *epre_ptr, *eq_ptr;
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);
            _sale->foe->get_double(&lid,&eq_ptr,_sale->e_q);
            double *edvel_ptr;
            _sale->foe->get_double(&lid,&edvel_ptr,_sale->e_div_vel);

            // strain and gradvel
            double *egradvel_ptr,*este_ptr;
            _sale->foe->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);
            _sale->foe->get_double(&lid,&este_ptr,_sale->e_ste);

            double *meng_ptr,*mden_ptr;
            _sale->foe->get_double(&lid,&meng_ptr,_sale->m_eng);
            _sale->foe->get_double(&lid,&mden_ptr,_sale->m_den);

            double *evel_ptr;
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);

            double *edebug_ptr;
            _sale->foe->get_double(&lid,&edebug_ptr,_sale->e_debug);

            double *esigma_ptr;
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);
            double local_pre = epre_ptr[0] + eq_ptr[0];
            double dE_pre = - local_pre * (*edvel_ptr);
            double dE_ste = este_ptr[0] * egradvel_ptr[0]
                    + este_ptr[1]*(egradvel_ptr[1] + egradvel_ptr[2])
                    + este_ptr[2]*egradvel_ptr[3] + este_ptr[3]*egradvel_ptr[4];
            double egamfac = 0.;
            for(int matid=1;matid<nmat;++matid)
            {
                // compute the average gamma/csd**2
                if(evof_ptr[matid] > TOLVOF)
                    egamfac += evof_ptr[matid] * _sale->etb[matid].ref->gamfrac/mden_ptr[matid];
            }
            double dcvs = (*edvel_ptr) - egamfac*(dE_pre + dE_ste); // the cold volume strain

#ifdef UPDATE_ENERGY_AVG_TC
            // may be a more precious scheme from SALE manual, se the u_tc and v_tc
            // use averge velocity during eulerian step
            // when courant number < 0.2, this scheme give same results as previous
            double t_egradvel_ptr[5] = {0.};
            t_egradvel_ptr[0] = vec_dot(egrad_ptr+0,vvelx,4); /*partial v_r / partial r */
            t_egradvel_ptr[2] = vec_dot(egrad_ptr+0,vvely,4); /*partial v_z / partial r */
            t_egradvel_ptr[1] = vec_dot(egrad_ptr+4,vvelx,4); /*partial v_r / partial z */
            t_egradvel_ptr[3] = vec_dot(egrad_ptr+4,vvely,4); /*partial v_z / partial z */
            t_egradvel_ptr[4] = vec_dot(egrad_ptr+8,vvelx,4); /*partial v_r / r */
            double t_divvel = t_egradvel_ptr[0] + t_egradvel_ptr[3] + t_egradvel_ptr[4];
            double t_dE_pre = - local_pre * (t_divvel);
            double t_dE_ste = este_ptr[0] * t_egradvel_ptr[0]
                            + este_ptr[1]*(t_egradvel_ptr[1] + t_egradvel_ptr[2])
                            + este_ptr[2]*t_egradvel_ptr[3] + este_ptr[3]*t_egradvel_ptr[4];
            dE_pre = 0.5*(dE_pre + t_dE_pre);
            dE_ste = 0.5*(dE_ste + t_dE_ste);
            dcvs = t_divvel - egamfac*(dE_pre + dE_ste);
#endif

            for(int matid=1;matid<nmat;++matid)
            {
                if(evof_ptr[matid] <= TOLVOF || mden_ptr[matid] <= TOLRHO)
                {
                    meng_ptr[matid] = 0.;
                    if(_sale->porosity_model) mpty_ptr[matid] = 1.;
                } else
                {
                    // double local_de_factor = mden_ptr[matid]/eden_ptr[0];
                    double local_de_factor = 1.0;
                    meng_ptr[matid] = Max(0., meng_ptr[matid] + (dE_pre + dE_ste)*dt*local_de_factor);
                    // internal energy should be always positive
                    if(_sale->porosity_model){
                        double (*porfunc)(double*,const state_reference *) = (double (*)(double*,const state_reference *))_sale->etb[matid].ref->porosity_increment_func;
                        double delta_pty = 0.;
#ifdef EXCLUDE_PML_POROSITY
                        if(esigma_ptr[X] + esigma_ptr[Y] <= 0.)
                        {
                            delta_pty = porfunc((double []){mpty_ptr[matid],*ecvs_ptr,*epre_ptr,dcvs, evel_ptr[0], evel_ptr[1]},_sale->etb[matid].ref);
                        }
#else
                        delta_pty = porfunc((double []){mpty_ptr[matid],*ecvs_ptr,*epre_ptr,dcvs, evel_ptr[0], evel_ptr[1]},_sale->etb[matid].ref);
#endif
                        mpty_ptr[matid] += delta_pty*dt;
                        mpty_ptr[matid] = Wind(mpty_ptr[matid],1.0,_sale->etb[matid].ref->alphamax2);
                    }
                }
            }
            *ecvs_ptr = Wind(*ecvs_ptr+dcvs*dt,-10.0,10.0);

        }
    }
}

void update_e_den(sale2d_var * _sale, node_type _nt)
{
    /*
     * update the density of whole element after every cycle
     * VARIABLES must be synchronous
     *      m_den, m_vof
     * VARIABLE overwritten
     *      e_den
     */
    proc_info_2d* proc_e_den = (proc_info_2d*) _sale->e_den->proc_self;
    int nx = proc_e_den->nx,ny=proc_e_den->ny,myrank=_sale->e_den->rank;
    int nmat = _sale->nmat;

    for(int k=0;k<nx;++k) {
        for (int j = 0; j < ny; ++j) {
            fv_id_2d lid = {.x=k, .y=j};
            node_type lnt = useless_section;
            _sale->foe->lid2nty(myrank, &lid, &lnt, _sale->e_vel);
            if (!(lnt & _nt)) continue;

            double *eden_ptr, *mden_ptr, *mvof_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&mden_ptr,_sale->m_den);
            _sale->foe->get_double(&lid,&mvof_ptr,_sale->m_vof);

            *eden_ptr = 0.0;
            for(int matid = 1;matid < nmat; ++matid)
            {
                eden_ptr[matid] += mden_ptr[matid] * mvof_ptr[matid];
            }

        }
    }

}


void update_v_vof(sale2d_var * _sale, node_type _nt)
{
    /*
     * copy the vof from cells, set v_vof accord to e_vof
     * VARIABLES need to be synchronized in advance:
     *    e_vof
     * VARIABLES overwritten AFTER:
     *    v_vof
     *
     * change loop sequence according to LID2index
     */
    proc_info_2d* proc_v_vof = (proc_info_2d*) _sale->v_vof->proc_self;
    int nx = proc_v_vof->nx,ny=proc_v_vof->ny,myrank=_sale->v_vof->rank;
    int nmat = _sale->nmat;

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _sale->fov->get_nty(&lid,&lnt,_sale->v_vof);

            if(!(_nt&lnt)) continue;

            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};
            char *vof_cptr;
            _sale->fov->get_data2(&lid,&vof_cptr,_sale->v_vof);
            vec_zero((double *) vof_cptr,nmat);

            for(int i=0;i<4;i++)
            {
                char * evof_cptr;
                _sale->foe->get_data2(elid[i],&evof_cptr,_sale->e_vof);
                scaler_add((double *)vof_cptr,nmat,(double *)(evof_cptr),0.25);
            }
        }
    }
}


void clean_fv(field_var * _var)
{
    int n_component = _var->bytes_per_unit / sizeof(double);
    double * d_data = (double *) _var->data;
    for(int k=0;k<_var->data_size;++k)
    {
        for(int j=0;j<n_component;++j)
        {
            d_data[n_component*k + j] = 0.;
        }
    }
}


void set_e_debug(sale2d_var * _sale)
{
    /*
     * copy the vof from cells, set v_vof accord to e_vof
     * VARIABLES need to be synchronized in advance:
     *    e_vof
     * VARIABLES overwritten AFTER:
     *    v_vof
     */
    proc_info_2d* proc_ebug = (proc_info_2d*) _sale->e_debug->proc_self;
    int nx = proc_ebug->nx,ny=proc_ebug->ny,myrank=_sale->e_debug->rank;
    int nmat = _sale->nmat;

    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;++j)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _sale->foe->lid2nty(myrank,&lid,&lnt,_sale->v_vof);

            double * edebug;
            _sale->foe->get_double(&lid,&edebug,_sale->e_debug);
            double * epos;
            _sale->foe->get_double(&lid,&epos,_sale->e_pos);
            double * evol;
            _sale->foe->get_double(&lid,&evol,_sale->e_vol);
        }
    }
}

void update_e_vel(sale2d_var * _sale, node_type _nt)
{

    /*
     * calculate the velocity of element, store results in e_vel
     * (*e_vel will also used to store increment of vertex velocity)
     * the increment of vertex velocity has been moved into e_dvel.
     */

    field_var * _var = _sale->e_vel;
    field_opt * _opt = _sale->foe;

    proc_info_2d* proc_this = (proc_info_2d*) _sale->e_vel->proc_self;
    int nx = proc_this->nx,ny=proc_this->ny,myrank=_var->rank;
    int nmat = _sale->nmat;


    for(int j = 0; j < ny; ++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k, .y=j};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vel);
            if (!(lnt & _nt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};
            double *vvel_ptr[4],*vmass_ptr[4];
            for(int kv=0;kv<4;++kv)
            {
                _sale->fov->get_double(vlid[kv],&(vvel_ptr[kv]),_sale->v_vel);
                _sale->fov->get_double(vlid[kv],&(vmass_ptr[kv]),_sale->v_mass);
            }

            double *evel_ptr,*ewpt_ptr,*epos_ptr;
            _sale->foe->get_double(&lid,&evel_ptr,_sale->e_vel);
            _sale->foe->get_double(&lid,&ewpt_ptr,_sale->e_wpt);
            _sale->foe->get_double(&lid,&epos_ptr,_sale->e_pos);

            double *eden_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            double *esubvol_ptr;
            _sale->foe->get_double(&lid,&esubvol_ptr,_sale->e_subvol);
            vec_zero(evel_ptr,2);

            // skip the cell is full of vacuum
            if(eden_ptr[0] < TOLRHO) continue;

            /*
             * first check weather this cell is an isolated cell
             * all of its neighbours is void cell
             *
             * in numerical experiment, surprised velocity is removed
             *
             * the code/algorithm of detect isolated cell is moved into update_v_vel_inc
             * in the update of cells variables should not change variables on vertexes!!
             * the results will change with loop sequence ~
             */


            int vacuum_vertex = 0;
            for(int kv=0;kv<4;++kv)
            {
                if(vmass_ptr[kv][1] >= 4.0 || (vmass_ptr[kv][1] >= 3.0 && ewpt_ptr[kv] < 0.)) vacuum_vertex++;
            }
            if(4 == vacuum_vertex) continue;

            int solid_vertex = 0;
            double solid_mass = 0.;
            for(int kv=0;kv<4;++kv)
            {
                // if vertex is far from the material/vacuum interface, skip it in calculate the cell average velocity
//                if(ewpt_ptr[kv] > 0.5*Max(epos_ptr[2],epos_ptr[3]) || eden_ptr[0] < TOLRHO) continue;
                if(ewpt_ptr[kv] > 0.0 || eden_ptr[0] < TOLRHO) continue;
                solid_vertex++;
                scaler_add(evel_ptr,2,vvel_ptr[kv],vmass_ptr[kv][0]);
                solid_mass += vmass_ptr[kv][0];
            }


            if(solid_mass > 0.)
                vec_scaler(evel_ptr,2,1.0/solid_mass);
            else
                vec_zero(evel_ptr,2);

        }
    }

}


void update_v_vel_inc(sale2d_var * _sale, node_type _nt)
{
    /*
     * e_vel is the increment of velocity in flux
     * update the velocity of vertex according its neighbour elements
     * VARIABLES need to be synchronized in advance:
     *      e_vel(last cycle), e_den, e_dvel
     * VARIABLES will be overwritten AFTER:
     *      v_vel, v_mass
     */
    proc_info_2d* proc_v_vel = (proc_info_2d*) _sale->v_vel->proc_self;
    int nx = proc_v_vel->nx,ny=proc_v_vel->ny,myrank=_sale->v_vel->rank;
    unsigned int n_component = _sale->v_vel->bytes_per_unit / sizeof(double);


    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _sale->fov->get_nty(&lid,&lnt,_sale->v_vel);

            if(!(_nt&lnt)) continue;

            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};
            int v_offset[4] = {2, 3, 0, 1};

            double *vel_ptr;
            _sale->fov->get_double(&lid,&vel_ptr,_sale->v_vel);
            double * vmass_ptr;
            _sale->fov->get_double(&lid,&vmass_ptr,_sale->v_mass);
            double * vvof_ptr;
            _sale->fov->get_double(&lid,&vvof_ptr,_sale->v_vof);
            double *vvol_ptr;
            _sale->fov->get_double(&lid,&vvol_ptr,_sale->v_vol);


            double vel_inc[2],vel_old[2];
            vec_zero(vel_inc,2);
            vec_copy(vel_ptr,vel_old,2);
            double vvol = 0., vmass = 0.;

            double *eden_ptr[4], *evol_ptr[4], *evel_ptr[4], *ewpt_ptr[4], *edvel_ptr[4];
            for(int i=0;i<4;i++)
            {
                _sale->foe->get_double(elid[i],&(evel_ptr[i]),_sale->e_vel);
                _sale->foe->get_double(elid[i],&(evol_ptr[i]),_sale->e_vol);
                _sale->foe->get_double(elid[i],&(eden_ptr[i]),_sale->e_den);
                _sale->foe->get_double(elid[i],&(ewpt_ptr[i]),_sale->e_wpt);
                _sale->foe->get_double(elid[i],&(edvel_ptr[i]),_sale->e_dvel);
            }


            int void_cell = 0;
            int void_corner = 0;
            for(int i=0;i<4;++i)
            {
                if(ewpt_ptr[i][v_offset[i]] > 0.) void_corner++;
                if(eden_ptr[i][0] < TOLRHO)
                {
                    void_cell ++;
                    continue;
                }

#ifndef MEANVERTEXVOL
                double local_vertex_vol = vvol_ptr[i+1];
#else
                double local_vertex_vol = 0.25*vvol_ptr[0];
#endif
                vvol  += local_vertex_vol;
                vmass += eden_ptr[i][0] * local_vertex_vol;
                scaler_add(vel_inc, 2, edvel_ptr[i],  eden_ptr[i][0] * local_vertex_vol);
            }


            // the increment of momentum is vec_inc
            // the original momentum in this vertex is vel_old * vmass_ptr[0]
            if(vmass > TOLRHO * vvol)
            {
                liner_op(vel_ptr,2,vel_inc,vel_old,1.0/vmass,vmass_ptr[0]/vmass);
                vmass_ptr[0] = vmass;
                vmass_ptr[1] = vvol;
            }
            else
            {
                vmass_ptr[0] = 0.;
                vmass_ptr[1] = 0.;
                vec_zero(vel_ptr,2);
            }

            // check void element
            // change the vamss_ptr[1] to number of void corners
            vmass_ptr[1] = 1.0*void_corner;


            // velocity cutoff
            double vel_norm = vec_len(vel_ptr,2);
            if(vel_norm > _sale->velocity_max)
            {
                vec_scaler(vel_ptr,2, _sale->velocity_max / vel_norm);
            }

#ifdef VELMIN_CUTOFF
            // delete the tiny velocity, keep the div(vel) stable
            // open this option with cautions ??
            double *vacc_ptr;
            _sale->fov->get_double(&lid,&vacc_ptr,_sale->v_acc);
            //if(vec_len(vacc_ptr,2) < 10.*_sale->acceleration_min &&  vel_norm < _sale->velocity_min)
            if(vel_norm < _sale->velocity_min * VELMIN_CUTOFF)
            {
                vec_zero(vel_ptr,2);
            }
#endif

        }
    }
}


void check_v_vel(sale2d_var * _sale, node_type _nt)
{
    /*
     * process some vertexes that does not have any materials
     * those vertexes will fail and arise surprise velocity in update_v_vel_inc
     */

    proc_info_2d *proc_v_vel = (proc_info_2d *) _sale->v_vel->proc_self;
    int nx = proc_v_vel->nx, ny = proc_v_vel->ny, myrank = _sale->v_vel->rank;
    unsigned int n_component = _sale->v_vel->bytes_per_unit / sizeof(double);


    for(int j = 0; j < ny; ++j)
    {
        for(int k = 0; k < nx; ++k)
        {
            fv_id_2d lid = {.x=k, .y=j};
            node_type lnt = useless_section;
            _sale->fov->get_nty(&lid, &lnt, _sale->v_vel);

            if(!(_nt & lnt)) continue;

            fv_id_2d lb = {k - 1, j - 1}, rb = {k, j - 1}, ru = {k, j}, lu = {k - 1, j};
            fv_id_2d *elid[4] = {&lb, &rb, &ru, &lu};
            int v_offset[4] = {2, 3, 0, 1};

            double *vel_ptr;
            _sale->fov->get_double(&lid, &vel_ptr, _sale->v_vel);
            double *vmass_ptr;
            _sale->fov->get_double(&lid, &vmass_ptr, _sale->v_mass);
            double *vvof_ptr;
            _sale->fov->get_double(&lid, &vvof_ptr, _sale->v_vof);
            double *vvol_ptr;
            _sale->fov->get_double(&lid, &vvol_ptr, _sale->v_vol);


            double *eden_ptr[4], *evol_ptr[4], *evel_ptr[4], *ewpt_ptr[4], *edvel_ptr[4];
            for(int i = 0; i < 4; i++)
            {
                _sale->foe->get_double(elid[i], &(evel_ptr[i]), _sale->e_vel);
                _sale->foe->get_double(elid[i], &(evol_ptr[i]), _sale->e_vol);
                _sale->foe->get_double(elid[i], &(eden_ptr[i]), _sale->e_den);
                _sale->foe->get_double(elid[i], &(ewpt_ptr[i]), _sale->e_wpt);
                _sale->foe->get_double(elid[i], &(edvel_ptr[i]), _sale->e_dvel);
            }


            int void_cell = 0;
            int void_corner = 0;
            for(int i = 0; i < 4; ++i)
            {
                if(eden_ptr[i][0] < TOLRHO) void_cell++;
                if(ewpt_ptr[i][v_offset[i]] > 0.) void_corner++;

            }

            if(void_cell >= 4) continue;

            if(void_corner >= 4)
            {
                // this vertex does not have any materials
                // update momentum fails here!

#ifdef VOID_VERTEX_AVG_VEL
                // interpolate the void vertex through averge velocity
                int nearest_cell = 0;
                double avg_vel[2] = {0.};
                for(int kv = 0; kv < 4; ++kv)
                {
                    if(eden_ptr[kv][0] < TOLRHO) continue;
                    scaler_add(avg_vel,2,evel_ptr[kv],1.0);
                    nearest_cell++;
                }
                scaler_move(vel_ptr,2,avg_vel,1.0/nearest_cell);
#else
                // get average velocity weighted by density
                double avg_vel[2] = {0.};
                double avg_den = 0.;
                for(int kv = 0; kv < 4; ++kv)
                {
                    if(eden_ptr[kv][0] < TOLRHO) continue;
                    avg_den += eden_ptr[kv][0];
                    scaler_add(avg_vel,2,evel_ptr[kv],eden_ptr[kv][0]);
                }
                scaler_move(vel_ptr,2,avg_vel,1.0/avg_den);
#endif
            }

        }
    }
}


void update_e_gradvel(sale2d_var * _sale, node_type _nt)
{
    /*
     * update the gradient of velocity in element
     * VARIABLES required to be synchronized in advance
     *      v_vel, e_grad
     * VARIABLES overwritten:
     *      e_grad_vel
     */

    field_opt * _opt = _sale->foe;
    field_var * _var = _sale->e_grad_vel;

    proc_info_2d * proc_gradvel = (proc_info_2d*) _var->proc_self;
    int nx = proc_gradvel->nx, ny=proc_gradvel->ny, myrank=_sale->e_grad_vel->rank;

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _opt->get_nty(&lid,&lnt,_sale->e_grad_vel);
            if(!(_nt&lnt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};

            double * vvel_ptr[4];
            double vvelx[4] = {0.},vvely[4]={0.};

            for(int i=0;i<4;++i)
            {
                _sale->fov->get_double(vlid[i],&(vvel_ptr[i]),_sale->v_vel);
                vvelx[i] = vvel_ptr[i][0];
                vvely[i] = vvel_ptr[i][1];
            }

            double * egrad_ptr;
            _sale->foe->get_double(&lid,&egrad_ptr,_sale->e_grad);

            double * egradvel_ptr;
            _sale->foe->get_double(&lid,&egradvel_ptr,_sale->e_grad_vel);

            egradvel_ptr[0] = vec_dot(egrad_ptr+0,vvelx,4); /*partial v_r / partial r */
            egradvel_ptr[2] = vec_dot(egrad_ptr+0,vvely,4); /*partial v_z / partial r */
            egradvel_ptr[1] = vec_dot(egrad_ptr+4,vvelx,4); /*partial v_r / partial z */
            egradvel_ptr[3] = vec_dot(egrad_ptr+4,vvely,4); /*partial v_z / partial z */
            egradvel_ptr[4] = vec_dot(egrad_ptr+8,vvelx,4); /*partial v_r / r */

            double * edivvel_ptr = NULL;
            _sale->foe->get_double(&lid,&edivvel_ptr,_sale->e_div_vel);
            edivvel_ptr[0] = egradvel_ptr[0] + egradvel_ptr[3] + egradvel_ptr[4];

            double * ecsd_ptr, *ecenter_ptr, *evol_ptr, *evof_ptr, *eden_ptr;
            _sale->foe->get_double(&lid,&ecsd_ptr,_sale->e_csd);
            _sale->foe->get_double(&lid,&ecenter_ptr,_sale->e_pos);
            _sale->foe->get_double(&lid,&evol_ptr,_sale->e_vol);
            _sale->foe->get_double(&lid,&evof_ptr,_sale->e_vof);
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);

            double * eavisp_ptr;
            _sale->foe->get_double(&lid,&eavisp_ptr,_sale->e_q);

            double * epre_ptr;
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);

            fv_id_2d gid = {0};
            _sale->foe->lid2gid(myrank,&lid,&gid,_sale->e_grad_vel);

            //generate artificial pressure
            //if((_sale->partpres == 0 && evof_ptr[VACUUM] > TOLVOF) || (_sale->partpres == 1 && evof_ptr[VACUUM] > TOLVOF))
            // the second condition is unstable:ignore the "partpres" parameters
            // the extension cells ignore the artifical pressure
            if(evof_ptr[VACUUM] > TOLVOF || *epre_ptr <= 0.0)
            {
                eavisp_ptr[0] = 0.;
            } else
            {
                double earea = (*evol_ptr)/ecenter_ptr[0];
                earea = fabs(earea);
                double q_linear = _sale->avis1*(*ecsd_ptr)* sqrt(earea);
                double q_quadratic = _sale->avis2 * (*edivvel_ptr)*earea;
                *eavisp_ptr = MIN(0.,*edivvel_ptr) * (*eden_ptr) * (q_quadratic - q_linear);
#ifdef TOLQPRE
                if(*eavisp_ptr < TOLQPRE) *eavisp_ptr = 0.;
#endif
            }

        }
    }
}


void calculate_v_dis(sale2d_var * _sale, node_type _nt, double dt)
{
    /*
     * update the displacement of vertex in every cycle
     * VARIABLES need to be synchronized in advance:
     *      v_vel, v_acc
     * VARIABLES will be overwritten after:
     *      v_dis
     */
    proc_info_2d* proc_v_dis = (proc_info_2d*) _sale->v_dis->proc_self;
    int nx = proc_v_dis->nx,ny=proc_v_dis->ny,myrank=_sale->v_dis->rank;
    size_t n_component = _sale->v_dis->bytes_per_unit/ sizeof(double);

    for(int j=0;j<ny;++j)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _sale->fov->get_nty(&lid,&lnt,_sale->v_dis);
            if(!(_nt&lnt)) continue;

            double *vdis_ptr, *vvel_ptr, *vacc_ptr;
            _sale->fov->get_double(&lid,&vdis_ptr,_sale->v_dis);
            _sale->fov->get_double(&lid,&vvel_ptr,_sale->v_vel);
            _sale->fov->get_double(&lid,&vacc_ptr,_sale->v_acc);
            liner_op(vdis_ptr,n_component,vvel_ptr,vacc_ptr,dt,0.5*dt*dt);
            scaler_add(vvel_ptr,n_component,vacc_ptr,dt);
        }
    }
}

void calculate_flux_average(sale2d_var * _sale, node_type _nt)
{
    /*
     * the flux data in bx & by must been clean before this function !
     * VARIABLES need to be synchronized in advance
     *    e_vel
     *    v_vof, v_dis
     *    m_den, m_eng
     * VARIABLES will be overwritten:
     *    bx_phi, by_phi
     */
    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

   /* int debugvalue = 1;
    while (myrank==1 && debugvalue)
    {
        sleep(2);
    }*/

    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;j++)
        {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->lid2nty(myrank, &lid, &lnt, _sale->e_vof);
            if(!(_nt&lnt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};
            char   *vvof_cptr[4], *vpos_cptr[4], *vdis_cptr[4];
            double *vvof_dptr[4], *vpos_dptr[4], *vdis_dptr[4];
            for(int i = 0; i < 4; ++i)
            {
                _sale->fov->get_data2(vlid[i], &vvof_cptr[i], _sale->v_vof);
                _sale->fov->get_data2(vlid[i], &vpos_cptr[i], _sale->v_pos);
                _sale->fov->get_data2(vlid[i], &vdis_cptr[i], _sale->v_dis);
                vvof_dptr[i] = (double *) vvof_cptr[i];
                vpos_dptr[i] = (double *) vpos_cptr[i];
                vdis_dptr[i] = (double *) vdis_cptr[i];
            }

            if(vec_dis(vpos_dptr[0],(double []){12000.0,12000.0},2) < 3000.0)
            {
                double l_vol = get_flux_volcyl(vpos_dptr[0],vpos_dptr[3],vdis_dptr[0],vdis_dptr[3]);
                double r_vol = get_flux_volcyl(vpos_dptr[2],vpos_dptr[1],vdis_dptr[2],vdis_dptr[1]);
                double b_vol = get_flux_volcyl(vpos_dptr[1],vpos_dptr[0],vdis_dptr[1],vdis_dptr[0]);
                double u_vol = get_flux_volcyl(vpos_dptr[3],vpos_dptr[2],vdis_dptr[3],vdis_dptr[2]);
            }
            double l_vol = get_flux_volcyl(vpos_dptr[0],vpos_dptr[3],vdis_dptr[0],vdis_dptr[3]);
            double r_vol = get_flux_volcyl(vpos_dptr[2],vpos_dptr[1],vdis_dptr[2],vdis_dptr[1]);
            double b_vol = get_flux_volcyl(vpos_dptr[1],vpos_dptr[0],vdis_dptr[1],vdis_dptr[0]);
            double u_vol = get_flux_volcyl(vpos_dptr[3],vpos_dptr[2],vdis_dptr[3],vdis_dptr[2]);

            if(l_vol <= 0 && r_vol<=0 && b_vol<=0 && u_vol<=0)
            {
                continue;
            }

            double *egrad_ptr = NULL;
            _sale->foe->get_double(&lid,&egrad_ptr,_sale->e_grad);


            char *epos_cptr;
            double *escale_dptr;
            _sale->foe->get_data2(&lid, &epos_cptr, _sale->e_pos);
            escale_dptr = (double *) (epos_cptr) + 2;

            char *evof_cptr;
            _sale->foe->get_data2(&lid, &evof_cptr, _sale->e_vof);

            double flux_vol[10][4] = {0};

            double *ewpt_ptr;
            _sale->foe->get_double(&lid,&ewpt_ptr,_sale->e_wpt);

            double *evof_ptr;
            _sale->foe->get_double(&lid,&evof_ptr,_sale->e_vof);

            for(int cutid=0;cutid<nmat;++cutid)
            {

                if(evof_ptr[cutid] < TOLVOF)
                {
                    if(l_vol > 0) flux_vol[cutid][0] = 0.0;
                    if(r_vol > 0) flux_vol[cutid][2] = 0.0;
                    if(b_vol > 0) flux_vol[cutid][1] = 0.0;
                    if(u_vol > 0) flux_vol[cutid][3] = 0.0;

                    if(0 == cutid)
                    {
                        for(int kv=0;kv<4;++kv)
                        {
                            ewpt_ptr[kv] = -0.25;
                        }
                    }
                    continue;
                } else if(evof_ptr[cutid] > 1.0 - TOLVOF)
                {
                    if(l_vol > 0) flux_vol[cutid][0] = l_vol*evof_ptr[cutid];
                    if(r_vol > 0) flux_vol[cutid][2] = r_vol*evof_ptr[cutid];
                    if(b_vol > 0) flux_vol[cutid][1] = b_vol*evof_ptr[cutid];
                    if(u_vol > 0) flux_vol[cutid][3] = u_vol*evof_ptr[cutid];
                    if(0 == cutid)
                    {
                        for(int kv=0;kv<4;++kv)
                        {
                            ewpt_ptr[kv] = 0.25;
                        }
                    }
                    continue;
                }


                double local_vof_grad[2];
                double local_vvof[4] = {vvof_dptr[0][cutid], vvof_dptr[1][cutid],
                                        vvof_dptr[2][cutid], vvof_dptr[3][cutid]};
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

                double dstar = cut_length(*((double *) (evof_cptr) + cutid), (double[]) {cut_a, cut_b});
                double local_wpt[4], sum_wpt = 0., max_wpt;
                for (int i = 0; i < 4; ++i) {
                    local_wpt[i] = vec_dot(local_vof_grad, (double *) (vpos_cptr[i]), 2);
                    sum_wpt += local_wpt[i];
                }

                max_wpt = MAX(MAX(local_wpt[0], local_wpt[1]), MAX(local_wpt[2], local_wpt[3]));
                for (int i = 0; i < 4; ++i) {
                    local_wpt[i] = local_wpt[i] - max_wpt + dstar;
                }

                if(l_vol > 0) flux_vol[cutid][0] = l_vol * geometric_flux_ratio(local_wpt[3], local_wpt[0]);
                if(r_vol > 0) flux_vol[cutid][2] = r_vol * geometric_flux_ratio(local_wpt[2], local_wpt[1]);
                if(b_vol > 0) flux_vol[cutid][1] = b_vol * geometric_flux_ratio(local_wpt[0], local_wpt[1]);
                if(u_vol > 0) flux_vol[cutid][3] = u_vol * geometric_flux_ratio(local_wpt[3], local_wpt[2]);

                /*if(l_vol > 0) flux_vol[cutid][0] = l_vol*evof_ptr[cutid];
                if(r_vol > 0) flux_vol[cutid][2] = r_vol*evof_ptr[cutid];
                if(b_vol > 0) flux_vol[cutid][1] = b_vol*evof_ptr[cutid];
                if(u_vol > 0) flux_vol[cutid][3] = u_vol*evof_ptr[cutid];*/


                for(int kb=0;kb<4;++kb)
                {
                    flux_vol[cutid][kb] = MIN(flux_vol[cutid][kb],evof_ptr[cutid]);
                }

                if(cutid == 0)
                {
                    vec_copy(local_wpt,ewpt_ptr,4);
                }

            }

            fv_id_2d bxl = {k,j}, bxr={k+1,j}, byb={k,j}, byu={k,j+1};
            char *bphi_cptr[4];
            _sale->fobx->get_data2(&bxl,&bphi_cptr[0],_sale->bx_phi);
            _sale->fobx->get_data2(&bxr,&bphi_cptr[2],_sale->bx_phi);
            _sale->foby->get_data2(&byb,&bphi_cptr[1],_sale->by_phi);
            _sale->foby->get_data2(&byu,&bphi_cptr[3],_sale->by_phi);
            double *bphi_dptr[4];
            double direction[4] = {1.0,1.0,-1.0,-1.0};
            double lbru_vol[4] = {l_vol,b_vol,r_vol,u_vol};

            char * eden_cptr;
            _sale->foe->get_data2(&lid,&eden_cptr,_sale->m_den);

            char * eeng_cptr;
            _sale->foe->get_data2(&lid,&eeng_cptr,_sale->m_eng);

            char * evel_cptr;
            _sale->foe->get_data2(&lid,&evel_cptr,_sale->e_vel);
            for(int i=0;i<4;i++)
            {
                if(lbru_vol[i] <=0) continue; //only process the outflow/upwind case, negative volume indicate inflow
                bphi_dptr[i] = (double *) bphi_cptr[i];
                for(int cutid=0;cutid<nmat;++cutid)
                {
                    // set vol
                    bphi_dptr[i][cutid + 0*nmat] = direction[i]*flux_vol[cutid][i];
                    // set density
                    bphi_dptr[i][cutid + 1*nmat] = ((double *)(eden_cptr))[cutid];
                    // set energy
                    bphi_dptr[i][cutid + 2*nmat] = ((double *)(eeng_cptr))[cutid];
                }
                vec_copy((double *)evel_cptr,bphi_dptr[i]+3*nmat,2);
            }
        }
    }
}


void calculate_flux_upwind(sale2d_var * _sale, node_type _nt)
{
    /*
     * use traditional upwind stencil
     * calculate flux between elements
     */
    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;


    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;++j)
        {
            fv_id_2d elid = {k,j};
            node_type lnt = useless_section;
            _sale->foe->lid2nty(myrank,&elid,&lnt,_sale->e_vof);
            if(!(lnt&_nt)) continue;

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};
            double *vpos_dptr[4], *vdis_dptr[4];
            for(int i=0;i<4;++i)
            {
                _sale->fov->get_double(vlid[i],&vpos_dptr[i],_sale->v_pos);
                _sale->fov->get_double(vlid[i],&vdis_dptr[i],_sale->v_dis);
            }

            /*if(vec_dis(vpos_dptr[0],(double []){12000.0,12000.0},2) < 3000.0)
            {
                double l_vol = get_flux_volcyl(vpos_dptr[0],vpos_dptr[3],vdis_dptr[0],vdis_dptr[3]);
                double r_vol = get_flux_volcyl(vpos_dptr[2],vpos_dptr[1],vdis_dptr[2],vdis_dptr[1]);
                double b_vol = get_flux_volcyl(vpos_dptr[1],vpos_dptr[0],vdis_dptr[1],vdis_dptr[0]);
                double u_vol = get_flux_volcyl(vpos_dptr[3],vpos_dptr[2],vdis_dptr[3],vdis_dptr[2]);
            }*/


            double l_vol = get_flux_volcyl(vpos_dptr[0],vpos_dptr[3],vdis_dptr[0],vdis_dptr[3]);
            double r_vol = get_flux_volcyl(vpos_dptr[2],vpos_dptr[1],vdis_dptr[2],vdis_dptr[1]);
            double b_vol = get_flux_volcyl(vpos_dptr[1],vpos_dptr[0],vdis_dptr[1],vdis_dptr[0]);
            double u_vol = get_flux_volcyl(vpos_dptr[3],vpos_dptr[2],vdis_dptr[3],vdis_dptr[2]);

            double * evof_ptr = NULL;
            _sale->foe->get_double(&elid,&evof_ptr,_sale->e_vof);

            double flux_vol[MAXMAT][4] = {0.};
            for(int cutid=0;cutid<nmat;++cutid)
            {
                flux_vol[cutid][0] =      MIN(l_vol,0)*evof_ptr[cutid];
                flux_vol[cutid][2] = -1.0*MIN(r_vol,0)*evof_ptr[cutid];
                flux_vol[cutid][1] =      MIN(b_vol,0)*evof_ptr[cutid];
                flux_vol[cutid][3] = -1.0*MIN(u_vol,0)*evof_ptr[cutid];
            }

            fv_id_2d bxl = {k,j}, bxr={k+1,j}, byb={k,j}, byu={k,j+1};
            double * bphi_ptr[4];
            _sale->fobx->get_double(&bxl,&bphi_ptr[0],_sale->bx_phi);
            _sale->fobx->get_double(&bxr,&bphi_ptr[2],_sale->bx_phi);
            _sale->foby->get_double(&byb,&bphi_ptr[1],_sale->by_phi);
            _sale->foby->get_double(&byu,&bphi_ptr[3],_sale->by_phi);

            double * mden = NULL, *meng = NULL, *evel = NULL;
            _sale->foe->get_double(&elid,&mden,_sale->m_den);
            _sale->foe->get_double(&elid,&meng,_sale->m_eng);
            _sale->foe->get_double(&elid,&evel,_sale->e_vel);

            double lbru_vol[4] = {l_vol,b_vol,r_vol,u_vol};
            for(int i=0;i<4;++i)
            {
                if(lbru_vol[i] >=0) continue;
                for(int cutid=0;cutid<nmat;++cutid)
                {
                    bphi_ptr[i][cutid + 0*nmat] = flux_vol[cutid][i];
                    bphi_ptr[i][cutid + 1*nmat] = mden[cutid];
                    bphi_ptr[i][cutid + 2*nmat] = meng[cutid];
                }
                vec_copy(evel,bphi_ptr[i]+3*nmat,2);
            }
        }
    }

}


double record_avg_value(sale2d_var * _sale)
{
    /*
     * record the initial value used for damping
     * the average is the stable solution for equation at far from interested section.
     * at this time all partial_t term should be zero
     * at the those section, the number of materials should be ONE!
     * variables store in dampref of inital state is
     * dampref[0] = density
     * dampref[1] = energy*density (if allocated)
     * dampref[2] = pressure       (if allocated)
     * dampref[3] = local sound speed (if allocated)
     */

    proc_info_2d * proc_self = _sale->e_den->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int nvar_in_dampref = _sale->e_dampref->bytes_per_unit/ sizeof(double);
    double local_damp_time = 0.;
    for(int j = 0; j < ny; j++)
    {
        for (int k = 0; k < nx; ++k)
        {
            fv_id_2d lid = {.x=k, .y=j};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if(lnt & useless_section) continue;
            double * eden_ptr, *edampref_ptr;
            _sale->foe->get_double(&lid,&eden_ptr,_sale->e_den);
            _sale->foe->get_double(&lid,&edampref_ptr,_sale->e_dampref);
            edampref_ptr[0] = eden_ptr[0];

            // estimate the time when impact wave arrive here (the impact point is 0,0)
            double * epos_ptr, * ecsd_ptr, *esigma_ptr;
            _sale->foe->get_double(&lid,&epos_ptr,_sale->e_pos);
            _sale->foe->get_double(&lid,&ecsd_ptr,_sale->e_csd);
            _sale->foe->get_double(&lid,&esigma_ptr,_sale->e_sigma);
            double local_distance = vec_len(epos_ptr,2);

            if((esigma_ptr[X] > 0. || esigma_ptr[Y] > 0.) && ecsd_ptr[0] > 2.0)
            {
                local_damp_time = MAX(local_damp_time,local_distance/ecsd_ptr[0]);
            }


            if(nvar_in_dampref == 1) continue;
            double * mvof_ptr, *meng_ptr;
            _sale->foe->get_double(&lid,&mvof_ptr,_sale->m_vof);
            _sale->foe->get_double(&lid,&meng_ptr,_sale->m_eng);
            edampref_ptr[1] = 0.0;
            for(int matid=1;matid < _sale->nmat;++matid)
            {
                edampref_ptr[1] += mvof_ptr[matid] * meng_ptr[matid];
            }
            if(nvar_in_dampref == 2) continue;

            double * epre_ptr;
            _sale->foe->get_double(&lid,&epre_ptr,_sale->e_pre);
            edampref_ptr[2] = epre_ptr[0];

            if(nvar_in_dampref == 3) continue;
            edampref_ptr[3] = ecsd_ptr[0];
        }
    }

    return local_damp_time;
}

double minimum_pressure(double dam, double minpint, double minpdam)
{
    return (1.0 - dam)*minpint + dam*minpdam;
}


void update_state(sale2d_var * _sale, node_type _nt)
{
    /*
     * after the advect_flux function, the (density,energy) of
     * every materials in a element is specified.
     * EOS calculate the (sound of speed, pressure) through (density,energy)
     * (sound of speed, pressure) <= EOS(density,energy)
     * rewrite using get_double()
     * insert functionality of update_e_den()
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    static int nan_check_counter = 0;

    for(int j = 0; j < ny; j++)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if (!(lnt&_nt)) continue;

            double * evof, *mden, *meng, *edam, *mpty;
            _sale->foe->get_double(&lid,&evof,_sale->e_vof);
            _sale->foe->get_double(&lid,&mden,_sale->m_den);
            _sale->foe->get_double(&lid,&meng,_sale->m_eng);
            _sale->foe->get_double(&lid,&edam,_sale->e_dam);
            _sale->foe->get_double(&lid,&mpty,_sale->m_pty);

            double mpre[MAXMAT] = {0.}, mcsd[MAXMAT] = {0.}, mtem[MAXMAT] = {0.};

            for(int matid=1;matid<nmat;++matid)
            {
                if(evof[matid] < TOLVOF)
                {
                    mpre[matid] = 0.;
                    mcsd[matid] = 1.;
                    mtem[matid] = 0.;
                }
                else
                {
                    double mpre_min = minimum_pressure(edam[0],_sale->etb[matid].ref->minPint,_sale->etb[matid].ref->minPdam);
                    // crop pressure of materials in calculate_state
                    // calculate_state(_sale->etb + matid,mden[matid],meng[matid],mtem+matid,mpre+matid,mcsd+matid);
                    double mcfactor = calculate_state_correct(_sale->etb + matid,mden[matid],meng[matid],mpre_min,mtem+matid,mpre+matid,mcsd+matid,mpty+matid);
                    // adjust variables related to volume of material
                    // vof, den, eng

                    if(mcfactor > 1.0)
                    {
                        mden[matid]  *= mcfactor;
                        meng[matid]  *= mcfactor;
                        evof[VACUUM] += evof[matid]*(1.0 - 1.0/mcfactor);
                        evof[matid]  *= 1.0/mcfactor;
                    }
                }
            }


            double * epre, *ecsd, *etem, *eden;
            _sale->foe->get_double(&lid,&epre,_sale->e_pre);
            _sale->foe->get_double(&lid,&ecsd,_sale->e_csd);
            _sale->foe->get_double(&lid,&etem,_sale->e_tem);
            _sale->foe->get_double(&lid,&eden,_sale->e_den);


            double * emu;
            _sale->foe->get_double(&lid,&emu,_sale->e_mu);
            *emu = 0.0;

            *epre = 0.0;
            *ecsd = 0.0;
            *etem = 0.0;
            *eden = 0.0;

            double sum_mat = 0.;
            for(int matid=1;matid<nmat;++matid)
            {
                *epre += evof[matid] * mpre[matid];
                *ecsd += evof[matid] / Max(1.0,mcsd[matid]);
                *etem += evof[matid] * mtem[matid];
                *eden += evof[matid] * mden[matid];
                *emu  += evof[matid] / Max(1.0,mden[matid]*mcsd[matid]*mcsd[matid]*_sale->etb[matid].ref->GaCoef);
            }

            if((*emu) > 0.) *emu = 1.0/(*emu);
            if((*ecsd) > 0.) *ecsd = 1.0/(*ecsd);


            if((_sale->partpres == 0 && evof[VACUUM] > TOLVOF) || (_sale->partpres == 1 && evof[VACUUM] > 1.0 - TOLVOF))
            {
                *epre = 0.;
                *ecsd = 1.0;
                *emu  = 0.;
            } else if(_sale->partpres == 1)
            {
                *ecsd = *ecsd*(1.0 - evof[VACUUM]);
                *etem = *etem/(1.0 - evof[VACUUM]);
                *emu  = *emu *(1.0 - evof[VACUUM]);
            }
        }
    }
}

void update_state_balance(sale2d_var * _sale, node_type _nt)
{
    /*
     * after the advect_flux function, the (density,energy) of
     * every materials in a element is specified.
     * EOS calculate the (sound of speed, pressure) through (density,energy)
     * (sound of speed, pressure) <= EOS(density,energy)
     * rewrite using get_double()
     * insert functionality of update_e_den()
     */

    proc_info_2d * proc_self = _sale->e_vof->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int myrank = _sale->e_vof->rank;
    int nmat = _sale->nmat;

    for(int j = 0; j < ny; j++)
    {
        for(int k=0;k<nx;++k)
        {
            fv_id_2d lid = {.x=k, .y=j, 0, 0};
            node_type lnt = useless_section;
            _sale->foe->get_nty(&lid, &lnt, _sale->e_vof);
            if (!(lnt&_nt)) continue;

            double * evof, *mden, *meng, *edam, *mpty;
            _sale->foe->get_double(&lid,&evof,_sale->e_vof);
            _sale->foe->get_double(&lid,&mden,_sale->m_den);
            _sale->foe->get_double(&lid,&meng,_sale->m_eng);
            _sale->foe->get_double(&lid,&edam,_sale->e_dam);
            _sale->foe->get_double(&lid,&mpty,_sale->m_pty);

            double mpre[MAXMAT] = {0.}, mcsd[MAXMAT] = {0.}, mtem[MAXMAT] = {0.};
            int cmat = 0; // the number of types materials in this cell

            for(int matid=1;matid<nmat;++matid)
            {
                if(evof[matid] > TOLVOF) cmat++;
            }


            if(cmat < 2 || evof[VACUUM] > TOLVOF)
            {
                for(int matid=1;matid<nmat;++matid)
                {
                    if(evof[matid] < TOLVOF)
                    {
                        mpre[matid] = 0.;
                        mcsd[matid] = 1.;
                        mtem[matid] = 0.;
                    }
                    else
                    {
                        double mpre_min = minimum_pressure(edam[0],_sale->etb[matid].ref->minPint,_sale->etb[matid].ref->minPdam);
                        // crop pressure of materials in calculate_state
                        // calculate_state(_sale->etb + matid,mden[matid],meng[matid],mtem+matid,mpre+matid,mcsd+matid);
                        double mcfactor = calculate_state_correct(_sale->etb + matid,mden[matid],meng[matid],mpre_min,mtem+matid,mpre+matid,mcsd+matid,mpty+matid);
                        // adjust variables related to volume of material
                        // vof, den, eng
                        if(mcfactor > 1.0)
                        {
                            mden[matid]  *= mcfactor;
                            meng[matid]  *= mcfactor;
                            evof[VACUUM] += evof[matid]*(1.0 - 1.0/mcfactor);
                            evof[matid]  *= 1.0/mcfactor;
                        }
                    }
                }
            }
            else
            {
                // only adjust the mixed cells, 2 iteration
                double adjust_vof[MAXMAT] = {0.}; // the vof used for adjust pressure
                for(int k_adjust=0;k_adjust < 2*cmat; ++k_adjust)
                {
                    double average_pre = 0.;
                    for(int matid=1;matid<nmat;++matid)
                    {
                        if(evof[matid] < TOLVOF)
                        {
                            mpre[matid] = 0.;
                            mcsd[matid] = 1.;
                            mtem[matid] = 0.;
                        }
                        else
                        {
                            double mpre_min = minimum_pressure(edam[0],_sale->etb[matid].ref->minPint,_sale->etb[matid].ref->minPdam);
                            double mcfactor = calculate_state_correct(_sale->etb + matid,mden[matid],meng[matid],mpre_min,mtem+matid,mpre+matid,mcsd+matid,mpty+matid);
                            if(mcfactor > 1.0) // notice: mfactor is used to adjust negative pressure
                            {
                                mden[matid]  *= mcfactor;
                                meng[matid]  *= mcfactor;
                                evof[VACUUM] += evof[matid]*(1.0 - 1.0/mcfactor);
                                evof[matid]  *= 1.0/mcfactor;
                            }
                        }
                        average_pre += evof[matid]*mpre[matid];
                    }
                    // after average_pre calculated, lets estimated the vof used for adjust


                }
            }




            double * epre, *ecsd, *etem, *eden;
            _sale->foe->get_double(&lid,&epre,_sale->e_pre);
            _sale->foe->get_double(&lid,&ecsd,_sale->e_csd);
            _sale->foe->get_double(&lid,&etem,_sale->e_tem);
            _sale->foe->get_double(&lid,&eden,_sale->e_den);


            double * emu;
            _sale->foe->get_double(&lid,&emu,_sale->e_mu);
            *emu = 0.0;

            *epre = 0.0;
            *ecsd = 0.0;
            *etem = 0.0;
            *eden = 0.0;

            double sum_mat = 0.;
            for(int matid=1;matid<nmat;++matid)
            {
                *epre += evof[matid] * mpre[matid];
                *ecsd += evof[matid] / Max(1.0,mcsd[matid]);
                *etem += evof[matid] * mtem[matid];
                *eden += evof[matid] * mden[matid];
                *emu  += evof[matid] / Max(1.0,mden[matid]*mcsd[matid]*mcsd[matid]*_sale->etb[matid].ref->GaCoef);
            }

            if((*emu) > 0.) *emu = 1.0/(*emu);
            if((*ecsd) > 0.) *ecsd = 1.0/(*ecsd);


            if((_sale->partpres == 0 && evof[VACUUM] > TOLVOF) || (_sale->partpres == 1 && evof[VACUUM] > 1.0 - TOLVOF))
            {
                *epre = 0.;
                *ecsd = 1.0;
                *emu  = 0.;
            } else if(_sale->partpres == 1)
            {
                *ecsd = *ecsd*(1.0 - evof[VACUUM]);
                *etem = *etem/(1.0 - evof[VACUUM]);
                *emu  = *emu *(1.0 - evof[VACUUM]);
            }
        }
    }
}


void init_v_field(sale2d_var * _sale)
{
    field_var * _var = _sale->v_vel;
    field_opt * _opt = _sale->fov;

    proc_info_2d * proc_self = (proc_info_2d *) _var->proc_self;
    int nx = proc_self->nx, ny=proc_self->ny;
    int nmat = _sale->nmat;

    for(int k=0;k<nx;++k)
    {
        for(int j=0;j<ny;++j)
        {
            fv_id_2d lid = {.x=k,.y=j};
            node_type lnt = useless_section;
            _opt->get_nty(&lid,&lnt,_var);
            if(useless_section & lnt) continue;

            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};
            int n_component = _var->bytes_per_unit/ sizeof(double);

            double *vel_ptr,*mass_ptr;
            _sale->fov->get_double(&lid,&vel_ptr,_sale->v_vel);
            vec_zero(vel_ptr,n_component);
            _sale->fov->get_double(&lid,&mass_ptr,_sale->v_mass);
            double *vsubvol_ptr,*vvol_ptr;
            _sale->fov->get_double(&lid,&vsubvol_ptr,_sale->v_subvol);
            _sale->fov->get_double(&lid,&vvol_ptr,_sale->v_vol);

            double *eden_ptr, *evol_ptr, *evel_ptr, emass =0., evol = 0.;
            double local_frac = 0.25;
            for(int i=0;i<4;i++)
            {
                _sale->foe->get_double(elid[i],&evel_ptr,_sale->e_vel);
                _sale->foe->get_double(elid[i],&evol_ptr,_sale->e_vol);
                _sale->foe->get_double(elid[i],&eden_ptr,_sale->e_den);

                double local_vertex_vol = vvol_ptr[i+1];
                // the old version is local_frac*(*evol_ptr)
                scaler_add(vel_ptr,n_component,evel_ptr,(*eden_ptr)*local_vertex_vol);
                emass += (*eden_ptr)*local_vertex_vol;
                evol  += local_vertex_vol;
            }

            if(emass > TOLRHO * evol)
            {
                vec_scaler(vel_ptr,n_component,1.0/emass);
                mass_ptr[0] = emass;
                mass_ptr[1] = evol;
            } else
            {
                vec_zero(vel_ptr,n_component);
                mass_ptr[0] = 0.;
                mass_ptr[1] = 0.;
            }
        }
    }

}


double calculate_dt(sale2d_var * _sale, double cur_dt)
{
    double rtn_dt = 1.1*cur_dt;
    proc_info_2d * e_proc = _sale->e_vof->proc_self;
    int e_nx = e_proc->nx, e_ny = e_proc->ny;

    for(int j=0;j<e_ny;++j)
    {
        for(int k=0;k<e_nx;++k)
        {
            fv_id_2d elid = {.x=k, .y=j};

            node_type nt;
            _sale->foe->get_nty(&elid,&nt,_sale->e_vel);

            if(nt&receive_section || nt&useless_section) continue;

            double * evof_ptr;
            _sale->foe->get_double(&elid,&evof_ptr,_sale->e_vof);
            if(evof_ptr[0] > 1.0 - TOLVOF) continue;

            double * evel, *ecsd, *epos;

            _sale->foe->get_double(&elid,&evel,_sale->e_vel);
            _sale->foe->get_double(&elid,&ecsd,_sale->e_csd);
            _sale->foe->get_double(&elid,&epos,_sale->e_pos);
            

            /*
             * when the cell is near r=0, the cfl condition is not same as cartesian
             * add cyl_ratio to reduce the dt of whole model and increase accuracy around r = 0
             * in most of the time this option is disabled
             */
            double cyl_ratio = 1.0 ;
#ifdef UNIFORM_XGRID_DT
            if(fabs(epos[0]) > 0.0)
            {
                cyl_ratio = fabs(epos[0])/(fabs(epos[0]) + 0.5*epos[2]);
            }
#endif

            // dt at r/z direction
            double dt_r = cyl_ratio*epos[2]/(fabs(evel[X]) + ecsd[0]);
            double dt_z = 1.0*epos[3]/(fabs(evel[Y]) + ecsd[0]);
            double dt_rz = DTF*MIN(dt_r,dt_z);
#ifdef CHECK_NAN
            if(isnan(dt_r) || isnan(dt_z))
            {
                fprintf(stderr,"Warning vel=(%f,%f), csd=(%f), pos=(%f,%f) generate Nan dt\n",evel[X],evel[Y],ecsd[0],epos[X],epos[Y]);
                continue;
            }
#endif
            rtn_dt = MIN(rtn_dt,dt_rz);
        }
    }

    return rtn_dt;
}

void init_sale2d_ghost(sale2d_var * _sale, mesh2d_info * _minfo)
{
    // allocate the ghost information
    init_field_ght_table(_sale->m_vof,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->m_vof,ght_var_always_copy);

    init_field_ght_table(_sale->e_vel,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_vel,ght_velocity);

    init_field_ght_table(_sale->e_dvel,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_dvel,ght_velocity_inc);

    init_field_ght_table(_sale->e_den,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_den,ght_density);

    init_field_ght_table(_sale->e_pre,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_pre,ght_var_always_copy);

    init_field_ght_table(_sale->e_ste,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_ste,ght_stress);

    init_field_ght_table(_sale->e_q,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_q,ght_var_always_copy);

    init_field_ght_table(_sale->m_den,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->m_den,ght_density);

    init_field_ght_table(_sale->m_eng,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->m_eng,ght_density);

    init_field_ght_table(_sale->m_pty,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->m_pty,ght_var_always_copy);

    init_field_ght_table(_sale->e_cvs,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_cvs,ght_var_always_copy);

    // allocate debug ghost value
    init_field_ght_table(_sale->e_debug,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_debug,ght_var_always_copy);

    // allocate debug vibration
    init_field_ght_table(_sale->e_vib,_sale->foe);
    set_ghost_pattern(_minfo->go_list,_sale->e_vib,ght_var_always_copy);

//    init_field_ght_table(_sale->v_vel,_sale->fov);
//    set_ghost_pattern(_minfo->go_list,_sale->v_vel,ght_velocity);

    /*
     * subroutine about restrict vertexes is also set here
     * the last args in set_restrict_pattern (var_type) is not used
     */

    init_field_rst_table(_sale->v_dis,_sale->fov);
    set_restrict_pattern(_minfo->go_list,_sale->v_dis,ght_normal);

    init_field_rst_table(_sale->v_acc,_sale->fov);
    set_restrict_pattern(_minfo->go_list,_sale->v_acc,ght_normal);

    init_field_rst_table(_sale->v_vel,_sale->fov);
    set_restrict_pattern(_minfo->go_list,_sale->v_vel,ght_normal);

}

int init_sale2d_tracer(sale2d_var * _sale, mesh2d_info * _mesh)
{
    _sale->tracer_list = fl_new();
    fl_mpi_init(_sale->tracer_list,_sale->e_pos,_sale->e_vof,_sale->foe,
                _sale->v_pos,_sale->v_dis,_sale->v_vel,_sale->v_bf,_sale->fov);
    fl_record_init(_sale->tracer_list,_sale->e_den,_sale->e_tem,_sale->e_pre,_sale->v_acc);

    // directly set the threshold for ejecta
    _sale->tracer_list->threshold = _sale->threshold_z;
    fl_tracer_init(_sale->tracer_list,_mesh->tracer_strip);
}

int export_init_tracers(sale2d_var * _sale, mesh2d_info * _mesh)
{
    int num_tracers = _sale->tracer_list->len - _sale->tracer_list->dirty_nodes;
    char init_txt_name[200];
    sprintf(init_txt_name,"./txt/mass.tracer.proc%04d.txt",_sale->e_pos->rank);
    FILE * fp = fopen(init_txt_name,"w");
    if(fp == NULL)
    {
        fprintf(stderr,"cannot create %s\n",init_txt_name);
        exit(0);
    }

    fprintf(fp,"%d\n",num_tracers);

    field_list_iterator  fit;
    fl_iterator_init(&fit,_sale->tracer_list,field_list_head);
    field_list_node * node = NULL;

    while(NULL != (node = fl_iterator_next(&fit)))
    {
        if(node->state != field_list_deleted)
        {
            // the data in recorder is updated, if need output
            double *node_evol_ptr, *node_eden_ptr;
            _sale->foe->get_double(&(node->id),&node_eden_ptr,_sale->e_den);
            _sale->foe->get_double(&(node->id),&node_evol_ptr,_sale->e_vol);

            fprintf(fp,"%d, %d, %d, %10.5e,%10.5e, %10.5e, %10.5e\n", node->data.tag, node->data.gx, node->data.gy,
                    node->data.position[X], node->data.position[Y], node_eden_ptr[0], node_evol_ptr[0]);
        }
    }

    fclose(fp);
}


