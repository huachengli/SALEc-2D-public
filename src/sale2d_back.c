//
// Created by huach on 30/11/2022.
//
// functions not used directly in the main
// unclassified functions, most of them are related to io
#include "sale2d.h"

double flux_ratio(const double * _p)
{
    double a = _p[0], b=_p[1];
    if(b > a)
    {
        a = _p[1]; b=_p[0];
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


void load_projectile_info(InputFile * ifp, projectile_info * _pinfo)
{
    int nmat = GetValueI(ifp,"material.nm","3");
    _pinfo->material = -1;
    char target_name[80];
    GetValueS(ifp,"projectile.material",target_name,"unknown");

    for(int matid=0;matid<nmat;++matid)
    {
        char tmp_name[80];
        GetValueSk(ifp,"material.name",tmp_name,matid,"vacuum_");
        if(0 == strcasecmp(target_name,tmp_name))
        {
            _pinfo->material = matid;
            break;
        }
    }

    if(_pinfo->material < 0)
    {
        fprintf(stderr,"can not get material id of [%s]\n",target_name);
        exit(0);
    }

    _pinfo->radiu = GetValueD(ifp,"projectile.radiu","1.0");
    _pinfo->center[0] = GetValueDk(ifp,"projectile.center",0,"0.0");
    _pinfo->center[1] = GetValueDk(ifp,"projectile.center",1,"300.0");
    _pinfo->velocity[0] = GetValueDk(ifp, "projectile.velocity", 0, "0.0");
    _pinfo->velocity[1] = GetValueDk(ifp, "projectile.velocity", 1, "-1000.0");
    _pinfo->damage = GetValueDk(ifp,"projectile.damage",1,"1.0");
    _pinfo->pressure = GetValueDk(ifp,"projectile.pressure",1,"1.0");
    _pinfo->temperature = GetValueDk(ifp,"projectile.temperature",1,"293.0");

}

void load_mesh_info(InputFile * ifp, mesh2d_info * _minfo, const char * eospath)
{

    // section processor
    _minfo->npgx = GetValueI(ifp,"processor.npgx","2");
    _minfo->npgy = GetValueI(ifp,"processor.npgy","2");

    // section mesh
    _minfo->dx = GetValueD(ifp,"mesh.dx","1.");
    _minfo->dy = GetValueD(ifp,"mesh.dy","1.");
    _minfo->oy = GetValueD(ifp,"mesh.oy","0.7");
    _minfo->ox = GetValueD(ifp,"mesh.ox","0.0");
    _minfo->nx = GetValueI(ifp,"mesh.npx","32");
    _minfo->ny = GetValueI(ifp,"mesh.npy","32");

    char GridExtOpt[MaxStrLen];
    GetValueS(ifp,"mesh.ext",GridExtOpt,"off");
    if(strcasecmp(GridExtOpt,"off")==0)
    {
        _minfo->x_ext = -1.0;
        _minfo->y_ext = -1.0;
        _minfo->nx_ext[0] = _minfo->nx_ext[1] = 0;
        _minfo->ny_ext[0] = _minfo->ny_ext[1] = 0;

    } else
    {
        _minfo->x_ext     = GetValueDk(ifp,"mesh.ex",0,"1.05");
        _minfo->nx_ext[0] = GetValueIk(ifp,"mesh.ex",1,"0");
        _minfo->nx_ext[1] = GetValueIk(ifp,"mesh.ex",2,"0");

        _minfo->y_ext     = GetValueDk(ifp,"mesh.ey",0,"1.05");
        _minfo->ny_ext[0] = GetValueIk(ifp,"mesh.ey",1,"0");
        _minfo->ny_ext[1] = GetValueIk(ifp,"mesh.ey",2,"0");
    }

    // set numerical parameters
    _minfo->avis2 = GetValueD(ifp,"numerical.AVISQ","1.2");
    _minfo->avis1 = GetValueD(ifp,"numerical.AVISL","0.5");
    _minfo->anc   = GetValueD(ifp,"numerical.ANC","0.1");
    _minfo->anc_depth = GetValueD(ifp, "numerical.ANCD", "2.");
    _minfo->anc_width = GetValueD(ifp, "numerical.ANCW", "3.");
    GetValueS(ifp,"numerical.ANCPAT",_minfo->anc_pattern, "pml");
    _minfo->ejecta_threshold = GetValueD(ifp, "ejecta.threshold", "2.");
    _minfo->ejecta_resistance = GetValueD(ifp, "ejecta.resistance", "-1.0");
    _minfo->ejecta_friction   = GetValueD(ifp, "ejecta.friction", "0.0");

    _minfo->density_cutoff = GetValueD(ifp,"numerical.MinDensity","1.0e-2");
    _minfo->velocity_max = GetValueD(ifp, "numerical.MaxVelocity", "-1.63");
    _minfo->velocity_min = GetValueD(ifp, "numerical.MinVelocity", "-1e-8");
    _minfo->acceleration_min = GetValueD(ifp, "numerical.MinAcceleration", "-1e-3");
    _minfo->partpres = GetValueI(ifp,"numerical.partpres","0");

    // damp parameter
    char DampOpt[MaxStrLen];
    GetValueS(ifp,"numerical.damping",DampOpt,"off");
    _minfo->damp_x = GetValueI(ifp,"numerical.DampX","0");
    _minfo->damp_y = GetValueI(ifp,"numerical.DampY","0");
    _minfo->lambda_x = GetValueD(ifp,"numerical.LambdaX","10");
    _minfo->lambda_y = GetValueD(ifp,"numerical.LambdaY","10");
    _minfo->damp_t   = GetValueD(ifp,"numerical.damp_time","-1.0");
    if(_minfo->damp_t < 0.){
        _minfo->damp_t = GetValueD(ifp,"cycle.time1","100.0");
    }
    _minfo->damp_x_left = GetValueI(ifp,"numerical.DampX_Left","0");
    _minfo->damp_y_top = GetValueI(ifp,"numerical.DampY_Right","0");
    if(strcasecmp(DampOpt,"on")==0 || strcasecmp(DampOpt,"pml")==0)
    {
        // perfect match layer
        _minfo->damp_type = pml;
    }
    else if(strcasecmp(DampOpt,"pml2")==0)
    {
        // perfect match layer option 2
        _minfo->damp_type = pml2;
    }
    else if(strcasecmp(DampOpt,"pml3")==0)
    {
        // perfect match layer option 2
        _minfo->damp_type = pml3;
    }
    else if(strcasecmp(DampOpt,"sbc")==0)
    {
        // sponge boundary condition
        _minfo->damp_type = sbc;

        _minfo->lambda_x = Wind(_minfo->lambda_x,0.0,0.5);
        _minfo->lambda_y = Wind(_minfo->lambda_y,0.0,0.5);
    } else if(strcasecmp(DampOpt,"ctlg") == 0 || strcasecmp(DampOpt,"rock")==0)
    {
        _minfo->damp_type = ctlg;
    }
    else
    {
        _minfo->damp_type = none;
        _minfo->damp_x = 0;
        _minfo->damp_y = 0;
        _minfo->lambda_x = 0.;
        _minfo->lambda_y = 0.;
        _minfo->damp_t = -1.0;
        _minfo->damp_x_left = 0;
        _minfo->damp_y_top = 0;
    }


    // some default values, not changed in inputfile
    _minfo->xpad = 2;
    _minfo->ypad = 2;

    _minfo->v_nx = _minfo->nx + 2*_minfo->xpad + 1;
    _minfo->v_ny = _minfo->ny + 2*_minfo->ypad + 1;
    _minfo->e_nx = _minfo->v_nx - 1;
    _minfo->e_ny = _minfo->v_ny - 1;

    // section of material
    _minfo->n_materials = GetValueI(ifp,"material.nm","2");
    char local_etype[80], local_ename[80];
    for(int matid=0;matid<_minfo->n_materials;++matid)
    {
        GetValueSk(ifp,"material.name",local_ename,matid,"vacuum_");
        GetValueSk(ifp,"material.postfix",local_etype,matid,"null");

        if(0==strcasecmp(local_etype,"aneos"))
        {
            _minfo->etype[matid] = ANEOS;
            sprintf(_minfo->ename[matid], "%s/%s.aneos", eospath,local_ename);
        } else if(0==strcasecmp(local_etype,"Tillotson"))
        {
            _minfo->etype[matid] = Tillotson;
            sprintf(_minfo->ename[matid], "%s/%s.tillotson", eospath,local_ename);
        } else if(0 == strcasecmp(local_etype,"aneos.sd"))
        {
            _minfo->etype[matid] = ANEOSSieDen;
            sprintf(_minfo->ename[matid],"%s/%s.aneos.sd", eospath,local_ename);
        } else if(0==strcasecmp(local_etype,"aneos.td"))
        {
            _minfo->etype[matid] = ANEOS;
            sprintf(_minfo->ename[matid], "%s/%s.aneos.td", eospath,local_ename);
        } else
        {
            _minfo->etype[matid] = Unknown;
            sprintf(_minfo->ename[matid], "%s/%s.Unknown", eospath,local_ename);
            if(0 == matid) continue;
            fprintf(stderr,"##[unknown eos type] encounter with %s!", _minfo->ename[matid]);
            exit(0);
        }

    }

    GetValueS(ifp,"output.format",_minfo->format,"binary");
    GetValueS(ifp,"output.prefix",_minfo->prefix,"ParaTest");

    _minfo->gravity[2] = 0;
    _minfo->gravity[1] = GetValueDk(ifp,"condition.gravity",1,"-1.62");
    _minfo->gravity[0] = GetValueDk(ifp,"condition.gravity",0,"0.0");

    // load boundary condition
    GetValueS(ifp,"boundary.left",_minfo->ghost_boundary[0],"freeslip");
    GetValueS(ifp,"boundary.right",_minfo->ghost_boundary[1],"outflow");
    GetValueS(ifp,"boundary.bottom",_minfo->ghost_boundary[2],"freeslip");
    GetValueS(ifp,"boundary.top",_minfo->ghost_boundary[3],"outflow");

    _minfo->local_boundary[0] = left_boundary;
    _minfo->local_boundary[1] = right_boundary;
    _minfo->local_boundary[2] = bottom_boundary;
    _minfo->local_boundary[3] = top_boundary;

    // tracer info
    _minfo->tracer_strip = GetValueI(ifp,"numerical.tarcer_strip","2");

    for(int k=0;k<4;++k)
    {
        _minfo->go_list[k] = boundary_str2opt(_minfo->ghost_boundary[k],_minfo->local_boundary[k]);
    }

    if(_minfo->ox < 1e-6) _minfo->go_list[0] = ght_mirror;


    // load state reference information
    _minfo->ref = (state_reference *) malloc(sizeof(state_reference)*_minfo->n_materials);
    load_state_reference(ifp,_minfo->ref,_minfo->n_materials);
    load_target_info(ifp,&(_minfo->tinfo));
    load_projectile_info(ifp,&(_minfo->pinfo));
    load_target_addition_info(ifp,&(_minfo->tainfo));
}


void load_state_reference(InputFile * ifp, state_reference * _sref, int nmat)
{
    if(NULL == _sref)
    {
        fprintf(stderr,"memory error in load_state_reference\n");
        exit(0);
    }

    for(int matid=1;matid<nmat;++matid)
    {
        state_reference * _ref = _sref + matid;

        _ref->MaterialId = matid;
        /*
         * VaporDen, MeltDen
         * is set by eos
         */

        double nu = GetValueDk(ifp,"material.poisson",matid,"0.25");
        _ref->GaCoef = 1.5* (1.0 - 2.0*nu)/(1.0 + nu);


        // (ROCK)
        _ref->Yint0    = GetValueDk(ifp,"material.yint0",matid,"1.0");
        _ref->Yfricint = GetValueDk(ifp,"material.yintfri",matid,"1.0");
        _ref->Ylimint  = GetValueDk(ifp,"material.yintlim",matid,"2.50e9");
        _ref->Ydam0    = GetValueDk(ifp,"material.ydam0",matid,"1.0");
        _ref->Yfricdam = GetValueDk(ifp,"material.ydamfri",matid,"1.0");
        _ref->Ylimdam  = GetValueDk(ifp,"material.ydamlim",matid,"2.50e9");

        // (strength function)
        GetValueSk(ifp,"material.yshear",_ref->yfunc,matid,"SimpleRock");

        // (tensile failure)
        _ref->Yten0 = GetValueDk(ifp,"material.yten0",matid,"3.0e8");
        GetValueSk(ifp,"material.ytens",_ref->tfunc,matid,"SimpleTensile");

        // (IVANOV)
        _ref->IvanA = GetValueDk(ifp,"material.IvanA",matid,"1.0e-4");
        _ref->IvanB = GetValueDk(ifp,"material.IvanB",matid,"1.0e-11");
        _ref->IvanC = GetValueDk(ifp,"material.IvanC",matid,"3.0e8");

        // (COLLINS)
        _ref->Pbd  = GetValueDk(ifp,"material.Pbd",matid,"2.7e9");
        _ref->Pbp  = GetValueDk(ifp,"material.Pbp",matid,"4.1e9");
        GetValueSk(ifp,"material.damage",_ref->dSDf,matid,"SimpleShear");


        // set Johson-Cook strength parameters
        _ref->JCA      = GetValueDk(ifp,"material.JcA",matid,"4.90e7");
        _ref->JCB      = GetValueDk(ifp,"material.JcB",matid,"1.57e8");
        _ref->JCN      = GetValueDk(ifp,"material.JcN",matid,"1.67e-1");
        _ref->JCC      = GetValueDk(ifp,"material.JcC",matid,"1.60e-2");
        _ref->JCM      = GetValueDk(ifp,"material.JcM",matid,"1.70");
        _ref->JCTREF   = GetValueDk(ifp,"material.JcTref",matid,"8.00e2");
        _ref->JCMinP   = GetValueDk(ifp,"material.JcminP",matid,"-2.44e9");
        _ref->JCMaxPlastic = GetValueDk(ifp,"material.JCMaxPlastic",matid,"10.0");


        if(strcasecmp(_ref->yfunc,"JohnsonCook2") == 0 || strcasecmp(_ref->yfunc,"JohnsonCook1") == 0)
        {
            // set minimum pressure for JNCK model
            _ref->minPdam = _ref->JCMinP;
            _ref->minPint = _ref->JCMinP;
        }
        else if(strcasecmp(_ref->yfunc,"VonMises1") == 0 || strcasecmp(_ref->yfunc,"VonMises2") == 0)
        {
            // set minimum pressure for Von Mises strength model
            _ref->minPdam = _ref->JCMinP;
            _ref->minPint = _ref->JCMinP;
        }
        else if(strstr(_ref->yfunc,"Rock") != NULL || strcasecmp(_ref->yfunc,"IceStrength") == 0)
        {
            // set the parameters related to pressure_min : ROCK
            _ref->minPint = _ref->Yint0/_ref->Yfricint * (_ref->Yint0/_ref->Ylimint - 1.0);
            _ref->minPdam = -1.0*_ref->Ydam0/_ref->Yfricdam;
        }
        else
        {
            _ref->minPdam = 0.;
            _ref->minPint = 0.;
        }


        // set thermal temperature parameters
        _ref->Asimon   = 1.0/GetValueDk(ifp,"material.SimonA",matid,"6.00e9");
        _ref->Csimon   = 1.0/GetValueDk(ifp,"material.SimonC",matid,"3.00");
        _ref->Tmelt0   = GetValueDk(ifp,"material.SimonT0",matid,"933.0");
        _ref->Tfrac    = GetValueDk(ifp,"material.OhnakaXi",matid,"1.2");
        _ref->Tdelta   = GetValueDk(ifp,"material.Tdelta",matid,"0.0");
        // liquids temperature parameters
        _ref->PolyLiq  = GetValueIk(ifp,"material.PolyLiq",matid,"0");
        _ref->Tlidc1   = GetValueDk(ifp,"material.Tlidc1",matid,"0.0");
        _ref->Tlidc2   = GetValueDk(ifp,"material.Tlidc2",matid,"0.0");
        _ref->Tlidc3   = GetValueDk(ifp,"material.Tlidc3",matid,"0.0");
        // solidus temperature parameters
        _ref->PolySol  = GetValueIk(ifp,"material.PolySol",matid,"0");
        _ref->Tsolc1   = GetValueDk(ifp,"material.Tsolc1",matid,"0.0");
        _ref->Tsolc2   = GetValueDk(ifp,"material.Tsolc2",matid,"0.0");
        _ref->Tsolc3   = GetValueDk(ifp,"material.Tsolc3",matid,"0.0");

        // set acfl parameters
        double BlockSize = GetValueDk(ifp,"material.BlockSize",matid,"-1.0");
        if(BlockSize < 0.)
        {
            BlockSize = GetValueD(ifp,"projectile.radiu","1.0");
        }
        _ref->GammaEta  = GetValueDk(ifp, "material.GammaEta", matid, "8.0e-3");
        _ref->GammaBeta = GetValueDk(ifp, "material.GammaBeta", matid, "1.15e2");
        _ref->Toff     = GetValueDk(ifp,"material.Toff",matid,"16.0");
        _ref->Cvib     = GetValueDk(ifp,"material.Cvib",matid,"0.1");
        _ref->VibMax   = GetValueDk(ifp,"material.VibMax",matid,"200.0");
        _ref->Pvlim    = GetValueDk(ifp,"material.Pvlim",matid,"2.50e10");
        _ref->Acvis    = _ref->GammaEta * BlockSize * 5000.0;
        _ref->Tdamp    = _ref->GammaBeta * BlockSize / 5000.0;

        // set for melt viscoity
        _ref->melt_viscosity = GetValueDk(ifp, "material.viscosityM", matid, "0.0");

        // get porosity
        GetValueSk(ifp,"material.Porosity",_ref->porosityfunc,matid,"None");
        if(strcasecmp(_ref->porosityfunc,"Wunnemann") == 0)
        {
            _ref->alpha0 = GetValueDk(ifp,"material.alpha0",matid,"1.0");
            _ref->alphax = GetValueDk(ifp,"material.alphax",matid,"1.0");
            _ref->epse0  = GetValueDk(ifp,"material.epse0",matid,"-1.0e-5");
            _ref->kappa  = GetValueDk(ifp,"material.kappa",matid,"0.98");
            _ref->chi    = GetValueDk(ifp,"material.chi",matid,"1.0");
            _ref->avel   = GetValueDk(ifp,"material.avel",matid,"1e-4");
            // derive other parameters
            _ref->alphamax = _ref->alpha0;
            _ref->epsex = _ref->epse0 + log(_ref->alphax/_ref->alpha0)/_ref->kappa;
            _ref->epsec = _ref->epsex + 2.0*(1.0-_ref->alphax)/(_ref->kappa*_ref->alphax);
            if(_ref->epsex < _ref->epsec) _ref->epsec = _ref->epsex;
            _ref->alphae = _ref->alpha0*(1.0+(1.0-_ref->chi*_ref->chi)*_ref->epse0);
            _ref->alphamax2 = _ref->alpha0;
        } else if(strcasecmp(_ref->porosityfunc,"None") == 0)
        {
            // no porosity set
            _ref->alpha0   = 1.0;
            _ref->epse0    = 0.;
            _ref->kappa    = 0.;
            _ref->chi      = 0.;
            _ref->alphax   = 1.0;
            _ref->alphamax = 1.0;
        }
        else
        {
            fprintf(stdout,"undefined porosity model:%s\n",_ref->porosityfunc);
            exit(0);
        }
    }
}


ght_opt boundary_str2opt(const char * boundary_str, boundary_type local_bt)
{
    // just support outflow, freeslip and fixed boundary
    if(strcasecmp(boundary_str,"outflow") == 0)
    {
        return ght_outflow;
    } else if(strcasecmp(boundary_str,"fixed") == 0)
    {
        return ght_fixed;
    } else if(strcasecmp(boundary_str,"freeslip") == 0)
    {
        if((local_bt&left_boundary)||(local_bt&right_boundary))
        {
            return ght_reverse_vel_x;
        } else
        {
            return ght_reverse_vel_y;
        }
    }
    else
    {
        return ght_zero;
    }
}


void set_corerank_map(mesh2d_info * _minfo)
{
    // physical processor number and rank size
    int p_nprcs,l_nprcs;
    MPI_Comm_size(MPI_COMM_WORLD,&p_nprcs);
    l_nprcs = _minfo->npgx * _minfo->npgy;
    if(p_nprcs >= l_nprcs) return;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    int x_rank = myrank/_minfo->npgy;
    int y_rank = myrank%_minfo->npgy;

    int num_full_core = (int) floor(_minfo->oy * _minfo->npgy);
    int num_res_core  = _minfo->npgy - num_full_core;
    int num_full_size = num_full_core*_minfo->npgx;

    int core_rank;
    if(y_rank < num_full_core)
    {
        core_rank = x_rank*num_full_core + y_rank;
    } else if(num_res_core%2 == 0)
    {
        int res_y_rank = (y_rank - num_full_core)/2;
        core_rank = x_rank*(num_res_core/2) + res_y_rank + num_full_size;
    } else
    {
        int res_y_rank = (y_rank - num_full_core + 1)/2;
        core_rank = x_rank*(num_res_core/2 + 1) + res_y_rank + num_full_size;
    }

    char env_string[256];
    sprintf(env_string, "OMP_PROC_BIND=true,OMP_PLACES=cores,OMP_NUM_THREADS=1,OMP_PROC_BIND_VERBOSE=1,OMP_DISPLAY_ENV=true,OMP_PLACES=%d", core_rank);
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "env", env_string);
    MPI_Comm_set_info(MPI_COMM_WORLD, info);

}

void make_empty_dir(const char * dir)
{
    if(access(dir, F_OK) == 0)
    {
        DIR * _dir = opendir(dir);
        struct dirent * _tmp;
        while(1)
        {
            _tmp = readdir(_dir);
            if(NULL == _tmp) break;
            if(strcasecmp(".",_tmp->d_name)==0 || strcasecmp("..",_tmp->d_name)==0) continue;
            if(_tmp->d_type == DT_REG)
            {
                char _tmp_name[4096];
                snprintf(_tmp_name,4096,"%s/%s",dir,_tmp->d_name);
                remove(_tmp_name);
            }

        }
        closedir(_dir);
    } else
    {
        mkdir(dir, 0700);
    }
}

void put_cpuinfo(FILE * fp)
{
    FILE* cpuinfo = fopen("/proc/cpuinfo", "r");
    if (cpuinfo == NULL) {
        fprintf(stderr, "Failed to open /proc/cpuinfo\n");
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), cpuinfo)) {
        fprintf(fp, "%s", line);
    }
    fclose(cpuinfo);
}