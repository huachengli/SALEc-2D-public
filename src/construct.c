//
// Created by huach on 3/7/2023.
//
// functions related to grid geometry
// calculate grid parameters of grad
//           cell/vertex volume
//
#include "sale2d.h"
#include "interpolate_2d9.h"

void init_mesh_pos(sale2d_var * _sale, mesh2d_info * _minfo)
{
    /*
     * set the position of vertex
     * and derive the geometry parameters
     */


    int my_rank = _sale->v_pos->rank;
    int num_y_cells = _minfo->npgy*_minfo->ny + 2*_minfo->ypad;
    int YeL = _minfo->ny_ext[0], YeR = _minfo->ny_ext[1];
    double Yext = _minfo->y_ext;

    int num_x_cells = _minfo->npgx*_minfo->nx + 2*_minfo->xpad;
    double Xext = _minfo->x_ext;
    int XeL = _minfo->nx_ext[0], XeR = _minfo->nx_ext[1];

    int XdR = _minfo->damp_x, XdL = _minfo->damp_x_left;
    int YdL = _minfo->damp_y, YdR = _minfo->damp_y_top;

    double Ly = SpacingY(_minfo->npgy*_minfo->ny + _minfo->ypad,Yext,YeL,YeR,YdL,YdR,num_y_cells)
                - SpacingY(_minfo->ypad,Yext,YeL,YeR,YdL,YdR,num_y_cells);
    int y0_i;
    for(y0_i = 0; y0_i < num_y_cells+1; ++y0_i)
    {
        if(SpacingY(y0_i + _minfo->ypad,Yext,YeL,YeR,YdL,YdR,num_y_cells) >= Ly*_minfo->oy) break;
    }

    double y0 =  - SpacingY(y0_i + _minfo->ypad,Yext,YeL,YeR,YdL,YdR,num_y_cells) * _minfo->dy;


    // fix x0 in vertex
    /*
     * x0 is must be put at high resolution area
     */
    double Lx = SpacingX(_minfo->npgx*_minfo->nx + _minfo->xpad,Xext,XeL,XeR,XdL,XdR,num_x_cells,_minfo->xpad)
                - SpacingX(_minfo->xpad,Xext,XeL,XeR,XdL,XdR,num_x_cells,_minfo->xpad);


    int x0_i;
    for(x0_i = 0;x0_i<num_x_cells+1;++x0_i)
    {
        if(SpacingX(x0_i,Xext,XeL,XeR,XdL,XdR,num_x_cells,_minfo->xpad) >= Lx*_minfo->ox) break;
    }

    double x0 = - SpacingX(x0_i + _minfo->xpad,Xext,XeL,XeR,XdL,XdR,num_x_cells,_minfo->xpad) * _minfo->dx;

    // first loop on vertex just set its coordinates
    for(int k=0;k<_minfo->v_nx;++k)
    {
        for(int j=0;j<_minfo->v_ny;++j)
        {
            fv_id_2d tmp_lid = {.x=k,.y=j,.u=0,.v=0};
            fv_id_2d tmp_gid = {0};
            double * vpos_ptr;

            _sale->fov->lid2gid(my_rank,&tmp_lid,&tmp_gid,_sale->v_pos);
            _sale->fov->get_double(&tmp_lid,&vpos_ptr,_sale->v_pos);

            vpos_ptr[0]   = SpacingX(tmp_gid.x,Xext,XeL,XeR,XdL,XdR,num_x_cells,_minfo->xpad) * _minfo->dx + x0;
            vpos_ptr[1]   = SpacingY(tmp_gid.y,Yext,YeL,YeR,YdL,YdR,num_y_cells) * _minfo->dy + y0;
        }
    }

    // loop on full range of gid.y, set the lub_y0/y1
    if(_sale->tracer_lubrication >= 1)
    {
        _sale->tracer_lub_gy0 = 0;
        _sale->tracer_lub_gy1 = num_y_cells;
        int j = 0;
        for(;j<num_y_cells+1;++j)
        {
            double tmpy = SpacingY(j,Yext,YeL,YeR,YdL,YdR,num_y_cells) * _minfo->dy + y0;
            _sale->tracer_lub_gy0 = j;
            if(tmpy >= _sale->tracer_lubrication*_minfo->dy*(-1.0)) break;
        }
        for(;j<num_y_cells+1;++j)
        {
            double tmpy = SpacingY(j,Yext,YeL,YeR,YdL,YdR,num_y_cells) * _minfo->dy + y0;
            _sale->tracer_lub_gy1 = j;
            if(tmpy >= 0.0) break;
        }
    }


    for(int k=0;k<_minfo->e_nx;++k)
    {
        for(int j=0;j<_minfo->e_ny;++j)
        {
            fv_id_2d tmp_lid = {.x=k,.y=j,.u=0,.v=0};
            fv_id_2d v_lb = {.x=k,.y=j,.u=0,.v=0};
            fv_id_2d v_rb = {.x=k+1,.y=j,.u=0,.v=0};
            fv_id_2d v_ru = {.x=k+1,.y=j+1,.u=0,.v=0};
            fv_id_2d v_lu = {.x=k,.y=j+1,.u=0,.v=0};

            char *d_lb, *d_rb, *d_lu, *d_ru;
            _sale->fov->get_data(&v_lb,&d_lb,_sale->v_pos);
            _sale->fov->get_data(&v_rb,&d_rb,_sale->v_pos);
            _sale->fov->get_data(&v_lu,&d_lu,_sale->v_pos);
            _sale->fov->get_data(&v_ru,&d_ru,_sale->v_pos);


            char *de_pos;
            _sale->foe->get_data(&tmp_lid,&de_pos,_sale->e_pos);

            double * de_pos2 = (double *) de_pos;
            double * d_lb2 = (double *) d_lb;
            double * d_rb2 = (double *) d_rb;
            double * d_ru2 = (double *) d_ru;
            double * d_lu2 = (double *) d_lu;

            // set the centre of element
            char * suf, *vol, *subvol;
            _sale->foe->get_data(&tmp_lid,&suf,_sale->e_suf);
            _sale->foe->get_data(&tmp_lid,&vol,_sale->e_vol);
            _sale->foe->get_data(&tmp_lid,&subvol,_sale->e_subvol);

            // set the surface and volume of element
            generate2dcyl_suf_and_vol((double []){
                                              d_lb2[0],d_lb2[1],d_rb2[0],d_rb2[1],
                                              d_ru2[0],d_ru2[1],d_lu2[0],d_lu2[1],},
                                      (double *)suf,(double *)subvol,(double *)vol);

            //set gradient parameter of variables
            char * grad_cptr;
            _sale->foe->get_data2(&tmp_lid,&grad_cptr,_sale->e_grad);
            generate2dcyl_grad((double []){
                    d_lb2[0],d_lb2[1],d_rb2[0],d_rb2[1],
                    d_ru2[0],d_ru2[1],d_lu2[0],d_lu2[1],},(double *) grad_cptr);

            //set the center and scale of element
            generate2d_center((double []){
                    d_lb2[0],d_lb2[1],d_rb2[0],d_rb2[1],
                    d_ru2[0],d_ru2[1],d_lu2[0],d_lu2[1],}, de_pos2);

            fv_id_2d tmp_gid = {0,0};
            _sale->foe->lid2gid(my_rank,&tmp_lid,&tmp_gid,_sale->e_vof);
            double * ealpha_ptr,*esigma_ptr;
            _sale->foe->get_double(&tmp_lid,&ealpha_ptr,_sale->e_alpha);
            _sale->foe->get_double(&tmp_lid,&esigma_ptr,_sale->e_sigma);
            vec_zero(ealpha_ptr,2);
            vec_zero(esigma_ptr,2);
            clean_fv(_sale->e_beta);

            double delta_len_y = Spacing2(_minfo->damp_y,Yext,YeL,YeR,YdL,YdR,num_y_cells) * _minfo->dy;
            double delta_len_x = Spacing2(num_x_cells,Xext,XeL,XeR,XdL,XdR,num_x_cells) * _minfo->dx -
                                 Spacing2(num_x_cells - _minfo->damp_x,Xext,XeL,XeR,XdL,XdR,num_x_cells) * _minfo->dx;
            double delta_y   = de_pos2[Y] - y0;
            double delta_x   = Spacing2(num_x_cells,Xext,XeL,XeR,XdL,XdR,num_x_cells) * _minfo->dx - (de_pos2[X] - x0);
            double csd_approx = 5.0e3;
            double sbc_b = 0.2, pml_b = 2.3;

            switch (_minfo->damp_type)
            {
                case pml:
                case pml2:
                    esigma_ptr[Y] = esigma_ptr[X] = 0.;

                    if(delta_y <= delta_len_y)
                    {
                        esigma_ptr[Y] = pow(1.0 - delta_y/delta_len_y,2.0)/(2.0*delta_len_y)*(3.0*csd_approx)*_minfo->lambda_y;
                    }

                    if(delta_x <= delta_len_x)
                    {
                        esigma_ptr[X] = pow(1.0 - delta_x/delta_len_x,2.0)/(2.0*delta_len_x)*(3.0*csd_approx)*_minfo->lambda_x;
                    }
                    _sale->alpha_damp_x = (csd_approx*2.0)/delta_len_x * pml_b;
                    _sale->alpha_damp_y = (csd_approx*2.0)/delta_len_y * pml_b;
                    break;
                case pml3:
                    esigma_ptr[Y] = esigma_ptr[X] = 0.;

                    if(delta_y <= delta_len_y)
                    {
                        esigma_ptr[Y] = pow(1.0 - delta_y/delta_len_y,2.0)/(2.0*delta_len_y)*(3.0*csd_approx)*_minfo->lambda_y;
                    }

                    if(delta_x <= delta_len_x)
                    {
                        esigma_ptr[X] = pow(1.0 - delta_x/delta_len_x,2.0)/(2.0*delta_len_x)*(3.0*csd_approx)*_minfo->lambda_x;
                    }
                    _sale->alpha_damp_x = (csd_approx*2.0)/delta_len_x * pml_b;
                    _sale->alpha_damp_y = (csd_approx*2.0)/delta_len_y * pml_b;
                    break;
                case sbc:
                    esigma_ptr[Y] = esigma_ptr[X] = 1.0;

                    if(delta_y <= delta_len_y)
                    {
                        esigma_ptr[Y] = hanning(0.5*((1.0 - _minfo->lambda_y)*delta_y/delta_len_y + _minfo->lambda_y));
                    }

                    if(delta_x <= delta_len_x)
                    {
                        esigma_ptr[X] = hanning(0.5*((1.0 - _minfo->lambda_x)*delta_x/delta_len_x + _minfo->lambda_x));
                    }

                    esigma_ptr[X] = Min(esigma_ptr[X],esigma_ptr[Y]);
                    break;
                default:
                    esigma_ptr[X] = esigma_ptr[Y] = 0.;
                    if(delta_y <= delta_len_y) esigma_ptr[Y]=1.0;
                    if(delta_x <= delta_len_x) esigma_ptr[X]=1.0;
                    break;
            }

        }
    }

    /*
     * make another loop on vertex to collect volume information of element
     * just reduce the memory loading in calculating
     */

    for(int k=0;k<_minfo->v_nx;++k)
    {
        for(int j=0;j<_minfo->v_ny;++j)
        {
            fv_id_2d tmp_lid = {.x=k,.y=j};
            node_type local_node_type = useless_section;
            _sale->fov->get_nty(&tmp_lid,&local_node_type,_sale->v_vol);
            // skip the useless vertex (they are on the boundary of domain, just have 2 element neighbours)
            if(useless_section & local_node_type) continue;
            fv_id_2d lb={k-1,j-1}, rb={k,j-1}, ru={k,j}, lu={k-1,j};
            fv_id_2d * elid[4] = {&lb,&rb,&ru,&lu};

            double * esubvol_ptr[4], *esuf_ptr[4];
            for(int i=0;i<4;++i)
            {
                _sale->foe->get_double(elid[i],&esubvol_ptr[i],_sale->e_subvol);
                _sale->foe->get_double(elid[i],&esuf_ptr[i],_sale->e_suf);
            }

            double * vvol_ptr, *vsubvol_ptr, *vsuf_ptr;
            _sale->fov->get_double(&tmp_lid,&vvol_ptr,_sale->v_vol);
            _sale->fov->get_double(&tmp_lid,&vsubvol_ptr,_sale->v_subvol);
            _sale->fov->get_double(&tmp_lid,&vsuf_ptr,_sale->v_suf);

            vvol_ptr[0] = 0.;

            int v_offset[4] = {2, 3, 0, 1};
            for(int i=0;i<4;++i)
            {
                vvol_ptr[0]    += esubvol_ptr[i][4 + v_offset[i]];
                vsubvol_ptr[i]  = esubvol_ptr[i][v_offset[i]];
                vvol_ptr[i + 1] = esubvol_ptr[i][v_offset[i]+4];
                vsuf_ptr[2*i    ]  = esuf_ptr[i][v_offset[i]*2];
                vsuf_ptr[2*i + 1]  = esuf_ptr[i][v_offset[i]*2 + 1];
            }

            double * valpha_ptr;
            _sale->fov->get_double(&tmp_lid,&valpha_ptr,_sale->v_alpha);
            vec_zero(valpha_ptr,2);
        }
    }

    /*
     * after position has been set, calculate the laplacian operator
     * of vertex.
     */

    for(int k=0;k<_minfo->v_nx;++k)
    {
        for(int j=0;j<_minfo->v_ny;++j)
        {
            fv_id_2d tmp_lid = {.x=k,.y=j};
            node_type local_node_type = useless_section;
            _sale->fov->get_nty(&tmp_lid,&local_node_type,_sale->v_vol);
            if(useless_section & local_node_type) continue; // skip the physical boundary vertexes

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

            double *vpos_ptr;
            double tmp_vpos[18] = {0.};
            for(int kv=0;kv<9;++kv)
            {
                _sale->fov->get_double(vlid+kv,&(vpos_ptr),_sale->v_pos);
                tmp_vpos[2*kv + 0] = vpos_ptr[0];
                tmp_vpos[2*kv + 1] = vpos_ptr[1];
            }
            double *vlaplace_ptr;
            _sale->fov->get_double(&tmp_lid,&vlaplace_ptr,_sale->v_laplace);
            generate2d_laplace(tmp_vpos,vlaplace_ptr);
        }
    }

    /*
     * set some info required on bx/by
     */
    proc_info_2d * proc_bx = _sale->bx_normal->proc_self;
    int bx_nx = proc_bx->nx, bx_ny = proc_bx->ny;

    for(int j=0;j<bx_ny;++j)
    {
        for(int k=0;k<bx_nx;k++)
        {
            fv_id_2d blid = {k,j};
            node_type local_node_type = useless_section;
            _sale->fobx->get_nty(&blid,&local_node_type,_sale->bx_normal);
            if(local_node_type & useless_section) continue;
            fv_id_2d u_elid = {k,j}, b_elid={k-1,j};
            double * u_epos_ptr, *b_epos_ptr;
            _sale->foe->get_double(&u_elid,&u_epos_ptr,_sale->e_pos);
            _sale->foe->get_double(&b_elid,&b_epos_ptr,_sale->e_pos);

            double *bxnormal_ptr;
            _sale->fobx->get_double(&blid,&bxnormal_ptr,_sale->bx_normal);
            liner_op(bxnormal_ptr,2,u_epos_ptr,b_epos_ptr,1.0,-1.0);
            vec_norm(bxnormal_ptr,2);
        }
    }

    proc_info_2d * proc_by = _sale->by_normal->proc_self;
    int by_nx = proc_by->nx, by_ny = proc_by->ny;
    for(int j=0;j<by_ny;++j)
    {
        for(int k=0;k<by_nx;++k)
        {
            fv_id_2d blid = {k,j};
            node_type local_node_type = useless_section;
            _sale->foby->get_nty(&blid,&local_node_type,_sale->by_normal);
            if(local_node_type & useless_section) continue;
            fv_id_2d r_elid = {k,j}, l_elid={k,j-1};
            double * r_epos_ptr, * l_epos_ptr;
            _sale->foe->get_double(&r_elid,&r_epos_ptr,_sale->e_pos);
            _sale->foe->get_double(&l_elid,&l_epos_ptr,_sale->e_pos);

            double *bynormal_ptr;
            _sale->foby->get_double(&blid,&bynormal_ptr,_sale->by_normal);
            liner_op(bynormal_ptr,2,r_epos_ptr,l_epos_ptr,1.0,-1.0);
            vec_norm(bynormal_ptr,2);
        }
    }

    /*
     * output info of grid to txt
     */
    if(_sale->e_pos->rank == 0)
    {
        proc_info_2d * proc_main = (proc_info_2d *) _sale->e_pos->domain;
        int ngx = proc_main->nx;
        int ngy = proc_main->ny;
        FILE * grid_fp = fopen("txt/grid.txt","w");
        fprintf(grid_fp,"%d %d\n",ngx,ngy);
        for(int k=0;k<ngx;++k)
        {
            fprintf(grid_fp,"%10.5e\n",SpacingX(k,Xext,XeL,XeR,XdL,XdR,num_x_cells,_minfo->xpad) * _minfo->dx + x0);
        }
        for(int j=0;j<ngy;++j)
        {
            fprintf(grid_fp,"%10.5e\n",SpacingY(j,Yext,YeL,YeR,YdL,YdR,num_y_cells) * _minfo->dy + y0);
        }
        fclose(grid_fp);
    }

}

double Sk(int k,double ext)
{
    return ext*(pow(ext,k) - 1.0)/(ext - 1.0);
}

double Spacing(int k,double ext, int eL, int eR, int N)
{
    double Xk = 0.0;
    if(k<eL)
        Xk = Sk(eL,ext) - Sk(eL-k,ext);
    else if(k<N-eR)
        Xk = Sk(eL,ext) + (k - eL);
    else
        Xk = Sk(eL,ext) + (N - eL - eR) + Sk(k-N+eR,ext);
    return Xk;
}

double Spacing2(int k,double ext, int eL,int eR, int dL,int dR,int N)
{
    double Xk = 0.;
    int rk = N - 1 - k;
    if(k < dL)
        Xk = pow(ext,eL)*k;
    else if(k < dL + eL)
        Xk = Sk(eL,ext) - Sk(dL+eL-k,ext) + pow(ext,eL)*dL;
    else if(k < N  - (dR + eR))
        Xk = Sk(eL,ext) + pow(ext,eL)*dL + (k - (dL + eL));
    else if(k < N  - dR)
        Xk = Sk(eL,ext) + pow(ext,eL)*dL + (N - (dL + eL) - (dR + eR)) + Sk(k - (N - dR - eR),ext);
    else
        Xk = Sk(eL,ext) + pow(ext,eL)*dL + (N  - (dL + eL) - (dR + eR)) + Sk(eR,ext) + pow(ext,eR)*(k - (N - dR));
    return Xk;
}


/*
 * Spacing3 is used to set coordinate at r-direction
 * and do not support extension grid at left ..
 */

double cyl_dh(double r1)
{
    /*
     * generate the coordinate of point at right of r1
     */
    return 1.0 + sqrt(1.0 + r1*r1);

}

double Spacing3(int k,double ext, int eR, int dR,int N, int rpad)
{
    const int num_approx_points = 30;
    const double cyl_approx_points[] = {
            0.000000000,    2.000000000,    3.236067977,    4.387054171,    5.499582680,
            6.589759356,    7.665202800,    8.730157435,    9.787243529,   10.838197797,
            11.884233160,   12.926231501,   13.964854832,   15.000613218,   16.033908238,
            17.065061886,   18.094336406,   19.121948294,   20.148078404,   21.172879402,
            22.196481363,   23.218996037,   24.240520153,   25.261138005,   26.280923505,
            27.299941831,   28.318250749,   29.335901705,   30.352940719,   31.369409120,
            32.385344168,   33.400779572,   34.415745930,   35.430271098,   36.444380515,
            37.458097470,   38.471443341,   39.484437799,   40.497098986,   41.509443668,
    };

    if(N - eR - dR < num_approx_points)
    {
        return Spacing2(k,ext,0,eR,0,dR,N);
    }

    if(k < rpad)
    {
        return 2.0*k;
    }
    else if(k - rpad < num_approx_points)
    {
        return cyl_approx_points[k - rpad] + 2.0*rpad;
    }
    else
    {
        return Spacing2(k ,ext,0,eR,0,dR,N) - Spacing2(num_approx_points+rpad-1,ext,0,eR,0,dR,N) + cyl_approx_points[num_approx_points - 1] + 2.0*rpad;
    }
}

double SpacingX(int k,double ext, int eL,int eR, int dL,int dR,int N,int rpad)
{
#ifndef UNIFORM_XGRID
    return Spacing3(k,ext,eR,dR,N,rpad);
#else
    return Spacing2(k,ext,eL,eR,dL,dR,N);
#endif
}

double SpacingY(int k,double ext, int eL,int eR, int dL,int dR,int N)
{
    return Spacing2(k,ext,eL,eR,dL,dR,N);
}

double hanning(double x)
{
    return 0.5*(1.0 - cos(2.0*M_PI*x));
}

void generate2dcyl_suf_and_vol(const double * _pos, double * _suf, double * _vol_sub, double * _vol_cyl)
{
    /*
     * _suf is the surface of element contributed to vertex, length of _suf is 4*2
     * _vol_cyl is the volume in Cylindrical coordinate, length of _vol_cyl is 1
     * _vol_sub is the volume calculated in cartesian coordinate, length of _vol_sub is 4 (number of sub-volume)
     * in regin r < 0, _vol_sub should be negative,
     * _vol_sub[0,1,2,3] represent the term Srz*e_r, in r < 0, direction of e_r is point to negative
     * _vol_sub[4,5,6,7] is the volume contributed to vertex
     */
    // local coordinate of sub element
    double local_xl[4][4][2] = {
            /*
             * left bottom, right bottom, right up, left up
             */
            {{-1.0,-1.0}, { 0.0,-1.0}, { 0.0, 0.0}, {-1.0, 0.0}},
            {{ 0.0,-1.0}, { 1.0,-1.0}, { 1.0, 0.0}, { 0.0, 0.0}},
            {{ 0.0, 0.0}, { 1.0, 0.0}, { 1.0, 1.0}, { 0.0, 1.0}},
            {{-1.0, 0.0}, { 0.0, 0.0}, { 0.0, 1.0}, {-1.0, 1.0}}
    };
    double local_xg[4][2];
    double sub_vol[4];
    double sub_rvol[4];
    double sign_er = _pos[0] + _pos[2] + _pos[4] + _pos[6];
    for(int k=0;k<4;++k)
    {
        // set local coordinate
        for(int j=0;j<4;++j) Interpolate2dvec(_pos,local_xl[k][j],local_xg[j]);

        // calculate volume in Cylindrical coordinate
        sub_rvol[k] = Integrate2dcyl(local_xg[0],(double [4]){1.,1.,1.,1.});
        // calculate volume in cartesian coordinate
        sub_vol[k] = Integrate2d(local_xg[0],(double [4]){1.,1.,1.,1.});
        if(sign_er < 0.0) sub_vol[k] *= -1.0;
        // calculate sub surface
        for(int j=0;j<(k+1);++j) Polygon2dRotate(local_xg[0],4);
        Polygon2dcyl(local_xg[0], 3, _suf + k * 2);
    }

    for(int k=0;k<4;++k) sub_rvol[k] = fabs(sub_rvol[k]);
    _vol_cyl[0] = fabs(sub_rvol[0]) + fabs(sub_rvol[1]) + fabs(sub_rvol[2]) + fabs(sub_rvol[3]);
    vec_copy(sub_vol  , _vol_sub, 4);
    vec_copy(sub_rvol ,_vol_sub+4,4);
}

void generate2dcyl_grad(const double * _pos, double * _grad)
{
    double local_corner[4][4] = {
            {1., 0., 0., 0.}, /* left bottom, right bottom, right up, left up */
            {0., 1., 0., 0.},
            {0., 0., 1., 0.},
            {0., 0., 0., 1.}
    };

    double tmp_vol = Integrate2dcyl(_pos,(double [4]){1.,1.,1.,1.});
    for(int j=0;j<4;++j)
    {
        for(int k=0;k<2;++k)
        {
            _grad[4*k + j] = Integrate2dcyl_Fdx(_pos,local_corner[j],k);
        }
    }

    /*
     * factor for v_r/r in divergence
     * in r<0, coordinate line is reversed
     */
    double sign_er = _pos[0] + _pos[2] + _pos[4] + _pos[6];

    for(int j=0;j<4;++j)
    {
        _grad[8+j] = Integrate2d(_pos,local_corner[j]);
        if(sign_er < 0.) _grad[8+j] *= -1.0;
    }

    vec_scaler(_grad,12,1.0/tmp_vol);
}

void generate2d_center(const double * _pos, double * _center)
{
    /*
     * _center =
     * { center of X direction, center of Y direction,
     *   scale  of X direction, scale  of Y direction }
     */
    vec_zero(_center,4);
    for(int k=0;k<4;++k)
    {
        //some case the center should be the integrated average of
        //corner coordinates
        scaler_add(_center,2,_pos +2*k,0.25);
    }

    // find the max scale of element
    for(int k=0;k<4;++k)
    {
        for(int j=0;j<4;++j)
        {
            _center[2] = MAX(fabs(_pos[2*k + 0]-_pos[2*j + 0]),_center[2]);
            _center[3] = MAX(fabs(_pos[2*k + 1]-_pos[2*j + 1]),_center[3]);
        }
    }
}

void generate2d_laplace(const double * _pos, double * _vl)
{
    double local_corner[9] = {0.};
    for(int k=0;k<9;++k)
    {
        vec_zero(local_corner,9);
        local_corner[k] = 1.;
        _vl[k] = Laplace2d9(_pos,local_corner,LAPLICAN_GAMMA);
    }

}
