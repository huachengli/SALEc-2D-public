//
// Created by huacheng on 5/30/23.
//

#include "write_vtk.h"


void sale_array_vec_binary(FILE * fp,const char *name,field_var * fv)
{
    vtk_dataarray_vec(fp,name,"binary",(double *)fv->data,fv->data_size,fv->bytes_per_unit/ sizeof(double ));
}

void sale_array_vec_ascii(FILE * fp,const char *name,field_var * fv)
{
    vtk_dataarray_vec(fp,name,"ascii",(double *)fv->data,fv->data_size,fv->bytes_per_unit/ sizeof(double ));
}


void pvtk_array_description(FILE * fp,const char *name,field_var * fv)
{
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n",
            name,fv->bytes_per_unit/ sizeof(double));
}

void sale2d_write_vtm(int nproc, int cycles, const char _prefix[])
{
    const char header[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <vtkMultiBlockDataSet>\n";


    // prepare vtm file for vts
    char fname[50];
    sprintf(fname,"%s.%04d.vtm",_prefix,cycles);
    FILE * fp = fopen(fname,"w");
    fputs(header, fp);

    for(int k=0; k<nproc; k++) {
        fprintf(fp, "    <DataSet index=\"%d\" file=\"./vts/%s.proc%04d.%04d.vts\"/>\n",
                k,_prefix, k, cycles);
    }
    fputs("  </vtkMultiBlockDataSet>\n",fp);
    fputs("</VTKFile>",fp);
    fclose(fp);

    // prepare vtm file for vtp
    sprintf(fname,"%s.tracer.%04d.vtm",_prefix,cycles);
    fp = fopen(fname,"w");
    fputs(header, fp);

    for(int k=0; k<nproc; k++) {
        fprintf(fp, "    <DataSet index=\"%d\" file=\"./vtp/%s.tracer.proc%04d.%04d.vtp\"/>\n",k,_prefix, k, cycles);
    }
    fputs("  </vtkMultiBlockDataSet>\n",fp);
    fputs("</VTKFile>",fp);
    fclose(fp);
}

void sale2d_write_pvtk(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl)
{
    /*
     * in most time we can use the vtm file
     * the sequence of variable written in pvts should same with vts
     */
    int cycles = _ccl->k_out_step - 1;
    char fname[256];
    sprintf(fname,"%s.%04d.pvts",_minfo->prefix,cycles);
    FILE * fp = fopen(fname,"w");

    const char header[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"PStructuredGrid\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n";
    fputs(header, fp);

    // construct information on PStructuredGrid
#ifdef OUTPUT_GHT
    fprintf(fp,"<PStructuredGrid WholeExtent=\"%d %d %d %d 0 0\" GhostLevel=\"%d\">\n",
            _minfo->xpad,_minfo->npgx*_minfo->nx+2*_minfo->xpad,0,_minfo->npgy*_minfo->ny+_minfo->ypad,0);
#else
    fprintf(fp,"<PStructuredGrid WholeExtent=\"%d %d %d %d 0 0\" GhostLevel=\"%d\">\n",
            _minfo->xpad,_minfo->npgx*_minfo->nx+_minfo->xpad,_minfo->ypad,_minfo->npgy*_minfo->ny+_minfo->ypad,0);
#endif
    pvtk_point_data_header(fp);
    pvtk_array_description(fp,"v_vel",_sale->v_vel);
    pvtk_array_description(fp,"v_acc",_sale->v_acc);
    pvtk_array_description(fp,"v_gra",_sale->v_bf);
    pvtk_array_description(fp,"v_vof",_sale->v_vof);

#ifdef DEBUG_MODE
    pvtk_array_description(fp,"v_mass",_sale->v_mass);
    pvtk_array_description(fp,"v_dis",_sale->v_dis);
    pvtk_array_description(fp,"v_vol",_sale->v_vol);
    pvtk_array_description(fp,"v_subvol",_sale->v_subvol);
    pvtk_array_description(fp,"v_suf",_sale->v_suf);
    pvtk_array_description(fp,"v_laplacian",_sale->v_laplace);
    pvtk_array_description(fp,"v_debug",_sale->v_debug);
#endif
    pvtk_point_data_trailer(fp);

    pvtk_cell_data_header(fp);
    pvtk_array_description(fp,"e_den",_sale->e_den);
    pvtk_array_description(fp,"e_pre",_sale->e_pre);
    pvtk_array_description(fp,"e_wpt",_sale->e_wpt);
    pvtk_array_description(fp,"e_vof",_sale->e_vof);
    pvtk_array_description(fp,"e_csd",_sale->e_csd);
    pvtk_array_description(fp,"e_dam",_sale->e_dam);
    pvtk_array_description(fp,"e_tps",_sale->e_tps);
    pvtk_array_description(fp,"e_q",_sale->e_q);
    pvtk_array_description(fp,"e_tem",_sale->e_tem);
    pvtk_array_description(fp,"e_eng",_sale->m_eng);
    pvtk_array_description(fp,"m_den",_sale->m_den);
    pvtk_array_description(fp,"e_ste2s",_sale->e_ste2s);
    pvtk_array_description(fp,"e_sta2s",_sale->e_sta2s);
    pvtk_array_description(fp,"e_ste",_sale->e_ste);
    pvtk_array_description(fp,"e_strength",_sale->e_strength);
    pvtk_array_description(fp,"e_vib",_sale->e_vib);
    pvtk_array_description(fp,"e_ejt",_sale->e_ejt);

#ifdef DEBUG_MODE
    pvtk_array_description(fp,"edebug",_sale->e_debug);
    pvtk_array_description(fp,"esuf",_sale->e_suf);
    pvtk_array_description(fp,"esubvol",_sale->e_subvol);
    pvtk_array_description(fp,"egrad",_sale->e_grad);
    pvtk_array_description(fp,"divvel",_sale->e_div_vel);
    pvtk_array_description(fp,"e_alpha",_sale->e_alpha);
    pvtk_array_description(fp,"e_sigma",_sale->e_sigma);
    pvtk_array_description(fp,"e_pos",_sale->e_pos);
    pvtk_array_description(fp,"e_vol",_sale->e_vol);
    pvtk_array_description(fp,"e_vel",_sale->e_vel);
    pvtk_array_description(fp,"e_dampref",_sale->e_dampref);
    pvtk_array_description(fp,"e_gradvel",_sale->e_grad_vel);
    pvtk_array_description(fp,"e_cnd",_sale->e_cnd);
#endif

    pvtk_cell_data_trailer(fp);

    // add coordinate tag
    fprintf(fp,"    <PPoints>\n"
               "      <PDataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\"/>\n"
               "    </PPoints>\n");

    int n_prcs = _minfo->npgx*_minfo->npgy;
    for(int k=0;k<n_prcs;++k)
    {
        int x_rank = k/_minfo->npgy;
        int y_rank = k%_minfo->npgy;

#ifdef OUTPUT_GHT
        int x_start = _minfo->xpad, x_end = _minfo->nx+_minfo->xpad;
        int y_start = _minfo->ypad, y_end = _minfo->ny+_minfo->ypad;
        if(x_rank <= 0)
        {
            x_start = 0;
        }
        else if(x_rank >= _minfo->npgx-1)
        {
            x_end = _minfo->nx+2*_minfo->xpad;
        }

        if(y_rank <= 0)
        {
            y_start = 0;

        } else if(y_rank >= _minfo->npgy-1)
        {
            y_end = _minfo->ny+2*_minfo->ypad;
        }

        proc_info_2d *proc_self = (proc_info_2d *) (_sale->v_pos->proc) + k;
        int stx = proc_self->stx;
        int sty = proc_self->sty;
        fprintf(fp,"      <Piece Extent=\"%d %d %d %d 0 0\" Source=\"./vts/%s.proc%04d.%04d.vts\"/>\n",
                x_start + stx,x_end + stx,y_start + sty,y_end + sty,
                _minfo->prefix,k,cycles);
#else
        fprintf(fp,"      <Piece Extent=\"%d %d %d %d 0 0\" Source=\"%s.proc%04d.%04d.vts\"/>\n",
                x_rank*_minfo->nx+_minfo->xpad,(x_rank+1)*_minfo->nx+_minfo->xpad,
                y_rank*_minfo->ny+_minfo->ypad,(y_rank+1)*_minfo->ny+_minfo->ypad,
                _minfo->prefix,k,cycles);
#endif
    }

    fputs("</PStructuredGrid>\n",fp);
    fputs("</VTKFile>",fp);
    fclose(fp);
}

void sale2d_write_vtk(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl)
{
    char output_name[200];
    int myrank = _ccl->rank;

    /*
     * get the start of ext distribution
     *
     */
    proc_info_2d *proc_self = (proc_info_2d *) _sale->v_pos->proc_self;
    int stx = proc_self->stx;
    int sty = proc_self->sty;

    sprintf(output_name,"./vts/%s.proc%04d.%04d.vts",_minfo->prefix,myrank,_ccl->k_out_step-1);
    FILE * fp = fopen(output_name,"w");
    if(NULL == fp)
    {
        fprintf(stdout,"cannot create %s\n",output_name);
        exit(0);
    }

    char piece_extent[64], whole_extent[64];
    sprintf(piece_extent,"%d %d %d %d 0 0",
            0 + stx,_minfo->v_nx-1 + stx,
            0 + sty,_minfo->v_ny-1 + sty);


    int x_start = _minfo->xpad, x_end = _minfo->nx+_minfo->xpad;
    int y_start = _minfo->ypad, y_end = _minfo->ny+_minfo->ypad;

#ifdef OUTPUT_GHT
    if(_ccl->x_rank <= 0)
    {
        x_start = 0;
    }
    else if(_ccl->x_rank >= _minfo->npgx-1)
    {
        x_end = _minfo->nx+2*_minfo->xpad;
    }

    if(_ccl->y_rank <= 0)
    {
        y_start = 0;

    } else if(_ccl->y_rank >= _minfo->npgy-1)
    {
        y_end = _minfo->ny+2*_minfo->ypad;
    }
#endif

    sprintf(whole_extent,"%d %d %d %d 0 0",
            x_start + stx,x_end + stx,
            y_start + sty,y_end + sty);

    vts_file_header(fp,piece_extent,whole_extent);
    vtk_point_data_header(fp);
    sale_array_vec_binary(fp,"v_vel",_sale->v_vel);
    sale_array_vec_binary(fp,"v_acc",_sale->v_acc);
    sale_array_vec_binary(fp,"v_gra",_sale->v_bf);
    sale_array_vec_binary(fp,"v_vof",_sale->v_vof);

#ifdef DEBUG_MODE
    sale_array_vec_binary(fp,"v_mass",_sale->v_mass);
    sale_array_vec_binary(fp,"v_dis",_sale->v_dis);
    sale_array_vec_binary(fp,"v_vol",_sale->v_vol);
    sale_array_vec_binary(fp,"v_subvol",_sale->v_subvol);
    sale_array_vec_binary(fp,"v_suf",_sale->v_suf);
    sale_array_vec_binary(fp,"v_laplacian",_sale->v_laplace);
    sale_array_vec_binary(fp,"v_debug",_sale->v_debug);
#endif
    vtk_point_data_trailer(fp);
    vtk_cell_data_header(fp);
    sale_array_vec_binary(fp,"e_den",_sale->e_den);
    sale_array_vec_binary(fp,"e_pre",_sale->e_pre);
    sale_array_vec_binary(fp,"e_wpt",_sale->e_wpt);
    sale_array_vec_binary(fp,"e_vof",_sale->e_vof);
    sale_array_vec_binary(fp,"e_csd",_sale->e_csd);
    sale_array_vec_binary(fp,"e_dam",_sale->e_dam);
    sale_array_vec_binary(fp,"e_tps",_sale->e_tps);
    sale_array_vec_binary(fp,"e_q",_sale->e_q);
    sale_array_vec_binary(fp,"e_tem",_sale->e_tem);
    sale_array_vec_binary(fp,"e_eng",_sale->m_eng);
    sale_array_vec_binary(fp,"m_den",_sale->m_den);
    sale_array_vec_binary(fp,"e_ste2s",_sale->e_ste2s);
    sale_array_vec_binary(fp,"e_sta2s",_sale->e_sta2s);
    sale_array_vec_binary(fp,"e_ste",_sale->e_ste);
    sale_array_vec_binary(fp,"e_strength",_sale->e_strength);
    sale_array_vec_binary(fp,"e_vib",_sale->e_vib);
    sale_array_vec_binary(fp,"e_ejt",_sale->e_ejt);
    sale_array_vec_binary(fp,"e_cvs",_sale->e_cvs);
    sale_array_vec_binary(fp,"m_pty",_sale->m_pty);
#ifdef DEBUG_MODE
    sale_array_vec_binary(fp,"edebug",_sale->e_debug);
    sale_array_vec_binary(fp,"esuf",_sale->e_suf);
    sale_array_vec_binary(fp,"esubvol",_sale->e_subvol);
    sale_array_vec_binary(fp,"egrad",_sale->e_grad);
    sale_array_vec_binary(fp,"divvel",_sale->e_div_vel);
    sale_array_vec_binary(fp,"e_alpha",_sale->e_alpha);
    sale_array_vec_binary(fp,"e_sigma",_sale->e_sigma);
    sale_array_vec_binary(fp,"e_pos",_sale->e_pos);
    sale_array_vec_binary(fp,"e_vol",_sale->e_vol);
    sale_array_vec_binary(fp,"e_vel",_sale->e_vel);
    sale_array_vec_binary(fp,"e_dampref",_sale->e_dampref);
    sale_array_vec_binary(fp,"e_gradvel",_sale->e_grad_vel);
    sale_array_vec_binary(fp,"e_cnd",_sale->e_cnd);
#endif

    vtk_cell_data_trailer(fp);
    vtk_output_coord(fp,_minfo->format,(double *)_sale->v_pos->data,_sale->v_pos->data_size);
    vts_file_trailer(fp);
    fclose(fp);

    // at last write the pvtk file
    if(_ccl->rank == 1)
    {
        sale2d_write_pvtk(_sale,_minfo,_ccl);
        sale2d_write_vtm(_minfo->npgy*_minfo->npgx,_ccl->k_out_step-1,_ccl->prefix);
    }
}


void sale2d_write_pvtp(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl)
{
    int cycles = _ccl->k_out_step - 1;
    char fname[256];
    sprintf(fname,"%s.tracer.%04d.pvtp",_minfo->prefix,cycles);
    FILE * fp = fopen(fname,"w");

    const char header[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    fputs(header, fp);

    // construct information on PPolyData
    fprintf(fp,"<PPolyData GhostLevel=\"%d\">\n",0);

    pvtk_point_data_header(fp);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","matid",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","tag",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","gx",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","gy",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","vel",2ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","den",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","tem",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","pre",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","mpre",1ul);
    fprintf(fp,"       <PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lu\"/>\n","mtem",1ul);
    pvtk_point_data_trailer(fp);

    // add coordinate tag
    fprintf(fp,"    <PPoints>\n"
               "      <PDataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\"/>\n"
               "    </PPoints>\n");

    int n_prcs = _minfo->npgx*_minfo->npgy;
    for(int k=0;k<n_prcs;++k)
    {
        fprintf(fp,"      <Piece Source=\"./vtp/%s.tracer.proc%04d.%04d.vtp\"/>\n",_minfo->prefix,k,cycles);
    }

    fputs("</PPolyData>\n",fp);
    fputs("</VTKFile>",fp);
    fclose(fp);
}

void sale2d_write_vtp(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl)
{
#ifndef UPDATE_TRACERS
    return;
#endif

    if(NULL == _sale->tracer_list) return;
    int num_tracers = _sale->tracer_list->len - _sale->tracer_list->dirty_nodes;
    if(num_tracers < 0) return;


    char output_name[200];
    int myrank = _ccl->rank;
    /*
     * get the start of ext distribution
     *
     */
    proc_info_2d *proc_self = (proc_info_2d *) _sale->v_pos->proc_self;
    sprintf(output_name,"./vtp/%s.tracer.proc%04d.%04d.vtp",_minfo->prefix,myrank,_ccl->k_out_step-1);

    FILE * fp = fopen(output_name,"w");
    if(NULL == fp)
    {
        fprintf(stderr,"cannot create %s\n",output_name);
        exit(0);
    }

    int write_dummpy_point = 0;

    if(num_tracers == 0)
    {
        write_dummpy_point = 1;
        num_tracers = 1;
    }

    // malloc memory for output
    float * tracer_pos = (float *) MV_MALLOC(sizeof(float)*num_tracers*3);
    float * tracer_mat = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_tag = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_gx  = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_gy  = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_vel = (float *) MV_MALLOC(sizeof(float)*num_tracers*2);
    float * tracer_den = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_tem = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_pre = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_mtem = (float *) MV_MALLOC(sizeof(float)*num_tracers);
    float * tracer_mpre = (float *) MV_MALLOC(sizeof(float)*num_tracers);

    int tracer_pos_k = 0;
    field_list_iterator  fit;
    fl_iterator_init(&fit,_sale->tracer_list,field_list_head);
    field_list_node * node = NULL;

    if(write_dummpy_point)
    {
        tracer_pos[tracer_pos_k*3 + X] = 0.0f;
        tracer_pos[tracer_pos_k*3 + Y] = 0.0f;
        tracer_pos[tracer_pos_k*3 + Z] = 0.0f;
        tracer_mat[tracer_pos_k]       = -1.0f;
        tracer_tag[tracer_pos_k]       = -1.0f;
        tracer_gx[tracer_pos_k]        = -1.0f;
        tracer_gy[tracer_pos_k]        = -1.0f;
        tracer_vel[tracer_pos_k*2 + X]        = 0.0f;
        tracer_vel[tracer_pos_k*2 + Y]        = 0.0f;
        tracer_den[tracer_pos_k] = 0.0f;
        tracer_tem[tracer_pos_k] = 0.0f;
        tracer_pre[tracer_pos_k] = 0.0f;
        tracer_mtem[tracer_pos_k] = 0.0f;
        tracer_mpre[tracer_pos_k] = 0.0f;
    }

    while(NULL != (node = fl_iterator_next(&fit)))
    {
        if(node->state != field_list_deleted)
        {
            tracer_pos[tracer_pos_k*3 + X] = (float ) node->data.position[X];
            tracer_pos[tracer_pos_k*3 + Y] = (float ) node->data.position[Y];
            tracer_pos[tracer_pos_k*3 + Z] = 0.0f;
            tracer_mat[tracer_pos_k]       = 1.0f * node->data.matid;
            tracer_tag[tracer_pos_k]       = 1.0f * node->data.tag;
            tracer_gx[tracer_pos_k]       = 1.0f * node->data.gx;
            tracer_gy[tracer_pos_k]       = 1.0f * node->data.gy;
            tracer_vel[tracer_pos_k*2 + X]    = (float ) node->data.velocity[X];
            tracer_vel[tracer_pos_k*2 + Y]    = (float ) node->data.velocity[Y];
            tracer_mtem[tracer_pos_k]     = (float ) node->data.mtem;
            tracer_mpre[tracer_pos_k]     = (float ) node->data.mpre;


            // the data in recorder is updated, if need output
            double *node_eden_ptr, *node_etem_ptr, *node_epre_ptr;
            _sale->foe->get_double(&(node->id),&node_eden_ptr,_sale->e_den);
            _sale->foe->get_double(&(node->id),&node_etem_ptr,_sale->e_tem);
            _sale->foe->get_double(&(node->id),&node_epre_ptr,_sale->e_pre);
            // after fetch the pointer, copy data to tracer
            tracer_den[tracer_pos_k] = (float) node_eden_ptr[0];
            tracer_tem[tracer_pos_k] = (float) node_etem_ptr[0];
            tracer_pre[tracer_pos_k] = (float) node_epre_ptr[0];
            tracer_pos_k++;
        }
    }


    if(num_tracers >= 1)
    {
        vtp_file_header(fp,num_tracers);
        vtk_point_data_header(fp);
        vtk_dataarray_vec_f(fp,"matid","binary",tracer_mat,num_tracers,1);
        vtk_dataarray_vec_f(fp,"tag","binary"  ,tracer_tag,num_tracers,1);
        vtk_dataarray_vec_f(fp,"gx","binary"   ,tracer_gx,num_tracers,1);
        vtk_dataarray_vec_f(fp,"gy","binary"   ,tracer_gy,num_tracers,1);
        vtk_dataarray_vec_f(fp,"vel","binary"  ,tracer_vel,num_tracers,2);
        vtk_dataarray_vec_f(fp,"den","binary"  ,tracer_den,num_tracers,1);
        vtk_dataarray_vec_f(fp,"tem","binary"  ,tracer_tem,num_tracers,1);
        vtk_dataarray_vec_f(fp,"pre","binary"  ,tracer_pre,num_tracers,1);
        vtk_dataarray_vec_f(fp,"mpre","binary"  ,tracer_mpre,num_tracers,1);
        vtk_dataarray_vec_f(fp,"mtem","binary"  ,tracer_mtem,num_tracers,1);

        vtk_point_data_trailer(fp);

        vtk_point_header(fp);
        vtk_dataarray_vec_f(fp,"coordinate","binary",tracer_pos,num_tracers,3);
        vtk_point_trailer(fp);
    } else
    {
        vtp_file_header(fp,1);
        vtk_point_data_header(fp);
        vtk_point_data_trailer(fp);

        vtk_point_header(fp);
        fprintf(fp,"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n"
                   "          0.0 -1.0 0.0\n"
                   "        </DataArray>\n");
        vtk_point_trailer(fp);
    }
    vtp_file_trailer(fp);
    fclose(fp);


    if(_ccl->out_tracer_txt)
    {
        // export tracer info in txt
        sprintf(output_name,"./txt/%s.tracer.proc%04d.%04d.txt",_minfo->prefix,myrank,_ccl->k_out_step-1);
        FILE * txt_fp = fopen(output_name,"w");
        if(NULL == txt_fp)
        {
            fprintf(stderr,"cannot open %s\n", output_name);
            exit(0);
        }
        fprintf(txt_fp,"   %d\n",num_tracers - write_dummpy_point);
        for(tracer_pos_k=0;tracer_pos_k<num_tracers;++tracer_pos_k)
        {
            fprintf(txt_fp,"%10f, %5f, %5f, %5f, %10.5e, %10.5e, %10.5e, %10.5e, %10.5e, %10.5e, %10.5e, %10.5e\n",
                    tracer_tag[tracer_pos_k], tracer_mat[tracer_pos_k],tracer_gx[tracer_pos_k],tracer_gy[tracer_pos_k],
                    tracer_pos[3*tracer_pos_k + X],tracer_pos[3*tracer_pos_k + Y],tracer_vel[2*tracer_pos_k + X],tracer_vel[2*tracer_pos_k + Y],
                    tracer_mpre[tracer_pos_k],tracer_mtem[tracer_pos_k],tracer_pre[tracer_pos_k],tracer_tem[tracer_pos_k]);
        }
        fclose(txt_fp);
    }

    MV_FREE(tracer_pos);
    MV_FREE(tracer_mat);
    MV_FREE(tracer_tag);
    MV_FREE(tracer_gx);
    MV_FREE(tracer_gy);
    MV_FREE(tracer_vel);
    MV_FREE(tracer_den);
    MV_FREE(tracer_tem);
    MV_FREE(tracer_pre);
    MV_FREE(tracer_mtem);
    MV_FREE(tracer_mpre);

    if(_ccl->rank == 1)
    {
        sale2d_write_pvtp(_sale,_minfo,_ccl);
    }
}

