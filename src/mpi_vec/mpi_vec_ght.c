//
// Created by huach on 19/5/2023.
//

#include "mpi_vec.h"

#define LEFT 0
#define RIGHT 1
#define BOTTOM 2
#define TOP 3

mpi_vec_state init_field_ght_table(field_var * fv, field_opt * fo)
{
    proc_info_2d * proc_self = (proc_info_2d *) (fv->proc) + fv->rank;
    proc_info_2d * proc_domain = (proc_info_2d *) (fv->domain);
    if(NULL == fv->data_type) return malloc_error;

    // check this cell is on a boundary
    int invalid_neighbour = 0;
    for(int nei_proc=0;nei_proc<NNP2;++nei_proc)
    {
        int nei_rank = proc_self->neighbour[nei_proc];
        if(nei_rank < 0) invalid_neighbour++;
    }

    if(0 == invalid_neighbour)
    {
        // this progress is not on the boundary, skip
        fv->gt_list_len = 0;
        fv->gt_list = NULL;
        return success;
    }

    unsigned int fv_gt_len = 0;
    // check number of ghost cell
    for(int index=0;index<fv->data_size;++index)
    {

//        if(fv->data_type[index] != receive_section)
//        {
//            // all the ghost cell is in receive section.
//            continue;
//        }

        // loop in the data_type
        fv_id_2d lid,gid;
        fo->index2lid(fv->rank,index,&lid,fv);
        fo->lid2gid(fv->rank,&lid,&gid,fv);

        boundary_type local_bt = undefined_boundary;
        fo->get_boundary(&gid,&local_bt,fv);
        if(local_bt & check_boundary)
        {
            fv_gt_len ++;
            // record the number of boundary cells
        } else
        {
            continue;
        }
    }

    ghost_table * tmp_gt_list = (ghost_table *) malloc(sizeof(ghost_table)*fv_gt_len);
    // loop on the elements and init ghost table

    unsigned int k_gt_table = 0;

    for(int index=0;index<fv->data_size;++index)
    {

//        if(fv->data_type[index] != receive_section)
//        {
//            continue;
//        }

        // loop in the data_type
        fv_id_2d lid,gid;
        fo->index2lid(fv->rank,index,&lid,fv);
        fo->lid2gid(fv->rank,&lid,&gid,fv);

        boundary_type local_bt = undefined_boundary;
        fo->get_boundary(&gid,&local_bt,fv);

        if(undefined_boundary == (local_bt & check_boundary)) continue;


        fv_id_2d ght_a = {.x = lid.x, .y= lid.y},
                ght_b = {.x = lid.x, .y= lid.y},
                ght_c = {.x = lid.x, .y= lid.y};


        if(left_boundary & local_bt)
        {
            ght_a.x = lid.x + 1;
            ght_b.x = lid.x + 2;
            ght_c.x = lid.x + 3;
        }

        if(right_boundary & local_bt)
        {
            ght_a.x = lid.x - 1;
            ght_b.x = lid.x - 2;
            ght_c.x = lid.x - 3;
        }

        if(bottom_boundary & local_bt)
        {
            ght_a.y = lid.y + 1;
            ght_b.y = lid.y + 2;
            ght_c.y = lid.y + 3;
        }

        if(top_boundary & local_bt)
        {
            ght_a.y = lid.y - 1;
            ght_b.y = lid.y - 2;
            ght_c.y = lid.y - 3;
        }

        /*
         * at corners of physical domain
         * some ghost reference cell is still ghost cells
         * just move them into send_section
         */
        if(MIN(gid.x, abs(proc_domain->nx-1-gid.x)) + MIN(gid.y, abs(proc_domain->ny-1-gid.y)) <= 1)
        {
            ghost_pad_move(&ght_b,fv);
            ghost_pad_move(&ght_c,fv);
        }


        // fetch the data char ptr, store it into ghost info
        fo->get_data(&lid,&(tmp_gt_list[k_gt_table].data),fv);
        fo->get_data(&ght_a,&(tmp_gt_list[k_gt_table].ref_a),fv);
        fo->get_data(&ght_b,&(tmp_gt_list[k_gt_table].ref_b),fv);
        fo->get_data(&ght_c,&(tmp_gt_list[k_gt_table].ref_c),fv);
        tmp_gt_list[k_gt_table].tag = local_bt;
        tmp_gt_list[k_gt_table].opt = ght_undefined;
        k_gt_table++;
    }

    fv->gt_list = tmp_gt_list;
    fv->gt_list_len = fv_gt_len;

    // check the number of ghost cells
    if(k_gt_table != fv_gt_len)
    {
        return malloc_error;
    }
    else
    {
        return success;
    }

}

mpi_vec_state set_ghost_pattern(ght_opt * go_list, field_var * fv, ght_var_type var_type)
{
    /*
     * var_type define conversion between go_list and ght_opt in field vars
     * var_type 0 set all ght_outflow
     *          1 just copy from go_lits
     *          2 set for vel/vector
     *          3 set foe ste/tensor
     *  go_list / left, right, top, bottom /
     */
    ght_opt local_go_list[4];

    switch (var_type)
    {
        case ght_var_always_copy:
            // scalar variables, only try copy from internal
            local_go_list[0] = ght_outflow;
            local_go_list[1] = ght_outflow;
            local_go_list[2] = ght_outflow;
            local_go_list[3] = ght_outflow;
            break;
        case ght_normal:
            // just copy go_list
            for(int k=0;k<4;++k) local_go_list[k] = go_list[k];
            break;
        case ght_velocity_inc:
            /*
             * this option is processed in boundary_str2opt,
             * check for ght_mirror, change it to ght_reverse_
             */
            for(int k=0;k<4;++k) local_go_list[k] = go_list[k];
            if(go_list[0] == ght_mirror) local_go_list[0]  = ght_reverse_vel_x;
            if(go_list[3] == ght_outflow) local_go_list[3] = ght_reverse_vel_absy;
            break;
        case ght_velocity:
            for(int k=0;k<4;++k) local_go_list[k] = ght_outflow_vertex;
            if(go_list[0] == ght_mirror) local_go_list[0]  = ght_reverse_vel_x;
            break;
        case ght_stress:
            /*
             * set the ghost of ste,
             * in most time, it should be ght_outflow
             * but with ght_mirror at left boundary it should be set at ght_reverse_ste_1
             */
            local_go_list[0] = ght_copy_ste;
            local_go_list[1] = ght_copy_ste;
            local_go_list[2] = ght_copy_ste;
            local_go_list[3] = ght_copy_ste;
            if(go_list[0] == ght_mirror) local_go_list[0] = ght_reverse_ste_1;
            break;
        case ght_density:
            local_go_list[LEFT] = ght_outflow;
            local_go_list[RIGHT] = ght_outflow;
            local_go_list[BOTTOM] = ght_outflow;
            local_go_list[TOP] = ght_clean;
            break;
        default:
            local_go_list[0] = ght_outflow;
            local_go_list[1] = ght_outflow;
            local_go_list[2] = ght_outflow;
            local_go_list[3] = ght_outflow;
            break;
    }

    if(fv->gt_list_len >0 && fv->gt_list[0].opt == ght_undefined )
    {
        for(int k=0;k<fv->gt_list_len;++k)
        {
            boundary_type current_bt = fv->gt_list[k].tag;
            if(current_bt & left_boundary)
            {
                fv->gt_list[k].opt = local_go_list[0];
            } else if(current_bt & right_boundary)
            {
                fv->gt_list[k].opt = local_go_list[1];
            } else if(current_bt & bottom_boundary)
            {
                fv->gt_list[k].opt = local_go_list[2];
            } else if(current_bt & top_boundary)
            {
                fv->gt_list[k].opt = local_go_list[3];
            } else
            {
                fv->gt_list[k].opt = ght_zero;
            }
        }
        return success;
    }
    else if(fv->gt_list_len > 0 && fv->gt2_list == NULL)
    {
        fv->gt2_list_len = fv->gt_list_len;
        fv->gt2_list = (ghost_table *) malloc(sizeof(ghost_table)*fv->gt_list_len);
        for(int k=0;k<fv->gt_list_len;++k)
        {
            boundary_type current_bt = fv->gt_list[k].tag;
            memcpy(fv->gt2_list+k,fv->gt_list+k, sizeof(ghost_table));
            if(current_bt & left_boundary)
            {
                fv->gt2_list[k].opt = local_go_list[0];
            } else if(current_bt & right_boundary)
            {
                fv->gt2_list[k].opt = local_go_list[1];
            } else if(current_bt & bottom_boundary)
            {
                fv->gt2_list[k].opt = local_go_list[2];
            } else if(current_bt & top_boundary)
            {
                fv->gt2_list[k].opt = local_go_list[3];
            } else
            {
                fv->gt2_list[k].opt = ght_zero;
            }

        }
        return success;
    }
    else
    {
        return duplicate_error;
    }

}

mpi_vec_state ghost_pad_fill(ghost_table * gt, field_var * fv)
{
    /*
     * fill ghost cell use element ref_b and ref c
     *  |<- gt ->|<- ref_a ->|<- ref_b ->|<- ref_c ->
     *  (calculate boundary)  is on the left of gt
     *  (physical boundary) is between ref_a anb ref_b
     */
    int n_double = (int) (fv->bytes_per_unit/ sizeof(double));
    double * _gt_data = (double *)gt->data;
    double * _ra_data = (double *)gt->ref_a;
    double * _rb_data = (double *)gt->ref_b;
    double * _rc_data = (double *)gt->ref_c;

    switch (gt->opt) {
        case ght_outflow:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            break;
        case ght_fixed:
            scaler_move((double *)gt->data , n_double, (double *)gt->ref_c, -1.0);
            scaler_move((double *)gt->ref_a, n_double, (double *)gt->ref_b,-1.0);
            break;
        case ght_zero:
            vec_zero(_gt_data,n_double);
            vec_zero(_ra_data,n_double);
            break;
        case ght_reverse_vel_y:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[1] = -1.0*_rc_data[1];
            _ra_data[1] = -1.0*_rb_data[1];
            break;
        case ght_reverse_vel_x:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[0] = -1.0*_rc_data[0];
            _ra_data[0] = -1.0*_rb_data[0];
            break;
        case ght_reverse_vel_absy:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[1] = fabs(_rc_data[1]);
            _ra_data[1] = fabs(_rb_data[1]);
            break;
        case ght_copy_ste:
//            liner_op(_ra_data,n_double,_rb_data,_rc_data,2.0,-1.0);
//            liner_op(_gt_data,n_double,_ra_data,_rb_data,2.0,-1.0);
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            break;
        case ght_reverse_ste_1:
//            vec_copy(_rc_data,_rb_data,n_double);
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[1] = -1.0*_rc_data[1];
            _ra_data[1] = -1.0*_rb_data[1];
            break;
        case ght_outflow_vertex:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rc_data,_ra_data,n_double);
            vec_copy(_rc_data,_rb_data,n_double);
            break;
        case ght_clean:
            /*
             * this boundary keep density/energy decreased (usually used in top boundary)
             * when materials flow back into domain.
             */
            vec_zero(_ra_data,n_double);
            vec_zero(_rb_data,n_double);
            vec_zero(_gt_data,n_double);
        default:
            // unknown operation !
            return error;
    }

    return success;
}


mpi_vec_state ghost_opt_fill(ghost_table * gt, field_var * fv, ght_opt _opt)
{
    /*
     * fill ghost cell use element ref_b and ref c
     *  |<- gt ->|<- ref_a ->|<- ref_b ->|<- ref_c ->
     *  (calculate boundary)  is on the left of gt
     *  (physical boundary) is between ref_a anb ref_b
     *
     *  this function is very like ghost_pad_fill, if the _opt is not ght_undefined
     *  the pad section will be filled with _opt, otherwise ghost_opt_fill do same thing with ghost_pad_fill
     *
     *  when ghost_pad_fill is used in e_vel, this variables need different fill pattern in different calls,
     *  1) in the start of cycles, e_vel is the increment velocity of vertex.
     *  2) but after the v_vel has been updated, e_vel is used as the velocity of cells
     *
     *  the ght_opt in 'the increment velocity of vertex' is related to the boundary type
     *  while in 'the velocity of cells', ght_opt is usually ght_outflow except for cells close to the axis of symmetry
     *  (r approx -0)
     */
    int n_double = (int) (fv->bytes_per_unit/ sizeof(double));
    double * _gt_data = (double *)gt->data;
    double * _ra_data = (double *)gt->ref_a;
    double * _rb_data = (double *)gt->ref_b;
    double * _rc_data = (double *)gt->ref_c;

    ght_opt case_opt = _opt;
    if(case_opt == ght_undefined)
    {
        case_opt = gt->opt;
    }

    switch (case_opt) {
        case ght_outflow:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            break;
        case ght_fixed:
            scaler_move((double *)gt->data , n_double, (double *)gt->ref_c, -1.0);
            scaler_move((double *)gt->ref_a, n_double, (double *)gt->ref_b,-1.0);
            break;
        case ght_zero:
            vec_zero(_gt_data,n_double);
            vec_zero(_ra_data,n_double);
            break;
        case ght_reverse_vel_y:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[1] = -1.0*_rc_data[1];
            _ra_data[1] = -1.0*_rb_data[1];
            break;
        case ght_reverse_vel_x:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[0] = -1.0*_rc_data[0];
            _ra_data[0] = -1.0*_rb_data[0];
            break;
        case ght_reverse_vel_absy:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[1] = fabs(_rc_data[1]);
            _ra_data[1] = fabs(_rb_data[1]);
            break;
        case ght_copy_ste:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            break;
        case ght_reverse_ste_1:
            vec_copy(_rc_data,_gt_data,n_double);
            vec_copy(_rb_data,_ra_data,n_double);
            _gt_data[1] = -1.0*_rc_data[1];
            _ra_data[1] = -1.0*_rb_data[1];
            break;
        default:
            // unknown operation !
            return error;
    }

    return success;
}


mpi_vec_state init_field_rst_table(field_var * fv, field_opt * fo)
{
    // like ght, the restrict is only applied on the physical boundary
    proc_info_2d * proc_self = (proc_info_2d *) (fv->proc) + fv->rank;
    proc_info_2d * proc_domain = (proc_info_2d *) (fv->domain);
    if(NULL == fv->data_type) return malloc_error;

    // check this cell is on a boundary
    int invalid_neighbour = 0;
    for(int nei_proc=0;nei_proc<NNP2;++nei_proc)
    {
        int nei_rank = proc_self->neighbour[nei_proc];
        if(nei_rank < 0) invalid_neighbour++;
    }

    if(0 == invalid_neighbour)
    {
        // this progress is not on the boundary, skip
        fv->rt_list_len = 0;
        fv->rt_list = NULL;
        return success;
    }

    // check number of ghost cell
    unsigned int fv_rt_len = 0;
    for(int index=0;index<fv->data_size;++index)
    {
        fv_id_2d lid,gid;
        fo->index2lid(fv->rank,index,&lid,fv);
        fo->lid2gid(fv->rank,&lid,&gid,fv);

        boundary_type local_bt = undefined_boundary;
        fo->get_boundary(&gid,&local_bt,fv);

        if(local_bt & check_pad)
        {
            fv_rt_len ++;
        }
    }

    restrict_table * tmp_rt_table = (restrict_table *) malloc(sizeof(restrict_table)*fv_rt_len);
    unsigned int k_rt_table = 0;

    for(int index=0;index<fv->data_size;++index)
    {
        fv_id_2d lid,gid;
        fo->index2lid(fv->rank,index,&lid,fv);
        fo->lid2gid(fv->rank,&lid,&gid,fv);

        boundary_type local_bt = undefined_boundary;
        fo->get_boundary(&gid,&local_bt,fv);

        if(undefined_boundary == (local_bt & check_pad))
        {
            continue;
        }

        fv_id_2d rst_a = {.x = lid.x, .y= lid.y};
        fv_id_2d rst_b = {.x = lid.x, .y= lid.y};

        if(left_pad & local_bt)
        {
            rst_a.x = lid.x + 1;
            rst_b.x = lid.x - 1;
        }

        if(right_pad & local_bt)
        {
            rst_a.x = lid.x - 1;
            rst_b.x = lid.x + 1;
        }

        if(bottom_pad & local_bt)
        {
            rst_a.y = lid.y + 1;
            rst_b.y = lid.y - 1;
        }

        if(top_pad & local_bt)
        {
            rst_a.y = lid.y - 1;
            rst_b.y = lid.y + 1;
        }

        // fetch data char ptr and store it into rt_table
        fo->get_data(&lid,&(tmp_rt_table[k_rt_table].data),fv);
        fo->get_data(&rst_a,&(tmp_rt_table[k_rt_table].ref_a),fv);
        fo->get_data(&rst_b,&(tmp_rt_table[k_rt_table].ref_b),fv);
        tmp_rt_table[k_rt_table].tag = local_bt;
        k_rt_table++;
    }

    fv->rt_list = tmp_rt_table;
    fv->rt_list_len = fv_rt_len;
    if(k_rt_table != fv_rt_len)
    {
        return malloc_error;
    }
    else
    {
        return success;
    }

}

mpi_vec_state set_restrict_pattern(ght_opt * go_list, field_var * fv, ght_var_type var_type)
{
    /*
     * the go_list is ghost type of go_list, those boundaries are usually
     * to restricted some variables on vertex like velocity or displacement
     * restrict don't need information about its neighbour vertexes
     * the ghost is usually used on cells and those ghost cells need info of its neighbour
     */

    ght_opt local_ro_list[4] = {rst_none,rst_none,rst_none,rst_none};

    /*
     * int most time just set restrict on bottom and right boundary,
     * but if bottom or right boundaries are out flow, the restrict will not be applied
     * NOTICE: boundary type stored in go_list is
     *       [0] = left, [1] = right, [2] = bottom, [3] = top
     */

    // right boundary
    switch (go_list[RIGHT])
    {
        case ght_outflow:
            local_ro_list[RIGHT] = rst_none;
            break;
        case ght_freeslip:
        case ght_reverse_vel_x:
            local_ro_list[RIGHT] = rst_x_zero;
            break;
        case ght_fixed:
            local_ro_list[RIGHT] = rst_zero;
            break;
        default:
            local_ro_list[RIGHT] = rst_zero;
    }

    // bottom boundary
    switch (go_list[BOTTOM])
    {
        case ght_outflow:
            local_ro_list[BOTTOM] = rst_none;
            break;
        case ght_freeslip:
        case ght_reverse_vel_y:
            local_ro_list[BOTTOM] = rst_y_zero;
            break;
        case ght_fixed:
            local_ro_list[BOTTOM] = rst_zero;
            break;
        default:
            local_ro_list[BOTTOM] = rst_zero;
    }

    switch (go_list[LEFT])
    {
        case ght_reverse_vel_x:
        case ght_mirror:
            local_ro_list[LEFT] = rst_x_zero;
            break;
        case ght_outflow:
        default:
            local_ro_list[LEFT] = rst_none;
    }

//    switch(go_list[TOP])
//    {
//        case ght_outflow:
//            local_ro_list[TOP] = rst_none;
//            break;
//        case ght_freeslip:
//            local_ro_list[TOP] = rst_y_zero;
//            break;
//        default:
//            local_ro_list[TOP] = rst_none;
//            break;
//    }


    for(int k=0;k<fv->rt_list_len;++k)
    {
        boundary_type current_bt = fv->rt_list[k].tag;

        fv->rt_list[k].opt = rst_none;
        if(current_bt & left_pad)
        {
            fv->rt_list[k].opt = local_ro_list[LEFT];
        }

        if(current_bt & right_pad)
        {
            fv->rt_list[k].opt = local_ro_list[RIGHT];
        }

        if(current_bt & top_pad)
        {
            fv->rt_list[k].opt = local_ro_list[TOP];
        }

        if(current_bt & bottom_pad)
        {
            fv->rt_list[k].opt = local_ro_list[BOTTOM];
        }

        // check the corner boundary
        if((current_bt & bottom_pad) && (current_bt & right_pad))
        {
            if((local_ro_list[BOTTOM] != rst_none) && (local_ro_list[RIGHT] != rst_none))
                fv->rt_list[k].opt = rst_zero;
        }
    }

}

mpi_vec_state apply_restrict_vertex(restrict_table * rt, field_var * fv)
{
    /*
     * different the ghost cells
     * ref_a is the vertexes in send_section
     * ref_b is in the ghost cells
     * data is vertexes between ghost and send cells
     */
    int n_double = (int) (fv->bytes_per_unit/ sizeof(double));
    double * _rt_data = (double *) rt->data;
    double * _ra_data = (double *) rt->ref_a;
    double * _rb_data = (double *) rt->ref_b;

    vec_copy(_ra_data,_rb_data,n_double);
    vec_copy(_ra_data,_rt_data,n_double);

    switch (rt->opt)
    {
        case rst_zero:
            vec_zero(_rt_data,n_double);
            break;
        case rst_x_zero:
            _rt_data[0] = 0.;
            _rt_data[1] = _ra_data[1];
            break;
        case rst_y_zero:
            _rt_data[0] = _ra_data[0];
            _rt_data[1] = 0.;
            break;
        case rst_y_positive:
            _rt_data[Y] = fabs(_rt_data[Y]);
            break;
        case rst_x_positive:
            _rt_data[X] = fabs(_rt_data[X]);
            break;
        case rst_none:
        default:
            break;
    }
}

mpi_vec_state sync_ghost_cells(field_var * fv)
{
    if(0 == fv->sync_gt_state)
    {
        for(int k=0;k<fv->gt_list_len;++k)
        {
            ghost_pad_fill(fv->gt_list + k,fv);
        }
    }
    else
    {
        for(int k=0;k<fv->gt2_list_len;++k)
        {
            ghost_pad_fill(fv->gt2_list + k,fv);
        }
    }

    return success;
}

/*
 * restrict apply the physical boundary on vertex
 * v_dis, v_acc, v_vel
 */

mpi_vec_state sync_restrict_vertexes(field_var * fv)
{
#ifndef RESTRICT_PHY_VERTEX
    return success;
#else
    for(int k=0;k<fv->rt_list_len;++k)
    {
        apply_restrict_vertex(fv->rt_list +k,fv);
    }
    return success;
#endif
}

mpi_vec_state flush_ght(field_var * fv)
{
    /*
     * there are two different list of ghost_table defined in field_var
     * in most time gt2_list is NULL, except for the variable is e_vel
     */
    if(fv->gt_list_len > 0 && fv->gt2_list == NULL)
    {
        sync_ghost_cells(fv);
    }
    else
    {
        sync_ghost_cells(fv);
        fv->sync_gt_state = (fv->sync_gt_state + 1)%2;
    }
    return success;
}
