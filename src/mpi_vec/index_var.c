//
// Created by huach on 11/7/2023.
//

#include "index_var.h"

mpi_vec_state init_index_2d(field_2d * fi, field_var * fv)
{
    /*
     * initialize the variables only related with lid.x or lid.y
     * the e_pos and v_pos can be represented by this type variables
     *
     * there are no complex info exchange between progress.
     * it is unnecessary to initialize route table or proc info.
     * in the future version index_var should be implemented by another struct
     */

    MPI_Comm _comm = MPI_COMM_WORLD;

    int rank;
    MPI_Comm_rank(_comm,&rank);
    int x_rank = rank/fi->npgy;
    int y_rank = rank%fi->npgy;

    if((fi->nx > 0) && (fi->ny <= 0))
    {
        MPI_Comm_split(_comm,x_rank,y_rank,&(fv->comm));
        fv->bytes_per_unit = fi->bytes_per_unit;
        fv->data_size = fi->nx + 2*fi->xpad;
        fv->padding = x_pad;
    }
    else if((fi->ny > 0) && (fi->nx <= 0))
    {
        MPI_Comm_split(_comm,y_rank,x_rank,&(fv->comm));
        fv->bytes_per_unit = fi->bytes_per_unit;
        fv->data_size = fi->ny + 2*fi->ypad;
        fv->padding = y_pad;
    }
    else
    {
        fprintf(stderr,"%s[%d]: nx=%d and ny=%d\n",__func__ ,__LINE__,fi->nx,fi->ny);
        exit(0);
    }

    // at last only the data is malloced, no buffer needed
    fv->rank = rank;
    fv->data = malloc(fv->data_size* sizeof(char)*fv->bytes_per_unit);

    // malloc request for reduce opt, only use the recv_req
    fv->recv_req = malloc(sizeof(MPI_Request)*2);
    fv->recv_rtb_size = 2;
    fv->send_req = NULL;
    fv->send_rtb_size = 0;

    // set proc as NULL to avoid undefined pointers
    fv->proc = fv->domain = fv->proc_self = NULL;

    if(NULL == fv->data) return malloc_error;
    return success;
}

mpi_vec_state sync_index_start(field_var * fv,MPI_Op sOP)
{
#ifdef COMM_OFF
    return success;
#endif

    MPI_Iallreduce(MPI_IN_PLACE,fv->data,fv->data_size,MPI_DOUBLE,sOP,fv->comm,fv->recv_req);

#ifdef DEBUG_RANK
    MPI_Barrier(fv->comm);
#endif
    return success;
}

mpi_vec_state sync_index_complete(field_var * fv)
{
#ifdef COMM_OFF
    return success;
#endif
    MPI_Wait(fv->recv_req,MPI_STATUSES_IGNORE);
    return success;
}

mpi_vec_state lid2double_index2d_default(void * _lid, double ** dest, field_var * var)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    char * dest_cptr = NULL;
    // the proc is undefined in index_var,
    // check the validity of lid using padding and data_size

    mpi_vec_state rtn_value = success;

    if(lid->x < 0 || lid->y < 0 || (x_pad == var->padding && lid->x > var->data_size) || (y_pad == var->padding && lid->y > var->data_size))
    {
        rtn_value = id_error;
        dest_cptr = NULL;
    } else
    {
        rtn_value = success;

        if(var->padding == x_pad)
        {
            dest_cptr = (char *)(var->data + lid->x*var->bytes_per_unit);
        }
        else if(var->padding == y_pad)
        {
            dest_cptr = (char *)(var->data + lid->y*var->bytes_per_unit);
        }
        else
        {
            rtn_value = error;
            dest_cptr = NULL;
        }
    }

    *dest = (double *) dest_cptr;
    return rtn_value;
}

mpi_vec_state lid2data_index2d_default(void * _lid, char ** dest, field_var * var)
{
    /*
     * same as lid2double_index2d_default, just change the type of dest pointer
     */
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + var->rank;
    char * dest_ptr = NULL;

    mpi_vec_state rtn_value = success;
    if(lid->x < 0 || lid->y < 0 || (x_pad == var->padding && lid->x > var->data_size) || (y_pad == var->padding && lid->y > var->data_size))
    {
        rtn_value = id_error;
        dest_ptr = NULL;
    } else
    {
        rtn_value = success;

        if(var->padding == x_pad)
        {
            dest_ptr = (char *)(var->data + lid->x * var->bytes_per_unit);
        }
        else if(var->padding == y_pad)
        {
            dest_ptr = (char *)(var->data + lid->y * var->bytes_per_unit);
        }
        else
        {
            rtn_value = error;
            dest_ptr = NULL;
        }
    }

    *dest =  dest_ptr;
    return rtn_value;
}

mpi_vec_state init_index_opt_2d_default(field_opt * fo)
{
    fo->gid2lid    = NULL;
    fo->lid2gid    = NULL;
    fo->lid2pri    = NULL;
    fo->gid2tag    = NULL;
    fo->lid2nty    = NULL;
    fo->index2lid  = NULL;
    fo->lid2index  = NULL;
    fo->get_data   = lid2data_index2d_default;
    fo->get_double = lid2double_index2d_default;
    fo->get_data2  = NULL;
    fo->get_nty    = NULL;
    fo->get_nty2   = NULL;
    fo->get_boundary = NULL;

    return success;
}