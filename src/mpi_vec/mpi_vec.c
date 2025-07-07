//
// Created by huacheng on 10/28/22.
//

#include "mpi_vec.h"

int cmp_route_table(const void * x, const void * y)
{
    route_table * a = (route_table *) x;
    route_table * b = (route_table *) y;

    int rtn_value = 0;
    if(a->rank != b->rank)
        rtn_value = a->rank - b->rank;
    else if(a->gx != b->gx)
        rtn_value = a->gx - b->gx;
    else if(a->gy != b->gy)
        rtn_value = a->gy - b->gy;
    else if(a->tag != b->tag)
        rtn_value = a->tag - b->tag;
    else
        rtn_value = a->priority - b->priority;

    return rtn_value;
}

mpi_vec_state gid2lid_2d_default(int rank, void * _gid, void * _lid, field_var * var)
{
    fv_id_2d * gid = (fv_id_2d *) _gid;
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_this = (proc_info_2d *) (var->proc) + rank;

    lid->x = gid->x - proc_this->stx;
    lid->y = gid->y - proc_this->sty;

    if(lid->x < 0 || lid->x >= proc_this->nx || lid->y < 0 || lid->y >= proc_this->ny)
        return id_error;
    else
        return success;
}

mpi_vec_state lid2gid_2d_default(int rank, void* _lid, void* _gid, field_var * var)
{
    fv_id_2d * gid = (fv_id_2d *) _gid;
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_list = (proc_info_2d *) var->proc;
    proc_info_2d * proc_this = proc_list + rank;
    gid->x = lid->x + proc_this->stx;
    gid->y = lid->y + proc_this->sty;

    proc_info_2d * proc_domain = (proc_info_2d *) var->domain;
    if(gid->x < 0 || gid->x >= proc_domain->nx || gid->y < 0 || gid->y >= proc_domain->ny)
        return id_error;
    else
    {
        return success;
    }

    /*
     * in sometime, gid->y + proc_domain->ny * gid->x > max of(int),
     * risk for overflow
     */

}

int lid2pri_2d_default(int rank,void* _lid, field_var * var)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_list = (proc_info_2d *) var->proc;
    proc_info_2d * proc_this = proc_list + rank;

    int rtn_value_x = MIN(lid->x,abs(proc_this->nx-1-lid->x)),
        rtn_value_y = MIN(lid->y,abs(proc_this->ny-1-lid->y));
    return MIN(rtn_value_x,rtn_value_y);
}


mpi_vec_state gid2boundary_2d_default(void * _gid, void* _bound, field_var * var)
{
    proc_info_2d * proc_domain = (proc_info_2d *) var->domain;
    int nx = proc_domain->nx, ny= proc_domain->ny, xpad = proc_domain->xpad, ypad = proc_domain->ypad;
    fv_id_2d * gid = (fv_id_2d *) _gid;
    boundary_type local_bound = undefined_boundary;

    if(gid->x == 0) local_bound |= left_boundary;
    if(gid->y == 0) local_bound |= bottom_boundary;
    if(gid->x == (nx -1)) local_bound |= right_boundary;
    if(gid->y == (ny -1)) local_bound |= top_boundary;

    if(gid->x == xpad) local_bound |= left_pad;
    if(gid->y == ypad) local_bound |= bottom_pad;
    if(gid->x == (nx - 1 - xpad)) local_bound |= right_pad;
    if(gid->y == (ny - 1 - ypad)) local_bound |= top_pad;


    boundary_type * bound = (boundary_type *) _bound;

    // check the validity of gid;
    if(gid->x < 0 || gid->y < 0 || gid->x > (nx -1) || gid->y > (ny-1))
    {
        *bound = undefined_boundary;
        return id_error;
    } else
    {
        *bound = local_bound;
        return success;
    }
}


mpi_vec_state node_type_2d_default_vertex(int rank, void* _lid, node_type * nty, field_var * var)
{
    /*
     * different from elements, variables in vertex are usually determined by variables elements.
     * in most time, synchronized operation is taken on elements, vertex value is updated passively or lazy evaluated.
     * sync_* and complete_* function is not suitable for variables initialized with this function.
     */
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + rank;
    fv_id_2d * lid = (fv_id_2d *) _lid;

    if(lid->x < 0 || lid->x > proc_self->nx || lid->y <0 || lid->y > proc_self->ny)
    {
        *nty = useless_section;
        return id_error;
    }

    int x_distance = MIN(abs(lid->x), abs(proc_self->nx - 1 - lid->x));
    int y_distance = MIN(abs(lid->y), abs(proc_self->ny - 1 - lid->y));

    if(x_distance <= 0 || y_distance <= 0)
    {
        *nty = useless_section;
    } else if(x_distance <= proc_self->xpad +1 || y_distance <= proc_self->ypad + 1)
    {
        *nty = update_section;
    } else
    {
        *nty = internal_section;
    }

    if(x_distance <= 1 || y_distance <= 1)
    {
        *nty |= pad_section;
    }

    /*
     * some variables is shared on boundary, set send/recv pattern if necessary
     * this section may be used in restricted process
     */
    if(x_distance >= proc_self->xpad && y_distance >= proc_self->ypad)
    {
        if((proc_self->nx - 1 - lid->x) == (proc_self->xpad) || (proc_self->ny - 1 - lid->y) == (proc_self->ypad))
        {
            *nty |= receive_section;
        } else if((lid->x == proc_self->xpad) || (lid->y == proc_self->ypad))
        {

            *nty |= send_section;
        }
    }

    /*
     * except for the vertexes on boundary between procs (x/y_distance == x/ypad)
     * neighbours of those vertexes also should be synchronized.
     * x/y_distance == x/pad + 1 (send from internal)
     * x/y_distance == x/pad - 1 (receive from update)
     * required by the calculation of v_laplacian operator which uses 9 points stencil
     */


    if(x_distance == (proc_self->xpad + 1) && y_distance >= (proc_self->ypad + 1))
    {
        *nty |= send_section;
    }

    if(x_distance == (proc_self->xpad - 1) && y_distance >= (proc_self->ypad - 1))
    {
        *nty |= receive_section;
    }

    if(y_distance == (proc_self->ypad + 1) && x_distance >= (proc_self->xpad + 1))
    {
        *nty |= send_section;
    }

    if(y_distance == (proc_self->ypad - 1) && x_distance >= (proc_self->xpad - 1))
    {
        *nty |= receive_section;
    }


    if(__builtin_popcount(*nty & (send_section|receive_section|internal_section)) >= 2)
    {
        /*
         * send/receive/internal should be mutual
         */
        fprintf(stderr,"%s:%d: error send/recv section at lid=(%d,%d,%d,%d,%d)\n",__func__,__LINE__,lid->x,lid->y,x_distance,y_distance,*nty);
        exit(0);
    }

    return success;
}

mpi_vec_state node_type_2d_default_element(int rank, void* _lid, node_type *nty, field_var * var)
{
    /*
     *  node_type of elements imply the information communicate pattern
     *  receive_section,  receive values from neighbor
     *  send_section, send values to neighbors
     *  internal_section, values do not involve any send/receive operation
     */
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + rank;
    fv_id_2d * lid = (fv_id_2d *) _lid;
    int x_distance = MIN(abs(lid->x), abs(proc_self->nx - 1 - lid->x));
    int y_distance = MIN(abs(lid->y), abs(proc_self->ny - 1 - lid->y));

    if(x_distance < proc_self->xpad || y_distance < proc_self->ypad)
    {
        *nty = receive_section;
    }
    else if(x_distance < 2*proc_self->xpad || y_distance < 2*proc_self->ypad)
    {
        *nty = send_section;
    } else
    {
        *nty = internal_section;
    }

    return success;
}

mpi_vec_state node_type_2d_default_bx(int rank, void* _lid, node_type * nty, field_var * var)
{
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + rank;
    fv_id_2d * lid = (fv_id_2d *) _lid;
    int x_distance = MIN(abs(lid->x), abs(proc_self->nx - 1 - lid->x));
    int y_distance = MIN(abs(lid->y), abs(proc_self->ny - 1 - lid->y));

    if(x_distance == 0 || y_distance == 0)
    {
        *nty = useless_section;
    } else if(x_distance <= proc_self->xpad | y_distance < proc_self->ypad)
    {
        *nty = update_section;
    } else
    {
        *nty = internal_section;
    }
    return success;

}

mpi_vec_state node_type_2d_default_by(int rank, void* _lid, node_type * nty, field_var * var)
{
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + rank;
    fv_id_2d * lid = (fv_id_2d *) _lid;
    int x_distance = MIN(abs(lid->x), abs(proc_self->nx - 1 - lid->x));
    int y_distance = MIN(abs(lid->y), abs(proc_self->ny - 1 - lid->y));

    if(x_distance == 0 || y_distance == 0)
    {
        *nty = useless_section;
    } else if(x_distance < proc_self->xpad | y_distance <= proc_self->ypad)
    {
        *nty = update_section;
    } else
    {
        *nty = internal_section;
    }
    return success;
}


int gid2tag_2d_default(int rank,void* _gid,field_var * fv)
{
    return fv->dup_attr;
}


mpi_vec_state init_field_opt_2d_default(field_opt * fo)
{
    fo->gid2lid = gid2lid_2d_default;
    fo->lid2gid = lid2gid_2d_default;
    fo->lid2pri = lid2pri_2d_default;
    fo->gid2tag = gid2tag_2d_default;
    fo->lid2nty = node_type_2d_default_element;
    fo->index2lid = index2lid_2d_default;
    fo->lid2index = lid2index_2d_default;
    fo->get_data = lid2data_2d_default;
    fo->get_double = lid2double_2d_default;
    fo->get_data2 = lid2data_nocheck_2d_default;
    fo->get_nty = lid2node_type_default;
    fo->get_nty2 = lid2node_type_nocheck_default;
    fo->get_boundary = gid2boundary_2d_default;

    if(NULL==fo->get_data|| NULL==fo->get_data2 || NULL==fo->get_double
        ||NULL==fo->gid2lid || NULL==fo->lid2gid || NULL==fo->lid2pri
        || NULL==fo->gid2tag || NULL==fo->lid2nty || NULL==fo->index2lid
        || NULL==fo->lid2index || NULL==fo->get_boundary)
        return nullptr_error;
    else
        return success;
}

mpi_vec_state init_field_route_table(field_var * fv, field_opt * fo)
{
    proc_info_2d * proc_self = (proc_info_2d *) (fv->proc) + fv->rank;
    fv->data_type = (node_type *) malloc(sizeof(node_type)*fv->data_size);
    if(NULL == fv->data_type) return malloc_error;

    unsigned int recv_route_size = 0;
    unsigned int send_route_size = 0;
    for(int index=0;index<fv->data_size;++index)
    {
        fv_id_2d lid,gid;
        fo->index2lid(fv->rank,index,&lid,fv);
        fo->lid2gid(fv->rank,&lid,&gid,fv);
        fo->lid2nty(fv->rank,&lid,fv->data_type+index,fv);

        if(fv->data_type[index] & internal_section || fv->data_type[index] & useless_section)
        {
            // internal section data is never send/recv
            // useless section will be skipped in communication
            continue;
        }

        for(int nei_proc=0;nei_proc<NNP2;++nei_proc)
        {
            int nei_rank = proc_self->neighbour[nei_proc];
            if(nei_rank == fv->rank || nei_rank < 0)
            {
                // skip the rank same with local processor
                continue;
            }

            fv_id_2d nei_lid;
            node_type nei_type;
            mpi_vec_state nei_id_state = fo->gid2lid(nei_rank, &gid, &nei_lid, fv);
            mpi_vec_state nei_nt_state = fo->lid2nty(nei_rank, &nei_lid, &nei_type, fv);
            if(success != nei_id_state || success != nei_nt_state)
            {
                // this data is not in neighbour
                continue;
            }

            // first scan just get the number of send_section/recv route map
            // send and recv operation must be paired
            if(nei_type & receive_section && fv->data_type[index] & send_section)
            {
                send_route_size ++;
            }
            else if(nei_type & send_section && fv->data_type[index] & receive_section)
            {
                recv_route_size ++;
            }

        }
    }

    int send_index = 0;
    int recv_index = 0;
    fv->send_rtb = (route_table *) malloc(sizeof(route_table)*(send_route_size));
    fv->recv_rtb = (route_table *) malloc(sizeof(route_table)*(recv_route_size));
    fv->send_rtb_size = send_route_size;
    fv->recv_rtb_size = recv_route_size;

    // loop same as the previous, after route map have been malloc
    for(int index=0;index<fv->data_size;++index)
    {
        fv_id_2d lid,gid;
        fo->index2lid(fv->rank,index,&lid,fv);
        fo->lid2gid(fv->rank,&lid,&gid,fv);

        if(fv->data_type[index] & internal_section || fv->data_type[index] & useless_section)
        {
            continue;
        }

        for(int nei_proc=0;nei_proc<NNP2;++nei_proc)
        {
            int nei_rank = proc_self->neighbour[nei_proc];
            if(nei_rank == fv->rank || nei_rank < 0)
            {
                continue;
            }

            fv_id_2d nei_lid;
            node_type nei_type;
            mpi_vec_state nei_id_state = fo->gid2lid(nei_rank, &gid, &nei_lid, fv);
            mpi_vec_state nei_nt_state = fo->lid2nty(nei_rank, &nei_lid, &nei_type, fv);
            if(success != nei_id_state || success != nei_nt_state)
            {
                continue;
            }


            if(nei_type & receive_section && fv->data_type[index] & send_section)
            {
                fv->send_rtb[send_index].rank = nei_rank;
                fv->send_rtb[send_index].priority = fo->lid2pri(fv->rank,&lid,fv);
                fv->send_rtb[send_index].tag = fo->gid2tag(fv->rank,&gid,fv);
                fo->get_data(&lid, &(fv->send_rtb[send_index].data), fv);
                fv->send_rtb[send_index].gx = gid.x;
                fv->send_rtb[send_index].gy = gid.y;

                send_index++;
            }
            else if(nei_type & send_section && fv->data_type[index] & receive_section)
            {
                fv->recv_rtb[recv_index].rank = nei_rank;
                fv->recv_rtb[recv_index].priority = fo->lid2pri(fv->rank,&lid,fv);
                fv->recv_rtb[recv_index].tag = fo->gid2tag(nei_rank,&gid,fv);
                fo->get_data(&lid, &(fv->recv_rtb[recv_index].data), fv);
                fv->recv_rtb[recv_index].gx = gid.x;
                fv->recv_rtb[recv_index].gy = gid.y;
                recv_index ++;
            }
        }
    }

    compress_route_table(fv->send_rtb,fv->send_rtb_size,&(fv->send_list),&(fv->send_list_len));
    compress_route_table(fv->recv_rtb,fv->recv_rtb_size,&(fv->recv_list),&(fv->recv_list_len));
}



void ghost_pad_move(fv_id_2d * lid,field_var * fv)
{
    proc_info_2d * proc_self = (proc_info_2d *) (fv->proc) + fv->rank;

    if(lid->x < proc_self->xpad) lid->x = proc_self->xpad;
    if(lid->x > proc_self->nx - 1 - proc_self->xpad) lid->x = proc_self->nx - 1 - proc_self->xpad;

    if(lid->y < proc_self->ypad) lid->y = proc_self->ypad;
    if(lid->y > proc_self->ny - 1 - proc_self->ypad) lid->y = proc_self->ny - 1 - proc_self->ypad;
}



mpi_vec_state compress_route_table(route_table * rtb, unsigned int rtb_size, route_table_continuous ** c_rtb, unsigned int * c_rtb_len)
{
    if(NULL == rtb || rtb_size == 0)
    {
        *c_rtb = NULL;
        *c_rtb_len = 0;
        return nullptr_error;
    }

    qsort(rtb,rtb_size, sizeof(route_table),cmp_route_table);

    unsigned int n_src = 1;
    int recent_src = rtb[0].rank;

    for(unsigned int k =0;k < rtb_size; ++k)
    {
        if(recent_src != rtb[k].rank)
        {
            n_src++;
            recent_src = rtb[k].rank;
        }
    }

    *c_rtb_len = n_src;
    *c_rtb = (route_table_continuous *) malloc(sizeof(route_table_continuous)*n_src);
    if(NULL == *c_rtb)
    {
        *c_rtb_len = 0;
        return malloc_error;
    }

    recent_src = rtb[0].rank;
    unsigned int c_rtb_index = 0;
    (*c_rtb)[0].head = 0;
    (*c_rtb)[0].len  = 1;
    for(int k =1;k < rtb_size; ++k)
    {
        if(recent_src != rtb[k].rank)
        {
            c_rtb_index ++;
            recent_src = rtb[k].rank;
            (*c_rtb)[c_rtb_index].head = k;
            (*c_rtb)[c_rtb_index].len  = 1;
        } else
        {
            (*c_rtb)[c_rtb_index].len += 1;
        }
    }

    if(c_rtb_index == (n_src-1))
        return success;
    else
        return error;
}

mpi_vec_state init_field_buffer(field_var * fv)
{
    fv->recv_buffer_size = fv->recv_rtb_size;
    fv->send_buffer_size = fv->send_rtb_size;
    fv->recv_buffer = (char *)malloc(sizeof(char)*fv->bytes_per_unit*fv->recv_rtb_size);
    fv->send_buffer = (char *)malloc(sizeof(char)*fv->bytes_per_unit*fv->send_rtb_size);
    if(NULL==fv->recv_buffer || NULL==fv->send_buffer)
        return malloc_error;

    fv->num_recv_req = fv->recv_list_len;
    fv->recv_req = (MPI_Request *) malloc(sizeof(MPI_Request)*fv->num_recv_req);
    fv->num_send_req = fv->send_list_len;
    fv->send_req = (MPI_Request *) malloc(sizeof(MPI_Request)*fv->num_send_req);
    if(NULL==fv->recv_req || NULL==fv->send_req)
        return malloc_error;

    for(int k=0;k<fv->num_recv_req;k++)
    {
        unsigned int local_buffer_offset = fv->recv_list[k].head*fv->bytes_per_unit;
        unsigned int local_buffer_size = fv->recv_list[k].len*fv->bytes_per_unit;
        int remote_src_rank = (fv->recv_rtb + fv->recv_list[k].head)->rank;
        int remote_src_tag  = (fv->recv_rtb + fv->recv_list[k].head)->tag;
        MPI_Recv_init(fv->recv_buffer+local_buffer_offset,(int)local_buffer_size,MPI_CHAR,
                      remote_src_rank,remote_src_tag,
                      fv->comm,fv->recv_req+k);
    }

    for(int k=0;k<fv->num_send_req;k++)
    {
        unsigned int local_buffer_offset = fv->send_list[k].head*fv->bytes_per_unit;
        unsigned int local_buffer_size = fv->send_list[k].len*fv->bytes_per_unit;
        int remote_dst_rank = (fv->send_rtb + fv->send_list[k].head)->rank;
        int remote_dst_tag = (fv->send_rtb + fv->send_list[k].head)->tag;
        MPI_Send_init(fv->send_buffer+local_buffer_offset,(int)local_buffer_size,MPI_CHAR,
                      remote_dst_rank,remote_dst_tag,
                      fv->comm,fv->send_req+k);
    }

    return success;
}

mpi_vec_state sync_start(field_var * fv)
{
    // pack the data into buffer
#ifdef COMM_OFF
    return success;
#endif
    for(unsigned int k=0;k<fv->send_rtb_size;++k)
    {
        memcpy(fv->send_buffer+k*fv->bytes_per_unit,fv->send_rtb[k].data,fv->bytes_per_unit);
    }

#ifdef DEBUG_RANK
    MPI_Barrier(fv->comm);
#endif

    MPI_Startall((int)fv->num_send_req,fv->send_req);
    MPI_Startall((int)fv->num_recv_req,fv->recv_req);
    return success;
}

mpi_vec_state sync_complete(field_var * fv)
{
#ifdef COMM_OFF
    return success;
#endif

    MPI_Waitall((int)fv->num_send_req,fv->send_req,MPI_STATUSES_IGNORE);
    MPI_Waitall((int)fv->num_recv_req,fv->recv_req,MPI_STATUSES_IGNORE);
    // unpack data from buffer to data array
    for(unsigned int k=0;k<fv->recv_rtb_size;++k)
    {
        memcpy(fv->recv_rtb[k].data,fv->recv_buffer+k*fv->bytes_per_unit,fv->bytes_per_unit);
    }

    flush_ght(fv);
#ifdef DEBUG_RANK
    MPI_Barrier(fv->comm);
#endif
    return success;
}


mpi_vec_state init_field_proc_2d_default(field_2d *fi, field_var *fv)
{
    /*
     * this initialize function is designed for elements,
     * adjacent process share 2*xpad elements (padding = element_pad)
     */
    static MPI_Comm _comm = MPI_COMM_NULL;
    static int _comm_npgx = 0;
    static int _comm_npgy = 0;

    if(_comm != MPI_COMM_NULL && _comm_npgx == fi->npgx && _comm_npgy == fi->npgy)
    {
        MPI_Comm_dup(_comm,&fv->comm);
    } else
    {
        MPI_Cart_create(MPI_COMM_WORLD,2,
                        (int[]){fi->npgx,fi->npgy},(int[]){0,0},
                        0,&fv->comm);
        _comm = fv->comm;
        _comm_npgx = fi->npgx;
        _comm_npgy = fi->npgy;
    }

    MPI_Comm_rank(fv->comm,&fv->rank);
    fv->proc_num = fi->npgx * fi->npgy;
    fv->proc = malloc(sizeof(proc_info_2d)*(fv->proc_num+1));
    if(fv->proc == NULL) return malloc_error;
    fv->domain = (proc_info_2d *)(fv->proc) + fv->proc_num;

    proc_info_2d * proc_self = (proc_info_2d *)(fv->proc) + fv->rank;
    proc_info_2d * proc_main = (proc_info_2d *)(fv->domain);

    fv->proc_self = proc_self;

    switch (fi->padding) {
        case element_pad:
            proc_self->stx = fi->nx * (fv->rank/fi->npgy);
            proc_self->sty = fi->ny * (fv->rank%fi->npgy);
            proc_self->nx = fi->nx + 2*fi->xpad;
            proc_self->ny = fi->ny + 2*fi->ypad;
            proc_self->xpad = fi->xpad;
            proc_self->ypad = fi->ypad;

            proc_main->stx = 0;
            proc_main->sty = 0;
            proc_main->nx = fi->nx * fi->npgx + 2*fi->xpad;
            proc_main->ny = fi->ny * fi->npgy + 2*fi->ypad;
            proc_main->xpad = fi->xpad;
            proc_main->ypad = fi->ypad;
            break;
        case vertex_pad:
            proc_self->stx = (fi->nx-1) * (fv->rank/fi->npgy);
            proc_self->sty = (fi->ny-1) * (fv->rank%fi->npgy);
            proc_self->nx = (fi->nx) + 2*fi->xpad;
            proc_self->ny = (fi->ny) + 2*fi->ypad;
            proc_self->xpad = fi->xpad;
            proc_self->ypad = fi->ypad;

            proc_main->stx = 0;
            proc_main->sty = 0;
            proc_main->nx = (fi->nx-1) * fi->npgx + 2*fi->xpad + 1;
            proc_main->ny = (fi->ny-1) * fi->npgy + 2*fi->ypad + 1;
            proc_main->xpad = fi->xpad;
            proc_main->ypad = fi->ypad;
            break;
        default:
            proc_self->stx = fi->nx * (fv->rank/fi->npgy);
            proc_self->sty = fi->ny * (fv->rank%fi->npgy);
            proc_self->nx = fi->nx + 2*fi->xpad;
            proc_self->ny = fi->ny + 2*fi->ypad;
            proc_self->xpad = fi->xpad;
            proc_self->ypad = fi->ypad;

            proc_main->stx = 0;
            proc_main->sty = 0;
            proc_main->nx = fi->nx * fi->npgx + 2*fi->xpad;
            proc_main->ny = fi->ny * fi->npgy + 2*fi->ypad;
            proc_main->xpad = fi->xpad;
            proc_main->ypad = fi->ypad;
            break;
    }

    for(int k=0;k<NNP2;++k)
    {
        int sx = k/3 - 1,sy = k%3 - 1;
        int neix = sx + fv->rank/fi->npgy;
        int neiy = sy + fv->rank%fi->npgy;
        if(neix<0||neix>=fi->npgx||neiy<0||neiy>=fi->npgy)
            proc_self->neighbour[k] = -1;
        else
            proc_self->neighbour[k] = neiy + fi->npgy*neix;
    }

    proc_info_2d * proc_self_copy = (proc_info_2d *) malloc(sizeof(proc_info_2d));
    memcpy(proc_self_copy,proc_self, sizeof(proc_info_2d)/sizeof(char));
    MPI_Allgather(proc_self_copy, sizeof(proc_info_2d)/sizeof(char),MPI_CHAR,
                  fv->proc,sizeof(proc_info_2d)/sizeof(char),MPI_CHAR,
                  fv->comm);
    free(proc_self_copy);

    // set space of data
    fv->bytes_per_unit = fi->bytes_per_unit;
    fv->data_size = proc_self->nx * proc_self->ny;
    fv->data = (char *) malloc(sizeof(char)*fi->bytes_per_unit*fv->data_size);
    if(NULL==fv->data) return malloc_error;
    return success;
}


mpi_vec_state set_null_ght(field_var * fv)
{
    fv->gt_list = NULL;
    fv->gt_list_len = 0;

    fv->rt_list = NULL;
    fv->rt_list_len = 0;

    fv->gt2_list = NULL;
    fv->gt2_list_len = 0;

    fv->sync_gt_state = 0;
    return success;
}

mpi_vec_state init_field_2d(field_2d *fi, field_var *fv, field_opt *fo)
{
    mpi_vec_state rtn_value = success;
    rtn_value &= init_field_proc_2d_default(fi,fv);
    if(NULL == fo->gid2lid)
    {
        rtn_value &= init_field_opt_2d_default(fo);
    }

    fv->dup_attr = 0;
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN,
                           &fv->dup_key, NULL);
    int *dup_attr_ptr = (int *) malloc(sizeof(int));
    *dup_attr_ptr = fv->dup_attr;
    MPI_Comm_set_attr(fv->comm,fv->dup_key,dup_attr_ptr);

    rtn_value &= init_field_route_table(fv,fo);
    rtn_value &= init_field_buffer(fv);
    rtn_value &= set_null_ght(fv);

    return rtn_value;
}


mpi_vec_state clean_field(field_var * fv)
{
    free(fv->data);
    free(fv->data_type);
    free(fv->send_rtb);
    free(fv->recv_rtb);
    free(fv->send_list);
    free(fv->recv_list);
    free(fv->send_buffer);
    free(fv->recv_buffer);

    if(fv->dup_attr == 0)
    {
        int * dup_attr_ptr; int flag;
        MPI_Comm_get_attr(fv->comm,fv->dup_key,&dup_attr_ptr,&flag);
        free(dup_attr_ptr);
        free(fv->proc);
        MPI_Comm_free(&fv->comm);
    }
}

mpi_vec_state duplicate_field_2d(field_var *fv_src, field_var * fv_dst, field_opt *fo)
{
    fv_dst->dup_key = fv_src->dup_key;
    fv_dst->comm = fv_src->comm;
    int *dup_attr_ptr;
    int dup_flag;
    MPI_Comm_get_attr(fv_src->comm,fv_src->dup_key,&dup_attr_ptr,&dup_flag);
    if(!dup_flag) return duplicate_error;
    *dup_attr_ptr += 1;
    fv_dst->dup_attr = *dup_attr_ptr;

    fv_dst->proc = fv_src->proc;
    fv_dst->domain = fv_src->domain;
    fv_dst->proc_num = fv_src->proc_num;
    fv_dst->proc_self = fv_src->proc_self;
    fv_dst->rank = fv_src->rank;

    fv_dst->bytes_per_unit = fv_src->bytes_per_unit;
    fv_dst->data_size = fv_src->data_size;
    fv_dst->data = (char *) malloc(sizeof(char)*fv_dst->bytes_per_unit*fv_dst->data_size);
    if(NULL==fv_dst->data) return malloc_error;

    mpi_vec_state rtn_value = success;
    rtn_value &= init_field_route_table(fv_dst,fo);
    rtn_value &= init_field_buffer(fv_dst);
    rtn_value &= set_null_ght(fv_dst);
    if(rtn_value != success)
        return duplicate_error;
    else
        return success;
}

mpi_vec_state duplicate_field_2d_n(field_var *fv_src, field_var * fv_dst, field_opt *fo, int _n)
{
    fv_dst->dup_key = fv_src->dup_key;
    fv_dst->comm = fv_src->comm;
    int *dup_attr_ptr;
    int dup_flag;
    MPI_Comm_get_attr(fv_src->comm,fv_src->dup_key,&dup_attr_ptr,&dup_flag);
    if(!dup_flag) return duplicate_error;
    *dup_attr_ptr += 1;
    fv_dst->dup_attr = *dup_attr_ptr;
    fv_dst->proc = fv_src->proc;
    fv_dst->domain = fv_src->domain;
    fv_dst->proc_num = fv_src->proc_num;
    fv_dst->proc_self = fv_src->proc_self;

    fv_dst->rank = fv_src->rank;
    if(_n < 0)
        fv_dst->bytes_per_unit = abs(_n) * sizeof(double);
    else
        fv_dst->bytes_per_unit = _n * fv_src->bytes_per_unit;

    fv_dst->data_size = fv_src->data_size;
    fv_dst->data = (char *) malloc(sizeof(char)*fv_dst->bytes_per_unit*fv_dst->data_size);
    if(NULL==fv_dst->data) return malloc_error;

    mpi_vec_state rtn_value = success;
    rtn_value &= init_field_route_table(fv_dst,fo);
    rtn_value &= init_field_buffer(fv_dst);
    rtn_value &= set_null_ght(fv_dst);
    if(rtn_value != success)
        return duplicate_error;
    else
        return success;
}


mpi_vec_state lid2index_2d_default(int rank, void * _lid, int * _index, field_var * fv)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(fv->proc) + rank;
    if(lid->x < 0 || lid->x >= proc_self->nx || lid->y < 0 || lid->y >= proc_self->ny)
        return id_error;
    *_index = LID2INDEX(lid->x,lid->y,proc_self->nx,proc_self->ny);
    return success;
}

mpi_vec_state index2lid_2d_default(int rank, int _index, void * _lid, field_var *fv)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(fv->proc) + rank;

#ifdef COLOUM_ARRAY
    lid->x = _index%proc_self->nx;
    lid->y = _index/proc_self->nx;
#else
    lid->x = _index/proc_self->ny;
    lid->y = _index%proc_self->ny;
#endif
    if(_index >=0 && _index<fv->data_size)
        return success;
    else
        return id_error;
}

mpi_vec_state lid2data_2d_default(void * _lid, char ** dest, field_var * var)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + var->rank;

    mpi_vec_state rtn_value = success;
    if(lid->x < 0 || lid->y < 0 || lid->x>= proc_self->nx || lid->y >= proc_self->ny)
    {
        rtn_value = id_error;
        *dest = NULL;
    } else
    {
        rtn_value = success;
        *dest = (char *)(var->data + LID2INDEX(lid->x,lid->y,proc_self->nx,proc_self->ny)*var->bytes_per_unit);
    }
    return rtn_value;
}

mpi_vec_state lid2double_2d_default(void * _lid, double ** dest, field_var * var)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + var->rank;
    char * dest_cptr = NULL;

    mpi_vec_state rtn_value = success;
    if(lid->x < 0 || lid->y < 0 || lid->x>= proc_self->nx || lid->y >= proc_self->ny)
    {
        rtn_value = id_error;
        dest_cptr = NULL;
    } else
    {
        rtn_value = success;
        dest_cptr = (char *)(var->data + LID2INDEX(lid->x,lid->y,proc_self->nx,proc_self->ny)*var->bytes_per_unit);
    }

    *dest = (double *) dest_cptr;
    return rtn_value;
}

mpi_vec_state lid2data_nocheck_2d_default(void * _lid, char ** dest, field_var * var)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(var->proc) + var->rank;
    *dest = (char *)(var->data + LID2INDEX(lid->x,lid->y,proc_self->nx,proc_self->ny)*var->bytes_per_unit);
    return success;
}

mpi_vec_state lid2node_type_nocheck_default(void * _lid, node_type * ntype, field_var * fv)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(fv->proc) + fv->rank;
    *ntype = *(fv->data_type + LID2INDEX(lid->x,lid->y,proc_self->nx,proc_self->ny));
    return success;
}

mpi_vec_state lid2node_type_default(void * _lid, node_type * ntype, field_var * fv)
{
    fv_id_2d * lid = (fv_id_2d *) _lid;
    proc_info_2d * proc_self = (proc_info_2d *)(fv->proc) + fv->rank;
    if(lid->x < 0 || lid->y < 0 || lid->x>= proc_self->nx || lid->y >= proc_self->ny)
    {
        *ntype = useless_section;
        return id_error;
    }
    *ntype = *(fv->data_type + LID2INDEX(lid->x,lid->y,proc_self->nx,proc_self->ny));
    return success;
}