//
// Created by huach on 25/5/2023.
//

#include "tracer.h"


int fl_node_nei(fv_id_2d * _id, proc_info_2d * _proc)
{

    int nei_x = GET_OFFSET4(_id->x,0,_proc->xpad, _proc->nx - 1 - _proc->xpad,_proc->nx - 1);
    int nei_y = GET_OFFSET4(_id->y,0,_proc->ypad, _proc->ny - 1 - _proc->ypad,_proc->ny - 1);
//    if (_id->x < 0 || _id->x > _proc->nx - 1)
//    {
//        nei_x = -1;
//    }
//    else if (_id->x < _proc->xpad)
//    {
//        nei_x = 0;
//    }
//    else if (_id->y > _proc->nx - _proc->xpad - 1)
//    {
//        nei_x = 2;
//    }
//    else
//    {
//        nei_x = 1;
//    }

//    if (_id->y < 0 || _id->y > _proc->ny - 1)
//    {
//        nei_y = -1;
//    }
//    else if(_id->y < _proc->ypad)
//    {
//        nei_y = 0;
//    }
//    else if(_id->y > _proc->ny - _proc->ypad - 1)
//    {
//        nei_y = 1
//    }
//    else
//    {
//        nei_y = 2;
//    }

    return (nei_y >= 0 && nei_x >= 0) ? nei_y + 3 * nei_x : -1;

}


field_list_state fl_update_node(field_list * self, field_list_node * node,double dt)
{
    if(node->state == field_list_deleted)
    {
        return node->state;
    }

    if(node->state == field_list_new)
    {
        // this tracer is received from other process, before calculate its displacement
        // the id must be reset
        node->id.x = bisearch(node->data.position[X], self->vx, self->nvx);
        node->id.y = bisearch(node->data.position[Y], self->vy, self->nvy);
        node->state = field_list_used;
    }

    proc_info_2d * proc_self = (proc_info_2d *) self->proc_self;
    int k = node->id.x, j = node->id.y;

    fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
    fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};

    double *vdis_ptr[4], *vpos_ptr[4], *vvel_ptr[4], *vbf_ptr[4], *vacc_ptr[4];
    field_opt * fo = self->fov;
    for(int vid=0;vid<4;++vid)
    {
        fo->get_double(vlid[vid],&vdis_ptr[vid],self->v_dis);
        fo->get_double(vlid[vid],&vpos_ptr[vid],self->v_pos);
        fo->get_double(vlid[vid],&vvel_ptr[vid],self->v_vel);
        fo->get_double(vlid[vid],&vbf_ptr[vid],self->v_bf);
        fo->get_double(vlid[vid],&vacc_ptr[vid],self->v_acc);
    }
    double *evof_ptr;
    self->foe->get_double(&lb,&evof_ptr,self->e_vof);

    // calculate the relative position in cells
    // relative coordinate (xi,mu) should in [-1,1]x[-1,1]
    double relative_pos[2];
    liner_op(relative_pos,2,node->data.position,vpos_ptr[0],1.0,-1.0);
    relative_pos[X] = relative_pos[X]/(self->vx[k+1] - self->vx[k])*2.0 - 1.0;
    relative_pos[Y] = relative_pos[Y]/(self->vy[j+1] - self->vy[j])*2.0 - 1.0;

    // interpolate the displacement of tracer
    double tracer_dis[2];
    double tracer_vel[2];
    double tracer_bf[2];
    if(evof_ptr[VACUUM] > TOLVOF)
    {
        // if the tracer is in void cells, interpolate its acceleration
        // and move according
        tracer_bf[X] = Interpolate2d_k(vacc_ptr,relative_pos,X);
        tracer_bf[Y] = Interpolate2d_k(vacc_ptr,relative_pos,Y);

        liner_op(tracer_dis,2,node->data.velocity,tracer_bf,dt,0.5*dt*dt);

        // move the tracer, update velocity
        scaler_add(node->data.velocity,2,tracer_bf,dt);
        scaler_add(node->data.position,2,tracer_dis,1.0);
    }
    else
    {
        tracer_dis[X] = Interpolate2d((double []){vdis_ptr[0][0],vdis_ptr[1][0],vdis_ptr[2][0],vdis_ptr[3][0]},relative_pos);
        tracer_dis[Y] = Interpolate2d((double []){vdis_ptr[0][1],vdis_ptr[1][1],vdis_ptr[2][1],vdis_ptr[3][1]},relative_pos);

        tracer_vel[X] = Interpolate2d_k(vvel_ptr,relative_pos,X);
        tracer_vel[Y] = Interpolate2d_k(vvel_ptr,relative_pos,Y);

        scaler_add(node->data.position,2,tracer_dis,1.0);
        scaler_move(node->data.velocity,2,tracer_vel,1.0);
    }

    // check whether this tracer is the ejecta, the ejecta is represented by negative matid
    if(node->data.matid > 0 && node->data.position[Y] > self->threshold)
    {
        node->data.matid = -node->data.matid;
    }

    // record the max pressure/temperature
    double *epre_ptr;
    double *etem_ptr;
    fv_id_2d elid = {k,j};
    self->foe->get_double(&elid,&epre_ptr,self->e_pre);
    self->foe->get_double(&elid,&etem_ptr,self->e_tem);

    node->data.mpre = MAX(node->data.mpre,*epre_ptr);
    node->data.mtem = MAX(node->data.mtem,*etem_ptr);
    node->record.pre = *epre_ptr;
    node->record.tmp = *etem_ptr;

    // update the id of tracer, this tracer should in cell(j-1,k-1) to cell(j+1,k+1)
    // because the displacement < min(delta_x,delta_y)
    int x_offset = GET_OFFSET2(node->data.position[X],self->vx[k],self->vx[k+1]);
    int y_offset = GET_OFFSET2(node->data.position[Y],self->vy[j],self->vy[j+1]);

    node->id.x += x_offset;
    node->id.y += y_offset;

    int nei_id = fl_node_nei(&(node->id),self->proc_self);

    if(nei_id < 0)
    {
        fprintf(stderr,"tracer position err [%d]\n",node->data.tag);
        fl_del(self,node);
        return field_list_deleted;
    }
    else if(proc_self->neighbour[nei_id] < 0)
    {
        fl_del(self,node);
        return field_list_deleted;
    }
    else if(proc_self->neighbour[nei_id] == self->rank)
    {
        node->state = field_list_used;
        node->target_proc  = nei_id;
        return field_list_used;
    }
    else
    {
        node->state = field_list_send;
        node->target_proc  = nei_id;
        return field_list_send;
    }

}

void fl_send_node_init(field_list * self, int * _list_len)
{
    field_list_node * node = NULL;
    field_list_iterator * fit = NULL;
    int send_node_list_len[NNP2] = {0};

    if(NULL == _list_len)
    {
        fit = fl_iterator_new(self,field_list_head);
        while(NULL != (node = fl_iterator_next(fit)))
        {
            if(node->state == field_list_send)
            {
                send_node_list_len[node->target_proc] ++;
            }
        }
        fl_iterator_del(fit);
    } else
    {
        for(int k=0;k<NNP2;++k) send_node_list_len[k] = _list_len[k];
    }


    // malloc data for send_list
    fl_data * send_node_list[NNP2];
    for(int node_i=0;node_i<NNP2;++node_i)
    {
        if(send_node_list_len[node_i] > 0)
        {
            send_node_list[node_i] = (fl_data *) MV_MALLOC(sizeof(fl_data)*send_node_list_len[node_i]);
        }
        else
        {
            send_node_list[node_i] = NULL;
        }
    }

    // collect tracer in send_node_list
    fit = fl_iterator_new(self,field_list_head);
    int send_node_list_k[NNP2] = {0};
    while(NULL != (node = fl_iterator_next(fit)))
    {
        if(node->state == field_list_send)
        {
            memcpy(send_node_list[node->target_proc] + send_node_list_k[node->target_proc],
                   &(node->data), sizeof(fl_data));
            send_node_list_k[node->target_proc] ++;
            fl_del(self, node);
        }
    }
    fl_iterator_del(fit);

    // at last copy send_node_list_len into send_rtb
    proc_info_2d * proc_self = (proc_info_2d *) self->proc_self;
    for(int k=0;k<self->send_rtb_size;++k)
    {
        self->send_rtb[k].priority = 0;
        self->send_rtb[k].data     = NULL;
        for(int j=0;j<NNP2;++j)
        {
            if(proc_self->neighbour[j] == self->send_rtb[k].rank)
            {
                self->send_rtb[k].priority = send_node_list_len[j];
                self->send_rtb[k].data = (char *) send_node_list[j];
                break;
            }
        }
    }
}

void fl_recv_node_init(field_list * self)
{
    for(int k=0;k<self->recv_rtb_size;++k)
    {
        if(self->recv_rtb[k].priority > 0)
        {
            self->recv_rtb[k].data = MV_MALLOC(sizeof(fl_data)*self->recv_rtb[k].priority);
        }
        else
        {
            self->recv_rtb[k].data = NULL;
        }
    }
}

void fl_tracer_update(field_list * self, double dt)
{
    /*
     * move/update tracer position
     * if the tracer need send/deleted, it should be marked as field_list_send/field_list_del
     * the del operation will be completed by fl_check
     */

    field_list_iterator * fit = fl_iterator_new(self,field_list_head);
    field_list_node * node = NULL;

    int send_node_list_len[NNP2] = {0};
    while(NULL != (node = fl_iterator_next(fit)))
    {
        field_list_state node_state = fl_update_node(self,node,dt);
        if(node_state == field_list_send)
        {
            send_node_list_len[node->target_proc] += 1;
        }
    }
    fl_iterator_del(fit);
    fl_check(self); // clean some dirty memory
    fl_send_node_init(self,send_node_list_len);
}

void fl_mpi_init(field_list * self, field_var * e_pos, field_var * e_vof,field_opt * foe,
                 field_var * v_pos, field_var * v_dis, field_var * v_vel,field_var * v_bf,field_opt * fo)
{
    /*
     * just copy some process info from field_var,
     * in most time the fv should be field_var on elements
     */
    MPI_Comm_dup(e_pos->comm, &self->comm);
    self->proc      = e_pos->proc;
    self->domain    = e_pos->domain;
    self->proc_num  = e_pos->proc_num;
    self->proc_self = e_pos->proc_self;
    self->rank      = e_pos->rank;
    MPI_Comm_rank(self->comm,&self->rank);

    // in field_list use route_table.priority as the number of nodes need to be exchange
    self->send_rtb_size = e_pos->send_list_len;
    self->send_rtb = (route_table *) MV_MALLOC(sizeof(route_table)*self->send_rtb_size);
    self->send_info_req = (MPI_Request *) MV_MALLOC(sizeof(MPI_Request) * self->send_rtb_size);
    for(int k=0; k < e_pos->send_list_len; ++k)
    {
        memcpy(self->send_rtb +k, e_pos->send_rtb + e_pos->send_list[k].head, sizeof(route_table));
        self->send_rtb[k].priority = 0;
        self->send_rtb[k].data = NULL;
        MPI_Send_init(&(self->send_rtb[k].priority),1,MPI_INT,self->send_rtb[k].rank,self->send_rtb[k].tag,self->comm, self->send_info_req + k);

    }

    self->recv_rtb_size = e_pos->recv_list_len;
    self->recv_rtb = (route_table *) MV_MALLOC(sizeof(route_table)*self->recv_rtb_size);
    self->recv_info_req = (MPI_Request *) MV_MALLOC(sizeof(MPI_Request) * self->recv_rtb_size);
    for(int k=0; k < e_pos->recv_list_len; ++k)
    {
        memcpy(self->recv_rtb +k, e_pos->recv_rtb + e_pos->recv_list[k].head, sizeof(route_table));
        self->recv_rtb[k].priority = 0;
        self->recv_rtb[k].data = NULL;
        MPI_Recv_init(&(self->recv_rtb[k].priority),1,MPI_INT,self->recv_rtb[k].rank, self->recv_rtb[k].tag,self->comm, self->recv_info_req + k);
    }

    // malloc mem for node exchange request
    self->send_node_req = (MPI_Request *) MV_MALLOC(sizeof(MPI_Request) * self->send_rtb_size);
    self->recv_node_req = (MPI_Request *) MV_MALLOC(sizeof(MPI_Request) * self->recv_rtb_size);

    // set velocity and velocity field in field_list
    self->v_pos = v_pos;
    self->v_dis = v_dis;
    self->v_vel = v_vel;
    self->v_bf  = v_bf;
    self->fov = fo;
    self->e_pos = e_pos;
    self->e_vof = e_vof;
    self->foe = foe;

    proc_info_2d * proc_pos = (proc_info_2d *) v_pos->proc_self;
    self->nvx = proc_pos->nx;
    self->nvy = proc_pos->ny;

    self->vx = (double *) MV_MALLOC(sizeof(double)*proc_pos->nx);
    self->vy = (double *) MV_MALLOC(sizeof(double)*proc_pos->ny);

    for(int k=0;k<proc_pos->nx;++k)
    {
        fv_id_2d lid = {.x=k,.y=1};
        double * vpos_ptr;
        fo->get_double(&lid,&vpos_ptr,v_pos);
        self->vx[k] = vpos_ptr[0];
    }

    for(int j=0;j<proc_pos->ny;++j)
    {
        fv_id_2d lid = {.x=1,.y=j};
        double *vpos_ptr;
        fo->get_double(&lid,&vpos_ptr,v_pos);
        self->vy[j] = vpos_ptr[1];
    }

}

void fl_record_init(field_list * self, field_var * e_den, field_var * e_tmp, field_var * e_pre, field_var * v_acc)
{
    // just copy the pointer
    // those filed var will be used in output
    // the field operator has been set in fl_mpi_init
    self->e_tem = e_tmp;
    self->e_pre = e_pre;
    self->e_den = e_den;

    // if other variables needed to recorded
    // insert in following .....
    self->v_acc = v_acc;

}

void fl_tracer_init(field_list * self, int strip)
{

    if(strip < 0) strip = -strip;

    proc_info_2d * proc_self = (proc_info_2d *) self->proc_self;
    proc_info_2d * proc_domain = (proc_info_2d *) self->domain;

    int nmat = self->e_vof->bytes_per_unit/ sizeof(double);

    for(int j=0;j< proc_self->ny;++j)
    {
        for(int k=0;k< proc_self->nx;++k)
        {
            fv_id_2d lid = {.x=k,.y=j};

            node_type local_nt = useless_section;
            self->foe->get_nty(&lid,&local_nt,self->e_pos);
            if(!(local_nt & (internal_section|send_section))) continue;

            fv_id_2d gid = {0};
            self->foe->lid2gid(self->rank,&lid,&gid,self->e_pos);

            if(gid.x % strip != 0 || gid.y % strip != 0) continue;

            fl_data tmp_data;

            // set global id
            tmp_data.gx    = gid.x;
            tmp_data.gy    = gid.y;
            tmp_data.tag   = LID2INDEX(gid.x,gid.y,proc_domain->nx,proc_domain->ny);

            double * epos_ptr = NULL;
            self->foe->get_double(&lid,&epos_ptr,self->e_pos);
            tmp_data.position[X] = epos_ptr[X];
            tmp_data.position[Y] = epos_ptr[Y];

            double * evof_ptr = NULL;
            self->foe->get_double(&lid,&evof_ptr,self->e_vof);
            tmp_data.matid = 0;
            for(int matid=1;matid<nmat;++matid)
            {
                if(evof_ptr[tmp_data.matid] < evof_ptr[matid])
                {
                    tmp_data.matid = matid;
                }
            }

            fv_id_2d lb = {k, j}, rb = {k + 1, j}, ru = {k + 1, j + 1}, lu = {k, j + 1};
            fv_id_2d *vlid[4] = {&lb, &rb, &ru, &lu};

            double *vvel_ptr[4];
            field_opt * fo = self->fov;
            for(int vid=0;vid<4;++vid)
            {
                fo->get_double(vlid[vid],&vvel_ptr[vid],self->v_vel);
            }
            tmp_data.velocity[X] = Interpolate2d_k(vvel_ptr,(double []){0.,0.},X);
            tmp_data.velocity[Y] = Interpolate2d_k(vvel_ptr,(double []){0.,0.},Y);

            if(0 != tmp_data.matid)
            {
                field_list_node * tmp_node = fl_node_new(&tmp_data);
                tmp_node->id.x = k;
                tmp_node->id.y = j;
                fl_rpush(self,tmp_node);
            }

        }
    }
}

