//
// Created by huach on 22/5/2023.
//

#include "linked_list.h"

field_list_node * fl_node_new(fl_data * _t)
{
    field_list_node * self = (field_list_node *) MV_MALLOC(sizeof(field_list_node));
    self->prev = NULL;
    self->next = NULL;
    memcpy(&(self->data), _t, sizeof(tracer_t));
    self->state = field_list_used;
    self->id.x = 0;
    self->id.y = 0;
    return self;
}

void fl_node_del(field_list_node * node)
{
    if(NULL != node)
    {
        MV_FREE(node);
    }
}

field_list * fl_new()
{
    field_list * self = (field_list *) MV_MALLOC(sizeof(field_list));
    self->head = NULL;
    self->tail = NULL;
    self->len  = 0;
    self->dirty_nodes = 0;
    self->match = NULL;
    self->_add  = NULL;
    self->proc_self = NULL;
    self->proc = NULL;
    self->domain = NULL;
    return self;
}

field_list_node * fl_rpush(field_list * self, field_list_node * node)
{
    if(node == NULL) return NULL;
    if(self->len >= 1)
    {
        node->prev = self->tail;
        node->next = NULL;
        self->tail->next = node;
        self->tail = node;
    }
    else
    {
        self->head = self->tail = node;
        node->prev = node->next = NULL;
    }

    ++self->len;
    return node;
}

field_list_node * fl_lpush(field_list * self, field_list_node * node)
{
    if(node == NULL) return node;

    if(self->len >= 1)
    {
        node->prev = NULL;
        node->next = self->head;
        self->head->prev = node;
        self->head = node;
    }
    else
    {
        self->head = self->tail = node;
        node->prev = node->next = NULL;
    }

    ++ self->len;
    return node;

}

field_list_node * fl_rpop(field_list * self)
{
    if(self->len <= 0) return NULL;
    field_list_node * node = self->tail;
    self->len --;
    if(node->state == field_list_deleted) self->dirty_nodes--;

    if(self->len >= 1)
    {
        self->tail = node->prev;
        self->tail->next = NULL;
    }
    else
    {
        self->tail = self->head = NULL;
    }

    node->next = node->prev = NULL;
    return node;
}

field_list_node * fl_lpop(field_list * self)
{
    if(self->len <= 0) return NULL;
    field_list_node * node = self->head;
    self->len --;
    if(node->state == field_list_deleted) self->dirty_nodes--;

    if(self->len >= 1)
    {
        self->head = node->next;
        self->head->prev = NULL;
    }
    else
    {
        self->tail = self->head = NULL;
    }

    node->next = node->prev = NULL;
    return node;
}


void fl_remove(field_list *self, field_list_node *node) {
    node->prev
    ? (node->prev->next = node->next)
    : (self->head = node->next);

    node->next
    ? (node->next->prev = node->prev)
    : (self->tail = node->prev);

    MV_FREE(node);
    --self->len;
}

void fl_del(field_list * self, field_list_node * node)
{
    node->state = field_list_deleted;
    self->dirty_nodes ++;
}

void fl_check(field_list * self)
{
    if(self->dirty_nodes*2 > self->len && self->dirty_nodes >= 16)
    {
        field_list_iterator * fit = fl_iterator_new(self,field_list_head);
        field_list_node * cur = NULL;
        while( NULL != (cur = fl_iterator_next(fit)))
        {
            if(self->dirty_nodes == 0) break;
            if(cur->state == field_list_deleted)
            {
                fl_remove(self,cur);
                self->dirty_nodes --;
            }
        }
        fl_iterator_del(fit);
    }
}

field_list_node * fl_add(field_list * self, fl_data * data)
{
    if(data == NULL)
    {
        return NULL;
    }
    else if(self->dirty_nodes == 0)
    {
        field_list_node * node = fl_node_new(data);
        node->state = field_list_new;
        return fl_rpush(self,node);
    }
    else
    {
        field_list_iterator * fit = fl_iterator_new(self,field_list_head);
        field_list_node * cur = NULL;
        while(NULL != (cur = fl_iterator_next(fit)))
        {
            if(cur->state == field_list_deleted)
            {
                memcpy(&(cur->data),data, sizeof(fl_data));
                cur->state = field_list_new;
                self->dirty_nodes--;
                break;
            }
        }
        fl_iterator_del(fit);
        return cur;
    }
}

field_list_node * fl_add_n(field_list * self, fl_data * data, int len)
{
    field_list_node * node = NULL;
    for(int k=0;k<len;++k)
    {
        node = fl_add(self,data + k);
    }
    return node;
}

field_list_node * fl_add_n2(field_list * self, fl_data * data, int len)
{
    field_list_node * node = NULL;
    field_list_iterator * fit = fl_iterator_new(self,field_list_head);

    for(int k=0;k<len;++k)
    {
        if(self->dirty_nodes <= 0)
        {
            field_list_node * tmp_node = fl_node_new(data);
            tmp_node->state = field_list_new;
            fl_rpush(self,tmp_node);
        }
        else
        {
            field_list_node * cur = NULL;
            while(NULL != (cur = fl_iterator_next(fit)))
            {
                if(cur->state == field_list_deleted)
                {
                    memcpy(&(cur->data),data, sizeof(fl_data));
                    cur->state = field_list_new;
                    self->dirty_nodes--;
                    break;
                }
            }
        }
    }
    fl_iterator_del(fit);
    return node;
}


field_list_iterator * fl_iterator_new(field_list * self, field_list_direction direction)
{
    field_list_iterator * self_iterator = (field_list_iterator *) MV_MALLOC(sizeof(field_list_iterator));
    self_iterator->direction = direction;
    self_iterator->next = (direction == field_list_head)? self->head : self->tail;
    return self_iterator;
}

field_list_iterator * fl_iterator_init(field_list_iterator* self_iterator,field_list * self, field_list_direction direction)
{
    self_iterator->direction = direction;
    self_iterator->next = (direction == field_list_head)? self->head : self->tail;
    return self_iterator;
}

void fl_iterator_del(field_list_iterator * self)
{
    MV_FREE(self);
}

field_list_node * fl_iterator_next(field_list_iterator * self)
{
    field_list_node * rtn_node = self->next;
    if(rtn_node != NULL)
    {
        self->next = (self->direction == field_list_head)? rtn_node->next : rtn_node->prev;
    }
    return rtn_node;
}

field_list_node * fl_find(field_list * self, void * val)
{
    field_list_iterator * it = fl_iterator_new(self, field_list_head);
    field_list_node * node;
    while ((node = fl_iterator_next(it))) {
        if (self->match) {
            if (self->match(val, &node->data)) {
                fl_iterator_del(it);
                return node;
            }
        } else {
            if (val == &node->data) {
                fl_iterator_del(it);
                return node;
            }
        }
    }

    fl_iterator_del(it);
    return NULL;
}

field_list_node * fl_at(field_list * self, int index)
{
    field_list_direction direction = field_list_head;
    if (index < 0) {
        direction = field_list_tail;
        index = ~index;
    }

    if ((unsigned int)index < self->len)
    {
        field_list_iterator *it = fl_iterator_new(self, direction);
        field_list_node *node = fl_iterator_next(it);
        while(index--) node = fl_iterator_next(it);
        fl_iterator_del(it);
        return node;
    }

    return NULL;
}

void fl_sync_node_start(field_list * self)
{
#ifdef DEBUG_RANK
    MPI_Barrier(self->comm);
#endif

    for(int k=0;k<self->send_rtb_size;++k)
    {
        MPI_Isend(self->send_rtb[k].data,self->send_rtb[k].priority* sizeof(fl_data),MPI_CHAR,
                  self->send_rtb[k].rank,self->send_rtb[k].tag,self->comm,self->send_node_req+k);
    }

    for(int k=0;k<self->recv_rtb_size;++k)
    {
        MPI_Irecv(self->recv_rtb[k].data,self->recv_rtb[k].priority* sizeof(fl_data),MPI_CHAR,
                  self->recv_rtb[k].rank,self->recv_rtb[k].tag,self->comm,self->recv_node_req+k);
    }
}

void fl_sync_node_complete(field_list * self)
{
    MPI_Waitall(self->send_rtb_size,self->send_node_req,MPI_STATUSES_IGNORE);
    MPI_Waitall(self->recv_rtb_size,self->recv_node_req,MPI_STATUSES_IGNORE);

#ifdef DEBUG_RANK
    MPI_Barrier(self->comm);
#endif

    // after sync_node, free memory allocated in send/recv_rtb
    // copy data in recv_rtb to fl

    // delete mem for send
    for(int k=0;k<self->send_rtb_size;++k)
    {
        if(self->send_rtb[k].priority > 0)
        {
            MV_FREE(self->send_rtb[k].data);
        }
    }

    // copy recv data into fl before del
    for(int k=0;k<self->recv_rtb_size;++k)
    {
        if(self->recv_rtb[k].priority > 0)
        {
            fl_data * _recv_data = (fl_data *)self->recv_rtb[k].data;
            fl_add_n(self,_recv_data,self->recv_rtb[k].priority);
            MV_FREE(_recv_data);
        }
    }

}

void fl_sync_info_start(field_list * self)
{
#ifdef DEBUG_RANK
    MPI_Barrier(self->comm);
#endif
    MPI_Startall(self->send_rtb_size,self->send_info_req);
    MPI_Startall(self->recv_rtb_size,self->recv_info_req);
}

void fl_sync_info_complete(field_list * self)
{
    MPI_Waitall(self->send_rtb_size,self->send_info_req,MPI_STATUSES_IGNORE);
    MPI_Waitall(self->recv_rtb_size,self->recv_info_req,MPI_STATUSES_IGNORE);

#ifdef DEBUG_RANK
    MPI_Barrier(self->comm);
#endif

    // after recv the num of tracers will be received,
    // allocate the memory
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
