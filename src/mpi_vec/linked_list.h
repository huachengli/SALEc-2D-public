//
// Created by huach on 22/5/2023.
//

#ifndef SALE_REBUILD_LINKED_LIST_H
#define SALE_REBUILD_LINKED_LIST_H

#include "mpi_vec_types.h"
#include "mpi_vec.h"
#include "stdlib.h"
#ifndef MV_MALLOC
#define MV_MALLOC malloc
#endif
#ifndef MV_FREE
#define MV_FREE free
#endif

typedef struct Tracer
{
    double position[2];
    double velocity[2];
    float mpre;
    float mtem;
    int gx;
    int gy;
    int tag;
    int matid;
} tracer_t;

typedef struct Recorder
{
    double pre;
    double tmp;
    double den;
} recorder_t;

typedef enum FieldListState
{
    field_list_used    = 0x01,
    field_list_unused  = 0x02,
    field_list_new     = 0x04,
    field_list_deleted = 0x08,
    field_list_send    = 0x10
} field_list_state;

#define fl_data tracer_t

typedef struct field_list_node_impl
{
    fv_id_2d id;
    fl_data data;
    recorder_t record;
    int target_proc;
    field_list_state state;
    struct field_list_node_impl * prev;
    struct field_list_node_impl * next;
} field_list_node;


typedef struct FieldList
{
    // section for basic list
    field_list_node * head;
    field_list_node * tail;
    unsigned int dirty_nodes;
    unsigned int len;
    int (*match)(fl_data *a, fl_data *b);
    field_list_state (*_add)(field_list_node * node);

    // section for mpi
    MPI_Comm comm;
    void * proc;
    void * proc_self;
    void * domain;
    int proc_num;
    int rank;

    route_table * send_rtb;
    route_table * recv_rtb;
    unsigned int send_rtb_size;
    unsigned int recv_rtb_size;

    MPI_Request * send_info_req;
    MPI_Request * recv_info_req;
    MPI_Request * send_node_req;
    MPI_Request * recv_node_req;

    // section for field var require when update fl_data
    field_var * v_pos;
    field_var * v_dis;
    field_var * v_vel;
    field_var * v_bf;
    field_var * v_acc;
    field_opt * fov;
    field_var * e_vof;
    field_var * e_pos;
    field_opt * foe;
    double * vx;
    double * vy;
    int nvx;
    int nvy;

    // data section for update record
    field_var * e_tem;
    field_var * e_den;
    field_var * e_pre;

    // threshold for ejecta
    double threshold;
} field_list;

typedef enum FieldListDirection
{
    field_list_head,
    field_list_tail
} field_list_direction;

typedef struct FieldListIterator
{
    field_list_node * next;
    field_list_direction direction;
} field_list_iterator;


field_list_node * fl_node_new(fl_data * _t);
void fl_node_del(field_list_node * node);
field_list * fl_new();
field_list_iterator * fl_iterator_new(field_list * self, field_list_direction direction);
field_list_iterator * fl_iterator_init(field_list_iterator* self_iterator,field_list * self, field_list_direction direction);
void fl_iterator_del(field_list_iterator * self);
field_list_node * fl_iterator_next(field_list_iterator * self);

// basic operation
field_list_node * fl_rpush(field_list * self, field_list_node * node);
field_list_node * fl_lpush(field_list * self, field_list_node * node);
field_list_node * fl_rpop(field_list * self);
field_list_node * fl_lpop(field_list * self);
field_list_node * fl_find(field_list * self, void * val);
field_list_node * fl_at(field_list * self, int index);
void fl_remove(field_list *self, field_list_node *node);
void fl_del(field_list * self, field_list_node * node);
void fl_check(field_list * self);
field_list_node * fl_add(field_list * self, fl_data * data);
field_list_node * fl_add_n(field_list * self, fl_data * data, int len);

void fl_sync_node_start(field_list * self);
void fl_sync_node_complete(field_list * self);
void fl_sync_info_start(field_list * self);
void fl_sync_info_complete(field_list * self);

#endif //SALE_REBUILD_LINKED_LIST_H
