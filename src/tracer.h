//
// Created by huach on 25/5/2023.
//

#ifndef SALE_REBUILD_TRACER_H
#define SALE_REBUILD_TRACER_H


#include <linked_list.h>
#include <interpolate.h>
#include "sale2d_constant.h"

#define GET_OFFSET2(t, x, y) (((t) < (x)) ? -1 : ((t) < (y) && (t) >= (x)) ? 0 : 1)
#define GET_OFFSET4(t, x, y, z, u) (((t) < (x) || (t) > (u)) ? -1 : ((t) >= (x) && (t) < (y)) ? 0 : ((t) >= (y) && (t) <= (z)) ? 1 : 2)

void fl_send_node_init(field_list * self, int * _list_len);
void fl_recv_node_init(field_list * self);
field_list_state fl_update_node(field_list * self, field_list_node * node,double dt);
void fl_tracer_update(field_list * self,double dt);

void fl_tracer_init(field_list * self, int strip);
void fl_mpi_init(field_list * self, field_var * e_pos, field_var * e_vof,field_opt * foe,
                 field_var * v_pos, field_var * v_dis, field_var * v_vel,field_var * v_bf,field_opt * fo);
void fl_record_init(field_list * self, field_var * e_den, field_var * e_tmp, field_var * e_pre, field_var * v_acc);
#endif //SALE_REBUILD_TRACER_H
