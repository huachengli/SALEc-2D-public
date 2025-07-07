//
// Created by huacheng on 10/28/22.
//

#ifndef MPI_VEC_MPI_VEC_H
#define MPI_VEC_MPI_VEC_H

#include "mpi_vec_types.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "linear.h"


#define RESTRICT_PHY_VERTEX 1
#define COLOUM_ARRAY 1
#ifdef COLOUM_ARRAY
#define LID2INDEX(ix,jy,nx,ny) ((ix) + (nx)*(jy))
#define INDEX2LIDX(index,nx,ny) ((index)%(nx))
#define INDEX2LIDY(index,nx,ny) ((index)/(nx))
#else
#define LID2INDEX(ix,jy,nx,ny) ((jy) + (ny)*(ix))
#define INDEX2LIDX(index,nx,ny) ((index)/(ny))
#define INDEX2LIDY(index,nx,ny) ((index)%(ny))
#endif

int cmp_route_table(const void * x, const void * y);
mpi_vec_state gid2lid_2d_default(int rank, void* _gid, void* _lid, field_var * var);
mpi_vec_state lid2gid_2d_default(int rank, void* _lid, void* _gid, field_var * var);
int lid2pri_2d_default(int rank,void* _lid, field_var * var);
int gid2tag_2d_default(int rank,void* _gid,field_var * fv);
mpi_vec_state node_type_2d_default_element(int rank, void* _lid, node_type *, field_var * var);
mpi_vec_state node_type_2d_default_vertex(int rank, void* _lid, node_type *nty, field_var * var);
mpi_vec_state lid2index_2d_default(int rank, void * _lid, int * _index, field_var * fv);
mpi_vec_state index2lid_2d_default(int rank, int _index, void * _lid, field_var *fv);
mpi_vec_state lid2data_2d_default(void * _lid, char ** dest, field_var * var);
mpi_vec_state lid2double_2d_default(void * _lid, double ** dest, field_var * var);
mpi_vec_state gid2boundary_2d_default(void * _gid, void* _bound, field_var * var);
mpi_vec_state lid2data_nocheck_2d_default(void * _lid, char ** dest, field_var * var);
mpi_vec_state lid2node_type_default(void * _lid, node_type * ntype, field_var * fv);
mpi_vec_state lid2node_type_nocheck_default(void * _lid, node_type * ntype, field_var * fv);

mpi_vec_state init_field_proc_2d_default(field_2d *fi, field_var *fv);
mpi_vec_state set_null_ght(field_var * fv);
mpi_vec_state init_field_opt_2d_default(field_opt * fo);
mpi_vec_state init_field_route_table(field_var * fv, field_opt * fo);
mpi_vec_state init_field_ght_table(field_var * fv, field_opt * fo);
mpi_vec_state init_field_rst_table(field_var * fv, field_opt * fo);
void ghost_pad_move(fv_id_2d * lid,field_var * fv);
mpi_vec_state ghost_pad_fill(ghost_table * gt, field_var * fv);
mpi_vec_state ghost_opt_fill(ghost_table * gt, field_var * fv, ght_opt _opt);
mpi_vec_state sync_ghost_cells(field_var * fv);
mpi_vec_state sync_restrict_vertexes(field_var * fv);
mpi_vec_state set_ghost_pattern(ght_opt * go_list, field_var * fv, ght_var_type var_type);
mpi_vec_state set_restrict_pattern(ght_opt * go_list, field_var * fv, ght_var_type var_type);
mpi_vec_state compress_route_table(route_table * rtb, unsigned int rtb_size, route_table_continuous ** c_rtb, unsigned int * c_rtb_len);
mpi_vec_state init_field_buffer(field_var * fv);

mpi_vec_state init_field_2d(field_2d *fi, field_var *fv, field_opt *fo);
mpi_vec_state clean_field(field_var * fv);
mpi_vec_state sync_start(field_var * fv);
mpi_vec_state sync_complete(field_var * fv);
mpi_vec_state flush_ght(field_var * fv);

mpi_vec_state duplicate_field_2d(field_var *fv_src, field_var * fv_dst, field_opt *fo);
mpi_vec_state duplicate_field_2d_n(field_var *fv_src, field_var * fv_dst, field_opt *fo, int _n);


#endif //MPI_VEC_MPI_VEC_H
