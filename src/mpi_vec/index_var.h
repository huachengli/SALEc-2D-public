//
// Created by huach on 11/7/2023.
//

#ifndef SALE_REBUILD_INDEX_VAR_H
#define SALE_REBUILD_INDEX_VAR_H

#include "mpi_vec.h"

typedef struct index_var_impl
{
    char * data;

    char * (*get)(fv_id_2d * lid);
    void (*sync_start)(struct index_var_impl * source);
    void (*sync_complete)(struct index_var_impl * source);
    MPI_Comm comm;
} index_var;

mpi_vec_state init_index_2d(field_2d * fi, field_var * fv);
mpi_vec_state sync_index_start(field_var * fv,MPI_Op sOP);
mpi_vec_state sync_index_complete(field_var * fv);
mpi_vec_state lid2double_index2d_default(void * _lid, double ** dest, field_var * var);
mpi_vec_state lid2data_index2d_default(void * _lid, char ** dest, field_var * var);
mpi_vec_state init_index_opt_2d_default(field_opt * fo);

#endif //SALE_REBUILD_INDEX_VAR_H
