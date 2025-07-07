//
// Created by huacheng on 10/28/22.
//

#ifndef MPI_VEC_MPI_VEC_TYPES_H
#define MPI_VEC_MPI_VEC_TYPES_H
#include <mpi.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

#ifndef MPI_VEC_MALLOC
#define MV_MALLOC malloc
#endif

#ifndef MPI_VEC_FREE
#define MV_FREE free
#endif


#define NNP2 9 // number of neighbour processor in 2d
typedef struct proc_info_2d_impl
{
    int neighbour[NNP2];
    int stx;
    int nx;
    int xpad;
    int sty;
    int ny;
    int ypad;
} proc_info_2d;


typedef enum route_table_type__
{
  send_route,
  recv_route,
  other_route
} rt_type;

typedef enum boundary_type_impl
{
    left_boundary   = 0x01,
    right_boundary  = 0x02,
    bottom_boundary = 0x04,
    top_boundary    = 0x08,
    check_boundary  = 0x0f,
    left_pad        = 0x10,
    right_pad       = 0x20,
    top_pad         = 0x40,
    bottom_pad      = 0x80,
    check_pad       = 0xf0,
    undefined_boundary = 0x00
} boundary_type;

typedef struct route_table_impl
{
    int rank;
    int tag;
    int priority;
    int gx;
    int gy;
    char * data;
} route_table;

typedef enum ght_table_opt_impl
{
    ght_outflow,
    ght_outflow_vertex,
    ght_fixed,
    ght_zero,
    ght_clean,
    ght_mirror,
    ght_freeslip,
    ght_reverse_vel_y,
    ght_reverse_vel_absy,
    ght_reverse_vel_x,
    ght_reverse_ste_1,
    ght_copy_ste,
    rst_x_zero,
    rst_y_zero,
    rst_y_positive,
    rst_x_positive,
    rst_zero,
    rst_none,
    ght_undefined
} ght_opt;

typedef enum ght_var_type_impl
{
    ght_var_always_copy,
    ght_normal,
    ght_velocity_inc,
    ght_velocity,
    ght_stress,
    ght_density
} ght_var_type;


typedef struct ght_table_impl
{
    char * data;
    char * ref_a;
    char * ref_b;
    char * ref_c;
    boundary_type tag;
    ght_opt opt;
} ghost_table;

typedef struct rst_table_impl
{
    char * data;
    char * ref_a;
    char * ref_b;
    boundary_type tag;
    ght_opt opt;
} restrict_table;


typedef struct route_table_continuous_impl
{
    int head;
    int len;
} route_table_continuous;

typedef enum field_var_pattern__
{
    Rect1,
    Rect2,
    Rect3
} fv_pattern;

typedef struct field_var_id_2d_impl
{
    int x;
    int y;
    int u;
    int v;
} fv_id_2d;

typedef enum pad_type_impl
{
    element_pad,
    vertex_pad,
    bx_pad,
    by_pad,
    x_pad,
    y_pad,
    unknown_pad
} pad_type;

typedef struct field_2d_impl
{
    int npgx;
    int npgy;
    int xpad;
    int ypad;
    int nx;
    int ny;
    int bytes_per_unit;
    pad_type padding;
    int * nx_list;
    int * ny_list;
} field_2d;

typedef enum node_type_impl
{
    receive_section = 0x01,
    send_section    = 0x02,
    internal_section= 0x04,
    update_section  = 0x08,
    useless_section = 0x10,
    share_section   = 0x20,
    pad_section     = 0x40
} node_type;

struct field_opt_impl;
typedef struct field_opt_impl field_opt;

typedef struct field_var_impl
{
    char * data;
    node_type * data_type;
    unsigned int data_size;
    unsigned int bytes_per_unit;

    route_table * send_rtb;
    route_table * recv_rtb;
    unsigned int send_rtb_size;
    unsigned int recv_rtb_size;

    route_table_continuous * send_list;
    route_table_continuous * recv_list;
    unsigned int send_list_len;
    unsigned int recv_list_len;

    ghost_table * gt_list;
    unsigned int gt_list_len;

    ghost_table * gt2_list;
    unsigned int gt2_list_len;

    int sync_gt_state;

    restrict_table * rt_list;
    unsigned int rt_list_len;

    char * send_buffer;
    unsigned int send_buffer_size;
    MPI_Request * send_req;
    unsigned int num_send_req;

    char * recv_buffer;
    unsigned int recv_buffer_size;
    MPI_Request * recv_req;
    unsigned int num_recv_req;

    void * proc;
    void * proc_self;
    int proc_num;
    void * domain;
    int rank;
    MPI_Comm comm;
    int dup_key;
    int dup_attr;

    pad_type padding;
} field_var;


typedef enum mpi_vec_check_impl
{
    success       = 0,
    error         = 1,
    id_error      = 2,
    malloc_error  = 4,
    nullptr_error = 8,
    unknown      = 16,
    duplicate_error = 32
} mpi_vec_state;

struct field_opt_impl
{
    mpi_vec_state (*gid2lid)(int, void*, void*, field_var *);
    mpi_vec_state (*lid2gid)(int, void*, void*, field_var *);
    int (*lid2pri)(int,void*, field_var *);
    int (*gid2tag)(int,void*, field_var *);
    mpi_vec_state (*lid2nty)(int, void*, node_type *, field_var *);
    mpi_vec_state (*get_data)(void*, char**, field_var *);
    mpi_vec_state (*get_data2)(void*, char**, field_var *);
    mpi_vec_state (*get_nty)(void*, node_type*, field_var*);
    mpi_vec_state (*get_nty2)(void*, node_type*, field_var*);
    mpi_vec_state (*index2lid)(int, int, void*, field_var*);
    mpi_vec_state (*lid2index)(int, void*, int*, field_var*);
    mpi_vec_state (*get_double)(void*, double **, field_var *);
    mpi_vec_state (*get_boundary)(void*, void*, field_var*);
};



#define MAX(a,b) ((a)>(b) ? (a):(b))
#define MIN(a,b) ((a)<(b) ? (a):(b))


#include <stdio.h>
#include <mpi.h>

#define PRINT_CONDITION(rank) (rank == DEBUG_RANK)  // Change 2 to the desired rank

// Define MPI_Printf macro
#define MPI_Printf(rank, ...) \
    do { \
        if (PRINT_CONDITION(rank)) { \
            fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__); \
            fprintf(stderr, __VA_ARGS__); \
        } \
    } while (0)


#ifdef COMM_OFF
#define MPI_Send(buf, count, datatype, dest, tag, comm) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Send called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Recv(buf, count, datatype, source, tag, comm, status) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Recv called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Allreduce called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Barrier(comm) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Barrier called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Startall(count, array_of_requests) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Startall called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Waitall(count, array_of_requests, array_of_statuses) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Waitall called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Isend(buf, count, datatype, dest, tag, comm, request) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Isend called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#define MPI_Irecv(buf, count, datatype, source, tag, comm, request) \
    do { \
        int my_rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
        if (PRINT_CONDITION(my_rank)) { \
            MPI_Printf(my_rank, "MPI_Irecv called, but disabled in single process debug mode.\n"); \
        } \
    } while (0)

#endif


#endif //MPI_VEC_MPI_VEC_TYPES_H
