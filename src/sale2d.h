//
// Created by huach on 14/11/2022.
//

#ifndef SALE_REBUILD_SALE2D_H
#define SALE_REBUILD_SALE2D_H
#include "interpolate.h"
#include "mpi_vec.h"
#include "index_var.h"
#include "InputParser.h"
#include "eos_state.h"
#include "output.h"
#include "tracer.h"
#include "sale2d_constant.h"
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

typedef enum FluxUpdateStrategy
{
    vaccum_priority,
    material_priority,
    mixed_strategy
} flux_strategy;

typedef enum StrengthMixed
{
    min_strength,
    max_strength,
    avg_strength,
    undefined_strength
} strength_mixed;

typedef struct sale2d_impl
{
    // variables on elements related to velocity
    field_var *v_vel;
    field_var *v_bf;
    field_var *v_acc;
    field_var *e_vel;
    field_var *v_dis;
    field_var *e_div_vel;
    field_var *e_grad_vel;
    field_var *e_ste; // stress in element
    field_var *e_mu;
    field_var *e_dvel; // the velocity increment of cells

    // state variables (one double)
    field_var *e_pre;
    field_var *e_tem;
    field_var *e_den;
    field_var *e_eng;
    field_var *e_csd; // sound of speed;
    field_var *e_dam;
    field_var *e_cvs; // cold volume strain.

    // state variables for different materials (nmat double)
    field_var *m_vof;
    field_var *m_den;
    field_var *m_eng;
    field_var *m_pty; // the porosity

    // variables defined before loop, set_mesh_pos
    field_var * v_pos;
    field_var * e_pos;
    field_var * e_vol; // 1 component, volume of element(cyl)
    field_var * e_suf;
    field_var * e_grad;
    field_var * v_laplace;


    field_var * e_subvol; // 8 components, 4 volume from vertexes (cart) and 4 volume (cyl)

    // VARIABLES ON BOUNDARY
    field_var *bx_phi;
    field_var *by_phi;
    field_opt *fobx;
    field_opt *foby;

    // variables for flux
    field_var * e_vof; // just copy of m_vof, in most of the time
    field_var * v_vof;
    field_var * e_wpt;
    field_var * v_mass;


    // debug variables
    field_var * e_debug;
    field_var * v_debug;

    // variables used for numerical stability
    field_var * e_q;

    // related to strength
    field_var * e_strength;
    field_var * e_tps;
    field_var * e_ste2s;
    field_var * e_sta2s;
    field_var * e_cnd; // some condition flag, updated in every cycle

    // some info collect from e_vol/e_subvol
    field_var * v_vol;    // 4 components, [0-3] volumes from cells(cart)
    field_var * v_subvol; // 5 components,[0] = sum of volumes from cells (cyl), [1-4] = volumes from cells(cyl)
    field_var * v_suf;

    // parameters from damping
    field_var * e_alpha; // pml auxiliary vars on cells
    field_var * v_alpha; // pml auxiliary vars on vertexes
    field_var * e_sigma; // damp factor
    field_var * v_beta; // reserved for damp on vertex
    field_var * e_dampref; // reference damp state (density, energy, pressure, ..)
    field_var * e_beta; // reserved for damp on element

    // normal vector required by algebraic reconstruction,
    field_var * bx_normal;
    field_var * by_normal;


    //variables used for acoustic fluid
    field_var * e_vib; // the square of vib velocity
    // RECORD THE VOLUME OF EJECTA MATERIALS
    field_var * e_ejt;
    int fv_counter; // number of field_var pointer declared in this struct
    int phi_length;
    // operation collection on vertex/element
    field_opt *foe;
    field_opt *fov;

    // variable to record profile dynamically
    field_var * y_profile;
    // operation for index_var
    field_opt * foi;

    int nmat;
    eos_table * etb;
    double avis1; // linear    artf. viscosity
    double avis2; // quadratic artf. viscosity
    double anc;
    double damp_time;

    double density_cutoff;
    double velocity_max;
    double velocity_min;
    double acceleration_min;

    // flags for some model
    int partpres;
    int acoustic_fluid;
    int porosity_model;
    strength_mixed strength_mod;
    int tensile_failure;
    flux_strategy flux_opt;

    // section for tracer
    field_list * tracer_list;

    // some args for damp
    // alpha_damp is only used in pml
    double alpha_damp_x;
    double alpha_damp_y;

    // height for trace ejecta materials
    double threshold_z;
    double flow_resistance;
    double flow_friction;
    double min_dx;

    // tracers lubricate
    int tracer_lubrication;
    int tracer_lub_gy0;
    int tracer_lub_gy1;
    double tracer_lub_y0;
    double tracer_lub_frac;
} sale2d_var;

typedef struct CycleControl
{
    double start_time;
    double end_time;
    double max_step;
    double step_time;
    double dt;
    double local_dt;
    double max_dt;


    double out_time;
    char prefix[100];
    int init_step;
    int k_out_step;
    int i_out_step;
    int i_step;

    int sync_vel;

    int out_tracer;
    int out_vts;
    int out_tracer_txt;
    int out_log;
    int out_init;
    int out_cycle;
    int cycle_continue;
    int rank;
    int x_rank;
    int y_rank;

    FILE * log;
} cycle_control;


// numerical center of element
/*
 * interval [a,b) f(x) = x, b = a + dx
 * integral f(x) = (b^3 - a^3)/3, volume of this interval (b^2 - a^2)/2
 * average = 2/3*(b^2+ab+a^2)/(a+b) = 2/3*(a + b^2/(a+b)) = a + 2/3*(a + dx)^2/(2a + dx) - a/3
 * = a + 1/3*((a+dx)^2/(a+dx/2) - a ) = a + 1/3*(a^2 + 2a*dx + dx^2 - a^2 - a*dx/2)/(a + dx/2)
 * = a + (a*dx/2 + dx^2/3)/(a + dx/2)
 */

typedef struct ProjectileDefition
{
    int material;
    double radiu;
    double center[3];
    double velocity[3];
    double damage;
    double pressure;
    double temperature;
} projectile_info;

void load_projectile_info(InputFile * ifp, projectile_info * _pinfo);


typedef enum TargetType
{
    sphere_target,
    plane_target,
    unknown_target
} target_type;

typedef enum TemProfType
{
    conductive_profile,
    convective_profile,
    condconv_profile,
    const_profile,
    lunar_profile,
    user_profile,
    unknown_profile
} tem_prof_type;

typedef enum DampOptType
{
    pml,
    pml2,
    pml3,
    sbc,
    ctlg,
    none
} damp_type_t;

#define MAXLAYERS 10



typedef struct TargetDefition
{
    int n_layer;
    int material[MAXLAYERS];
    double depth[MAXLAYERS];

    int    n_profile;
    double max_depth;

    double t_lith;
    double d_lith;
    double d_mantle;
    double t_mantle;
    double surf_gravity;
    double surf_tem;
    double surf_dtdz;
    double surf_pre;
    double r_planet;
    double center[3];
    target_type layer_type;
    tem_prof_type tem_type;

} target_info;

void load_target_info(InputFile * ifp, target_info * _tinfo);

typedef struct mesh2d_info_impl
{
    int npgx;
    int npgy;
    int xpad;
    int ypad;
    int nx;
    int ny;
    int v_nx;
    int v_ny;
    int e_nx;
    int e_ny;
    int nx_ext[2];
    int ny_ext[2];

    double oy;
    double ox;
    double dx;
    double dy;
    double x_ext; // the factor use to expand grid
    double y_ext;

    int n_materials;
    char ename[MAXMAT][200];
    eos_type etype[MAXMAT];

    char prefix[64];
    char format[64];

    target_info tinfo;
    projectile_info pinfo;
    double gravity[3];


    char ghost_boundary[4][20];
    boundary_type local_boundary[4];
    ght_opt go_list[4];
    int flux_restriction[4];

    double avis1;
    double avis2;
    double anc;
    double anc_depth;
    double anc_width;
    char anc_pattern[64];
    double ejecta_threshold;
    double ejecta_resistance;
    double ejecta_friction;
    double density_cutoff;
    double velocity_max;
    double velocity_min;
    double acceleration_min;
    int partpres;

    damp_type_t damp_type;
    double damp_t;
    int damp_x;
    int damp_y;
    int damp_x_left; // Warning in most time, this value should be 0, no damping zone around r=0
    int damp_y_top;
    double lambda_x;
    double lambda_y;

    int tracer_strip;

    state_reference * ref;
} mesh2d_info;

void init_sale2d_var(sale2d_var * _sale, mesh2d_info * _minfo);
void init_sale2d_ghost(sale2d_var * _sale, mesh2d_info * _minfo);
void clean_sale_var(sale2d_var * _sale);

void generate2dcyl_suf_and_vol(const double * _pos, double * _suf, double * _vol_sub, double * _vol_cyl);
void generate2dcyl_grad(const double * _pos, double * _grad);
void generate2d_center(const double * _pos, double * _center);
void generate2d_laplace(const double * _pos, double * _vl);

void init_mesh_pos(sale2d_var * _sale, mesh2d_info * _minfo);
void update_velocity(sale2d_var * _sale, node_type _nt,cycle_control * _ccl);
void update_energy(sale2d_var * _sale, node_type _nt, double dt);
void update_e_den(sale2d_var * _sale, node_type _nt);
void calculate_flux_average(sale2d_var * _sale, node_type _nt);
void calculate_flux_upwind(sale2d_var * _sale, node_type _nt);

void update_v_vof(sale2d_var * _sale, node_type _nt);
void calculate_v_dis(sale2d_var * _sale, node_type _nt, double dt);
void update_e_gradvel(sale2d_var * _sale, node_type _nt);

double cut_length(double _vf, double * _p);
double get_flux_volcyl(const double * base_a,const double * base_b, const double * dis_a, const double * dis_b);
double geometric_flux_ratio(double  _p1, double _p2);

void load_mesh_info(InputFile * ifp, mesh2d_info * _minfo, const char * eospath);
void load_state_reference(InputFile * ifp, state_reference * _sref, int nmat);
void init_sale_mat(sale2d_var * _sale, mesh2d_info * _minfo);
void init_sale2d_other(sale2d_var * _sale, InputFile * ifp);
void init_cycle_control(InputFile * ifp, cycle_control * _cyl);
void clean_cycle_control(cycle_control * _ccl);
void update_cycle_control(cycle_control * _ccl, FILE * fp);
void sale2d_write_log(FILE * fp, cycle_control * _ccl);
void sale2d_write_vtm(int nproc, int cycles, const char _prefix[]);

ght_opt boundary_str2opt(const char * boundary_str, boundary_type local_bt);
double calculate_dt(sale2d_var * _sale, double cur_dt);
void update_state(sale2d_var * _sale, node_type _nt);
void set_e_debug(sale2d_var * _sale);
double minimum_pressure(double dam, double minpint, double minpdam);

// content of advect.c
void calculate_flux_primary(sale2d_var * _sale, node_type _nt);
void advect_flux(damp_type_t dampType,sale2d_var * _sale, node_type _nt, double dt);
void advect_flux_default(sale2d_var * _sale, node_type _nt);
void advect_flux_pml3(sale2d_var * _sale, node_type _nt, double dt);
void clean_flux(sale2d_var * _sale);
//#define calculate_flux calculate_flux_primary
void _calculate_flux_set_bphi(sale2d_var * _sale, fv_id_2d * _elid, double lbru_vol[],double flux_vol[][4]);
double algebraic_flux_ratio(double VOFU, double VOFD, double VOFA, double CourantD,double bnormal[],double grad_vof[],
                            double (*vofAlgorithm)(double,double,double));
double vofMSTACS(double nVOFD, double Gf, double CourantD);
double vofCICSAM(double nVOFD, double Gf, double CourantD);
double vofSTACS(double nVOFD, double Gf);
double vofMHRIC(double nVOFD, double Gf);
void calculate_flux_mixed(sale2d_var * _sale, node_type _nt);
void calculate_flux(sale2d_var * _sale, node_type _nt);

// function: velocity increment
void update_v_vel_inc(sale2d_var * _sale, node_type _nt);
void check_v_vel(sale2d_var * _sale, node_type _nt);
void update_e_vel(sale2d_var * _sale, node_type _nt);
void clean_fv(field_var * _var);
void init_v_field(sale2d_var * _sale);
mpi_vec_state sync_all(sale2d_var * _sale);

// velocity/acceleration update with  pml boundary
void calculate_acc(damp_type_t dampType,sale2d_var * _sale, node_type _nt, double dt);
// pml boundary
double record_avg_value(sale2d_var * _sale);
void damp_on_cell(damp_type_t dampType,sale2d_var * _sale,double dt,double t);
void damp_on_cell_pml(sale2d_var * _sale,double dt,double t);
void damp_on_cell_pml2(sale2d_var * _sale,double dt,double t);
void damp_on_cell_pml3(sale2d_var * _sale,double dt,double t);
void damp_on_cell_sbc(sale2d_var * _sale,double dt,double t);
// tracer
int init_sale2d_tracer(sale2d_var * _sale,mesh2d_info * _mesh);
int export_init_tracers(sale2d_var * _sale, mesh2d_info * _mesh);

// other function
void make_empty_dir(const char * dir);
void put_cpuinfo(FILE * fp);

/*
 * some function used for grid
 */
double Spacing(int k,double ext, int eL, int eR, int N);
double Spacing2(int k,double ext, int eL,int eR, int dL,int dR,int N);
double Spacing3(int k,double ext, int eR, int dR,int N, int rpad);
double SpacingX(int k,double ext, int eL,int eR, int dL,int dR,int N,int rpad);
double SpacingY(int k,double ext, int eL,int eR, int dL,int dR,int N);
double hanning(double x);
#endif //SALE_REBUILD_SALE2D_H