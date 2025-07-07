//
// Created by huacheng on 2/23/23.
//

#include "eos_state.h"
#include "sale2d.h"
#include "stdbool.h"
#include "time.h"
#include <sys/stat.h>
#include "interpolate_2d9.h"
#include "sale2d.h"

void compare_binary_ascii(int test_size);
double test_ablinear(double scale);
void test_aneossieden();
int test_spacing2(double ext, int eL,int eR, int dL,int dR,int N);
int test_spacing3(double ext, int eR, int dR, int N, int rpad);

int test_Int2d9();
int test_laplace2d9(double gamma);
void test_aneossimple();


int main()
{
    // test the anoes
    test_aneossimple();

    // compare_binary_ascii(6*1024);
    // compare_binary_ascii(6*1024*1024);

    /*srand(time(NULL));
    double sum_delta = 0.0;
    for(int k=0;k< 1<<16 ;++k)
    {
        sum_delta += test_ablinear(10.0);
    }

    fprintf(stdout,"delta is %f\n",sum_delta);*/

//    test_aneossieden();
//    test_spacing2(1.05,8,3,0,0,60);
//    test_spacing3(1.05,3,4,50,2);
//    test_Int2d9();
    test_laplace2d9(1.0);
    test_laplace2d9(0.5);
    test_laplace2d9(sqrt(1./2.));
    test_laplace2d9(sqrt(3./4.));
//
//    test_laplace2d9(1.0);
//    test_laplace2d9(0.5);
//    test_laplace2d9(sqrt(1./2.));
//    test_laplace2d9(sqrt(3./4.));

    double Yi[9] = {
            -1., 1., 2.,
            0., 1., 2.,
            0., 1., 2.

    };

    double Yi2[9] = {0.};

//    for(int k=0;k<9;++k)
//    fprintf(stdout,"%f\n", Interpolate2d9_Fdx(NULL,Yi,GIPS2d9[k],0));

//    for(int k=0;k<20;++k)
//    {
//        double xi = 0.1*k - 1.0;
//        double eta = 0.0;
//        double local_xl[2] = {xi,eta};
//        int j0 = 3;
//        fprintf(stdout,"%f %f %f %f\n",xi,Ip2d9_X[j0+0][0](local_xl),Ip2d9_X[j0+1][0](local_xl),Ip2d9_X[j0+2][0](local_xl));
//    }


    return 1;
}

int test_Int2d9()
{
    double tYi[9] = {0.};
    for(int k=0;k<9;++k)
    {
        tYi[k] = 1.;
    }


    fprintf(stdout,"%s: %f\n\n",__func__, Int2d9(NULL,tYi));
    return 0;
}

int test_laplace2d9(double gamma)
{
    /*
     * when gamma = sqrt(3/4) = Oono-Puri nine point
     * gamma = sqrt(1.0/2) = Patra-Karttunen
     * check at wiki/Nine-point stencil
     */
    const double xlist[3] = {-1.0, 0.0, 1.0};
    const double ylist[3] = {0.0,1.0,2.0};
    double tXi[9][2] = {0.};
    for(int k=0;k<3;++k)
    {
        for(int j=0;j<3;++j)
        {
            tXi[3*j + k][0] = xlist[k];
            tXi[3*j + k][1] = ylist[j];
        }

    }
    double tYi[9] = {0.};
    double p[9] = {0.};
    double c[9] = {0.};
    double sum_p = 0.,sum_c = 0.;
    fprintf(stdout,"%s: gamma = %10.5f\n", __func__ ,gamma);
    for(int k=0;k<9;++k)
    {
        tYi[k] = 1.0;
        p[k] = Laplace2d9(tXi[0],tYi,gamma);
        c[k] = Laplace2d9cyl(tXi[0],tYi,gamma);
        tYi[k] = 0.;

        fprintf(stdout,"%7.4f;%7.4f ",p[k],c[k]);
        if(k%3 == 2) fprintf(stdout,"\n");
        sum_p += p[k];
        sum_c += c[k];

    }
    fprintf(stdout,"    sum = %f,%f\n\n",sum_p,sum_c);
    return 0.;
}

int test_spacing3(double ext, int eR, int dR, int N, int rpad)
{
    for(int k =0; k< N;++k)
    {
        if(k%10 == 0) fprintf(stdout,"\n");
        fprintf(stdout,"%10.5f  ", Spacing3(k+1,ext,eR,dR,N,rpad) - Spacing3(k,ext,eR,dR,N,rpad));
    }
    return 0;
}

int test_spacing2(double ext, int eL,int eR, int dL,int dR,int N)
{
    for(int k =0; k< N;++k)
    {
        if(k%10 == 0) fprintf(stdout,"\n");
        fprintf(stdout,"%10.5f", Spacing2(k+1,ext,eL,eR,dL,dR,N) - Spacing2(k,ext,eL,eR,dL,dR,N));
    }

    for(int k =0; k< N;++k)
    {
        if(k%10 == 0) fprintf(stdout,"\n");
        fprintf(stdout,"%10.5f", Spacing(k+1,ext,eL,eR,N) - Spacing(k,ext,eL,eR,N));
    }
    return 0;
}


double test_ablinear(double scale)
{
    double Xi[4];
    double Yi[4];
    double xl[2];
    double xg[2];
    double xlt[2];
    double xgt[2];
    for(int k=0;k<4;++k)
    {
        Xi[k] = (double) scale * rand()/RAND_MAX - scale/2.0;
        Yi[k] = (double) scale * rand()/RAND_MAX - scale/2.0;
        if(k < 2)
        {
            xl[k] = (double) 2.0 * rand()/RAND_MAX - 1.0;
        }
    }

//    Xi[1] = (double) (Xi[1] - Xi[0]) * rand()/RAND_MAX + Xi[0];
//    Xi[3] = (double) (Xi[1] - Xi[0]) * rand()/RAND_MAX + Xi[0];

    xg[0] = Interpolate2d(Xi,xl);
    xg[1] = Interpolate2d(Yi,xl);

    abilinear((double []){Xi[0],Xi[1],Xi[2],Xi[3],xg[0]},(double []){Yi[0],Yi[1],Yi[2],Yi[3],xg[1]},xlt);

    xgt[0] = Interpolate2d(Xi,xlt);
    xgt[1] = Interpolate2d(Yi,xlt);

    double xl_dis = vec_dis(xl,xlt,2);
    double xg_dis = vec_dis(xg,xgt,2);
    return xg_dis;
}

int test_compress(int argc, char ** argv)
{
    const int test_size = 4096*8;
    const int test_loop = 512;
    double compress_factor = 0.;

    srand(time(NULL));
    unsigned char * test_array = (unsigned char *) malloc(sizeof(unsigned char)*test_size);
    for(int k=0;k<test_loop;++k)
    {
        // generate data
        for(int j=0;j<test_size;++j)
        {
            test_array[j] = (unsigned char )(rand()%255);
        }

        int compressedlength = 0;
        unsigned char * compressedpointer = NULL;
        zlibcompress(test_array,test_size,&compressedpointer,&compressedlength);
        compress_factor += compressedlength;
        if(compressedpointer) free(compressedpointer);
    }

    free(test_array);
    fprintf(stdout,"test_size = %d, mean compressed ratio is %f\n",test_size,compress_factor/(test_size*test_loop));
    compare_binary_ascii(12*1024*1024);
    return 1;
}

double get_file_size_mb(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
        double fileSizeInBytes = (double)st.st_size;
        double fileSizeInMB = fileSizeInBytes / (1024 * 1024);
        return fileSizeInMB;
    }
    return -1.0;  // Return -1 on error
}

void compare_binary_ascii(int test_size)
{
    /*
     * compare file size written in bin/txt
     */
    if(test_size < 0) test_size = -test_size;
    const char fname_bin[] = "compare.bin";
    const char fname_txt[] = "compare.txt";
    const char fname_raw[] = "compare.raw";
    const int vnumber = RAND_MAX;

    float * test_array = malloc(sizeof(float)*test_size);
    srand(time(NULL));
    for(int j=0;j<test_size;++j)
    {
        test_array[j] = (float) (rand()%vnumber) /vnumber * 2e6 - 1.5e6;
    }

    FILE * fp_bin = fopen(fname_bin,"w");
    FILE * fp_txt = fopen(fname_txt,"w");
    FILE * fp_raw = fopen(fname_raw,"w");

    write_binary_array(test_size,test_array,fp_bin);
    write_ascii_array(test_size,6,test_array,fp_txt);
    fwrite(test_array, sizeof(float),test_size,fp_raw);
    fclose(fp_bin);
    fclose(fp_txt);
    fclose(fp_raw);
    free(test_array);

    double size_bin = get_file_size_mb(fname_bin);
    double size_txt = get_file_size_mb(fname_txt);
    double size_raw = get_file_size_mb(fname_raw);

    fprintf(stdout,"******test report******\n");
    fprintf(stdout,"%d float written in %s and %s\n", test_size,fname_bin,fname_txt);
    fprintf(stdout," %s is %f MB\n %s is %f MB\n %s is %f MB\n",fname_bin,size_bin,fname_txt,size_txt,fname_raw,size_raw);
    fprintf(stdout,"ratio : %s: %f %%\n",fname_bin,100.0 * size_bin/size_raw);
    fprintf(stdout,"ratio : %s: %f %%\n",fname_txt,100.0 * size_txt/size_raw);
    fprintf(stdout,"***********************\n");
}


int test_write_binary(int argc, char ** argv)
{
    /*
     * this function is used for test write_binary_array
     */

    FILE * fp = fopen("test_functions.txt","w");
    if(NULL == fp)
    {
        fprintf(stderr,"cannot create/open test files\n");
    }

//    FILE * fp = stdout;
    const int test_size = 256;
    const int min_test_size = 1;
    float * test_array = malloc(sizeof(float) * test_size);
    for(int k=0;k<test_size;++k)
    {
        test_array[k] = (float) rand() / RAND_MAX * 15.0f;
    }


    for(int k=test_size;k>min_test_size;k--)
    {
        fprintf(fp,"head of test for write binary array (%d)\n",k);
        fputs("\n",fp);
        write_binary_array(k,test_array,fp);
        fputs("\n",fp);
        fprintf(fp,"tail of test for write binary array (%d)\n",k);
    }


    free(test_array);
    if(fp != stdout) fclose(fp);
}


int test_sync(int argc, char ** argv)
{
    MPI_Init(&argc,&argv);
    field_var * local_fv = (field_var *)malloc(sizeof(field_var)*10);
    field_var * pfv = local_fv;
    field_opt * pfo = (field_opt *) malloc(sizeof(field_opt));
    init_field_opt_2d_default(pfo);

    pfo->lid2nty = node_type_2d_default_vertex;
    init_field_2d(&(field_2d){.npgx=2,.npgy=2,.xpad=3,.ypad=3,
                          .nx=10,.ny=15,.bytes_per_unit=2*sizeof(double)/sizeof(char),
                          .padding=vertex_pad},
                  pfv,pfo);
    sync_start(pfv);
    sync_complete(pfv);
    MPI_Finalize();
    return 0;
}


int test_tillpre()
{
    tillotson_table  tilltable;
    FILE * fp = fopen("../eos/aluminu.tillotson","r");
    LoadTillEOS(&tilltable,fp);

    double init_tem = 293.0;
    double init_pre = -0.25;
    double init_eng = 1006.3;
    double init_den = 8.65622e-5;

    double csd = 0., den = 0., eng = 0.;

    double pre_cold,csd_cold,pre_hot,csd_hot;

    for(int k=-20;k<=20;++k)
    {
        tillotson_calculate_state(&tilltable,init_den,init_eng,&init_tem,&init_pre,&csd);
        TillHot(&tilltable,init_eng/init_den,init_den,&pre_hot,&csd_hot);
        TillCold(&tilltable,init_eng/init_den,init_den,&pre_cold,&csd_cold);

        fprintf(stdout,"%20e,%20e,%20e\n",init_eng,init_den,init_pre);
    }
    UnAllocateTillEOS(&tilltable);
}

void test_aneossieden()
{
    struct ANEOSTable atable;
    FILE * fp = fopen("../eos/h2o_ice.aneos.sd","r");
    LoadANEOSSieDen(&atable,fp);
    struct StateReference aref;
    ANEOSSieDenInitStateRef(&atable,&aref);

    double init_pre = 1.0e5;
    double init_tem = 273.0;
    double init_eng = 0.0;
    double init_den = 0.;

    double test_pre = 0.;
    double test_eng = 0.;

    init_eng = ANEOSSieDenInterpolateTP(&atable,init_tem,init_pre,EOSENG);
    init_den = ANEOSSieDenInterpolateTP(&atable,init_tem,init_pre,EOSDEN);
    test_pre = ANEOSSieDenInterpolateTD(&atable,init_tem,init_den,EOSPRE);
    test_eng = ANEOSSieDenInterpolateTD(&atable,init_tem,init_den,EOSENG);

    double cal_pre, cal_csd, cal_tem;
    aneossieden_calculate_state(&atable,init_den,init_eng,&cal_tem,&cal_pre,&cal_csd);

    UnAllocateANEOS(&atable);
}

void test_aneossimple()
{
    struct ANEOSTable atable;
    FILE * fp = fopen("../eos/basalt_.aneos.td","r");
    LoadANEOS(&atable,fp);
    struct StateReference aref;
    ANEOSInitStateRef(&atable,&aref);

    double init_pre = 5.0e6;
    double init_tem = 293.0;
    double init_eng = 0.0;
    double init_den = 0.;
    double init_pty = 1.05;

    double test_pre = 0.;
    double test_eng = 0.;

    init_eng = ANEOSInterpolateTP(&atable,init_tem,init_pre*init_pty,EOSENG)/init_pty;
    init_den = ANEOSInterpolateTP(&atable,init_tem,init_pre*init_pty,EOSDEN)/init_pty;
    // test_pre = ANEOSInterpolateTD(&atable,init_tem,init_den*init_pty,EOSPRE);
    // test_eng = ANEOSInterpolateTD(&atable,init_tem,init_den*init_pty,EOSENG);
    // test_pre = ANEOSInterpolateTP(&atable,init_tem,init_pre,EOSDEN);

    double max_pty = init_pty, min_pty = 1.00;

    double cal_pre, cal_csd, cal_tem;
    aneos_calculate_state(&atable,init_den*init_pty,init_eng*init_pty,&cal_tem,&cal_pre,&cal_csd);
    cal_pre = cal_pre/init_pty;
    fprintf(stdout,"p=%e, dp=%e,ip=%e\n",init_pty,cal_pre - init_pre,cal_pre);
    for(double i_pty = max_pty; i_pty > min_pty; i_pty -= 1e-3 )
    {
        aneos_calculate_state(&atable,init_den*i_pty,init_eng*i_pty,&cal_tem,&cal_pre,&cal_csd);
        cal_pre = cal_pre/i_pty;
        fprintf(stdout,"p=%e, delta=%e,d=%e,\n",i_pty,cal_pre - init_pre,cal_pre/i_pty);
    }
    aneos_calculate_state(&atable,init_den,init_eng,&cal_tem,&cal_pre,&cal_csd);

    UnAllocateANEOS(&atable);
}