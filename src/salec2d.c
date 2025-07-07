//
// Created by huacheng on 11/3/22.
//

#include "sale2d.h"
#include "strength.h"
#include "ejecta.h"
#include "target_profile.h"
#include "write_vtk.h"

int main(int argc, char ** argv)
{
    MPI_Init(&argc,&argv);

    char inpfile[256];
    char eospath[256];

    if(argc == 1)
    {
        sprintf(inpfile,"../sale2d.inp");
        sprintf(eospath,"../eos");
    }
    else if(argc == 2)
    {
        sprintf(inpfile,"%s",argv[1]);
        sprintf(eospath,"../eos");
    } else if(argc == 3)
    {
        sprintf(inpfile,"%s",argv[1]);
        sprintf(eospath,"%s",argv[2]);
    } else
    {
        fprintf(stderr,"invalid args for salec2d\n");
        exit(0);
    }

    sale2d_var sale2d;
    mesh2d_info minfo;
    cycle_control ccl;

    double time_record[3] = {0.};
    // record time of start init
    time_record[0] = MPI_Wtime();

    InputFile * ifp = OpenInputFile(inpfile);

    init_cycle_control(ifp,&ccl);

    load_mesh_info(ifp,&minfo,eospath);
    init_sale_mat(&sale2d, &minfo);
    select_strength_model(&sale2d);
    init_sale2d_var(&sale2d,&minfo);
    init_sale2d_other(&sale2d,ifp);
    init_sale2d_ghost(&sale2d,&minfo);

    CloseInputFile(ifp);

    init_mesh_pos(&sale2d,&minfo);

    // debug entry for mpi/ it should be complied with -O0
#ifdef DEBUG_RANK
    volatile int debugvalue = 1;
    if(ccl.rank==DEBUG_RANK)
    {
        fprintf(stdout,"(rank = %d)reach the position mpi stucked\n",ccl.rank);
        fflush(stdout);
    }
    while (debugvalue &&  ccl.rank==DEBUG_RANK)
    {
        sleep(3);
    }
#endif

    init_material_distribution(&minfo, &sale2d);

    update_state(&sale2d,internal_section|send_section);

    sync_start(sale2d.e_den);
    sync_complete(sale2d.e_den);
    flush_ght(sale2d.e_pre);
    flush_ght(sale2d.e_ste);
    flush_ght(sale2d.m_eng);
    flush_ght(sale2d.m_pty);
    flush_ght(sale2d.e_cvs);

    double local_damp_time = record_avg_value(&sale2d);
    MPI_Allreduce(&(local_damp_time),&(minfo.damp_t),1,MPI_DOUBLE,MPI_MAX,sale2d.e_pre->comm);
    minfo.damp_t *= 2.0;
    sale2d.damp_time = Max(minfo.damp_t,sale2d.damp_time);

    init_v_field(&sale2d);
    clean_fv(sale2d.e_vel);
    clean_fv(sale2d.e_ste);
    clean_fv(sale2d.e_dvel);

    // record time of start cycle
    time_record[1] = MPI_Wtime();

    ccl.out_vts = 1;
    ccl.k_out_step = 1;

#ifdef UPDATE_TRACERS
    init_sale2d_tracer(&sale2d,&minfo);
    export_init_tracers(&sale2d,&minfo);
#endif

#if  defined(DEBUG_RANK) && defined(DEBUG_RANK_ONLY)
    debugvalue = 1;
    while (debugvalue)
    {
        sleep(3);
    }
#endif

    while(ccl.cycle_continue)
    {
        // dynamic track the crater profile
#ifdef UPDATE_EJECTA
        trace_profile(&sale2d);
        sync_index_start(sale2d.y_profile,MPI_MIN);
        sync_start(sale2d.e_ejt);
#endif
        /*
         * update the vof on vertex, according to values in elements
         * update_v_vof requires synchronous <m_vof>
         */
        sync_start(sale2d.m_vof);
        sync_start(sale2d.m_den);
        sync_start(sale2d.e_dvel);
        sync_start(sale2d.e_den);
        sync_start(sale2d.e_wpt);
        sync_start(sale2d.e_pre);

        update_v_vof(&sale2d, internal_section);
        update_v_vel_inc(&sale2d,internal_section);


        sync_complete(sale2d.m_vof);
        update_v_vof(&sale2d, update_section);


        /*
         * update velocity on vertex, ...
         * update_v_vel requires synchronous <e_den>, <e_vel>
         */

        sync_complete(sale2d.e_dvel);
        sync_complete(sale2d.e_den);
        sync_complete(sale2d.e_wpt);
        update_v_vel_inc(&sale2d,update_section);
        sync_restrict_vertexes(sale2d.v_vel);

        // calculate the velocity of element, store results in e_vel
        update_e_vel(&sale2d,send_section);
        sync_start(sale2d.e_vel);
        update_e_vel(&sale2d,internal_section);
        sync_complete(sale2d.e_vel);
        check_v_vel(&sale2d,update_section|internal_section);

        // calculate div and grad of velocity
        // artificial pressure is set
        update_e_gradvel(&sale2d, internal_section | send_section);
        sync_start(sale2d.e_q);


        // update time step
        ccl.local_dt = Min(ccl.max_dt,calculate_dt(&sale2d,ccl.dt));
        MPI_Allreduce(&(ccl.local_dt),&(ccl.dt),1,MPI_DOUBLE,MPI_MIN,sale2d.e_pre->comm);

        /*
         * calculate acfl vibration, only the block model is implemented
         * after the strength is determined by acfl, update the stress and strain
         */
        block_vibration(&sale2d,internal_section|send_section,ccl.dt,ccl.step_time);
        update_e_stress(&sale2d,internal_section|send_section,ccl.dt);
        sync_start(sale2d.e_tps);
        sync_start(sale2d.e_dam);
        sync_start(sale2d.e_vib);
        sync_start(sale2d.e_ste);

        /*
         *  calculate acceleration of vertex,
         *  calculate_acc requires synchronous <e_pre>, <e_stress>, <e_den>
         */

        calculate_acc(minfo.damp_type,&sale2d,internal_section,ccl.dt);
        sync_complete(sale2d.e_pre);
        sync_complete(sale2d.e_ste);
        sync_complete(sale2d.e_q);
        calculate_acc(minfo.damp_type,&sale2d,update_section,ccl.dt);
        sync_restrict_vertexes(sale2d.v_acc);

        // calculate displacement in every vertex, and energy change in every element
        calculate_v_dis(&sale2d, update_section | internal_section, ccl.dt);

        // after calculate_vis, the anc acceleration has been calculated
        sync_start(sale2d.v_vel);

        sync_restrict_vertexes(sale2d.v_dis);
        update_energy(&sale2d,internal_section|send_section,ccl.dt);
        sync_start(sale2d.m_eng);
        if(sale2d.porosity_model) sync_start(sale2d.m_pty);
        if(sale2d.porosity_model) sync_start(sale2d.e_cvs);


        /*
         * update tracer info
         */
#ifdef UPDATE_TRACERS
        fl_tracer_update(sale2d.tracer_list,ccl.dt);
        fl_sync_info_start(sale2d.tracer_list);
        fl_sync_info_complete(sale2d.tracer_list);
        fl_sync_node_start(sale2d.tracer_list);
#endif
        /*
         * calculate the flux across boundary
         * calculate_flux_average requires synchronous <e_vel>, <m_vof>, <m_den>, <m_eng>
         */

        clean_flux(&sale2d);
        calculate_flux(&sale2d,internal_section|send_section);

        sync_complete(sale2d.m_den);
        sync_complete(sale2d.e_tps);
        sync_complete(sale2d.e_dam);
        sync_complete(sale2d.e_vib);
        sync_complete(sale2d.m_eng);
        if(sale2d.porosity_model) sync_complete(sale2d.m_pty);
        if(sale2d.porosity_model) sync_complete(sale2d.e_cvs);
#ifdef UPDATE_EJECTA
        sync_complete(sale2d.e_ejt);
        sync_index_complete(sale2d.y_profile);
#endif
        // some ghost information need inserted, sync_ghost is called automatically in sync_complete
        calculate_flux(&sale2d,receive_section);

#ifdef UPDATE_TRACERS
        fl_sync_node_complete(sale2d.tracer_list); // complete of tracers updating
        update_strength_condition(&sale2d);
#endif
        if(ccl.out_vts)
        {
            sale2d_write_vtk(&sale2d, &minfo, &ccl);
            sale2d_write_vtp(&sale2d, &minfo, &ccl);
        }

        sale2d_write_log(stdout,&ccl);

        /*
         * just update state of internal/send section elements,
         * elements in receive section will receive values in next cycle through sync_*
         */

        advect_flux(minfo.damp_type,&sale2d,internal_section|send_section,ccl.dt);
        damp_on_cell(minfo.damp_type,&sale2d,ccl.dt,ccl.step_time);
#ifdef UPDATE_EJECTA
        ejecta_flow(&sale2d,internal_section|send_section,ccl.dt); // decrease velocity/momentum in cells if necessary
#endif
        update_state(&sale2d,internal_section|send_section);

        /*
         * check the stop condition and other control tags
         */
        update_cycle_control(&ccl,stdout);

        /*
         * velocity on vertex is updated through its increment (flux)
         * and local acceleration (calculate_v_dis). before add the effect flux
         * synchronize it through sync_
         */
        sync_complete(sale2d.v_vel);
    }

    /*
     * many bugs in clean_sale_var()
     * left free of memory to system
     */
    // clean_sale_var(&sale2d);

    // record time of end cycle
    time_record[2] = MPI_Wtime();
    // report time costed
    if(ccl.rank == 0)
    {
        fprintf(stdout,"--*****************************--\n");
        fprintf(stdout,"  initialized  -- %10.3f s \n",  time_record[1] - time_record[0]);
        fprintf(stdout,"  all  cycles  -- %10.3f s \n",  time_record[2] - time_record[1]);
        fprintf(stdout,"  1k   cycles  -- %10.3f s \n",  (time_record[2] - time_record[1])/ccl.i_step*1000.0);
        fprintf(stdout,"  100k cycles  -- %10.3f h \n",  (time_record[2] - time_record[1])/ccl.i_step*1.0e5/3600.);
        fprintf(stdout,"--*****************************--\n");
    }

    MPI_Finalize();
    return 0;
}

