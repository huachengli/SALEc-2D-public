//
// Created by huacheng on 5/30/23.
// use functions in io/output.*
// write grid/tracer information into vts/vtp file
//

#ifndef SALE_REBUILD_WRITE_VTK_H
#define SALE_REBUILD_WRITE_VTK_H

#include "tracer.h"
#include "sale2d.h"

void sale2d_write_vtp(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl);
void sale2d_write_vtk(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl);
void sale2d_write_pvtk(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl);
void sale2d_write_pvtp(sale2d_var * _sale, mesh2d_info * _minfo, cycle_control * _ccl);


#endif //SALE_REBUILD_WRITE_VTK_H
