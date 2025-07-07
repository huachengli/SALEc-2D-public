//
// Created by huach on 12/7/2023.
//

#ifndef SALE_REBUILD_EJECTA_H
#define SALE_REBUILD_EJECTA_H

#include "sale2d.h"

void trace_profile(sale2d_var * _sale);
void ejecta_flow(sale2d_var * _sale,node_type _nt, double dt);

#endif //SALE_REBUILD_EJECTA_H
