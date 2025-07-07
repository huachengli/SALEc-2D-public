//
// Created by huach on 1/6/2023.
//

#ifndef SALE_REBUILD_SALE2D_CONSTANT_H
#define SALE_REBUILD_SALE2D_CONSTANT_H


#define MAXMAT 10
#define MAXPHI 42
#define TOLVOF 1.0e-6
#define TOLRHO 1.0e-4
#define VACUUM 0
#define TOLQPRE 1e-4
#define VELMIN_CUTOFF 1  // delete velocity less than VELMIN
#define ACCMIN_CUTOFF 1
#define VACUUM_COMPRESSED 0
#define DTF 0.2
#define CISSAMK1 1.0
#define WIKINS_RORARTION 1
#define LAPLICAN_GAMMA 0.5
#define UPDATE_TRACERS 1
#define UPDATE_EJECTA 1
#define UNIFORM_XGRID 1
#define FLUXSTRATEGY mixed_strategy
#define AVOFALG vofMSTACS
/*
 * select from vofMSTACS,vofCICSAM,vofSTACS,vofMHRIC
 */
#define STRENGTHMOD avg_strength
//#define MEANVERTEXVOL 1
#define EXCLUDE_PML_POROSITY 1
#endif //SALE_REBUILD_SALE2D_CONSTANT_H
