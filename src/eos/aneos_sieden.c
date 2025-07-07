//
// Created by huacheng on 6/7/23.
// add some function to process ANEOSTable at SieDen
// this type use struct ANEOSTable to store data same as ANEOSTable at TmpDen
//

#include "eos_state.h"
#include "interpolate.h"

// the allocate/free is same
void LoadANEOSSieDen(struct ANEOSTable *_t, FILE *fp)
{
    // the comment line should be processed before this function
    // process them here is not a right solution
    Over(fp,3);
    if(2!=fscanf(fp, "%d %d", &_t->nDen, &_t->nTem))
    {
        fprintf(stdout,"Can not get number of indexes for density or temperature!\n");
        exit(0);
    }
    AllocateANEOS(_t);

    if(2!=fscanf(fp, "%le %le", &_t->Pnorm, &_t->Tnorm))
    {
        fprintf(stdout, "Can not get the standary state from EOS file!\n");
        exit(0);
    }


    for (int i = 0; i < _t->nTem; i++)
    {
        if(1!=fscanf(fp, "%le", &_t->xTem[i]))
        {
            fprintf(stdout,"Can not get table value!!\n");
            exit(0);
        }
    }

    for (int i = 0; i < _t->nDen; i++)
    {
        if(1!=fscanf(fp, "%le", &_t->yDen[i]))
            fprintf(stdout, "Can not get index values of density!\n");
    }


    for (int j = 0; j < _t->nDen; j++)
    {
        for (int i = 0; i < _t->nTem; i++)
        {
            if(3!=fscanf(fp, "%le %le %le\n", &_t->Data[i][j][0],&_t->Data[i][j][1],&_t->Data[i][j][2]))
            {
                fprintf(stdout, "Can not get table value from EOS table!\n");
                exit(0);
            }
        }
    }

    // at the last red entropy into memory, this section data is not included in TmpDen type
    double dummy_entropy,dummy_strength;
    for (int j = 0; j < _t->nDen; j++)
    {
        for (int i = 0; i < _t->nTem; i++)
        {
            if(2!=fscanf(fp, "%le %le\n", &dummy_entropy,&dummy_strength))
            {
                fprintf(stdout, "Can not get table value from EOS table!\n");
                exit(0);
            }
        }
    }

    fclose(fp);
}


double ANEOSSieDenInterpolateTD(struct ANEOSTable *_t,double tTem, double tDen, int DataId)
{
    // the first step locate the density relative pos
    int DI = bisearch(tDen,_t->yDen,_t->nDen);
    int DJ = DI + 1;

    double Dxr = (tDen - _t->yDen[DI])/(_t->yDen[DJ] - _t->yDen[DI]);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;

    int EI = 0;
    int EJ = _t->nTem - 1;

    double IRT = _t->Data[EI][DI][0]*Dxl + _t->Data[EI][DJ][0]*Dxr;
    double JRT = _t->Data[EJ][DI][0]*Dxl + _t->Data[EJ][DJ][0]*Dxr;

    while (EI + 1 < EJ)
    {
        int EM = (EI + EJ)/2;
        double MRT = _t->Data[EM][DI][0]*Dxl + _t->Data[EM][DJ][0]*Dxr;
        if(MRT > tTem)
        {
            EJ = EM;
            JRT = MRT;
        } else
        {
            EI = EM;
            IRT = MRT;
        }
    }

    double Exr = (tTem - IRT)/(JRT - IRT);
    double Exl = Wind(1.0 - Exr,0.0,1.0);

    double ret_val = 0.;
    switch (DataId)
    {
        case EOSENG:
            ret_val = _t->xTem[EI]*Exl + _t->xTem[EJ]*Exr;
            break;
        case EOSCSD:
        case EOSPRE:
            ret_val = (_t->Data[EI][DI][DataId] * Dxl + _t->Data[EI][DJ][DataId] * Dxr) * Exl
                      + (_t->Data[EJ][DI][DataId] * Dxl + _t->Data[EJ][DJ][DataId] * Dxr) * Exr;
            break;
        default:
            fprintf(stderr,"Invalid data id in %s\n",__func__);
            exit(0);
            break;
    }

    if(DataId == EOSENG)
    {
        ret_val *= tDen;
    }

    return ret_val;
}


double ANEOSSieDenInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId)
{
    // in the SieDen type ANEOS
    // x is <Energy>， y is <Density>
    // data represent <temperature, pressure, csound>
    int TI = 0;
    int DI = 0;
    double xl[2] = {0.};
    int sl_state = 0;

    for(TI=0;TI<_t->nTem-1;++TI)
    {
        for(DI=0;DI<_t->nDen-1;++DI)
        {
            double pre_in[5] = {_t->Data[TI][DI][EOSPRE],_t->Data[TI+1][DI][EOSPRE],
                                _t->Data[TI+1][DI+1][EOSPRE],_t->Data[TI][DI+1][EOSPRE],tPre};
            double tem_in[5] = {_t->Data[TI][DI][EOSENG],_t->Data[TI+1][DI][EOSENG],
                                _t->Data[TI+1][DI+1][EOSENG],_t->Data[TI][DI+1][EOSENG],tTem};

            int pre_inter = vec_cmp(pre_in,tPre,4);
            int tem_inter = vec_cmp(tem_in,tTem,4);

            if(pre_inter <=0 || pre_inter>= 4) continue;
            if(tem_inter <=0 || tem_inter>= 4) continue;
            sl_state = abilinear(tem_in,pre_in,xl);
            if(sl_state != 0) break;
        }
        if(sl_state != 0) break;
    }

    xl[0] = (xl[0] + 1.0)/2.0;
    xl[1] = (xl[1] + 1.0)/2.0;

    double ret_val = 0.;
    switch(DataId)
    {
        case EOSENG:
            ret_val = _t->xTem[TI]*(1.0 - xl[0]) + _t->xTem[TI+1]*xl[0];
            break;
        case EOSDEN:
            ret_val = _t->yDen[DI]*(1.0 - xl[1]) + _t->yDen[DI+1]*xl[1];
            break;
        case EOSCSD:
            ret_val = _t->Data[TI][DI][EOSCSD]*(1.0 - xl[0])*(1.0 - xl[1]) + _t->Data[TI+1][DI][EOSCSD]*xl[0]*(1.0 - xl[1])
                    + _t->Data[TI+1][DI+1][EOSCSD]*xl[0]*xl[0] + _t->Data[TI][DI+1][EOSCSD]*(1.0 - xl[0])*xl[1];
            break;
        case EOSPRE:
            ret_val = tPre;
        case EOSTEM:
            ret_val = tTem;
        default:
            fprintf(stderr,"Invalid data id in %s\n",__func__);
            exit(0);
            break;
    }

    if(DataId == EOSENG)
    {
        ret_val *= _t->yDen[DI]*(1.0 - xl[1]) + _t->yDen[DI+1]*xl[1];
    }

    return ret_val;
}

void ANEOSSieDenInitStateRef(struct ANEOSTable *_t, struct StateReference *_s)
{
    /*
     * Init the state reference of material
     * (T norm, P norm) -> (D norm)
     */

    int TI = 0;
    int DI = 0;
    double xl[2] = {0.};
    int sl_state = 0;

    for(TI=0;TI<_t->nTem-1;++TI)
    {
        for(DI=0;DI<_t->nDen-1;++DI)
        {
            double pre_in[5] = {_t->Data[TI][DI][EOSPRE],_t->Data[TI+1][DI][EOSPRE],
                                _t->Data[TI+1][DI+1][EOSPRE],_t->Data[TI][DI+1][EOSPRE],_t->Pnorm};
            double tem_in[5] = {_t->Data[TI][DI][EOSENG],_t->Data[TI+1][DI][EOSENG],
                                _t->Data[TI+1][DI+1][EOSENG],_t->Data[TI][DI+1][EOSENG],_t->Tnorm};

            int pre_inter = vec_cmp(pre_in,_t->Pnorm,4);
            int tem_inter = vec_cmp(tem_in,_t->Tnorm,4);

            if(pre_inter <=0 || pre_inter>= 4) continue;
            if(tem_inter <=0 || tem_inter>= 4) continue;
            sl_state = abilinear(tem_in,pre_in,xl);
            if(sl_state != 0) break;
        }
        if(sl_state != 0) break;
    }

    int TJ = TI+1, DJ=DI+1;

    _t->Dnorm = _t->yDen[DI]*(1.0 - xl[1])/2.0 + _t->yDen[DJ]*(1.0 + xl[1])/2.0;
    double Cs = _t->Data[TI][DI][EOSCSD]*IpB[0](xl) + _t->Data[TI+1][DI][EOSCSD]*IpB[1](xl);
                + _t->Data[TJ][DJ][EOSCSD]*IpB[2](xl) + _t->Data[TI][DI+1][EOSCSD]*IpB[3](xl);

    // Get the density of vapor and melt accord to the Dnorm
    _s->MeltDen =  0.85 * _t->Dnorm;
    _s->VaporDen = 0.04 * _t->Dnorm;
    _s->Viscosity = 0.0;

    // calculate dpdt and dedt
    double dpde = IpB_X[0][0](xl)*_t->Data[TI][DI][EOSPRE]
                  + IpB_X[1][0](xl)*_t->Data[TJ][DI][EOSPRE]
                  + IpB_X[2][0](xl)*_t->Data[TJ][DJ][EOSPRE]
                  + IpB_X[3][0](xl)*_t->Data[TI][DJ][EOSPRE];

    double dtde = IpB_X[0][0](xl)*_t->Data[TI][DI][EOSENG]
                  + IpB_X[1][0](xl)*_t->Data[TJ][DI][EOSENG]
                  + IpB_X[2][0](xl)*_t->Data[TJ][DJ][EOSENG]
                  + IpB_X[3][0](xl)*_t->Data[TI][DJ][EOSENG];

    dpde *= 2.0/(_t->xTem[TJ] - _t->xTem[TI]);
    dtde *= 2.0/(_t->xTem[TJ] - _t->xTem[TI]);

    _s->den0 = _t->Dnorm;
    _s->thermexp = (dpde/dtde)/(Cs*Cs*_t->Dnorm);
    _s->c_heat = 1.0/dtde;
    _s->csd0 = Cs;
    _s->gamfrac = _s->thermexp/_s->c_heat;
}

double aneossieden_calculate_state(aneos_table * _t,double den,double eng,double * tem,double * pre,double * csd)
{
    double tDen = den;
    double tEng = eng/den;

    int DI = bisearch(tDen,_t->yDen,_t->nDen);
    int DJ = DI + 1;

    int TI = bisearch(tEng,_t->xTem,_t->nTem);
    int TJ = TI + 1;

    double Dxr = (tDen - _t->yDen[DI])/(_t->yDen[DJ] - _t->yDen[DI]);
    double Dxl = Wind(1.0 - Dxr,0.0,1.0);

    double Exr = (tEng - _t->xTem[TI])/(_t->xTem[TJ] - _t->xTem[TI]);
    double Exl = Wind(1.0 - Exr,0.,1.0);

    *tem = (_t->Data[TI][DI][EOSENG] * Dxl + _t->Data[TI][DJ][EOSENG] * Dxr) * Exl
           + (_t->Data[TJ][DI][EOSENG] * Dxl + _t->Data[TJ][DJ][EOSENG] * Dxr) * Exr;

    *pre = (_t->Data[TI][DI][EOSPRE] * Dxl + _t->Data[TI][DJ][EOSPRE] * Dxr) * Exl
           + (_t->Data[TJ][DI][EOSPRE] * Dxl + _t->Data[TJ][DJ][EOSPRE] * Dxr) * Exr;

    *csd = (_t->Data[TI][DI][EOSCSD] * Dxl + _t->Data[TI][DJ][EOSCSD] * Dxr) * Exl
           + (_t->Data[TJ][DI][EOSCSD] * Dxl + _t->Data[TJ][DJ][EOSCSD] * Dxr) * Exr;

}