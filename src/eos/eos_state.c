//
// Created by huacheng on 1/11/22.
//

#include "eos_state.h"
#include "interpolate.h"
inline double Wind(double x, double limitL, double limitR)
{
    if(x < limitL)
        return limitL;
    else if(x > limitR)
        return limitR;
    else
        return x;
}

inline void Over(FILE * fp,int n)
{
    /*
     * To process the head of ANEOS table
     * The first version is only able to process some comment line
     * in the head.
     */
    char buf[1024];
    for(int k=0;k<n;k++)
    {
        if(NULL==fgets(buf,sizeof(buf),fp))
            fprintf(stdout,"error in fgets(). skip over some lines is not valid.\n");
    }
}

inline double Max(double a, double b)
{
    if (a > b)
        return a;
    else
        return b;
}

inline double Min(double a, double b)
{
    if (a < b)
        return a;
    else
        return b;
}

inline double Sqrt_trunc(double x)
{
    if(x < 1.0e-2)
        return -sqrt(fabs(x));
    else
        return sqrt(x);
}


void AllocateANEOS(struct ANEOSTable *_t)
{
    _t->yDen = (double *) malloc(_t->nDen * sizeof(double));
    _t->xTem = (double *) malloc(_t->nTem * sizeof(double));
    _t->Data = (double ***) malloc(_t->nTem * sizeof(double **));

    for (int k = 0; k < _t->nTem; k++)
    {
        _t->Data[k] = (double **) malloc(_t->nDen * sizeof(double *));
        for (int i = 0; i < _t->nDen; i++)
        {
            _t->Data[k][i] = (double *) malloc(3 * sizeof(double));
        }
    }
}

void UnAllocateANEOS(struct ANEOSTable *_t)
{

    free(_t->xTem);
    free(_t->yDen);
    for (int k = 0; k < 3; k++) // need to be fixed
    {
        for (int i = 0; i < _t->nTem; i++)
        {
            free(_t->Data[k][i]);
        }
        free(_t->Data[k]);
    }
    free(_t->Data);
}

void LoadANEOS(struct ANEOSTable *_t, FILE *fp)
{
    // the comment line should be processed before this function
    // process them here is not a right solution
    Over(fp,3);
    if(2!=fscanf(fp, "%d %d", &_t->nTem, &_t->nDen))
    {
        fprintf(stdout,"Can not get number of indexes for density or temperature!\n");
        exit(0);
    }
    AllocateANEOS(_t);

    if(3!=fscanf(fp, "%le %le %le", &_t->Tnorm, &_t->Dnorm, &_t->Pnorm))
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

    for (int i = 0; i < _t->nTem; i++)
    {
        for (int j = 0; j < _t->nDen; j++)
        {
            if(3!=fscanf(fp, "%le %le %le\n", &_t->Data[i][j][0],&_t->Data[i][j][1],&_t->Data[i][j][2]))
            {
                fprintf(stdout, "Can not get table value from EOS table!\n");
                exit(0);
            }
        }
    }

    for (int i = 0; i < _t->nDen; i++)
    {
        if(1!=fscanf(fp, "%le", &_t->yDen[i]))
            fprintf(stdout, "Can not get index values of density!\n");
    }

    // in the raw table , the first column is energy
    // but in calculation energy*density is used more usually
    // this will make a difference to the linear interpolation of energy in
    // (temperature,density). (a linear interpolation of energy*density)
    for (int i = 0; i < _t->nTem; i++)
    {
        for (int j = 0; j < _t->nDen; j++)
        {
            _t->Data[i][j][0] *= _t->yDen[j];
        }
    }

    fclose(fp);

}

void ANEOSInitStateRef(struct ANEOSTable *_t, struct StateReference *_s)
{
    /*
     * change the Dnorm accord to EOS
     * make (Tnorm,Dnorm,Pnorm) is consistent with ANEOSTable
     * this segment is copied from ANEOSCalL3t
     */
    int TI = 0;
    int TJ = _t->nTem;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        if (_t->xTem[TM] > _t->Tnorm)
            TJ = TM;
        else
            TI = TM;
    }

    double Txr = (_t->Tnorm - _t->xTem[TI]) / (_t->xTem[TJ] - _t->xTem[TI]);
    double Txl = 1.0 - Txr;

    int DI = 0;
    int DJ = _t->nDen - 1;
    double IRP = _t->Data[TI][DI][1] * Txl + _t->Data[TJ][DI][1] * Txr;
    double JRP = _t->Data[TI][DJ][1] * Txl + _t->Data[TJ][DJ][1] * Txr;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        double MRP = _t->Data[TI][DM][1] * Txl + _t->Data[TJ][DM][1] * Txr;
        if (MRP > _t->Pnorm)
        {
            DJ = DM;
            JRP = MRP;
        }
        else
        {
            DI = DM;
            IRP = MRP;
        }
    }

    double Dxr = (_t->Pnorm - IRP) / (JRP - IRP);
    double Dxl = 1.0 - Dxr;
    _t->Dnorm = _t->yDen[DI] * Dxl + _t->yDen[DJ] * Dxr;

    // Get the density of vapor and melt accord to the Dnorm
    _s->MeltDen =  0.85 * _t->Dnorm;
    _s->VaporDen = 0.04 * _t->Dnorm;
    _s->Viscosity = 0.0;

    double lcoord[2] = {Txr*2.0 - 1.0,Dxr*2.0 - 1.0};
    double dpdt = IpB_X[0][0](lcoord)*_t->Data[TI][DI][1]
                  + IpB_X[1][0](lcoord)*_t->Data[TJ][DI][1]
                  + IpB_X[2][0](lcoord)*_t->Data[TJ][DJ][1]
                  + IpB_X[3][0](lcoord)*_t->Data[TI][DJ][1];

    dpdt *= 2.0/(_t->xTem[TJ] - _t->xTem[TI]);

    double Cs = (_t->Data[TI][DI][2] * Dxl + _t->Data[TI][DJ][2] * Dxr) * Txl
                + (_t->Data[TJ][DI][2] * Dxl + _t->Data[TJ][DJ][2] * Dxr) * Txr;

    _s->thermexp = dpdt/(Cs*Cs*_t->Dnorm);
    double dedt = IpB_X[0][0](lcoord)*_t->Data[TI][DI][0]/_t->yDen[DI]
                  + IpB_X[3][0](lcoord)*_t->Data[TI][DJ][0]/_t->yDen[DJ]
                  + IpB_X[1][0](lcoord)*_t->Data[TJ][DI][0]/_t->yDen[DI]
                  + IpB_X[2][0](lcoord)*_t->Data[TJ][DJ][0]/_t->yDen[DJ];

    dedt *= 2.0/(_t->xTem[TJ] - _t->xTem[TI]);
    _s->c_heat = dedt;
    _s->den0 = _t->Dnorm;
    _s->csd0 = Cs;
    _s->gamfrac = _s->thermexp/_s->c_heat;
}

double ANEOSInterpolateTP(struct ANEOSTable *_t, double tTem, double tPre, int DataId)
{
    /*
     * get data using (Ten,Pre)
     * copy from ANEOSInitStateRef
     */

    int TI = 0;
    int TJ = _t->nTem;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        if (_t->xTem[TM] > tTem)
            TJ = TM;
        else
            TI = TM;
    }

    double Txr = (tTem - _t->xTem[TI]) / (_t->xTem[TJ] - _t->xTem[TI]);
    double Txl = 1.0 - Txr;

    int DI = 0;
    int DJ = _t->nDen - 1;
    double IRP = _t->Data[TI][DI][1] * Txl + _t->Data[TJ][DI][1] * Txr;
    double JRP = _t->Data[TI][DJ][1] * Txl + _t->Data[TJ][DJ][1] * Txr;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        double MRP = _t->Data[TI][DM][1] * Txl + _t->Data[TJ][DM][1] * Txr;
        if (MRP > tPre)
        {
            DJ = DM;
            JRP = MRP;
        }
        else
        {
            DI = DM;
            IRP = MRP;
        }
    }

    double Dxr = (tPre - IRP) / (JRP - IRP);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;
    if(DataId == -1)
        return _t->yDen[DI] * Dxl + _t->yDen[DJ] * Dxr;
    else
        return (_t->Data[TI][DI][DataId] * Dxl + _t->Data[TI][DJ][DataId] * Dxr) * Txl
               + (_t->Data[TJ][DI][DataId] * Dxl + _t->Data[TJ][DJ][DataId] * Dxr) * Txr;

}

double ANEOSInterpolateTD(struct ANEOSTable *_t,double tTem, double tDen, int DataId)
{
    /*
     * Interpolate the data using (Tem,Den)
     */

    int TI = 0;
    int TJ = _t->nTem - 1;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        if (_t->xTem[TM] > tTem)
            TJ = TM;
        else
            TI = TM;
    }

    double Txr = (tTem - _t->xTem[TI]) / (_t->xTem[TJ] - _t->xTem[TI]);
    Txr = Wind(Txr,0.0,1.0);
    double Txl = 1.0 - Txr;

    int DI = 0;
    int DJ = _t->nDen - 1;
    while(DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        if (_t->yDen[DM] > tDen)
            DJ = DM;
        else
            DI = DM;
    }

    double Dxr = (tDen - _t->yDen[DI])/(_t->yDen[DJ] - _t->yDen[DI]);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;

    if(DataId <3 && DataId >=0)
    {
        return (_t->Data[TI][DI][DataId] * Dxl + _t->Data[TI][DJ][DataId] * Dxr) * Txl
               + (_t->Data[TJ][DI][DataId] * Dxl + _t->Data[TJ][DJ][DataId] * Dxr) * Txr;
    }

}


double TillCold(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{
    /*
     * Tillotson equation of state
     * cold compressed states, Eta >= 1 or EngIn <= Eiv
     */
    double TillEta = RhoIn/_t->TLRho0;
    double TillMu  = TillEta - 1.0;
    double TillPhi0 = EngIn/(_t->TLE0*TillEta*TillEta);
    double TillPhi1 = 1.0/(TillPhi0 + 1.0);
    double temp_B = TillEta < 1.0 ? _t->TLA : _t->TLB;

    *Pres   = (_t->TLa + _t->TLb*TillPhi1)*EngIn*RhoIn + _t->TLA*TillMu + temp_B*TillMu*TillMu;

    double CsSquare1 = (1.0/_t->TLRho0)*(_t->TLA + 2.0*temp_B*TillMu) + EngIn*(_t->TLa + _t->TLb*TillPhi1*(3.0 - 2.0*TillPhi1));
    double CsSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1)*(*Pres)/RhoIn;
    *Cs = Sqrt_trunc(CsSquare1 + CsSquare2);

    return *Cs;
}

double TillHot(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{
    /*
     * Tillotson equation of state
     * expanded states, Eta <= 1 && EngIn >= Ecv
     */
    double TillEta = RhoIn/_t->TLRho0;
    double TillMu  = TillEta - 1.0;
    double TillXi  = _t->TLRho0/RhoIn - 1.0;
    double TillPhi0 = EngIn/(_t->TLE0*TillEta*TillEta);
    double TillPhi1 = 1.0/(TillPhi0 + 1.0);

    double Exp_A2   = 0.0;
    if(_t->TLAlpha*TillXi*TillXi < 100.0)
        Exp_A2 = exp(-1.0*_t->TLAlpha*TillXi*TillXi);
    double Exp_B1  = 0.0;
    if(_t->TLBeta*TillXi < 100.0)
        Exp_B1 = exp(-1.0*_t->TLBeta*TillXi);
    *Pres = _t->TLa*RhoIn*EngIn + (_t->TLb*EngIn*RhoIn*TillPhi1 + _t->TLA*TillMu*Exp_B1)*Exp_A2;

    double CsSquare1 = _t->TLa*EngIn + _t->TLb*EngIn*TillPhi1*(3.0 + 2.0*_t->TLAlpha*TillXi/TillEta - 2.0*TillPhi1)*Exp_A2
                       + _t->TLA/_t->TLRho0*(1.0 - (TillXi/TillEta)*(_t->TLBeta+ 2.0*_t->TLAlpha*TillXi))*Exp_A2*Exp_B1;
    double CsSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1*Exp_A2)*(*Pres)/RhoIn;
    *Cs = Sqrt_trunc(CsSquare1 + CsSquare2);

    return *Cs;
}

double TillTran(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{
    /*
     * Tillotson equation of state
     * transition between Cold and Hot
     */

    double TillETran = _t->TLEcv - _t->TLEiv;
    double TilldEr   = EngIn - _t->TLEiv;
    double TilldEl   = _t->TLEcv - EngIn;

    double TillEta = RhoIn/_t->TLRho0;
    double TillMu  = TillEta - 1.0;
    double TillPhi0 = EngIn/(_t->TLE0*TillEta*TillEta);
    double TillPhi1 = 1.0/(TillPhi0 + 1.0);
    double TillXi  = _t->TLRho0/RhoIn - 1.0;

    double temp_B = TillEta < 1.0 ? _t->TLA : _t->TLB;

    double PresCold = (_t->TLa + _t->TLb*TillPhi1)*EngIn*RhoIn + _t->TLA*TillMu + temp_B*TillMu*TillMu;

    double Exp_A2   = 0.0;
    if(_t->TLAlpha*TillXi*TillXi < 100.0)
    {
        Exp_A2 = exp(-1.0*_t->TLAlpha*TillXi*TillXi);
    }
    double Exp_B1  = 0.0;
    if(_t->TLBeta*TillXi < 100.0)
    {
        Exp_B1 = exp(-1.0*_t->TLBeta*TillXi);
    }

    double PresHot = _t->TLa*RhoIn*EngIn + (_t->TLb*EngIn*RhoIn*TillPhi1 + _t->TLA*TillMu*Exp_B1)*Exp_A2;

    *Pres = (TilldEl*PresCold + TilldEr*PresHot)/TillETran;

    double CsColdSquare1 = (1.0/_t->TLRho0)*(_t->TLA + 2.0*temp_B*TillMu) + EngIn*(_t->TLa + _t->TLb*TillPhi1*(3.0 - 2.0*TillPhi1));
    double CsColdSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1)*(*Pres)/RhoIn;

    double CsHotSquare1 = _t->TLa*EngIn + _t->TLb*EngIn*TillPhi1*(3.0 + 2.0*_t->TLAlpha*TillXi/TillEta - 2.0*TillPhi1)*Exp_A2
                          + _t->TLA/_t->TLRho0*(1.0 - (TillXi/TillEta)*(_t->TLBeta+ 2.0*_t->TLAlpha*TillXi))*Exp_A2*Exp_B1;
    double CsHotSquare2 = (_t->TLa + _t->TLb*TillPhi1*TillPhi1*Exp_A2)*(*Pres)/RhoIn;

    double Cs2Cold = CsColdSquare1 + CsColdSquare2;
    double Cs2Hot  = CsHotSquare1  + CsHotSquare2 ;

    double CsSquare1 = (TilldEl*Cs2Cold + TilldEr*Cs2Hot)/TillETran;
    double CsSquare2 = (PresHot - PresCold)/TillETran * (*Pres)/(RhoIn * RhoIn);
    *Cs = Sqrt_trunc(CsSquare1 + CsSquare2);

    return *Cs;
}


double TillPres(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{

    if((RhoIn >= _t->TLRho0) || (EngIn < _t->TLEiv))
        TillCold(_t,EngIn,RhoIn,Pres,Cs);
    else if(EngIn > _t->TLEcv)
        TillHot(_t,EngIn,RhoIn,Pres,Cs);
    else
        TillTran(_t,EngIn,RhoIn,Pres,Cs);

    if(*Pres <= 1.0e-2)
        *Cs = Sqrt_trunc(_t->TLA/_t->TLRho0);
    return *Cs;
}

double TillPres_split(const struct TillotsonTable * _t, double EngIn, double RhoIn, double * Pres, double * Cs)
{

    if((RhoIn >= _t->TLRho0) || (EngIn < _t->TLEiv))
        TillCold(_t,EngIn,RhoIn,Pres,Cs);
    else if(EngIn > _t->TLEcv)
        TillHot(_t,EngIn,RhoIn,Pres,Cs);
    else
        TillTran(_t,EngIn,RhoIn,Pres,Cs);

    if(*Pres <= 1.0e-2)
        *Cs = Sqrt_trunc(_t->TLA/_t->TLRho0);
    return *Cs;
}


double TilldEdr(const struct TillotsonTable * _t, double EngIn, double RhoIn)
{
    double tCs,tPres;
    TillPres_split(_t,EngIn,RhoIn,&tPres,&tCs);
    return tPres/(RhoIn*RhoIn);
}

void TillColdEnergy(struct TillotsonTable * _t)
{
    double RefDen = _t->TLRho0;
    double MaxDen = _t->TLRho0*10.0;
    int nDen = _t->nDen;
    double StepDen = MaxDen/(nDen-1);
    double * xDen = (double *) malloc(sizeof(double)*nDen);
    double * yEng = (double *) malloc(sizeof(double)*nDen);
    for(int k=0;k<nDen;k++)
    {
        xDen[k] = k*StepDen;
        yEng[k] = 0.0;
    }

    int IndexM = (int)(floor(RefDen/MaxDen*nDen));

    xDen[IndexM] = _t->TLRho0;
    yEng[IndexM] = 0.0;

    double K1,K2,K3,K4;
    for(int k=IndexM+1;k<nDen;k++)
    {
        StepDen = xDen[k] - xDen[k-1];
//        K1 = TilldEdr(_t,yEng[k-1]                  ,xDen[k-1]);
//        K2 = TilldEdr(_t,yEng[k-1] + 0.50*K1*StepDen,xDen[k-1]+0.50*StepDen);
//        K3 = TilldEdr(_t,yEng[k-1] + 0.75*K2*StepDen,xDen[k-1]+0.75*StepDen);
//        yEng[k] = yEng[k-1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*StepDen/9.0;

        K1 = TilldEdr(_t,yEng[k-1],xDen[k-1]);
        K2 = TilldEdr(_t,yEng[k-1] + 0.50*K1*StepDen,xDen[k-1]+0.50*StepDen);
        K3 = TilldEdr(_t,yEng[k-1] + 0.50*K2*StepDen,xDen[k-1]+0.50*StepDen);
        K4 = TilldEdr(_t,yEng[k-1] +      K3*StepDen,xDen[k-1]+     StepDen);
        yEng[k] = yEng[k-1] + (K1 + 2.0*K2 + 2.0*K3 + K4)*StepDen/6.0;
    }

    for(int k=IndexM-1;k>=1;k--)
    {
        StepDen = xDen[k] - xDen[k+1];
//        K1 = TilldEdr(_t,yEng[k+1],xDen[k+1]);
//        K2 = TilldEdr(_t,yEng[k+1] + 0.50*K1*StepDen,xDen[k+1]+0.50*StepDen);
//        K3 = TilldEdr(_t,yEng[k+1] + 0.75*K2*StepDen,xDen[k+1]+0.75*StepDen);
//        yEng[k] = yEng[k+1] + (2.0*K1 + 3.0*K2 + 4.0*K3)*StepDen/9.0;
        K1 = TilldEdr(_t,yEng[k+1],xDen[k+1]);
        K2 = TilldEdr(_t,yEng[k+1] + 0.50*K1*StepDen,xDen[k+1]+0.50*StepDen);
        K3 = TilldEdr(_t,yEng[k+1] + 0.50*K2*StepDen,xDen[k+1]+0.50*StepDen);
        K4 = TilldEdr(_t,yEng[k+1] +      K3*StepDen,xDen[k+1]+     StepDen);
        yEng[k] = yEng[k+1] + (K1 + 2.0*K2 + 2.0*K3 + K4)*StepDen/6.0;
    }

    yEng[0] = yEng[1];

    _t->xDen = xDen;
    _t->yEng = yEng;
    _t->StepDen = MaxDen/(nDen-1);
}

double TillTemp(const struct TillotsonTable * _t,double EngIn, double RhoIn)
{
    int DL = (int)(Wind(floor(RhoIn/_t->StepDen),0,_t->nDen-2));
    double Dxl = DL + 1.0 - RhoIn/_t->StepDen;
    double Dxr = 1.0 - Dxl;
    double rEngIn = Dxl*_t->yEng[DL] + Dxr*_t->yEng[DL+1];
    return Max((EngIn - rEngIn)/_t->TLCv + _t->TLTref,0.0);
}

double TillEOSInterpolateTP(struct TillotsonTable *_t, double tTem, double tPre, int DataId)
{
    double ThermalE = (tTem - _t->TLTref)*_t->TLCv;
    int DI = 0;
    int DJ = _t->nDen - 1;
    while(DI + 1 < DJ)
    {
        int DM = (DI + DJ)/2;
        double DEk = _t->yEng[DM] + ThermalE;
        double DPk,Csk;
        TillPres(_t,DEk,_t->xDen[DM],&DPk,&Csk);
        if(DPk > tPre)
            DJ = DM;
        else
            DI = DM;
    }

    double LIPre,LCs,RJPre,RCs;
    TillPres(_t,_t->yEng[DI]+ThermalE,_t->xDen[DI],&LIPre,&LCs);
    TillPres(_t,_t->yEng[DJ]+ThermalE,_t->xDen[DJ],&RJPre,&RCs);

    double LIDen = _t->xDen[DI], RJDen = _t->xDen[DJ], LMDen, LMEng, LMPre, LMCs;
    for(int k=0;k<50;k++)
    {
        double DenError = (RJDen - LIDen)/RJDen;
        LMDen = (LIDen + RJDen)/2.0;
        if(DenError < 1e-12)
            break;
        double Dxl = (_t->xDen[DJ] - LMDen)/_t->StepDen;
        double Dxr = 1.0 - Dxl;
        LMEng = Dxl*_t->yEng[DI] + Dxr*_t->yEng[DJ] + ThermalE;
        TillPres(_t,LMEng,LMDen,&LMPre,&LMCs);
        if(LMPre > tPre)
            RJDen = LMDen;
        else
            LIDen = LMDen;
    }

    if(-1==DataId)
        return LMDen;
    else if(EOSENG == DataId)
        return LMEng*LMDen;
    else if(EOSPRE == DataId)
        return tPre;
    else if(EOSCSD == DataId)
        return LMCs;
    else
    {
        fprintf(stdout,"unsupported data id in InterpolateTP\n");
        exit(0);
    }
}

void TillInitStateRef(struct TillotsonTable *_t, struct StateReference *_s)
{
    _s->MeltDen   = 0.85*_t->TLRho0;
    double GasMin = 1.0e-4*_t->TLRho0;
    _s->VaporDen = GasMin*10.0;
    _s->Viscosity = 0.0;

    _s->c_heat = _t->TLCv;
    _s->thermexp = _t->TLCv*(_t->TLb + _t->TLa)/(_t->TLA/_t->TLRho0);
    _s->den0 = _t->TLRho0;

    _s->csd0 = sqrt(_t->TLA/_t->TLRho0);
    _s->gamfrac = (_t->TLb + _t->TLa)/(_t->TLA/_t->TLRho0);
}


double TillEOSInterpolateTD(struct TillotsonTable *_t, double tTem, double tDen, int DataId)
{
    double ThermalE = (tTem - _t->TLTref)*_t->TLCv;
    int DI = 0;
    int DJ = _t->nDen - 1;
    while(DI + 1 < DJ)
    {
        int DM = (DI + DJ)/2;
        if(_t->xDen[DM] > tDen)
            DJ = DM;
        else
            DI = DM;
    }

    double Dxl = (_t->xDen[DJ]-tDen)/(_t->xDen[DJ]-_t->xDen[DI]);
    Dxl = Wind(Dxl,0.0,1.0);
    double Dxr = 1.0 - Dxl;
    double LMEng = _t->yEng[DI]*Dxl + _t->yEng[DJ]*Dxr + ThermalE;


    double LMPre,LMCs;
    TillPres(_t,LMEng,tDen,&LMPre,&LMCs);

    if(-1==DataId)
        return tDen;
    else if(0 == DataId)
        return LMEng*tDen;
    else if(1 == DataId)
        return LMPre;
    else if(2 == DataId)
        return LMCs;
    else
    {
        fprintf(stdout,"unsupported data id in InterpolateTP\n");
        exit(0);
    }
}

double LoadTillEOS(struct TillotsonTable * _t, FILE * fp)
{
    InputFile * ifp = ParseInputFile(fp);
    _t->TLRho0 = GetValueD(ifp,"Tillostson.Rho0","2.7e3");
    _t->TLCv   = GetValueD(ifp,"ColdEnergy.Cv","8.96e2");
    _t->TLA    = GetValueD(ifp,"Tillostson.BulkA","7.52e10");
    _t->TLB    = GetValueD(ifp,"Tillostson.BulkB","6.50e10");
    _t->TLE0   = GetValueD(ifp,"Tillostson.E0","5.0e6");
    _t->TLa    = GetValueD(ifp,"Tillostson.Ta","0.50");
    _t->TLb    = GetValueD(ifp,"Tillostson.Tb","1.63");
    _t->TLAlpha = GetValueD(ifp,"Tillostson.Alpha","5.0");
    _t->TLBeta = GetValueD(ifp,"Tillostson.Beta","5.0");
    _t->TLEiv = GetValueD(ifp,"Tillostson.Eiv","3.00e6");
    _t->TLEcv = GetValueD(ifp,"Tillostson.Ecv","1.39e7");
    _t->TLTref = GetValueD(ifp,"ColdEnergy.Tref","293.0");
    _t->nDen = GetValueI(ifp,"ColdEnergy.RhoStep","1501");
    CloseInputFile(ifp);
    TillColdEnergy(_t);

#ifdef F90_DEBUG
    FILE * efp = fopen("salec_ec.txt","w");
    fprintf(efp,"%d\n",_t->nDen);
    for(int j=0;j<_t->nDen;++j)
    {
        fprintf(efp,"%27.9f,%27.9f\n",_t->xDen[j],_t->yEng[j]);
    }
#endif
}

void UnAllocateTillEOS(struct TillotsonTable *_t)
{
    free(_t->xDen);
    free(_t->yEng);
}

void ANEOSLowDenCorrect(struct ANEOSTable * _t, double plimit)
{
    // make the pressure at low density& low temperature is negative
    for(int k=0;k<_t->nTem;++k)
    {
        for(int j=0;j<_t->nDen;++j)
        {
            _t->Data[k][j][EOSPRE] -= plimit;
        }
    }
}

void ANEOSWrite(struct ANEOSTable * _t,const char fname[], const char comment[] )
{
    // remove "\n" in comment
    // ...
    //
    FILE * fp = fopen(fname,"w");
    fprintf(fp,"# ANEOS table for granite\n"
               "# generated by EoSTool\n"
               "# %s\n",comment);
    // write the rows and columns of aneos table
    fprintf(fp,"%15d %15d\n",_t->nTem,_t->nDen);
    // set the normal state (T0,RHO0,P0)
    fprintf(fp,"%22.15e %22.15e %22.15e\n",_t->Tnorm,_t->Dnorm,_t->Pnorm);
    // write the rows (Temperature)
    int WriteCounter = 0;
    for(int k=0;k<_t->nTem;++k)
    {
        fprintf(fp,"%22.15e ",_t->xTem[k]);
        if(++WriteCounter%4 == 0) fprintf(fp,"\n");
    }

    // write the data (Internal energy, pressure, speed of sound)
    for(int k=0;k<_t->nTem;++k)
    {
        for(int j=0;j<_t->nDen;++j)
        {
            fprintf(fp,"%22.15e ",_t->Data[k][j][0]/_t->yDen[j]);
            if(++WriteCounter%4 == 0) fprintf(fp,"\n");
            fprintf(fp,"%22.15e ", _t->Data[k][j][1]);
            if(++WriteCounter%4 == 0) fprintf(fp,"\n");
            fprintf(fp,"%22.15e ",_t->Data[k][j][2]);
            if(++WriteCounter%4 == 0) fprintf(fp,"\n");
        }
    }

    // write the columns (Density)
    for(int k=0;k<_t->nDen;++k)
    {
        fprintf(fp,"%22.15e",_t->yDen[k]);
        if(++WriteCounter%4 == 0) fprintf(fp,"\n");
    }

    fclose(fp);
}

double calculate_state(eos_table * _e,double den,double eng,double * tem,double * pre,double * csd)
{
    switch (_e->type)
    {
        case ANEOS:
            aneos_calculate_state((aneos_table *) _e->table,den,eng,tem,pre,csd);
            break;
        case ANEOSSieDen:
            aneossieden_calculate_state((aneos_table *) _e->table,den,eng,tem,pre,csd);
            break;
        case Tillotson:
            tillotson_calculate_state((tillotson_table *)_e->table,den,eng,tem,pre,csd);
            break;
        default:
            fprintf(stderr,"unimplemented eos type in %s\n",__func__ );
    }

    return -1.0;
}

double calculate_state_correct(eos_table * _e,double den,double eng,double plim,double * tem,double * pre,double * csd, double * pty)
{
    /*
     * same as the calculate_state in most the time
     * if the pre of output < plim, this function will return the compression factor of this material
     * or cap the output pressure
     */
    double tmp_tem = 0., tmp_pre = 0., tmp_csd = 0., max_cfactor = 2., min_cfactor = 1.0, tmp_cfactor = 1.0;
    calculate_state(_e,den*pty[0],eng*pty[0],&tmp_tem,&tmp_pre,&tmp_csd);
    tmp_pre = tmp_pre/pty[0];
    if(pty[0] > 1.0)
    {
        double chi = _e->ref->chi;
        double alpha0 = _e->ref->alpha0;
        tmp_csd *= (1.0 + (chi - 1.0) * fmax(0.0, fmin(1.0, (*pty - 1.0) / (alpha0 - 1.0))));
    }

    if(tmp_pre < plim)
    {
#ifndef CAP_MINPRESSURE
        int k_cycles = 50; // consistent with TOLVOF
        while(fabs(plim - tmp_pre) > 0.01*fabs(plim)  && k_cycles > 0)
        {
            k_cycles--;
            tmp_cfactor = 0.5*(max_cfactor + min_cfactor);
            calculate_state(_e,den*tmp_cfactor,eng*tmp_cfactor,&tmp_tem,&tmp_pre,&tmp_csd);
            if(tmp_pre < plim)
            {
                min_cfactor = tmp_cfactor;
            }
            else
            {
                max_cfactor = tmp_cfactor;
            }
        }
#else
        tmp_pre = Max(tmp_pre,plim);
#endif
    }

    *tem = tmp_tem;
    *pre = tmp_pre;
    *csd = tmp_csd;
    return tmp_cfactor;
}

void aneos_calculate_state(aneos_table * _t,double den,double eng,double * tem,double * pre,double * csd)
{
    int DI = 0;
    int DJ = _t->nDen - 1;
    while (DI + 1 < DJ)
    {
        int DM = (DI + DJ) / 2;
        if (_t->yDen[DM] > den)
            DJ = DM;
        else
            DI = DM;
    }

    // the local coordinate for Density
    double Dxr = (den - _t->yDen[DI]) / (_t->yDen[DJ] - _t->yDen[DI]);
    Dxr = Wind(Dxr,0.0,1.0);
    double Dxl = 1.0 - Dxr;

    double Energy = eng;
    int TI = 0;
    double IRE = _t->Data[TI][DI][0] * Dxl + _t->Data[TI][DJ][0] * Dxr;
    int TJ = _t->nTem - 1;
    double JRE = _t->Data[TJ][DI][0] * Dxl + _t->Data[TJ][DJ][0] * Dxr;
    while (TI + 1 < TJ)
    {
        int TM = (TI + TJ) / 2;
        double MRE = _t->Data[TM][DI][0] * Dxl + _t->Data[TM][DJ][0] * Dxr;
        if (MRE > Energy)
        {
            TJ = TM;
            JRE = MRE;
        }
        else
        {
            TI = TM;
            IRE = MRE;
        }
    }

    // the local coordinate for temperature
    double Txr = (Energy - IRE) / (JRE - IRE);
    Txr = Wind(Txr,0.0,1.0);
    double Txl = 1.0 - Txr;
    *tem = _t->xTem[TI] * Txl + _t->xTem[TJ] * Txr;
    *pre = (_t->Data[TI][DI][1] * Dxl + _t->Data[TI][DJ][1] * Dxr) * Txl
           + (_t->Data[TJ][DI][1] * Dxl + _t->Data[TJ][DJ][1] * Dxr) * Txr;
    *csd = (_t->Data[TI][DI][2] * Dxl + _t->Data[TI][DJ][2] * Dxr) * Txl
           + (_t->Data[TJ][DI][2] * Dxl + _t->Data[TJ][DJ][2] * Dxr) * Txr;
}


void tillotson_calculate_state(tillotson_table * _t,double den,double eng,double * tem,double * pre,double * csd)
{
    double EngIn = eng/den;
    double RhoIn = den;
    *tem = TillTemp(_t,EngIn,RhoIn);
    TillPres(_t,EngIn,RhoIn,pre,csd);
}

void load_eos_table(eos_table * _e, FILE * fp)
{
    /*
     * this subroutine will allocate memory automatically
     */

    if(_e->type == ANEOS)
    {
        if(NULL == _e->table)
        {
            _e->table = malloc(sizeof(aneos_table));
        }
        LoadANEOS((aneos_table *)_e->table,fp);
        if(NULL == _e->ref)
        {
             _e->ref = (state_reference *) malloc(sizeof(state_reference));
        }
        ANEOSInitStateRef(_e->table,_e->ref);
    } else if(_e->type == Tillotson)
    {
        if(NULL == _e->table)
        {
            _e->table = malloc(sizeof(tillotson_table));
        }
        LoadTillEOS((tillotson_table *)_e->table,fp);
        if(NULL == _e->ref)
        {
            _e->ref = (state_reference *) malloc(sizeof(state_reference));
        }
        TillInitStateRef(_e->table,_e->ref);
    } else if(_e->type == ANEOSSieDen)
    {
        if(NULL == _e->table)
        {
            _e->table = malloc(sizeof(aneos_table));
        }
        LoadANEOSSieDen((aneos_table *)_e->table,fp);
        if(NULL == _e->ref)
        {
            _e->ref = (state_reference *) malloc(sizeof(state_reference));
        }
        ANEOSSieDenInitStateRef(_e->table,_e->ref);
    }
    else
    {
        fprintf(stdout,"cannot load unknown type eos file in %s\n",__func__);
        exit(0);
    }
}

double interpolate_eos_tp(eos_table * _e,double tem,double pre, int DataId)
{
    double ret_value = 0.;
    switch (_e->type)
    {
        case ANEOS:
            ret_value =  ANEOSInterpolateTP(_e->table,tem,pre,DataId);
            break;
        case Tillotson:
            ret_value = TillEOSInterpolateTP(_e->table,tem,pre,DataId);
            break;
        case ANEOSSieDen:
            ret_value = ANEOSSieDenInterpolateTP(_e->table,tem,pre,DataId);
            break;
        default:
            fprintf(stderr,"unimplemented eos type in %s\n",__func__ );
            exit(0);
    }
    return ret_value;
}

double interpolate_eos_td(eos_table * _e,double tem,double den, int DataId)
{
    double ret_value = 0.;

    switch (_e->type)
    {
        case ANEOS:
            ret_value = ANEOSInterpolateTD(_e->table,tem,den,DataId);
            break;
        case Tillotson:
            ret_value = TillEOSInterpolateTD(_e->table,tem,den,DataId);;
            break;
        case ANEOSSieDen:
            ret_value = ANEOSSieDenInterpolateTD(_e->table,tem,den,DataId);
            break;
        default:
            fprintf(stderr,"unimplemented eos type in %s\n",__func__ );
            exit(0);
    }
    return ret_value;
}


