//
// Created by huach on 9/11/2023.
//

#include "iSALEInputParser.h"

int IsiSALEComment(const char c[])
{
    // check the first non-space character
    // if in "-#<", this line is comment
    if(strlen(c)==0)
        return 1;
    char CommentHead[5]   = "-#<";
    int IsComm = 0;
    for(int k=0;k<5;k++)
    {
        if(c[0] == CommentHead[k])
        {
            IsComm = 1;
            break;
        }
    }
    return IsComm;
}

int IsiSALEList(const char vname[])
{
    const char VnameList[][MaxStrLen] = {
          "GRIDH","GRIDV", "GRIDSPCM",
          "LAYMAT", "LAYPOS","LAYTPROF",
          "TR_SPCH","TR_SPCV"
          // other name of list adding here
    };

    int VnameListLen = sizeof(VnameList) / sizeof(VnameList[0]);
    for(int k=0;k<VnameListLen;++k)
    {
        if(strcasecmp(vname, VnameList[k]) == 0) return k;
    }
    return -1;
}


InputFile * OpeniSALEInputFile(const char fname[])
{
    FILE * fp = fopen(fname,"r");
    if(NULL == fp)
    {
        fprintf(stdout,"Cannot open input file * %s * \n",fname);
        exit(0);
    }
    InputFile * ifp = (InputFile *) malloc(sizeof(InputFile));
    ifp->Len = 0;

    char LineBuffer[300];
    char MainDelimiter[] = ":";
    char SubDelimiter[]  = "\"";
    int is_asteroid = 0;

    while(fscanf(fp,"%[^\n]",LineBuffer)!=EOF)
    {
        fgetc(fp);
        trim(LineBuffer);
        int LineBufferLen = strlen(LineBuffer);
        if(strcasecmp(LineBuffer,"#ISINP") == 0)
        {
            is_asteroid = 1;
            continue;
        }
        if(strcasecmp(LineBuffer,"#ISINP") == 0)
        {
            is_asteroid = 0;
            continue;
        }

        if((LineBufferLen==0) || IsiSALEComment(LineBuffer)) continue;


        char tKey[MaxStrLen], tValue[MaxStrLen], tKeyShort[MaxStrLen];
        int r = Strok(LineBuffer,MainDelimiter,tKey);
        Strok(LineBuffer+r,SubDelimiter,tValue);
        if((0== strlen(tKey)) || (0== strlen(tValue)))
        {
            LineBuffer[0] = '\0';
            continue;
        }
        Strok(tKey," ",tKeyShort); // remove the comment after VARNAME

        strcpy(ifp->Key[ifp->Len], tKeyShort);
        if(IsiSALEList(tKeyShort) >= 0 || is_asteroid == 0)
        {
            // all the variables in material.inp is list
            sprintf(ifp->Value[ifp->Len],"[%s]",tValue);
        }
        else
        {
            strcpy(ifp->Value[ifp->Len],tValue);
        }


        trim(ifp->Key[ifp->Len]);
        trim(ifp->Value[ifp->Len]);
        int ncomp = Strrpl(ifp->Value[ifp->Len],MainDelimiter,',');
        ifp->Len++;

        // get the number of materials
        if(strcasecmp("MATNAME",tKeyShort)==0)
        {
            sprintf(ifp->Key[ifp->Len],"NM");
            sprintf(ifp->Value[ifp->Len],"%d",ncomp+1);
            ifp->Len++;
        }
        LineBuffer[0] = '\0';
    }
    fclose(fp);
    SortInputFile(ifp,0,ifp->Len-1);
    return ifp;
}

void * CloseiSALEInputFile(InputFile * ifp)
{
    free(ifp);
}

int ViewInputFile2(InputFile * ifp)
{
    fprintf(stdout,"#### %d variables IN InputFile ####\n",ifp->Len);
    fprintf(stdout,"[ Id]:       Key = Value\n");
    for(int k=0;k<ifp->Len;++k)
    {
        fprintf(stdout,"[%4d]%-25s = %s\n",k,ifp->Key[k],ifp->Value[k]);
    }
    fprintf(stdout,"#### End of InputFile ####\n");
    return ifp->Len;
}

int ViewInputFile(InputFile * ifp)
{
    fprintf(stdout,"#### %d variables IN InputFile ####\n",ifp->Len);
    for(int k=0;k<ifp->Len;++k)
    {
        fprintf(stdout,"%-25s = %s\n",ifp->Key[k],ifp->Value[k]);
    }
    fprintf(stdout,"#### End of InputFile ####\n");
    return ifp->Len;
}


int GenerateMesh(InputFile * ifp, InputFile * new_ifp)
{
    int nx = GetValueIk(ifp,"GRIDH",1,"0");
    int ex[2] = {
            GetValueIk(ifp,"GRIDH",0,"0"),
            GetValueIk(ifp,"GRIDH",2,"0")
    };
    int ny = GetValueIk(ifp,"GRIDV",1,"0");
    int ey[2] = {
            GetValueIk(ifp,"GRIDV",0,"0"),
            GetValueIk(ifp,"GRIDV",2,"0")
    };
    int Nx = nx + ex[0] + ex[1];
    int Ny = ny + ey[0] + ey[1];
    double extent_factor = GetValueDf90(ifp,"GRIDEXT","1.");
    double grid_spec = GetValueDf90(ifp,"GRIDSPC","1.");
    double dx_max[2] = {
            GetValueDkf90(ifp,"GRIDSPCM",0,"0.0"),
            GetValueDkf90(ifp,"GRIDSPCM",1,"0.0")
    };

    if(dx_max[1] >= 0. && dx_max[1] <= 1.0e-2)
        dx_max[1] = dx_max[0];

    if(dx_max[0] > 0)
        dx_max[0] = dx_max[0]/grid_spec;
    else
        dx_max[0] = -1.0*dx_max[0];

    if(dx_max[1] > 0)
        dx_max[1] = dx_max[1]/grid_spec;
    else
        dx_max[1] = -1.0*dx_max[1];

    double Oy = 0.;
    double Ox = 0.;

    int obj_num = GetValueI(ifp,"OBJNUM","1");
    int lay_num = GetValueI(ifp,"LAYNUM","1");
    int lay_pos[MaxObjNum] = {0};
    for(int k=0;k<lay_num;++k)
    {
        lay_pos[k] = GetValueIk(ifp,"LAYPOS",k,"1");
    }

    // build the dy and dx distribution
    double * dx, *xpos;
    int xpml[2] = {0};
    dx = (double *) malloc(sizeof(double)*(Nx));
    xpos = (double *) malloc(sizeof(double)*(Nx+1));

    for(int k=0;k<Nx;++k)
    {
        dx[k] = 1.0;
    }
    for(int k= ex[0] - 1;k>=0;--k)
    {
        dx[k] = extent_factor*dx[k+1];
        if(dx[k] > dx_max[0])
        {
            dx[k] = dx_max[0];
            xpml[0]++;
        }
    }
    for(int k=Nx - ex[1];k<Nx;++k)
    {
        dx[k] = extent_factor*dx[k-1];
        if(dx[k] > dx_max[0])
        {
            dx[k] =  dx_max[0];
            xpml[1] ++;
        }
    }
    xpos[0] = 0.;
    for(int k=0;k<Nx;++k)
    {
        xpos[k+1] = xpos[k] + dx[k];
    }

    double * dy, *ypos;
    int ypml[2] = {0};
    dy = (double *) malloc(sizeof(double)*(Ny));
    ypos = (double *) malloc(sizeof(double)*(Ny+1));
    for(int k=0;k<Ny;++k)
    {
        dy[k] = 1.0;
    }
    for(int k=ey[0] - 1;k>=0;--k)
    {
        dy[k] = extent_factor*dy[k+1];
        if(dy[k] > dx_max[1])
        {
            dy[k] = dx_max[1];
            ypml[0]++;
        }
    }
    for(int k=Ny - ey[1];k<Ny;++k)
    {
        dy[k] = extent_factor*dy[k-1];
        if(dy[k] > dx_max[1])
        {
            dy[k] = dx_max[1];
            ypml[1]++;
        }
    }
    ypos[0] = 0.;
    for(int k=0;k<Ny;++k)
    {
        ypos[k+1] = ypos[k] + dy[k];
    }

    for(int k=0;k<lay_num;++k)
    {
        Oy = Oy>ypos[lay_pos[k]]/ypos[Ny]?Oy:ypos[lay_pos[k]]/ypos[Ny];
    }

    char tKey[MaxStrLen], tValue[MaxStrLen];
    snprintf(tKey,MaxStrLen,"mesh.dx");
    snprintf(tValue,MaxStrLen,"%f",grid_spec);
    AddKeyValue(new_ifp,tKey,tValue);

    snprintf(tKey,MaxStrLen,"mesh.dy");
    AddKeyValue(new_ifp,tKey,tValue);

    snprintf(tKey,MaxStrLen,"mesh.npx");
    snprintf(tValue,MaxStrLen,"%d",Nx);
    AddKeyValue(new_ifp,tKey,tValue);

    snprintf(tKey,MaxStrLen,"mesh.npy");
    snprintf(tValue,MaxStrLen,"%d",Ny);
    AddKeyValue(new_ifp,tKey,tValue);

    snprintf(tKey,MaxStrLen,"mesh.Ox");
    snprintf(tValue,MaxStrLen,"%f",Ox);
    AddKeyValue(new_ifp,tKey,tValue);

    snprintf(tKey,MaxStrLen,"mesh.Oy");
    snprintf(tValue,MaxStrLen,"%f",Oy);
    AddKeyValue(new_ifp,tKey,tValue);

    AddKeyValue(new_ifp,"mesh.ext","on");
    snprintf(tKey,MaxStrLen,"mesh.ex");
    snprintf(tValue,MaxStrLen,"[%f,%d,%d]",extent_factor,ex[0]-xpml[0],ex[1]-xpml[1]);
    AddKeyValue(new_ifp,tKey,tValue);

    snprintf(tKey,MaxStrLen,"mesh.ey");
    snprintf(tValue,MaxStrLen,"[%f,%d,%d]",extent_factor,ey[0]-ypml[0],ey[1]-ypml[1]);
    AddKeyValue(new_ifp,tKey,tValue);


    if(xpml[0]+xpml[1]+ypml[1]+ypml[1] > 0)
    {
        AddKeyValue(new_ifp,"numerical.damping","ctlg");

        snprintf(tKey,MaxStrLen,"numerical.dampX");
        snprintf(tValue,MaxStrLen,"%d",xpml[1]);
        AddKeyValue(new_ifp,tKey,tValue);

        snprintf(tKey,MaxStrLen,"numerical.dampY");
        snprintf(tValue,MaxStrLen,"%d",ypml[0]);
        AddKeyValue(new_ifp,tKey,tValue);

        snprintf(tKey,MaxStrLen,"numerical.dampX_LEFT");
        snprintf(tValue,MaxStrLen,"%d",xpml[0]);
        AddKeyValue(new_ifp,tKey,tValue);

        snprintf(tKey,MaxStrLen,"numerical.dampY_TOP");
        snprintf(tValue,MaxStrLen,"%d",ypml[1]);
        AddKeyValue(new_ifp,tKey,tValue);
    }
    else
    {
        AddKeyValue(new_ifp,"numerical.damping","off");
    }

    // add depth of layers
    snprintf(tValue,MaxStrLen,"%d",lay_num);
    AddKeyValue(new_ifp,"target.number",tValue);
    tValue[0] = '\0';
    strcat(tValue,"[");
    for(int k=0;k<lay_num;++k)
    {
        double layer_top = ypos[lay_pos[k]];
        double layer_bottom = 0.;
        if(k + 1< lay_num)
            layer_bottom = ypos[lay_pos[k]];
        char tmp[MaxStrLen];
        snprintf(tmp,MaxStrLen,"%.8e,",layer_top - layer_bottom);
        strcat(tValue,tmp);
    }
    int tValueLen = strlen(tValue);
    tValue[tValueLen-1] = ']';
    AddKeyValue(new_ifp,"target.depth",tValue);

}

int GetStrkf90(char strlist[],int pos, char value[])
{
    int len = strlen(strlist);
    int lpos = 0,rpos=0,k=0;
    for(int i=0;i<len;i++)
    {
        if(':'==strlist[i])
        {
            lpos++;
            continue;
        }
        if(lpos>pos) break;
        if(lpos==pos)
        {
            value[k++] = strlist[i];
        }
    }
    value[k]='\0';
    return k;
}


double GetValueDf90(InputFile * ifp,const char key[], char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
    {
        char tmp[MaxStrLen];
        strcpy(tmp,ifp->Value[r]);
        Strrpl(tmp,"DE",'e');
        return strtod(tmp,NULL);
    }
}

double GetValueDkf90(InputFile * ifp,const char key[], int k,char dvalue[])
{
    int r = SearchInput(ifp,key);
    if(r < 0)
        return strtod(dvalue,NULL);
    else
    {
        char tmp[MaxStrLen];
        if(GetStrk(ifp->Value[r],k,tmp))
        {
            Strrpl(tmp,"DE",'e');
            return strtod(tmp,NULL);
        } else
        {
            return strtod(dvalue,NULL);
        }
    }
}

double GenerateMaterial(InputFile * ifp, InputFile * new_ifp)
{
    char tValue[MaxStrLen],tmp[MaxStrLen];
    int tlen = 0;
    // check for material name and eos name
    int nm = GetValueI(ifp,"NM","1"); // add vacuum
    snprintf(tValue,MaxStrLen,"%d",nm+1);
    AddKeyValue(new_ifp,"material.nm",tValue);

    GetValueS(ifp,"MATNAME",tmp,"dunite_");
    snprintf(tValue,MaxStrLen,"[vacuum_,%s]",tmp+1);
    tlen = strlen(tValue); tValue[tlen-1] = '\0';
    AddKeyValue(new_ifp,"material.name",tValue);

    tValue[0] = '['; tValue[1] = '\0';
    strcat(tValue,"null   ");
    for(int k=0;k<nm;++k)
    {
        strcat(tValue,",");
        GetValueSk(ifp,"EOSTYPE",tmp,k,"unknown");
        if(strcasecmp("aneos",tmp) == 0)
        {
            strcat(tValue,"aneos");
        }
        else if(strcasecmp("tillo",tmp) == 0)
        {
            strcat(tValue,"tillotson");
        }
        else
        {
            strcat(tValue,"unknown");
        }
    }
    strcat(tValue,"]");
    AddKeyValue(new_ifp,"material.postfix",tValue);

    // poisson number
    ConvertInputVAR(ifp,new_ifp,"POIS","material.poisson",'f');

    // Simon melt model and ohnaka soft
    ConvertInputVAR(ifp,new_ifp,"TMELT0","material.SimonT0",'f');
    ConvertInputVAR(ifp,new_ifp,"TFRAC","material.OhnakaXi",'f');
    ConvertInputVAR(ifp,new_ifp,"ASIMON","material.SimonA",'f');
    ConvertInputVAR(ifp,new_ifp,"CSIMON","material.SimonC",'f');

    // Rock strength
    ConvertInputVAR(ifp,new_ifp,"YDAM0","material.ydam0",'f');
    ConvertInputVAR(ifp,new_ifp,"FRICDAM","material.ydamfri",'f');
    ConvertInputVAR(ifp,new_ifp,"YLIMDAM","material.ydamlim",'f');
    ConvertInputVAR(ifp,new_ifp,"YINT0","material.yint0",'f');
    ConvertInputVAR(ifp,new_ifp,"FRICINT","material.yintfri",'f');
    ConvertInputVAR(ifp,new_ifp,"YLIMINT","material.yintlim",'f');

    //Tillotson strength
    ConvertInputVAR(ifp,new_ifp,"JC_A","material.JcA",'f');
    ConvertInputVAR(ifp,new_ifp,"JC_B","material.JcB",'f');
    ConvertInputVAR(ifp,new_ifp,"JC_N","material.JcN",'f');
    ConvertInputVAR(ifp,new_ifp,"JC_C","material.JcC",'f');
    ConvertInputVAR(ifp,new_ifp,"JC_M","material.JcM",'f');
    ConvertInputVAR(ifp,new_ifp,"JC_TREF","material.JcTref",'f');

    // Ivanov damage
    ConvertInputVAR(ifp,new_ifp,"IVANOV_A","material.IvanA",'f');
    ConvertInputVAR(ifp,new_ifp,"IVANOV_B","material.IvanB",'f');
    ConvertInputVAR(ifp,new_ifp,"IVANOV_C","material.IvanC",'f');

    // Collins model

    // ACFL damage model
    ConvertInputVAR(ifp,new_ifp,"GAMETA","material.GammaEta",'f');
    ConvertInputVAR(ifp,new_ifp,"GAMBETA","material.GammaBeta",'f');
}

int ConvertInputVAR(InputFile * ifp,InputFile * new_ifp, const char _varname[], const char _newname[],char _type)
{
    char tmp[MaxStrLen],tValue[MaxStrLen];
    GetValueS(ifp,_varname,tmp,"[ ]");
    snprintf(tValue,MaxStrLen,"[,%s]",tmp+1);
    int tlen = strlen(tValue); tValue[tlen-1] = '\0';
    if(_type == 'f' || _type == 'D')
    {
        Strrpl(tValue,"DE",'e');
    }
    return AddKeyValue(new_ifp,_newname,tValue);
}