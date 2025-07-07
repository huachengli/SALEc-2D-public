//
// Created by huach on 9/11/2023.
//

#ifndef SALE_REBUILD_ASTEROIDINPUTPARSER_H
#define SALE_REBUILD_ASTEROIDINPUTPARSER_H

#include "InputParser.h"
#define MaxObjNum 20 // this value in iSALE is 10
InputFile * OpeniSALEInputFile(const char fname[]);
void * CloseiSALEInputFile(InputFile * ifp);
int ViewInputFile(InputFile * ifp);
int GenerateMesh(InputFile * ifp,InputFile * new_ifp);
double GenerateMaterial(InputFile * ifp, InputFile * new_ifp);
double GetValueDf90(InputFile * ifp,const char key[], char dvalue[]);
double GetValueDkf90(InputFile * ifp,const char key[], int k,char dvalue[]);
int ConvertInputVAR(InputFile * ifp,InputFile * new_ifp, const char _varname[], const char _newname[],char _type);

#endif //SALE_REBUILD_ASTEROIDINPUTPARSER_H
