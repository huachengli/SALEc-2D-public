//
// Created by huach on 9/11/2023.
//

/*
 * read the asteroid file and convert it to salec2d.inp
 *
 * generate some terms need copy to salec2d.inp
 */

#include "iSALEInputParser.h"

int main()
{
    fprintf(stdout,"####Read asteroid.inp&material.inp from iSALE2D:\n");
    fprintf(stdout,"####Following LINES can be put int to *.inp for salec2d\n");
    // load asteroid.inp
    InputFile * isale_ifp = OpeniSALEInputFile("./asteroid.inp");
    InputFile * salec_ifp = NewInputFile();
    GenerateMesh(isale_ifp,salec_ifp);
    CloseInputFile(isale_ifp);

    // load material.inp
    InputFile * isale_ifp_material = OpeniSALEInputFile("./material.inp");
    GenerateMaterial(isale_ifp_material,salec_ifp);
//    ViewInputFile(isale_ifp_material);
    ViewInputFile(salec_ifp);
    CloseInputFile(isale_ifp_material);
    CloseInputFile(salec_ifp);
    return 0;
}