#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>
#include <armadillo>
#include <time.h>
#include "pp.h"


//Parameters of the system

using namespace arma;


    extern int N; /*number of particles*/
    extern double L; /*lenght of the box*/
    extern double PASSI; /*number of steps*/
    extern double DELTA; /*evolution parameter*/
    extern double T; /*temperature*/
    extern double BINS; /*number of bins,better if submultiple of PASSI*/


class Tools : public PP
{
    public:
    mat current;
	mat evolved;

    void pbcize(int q);
    double distance(mat v, int i, int j); /*applies pbc and calculates squared distance between two particles*/
    double energycounter(); /*calculates lennard-jones total energy*/
    double virialcounter(); /*calculates total virial*/
    param evolvevector (); /*evolves configuration, energy and virial*/
};
#endif // TOOLS_H
