#ifndef PP_H
#define PP_H

#include <stdlib.h>
#include <armadillo>
#include <time.h>

using namespace arma;

    extern int N; /*number of particles*/
    extern double L; /*lenght of the box*/
    extern double PASSI; /*number of steps*/
    extern double DELTA; /*evolution parameter*/
    extern double T; /*temperature*/
    extern double BINS; /*number of bins,better if submultiple of PASSI*/

struct param                        //a type consisting in three variables, which will be useful
{
	double first;
	double second;
};

class PP
{
    public:
    param jackknife (vec v); //calculates mean values and errors
};

#endif // PP_H
