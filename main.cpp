#include <iostream>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <stdlib.h>
#include <armadillo>
#include <time.h>
#include "tools.h"
#include "pp.h"



using namespace std;
using namespace arma;

int N; /*number of particles*/
double L; /*lenght of the box*/
double PASSI; /*number of steps*/
double DELTA; /*evolution parameter*/
double T; /*temperature*/
double BINS; /*number of bins,better if submultiple of PASSI*/


int main(int argc, char *argv[])
{
    srand (time(NULL));

    N=100;
    L=10;
    PASSI=10000;
    DELTA=0.5;
    T=1.0;
    BINS=PASSI/100;

	int a=0,i=0,j=0;

	param deltas={0,0,0};
	param energy={0,0,0};
	param pressure={0,0,0};

    Tools suzzu;
    PP postproc;

    suzzu.current = zeros<mat>(N,3);
    suzzu.evolved = zeros<mat>(N,3);
	vec energyvector = zeros<vec>(PASSI);
	vec pressurevector = zeros<vec>(PASSI);

	for(i=0;i<N;i++)
	{
        for(j=0;j<3;j++)
        {
            suzzu.current(i,j)=((double) rand() / (RAND_MAX))*L-L/2;
        }
	}

	energyvector(0)=suzzu.energycounter();
	pressurevector(0)=suzzu.pressurecounter();

    i=1;
    do
	{
		deltas=suzzu.evolvevector();

		if (deltas.error!=1)
		{
			energyvector(i)=energyvector(i-1)+deltas.first;
			pressurevector(i)=pressurevector(i-1)+deltas.second;
			i++;
			a++;
		}
		else
		{
			energyvector(i)=energyvector(i-1);
			pressurevector(i)=pressurevector(i-1);
			i++;
		}

	}while (i<500); //shakedown, discard first elements which are not significant

    energyvector(0)=energyvector(499);
    pressurevector(0)=pressurevector(499);
    i=1;
	do
	{
		deltas=suzzu.evolvevector();

		if (deltas.error!=1)
		{
			energyvector(i)=energyvector(i-1)+deltas.first;
			pressurevector(i)=pressurevector(i-1)+deltas.second;
			i++;
			a++;
		}
		else
		{
			energyvector(i)=energyvector(i-1);
			pressurevector(i)=pressurevector(i-1);
			i++;
		}

	}while (i<PASSI);

	energy=postproc.jackknife(energyvector);
	pressure=postproc.jackknife(pressurevector);
    //cout<<energyvector<<endl;
    cerr<<"#! rate of acceptance is "<<(double)a/i<<endl;
    cerr<<"#! energy is "<<energy.first<<" with error "<<fabs(energy.second/energy.first)*100<<"%"<<endl;
	cerr<<"#! pressure is "<<(pressure.first+N*T)/(L*L*L)<<" with error "<<fabs(pressure.second/pressure.first)*100<<"%"<<endl;

    return 0;
}
