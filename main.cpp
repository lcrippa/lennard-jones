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
    PASSI=10000;
    DELTA=0.5;
    T=1;
    BINS=PASSI/50;

    double rho,random,e,vir;
	int i=0,j=0,shakedown=2000;

	param deltas={0,0};
	param energy={0,0};
	param virial={0,0};

    Tools stuff;
    PP postproc;

    stuff.current = zeros<mat>(N,3);
    stuff.evolved = zeros<mat>(N,3);
	vec energyvector = zeros<vec>(PASSI);
	vec virialvector = zeros<vec>(PASSI);

    for(rho=0.001;rho<0.9;rho=rho+0.05)
    {

        L=pow(N/rho,0.333333);

        do //initialize positions, repeat until initial randomization is not absurd
        {
            for(i=0;i<N;i++)
            {
                for(j=0;j<3;j++)
                {
                    random=(double) rand() / RAND_MAX;
                    stuff.current(i,j)=(random-0.5)*L;
                }
            }

            e=stuff.energycounter();
            vir=stuff.virialcounter();
        }while(vir>10000000000);

        for(i=0;i<shakedown;i++)
        {
            deltas=stuff.evolvevector();
            e=e+deltas.first;
            vir=vir+deltas.second;
        }

        energyvector(0)=e;
        virialvector(0)=vir;

        for(i=1;i<PASSI;i++)
        {
            deltas=stuff.evolvevector();
            energyvector(i)=energyvector(i-1)+deltas.first;
            virialvector(i)=virialvector(i-1)+deltas.second;

        }

        energy=postproc.jackknife(energyvector);
        virial=postproc.jackknife(virialvector);

        cout<<rho<<"    "<<rho*T+virial.first/(L*L*L)<<"    "<<virial.second/(L*L*L)<<endl;
    }
    return 0;
}
