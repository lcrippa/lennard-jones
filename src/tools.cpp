#include "tools.h"
#include <math.h>
#include <time.h>
#include "pp.h"

void Tools::pbcize(int q)
{
    int j;
    double hL=L*0.5;

    for(j=0;j<3;j++)
    {
        if (evolved(q,j)>=hL)
        {
            evolved(q,j)=evolved(q,j)-L;
        }
        else if (evolved(q,j)<-hL)
        {
            evolved(q,j)=evolved(q,j)+L;
        }
    }
}

double Tools::distance(mat v, int i, int j)
{
	double dx,dy,dz,hL,dist=0;
	hL=L*0.5;

	dx  = v(i,0)-v(j,0);
	dy  = v(i,1)-v(j,1);
	dz  = v(i,2)-v(j,2);

	if (dx>hL)       dx=dx-L;
	else if (dx<-hL) dx=dx+L;
	if (dy>hL)       dx=dx-L;
	else if (dy<-hL) dx=dx+L;
	if (dz>hL)       dx=dx-L;
	else if (dz<-hL) dx=dx+L;

	dist = dx*dx + dy*dy + dz*dz;

	return dist;
}

double Tools::energycounter()
{
	int i,j;
	double U,dist;
	U=0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<i;j++)
		{
            dist=distance(current,i,j);
            U=U+4.0*(pow(dist,-6)-pow(dist,-3));
		}
	}
	return U;
}

double Tools::pressurecounter()
{
	int i,j;
	double p=0,dist;

	for(i=0;i<N;i++)
	{
		for(j=0;j<i;j++)
		{
			dist=distance(current,i,j);
			p=p+(4.0/3.0)*(12.0*pow(dist,-6)-6.0*pow(dist,-3));
		}
	}
	return p;
}

param Tools::evolvevector()
{
	int j,z;
	double random,dold,dnew;
	param deltas={0,0,0};

	z=rand() % N;
	evolved=current;

    random=(double) rand() / (RAND_MAX);
	evolved(z,0)=evolved(z,0)+2*DELTA*(random-0.5);
	random=(double) rand() / (RAND_MAX);
	evolved(z,1)=evolved(z,1)+2*DELTA*(random-0.5);
	random=(double) rand() / (RAND_MAX);
    evolved(z,2)=evolved(z,2)+2*DELTA*(random-0.5);

    pbcize(z);

	random=(double) rand() / (RAND_MAX);

    for(j=0;j<N;j++)
    {
        if (j!=z)
        {
            dold=distance(current,z,j);
            dnew=distance(evolved,z,j);
            deltas.first=deltas.first+4.0*(pow(dnew,-6)-pow(dnew,-3))-4.0*(pow(dold,-6)-pow(dold,-3));
        }
    }
	if(random<exp(-deltas.first/T))
	{
		for(j=0;j<N;j++)
		{
			if (j!=z)
			{
                dold=distance(current,z,j);
                dnew=distance(evolved,z,j);
				deltas.second+=(4.0/3.0)*(12.0/pow(dnew,6)-6.0/pow(dnew,3))-(4.0/3.0)*(12.0/pow(dold,6)-6.0/pow(dold,3));
			}
		}
		current=evolved;
		deltas.error=0;
	}
	else
	{
		deltas.error=1;
	}
    return deltas;
}
