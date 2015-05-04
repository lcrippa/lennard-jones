#include "tools.h"
#include <math.h>
#include <time.h>
#include "pp.h"

void Tools::pbcize(int q) //applies pbc (only on one particle because I evolve only one each time)
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

double Tools::distance(mat v, int i, int j) //returns squared distance
{
	double dx,dy,dz,hL,dist=0;
	hL=L*0.5;

	dx  = v(i,0)-v(j,0);
	dy  = v(i,1)-v(j,1);
	dz  = v(i,2)-v(j,2);

	if (dx>hL)       {dx=dx-L;}
	else if (dx<-hL) {dx=dx+L;}
	if (dy>hL)       {dy=dy-L;}
	else if (dy<-hL) {dy=dy+L;}
	if (dz>hL)       {dz=dz-L;}
	else if (dz<-hL) {dz=dz+L;}

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
            dist=pow(distance(current,i,j),-3);
            U=U+4.0*dist*(dist-1.0);
		}
	}
	return U;
}

double Tools::virialcounter()
{
	int i,j;
	double p=0,dist;

	for(i=0;i<N;i++)
	{
		for(j=0;j<i;j++)
		{
			dist=pow(distance(current,i,j),-3);
			p=p+16.0*dist*dist-8.0*dist;
		}
	}
	return p;
}

param Tools::evolvevector()
{
	int i,j,z;
	double random,dold,dnew;
	param deltas={0,0};

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
            dold=pow(distance(current,z,j),-3);
            dnew=pow(distance(evolved,z,j),-3);
            deltas.first=deltas.first+4.0*dnew*(dnew-1.0)-4.0*dold*(dold-1.0);
        }
    }
	if(random<exp(-deltas.first/T))
	{
		for(j=0;j<N;j++)
		{
			if (j!=z)
			{
                dold=pow(distance(current,z,j),-3);
                dnew=pow(distance(evolved,z,j),-3);
				deltas.second=deltas.second+16.0*(dnew*dnew-dold*dold)-8.0*(dnew-dold);
			}
		}
		current=evolved;
	}
	else
	{
		deltas.first=0.0;
		deltas.second=0.0;
	}
    return deltas;
}
