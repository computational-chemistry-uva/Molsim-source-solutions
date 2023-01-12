#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

// setup simple cubic lattice
void Lattice(void)
{
  int Ix,Iy,Iz,i,Ip,In;
  double eps,tol,tol1,a,precision;
  int k;

  i=0;
  Ip=0;
  In=0;

  for(Ix=0;Ix<NumberOfCells;Ix++)
  {
    for(Iy=0;Iy<NumberOfCells;Iy++)
    {
      for(Iz=0;Iz<NumberOfCells;Iz++)
      {
        i=Iz+Iy*NumberOfCells+Ix*NumberOfCells*NumberOfCells;

        if(i>MAX_PARTICLES)
        {
          printf("Increase Maxpart\n");
          exit(0);
        }

        Positions[i].x=Ix+0.5;
        Positions[i].y=Iy+0.5;
	Positions[i].z=Iz+0.5;
	//MODIFICATION FOR PART 4: random displacements on the lattice
	   /*
	Positions[i].x=Ix+0.5+0.3*(RandomNumber()-0.5);
        Positions[i].y=Iy+0.5+0.3*(RandomNumber()-0.5);
        Positions[i].z=Iz+0.5+0.3*(RandomNumber()-0.5);
	*/
        if(((Ix+Iy+Iz)%2)==0)
        {
          Charges[i]=1.0;
          Ip++;
        }
        else
        {
          Charges[i]=-1.0;
          In++;
        }
      }
    }
  }

  NumberOfParticles=CUBE(NumberOfCells);
  Box=NumberOfCells;

  // write info
  printf("Box                  : %f\n",Box);
  printf("Number Of Cells      : %d\n",NumberOfCells);
  printf("Number Of Ions       : %d\n",NumberOfParticles);
  printf("Number Of Pos. Ions  : %d\n",Ip);
  printf("Number Of Neg. Ions  : %d\n",In);
  printf("Alpha                : %f\n",Alpha);
  printf("Kmax                 : %d\n",NumberOfWavevectors);
  printf("\n");
  for(precision=1e-3;precision>1e-10;precision/=10.0)
  {
    eps=MIN(precision>0?precision:-precision,0.5);
    tol=sqrt(fabs(log(eps*0.5*Box)));
    a=sqrt(fabs(log(eps*0.5*Box*tol)))/(0.5*Box);
    tol1=sqrt(-log(eps*0.5*Box*SQR(2.0*tol*a)));
    k=rint(0.25+Box*a*tol1/M_PI);
    printf("precision: %lg recommended Alpha: %lg recommended: %d\n",precision,a,k);
  }
  printf("\n");
}
