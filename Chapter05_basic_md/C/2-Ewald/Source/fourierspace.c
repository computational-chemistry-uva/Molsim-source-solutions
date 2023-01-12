#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

COMPLEX lm[MAX_PARTICLES];
COMPLEX Ck[MAX_PARTICLES];

COMPLEX Em[MAX_PARTICLES][MAX_WAVEVECTOR+1];
COMPLEX En[MAX_PARTICLES][MAX_WAVEVECTOR+1];
COMPLEX El[MAX_PARTICLES][2];

// Fourier Part Of The Ewald Summation
void Fourierspace(double *Ukappa,double *Uself)
{
  int i,j,m,n;
  int Mmin,Nmin,Ll,Mm,Nn;
  double Cs,Rksq,temp,fac,energy_sum,exp_term;
  VECTOR Ss,Rk;
  double Uself_sum,UEwaldFourier;
  COMPLEX Cksum;
  MATRIX StressEwaldFourier;

  fac=0.0;
  StressEwaldFourier.ax=StressEwaldFourier.ay=StressEwaldFourier.az=0.0;
  StressEwaldFourier.bx=StressEwaldFourier.by=StressEwaldFourier.bz=0.0;
  StressEwaldFourier.cx=StressEwaldFourier.cy=StressEwaldFourier.cz=0.0;
  Ss.x=Ss.y=Ss.z=0.0;

  energy_sum=0.0;
  for(i=0;i<NumberOfParticles;i++)
  {
    Em[i][0].re=1.0;
    Em[i][0].im=0.0;
    El[i][0].re=1.0;
    El[i][0].im=0.0;
    En[i][0].re=1.0;
    En[i][0].im=0.0;
    Ss.x=2.0*M_PI*Positions[i].x/Box;
    Ss.y=2.0*M_PI*Positions[i].y/Box;
    Ss.z=2.0*M_PI*Positions[i].z/Box;
    El[i][1].re=cos(Ss.x);
    El[i][1].im=sin(Ss.x);
    Em[i][1].re=cos(Ss.y);
    Em[i][1].im=sin(Ss.y);
    En[i][1].re=cos(Ss.z);
    En[i][1].im=sin(Ss.z);
  }

  for(j=2;j<=NumberOfWavevectors;j++)
    for(i=0;i<NumberOfParticles;i++)
    {
      Em[i][j].re=Em[i][j-1].re*Em[i][1].re-Em[i][j-1].im*Em[i][1].im;
      Em[i][j].im=Em[i][j-1].im*Em[i][1].re+Em[i][j-1].re*Em[i][1].im;
    }
  for(j=2;j<=NumberOfWavevectors;j++)
    for(i=0;i<NumberOfParticles;i++)
    {
      En[i][j].re=En[i][j-1].re*En[i][1].re-En[i][j-1].im*En[i][1].im;
      En[i][j].im=En[i][j-1].im*En[i][1].re+En[i][j-1].re*En[i][1].im;
    }

  Mmin=0;
  Nmin=1;
  for(Ll=0;Ll<=NumberOfWavevectors;Ll++)     // Loop Over All K-Vectors
  {
    Rk.x=2.0*M_PI*Ll/Box;
    if(Ll==1)
      for(i=0;i<NumberOfParticles;i++)
      {
        El[i][0].re=El[i][1].re;
        El[i][0].im=El[i][1].im;
      }
    else if(Ll>1)
      for(i=0;i<NumberOfParticles;i++)
      {
        Cs=El[i][0].re;
        El[i][0].re=El[i][1].re*Cs-El[i][1].im*El[i][0].im;
        El[i][0].im=El[i][1].re*El[i][0].im+El[i][1].im*Cs;
      }

    for(Mm=Mmin;Mm<=NumberOfWavevectors;Mm++)
    {
      Rk.y=2.0*M_PI*Mm/Box;
      m=abs(Mm);
      if(Mm>=0)
        for(i=0;i<NumberOfParticles;i++)
        {
          lm[i].re=El[i][0].re*Em[i][m].re-El[i][0].im*Em[i][m].im;
          lm[i].im=El[i][0].im*Em[i][m].re+Em[i][m].im*El[i][0].re;
        }
      else
        for(i=0;i<NumberOfParticles;i++)
        {
          lm[i].re=El[i][0].re*Em[i][m].re+El[i][0].im*Em[i][m].im;
          lm[i].im=El[i][0].im*Em[i][m].re-Em[i][m].im*El[i][0].re;
        }

      for(Nn=Nmin;Nn<=NumberOfWavevectors;Nn++)
      {
        n=abs(Nn);
        Rk.z=2.0*M_PI*Nn/Box;
        Rksq=SQR(Rk.x)+SQR(Rk.y)+SQR(Rk.z);

        Cksum.re=0.0;
        Cksum.im=0.0;

        if(Nn>=0)
        {
          for(i=0;i<NumberOfParticles;i++)
          {
            Ck[i].re=Charges[i]*(lm[i].re*En[i][n].re-lm[i].im*En[i][n].im);
            Ck[i].im=Charges[i]*(lm[i].im*En[i][n].re+lm[i].re*En[i][n].im);
            Cksum.re+=Ck[i].re;
            Cksum.im+=Ck[i].im;
          }
        }
        else
        {
          for(i=0;i<NumberOfParticles;i++)
          {
            Ck[i].re=Charges[i]*(lm[i].re*En[i][n].re+lm[i].im*En[i][n].im);
            Ck[i].im=Charges[i]*(lm[i].im*En[i][n].re-lm[i].re*En[i][n].im);
            Cksum.re+=Ck[i].re;
            Cksum.im+=Ck[i].im;
          }
        }
        exp_term=exp((-0.25/SQR(Alpha))*Rksq)/Rksq;
        temp=(SQR(Cksum.re)+SQR(Cksum.im))*exp_term;

        energy_sum+=temp;

        // Stress tensor
        fac=2.0*(1.0/Rksq-(-0.25/SQR(Alpha)))*temp;
        StressEwaldFourier.ax-=fac*Rk.x*Rk.x;
        StressEwaldFourier.ay-=fac*Rk.x*Rk.x;
        StressEwaldFourier.az-=fac*Rk.x*Rk.x;

        StressEwaldFourier.bx-=fac*Rk.y*Rk.y;
        StressEwaldFourier.by-=fac*Rk.y*Rk.y;
        StressEwaldFourier.bz-=fac*Rk.y*Rk.y;

        StressEwaldFourier.cx-=fac*Rk.z*Rk.z;
        StressEwaldFourier.cy-=fac*Rk.z*Rk.z;
        StressEwaldFourier.cz-=fac*Rk.z*Rk.z;

        for(i=0;i<NumberOfParticles;i++)
        {
          // Forces
          fac=(8.0*M_PI/CUBE(Box))*(Ck[i].im*Cksum.re-Ck[i].re*Cksum.im)*exp_term;

          Forces[i].x+=fac*Rk.z;
          Forces[i].y+=fac*Rk.z;
          Forces[i].z+=fac*Rk.z;
        }
      }
      Nmin=-NumberOfWavevectors;
    }
    Mmin=-NumberOfWavevectors;
  }

  // energy
  UEwaldFourier=(4.0*M_PI/CUBE(Box))*energy_sum;

  Uself_sum=0.0;
  for(i=0;i<NumberOfParticles;i++)
    Uself_sum-=SQR(Charges[i]);
  Uself_sum*=Alpha/sqrt(M_PI);

  // stress tensor
  StressEwaldFourier.ax+=energy_sum;
  StressEwaldFourier.by+=energy_sum;
  StressEwaldFourier.cz+=energy_sum;

  fac=4.0*M_PI/CUBE(Box);
  StressEwaldFourier.ax*=fac;  StressEwaldFourier.bx*=fac;  StressEwaldFourier.cx*=fac;
  StressEwaldFourier.ay*=fac;  StressEwaldFourier.by*=fac;  StressEwaldFourier.cy*=fac;
  StressEwaldFourier.az*=fac;  StressEwaldFourier.bz*=fac;  StressEwaldFourier.cz*=fac;

  *Ukappa=UEwaldFourier;
  *Uself=Uself_sum;
}
