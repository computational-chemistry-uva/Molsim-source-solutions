#include <math.h>
#include "system.h"

VECTOR Positions[MAX_PARTICLES];
VECTOR Forces[MAX_PARTICLES];
double Charges[MAX_PARTICLES];
double Alpha;
double Box;

int NumberOfParticles;
int NumberOfCells;
int NumberOfWavevectors;

// COMPLEMENTARY ERROR FUNCTION; THIS IS A REALLY GOOD ONE !!!
double ErrorFunctionComplement(double x)
{
  double t,u,y;
  const double PA=3.97886080735226000,P0=2.75374741597376782e-1,P1=4.90165080585318424e-1;
  const double P2=7.74368199119538609e-1,P3=1.07925515155856677,P4=1.31314653831023098;
  const double P5=1.37040217682338167,P6=1.18902982909273333,P7=8.05276408752910567e-1;
  const double P8=3.57524274449531043e-1,P9=1.66207924969367356e-2,P10=-1.19463959964325415e-1;
  const double P11=-8.38864557023001992e-2,P12=2.49367200053503304e-3,P13=3.90976845588484035e-2;
  const double P14=1.61315329733252248e-2,P15=-1.33823644533460069e-2,P16=-1.27223813782122755e-2;
  const double P17=3.83335126264887303e-3,P18=7.73672528313526668e-3,P19=-8.70779635317295828e-4;
  const double P20=-3.96385097360513500e-3,P21=1.19314022838340944e-4,P22=1.27109764952614092e-3;

  t=PA/(PA+fabs(x));
  u=t-0.5;
  y=(((((((((P22*u+P21)*u+P20)*u+P19)*u+P18)*u+P17)*u+P16)*u+P15)*u+P14)*u+P13)*u+P12;
  y=((((((((((((y*u+P11)*u+P10)*u+P9)*u+P8)*u+P7)*u+P6)*u+P5)*u+P4)*u+P3)*u+P2)*u+P1)*u+P0)*t*exp(-x*x);
  if(x<0) y=2.0-y;
  return y;
}

double ErrorFunction(double x)
{
  return 1.0-ErrorFunctionComplement(x);
}

