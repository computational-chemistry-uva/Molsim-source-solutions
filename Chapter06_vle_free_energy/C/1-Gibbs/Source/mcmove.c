#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"
#include "conf.h"  
#include "system.h"  

int Mcmove(double En[2], double Vir[2], int *Attempt, int *Nacc, double Dr)
{
  // Attempts To Displace A Randomly Selected Particle
  
  double Enn, Eno, Xn, Yn, Zn, Viro, Virn;
  
  // WHAT IS THIS?
  //   Dimension En(*), Vir(*)
  int Idmove, Jb, Ido;
  
  (*Attempt)++;
  Jb = 0;
  
  // ---Select A Particle At Random
  Idmove = (int)((double)(Npart)*RandomNumber());
  Ido = Id[Idmove];
  
  // ---Calculate Energy Old Configuration
  Eneri(X[Idmove], Y[Idmove], Z[Idmove], Idmove, Jb, &Eno, &Viro, Ido);
  
  // ---Give Particle A Random Displacement
  Xn = X[Idmove] + (RandomNumber()-0.5)*Dr;
  Yn = Y[Idmove] + (RandomNumber()-0.5)*Dr;
  Zn = Z[Idmove] + (RandomNumber()-0.5)*Dr;
  
  // ---Calculate Energy New Configuration:
  Eneri(Xn, Yn, Zn, Idmove, Jb, &Enn, &Virn, Ido);
  
  // ---Acceptance Test
  if(RandomNumber()<exp(-Beta*(Enn-Eno)))
    {
      //    --Accepted
      (*Nacc)++;
      En[Ido]+=(Enn-Eno);
      Vir[Ido]+=(Virn-Viro);
      
      //    ---Put Particle In Simulation Box
      if(Xn<0)        Xn = Xn + Box[Ido];
      if(Xn>Box[Ido]) Xn = Xn - Box[Ido];
      if(Yn<0)        Yn = Yn + Box[Ido];
      if(Yn>Box[Ido]) Yn = Yn - Box[Ido];
      if(Zn<0)        Zn = Zn + Box[Ido];
      if(Zn>Box[Ido]) Zn = Zn - Box[Ido];
      X[Idmove] = Xn;
      Y[Idmove] = Yn;
      Z[Idmove] = Zn;
    }
  return(0);
}
