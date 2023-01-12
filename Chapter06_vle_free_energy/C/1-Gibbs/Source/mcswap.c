#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"
#include "conf.h"  
#include "chem.h"  

int Mcswap(double En[2], double Vir[2],int *Attempt, int *Acc)
{
  // ---Exchange A Particle Bewteen The Two Boxes
  
  double Xn, Yn, Zn, Enn, Virn, Eno, Viro; 
  double Arg, Vola, Vold; 
  double Xo, Yo, Zo, Dele;
  int IndexAdd, Iadd, Idel, Jb, Idi;
  
  (*Attempt)++;
  
  // ===Select A Box At Random
  
  if(RandomNumber()<0.5)
    {
      Iadd = 0;
      Idel = 1;
   }
  else
   {
     Iadd = 1;
     Idel = 0;
   }
  
  Vola = pow(Box[Iadd],3);
  Vold = pow(Box[Idel],3);
  
  // ---Add A Particle To Box Iadd
  
  Xn = Box[Iadd]*RandomNumber();
  Yn = Box[Iadd]*RandomNumber();
  Zn = Box[Iadd]*RandomNumber();
  
  // ---Calculate Energy Of This Particle
  
  Jb = 0;
  IndexAdd = Npart;
  
  Eneri(Xn, Yn, Zn, IndexAdd, Jb, &Enn, &Virn, Iadd);
  
  // ---Calculate Contibution To The Chemical Potential:
  Arg = -Beta*Enn;
  Chp[Iadd] = Chp[Iadd] + Vola*exp(Arg)/(double)(Npbox[Iadd]+1);
  if(Npbox[Iadd]==Npart) 
    {
      Chp[Iadd]+=Vola*exp(Arg)/(double)(Npbox[Iadd]+1);
    }
  Ichp[Iadd]++;
  
  // ---Delete Particle From Box B:  
  if(Npbox[Idel]==0)
    { 
      return(0);
    }     
  Idi = -1;
  while(Idi!=Idel)
    {
      IndexAdd = (int)((double)(Npart)*RandomNumber());
      Idi = Id[IndexAdd];
    }
  Xo = X[IndexAdd];
  Yo = Y[IndexAdd];
  Zo = Z[IndexAdd];
  Eneri(Xo, Yo, Zo, IndexAdd, Jb, &Eno, &Viro, Idel);
  
  // ---Acceptance Test:
  Dele = Enn - Eno + log(Vold*(double)((Npbox[Iadd]+1))/
			 (Vola*(double)(Npbox[Idel])))/Beta;
  
  if(RandomNumber()<exp(-Beta*Dele))
    { 
      //    ---Accepted:
      (*Acc)++;
      Npbox[Iadd]++;
      X[IndexAdd] = Xn;
      Y[IndexAdd] = Yn;
      Z[IndexAdd] = Zn;
      Id[IndexAdd] = Iadd;
      En[Iadd]+= Enn;
      Vir[Iadd]+= Virn;
      Npbox[Idel]--;
      En[Idel]-= Eno;
      Vir[Idel]-= Viro;
    }
  return(0);
}
