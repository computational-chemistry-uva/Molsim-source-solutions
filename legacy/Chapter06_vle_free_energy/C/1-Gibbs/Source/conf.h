#include "parameter.h"
extern double X[Npmax],Y[Npmax],Z[Npmax];
extern int Id[Npmax],Npart;
extern int Npbox[2];

//      Common /Conf1/ X(Npmax),Y(Npmax),Z(Npmax),Id(Npmax),Npart
//      Common /Conf2/ Npbox(2)
//    X(I),Y(I),Z(I)    : Position Particle I
//    Id(I)             : Box 1 Or 2
//    Npart             : Actual Number Of Particles
//    Npbox(I)          : Number Of Particles In Box I
