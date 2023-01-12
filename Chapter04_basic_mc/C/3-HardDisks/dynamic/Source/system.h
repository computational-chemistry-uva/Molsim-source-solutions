#include <stdio.h>

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define SQR(x) ((x)*(x))

#define FALSE 0
#define TRUE 1

#define BOXSIZE 10.0
#define MAX_PARTICLES 100

enum{INITIALIZE,SAMPLE,WRITE_RESULTS};

typedef struct
{
  double x;
  double y;
} VECTOR;

extern VECTOR Positions[MAX_PARTICLES];
extern int NumberOfParticles;
extern int PBC;

void WritePdb(FILE *FilePtr);
void Sample(int Option);

