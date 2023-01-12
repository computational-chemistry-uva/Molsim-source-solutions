#define MAX_PARTICLES 10000
#define MAX_WAVEVECTOR 22

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN(x,y) ((x)<(y)?(x):(y))

typedef struct
{
  double x;
  double y;
  double z;
} VECTOR;

typedef struct
{
  long double ax;
  long double ay;
  long double az;

  long double bx;
  long double by;
  long double bz;

  long double cx;
  long double cy;
  long double cz;

} MATRIX;

typedef struct
{
  double re;
  double im;
} COMPLEX;


extern VECTOR Positions[MAX_PARTICLES];
extern VECTOR Forces[MAX_PARTICLES];
extern double Charges[MAX_PARTICLES];

extern double Alpha;
extern double Box;

extern int NumberOfParticles;
extern int NumberOfCells;
extern int NumberOfWavevectors;

double ErrorFunctionComplement(double x);
double ErrorFunction(double x);

void Lattice(void);
void Fourierspace(double *Ukappa,double *Uself);
void Realspace(double *Ureal);
