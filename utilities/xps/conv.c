#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#define N 100000

// Convolution with Lorentz and Fermi function,
// compile with:  gcc conv.c -oconv -lm
// modify the N above if needed.
// Use at your own risk.
// Jerzy Goraus (2003)




typedef struct { 
double x,y;
} pT;


pT p1[N],p2[N];
double a_lor,Ef;

int cmp1(pT *t1, pT *t2)
{
double t=t1->x-t2->x;
return (int)(2*t/fabs(t));
}


inline double lor(double x)
{
return (1/(1+x*x*a_lor));

}

inline double fermi(double x)
{
return 1/(1+exp((x-Ef)/0.02569));

}


main (int argc, char **argv)
{
  int m=4,i1;
  double *tabl=malloc(1600);
 FILE *f;
 double dE=0.4,Ef=0,DE=0.2;
  char *buffer=malloc(256);
  if (!((argc==2)||(argc==4)))
  {
  printf("\nconv: convolution with Lorentz and Fermi function\n\
  conv <filename > [{ FWHM Ef }]\nfilename is xy ascii data file,\
   FWHM - Full Width at Half Maximum  default : 0.4eV     \n \
   Ef - Fermi Energy default : 0 eV n\n");exit(0);
  }
  if (argc==4)
  {
  sscanf(argv[2],"%lf",&dE);
  sscanf(argv[3],"%lf",&Ef);
  if (dE>20) {printf("%i value too high\n",m);exit(1);}
  };
  a_lor=4/(dE*dE);
  srand (time (NULL));
  double y,sum;
  
  int i=0,k,n1=0,n2=0;
  f = fopen (argv[1], "r");
  if (f==NULL) { printf("can't open for reading %s\n",argv[1]); exit(1);}
  while (!feof(f))
  {
  fgets(buffer,255,f);
  sscanf(buffer,"%lf %lf",&(p1[n1].x),&(p1[n1].y));
  n1++;
  assert(n1<N);
}
qsort(p1,n1,sizeof(pT),cmp1);
DE=p1[1].x-p1[0].x;

for (i=0; i<n1; i++)
{
for (sum=0,k=0; k<n1; k+=2) // simple minded simpson rule
sum+=DE/3*((p1[k].y*lor(p1[k].x-p1[i].x)*fermi(p1[i].x))+4*(p1[k+1].y*lor(p1[k+1].x-p1[i].x)*fermi(p1[i].x))\
+(p1[k+2].y*lor(p1[k+2].x-p1[i].x)*fermi(p1[i].x)));

printf("%lg %lg\n",p1[i].x,sum);

}
}