#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

#define re  2.740970// A
#define fe  3.542906// (ev/A)
#define a  10.068772
#define b  8.586325
#define A 2.833515//eV
#define B 3.302170//eV
#define k  0.172486 //kappa
#define l  0.172486// lambda
#define pe  33.918623//(eV/A)
#define ps  33.918623//(eV/A)
#define Fn0 -5.223983//eV
#define Fn1 -.071260//eV
#define Fn2  0.227405 //eV
#define Fn3  -4.925318 //eV
#define F0  -5.230007 //eV
#define F1  0 //eV
#define F2  0.244155 //eV
#define F3 -0.156871 //eV
#define n  0.307301 // a fitted paramteter
#define Fe  -5.230007 //eV


#define EOK 0
#define ERROR -1


double Z1 (double r) {
   return exp (-a*((r-re)/re));
}

double Z2 (double r) {
   return exp (-b*((r-re)/re));
}

double Z3 (double r) {
   return (1+pow(((r/re-k),20));
}

double Z4 (double r) {
   return (1+pow((r-l*re)/re),20));
}

double phi (double r) {
return ((A*Z1(r)*Z4(r)-B*Z2(r)*Z3(r))/(Z3(r)*Z4(r)));
}

double rho (double r) {
   return fe*exp (-b*((r-re)/re)))/(1+pow((r-re*l)/re),20));
}

double F (double rho_) {
   if (rho_< (.85*pe))     return Fn0
         + Fn1*pow (((rho_-(.85*pe))/(.85*pe)),1)
         + Fn2*pow (((rho_-(.85*pe))/(.85*pe)),2)
         + Fn3*pow (((rho_-(.85*pe))/(.85*pe)),1);
   else if (rho_>= (1.15*pe))
      return (Fe*(1-n*log((rho_/(1.15*pe))))*pow((rho_/(1.15*pe)),n)) ;
   else return Fn0
         + Fn1*pow (((rho_-(.85*pe))/(pe)),1)
         + Fn2*pow (((rho_-(.85*pe))/(pe)),2)
         + Fn3*pow (((rho_-(.85*pe))/(pe)),1);

}

int main (void) {
   int Nr = 10001;
   double rmax = 9.;
   double dr = rmax/(double)Nr;
   int Nrho = 10001;
   double rhomax = rho (0.);
   double drho = rhomax/(double)Nrho;

   int atomic_number = 1;
   double mass = 63.55;
   double lattice_constant = 3.615;
   char lattice_type[] = "FCC";

   int i;

   char LAMMPSFilename[] = "Cu_potentialbyAjoy.eam";
   FILE *LAMMPSFile = fopen (LAMMPSFilename, "w");
   if (!LAMMPSFile) exit (ERROR);

   // Header for setfl format
   fprintf (LAMMPSFile, \
         "#-> LAMMPS Potential File in  <-#\n"\
         "# Mishin Cu EAM1 PRB(2001)63:224106\n"\
         "# Implemented by G. Ziegenhain (2007) gerolf@ziegenhain.com \n"\
         "%d Cu\n"\
         "%d %20.20f %d %20.20f %20.20f\n"\
         "%d %20.20f %20.20f %s\n",
         atomic_number,
         Nrho, drho, Nr, dr, rc,
         atomic_number, mass, lattice_constant, lattice_type);

   // Embedding function
   for (i = 0; i < Nrho; i++)
      fprintf (LAMMPSFile, "%20.20f\n", F ((double)i*drho));
   // Density function
   for (i = 0; i < Nr; i++)
      fprintf (LAMMPSFile, "%20.20f\n", rho ((double)i*dr));
   // Pair potential
   for (i = 0; i < Nr; i++)
      fprintf (LAMMPSFile, "%20.20f\n", phi ((double)i*dr);

   fclose (LAMMPSFile);
   return (EOK);
}
