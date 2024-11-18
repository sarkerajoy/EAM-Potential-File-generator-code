#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
//#include <unistd.h>
//#include <stdarg.h>
#include <string.h>


#define re  2.556191// A
#define fe  1.554467// (ev/A)
#define a  8.163915
#define b  4.354088
#define A 0.349157//eV
#define B 0.484301//eV
#define k  0.297849 //kappa
#define l  0.768450// lambda
#define pe  21.236488//(eV/A)
#define ps  21.236488//(eV/A)
#define Fn0 -2.314969//eV
#define Fn1 -.158956//eV
#define Fn2  0.491668 //eV
#define Fn3  -1.664346 //eV
#define F0  -2.328566 //eV
#define F1  0 //eV
#define F2  0.566203 //eV
#define F3 -0.2540413 //eV
#define n  0.699553 // a fitted paramteter
#define Fe  -2.328566//eV


#define EOK 0
#define ERROR -1


double Z1 (double r) {
   return exp (-a*((r-re)/re));
}

double Z2 (double r) {
   return exp (-b*((r-re)/re));
}

double Z3 (double r) {
   //return 1+ pow((r/re-k),20);
   double temp = r / re - k;
   return 1+pow(temp,20);
}

double Z4 (double r) {
 double temp=r/re-l ;
 return 1+ pow(temp,20);
  
}

double phi (double r) {
return ((A*Z1(r)*Z4(r)-B*Z2(r)*Z3(r))/(Z3(r)*Z4(r)));
}

double rho (double r) {
   double temp1 = (r -re*l) / re;

   return fe* exp(-b * ( (r-re) /re )) / (1 + pow(temp1,20));
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
   int cutoff=5;
   double rmax = 5;
   double dr = rmax/(double)Nr;
   int Nrho = 10001;
   double rhomax =120.3107496;
   double drho = rhomax/(double)Nrho;

   int atomic_number = 2;
   double mass = 63.55;
   double lattice_constant = 3.615;
   char lattice_type[] = "FCC";

   int i;

   char LAMMPSFilename[] = "C-C_potentialbyAjoy.txt";
   FILE *LAMMPSFile = fopen (LAMMPSFilename, "w");
   if (!LAMMPSFile) return 1;

   // Header for setfl format
   fprintf (LAMMPSFile, \
         "#-> LAMMPS Potential File in  <-#\n"\
         "%d W Cu\n"\
         "%d %20.20e %d %20.20e  %d \n"\
         "%d %20.20f %20.20f %s\n",
         atomic_number, 
         Nrho, drho, Nr, dr, cutoff,
         atomic_number , mass , lattice_constant, lattice_type);

   // Embedding function
   //for (i = 0; i < Nrho; i++) 
      //fprintf (LAMMPSFile, "%20.20e\n", F ((double)i*drho));
      
   // Density function
   //for (i = 0; i < Nr; i++) 
     // fprintf (LAMMPSFile, "%20.20e\n", rho ((double)i*dr));
    
    //Pair potential
   for (i = 0; i < Nr; i++)   
      fprintf (LAMMPSFile, "%20.20e\n", phi(i*dr));

   fclose (LAMMPSFile);
   return 0;
}