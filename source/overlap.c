#include "global.h"
#include "prototypes.h"

/*..............................................................................*/
 
/*
Determines whether two "2D discs" overlap.  The distance between the axes
below which overlap is defined is supplied as "crit".  For overlap of the cores,
this is 1, but for overlap of conducting shells it may be any non-negative
number.  The vector from the centre of disc 2 to 1, r_cm, (in real,
not scaled coordinates) must be sent as an argument in this version.

Returns 1 if there is an overlap, 0 if not.
*/

int overlap(struct vector r_cm, double diameter1, double diameter2, double shell)
{
   double sep; //scalar separation of two discs
   double thresh; //threshold value for contact of two discs

   sep = DOT(r_cm, r_cm);
   thresh = (diameter1 * shell)/2 + (diameter2 * shell)/2;
   
   if (sep <= SQ(thresh)) {
      //printf("overlap\n");
      return 1;
   } else {
      //printf("no overlap\n");
      return 0;
   }

}

/*..............................................................................*/

/*
Returns the vector pointing from the centre of mass of particle 2 to the
centre of mass of the closest image of particle 1 in 2D.
*/

struct vector image(struct vector r1, struct vector r2, struct vector box)
{
   struct vector r12;

   r12.x = r1.x - r2.x;
   r12.y = r1.y - r2.y;

   r12.x = box.x * (r12.x - anint(r12.x));
   r12.y = box.y * (r12.y - anint(r12.y));

   return r12;
}

/*..............................................................................*/
 
/*
Returns the nearest integer to its argument as a double precision number. e.g.
anint(-0.49) = 0.0 and anint(-0.51) = -1.0. Equivalent to the Fortran intrinsic
ANINT.
*/
 
double anint(double arg)
{
   if (arg < 0) {
      return (double)( (long)(arg-0.5) );
   } else {
      return (double)( (long)(arg+0.5) );
   }
}
                                
/*..............................................................................*/

/*
Determines whether the core of a specified "2D disc" (testp) overlaps
with the core of any other particle.  Returns 1 if an overlap is detected and 0
if not.
*/

int pairo(long ntot, long testp, struct disc *particle,
          struct vector box, struct disc **cfirst, long **neighbour)
{
   long *cell;
   struct disc *test;
   struct vector r_cm;

   /* Loop over all cells adjacent to particle */
   cell = &neighbour[particle[testp].cell][0];
   while (*cell >= 0) {
      /* Loop over all particles in cell */
      test = cfirst[*cell];
      while (test) {

         if (testp != test->idx) {
            //printf("a\n");
            r_cm = image(particle[testp].pos, test->pos, box);
            //compares on hard core overlap - i.e. shell is set to 1
            if ( overlap(r_cm, particle[testp].diameter, test->diameter, 1.0) )
            { return 1; }
         }

      test = test->next;
      }  /* End of loop over particles in adjacent cell */

      cell++;
   }  /* End of loop of adjacent cells */

   return 0;
}

/*................................................................................*/
//return 1 for reject, 0 for accept
int decision(struct vector vold, long oldcell, struct disc *particle,
                    int potential, struct vector box, long ntot, long testp,
                    struct disc **cfirst, long **neighbour, double kt,
                    double wellrange, double welldepth,
                    double **extPotential, long extPotentialLength,
                    double prefactor, double expprefactor,
                    double extPotentialR0, double extPotentialSpacing,
                    double dipole_strength, double dipole_cutoff, double ljTrunc)
{
   long newcell;
   double deltaE;
   double diameter;
   struct vector vnew;

   newcell = particle[testp].cell;
   vnew = particle[testp].pos;
   diameter = particle[testp].diameter;
   deltaE = 0;
   
    
   //for hardsphere only, always move if there's no overlap
   if (potential == 0) {
      //overlap
      if (!pairo(ntot, testp, particle, box, cfirst, neighbour)) {
         //printf("acc\n");
         return 0;
      } else {
         //printf("rej\n");
         return 1;
      }

   //PHS
   } else if (potential == 1) {
      deltaE = phs(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                   box, particle, testp);
      return metropolis(deltaE, kt);

   //WCA
   } else if (potential == 2) {
      deltaE = wca(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                   box, particle, testp);
      return metropolis(deltaE, kt);

   //SQUARE WELL
   } else if (potential == 3) {
      deltaE = squarewell(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                          box, particle, testp, wellrange, welldepth);
      if (isnan(deltaE)) {
         return 1;
      }
      return metropolis(deltaE, kt);
   
   //CUSTOM SQUARE WELL-PHS COMBINATION
   } else if (potential == 4){
      deltaE = custompotential(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                          box, particle, testp, wellrange, welldepth);
      return metropolis(deltaE, kt);

   //LOAD FROM FILE
   } else if (potential==5) {
      deltaE = externalpotential(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
         box, particle, testp, extPotential, extPotentialLength, extPotentialR0, extPotentialSpacing);

      return metropolis(deltaE, kt);

   //LENNARD-JONES
   } else if (potential == 6) {
      deltaE = lj(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                   box, particle, testp, ljTrunc);
      return metropolis(deltaE, kt);

   //YUKAWA
   } else if (potential == 7) {
      deltaE = yukawa(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                   box, particle, testp, prefactor, expprefactor);
      return metropolis(deltaE, kt);

   // SPLIT POTENTIAL
   } else if (potential == 8) {
      deltaE = split_potential(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
         box, particle, testp, extPotential, extPotentialLength, extPotentialR0, extPotentialSpacing);

      return metropolis(deltaE, kt);

   } else if (potential == 9) {
      deltaE = aligned_dipole(vold, oldcell, vnew, newcell, diameter, cfirst, neighbour,
                          box, particle, testp, dipole_strength, dipole_cutoff);
      if (isnan(deltaE)) {
         return 1;
      }

      return metropolis(deltaE, kt);

   } else {
      die("Decision() not made. Potential not selected.");
      return 1;
   }
   return 1;
}


//return 0 for accept, 1 for reject
int metropolis(double deltaE, double kt) 
{
   double accprob;
   if (deltaE <=0) {
      //printf("acc\n");
      return 0;
   } else {
      accprob = exp(-deltaE/kt);
      if (ran2(&seed) <= accprob) {
         //printf("acc: %lf\n", accprob);
         return 0;
      }
      //printf("rej\n");
      return 1;
   }
}
