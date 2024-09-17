#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Returns the cell number for a position in reduced coordinated
(on the unit square -0.5<=x<0.5, -0.5<=y<0.5).
*/

long getcell(struct vector pos, long ncellx, long ncelly)
{
   return (long)( (pos.y+0.5) * ncelly) * ncellx + (long)( (pos.x+0.5) * ncellx);
}

/*..............................................................................*/

/*
Normalise a vector to have unit length.  For speed during heavy use, it is
not checked that the supplied vector has non-zero length.
*/

void normalise(struct vector *u)
{
   double tot;

   tot = sqrt( DOT(*u,*u) );

   (*u).x /= tot;
   (*u).y /= tot;
}

/*..............................................................................*/
/*
Returns angle from x axis of given vector
*/
double get_angle(struct vector r_cm) {
    double tmp_angle;
    if (r_cm.x == 0.0) {
        tmp_angle = atan(r_cm.y);
    } else {
        tmp_angle = atan(r_cm.y/r_cm.x);
    }

    // Adjust angle based on vector to give correct value
    if (r_cm.x < 0.0) {tmp_angle += M_PI;}
    else if( r_cm.y < 0.0) {tmp_angle += 2 * M_PI;}
    else {tmp_angle = tmp_angle;}

    return tmp_angle;
}

/*..............................................................................*/
/* 
Return area of disc
*/

double discarea(double length, double shell) 
{
   return M_PI * SQ(length * shell);
}

/*..............................................................................*/

/* 
Return value of potential via linear interpolation 
*/

double interpolate(double r1, double r2, double v1, double v2, double r_target)
{
   double grad, val;

   //printf("x1: %lf x2: %lf y1: %lf y2: %lf\n", r1, r2, v1, v2);

   grad = (v2-v1)/(r2-r1);

   val = v1 + (r_target-r1) * grad;
   //printf("target: %lf pot: %lf\n",r_target,val);
   return val;
}

/*..............................................................................*/

/* 
Calculate the energy at a given point from external potential.
Assume that the potential is equally spaced and the first value for r is equal to the spacing.

r
*/

double get_ext_energy(double **extPotential, long extPotentialLength, double r_squared, double extPotentialR0, double extPotentialSpacing)
{
   long i;
   int success;
   double r;
   double r0, r1, r2, v1, v2;
   success = 0;
   r = sqrt(r_squared); //sqrt(r);
   v1 = v2 = r1 =r2 = 0.0;

   if (r < extPotential[0][0]) return 1000000.0;
   if (r > extPotential[extPotentialLength-1][0]) return 0.0;

   r0 = extPotential[0][0];

   // i = floor(r/r0);
   i = floor((r - extPotentialR0)/extPotentialSpacing);
   
   //deal with special case where i=0 because r=r0
      r1 = extPotential[i][0];
      v1 = extPotential[i][1];
   
   r2 = extPotential[i+1][0];
   v2 = extPotential[i+1][1];

   return interpolate(r1, r2, v1, v2, r);

   // for (i = 0; i < extPotentialLength - 1; i++){
   //    if (r_sqrt >= extPotential[i][0] && r_sqrt < extPotential[i+1][0]) {
   //       r1 = extPotential[i][0];
   //       r2 = extPotential[i+1][0];
   //       v1 = extPotential[i][1];
   //       v2 = extPotential[i+1][1];
   //       success=1;
   //       break;
   //    }
   // }

// if (!success) die ("Separation of particles is not within range of external potential.");

//    return interpolate(r1, r2, v1, v2, r_sqrt);
}
