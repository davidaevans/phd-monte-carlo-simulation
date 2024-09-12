#include "global.h"
#include "prototypes.h"

/*..............................................................................*/
 
/*
Determines the energy change in trial move for the PHS potential.
Set e=1.
Returns the energy change explicitly
*/

double phs(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, deltaE;
    double diameter2;
    double r;
    double critval, critval2, tmp;

    energyold = energynew = deltaE = 0;

    diameter2 = SQ(diameter);
    tmp = (double) 50/49;
    critval = tmp * diameter;
    critval2 = SQ(critval);
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);
                if (r <= critval2) {
                    
                   energyold += 50 * pow(tmp, 49) * (pow((diameter2/r), 25) - pow((diameter2/r), 24.5)) + 1;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                if (r <= critval2) {
                    energynew += 50 * pow(tmp, 49) * (pow((diameter2/r), 25) - pow((diameter2/r), 24.5)) + 1;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}

double wca(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, deltaE;
    double diameter2;
    double r;
    double critval, critval2, tmp;

    energyold = energynew = deltaE = 0;

    diameter2 = SQ(diameter);
    tmp = pow(2.0,(1.0/6.0));
    critval = tmp * diameter;
    critval2 = SQ(critval);
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);
                if (r <= critval2) {
                    
                    energyold += 4 * (pow((diameter2/r), 6) - pow((diameter2/r), 3)) + 1;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                if (r <= critval2) {
                    energynew += 4 * (pow((diameter2/r), 6) - pow((diameter2/r), 3)) + 1;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}


double squarewell(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double wellrange, double welldepth)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew;
    double r;
    double critval, critval2;

    energyold = energynew = 0;

    critval = wellrange;
    critval2 = SQ(critval);
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);

                //if there is overlap between hard cores, can immediately reject
                //although there shouldn't be overlap in the old configuration,
                //double check just in case
                if (overlap(r_cm, diameter, test->diameter,1)) return NAN;
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);
                if (r <= critval2) {
                    
                    energyold -= welldepth;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //immediately reject if there is hardcore overlap
                if (overlap(r_cm, diameter, test->diameter,1)) {
                    return NAN;
                }
                
                if (r <= critval2) {
                    energynew -= welldepth;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}

double custompotential(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double wellrange, double welldepth)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, diameter2;
    double r;
    double critval, critval2, wellrange2;
    double tmp;

    tmp = 50/49;

    energyold = energynew = 0;

    critval = tmp * diameter; //cutoff between PHS and SW
    critval2 = SQ(critval);
    wellrange2 = SQ(wellrange);
    diameter2 = SQ(diameter);

    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);

                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);

                //printf("r: %lf\n", r);
                // if less than critical value, it's shifted PHS potential
                if (r <= critval2) {
                    energyold += 50 * pow(tmp, 49) * (pow((diameter2/r), 25) - pow((diameter2/r), 24.5)) + 1 - welldepth;
                // otherwise it is the square well potential
                } else {
                    if (r <= wellrange2){
                        energyold -= welldepth;
                    } else {
                        energyold += 0;
                    }
                }

            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //immediately reject if there is hardcore overlap
                // if less than critical value, it's shifted PHS potential
                if (r <= critval2) {
                    energynew += 50 * pow(tmp, 49) * (pow((diameter2/r), 25) - pow((diameter2/r), 24.5)) + 1 - welldepth;
                // otherwise it is the square well potential
                } else {
                    if (r <= wellrange2){
                        energynew -= welldepth;
                    } else {
                        energynew += 0;
                    }
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */ 
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}

double externalpotential(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double **extPotential, long extPotentialLength, double extPotentialR0, double extPotentialSpacing) {

    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, deltaE;
    double r;
    double critval, critval2;
    //int boolflag = 0;

    energyold = energynew = deltaE = 0;

    // Critical value is given by the maximum distance in the loaded potential
    critval = extPotential[extPotentialLength-1][0];
    critval2 = SQ(critval);

    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);
                if (r <= critval2) {
                    //if (r <= 0.09) {printf("%lf\n",sqrt(r));}
                    energyold += get_ext_energy(extPotential, extPotentialLength, r, extPotentialR0, extPotentialSpacing);
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                if (r <= critval2) {
                    //if (r <= 0.09) {printf("%lf\n",sqrt(r));}
                    energynew += get_ext_energy(extPotential, extPotentialLength, r, extPotentialR0, extPotentialSpacing);
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    //if (boolflag==1) printf("DELTA E: %lf\n", energynew-energyold);
    return energynew - energyold;
}


double lj(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double ljTrunc)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, deltaE;
    double r, shift_r;
    double critval2;
    static double shift;
    
    // Set cutoff to be 5 times diameter to match cell size
    critval2 = SQ(5.0 * diameter);
    if (!shift) {
        // If using a truncated version, calculate shift
        if (ljTrunc > 0.0) {
            shift_r = SQ(ljTrunc);
            shift = 4 * ( pow((1/shift_r),6) - pow((1/shift_r),3) );
            critval2 = shift_r;
        } else {
            shift = 0.0;
        }
    }
    // printf("shift: %lf\n", shift);    
    energyold = energynew = deltaE = 0;
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);

                if (r <= critval2) {
                    energyold += 4 * ( pow((1/r),6) - pow((1/r),3) ) - shift;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);

                if (r <= critval2) {
                    energynew += 4 * ( pow((1/r),6) - pow((1/r),3) ) - shift;
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}


double yukawa(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double prefactor, double expprefactor)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, deltaE;
    double r;

    energyold = energynew = deltaE = 0;
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                // Hard core repulsion for r < 1
                if (r < 1){
                    energyold += 1000000.0;
                } else {
                    //printf("r: %lf\n", r);
                    energyold += prefactor * exp(-expprefactor * sqrt(r));
                }

            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                // Hard core repulsion for r < 1
                if (r < 1){
                    energynew += 1000000.0;
                } else {
                    //printf("r: %lf\n", r);
                    energynew += prefactor * exp(-expprefactor * sqrt(r));
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}

/*..............................................................................*/
 
/*
Determines the energy change in trial move for the PHS potential.
Set e=1.
Returns the energy change explicitly
*/

double split_potential(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double **extPotential, long extPotentialLength, double extPotentialR0, double extPotentialSpacing)
{
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew, deltaE;
    double diameter2;
    double r;
    double critvalCont, critvalCont2, tmp;
    double critvalExt, critvalExt2;

    energyold = energynew = deltaE = 0;

    diameter2 = SQ(diameter);
    tmp = (double) 50.0/49.0;
    critvalCont = tmp * diameter;
    critvalCont2 = SQ(critvalCont);

    // Critical value is given by the maximum distance in the loaded potential
    critvalExt = extPotential[extPotentialLength-1][0];
    critvalExt2 = SQ(critvalExt);
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);

                //Continuous potential contribution
                if (r <= critvalCont2) {
                    
                   energyold += 50 * pow(tmp, 49) * (pow((diameter2/r), 25) - pow((diameter2/r), 24.5)) + 1;
                }

                //Discrete potential contribution
                if (r <= critvalExt2) {
                    //if (r <= 0.09) {printf("%lf\n",sqrt(r));}
                    energyold += get_ext_energy(extPotential, extPotentialLength, r, extPotentialR0, extPotentialSpacing);
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);

                //Continuous potential contribution
                if (r <= critvalCont2) {
                    energynew += 50 * pow(tmp, 49) * (pow((diameter2/r), 25) - pow((diameter2/r), 24.5)) + 1;
                }

                //External/Discrete potential contribution
                if (r <= critvalExt2) {
                    //if (r <= 0.09) {printf("%lf\n",sqrt(r));}
                    energynew += get_ext_energy(extPotential, extPotentialLength, r, extPotentialR0, extPotentialSpacing);
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;


}



double aligned_dipole(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double dipole_strength, double dipole_cutoff) {

    // aligned dipoles with wca repulsion
    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew;
    double r;
    double critval, critval2;
    double wca_critval, wca_critval2;
    double diameter2;

    energyold = energynew = 0;

    critval = dipole_cutoff;
    critval2 = SQ(critval);

    wca_critval = pow(2.0,(1.0/6.0)) * diameter;
    wca_critval2 = SQ(wca_critval);
    
    diameter2 =  SQ(diameter);

    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);

                //if there is overlap between hard cores, can immediately reject
                //although there shouldn't be overlap in the old configuration,
                //double check just in case
                if (overlap(r_cm, diameter, test->diameter,1)) return NAN;
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);
                if (r <= critval2) {
                    // WCA component
                    if (r <= wca_critval2){
                        energyold += 4 * (pow((diameter2/r), 6) - pow((diameter2/r), 3)) + 1;
                    }
                    // Dipole component
                    energyold += SQ(dipole_strength) * CUBE(diameter) / CUBE(sqrt(r)) * ( DOT(particle[testp].dir, test->dir) - 3.0/r * DOT(particle[testp].dir, r_cm) * DOT(test->dir, r_cm));
                    // tmp_energy = SQ(dipole_strength) * CUBE(diameter) / CUBE(sqrt(r)) * ( DOT(particle[testp].dir, test->dir) - 3.0/r * DOT(particle[testp].dir, r_cm) * DOT(test->dir, r_cm));
                    // printf("strength: %lf diameter: %lf r: %lf\n", dipole_strength, diameter, sqrt(r));
                    // printf("testp -> x: %lf y: %lf dir_x: %lf dir_y: %lf\n", particle[testp].pos.x, particle[testp].pos.y, particle[testp].dir.x, particle[testp].dir.y);
                    // printf("insrt -> x: %lf y: %lf dir_x: %lf dir_y: %lf\n", test->pos.x, test->pos.y, test->dir.x, test->dir.y);
                    // printf("r_cm  -> x: %lf y: %lf\n", r_cm.x, r_cm.y);
                    // printf("energy: %lf\n", tmp_energy);
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //immediately reject if there is hardcore overlap
                if (overlap(r_cm, diameter, test->diameter,1)) {
                    return NAN;
                }
                
                if (r <= critval2) {

                    // WCA component
                    if (r <= wca_critval2){
                        energynew += 4 * (pow((diameter2/r), 6) - pow((diameter2/r), 3)) + 1;
                    }
                    // Dipole component
                    energynew += SQ(dipole_strength) * CUBE(diameter) / CUBE(sqrt(r)) * ( DOT(particle[testp].dir, test->dir) - 3.0/r * DOT(particle[testp].dir, r_cm) * DOT(test->dir, r_cm));
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;
}

double stockmayer(struct vector vold, long cellold, struct vector vnew, long cellnew, double diameter,
           struct disc **cfirst, long **neighbour, struct vector box, struct disc *particle, 
           long testp, double dipole_strength, double dipole_cutoff) {


    long *cell;
    struct disc *test;
    struct vector r_cm;

    double energyold, energynew;
    double r;
    double critval, critval2;

    energyold = energynew = 0;

    critval = dipole_cutoff;
    critval2 = SQ(critval);
    
    //absolute energy of old config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellold][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                //printf("test: %ld\ntestp: %ld\n", test->idx, testp);
                r_cm = image(vold, test->pos, box);
                //printf("tx: %lf, ty: %lf\n", test->pos.x, test->pos.y);
                //printf("vx: %lf, vy: %lf\n", vold.x, vold.y);
                //printf("px: %lf, py: %lf\n", particle[testp].pos.x, particle[testp].pos.y);
                //printf("dx: %lf, dy: %lf\n", r_cm.x, r_cm.y);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                //printf("r: %lf\n", r);
                if (r <= critval2) {
                    
                    // Lennard Jones contribution
                    energyold += 4 * ( pow((1/r),6) - pow((1/r),3) );

                    // Dipole contribution
                    energyold += SQ(dipole_strength) * CUBE(diameter) / CUBE(sqrt(r)) * ( DOT(particle[testp].dir, test->dir) - 3.0/r * DOT(particle[testp].dir, r_cm) * DOT(test->dir, r_cm));
                    // tmp_energy = SQ(dipole_strength) * CUBE(diameter) / CUBE(sqrt(r)) * ( DOT(particle[testp].dir, test->dir) - 3.0/r * DOT(particle[testp].dir, r_cm) * DOT(test->dir, r_cm));
                    // printf("strength: %lf diameter: %lf r: %lf\n", dipole_strength, diameter, sqrt(r));
                    // printf("testp -> x: %lf y: %lf dir_x: %lf dir_y: %lf\n", particle[testp].pos.x, particle[testp].pos.y, particle[testp].dir.x, particle[testp].dir.y);
                    // printf("insrt -> x: %lf y: %lf dir_x: %lf dir_y: %lf\n", test->pos.x, test->pos.y, test->dir.x, test->dir.y);
                    // printf("r_cm  -> x: %lf y: %lf\n", r_cm.x, r_cm.y);
                    // printf("energy: %lf\n", tmp_energy);
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */    
    //printf("eo: %lf\n", energyold);
    //absolute energy of new config
    /* Loop over all cells adjacent to particle */
    cell = &neighbour[cellnew][0];
    while (*cell >= 0) {
        /* Loop over all particles in cell */
        test = cfirst[*cell];
        while (test) {

            if (testp != test->idx) {
                r_cm = image(vnew, test->pos, box);
                //calculate pair energy contribution
                r = DOT(r_cm, r_cm);
                
                if (r <= critval2) {
                    //Lennard-Jones contribution
                    energynew += 4 * ( pow((1/r),6) - pow((1/r),3) );
                    // Dipole contribution
                    energynew += SQ(dipole_strength) * CUBE(diameter) / CUBE(sqrt(r)) * ( DOT(particle[testp].dir, test->dir) - 3.0/r * DOT(particle[testp].dir, r_cm) * DOT(test->dir, r_cm));
                }
            }

        test = test->next;
        }  /* End of loop over particles in adjacent cell */

        cell++;
    }  /* End of loop of adjacent cells */   
    //printf("en: %lf\n", energynew);
    //return energy difference
    return energynew - energyold;
}




void checkWCA() {
	double r, r2, tmpV, diameter2, tmp, critval, diameter;
	FILE *outfile = NULL;
	
	outfile = fopen("wca-potential.dat", "w");
	if (outfile == NULL) die ("Could not open wca-potential.dat");
	
	diameter = 1.0;
    	tmp = pow(2,(1.0/6.0));
	critval = tmp * diameter;
	r = 0.025;
	diameter2 = SQ(diameter);
	printf("Diameter squared: %lf\n",diameter2);
	while (r < 5) {
		r2=SQ(r);
		if (r < critval) {
        		tmpV = 4 * (pow((diameter2/r2), 6) - pow((diameter2/r2), 3)) + 1;
		} else {
			tmpV = 0.0;
		}
		fprintf(outfile, "%lf %lf\n", r, tmpV);
		r += 0.025;	
	}
	fclose(outfile);
}
