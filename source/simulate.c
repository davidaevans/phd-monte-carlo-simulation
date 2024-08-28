#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

void simulate(long ntot, long *npart, struct vector box, double *length,
   double *response, long nsweeps, long report, long movie, long adjust,
   long record, double shell, struct disp *trans, struct disp *rot, int periodic,
   struct disc *particle, int potential, double kt, double wellrange, double welldepth,
   double **extPotential, long extPotentialLength, double prefactor, double expprefactor,
   double extPotentialR0, double extPotentialSpacing, struct vector fixed_dipole, 
   double dipole_strength, double dipole_cutoff, double ljTrunc)
{
   static double denom;
   static double longer;           /* Length of the longer of the two disc species */
   double shorter;
   static double cutoff_phs, cutoff_wca;
   long cell;           /* Cell-list cell of a given particle */
   long i, j;
   static long ncellx, ncelly;     /* Number of cell-list cells in each direction */
   static long ncells;  /* Total number of cell-list cells */
   static long **neighbour=NULL;   /* List of neighbouring cells for each cell */
   long next_adjust;    /* Next sweep number for adjusting step sizes */
   long next_frame;     /* Next sweep number for dumping a movie frame */
   long oldcell;        /* For saving original cell of test particle */
   long step;           /* Current step number within sweep */
   long sweep;          /* Current sweep number */
   long testp;          /* Particle for trial move */
   static struct disc **cfirst=NULL;    /* Array of pointers to first particle in cell list */
   struct vector vold;  /* For saving original position or direction vector */
   FILE *mf;            /* Handle for movie file */
   FILE *dumpfile;       /* Handle for configs file */
   char dumpfilename[50]; /* Buffer for name of configs file */
   char *potshortname[3];


   *potshortname = "hs";
   if (potential == 1) *potshortname = "phs";
   if (potential == 2) {
	   *potshortname = "wca";
	   checkWCA();
   }
   if (potential == 3) *potshortname = "sw";
   if (potential == 4) *potshortname = "psw";
   if (potential == 5) *potshortname = "ext";
   if (potential == 6) {
      *potshortname = "lj";
      checkLJPotential();
   }
   if (potential == 7) *potshortname = "yuk";
   if (potential == 8) *potshortname = "mix";
   if (potential == 9) *potshortname = "dip";
   if (potential == 10) *potshortname = "sto";
   sprintf(dumpfilename, "configs-%s-%.1lfx%.1lf-npart_%ld.dat", *potshortname, box.x, box.y, ntot);

   /*=== Initialise cell list on first call ===*/
   if (!cfirst) {
      /* Find cell list dimensions and allocate memory */
      longer = length[0];
      if (length[1] > length[0]) { longer = length[1]; }
      cutoff_phs = 50/49 * longer;
      cutoff_wca = pow(2, 1.0/6.0) * longer;
      denom = longer * shell;
      if (cutoff_phs > denom && potential == 1) denom = cutoff_phs;
      if (cutoff_wca > denom && potential == 2) denom = cutoff_wca;
      if (wellrange > denom  && potential == 3) denom = wellrange;
      if (wellrange > denom && potential == 4) denom = wellrange;
      if (potential == 5){
         //max range of potential
         denom = extPotential[extPotentialLength-1][0];
         checkPointential(extPotential, extPotentialLength, extPotentialR0, extPotentialSpacing);
         //die("potential checked");
      }
      if (potential == 6) {
         if (ljTrunc > 0.0){
            denom = ljTrunc;
         } else {
            denom = 5.0*longer;
         }
      }
      if (potential == 7) denom = 3*longer;
      if (potential == 8) {
         if (cutoff_phs > denom) denom = cutoff_phs;
         if (extPotential[extPotentialLength-1][0] > denom) denom = extPotential[extPotentialLength-1][0];
      }
      if (potential == 9 || potential == 10) denom = dipole_cutoff;
      ncellx = (long)(box.x / (denom));
      ncelly = (long)(box.y / (denom));
      ncells = ncellx * ncelly;
      cfirst = (struct disc **)malloc(sizeof(struct disc *) * ncells);
      neighbour = (long **)malloc(sizeof(long *) * ncells);
      for (i=0; i<ncells; i++) {
         neighbour[i] = (long *)malloc(sizeof(long) * 10);
      }
      printf("Creating neighbour list\n");
      fflush(stdout);
      /* Work out neighbouring cells for each cell by pointer */
      /* Interior of box */
      for (i=0; i<ncellx; i++) {
         for (j=0; j<ncelly; j++) {
            neighbour[j*ncellx+i][0] = ((j-1) + (j==0?ncelly:0))*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][1] = ((j-1) + (j==0?ncelly:0))*ncellx + i;
            neighbour[j*ncellx+i][2] = ((j-1) + (j==0?ncelly:0))*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][3] = j*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][4] = j*ncellx + i;
            neighbour[j*ncellx+i][5] = j*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][6] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + ((i-1) + (i==0?ncellx:0));
            neighbour[j*ncellx+i][7] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + i;
            neighbour[j*ncellx+i][8] = ((j+1) - (j==ncelly-1?ncelly:0))*ncellx + ((i+1) - (i==ncellx-1?ncellx:0));
            neighbour[j*ncellx+i][9] = -1;  /* end token */
         }
      }
      
      /* Put the particles in the cells */
      for (i=0; i<ntot; i++) { particle[i].next = particle[i].prev = NULL; }
      for (i=0; i<ncells; i++) { cfirst[i] = NULL; }
      for (i=0; i<ntot; i++) {
         cell = getcell(particle[i].pos, ncellx, ncelly);
         //printf("x: %lf y: %lf\n", particle[i].pos.x, particle[i].pos.y); 
         //printf("cell: %ld\n",cell); 
         particle[i].cell = cell;
         if (cfirst[cell]) cfirst[cell]->prev = &particle[i];
         particle[i].next = cfirst[cell];
         cfirst[cell] = &particle[i];
      }

      printf ("Cell list grid: %ld x %ld\n\n", ncellx, ncelly);
   }

/*
   for (i=0;i<ncells;i++) {
      printf ("cell %ld:", i);
      for (j=0;j<10;j++) {
         printf (" %ld", neighbour[i][j]);
      }
      printf ("\n");
   }
*/

   /*=== Initialise counters etc. ===*/

   next_adjust = adjust;
   next_frame = movie;

   for (i=0; i<=1; i++) {
      trans[i].acc = trans[i].rej = 0;
      rot[i].acc = rot[i].rej = 0;
   }

   if (movie > 0) {
      //mf = fopen("movie", "w");
      dumpfile = fopen(dumpfilename, "a");
   } else {
      mf = NULL;
      dumpfile = NULL;
   }

   //for (i = 0; i < ntot; i++) printf("%ld  x: %lf y: %lf\n", particle[i].idx, particle[i].pos.x, particle[i].pos.y);

   /*=== The simulation ===*/
   for (sweep=1; sweep<=nsweeps; sweep++) {
      // printf("sweep: %ld\n", sweep);
      if(sweep%1000 == 0) printf("sweep: %ld\n",sweep);
      fflush(stdout);
      /*=== One translation and one rotation per particle on average ===*/
      for (step=0; step<ntot; step++) {
         testp = (long)(ran2(&seed) * ntot);

   
         /* Translation step */
         //printf("%ld   1x: %lf 1y: %lf\n", particle[testp].idx, particle[testp].pos.x, particle[testp].pos.y);
         vold = particle[testp].pos;
         //printf("%ld   2x: %lf 2y: %lf\n", testp, vold.x, vold.y);
         oldcell = particle[testp].cell;
         particle[testp].pos.x += (ran2(&seed) - 0.5) * 2.0 * trans[particle[testp].species].mx/box.x;
         particle[testp].pos.y += (ran2(&seed) - 0.5) * 2.0 * trans[particle[testp].species].mx/box.y;
         if (particle[testp].pos.x > 0.5) particle[testp].pos.x-=1.0;
         if (particle[testp].pos.x < -0.5) particle[testp].pos.x+=1.0;
         if (particle[testp].pos.y > 0.5) particle[testp].pos.y-=1.0;
         if (particle[testp].pos.y < -0.5) particle[testp].pos.y+=1.0;
         cell = getcell(particle[testp].pos, ncellx, ncelly);
         if (cell != oldcell) {
            /* Remove particle from original cell */
            if (particle[testp].prev) {
               (particle[testp].prev)->next = particle[testp].next;
            } else {
               cfirst[oldcell] = particle[testp].next;
               if (cfirst[oldcell]) cfirst[oldcell]->prev = NULL;
            }
            if (particle[testp].next) {
               (particle[testp].next)->prev = particle[testp].prev;
            }
            /* Insert particle into new cell */
            particle[testp].cell = cell;
            if (cfirst[cell]) cfirst[cell]->prev = &particle[testp];
            particle[testp].next = cfirst[cell];
            particle[testp].prev = NULL;
            cfirst[cell] = &particle[testp];
         }
         //printf("predecision\n");
	      fflush(stdout);
         if ( decision(vold, oldcell, particle, potential, box, ntot, testp, cfirst, neighbour, kt, wellrange, welldepth, extPotential, extPotentialLength, prefactor, expprefactor, extPotentialR0, extPotentialSpacing, dipole_strength, dipole_cutoff, ljTrunc) ) {
            //printf("postdecision-acc\n");
            /* Reject due to overlap */
            particle[testp].pos = vold;
            trans[particle[testp].species].rej++;
            if (cell != oldcell) {
               /* Remove particle from trial cell, given that it must be the first one in the cell */
               if (particle[testp].next) (particle[testp].next)->prev = NULL;
               cfirst[cell] = particle[testp].next;
               /* Reinsert particle into original cell */
               particle[testp].cell = oldcell;
               if (cfirst[oldcell]) cfirst[oldcell]->prev = &particle[testp];
               particle[testp].next = cfirst[oldcell];
               particle[testp].prev = NULL;
               cfirst[oldcell] = &particle[testp];
            }
         } else {
            //printf("postdecision-acc\n");
            /* Accept */
            trans[particle[testp].species].acc++;
         }
         

      }

      /*=== Adjust step sizes during equilibration ===*/
      if (sweep == next_adjust) {
         maxstep(&trans[0], longer, 0.001);
         maxstep(&trans[1], longer, 0.001);
         next_adjust += adjust;
      }


      /* Writing of movie frame */
      if (sweep == next_frame) {
         //fprintf (mf, "%ld\n", ntot);
         //fprintf (mf, "sweep %ld;  box %.10lf %.10lf 1.0\n", sweep, box.x, box.y);
         //draw (mf, box, ntot, particle);
         dumpconfigs_LAMMPS(dumpfile, particle, box, ntot, sweep);
         //fflush (mf);
         next_frame += movie;
      }

   }
   /**** End of sweeps loop ****/

   //if (movie > 0) fclose (mf);
}

/*..............................................................................*/

/*
Adjusts the maximum displacement within the specified limits and resets the
acceptance counters to zero.
*/
 
void maxstep(struct disp *x, double hi, double lo)
{
   if (RATIO(*x) < 0.5) {
      (*x).mx *= 0.95;
   } else {
      (*x).mx *= 1.05;
   }
 
   if ( (*x).mx > hi ) (*x).mx = hi;
   if ( (*x).mx < lo ) (*x).mx = lo;
 
   (*x).acc = (*x).rej = 0;
}

/*..............................................................................*/

/*
Accumulate a value into the statistics and update the mean and rms values.
*/
void accumulate(struct stat *q, double x)
{
   (*q).sum += x;
   (*q).sum2 += x*x;
   (*q).samples++;
   (*q).mean = (*q).sum / (*q).samples;
   (*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
                        (*q).sum * (*q).sum / (*q).samples / (*q).samples));
}

/*................................................................................*/

/*
Accumulate a batch of values into the statistics and update the mean and rms values.
*/
void baccumulate(struct stat *q, double x, double x2, long n)
{
   if (LONG_MAX - (*q).samples < n) {
      fprintf (stderr, "ERROR: Overflow of sample number in baccumulate\n");
      exit (1);
   }
   (*q).sum += x;
   (*q).sum2 += x2;
   (*q).samples += n;
   if ((*q).samples > 0) {
      (*q).mean = (*q).sum / (*q).samples;
      (*q).rms = sqrt(fabs((*q).sum2 / (*q).samples -
                        (*q).sum * (*q).sum / (*q).samples / (*q).samples));
   }
}

/*................................................................................*/

/*
Debugging tool to check that cell lists are consistent.
*/

void checkcells(long ntot, long ncells, struct disc **cfirst,
   long **neighbour, struct disc *particle)
{
   long i;
   static int *seen=NULL;
   struct disc *test;

   if (!seen) seen = (int *)malloc(ntot * sizeof(int));

   /* Check that all particles appear exactly once in the cell list */
   for (i=0; i<ntot; i++) seen[i]=0;
   for (i=0; i<ncells; i++) {
      test = cfirst[i];
      while (test) {
         seen[test->idx]++;
         test = test->next;
      }
   }
   for (i=0; i<ntot; i++) {
      if (seen[i] == 0) {
         printf ("CELL LIST ERROR: Particle %ld does not appear in any cell list\n", i);
         exit(1);
      }
      if (seen[i] > 1) {
         printf ("CELL LIST ERROR: Particle %ld appears in more than one cell list\n", i);
         exit(1);
      }
   }

   /* Check consistency of prev<->next relations */
   for (i=0; i<ntot; i++) {
      if (particle[i].next) {
         if ( (particle[i].next)->prev != &particle[i] ) {
            printf ("CELL LIST ERROR: i->next->prev != i for i=%ld\n", i);
            exit(1);
         }
      }
      if (particle[i].prev) {
         if ( (particle[i].prev)->next != &particle[i] ) {
            printf ("CELL LIST ERROR: i->prev->next != i for i=%ld\n", i);
            exit(1);
         }
      }
   }
}

/*................................................................................*/


void checkPointential(double **extPotential, long extPotentialLength, double extPotentialR0, double extPotentialSpacing) {
   long i;
   double r, pot;
   FILE *outfile = NULL;

   printf("Checking external potential via interpolation\n");
   fflush(stdout);
   outfile = fopen("interpolated-potential.dat", "w");
   if (!outfile) die ("Could not open file iterpolated-potential.dat");

   r = 0.0;
   while (r < 5.0) {
      
      pot = get_ext_energy(extPotential, extPotentialLength, SQ(r),extPotentialR0, extPotentialSpacing);
      //printf("r: %lf pot: %lf\n", r, pot);
      fprintf(outfile, "%lf %lf\n", r, pot);

      r += 0.001;
   }

   fclose(outfile);
}

void checkLJPotential() {
   double shift, critval, critval2, diameter, energy, r, i;
   FILE *outfile = NULL;

   
   shift = 4 * ( pow((1/SQ(2.5)),6) - pow((1/SQ(2.5)),3) );
   diameter = 1.0;

   outfile = fopen("lj-potential.dat", "w");
   if (!outfile) die ("Could not open lj-potential.dat");

   critval2 = SQ(2.5 * diameter);
   i = 0.01;
   while(i<5.0){
      r = SQ(i);
      if (r <= critval2) {
         energy = 4 * ( pow((1/r),6) - pow((1/r),3) ) - shift;
      } else {
         energy = 0.0;
      }
      fprintf(outfile, "%lf %lf\n", sqrt(r), energy);
      i+= 0.01;
   }
   fclose(outfile);
}