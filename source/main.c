#include "global.h"
#include "prototypes.h"
long seed;
/*..............................................................................*/

int main()
{
   double fieldk;              /* External field strength */
   double *length;             /* Length of the cylinders */
   double power;               /* Power of L/D in energy function */
   double *response;           /* Field prefactor for each species */
   double shell;               /* Distance between axes for "connected" discs */
   double kt;                  /* kt for use in metropolist algorithm */
   double wellrange;           /* Range of square well potential */
   double welldepth;           /* Depth of square well potential */
   double prefactor;           /* prefactor for Yukawa potential a in V = -a * exp(-bx) */
   double expprefactor;        /* exponential prefactor for Yukawa potential b in V = -a * exp(-bx) */
   double **extPotential;       /* Array of external potential */
   double extPotentialSpacing; /* Spacing of r in external potential */
   double extPotentialR0;      /* Initial value of r in the external potential */
   double ljTrunc;             /* Value to truncate and shift LJ at. If negative, don't truncate and shift - cutoff becomes 5 sigma */
   double r;
   double r2;
   double tmp_engy;
   double dipole_strength, dipole_cutoff;
   int periodic;               /* 1=wrapping criterion, 0=spanning */
   int potential;              /* 0=HS, 1=PHS, 2=WCA, 3=SW, 4=CUS*/
   int loadconfig;             /* 0=generate initial config, 1=load from config.equil file */
   long adjust;                /* Number of sweeps between step size adjustments */
   long equilibrate;           /* Number of equilibration sweeps */
   long movie;                 /* Number of sweeps between movie frames */
   long *npart;                /* Number of particles */
   long nsweeps;               /* Number of pdiscuction sweeps */
   long ntot;                  /* Totl number of particles */
   long record;                /* Number of sweeps between percolation tests */
   long report;                /* Number of sweeps between statistics reports */
   long extPotentialLength;    /* If reading from external potential, how many points to read in */
   long i;                     /* Iterator */
   struct disp *rot;           /* Maximum step in theta */
   struct disp *trans;         /* Maximum (unscaled) translation step */
   struct disc *particle;      /* Configuration of entire system */
   struct vector box;          /* Simulation cell dimensions */
   struct vector fixed_dipole; /* Direction of fixed dipole */
   FILE *equilfile;            /* File to write equilibrated config too */
   FILE *outpot;

   equilfile = fopen("config.equil", "w+");
   if (equilfile == NULL) die ("Could not open config.equil");

   printf ("\nMC Simulation of Discs");
   printf ("\n--------------------------\n\n");


   npart = (long *)malloc( 2 * sizeof(long) );
   length = (double *)malloc( 2 * sizeof(double) );
   response = (double *)malloc( 2 * sizeof(double) );
   trans = (struct disp *)malloc( 2 * sizeof(struct disp) );
   rot = (struct disp *)malloc( 2 * sizeof(struct disp) );

   /* Get user parameters */
   read_options(npart, &ntot, &box, length, &nsweeps, &report, &record, &movie,
      &shell, &fieldk, &equilibrate, &power, response, &periodic, &adjust, trans,
      rot, &potential, &kt, &wellrange, &welldepth, &loadconfig, &extPotentialLength, &prefactor, &expprefactor, &fixed_dipole, &dipole_strength, &dipole_cutoff, &ljTrunc);
   
   /* Set aside memory for the configuration */
   particle = (struct disc *)malloc(ntot * sizeof(struct disc));

   // if using external potential to be read in
   printf("Setting up arrays for external potential\n");   
   fflush(stdout);
   extPotential = (double **) malloc(sizeof(double *) * extPotentialLength);
   for (i=0;i<extPotentialLength; i++){
      // two elements: r V(r)
      extPotential[i] = (double *) malloc(sizeof(double) * 2);
      extPotential[i][0] = extPotential[i][1] = 0.0;
   }
   printf("Loading external potential\n");
   fflush(stdout);
   if (potential==5 || potential ==8) {
      load_potential(extPotential, extPotentialLength, &extPotentialSpacing, &extPotentialR0);
   }
   printf("External Potential read. Initial r: %lf, spacing: %lf\n", extPotentialR0, extPotentialSpacing);

   //if flag is set to generate configuration
   if (!loadconfig){

      // equilfile = fopen("config-equil", "w");
      // if (equilfile == NULL) die ("Could not open config-equil");

      generateconfig(box, ntot, length[0]);

      init_config(npart, box, length, response, particle, fixed_dipole);
      printf ("Read %ld positions and orientations from \"config.init\"\n\n", ntot);

      printf ("Beginning cell set-up and equilibration\n\n");
      fflush(stdout);
      simulate(ntot, npart, box, length, response, equilibrate, 0, 0, adjust, 0,
         shell, trans, rot, periodic, particle, potential, kt, wellrange, welldepth, extPotential, extPotentialLength, 
         prefactor, expprefactor, extPotentialR0, extPotentialSpacing, fixed_dipole, dipole_strength, dipole_cutoff,
         ljTrunc);
      printf ("Equilibrated step sizes for each species:\n");
      printf ("   Translation: %.6le, %.6le\n", trans[0].mx, trans[1].mx);
      printf ("Writing equilibrated configuration to config.equil\n");
      dumpconfigs_equil(equilfile, particle, box, ntot, 0, trans);
      printf ("\n");
      fclose(equilfile);
   // otherwise just read configuration from config.equil
   } else {
      
      printf("Reading equilibrated configuration from \"config-equil\"\n");
      loadconfig_equil(equilfile, particle, *npart, trans, box);
      printf("Equilibrated maximum step sizes: %lf %lf\n\n", trans[0].mx, trans[1].mx);
      
   }

   printf ("Beginning main run\n\n");
   fflush (stdout);
   simulate(ntot, npart, box, length, response, nsweeps, report, movie, 0, record,
      shell, trans, rot, periodic, particle, potential, kt, wellrange, welldepth, extPotential, extPotentialLength, 
      prefactor, expprefactor, extPotentialR0, extPotentialSpacing, fixed_dipole, dipole_strength, dipole_cutoff,
      ljTrunc);
   printf("Acceptance Ratios:\n");
   printf("acc - species 1: %lf\n", RATIO(trans[0]));
   printf("acc - species 2: %lf\n", RATIO(trans[1]));
   printf ("\nDone\n\n");

   return 0;
}

/*..............................................................................*/

/*
Print error message and exit.
*/

void die(char string[])
{
   fprintf (stderr, "\nERROR: %s\n\n", string);
   exit (1);
}

/*................................................................................*/
