#include "global.h"
#include "prototypes.h"

/*..............................................................................*/

/*
Reads the run parameters from the external file "options".  See INSTRUCTIONS
for a list of keywords.
*/

void read_options(long *npart, long *ntot, struct vector *box, double *length,
   long *nsweeps, long *report, long *record, long *movie, double *shell,
   double *fieldk, long *equilibrate, double *power, double *response,
   int *periodic, long *adjust, struct disp *trans, struct disp *rot, 
   int *potential, double *kt, double *wellrange, double *welldepth, int *loadconfig,
   long *extPotentialLength, double *prefactor, double *expprefactor, struct vector *fixed_dipole, 
   double *dipole_strength, double *dipole_cutoff, double *ljTrunc)
{
   char command[20];
   char option[20];
   char error[200];
   char initoption[9];
   double boxarea;    /* Area of simulation cell */
   double minbox;
   double cpf[2];     /* Core packing fraction of each species */
   double spf[2];     /* Shell packing fraction of each species */
   double dipole_x, dipole_y, dipole_norm; /* read in x/y components of fixed dipole*/
   FILE *infile;


   // /* Prototypes for keywords library */
   // int  read_line(FILE *);
   // int  get_string(char [], int);
   // void upper_case(char []);
   // int  get_int(long int *);
   // int  get_double(double *);


   /*--- 0. Defaults ---*/
   *adjust = 50;                      /* Sweeps between step size adjustments (equilibration) */
   box->x = box->y -1.0;              /* Box dimensions */
   *equilibrate = 100;                /* Number of equilibration sweeps */
   *fieldk = 0.0;                     /* Field strength (see README) */
   length[0] = 0.0;                   /* L/D aspect ratio of cylindrical part of disc species 1 */
   length[1] = 0.0;                   /* L/D aspect ratio of cylindrical part of disc species 2 */
   npart[0] = 0;                      /* Number of discs species 1 */
   npart[1] = 0;                      /* Number of discs species 2 */
   *nsweeps = 100;                    /* Number of accumulation sweeps */
   *periodic = 1;                     /* Wrapping rather than spanning percolation criterion */
   *potential = 0;                    /* Use hard sphere potential by default */
   *power = 0.0;                      /* Power of aspect ratio in energy function */
   *record = 10;                      /* Number of sweeps between percolation tests */
   *report = 100;                     /* Number of sweeps between statistics reports */
   seed = -1;                         /* Random number seed */
   *shell = 1.2;                      /* Distance between axes of connected discs */
   *kt = 1.0;                         /* Default kt = 1 for use in metropolis algorithm */
   trans[0].mx = trans[1].mx = 0.1;   /* Initial maximum translation step sizes */
   rot[0].mx = rot[1].mx = 0.1;       /* Initial maximum rotation step sizes */
   *wellrange = 0;                    /* Range of square well potential */
   *welldepth = 0;                    /* Depth of square well potential */
   *loadconfig = 0;                   /* Default is to generate configurations (0) rather than load from config.equil (1) */
   *extPotentialLength = 0;         /* Number of lines to read in for external potential file */
   *prefactor = 0;
   *expprefactor = 0;
   dipole_x = dipole_y = *dipole_strength = *dipole_cutoff = 0.0;
   *ljTrunc = -1.0;                   /* Value to truncate and shift LJ at. If negative, don't truncate and shift */
   


   /*--- 1. Read in values ---*/

   infile = fopen("options-mc", "r");
   if (infile == NULL) die ("Could not open \"options\" file");
   printf ("Reading run options from the \"options\" file\n\n");

   while ( read_line(infile) ) {
      get_string(command, sizeof(command));
      upper_case(command);

      if (*command == '#') {
         continue;

      } else if (strcmp(command, "ADJUST") == 0) {
         if (!get_int(adjust)) die ("Could not read number of sweeps between step size adjustments after ADJUST");

      } else if (strcmp(command, "BOX") == 0) {
         if (!get_double(&(box->x))) die ("Could not read x box dimension after BOX");
         if (!get_double(&(box->y))) die ("Could not read y box dimension after BOX");

      } else if (strcmp(command, "DEFINITION") == 0) {
         get_string(option, sizeof(option));
         upper_case(option);
         if (strcmp(option, "SPANNING") == 0) {
            *periodic = 0;
         } else if (strcmp(option, "WRAPPING") == 0) {
            *periodic = 1;
         } else {
            sprintf (error, "Unrecognised option after DEFINITION keyword: %s", option);
            die (error);
         }

      } else if (strcmp(command, "EQUILIBRATE") == 0) {
         if (!get_int(equilibrate)) die ("Could not read number of equilibration sweeps after EQUILIBRATE");
      
      } else if (strcmp(command, "INITIAL") == 0) {
         get_string(initoption, sizeof(initoption));
         upper_case(initoption);
         if (strcmp(initoption, "GENERATE") == 0) {
            *loadconfig = 0;
         } else if (strcmp(initoption, "LOAD") == 0) {
            *loadconfig = 1;
         } else {
            sprintf(error, "Unrecognised option after INITIAL keyword: %s", initoption);
            die (error);
         }

      } else if (strcmp(command, "FIELD") == 0) {
         if (!get_double(fieldk)) die ("Could not read field strength after FIELD");

      } else if (strcmp(command, "LENGTH") == 0) {
         if (!get_double(&length[0])) die ("Could not read L/D ratio of species 1 after LENGTH");
         if (!get_double(&length[1])) die ("Could not read L/D ratio of species 2 after LENGTH");

      } else if (strcmp(command, "MAXSTEP") == 0) {
         while ( get_string(option, sizeof(option)) ) {
            upper_case(option);
            if (strcmp(option, "ROT") == 0) {
               if (!get_double(&rot[0].mx)) die ("Could not read max rotation of species 1 after MAXSTEP ROT");
               if (!get_double(&rot[1].mx)) die ("Could not read max rotation of species 2 after MAXSTEP ROT");
            } else if (strcmp(option, "TRANS") == 0) {
               if (!get_double(&trans[0].mx)) die ("Could not read max translation of species 1 after MAXSTEP TRANS");
               if (!get_double(&trans[1].mx)) die ("Could not read max translation of species 2 after MAXSTEP TRANS");
            } else {
               sprintf (error, "Unrecognised option after MAXSTEP keyword: %s", option);
               die (error);
            }
         }
      } else if (strcmp(command, "POTENTIAL") == 0) {
         get_string(option, sizeof(option));
         upper_case(option);
         if (strcmp(option, "HS") == 0) {
            *potential = 0;
         } else if (strcmp(option, "PHS") == 0) {
            *potential = 1;
         } else if (strcmp(option, "WCA") == 0) {
            *potential = 2;
         } else if (strcmp(option, "SW") == 0) {
            *potential = 3;
            if (!get_double(wellrange)) die ("Could not read potential well range after POTENTIAL SW");
            if (!get_double(welldepth)) die ("Could not read potential well depth after POTENTIAL SW"); 
         } else if (strcmp(option, "PSW") == 0) {
            *potential = 4;
            if (!get_double(wellrange)) die ("Could not read potential well range after POTENTIAL PSW");
            if (!get_double(welldepth)) die ("Could not read potential well depth after POTENTIAL PSW");
         } else if (strcmp(option, "EXT") == 0){
            *potential = 5;
            if (!get_int(extPotentialLength)) die ("Could not read external potential length after POTENTIAL EXT.");
         } else if (strcmp(option, "LJ") == 0) {
            *potential = 6;
            if (!get_double(ljTrunc)) die ("Could not read truncation point after POTENTIAL LJ");
         } else if (strcmp(option, "YUK") == 0) {
            *potential = 7;
            if (!get_double(prefactor)) die ("Could not read prefactor after POTENTIAL YUK");
            if (!get_double(expprefactor)) die ("Could not read exponential prefact after POTENTIAL YUK");
         } else if (strcmp(option, "SPLIT") == 0){
            *potential = 8;
            if (!get_int(extPotentialLength)) die ("Could not read external potential length after POTENTIAL EXT.");
         } else if (strcmp(option, "DIPOLE") == 0) {
            *potential = 9;
            if (!get_double(&dipole_x)) die ("Could not read x component of dipole direction after DIPOLE");
            if (!get_double(&dipole_y)) die ("Could not read y component of dipole direction after DIPOLE");
            if (!get_double(dipole_strength)) die ("Could not read dipole strength after DIPOLE");
            if (!get_double(dipole_cutoff)) die ("Could not read dipole cutoff after DIPOLE");
         } else if (strcmp(option, "STOCKMAYER") == 0) {
            *potential = 10;
            if (!get_double(&dipole_x)) die ("Could not read x component of dipole direction after STOCKMAYER");
            if (!get_double(&dipole_y)) die ("Could not read y component of dipole direction after STOCKMAYER");
            if (!get_double(dipole_strength)) die ("Could not read dipole strength after STOCKMAYER");
            if (!get_double(dipole_cutoff)) die ("Could not read dipole cutoff after STOCKMAYER");
         } else{
            sprintf (error, "Unrecognised option after POTENTIAL keyword: %s", option);
            die (error);
         }
      
      } else if (strcmp(command, "MOVIE") == 0) {
         if (!get_int(movie)) die ("Could not read number sweeps between movie frames after MOVIE");

      } else if (strcmp(command, "POWER") == 0) {
         if (!get_double(power)) die ("Could not read power of (L/D) after POWER");

      } else if (strcmp(command, "RECORD") == 0) {
         if (!get_int(record)) die ("Could not read number sweeps between percolation tests after RECORD");

      } else if (strcmp(command, "DISCS") == 0) {
         if (!get_int(&npart[0])) die ("Could not read number of discs of species 1 after DISCS");
         if (!get_int(&npart[1])) die ("Could not read number of discs of species 2 after DISCS");

      } else if (strcmp(command, "SEED") == 0) {
         if (!get_int(&seed)) die ("Could not read random number seed after SEED");
         seed = -abs(seed);
      } else if (strcmp(command, "KT") == 0 ) {
         if (!get_double(kt)) die ("Could not read value of kt after KT");      
      } else if (strcmp(command, "SHELL") == 0) {
         if (!get_double(shell)) die ("Could not read connection distance after SHELL");

      } else if (strcmp(command, "STATISTICS") == 0) {
         if (!get_int(report)) die ("Could not read number sweeps between statistics reports after STATISTICS");

      } else if (strcmp(command, "SWEEPS") == 0) {
         if (!get_int(nsweeps)) die ("Could not read number of sweeps after SWEEPS");

      } else {
         sprintf (error, "Unrecognised keyword: %s", command);
         die (error);
      }
   }

   fclose(infile);


   /*--- 2. Validity checks ---*/

   *ntot = npart[0] + npart[1];
   if (*ntot < 1) {
      die ("The number of discs must be at least 1.");
   }

   if (length[0] < 0.0 || length[1] < 0.0) {
      die ("The L/D ratio for cylinders cannot be negative.");
   }

   if (seed == 0) {
      die ("The random seed must be a negative integer (not zero).");
   }

   if (*shell < 0.0) {
      die ("The shell thickness cannot be negative (but can be less than 1).");
   }

   if (*kt <= 0.0) {
      die ("The value of kt must be greater than 0.0");
   }

   if (*wellrange < 0) {
      die ("Potential well range must be 0 or greater.");
   }

   minbox = length[0] * 2.0 + *shell * 2.0;
   if (length[1] > length[0]) {
      minbox = length[1] * 2.0 + *shell * 2.0;
   }
   if ( (box->x < minbox) || (box->y < minbox) ) {
      die ("Both box lengths must be at least two full discs long for the longer species.");
   }

   if (*report > *nsweeps) *report=*nsweeps;

   response[0] = pow(length[0], *power) * *fieldk / 2.0;
   response[1] = pow(length[1], *power) * *fieldk / 2.0;

   boxarea = (box->x) * (box->y) ;

   cpf[0] = npart[0] * discarea(length[0], 1.0) / boxarea;
   cpf[1] = npart[1] * discarea(length[1], 1.0) / boxarea;
   
   spf[0] = npart[0] * discarea(length[0], *shell) / boxarea;
   spf[1] = npart[1] * discarea(length[1], *shell) / boxarea;


   // Normalise dipole direction
   if (dipole_x == 0 && dipole_y == 0) {
      printf("DIPOLE X: %lf DIPOLE Y: %lf\n", dipole_x, dipole_y);
      // dipole_norm = sqrt(SQ(dipole_x) + SQ(dipole_y));
      // fixed_dipole->x = dipole_x/dipole_norm;
      // fixed_dipole->y = dipole_y/dipole_norm;
      fixed_dipole->x = fixed_dipole->y = 0.0;
   } else {
      printf("DIPOLE X: %lf DIPOLE Y: %lf\n", dipole_x, dipole_y);
      // fixed_dipole->x = fixed_dipole->y = 0.0;
      dipole_norm = sqrt(SQ(dipole_x) + SQ(dipole_y));
      fixed_dipole->x = dipole_x/dipole_norm;
      fixed_dipole->y = dipole_y/dipole_norm;
   }


   /*--- 3. Summarize results on standard output ---*/
   if (*potential == 0) {
      printf (" Potential:                                %s\n", "HARD SPHERE");
   } else if (*potential == 1) {
      printf (" Potential:                                %s\n", "PSEUDO-HARD SPHERE");
   } else if (*potential == 2) {
      printf (" Potential:                                %s\n", "WCA");
   } else if (*potential == 3) {
      printf (" Potential:                                %s\n", "SQUARE WELL");
      printf ("    Potential range:                       %lf\n", *wellrange);
      printf ("    Potential depth:                       %lf\n", *welldepth);
   } else if (*potential == 4) {
      printf (" Potential:                                %s\n", "CUSTOM PHS-SW COMBINATION");
      printf ("    Potential range:                       %lf\n", *wellrange);
      printf ("    Potential depth:                       %lf\n", *welldepth);
   } else if (*potential == 5) {
      printf (" Potential:                                %s\n", "EXTERNAL reading from file potential.dat");
   } else if (*potential == 6) {
      printf (" Potential:                                %s\n", "LENNARD-JONES");
      printf ("    Truncated and shifted at:              %lf\n", *ljTrunc);
   } else if (*potential == 7) {
      printf (" Potential:                                %s\n", "YUKAWA");
      printf ("    Form  : %.5lf * exp (-%.5lf*x)          \n", *prefactor, *expprefactor);
   } else if (*potential == 8) {
      printf (" Potential:                                %s\n", "SPLIT (active) reading from file potential.dat");
   } else if (*potential == 9) {
      printf (" Potential:                                %s\n", "ALIGNED DIPOLE");
      printf ("    Dipole direction:                      (%.6lf, %.6lf)\n", fixed_dipole->x, fixed_dipole->y);
      printf ("    Dipole strength:                       %lf\n", *dipole_strength);
      printf ("    Dipole cutoff:                         %lf\n", *dipole_cutoff);
   } else if (*potential == 9) {
      printf (" Potential:                                %s\n", "ALIGNED STOCKMAYER");
      printf ("    Dipole direction:                      (%.6lf, %.6lf)\n", fixed_dipole->x, fixed_dipole->y);
      printf ("    Dipole strength:                       %lf\n", *dipole_strength);
      printf ("    Dipole cutoff:                         %lf\n", *dipole_cutoff);
   } else {
      printf (" Potential:                                %s\n", "ERROR");
   }
   
   printf (" Simulation cell dimensions:               %.8lf, %.8lf\n",
      box->x, box->y);
   printf (" KT:                                       %.8lf\n", *kt);
   printf (" External field:                           %.10le\n", *fieldk);
   printf (" Power of L/D in energy function:          %.8lf\n", *power);
   printf (" Number of discs:\n");
   printf ("    Species 1:                             %ld\n", npart[0]);
   printf ("    Species 2:                             %ld\n", npart[1]);
   printf ("    Total:                                 %ld\n", *ntot);
   printf (" Cylinder length:\n");
   printf ("    Species 1:                             %.8lf\n", length[0]);
   printf ("    Species 1:                             %.8lf\n", length[1]);
   printf (" Shell thickness for connectivity:         %.10lf\n", *shell);

   if (*shell < 1.0e-6) {
      printf (" ***\n");
      printf (" *** WARNING: The shell thickness is unhygienically low.\n");
      printf (" ***          Recompile with STICKS=0 switch for thin stick limit.\n");
      printf (" ***\n");
   }

   printf (" Core packing fractions:\n");
   printf ("    Species 1:                             %.10lf\n", cpf[0]);
   printf ("    Species 2:                             %.10lf\n", cpf[1]);
   printf ("    Total:                                 %.10lf\n", cpf[0]+cpf[1]);

   printf (" Core + shell packing fractions:\n");
   printf ("    Species 1:                             %.10lf\n", spf[0]);
   printf ("    Species 2:                             %.10lf\n", spf[1]);
   printf ("    Total:                                 %.10lf\n", spf[0]+spf[1]);
   printf (" Number of sweeps:\n");
   printf ("    Equilibration:                         %ld\n", *equilibrate);
   printf ("    Main run:                              %ld\n", *nsweeps);
   printf ("    Between percolation tests:             %ld\n", *record);
   printf ("    Between statistics reports:            %ld\n", *report);
   printf ("    Between step size adjustments (equil): %ld\n", *adjust);
   if (*movie > 0) {
      printf ("    Between movie frames:                  %ld\n", *movie);
   } else {
      printf (" No movie\n");
   }
   printf (" Initial maximum step sizes for each species:\n");
   printf ("    Translation:                           %.6le, %.6le\n", trans[0].mx, trans[1].mx);
   printf ("    Rotation:                              %.6le, %.6le\n", rot[0].mx, rot[1].mx);
   printf (" Random number seed:                       %ld\n", seed);
   printf (" Percolation definition:                   %s\n", *periodic?"WRAPPING":"SPANNING");
   printf ("\n");

}

/*..............................................................................*/

/*
Reads in the initial configuration from the file "config.init".  Each line
contains the two components of the position vector and two components of
the direction vector for a disc.  The direction vector is normalised
after being read in.  There is no check for overlaps at this stage.
*/

void init_config(long *npart, struct vector box, double *length, double *response,
   struct disc *particle, struct vector fixed_dipole)
{
   long i;
   FILE *infile;

   infile = fopen("config.init", "r");
   if (infile == NULL) {
      fprintf (stderr, "\nERROR: Could not open config.init file.\n\n");
      exit (1);
   }

   for (i=0; i<npart[0]+npart[1]; i++) {
      if (fscanf(infile, "%le %le",
         &particle[i].pos.x, &particle[i].pos.y)
         != 2) {
         fprintf (stderr,
            "ERROR: Could not read coordinates for particle %ld.\n\n", i+1);
         exit (1);
      }
      /* Scale position vector to the unit cube */
      particle[i].pos.x -= 0.5;
      particle[i].pos.y -= 0.5;

      /* Species 1 discs must be listed first */
      particle[i].idx = i;
      if (i<npart[0]) {
         particle[i].species = 0;
         particle[i].diameter = length[0];
         particle[i].response = response[0];
      } else {
         particle[i].species = 1;
         particle[i].diameter = length[1];
         particle[i].response = response[1];
      }

      particle[i].dir = fixed_dipole;
      // printf("dir_x: %lf dir_y: %lf\n", particle[i].dir.x, particle[i].dir.y);
   }

   fclose (infile);
   fflush (stdout);
}

/*..............................................................................*/
