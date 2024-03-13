#include "global.h"
#include "prototypes.h"

/*................................................................................*/

/*
Tests for percolation in the x and y directions seperately, returning 1 or 0
(as double precision values) in each component of a vector for percolation and no
percolation, respectively.  Also returns the total number of contacts in the
system.
*/

struct vector percolate(long npart, struct disc *particle, double shell,
   struct vector box, long *ctot,
   struct stat *npc1, struct stat *npc2, struct disc **cfirst,
   long **neighbour)
{
   int analyse;
   long build;
   long explore;
   long n[2];
   long placed;
   long first;
   long i, j, k;
   long ind;
   struct vector link;
   struct vector psep;
   struct vector pcomp;
   static int *done;
   static long **conn=NULL;
   static long *members=NULL;
   static long *nc=NULL;
   static struct vector *phys=NULL;

   /* Allocate memory on first call, then keep it */
   if (!nc) {
      conn = (long **)malloc(sizeof(long *) * npart);
      for (i=0; i<npart; i++) conn[i]=(long *)malloc(sizeof(long) * MAXO);
      done = (int *)malloc(sizeof(int) * npart);
      members = (long *)malloc(sizeof(long) * npart);
      nc = (long *)malloc(sizeof(long) * npart);
      phys = (struct vector *)malloc(sizeof(struct vector) * npart);
   }

   /* Nonpercolating until proven percolating */
   pcomp.x = pcomp.y = 0.0;

   /* Compile lists of contacts for each disc */
   touching(npart, particle, box, shell, nc, conn, ctot, cfirst, neighbour);

   /* Test the cluster to which each particle belongs, unless it has already
      been tested. */

      for (i=0; i<npart; i++) done[i]=0;

      for (first=0; first<npart; first++) {
         if (done[first]) continue;

         /* Construct a unified instance of the cluster */
         /* Put the first particle at the origin */
         members[0] = first;
         phys[first].x = phys[first].y = 0.0;
         done[first] = 1;
         placed = 1;
         n[0] = n[1] = 0;
         n[particle[first].species] = 1;
         explore = 0;
         analyse = 0;
         /* Place touching particles relative to those already placed */
         while (explore < placed) {
            ind = members[explore];
            for (j=0; j<nc[ind]; j++) {
               build = conn[ind][j];
               if (!done[build]) {
                  members[placed] = build;
                  done[build] = 1;
                  link = image(particle[build].pos, particle[ind].pos, box);
                  phys[build].x = phys[ind].x + link.x;
                  phys[build].y = phys[ind].y + link.y;
                  placed++;
                  n[particle[build].species]++;
               }
            }
            explore++;
         }

         /* Now check distances between particles that are bonded when the periodic
            boundary conditions are switched off */

         /* x percolation (turn off just x boundary conditions) */
         if (pcomp.x == 0.0) {  /* Skip test if a previous cluster percolated */
            for (j=0; j<placed; j++) {
               ind = members[j];
               for (i=0; i<nc[ind]; i++) {
                  k = conn[ind][i];
                  psep.x = phys[ind].x - phys[k].x;
                  psep.y = phys[ind].y - phys[k].y;
                  psep.y = psep.y - box.y * anint(psep.y/box.y);
                  if (!overlap(psep, particle[ind].diameter, particle[k].diameter, shell)) {
                     pcomp.x = 1.0;
                     analyse = 1;
                     goto ytest;   /* Skip to test for y-percolation */
                  }
               }
            }
         }

         /* y percolation (turn off just y boundary conditions) */
         ytest:;
         if (pcomp.y == 0.0) {  /* Skip test if a previous cluster percolated */
            for (j=0; j<placed; j++) {
               ind = members[j];
               for (i=0; i<nc[ind]; i++) {
                  k = conn[ind][i];
                  psep.x = phys[ind].x - phys[k].x;
                  psep.y = phys[ind].y - phys[k].y;
                  psep.x = psep.x - box.x * anint(psep.x/box.x);
                  if (!overlap(psep, particle[ind].diameter, particle[k].diameter, shell)) {
                     pcomp.y = 1.0;
                     analyse = 1;
                     goto stoptest;
                  }
               }
            }
         }

         stoptest:;

         if (analyse == 1) {
            /* If a percolating cluster has been found then gather some statistics on it */
            accumulate(npc1, (double)n[0]);
            accumulate(npc2, (double)n[1]);
            if (pcomp.x + pcomp.y > 1.99) goto escape;
         }

      }  /* End of loop over particles as potential starting points for a chain */

      escape:;
      return pcomp;
}

/*................................................................................*/

/*
For each particle, makes a list of other discs connected to it
(i.e., with axes within a distance of "shell").  The list is returned in "conn"
and the number of connections for each disc is returned in "nc".
The total number of connections is returned as ctot.
*/

void touching(long npart, struct disc *particle, struct vector box,
              double shell, long *nc, long **conn, long *ctot,
              struct disc **cfirst, long **neighbour)
{
   long i;
   long *cell;
   struct disc *test;
   struct vector r_cm;

   *ctot = 0;
   for (i=0; i<npart; i++) nc[i]=0;

   /* Loop over all particles */
   for (i=0; i<npart; i++) {
      /* Loop over all cells adjacent to particle */
      cell = &neighbour[particle[i].cell][0];
      while (*cell >= 0) {
         /* Loop over all particles in cell */
         test = cfirst[*cell];
         while (test) {

            if (i != test->idx) {
               r_cm = image(particle[i].pos, test->pos, box);
               if ( overlap(r_cm, particle[i].diameter, test->diameter, shell) ) {
                  if (nc[i] >= MAXO || nc[test->idx] >= MAXO) {
                     fprintf (stderr,
                        "ERROR: Compiled maximum number of overlaps per disc exceeded\n");
                     fprintf (stderr, "MAXO = %d\n", MAXO);
                     exit (99);
                  }
                  conn[i][nc[i]] = test->idx;
                  nc[i]++;
                  (*ctot)++;
               }  /* End of pair overlap test */
            }

         test = test->next;
         }  /* End of loop over particles in adjacent cell */

         cell++;
      }  /* End of loop of adjacent cells */
   }  /* End of loop over all particles */

   /* Correct for fact that each pair is treated twice. */
   *ctot /= 2;

}

/*................................................................................*/
