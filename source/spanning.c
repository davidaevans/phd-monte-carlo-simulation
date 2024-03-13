#include "global.h"
#include "prototypes.h"

/*................................................................................*/

/*
Tests for percolation in the x and y directions seperately, returning 1 or 0
(as double precision values) in each component of a vector for percolation and no
percolation, respectively.  Also returns the total number of contacts in the
system.
*/

struct vector spantest(long npart, struct disc *particle, double shell,
   struct vector box, long *ctot,
   struct stat *npc1, struct stat *npc2, struct disc **cfirst,
   long **neighbour)
{
   int analyse;
   int top, bottom, left, right;   /* Logical variables for touching boundaries */
   long build;
   long explore;
   long n[2];
   long placed;
   long first;
   long i, j;
   long ind;
   struct vector pcomp;
   static int *done;
   static long **conn=NULL;
   static long *members=NULL;
   static long *nc=NULL;

   /* Allocate memory on first call, then keep it */
   if (!nc) {
      conn = (long **)malloc(sizeof(long *) * npart);
      for (i=0; i<npart; i++) conn[i]=(long *)malloc(sizeof(long) * MAXO);
      done = (int *)malloc(sizeof(int) * npart);
      members = (long *)malloc(sizeof(long) * npart);
      nc = (long *)malloc(sizeof(long) * npart);
   }

   /* Nonpercolating until proven percolating */
   pcomp.x = pcomp.y = 0.0;

   /* Compile lists of contacts for each disc */
   touching_nopbc(npart, particle, box, shell, nc, conn, ctot, cfirst, neighbour);

   /* Test the cluster to which each particle belongs, unless it has already
      been tested. */

      for (i=0; i<npart; i++) done[i]=0;

      for (first=0; first<npart; first++) {
         if (done[first]) continue;

         /* Keep track of particles already considered */
         members[0] = first;
         done[first] = 1;
         placed = 1;
         n[0] = n[1] = 0;
         n[particle[first].species] = 1;
         explore = 0;
         analyse = 0;
         top = bottom = left = right = 0;
         /* Place touching particles relative to those already placed */
         while (explore < placed) {
            ind = members[explore];
            for (j=0; j<nc[ind]; j++) {
               build = conn[ind][j];
               if (!done[build]) {
                  members[placed] = build;
                  done[build] = 1;
                  placed++;
                  n[particle[build].species]++;
                  /* Treat box edges as long discs of zero width */
                  /* Test contact with top boundary */
                  if (particle[build].pos.x * box.x 
                     + (particle[build].diameter*shell)/2 >= box.x/2) {right=1;}
                  if (particle[build].pos.x * box.y 
                     - (particle[build].diameter*shell)/2 <= box.x/2) {left=1;}
                  if (particle[build].pos.y * box.y
                     + (particle[build].diameter*shell)/2 >= box.y/2) {top=1;}
                  if (particle[build].pos.y * box.y
                     - (particle[build].diameter*shell)/2 <= box.y/2) {bottom=1;}
               }
            }
            explore++;
         }

         if (left && right) {
            pcomp.x = 1.0;
            analyse = 1;
         }
         if (top && bottom) {
            pcomp.y = 1.0;
            analyse = 1;
         }

         if (analyse == 1) {
            /* If a percolating cluster has been found then gather some statistics on it */
            accumulate(npc1, (double)n[0]);
            accumulate(npc2, (double)n[1]);

            if ( pcomp.x + pcomp.y > 1.99 ) goto escape;
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
In this version, periodic boundary conditions are not applied.
*/

void touching_nopbc(long npart, struct disc *particle,
   struct vector box,  double shell, long *nc, long **conn, long *ctot,
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
               r_cm.x = (particle[i].pos.x - test->pos.x) * box.x;
               r_cm.y = (particle[i].pos.y - test->pos.y) * box.y;
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
