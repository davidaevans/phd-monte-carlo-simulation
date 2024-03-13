#include "global.h"
#include "prototypes.h"

void generateconfig(struct vector box, long ntot, double diameter) 
{
    double sx, sy;
    long squaresize;
    long i, j, count;
    FILE *configinit;


    configinit = NULL;
    configinit = fopen("config.init", "w");
    if (configinit == NULL) die("Could not open config.init in generate.c");

    count = 0;
    squaresize = ceil(sqrt(ntot));

    sx = box.x/(squaresize); 
    sy = box.y/(squaresize);

    if (sx < diameter || sy < diameter) die ("Too many particles for square lattice in generate.c");

    for (i = 0; i < squaresize; i++){
        for (j = 0; j < squaresize; j++) {
            if (count < ntot){
                fprintf(configinit,
                "%le %le\n",
                ((i * sx)/box.x),
                ((j * sy)/box.y));
                count++;
            }
        }
    }

    fclose(configinit);

}