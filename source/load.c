#include "global.h"
#include "prototypes.h"

void loadconfig_equil(FILE *equilfile, struct disc *particle, long npart, struct disp *trans,
        struct vector box) {
    long i;
    double xpos, xpos_scaled;
    double ypos, ypos_scaled;
    char optionname[30];
    long tmp_id;
    long type;
    int headerln;


    // equilfile = fopen("config-equil", "r");
    // if (equilfile == NULL) die ("Could not open config-equil");
    //get equilibrated translational step size
    read_line(equilfile);
    if(!get_string(optionname, sizeof(optionname))) die ("Error reading equilibrated step size");
    read_line(equilfile);
    if(!get_double(&trans[0].mx)) die ("Error reading equilibrated step size value");
    if(!get_double(&trans[1].mx)) die ("Error reading equilibrated step size value");
    // fscanf(configfile, "%s %lf", optionname, &trans[0].mx);
    // printf("%s %lf\n",optionname, trans[0].mx);
    for (headerln = 0; headerln < 9; headerln++) {
        if(!read_line(equilfile)){
            die ("Could not read header line.");
        }
    }

    for(i = 0; i<npart;i++) {
        
        if(read_line(equilfile)){
            //id type x y....
            if(!get_int(&tmp_id)) die ("Could not read particle ID");
            if(!get_int(&type)) die ("Could not read particle type");
            if(!get_double(&xpos)) die ("Could not read particle x position");
            if(!get_double(&ypos)) die ("Could not read particle y position");

            xpos_scaled = xpos/box.x;
            ypos_scaled = ypos/box.y;

            //Some LAMMPS coorinates are out of scaled bounds so reset them to be at boundary
            if (xpos_scaled >= 1.0) {xpos_scaled = 1.0;}
            if (xpos_scaled < 0.0) {xpos_scaled = 0.0;}
            if (ypos_scaled >= 1.0) {ypos_scaled = 1.0;}
            if (ypos_scaled < 0.0) {ypos_scaled = 0.0;}

            particle[i].idx = tmp_id;
            particle[i].pos.x = xpos_scaled;
            particle[i].pos.y = ypos_scaled;
            particle[i].species = 0;
            particle[i].diameter = 1;

            //if particle positions are still out of bounds, throw error
            if (particle[i].pos.x > 0.5 || particle[i].pos.x < -0.5){
		        printf("Particle %ld: x position out of bounds\n", i+1);
		        printf("Particle: %ld, x: %lf, y: %lf\n", i+1, particle[i].pos.x, particle[i].pos.y);
		        die("Particle x position outside of scaled coordinates bounds");
            }
	        if (particle[i].pos.y > 0.5 || particle[i].pos.y < -0.5) {
		        printf("Particle %ld: y position out of bounds\n", i+1); 
		        printf("Particle: %ld, x: %lf, y: %lf\n", i+1, particle[i].pos.x, particle[i].pos.y);
		        die("Particle y position outside of scaled coordinates bounds");
            }

        } else {
            die("Configuration could not be read in. Load.c: ln62");
        }
        
    }
    // fclose(configfile);
    fclose(equilfile); 
}


void load_potential(double **extPotential, long extPotentialLength, double *extPotentialSpacing, double *extPotentialR0) {
    long i;

    FILE *infile;

    infile = fopen("potential.dat", "r");
    if (!infile) die ("Could not open file potential.dat.");

    for (i = 0; i < extPotentialLength; i++) {
        if (read_line(infile)) {

            if(!get_double(&extPotential[i][0])) die ("Could not read position in potential.dat");
            if(!get_double(&extPotential[i][1])) die ("Could not read potential in potential.dat");
        } else {
            die ("Could not read in potential from potential.dat");
        }
    }

    *extPotentialR0 = extPotential[0][0];
    *extPotentialSpacing = extPotential[1][0] - *extPotentialR0;

}