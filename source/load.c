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

void load_2D_potential(FILE *potential_file, double **potential, long nbins_r, long nbins_theta, double *dr, double *dtheta) {
    
    long i, j; // Iterators
    double tmp_r, tmp_theta, tmp_V; // Temporary variables for reading in values
    double r0, r1, t0, t1;
    
    *dr = *dtheta = 0.0;
    r0 = r1 = t0 = t1 = 0.0;

    for (i = 0; i < nbins_r; i++) {
        for (j = 0; j < nbins_theta; j++) {
            if (read_line(potential_file)) {
                if (!get_double(&tmp_r)) die ("Could not read r value of potential");
                if (!get_double(&tmp_theta)) die ("Could not read r value of potential");
                if (!get_double(&tmp_V)) die ("Could not read r value of potential");

                potential[i][j] = tmp_V;

                if (i == 0) {r0 = tmp_r;}
                if (i == 1) {r1 = tmp_r;}
                if (j == 0) {t0 = tmp_theta;}
                if (j == 1) {t1 = tmp_theta;}

            } else {
                die("Error reading line in potential file - load_2D_potential\n");
            }
        }
    }

    *dr = r1 - r0;
    *dtheta = t1 - t0;
    if (*dr < 0.0001) {die ("dr too small");}
    if (*dtheta < 0.1) {die ("dtheta too small");}

    printf("dr: %lf dtheta: %lf\n", *dr, *dtheta);

}