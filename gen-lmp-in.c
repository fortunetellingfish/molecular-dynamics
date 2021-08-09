#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]){

    if (argc!=9){
        printf("Require 8 args: name of dest. file, nSteps, pair_style, pair_coeff, temperature, configuration file, and dump file.\n");
        printf("Make sure any multi-coefficient command arguments are surrounded by quotation marks.\n");
        printf("Exiting.\n");
        return 1;
    }


    char infileName[32];
    strcpy(argv[1], infileName);

    int nSteps = atoi(argv[2]);
    double dt = atof(argv[3]);

    char pair_style[64];
    strcpy(argv[4], pair_style);

    char pair_coeff[32];
    strcpy(argv[5], pair_coeff);

    double temp = atof(argv[6]);

    char config_file[32];
    strcpy(argv[7], config_file);

    char dump_file[32];
    strcpy(argv[8], dump_file);

    FILE* fp = fopen(infileName, "w");

    fprintf(fp, "variable t index %i\n", nSteps);
    fprintf(fp, "units \t lj\n");
    fprintf(fp, "atom_style \t atomic\n");

    fprintf(fp, "\n");

    fprintf(fp, "read_data \t %s\n", config_file);

    fprintf(fp, "timestep \t %lf\n", dt);

    fprintf(fp, "\n");

    fprintf(fp, "mass \t 1 1.0\n");
    fprintf(fp, "pair_style \t %s\n", pair_style);
    fprintf(fp, "pair_coeff \t %s\n", pair_coeff);

    fprintf(fp, "\n");

    fprintf(fp, "neighbour \t 0.3 bin\n");
    fprintf(fp, "neigh_modify \t delay 0 every 40 check no\n");
    fprintf(fp, "\n\n");

    fprintf(fp, "#OUTPUT\n");
    fprintf(fp, "dump \t %s all xyz 50 %s\n", dump_file, dump_file);
    fprintf(fp, "fix \t 1 all nvt temp %lf %lf 10\n", temp, temp);
    fprintf(fp, "thermo \t 1\n");
    fprintf(fp, "thermo_style \t custom elapsed pe ke etotal press temp\n");
    fprintf(fp, "run \t $t\n");

    fclose(fp);

    return 0;
}
