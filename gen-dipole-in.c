#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{

    if(argc!=8){
        printf("Require 7 args: Filename, dt, L, N, dipoleInt, nSteps, temperature. Exiting.\n");
        return 1;
    }

    char filename[32];
    strcpy(filename, argv[1]);

    double dt = atof(argv[2]);

    double L = atof(argv[3]);
    int N = atoi(argv[4]);

    double dipoleInt = atof(argv[5]);

    int nSteps = atoi(argv[6]);
    double temp = atof(argv[7]);

    double density = N/(L*L*L);


    strcat(filename, ".in");

    FILE* fp = fopen(filename, "w");

    fprintf(fp, "alpha\t\t0.20\n");
    fprintf(fp, "deltaT\t\t%lf\n", dt);
    fprintf(fp, "density\t\t%lf\n", density);
    fprintf(fp, "dipoleInt\t%lf\n", dipoleInt);
    fprintf(fp, "fSpaceLimit\t7\n");
    fprintf(fp, "initUcell\t%lf %lf %lf\n", L, L, L);
    fprintf(fp, "limitRdf\t200\n");
    fprintf(fp, "mInert\t\t0.025\n");
    fprintf(fp, "randSeed\t17\n");
    fprintf(fp, "rangeRdf\t2.5\n");
    fprintf(fp, "sizeHistRdf\t125\n");
    fprintf(fp, "stepAdjustTemp\t250\n");
    fprintf(fp, "stepAvg\t\t1\n");
    fprintf(fp, "stepEquil\t10\n");
    fprintf(fp, "stepLimit\t%i\n", nSteps);
    fprintf(fp, "stepRdf\t\t200\n");
    fprintf(fp, "temperature\t%lf\n", temp);
    fprintf(fp, "initConf\t xyzv\n");
    fprintf(fp, "potSR\t\tcolloid\n");


    fclose(fp);
    return 0;
}
