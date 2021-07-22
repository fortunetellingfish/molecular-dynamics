#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct Tuple{
    double pe, f;
} Tuple;

Tuple compute_row(double r, double rc, double alpha){
    //calculate potential energy and force as a fn of r due to the short range ptl proposed by Wang
    //et al. (2019). sigma=1.

    double ri = 1./r;
    double ri2 = ri*ri;
    double rc2 = rc*rc;
    double rcri2 = rc2*ri2;

    Tuple row;

    row.pe = alpha * (ri2 - 1.) * (rcri2 - 1.) * (rcri2 - 1.); // Potential energy
    row.f = 2.*alpha*ri2*ri2*(rcri2 - 1.)*(3.*rcri2 - 2.*rc2 -1.); // force

    return row;
}

int main(int argc, char* argv[]){

    if(argc!=5){
        printf("Require 4 args: L (box length), resolution, rc (cutoff),  filename. Exiting.\n");
        return 1;
    }

    double L = atof(argv[1]); // box length
    double res = atof(argv[2]); // resolution of potential table
    double rc = atof(argv[3]);

    int N = L/(2*res) + 1;

    // get filename from argv[4] directly
    FILE* fp = fopen(argv[4], "w");
    fprintf(fp, "#UNITS: real\n\n");
    fprintf(fp, "colloid\n");
    fprintf(fp, "N %i\n\n", N);

    double alpha = 2 * rc*rc * pow( 3./(2.*(rc*rc - 1.)) , 3);

    double r = 0.;
    int i=0;

    while(r<L/2.){
        r += res;
        i++;
        if (r>rc){
            //write 0 for pe and force for all r beyond rc
            fprintf(fp, "%i %lf 0 0\n", i, r);
        }

        else{
            // calculate potential
            Tuple row = compute_row(r, rc, alpha);
            fprintf(fp, "%i %lf %lf %lf\n", i, r, row.pe, row.f);
        }
    }

    fclose(fp);
    return 0;
}
