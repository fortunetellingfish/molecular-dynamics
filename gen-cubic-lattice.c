#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void setCubicLattice(double x[], double y[], double z[], double Lx, double Ly, double Lz, int nx, int ny, int nz){
    double dx = Lx/((double) nx); //distance btw columns
    double dy = Ly/((double) ny); //distance btw rows
    double dz = Lz/((double) nz); //distance btw layers

    for(int ix=0; ix<nx; ++ix){
        for(int iy=0; iy<ny; ++iy){
            for(int iz=0; iz<nz; ++iz){
                int i = ix + ny*(iy + iz*nz);
                printf("%i\n", i);
                x[i] = dx * (ix+0.5);
                y[i] = dy * (iy+0.5);
                z[i] = dz * (iz+0.5);
            }
        }
    }
}

void setVelocities(double vx[], double vy[], double vz[], double initialKE, int N){
    double vxSum = 0.0;
    double vySum = 0.0;
    double vzSum = 0.0;

    //srand48(time(NULL)); //seed the RNG

    for(int i=0; i<N; ++i){//assign random initial velocities
        vx[i] = drand48() - 0.5;
        vy[i] = drand48() - 0.5;
        vz[i] = drand48() - 0.5;

        vxSum += vx[i];
        vySum += vy[i];
        vzSum += vz[i];
    }
    double vxcm = vxSum/((double) N); // centre of mass momentum (velocity)
    double vycm = vySum/((double) N);
    double vzcm = vzSum/((double) N);

    for (int i=0; i<N; ++i){
        vx[i] -= vxcm;
        vy[i] -= vycm;
        vz[i] -= vzcm;
    }

    //rescale velocities to get desired initial KE
    double v2sum = 0;
    for(int i=0; i<N; ++i){
        v2sum += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    double kePerParticle = 0.5*v2sum/((double) N);
    double rescale = sqrt(initialKE/kePerParticle);

    for(int i=0; i<N; ++i){
        vx[i] *= rescale;
        vy[i] *= rescale;
        vz[i] *= rescale;
    }
}

void printFile(char name[], double x[], double y[], double z[], double vx[], double vy[], double vz[], double density, int N){
    strcat(name, ".cfg");
    FILE* fp = fopen(name, "w");
    fprintf(fp, "%lf \t %i \t 0\n", density, N); // write the file header

    for (int i=0; i<N; i++){
        fprintf(fp, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf\n", x[i], y[i], z[i], vx[i], vy[i], vz[i]);
    }
    fclose(fp);
}

int main(int argc, char **argv){
    double Lx = 5;
    double Ly = 5;
    double Lz = 5;

    int nx = 5;
    int ny = 5;
    int nz = 5;

    int N = nx*ny*nz;

    double x[N];
    double y[N];
    double z[N];

    double vx[N];
    double vy[N];
    double vz[N];

    char name[32] = "test";

    setCubicLattice(x, y, z, Lx, Ly, Lz, nx, ny, nz);
    setVelocities(vx, vy, vz, 0.4, N);
    printFile(name, x, y, z, vx, vy, vz, 0.109, N);

    return 0;
}
