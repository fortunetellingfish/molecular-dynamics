#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

typedef struct ETuple{
    double k;
    double u;
    double e;
} ETuple;

double pbcPosition(double s, double L){
    //applies the periodic BC's to a position
    return s<0 ? fmod(s, L)+L : fmod(s,L);
}

double pbcSeparation(double ds, double L){
    // finds shortest distance between particles by applying PBC's
    if (ds <= -0.5*L){ds += L;}
    else if (ds >= 0.5*L){ds -= L;}
    return ds;
}

void ljp(double uf[], double r2){
    //calculates potential and force per unit distance assoc. w/ Lennard-Jones ptl
    //uf is an array with ptl energy at 0 and f over r at 1
    //r2 is r squared
    // eps = 1, sig = 1
    double oneOverR2 = 1.0/(r2);
    double oneOverR6 = oneOverR2 * oneOverR2 * oneOverR2;
    double oneOverR12 = oneOverR6 * oneOverR6;
    uf[0] = 4.0 * ( oneOverR12 - oneOverR6); //potential
    uf[1] = 24.0/r2 * ( 2.0 * oneOverR12 - oneOverR6); //force over r
}

double computeAcceleration(double ax[], double ay[], double az[], double x[], double y[], double z[], double uf[], double Lx, double Ly, double Lz, int N){
    //compute acceleration
    for (int k=0; k<N; k++){
        ax[k] = 0;
        ay[k] = 0;
        az[k] = 0;
    }

    double pe = 0.0;

    for (int k=0; k<N-1; k++){
        for (int j=k+1; j<N; j++){
            double dx = pbcSeparation(x[k] - x[j], Lx);
            double dy = pbcSeparation(y[k] - y[j], Ly);
            double dz = pbcSeparation(z[k] - z[j], Lz);
            double r2 = dx*dx + dy*dy + dz*dz;
            ljp(uf, r2);
            double fx = uf[1] * dx;
            double fy = uf[1] * dy;
            double fz = uf[1] * dz;

            ax[k] += fx;
            ay[k] += fy;
            az[k] += fz;
            ax[j] -= fx;
            ay[j] -= fy;
            az[j] -= fz;

            pe += uf[0];
           }
    }
    return pe; //returns the potential energy of the system at this time step
}

ETuple verlet_step(double x[], double y[], double z[], double vx[], double vy[], double vz[], double ax[], double ay[], double az[], double uf[], double Lx, double Ly, double Lz, double dt, int N, FILE *pos, int write){
    //take one step in the Verlet Algorithm
    //each x[i], vx[i], etc is associated with ONE particle
    double halfdt = 0.5 * dt;
    double halfdt2 = halfdt * dt;

    double ke = 0.0;

    //loop over all particles
    for (int i=0; i<N; i++){
        x[i] = x[i] + vx[i] * dt + ax[i] * halfdt2;
        y[i] = y[i] + vy[i] * dt + ay[i] * halfdt2;
        z[i] = z[i] + vz[i] * dt + az[i] * halfdt2;

        vx[i] += ax[i] * halfdt;
        vy[i] += ay[i] * halfdt;
        vz[i] += az[i] * halfdt;
    }

    double pe = computeAcceleration(ax, ay, az, x, y, z, uf, Lx, Ly, Lz, N);

    //add new acceleration terms
    for(int i=0; i<N; i++){
        vx[i] += ax[i]*halfdt;
        vy[i] += ay[i]*halfdt;
        vz[i] += az[i]*halfdt;

      	ke += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        x[i] = pbcPosition(x[i], Lx);
        y[i] = pbcPosition(y[i], Ly);
        z[i] = pbcPosition(z[i], Lz);
        if(write==0){
            fprintf(pos, "1 %lf %lf %lf\n", x[i], y[i], z[i]);
        }
    }
    ke = ke * 0.5;

    ETuple energies = {ke, pe, ke+pe};

    return energies; //return the kinetic and potential energies of this step
}

void setVelocities(double vx[], double vy[], double vz[], double initialKE, int N){
    srand48(time(NULL));
    double vxSum = 0.0;
    double vySum = 0.0;
    double vzSum = 0.0;

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
        vx[i] = vx[i]*rescale;
        vy[i] = vy[i]*rescale;
        vz[i] = vz[i]*rescale;
    }
}

void setRectangularLattice(double x[], double y[], double z[], double Lx, double Ly, double Lz, int nx, int ny, int nz){
    double dx = Lx/((double) nx); //distance btw columns
    double dy = Ly/((double) ny); //distance btw rows
    double dz = Lz/((double) nz); //distance btw layers

    for(int ix=0; ix<nx; ++ix){
        for(int iy=0; iy<ny; ++iy){
            for(int iz=0; iz<nz; ++iz){
                int i = ix + ny*(iy + iz*nz);
                x[i] = dx * ((double) ix+0.5);
                y[i] = dy * ((double) iy+0.5);
                z[i] = dz * ((double) iz+0.5);
            }
        }
    }
}

double measureTemp(double ke, int N){
    return (2.0 * ke)/(3.0*((double) N));
}

double calcScaleFactor(double currentTemp, double desiredTemp){
    return sqrt(desiredTemp/currentTemp);
}

void rescaleVelocities(double vx[], double vy[], double vz[], double ke, double temp, int N){
    double tk = measureTemp(ke, N);
    double lambda = calcScaleFactor(tk, temp);
    for (int i=0; i<N; i++){
        vx[i] = vx[i]*lambda;
        vy[i] = vy[i]*lambda;
        vz[i] = vz[i]*lambda;
    }
}

void computeNeighbourList(double x[], double y[], double z[], double Lx, double Ly, double Lz, int numberInList[], int **list, double r2ListCutoff, int N){
    for(int i=0; i<N-1; i++){
        numberInList[i] = 0;
        for(int j=i+1; j<N; j++){
            double dx = pbcSeparation(x[i] - x[j], Lx);
            double dy = pbcSeparation(y[i] - y[j], Ly);
            double dz = pbcSeparation(z[i] - y[j], Lz);
            double r2 = dx*dx + dy*dy + dz*dz;
            if(r2<r2ListCutoff){
                list[i][numberInList[i]] = j;
                numberInList[i]++;
            }
        }
    }
}

void write_lmp_config(double x[], double y[], double z[], double vx[], double vy[], double vz[], int N, double Lx, double Ly, double Lz){
    FILE* fp = fopen("lmp_config.data", "w");
    fprintf(fp, "#initial configuration suitable for lammps\n\n");

    fprintf(fp, "%i atoms\n", N);
    fprintf(fp, "1 atom types\n");
    fprintf(fp, "0 %lf xlo xhi\n", Lx);
    fprintf(fp, "0 %lf ylo yhi\n", Ly);
    fprintf(fp, "0 %lf zlo zhi\n", Lz);
    fprintf(fp, "\n");

    fprintf(fp, "Atoms\n\n");
    for (int i=0; i<N; i++){
        fprintf(fp, "%i 1 %lf %lf %lf\n", i, x[i], y[i], z[i]);
    }

    fprintf(fp, "Velocities\n\n");
    for (int i=0; i<N; i++){
        fprintf(fp, "%i %lf %lf %lf\n", i, vx[i], vy[i], vz[i]);
    }

}

int main(int argc, char **argv){

    if(argc!=10){
        fprintf(stderr, "mdverlet requires 9 args: nx, ny, nz, Lx, Ly, Lz, dt, tmax, initialKE.\n");
        return 1;
    }

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int nz = atoi(argv[3]);
    double Lx = atof(argv[4]);
    double Ly = atof(argv[5]);
    double Lz = atof(argv[6]);
    double dt = atof(argv[7]);
    double tmax = atof(argv[8]);
    double initialKE = atof(argv[9]);

    double *x, *y, *z, *vx, *vy, *vz, *ax, *ay, *az;

    int N = nx*ny*nz;

    x = malloc(N*sizeof(double));
    y = malloc(N*sizeof(double));
    z = malloc(N*sizeof(double));
    vx = malloc(N*sizeof(double));
    vy = malloc(N*sizeof(double));
    vz = malloc(N*sizeof(double));
    ax = malloc(N*sizeof(double));
    ay = malloc(N*sizeof(double));
    az = malloc(N*sizeof(double));

    double uf[2] = {0};

    double t=0.0;

    double totalPEAccumulator=0.0;

    double radius = 0.5;

    setRectangularLattice(x, y, z, Lx, Ly, Lz, nx, ny, nz);
    setVelocities(vx, vy, vz, initialKE, N);
    computeAcceleration(ax, ay, az, x, y, z, uf, Lx, Ly, Lz, N);

    write_lmp_config(x, y, z, vx, vy, vz, N, Lx, Ly, Lz);

    FILE* fp = fopen("test.dat", "w");
    FILE* pos = fopen("pos.xyz", "w");

    fprintf(fp, "# t \t KE \t PE\n");
    fprintf(fp, "%lf \t %lf \t %lf \t %lf\n", t, ((double) N)*initialKE, uf[0], ((double) N)*initialKE+uf[0]);

    fprintf(pos, "%i\n\n", N);

    int steps=0;
    int write;

    while(t<tmax){
        t+=dt;
        steps++;
        write = steps % 10;
        ETuple energies = verlet_step(x, y, z, vx, vy, vz, ax, ay, az, uf, Lx, Ly, Lz, dt, N, pos, write);
        fprintf(fp, "%lf \t %lf \t %lf \t %lf\n", t, energies.k, energies.u, energies.e);
        if (write == 0){
            fprintf(pos, "%i\n\n", N);
        }
    }

    fclose(fp);
    fclose(pos);
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(ax);
    free(ay);
    free(az);

    return 0;
}
