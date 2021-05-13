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
    //applies the periodic BC's
    return s<0 ? fmod(s, L)+L : fmod(s,L);
}

double pbcSeparation(double ds, double L){
    return fmod(ds, L);
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

ETuple verlet_step(double x[], double y[], double vx[], double vy[], double ax[], double ay[], double uf[], double Lx, double Ly, double dt, int N){
    //take one step in the Verlet Algorithm
    //each x[i], vx[i], etc is associated with ONE particle
    double halfdt = 0.5 * dt;
    double halfdt2 = halfdt * dt;

    double ke;

    //loop over all particles
    for (i=0; i<N; i++){
        x[i] = x[i] + vx[i] * dt + ax[i] * halfdt2;
        y[i] = y[i] + vy[i] * dt + ay[i] * halfdt2;

        vx[i] += ax[i] * halfdt;
        vy[i] += ay[i] * halfdt;
    }
    double pe = computeAcceleration(ax, ay, x, y, uf, Lx, Ly, N);

    //add new acceleration terms
    for(int i=0; i<N; i++){
        vx[i] += ax[i]*halfdt;
        vy[i] += ay[i]*halfdt;

        ke += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i]);
    }

    ETuple energies = {ke, pe, ke+pe};

    return energies; //return the kinetic and potential energies of this step
}

double computeAcceleration(double ax[], double ay[], double x[], double y[], double uf[], double Lx, double Ly, int N){    
    //compute acceleration
    for (int k=0; k<N; k++){
        ax[k] = 0;
        ay[k] = 0;
    }

    double pe = 0.0;

    for (int k=0; k<N-1; k++){
        for (int j=k+1; j<N; j++){
            double dx = pbcSeparation(x[k] - x[j], Lx);
            double dy = pbcSeparation(y[k] - y[j], Ly);
            double r2 = dx*dx + dy*dy;
            ljp(uf, r2);
            double fx = uf[1] * dx;
            double fy = uf[1] * dy;

            ax[k] += fx;
            ay[k] += fy;
            ax[j] -= fx;
            ay[j] -= fy;

            pe += uf[0];
        }
    }

    return pe; //returns the potential energy of the system at this time step
}

void setVelocities(double vx[], double vy[], double initialKE, int N){
    double vxSum = 0.0;
    double vySum = 0.0;

    srand48(time(NULL)); //seed the RNG

    for(int i=0; i<N; ++i){//assign random initial velocities
        vx[i] = drand48() - 0.5;
        vy[i] = drand48() - 0.5;

        vxSum += vx[i];
        vySum += vy[i];
    }
    double vxcm = vxSum/(double) N; // centre of mass momentum (velocity)
    double vycm = vySum/(double) N;
    for (int i=0; i<N; ++i){
        vx[i] -= vxcm;
        vy[i] -= vycm;
    }
    //rescale velocities to get desired initial KE
    double v2sum = 0;
    for(int i=0; i<N; ++i){
        v2sum += vx[i] * vx[i] + vy[i] * vy[i];
    }
    double kePerParticle = 0.5*v2sum/ (double) N;
    double rescale = sqrt(initialKE/kePerParticle);
    for(int i=0; i<N; ++i){
        vx[i] *= rescale;
        vy[i] *= rescale;
    }
}

void setPositions(double x[], double y[], double Lx, double Ly, int N){
    double rMin2 = pow(2.0, 1.0/3.0); //minimum separation distance
    bool overlap;
    srand48(time(NULL));
    for(int i=0; i<N; ++i){
        do{
            overlap = 0;
            x[i] = Lx * drand48();
            y[i] = Ly * drand48();
            int j=0;
            while(j<i && !overlap){
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                if (dx*dx + dy*dy < rMin2){
                    overlap = 1;
                }
                j++;
            }
        } while(overlap);
    }
}

void setRectangularLattice(double x[], double y[], double Lx, double Ly, int nx, int ny){
    double dx = Lx/(double) nx; //distance btw columns
    double dy = Ly/(double) ny; //distance btw rows

    for(int ix=0; ix<nx; ++ix){
        for(int iy=0; iy<ny; ++iy){
            int i = ix + iy*ny;
            x[i] = dx * (ix+0.5);
            y[i] = dy * (iy+0.5);
        }
    }
}

int main(int argc, char *argv[]){

    if(argc!=9){
        fprintf(stderr, "mdverlet requires 8 args: nx, ny, Lx, Ly, dt, tmax, initialKE, and init.\n");
         return 1;
    }

    int nx = argv[1];
    int ny = argv[2];
    double Lx = argv[3];
    double Ly = argv[4];
    double dt = argv[5];
    double tmax = argv[6];
    double initialKE = argv[7];
    char init[5] = argv[8];

    double *x, *y, *vx, *vy, *ax, *ay;
    x = calloc(N,sizeof(double));
    y = calloc(N,sizeof(double));
    vx = calloc(N,sizeof(double));
    vy = calloc(N,sizeof(double));
    ax = calloc(N,sizeof(double));
    ay = calloc(N,sizeof(double));

    double uf[2] = {0};

    double rho = (double) N/(Lx*Ly);

    double t=0.0;

    double totalPEAccumulator=0.0;

    double radius = 0.5;

    N = nx*ny;

    if (strcmp(init, "rect"){
        setRectangularLattice(x, y, Lx, Ly, nx, ny);
    }
    else if(strcmp(init, "rand")){
        setPositions(x, y, Lx, Ly, N);
    }
    else {
        fprintf(stderr, "arg init can have values 'rect' or 'rand' only.\n");
        return 1;
    }

    setVelocities(vx, vy, initialKE, N);
    computeAcceleration(ax, ay, x, y, uf, Lx, Ly, N);

    FILE* fp = fopen("test.dat", w);

    fprintf(fp, "# t \t KE \t PE\n");
    fprintf(fp, "%lf \t %lf \t %lf \t %lf\n", t, initialKE, uf[0], initialKE+uf[0]);

    while(t<tmax){
        t+=dt;
        ETuple energies = verlet_step(x, y, vx, vy, ax, ay, uf, Lx, Ly, i, dt, N);
        fprintf(fp, "%lf \t %lf \t %lf \n", t, energies.k, energies.u, energies.e);
    }
}











