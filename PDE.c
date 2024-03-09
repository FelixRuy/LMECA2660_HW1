#include "PDE.h"
#include "thomas.h"
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"

#define MAG   "\x1B[35m"
#define GRN   "\x1B[32m"
#define RESET "\x1B[0m"


double g(double xi, double L, double a){
    return xi - a * L / (2* M_PI) * sin(2*M_PI*xi/L);
}

double dg(double xi, double L, double a){
    return 1 - a * cos(2*M_PI*xi/L);
}

int comp_dxu(double * u, double * dxu, int N, double h, char * method, int mode){

    if (strcmp(method, "E2")==0){
        //printf(GRN "E2 method selected.\n" RESET);
        double den = 2*h;
        dxu[0] = (u[1] - u[N-1])/den;  // Periodic domain
        for (int i = 1; i < N-1; i++){
            dxu[i] = (u[i+1] - u[i-1])/den;
        }
        dxu[N-1] = (u[0] - u[N-2])/den; // Periodic domain
        return 0;
    }
    if (strcmp(method, "E4")==0){
        //printf(GRN "E4 method selected.\n" RESET);
        double den = 12*h;
        dxu[0] = (u[N-2] + 8*u[1] - 8*u[N-1] - u[2])/den;  // Periodic domain
        dxu[1] = (u[N-1] + 8*u[2] - 8*u[0] - u[3])/den;  // Periodic domain
        for (int i = 2; i < N-2; i++){
            dxu[i] = (u[i-2] + 8*u[i+1] - 8*u[i-1] - u[i+2])/den;
        }
        dxu[N-2] = (u[N-4] + 8*u[N-1] - 8*u[N-3] - u[0])/den;  // Periodic domain
        dxu[N-1] = (u[N-3] + 8*u[0] - 8*u[N-2] - u[1])/den;  // Periodic domain
        return 0;
    }
    if (strcmp(method, "ED")==0){
        double den = 6*h;
        dxu[0] = (u[N-2] - 6*u[N-1] + 3*u[0] + 2*u[1])/den;  // Periodic domain
        dxu[1] = (u[N-1] - 6*u[0] + 3*u[1] + 2*u[2])/den;  // Periodic domain
        for (int i = 2; i < N-1; i++){
            dxu[i] = (u[i-2] - 6*u[i-1] + 3*u[i] + 2*u[i+1])/den;
        }
        dxu[N-1] = (u[N-3] - 6*u[N-2] + 3*u[N-1] + 2*u[0])/den;  // Periodic domain
        return 0;
    }
    if (strcmp(method, "I4")==0){
        //printf(GRN "I4 method selected.\n" RESET);
        double * temp = (double *) calloc(N, sizeof(double));
        double cons = 3./(4*h);
        temp[0] = (u[1]-u[N-1])*cons;  // Periodic domain
        temp[N-1] = (u[0]-u[N-2])*cons;  // Periodic domain
        for (int i = 1; i < N-1; i++) temp[i] = (u[i+1]-u[i-1])*cons; // Inside domain
        solve_period_3diag(N, 1.0, 1.0/4.0, 1.0/4.0, dxu, temp); // Solves implicit from thomas.h
        free(temp);
        return 0;
    }
    printf(MAG "Unknown method : select E2, E4, I4, ED. No update done.\n" RESET);
    return -1;
}

// Update u via RK4 scheme
int update_u(double * u, double * grid, double ** k, int N, double h, double L, double c, double dt, double a, char * method, int mode){

    // temp vector to store the derivatives of rk4 steps
    double * temp = (double *) malloc(N*sizeof(double));

    // computes k1 = dt * -c * du/dx
    if(mode==0) for (int i = 0; i < N; i++) temp[i] = -c*dt*u[i];
    if(mode==1) for (int i = 0; i < N; i++) temp[i] = -dt*u[i]*c/dg(grid[i], L, a);
    comp_dxu(temp, k[0], N, h, method, mode);

    // computes k2 = dt * -c * d(u+ k1/2)/dx
    if(mode==0) for(int i=0; i<N; i++) temp[i] = -c*dt*(u[i] + 0.5*k[0][i]);
    if(mode==1) for(int i=0; i<N; i++) temp[i] = -dt*(u[i] + 0.5*k[0][i])*c/dg(grid[i], L, a);
    comp_dxu(temp, k[1], N, h, method, mode);

    // computes k3 = dt * -c * d(u+ k2/2)/dx
    if(mode==0) for(int i=0; i<N; i++) temp[i] = -c*dt*(u[i] + 0.5*k[1][i]);
    if(mode==1) for(int i=0; i<N; i++) temp[i] = -dt*(u[i] + 0.5*k[1][i])*c/dg(grid[i], L, a);
    comp_dxu(temp, k[2], N, h, method, mode);

    // computes k4 = dt * -c * d(u+ k3)/dx
    if(mode==0) for(int i=0; i<N; i++) temp[i] = -c*dt*(u[i] + k[2][i]);
    if(mode==1) for(int i=0; i<N; i++) temp[i] = -dt*(u[i] + k[2][i])*c/dg(grid[i], L, a);
    comp_dxu(temp, k[3], N, h, method, mode);

    // updates u
    for(int i=0; i<N; i++){
        u[i] = u[i] + (k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i])/6;
    }
    
    return 0;
}