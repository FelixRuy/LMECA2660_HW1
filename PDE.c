#include "PDE.h"
#include "thomas.h"
#include "vector.h"
#include <stdio.h>
#include <string.h>

#define MAG   "\x1B[35m"
#define GRN   "\x1B[32m"
#define RESET "\x1B[0m"

int comp_dxu(double * u, double * dxu, int N, double h, char * method){

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
        double cons = 3/h;
        temp[0] = (u[1]-u[N-1])*cons;  // Periodic domain
        for (int i = 1; i < N-1; i++){
            temp[i] = (u[i+1]-u[i-1])*cons;
        }
        temp[N-1] = (u[0]-u[N-2])*cons;  // Periodic domain
        solve_period_3diag(N, 4., 1., 1., dxu, temp); // Solves implicit from thomas.h
        free(temp);
        return 0;
    }
    printf(MAG "Unknown method : select E2, E4, I4, ED. No update done.\n" RESET);
    return -1;
}

// Update u via RK4 scheme
int update_u(double * u, double ** k, int N, double h, double c, double dt, char * method){

    // temp vector to store the derivatives of rk4 steps
    double * temp = (double *) calloc(N, sizeof(double));

    // computes k1 = dt * -c * du/dx
    comp_dxu(u, k[0], N, h, method);
    for(int i=0; i<N; i++) k[0][i] = -c*dt*k[0][i];

    // computes k2 = dt * -c * d(u+ k1/2)/dx
    for(int i=0; i<N; i++) temp[i] = -c*dt*(u[i] + 0.5*k[0][i]);
    comp_dxu(temp, k[1], N, h, method);

    // computes k3 = dt * -c * d(u+ k2/2)/dx
    for(int i=0; i<N; i++) temp[i] = -c*dt*(u[i] + 0.5*k[1][i]);
    comp_dxu(temp, k[2], N, h, method);

    // computes k4 = dt * -c * d(u+ k3)/dx
    for(int i=0; i<N; i++) temp[i] = -c*dt*(u[i] + k[2][i]);
    comp_dxu(temp, k[3], N, h, method);

    // updates u
    for(int i=0; i<N; i++){
        u[i] = u[i] + (k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i])/6;
    }
    return 0;
}


// Update u via RK4 scheme
int update_u_bad(double * u, double * dx1u, double * dx2u, double * dx3u, double * dx4u, int N, double h, double c, double dt){
    double k1, k2, k3, k4;
    double cdt, cdt2, cdt3, cdt4;
    double * temp = (double *) calloc(N, sizeof(double));
    // Computes the powers of multiplicative const to avoid to do it in the for loop
    cdt = -c*dt; cdt2 = cdt*cdt; cdt3 = cdt2*cdt; cdt4 = cdt3*cdt;

    for (int i=0; i<N; i++){
        k1 = cdt * dx1u[i];
        k2 = k1 + cdt2/2 * dx2u[i];
        k3 = k2 + cdt3/4 * dx3u[i];
        k4 = 2*k3 - k1 + cdt4/4 * dx4u[i];
        temp[i] = u[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
    }
    memcpy(u, temp, N*sizeof(double)); // copy the new u^n+1 to previous u^n
    free(temp);
    return 0;

    // BAD BAD BAD BAD BAD BAD 

    /* if (strcmp(method, "E2")==0){
        for(int i = 0; i<N; i++){
            k1 = cdt*dxu[i];
            k2 = k1 + cdt2/(4*h) * (dxu[(i+1)%N] - dxu[(i-1+N)%N]);
            k3 = k2 + cdt3/(16*h2) * (dxu[(i+2)%N] - 2*dxu[i] + dxu[(i-2+N)%N]);
            k4 = 2*k3 - k1 + cdt4/(32*h3) * (dxu[(i+3)%N] - 3*dxu[(i+1)%N] + 3*dxu[(i-1+N)%N] - dxu[(i-3+N)%N]);
            temp[i] = u[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
        }
        memcpy(u, temp, N*sizeof(double)); // copy the new u^n+1 to u
        free(temp);
        return 0;
    }
    if (strcmp(method, "E4")==0){
        for(int i=0; i<N; i++){
            k1 = cdt*(dxu[i]);
            k2 = k1 + cdt2/(24*h) * (dxu[(i-2+N)%N] - dxu[(i+2)%N] + 8*dxu[(i+1)%N] - 8*dxu[(i-1+N)%N]);
            k3 = k2 + cdt3/(576*h2) * (dxu[(i-4+N)%N] - 16*dxu[(i-3+N)%N] + 64*dxu[(i-2+N)%N] + 16*dxu[(i-1+N)%N] - 130*dxu[i] + 16*dxu[(i+1)%N] + 64*dxu[(i+2)%N] - 16*dxu[(i+3)%N] + dxu[(i+4)%N]);
            k4 = 2*k3 - k1 + cdt4/(6912*h3) * (dxu[(i-6+N)%N] - 24*dxu[(i-5+N)%N] + 192*dxu[(i-4+N)%N] - 488*dxu[(i-3+N)%N] - 387*dxu[(i-2+N)%N] + 1584*dxu[(i-1+N)%N] - 256*dxu[i] - 1584*dxu[(i+1)%N] + 387*dxu[(i+2)%N] + 488*dxu[(i+3)%N] - 192*dxu[(i+4)%N] + 24*dxu[(i+5)%N] - dxu[(i+6)%N]);
            temp[i] = u[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
        }
        memcpy(u, temp, N*sizeof(double)); // copy the new u^n+1 to u
        free(temp);
        return 0; 
    }
    if (strcmp(method, "ED")==0){
        for(int i=0; i<N; i++){
            k1 = cdt*(dxu[i]);
            k2 = k1 + cdt2/(12*h) * (dxu[(i-2+N)%N] - 6*dxu[(i-1+N)%N] + 3*dxu[i] + 2*dxu[(i+1)%N]);
            k3 = k2 + cdt3/(144*h2) * (dxu[(i-4+N)%N] - 12*dxu[(i-3+N)%N] + 42*dxu[(i-2+N)%N] - 32*dxu[(i-1+N)%N] - 15*dxu[i] + 12*dxu[(i+1)%N] + 4*dxu[(i+2)%N]);
            k4 = 2*k3 - k1 + cdt4/(864*h3) * (dxu[(i-6+N)%N] -18*dxu[(i-5+N)%N] + 117*dxu[(i-4+N)%N] - 318*dxu[(i-3+N)%N] + 279*dxu[(i-2+N)%N] + 90*dxu[(i-1+N)%N] - 177*dxu[i] - 18*dxu[(i+1)%N] + 36*dxu[(i+2)%N] + 8*dxu[(i+3)%N]);
            temp[i] = u[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
        }
        memcpy(u, temp, N*sizeof(double)); // copy the new u^n+1 to u
        free(temp);
        return 0;
    }
    if(strcmp(method, "I4")==0){

        NULL;

        memcpy(u, temp, N*sizeof(double)); // copy the new u^n+1 to u
        free(temp);
        return 0;
    }
    printf(MAG "Unknown method : select E2, E4, I4, ED. No update done.\n" RESET);
    return -1; */
}