#ifndef PDE_H
#define PDE_H

// Computes dxu  using E2,E4,I4,ED methods
int comp_dxu(double * u, double * dxu, int N, double h, char * method, int mode);
// Updates u using a RK4 scheme
int update_u(double * u, double * grid, double ** k, int N, double h, double L, double c, double dt, double a, char * method, int mode);
// mapping
double g(double xi, double L, double a);
// derivative of mapping
double dg(double xi, double L, double a);


#endif