#ifndef PDE_H
#define PDE_H

// Computes dxu  using E2,E4,I4,ED methods
int comp_dxu(double * u, double * dxu, int N, double h, char * method);
// Updates u using a RK4 scheme
int update_u(double * u, double ** k, int N, double h, double c, double dt, char * method);
int update_u_bad(double * u, double * dx1u, double * dx2u, double * dx3u, double * dx4u, int N, double h, double c, double dt);

#endif