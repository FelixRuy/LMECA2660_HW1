#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"

#include "vector.h"
#include "PDE.h"

#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define RESET "\x1B[0m" 

double true_u(double x, double t, double sigma, double c, double U, double L){
    return U * exp(-pow((x-c*t), 2) / (2 * pow(sigma, 2)));
}

void write_diagnostic(double * u, double * u_sol, double h, double sig, double U, double time, int N, char * filename){
    FILE * file = fopen
    (filename, "a");
    double Ih = 0.0;
    double Eh = 0.0;
    double Rh = 0.0;
    for (int i = 0; i < N; i++){
        Ih += u[i];
        Eh += u[i]*u[i];
        Rh += pow(u[i] - u_sol[i], 2);
    }
    Ih = Ih * h / (sig*U);
    Eh = Eh * h / (sig*U*U);
    Rh = Rh * h / (sig*U*U);
    fprintf(file, "(time, Ih, Eh, Rh) : %f, %f, %f, %f", time, Ih, Eh, Rh);
    fprintf(file, "\n");
    fclose(file);
}


int main(int argc, char ** argv){
    if (argc < 6){
    printf("Usage: \n"
			"./exec <CFL> <RESOL> <method> <filename_solution> <filename_diagnostic>\n" 
			"---------------------------- \n\n"
			"- CFL : time step condition \n"
			"- resol : h/sigma for the periodic domain. \n"
            "- method : 'E2, E4', 'I4' or 'ED' \n"
            "- filename_solution : output file with computed solution. \n"
            "- filename_diagnostic : output file with global diagnostic. \n"
      "\n");
		return -1;
  } 
    //_____________________Uniform grid size_____________________
    int res; // To ensure that there is no error
    //_____________________Parameters____________________________
    char * method = argv[3]; // Method for du/dx
    char * file_comp_sol = argv[4]; // Output file with computed solution
    char * file_glob_diag = argv[5]; // Output file with global diagnostic
    double CFL = atof(argv[1]); // CFL number for stability
    double resol_param = atof(argv[2]); // Resolution parameter
    double dom_s = 16.; // L/sigma for the periodic domain
    double L = 1.; // Length of the domain
    double U = 1.; // Constant for the gaussian 
    double c = 1.; // propagation velocity
    double sigma = 1.0/dom_s; // Standard deviation of the gaussian
    double h = sigma/resol_param ; // Grid size
    int N = (int) L / h; // Number of grid points
    double dt = CFL * h / c; // Time step
    printf(BLU "Parameters :\n–––––––––––––––––––––––\nCFL = %f\ndom_s = %f\nL = %f\nU = %f\nc = %f\nsigma = %f\nh = %f\nN = %d\ndt = %f \n\n" RESET, CFL, dom_s, L, U, c, sigma, h, N, dt);
    printf(BLU "Method Selected :\n–––––––––––––––––––––––\n%s\n\n" RESET, method);
    //_____________________Initial condition______________________
    double * xi = (double *)calloc(N, sizeof(double)); // set grid points
    double * u = (double *)calloc(N, sizeof(double)); // set u0
    for (int i = 0; i < N; i++){
        xi[i] = -L/2 + i*h;
        u[i] = U * exp(-pow(xi[i], 2)/(2*pow(sigma, 2)));
    }
    //_____________________Allocations_____________________________
    /* double * dx1u = (double *)calloc(N, sizeof(double)); // alloc du/dx
    double * dx2u = (double *)calloc(N, sizeof(double)); // alloc ddu/dx2
    double * dx3u = (double *)calloc(N, sizeof(double)); // alloc dddu/dx3
    double * dx4u = (double *)calloc(N, sizeof(double)); // alloc ddddu/dx4 */
    double ** k = malloc(4 * sizeof(double *));
    for (int i = 0; i<4; i++) k[i] = (double *)calloc(N, sizeof(double));
    //_____________________Update phase___________________________
    double time = 0;
    double * u_sol = (double *)calloc(N, sizeof(double));
    while(time*c/L < 1){
        //____________________________________________
        res = update_u(u, k, N, h, c, dt, method); // Update u with derivatives via RK4 scheme

        /* res = comp_dxu(u, dx1u, N, h, method); // Compute du/dx
        if(time == 0){
            printf(YEL "First derivative computed (over u0) :\n" RESET);
            print_vector(dx1u, N);
        }
        res = comp_dxu(dx1u, dx2u, N, h, method); // Compute ddu/dx2
        res = comp_dxu(dx2u, dx3u, N, h, method); // Compute dddu/dx3
        res = comp_dxu(dx3u, dx4u, N, h, method); // Compute ddddu/dx4
        res = update_u(u, dx1u, dx2u, dx3u, dx4u, N, h, c, dt); // Update u with derivatives via RK4 scheme */
        //____________________________________________
        write_vector_file(u, N, time, file_comp_sol); // Write to file computed solution
        //____________________________________________
        for (int i = 0; i < N; i++){
            u_sol[i] = true_u(xi[i], time, sigma, c, U, L); // Compute analytical solution
        }
        write_vector_file(u_sol, N, time, "data/analytical_solution.txt"); // Write to file analytical solution
        //____________________________________________
        write_diagnostic(u, u_sol, h, sigma, U, time, N, file_glob_diag); // Write to file global diagnostic
        //____________________________________________
        time = time + dt; // Update time with time-step
    }
    write_vector_file(xi, N, -1., file_comp_sol); // Write grid points to file
    write_vector_file(xi, N, -1., "data/analytical_solution.txt"); // Write grid points to file solution

    free(xi);
    for(int i = 0; i<4; i++) free(k[i]);
    free(u_sol);
    return 0;
}
