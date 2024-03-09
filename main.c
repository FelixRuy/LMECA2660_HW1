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


double true_u(double x, double t, double sigma, double c, double U, double L, double INIT_COND) {
    double w_x = fmod(x - c * t, L);
    w_x -= (w_x > 0.5 * L) ? L : 0;  // Ensures w_x is within [-L/2, L/2]
    w_x += (w_x < -0.5 * L) ? L : 0; // Ensures w_x is within [-L/2, L/2]

    if(INIT_COND == 0) return U * exp(-pow(w_x, 2) / (pow(sigma, 2)));
    if(INIT_COND == 1){
        double kp = 2*M_PI/L * 12;
        return U * cos(kp*w_x) * exp(-pow(w_x, 2) / (pow(sigma, 2)));
    }
    return 0.;
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
        Rh += (u[i] - u_sol[i])* (u[i] - u_sol[i]);
    }
    Ih = Ih * h / (sig*U);
    Eh = Eh * h / (sig*U*U);
    Rh = Rh * h / (sig*U*U);

    if(time == 0.25) printf(GRN "Diagnostic for ct/L = 1/4 : Ih = %f, Eh = %f, Rh = %.16f\n" RESET, Ih, Eh, Rh);
    if(time == 0.5) printf(GRN "Diagnostic for ct/L = 1/2 : Ih = %f, Eh = %f, Rh = %.16f\n" RESET, Ih, Eh, Rh);
    if(time == 1) printf(GRN "Diagnostic for ct/L = 1 : Ih = %f, Eh = %f, Rh = %.16f\n" RESET, Ih, Eh, Rh);

    fprintf(file, "(time, Ih, Eh, Rh) : %.10f, %.10f, %.10f, %.10f", time, Ih, Eh, Rh);
    fprintf(file, "\n");
    fclose(file);
}


int main(int argc, char ** argv){
    if (argc < 9){
        printf("Usage: \n"
			"./exec <CFL> <RESOL> <method> <filename_solution> <filename_diagnostic>\n" 
			"---------------------------- \n\n"
			"- CFL : time step condition \n"
			"- resol : h/sigma for the periodic domain. \n"
            "- method : 'E2, E4', 'I4' or 'ED' \n"
            "- MODE : 0 for uniform grid, 1 for non-uniform grid. \n"
            "- INIT_COND : 0 for gaussian, 1 for wave_packet. \n"
            "- filename_solution : output file with computed solution. \n"
            "- filename_diagnostic : output file with global diagnostic. \n"
            "- filename_solution_phys : output file with computed solution in physical space (if mode == 1). \n"
            "---------------------------- \n\n");
		return -1;
    } 
    //_____________________Uniform grid size_____________________
    int res; // To ensure that there is no error
    //_____________________Parameters____________________________
    char * method = argv[3]; // Method for du/dx
    char * file_comp_sol = argv[6]; // Output file with computed solution
    char * file_comp_sol_phys = argv[8]; // Output file with computed solution in physical space
    char * file_glob_diag = argv[7]; // Output file with global diagnostic
    int MODE = atoi(argv[4]); // uniform (0) or non-uniform (1) grid
    int INIT_COND = atoi(argv[5]); // (0) for gaussian, (1) for wave_packet
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
    double a = 3./5.;

    printf(BLU "Mode Selected :\n–––––––––––––––––––––––\n%s\n\n" RESET, MODE ? "Non-uniform grid" : "Uniform grid");
    printf(BLU "Initial condition Selected :\n–––––––––––––––––––––––\n%s\n\n" RESET, INIT_COND ? "Wave packet" : "Gaussian");
    printf(BLU "Parameters :\n–––––––––––––––––––––––\nCFL = %f\ndom_s = %f\nresol = %f\nL = %f\nU = %f\nc = %f\nsigma = %f\nh = %f\nN = %d\ndt = %f\na = %f \n\n" RESET, CFL, dom_s, resol_param, L, U, c, sigma, h, N, dt, a);
    printf(BLU "Method Selected :\n–––––––––––––––––––––––\n%s\n\n" RESET, method);


    if(MODE==0){
        //_____________________Initial condition______________________
        double * xi = (double *)calloc(N, sizeof(double)); // set grid points (uniform spacing)
        double * u = (double *)calloc(N, sizeof(double)); // set u0
        for (int i = 0; i < N; i++){
            xi[i] = -L/2 + i*h;
            u[i] = true_u(xi[i], 0, sigma, c, U, L, INIT_COND);
        }
        //_____________________Allocations_____________________________
        double ** k = malloc(4 * sizeof(double *));
        for (int i = 0; i<4; i++) k[i] = (double *)calloc(N, sizeof(double));
        //_____________________Update phase___________________________
        double time = 0;
        double * u_sol = (double *)calloc(N, sizeof(double));
        //_____________________Write initial condition_________________
        write_vector_file(u, N, U, time, file_comp_sol);
        for(int i = 0; i < N; i++){
            u_sol[i] = true_u(xi[i], time, sigma, c, U, L, INIT_COND);
        }
        write_vector_file(u_sol, N, U, time, "data/analytical_solution.txt");
        write_diagnostic(u, u_sol, h, sigma, U, time, N, file_glob_diag);
        //_____________________Main loop______________________________
        while(time*c/L < 1){
            //____________________________________________
            res = update_u(u, xi, k, N, h, L, c, dt, a, method, MODE); // Update u with derivatives via RK4 scheme
            //____________________________________________
            time = time + dt; // Update time with time-step
            //____________________________________________
            write_vector_file(u, N, U, time, file_comp_sol); // Write to file computed solution
            //____________________________________________
            for (int i = 0; i < N; i++){
                u_sol[i] = true_u(xi[i], time, sigma, c, U, L, INIT_COND); // Compute analytical solution
            }
            write_vector_file(u_sol, N, U, time, "data/analytical_solution.txt"); // Write to file analytical solution
            //____________________________________________
            write_diagnostic(u, u_sol, h, sigma, U, time, N, file_glob_diag); // Write to file global diagnostic
        }
        write_vector_file(xi, N, L, -1., file_comp_sol); // Write grid points to file
        write_vector_file(xi, N, L, -1., "data/analytical_solution.txt"); // Write grid points to file solution

        free(xi);
        for(int i = 0; i<4; i++) free(k[i]);
        free(u_sol);
        return 0;
    }


    if(MODE==1){
        //_____________________Initial condition______________________
        double * xi = (double *)calloc(N, sizeof(double)); // set grid points (uniform spacing)
        double * x = (double *)calloc(N, sizeof(double)); // set grid points (non-uniform spacing)
        double * v = (double *)calloc(N, sizeof(double)); // set v
        double * u = (double *)calloc(N, sizeof(double)); // set u
        for (int i = 0; i < N; i++){
            xi[i] = -L/2 + i*h;
            x[i] = g(xi[i], L, a);
            v[i] = dg(xi[i], L, a)*true_u(x[i], 0, sigma, c, U, L, INIT_COND);
            u[i] = v[i]/dg(xi[i], L, a);
        }
        //_____________________Allocations_____________________________
        double ** k = malloc(4 * sizeof(double *));
        for (int i = 0; i<4; i++) k[i] = (double *)calloc(N, sizeof(double));
        //_____________________Update phase___________________________
        double time = 0;
        double * v_sol = (double *)calloc(N, sizeof(double));
        double * u_sol = (double *)calloc(N, sizeof(double));
        //_____________________Write initial condition_________________
        write_vector_file(v, N, U, time, file_comp_sol);
        write_vector_file(u, N, U, time, file_comp_sol_phys);
        for(int i = 0; i < N; i++){
            v_sol[i] = dg(xi[i], L, a)*true_u(x[i], time, sigma, c, U, L, INIT_COND);
            u_sol[i] = true_u(x[i], time, sigma, c, U, L, INIT_COND);
        }
        write_vector_file(v_sol, N, U, time, "data/analytical_solution.txt");
        write_vector_file(u_sol, N, U, time, "data/analytical_solution_phys.txt");
        write_diagnostic(v, v_sol, h, sigma, U, time, N, file_glob_diag);

        //_____________________Main loop______________________________
        while(time*c/L < 1){
            //____________________________________________
            res = update_u(v, xi, k, N, h, L, c, dt, a, method, MODE); // Update u with derivatives via RK4 scheme
            //____________________________________________
            time = time + dt; // Update time with time-step
            //____________________________________________
            write_vector_file(v, N, U, time, file_comp_sol); // Write to file computed solution
            write_vector_file(u, N, U, time, file_comp_sol_phys);
            //____________________________________________
            for (int i = 0; i < N; i++){
                u[i] = v[i]/dg(xi[i], L, a);
                v_sol[i] = dg(xi[i], L, a)*true_u(x[i], time, sigma, c, U, L, INIT_COND); // Compute analytical solution
                u_sol[i] = true_u(x[i], time, sigma, c, U, L, INIT_COND);
            }
            write_vector_file(v_sol, N, U, time, "data/analytical_solution.txt"); // Write to file analytical solution
            write_vector_file(u_sol, N, U, time, "data/analytical_solution_phys.txt");
            //____________________________________________
            write_diagnostic(v, v_sol, h, sigma, U, time, N, file_glob_diag); // Write to file global diagnostic for numerical solution
        }
        write_vector_file(xi, N, L, -1., file_comp_sol); // Write grid points to file
        write_vector_file(xi, N, L, -1., "data/analytical_solution.txt"); // Write grid points to file solution
        write_vector_file(x, N, L, -1., "data/analytical_solution_phys.txt"); // Write grid points to file solution
        write_vector_file(x, N, L, -1., file_comp_sol_phys); // Write grid points to file solution
        
        free(xi);
        free(x);
        free(v);
        for(int i = 0; i<4; i++) free(k[i]);
        free(v_sol);
        return 0;
    }
    return -1;
}
