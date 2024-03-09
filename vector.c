#include "vector.h"
#include <stdio.h>
#include <stdlib.h>

void mult_vect_cons(double * v, double c, double * res, int N){
    for (int i = 0; i < N; i++) res[i] = c*v[i];
}

void add_vect(double * v1, double * v2, double * res, int N){
    for (int i = 0; i < N; i++) res[i] = v1[i] + v2[i];
}

void print_vector(double * v, int N){
    printf("[");
    for (int i = 0; i < N-1; i++) printf("%f, ", v[i]);
    printf("%f]\n", v[N-1]);
}

void write_vector_file(double * v, int N, double scale, double time, char * filename){
    FILE * file = fopen
    (filename, "a");
    fprintf(file, "Time = %f : ", time);
    for (int i = 0; i < N; i++){
        fprintf(file, "%f, ", v[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}
