#ifndef VECTOR_H
#define VECTOR_H
void mult_vect_cons(double * v, double c, double * res, int N);
void add_vect(double * v1, double * v2, double * res, int N);
void print_vector(double * v, int N);
void write_vector_file(double * v, int N, double time, char * filename);

#endif