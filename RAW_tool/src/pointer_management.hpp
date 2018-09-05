#ifndef POINTER_MANAGEMENT_H
#define POINTER_MANAGEMENT_H

/*----
Author: Jérôme Dubois
Contributors:
Version: 1.1
Date: 07/2018
Descritption:
A simple C librairy made to easily manage "vector" and "matrix"
pointers. Just need to call create, resize and delete function.
Usage exemple for int:

int *vec=NULL;
int numel=12;
create_vector_int(&vec, numel);
resize_vector_int(&vec, numel+4);
delete_vector_int(&vec);

successfully pass varied valgrind tests
----*/


//complex typedef
typedef double _Complex cplxd;
typedef float _Complex cplxf;

// char vector management
void create_vector_char (char **data_vec, int N_data);
void create_0_vector_char (char **data_vec, int N_data);
void resize_vector_char (char **data_vec, int N_data);
void delete_vector_char (char **data_vec);
// char matrix management
void create_matrix_char (char ***data_mat, int Nx_data, int Ny_data);
void create_0_matrix_char (char ***data_mat, int Nx_data, int Ny_data);
void resize_matrix_char (char ***data_mat, int Nx_data, int Ny_data, int Nx_old_data, int Ny_old_data);
void delete_matrix_char (char ***data_mat, int Nx_data);

// int vector management
void create_vector_int (int **data_vec, int N_data);
void create_0_vector_int (int **data_vec, int N_data);
void resize_vector_int (int **data_vec, int N_data);
void delete_vector_int (int **data_vec);
// int matrix management
void create_matrix_int (int ***data_mat, int Nx_data, int Ny_data);
void create_0_matrix_int (int ***data_mat, int Nx_data, int Ny_data);
void resize_matrix_int (int ***data_mat, int Nx_data, int Ny_data, int Nx_old_data, int Ny_old_data);
void delete_matrix_int (int ***data_mat, int Nx_data);

// int16_t vector management
void create_vector_int16_t (int16_t **data_vec, int N_data);
void create_0_vector_int16_t (int16_t **data_vec, int N_data);
void resize_vector_int16_t (int16_t **data_vec, int N_data);
// int matrix management
void delete_vector_int16_t (int16_t **data_vec);
void create_matrix_int16_t (int16_t ***data_mat, int Nx_data, int Ny_data);
void create_0_matrix_int16_t (int16_t ***data_mat, int Nx_data, int Ny_data);
void resize_matrix_int16_t (int16_t ***data_mat, int Nx_data, int Ny_data, int Nx_old_data, int Ny_old_data);
void delete_matrix_int16_t (int16_t ***data_mat, int Nx_data);

// float vector management
void create_vector_float (float **data_vec, int N_data);
void create_0_vector_float (float **data_vec, int N_data);
void resize_vector_float (float **data_vec, int N_data);
void delete_vector_float (float **data_vec);
// float matrix management
void create_matrix_float (float ***data_mat, int Nx_data, int Ny_data);
void create_0_matrix_float (float ***data_mat, int Nx_data, int Ny_data);
void resize_matrix_float (float ***data_mat, int Nx_data, int Ny_data, int Nx_old_data, int Ny_old_data);
void delete_matrix_float (float ***data_mat, int Nx_data);

// double vector management
void create_vector_double (double **data_vec, int N_data);
void create_0_vector_double (double **data_vec, int N_data);
void resize_vector_double (double **data_vec, int N_data);
void delete_vector_double (double **data_vec);
// float matrix management
void create_matrix_double (double ***data_mat, int Nx_data, int Ny_data);
void create_0_matrix_double (double ***data_mat, int Nx_data, int Ny_data);
void resize_matrix_double (double ***data_mat, int Nx_data, int Ny_data, int Nx_old_data, int Ny_old_data);
void delete_matrix_double (double ***data_mat, int Nx_data);

// complex float vector management
void create_vector_cplxf (cplxf **data_vec, int N_data);
void create_0_vector_cplxf (cplxf **data_vec, int N_data);
void resize_vector_cplxf (cplxf **data_vec, int N_data);
void delete_vector_cplxf (cplxf **data_vec);

// complex double vector management
void create_vector_cplxd (cplxd **data_vec, int N_data);
void create_0_vector_cplxd (cplxd **data_vec, int N_data);
void resize_vector_cplxd (cplxd **data_vec, int N_data);
void delete_vector_cplxd (cplxd **data_vec);

#endif
