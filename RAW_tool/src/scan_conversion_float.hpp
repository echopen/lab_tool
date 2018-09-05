#ifndef SCAN_CONVERSION_H
#define SCAN_CONVERSION_H

/*----
Author: Jérôme Dubois
Contributors:
Version: 1.1
Date: 08/2018
Descritption:
A C librairy create to make scan conversion using
linear or weight interpolation. Weight interpolation often use
in acoustic imaging, but linear interpolation give better results.

This librairy use structures to facilitate the use in C code. These structures can
be resize. Need the pointer_management.h librairy.

Usage:

scan_conv scan_conv_struct;
float **dataf=NULL;
create_matrix_float(&dataf, Nline, Npoint);
create_scan_conv_struct (&scan_conv_struct, Npoint, Nline, sector, r0, rf, nx, ny, 0); //initiate the structure and detemined the wieght matrix to make scan conversion
image_scan_conversion (&scan_conv_struct, dataf);

----*/

#define PIf 3.14159265358979323846

typedef struct x_y_tensor x_y_tensor;
typedef struct r_theta_tensor r_theta_tensor;
typedef struct scan_conv scan_conv;

void create_x_y_tensor (x_y_tensor *xyt, int Nxt, int Nyt);
void fill_x_y_vector (scan_conv *scan_conv_struct);
void resize_x_y_tensor (x_y_tensor *xyt, int Nxt, int Nyt);
void clear_x_y_tensor (x_y_tensor *xyt);

void create_r_theta_tensor (r_theta_tensor *rtt, float R0t, float Rft, int Nrt, int Nlt, int Nxt, int Nyt);
void fill_r_theta_vectors (scan_conv *scan_conv_struct);
void fill_r_theta_matrices (scan_conv *scan_conv_struct); //also fill matrix of x_y_tensor
void resize_r_theta_tensor (r_theta_tensor *rtt, float R0t, float Rft, int Nrt, int Nlt, int Nxt, int Nyt);
void clear_r_theta_tensor(r_theta_tensor *rtt);

void create_image (scan_conv *scan_conv_struct, int Nxt, int Nyt);
void resize_image (scan_conv *scan_conv_struct, int Nxnew, int Nynew, int Nxold, int Nyold);
void delete_image (scan_conv *scan_conv_struct);

void create_indicial_x_y_buffers (scan_conv *scan_conv_struct);
void fill_indicial_x_y_buffers (scan_conv *scan_conv_struct);
void resize_indicial_x_y_buffers (scan_conv *scan_conv_struct, int Nin_old, int Nout_old);
void delete_indicial_x_y_buffers (scan_conv *scan_conv_struct);

void create_indicial_weight_matrix(scan_conv *scan_conv_struct);
void resize_indicial_weight_matrix(scan_conv *scan_conv_struct, int Nin_old);
void delete_indicial_weight_matrix(scan_conv *scan_conv_struct);
void indicial_matrix_calculation(scan_conv *scan_conv_struct);
void weight_matrix_calculation(scan_conv *scan_conv_struct);

void create_scan_conv_struct (scan_conv *scan_conv_struct, int Nr_probe, int Nline_probe, float sector_probe, float R0_probe, float Rf_probe, int Nx_im, int Ny_im, int option_selection);
void resize_scan_conv_struct (scan_conv *scan_conv_struct, int Nr_probe, int Nline_probe, float sector_probe, float R0_probe, float Rf_probe, int Nx_im, int Ny_im, int option_selection);
void delete_scan_conv_struct (scan_conv *scan_conv_struct);

void image_scan_conversion (scan_conv *scan_conv_struct, float **image_r_theta);
void change_image_background (scan_conv *scan_conv_struct, float background);

//structure definition
struct x_y_tensor
{
	int Nx; //number of point of the x vector
	int Ny; //number of point of the y vector

	float *x_vector;
	float *y_vector;
	int **indicial_x_y_tmp;
	int **indicial_x_y_out_tmp;
};

struct r_theta_tensor
{
	float R0; //begin of measuring line
	float Rf; //end of measuring line
	int Nr;
	int Nl;
	int Nx;
	int Ny;

	float *r_vector;
	float *theta_vector;
	float **r_matrix;
	float **theta_matrix;
};

struct scan_conv
{
	int N_point_to_change; //number of pixel to calcul
	int N_out; //number of pixel outside the sector
	int **indicial_x_y; //2*N_point_to_change matrix, indicial in the image of the N_point_to_changepixel that have to be calculated in the image
	int **indicial_x_y_out; //2*N_point_to_change matrix, gather the indices of the point of the image out the imaging sector
	int **indicial_r_theta; //size 4*N_point_to_change, indicial of the 4 points in the r theta plan that are used in the calculus of the Nth pixel, pixel listed in indicial_x_y
	float **weight; //size 4*N_point_to_change, weight of the 4 points listed in indicial_r_theta for the calcul of the pixel listed in indicial_x_y
	int Nr; //number of point per line from the probe
	int Nline; //number of line per image from the probe
	float sector; //sector of the image
	int Nx; //number of pixel along x in the image
	int Ny; //number of pixel along y in the image
	int **image; //image
	int option; //type of scan conversion, 0 linear interpolation, 1 weight/distance interpolation

	x_y_tensor *cartesian_tensor;
	r_theta_tensor *polar_tensor;
};

#endif
