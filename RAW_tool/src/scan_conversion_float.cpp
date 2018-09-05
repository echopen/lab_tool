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
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"pointer_management.hpp"

#include"scan_conversion_float.hpp"

//x_y_tensor
void create_x_y_tensor (x_y_tensor *xyt, int Nxt, int Nyt)
{
	xyt->Nx=Nxt;
	xyt->Ny=Nyt;

	create_vector_float(&xyt->x_vector,Nxt);
	create_vector_float(&xyt->y_vector,Nyt);
	create_matrix_int(&xyt->indicial_x_y_tmp,2,Nxt*Nyt);
	create_matrix_int(&xyt->indicial_x_y_out_tmp,2,Nxt*Nyt);
}

void fill_x_y_vector (scan_conv *scan_conv_struct)
{
  int i=0;
  float half_sec=scan_conv_struct->sector/2.0;
  float x0=-scan_conv_struct->polar_tensor->Rf*sin(half_sec);
	float xf= -x0;
	float delta_x=(2.0*xf)/(((float)scan_conv_struct->Nx)-1.0);
	float y0=scan_conv_struct->polar_tensor->R0*cos(half_sec);
	float yf=scan_conv_struct->polar_tensor->Rf;
	float delta_y=(yf-y0)/(((float)scan_conv_struct->Ny)-1.0);

	scan_conv_struct->cartesian_tensor->x_vector[0]=x0;
	scan_conv_struct->cartesian_tensor->y_vector[0]=y0;
	if (scan_conv_struct->Nx<=scan_conv_struct->Ny)
	{
		for (i=1 ; i<scan_conv_struct->Nx-1 ; i++)
		{
			scan_conv_struct->cartesian_tensor->x_vector[i]=scan_conv_struct->cartesian_tensor->x_vector[i-1]+delta_x;
			scan_conv_struct->cartesian_tensor->y_vector[i]=scan_conv_struct->cartesian_tensor->y_vector[i-1]+delta_y;
		}
		for (i=scan_conv_struct->Nx-1 ; i<scan_conv_struct->Ny-1 ; i++)
		{
			scan_conv_struct->cartesian_tensor->y_vector[i]=scan_conv_struct->cartesian_tensor->y_vector[i-1]+delta_y;
		}
	}
	else
	{
		for (i=1 ; i<scan_conv_struct->Ny-1 ; i++)
                {
                        scan_conv_struct->cartesian_tensor->x_vector[i]=scan_conv_struct->cartesian_tensor->x_vector[i-1]+delta_x;
                        scan_conv_struct->cartesian_tensor->y_vector[i]=scan_conv_struct->cartesian_tensor->y_vector[i-1]+delta_y;
                }
                for (i=scan_conv_struct->Ny-1 ; i<scan_conv_struct->Nx-1 ; i++)
                {
                        scan_conv_struct->cartesian_tensor->x_vector[i]=scan_conv_struct->cartesian_tensor->x_vector[i-1]+delta_x;
                }
	}
	scan_conv_struct->cartesian_tensor->x_vector[scan_conv_struct->Nx-1]=xf;
	scan_conv_struct->cartesian_tensor->y_vector[scan_conv_struct->Ny-1]=yf;
}

void resize_x_y_tensor (x_y_tensor *xyt, int Nxt, int Nyt)
{
	int Nyold=xyt->Nx*xyt->Ny;
	xyt->Nx=Nxt;
	xyt->Ny=Nyt;

	resize_vector_float(&xyt->x_vector,Nxt);
	resize_vector_float(&xyt->y_vector,Nyt);
	resize_matrix_int(&xyt->indicial_x_y_tmp,2,Nxt*Nyt,2,Nyold);
	resize_matrix_int(&xyt->indicial_x_y_out_tmp,2,Nxt*Nyt,2,Nyold);
}

void clear_x_y_tensor (x_y_tensor *xyt)
{
	delete_vector_float(&xyt->x_vector);
	delete_vector_float(&xyt->y_vector);
	delete_matrix_int(&xyt->indicial_x_y_tmp,2);
	delete_matrix_int(&xyt->indicial_x_y_out_tmp,2);
}

//r_theta_tensor
void create_r_theta_tensor (r_theta_tensor *rtt, float R0t, float Rft, int Nrt, int Nlt, int Nxt, int Nyt)
{
	rtt->R0=R0t;
	rtt->Rf=Rft;
	rtt->Nr=Nrt;
	rtt->Nl=Nlt;
	rtt->Nx=Nxt;
	rtt->Ny=Nyt;
	create_vector_float(&rtt->r_vector,rtt->Nr);
	create_vector_float(&rtt->theta_vector,rtt->Nl);
	create_matrix_float(&rtt->r_matrix,rtt->Nx,rtt->Ny);
	create_matrix_float(&rtt->theta_matrix,rtt->Nx,rtt->Ny);
}

void fill_r_theta_vectors (scan_conv *scan_conv_struct)
{
  int i=0;
  float delta_r=(scan_conv_struct->polar_tensor->Rf-scan_conv_struct->polar_tensor->R0)/(((float)scan_conv_struct->Nr)-1.0);
	float delta_theta=scan_conv_struct->sector/(((float)scan_conv_struct->Nline)-1.0);
  float half_sec=scan_conv_struct->sector/2;

	scan_conv_struct->polar_tensor->r_vector[0]=scan_conv_struct->polar_tensor->R0;
	for (i=1 ; i<scan_conv_struct->Nr-1 ; i++) {scan_conv_struct->polar_tensor->r_vector[i]=scan_conv_struct->polar_tensor->r_vector[i-1]+delta_r;}
	scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->Nr-1]=scan_conv_struct->polar_tensor->Rf;

  scan_conv_struct->polar_tensor->theta_vector[0]=-half_sec;
	for (i=1 ; i<scan_conv_struct->Nline-1 ; i++) {scan_conv_struct->polar_tensor->theta_vector[i]=scan_conv_struct->polar_tensor->theta_vector[i-1]+delta_theta;}
	scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->Nline-1]=half_sec;
}

void fill_r_theta_matrices (scan_conv *scan_conv_struct)
{
  int i=0, j=0, N_tmp=0, N_out_tmp=0;
  float half_sec=scan_conv_struct->sector/2;

	for (j=0 ; j<scan_conv_struct->Ny ; j++)
	{
		for (i=0 ; i<scan_conv_struct->Nx ; i++)
		{
			scan_conv_struct->polar_tensor->r_matrix[i][j]=sqrt(scan_conv_struct->cartesian_tensor->x_vector[i]*scan_conv_struct->cartesian_tensor->x_vector[i]+scan_conv_struct->cartesian_tensor->y_vector[j]*scan_conv_struct->cartesian_tensor->y_vector[j]);
			scan_conv_struct->polar_tensor->theta_matrix[i][j]=atan(scan_conv_struct->cartesian_tensor->x_vector[i]/scan_conv_struct->cartesian_tensor->y_vector[j]);
			if (scan_conv_struct->polar_tensor->r_matrix[i][j]<scan_conv_struct->polar_tensor->R0 || scan_conv_struct->polar_tensor->r_matrix[i][j]>scan_conv_struct->polar_tensor->Rf || scan_conv_struct->polar_tensor->theta_matrix[i][j]<-half_sec || scan_conv_struct->polar_tensor->theta_matrix[i][j]>half_sec)
			{
				scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[0][N_out_tmp]=i;
				scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[1][N_out_tmp]=j;
				N_out_tmp++;
			}
			else
			{
				scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[0][N_tmp]=i;
				scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[1][N_tmp]=j;
				N_tmp++;
			}
		}
	}
	scan_conv_struct->N_point_to_change=N_tmp;
  scan_conv_struct->N_out=scan_conv_struct->Nx*scan_conv_struct->Ny-N_tmp;
}

void resize_r_theta_tensor (r_theta_tensor *rtt, float R0t, float Rft, int Nrt, int Nlt, int Nxt, int Nyt)
{
	int Nxold=rtt->Nx, Nyold=rtt->Ny;

  rtt->R0=R0t;
  rtt->Rf=Rft;
  rtt->Nr=Nrt;
  rtt->Nl=Nlt;
  rtt->Nx=Nxt;
  rtt->Ny=Nyt;
	resize_vector_float(&rtt->r_vector,rtt->Nr);
	resize_vector_float(&rtt->theta_vector,rtt->Nl);
	resize_matrix_float(&rtt->r_matrix,rtt->Nx,rtt->Ny, Nxold, Nyold);
	resize_matrix_float(&rtt->theta_matrix,rtt->Nx,rtt->Ny, Nxold, Nyold);
}

void clear_r_theta_tensor(r_theta_tensor *rtt)
{
	delete_vector_float(&rtt->r_vector);
	delete_vector_float(&rtt->theta_vector);
	delete_matrix_float(&rtt->r_matrix,rtt->Nx);
	delete_matrix_float(&rtt->theta_matrix,rtt->Nx);
}

//image management
void create_image (scan_conv *scan_conv_struct, int Nxt, int Nyt)
{
	create_matrix_int(&scan_conv_struct->image, Nxt, Nyt);
}

void resize_image (scan_conv *scan_conv_struct, int Nxnew, int Nynew, int Nxold, int Nyold)
{
	resize_matrix_int(&scan_conv_struct->image, Nxnew, Nynew, Nxold, Nyold);
}

void delete_image (scan_conv *scan_conv_struct)
{
	delete_matrix_int(&scan_conv_struct->image,scan_conv_struct->Nx);
}


void create_indicial_x_y_buffers (scan_conv *scan_conv_struct)
{
	create_matrix_int(&scan_conv_struct->indicial_x_y, 2, scan_conv_struct->N_point_to_change);
  create_matrix_int(&scan_conv_struct->indicial_x_y_out, 2, scan_conv_struct->N_out);
}

void fill_indicial_x_y_buffers (scan_conv *scan_conv_struct)
{
  int i=0;

	if (scan_conv_struct->N_point_to_change<=scan_conv_struct->N_out)
	{
		for (i=0 ; i<scan_conv_struct->N_point_to_change ; i++)
		{
                        //printf("first loop i =%i\n",i);
			scan_conv_struct->indicial_x_y[0][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[0][i];
			scan_conv_struct->indicial_x_y[1][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[1][i];
			scan_conv_struct->indicial_x_y_out[0][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[0][i];
			scan_conv_struct->indicial_x_y_out[1][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[1][i];
		}
		for (i=scan_conv_struct->N_point_to_change ; i<scan_conv_struct->N_out ; i++)
		{
                        //printf("second loop i =%i\n",i);
			scan_conv_struct->indicial_x_y_out[0][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[0][i];
			scan_conv_struct->indicial_x_y_out[1][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[1][i];
		}
	}
	else
	{
    for (i=0 ; i<scan_conv_struct->N_out ; i++)
    {
      scan_conv_struct->indicial_x_y[0][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[0][i];
      scan_conv_struct->indicial_x_y[1][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[1][i];
      scan_conv_struct->indicial_x_y_out[0][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[0][i];
      scan_conv_struct->indicial_x_y_out[1][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_out_tmp[1][i];
    }
    for (i=scan_conv_struct->N_out ; i<scan_conv_struct->N_point_to_change ; i++)
    {
      scan_conv_struct->indicial_x_y[0][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[0][i];
      scan_conv_struct->indicial_x_y[1][i]=scan_conv_struct->cartesian_tensor->indicial_x_y_tmp[1][i];
    }
	}
}

void resize_indicial_x_y_buffers (scan_conv *scan_conv_struct, int Nin_old, int Nout_old)
{
	resize_matrix_int(&scan_conv_struct->indicial_x_y, 2, scan_conv_struct->N_point_to_change, 2, Nin_old);
  resize_matrix_int(&scan_conv_struct->indicial_x_y_out, 2, scan_conv_struct->N_out, 2, Nout_old);
}

void delete_indicial_x_y_buffers (scan_conv *scan_conv_struct)
{
	delete_matrix_int(&scan_conv_struct->indicial_x_y, 2);
  delete_matrix_int(&scan_conv_struct->indicial_x_y_out, 2);
}

void create_indicial_weight_matrix(scan_conv *scan_conv_struct)
{
  create_matrix_int(&scan_conv_struct->indicial_r_theta, 4, scan_conv_struct->N_point_to_change);
  create_matrix_float(&scan_conv_struct->weight, 4, scan_conv_struct->N_point_to_change);
}

void resize_indicial_weight_matrix(scan_conv *scan_conv_struct, int Nin_old)
{
  resize_matrix_int(&scan_conv_struct->indicial_r_theta, 4, scan_conv_struct->N_point_to_change, 4, Nin_old);
  resize_matrix_float(&scan_conv_struct->weight, 4, scan_conv_struct->N_point_to_change, 4, Nin_old);
}

void delete_indicial_weight_matrix(scan_conv *scan_conv_struct)
{
  delete_matrix_int(&scan_conv_struct->indicial_r_theta, 4);
  delete_matrix_float(&scan_conv_struct->weight, 4);
}

void indicial_matrix_calculation(scan_conv *scan_conv_struct)
{
	int i=0, k=0, l=0;
	float r_old=0.0, r_new=0.0, theta_old=0.0, theta_new=0.0;

	for (i=0 ; i<scan_conv_struct->N_point_to_change ; i++)
	{
		r_new=scan_conv_struct->polar_tensor->r_matrix[scan_conv_struct->indicial_x_y[0][i]][scan_conv_struct->indicial_x_y[1][i]];
		theta_new=scan_conv_struct->polar_tensor->theta_matrix[scan_conv_struct->indicial_x_y[0][i]][scan_conv_struct->indicial_x_y[1][i]];
		if (r_new>r_old) //r increase and decrease continuously
		{
			while (scan_conv_struct->polar_tensor->r_vector[k+1]<r_new && k<scan_conv_struct->Nr-1) {k++;}
		}
		else if (r_new<r_old)
		{
			while (scan_conv_struct->polar_tensor->r_vector[k]>r_new && k>0) {k--;}
		}
		r_old=r_new;
		if (theta_new<theta_old) //theta increase continuously but decrease from theta to -theta
		{
			l=0;
			while (scan_conv_struct->polar_tensor->theta_vector[l+1]<theta_new && l<scan_conv_struct->Nline-1) {l++;}
		}
		else if (theta_new>theta_old)
		{
			while (scan_conv_struct->polar_tensor->theta_vector[l+1]<theta_new && l<scan_conv_struct->Nline-1) {l++;}
		}
		theta_old=theta_new;

		scan_conv_struct->indicial_r_theta[0][i]=k;
		scan_conv_struct->indicial_r_theta[1][i]=k+1;
		scan_conv_struct->indicial_r_theta[2][i]=l;
		scan_conv_struct->indicial_r_theta[3][i]=l+1;
	}
}

void weight_matrix_calculation(scan_conv *scan_conv_struct)
{
	int i=0;
	float ri=0.0, ti=0.0, r1=0.0, r2=0.0, t1=0.0, t2=0.0, P1=0.0, P2=0.0, P3=0.0, P4=0.0;

	if (scan_conv_struct->option==0) //linear interpolation
	{
		for (i=0 ; i<scan_conv_struct->N_point_to_change ; i++)
		{
			ri=scan_conv_struct->polar_tensor->r_matrix[scan_conv_struct->indicial_x_y[0][i]][scan_conv_struct->indicial_x_y[1][i]];
			ti=scan_conv_struct->polar_tensor->theta_matrix[scan_conv_struct->indicial_x_y[0][i]][scan_conv_struct->indicial_x_y[1][i]];
			r1=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[0][i]];
			r2=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[1][i]];
			t1=scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[2][i]];
			t2=scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[3][i]];
			P1=(r2-ri)/(r2-r1);
			P2=(ri-r1)/(r2-r1);
			P3=(t2-ti)/(t2-t1);
			P4=(ti-t1)/(t2-t1);
			scan_conv_struct->weight[0][i]=P1*P3;
			scan_conv_struct->weight[1][i]=P1*P4;
			scan_conv_struct->weight[2][i]=P2*P3;
			scan_conv_struct->weight[3][i]=P2*P4;
		}
	}
	else if (scan_conv_struct->option==1) //distance interpolation
	{
		float r3=0.0, r4=0.0;
		float xi=0.0, xt=0.0, yi=0.0, yt=0.0;
		float norm=0.0;
		for (i=0 ; i<scan_conv_struct->N_point_to_change ; i++)
		{
			xi=scan_conv_struct->cartesian_tensor->x_vector[scan_conv_struct->indicial_x_y[0][i]];
			yi=scan_conv_struct->cartesian_tensor->y_vector[scan_conv_struct->indicial_x_y[1][i]];
			xt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[0][i]]*sin(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[2][i]]);
			yt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[0][i]]*cos(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[2][i]]);
			r1=sqrt((xi-xt)*(xi-xt)+(yi-yt)*(yi-yt));
			xt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[0][i]]*sin(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[3][i]]);
      yt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[0][i]]*cos(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[3][i]]);
      r2=sqrt((xi-xt)*(xi-xt)+(yi-yt)*(yi-yt));
			xt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[1][i]]*sin(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[2][i]]);
      yt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[1][i]]*cos(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[2][i]]);
      r3=sqrt((xi-xt)*(xi-xt)+(yi-yt)*(yi-yt));
			xt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[1][i]]*sin(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[3][i]]);
      yt=scan_conv_struct->polar_tensor->r_vector[scan_conv_struct->indicial_r_theta[1][i]]*cos(scan_conv_struct->polar_tensor->theta_vector[scan_conv_struct->indicial_r_theta[3][i]]);
      r4=sqrt((xi-xt)*(xi-xt)+(yi-yt)*(yi-yt));
			norm=1/r1+1/r2+1/r3+1/r4;
			norm=1/norm;
			scan_conv_struct->weight[0][i]=norm/r1;
      scan_conv_struct->weight[1][i]=norm/r2;
      scan_conv_struct->weight[2][i]=norm/r3;
      scan_conv_struct->weight[3][i]=norm/r4;
		}
	}
}

void create_scan_conv_struct (scan_conv *scan_conv_struct, int Nr_probe, int Nline_probe, float sector_probe, float R0_probe, float Rf_probe, int Nx_im, int Ny_im, int option_selection)
{
  //init attribuite of the structure
	sector_probe=sector_probe*PIf/180.0;
	scan_conv_struct->cartesian_tensor=(x_y_tensor *)malloc(sizeof(x_y_tensor));
	scan_conv_struct->polar_tensor=(r_theta_tensor *)malloc(sizeof(r_theta_tensor));
	scan_conv_struct->Nr=Nr_probe;
	scan_conv_struct->Nline=Nline_probe;
	scan_conv_struct->sector=sector_probe;
  scan_conv_struct->polar_tensor->R0=R0_probe;
  scan_conv_struct->polar_tensor->Rf=Rf_probe;
	scan_conv_struct->Nx=Nx_im;
	scan_conv_struct->Ny=Ny_im;
	scan_conv_struct->option=option_selection;
  create_x_y_tensor (scan_conv_struct->cartesian_tensor, Nx_im, Ny_im);
  fill_x_y_vector(scan_conv_struct);
  create_r_theta_tensor (scan_conv_struct->polar_tensor, R0_probe, Rf_probe, Nr_probe, Nline_probe, Nx_im, Ny_im);
  fill_r_theta_vectors(scan_conv_struct);
  fill_r_theta_matrices(scan_conv_struct);
  create_image(scan_conv_struct, Nx_im, Ny_im);
  create_indicial_x_y_buffers(scan_conv_struct);
  fill_indicial_x_y_buffers(scan_conv_struct);
  create_indicial_weight_matrix(scan_conv_struct);
  indicial_matrix_calculation(scan_conv_struct);
  weight_matrix_calculation(scan_conv_struct);

  //clear all temporary buffer
  clear_x_y_tensor(scan_conv_struct->cartesian_tensor);
  clear_r_theta_tensor(scan_conv_struct->polar_tensor);
	printf("scan conv struct initiated\n");
}

void resize_scan_conv_struct (scan_conv *scan_conv_struct, int Nr_probe, int Nline_probe, float sector_probe, float R0_probe, float Rf_probe, int Nx_im, int Ny_im, int option_selection)
{
  //take old value for resizing attribuite of the structure
	sector_probe=sector_probe*PIf/180.0;
  int Nx_old, Ny_old, Nin_old, Nout_old;//, option_old, Nr_old, Nline_old;
  //float sector_old, R0_old, Rf_old;
	Nx_old=scan_conv_struct->Nx;
	Ny_old=scan_conv_struct->Ny;
  Nin_old=scan_conv_struct->N_point_to_change;
  Nout_old=scan_conv_struct->N_out;
  /*Nr_old=scan_conv_struct->Nr;
	Nline_old=scan_conv_struct->Nline;
	sector_old=scan_conv_struct->sector;
  R0_old=scan_conv_struct->polar_tensor->R0;
  Rf_old=scan_conv_struct->polar_tensor->Rf;
	option_old=scan_conv_struct->option;*/

	scan_conv_struct->Nr=Nr_probe;
	scan_conv_struct->Nline=Nline_probe;
	scan_conv_struct->sector=sector_probe;
  scan_conv_struct->polar_tensor->R0=R0_probe;
  scan_conv_struct->polar_tensor->Rf=Rf_probe;
	scan_conv_struct->Nx=Nx_im;
	scan_conv_struct->Ny=Ny_im;
	scan_conv_struct->option=option_selection;

  create_x_y_tensor (scan_conv_struct->cartesian_tensor, Nx_im, Ny_im);
  fill_x_y_vector(scan_conv_struct);
  create_r_theta_tensor (scan_conv_struct->polar_tensor, R0_probe, Rf_probe, Nr_probe, Nline_probe, Nx_im, Ny_im);
  fill_r_theta_vectors(scan_conv_struct);
  fill_r_theta_matrices(scan_conv_struct);
	if (Nx_im!=Nx_old || Ny_im!=Ny_old)  {resize_image(scan_conv_struct, Nx_im, Ny_im, Nx_old, Ny_old);}
  resize_indicial_x_y_buffers(scan_conv_struct, Nin_old, Nout_old);
  fill_indicial_x_y_buffers(scan_conv_struct);
  resize_indicial_weight_matrix(scan_conv_struct, Nin_old);
  indicial_matrix_calculation(scan_conv_struct);
  weight_matrix_calculation(scan_conv_struct);

  //clear all temporary buffer
  clear_x_y_tensor(scan_conv_struct->cartesian_tensor);
  clear_r_theta_tensor(scan_conv_struct->polar_tensor);
}

void delete_scan_conv_struct (scan_conv *scan_conv_struct)
{
  delete_image(scan_conv_struct);
  delete_indicial_x_y_buffers(scan_conv_struct);
  delete_indicial_weight_matrix(scan_conv_struct);
	free(scan_conv_struct->cartesian_tensor);
	free(scan_conv_struct->polar_tensor);
}

void image_scan_conversion (scan_conv *scan_conv_struct, float **image_r_theta)
{
  int i=0;
	float P1=0.0, P2=0.0, P3=0.0, P4=0.0;
  for (i=0 ; i<scan_conv_struct->N_point_to_change ; i++)
  {
		P1=scan_conv_struct->weight[0][i]*image_r_theta[scan_conv_struct->indicial_r_theta[2][i]][scan_conv_struct->indicial_r_theta[0][i]];
		P2=scan_conv_struct->weight[1][i]*image_r_theta[scan_conv_struct->indicial_r_theta[3][i]][scan_conv_struct->indicial_r_theta[0][i]];
		P3=scan_conv_struct->weight[2][i]*image_r_theta[scan_conv_struct->indicial_r_theta[2][i]][scan_conv_struct->indicial_r_theta[1][i]];
		P4=scan_conv_struct->weight[3][i]*image_r_theta[scan_conv_struct->indicial_r_theta[3][i]][scan_conv_struct->indicial_r_theta[1][i]];
    scan_conv_struct->image[scan_conv_struct->indicial_x_y[0][i]][scan_conv_struct->indicial_x_y[1][i]]=(int)(P1+P2+P3+P4);
  }
}

void change_image_background (scan_conv *scan_conv_struct, float background)
{
	int i=0;
	for (i=0 ; i<scan_conv_struct->N_out ; i++)
	{
		scan_conv_struct->image[scan_conv_struct->indicial_x_y_out[0][i]][scan_conv_struct->indicial_x_y_out[1][i]]=background;
	}
}
