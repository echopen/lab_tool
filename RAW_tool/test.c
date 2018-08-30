#include<stdio.h>
#include<stdlib.h>

#include"src/TCP_API.h"
#include"src/pointer_management.h"
#include"src/sig_proc_float.h"
#include"src/scan_conversion_float.h"

#define Port 7538

float r0=0.0, rf=0.0;
int Nline=0, dec=0, Npoint=0;
double sector=0.0;
int mode_RP;

void loadb(char *name, float *buff, int N)
{
  FILE *h=NULL;
  float tmp=0.0;
  int i=0;

  h=fopen(name,"r");
  if (h==NULL)
  {
    printf("file not found\n");
  }
  else
  {
    printf("file opened\n");
    for (i=0 ; i<N ; i++)
    {
      fscanf(h, "%f", &tmp);
      buff[i]=tmp;
    }
    printf("close file\n");
    fclose(h);
  }
}

void writeb(char *name, cplxf *buff, int N)
{
  FILE *h=NULL;
  float tmp=0.0;
  int i=0;

  h=fopen(name,"w+");
  if (h==NULL)
  {
    printf("file not found\n");
  }
  else
  {
    for (i=0 ; i<N ; i++)
    {
      fprintf(h, "%f %f\n", creal(buff[i]), cimag(buff[i]));
    }
    fclose(h);
  }
}

void writef(char *name, float *buff, int N)
{
  FILE *h=NULL;
  float tmp=0.0;
  int i=0;

  h=fopen(name,"w+");
  if (h==NULL)
  {
    printf("file not found\n");
  }
  else
  {
    for (i=0 ; i<N ; i++)
    {
      fprintf(h, "%f\n", buff[i]);
    }
    fclose(h);
  }
}

void writeim(char *name, int16_t **buff, int Nline, int col)
{
  FILE *h=NULL;
  int i=0, j=0;

  h=fopen(name,"w+");
  if (h==NULL)
  {
    printf("file not found\n");
  }
  else
  {
    for (i=0 ; i<Nline ; i++)
    {
      for (j=0 ; j<col ; j++)
      {
        fprintf(h, "%i ", buff[i][j]);
      }
      fprintf(h,"\n");
    }
    fclose(h);
  }
}

void writeimfloat(char *name, float **buff, int Nline, int col)
{
  FILE *h=NULL;
  int i=0, j=0;

  h=fopen(name,"w+");
  if (h==NULL)
  {
    printf("file not found\n");
  }
  else
  {
    for (i=0 ; i<Nline ; i++)
    {
      for (j=0 ; j<col ; j++)
      {
        fprintf(h, "%f ", buff[i][j]);
      }
      fprintf(h,"\n");
    }
    fclose(h);
  }
}

void init_line(SOCKET *sock, int Npoint)
{
  int16_t *data=NULL;
  int old=0, new=0;
  create_vector_int16_t(&data, Npoint+1);
  while(1)
  {
    receive_int16_TCP_client(sock, data, Npoint+1);
    new=data[0];
    //printf("line = %i\n",new);
    if (new==1 & old==2) {break;}
    old=new;
  }
  delete_vector_int16_t(&data);
}

void get_image(SOCKET *sock, int16_t **data, int Nl, int Np)
{
  int16_t line=0;
  int k=0;
  for (k=0 ; k<Nl ; k++)
  {
    receive_int16_TCP_client(sock, &line, 1);
    //printf("nline = %i\n",line);
    receive_int16_TCP_client(sock, data[line], Np);
  }
}

int main(int argc, char const *argv[])
{
  //socket variable
  SOCKET sock;
  //const char *IP="192.168.43.223";
  const char *IP="0.0.0.0";

  init_TCP_client(&sock, IP, Port);
  get_RP_settings(&sock);
  printf("r0=%f\n",r0);
  printf("rf=%f\n",rf);
  printf("dec=%i\n",dec);
  printf("Nline=%i\n",Nline);
  printf("sector=%f\n",sector);
  printf("mode_RP=%i\n",mode_RP);
  float fech=0.0;
  fech=125.0/((float)dec);

  int Npoint=0, loop=0;
  Npoint=(int)(2.0*(rf-r0)*125.0/1.48/((double)dec));
	if (Npoint>16384){Npoint=16384;}
	printf("Npoint = %i\n",Npoint);

  //init struc fft
  sig_proc_float * spf1;
  init_struct_spf(&spf1,Npoint);


  int16_t **image0=NULL;
  float **imagee=NULL;
  create_matrix_int16_t(&image0, Nline, Npoint);
  create_matrix_float(&imagee, Nline, Npoint);

  //wait for line 1 to start image at beggining
  init_line(&sock, Npoint);
  //get image
  get_image(&sock, image0, Nline, Npoint);
  //char toto[20]="image_RAW.txt";
  //writeim(toto, image0, Nline, Npoint);

  int k=0, line=0;
  for (k=0 ; k<Npoint ; k++)
  {
    spf1->sig[k]=(float)image0[3][k];
  }
  //average of the signal is determined for "zero" padding
  mean_fsignal_calculus(spf1);
  fftf_preparation(spf1);
  paddingf(spf1, spf1->mean_sig);
  init_hilberttransformf(spf1, fech, 0.5, 8.0, 0);//determination of the pass band filer to apply when using hilbert transform

  for(line=0 ; line<Nline ; line++)
  {
    for (k=0 ; k<Npoint ; k++)
    {
      spf1->sig[k]=(float)image0[line][k];
    }
    fftf_preparation(spf1); //fft preparation must be done each time due to fft algorithm that modify input "vector"
    paddingf(spf1, spf1->mean_sig);
    envelope_detectionf(spf1);
    for (k=0 ; k<Npoint ; k++)
    {
      imagee[line][k]=spf1->envelope[k];
    }
  }

  //char toto1[50]="image_env.txt";
  //writeimfloat(toto1, imagee, Nline, Npoint);

  int nx=256, ny=256;
  scan_conv scan_conv_struct={0};
  create_scan_conv_struct (&scan_conv_struct, Npoint, Nline, (float)sector, r0, rf, nx, ny, 0);
  image_scan_conversion (&scan_conv_struct, imagee);
  char toto2[50]="image_sc.txt";
  writeimfloat(toto2, scan_conv_struct.image, nx, ny);



  delete_matrix_int16_t(&image0, Nline, Npoint);
  delete_matrix_float(&imagee, Nline, Npoint);
  delete_scan_conv_struct(&scan_conv_struct);
  return 0;
}
