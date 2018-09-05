#include<stdio.h>
#include<stdlib.h>
#include<string>


#include"src/pointer_management.hpp"
#include"src/TCP_API.hpp"
#include"src/sig_proc_float.hpp"
#include"src/scan_conversion_float.hpp"
#include"src/methods.hpp"

typedef int SOCKET;
typedef float _Complex cplxf;

#define Port 7538

float r0=0.0, rf=0.0;
int Nline=0, decimation=0, Npoint=0;
double sector=0.0;
int mode_RP;


int main(void)
{
  //socket variable
  SOCKET sock;
  //const char *IP="192.168.43.223";
  const char *IP="0.0.0.0";

  init_TCP_client(&sock, IP, Port);
  get_RP_settings(&sock);
  float fech=0.0;
  fech=125.0/((float)decimation);

  int Npoint=0, loop=0;
  Npoint=(int)(2.0*(rf-r0)*125.0/1.48/((double)decimation));
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

  image_envelop(spf1, image0, imagee, fech, Nline, Npoint);

  unsigned int nx=256, ny=256;
  scan_conv scan_conv_struct={0};
  create_scan_conv_struct (&scan_conv_struct, Npoint, Nline, (float)sector, r0, rf, nx, ny, 0);
  image_scan_conversion (&scan_conv_struct, imagee);
  std::string toto2="./image/image_sc.txt";
  writeim(toto2, scan_conv_struct.image, nx, ny);

  delete_matrix_int16_t(&image0, Nline);
  delete_matrix_float(&imagee, Nline);
  delete_scan_conv_struct(&scan_conv_struct);
  return 0;
}
