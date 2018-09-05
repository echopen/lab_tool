/*-----------------------------------------
Author: Jérôme Dubois
Contributors:
Version: 1.0
Date: 09/2018
Descritption:
-----------------------------------------*/


#include<iostream>
#include<fstream>
#include<string>

#include"pointer_management.hpp"
#include"TCP_API.hpp"
#include"sig_proc_float.hpp"
//typedef int SOCKET;
//typedef float _Complex cplxf;
#include"methods.hpp"

void loadvec(std::string name, float *buff, const unsigned int N)
{
  std::ifstream txtfile;
  txtfile.open(name);
  if (!txtfile)
  {
    std::cout << "error opening file to load" << std::endl;
    exit(1);
  }

  for (unsigned int i=0 ; i<N ; i++) {txtfile >> buff[i];}
  txtfile.close();
}


void writevec(std::string name, cplxf *buff, const unsigned int N)
{
  std::ofstream txtfile;
  txtfile.open(name);
  if (!txtfile)
  {
    std::cout << "error opening file to write" << std::endl;
    exit(1);
  }

  for (unsigned int i=0 ; i<N ; i++) {txtfile << buff[i] << std::endl;}
  txtfile.close();
}

void writevec(std::string name, float *buff, const unsigned int N)
{
  std::ofstream txtfile;
  txtfile.open(name);
  if (!txtfile)
  {
    std::cout << "error opening file to write" << std::endl;
    exit(1);
  }

  for (unsigned int i=0 ; i<N ; i++) {txtfile << buff[i] << std::endl;}
  txtfile.close();
}


void writeim(std::string name, int16_t **buff, const unsigned int Nline, const unsigned int col)
{
  std::ofstream txtfile;
  txtfile.open(name);
  if (!txtfile)
  {
    std::cout << "error opening file to write" << std::endl;
    exit(1);
  }

  for (unsigned int i=0 ; i<Nline ; i++)
  {
    for (unsigned int j=0 ; j<col ; j++)
    {
      txtfile << buff[i][j] << " ";
    }
    txtfile << std::endl;
  }
  txtfile.close();
}

void writeim(std::string name, int **buff, const unsigned int Nline, const unsigned int col)
{
  std::ofstream txtfile;
  txtfile.open(name);
  if (!txtfile)
  {
    std::cout << "error opening file to write" << std::endl;
    exit(1);
  }

  for (unsigned int i=0 ; i<Nline ; i++)
  {
    for (unsigned int j=0 ; j<col ; j++)
    {
      txtfile << buff[i][j] << " ";
    }
    txtfile << std::endl;
  }
  txtfile.close();
}

void writeim(std::string name, float **buff, const unsigned int Nline, const unsigned int col)
{
  std::ofstream txtfile;
  txtfile.open(name);
  if (!txtfile)
  {
    std::cout << "error opening file to write" << std::endl;
    exit(1);
  }

  for (unsigned int i=0 ; i<Nline ; i++)
  {
    for (unsigned int j=0 ; j<col ; j++)
    {
      txtfile << buff[i][j] << " ";
    }
    txtfile << std::endl;
  }
  txtfile.close();
}

void init_line(SOCKET *sock, int Npoint)
{
  int16_t *data=NULL;
  int old_line=0, new_line=0;
  create_vector_int16_t(&data, Npoint+1);
  while(1)
  {
    receive_int16_TCP_client(sock, data, Npoint+1);
    new_line=data[0];
    //printf("line = %i\n",new);
    if (new_line==1 & old_line==2) {break;}
    old_line=new_line;
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

void image_envelop (sig_proc_float *spf1, int16_t **datai, float **dataf, float fech, const int Nline, const int Npoint)
{
    int k=0, line=0;
    for (k=0 ; k<Npoint ; k++)
    {
      spf1->sig[k]=(float)datai[3][k];
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
        spf1->sig[k]=(float)datai[line][k];
      }
      fftf_preparation(spf1); //fft preparation must be done each time due to fft algorithm that modify input "vector"
      paddingf(spf1, spf1->mean_sig);
      envelope_detectionf(spf1);
      for (k=0 ; k<Npoint ; k++)
      {
        dataf[line][k]=spf1->envelope[k];
      }
    }
}
