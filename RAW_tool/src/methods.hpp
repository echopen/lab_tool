#ifndef METHODS_HPP
#define METHODS_HPP

/*-----------------------------------------
Author: Jérôme Dubois
Contributors:
Version: 1.1
Date: 09/2018
Descritption:
-----------------------------------------*/


void loadvec(std::string name, float *buff, const unsigned int N);
void writevec(std::string name, cplxf *buff, const unsigned int N);
void writevec(std::string name, float *buff, const unsigned int N);
void writeim(std::string name, int16_t **buff, const unsigned int Nline, const unsigned int col);
void writeim(std::string name, int **buff, const unsigned int Nline, const unsigned int col);
void writeim(std::string name, float **buff, const unsigned int Nline, const unsigned int col);
void init_line(SOCKET *sock, int Npoint);
void get_image(SOCKET *sock, int16_t **data, int Nl, int Np);
void image_envelop (sig_proc_float *spf1, int16_t **datai, float **dataf, float fech, const int Nline, const int Npoint);

#endif
