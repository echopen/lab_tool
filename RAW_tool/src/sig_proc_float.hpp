#ifndef SIG_PROC_H
#define SIG_PROC_H

/*-----------------------------------------
Author: Jérôme Dubois
Contributors:
Version: 1.1
Date: 07/2018
Descritption:
C scrit for signal processing based on the fft algorithme from rosettacode
source from https://rosettacode.org/wiki/Fast_Fourier_transform
all other algorithmes are made by Jérôme Dubois
algorithmes are made for both float and double precision
-----------------------------------------*/

//Pi
const double PI = 3.141592653589793238460;

//Structures
typedef struct fftf_table fftf_table;
typedef struct sig_proc_float sig_proc_float;

//Functions prototype
int int_pow(int x, int power);//power>= 0
int power_two(int x, int *pow);//determine the power of 2 superior and nearest to x return 0 if length(x) is a power of 2, 1 else
void fill_ftable(fftf_table *tbf);
void init_struct_fftf_table (fftf_table **tbf, int N); //N number of point for initialisation
void resize_struct_fftf_table(fftf_table *tbf, int N);
void delete_struct_fftf_table(fftf_table *tbf);
void init_struct_spf (sig_proc_float **spf, int N);
void resize_struct_spf (sig_proc_float *spf, int N);
void delete_struct_spf (sig_proc_float *spf);

void fftf_preparation(sig_proc_float *spf); //fill sigc with sig (if user want to use sig buffer)
void mean_fsignal_calculus(sig_proc_float *spf);
void paddingf (sig_proc_float *spf, float value);// value_padding, pad is filled with value from Npoint Npad
void rec_fftf(cplxf *buf, cplxf *out, int step, fftf_table *tbf);
void fftf(sig_proc_float *spf);//need a signal with length equal to a power of 2
void rec_ifftf(cplxf *buf, cplxf *out, int step, fftf_table *tbf);
void ifftf(sig_proc_float *spf);//need a signal with length equal to a power of 2
void Nifftf(sig_proc_float *spf);//normalized ifft: Nifft(fft(s))=s
void ifftf_filt(sig_proc_float *spf);//same than ifftf but use fft_filt not fft
void Nifftf_filt(sig_proc_float *spf);
void filterf_parameter(sig_proc_float *spf, float fech , float f0 , float fm , int methode);
void reset_filtf_buffer(sig_proc_float *spf);
void apply_filterf(sig_proc_float *spf);
void init_hilberttransformf(sig_proc_float *spf, float fech , float f0 , float fm , int methode);
void hilberttransformf(sig_proc_float *spf);
void Nhilberttransformf(sig_proc_float *spf); //normalized Hilbert transform
void envelope_detectionf(sig_proc_float *spf); //proportional to envelope
void Nenvelope_detectionf(sig_proc_float *spf); //true envelope

//structure definition
struct fftf_table
{
	int Npoint; //number of point of the signal
	int Npad; //number of point of the signal with zeros padding

	cplxf *fftf_table;
  cplxf *ifftf_table;
};

struct sig_proc_float
{
  float mean_sig; //used when padding with mean of the signal

  fftf_table *tbf;
  float *sig;
  float *envelope;
  cplxf *sigc;
  cplxf *fft;
  cplxf *ifft;

  float *window_hilb; //window for filtering size not Npad but size of the window
  int n0; //position of the first point of the window
  int nf; //position of the last point
  int N_win; //number of point of window_hilb
  cplxf *fft_filt;
};

#endif
