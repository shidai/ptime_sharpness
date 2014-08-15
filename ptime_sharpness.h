#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf.h>
#include <fftw3.h>

struct my_params {double a; double b;};

//////////////////////////////////////////////////////////////////////
double find_peak_value (int n, double *s);
int remove_baseline (double *in, double frac_off, int n, double *out);
int off_pulse (int nphase, double *in, double *out, double frac);
int dft_profiles (int N, double *in, fftw_complex *out);
int pre (double *s, int nphase, double *real_s, double *ima_s);
int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new);

/////////////////////////////////////////////////////////////////////
int getNint (char *filename, int *n);
int readResi (char *filename, double *resi, double *error);
double read_psrfreq ( char *name );
int get_nchan ( char *name );
int get_npol ( char *name );
int get_nphase ( char *name );
int get_subint ( char *name );
int read_scl ( char *name, int subint, double *scl, int nchan, int npol);
int read_offs ( char *name, int subint, double *offs, int nchan, int npol);
int read_value ( char *name, int subint, double *value, int nphase, int nchan, int npol);
int read_prof ( char *name, int subint, double *profile, int nphase, int npol, int nchan);
int check_std ( char *name, int subint, int mode, int nchn, int nphase);
int read_std_pt ( char *name, double *profile, int nphase, int nchn);
int read_std ( char *name, int subint, double *profile, int nphase, int mode, int nchn);

///////////////////////////////////////////////////////////////////////
int sharpness (double *y, int nphase, double psrfreq, double frac);

