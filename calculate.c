// calculate the sharpness of profiles
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
#include "ptime_sharpness.h"

#define ITMAX 10000  // Maximum allowed number of iterations.
#define EPS 3.0e-16 // Machine double floating-point precision.
//#define EPS 3.0e-8 // Machine floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}

	return temp[n-1];
}

int remove_baseline (double *in, double frac_off, int n, double *out)
{
	// define the off_pulse range of std, frac_off is the fraction of the phase
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse(n, in, off_0, frac_off);

	// normalize std
	int i;
	double baseline = 0.0;
  for (i = 0; i < num_off; i++)
  {
      baseline += off_0[i];
	//printf ("%lf \n", off_0[i]);
      //average_s += s_off[i];
  }
	baseline = baseline/num_off;

  printf ("the baseline of std is: %lf \n", baseline);
  //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int off_pulse (int nphase, double *in, double *out, double frac)
// define the on_pulse range, 10% of the phase
{
	int n = nphase;
	int num = (int)(n*frac);
	int i,j;
	double small;
	double ave;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num; j++)
			{
				small += in[j];
			}
			small = small/num;
		}
			
		ave = 0.0;
		for(j = 0; j < num; j++)
		{
			if ((i+j) > n-1)
			{
				ave += in[(i+j)-(n-1)];
			}
			else 
			{
				ave += in[i+j];
			}
		}
		ave = ave/num;

		small = (ave <= small ? ave : small);
		index = (ave <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	for (i = 0; i < num; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

int pre (double *s, int nphase, double *real_s, double *ima_s)
// dft std, output real_s, ima_s  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i = 0; i < nphase; i++)
	{
		test[i] = s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

  fftw_complex *out_s;
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	

	{
	  dft_profiles(nphase,s,out_s);
	  //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	  //double amp_s[N/2],phi_s[N/2];
	  //double amp_p[N/2],phi_p[N/2];
		/*
		for (i = 0; i < nphase; i++)                                
		{                                                        
			printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		}
		*/

		for (i = 0; i < nphase/2+1; i++)                                
		{                                                        
			real_s[i]=out_s[i][0];                                       
			ima_s[i]=out_s[i][1];                                        
		}
										
	}

	fftw_free(out_s); 
	fftw_free(out_t); 

	return 0;
}

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new)
// inverse dft, and add pi/2 phase
{
	double *dp;
  fftw_plan plan;
	fftw_complex *cp;

  dp = (double *)malloc(sizeof (double) * ncount);
	cp = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp, 0, sizeof (double) * ncount);
	memset(cp, 0, sizeof (fftw_complex) * ncount);

	// initialize the dft...
	double *dp_t;
  fftw_plan plan_t;
	fftw_complex *cp_t;

  dp_t = (double *)malloc(sizeof (double) * ncount);
	cp_t = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp_t, 0, sizeof (double) * ncount);
	memset(cp_t, 0, sizeof (fftw_complex) * ncount);

	int i;
  double real,ima,amp,cosina,sina;

	for (i = 0; i < ncount; i++)
	{
		if (i < ncount/2+1)
		//if (i < ncount/2)
		{
      real = real_p[i];
      ima = ima_p[i];
			amp = sqrt(real*real+ima*ima);

			if ( amp != 0.0 )
			{
				cosina = real/amp;
				sina = ima/amp;

				cp[i][0] = amp*(cosina*cos(-i*3.1415926)-sina*sin(-i*3.1415926));
				cp[i][1] = amp*(sina*cos(-i*3.1415926)+cosina*sin(-i*3.1415926));

				cp_t[i][0] = real_p[i];
				cp_t[i][1] = ima_p[i];
			}
			else
			{
				cp[i][0]=0.0;
				cp[i][1]=0.0;
				cp_t[i][0]=0.0;
				cp_t[i][1]=0.0;
			}
		}
		else
		{
			cp[i][0]=0.0;
			cp[i][1]=0.0;
			cp_t[i][0]=0.0;
			cp_t[i][1]=0.0;
		}
	}

  plan_t = fftw_plan_dft_c2r_1d(ncount, cp_t, dp_t, FFTW_MEASURE);

  fftw_execute(plan_t);

  fftw_destroy_plan(plan_t);

	/////////////////////////////////////////////////////////////////

  plan = fftw_plan_dft_c2r_1d(ncount, cp, dp, FFTW_MEASURE);

  fftw_execute(plan);

  fftw_destroy_plan(plan);

	for (i = 0; i < ncount; i++)
	{
		p_new[i] = dp[i]/ncount;  // normalized by the ncount
		//printf ("%lf\n", p_new[i]);
	}

	return 0;
}

int sharpness (double *y, int nphase, double psrfreq, double frac)
{
	// move the std to the center
	double real_s[nphase/2+1],ima_s[nphase/2+1];
	pre(y, nphase, real_s, ima_s);

	/*
	int j;
	for (j = 0; j < nphase/2+1; j++)
	{
		printf ("%d %lf\n", j, real_s[j]);
	}
	*/

	double std_0[nphase];
	inverse_dft (real_s, ima_s, nphase, std_0);
	
	// remove the baseline
	double std_1[nphase];
	remove_baseline (std_0, frac, nphase, std_1);

	// normalize the std
	double peak;
	peak = find_peak_value(nphase, std_1);
	printf ("peak is: %lf\n", peak);

  double u_0[nphase];

  double t_step=(1.0/psrfreq)/nphase;
  //double t_step=(1.0/173.687946184761)/1024;
  
  int i;
  for (i = 0; i < nphase; i++)
  {
		u_0[i]=std_1[i]/peak;
  }

	/*
	int j;
	for (j = 0; j < nphase; j++)
	{
		printf ("%d %lf\n", j, std_0[j]);
	}
	*/

  double u_1=0.0;
  
  for (i = 0; i < nphase-1; i++)
  //for (i = 0; i < nphase; i++)
  {
		u_1 += (u_0[i+1]-u_0[i])*(u_0[i+1]-u_0[i]);
	}

  double weff;

  weff = 1000000*t_step/u_1;

	printf ("W_eff: %.10lf \n", weff);

	return 0;
}

