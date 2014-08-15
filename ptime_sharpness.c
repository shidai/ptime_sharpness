//"Usage: ptime_sharpness -s tname -on frac\n"
//"Calculate the sharpness of profile\n"
//tname: templates; -s: standard template format (ptime template not supported);
//frac: on-pulse fraction; 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptime_sharpness.h"

int main (int argc, char *argv[])
{
	//int h,i,j,k;

	//////////////////////////////////////////////////////
	char tname[128];   // name of template
	double frac; // on-pulse fraction

	int i;
	for (i=0;i<argc;i++)
  {
		if (strcmp(argv[i],"-s")==0)
		{
			strcpy(tname,argv[++i]);
			printf ("standard template format\n");
			//sscanf(argv[++i],"%d",&nbin);
		}
		else if (strcmp(argv[i],"-on")==0)
		{
			frac = atof(argv[++i]);
		}
		//else 
		//{
		//	printf ("Wrong option!\n");
		//	exit(1);
		//}
  }
	//printf ("%d\n", smode);

	/////////////////////////////////////////////////////////
	// read std
	char data[] = "[SUBINT]";
	char psrparam[] = "[PSRPARAM]";

	char std[50];
	char std_para[50];

	strcpy(std,tname);
	strcpy(std_para,tname);

	strcat(std, data);
	strcat(std_para, psrparam);

	double psrfreq;
	psrfreq = read_psrfreq(std_para);
	printf ("PSR frequency: %.15lf\n", psrfreq);

	int nphase, nsub, nchn;
	nchn = get_nchan(std);	
	nsub = get_subint(std);	
	nphase = get_nphase(std);	

	double s[nphase];
	int mode = 0; // only std supported
	read_std(std,nsub,s,nphase,mode,nchn);
	/////////////////////////////////////////////////////////

	// calculate the sharpness of profiles
	sharpness (s, nphase, psrfreq, frac);

	return 0;
}
