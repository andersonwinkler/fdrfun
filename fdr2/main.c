#include <stdio.h>
#include <stdlib.h>
#include <nifti2_io.h>
#include "bky.h"
#include "core.h"
#include "nifti1.h"

#define flt double
#define DT_CALC DT_FLOAT64
#define DT64

#if defined(__APPLE__)
#define kOS "MacOS"
#elif (defined(__linux) || defined(__linux__))
#define kOS "Linux"
#else
#define kOS "Windows"
#endif

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#if defined(__ICC) || defined(__INTEL_COMPILER)
	#define kCCsuf  " IntelCC" STR(__INTEL_COMPILER)
#elif defined(_MSC_VER)
	#define kCCsuf  " MSC" STR(_MSC_VER)
#elif defined(__clang__)
	#define kCCsuf  " Clang" STR(__clang_major__) "." STR(__clang_minor__) "." STR(__clang_patchlevel__)
#elif defined(__GNUC__) || defined(__GNUG__)
    #define kCCsuf  " GCC" STR(__GNUC__) "." STR(__GNUC_MINOR__) "." STR(__GNUC_PATCHLEVEL__)
#else
	#define kCCsuf " CompilerNA" //unknown compiler!
#endif

#if defined(__arm__) || defined(__ARM_ARCH)
    #define kCPUsuf " ARM"
#elif defined(__x86_64)
    #define kCPUsuf " x86-64"
#else
    #define kCPUsuf " " //unknown CPU
#endif

#define kVersDate "v1.0.20230609 FDR2" 
#define kFdrVers kVersDate " " kCCsuf kCPUsuf


/*---------------------------------------------------------------
  compute log of complete beta function, using the
  Unix math library's log gamma function.  If this is
  not available, see the end of this file.
-----------------------------------------------------------------*/



/*---------------------------------------------------------------
     TRANSLATED FROM THE ORIGINAL FORTRAN:
     algorithm as 63  appl. statist. (1973), vol.22, no.3

     computes incomplete beta function ratio for arguments
     x between zero and one, p and q positive.
     log of complete beta function, beta, is assumed to be known
-----------------------------------------------------------------*/

#define ZERO 0.0
#define ONE  1.0
#define ACU  1.0e-15

double incbeta( double x , double p , double q , double beta )
{
   double betain , psq , cx , xx,pp,qq , term,ai , temp , rx ;
   int indx , ns ;

   if( p <= ZERO || q <= ZERO ) return -1.0 ;  /* error! */

   if( x <= ZERO ) return ZERO ;
   if( x >= ONE  ) return ONE ;

   /**  change tail if necessary and determine s **/

   psq = p+q ;
   cx  = ONE-x ;
   if(  p < psq*x ){
      xx   = cx ;
      cx   = x ;
      pp   = q ;
      qq   = p ;
      indx = 1 ;
   } else {
      xx   = x ;
      pp   = p ;
      qq   = q ;
      indx = 0 ;
   }

   term   = ONE ;
   ai     = ONE ;
   betain = ONE ;
   ns     = qq + cx*psq ;

   /** use soper's reduction formulae **/

      rx = xx/cx ;

lab3:
      temp = qq-ai ;
      if(ns == 0) rx = xx ;

lab4:
      term   = term*temp*rx/(pp+ai) ;
      betain = betain+term ;
      temp   = fabs(term) ;
      if(temp <= ACU && temp <= ACU*betain) goto lab5 ;

      ai = ai+ONE ;
      ns = ns-1 ;
      if(ns >= 0) goto lab3 ;
      temp = psq ;
      psq  = psq+ONE ;
      goto lab4 ;

lab5:
      betain = betain*exp(pp*log(xx)+(qq-ONE)*log(cx)-beta)/pp ;
      if(indx) betain=ONE-betain ;

   return betain ;
}

double lnbeta( double p , double q )
{
   return (lgamma(p) + lgamma(q) - lgamma(p+q)) ;
}

double student_t2p( double tt , double dof )
{
   double bb , xx , pp ;

   if( tt <= 0.0 || dof < 1.0 ) return 1.0 ;

   bb = lnbeta( 0.5*dof , 0.5 ) ;
   xx = dof/(dof + tt*tt) ;
   pp = incbeta( xx , 0.5*dof , 0.5 , bb ) ;
   return pp ;
}

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

void showHelp(int argc, char *argv[]) {
	printf("%s (%llu-bit %s)\n", kFdrVers, (unsigned long long)sizeof(size_t) * 8, kOS);
	printf("Anderson Winkler, Thomas Nichols & Chris Rorden (2023) 2-Clause BSD License\n");
	printf("Usage:\n");
	printf(" %s -i <pvalimage> [options]\n", argv[0]);
	printf("Examples:\n");
	printf(" %s -i <pvalimage> -m <maskimage> -q 0.05\n", argv[0]);
	printf(" %s -i <pvalimage> -a <adjustedimage>\n", argv[0]);
	printf("Required Arguments:\n");
	printf("  -i  input NIfTI p-value, t-statistic of z-score image.\n");
	printf("Options:\n");
	printf("  -m  mask filename (default: none)\n");
	printf("  -q  q-value (FDR) threshold (0..1, default: 0.05)\n");
	printf("  -a  output NIfTI image with FDR-adjusted p-values (default: none)\n");
	printf("  -v  switch on verbose messages\n");
	printf("  -h  display help\n");
}

int main(int argc, char *argv[]) {
	//set defaults
	double qval = 0.05;
	char *pvalFnm = NULL;
	char *adjFnm = NULL;
	char *maskFnm = NULL;
	bool isVerbose = true;
	//read arguments
	int i = 1;
	char *ptr;
	while (i < (argc)) {
		if ((strlen(argv[i]) > 1) && (argv[i][0] == '-')) { //command
			if (argv[i][1] == 'a') {
				i++;
				adjFnm = argv[i];
			} //adj filename
			if ((argv[i][1] == 'h') || (strstr(argv[i], "--help") != NULL)) {
				showHelp(argc, argv);
				return 2;
			} //help
			if ((argv[i][1] == 'i') || (strstr(argv[i], "--in") != NULL)) {
				i++;
				pvalFnm = argv[i];
			} //pval filename
			if (argv[i][1] == 'm') {
				i++;
				maskFnm = argv[i];
			} //mask filename
			if (argv[i][1] == 'q') {
				i++;
				qval = strtod(argv[i], &ptr);
			} //help
			if (argv[i][1] == 'v') {
				isVerbose = true;
			} //verbose
		} //commands start with '-'
		i++; //read next parameter
	} //for each argument
	if (pvalFnm == NULL) {
		printf("compulsory p-value filename not provided\n");
		showHelp(argc, argv);
		return 2;
	}
	nifti_image *nim = nifti_image_read2(pvalFnm, 1);
	if (nim->data == NULL)
		exit(33);
	if (nim->nvox != (nim->nx * nim->ny * nim->nz)) {
		printf("Only designed for 3D images\n");
		exit(99);
	}
	//convert image data to float64
	in_hdr ihdr = set_input_hdr(nim);
	if (nifti_image_change_datatype(nim, DT_CALC, &ihdr) != 0) {
		printf("Error changing datatype.\n");
		nifti_image_free(nim);
		exit(1);
	}
	flt *img = (flt *)nim->data;
	//first step: convert data to pvalues
	if ((isVerbose) && (nim->intent_code != NIFTI_INTENT_PVAL)) {
		flt mn = img[0];
		flt mx = img[0];
		for (size_t i = 0; i < nim->nvox; i++) {
			mn = MIN(mn, img[i]);
			mx = MAX(mx, img[i]);
		}
		printf("Minimum and maximum intensities are: %g and %g\n", mn, mx);
	}
	if (nim->intent_code == NIFTI_INTENT_PVAL) {
		//printf("input already p-values!\n");
	} else if (nim->intent_code == NIFTI_INTENT_LOGPVAL) {
		for (size_t i = 0; i < nim->nvox; i++)
			img[i] = exp(-fabs(img[i]));
	} else if (nim->intent_code == NIFTI_INTENT_LOG10PVAL) {
		for (size_t i = 0; i < nim->nvox; i++)
			img[i] = pow(10.,-fabs(img[i]));
	} else if (nim->intent_code == NIFTI_INTENT_ZSCORE) {
		#ifdef DT32 //issue8
		flt mn = -5.41;
		#else
		flt mn = -8.29;
		#endif
		flt mx = 13.0;
		size_t nClamp = 0;
		for (size_t i = 0; i < nim->nvox; i++) {
			if ((img[i] < mn) || (img[i] >= mx)) nClamp++;
			img[i] = qg(MAX(img[i], mn));
		}
		if (nClamp > 0) printf("ztop clamped %zu extreme z-scores\n", nClamp);
	} else if (nim->intent_code == NIFTI_INTENT_TTEST) {
		for (size_t i = 0; i < nim->nvox; i++)
			img[i] = student_t2p(img[i], nim->intent_p1);
	} else if (nim->intent_code == NIFTI_INTENT_FTEST) {
		printf("Please convert data to p-values (F-Test not supported, intent code %d).\n", nim->intent_code);
		exit(1);
	} else {
		printf("Unknown intention (%d): assuming voxels are p-values.\n", nim->intent_code);
	}
	size_t nvox = nim->nvox;
	int* mask = (int*)malloc(nvox*sizeof(int));
	for (size_t i = 0; i < nvox; i++)
		mask[i] = -1;
	size_t nmaskvox = 0;
	if (maskFnm != NULL) {
		if (isVerbose)
			printf("Loading mask '%s'\n", maskFnm);
		nifti_image *mnim = nifti_image_read2(maskFnm, 1);
		
		if (mnim->data == NULL)
			exit(33);
		if (nim->nvox != mnim->nvox) {
			printf("Mask must have same dimensions as statistics.");
			exit(21);
		}
		in_hdr mhdr = set_input_hdr(mnim);
		nifti_image_change_datatype(mnim, DT_CALC, &mhdr);
		flt *mimg = (flt *)mnim->data;
		for (size_t i = 0; i < nvox; i++) {
			if (mimg[i] > 0.0001) {//check inclusive vs inclusive
				mask[i] = nmaskvox;
				nmaskvox ++;
			}
		}
		nifti_image_free(mnim);
	} else {
		//https://github.com/fithisux/FSL/blob/7aa2932949129f5c61af912ea677d4dbda843895/src/randomise/fdr.cc#L155
		for (size_t i = 0; i < nvox; i++) {
			if (img[i] < 0.9999) {//check inclusive vs inclusive
				mask[i] = nmaskvox;
				nmaskvox ++;
			}
		}
	}
	if (nmaskvox < 1) {
		printf("No voxels in mask\n");
		exit(23);
	}
	if (isVerbose) {
		printf("Size = (%lld %lld %lld)\n", nim->nx, nim->ny, nim->nz);
		printf("Dims = (%g,%g,%g)\n", nim->dx, nim->dy, nim->dz);
		printf("%zu of %lld voxels are in the mask\n", nmaskvox, nim->nvox);
	}
	double* padj = NULL;
	if (adjFnm != NULL)
		padj = (double*)malloc(nmaskvox*sizeof(double));
	double* pvals = (double*)malloc(nmaskvox*sizeof(double));
	for (size_t i = 0; i < nvox; i++) {
		if (mask[i] >= 0)
			pvals[mask[i]] = img[i];
	}
	double thresh = bky(pvals, padj, nmaskvox, qval);
	if (adjFnm != NULL) {
		for (size_t i = 0; i < nvox; i++) {
			img[i] = 1.0; //masked regions have a p-value of 1.0
			if (mask[i] >= 0)
				img[i] = pvals[mask[i]];
		}
	
		if (nifti_set_filenames(nim, adjFnm, 0, 1)) {
			printf("Unable to create a new file named '%s'\n", adjFnm);
			return 12;
		}
		//set NIfTI header
		char blank_string[128];
		memset(&blank_string[0], 0, sizeof(blank_string));
		memcpy(nim->descrip, blank_string, 79);
		nim->descrip[79] = '\0';
		strcat(nim->descrip, kVersDate); //target fslmaths version
		memcpy(nim->aux_file, blank_string, 23);
		nim->aux_file[23] = '\0';
		memcpy(nim->intent_name, blank_string, 15);
		strcat(nim->intent_name, "FDR-adjusted p");
		nim->intent_code = NIFTI_INTENT_ESTIMATE;
		nim->cal_min = 0.0;
		nim->cal_max = 0.0;
		nim->scl_slope = 1.0;
		nim->scl_inter = 0.0;
		nifti_save(nim, "");
	}
	nifti_image_free(nim);
	free(padj);
	free(pvals);
	free(mask);
	printf("Probability Threshold of %zu voxels with q-value %g\n", nvox, qval);
	printf("%.20f\n", thresh);
	return 0;
}
