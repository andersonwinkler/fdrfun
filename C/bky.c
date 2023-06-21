#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "radixsort.h"
#ifndef _MSC_VER
 #include <unistd.h>
#endif
#include <stdint.h>

/*
Note that the corrected and adjusted p-values do **not** depend
on the supplied q-value, but they do depend on the choice of c(V).

References:
* Benjamini & Hochberg. Controlling the false discovery
  rate: a practical and powerful approach to multiple testing.
  J. R. Statist. Soc. B (1995) 57(1):289-300.
* Yekutieli & Benjamini. Resampling-based false discovery rate
  controlling multiple test procedures for multiple testing
  procedures. J. Stat. Plan. Inf. (1999) 82:171-96.
* Benjamini, Krieger, and Yekutieli. Adaptive linear step-up
  procedures that control the false discovery rate.
  Biometrika. (2006) 93(3): 491–507.
*/

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

void sortFloats(float *input, float *sorted, size_t n) {
	uint32_t* idx_in = (uint32_t*)malloc(n*sizeof(uint32_t));
	uint32_t* idx_out = (uint32_t*)malloc(n*sizeof(uint32_t));
	radix11sort_f32(input, sorted, idx_in, idx_out, n);
	free(idx_in);
	free(idx_out);
}

float bky(float *pvals, size_t nvox, double qval) {
	double cV = 1.0;
	if ((qval < 0.0) || (qval > 1.0)) {
		printf("qval out of range (0-1).\n");
		return NAN;
	}
	//========[PART 1: FDR THRESHOLD]========================================
	//Sort p-values
	float* pvalsLo2Hi = (float*)malloc(nvox*sizeof(float));
	sortFloats(pvals, pvalsLo2Hi, nvox);
	if ((pvalsLo2Hi[0] < 0.0) || (pvalsLo2Hi[nvox-1] > 1.0)) {
		printf("Values out of range (0-1).\n");
		return NAN;
	}
	// indices that survive FDR, in the same size as the pvalues
	bool* idxthr = (bool*)malloc(nvox*sizeof(bool));
	for (int v = 0; v < nvox; v++)
		idxthr[v] = false;
	int v = 0;
	do {
		int v2 = v;
		bool pSurvives = false;
		do {
			// Line to be used as cutoff
			float thrline = (v2+1)*qval/(nvox+1-(v+1)*(1-qval))/cV;
			//detect if p-vals that survive the cutoff
			if (pvalsLo2Hi[v2] <= thrline) {
				pSurvives = true;
				idxthr[v2] = true;
			}
			v2++;
		} while ( v2 < nvox );
		if (!pSurvives)
			break;
		v++;
	} while ( v < nvox );
	float thr = -1.0;
	for (int v = 0; v < nvox; v++)
		if (idxthr[v])
			thr = pvalsLo2Hi[v];
	//Case when it does not cross: no pvalues survive
	if (thr < 0.0)
		thr = 0.0;
	else if (thr == 0.0) {
		//Deal with the case when all the points under the line
		//are equal to zero, and other points are above the line
		for (int v = 0; v < nvox; v++)
			if (idxthr[v])
				thr = (v+1)*qval/(nvox+1-(v+1)*(1-qval))/cV;
	}
	free(idxthr);
	free(pvalsLo2Hi);
	return thr;
}
