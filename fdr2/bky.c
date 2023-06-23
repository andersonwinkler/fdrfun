#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifndef _MSC_VER
 #include <unistd.h>
#endif
#include <stdint.h>
#include "bky.h"

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

struct fltidx
{
    flote value;
    int idx;
};

int cmp(const void *a, const void *b)
{
    struct fltidx *a1 = (struct fltidx *)a;
    struct fltidx *a2 = (struct fltidx *)b;
    if ((*a1).value > (*a2).value)
        return 1;
    else if ((*a1).value < (*a2).value)
        return -1;
    else
        return 0;
}

//The arrays pvals and padj should both have length nvox
// pvals provides the unsorted p-values
// padj will be filled with the adjusted p-values
//If padj is not required, pass a NULL pointer:
// bky(pvals, NULL, nvox, 0.05)

flote bky(flote *pvals, flote *padj, size_t nvox, double qval) {
	double cV = 1.0;
	if ((qval < 0.0) || (qval > 1.0)) {
		printf("qval out of range (0-1).\n");
		return NAN;
	}
	//========[PART 1: FDR THRESHOLD]========================================
	//Sort p-values
	struct fltidx* pvalsLo2Hi = malloc(nvox * sizeof(struct fltidx));
	 for (int i = 0; i < nvox; i++) {
		pvalsLo2Hi[i].value = pvals[i];
		pvalsLo2Hi[i].idx = i;
	}
	qsort (pvalsLo2Hi, nvox, sizeof(struct fltidx), cmp);
	if ((pvalsLo2Hi[0].value < 0.0) || (pvalsLo2Hi[nvox-1].value > 1.0)) {
		printf("Values should be in range 0..1 not %g..%g.\n", pvalsLo2Hi[0].value, pvalsLo2Hi[nvox-1].value);
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
			flote thrline = (v2+1)*qval/(nvox+1-(v+1)*(1-qval))/cV;
			//detect if p-vals that survive the cutoff
			if (pvalsLo2Hi[v2].value <= thrline) {
				pSurvives = true;
				idxthr[v2] = true;
			}
			v2++;
		} while ( v2 < nvox );
		if (!pSurvives)
			break;
		v++;
	} while ( v < nvox );
	flote thr = -1.0;
	for (int v = 0; v < nvox; v++)
		if (idxthr[v])
			thr = pvalsLo2Hi[v].value;
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
	if (padj != NULL) {
		//========[PART 2: FDR CORRECTED]========================================
		flote* pcor = (flote*)malloc(nvox*sizeof(flote));
		for (int v = 0; v < nvox; v++)
			pcor[v] = 1.0;
		for (int v = 0; v < nvox; v++) {
			//% For every p-value, this is the minimum q that will eventually
			// satisfy Definition #7 of Benjamini, Krieger and Yekutieli (2006).
			flote mn = INFINITY;
			for (int v2 = v; v2 < nvox; v2++) {
				flote q = pvalsLo2Hi[v2].value*(nvox+1-(v+1))*cV/((v2+1)-(v+1)*pvalsLo2Hi[v2].value);
				mn = MIN(mn, q);
			}
			if (mn > 1.0) break;
			pcor[v] = mn;
		}
		//========[PART 3: FDR ADJUSTED ]========================================
		// The p-adjusted for the current p-value is the cummulative maximum
		// up to it. Note that this is different from equation #3 of
		// Yekutieli & Benjamini (1999) that is used for p-value adjustment of
		// the usual FDR (BH).
		//padjs: sorted p-adjusted
		flote* padjs = (flote*)malloc(nvox*sizeof(flote));
		padjs[0] = pcor[0];
		for (int v = 1; v < nvox; v++)
			padjs[v] = MAX(padjs[v-1], pcor[v]);
		for (int v = 0; v < nvox; v++)
			padj[pvalsLo2Hi[v].idx] = padjs[v];
		free(padjs);
		free(pcor);
	}
	free(pvalsLo2Hi);
	return thr;
}
