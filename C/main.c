#include <stdio.h>
#include <stdlib.h>
#include "bky.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		printf("Usage:");
		printf(" %s filename.bin [qvalue] [verbose] [outfile]\n", argv[0]);
		printf("  - filename is binary float64 data with native endian\n");
		printf("  - qvalue is fdr threshold (0..1, default 0.05)\n");
		printf("  - verbose reporting (y/n default y)\n");
		printf("  - output filename for adjusted values (string, default none)\n");
		printf("Examples:\n");
		printf("  %s filename.bin 0.05\n", argv[0]);
		printf("  %s filename.bin 0.05 y padj.bin\n", argv[0]);
		return 2;
	}
	char *ptr;
	double qval = 0.05;
	if (argc > 2)	
		qval = strtod(argv[2], &ptr);
	if (argc > 2)
		qval = strtod(argv[2], &ptr);
	bool verbose = true;
	if ((argc > 3) && ((argv[3][0] == 'n') || (argv[3][0] == 'N')))
		verbose = false;
	FILE *fp;
	if( (fp = fopen(argv[1], "rb") ) == NULL ) {
		printf("Cannot open file: %s\n", argv[1]);
		return -1;
	}
	fseek(fp, 0L, SEEK_END);
	size_t sz = ftell(fp);
	rewind(fp);
	int nvox = sz >> 3; //8 bytes per float64
	double* pvals = (double*)malloc(nvox*sizeof(double));
	size_t n = fread(pvals, sizeof(double), nvox, fp);
	fclose(fp);
	double* padj = (double*)malloc(nvox*sizeof(double));
	double thresh = bky(pvals, padj, nvox, qval);
	if (argc > 4) {
		if( (fp = fopen(argv[4], "wb") ) == NULL ) {
			printf("Cannot open file: %s\n", argv[4]);
			return -1;
		}
		fwrite(padj, sizeof(double), nvox, fp);
		fclose(fp);
	}
	free(padj);
	free(pvals);
	if (verbose)
		printf("%d pvalues with FDR qvalue %g threshold %g\n", nvox, qval, thresh);
	else
		printf("%g\n", thresh);
	return 0;
}
