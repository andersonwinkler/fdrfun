#include <stdio.h>
#include <stdlib.h>
#include "bky.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		printf("Usage:");
		printf(" %s filename.bin [qvalue] [verbose]\n", argv[0]);
		printf("  - filename is binary float32 data with native endian\n");
		printf("  - qvalue is fdr threshold (0..1, default 0.05)\n");
		printf("  - verbose reporting (y/n default y)\n");
		printf("Example with 5%% false discovery rate:\n");
		printf("  %s filename.bin 0.05\n", argv[0]);
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
	if( (fp = fopen(argv[1], "r") ) == NULL ) {
		printf("Cannot open file: %s\n", argv[1]);
		return -1;
	}
	fseek(fp, 0L, SEEK_END);
	size_t sz = ftell(fp);
	rewind(fp);
	int nvox = sz >> 2; //4 bytes per float
	float* pvals = (float*)malloc(nvox*sizeof(float));
	size_t n = fread(pvals, sizeof(float), nvox, fp);
	fclose(fp);
	float thresh = bky(pvals, nvox, qval);
	free(pvals);
	if (verbose)
		printf("%d pvalues with FDR qvalue %g threshold %g\n", nvox, qval, thresh);
	else
		printf("%g\n", thresh);
	return 0;
}
