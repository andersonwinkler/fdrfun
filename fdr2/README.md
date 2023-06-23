### Introdcution

FDR2 is a simple program which takes in a p-value, z-score or t-statistic image and uses FDR ( False Discovery Rate ) theory to carry out multiple comparison correction. It extends the [FSL fdr](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDR) program being truly open source (BSD license), uses the more recent BKY (2006) method and has features for signal in both tails. 

You will first need to download the repository and compile the C Code:

```
git clone --branch AddCCode git@github.com:andersonwinkler/fdrfun.git
cd fdrfun/fdr2
make
./fdr2  -i ../narps-4735_50GV-hypo1_unthresh -m ../narps-mask -v -q 0.05 -a padj
```


### Links

 - Benjamini & Hochberg. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. R. Statist. Soc. B ([1995](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x)) 57(1):289-300.
 - Yekutieli & Benjamini. Resampling-based false discovery rate controlling multiple test procedures for multiple testin procedures. J. Stat. Plan. Inf. ([1999](https://www.sciencedirect.com/science/article/abs/pii/S0378375899000415?via%3Dihub)) 82:171-96.
 - Benjamini, Krieger, and Yekutieli. Adaptive linear step-up procedures that control the false discovery rate. Biometrika. ([2006](https://academic.oup.com/biomet/article-abstract/93/3/491/380683?redirectedFrom=fulltext)) 93(3): 491–507.

