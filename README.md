# CUWCA
Codon usage analysis using within-group CA (correspondence analysis)

This repository includes scripts that implement codon usage analysis that were carried out in previous studies ([Qiu et al 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/jpy.12294), [Qiu et al 2011](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-325)). The codon frequency data is analyzed using within-group correspondence analysis (WCA) available from 'dudi.coa' function implemented in the ADE-4 R package ([Thioulouse et al 1997](https://link.springer.com/article/10.1023/A:1018513530268)). WCA takes into consideration both amino acid composition and codon degeneracy and is less biased than relative synonymous codon usage (CA-RSCU) ([Suzuki et al 2008](https://academic.oup.com/dnaresearch/article/15/6/357/513030)).

## Dependencies
ADE-4 (https://cran.r-project.org/web/packages/ade4/index.html)

Perl (https://www.python.org/)

R (https://github.com/RomelTorres/alpha_vantage)

##  Quick start

1. Downdown this repository and move all files in 'scripts' directory to your path.

2. Change directory to 'exampleData' directory and execute the following command to analysis example input data:
```
codonUsageWithinCA.sh exampleCds.txt 300
```
It results in self-explanatory filtered sequence file, files containing analysis statstics, and data plot in PDF files.
 
## Usage
codonUsageWithinCA.sh takes two arugments

1. A file containing coding sequences (excluding non-coding regions such as 5'-UTR) in fasta format

2. A cutoff value to defined the minmal length for a sequences to be included in analysis. Short sequences will be discarded.