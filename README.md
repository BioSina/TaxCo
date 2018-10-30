# TaxCo
#### TaxCo - Correlation Analysis for Taxonomic Data
TaxCo computes correlations between taxonomic abundances and numeric metadata and presents the results in tabular and graphic form.

## Dependencies:
[NetworkX](https://pypi.org/project/networkx/2.2/)

[matplotlib](https://matplotlib.org)

## Usage:

```
usage: TaxCo.py [-h]
                (-qiime <L2> <L3> <L4> <L5> <L6> | -mothur <tax.summary file> | -megan <L2> <L3> <L4> <L5> <L6>)
                [-out <output directory>] [-pval <p-value>]
                [-cutoffs <decimal number> [<decimal number> ...]]
                [-coefficients <correlation coefficients> [<correlation coefficients> ...]]
                [-metadata <metadata>]

optional arguments:
  -h, --help            show this help message and exit
  -out <output directory>
                        Path to output directory
  -pval <p-value>       P-value cutoff
  -cutoffs <decimal number> [<decimal number> ...]
                        Cutoffs for correlation strength, with a dot as
                        decimal separator
  -coefficients <correlation coefficients> [<correlation coefficients> ...]
                        One or more correlation coefficients (kendall,
                        spearman, pearson)
  -metadata <metadata>  Metadata in the format MEGAN exports it in (tab
                        separated, samples in rows, data in columns)

Input files (mutually exclusive):
  -qiime <L2> <L3> <L4> <L5> <L6>
                        Paths to QIIME taxon summaries, sorted from L2 to L6.
                        Samples have to be in the same order in each file.
  -mothur <tax.summary file>
                        Path to mothur tax.summary file (including all
                        samples)
  -megan <L2> <L3> <L4> <L5> <L6>
                        Paths to MEGAN CSV files, sorted from L2 to L6.
                        Collapse tree to matching level and export the
                        summarized counts.
```
