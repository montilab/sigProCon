
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigProCon

<!-- badges: start -->

<!-- badges: end -->

It’s common to have a signature of genes, then to give each sample in a
bulk rna-seq dataset a score reflecting the activity of that signature,
such as by GSVA.

The goal of sigProCon is to identify which genes from a signature are
contributing most to the sample signature activity scores across all
samples.

For example, one may have a signature of genes that are upregulated
during the process of epithelial to mesencymal transition (EMT), which
is associated with tumor metastasis. One can use a tool (such as GSVA)
to assign a score to each patient in the TCGA dataset based on each
patient’s gene expression profile, that corresponds to how ‘active’ the
EMT signature is in each patient. It would suggest that a patient with a
higher score has more activity of the genes in the EMT geneset, and thus
may have a more aggressive tumor. Some genes will be very correlated
with the score, suggesting their larger role in EMT in this particular
set of samples, while some genes in the EMT signature may have very
little correlation.

This package simply takes the correlation of each gene in the signature
with the signature score, to identify the genes that are contributing
most to the score.

## Installation

You can install the development version of sigProCon from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("montilab/sigProCon@dev")
```

## Example with toy sample data:

``` r
data(toy_signature)
data(toy_eset)
print(length(toy_signature))
#> [1] 13
print(dim(toy_eset))
#> Features  Samples 
#>      100       20
```

### View toy data tables:

``` r
toy_output <- sigProCon::signature_projection_contributors(
  eset = toy_eset,
  signature = list(signature1 = toy_signature)
)

## Show some of the genes' correlation to the signature score
print(head(toy_output$score_cor))
#>          score_cor     pval_cor     insig leading_edge
#> gene_009 0.9117601 2.228287e-08 signature          yes
#> gene_010 0.8923537 1.239783e-07 signature          yes
#> gene_003 0.7753862 5.913725e-05 signature          yes
#> gene_008 0.7645522 8.655076e-05 signature          yes
#> gene_007 0.7428196 1.756087e-04 signature          yes
#> gene_001 0.7110160 4.408161e-04 signature          yes

## Show some of the signature scores
print(head(toy_output$sig_score))
#>           sig_score
#> sample_06 0.5496272
#> sample_11 0.5329154
#> sample_03 0.5251103
#> sample_16 0.4113453
#> sample_19 0.3909520
#> sample_07 0.2712722
```

### View toy heatmap that plots all signature and non-signature genes:

``` r
print(toy_output$heatmap_all_genes)
```

<img src="man/figures/README-toy-full.heatmap-1.png" width="100%" />

### View toy heatmap that plots only signature genes:

``` r
print(toy_output$heatmap_sig_genes)
```

<img src="man/figures/README-toy-sig.heatmap-1.png" width="100%" />

### View toy k.s. plot of signature genes and scores:

``` r
print(toy_output$ks$plot)
```

<img src="man/figures/README-toy-ksplot-1.png" width="100%" />

#### Passing a pre-computed score

``` r
toy_score <- sigProCon::omics_signature_score(
  eset = toy_eset,
  signature = list(signature1 = toy_signature),
  method = "GSVA"
)
## call w/ pre-computed score 
toy_output_precomp <- sigProCon::signature_projection_contributors(
  eset = toy_eset,
  signature = list(signature1 = toy_signature),
  sig_score = toy_score
)
## same results
isTRUE(all.equal(toy_output$score_cor, toy_output_precomp$score_cor))
#> [1] TRUE
isTRUE(all.equal(toy_output$sig_score, toy_output_precomp$sig_score))
#> [1] TRUE
```

## Example with package sample data (`eset` and `signature`):

``` r
#with our sample data
data(signature)
data(eset)

output <- sigProCon::signature_projection_contributors(
  eset = eset, 
  signature = signature)
```

### View tables:

``` r
## Show some of the genes' correlation to the signature score
print(head(output$score_cor))
#>          score_cor     pval_cor      insig leading_edge
#> Gene329  0.6039086 0.0004099159 background           no
#> Gene7262 0.5914609 0.0005769769 background           no
#> Gene637  0.5795738 0.0007896397 background           no
#> Gene7822 0.5793168 0.0007949097 background           no
#> Gene2907 0.5775437 0.0008321167 background           no
#> Gene5511 0.5719439 0.0009598320 background           no
## Show some of the signature scores
print(head(output$sig_score))
#>           sig_score
#> Sample5  0.20783161
#> Sample9  0.20495746
#> Sample23 0.17228859
#> Sample12 0.09991516
#> Sample11 0.09654342
#> Sample6  0.09226228
```

### View heatmap that plots all signature and non-signature genes:

``` r
print(output$heatmap_all_genes)
```

<img src="man/figures/README-full.heatmap-1.png" width="100%" />

### View heatmap that plots only signature genes:

``` r
print(output$heatmap_sig_genes)
```

<img src="man/figures/README-sig.heatmap-1.png" width="100%" />

### View k.s. plot of signature genes and scores:

``` r
print(output$ks$plot)
```

<img src="man/figures/README-ksplot-1.png" width="100%" />
