# Differential test with PACS on Seurat object

Run differential test using the similar style as in Seurat
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
function. Note 1: This function is modified from Seurat package v5.1.0
Note 2: Seurat must be installed to use this function.

## Usage

``` r
FindMarkersPACS(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  latent.vars = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = "count",
  reduction = NULL,
  ...
)
```

## Arguments

- object:

  A `Seurat` object

- ident.1:

  Identity class to define markers for; pass an object of class `phylo`
  or 'clustertree' to find markers for a node in a cluster tree; passing
  'clustertree' requires `BuildClusterTree` to have been run

- ident.2:

  A second identity class for comparison; if `NULL`, use all other cells
  for comparison; if an object of class `phylo` or 'clustertree' is
  passed to `ident.1`, must pass a node to find markers for

- latent.vars:

  Latent variables to be controlled for from the meta.data column of the
  Seurat object

- group.by:

  Regroup cells into a different identity class prior to performing
  differential expression (see example)

- subset.ident:

  Subset a particular identity class prior to regrouping. Only relevant
  if group.by is set (see example)

- assay:

  Assay to use in differential expression testing, for PACS, we use
  `count` as assay

- reduction:

  Reduction to use in differential expression testing - will test for DE
  on cell embeddings

- ...:

  Additional arguments for identifying differential peaks
