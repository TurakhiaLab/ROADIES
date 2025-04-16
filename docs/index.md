 <div align="center">

<img src="images/ROADIES_logo.png" width="300" height="300"/>

</div>

# Reference-free Orthology-free Annotation-free DIscordance aware Estimation of Species tree (ROADIES)

## ROADIES Video Tutorial

<iframe width="1000" height="600" src="https://www.youtube.com/embed/1sR741TvZnM?si=xfktnTaQj4LUsNp0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

## Introduction

Welcome to the official wiki of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES offers a fully automated, easy-to-use, scalable solution, eliminating any manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>

## Key Features
- **Reference-free**: ROADIES ensures unbiased results by eliminating reference bias, enabling accurate species tree inference by randomly sampling genes from raw genome assemblies.
- **Orthology-free**: ROADIES automates the process of species tree inference from their raw genome assemblies without requiring any intermediate gene annotations or orthologous groups. It allows multi-copy gene trees (inferred from homologous regions) and does not require the challenging step of orthologous detection prior to gene tree inference. 
- **Annotation-free**: ROADIES does not require any input genome annotations, as ROADIES randomly samples genes within the pipeline itself.
- **Discordance-aware**: Instead of single-copy genes, ROADIES considers multi-copy genes while analyzing species trees and takes care of the possible gene discordances such as paralogs, horizontal gene transfer, incomplete lineage sorting. It uses a state-of-the-art and statistically consistent discordance-aware method to combine gene trees into a species tree.
- **Scalability**: ROADIES handles both small-scale and large-scale datasets efficiently, including diverse life forms such as mammals, flies, and birds. ROADIES also scales efficiently with multiple cores and produces faster results.
- **Flexibility**: ROADIES allows users to tune the tradeoff between accuracy and runtime by configuring the parameters and tailoring the pipeline to their specific needs.
- **Debugging options**: ROADIES provides multiple plots as output for graphical analysis, making it easier for the user to debug. 

## ROADIES Pipeline Overview
ROADIES pipeline consists of multiple stages, from raw genome assemblies to species tree estimation, with several user-configurable parameters in each stage. 

- **Stage 1: Random Sampling**: ROADIES randomly samples a configured number of subsequences from input genomic assemblies. Each of the subsequences is treated as a gene.
- **Stage 2: Pairwise alignment**: All sampled subsequences or genes are aligned with all input assemblies individually to find the homologous regions using the pairwise alignment tool [LASTZ](https://lastz.github.io/lastz/). 
- **Stage 3: Filtering of alignments**: ROADIES filters the low-quality alignments to reduce further redundant computation and limits the number of repetitive alignments from small genomic regions. 
- **Stage 4: Multiple sequence alignment**: ROADIES gathers the homologous regions for all genes across species and performs multiple sequence alignments for each of them using [PASTA](https://github.com/smirarab/pasta). 
- **Stage 5: Gene tree estimation**: ROADIES estimates gene trees from multiple sequence alignments of each of the genes using the maximum-likelihood based tree estimation tool [RAxML-NG](https://github.com/amkozlov/raxml-ng).
- **Stage 6: Species tree estimation**: After getting all the gene trees, ROADIES performs species tree estimation using a widely adopted discordance aware species tree estimation tool [ASTRAL-Pro](https://github.com/chaoszhang/A-pro). 

<div align="center">

<img src="images/drawing_github.png"width="1000" height="300" />

</div>

## Modes of operation

ROADIES supports multiple modes of operation based on various user requirements considering the tradeoff between accuracy and runtime. 

- **Accurate-Mode**: This is the default mode of operation and is preferred for accuracy-critical use cases. Here, the multiple sequence alignment stage is performed by [PASTA](https://github.com/smirarab/pasta) and the tree-building stage is governed by [RAxML-NG](https://github.com/amkozlov/raxml-ng).
- **Fast-Mode**: This mode of operation is preferred for achieving faster results, for runtime-critical use cases. Here, the multiple sequence alignment and tree-building stage is performed by [MashTree](https://github.com/lskatz/mashtree).
- **Balanced-Mode**: This mode of operation is preferred where the user wants an optimal runtime vs accuracy tradeoff. Here, the multiple sequence alignment stage is performed by [PASTA](https://github.com/smirarab/pasta), and the tree-building stage is performed using [FastTree](http://www.microbesonline.org/fasttree/). 

!!! Note
    These modes of operation can be modified using command line argument `--mode` (details mentioned in the [Usage](index.md#other-command-line-arguments) section).

!!! Note
    ROADIES also supports another mode of operation for datasets from deeper evolutionary timescales. If your datasets are from distant timescales, we recommend to switch on this mode (by adding extra argument `--deep True`; details mentioned in the [Usage](index.md#other-command-line-arguments) section).

## Convergence Mechanism

The initial count of the genes is crucial to get the accurate species tree at the end. The number of genes sufficient for getting the accurate tree also varies with datasets. Hence, ROADIES incorporates an adaptive algorithm for establishing accurate trees by tracking its confident scores. It performs multiple iterations of the entire pipeline and stops if it gets the confident tree, otherwise it continues with more gene counts. The confidence of the tree is evaluated by the confidence of its branches (or local posterior probability). The tree having most of the confident branches with high posterior probability are considered to be confident and stable. 

!!! Note
    [ASTRAL-Pro3](https://github.com/chaoszhang/A-pro) provides the information of all the internal nodes in the form of quartets (and its support values, such as local posterior probability) for every species tree per iteration. ROADIES gathers this information and keeps track of all the nodes with high support values. If the percentage change in the number of highly supported nodes gets minimal with a given number of iterations, then we say that the species tree is now converged.

!!! Note
    Users have the option to run ROADIES in non converge mode (for only one iteration) using `--noconverge` argument (details mentioned in [Usage](index.md#other-command-line-arguments) section).

<img src="images/converge_manuscript.png"width="1000" height="300" />

