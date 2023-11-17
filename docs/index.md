 <div align="center">

<img src="images/ROADIES_logo.png" width="300" height="300"/>

</div>

# Reference-free Orthology-free Alignment-free DIscordance aware Estimation of Species tree (ROADIES)

## Introduction

Welcome to the official wiki of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES pipeline offers a fully automated, easy-to-use, scalable solution, eliminating any error-prone manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>

## Key Features
- **Orthology-free**: ROADIES automates the process of species tree inference from their raw genome assemblies without requiring any intermediate gene annotations or orthologous groups, making it effortless for users to generate accurate species trees.
- **Reference-free**: ROADIES ensures unbiased results by eliminating reference bias, enabling accurate species tree inference by randomly sampling genes from raw genome assemblies.
- **Discordance-aware**: Instead of single-copy genes, ROADIES considers multi-copy genes while analyzing species tree and takes care of the possible gene discordances such as paralogs, horizonal gene transfer, incomplete lineage sorting. 
- **Scalability**: ROADIES handles both small-scale and large-scale datasets efficiently, including diverse life forms such as mammals, flies, and birds. ROADIES also scales efficiently with multiple cores and produces faster results.
- **Flexibility**: ROADIES allows users to tune the tradeoff between accuracy and runtime by configuring the parameters and tailoring the pipeline to their specific needs.
- **Debugging options**: ROADIES provides multiple plots as output for graphical analysis, making it easier for the user to debug. 

## ROADIES Pipeline Overview
ROADIES pipeline consists of multiple stages, from raw genome assemblies to species tree estimation, with several user-configurable parameters in each stage. 

- **Stage 1: Random aampling**: ROADIES randomly samples a configured number of subsequences from input genomic assemblies. Each of the subsequences is treated as a gene.
- **Stage 2: Pairwise alignment**: All sampled subsequences are aligned with all input assemblies individually using [LASTZ](https://lastz.github.io/lastz/). 
- **Stage 3: Filtering of alignments**: ROADIES filters the low-quality alignments to reduce further redundant computation, gathers all homology data per gene, and limits the number of homologs per gene. 
- **Stage 4: Multiple sequence alignment**: After filtering, ROADIES gathers all genes from different species and performs multiple sequence alignments for every gene using [PASTA](https://github.com/smirarab/pasta). 
- **Stage 5: Gene tree estimation**: ROADIES estimates gene trees from multiple sequence alignments of each of the genes using [IQTREE](http://www.iqtree.org/).
- **Stage 6: Species tree estimation**: After gene tree estimation, ROADIES concatenates all gene trees in a single list and performs species tree estimation using [ASTRAL-Pro](https://github.com/chaoszhang/A-pro). 

<div align="center">

<img src="images/drawing_github.png"width="1000" height="300" />

</div>

## Modes of operation

ROADIES supports multiple modes of operation based on various user requirements considering the tradeoff between accuracy and runtime. 

- **Accurate-Mode**: This is the default mode of operation and is preferred for accuracy-critical use cases. Here, the multiple sequence alignment stage is performed by [PASTA](https://github.com/smirarab/pasta) and the tree building stage is governed by [IQTREE](http://www.iqtree.org/).
- **Fast-Mode**: This mode of operation is preferred for achieving faster results, for runtime-critical use cases. Here, the multiple sequence alignment and tree building stage is performed by [MashTree](https://github.com/lskatz/mashtree).
- **Balanced-Mode**: This mode of operation is preferred where the user wants an optimal runtime vs accuracy tradeoff. Here, the multiple sequence alignment stage is performed by [PASTA](https://github.com/smirarab/pasta), and the tree building stage is performed using [FastTree](http://www.microbesonline.org/fasttree/). 

!!! Note
    These modes of operation can be modified using command line arguments, mentioned in the [usage](usage.md#various-command-line-arguments) tab.

## Convergence Mechanism

ROADIES incorporates a method for establishing accurate and stable species trees. It performs multiple iterations of the entire pipeline and collects information on the percentage of highly supported nodes in the species tree for every iteration. 

[ASTRAL-Pro](https://github.com/chaoszhang/A-pro) provides the information of all the internal nodes in the form of quartets (and its support values such as local posterior probability) for every species tree per iteration. ROADIES gathers this information and keeps track of all the nodes with high support values. If the percentage change in the number of highly supported nodes gets minimal with a given number of iterations, then we say that the species tree is now converged.

!!! Note
    Users have the option to run ROADIES with both converge and no-converge options using command line arguments mentioned in [usage](usage.md#various-command-line-arguments) section.

<img src="images/converge_manuscript.png"width="1000" height="300" />
