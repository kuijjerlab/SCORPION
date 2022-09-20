# <img src="https://raw.githubusercontent.com/kuijjerlab/SCORPION/main/inst/logoSCORPION.png" width="30" title="SCORPION"> SCORPION

**SCORPION** (<ins>S</ins>ingle-<ins>C</ins>ell <ins>O</ins>riented <ins>R</ins>econstruction of <ins>P</ins>ANDA <ins>I</ins>ndividually <ins>O</ins>ptimized Gene Regulatory <ins>N</ins>etworks), is an R package that uses coarse-graining of single-cell/nuclei RNA-seq data to reduce sparsity and improve the ability to detect the gene regulatory network's underlying correlation structure. The coarse-grained data generated is then used to reconstruct the gene regulatory network using a network refinement strategy through the **PANDA** (<ins>P</ins>assing <ins>A</ins>ttributes between <ins>N</ins>etworks for <ins>D</ins>ata <ins>A</ins>ssimilation) message passing algorithm. This algorithm is designed to integrate multiple sources of information such as protein-protein interaction, gene expression, and sequence motif data to predict accurate regulatory relationships. Thanks to the use of the same baseline priors in each instance, this approach can reconstruct comparable, fully-connected, weighted, and directed transcriptome-wide single-cell gene regulatory networks suitable for use in population-level studies.
## Method
![method](https://raw.githubusercontent.com/kuijjerlab/SCORPION/main/inst/methodSCORPION.png)

## Usage
**SCORPION** is under active development, you can install it, using the following command:
```{r}
library(remotes)
install_github('kuijjerlab/SCORPION')
library(SCORPION)
```
