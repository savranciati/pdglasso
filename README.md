# Graphical Lasso for Coloured Gaussian Graphical Models for Paired Data

<img src="pdglasso_logo.png" height=150>

## How to install the package:

```{r message = FALSE, warning = FALSE}
library(devtools)
install_github("savranciati/pdglasso", build_vignettes=FALSE)
library(pdglasso)
```

## Relevant paper for the package

* Ranciati, S., Roverato, A., (2024). "On the application of Gaussian graphical models to paired data problems", **Statistics and Computing**, doi: https://doi.org/10.1007/s11222-024-10513-6 
* Ranciati, S., Roverato, A., (2024). "On the application of Gaussian graphical models to paired data problems", **arXiv pre-print**, https://arxiv.org/abs/2307.14160 [please check latest version]

### Code reproducibility and data availability
#### Updated version
The .zip folder `additional simulations' contains ten subfolders with further results concerning 10 different independent replicated runs from the same simulation scheme. These are related to an ongoing review process for a submitted version of the paper so they might slightly differ from simulations of the arXiv version.

#### arXiv version
Inside this repository we made publicly available a .zip file containing:
* *raw_data*, a folder containing the original data and some metadata information;
* *breast_cancer.RData*, the dataset used in Section 6 of the manuscript
* *breast_cancer_analysis.r*, the R script to reproduce the analysis of Section 6 of the manuscript.

Also, for replicability, we deposit a .zip version of the entire package and environment 'freezed' at the submission of the arXiv version of the paper.

## Package overview

An RCON model for paired data (pdRCON model) is a coloured Gaussian Graphical Model (GGM) where the $p$ variables are partitioned into a Left block $L$ and a right block $R$. Every variable in the left block has an homologous variable in the right block and certain types of equality *R*estrictions on the entries of the *CON*centration matrix $K$ are allowed. Every pdRCON model is uniquely represented by a Coloured Graph for Paired Data (pdColG) implemented in the form of a $p\times p$ symmetric matrix, where every entry is one of the values 0, 1 or 2, as follows: 

* The diagonal entries of the pdColG matrix are all equal to either 1 or 2. More specificaly, every variable in $L$ has an homologous variable in $R$ and the corresponding diagonal entries of $K$ can be constrained to have equal value, and in this case the corresponding entries of the pdColG matrix are set to 2. These are referred to as *coloured vertices* whereas the unconstrained diagonal entries of $K$ are encoded by the value 1 and are referred to as *uncoloured vertices*. 

* The off-diagonal entries of the pdColG matrix are equal to 0 for missing edges and either 1 or 2 for present edges, where the value 2 is used to encode *coloured edges* as detailed below. 

* For every pair of variables in $L$ there exists an homologous pair of variables in $R$, thereby identifying a pair of homologous edges. If both edges are present in the graph the corresponding off-diagonal entries of $K$ can be constrained to have equal value. These type of edges are referred to as *coloured symmetric inside-block edges* and the corresponding entries of the pdColG matrix are set to 2. On the other hand, if both edges are present but the corresponding parameters are unconstrained then they form a pair of *uncoloured symmetric inside-block edges*.

* We say that two variables are *across-block* if one variable belongs to $L$ and the other to $R$. For every pair of non-homologous across-block variables there exists an homologous pair across-block variables, thereby identifying a pair of homologous edges. If both edges are present in the graph the corresponding off-diagonal entries of $K$ can be constrained to have equal value. These type of edges are referred to as *coloured symmetric across-block edges* and the corresponding entries of the pdColG matrix are set to 2. On the other hand, if both edges are present but the corresponding parameters are unconstrained then they form a pair of *uncoloured symmetric across-block edges*.

The functions of this package make it possible to specify different types of pdRCON submodels of interest through the arguments *type* and *force.symm* which can both take as value any subvector of the character vector *c("vertex", "inside.block.edge", "across.block.edge")*; note that the names of the components can be abbreviated down, up to the first letter only, and are not case-sensitive. The argument *type* cannot be *NULL* and:

* if *type* contains the string "vertex" then coloured vertex symmetries are allowed and, if in addition also *force.symm* contains the string "vertex", then all vertices are coloured. 

* if *type* contains the string "inside.block.edge" then coloured inside- block edge symmetries are allowed and, if in addition also *force.symm* contains the string "inside.block.edge", then only coloured edges are allowed inside blocks. 

* if *type* contains the string "across.block.edge" then coloured across-block edge symmetries are allowed and, if in addition also *force.symm* contains the string "across.block.edge", then only coloured edges are allowed across blocks, with the exception of edges joining across-block homologous variables.

For the in-package documentation please check:
```{r message = FALSE, warning = FALSE}
help("pdglasso-package")
```
