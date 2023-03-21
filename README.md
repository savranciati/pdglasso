## Symmetric Graphical Lasso for Paired Data

How to install the package:

```{r message = FALSE, warning = FALSE}
library(devtools)
install_github("savranciati/pdglasso", build_vignettes=FALSE)
library(pdglasso)
```
The package contains function to fit and simulate from a coloured graphical model for paired data, using a graphical lasso approach with a set of custom penalties inducing symmetries in the estimated concentration matrix. The estimation is performed via the alternating directions method of multipliers (ADMM).
