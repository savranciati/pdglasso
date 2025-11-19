#' @details
#'
#' An RCON model for paired data (pdRCON model) is a coloured Gaussian Graphical
#' Model (GGM) where the \eqn{p} variables are partitioned into a Left block
#' \eqn{L} and a Right block \eqn{R}. The two blocks are not independent, every
#' variable in the left block has an homologous variable in the right block and
#' certain types of equality `R`estrictions on the entries of the
#' `CON`centration matrix \eqn{K} are allowed. A pdRCON model is represented by
#' a Coloured Graph for Paired Data (pdColG) with a vertex for every variable
#' and where every vertex and edge is either *coloured* or *uncoloured*. More
#' details on the equality constraints of pdRCON models, submodel classes
#' of interest and on the usage of the package are given in the following.
#'
#' @section pdRCON models - terminology and relevant submodel classes:
#'
#' A pdRCON model is a Gaussian Graphical Model with additional equality
#' restrictions on the entries of the concentration matrix. In the paired data
#' framework, there are three different types of equality restrictions of
#' interest, identified with the names *vertex*, *inside-block edge* and
#' *across-block edge*, respectively. Relevant submodel classes can be specified both
#' by allowing different combinations of restriction types and by forcing
#' different types of fully symmetric structures. In this package, different
#' submodel classes are identified by the arguments `type` and `force.symm` and
#' models are represented by coloured graphs for paired data encoded in the form
#' of a `pdColG` matrix.
#'
#'   * Every variable in \eqn{L} has an homologous variable in \eqn{R} and the
#' corresponding diagonal entries of \eqn{K} can be constrained to have equal
#' value. Such entries of  \eqn{K} are represented by *coloured vertices* of the
#' independence graph whereas the unconstrained diagonal entries are represented
#' by  *uncoloured vertices*. Different types of submodel classes of interest
#' may be obtained by (i) not allowing coloured vertices, (ii) allowing both
#' coloured and uncoloured vertices and (iii) allowing only coloured vertices.
#'
#'   * For every pair of variables in \eqn{L} there exists an homologous pair of
#' variables in \eqn{R}, thereby identifying a pair of homologous edges. If both
#' edges are present in the graph the corresponding off-diagonal entries of
#' \eqn{K} can be constrained to have equal value. These type of edges are
#' referred to as *coloured symmetric inside-block edges*. Different types of
#' submodel classes of interest may be obtained by (i) not allowing coloured
#' inside-block edges, (ii) allowing both coloured and uncoloured inside-block
#' edges and (iii) allowing only coloured inside-block edges.
#'
#'   * We say that two variables are
#'   *across-block* if one variable belongs to \eqn{L}
#' and the other to \eqn{R}. For every pair of non-homologous across-block
#' variables there exists an homologous pair across-block variables, thereby
#' identifying a pair of homologous edges. If both edges are present in the
#' graph the corresponding off-diagonal entries of \eqn{K} can be constrained to
#' have equal value. These type of edges are referred to as *coloured symmetric
#' across-block edges*. Different types of submodel classes of interest may be
#' obtained by (i) not allowing coloured across-block edges, (ii) allowing both
#' coloured and uncoloured across block edges and (iii) allowing only coloured
#' across-block edges, with the exception of edges joining a variable in \eqn{L}
#' with its homologous in \eqn{R}.
#'
#' * We remark that every coloured edge belongs to a pair of coloured symmetric
#' edges, either inside- or across- blocks. On the other hand, for an uncolored
#' edge its homologous edge may or may  not be present in the graph. In the case
#' where an uncoloured edge and its homologous are both present we say that they
#' form a pair of *uncoloured symmetric edges*, either inside- or across- blocks.
#'
#'
#' @section Use of the arguments `type` and `force.symm` to specify a submodel class:
#'
#'  The functions of this package make it possible to specify different types of
#'  pdRCON submodel classes of interest through the arguments `type` and
#'  `force.symm` which can both take as value any subvector of the character
#'  vector `c("vertex", "inside.block.edge", "across.block.edge")`; note that
#'  the names of the components can be abbreviated down, up to the first letter
#'  only, and are not case-sensitive. The argument `type` cannot be `NULL` and:
#'
#'  * If `type` contains the string `"vertex"` then coloured vertex symmetries are
#'  allowed and, if in addition also `force.symm` contains the string
#'  `"vertex"`, then all vertices are coloured.
#'
#'  * If `type` contains the string `"inside.block.edge"` then coloured inside-block edge symmetries are allowed and, if in addition also `force.symm`
#'  contains the string `"inside.block.edge"`, then only coloured edges are
#'  allowed inside blocks.
#'
#'  * If `type` contains the string `"across.block.edge"` then coloured across-block edge symmetries are allowed and, if in addition also `force.symm`
#'  contains the string `"across.block.edge"`, then only coloured edges are
#'  allowed across blocks,  with the exception of edges joining a variable in
#'  \eqn{L} with its homologous in \eqn{R}.
#'
#'
#' Note that  `force.symm` is a, possibly  `NULL`, subvector of  `type`.
#' Elements of `force.symm` which are not elements of `type` are ignored.
#'
#' @section Variables position and block structure of sample covariance matrices:
#'
#'   The functions of this package assume that the positions occupied by
#'   variables follow certain rules. More specifically, the positions from
#'   \eqn{1} to \eqn{q=p/2} correspond to the first group of variables, say
#'   \eqn{L}, whereas the positions from \eqn{q+1} to \eqn{p} are associated
#'   with the second group of variables, \eqn{R}. Furthermore, the variables of
#'   the first group are ordered in the same way as those of the second group,
#'   in the sense that for every \eqn{i=1,\ldots, q} the variable in position
#'   \eqn{i} is homologous to the variable in position \eqn{q+i}. Hence, for
#'   instance, the functions that receive in input a sample covariance matrix
#'   assume that the rows and columns of this matrix are ordered according to
#'   these rules, so that it can be partitioned into four \eqn{q\times q}
#'   submatrices naturally associated with the inside- and across- block
#'   components.
#'
#' @section Model representation through the `pdColG` matrix:
#'
#'  Every pdRCON model is uniquely represented by a Coloured Graph for Paired
#'  Data (pdColG) implemented in the form of a \eqn{p\times p} symmetric matrix
#'  where every entry is one of the values 0, 1 or 2, as follows:
#'
#'  * The diagonal entries of the pdColG matrix are all equal to either 1,
#'  for uncoloured vertices, or 2, for coloured vertices.
#'
#'  * The off-diagonal entries of the pdColG matrix are equal to 0 for missing
#'  edges and either 1 or 2 for present edges, where the value 2 is used to
#'  encode coloured edges.
#'
#'
#'
#' @author Saverio Ranciati (ORCID: 0000-0001-7880-9465) and Alberto Roverato (ORCID: 0000-0001-7984-3593)
#' \cr **Maintainer**: Saverio Ranciati <saverio.ranciati2@unibo.it>
#'
#' @references Ranciati, S., Roverato, A., (2024). "On the application of Gaussian graphical models to paired data problems", *Statistics and Computing*, *34*(6), 1-19. \url{https://doi.org/10.1007/s11222-024-10513-6}
#' @references Ranciati, S., Roverato, A., (2023). On the application of Gaussian graphical models to paired data problems. *arXiv pre-print*. \url{https://arxiv.org/abs/2307.14160}
#' @references Ranciati, S., Roverato, A., Luati, A. (2021). Fused graphical lasso for brain networks with symmetries. *Journal of the Royal Statistical Society Series C: Applied Statistics*, *70*(5), 1299-1322. \url{https://doi.org/10.1111/rssc.12514}
#' @references Højsgaard, S., Lauritzen, S. L. (2008). Graphical Gaussian models with edge and vertex symmetries. *Journal of the Royal Statistical Society Series B: Statistical Methodology*, *70*(5), 1005-1027. \url{https://doi.org/10.1111/j.1467-9868.2008.00666.x}
#'
#' @docType package
#' @name pdglasso-package
#'
# @keywords internal
#' @aliases pdglasso-package NULL
"_PACKAGE"

## usethis namespace: start
#' @importFrom MASS mvrnorm
## usethis namespace: end
NULL
