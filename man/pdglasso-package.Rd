% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdglasso-package.R
\docType{package}
\name{pdglasso-package}
\alias{pdglasso-package}
\title{pdglasso: Graphical Lasso for Coloured Gaussian Graphical Models for Paired Data}
\description{
This package deals with RCON models for paired data and implements an Alternating Directions Method of Multipliers (ADMM) algorithm to solve the penalized likelihood method introduced by Ranciati and Roverato (2023). Also functions for the computation of maximum likelihood estimates and for the generation of simulated pdRCON models and data are provided.
}
\details{
An RCON model for paired data (pdRCON model) is a coloured Gaussian Graphical
Model (GGM) where the \eqn{p} variables are partitioned into a Left block
\eqn{L} and a Right block \eqn{R}. The two blocks are not independent, every
variable in the left block has an homologous variable in the right block and
certain types of equality \code{R}estrictions on the entries of the
\code{CON}centration matrix \eqn{K} are allowed. A pdRCON model is represented by
a Coloured Graph for Paired Data (pdColG) with a vertex for every variable
and where every vertex and edge is either \emph{coloured} or \emph{uncoloured}. More
details on the equality constraints of pdRCON models, submodel classes
of interest and on the usage of the package are given in the following.
}
\section{pdRCON models - terminology and relevant submodel classes}{


A pdRCON model is a Gaussian Graphical Model with additional equality
restrictions on the entries of the concentration matrix. In the paired data
framework, there are three different types of equality restrictions of
interest, identified with the names \emph{vertex}, \emph{inside-block edge} and
\emph{across-block edge}, respectively. Relevant submodel classes can be specified both
by allowing different combinations of restriction types and by forcing
different types of fully symmetric structures. In this package, different
submodel classes are identified by the arguments \code{type} and \code{force.symm} and
models are represented by coloured graphs for paired data encoded in the form
of a \code{pdColG} matrix.
\itemize{
\item Every variable in \eqn{L} has an homologous variable in \eqn{R} and the
corresponding diagonal entries of \eqn{K} can be constrained to have equal
value. Such entries of  \eqn{K} are represented by \emph{coloured vertices} of the
independence graph whereas the unconstrained diagonal entries are represented
by  \emph{uncoloured vertices}. Different types of submodel classes of interest
may be obtained by (i) not allowing coloured vertices, (ii) allowing both
coloured and uncoloured vertices and (iii) allowing only coloured vertices.
\item For every pair of variables in \eqn{L} there exists an homologous pair of
variables in \eqn{R}, thereby identifying a pair of homologous edges. If both
edges are present in the graph the corresponding off-diagonal entries of
\eqn{K} can be constrained to have equal value. These type of edges are
referred to as \emph{coloured symmetric inside-block edges}. Different types of
submodel classes of interest may be obtained by (i) not allowing coloured
inside-block edges, (ii) allowing both coloured and uncoloured inside-block
edges and (iii) allowing only coloured inside-block edges.
\item We say that two variables are
\emph{across-block} if one variable belongs to \eqn{L}
and the other to \eqn{R}. For every pair of non-homologous across-block
variables there exists an homologous pair across-block variables, thereby
identifying a pair of homologous edges. If both edges are present in the
graph the corresponding off-diagonal entries of \eqn{K} can be constrained to
have equal value. These type of edges are referred to as \emph{coloured symmetric
across-block edges}. Different types of submodel classes of interest may be
obtained by (i) not allowing coloured across-block edges, (ii) allowing both
coloured and uncoloured across block edges and (iii) allowing only coloured
across-block edges, with the exception of edges joining a variable in \eqn{L}
with its homologous in \eqn{R}.
\item We remark that every coloured edge belongs to a pair of coloured symmetric
edges, either inside- or across- blocks. On the other hand, for an uncolored
edge its homologous edge may or may  not be present in the graph. In the case
where an uncoloured edge and its homologous are both present we say that they
form a pair of \emph{uncoloured symmetric edges}, either inside- or across- blocks.
}
}

\section{Use of the arguments \code{type} and \code{force.symm} to specify a submodel class}{


The functions of this package make it possible to specify different types of
pdRCON submodel classes of interest through the arguments \code{type} and
\code{force.symm} which can both take as value any subvector of the character
vector \code{c("vertex", "inside.block.edge", "across.block.edge")}; note that
the names of the components can be abbreviated down, up to the first letter
only, and are not case-sensitive. The argument \code{type} cannot be \code{NULL} and:
\itemize{
\item If \code{type} contains the string \code{"vertex"} then coloured vertex symmetries are
allowed and, if in addition also \code{force.symm} contains the string
\code{"vertex"}, then all vertices are coloured.
\item If \code{type} contains the string \code{"inside.block.edge"} then coloured inside-block edge symmetries are allowed and, if in addition also \code{force.symm}
contains the string \code{"inside.block.edge"}, then only coloured edges are
allowed inside blocks.
\item If \code{type} contains the string \code{"across.block.edge"} then coloured across-block edge symmetries are allowed and, if in addition also \code{force.symm}
contains the string \code{"across.block.edge"}, then only coloured edges are
allowed across blocks,  with the exception of edges joining a variable in
\eqn{L} with its homologous in \eqn{R}.
}

Note that  \code{force.symm} is a, possibly  \code{NULL}, subvector of  \code{type}.
Elements of \code{force.symm} which are not elements of \code{type} are ignored.
}

\section{Variables position and block structure of sample covariance matrices}{


The functions of this package assume that the positions occupied by
variables follow certain rules. More specifically, the positions from
\eqn{1} to \eqn{q=p/2} correspond to the first group of variables, say
\eqn{L}, whereas the positions from \eqn{q+1} to \eqn{p} are associated
with the second group of variables, \eqn{R}. Furthermore, the variables of
the first group are ordered in the same way as those of the second group,
in the sense that for every \eqn{i=1,\ldots, q} the variable in position
\eqn{i} is homologous to the variable in position \eqn{q+i}. Hence, for
instance, the functions that receive in input a sample covariance matrix
assume that the rows and columns of this matrix are ordered according to
these rules, so that it can be partitioned into four \eqn{q\times q}
submatrices naturally associated with the inside- and across- block
components.
}

\section{Model representation through the \code{pdColG} matrix}{


Every pdRCON model is uniquely represented by a Coloured Graph for Paired
Data (pdColG) implemented in the form of a \eqn{p\times p} symmetric matrix
where every entry is one of the values 0, 1 or 2, as follows:
\itemize{
\item The diagonal entries of the pdColG matrix are all equal to either 1,
for uncoloured vertices, or 2, for coloured vertices.
\item The off-diagonal entries of the pdColG matrix are equal to 0 for missing
edges and either 1 or 2 for present edges, where the value 2 is used to
encode coloured edges.
}
}

\references{
Ranciati, S., Roverato, A., (2023). On the application of Gaussian graphical models to paired data problems. \emph{arXiv pre-print}. \url{https://arxiv.org/abs/2307.14160}

Ranciati, S., Roverato, A., Luati, A. (2021). Fused graphical lasso for brain networks with symmetries. \emph{Journal of the Royal Statistical Society Series C: Applied Statistics}, \emph{70}(5), 1299-1322. \url{https://doi.org/10.1111/rssc.12514}

Højsgaard, S., Lauritzen, S. L. (2008). Graphical Gaussian models with edge and vertex symmetries. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, \emph{70}(5), 1005-1027 \url{https://doi.org/10.1111/j.1467-9868.2008.00666.x}
}
\author{
Saverio Ranciati (ORCID: 0000-0001-7880-9465) and Alberto Roverato (ORCID: 0000-0001-7984-3593)
\cr \strong{Maintainer}: Saverio Ranciati \href{mailto:saverio.ranciati2@unibo.it}{saverio.ranciati2@unibo.it}
}
