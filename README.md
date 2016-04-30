# blossomTree

This is blossomTree R package for blossom tree graph structure estimation and density estimation. 

For more information please contact zhe.liu.uchicago@gmail.com

#### Description

Blossom tree graphical models combine the ideas behind trees and Gaussian graphical models to form a new nonparametric family of graphical models. The approach is to attach nonparanormal blossoms, with arbitrary graphs, to a collection of nonparametric trees. The tree edges are chosen to connect variables that most violate joint Gaussianity. The non-tree edges are partitioned into disjoint groups, and assigned to tree nodes using a nonparametric partial correlation statistic. A nonparanormal blossom is then grown for each group using established methods based on the graphical lasso. The result is a factorization with respect to the union of the tree branches and blossoms, defining a high-dimensional joint density that can be efficiently estimated and evaluated on test points.

#### Installation

To install the [devtools](https://cran.r-project.org/package=devtools) package:

    install.packages("devtools")
    library(devtools)
    install_github("zhejosephliu/blossomTree")

#### Usage

Please read the documentation for the main function blossomTree.

#### Reference

Zhe Liu and John Lafferty. Blossom tree graphical models. NIPS, 2014.
