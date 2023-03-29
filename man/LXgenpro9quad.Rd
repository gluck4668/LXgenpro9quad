\name{LXgenpro9quad}
\alias{LXgenpro9quad}
\title{LXgenpro9quad}
\usage{
LXgenpro9quad_p1(gene_data,protein_data,species)
LXgenpro9quad_p2(xmin,xmax,ymin,ymax)

}
\description{
LXgenpro9quad can be used to analysis the nine quad between genes and proteins.

}
\examples{

install.packages("devtools")
library(devtools)

install_github("gluck4668/LXgenpro9quad")

library(LXgenpro9quad)
??LXgenpro9quad
#---------------------------------------
data(gene_data_example)
data(protein_data_example)

#---------------------------------------

gene_data <- "gene_data.xlsx"

protein_data <- "protein_data.xlsx"

species <- " mouse" # Here, it should be "human", "mouse", or "rat"

LXgenpro9quad_p1(gene_data,protein_data,species)

xmin <- -3
xmax <- 3
ymin <- -2
ymax <- 3

LXgenpro9quad_p2(xmin,xmax,ymin,ymax)

}
