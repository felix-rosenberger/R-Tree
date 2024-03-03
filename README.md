# Building and Querying an R-Tree in Python

This repository contains the input files and program for building up an R-Tree from scratch which is then further queried with a number of range queries to output the number of data points found per query. All code is run with Python 3.9. To run the program, the input files must be in the same repository than the program when it is run in the console.

## Goals and Outcomes
The goal of an R-Tree implementation is to efficiently query large amounts of data, compared to naively querying the data where each data point has to be accessed by a search algorithm. In this program, both options are tested and the R-Tree approach achieves the result 75 times faster than the naive approach.
