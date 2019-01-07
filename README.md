# CREx2-APP
*CREx2* is an dynamic programming algorithm that aims to compute a weight-minimum sequence of rearrangements for arbitrary mitochondrial gene orders and the following types of weighted rearrangement operations: inversions, transpositions, inverse transpositions, and tandem duplication random loss. *CREx2* considers only rearrangement operations that preserve common intervals, i.e., groups of genes that form an interval in both given gene orders. If the common intervals of a problem instance are organized in a linear structure, then *CREx2* has a linear runtime. Otherwise, there are two modes for *CREx2*: The first mode *CREx2-ILP* has been proposed in [Hartmann et al. 2018]. It computes an exact weight-minimum rearrangement scenario within an exponential runtime (in worst case). The second mode *CREx2-APP* has been proposed in [Hartmann 2019] it computes approximated solutions efficiently. Both versions of *CREx2* as well as *CREx* proposed in [Bernt et al. 2007] are implemented in the current version of *CREx2*. 

###Caution: This implementation does not contain the first mode of CREx2! Therefore, only approximated solutions can be computed by either CREx or CREx2-APP.###

Installation:
1. Download CREx2
2. Extract crex2.tar.gz (tar -zxvf crex2.tar.gz)
3. If no Linux distribution is used, compile the *CREx2* source code using the given Makefile. Otherwise, no further installation is needed.

*CREx2* is executed from the command line using:
@./crex2 -f <input file> [options]@
use -h for more information. 

For testing whether or not *CREx2* is installed correctly enter:
i) @./crex2 -f example.fas@ or
ii) @./crex2 -f example.fas -d@.

The first command produces two scenarios of rearrangements that transform the given gene orders gene_order_1 and gene_order_2 (see example.fas) into each other. Thereby, equally weighted rearrangements are used. The second command produces a table that gives only the lengths of these scenarios.

Feel free to send bug reports or any other kind of impressions or suggestions on *CREx2* to thartmann@informatik.uni-leipzig.de.

References:
[Bernt et al. 2007]\n
M. Bernt, D. Merkle, K. Ramsch, G. Fritzsch, M. Perserke, D. Bernhard, M. Schlegel, P. Stadler, M. Middendorf
"CREx: Inferring Genomic Rearrangements Based on Common Intervals"
Bioinformatics, 23(21): 2957-2958, 2007. 

[Hartmann et al. 2018]
Tom Hartmann, Matthias Bernt, Martin Middendorf
"An Exact Algorithm for Sorting by Weighted Preserving Genome Rearrangements"
IEEE/ACM Transactions on Computational Biology and Bioinformatics (in press)
DOI 10.1109/TCBB.2018.2831661

[Hartmann 2019]
Tom Hartmann
"Models and Algorithms for Sorting Permutations with Tandem Duplication and Random Loss"
PhD Thesis
University Leipzig 2018
