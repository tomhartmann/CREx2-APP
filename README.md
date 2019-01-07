# CREx2-APP
*CREx2* is an dynamic programming algorithm that aims to compute a weight-minimum sequence of rearrangements for arbitrary mitochondrial gene orders and the following types of weighted rearrangement operations: inversions, transpositions, inverse transpositions, and tandem duplication random loss. *CREx2* considers only rearrangement operations that preserve common intervals, i.e., groups of genes that form an interval in both given gene orders. If the common intervals of a problem instance are organized in a linear structure, then *CREx2* has a linear runtime. Otherwise, there are two modes for *CREx2*: The first mode *CREx2-ILP* has been proposed in [Hartmann et al. 2018]. It computes an exact weight-minimum rearrangement scenario within an exponential runtime (in worst case). The second mode *CREx2-APP* has been proposed in [Hartmann 2019] it computes approximated solutions efficiently. Both versions of *CREx2* as well as *CREx* proposed in [Bernt et al. 2007] are implemented in the current version of *CREx2*. Caution: This implementation does not contain the first mode of CREx2! Therefore, only approximated solutions can be computed by either CREx or CREx2-APP. 
