/** @file conserved.hpp
 * conserved interval functions
 */

#ifndef _CONSERVED_HPP_
#define _CONSERVED_HPP_

#include "genom.hpp"
//#include "structs.h" // GRAPPA (genome_struct)


using namespace std;




	// type of a pq-tree node
	// P : order of the childs is rearrangeable (round node)
	// Q : order of the childs is fix (square node)
#define Q 0
#define P 1

// orientation of a node of a pq-tree
#define PQORIENTED 1
#define PQUNORIENTED 0

using namespace std; 

	// a node of a pq-tree
typedef struct pqn{
	pqn* parent;		// the parent
	vector<pqn* > children;	// the child nodes 
	
	int idx;
	pair<int,int> i;	// the interval
	int type, 			// P, Q
		orientation;
} pqnode;


/**
 * initialise a pqnode
 * @param[in] parent the parent node
 * @param[out] n the node which should be initialised
 * @param[in] type the type of the node P/Q
 * @param[in] interval the corresponding interval
 * @param[in] orientation the orientation of the node/interval ORIENTED/UNORIENTED
 */
void init_pqnode(pqnode* parent, pqnode** n, int type, pair<int,int> interval, int orientation);

/**
 * compute the pq-tree from given irreducible conserved intervals, i.e. the components
 * @param[in] comp the given irreducible intervals = components
 * @param[in] comp_or the orientation of the components
 * @param[in] n the lenght of the underlying permutation(s)
 * @param[out] pqroot the root node of the intervaltree
 */
void pqtree(const vector<pair<int, int> > &comp, const vector<int> &comp_or, int n, pqnode **pqroot);

/**
 * determine the length of the branches in the pq-tree, where a branch is defined as 
 * the set of nodes up to, but excluding, the next node of degree > 2. the length of 
 * a branch is the number of included unoriented components. 
 * @param[in] pqroot the root of the pq-tree
 * @param[in] len used for counting -> set to 0
 * @param[out] the length of the branches 
 */
void pqtree_branches(pqnode *pqroot, int len, vector<int> &branch_lens);

/** 
 * removes all dangling oriented components and square nodes from the tree. the tree os transformed
 * to an unrooted tree. 
 * @param[in] p the root of the tree
 * @return true if the child should be removed (for internal use)
 */
bool pqtree_defoliate(pqnode **p);

/**
 * remove a pqtree from the memory
 * @param[in] p the root of the pqtree
 */
void pqtree_free(pqnode *p);

/**
 * print a pq-tree
 *@param[in] p the node where the output should start (usually the root)
 */
void pqtree_print(pqnode *p);

/**
 * move the root of the pqtree to the next node of degree > 2, therefor the tree is transformed
 * to an unrooted tree
 * @param[in] root the root of the pqtree
 */
void pqtree_reroot(pqnode **root);



/**
 * determines how much conserved intervals of genomes are destroyed
 * by the genom g
 *@param[in] genomes the genomes
 *@param[in] g the added genom
 *@param[in,out] ci_pre the conserved intervals of genomes (will be computed if not given, i.e. size = 0)
 *@param[in] circular treat genomes as circular
 *@return the number of destroyed conserved intervals
 */
int changes_conserved_intervals(
		const vector<genom> &genomes, const genom &g, 
		vector<pair<unsigned, unsigned> > &ci_pre, 
		vector<pair<unsigned, unsigned> > &ci_dif, 
		int circular);

/**
 * C function determines how much conserved intervals of genomes are destroyed
 * by the genom g
 *@param[in] genomes_ptr the genomes
 *@param[in] m number of genomes in genomes_ptr
 *@param[in] g_ptr the added genom
 *@param[in] n length of the genomes
 *@param[in,out] ci_pre_ptr the conserved intervals of genomes (will be computed if not given, i.e. size = 0)
 *@param[in,out] ci_pre_size number of conserved intervals in ci_pre_ptr
 *@param[in] circular treat genomes as circular
 *@return  the number of destroyed conserved intervals
 */
/*extern "C" int changes_conserved_intervals(
	struct genome_struct **genomes_ptr, int m,
	int *g_ptr, int n, 
	unsigned ***ci_pre_ptr, unsigned *ci_pre_size, 
	unsigned ***ci_dif_ptr, unsigned *ci_dif_size, 
	int circular);*/

/**
 * get the conserved intervals given the maximal chains
 *@param[in] mc maximal chains
 *@param[in] circular treat genomes as circular
 *@return the conserved intervals
 */
vector<pair<unsigned, unsigned> > conserved_intervals(
	const vector< vector<int> > &mc, int circular, int n);

/**
 * get the conserved intervals of two given genomes
 *@param[in] g1 input genome one
 *@param[in] g2 input genome two
 *@param[in] circular treat genomes as circular
 *@return the conserved intervals
 */
vector<pair<unsigned, unsigned> > conserved_intervals(
		const genom &g1, const genom &g2, int circular);

/**
 * get the conserved intervals of the set of given genomes
 *@param[in] genomes the input genomes
 *@param[in] circular treat genomes as circular
 *@return the conserved intervals
 */
vector<pair<unsigned, unsigned> > conserved_intervals(
		const vector<genom> &genomes, int circular);

/**
 * computes the irreducible conserved intervals of the source genome and the genome 'target'  
 * the code comes directly from the program from Jens Stoye. i only extracted things
 * i don't need here like output and computation of the breakpoints ....
 * see: On the Similarity of Sets of Permutations and its Applications to Genome
 *   Comparison. (A. Bergeron, J. Stoye)
 *@param[in] source the source genom
 *@param[in] target the target genom
 *@param[in,out] d memory
 *@return the irreducible conserved intervals
 */
vector< pair<int, int> > getIrreducibleConservedIntervals(
	const genom &source, const genom &target, vector<int> &cior);

/**
 * computes the maximal chains of irreducible conserved intervals of 
 * the source genome and the genome 'target'  
 *@param[in] source the source genom
 *@param[in] target the second genom
 *@return the maximal chains of conserved intervals
 */
vector< vector<int> > getMaximalChains(const genom &source, 
	const genom &target);

/**
 * computes the maximal chains of irreducible conserved intervals of the genome to the 
 * genomes in the target vector
 *@param targets the target genomes
 *@return the maximal chains of conserved intervals
 */
vector< vector<int> > getMaximalChains(const vector<genom> &targets);

/**
 * computes the maximal chains from the given irreducible intervals
 *@param[in] n length of the genomes
 *@param[in] ici irreducible intervals
 *@return maximal chains
 */
vector< vector<int> > getMaximalChains( unsigned n, 
	const vector< pair<int,int> > &ici);

/**
 * computes the reversals that preserve the conserved intervals
 * structure of the calling genome and the genomes in G
 *@param G the other genomes
 *@param[in,out] d memory
 *@return the vector of indices of the preserving reversals
 */
vector< pair<int,int> > getPreservingReversals(
	const vector<genom> &G);

/**
 * computes the reversals that preserve the conserved intervals 
 * structure of the source genome and the genom target
 *@param[in] source the source genome
 *@param[in] target the target genome
 *@param[in,out] d memory
 *@return the vector of indices of the preserving reversals
 */

vector< pair<int,int> > getPreservingReversals(
	const genom &source, const genom &target);
/**
 * compute the preserving reversals from the maximal chains of irreducible conserved intervals
 *@param mc the maximal chains
 *- each vector contains one maximal chain in such way that the first elements of the 
 * maximal chains are sorted in decreasing order (\f$ \forall i<j : mc[i][0] > mc[j][0]\f$)
 *- and the indices stored in each maximal chain are sorted in decreasing order (\f$ \forall k, i<j : mc[k][i] > mc[k][j]\f$)
 *.
 *@return the vector of indices of the preserving reversals
 */
vector< pair<int,int> > getPreservingReversalsFromMaxChains(
	const vector< vector<int> > &mc);

/**
 * computes the interval distance between the two given genomes
 * with the N1 & N2 conserved intervals and N conserved intervals in the union
 * = N1 * N2 - 2*N 
 *@param[in] source the source genome
 *@param[out] target the target genome
 *@param[in,out] d memory
 *@return the interval distance
 */
int intervalDistance(const genom &source, const genom &target);

/**
 * computes the sum of the intervaldistance from source to every genome in genomes
 *@param[in] source the source genome
 *@param[in] targets the target genomes
 *@param[in,out] d memory
 *@return the intervaldistance
 */
int intervalDistance(const genom &source, const vector<genom> &targets);


/**
 * compute the preserving reversal distance, i.e. the minimal number
 * of preserving reversals to transform one genome into the other
 * linear \f$ d = n + 1 - c + u \f$
 * circular .. the same therefor the genomes have to be normalised
 *@param[in] src a genome
 *@param[in] tgt a genome 
 *@param[in] circular circularity of the genomes
 *@param[in,out] d memory
 *@return the preserving reversal distance
 */
int preserving_reversal_distance(
	const genom &src, const genom &tgt, int circular);

/**
 * C function to get the preserving reversal distance
 *@param[in] src_p source genom (int array)
 *@param[in] tgt_p target genom (int array)
 *@param[in] n length of the genomes
 *@param[in] circular circularity of the genomes
 *@return the preserving reversal distance
 */
/*extern "C" int preserving_reversal_distance(
	int *src_p, int *tgt_p,
	int n, char circular);*/

/**
 * get the valid adjacencies of the conserved intervals of the genome set G
 *@param[in] G genom set
 *@param[out] dependencies stores for each signed element which signed element are implied by it
 *@return a matrix (2n+2 x 2n+2) where i,j is > 0 iff \f$(i,j)\in E \f$ corresponds to a valid adjacency
 */
vector< vector<int> > valid_adjacencies(const vector<genom> &G, 
		vector< vector< int > > &dependencies);

/**
 *C function to get the valid adjacencies of the conserved intervals of the genome set G
 *@param[in] passgenome genomes for which the computation should be done
 *@param[in] ngenomes number of genomes in passgenome
 *@param[in] ngenes number of genes in the genomes of passgenome
  *@param[in] circular circularity of the genomes
 *@param[out] dep stores for each signed element which signed element are implied by it
 *@param[out] dep_len length of each dep-list
 *@return a matrix (2n+2 x 2n+2) where i,j is > 0 iff \f$(i,j)\in E \f$ corresponds to a valid adjacency
 */
/*extern "C" int ** valid_adjacencies(
		struct genome_struct **passgenome, 
		int ngenomes, int ngenes, int circular,
		int ***dep, int **dep_len);*/

#endif
