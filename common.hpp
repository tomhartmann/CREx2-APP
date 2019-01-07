/** @file common.hpp
 * common interval functions
 */

#ifndef _COMMON_HPP_
#define _COMMON_HPP_

#include <list>

//#include "dstnc_inv.hpp"
#include "genom.hpp"
//#include "structs.h"	// GRAPPA (genome_struct)

// output trivial common intervals

#define WITHTRIV 1
#define WITHOUTTRIV 0

// use sign information
#define SIGNED 1
#define UNSIGNED 0

// properties of the interval tree nodes
#define UNK -1

#define LIN 0
#define PRI 1

#define INC 0
#define DEC 1


/**
 * the functions for computing (strong) common intervals, interval trees
 * as described in : A. Bergeron, C. Chauve, F. de Montgolfier and M. Raffinot.
 * Computing common intervals of K permutations, with applications to modular decomposition of graphs
 * In 13th Annual European Symposium on Algorithms (ESA), Oct. 2005.
 *
 * canonical_generator
 * common_intervals
 * generator
 * interval_tree
 * strong_intervals
 * support
 *
 * funtions for computation of the perfect distance
 * as described in : Severine Berard, Anne Bergeron, Cedric Chauve, Christophe Paul:
 * Perfect sorting by reversals is not always difficult (extended abstract).
 * WABI 2005, volume 3692 des Lecture Notes in Bioinformatics, 2005.
 *
 * interval_tree_properties
 * interval_tree_sign
 * interval_tree_type
 * perfect_distance
 * prime_distance
 * quotient_permutation
*/

	// a node of an interval tree
typedef struct itnode_struct{
	itnode_struct* parent;		// the parent
	vector<itnode_struct* > children;	// the child nodes

	pair<int,int> i;	// the interval
	int type, 			// LIN, PRI
		sign;			// POS, NEG
} itnode;


/**
 * compute the canonical generator from a given generator
 *@param[in,out] r the r vector of the generator, afterwards the r part of the canonical generator
 *@param[in,out] l the l vector of the generator, afterwards the l part of the canonical generator
 */
void canonical_generator(vector<int> &r, vector<int> &l);


/**
 * compare two intervals by their first index
 * @param[in] a the 1st interval
 * @param[in] b the 2nd interval
 * @return true if a < b
 */
bool cmp_interval_first( pair<int, int> a, pair<int, int> b );

/**
 * compare two intervals by their second index
 * @param[in] a the 1st interval
 * @param[in] b the 2nd interval
 * @return true if a < b
 */
bool cmp_interval_second( pair<int, int> a, pair<int, int> b );

/**
 * combines two generators. asumption: the generators are computed relative to the
 * same permutation (otherwise it makes no sense)
 *@param[in,out] r the r part of the first generator, afterwards the r part of the combined generator
 *@param[in,out] l the l part of the first generator, afterwards the l part of the combined generator
 *@param[in,out] radd the r part of the 2nd generator
 *@param[in,out] ladd the l part of the 2nd generator
 */
void combine_generator(vector<int> &r, vector<int> &l,
	const vector<int> &radd, const vector<int> &ladd);

/**
 * compute the generator of the common intervals of the given genomes
 * @param[in] genomes the genomes (the first genome in the vector has to be the identity)
 * @param[in] n length of the genomes
 * @param[out] r the r part of the generator (must be of size n+1)
 * @param[out] l the l part of the generator (must be of size n+1)
 */

void combined_generator(const vector<genom> &genomes, int n, vector<int> &r, vector<int> &l);

/**
 * common intervals are computed given a generator
 * @param[in] r the r part of the generator
 * @param[in] l the l part of the generator
 * @param[in] circular true if circular genome
 * @param[in] trivial output the trivial common intervals? WITHOUTTRIV / WITHTRIV
 * 	note intervals of size 1 are included (contrary to the definition of Stoye,
 *  and contrary to common_old)
 * @param[in] sign output only signed common intervals
 * @param[out] comint the common intervals
 */
void common_intervals(const vector<int> &r, const vector<int> &l,
	int circular, int trivial, int sign, vector<pair<int, int> > &comint);

/**
 * compute the common intervals of g1 and g2; the indices are computed relative to g1
 * note: to get correct results for circular genomes both have to be normalised to 1
 * @param[in] g1 the first genome
 * @param[in] g2 the second genome
 * @param[in] circular true if circular genome
 * @param[in] trivial output the trivial common intervals? WITHOUTTRIV / WITHTRIV
 * 	note intervals of size 1 are included (contrary to the definition of Stoye,
 *  and contrary to common_old
 * @param[in] sign output only signed common intervals
 * @param[out] comint the common intervals
 */
void common_intervals(const genom &g1, const genom &g2, int circular,
	int trivial, int sign, vector<pair<int, int> > &comint);

/**
 * compute the common intervals of a vector of genomes; the indices are computed relative
 * to the first genome
 * note: to get correct results for circular genomes all have to be normalised to 1
 * @param[in] genomes the first genome
 * @param[in] circular true if circular genome
 * @param[in] trivial output the trivial common intervals? WITHOUTTRIV / WITHTRIV
 * 	note intervals of size 1 are included (contrary to the definition of Stoye,
 *  and contrary to common_old
 * @param[in] sign output only signed common intervals
 * @param[out] comint the common intervals
 */
void common_intervals(const vector<genom> &genomes, int circular, int trivial, int sign, vector<pair<int, int> > &comint);

///**
// * this function determines the elements which destroy the common intervals of the permutations
// * in 'genomes' if permutation 'gadd' is added.
// * @param[in] genomes a vector of genomes
// * @param[in] gadd a genome to add
// * @param[in] n length of the genomes
// * @param[in] circular circular genomes ??
// * @param[out] dest_elements which elements destroyed the cis
// * @param[out] comint_dif the destroyed common intervals
// */
//void common_intervals_diff(const vector<genom> &genomes, const genom &gadd,
//	int n, int circular, int sign);

/**
 * a c wrapper for the function just above. look there for a description
 * @param[in] genomes_ptr some genomes (m of length n)
 * @param[in] m number of genomes in genomes_ptr
 * @param[in] g_ptr another genome
 * @param[in] n length of the genomes
 * @param[in] circular circular genomes
 * @param[in] sign .. not implemented
 * @param[in] add see above
 * @return the number of destroyed intervals
 */
/*
extern "C" int common_intervals_diff_c(
	struct genome_struct **genomes_ptr, int m,
	int *g_ptr, int n, int circular, int sign, int add);*/


/**
 * computes the number of common intervals of 'genomes' which are destroyed if an
 * additional genome is added. the function should be called:
 * 1. with the genomes vector and an genome compliant with the common intervals of the
 *	genomes (just take an element of genomes) and the add parameter set to zero;
 *	the function computes the generator and the common intervals of genomes, which
 *	are stored for further usage until the function is again called with add = 0
 * 2. add = 1, an additional genome 'gadd' and the genomes vector (wich is not needed internally)
 * 	than the function determines how many common intervals are destroyed
 * @param[in] genomes a vector of genomes
 * @param[in] gadd a genome to add
 * @param[in] n length of the genomes
 * @param[in] circular circular genomes ??
 * @param[in] sign output only signed common intervals
 * @param[in] add control the behaviour of the function (.. haha, therefor its a parameter)
 *	.. if 0: build common intervals and generator of genomes and store for later usage
 *	.. if 1: add the genome and compute the difference
 * @return the number of destroyed (nontrival) common intervals intervals if add = 1; and
 *	if add= 0 just the number of common intervals of genomes
 */
int common_intervals_diff(const vector<genom> &genomes, const genom &gadd,
	int n, int circular, int sign, int add);

/**
 * get the genome (permutation) defined by the traversal of the leaves
 * @param[in] n the interval tree
 * @param[in] g the genome
 * @param[out] f the order of the leaves
 */
void front(itnode *n, const genom &g, genom &f);
/**
 * compute the generator of the common intervals of a genome (and the identity)
 *@param[in] p the given genome
 *@param[out] sup the r part of the generator
 *@param[out] inf the r part of the generator
 */
void generator(const genom &p, vector<int> &sup, vector<int> &inf);

/**
 * get the generator for a vector of genomes
 *@param[in] g the genomes
 *@param[out] r the r part of the generator
 *@param[out] l the l part of the generator
 */
//~ void generator(const vector<genom> &g, vector<int> &r, vector<int> &l);


/**
 * init an node of the interval tree
 * @param[in] i the interval to set
 * @param[in] parent the parent to set
 * @param[in,out] n the node to allocate
 */
void init_itnode(pair<int, int> &i, itnode* parent, itnode** n);

/**
 * compute the interval tree from strong intervals
 *@param[in] strong_int the given strong intervals
 *@param[in] n the lenght of the underlying permutation(s)
 *@param[in] g the underlying genome
 *@param[out] itroot the root node of the intervaltree
 */
void interval_tree(vector<pair<int, int> > &strong_int, int n,
	const genom &g, itnode **itroot);

/**
 * get the interval tree for two permutations; the intervals are given relative to g1
 * @param[in] g1 1st genome
 * @param[in] g2 2nd genome
 * @param[in] n length of the genomes
 * @param[out] tree_root the root of the tree
 */
void interval_tree(const genom &g1, const genom &g2, int n, itnode **tree_root);


/**
 * get the interval tree for a set of permutations; the intervals are given relative to the first genome
 * @param[in] genomes the set of genomes
 * @param[in] n length of the genomes
 * @param[out] tree_root the root of the tree
 */
void interval_tree(vector<genom> genomes, int n, itnode **tree_root );

/**
 * recursively copy an interval tree
 * @param[in] orig the original tree
 * @param[out] copy the copy
 */
void interval_tree_copy( itnode *orig, itnode **copy );

/**
 * get the depth of an interval tree; leaves have depth 0,
 * @param[in] p the node to start with
 * @return the depth
 */

int interval_tree_depth( itnode *p );

/**
 * free an interval tree
 *@param[in] itroot the root of the interval tree
 */
void interval_tree_free(itnode *itroot);

/**
 * get the siblings of a node
 * @param[in] p the node
 * @param[out] siblings the siblings
 */
void interval_tree_get_siblings( itnode *p, vector<itnode *> &siblings );

/**
 * check if the strong interval tree for g1,g2 is prime
 * @param[in] g1 a genome
 * @param[in] g2 another genome
 * @return true iff the tree contains a prime node
 */
bool interval_tree_is_prime( const genom &g1, const genom &g2 );

/**
 * get the leafs of an interval tree
 * @param[in] p the root of the tree
 * @param[out] node the nodes
 */
void interval_tree_leafs(itnode *p, vector<itnode *> &node);

/**
 * get the linear node of an interval tree
 * @param[in] p the start node
 * @param[out] lnode the linear nodes
 */
void interval_tree_linearnodes(itnode *p, list<itnode *> &lnode);

/**
 * get the nodes of an interval tree in postorder
 * @param[in] p the root of the tree
 * @param[out] node the nodes
 * @param[in] leaves get the leaf nodes iff True
 */
//void interval_tree_nodes(itnode *p, list<itnode *> &node, bool leaves=true);
template<typename _OutputIterator>
void interval_tree_nodes(itnode *p, _OutputIterator &__result, bool leaves){
	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_nodes( p->children[i], __result, leaves);
	}
	if( p->children.size() > 0 || leaves == true ){
	    *__result = p;
	    ++__result;
	}
}

/**
 * get the components of prime nodes of a interval tree
 * @param[in] p the tree (root)
 * @param[out] pcomp the components
 */
void interval_tree_pnode_components(itnode *p, vector<vector<itnode *> > &pcomp);

/**
 * given a top node of a prime component -> get the pnodes connected to it
 * @param[in] p the top node of a component
 * @param[out] pcomp the component
 */
void interval_tree_pnode_fill(itnode *p, vector<itnode *> &pcomp);

/**
 * get the top nodes of components
 * @param[in] p the root of an interval tree
 * @param[out] ptop the top nodes
 */
void interval_tree_pnode_top(itnode *p, vector<itnode *> &ptop);

/**
 * get the primenodes of an interval tree
 * @param[in] p the root of the tree
 * @param[out] pnode the prime nodes
 */
 void interval_tree_primenodes(itnode *p, list<itnode *> &pnode);

/**
 * print an interval tree
 * @param[in] p the node where the output should start (usually the root)
 * @param[in] g the source genom of sorting
 * @param[out] o the stream for the output
 * @param[in] prefix an optional prefix (applies only for the topmost node, i.e. not used in the recursion)
 */
void interval_tree_print(itnode *p, const genom &g, ostream &o, string prefix="");

/**
 * print the prefix for a node, i.e. the opening paranthesis etc.
 * @param[in] p the node
 * @param[out] o the stream for writing
 * @param[in] prefix an optional prefix for the prefix :-)
 */
void interval_tree_print_prefix( itnode *p, ostream &o, string prefix="");

/**
 * print the suffix for a node, i.e. the closing paranthesis
 * @param[in] p the node
 * @param[out] o the stream for writing
 */
void interval_tree_print_suffix( itnode *p, ostream &o );

///**
// * print an interval tree in dot format
// * @param[in] p the node where the output should start (usually the root)
// * @param[in] g the source genom of sorting
// * @param[in] fname filename
// * @param[in] ext the extension of the file
// * @param[in] namemap mapping from the integers representing the genes to strings
// */
//void interval_tree_print_dot(itnode *p, const genom &g, string fname, string ext, const vector<string> &namemap);

/**
 * print an interval tree in dot format
 * @param[in] p the node where the output should start (usually the root)
 * @param[in] g the source genom of sorting
 * @param[in] o the stream for the output
 * @param[in] namemap mapping from the integers representing the genes to strings
 */
void interval_tree_print_dot(itnode *p, const genom &g, ostream &o, const vector<string> &namemap=vector<string>());

/**
 * print an interval as a html table
 * @param[in] p the node where the output should start (usually the root)
 * @param[in] g the source genom of sorting
 * @param[in] namemap assings a name (string) to each element of the permutation
 * @param[in] prefix each table gets an id: prefix_istart_iend (where istart and iend are the start and end indices) of the corresponding interval
 * @param[in] o the stream for the output
 * @param[in] inv print the children of each node inverted (defaults to false)
 */
void interval_tree_print_table(itnode *p, const genom &g, const vector<string> &namemap, const string prefix, ostream &o, bool inv=false );

///**
// * get the perfect distance, except for 'hard to sort' nodes (prime nodes without linear parent)
// *@param[in] p the node of the tree where the computation schould start (root)
// *@param[in] g the underlying genome
// *@param[in] n the length of g
// *@param[out] dist the computed distance
// * @param[in] rd reversal distance
// */
//void interval_tree_properties(itnode *p, const genom &g, unsigned n, int &dist, dstnc_inv *rd);

/**
 * init a random strong interval tree consisting of only linear nodes
 * number of children is chosen at random from a given range
 * @param[in,out] r the root node (is allocated inside of the function)
 * @param[in,out] p the parent (init usually with NULL)
 * @param[in] i interval that should be represented by the root (init with (0,n-1))
 * @param[in] dm min out degree
 * @param[in] dM max out degree
 */
void interval_tree_random( itnode **r, itnode *p, pair<int,int> i, unsigned dm, unsigned dM  );

/**
 * simulate the effects of a reversal of the interval corresponding to the node
 * i.e. recursively 1) reverse the order of the children, switch sign : inc<->dec
 * note. because the leafs store only the interval the correct output is
 * handled by the print foo
 * @param[in] n the node
 */
void interval_tree_reverse(itnode *n);

/**
 * get the sign of the prime nodes of an interval tree, if possible
 *@param[in,out] n the starting node (usually the root)
 *@param[in] g the genome(s) underlying the interval tree
 */
void interval_tree_sign( itnode *n, const genom &g);

/**
 * get all sorting preserving reversals
 * @param[in] p a node of the strong interval tree
 * @param[in] g a genome
 * @param[in] n the length of the genome
 * @param[out] reversals the sorting reversals
 * @param[in] rd reversal distance
 */
//void interval_tree_sorting_reversals(itnode *p, const genom &g, int n, vector<pair<int,int> > &reversals, dstnc_inv *rd );

/**
 * get all sorting preserving reversals with the pqtree data structure
 * @param[in] g1 a genome
 * @param[in] g2 another genome
 * @param[out] reversals the sorting reversals
 */
void interval_tree_sorting_reversals(const genom g1, const genom g2, vector<pair<int,int> > &reversals);

/**
 * recursively switch all signs inc->dec, dec->inc
 * @param[in] p the node where the switching should start
 */
void interval_tree_switch_sign( itnode *p );

/**
 * recursively calculates the subpermutation of a given node
 * @param[in] p the node where subpermutation shall be calculated
 * @param[in] complete starting genome
 * @param[in/out] g the subpermutation
 */
void interval_tree_get_subpermutation(itnode *p, const genom start, genom &g, const unsigned child_numb);

/**
 * get the non trivial irreducible common intervals given a strong interval tree:
 * each prime node is a irreducible common interval + the union of two
 * adjacent childnodes of a linear node is a irreducible common interval
 *
 * @param[in] nd the root of the tree
 * @param[out] ii the non trivial irreducible intervals
 */
void irreducible_intervals(itnode *nd, vector<pair<int,int> > &ii);
/**
 * get the non trivial irreducible common intervals of two genomes
 * @param[in] g1 a genome
 * @param[in] g2 another genome
 * @param[in] n the length of the genomes
 * @param[out] ii the non trivial irreducible intervals
 */
void irreducible_intervals(const genom &g1, const genom &g2, int n, vector<pair<int,int> > &ii);

/**
 * get the non trivial irreducible common intervals of a set of genomes
 * @param[in] genomes the set of genomes
 * @param[in] n the length of the genomes
 * @param[out] ii the non trivial irreducible intervals
 */
void irreducible_intervals(vector<genom> genomes, int n, vector<pair<int,int> > &ii);


//void prime_sorting_reversals(itnode *p, int n, const genom &g, vector<pair<int,int> > &reversals, dstnc_inv *rdist);

/**
 * get the type (linear or prime) of each node of an interval tree
 * and assign signs to linear nodes
 * @param[in,out] pqroot the root of the interval tree
 * @param[in] r the r part of the generator
 * @param[in] l the l part of the generator
 * @param[in] g the underlying genome
 */
void interval_tree_type( itnode *pqroot, const vector<int> &r, const vector<int> &l, const genom &g );

///**
// * compute the perfect distance between g1 and g2, i.e. the minimal number
// * of reversals, which do not destroy the common intervals of g1 and g1,
// * needed to transform g1 into g2
// * note: to get correct results for circular genomes both have to be normalised to 1
// *@param[in] g1 a genome
// *@param[in] g2 anoter genome
// *@return the perfect distance
// */
//int perfect_distance(const genom g1, const genom g2);


///**
// * compute the perfect distance between g1 and g2 where the common intervals
// * of genomes have to be preserved
// * @param[in] g1 one genome
// * @param[in] g2 another genome
// * @param[in] genomes a vector of genomes
// * @param[in] n the length of all genomes
// * @return the perfect distance between g1 and g2, where the common intervals of genomes are preserved
// */
//int perfect_distance_glob(const genom &g1, const genom &g2, const vector<genom> &genomes, int n);


///**
// * computes the sum of the perfect distances between g1 and the genomes in
// * genomes, where the common intervals of \f${g_1, genomes}$\f are
// * preserved, the normal use of this function is to compute the distance sum
// * of some permutation to their (common interval preserving) median
// * @param[in] g1 as genome (usualy the median of genomes, which preserves the common intervals of genomes)
// * @param[in] genomes a set of genomes
// * @param[in] n the length of the genomes
// * @return the perfect distance
// */
//int perfect_distance_glob_sum(const genom &g1, const vector<genom> &genomes, int n );


/**
 * get the quotient permutation of a subtree
 * (note that the resulting permutation is unsigned, the signs have to be added elsewhere)
 *
 * @param[in] p root of subtree
 * @param[in] g underlying permutation
 * @param[in,out] quotient permutation to be filled
 * @param[in] circular consider the genome as circular (i.e. the tree as unrooted, so one of the siblings of the node
 * contributes to the quotient permutation)
 * @param[in] sign get the signed quotient permutation
 */
void quotient_permutation(itnode *p, const genom &g, genom &quotperm, int circular, int sign);

/**
 * get the quotient permutation of a permutation, i.e. given a permutation \f$ \pi \f$ over
 * \f$ k \f$ elements (e.g.: \f$ (3 6 1 4 )\f$ ) it returns  a permutation \f$ \pi' \f$ over the elements
 * \f$ \{1,\ldots, k \} \f$, so that \f$ |\pi_i| < |\pi_j| \rightarrow \pi'_i < \pi'_j \f$ (e.g.: \f$ (2 4 1 3 )\f$ ).
 * (note that the resulting permutation is unsigned, the signs have to be added elsewhere)
 *
 * @param[in,out] qout in: a permutation of n elements; out: a permutation over \f$ \{1...n\}\f$
 * @param[in] n the legth of the underlying permutation
 */
void quotient_permutation(genom &quot, int n);

/**
 * sorts a vector of intervals. the sorting is done with bucket sort in O(n)
 * uniq strong common intervals are removed (they have to be neighbours in the vector)
 *
 *@param[in] n the length of the underlying permutations, in subsequent calls should
 * 	never be an interval with one of the indices greater than n
 *@param[in,out] intervals the intervals to sort
 */
void sort_intervals(int n, vector<pair<int, int> > &intervals);

/**
 * get the strong intervals from a given (canonical) generator
 *@param[in] r the r part of the generator
 *@param[in] l the l part of the generator
 *@param[out] strong_int the strong intervals (has to be of size 2n), contains duplicates
 */
void strong_intervals(const vector<int> &r, const vector<int> &l, vector<pair<int,int> > &strong_int, bool unique=false, bool trivial = true);

/**
 * get the strong intervals of two permutations of length n
 * @param[in] g1 the first permutation
 * @param[in] g2 the second permutation
 * @param[in] n the length of the genomes
 * @param[out] strong_int the strong intervals (has to be of size 2n), no duplicates
 * @param[in] unique remove duplicate strong intervals
 */
void strong_intervals(const genom &g1, const genom &g2, int n, vector<pair<int,int> > &strong_int, bool unique=true, bool trivial = true);

/**
 * get the strong intervals of two permutations of length n
 * @param[in] g1 the first permutation
 * @param[in] g2 the second permutation
 * @param[in] n the length of the genomes
 * @param[out] r the r of the canonical generator (has to be of size n+1)
 * @param[out] l the l of the canonical generator (has to be of size n+1)
 * @param[out] strong_int the strong intervals (has to be of size 2n), contains duplicates
 * @param[in] uniq remove duplicate strong intervals
 */
void strong_intervals(const genom &g1, const genom &g2, int n, vector<int> &r, vector<int> &l,
	vector<pair<int,int> > &strong_int, bool uniq=false, bool trivial = true);

/**
 * get the strong intervals of a set of permutations of length n
 * @param[in] genomes the set of genomes
 * @param[in] n the length of the genomes
 * @param[out] strong_int the strong intervals (has to be of size 2n), contains duplicates
 * @param[in] uniq remove duplicate strong intervals
 */
void strong_intervals(const vector<genom> &genomes, int n, vector<pair<int,int> > &strong_int, bool uniq=true, bool trivial = true);

/**
 * get the strong intervals of a set of permutations of length n, additionally
 * the generator vectors are returned (can be used later, e.g. for interval tree
 * construction (node type recognition))
 * @param[in] genomes the set of genomes
 * @param[in] n the length of the genomes
 * @param[out] r the r of the canonical generator (has to be of size n+1)
 * @param[out] l the l of the canonical generator (has to be of size n+1)
 * @param[out] strong_int the strong intervals (has to be of size 2n), contains duplicates
 * @param[in] uniq remove duplicate strong intervals
 */
void strong_intervals(const vector<genom> &genomes, int n, vector<int> &r, vector<int> &l,
	vector<pair<int,int> > &strong_int, bool uniq=false, bool trivial = true);

/**
 * get the supports of r and l of a generator
 *@param[in] r the r part of the generator
 *@param[in] l the l part of the generator
 *@param[in,out] rsup the support of r
 *@param[in,out] lsup the support of l
 */
void support(const vector<int> &r, const vector<int> &l, vector<int> &rsup, vector<int> &lsup);


// these functions are not public
///**
// * normalises circular chromosomes in such a way that the linear
// * find_common_intervals finds all needed irreducible common intervals
// * in order to execute valid_adjacencies_common and srac correctly.
// * the function assumes that the chromosomes are normalised as circular,
// * i.e. are cutted before an arbitrary element
// *@param[in] g1 unnormalised genome one
// *@param[in] g2 unnormalised genome one
// *@param[in] g1n normalised genome one
// *@param[in] g2n normalised genome one
// *@param[in] sign use sign information to get the common intervals
// */
//void normalize_common(const genom &g1, const genom &g2,
//	genom &g1n, genom &g2n, int sign);
///**
// * normalises circular chromosomes in such a way that the linear
// * find_common_intervals finds all needed irreducible common intervals
// * in order to execute valid_adjacencies_common and srac correctly.
// * the function assumes that the chromosomes are normalised as circular,
// * i.e. are cutted before an arbitrary element
// *@param[in] genomes the genomes to normalize
// *@param[in] sign use sign information to get the common intervals
// *@return the normalised genomes
// */
//vector<genom> normalize_common(vector<genom> &genomes, int sign);
/**
 * find overlapping irreducible common intervals
 *@param[in] com_int the irreducible common intervals
 *@param[out] overlapping_int for each overlap-component a vector of intervals indices
 *@return always one :-)
 */
int overlapping_ici(vector<pair<int,int> > &com_int,
	vector<vector<unsigned> > &overlapping_int);

/**
 * find the valid adjacencies from the common intervals of the genomes
 * in g;
 * ATTENTION: the genomes in g may be returned shifted !! This
 * should only happen for circular genomes
 *@param[in,out] genomes the genomes
 *@param[in] sig if true use sign information for the common intervals
 *@param[in] circular treat genomes as circular
 *@param[out] element_range_map stores for each element the index of the range it belongs to
 *@param[out] range_size stores for each range how many elements it has
 *@param[out] range_max stores the max range index per overlap component
 *@return a matrix (2n+2 x 2n+2) where i,j is > 0 iff \f$(i,j)\in E \f$ corresponds to a valid adjacency
 */
vector< vector<int> > valid_adjacencies_common(
	vector< genom > &genomes,
	int sig, int circular,
	vector<vector<int> > &element_range_map,
	vector<vector<int> > &range_size,
	vector<int> &range_max);

/**
 * C function to find the valid adjacencies from the common intervals of the
 * genomes in g; especially (only) for calls from GRAPPA;
 *@param[in] passgenome genomes for which the computation should be done
 *@param[in] ngenomes number of genomes in passgenome
 *@param[in] ngenes number of genes in the genomes of passgenome
 *@param[in] sig if 1 then the sign information is used for common interval computation
 *@param[in] circular treat genomes as circular
 *@param[out] element_range_map stores for each element the index of the range it belongs to
 *@param[out] range_size stores for each range how many elements it has
 *@param[out] range_max stores the max range index per overlap component
 *@param[out] overlap_cmp_cnt number of components
 *@return a matrix (2n+2 x 2n+2) where i,j is > 0 iff \f$(i,j)\in E \f$ corresponds to a valid adjacency
 */
/*extern "C" int ** valid_adjacencies_common(
		struct genome_struct **passgenome,
		int ngenomes,
		int ngenes,
		int sig, int circular,
		int ***element_range_map,
		int ***range_size,
		int **range_max,
		int *overlap_cmp_cnt);*/
#endif
