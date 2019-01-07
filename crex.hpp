#ifndef DUPLO_HPP_
#define DUPLO_HPP_

#include <set>

#include "io.hpp"
#include "costfoo.hpp"
#include "common.hpp"
#include "helpers.hpp"
#include "rearrangements.hpp"
#include "tdl.hpp"

using namespace std;

class bpscenario : public alternative{
	private:
	/**
	 * recursively compute the rearrangement scenarios induced by breakpoints
	 * - apply the last element of current
	 * - find a reversal/transposition with the methods described in Zhao07
	 * - append the rearrangement to current
	 * - call bpscenario_add
	 * - pop current
	 * - continue finding other rearrangements
	 * @param[in] g1 origin genome
	 * @param[in] g2 target genome
	 * @param[in] current current solution of the recursion
	 * @param[in] maxalt
	 * @param[in] check each bpscenario is stored as ordered scenario we need an possibility
	 * to check if a scenario was already added, so we use a set of unordered scenarios
	 * @param[in] complete save only complete scenarios iff true
	 */
	void bpscenario_add( genom g1, const genom &g2, vector<rrrmt*> &current,
		unsigned maxalt, unordered *check, bool complete );
	public:

	/**
	 * construct an emptyscen scenario
	 */
	bpscenario();

	/**
	 * constructor: get a scenario from g1 to g2
	 * @param[in] g1 a genome
	 * @param[in] g2 another genome
	 * @param[in] maxalt maximal number of alternatives
	 * @param[in] complete save only complete scenarios iff true
	 */
	bpscenario(const genom &g1, const genom &g2, unsigned maxalt, bool complete);

	/**
	 * copy constructor
	 * @param[in] bps the original
	 */
	bpscenario(const bpscenario &bps );

	/**
	 * destructor
	 */
	~bpscenario();

//	void bpscenario_rec( const genom &g1, const genom &g2, hdata &hd );

	/**
	 * return a pointer to a new bpscenario
	 * @return the pointer
	 */
//	bpscenario* clone() const;
//
//
//	bpscenario* create() const;

//	/**
//	 * output function
//	 * @param[in] out stream to write into
//	 */
//	ostream &output(ostream &out) const;

	string typestr() const;
};

class crex2 : public alternative{
private:

	/**
	 * compute the cost and rearrangement scenario for a leaf node n
	 *
	 * @param[in] cst costs of the rearrangements
	 */
	void _crex2_leaf( itnode *n, const genom &g1, vector<float> &dps, vector<rrrmt*> &dpr,
			const costfoo *cst  );

	void _crex2_lin( itnode *n, const genom &g1, vector<float> &dps, vector<rrrmt*> &dpr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst, bool all);

	void _crex2_lin_k0( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr);
	/**
	 * determine the score and rearrangement scenario for inverting the sign of
	 * node n with one rearrangement
	 * @param[in] n root of the subtree
	 * @param[in] g1 origin genome
	 * @param[out] ms the score
	 * @param[out] mr the corresponding rearrangement
	 * @param[in] cdps scores of the children
	 * @param[in] cdpr corresponding rearrangement scenarios of the children
	 * @param[in] cst costs of the rearrangements
	 * @param[in] all get all alternatives
	 */
	void _crex2_lin_k1( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst, bool all);

	void _crex2_lin_k2( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst, bool all);

	void _crex2_lin_k3( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst, bool all);

	void _crex2_lin_all(itnode *n, const genom &g1, vector<float> &ms, vector<rrrmt*> &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst, bool all);

	/**
	 * function for handling prime nodes
	 *
	 * @param[in] n the root
	 * @param[in] g1 origin genome
	 * @param[out] ms the score
	 * @param[out] mr the rearrangement
	 * @param[out] cdps the scores of the child nodes
	 * @param[out] cdpr the rearrangement scenarios of the children
	 * @param[in] cst costs of the rearrangements
	 * @param[in] all get all scenarios
	 * @param[in] time_bound: time bound for ILP calculations in seconds
	 * @param[out] is_optimal: true iff scenarios are optimal
	 * @param[in] distance: only return scenarios with this number of rrrmts
	 */
	void _crex2_pri( itnode *n, const genom &g1, const genom &g, vector<float> &ms, vector<rrrmt*> &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst, bool all, const unsigned time_bound, bool &is_optimal,
			const unsigned distance);

	/**
	 * alternative function for handling prime nodes, using approximation algorithm for Sorting by I,T,iT,and TDRL
	 *
	 * @param[in] n the root
	 * @param[in] g1 origin genome
	 * @param[in] g identified genome
	 * @param[out] ms the score
	 * @param[out] mr the rearrangement
	 * @param[in/out] cdps the scores of the child nodes
	 * @param[in/out] cdpr the rearrangement scenarios of the children
	 * @param[in] cst costs of the rearrangements
	 */
	void _crex2_pri_approx(itnode *n, const genom &g1, const genom &g, vector<float> &ms, vector<rrrmt*> &mr,
			const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
			const costfoo *cst);

	/**
	 * determine parsimonious crex scenario for subtree
	 *
	 * @param[in] n root of the subtree
	 * @param[in] g1 origin genome
	 * @param[in] g is renamed genome (recall g1->g2 can be expressed as g2^{-1}*g1->iota, hence g:=g2^{-1}*g1)
	 * @param[out] dps dp scores for +/-
	 * @param[out] dpr corresponding rearrangement scenarios
	 * @param[in] cst costs of the rearrangements
	 * @param[in] complete solve linear nodes by enumerating all possible subsets
	 * @param[in] all get all scenarios
	 * @param[in] time_bound: time bound for ILP calculations in seconds
	 * @param[out] is_optimal: true iff scenarios are optimal
	 * @param[in] distance: only return scenarios with this number of rrrmts
	 */
	void _crex2_main( itnode *n, const genom &g1, const genom &g, vector<float> &dps, vector<rrrmt*> &dpr,
			const costfoo *cst, bool complete, bool all, const unsigned time_bound, bool &is_optimal, const unsigned distance, const bool approx);
public:
	/**
	 * construct an empty szenario
	 */
	crex2();

	/**
	 * constructor: get the scenario from g1 to g2
	 * @param[in] g1 a genome
	 * @param[in] g2 another genome
	 * @param[in] oriented treat genome as oriented
	 * @param[in] cst costs of the rearrangements
	 * @param[in] complete solve linear nodes by enumerating all possible subsets
	 * 		necessary e.g. for cost functions that do not score by type
	 * @param[in] all get all alternatives
	 * @param[in] time_bound: time bound for ILP calculations in seconds
	 * @param[in] distance: only return scenarios with this number of rrrmts
	 * @param[in] approx: use approximation algorithm for prime nodes
	 */
	crex2( const genom &g1, const genom &g2, bool oriented, const costfoo *cst, bool &is_optimal,
			bool complete = false, bool all=true, const unsigned time_bound=std::numeric_limits<unsigned>::max(),
			const unsigned distance=0, const bool approx=true);

	/**
	 * copy constructor
	 * @param[in] c the original
	 */
	crex2(const crex2 &c );

	/**
     * destructor
     */
	~crex2();

	/**
	 * @return "crex"
	 */
	string typestr() const;

};

class crex : public unordered{
	private:
	/**
	 * find single gene indels from src to tgt, i.e.:
	 * - genes not present in src that are present in tgt -> insertion
	 * - genes present in src that are not present in tgt -> deletion
	 * .
	 * additionally the elements of the genomes are renames such that we have
	 * permutations. the mapping in rnm gives the necessary renamings
	 * for getting the original elements.
	 *
	 * @param[in] src input genome
	 * @param[in] tgt input genome
	 * @param[out] srcp output genome
	 * @param[out] tgtp output genome
	 * @param[out] indels the list of insertions and deletions
	 * @param[out] rnm element re-renaming map
	 * @param[out] orgnmap the original the name map for the genomes with unequalized gene content
	 */
	void _crex_indels( const genom &src, const genom &tgt,
			genom &srcp, genom &tgtp, pair<vector<rrrmt*>,vector<rrrmt*> > &indels,
			vector<vector<int> > &rnm, vector<string> &orgnmap);

	/**
	 * find reversal patterns in the interval tree
	 * @param[in] node the remaining nodes of the interval tree
	 * @param[out] g1 origin genome
	 */
	bool _crex_reversals( itnode* node, const genom &g1);

	/**
	 * find reverse transposition patterns in the interval tree
	 * @param[in] node the remaining nodes of the interval tree
	 * @param[in] g1 origin genome
	 * @return true iff the reverse transposition pattern matched
	 */
	bool _crex_reversetranspositions(itnode* node, const genom &g1);

	/**
	 * find transposition patterns in the interval tree
	 * @param[in] node the remaining nodes of the interval tree
	 * @param[in] g1 origin genome
	 * @return true iff the transposition pattern matched
	 */
	bool _crex_transpositions( itnode* node, const genom &g1);

	/**
	 * find tdrl patterns in the interval tree
	 * @param[in] node the remaining nodes of the interval tree
	 * @param[in] g1 origin genome
	 * @param[in] g2 target genome
	 * @param[in] bps add an alternative scenario computed with the breakpoint method of ZhaoBourque07
	 * @param[in] maxszen maximal number of combined reversal tdl szenarios, no restriction if 0
	 */
	bool _crex_tdls(itnode* node, const genom &g1, const genom &g2, bool bps, unsigned maxszen);

	public:
	/**
	 * construct an empty szenario
	 */
	crex();
	/**
	 * constructor: get the scenario from g1 to g2
	 * @param[in] g1 a genome
	 * @param[in] g2 another genome
	 * @param[in] bps use the breakpoint method [ZhaoBourque07] to find rearrangements in prime nodes
	 * @param[in] maxszen maximal number of combined reversal + tdrls (default 0, i.e. no restriction)
	 * @param[in] mkalt add alternatives for transpositions and reverse transpositions
	 */
	crex(const genom &g1, const genom &g2, bool bps, unsigned maxszen = 0, bool mkalt = true);

	/**
	 * copy constructor
	 * @param[in] c the original
	 */
	crex(const crex &c );

	/**
	 * destructor
	 */
	~crex();
	/**
	 * return a pointer to a new crex scenario
	 * @return the pointer
	 */
	crex* clone() const;

	/**
	 * @return "crex"
	 */
	string typestr() const;
};

//ordered *adjust_ranges( vector<vector<bool> > &ts, itnode *nd, const genom &g, int n);
//
//void adjust_ranges( vector<vector<vector<bool> > > &ts, vector<vector<pair<int, int> > > &rs, itnode *nd, bool revfirst, const genom &g, int n, alternative *alts);


/**
 * add rearrangements (here tdrls) computed for a quotient permutation as alternative scenario
 * @param[in] ts the tdrls
 * @param[in] q the quotient permutation
 * @param[in] tgt the target quotiet permutation (for checking if the scenario is 'complete')
 * @param[in,out] the alternative scenario
 */
void add_quot_rrrmt( vector<vector<bool> > &ts, genom q, const genom &tgt, alternative *alt );

/**
 * add rearrangements (here reversal+tdrl) computed for a quotient permutation
 * as alternative scenarios,
 * @param[in] ts the tdrls
 * @param[in] rs the reversals
 * @param[in] q the quotient permutation
 * @param[in] tgt the target quotiet permutation (for checking if the scenario is 'complete')
 * @param[in] revfirst true: reversal first scenario
 * @param[in,out] the alternative scenario
 */
void add_quot_rrrmt( vector<vector<vector<bool> > > &ts,
	vector<vector<pair<int, int> > > &rs, const genom &q, const genom &tgt, bool revfirst, alternative *alt );



///**
// * get the number of common intervals of two given permutations
// * @param[in] p1 one permutation
// * @param[in] p2 the second permutation
// * @param[in] circular circular genomes ?
// * @return the number of common intervals
// */
//int cnt_intervals( vector<int> p1, vector<int> p2, int circular );
//
///**
// * get the number of breakpoints of two given permutations
// * @param[in] p1 one permutation
// * @param[in] p2 the second permutation
// * @param[in] circular circular genomes ?
// * @return the number of breakpoints
// */
//int cnt_breakpoints( vector<int> p1, vector<int> p2, int circular );
//
///**
// * get the reversal distance for two given permutations
// * @param[in] p1 one permutation
// * @param[in] p2 the second permutation
// * @param[in] circular circular genomes ?
// * @return the reversal distance
// */
//int cnt_reversals( vector<int> p1, vector<int> p2, int circular );
//
///**
// * @param[in] maxprev max number of reversals to try in a pnode (default = 0 -> all)
// */
//bool compare( genom g1, genom g2, vector<string> namemap, int circular, vector<string> &szenarios, int maxprev=0 );
//
///**
// * compare two genomes (given as permutations)
// * @param[in] p1 the first permutation
// * @param[in] p2 the second permutation
// * @param[in] namemap the mapping from the elements of the permutations to their real names
// * @param[in] circular circular genome ?
// * @param[in] norm_to normalize to element if circular
// * @param[in] maxprev max number of reversals to try in a pnode (default = 0 -> all)
// * @return the result of the comparison as html
// */
//vector<string> compare( vector<int> p1, vector<int> p2, vector<string> namemap, int circular, int norm_to, int maxprev = 0);
//
///**
// * check if two gene orders (of same length) are equal (wrapper for python)
// * @param[in] p1 gene order 1
// * @param[in] p2 gene order 2
// * @param[in] circular specify if the genomes are circular
// */
//bool equal(vector<int> p1, vector<int> p2, int circular);
//
//void find_reversals( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios );
//
////void find_reversals( list<itnode *> node, const genom &g1, const genom &g2, vector<pair<int,int> > &r);
//
//void find_reversetranspositions( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios );
//
////void find_reversetranspositions( list<itnode *> node, const genom &g1, const genom &g2, vector<vector<int> > &rt );
//
///**
// * get the tdls implied by the pnodes
// * @param[out] szenarios for each pnode: a list of szenarios (each szenario is a list of operations, specified by 3 strings type, orig, target )
// * @param[in] maxtdlrev max number of reversals to try in a pnode
// */
//void find_tdls( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios, int maxprev );
//
///**
// * find possible transpositions in an interval tree
// * - nodes with two children
// * - which have a sign different from the parent and the two children
// * @param[in] node the nodes to check
// */
//void find_transpositions( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios );
//
////void find_transpositions( list<itnode *> node, const genom &g1, const genom &g2, vector<vector<int> > &t );
//
///**
// * get a reversal / TDL szenario
// * @param[in] szenario the tdl szenario
// * @param[in] tdls the tdls to apply
// * @return the szenario vector of: operationname, origintree, resulttree
// */
//void get_tdl_szenario( genom g, const vector<vector<bool> > &tdls, const vector<pair<int, int> > &revs, itnode *node,
//	const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios);
//
/////**
//// * print a tdl html formatted
//// * @param[in] szenario the tdl szenario
//// * @param[in] tdls the tdls to apply
//// */
//////void print_tdl_szenario_html( genom g, const vector<vector<bool> > &tdls, const vector<pair<int, int> > &revs, itnode *node,
//////	const genom &g1, const vector<string> &nmap, stringstream &out );
////
//////void print_tdl( const genom &source, const genom &target, const genom &src, const genom &tgt, itnode *pnode, const vector<string> &namemap, ostream &o );
//
//pair<int, int> translate_interval(const genom &g1, const pair<int,int> &i1, const genom &g2 );
//
///**
// * get the text tree (parenthesised)
// * @param[in] p1 one permutation
// * @param[in] p2 the second permutation
// * @param[in] namemap mapping from indices to gene names (size n+1 !!)
// * @param[in] circular circular genomes ?
// * @param[in] normalize to ?
// * @return the dot description
// */
//string tree(vector<int> p1, vector<int> p2, vector<string> namemap, int circular, int norm_to);

/////**
//// * get the dot description of the strong intervals tree of two given permutations
//// * @param[in] p1 one permutation
//// * @param[in] p2 the second permutation
//// * @param[in] namemap mapping from indices to gene names (size n+1 !!)
//// * @param[in] circular circular genomes ?
//// * @param[in] normalize to ?
//// * @return the dot description
//// */
//////string tree_dot( vector<int> p1, vector<int> p2, vector<string> namemap, int circular, int norm_to );
////
/////**
//// * get the the strong intervals tree of two given permutations as html table
//// * @param[in] p1 one permutation
//// * @param[in] p2 the second permutation
//// * @param[in] namemap mapping from indices to gene names (size n+1 !!)
//// * @param[in] circular circular genomes ?
//// * @param[in] normalize to ?
//// * @return the dot description
//// */
//////string tree_table( vector<int> p1, vector<int> p2, vector<string> namemap, int circular, string prefix, int norm_to );
////
#endif/*DUPLO_HPP_*/
