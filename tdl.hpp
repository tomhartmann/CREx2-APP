#ifndef TDL_HPP_
#define TDL_HPP_

#include <limits>


/**
 * compute the smaller of:
 * a) tdl distance from origin to goal
 * b) tdl distance from goal to origin
 * @param[in] origin a genom
 * @param[in] goal a genom
 * @param[out] dir the smaller direction: 1 = a, -1 = b, 0 if equal
 */
int min_tdrl_distance(const genom &origin,const genom &goal, int &dir);

/**
 * print chains
 * @param[in] chains the chains to print
 */
void print_chains( const vector<vector<int> > &chains );


void print_tdrl( const vector<int> &tdrl, vector<string> *nmap = NULL );

void print_tdrl_scenario( const vector<vector<int> > &scenario, const genom &g, vector<string> *nmap);

/**
 * apply a tdrl; the tdrl vector specifies for each position in g if it is
 * kept in the first copy (0), or in the second copy (1)
 * @param[in] g the genome where the tdrl should be applied
 * @param[in] tdrl vector specifying the copyindex for each position
 */
void tdrl_apply( genom &g, const vector<int> &tdrl );

/**
 * compute the number of chains in g relative to id
 * @param[in] g a genome
 * @return number of chains or std::numeric_limits< unsigned >::max() in the case of a sign difference in
 * one element
 */
unsigned tdrl_chaincnt( const genom &g );

/**
 * compute the number of chains in g relative to h
 * @param[in] g a genome
 * @param[in] h another genome
 * @return number of chains or std::numeric_limits< unsigned >::max() in the case of a sign difference in
 * one element
 */
unsigned tdrl_chaincnt( const genom &g, const genom &h );

/**
 * compute chains in g wrt id
 * @param[in] g g
 * @param[out] chains the chains
 */
void tdrl_chains( const genom &g, vector<vector<int> > &chains );

/**
 * compute chains in g wrt h
 * @param[in] g g
 * @param[in] h h
 * @param[out] chains the chains
 */
void tdrl_chains( const genom &g, const genom &h, vector<vector<int> > &chains );

/**
 * compute tdrl distance from g to h
 * @param[in] g a genome
 * @param[in] g another genome
 * @return tdrl distance
 */
unsigned tdrl_distance( const genom &g, const genom &h );

/**
 * compute tdrl distance from g to id
 * @param[in] g a genome
 * @return tdrl distance
 */
unsigned tdrl_distance( const genom &g );

/**
 * compute tdrl distance from the number of chains
 * @param[in] cc number of chains
 * @return tdrl distance
 */
unsigned tdrl_distance( unsigned cc );


/**
 * compute the number of sorting tdrls given a number of chains (wo. breaking chains)
 * @param[in] c number of chains
 * @return the number of sorting of reversals
 */
unsigned tdrl_sorting_cnt( unsigned c );

/**
 * compute the number of sorting tdrls for g wrt. id (wo. breaking chains)
 * @param[in] g a genome
 * @return number of sorting tdrls
 */
unsigned tdrl_sorting_cnt( const genom &g );

/**
 * compute the number of sorting tdrls for g wrt. h (wo. breaking chains)
 * @param[in] g a genome
 * @param[in] h another genome
 * @return number of sorting tdrls
 */
unsigned tdrl_sorting_cnt( const genom &g, const genom &h );

/**
 * get all non chain breaking sorting tdrls for g wrt. id
 * @param[in] g a genome
 * @return the tdrls
 */
vector<vector<int> > tdrl_sorting( const genom &g );

/**
 * get all non chain breaking sorting tdrls for g wrt. h
 * @param[in] g a genome
 * @param[in] h another genome
 * @return the tdrls
 */
vector<vector<int> > tdrl_sorting( const genom &g, const genom &h );

/**
 * compute the number of ALL sorting tdrls for a permutation of length n with c chains;
 * if n=c the number of sorting restricted tdrls is returned
 * @garam[in] n length of the genome
 * @param[in] c number of chains
 * @return the number of sorting tdrls
 */
unsigned tdrl_sorting_cnt_all( unsigned n, unsigned c );

/**
 * compute ALL sorting tdrls for g wrt h
 * @garam[in] g a genome
 * @param[in] h another genome
 * @return the number of sorting tdrls
 */
unsigned tdrl_sorting_cnt_all( const genom &g, const genom &h );

/**
 * get all sorting tdrls for g wrt h
 * @garam[in] g a genome
 * @param[in] h another genome
 * @return the sorting tdrls
 */
vector<vector<int> > tdrl_sorting_all( const genom &g, const genom &h );


/**
 * get all sorting tdrl scenarios
 * @param[in] g a genome
 * @param[in] h another genome
 * @param[out] the scenarios
 * @param[in] bool restricted .. only restricted reversals (default: false)
 * @param[in] maxscen maximal number of scenarios to return (default: 0 = all)
 */
void tdrl_sort( const genom &g, const genom &h, vector<vector<vector<int> > > &scenarios, bool restricted = false, unsigned maxscen = 0 );







/**
 * count the number of negative blocks
 * @param quot the permuation
 * @return the number of negative blocks
 */
//int cnt_blocks(const genom &quot);

/**
 * return blocks of elements in quot with signs unequal to the signs of the corresponding elements in tgt
 * @param[in] quot the permuation
 * @param[in] tgt another permutation
 * @param[out] blocks the negative blocks
 */
void get_blocks(const genom &quot, const genom &tgt, vector<pair<int,int> > &blocks);

///**
// * get the blocks .. for tdl 1st szenarios
// * @param[in] quot the permuation
// * @param[in] pos get positive blocks if != 0; else get negative blocks
// * @param[out] blocks the negative/positive blocks
// */
//void get_blocks_inv(const genom &quot, int pos, vector<pair<int,int> > &blocks);

/**
 * check if a given tdl is a real one (asymmetric), or a transposition (symmetric).
 * the function can also consider only a sub string of the tdrl
 * @param[in] tdl the tdl (true: keep in first copy, false: keep in 2nd copy)
 * @param[in] s the start of the substring (default: 0)
 * @param[in] e the end of the substring (default: size of the tdrl)
 * @return true if it is a real tdl, false else
 */
bool is_tdl( const vector<bool> &tdl, unsigned s=0, unsigned e=std::numeric_limits< int >::max() );

/**
 * mark the elements to take in the original and the copy
 * @param[in] g a genome (the second is the id)
 * @param[out] marks the markinks
 */
void mark_tdl(const genom &g, vector<int> &marks);

/**
 * mark the elements to take in the original and the copy
 * @param[in] origin a genome
 * @param[in] goal a second genome
 * @param[out] marks the markinks
 */
void mark_tdl(const genom &origin,const genom &goal, vector<int> &marks);

/**
 * get the start and the end of the substring of the tdrl which is in normal form,
 * i.e. remove leading trues and trailing falses
 * @param[in] tdrl a tdrl
 * @param[out] start the start
 * @param[out] end the end
 */
void tdl_normalform( const vector<bool> &tdrl, int &start, int &end);

/**
 * get a random tdrl for a genome of length n
 * @param[in] n the length of the genome
 * @return the tdrl
 */
vector<bool> random_tdl( unsigned n);

/**
 * reverse the blocks and get all minimal tdrl scenarios,
 * i.e. reversal first szenarios
 * @param[in] g the genom with the negative blocks
 * @param[in] tgt the target permutation
 * @param[in] blocks the positive/negative blocks
 * @param[out] tdlsce the minimal tdl scenarios
 * @param[out] revszen the corresponding reversal scenarios
 * @param[in] used only for internal usage, to avoid recomputations
 * @param[in] maxszen maximal number of szenarios to enumerate (0 = unlimited)
 * @param[in] currz for internal usage
 * @param[out] bestd minimal number of tdrls
 */
void reverse_blocks(genom g, const genom &tgt, vector<pair<int, int> > blocks,
	vector<vector<vector<bool> > > &tdlsce, vector<vector<pair<int, int> > > &revszen,
	map<genom, bool> &used, unsigned maxszen, vector<pair<int, int> > &currz, int &bestd );

/**
 * @see reverse_blocks
 */
void reverse_blocks_inv(const genom &g, genom tgt, vector<pair<int, int> > blocks,
	vector<vector<vector<bool> > > &tdlsce,	vector<vector<pair<int, int> > > &revszen,
	map<genom, bool> &used, unsigned maxszen, vector<pair<int, int> > &currz, int &bestd );

/**
 * apply a tdl step; the keep vector specifies which genes are kept in the
 * first copy (1), and which are kept in the second copy (0)
 * @param[in,out] g the genome which should be rearranged
 * @param[in] a vector specifying if genes should be kept/deleted in the 1st/2nd copy
 */
void tdl( genom &g, const vector<bool> &keep );

///**
// * compute the distance from id to genom
// * @param[in] _genom the genome
// * @return the tdl distance
// */
//int tdl_distance(const genom &_genom);

///**
// * compute the distance from origin to goal
// * @param[in] origin the origin
// * @param[in] goal the goal
// * @return the diatance
// */
//int tdl_distance(const genom &origin,const genom &goal);


///**
// * get the tdlmedians in distance 1 if some exist
// * @param[in] g a genome
// * @param[in] h another genome
// * @param[out] median the medians (old content is deleted)
// */
//void tdl_bfmedian(const genom &g, const genom &h, vector<genom> &median );

/**
 * compute the tdl sorting szenario from the id to the permutation
 * see: Rao et al. 2005, a tdl specifies for each index!!
 * in which copy it is used .. remark that for tdrl steps different
 * from the 1st step the elemets at the indices are not the same as
 * in g
 * @param[in] g the goal
 * @param[out] the tdl szenario
 * @param[out] tdls the tdls to apply
 */
void tdl_sort(const genom &g, vector<vector<bool> > &tdls);

/**
 * compute the tdl sorting szenario from the origin to the goal
 * see: Rao et al. 2005
 * @param[in] origin the origin
 * @param[in] goal the goal
 * @param[out] szenario the tdl szenario
 * @param[out] tdls the tdls to apply
 */
void tdl_sort(const genom &origin, const genom &goal, vector<vector<bool> > &tdls);

/**
 * count the number of transpositions and tdls in a 'tdl szenario'
 * @param[in] tdls the tdl szenario
 * @return pair of the number of tdls and transpositions
 */
pair<int, int> tdl_transposition_count( const vector<vector<bool> > &tdls );


#endif /*TDL_HPP_*/
