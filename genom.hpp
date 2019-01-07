/**
 * @file genom.hpp
 * genom header file
 * @author M. Bernt
 */

#ifndef _GENOM_HPP_
#define _GENOM_HPP_

#include <iterator>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <stdexcept>

//// necessary GRAPPA includes for reversal distance
//#include "structs.h"
//#include "invdist.h"
//#include "med_util.h"
//
//// defines for common interval functions
//#define GENOMES_NO_END 0
//#define GENOMES_LINEAR +1
//#define GENOMES_CIRCULAR -1
//#define NO_NEW_BREAKPOINT -1

using namespace std;

struct hdata;
class genom;


/**
 * for treating signed gene orders .. unsigned / doubled
 */
enum SgnHandType { DOUBLE, UNSIGN };


/**
 * class for a genom.
 */
class genom{
	private:
		vector<string> *nmap;
	public:
		/*! name of the genom */
		//~ string name;
		/*! chromosom of the genom*/
		vector<int> chromosom;
		/*! flag thats true if the genom is circular */
		char circular;

		/**
		 * constructs an empty new genom
		 */
		genom();

		/**
		 * constructs a new idendity genom with length cnt
		 *@param cnt length of the genom
		 *@param c circularity of the genom
		 *@param[in] inv true -> construct inverse of the id
		 */
		genom(int cnt, char c, bool inv=false);

		//~ genom(string n, int cnt, bool c);
		/**
		 * contructs a new genome with the data from a c-style array
		 * @param genes c-style array of genes (length must be cnt)
		 * @param cnt length of the genom
		 * @param c circularity of the genom
		 * @param[in] nm namemap pointer to the name map or NULL (default: NULL)
		 */
		//~ genom(string n, int *genes, int cnt, bool c);
		genom(int *genes, int cnt, char c, vector<string> *nm = NULL);

		/**
		 * contructs a new genome with the permutation given in an int vector
		 * @param[in] p permuation
		 * @param[in] c circularity of the genom
		 * @param[in] nm pointer to a name map
		 * @param[in] norm_to normalize circular genomes to this element (default 1)
		 */
		genom(const vector<int> &p, int c, vector<string> *nm = NULL, int norm_to = 1);

		/**
		 * constructs a new genom with length cnt and distance or steps dist
		 * to the identity dist; the length of the reversals used is between
		 * start_len and end_len
		 * its important that the random number genarator is initialized (srand(time(NULL));)
		 * @param[in] cnt length of the genom
		 * @param[in] dist distance to the identity
		 * @param[in] start_len the start length of the reversals
		 * @param[in] end_len the end length of the reversals
		 * @param[in] steps (true) distance in steps or (false) reversal distance
		 * @param[in] random (true) real random (false) 'evolve' the identity
		 * @param[in] c circularity of the genom
		 * @param[in] unsig construct a permutation consisting of positive elements
		 * @param[in] hd helping memory
		 */
//		genom(int cnt, int dist, int start_len, int end_len, bool steps, bool random, char c, bool unsig, hdata &hd);

		/**
		 * returns the iterator to begin() of the chromosom
		 *@return iterator
		 */
		vector<int>::iterator begin();

		/**
		 * resets the name and the chromosom
		 */
		void clear();

		/**
		 * reverse the effect of origin.identify_g(goal)
		 * @param[in] origin the origin genom used in identify
		 */
		void de_identify_g( const genom &origin );

		/**
		 * calculate decomposition of genome into maximum increasing substrings and return in vector
		 * @param[in,out] deco: ordered vector of sets each containing a maximum increasing substring
		 */
		void decomposition(vector< set<int> > &deco);

		/**
		 * calculates decomposition of genome into strict maximum increasing substrings (consecutive
		 * 		elements differ by one) and return in vector
		 * @param[in,out] deco: ordered vector of sets each containing a maximum increasing substring
		 */
		void decomposition_strict(vector< set<int> > &deco);

		/**
		 * calculates the number of maximal increasing substrings of g
		 * @return: number of maximal increasing substrings
		 */
		unsigned int max_inc_substrings();

//		/**
//		 * calculates (grappa) the reversal distance of the calling genom an the genom g
//		 * @param g genom
//		 * @param[in] hd helping memory
//		 * @return the reversal distance
//		 */
//		int distance(const genom &g, hdata &hd) const;
//
//		/**
//		 * calculates the sum of the distances between the calling genom an the genomes in the vector
//		 * @param[in] respect vector of genomes
//		 * @param[in] hd helping memory
//		 * @return the sum of the reversal distances
//		 */
//		int distance(const vector<genom> &respect, hdata &hd) const;

		/**
		 * get the set of genes of the genomes
		 * as sorted vector to treat duplicates
		 * @return the set of genes
		 */
		vector<unsigned> genset() const;

		/**
		 * return the pointer to the namemap
		 * @return the pointer
		 */
		vector<string> *get_nmap() const;

		/**
		 * checks if the genome has duplicate elements
		 * @return true iff it has duplicates
		 */
		bool hasduplicates( ) const;

		/**
		 * returns the iterator to end() of the chromosom
		 *@return iterator
		 */
		vector<int>::iterator end();

		/**
		 * delete gene i from the genome and decrease all genes greater than i by one
		 * @param[in] i the gene
		 */
		void erase(int i);

		/**
		 * evolve : do r random reversals (/transpositions)
		 *@param[in] r number of reversals
		 *@param[in] t do every t-th step (randomly) a transposition,
		 *-t=0 never do a transposion
		 *-t=1 always do a transposition
		 *-t=2 do a transp. every 2nd step (randomly)
		 *.
		 *@param[in] l length of the reversals of transpositions
		 */
		void evolve(int r, int t, int l = 0);

		/**
		 * returns an integer pointer to the begin of the chromosom
		 * vector. that can be used as a standard C array.
		 */
		int * get_pointer() const;

		/**
		 * getter for the chromosom
		 *@return the chromosom vector
		 */
		vector<int> getChromosom() const;

		/**
		 * getter for the private member circular
		 *@return circular
		 */
		char getCircular() const;


		/**
		 * get the components of a genome
		 * - every component gets an index \f$i\f$ (\f$ i\in 1,\ldots ,c\f$ with \f$c\f$ : number of components)
		 * - every point belongs to exactly one component
		 * - in the resulting vector the indices of each component is stored
		 * - NOTE: components of length 1 get all the same component index : 0 and are stored as oriented
		 *.
		 *@param[in] g the genome
		 *@param[out] component_orientation the orientation of the components true if oriented / false else
		 *@param[out] component_cnt number of components
		 *@param[in,out] d memory
		 *@return the indice of the component for each point
		 */
//		vector<int> getComponents(const genom &g, vector<int> &component_orientation, unsigned &component_cnt) const;

		/**
		 * getter for the name of the node
		 *@return the name
		 */
		//~ string getName();

		/**
		 * get the orientation of the points
		 * the orientation of a point is:
		 * - positive (+1) if the neighboured elements of the permutaion are both positive
		 * - negative (-1) if the neighboured elements of the permutaion are both negative
		 * - 0 else
		 *.
		 *@param[in] g a genom
		 *@return the orientation of the points
		 */
		vector<int> getPointOrientations(const genom &g) const;


		/**
		 * return all possible reversals for the genome (only size interrests)
		 * the linear and circular case ist distinguished
		 *@return a vector of pairs (start, stop)
		 */
		vector< pair<int,int> > getReversals() const;

		/**
		 * return all reversals which are on the same cycle
		 * the array is sorted in that way that for all i<j:
		 * r[i]_2 < r[j]_2 and if r[i]_2 = r[j]_2 then
		 * r[i]_1 < r[j]_1.
		 * this behaviour could be changed be changing the
		 * order of the elements of the reversal, so that the
		 * ordering is first by index 1 and than index 2.
		 * (I have tried this, but it changes the results
		 * slightly. So I've decided to leave it as is for the
		 * sake of reproducibility.)
		 * @param g the second genome
		 * @return the reversals
		 */
		vector< pair<int,int> > getReversalsSameCycle(const genom &g) const;


		/**
		 * return all reversals wich act on different cycles of the same unoriented component
		 *@param[in] g a genome
		 *@return the reversals acting on different cycles of the same unoriented component
		 */
		vector< pair<int, int> > getReversalsSameUnorientedComponent(const genom &g) const;

		/**
		 * return all reversals wich act on different cycles of the same unoriented component and all
		 * reversals acting on the same cycle
		 * @param[in] g a genome
		 * @return the reversals
		 */
		vector< pair<int, int> > getReversalsSameUnorientedComponentAndSameCycle(const genom &g) const;

		/**
		 * get all reversals which dont break adjacencies of the calling genome and respect
		 * @param respect the 2nd genome
		 * @param[in] helping memory
		 */
		vector< pair<int,int> > getReversals_nonAdjacencyBreaking(const genom &respect) const;

		/**
		 * compute all possible transpositions. the returned vector consists of vectors of size 3
		 *@return the transpositions
		 */
		vector<vector<unsigned> > getTranspositions();

		/**
		 * computes the permutation p = p(q) (where p is this)
		 *@param[in] q the 2nd genome
		 *@return the 'identified' genome
		 */
		genom identify_g(const genom &q) const;

		/**
		 * computes the permutation p = p(q) (where p is this)
		 *@param[in] q the 2nd genome
		 *@return the 'identified' genome
		 */
		vector<int> identify(const genom &q) const;

		/**
		 * get the inverse permutation of a signed permutation
		 * (index i of the inverse stores the position of the element i
		 * in permutation g, if element i is negative then the element of
		 * the inverse is also negative)
		 * @param[in] sign add signs to the inverse (default: true)
		 * @return the inverse permutation
		 */
		vector<int> inverse( bool sign = true) const;

		/**
		 * check if the complete permutation consists of nagative elements only
		 * @return true iff there are only negative elements
		 */
		bool isneg() const;

		/**
		 * check if the complete permutation consists of positive elements only
		 * @return true iff there are only positive elements
		 */
		bool ispos() const;

		/**
		 * normalizes a genome:
		 * - if circular: rotate such that the element to is at position 0
		 *@param[in] to normalize to this element (default 1)
		 */
		void normalize(int to=1);

		/**
		 * return a reference to the element element at index i
		 * of the chromosom. This can be used to get the element
		 * or write it.
		 *@param[in] i the index
		 *@return the element at index i
		 */
		int& operator[](unsigned i);

		/**
		 * return the element element at index i of the chromosom.
		 * This can be used to get the element.
		 *@param[in] i the index
		 *@return the element at index i
		 */
		int operator[](unsigned i) const;

		/**
		 * checks if the chromosomes of the calling genom and the genom g are equal
		 *@param g genom
		 *@return true if equal false else
		 */
		bool operator==(const genom& g) const;

		/**
		 * checks if the chromosomes of the calling genom and the genom g are unequal
		 *@param g genom
		 *@return true if unequal false else
		 */
		bool operator!=(const genom& g) const;

		/**
		 * checks if the chromosomes are smaller (needed for saving in a map)
		 * the test is like a string < test; it's important that the genomes have
		 * equal length!!
		 *@param g genom
		 *@return true if < false else
		 */
		bool operator<(const genom& g) const;

		/**
		 * output operator for the genom (outputs the genom and then the name)
		 *@param out the stream to output in
		 *@param g the genom to output
		 *@return output
		 */
		friend ostream & operator<<(ostream &out, const genom &g);

		/**
		 * output operator for a vector of genomes
		 * it prints the Genomes and the distancematrix
		 *@param out the stream to output in
		 *@param G the genomes to output
		 *@return output
		 */
		friend ostream & operator<<(ostream &out, const vector<genom> &genomes);

		/**
		 * print an single element of the genome
		 * @param[in] idx the element at the specified position will be printed
		 * @param[in,out] out the stream to write into
		 * @param[in] inv print inv*chromosom[i] .. use for print an inverted element by setting to -1 (default 1)
		 * @param[in] plus the string which should be used for '+' (default '')
		 */
//		void print_element( unsigned idx, ostream &out, int inv = 1, string plus = "") const;

		/**
		 * append a new gene to the end of the chromosom vector
		 * @param gene the gene
		 */
		void push_back(int gene);

		/**
		 * randomise the order and orientation of the genes
		 */
		void randomise();

		/**
		 * get the interval of the given element set
		 * for a set X the result is [start,end]
		 * - i\in start..end -> c[i] \in X
		 * - i\in 0..i-start -> c[i] \not\in X
		 * - i\in end+1..n -> c[i] \not\in X
		 * .
		 * i.e., start and end included, counting starts with 0
		 * @param[in] elements the element set
		 */
		pair<int,int> interval(const set<int> &elements) const;

//		/**
//		 * Generates a 'unique' pair of number for a permutation.
//		 * Because of the limits of unsigned long the pair is not
//		 * really unique. The first element of the pair is generated
//		 * from the unsigned permutation generated by taking the
//		 * absolute values of the elements of the signed permutation
//		 * and computing the rank of the resulting permutation with
//		 * the algorithm given in Myrvold, Ruskey 2000 "Ranking and
//		 * Unranking Permutations in Linear Time".
//		 * And the second element of the pair is computed by computing
//		 * the decimal representation of the binary number resulting
//		 * out of the vector of signs.
//		 *@return the rank pair
//		 */
//		pair<unsigned, unsigned> rank(hdata &d);

		/**
		 *@todo comment
		 */
		vector<genom> reverseHamiltonian(vector<pair<int,int> > &r);

		/**
		 * reverse the genom by all the reversals given in the tupel
		 * keep track of the other index pairs if they are disturbed
		 * and recursively do the same for the remaining
		 * and return the resulting genomes
		 *@param r the reversals reverse
		 *@param used the genomes already used
		 *@param dist the reversal distance between the initial genomes (max. recursion depth)
		 *@return the genomes which are reached
		 */
		//vector<genom> reverse(vector< pair<int, int> > &r, map<genom, bool> &used, int dist);

		/**
		 * setter for the chromosom
		 *@param c the new chromosom vector
		 */
		void setChromosom(vector<int> c);

		/**
		 * setter for the circular property
		 * the normalize function is also called
		 *@param c the circularity of the genome
		 */
		void setCircular(char c);

		/**
		 * setter for the name of the node
		 *@param n the new name
		 */
		//~ void setName(string n);

		/**
		 * set the namemap
		 * @param[in] nm the new name map
		 */
		void set_nmap( vector<string> *nm );

		/**
		 * getter for the size (lenght) of the chromosom
		 *@return the size
		 */
		unsigned int size() const;
};

/**
 * eception thrown if the interval function fails
 */
class NoIntervalException : public std::runtime_error {
public:
	NoIntervalException( const genom &g, const set<int> &elements ) : std::runtime_error("NoIntervalException") {
//		cerr <<g<<endl;
//		copy(elements.begin(), elements.end(), ostream_iterator<int>(cerr, " ")); cerr <<endl;
	}
};

/**
 * eception thrown if the interval function fails
 */
class NonConsecutiveIntervalsException : public std::runtime_error {
public:
	NonConsecutiveIntervalsException( const genom &g, const set<set<int> > &t ) : std::runtime_error("NonConsecutiveIntervalsException") {
//		cerr << "error: transpose non tandem intervals"<<endl;
//		cerr <<g<<endl;
//		for(set<set<int> >::iterator it=t.begin(); it!=t.end(); it++){
//			cerr << "{"; copy(it->begin(), it->end(), ostream_iterator<int>(cerr, " "));cerr << "}"<<endl;
//		}
	}
};

///**
// * get the number of breakpoints
// * @param[in] g1 genome one
// * @param[in] g2 genomes two
// * @param[in] hd helping memory
// * @return number of breakpoints
// */
//int breakpoints(const genom &g1, const genom &g2, hdata &hd);
//
///**
// * get the number of breakpoints in the given genomes. to get correct results for
// * circular genomes, the genomes have to be normalised.
// *@param[in] genomes the genomes
// *@param[in] hd helping memory
// *@return number of breakpoints
// */
//int breakpoints(const vector<genom> &genomes, hdata &hd);

///**
// *@param[in] sign if 1 -> do not regard adjacencies with different signs
// */
//void condense (vector<genom> &genomes, char circ, int *condense_succ, int *condense_decode, int sign );

/**
 * get the cycles of the permutation p
 * - permutation: p =  p_0   p_1   p_2   p_3  ....   p_n
 *   points          0     1     2     3     4     n
 * - every cycle gets an index
 * - every point belongs to exact one cycle \f$i\f$ (\f$i \in 1,\ldots ,c\f$ with \f$c\f$ : number of cycles  )
 * - in the cycle vector are the cycle-indice of the points stored
 *.
 * @param[in] pi the gene order
 * @param[in] n the length of the genomes
 * @param[out] cycs the cycle-indices for each point-index (has to be of size n+1 !!)
 * @param[in,out] d memory
 * @return number of cycles
 */
int cycles(const vector<int> &pi, int n, vector<int> &cycs);

/**
 * de-identify a set of genomes. this is the opposite operation of identify.
 * basically deidentify renames the elements of given genomes with the given mapping.
 * @param[in] map the mapping
 * @param[in,out] genomes the genomes to rename
 */
void deidentify( const genom &map, vector<genom> &genomes );

//int distance(const genom &g1, const genom &g2, int n, hdata &hd);

/**
 * construct unsigned genome by replacing each element by 2 ordered
 * e -> 2e-1 2e; -e -> -2e -2e-1, additionally directed gene orders
 * are framed by 0 and 2n+1
 *
 * @param[in,out] g a signed gene order, replaced by the unsigned ones
 * @param[in,out] nmap the name map for the element (should have been adapted
 *   accordingly by double_nmap)
 * @param[in] directed directed: true undirected: false
 */
void double_genome( genom &g, vector<string> *nmap, bool directed );

/**
 * construct unsigned genomes by replacing each element by 2 ordered
 * e -> 2e-1 2e; -e -> -2e -2e-1, additionally directed gene orders
 * are framed by 0 and 2n+1
 *
 * @param[in,out] genomes the signed gene orders, replaced by the unsigned ones
 * @param[in,out] nmap the name map for the element (should have been adapted
 *   accordingly by double_nmap)
 * @param[in] directed directed: true undirected: false
 */
void double_genomes( vector<genom> &genomes, vector<string> *nmap,
		bool directed );

/**
 * construct the element to gene name map for doubled genomes
 * the given nmap is adapted accordingly:
 * - a gene X in the name map is replaced by +X and -X
 * - in the case of directed genomes additional elements 0 and 2n+1 are given
 *   names STA and END
 *
 * @param[in,out] nmap the name map
 * @param[in] directed directed: true undirected: false
 */
vector<string> double_nmap( const vector<string> &nmap, bool directed );

/**
 * get the elementary reversals of the given permutations (and id)
 * - permutation: p =  p_0   p_1   p_2   p_3  ....   p_n
 *   points          0     1     2     3     4     n
 * - for every element i in the permution there is an elementary reversal
 * 	which start at the point and
 *			right of i if i is positive
 *			left of i otherwise
 * 	and ends at the point
 *			left of the element (i+1) if it's positive and
 *			right of it else
 * - at every point there starts one elementary interval and one ends
 * - the endpoint of the interval starting at point i is stored in ei[i]
 * 		(note: 'starting' is choosen arbitrary : it's just one of the two intervals)
 * @param[in] pi the permutation genom
 * @param[in] n the length of the permutations
 * @param[out] ei the elementary intervals. all elements must be initialised with -1
 */
void elementary_intervals(const vector<int> &pi, int n, vector<pair<int,int> > &ei);

/**
 * construct a negative identity of length n
 * @param n the length of the genom
 * @return the negative id
 */
genom genom_nid(int n);

/**
 * construct all permutations of length n
 * @param[in] n length
 * @param[in] sig signed permutations
 * @param[out] all permutations
 */
void get_all_permutations(unsigned n, int circular, int sig, vector<genom> &all);

void identify(vector<genom> &genomes);

/**
 * apply a indel to a genome
 * the indel is given by three sets in a vector id
 * - id[0], id[2] are the framing sets of elements
 * - id[1] the elements inbetween id[0], id[2] that should be deleted
 *   or inserted between id[0], id[2]
 * - del true: deletion, false: insertion
 * @param[in] inv if true  then the elements are inserted as inverse
 */
void insertdelete(genom &g, const vector<set<int> > &id, bool del, bool inv);

/**
 * print an element. if the namemap is provided than print the name
 * @param[in] e the element
 * @param[in,out] out the stream to write in
 * @param[in] inv invert the output
 * @param[in] plus the symbol to use for '+' (default '')
 * @param[in] nmap the namemap (size = n+1)
 */
void print_element( int e, ostream &out, int inv, string plus, const vector<string> *nmap );

/**
 * reverses the interval between start and end
 * @param[in,out] g the genome
 * @param[in] start start of the interval
 * @param[in] end end of the interval
 */
void reverse(genom &g, int start, int end);

/**
 * reverses the interval given by the tupel
 * @param[in,out] g the genome
 * @param[in] r tupel(start, end)
 */
void reverse(genom &g, pair<int, int> r);

/**
 * do a reversal specified by the element set
 * @param[in,out] g the genom
 * @param[in] the element set
 */
void reverse(genom &g, const set<int> &elements);


/**
 * reverse an interval in all genomes
 * @param[in,out] genomes the genomes
 * @param[in] s start of the interval
 * @param[in] e end of the interval
 */
void reverse(vector<genom> &genomes, int s, int e);

/**
 * reverse the genom by all the reversals given in the tupel
 * keep track of the other index pairs if they are disturbed
 * and recursively do the same for the remaining
 * and return the resulting genomes
 *@param r the reversals reverse
 *@param used the genomes already used
 *@param dist the reversal distance between the initial genomes (max. recursion depth)
 *@param out ???
 */
//void reverse(vector< pair<int, int> > &r, map<genom, bool> &used, int dist, unsigned &out);

/**
 * transposition of a region. the segments \f$ \pi_i,\ldots,\pi_{j-1}\f$ and \f$ \pi_j,\ldots,\pi_{k-1}\f$
 * are exchanged. \f$k \leq n\f$, \f$i< j < k \f$ and \f$k \leq n \f$
 * @param[in,out] g the genome to modify with the transposition
 * @param[in] i begin of the region
 * @param[in] j end of the region
 * @param[in] k position where the region should be inserted
 */
void transpose(genom &g, unsigned i, unsigned j, unsigned k);

/**
 * transposition of a region.
 * @param[in,out] g the genome to modify with the transposition
 * @param[in] t and vector of the 3 needed indices
 */
void transpose(genom &g, vector<unsigned> &t);

/**
 * apply a transposition given as element sets
 */
void transpose(genom &g, const set<set<int> > &t);

void transpose(genom &g, const vector<set<int> > &t);

///**
// *
// */
//void uncondense ( vector<genom> &genomes, char circ, int *condense_succ, int *condense_decode, int orig_num_genes );

/**
 * the inverse to double_genomes applied to one genome
 *
 * @param[in] g an unsigned gene order
 * @param[in,out] nmap the name map for the (undoubled) elements
 * @param[in] directed directed: true undirected: false
 * @return the signed gene order
 */
genom undouble_genom( const genom &g, vector<string> *snmap, bool directed);

/**
 * the inverse to double_genomes
 *
 * @param[in,out] genomes the unsigned gene orders, replaced by the signed ones
 * @param[in,out] nmap the name map for the (undoubled) elements
 * @param[in] directed directed: true undirected: false
 */
void undouble_genomes( vector<genom> &genomes, vector<string> *nmap,
		bool directed );

/**
 * remove duplicates from a vector of genomes
 * @param[in,out] genomes the vector
 */
void unique(vector<genom> &genomes);

/**
 * remove duplicate genomes, also delete from names vector
 * @param[in,out] genomes the genomes
 * @param[in,out] names the names of the genomes
 * @param[in,out] tax taxomomy strings (of equal genomes are shortened to the most common part)
 * @param[in,out] cnt gives the multiplicity of each unique gene order
 * @param[in] report if true: write equalities to cerr; otherwise combine the names of equal genomes
 */
void unique( vector<genom> &genomes, vector<string> &names,
		vector<vector<string> > &tax,
		vector<unsigned> &cnt, bool report );

/**
 * make genomes unsigned, i.e. just remove all '-'
 * if linear: the permutations are framed by 0 and n+1
 * and the name map is adjusted accordingly
 *
 * @param[in,out] genomes the genomes
 * @param[in,out] nmap the element name map, is modified iff circular is false
 * @param[in] circular circularity of the genomes
 */
void unsign_genomes( vector<genom> &genomes, vector<string> *nmap, bool circular );

/**
 * get the other element of 2i 2i-1
 * @param an integer
 * @return 2i for 2i-1 and 2i-1 for 2i
 */
unsigned unsign_other( unsigned e );

/**
 * construct an edge of the permutation matching \f$ (2|\pi_i|-\nu(\pi_i), 2|\pi_{i+1}-1+\nu(\pi_{i+1})| ) \f$ where \f$ \nu(pi_i) = 1 \f$ if \f$ pi_i < 0\f$ else \f$ \nu(pi_i) = 0 \f$.
 * if \f$ |pi_i| \geq len \rightarrow (2*len+1,?)\f$ and if \f$ |\pi_{i+1}| > len \rightarrow (?, 2*len+1) \f$ is returned.
 * additionaly \f$ (?,0) \f$ is returned if \f$ \pi_{i+1} = 0 \f$.
 *@param pii \f$ \pi_i \f$
 *@param piip1 \f$ \pi_{i+1} \f$
 *@param len length of the underlying genome
 *@return the edge
 */
pair<unsigned, unsigned> matching_edge(int pii, int piip1, int len);

//struct hdata{
//		// indicate if the datastructure is already initialised
//	bool initialised;
//		// memory for identify and identify_g
//	genom g;
//	vector<int> pi_inv,
//		identified_chr;
//		// memory for rank
//	vector<unsigned> u_pi,
//		inv_u_pi;
//	/*@todo remove*/
//		// memory for distance
//	distmem_t *distmem;
//		// memory for capraras median solver (my wrapper)
//	int *genes;
//	struct genome_struct *gs[3];
//};
//
///**
// * checks if the help data structure is already initialised
// *@param[in] hd the help data structure
// *@return true if initialised; false else
// */
//bool is_initialised(hdata &hd);
//
///**
// * initialises the hdata help data structure
// *@param[in,out] d the help data structure
// *@param[in] len len of the genomes (for wich the data structure should be used)
// *@param[in] circ circularity of the genomes (for wich the data structure should be used)
// */
//void init_data(hdata &d, unsigned len, char circ);
//
///**
// * wrapper for hdata help data structure initialisation function; for usage from C
// *@param[in,out] d the help data structure
// *@param[in] len len of the genomes (for wich the data structure should be used)
// *@param[in] circ circularity of the genomes (for wich the data structure should be used)
// */
//extern "C" void init_data_c(hdata &d, unsigned len, char circ);
//
///**
// * frees the help data structure d
// *@param[in] d the help data structure
// */
//void free_data(hdata &d);
//
///**
// * wrapper for hdata help data structure free function; for usage from C
// *@param[in] d the help data structure
// */
//extern "C" void free_data_c(hdata &d);

#endif
