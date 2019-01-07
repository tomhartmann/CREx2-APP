/**
 *@file helpers.hpp some help functions
 */
#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_


#include <cctype>
#include <cmath>
#include <iostream>
#include <iterator>
#include <fstream>
#ifdef USEMPI
#include <mpi.h>
#endif//USEMPI
#include <sstream>
#include <string>
#include <stack>
#include <sys/types.h>
#include <unistd.h>
#include <string>
#include <vector>

//#define DEBUG_MARRAY

#include "genom.hpp"

using namespace std;

/**
 * comparison class for pointers to objects in a set
 */
struct DereferenceLess{
	template <typename PtrType>
	bool operator()(PtrType pT1, PtrType pT2) const{
    	return (*pT1 < *pT2);
    }
};


void error(char *, const char *);

/**
 * convert pair vector to string
 *@param[in] v the vector
 *@return string representation of v
 */
string output(vector<pair<unsigned, unsigned> > &v);

template<class Typ>
void marray_init( vector<Typ> &marray, vector<int> &lengths ){
	int size = 1;

#ifdef DEBUG_MARRAY
	cout << "initing (";
	for(unsigned i=0; i<lengths.size(); i++)
		cout << lengths[i]<<" ";
	cout << ")"<<endl;
#endif//DEBUG_MARRAY
	for(unsigned i=0; i<lengths.size(); i++){
		size *= lengths[i];
	}
	marray = vector<Typ>( size );
}

// gets the element at position[pos_0]...[pos_n] in the multidimensional array marray of dimension
// n=length.size(), where the i-th dimension has length length[i]
template<class Typ>
Typ marray_get( const vector<Typ> &marray, const vector<int> &lengths, const vector<int> &pos ){
	int p = 0,
		f = 1;

	for(unsigned i=0; i<lengths.size(); i++){
		p += ( pos[i] * f );
		f *= lengths[i];
	}
#ifdef DEBUG_MARRAY
	cout << "getting (";
	for(unsigned i=0; i<lengths.size(); i++)
		cout << lengths[i]<<" ";
	cout << ") at (";
	for(unsigned i=0; i<pos.size(); i++)
		cout << pos[i]<<" ";
	cout << ") = " << p <<endl;
#endif//DEBUG
	return marray[ p ];
}

// sets the element at position[pos_0]...[pos_n] in the multidimensional array marray of dimension
// n=length.size(), where the i-th dimension has length length[i]
template<class Typ>
void marray_set( vector<Typ> &marray, const vector<int> &lengths, const vector<int> &pos, Typ elem ){
	int p = 0,
		f = 1;

	for(unsigned i=0; i<lengths.size(); i++){
		p += ( pos[i] * f );
		f *= lengths[i];
	}
#ifdef DEBUG_MARRAY
	cout << "setting (";
	for(unsigned i=0; i<lengths.size(); i++)
		cout << lengths[i]<<" ";
	cout << ") at (";
	for(unsigned i=0; i<pos.size(); i++)
		cout << pos[i]<<" ";
	cout << ") = "<<p<<endl;
#endif//DEBUG_MARRAY
	marray[p] = elem;
}


/**
 * initialise a counter vector
 * @param[in] pos the vector
 * @param[in] size the desired length of the vector
 * @param[in] min the minimal value of the integers in the vector (default = 0)
 * @param[in] inc let the values be increasing \f$ (min, min+1,\ldots,min+size-1) \f$ otherwise
 * \f$ (min, \ldots, min ) \f$ (default: false)
 */
void counter_init( vector<int> &pos, unsigned size, unsigned min = 0, bool inc = false );

/**
 * invalidate a counter
 * @param[in,out] pos the counter to invalidate (everything is set to INT_MAX)
 */
void counter_invalidate( vector<int> &pos );

/**
 * check if a counter is valid
 * @param[in] pos the counter
 * @return true iff the counter is valid
 */
bool counter_valid( const vector<int> &pos );
/**
 * counter: given a vector of integers <max ..  compute another
 * vector of integers of the same length where all elements are also < max
 *
 * @param[in] max the maximal value of the integers in the vector is max-1
 * @param[in,out] pos the vector, if there is no next vector all elements
 *  are set to INT_MAX (or if inject > pos.size())
 * @param[in] inject do not change elements at positions < inject (default = 0)

 */
void counteradd( vector<int> &pos, int min, int max, unsigned inject = 0, bool inc = false );


/**
 * counting .. count up the elements in pos (each up to length)
 * if we get out of the range -> set all elements in pos to INT_MAX
 * @param[in] inject the inject position
 * @param[in] length the length of each element
 * @param[in, out] pos the ...
 */
void countt( unsigned inject, const vector<int> &lengths, vector<int> &pos);

///**
// * allocate an multidimensional array with length.size() dimensions
// * the length of each dimension is stored in the lengths array
// * @param[in] lengths the length of each dimension
// * @return the multidimensional array
// */
//template<class Typ>
//void marray_init( vector<Typ> &marray,  vector<int> &lengths );
//
///**
// * gets the element at position[pos_0]...[pos_n] in the multidimensional array marray of dimension
// * n=length.size(), where the i-th dimension has length length[i]
// * @param[in] marray the multidimensional array
// * @param[in] lengths the length of the dimensions
// * @param[in] pos the position to access
// * @return the element
// */
//template<class Typ>
//Typ marray_get( const vector<Typ> &marray, const vector<int> &lengths, const vector<int> &pos );
//
///**
// * sets the element at position[pos_0]...[pos_n] in the multidimensional array marray of dimension
// * n=length.size(), where the i-th dimension has length length[i]
// * @param[in] marray the multidimensional array
// * @param[in] lengths the length of the dimensions
// * @param[in] pos the position to access
// * @param[in] elem the value to set
// */
//template<class Typ>
//void marray_set( const vector<Typ> &marray, const vector<int> &lengths, const vector<int> &pos, Typ elem );

/**
 * ranks a permutation with the algorithm from
 * Myrvold, Ruskey 2000 "Ranking and Unranking Permutations
 * in Linear Time"
 *@param[in] n the length of the permutation
 *@param[in] pi the (unsigned) permutation
 *@param[in] pi_inv the inverse of pi;
 *@return the rank
 */
unsigned rank_perm(unsigned n, vector<unsigned> &pi, vector<unsigned> &pi_inv);

int string2int(string );

/**
 * transforms a given integer to a string and fills with 0s
 *@param i the integer to transform
 *@param m a maximum value (if digits(m)>digits(i) the string
 *	begins wit digits(m)-digits(i) additional 0s) m defaults to 0
 */
string int2string(int i, int m=0);

string float2string(float i, int m);

/**
 * string tokenizer
 * @param[in] str the string to split
 * @param[out] tokens the 'splits'
 * @param[in] delimiters the characters at which the string should be splitted (default ' ')
 */
void split_string(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ");

/**
 * remove leading whitespace from a string
 * @param[in,out] str the string
 */
void strip_leading_ws(string &str);

/**
 * remove trailing whitespace from a string
 * @param[in,out] str the string
 */
void strip_trailing_ws(string &str);

/**
 * return a lower case copy of the string
 * @param[in] s the string
 * @return a lowercase copy
 */
string lower( const string &s);

/**
 * return a upper case copy of the string
 * @param[in] s the string
 * @return a uppercase copy
 */
string upper( const string &s);


/**
 * binomial coefficient aka. compute n choose k
 * @param[in] n
 * @param[in] k
 * @return binomial coefficient
 */
unsigned binom(unsigned n, unsigned k);

/**
 * compute the n-th row of pascals triangle
 * @param[in] n the row number
 * @return an array containing the n-th row, the returned vector has length n+1: element for i in [1:n] i contains \binom{n}{i} and element 0 is 0
 */
//vector<unsigned> pascal_row(unsigned n);
vector<double> pascal_row(unsigned n);

/**
 * determine how many digits a number has
 *@param[in] number the number
 *@return the number of digits
 */
unsigned digits(const unsigned number);

/**
 * factorial
 *
 * ATTENTION IN THE CASE OF AN OVERFLOW 0 IS RETURNED
 *
 * @param a fact of
 * @return the factorial (or 0 when overflow)
 */
unsigned fact(unsigned a);

/**
 * binary function determining the minimum of two values
 */
template<typename T>
struct minfunc : public binary_function<T,T,T> {
	/**
	 * the functor
	 * @param[in] a a value
	 * @param[in] b another value
	 * @return the minimum of a and b
	 */
	unsigned operator() (T a, T b) {return min(a,b);}
};

/**
 * power of ints
 *@param b base
 *@param p power
 *@return power
 */
int pow(int b, int p);

/**
 * compute 2^p
 * @param p power 0
 * @return 2^p
 */
int ppow(int p);

/**
 * determine the sign of an integer
 * @param[in] a integer
 * @return +1 for positive, -1 for negative, 0 else
 */
int sign(const int &a);

/**
 * prints a interger stack
 *@param s the stack
 */
void outputStack(stack<int> &s);

/**
 * combine the elements of a bitvector with or
 *@param[in] b the bitvector
 *@return the result
 */
bool disjunction(const vector<bool> &b);

/**
 * combine the elements of a bitvector with and
 *@param[in] b the bitvector
 *@return the result
 */
bool conjunction(const vector<bool> &b);

/**
 * marks a given position in the bit_vector with a given value
 *@param[in,out] b the bitvector
 *@param[in] pos the position
 *@param[in] val the new value
 */
void mark(vector<bool> &b, unsigned pos, bool val);

/**
 * checks if the element at the given position is true
 *@param b the bitvector
 *@param pos the position
 *@return true (position is true) /false (else)
 */
bool isMarked(const vector<bool> &b, unsigned pos);

/**
 * count the number of true elements
 *@param b the bitvector
 *@return the number of true elements
 */
unsigned sum(const vector<bool> &b);

/**
 * prints the bitvector
 */
void print(const vector<bool> &b);

/**
 * initialize a random number generator
 *@param[in] init initvalue for the rng. defaults to time(NULL)+getpid()
 */
void init_rng(unsigned init=time(NULL)+getpid());

/**
 * get a random number from the generator
 *@return a 'random' number
 */
int ask_rng();

/**
 * get a float random number [0..1]
 */
float ask_rng_f();

/**
 * get a random number from a range
 * @param[in] s range start
 * @param[in] e range end
 *@return a 'random' number
 */
int ask_rng( int s, int e );

/**
 * get a random l-tuple (s_1, ..., s_l) such that
 * - \forall i\leq l: 0\leq s_i<m
 * - \forall i<l: s_i<s_{i+1}
 * - s_l - s_1 \leq d
 * .
 * the choice is made such that every tuple that fulfills the properties
 * is equally likely
 * @param[in] l the length of the tuple
 * @param[in] m 1+maximum value
 * @return a random l-tuple with increasing values smaller then m
 */
vector<unsigned> rng_inc_seq( unsigned l, unsigned m, unsigned d );

/**
 * compute the prefix sums of the probabilities given
 * @param[in] prob some probabilities
 */
void weighted_choice_init( vector<float> &prob);

/**
 * choose a random integer between 0 and prob.size()-1 given different 'probabilities'
 * of each integer
 * @param[in] prob probabilities of each integer (index)
  */
int weighted_choice( const vector<float> &prob);

/**
 * compares two scoretable entries (compares the first element)
 *@param[in] a 1st scoretable
 *@param[in] b 2nd scoretable
 *@return a[0]<b[0]
 */
bool scoreTableCmp(const vector<unsigned> &a, const vector<unsigned> &b);

/**
 * gets the current time. if mpi is used then it returns the time
 * from Wtime; else it uses clock_gettime
 *@return the time in seconds
 */
double get_time();

#endif//_HELPERS_HPP_
