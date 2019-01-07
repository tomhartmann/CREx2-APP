#ifndef REARRANGEMENTS_HPP_
#define REARRANGEMENTS_HPP_

#include <iostream>
#include <set>
#include <vector>

#include "genom.hpp"
#include "helpers.hpp"
#include "tdl.hpp"
class costfoo;

/**
 * rearrangement type codes, they are used for sorting the rearrangements in
 * unordered/alternative scenarios
 *
 * important note: placing reversals after TDRLs is important, imagine a
 * reversal which in includes all elements of a tdrl, then if the reversal is
 * applied before the tdrl the meaning of "keep elements in the first copy"
 * respectively "second copy" should also be reversed; with the given sorting
 * reversals are applied after tdrls and the problem never occurs
 */
#define EMPTYSCE 0
#define INDEL 1
#define TRA 2
#define TDRL 3
#define UORD 4
#define ORD 5
#define ALT 6
#define REV 7
#define RTRA 8

#define TRANM "T"
#define REVNM "I"
#define RTRANM "iT"
#define TDRLNM "TDRL"

#define ALTMIN 1
#define ALTMAX 2
#define ALTALL 3

using namespace std;

/**
 * comparison class for storing derived objects in a set, i.e. pointers to the base class
 * the derived classes have to implement 'bool operator(const base *) const'
 */
struct HDereferenceLess{
	template <typename PtrType>
	bool operator()(PtrType pT1, PtrType pT2) const{
    	return *pT1 < pT2;
    }
};

class emptyscen;
class indel;
class rev;
class transp;
class revtransp;
class tdrl;
class unordered;
class ordered;
class alternative;

/**
 * @param[in] n length of the genomes
 * @param[in] x probability of a breakpoint to be fragile
 * @param[in] y probability of choosing a fragile breakpoint
 * @param[out] bp the fragile and non fragile breakpoints (are initialized on
 * 		the first call, and reinitialized if size != n
 */
unsigned random_breakpoint( unsigned n, float x, float y, vector<vector<int> > &bp);

/**
 * class for rearrangement .. base class for the specific rearrangements
 */
class rrrmt{
	public:

	/**
	 * destructor
	 */
	virtual ~rrrmt();

	/**
	 * append a copy of the rearrangement(s) to the vector
	 * @param[in,out] v the vector to append
	 */
	virtual void append( vector<vector<rrrmt*> > &v ) const;

	/**
	 * apply a rearrangement on a given genome
	 * @param[in,out] g the genome
	 */
	virtual void apply( genom &g )const = 0;

	/**
	 * apply a rearrangement scenario on a given genome
	 * and store all possible intermediate gener oders in a set
	 *
	 * default behaviour apply rearrangement
	 *
	 * @param[in] g the genome
	 * @param[out] gs the genome set
	 */
	virtual void apply( genom &g, set<genom> &gs )const;

	/**
	 * return a pointer to a copy of a rearrangement
	 * this is a virtual copy constructor
	 * @return the pointer
	 */
	virtual rrrmt* clone() const = 0;

	/**
	 * virtual constructor
	 * @return an empty instance of an object, the default is emptyscen
	 */
	virtual rrrmt* create() const;

	/**
	 * deapply a rearrangement
	 * @param[in,out] g a genom
	 */
	virtual void deapply(genom &g) const = 0;

	/**
	 * check if two rearrangement scenarios are equal. derived
	 * classes must implement double dispatch. default implemententations
	 * return false. so, the foos needs to be overwritten in derived classes
	 * if it might potentially be true
	 *
	 * @param[in] r another rearrangement
	 */
	virtual bool equal( const rrrmt *r ) const  = 0;
	virtual bool equal( const emptyscen *r ) const;
	virtual bool equal( const indel *r ) const;
	virtual bool equal( const rev *r ) const;
	virtual bool equal( const transp *r) const;
	virtual bool equal( const revtransp *r) const;
	virtual bool equal( const tdrl *r) const;
	virtual bool equal( const unordered *r) const;
	virtual bool equal( const ordered *r) const;
	virtual bool equal( const alternative *r) const;

	/**
	 * get the set of all atomic rearrangements in the scenario
	 * @param[out] out the rearrangements should be stored here
	 */
	virtual void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL ) = 0;
	virtual void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL ) = 0;

	/**
	 * return minimum rrrmt scenario for a given set of costs
	 * @param[out] scen: rrrmt that contains parsimonious scenario
	 * @param[in] cst: costs of rrrmts
	 **/
	virtual void getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	virtual void elements( vector<set<int> > &elms, bool ab ) const = 0;

	/**
	 * check if rearrangement has composited subrearrangements
	 * @return true (overwrite for composited rearrangements)
	 */
	virtual bool hascomposite() const;

	/**
	 * compute the intersection, this function has to be overwritten by derived classes
	 * and implement a double dispatch
	 * @param r pointer to another rearrangement
	 * @param prefix use prefix of ordered scenarios, otherwise suffix
	 * @return pointer to the intersection
	 */
	virtual rrrmt* intersect( const rrrmt *r, bool prefix=false ) const = 0;
	/**
	 *
	 */
	virtual rrrmt* intersect( const emptyscen *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const indel *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const rev *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const transp *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const revtransp *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const tdrl *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const unordered *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const ordered *rr, bool prefix=false )const;
	/**
	 *
	 */
	virtual rrrmt* intersect( const alternative *rr, bool prefix=false )const;

	/**
	 * checks if a rearrangement is complete
	 * the default behauviour is true (usefull for atomic operations)
	 * @return true
	 */
	virtual bool iscomplete() const;

	/**
	 * return false per default
	 * overwrite for empty rearrangement
	 * @return true iff the rearrangement is empty
	 */
	virtual bool isempty() const;

	/**
	 * get the number of atomic events in the scenario
	 * rrrmt::length() returns 1, overwrite for complex scenarios
	 * @param[in] mode take max (ALTMAX) / min (ALTMIN) length of alternatives
	 * @return the number of events
	 */
	virtual unsigned length( int mode ) const;

	/**
	 * make alternative scenarios, the default behavior is
	 * to return a clone() .. overwrite if wished
	 * @return alternative scenarios
	 */
	virtual rrrmt *mkealt();

	/**
	 * get the number of alternative scenarios
	 * @param return number of alternative scenarios
	 */
	virtual unsigned nralt() const;

	/**
	 * get the minimum cost of all alternative rrrmts
	 * @param return minimum cost of alternative scenarios
	 */
	virtual float getmincost(const costfoo *cst) const;

	/**
	 * make a scenario parsimonious .. remove alternatives
	 * implemented behavior is .. return empty scenario
	 * this is overwritten in alternative
	 *
	 * @todo weighting scheme, currently only size
	 * @return parsimonious copy of the input scenario
	 */
//	virtual rrrmt* mkepar() const;

	/**
	 * comparison operator for rearrangements. so they can be stored in sets
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else directly compare the class data
	 * this must be overwritten in all derived classes
	 * @param[in] rrrmtp the second
	 */
	virtual bool operator<(const rrrmt *rrrmtp) const = 0;

	/**
	 * output operator, dispatches to the output function of derived classes
	 */
	friend ostream &operator<<(ostream &os, const rrrmt &t);

	/**
	 * call the output function with cout
	 * @param[in] l indentation level
	 * @param[in] d detail level
	 * @param[in] pquot character to use to indicate the sign of positive elements
	 */
	virtual void output(unsigned l=0, unsigned d=0, string pquot="");

	/**
	 * output foo
	 * @param[in,out] out the stream to write into
	 * @param[in] l intendation level
	 * @param[in] d level of detail
	 * - 0 print everything
	 * - 1 only print shortest alternative
	 * - 2 only print a very short description (one letter for atomic rearrangements and some parantheses for complex rearrangements)
	 * .
	 * @param[in] pquot string to put for parantheses (eg. for quoting in dot output)
	 * @param[in] ridx rearrangement index, if given the rearrangements index is printed
	 * 	instead of the rearrangement
	 */
	virtual ostream &output( ostream &out, unsigned l=0, unsigned d=0, string pquot="") const = 0;

	/**
	 * output foo, rearrangement index is updated dynamically
	 *
	 * default behaviour is to search the index for the given rearrangement
	 * if not included determine name and insert
	 * then output the name
	 *
	 * should be overwritten for the compex types
	 *
	 * @param[in,out] out the stream to write into
	 * @param[in,out] ridx the rearrangement index
	 */
	virtual ostream &output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx );

	/**
	 * replace an element i which is part of a rrrmt by the elements
	 * stored in the mapping at position i, this function is intended
	 * to be used for renaming rearrangements identified for a quotient
	 * permutation to the real permutation
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	virtual void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx) = 0;

	/**
	 * get a very short description of a rrrmt. usually only a letter
	 * or some brackets for combined scenarios
	 * @param[in,out] d the string where to append the description
	 */
//	virtual void shortdesc( string &d ) const = 0;

	/**
	 * return a pointer to a simplified version of the rearrangement
	 * per default: a clone. overwrite for more complex behaviour
	 * @return simplified rearrangement
	 */
	virtual rrrmt *simplify();

	/**
	 * get the rearrangement type
	 * @return the code for the rearrangement type
	 */
	virtual int type() const = 0;

	/**
	 * @return the name of the rearrangement as string
	 * @param[in] d level of detail
	 */
	virtual string typestrg(unsigned d=0) const = 0;
};

/**
 * class for an emptyscen scenario.
 * @todo implement emtyscen as singleton, problem: should be used as other rearrangements .. also the delete -> but delete should only be called once for a singleton
 */
class emptyscen : public rrrmt{
	public:

	/**
	 * apply nothing
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	/**
	 * return a pointer to a new emptyscen
	 * @return the pointer
	 */
	emptyscen* clone() const;

	/**
	 * deapply nothing
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for emptyscen
	 * @param[in] r another emptyscen
	 * @return true (always)
	 */
	virtual bool equal( const emptyscen *r ) const;

	void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL );
	/**
	 * intersect 'something' with the empty scenario
	 * @param[in] r another rearrangement
	 * @return a new empty scenario
	 */
	rrrmt *intersect(const rrrmt *r, bool prefix=false) const;

	/**
	 * @return true
	 */
	bool isempty() const;

	/**
	 * @param[in] mode take max (ALTMAX) / min (ALTMAX) length of alternatives; does not matter here
	 * @return 0
	 */
	unsigned length( int mode ) const;

	/**
	 * comparison operator for tdrls
	 * just return the result of the type index comparison,
	 * therefore there will be no more than one empty scenario in a set
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output( ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;


	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){}

//	void shortdesc( string &d ) const;

	/**
	 * @return EMPTYSCE
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

class indel : public rrrmt{
private:
	vector<string> *_nmap;	// pointer to a namemap
	bool _del; 		// true: deletion, false: insertion
	bool _inv; 		// if the inserted/deleted eleemnt is inverse
	vector<set<int> > _d;	// a vector of length 3 (left neighbors, elements, right neighbors); all unsigned
public:
	/**
	 * construct insertion or deletion of element e in g
	 * @param[in] e the element in g
	 * @param[in] d true: deletion, false: insertion
	 * @param[in] g the genome
	 */
	indel( int i, bool d, const genom &g );

	/**
	 * destructor
	 */
	virtual ~indel();

	/**
	 * apply a indel on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	/**
	 * return a pointer to a new indel
	 * @return the pointer
	 */
	indel* clone() const;

	/**
	 * deapply a indel
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for indels
	 * @param[in] r another indel
	 * @return true iff delete elements and context are equal
	 */
	virtual bool equal( const indel *r ) const;

	void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL );
	/**
	 * @todo documentation
	 */
	rrrmt *intersect(const rrrmt *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const indel *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const unordered *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const ordered *rr, bool prefix=false) const;

	/**
	 * comparison operator for indel
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else directly compare the class data
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output( ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);

	/**
	 * @return DEL
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

/**
 * class for a reversal
 */
class rev : public rrrmt{
	private:
		vector<string> *_nmap;	// pointer to a namemap
		set<int> r;
	public:

	/**
	 * construct a reversal of the range i..j
	 * @param[in] i the start index
	 * @param[in] j the end index
	 * @param[in] g the genomes for which the indices are computed
	 */
	rev(int i, int j, const genom &g);

	/**
	 * construct a set given the set of elements to reverse
	 * @param[in] rset the set
	 * @param[in] nm the namemap to set
	 */
	rev(const set<int> &rset, vector<string> *nm);

	/**
	 * contruct a random reversal applied to a given genome
	 * @param[in,out] g the genome
	 * @param[in] ml maximal number of genes involved in the reversal,
	 * 	not restricted if omitted
	 */
	rev(const genom &g, unsigned ml=std::numeric_limits< unsigned >::max());

	/**
	 * destructor
	 */
	virtual ~rev();

	/**
	 * apply a reversal on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	/**
	 * return a pointer to a new reversal
	 * @return the pointer
	 */
	rev* clone() const;

	/**
	 * deapply a reversal
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for inversions
	 * @param[in] r another inversion
	 * @return true iff the sets of inverted elements are equal
	 */
	virtual bool equal( const rev *r ) const;

	void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL );
	/**
	 * @todo documentation
	 */
	rrrmt *intersect(const rrrmt *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const rev *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const unordered *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const ordered *r, bool prefix=false) const;
	/**
	 *
	 */
//	rrrmt *intersect(const alternative &r) const;

	/**
	 * comparison operator for reversals
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else directly compare the class data
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output( ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);



//	void shortdesc( string &d ) const;

	/**
	 * @return REV
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

/**
 * class for a transposition
 */
class transp : public rrrmt{
	private:
		vector<string> *_nmap;	// pointer to a namemap
		set<set<int> > t;
	public:
	/**
	 * construct a transposition, i.e. move the range i..j-1 to position before
	 * position k. in other word the intervals i..j-1 and j..k-1 are exchanged
	 *
	 * @param[in] i i
	 * @param[in] j j
	 * @param[in] k k
	 * @param[in] g the genome for which the indicees are computed
	 */
	transp(int i, int j, int k, const genom &g);

	/**
	 * contruct a random transposition applied to a given genome
	 * @param[in,out] g the genome
	 * @param[in] ml maximal number of genes involved in the transposition,
	 * 	not restricted if omitted
	 */
	transp( const genom &g, unsigned ml = std::numeric_limits< unsigned >::max() );

	/**
	 * contruct a transposition from 2! sets
	 * @param[in] tra two sets specifying the transposition
	 * @param[in] nm the namemap to set
	 */
	transp( const vector<set<int> > &tra, vector<string> *nm);

	/**
	 * destructor
	 */
	virtual ~transp();

	/**
	 * apply a transposition on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	/**
	 * return a pointer to a new transposition
	 * @return the pointer
	 */
	transp* clone() const;

	/**
	 * deapply a transposition
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for transpositions
	 * @param[in] r another transposition
	 * @return true iff the transposed sets are equal
	 */
	virtual bool equal( const transp *r ) const;

	void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL );
	/**
	 * @todo documentation
	 */
	rrrmt *intersect(const rrrmt *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const transp *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const unordered *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const ordered *rr, bool prefix=false) const;

	/**
	 *
	 */
//	rrrmt *intersect(const alternative &r) const;
	/**
	 * construct atrenatives for a transposition, i.e.
	 * - 3 reversals
	 * - 1 reverse transposition + 1 reversal
	 * @return alternative scenario
	 */
	rrrmt *mkealt();

	/**
	 * comparison operator for transpositions
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else directly compare the class data
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output(ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;


	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);

//	void shortdesc( string &d ) const;

	/**
	 * @return TRA
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

/**
 * a reverse transposition
 */
class revtransp : public rrrmt{
	private:
		vector<string> *_nmap;	// pointer to a namemap
		vector<set<int> > rt;/*! 0: the reverse transposed elements .. 1 the transposed elements */
	public:

	/**
	 * construct a reverse transposition, i.e. inverse transpose the range i..j before position k
	 * @param[in] i start of the inverse transposed range
	 * @param[in] j end of the inverse transposed range
	 * @param[in] k position of insert
	 * @param[in] g genom
	 */
	revtransp(int i, int j, int k, const genom &g);

	/**
	 * construct a reverse transposition, i.e. transpose the ranges i..j and k..l,
	 * and reverse the range i..j
	 */
	revtransp(int i, int j, int k, int l, const genom &g);

	/**
	 * contruct a random reverse transposition applied to a given genome
	 * @param[in,out] g the genome
	 * @param[in] ml maximal number of genes involved in the transposition,
	 * 	not restricted if omitted
	 */
	revtransp( const genom &g, unsigned ml=std::numeric_limits< unsigned >::max() );

	/**
	 * construct a reverse transposition given the two sets
	 * @param[in] rtra the reverse transposed elements
	 * @param[in] tra the transposed elements
	 * @param[in] nm the namemap to set
	 */
	revtransp( const set<int> &rtra, const set<int> &tra, vector<string> *nm );

	/**
	 * destructor
	 */
	virtual ~revtransp();

	/**
	 * apply a reverse transposition on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	/**
	 * return a pointer to a new reverse transposition
	 * @return the pointer
	 */
	revtransp* clone() const;

	/**
	 * deapply a reverse transposition
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for inverse transpositions
	 * @param[in] r another inverse transpositions
	 * @return true iff the transposed sets are equal and the same is inverted
	 */
	virtual bool equal( const revtransp *r ) const;

	void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL );
	/**
	 * @todo documentation
	 */
	rrrmt *intersect(const rrrmt *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const revtransp *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const unordered *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const ordered *rr, bool prefix=false) const;
	/**
	 *
	 */
//	rrrmt *intersect(const alternative &r) const;

	/**
	 * get alternatives for a reverse transposition AB -> -BA (resp. AB -> B-A), i.e.
	 * - transposition + reversal AB -> BA -> -BA (resp. AB -> BA -> B-A differs only in the reversal)
	 * - reversal + reversal AB -> -B-A -> -BA (resp. AB -> -B-A -> B-A s.a.)
	 * @return alternative scenarios
	 */
	rrrmt *mkealt();

	/**
	 * comparison operator for reverse transpositions
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else directly compare the class data
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output(ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
     * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);

//	void shortdesc( string &d ) const;

	/**
	 * @return RTRA
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

/**
 * a single tdrl
 */
class tdrl : public rrrmt{
	private:
	vector<string> *_nmap;	// pointer to a namemap
	/**
	 * 0 the elements to keep in the first copy,
	 * 1 the elements to keep in the second copy
	 */
	vector<set<int> > _t;

	/**
	 * some information of the original order of the elements necessary for
	 * deapplying! imagine each 0/1-block gets an index, then the vector
	 * stores the abs(element) index mapping
	 */
	vector<int> _orgord;

	/**
	 *
	 */
	bool realtdl;

	/**
	 * the init function that is called by all constructors
	 * @param[in] tdl bool vector specifying for each element of the genome if it is in the 1st (true) or 2nd copy (false)
	 * @param[in] g the genome
	 */
	void _init( const vector<bool> &tdl, const genom &g );

	public:
	/**
	 * construct a tdrl
	 * @param[in] tdl bool vector specifying for each element of the genome if it is in the 1st (true) or 2nd copy (false)
	 * @param[in] g the genome
	 */
	tdrl(const vector<bool> &tdl, const genom &g);


	/**
	 * construct a set given the set of elements to rearrange
	 * @param[in] set: the set of rrrmt elements; set[0]=first copy; set[i]=second copy
	 * @param[in] nm: the namemap to set
	 * @param[in] rtdrl: a real TDRL ?
	 */
	tdrl(const vector<set<int> > &set, vector<string> *nm, bool rtdrl);

	/**
	 * contruct a random tdrl applied to a given genome
	 * @param[in,out] g the genome
	 * @param[in] randrange
	 * - true: generate a random tdrl in a randomly choosen interval
	 * - false: choose the copy index randomly for all elements of the perm
	 * .
	 * @param[in] ml maximum number of elements in the randomly choosen interval
	 */
	tdrl( const genom &g, bool randrange=false, unsigned ml=std::numeric_limits< unsigned >::max() );

	/**
	 * destructor
	 */
	virtual ~tdrl();

	/**
	 * apply a tdrl on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	/**
	 * return a pointer to a new tdrl
	 * @return the pointer
	 */
	tdrl* clone() const;

	/**
	 * deapply a tdrl
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for tdrls
	 * @param[in] r another tdrl
	 * @return true iff the set of elements kept in the 1st and 2nd are equal
	 */
	virtual bool equal( const tdrl *r ) const;

	void getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &out, int mode=ALTALL );
	/**
	 * @todo documentation
	 */
	rrrmt *intersect(const rrrmt *r, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const tdrl *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const unordered *rr, bool prefix=false) const;
	/**
	 *
	 */
	rrrmt *intersect(const ordered *rr, bool prefix=false) const;
	/**
	 *
	 */
//	rrrmt *intersect(const alternative &r) const;

	/**
	 * comparison operator for tdrls
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else directly compare the class data
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output(ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);

//	void shortdesc( string &d ) const;

	/*
	 * set orgord vector
	 */
	void set_orgord(vector<int > &vec);

	/**
	 * if the tdrl is a transposition .. return the transposition
	 * otherwise a clone
	 * @return a simplified rearrangement (transposition) or a clone
	 */
	rrrmt *simplify();

	/**
	 * @return TDRL
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

/**
 * a unordered scenario, the order of the rearrangements in the scenario is unimportant
 * the difference between unordered and alternative (both are sets) is: that all rearrangements
 * in unordered have to be applied, whereas in alternative only one
 */
class unordered : public rrrmt{
friend class ordered;
friend class alternative;
protected:
	set<rrrmt*, HDereferenceLess> _sce; 	/*! the (unordered) scenario */
	bool _complete;
public:
	/**
	 * construct a ordered scenario
	 */
	unordered();

	/**
	 * copy constructor
	 * @param[in] c the original
	 */
	unordered( const unordered &c );

	/**
	 * construct a unordered scenario from a vector, the order will be lost
	 * @param[in] c the vector
	 */
	unordered( const vector<rrrmt*> &c );

	/**
	 * destructor
	 */
	virtual ~unordered();

	/**
	 * @see rrrmt::append
	 * @param[in] v the vector to append to
	 */
	void append(  vector<vector<rrrmt*> > &v ) const;

	/**
	 * apply a unordered scenario on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;
	void apply( genom &g, set<genom> &gs )const;
	/**
	 * return a pointer to a new ordered scenario
	 * @return the pointer
	 */
	virtual unordered* clone() const;

	/**
	 * @see rrrmt::create
	 */
	virtual unordered* create() const;

	/**
	 * deapply an unordered scenario,
	 * i.e. deapply all rearrangements in the unordered scenario (in some order)
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * deapply rearrangements in the scenario (handling of alternative
	 * scenarios differs from the deapply method)
	 * @see rrrmt::deapplyall
	 */
//	vector<genom> deapplyall( const genom &g ) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	virtual void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for unordered scenarios
	 * @param[in] r another unordered scenario
	 * @return true iff each subscenario is equal
	 */
	virtual bool equal( const unordered *r ) const;

	/**
	 * get the set of all atomic rearrangements in the scenario
	 * @param[out] out the rearrangements should be stored here
	 */
	void getrrrmt( set<rrrmt*,HDereferenceLess> &rset, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &rset, int mode=ALTALL );

	/**
	 * return minimum rrrmt scenario for a given set of costs
	 * @param[out] scen: rrrmt that contains parsimonious scen
	 * @param[in] cst: costs of rrrmts
	 **/
	virtual void getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const;

	/**
	 * check if the unordered scenario has a composite (ordered, unordered ...) child
	 * @return true iff there is has at least one composite child
	 */
	bool hascomposite() const;

	/**
	 * add a rearrangement to the end
	 * @param[in] r the rrrmt
	 */
	bool insert(rrrmt *r);

	/**
	 * @todo documentation
	 */
	virtual rrrmt *intersect( const rrrmt *r, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const indel *r, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const rev *r, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const transp *r, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const revtransp *r, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const tdrl *r, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const ordered *r, bool prefix=false ) const;

	/**
	 * computes the intersection of 2 unordered scenarios. from the following
	 * three sets of rearrangements one with maximum length is chosen:
	 * 1. this \f$\cap\f$ every sub-scenario of rr
	 * 2. rr \f$\cap\f$ every subscenarion of this
	 * 3. the unordered scenario consisting of the results of the intersection
	 *    of all-vs-all sub-scenarios of rr and this
	 * @param[in] r another unordered scenario
	 * @param[in] prefix use prefix of ordered sub-scenarios
	 * @return the intersection
	 */
	virtual rrrmt *intersect( const unordered *r, bool prefix=false ) const;

	/**
	 *
	 */
	virtual rrrmt *intersect( const alternative *r, bool prefix=false ) const;

	/**
	 * @return true iff all rearrangements in the scenario are complete
	 *  and the scenario is not explicitely marked as incomplete
	 */
	bool iscomplete() const;

	/**
	 * @return return true iff the scenario is empty
	 */
	virtual bool isempty() const;

	/**
	 * get the number of atomic events in the scenario
	 * @param[in] mode take max (ALTMAX) / min (ALTMAX) length of alternatives
	 * @return the number of events
	 */
	virtual unsigned length( int mode ) const;

	/**
	 * recursively make alternatives
	 */
	virtual rrrmt *mkealt();

	/**
	 * get the number of alternative scenarios
	 * for an unordered scenario this is the product of the numbers from the children
	 *
	 * @param return number of alternative scenarios
	 */
	virtual unsigned nralt() const;

	/**
	 * get the minimum cost of all alternative rrrmts
	 * @param return minimum cost of alternative scenarios
	 */
	virtual float getmincost(const costfoo *cst) const;

	/**
	 * make parsimonious
	 * @see rrrmt::mkepar
	 */
//	virtual rrrmt *mkepar() const;

	/**
	 * comparison operator for unordered
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else lexicographically compare the sets (ordered)
	 * @param[in] rrrmtp the second
	 * @return true iff <
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	virtual ostream &output(ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx );

	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);

//	void shortdesc( string &d ) const;

	/**
	 * if the ordered scenario consists of only one rearrangements return
	 * this rearrangement, otherwise a clone
	 * @return simplified rearrangement / a clone
	 */
	rrrmt *simplify();

	/**
	 * @return COMB
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};


/**
 * a ordered scenario, the order of the rearrangements in the scenario is important
 */
class ordered : public rrrmt{
friend class alternative;
protected:
	vector<rrrmt*> sce; 	/*! the (ordered) scenario */
	bool _complete;		/*! is the scenario complete */
public:
	/**
	 * construct a ordered scenario
	 */
	ordered();

	/**
	 * copy constructor
	 * @param[in] c the original
	 */
	ordered( const ordered &c );

	/**
	 * construct a ordered scenario given a vector of rearrangements
	 * @param[in] c vector of rearrangements
	 */
	ordered( const vector<rrrmt*> &c, bool comp );

	/**
	 * construct a random scenario of k atomic operations (rev, tra, revtra,
	 * tdrl) for a genome of length n given probabilities for each operation
	 *
	 * @param[in] k the number of operations
	 * @param[in] prob the probabilities/frequencies (reversals,transpositions,
	 * reversetransp,tdrl)
	 * @param[in,out] g the genome on which the random rearrangements should be
	 * applied, this will be changed during the construction
	 * @param[in] randrange (applies only to tdrls)
	 * - true: generate a random tdrl in a randomly choosen interval
	 * - false: choose the copy index randomly for all elements of the perm
	 * @param[in] ml maximum number of elements affected by EACH of the random
	 * rearrangements
	 */
	ordered( int k, vector<float> prob, genom &g, vector<genom> &trace,
			bool randrange=false, unsigned ml=std::numeric_limits< unsigned >::max() );

	/**
	 * destructor
	 */
	virtual ~ordered();

	/**
	 * @see rrrmt::append
	 * here the ordered is appended as a vector
	 * @param[in] v the vector to append to
	 */
	void append( vector<vector<rrrmt*> > &v ) const;

	/**
	 * apply a ordered scenario on a given genome
	 * @param[in,out] g the genome
	 */
	void apply( genom &g )const;

	void apply( genom &g, set<genom> &gs )const;

	/**
	 * return a pointer to a new ordered scenario
	 * @return the pointer
	 */
	ordered* clone() const;

	/**
	 * @see rrrmt::create
	 */
	ordered* create() const;

	/**
	 * deapply a the ordered scenario on a given genome,
	 * i.e. deapply all rearrangements in the scenario in the stored order
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * deapply rearrangements in the scenario (handling of alternative
	 * scenarios differs from the deapply method)
	 * @see rrrmt::deapplyall
	 */
//	vector<genom> deapplyall( const genom &g ) const;

	/**
	 * get the elements which are rearranged
	 * @param[out] elms the elements
	 * @param[out] abs take the absolut values
	 */
	void elements( vector<set<int> > &elms, bool ab ) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * test for equality for ordered scenarios
	 * @param[in] r another ordered scenario
	 * @return true iff all subscenarios are equal
	 */
	virtual bool equal( const ordered *r ) const;

	/**
	 * get the set of all atomic rearrangements in the scenario
	 * @param[out] out the rearrangements should be stored here
	 */
	void getrrrmt( set<rrrmt*,HDereferenceLess> &rset, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &rset, int mode=ALTALL );

	/**
	 * return minimum rrrmt scenario for a given set of costs
	 * @param[out] scen: rrrmt that contains parsimonious scen
	 * @param[in] cst: costs of rrrmts
	 **/
	virtual void getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const;

	/**
	 * check if the ordered scenario has a composite (ordered, unordered ...) child
	 * @return true iff there is has at least one composite child
	 */
	bool hascomposite() const;

	/**
	 * add a rearrangement
	 * @param[in] r the rrrmt
	 */
//	bool insert(rrrmt *r);

	/**
	 * @todo documentation
	 */
	virtual rrrmt *intersect( const rrrmt *rr, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const indel *rr, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const rev *rr, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const transp *rr, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const revtransp *rr, bool prefix=false ) const;
	/**
	 *
	 */
	virtual rrrmt *intersect( const tdrl *rr, bool prefix=false ) const;

	/**
	 * get the intersection of 2 ordered scenarios as the longest of the
	 * following three options:
	 * 1. the last element of this \f$\cap\f$ rr
	 * 2. the last element of rr \f$\cap\f$ this
	 * 3. the longest ordered scenario consisting of \f$this_i \cap r_i\f$
	 *    where \f$|this_i \cap r_i| = |this_i| = |r_i| \f$ holds for all
	 *    but the last element
	 *
	 * if prefix is true then the prefix instead of the suffix is used
	 *
	 * @param[in] rr another ordered scenario
	 * @param[in] prefix take prefix
	 * @return the intersection
	 */
	virtual rrrmt *intersect( const ordered *rr, bool prefix=false ) const;

	/**
	 * compute intersection of an ordered and an unordered scenario as the
	 * longest of the following three options
	 * 1. 1st of this \f$\cap\f$ complete rr
	 * 2. complete this \f$ \cap \f$ each sub-scenario of rr
	 * 3. the ordered scenario consisting of the maximum intersections
	 *    of \f$ argmax_j this_i \cap rr_j\f$, such that all of these
	 *    intersections have the same length as \f$this_i\f$ (the last
	 *    is allowed to be incomplete)
	 * if prefix is false the suffix of the ordered scenario is used
	 *
	 *@param[in] rr an unordered scenario
	 *@param[in] prefix use prefix / suffix
	 *@return the intersection
	 */
	virtual rrrmt *intersect( const unordered *rr, bool prefix=false ) const;

	/**
	 *
	 */
	virtual rrrmt *intersect( const alternative *rr, bool prefix=false ) const;

	/**
	 * @return true iff all rearrangements in the scenario are complete
	 * and the scenario is not explicitely marked as incomplete
	 */
	bool iscomplete() const;

	/**
	 * @return true iff the scenario is empty
	 */
	bool isempty() const;

	/**
	 * get the number of atomic events in the scenario
	 * @param[in] mode take max (ALTMAX) / min (ALTMAX) length of alternatives
	 * @return the number of events
	 */
	unsigned length( int mode ) const;

	/**
	 * recursively make alternatives
	 */
	rrrmt *mkealt();

	/**
	 * make parsimonious
	 * @see rrrmt::mkepar
	 */
//	rrrmt *mkepar() const;

	/**
	 * get the number of alternative scenarios
	 * for an ordered scenario this is the product of the numbers from the children
	 * @param return number of alternative scenarios
	 */
	virtual unsigned nralt() const;

	/**
	 * get the minimum cost of all alternative rrrmts
	 * @param return minimum cost of alternative scenarios
	 */
	virtual float getmincost(const costfoo *cst) const;

	/**
	 * comparison operator for ordered
	 * - if the other rearrangement has a different type, then return the result of the type index comparison
	 * - else lexicographically compare the vectors
	 * @param[in] rrrmtp the second
	 */
	bool operator<(const rrrmt *rrrmtp) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output(ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx );

	/**
	 * add a rearrangement to the end
	 * @param[in] r the rrrmt
	 */
	void push_back(rrrmt *r);

	/**
	 * @see rrrmt::rename
	 * @param[in] mapping the mapping
	 * @param[in] nmap new name mapping
	 * @param[in] mx maximum element in the new genome
	 */
	void rename( const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx);

//	void shortdesc( string &d ) const;

	/**
	 * if the ordered scenario consists of only one rearrangements return
	 * this rearrangement, otherwise a clone
	 * @return simplified rearrangement / a clone
	 */
	rrrmt *simplify();

	/**
	 * @return UORD
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};

/**
 * alternative scenarios this is mainly like an unordered scenario.
 * the deapply, and type functions are overwritten
 */
class alternative : public unordered{
	private:
	public:

	/**
	 * construct a alterative scenario
	 */
	alternative() : unordered(){};

	/**
	 * copy constructor
	 * @param[in] c the original
	 */
	alternative( const alternative &c ): unordered(c){};

	/**
	 * destructor
	 */
	virtual ~alternative(){};

	/**
	 * append the smallest alternative (smallest complete or incomplete just the smallest)
	 * @param[in,out] v the vector to append
	 */
	void append( vector<vector<rrrmt*> > &v) const;

	/**
	 * apply alternative scenarios,
	 * i.e. deapply one (the "random" first) of the rearrangements
	 * @param[in,out] g a genom
	 */
	void apply( genom &g ) const;
	void apply( genom &g, set<genom> &gs )const;
	/**
	 * return a pointer to a new alternative scenario
	 * @return the pointer
	 */
	alternative* clone() const;

	/**
	 * @see rrrmt::create
	 */
	alternative* create() const;

	/**
	 * deapply alternative scenarios,
	 * i.e. deapply one (the "random" first) of the rearrangements
	 * @param[in,out] g a genom
	 */
	void deapply(genom &g) const;

	/**
	 * deapply the rearrangements, the difference to deapply is
	 * that for alternative scenarios all alternatives are
	 * deapplied; this function is implemented in rrrmt to
	 * call the deapply function on all elements in dg;
	 * and overwritten in alternative
	 * @param[in] g a genom
	 * @return the results od deapplying all rearrangements
	 */
//	vector<genom> deapplyall( const genom &g ) const;

//	/**
//	 * get the elements which are rearranged
//	 * @param[out] elms the elements
//	 * @param[out] abs take the absolut values
//	 */
//	void elements( vector<set<int> > &elms, bool ab ) const;

	set<rrrmt*, HDereferenceLess> get_alternatives() const;

	/**
	 * get one alternative at random
	 * if no alternatives -> emptyscen
	 * !return a random alternative
	 */
	rrrmt * get_random_alternative() const;

	/**
	 * get the set of all atomic rearrangements in the scenario
	 * @param[out] out the rearrangements should be stored here
	 */
	void getrrrmt( set<rrrmt*,HDereferenceLess> &rset, int mode=ALTALL );
	void getrrrmt( vector<rrrmt*> &rset, int mode=ALTALL );

	/**
	 * return minimum rrrmt scenario for a given set of costs
	 * @param[out] scen: rrrmt that contains parsimonious scen
	 * @param[in] cst: costs of rrrmts
	 **/
	virtual void getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const;

	/**
	 * test for equality: double dispatch
	 * @param[in] r another rearrangement scenario
	 * @return true iff equal
	 */
	virtual bool equal( const rrrmt *r ) const;

	/**
	 * @todo documentation
	 */
	virtual rrrmt *intersect( const rrrmt *rr, bool prefix=false ) const;

	/**
	 * intersection of an alternative and an unordered scenario as the larges
	 * of the following options
	 * 1. intersection of every sub-scenario of rr against this
	 * 2. intersection of every sub-scenario of this against rr
	 *@param[in] rr an unordered scenario
	 *@param[in] prefix use prefix of included ordered scenarion
	 *@return the intersection
	 */
	virtual rrrmt *intersect( const unordered *rr, bool prefix=false ) const;

	/**
	 * intersection of an alternative and an ordered scenario as longest of
	 * 1. intersection of every subscenario of this against rr
	 * 2. intersection of the first subscenario rr against this
	 *@param[in] rr an ordered scenario
	 *@param[in] prefix use prefix of included ordered scenarions
	 *@return the intersection
	 */
	virtual rrrmt *intersect( const ordered *rr, bool prefix=false ) const;

	/**
	 * intersection of two alternative scenarios as the longest of
	 * 1. intersection of every sub-scenario of this with rr
	 * 2. intersection of every sub-scenario of rr with this
	 * 3. the alternative consisting of the maximal intersections of each
	 *    sub-scenario of this with each subscenario of rr
	 *    if it is not 1:1 bi directional a warning is printed
	 *
	 *@param[in] rr an alternative scenario
	 *@param[in] prefix use prefix of included ordered scenarions
	 *@return the intersection
	 */
	virtual rrrmt *intersect( const alternative *rr, bool prefix=false ) const;

	/**
	 * @return true iff any rearrangement in the scenario is complete
	 *  and the scenario is not explicitely marked as incomplete
	 */
	bool iscomplete() const;

	/**
	 * @return true iff the scenario is empty, and false if there is one nonempty alternative
	 */
	bool isempty() const;

	/**
	 * return the smallest length of the alternatives
	 * @param[in] mode take max (ALTMAX) / min (ALTMAX) length of alternatives
	 * @return the min/max length
	 */
	unsigned length( int mode ) const;

	/**
	 * get the number of alternative scenarios
	 * for an alternative scenario this is the sum of the numbers from the children
	 * @param return number of alternative scenarios
	 */
	virtual unsigned nralt() const;

	/**
	 * get the minimum cost of all alternative rrrmts
	 * @param return minimum cost of alternative scenarios
	 */
	virtual float getmincost(const costfoo *cst) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output(ostream &out, unsigned l=0, unsigned d=0, string pquot=""  ) const;

	/**
	 * output operator
	 * @see rrrmt::output
	 */
	ostream &output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx );

//	/**
//	 * make parsimonious
//	 * @see rrrmt::mkepar
//	 * @todo introduce weighting scheme
//	 * @todo deal correctly with incomplete (bp)scenarios
//	 * @return a parsimonious copy of the rearrangement
//	 */
//	rrrmt *mkepar() const;

//	void shortdesc( string &d ) const;

	/**
	 * @return ALT
	 */
	int type() const;

	/**
	 * @return the name of the rearrangement as string
	 */
	string typestrg(unsigned d=0) const;
};


/**
 * get the rearrangement index of an rearrangement
 * - from the given map
 * - or insert a new one if not already present
 *
 * @param[in] r the rearrangement
 * @param[in,out] ridx rearrangement index
 * @return the index
 */
unsigned dynamic_ridx( rrrmt *r, map<rrrmt*, unsigned, HDereferenceLess> &ridx );

/**
 * get an index of the rearrangements in the set
 * @param[in] rset a set of rearrangements
 * @return a map containing a unique name for each unique rrrmt
 */
map<rrrmt*, unsigned, HDereferenceLess> getrrrmtidx( const set<rrrmt*, HDereferenceLess> &rset );
//
///**
// * function to determine the rearrangement scenario(s) with maximal size (>0)
// * @param[in] its a vector of rearrangement scenarios
// * @param[out] max the maximum (or 0 if none)
// * @param[out] the indices of the scenarios with maximal q
// */
//void getmax( const vector<rrrmt*> &its, unsigned &max, vector<unsigned> &maxidx);

/**
 * function to determine the rearrangement scenario(s) with maximal size (>0)
 * and maximal number of alternatives
 * @param[in,out] a set of rearrangements non maximal or zero-length will be deleted
 * @return the maximum size found
 */
unsigned getmax( set<rrrmt*,HDereferenceLess> &its );
#endif/*REARRANGEMENTS_HPP_*/
