/*
 * distance.h
 *
 *  Created on: Jun 12, 2009
 *      Author: maze
 */

#ifndef DISTANCE_HPP_
#define DISTANCE_HPP_


#include <iostream>

#include "genom.hpp"


/**
 * base class for computing distances between gene orders
 * derived classes must implement:
 * - adapt( n )
 * - calc( g1, g2 )
 * - clone( )
 * - output( )
 *
 * one might also overwrite calc(g1,g2,genomes)
 */
class dstnc {
private:
public:
	/**
	 * constructor: does nothing
	 */
	dstnc();

	/**
	 * destructor: does nothing
	 */
	virtual ~dstnc();

	/**
	 * adapt the distance calculator to a new genome length
	 * @param[in] n the new length
	 */
	virtual void adapt( unsigned n ) = 0;
	/**
	 * get the distance from src to tgt
	 * @param[in] src the source genome
	 * @param[in] tgt the target genome
	 * return the distance
	 */
	virtual unsigned calc( const genom &src, const genom &tgt ) = 0;

	/**
	 * calculate the distance from g1 to g2 some property of the
	 * genomes in a set of genomes
	 *
	 * note: if this function is not overwritten then the result of
	 * calc(g1,g2) is returned
	 *
	 * @param[in] src the source genome
	 * @param[in] tgt the traget genome
	 * @param[in] genomes the set of genomes
	 * @return the distance
	 */
	virtual unsigned calc( const genom &src, const genom &tgt, const vector<genom> &genomes );

	/**
	 * calculate the sum of the distances from src to the genomes in tgts
	 * @param[in] src a genomes
	 * @param[in] tgts genomes
	 * @return sum of distances
	 */
	unsigned calcsum( const genom &src, const vector<genom> &tgts );

	/**
	 * calculate the sum of the distances from src to the genomes in tgts
	 * preserving properties in genomes set
	 *
	 * @param[in] src a genomes
	 * @param[in] tgts genomes
	 * @param[in] genomes a list of genomes
	 * @return sum of distances
	 */
	unsigned calcsum( const genom &src, const vector<genom> &tgts, const vector<genom> &genomes );

	/**
	 * calculate the minimum of the distances from src to the genomes in tgts
	 * @param[in] src a genomes
	 * @param[in] tgts genomes
	 * @return minimum of distances
	 */
	unsigned calcmin( const genom &src, const vector<genom> &tgts );

	/**
	 * calculate the minimum of the distances from src to the genomes in tgts
	 * preserving properties in genomes set
	 * @param[in] src a genomes
	 * @param[in] tgts genomes
	 * @param[in] genomes a list of genomes
	 * @return minimum of distances
	 */
	unsigned calcmin( const genom &src, const vector<genom> &tgts, const vector<genom> &genomes  );

	/**
	 * calculate the maximum of the distances from src to the genomes in tgts
	 * @param[in] src a genomes
	 * @param[in] tgts genomes
	 * @return minimum of distances
	 */
	unsigned calcmax( const genom &src, const vector<genom> &tgts );

	/**
	 * calculate the maximum of the distances from src to the genomes in tgts
	 * preserving properties in genomes set
	 * @param[in] src a genomes
	 * @param[in] tgts genomes
	 * @param[in] genomes a list of genomes
	 * @return minimum of distances
	 */
	unsigned calcmax( const genom &src, const vector<genom> &tgts, const vector<genom> &genomes );

	/**
	 * return a pointer to a copy of a dstnc
	 * this is a virtual copy constructor
	 * @return the pointer
	 */
	virtual dstnc* clone() const = 0;

//	/**
//	 * virtual constructor
//	 * @return an empty instance of an object
//	 */
//	virtual dstnc* create() const = 0;

	/**
	 * output operator, dispatches to the output function of derived classes
	 * @param[in] os the stream to write in
	 * @param[in] the distance method to output
	 */
	friend ostream &operator<<(ostream &os, const dstnc &d){
		return d.output(os);
	};

	/**
	 * the output function
	 * @param[in] os the stream to write into
	 * @return the stream
	 */
	virtual ostream & output(ostream &os) const = 0;
};


/**
 * a distance class for combined distance measures
 * i.e. you can use it to get the sum, min, max of some distance measures
 */
template <typename Functor>
class dstnc_comb : public dstnc {
private:
	vector< dstnc* > _df;	/* distance functions */
	Functor _foo;
public:
	/**
	 * constructor
	 * @param[in] df a vector of distance functions
	 */

	dstnc_comb( const vector<dstnc *> &df, Functor foo  );

	/**
	 * destructor
	 */
	~dstnc_comb();

	/**
	 * adapt the distance function(s) to a new n
	 * @param[in] n the new permutation length
	 */
	void adapt( unsigned n );

	/**
	 * calculate the minimum of the distances from src to tgt
	 * @param[in] src the source
	 * @param[in] tgt the target
	 * @param the minimum distance
	 */
	unsigned calc( const genom &src, const genom &tgt );

	/**
	 * calculate the distance from g1 to g2 some property of the
	 * genomes in a set of genomes
     *
	 * @param[in] src the source genome
	 * @param[in] tgt the traget genome
	 * @param[in] genomes the set of genomes
	 * @return the distance
	 */
	unsigned calc( const genom &src, const genom &tgt, const vector<genom> &genomes );

	/**
	 * return a pointer to a copy of a dstnc this is a virtual copy constructor
	 * @return the pointer
	 */
	dstnc* clone() const;

	/**
	 * output function
	 * @param[in,out] os stream to write into
	 * @return the stream
	 */
	ostream & output(ostream &os) const;
};

/**
 * a distance class for modified distance measures
 * i.e. you can use it to get the foo(d)... of some distance measures d
 */
template <typename Functor>
class dstnc_mod : public dstnc {
private:
	dstnc* _df;	/* distance function */
	Functor _foo;
public:
	/**
	 * constructor
	 * @param[in] df a vector of distance functions
	 */

	dstnc_mod( const dstnc * &df, Functor foo  );

	/**
	 * destructor
	 */
	~dstnc_mod();

	/**
	 * adapt the distance function(s) to a new n
	 * @param[in] n the new permutation length
	 */
	void adapt( unsigned n );

	/**
	 * calculate the minimum of the distances from src to tgt
	 * @param[in] src the source
	 * @param[in] tgt the target
	 * @param the minimum distance
	 */
	unsigned calc( const genom &src, const genom &tgt );

	/**
	 * calculate the distance from g1 to g2 some property of the
	 * genomes in a set of genomes
     *
	 * @param[in] src the source genome
	 * @param[in] tgt the traget genome
	 * @param[in] genomes the set of genomes
	 * @return the distance
	 */
	unsigned calc( const genom &src, const genom &tgt, const vector<genom> &genomes );

	/**
	 * return a pointer to a copy of a dstnc this is a virtual copy constructor
	 * @return the pointer
	 */
	dstnc* clone() const;

	/**
	 * output function
	 * @param[in,out] os stream to write into
	 * @return the stream
	 */
	ostream & output(ostream &os) const;
};



template <typename Functor>
dstnc_comb<Functor>::dstnc_comb( const vector< dstnc * > &df, Functor foo   ){

	_foo = foo;
	for( unsigned i=0; i<df.size(); i++ ){
		_df.push_back( df[i]->clone() );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
dstnc_comb<Functor>::~dstnc_comb( ){
	for( unsigned i=0; i<_df.size(); i++ ){
		delete _df[i];
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
void dstnc_comb<Functor>::adapt( unsigned n ){
	for( unsigned i=0; i<_df.size(); i++ ){
		_df[i]->adapt( n );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
unsigned dstnc_comb<Functor>::calc( const genom &src, const genom &tgt ){

	unsigned r;

	r = _foo( _df[0]->calc(src, tgt), _df[1]->calc(src, tgt) );
	for( unsigned i = 2; i<_df.size(); i++ ){
		r = _foo( r, _df[i]->calc(src, tgt) );
	}
	return r;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
unsigned dstnc_comb<Functor>::calc( const genom &src, const genom &tgt,
		const vector<genom> &genomes ){

	unsigned r;

	r = _foo( _df[0]->calc(src, tgt, genomes),
			_df[1]->calc(src, tgt, genomes) );

	for( unsigned i = 2; i<_df.size(); i++ ){
		r = _foo( r, _df[i]->calc(src, tgt, genomes) );
	}
	return r;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
dstnc* dstnc_comb<Functor>::clone() const{
	return new dstnc_comb( _df, _foo );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
ostream & dstnc_comb<Functor>::output(ostream &os) const{
	os << "comb dist(";
	for( unsigned i=0; i<_df.size(); i++ ){
		os << *_df[i];
		if( i!=_df.size()-1 )
			os<<",";
	}
	os<<")";
	return os;
}




template <typename Functor>
dstnc_mod<Functor>::dstnc_mod( const dstnc * &df, Functor foo   ){
	_foo = foo;
	_df = df->clone();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
dstnc_mod<Functor>::~dstnc_mod( ){
	delete _df;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
void dstnc_mod<Functor>::adapt( unsigned n ){
	_df->adapt( n );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
unsigned dstnc_mod<Functor>::calc( const genom &src, const genom &tgt ){

	return _foo( _df->calc(src, tgt) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
unsigned dstnc_mod<Functor>::calc( const genom &src, const genom &tgt,
		const vector<genom> &genomes ){

	return _foo( _df->calc(src, tgt, genomes), _df->calc(src, tgt, genomes) );;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
dstnc* dstnc_mod<Functor>::clone() const{
	return new dstnc_mod( _df, _foo );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Functor>
ostream & dstnc_mod<Functor>::output(ostream &os) const{
	os << "mod dist("<< *_df <<")";
	return os;
}



#endif /* DISTANCE_HPP_ */
