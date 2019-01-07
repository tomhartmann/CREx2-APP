/*
 * distance.cpp
 *
 *  Created on: Feb 1, 2010
 *      Author: maze
 */

#include "dstnc.hpp"

#include <iostream>
#include <limits>

using namespace std;

dstnc::dstnc(){}

dstnc::~dstnc(){}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calc( const genom &src, const genom &tgt,
		const vector<genom> &genomes ){

	return calc(src, tgt);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calcsum( const genom &src, const vector<genom> &tgts ){
	unsigned d = 0;
	for (unsigned int i = 0; i < tgts.size(); i++)
		d += this->calc( src, tgts[i] );

	return d;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calcsum( const genom &src, const vector<genom> &tgts,
		const vector<genom> &genomes ){

	unsigned d = 0;
	for (unsigned int i = 0; i < tgts.size(); i++)
		d += this->calc( src, tgts[i], genomes );

	return d;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calcmin( const genom &src, const vector<genom> &tgts ){
	unsigned md = std::numeric_limits< unsigned >::max(),
		d;

	for (unsigned int i = 0; i < tgts.size(); i++){
		d = this->calc( src, tgts[i] );
		if( d < md ){
			md = d;
		}
	}

	return md;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calcmin( const genom &src, const vector<genom> &tgts,
		const vector<genom> &genomes ){

	unsigned md = std::numeric_limits< unsigned >::max(),
		d;

	for (unsigned int i = 0; i < tgts.size(); i++){
		d = this->calc( src, tgts[i], genomes );
		if( d < md ){
			md = d;
		}
	}

	return md;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calcmax( const genom &src, const vector<genom> &tgts ){
	unsigned md = 0,
		d;

	for (unsigned int i = 0; i < tgts.size(); i++){
		d = this->calc( src, tgts[i] );
		if( d > md ){
			md = d;
		}
	}

	return md;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dstnc::calcmax( const genom &src, const vector<genom> &tgts,
		const vector<genom> &genomes ){

	unsigned md = 0,
		d;

	for (unsigned int i = 0; i < tgts.size(); i++){
		d = this->calc( src, tgts[i], genomes );
		if( d > md ){
			md = d;
		}
	}

	return md;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

