/*
 * revdist.cpp
 *
 *  Created on: Jan 12, 2010
 *      Author: maze
 */

#include "dstnc_inv.hpp"

dstnc_inv::dstnc_inv( unsigned n, bool circular ) {
	_n = n;
	_distmem = new_distmem ( n );
	_circular = circular;
}

dstnc_inv::~dstnc_inv() {
	free_distmem ( _distmem );
}

void dstnc_inv::adapt(unsigned n){
	free_distmem ( _distmem );
	_n = n;
	_distmem = new_distmem ( n );
}


unsigned dstnc_inv::calc( const genom &src, const genom &tgt ){

	// Test for valid chromosom
	if ( _n < src.size() || _n < tgt.size() ) {
		cerr << "revdstnc: improper initialisation"<<endl;
		exit( EXIT_FAILURE );
	}

		// it seem that the vector is direct usable as c-array
		// it's ugly to cast to int*, but the genes should not
		// get altered during the distance calculation
	_g1.genes = src.get_pointer();
	_g2.genes = tgt.get_pointer();

		// (let) calculate the distance
	if (_circular) {
		return invdist_circular_mb ( &_g1, &_g2, src.size(), _distmem );
	} else{
		return invdist_noncircular_mb ( &_g1, &_g2, 0, src.size(), _distmem );
	}
}

dstnc* dstnc_inv::clone() const{
	return new dstnc_inv( _n, _circular );
}


//dstnc* dstnc_inv::create(){
//	return new dstnc_inv();
//}

ostream & dstnc_inv::output(ostream &os) const{
	os << "inversion distance (GRAPPA) for n="<<_n;
	return os;
}
