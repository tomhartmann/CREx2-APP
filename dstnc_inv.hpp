/*
 * revdist.h
 *
 *  Created on: Jan 12, 2010
 *      Author: maze
 */

#ifndef REVDIST_HPP_
#define REVDIST_HPP_

#include <iostream>

#include "dstnc.hpp"

// necessary GRAPPA includes for reversal distance
#include "structs.h"
#include "invdist.h"
#include "med_util.h"


class dstnc_inv : public dstnc {
private:
	distmem_t *_distmem;         /* helping memory for GRAPPA code */
	struct genome_struct _g1, _g2; /* GRAPPA code data structures for the genomes */
	unsigned _n;			     /* max. length of the genomes which can be handled */
	bool _circular;              /* circularity */
public:

	/**
	 * constructor
	 * @param[in] n the maximum length of the genomes
	 * @param[in] circular true iff reversal distance is to be computed for circular genomes
	 */
	dstnc_inv( unsigned n, bool circular );

	/**
	 * destructor
	 */
	virtual ~dstnc_inv();

	/**
	 * adapt the reversal distance calculator to a new genome length
	 * @param[in] n the new length
	 */
	void adapt(unsigned n);

	/**
	 * get the reversal distance between two gene orders
	 * @param[in] src a genome
	 * @param[in] tgt another genome
	 * @param[in] context @see distance::calc
	 * @return the reversal distance
	 */
	unsigned calc( const genom &src, const genom &tgt );

	/**
	 * return a pointer to a copy of a dstnc
	 * this is a virtual copy constructor
	 * @return the pointer
	 */
	dstnc* clone() const;

//	/**
//	 * virtual constructor
//	 * @return an empty instance of an object
//	 */
//	dstnc* create() const;

	/**
	 * output function
	 * @param[in] stream to write into
	 * @return the stream
	 */
	ostream & output(ostream &os) const;
};

#endif /* REVDIST_HPP_ */
