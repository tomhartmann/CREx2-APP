/*
 * rrrmtrandsit.hpp
 *
 * create a random rearrangement scenario from a SIT
 *
 *  Created on: Nov 28, 2013
 *      Author: maze
 */

#ifndef RRRMTRAND_HPP_
#define RRRMTRAND_HPP_

#include "rearrangements.hpp"

namespace std {

class rrrmt_rand: public unordered {
public:

	/**
	 * construct a random scenario of k atomic operations (rev, tra, revtra,
	 * tdrl) for a genome of length n given probabilities for each operation
	 *
	 * it may seem unlogical that an unordered scenario is constructed for
	 * simulations (instead of an ordered scenario), but for the comparison of
	 * the prediction and the simulation its convenient as this is done with
	 * the intersection.
	 *
	 * note: the set of rrrmts stored internally may be unequal to the set
	 * of actually applied rrrmts. this happens if the same rrrmt is applied more
	 * than once. in this case the rrrmt is only stored once in the set.
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
	rrrmt_rand( int k, vector<float> prob, genom &g, vector<genom> &trace, bool randrange, unsigned ml);

	virtual ~rrrmt_rand( );
};


class rrrmt_rand_sit: public unordered {
public:

	/**
	 * generate a random rearrangement scenario for a genome from a random SIT
	 * @param k number of rearrangements to apply
	 *
	 * @param[in,out] prob probabilities for the rearrangements (will be normalized internaly)
	 * @param[in,out] g the genome
	 * @param[out] trace the genomes on the random rearrangement scenario
	 * @param dm min degree of the SIT
	 * @param dm min degree of the SIT
	 *
	 * int k, vector<float> prob, genom &g, vector<genom> &trace, bool randrange, unsigned ml
	 */
	rrrmt_rand_sit( unsigned k, vector<float> prob, genom &g, vector<genom> &trace, unsigned dm, unsigned dM );

	virtual ~rrrmt_rand_sit( );
};

} /* namespace std */

#endif /* RRRMTRANDSIT_HPP_ */
