/*
 * counter.h
 *
 *  Created on: May 18, 2010
 *      Author: maze
 */

#ifndef COUNTER_HPP_
#define COUNTER_HPP_

#include <vector>
#include <iostream>


using namespace std;

class counter{
private:
	vector<unsigned> _cnt; // the counter
	unsigned _min,
		_max;
	bool _inc,
		_valid;
public:
	/**
	 * constructor
	 * @param[in] size the desired length of the vector
	 * @param[in] max the maximum values
	 * @param[in] min the minimal value of the integers in the vector (default = 0)
	 * @param[in] inc let the values be increasing, i.e. after initialisation:
	 * \f$ (min, min+1,\ldots,min+size-1) \f$ otherwise \f$ (min, \ldots, min ) \f$ (default: false)
	 */
	counter( unsigned size, unsigned max, unsigned min = 0, bool inc = false );

	/**
	 * default constructor
	 */
	counter();

	/**
	 * get the iterator to the beginning of the counter
	 */
	vector<unsigned>::iterator begin();

	/**
	 * get the iterator to the beginning of the counter
	 */
	vector<unsigned>::iterator end();

	/**
	 * get the internal representation
	 * @return the vector
	 */
	vector<unsigned> get_counter( );

	/**
	 * manually invalidate the counter
	 */
	void invalidate();

	/**
	 * check if the counter is valid
	 * @return true iff the counter is valid
	 */
	bool isvalid();

	counter& operator++();

	counter operator++(int);

	unsigned &operator[](unsigned i);
	const unsigned &operator[](unsigned i) const;

	friend ostream &operator<<(ostream &stream, counter c);
};

#endif /* COUNTER_H_ */
