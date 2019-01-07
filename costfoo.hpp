/*
 * costfoo.hpp
 *
 *  Created on: Nov 24, 2013
 *      Author: maze
 */

#ifndef COSTFOO_HPP_
#define COSTFOO_HPP_

#include <map>
#include <string>

#include "rearrangements.hpp"

using namespace std;

/**
 * abstract base class for cost functions for rearrangements
 *
 * subclasses must define output(), operator[]
 */
class costfoo{
public:
	/**
	 *
	 */
	virtual ~costfoo();

	/**
	 * output operator. redirects to output().
	 * @param[in,out] out the stream to write into
	 * @param[in] c the costfoo
	 * @return the stream written into
	 */
	friend ostream & operator<<(ostream &out, const costfoo &c);

	/**
	 * output function for the cost functions
	 * @param[in,out] out the stream to write into
	 * @param[in] c the costfoo
	 * @return the stream written into
	 */
	virtual ostream & output( ostream &out ) const = 0;

	/**
	 * get the cost of a given rearrangement
	 * @param[in] r the rearrangement
	 * @return the cost
	 */
	virtual float operator [] (const rrrmt *r) const = 0;
};

/**
 * cost function returning 1 for everything
 */
class costfoo_equi : public costfoo{
public:
	costfoo_equi();
	~costfoo_equi();

	/**
	 * output operator
	 * @param[in,out] out the stream to write to
	 * @param[in] c the cost foo to output
	 * @return the stream written into
	 */
	ostream & output(ostream &out) const;

	/**
	 * get the cost of rearrangement r
	 * @param[in] r the rearrangement
	 * @return the cost
	 */
	float operator [] (const rrrmt *r) const;
};

/**
 * cost function that allows a cost per type
 */
class costfoo_by_type : public costfoo{
private:
	map<string,float> _cost;
public:

	/**
	 * init makes basically nothing
	 */
	costfoo_by_type();

	/**
	 * destructor makes nothing except freeing the mem
	 */
	virtual ~costfoo_by_type();

	/**
	 * get the cost of a rearrangement of type tpe
	 * @param[in] tpe rearrangement type
	 * @return the cost
	 */
	float get( const string &tpe ) const;

	/**
	 * define the cost for a rearrangement of type r as c
	 * @param[in] tpe the type
	 * @param[in] c the cost
	 */
	void set( const string &r, const float c );

	/**
	 * output operator
	 * @param[in,out] out the stream to write to
	 * @param[in] c the cost foo to output
	 * @return the stream written into
	 */
	ostream & output(ostream &out ) const;

	/**
	 * get the cost of rearrangement r
	 * @param[in] r the rearrangement
	 * @return the cost
	 */
	float operator [] (const rrrmt *r) const;
};

#endif /* COSTFOO_HPP_ */
