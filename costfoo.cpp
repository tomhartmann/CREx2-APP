/*
 * costfoo.cpp
 *
 *  Created on: Nov 24, 2013
 *      Author: maze
 */
#include <cstdlib>
#include <iostream>
#include "costfoo.hpp"

costfoo::~costfoo(){}

ostream & operator<<(ostream &out, const costfoo &c){
	return c.output( out );
}

/*****************************************************************************/

costfoo_equi::costfoo_equi(){}

costfoo_equi::~costfoo_equi(){}

ostream & costfoo_equi::output( ostream &out ) const{
	out << "costs"<<endl;
	out << "\t 1 for everything"<<endl;
	return out;
}

float costfoo_equi::operator [] (const rrrmt *r) const{
	return 1.0;
}


/*****************************************************************************/

costfoo_by_type::costfoo_by_type() {
}

costfoo_by_type::~costfoo_by_type() {
}

void costfoo_by_type::set( const string &tpe, const float c ){
	_cost[tpe] = c;
}

float costfoo_by_type::get( const string &tpe ) const{
	if(_cost.find( tpe ) == _cost.end()){
		cerr << "error: no cost defined for "<< tpe << endl;
		exit(EXIT_FAILURE);
	}
	return _cost.at( tpe );
}

ostream & costfoo_by_type::output(ostream &out) const{
	out << "costs"<<endl;
	for( map<string,float>::const_iterator it=_cost.begin(); it!=_cost.end(); it++ ){
		out <<"\t" <<it->first << " -> "<< it->second<<endl;
	}
	return out;
}

float costfoo_by_type::operator [] (const rrrmt *r) const{
	return get( r->typestrg(1) );
}

