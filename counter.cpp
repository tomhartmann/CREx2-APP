/*
 * counter.cpp
 *
 *  Created on: May 18, 2010
 *      Author: maze
 */

#include <cstdlib>
#include <iterator>

#include "counter.hpp"

using namespace std;

counter::counter(){
	_valid = false;
}

counter::counter( unsigned size, unsigned max, unsigned min, bool inc ){
	_cnt = vector<unsigned>(size, 0);
	_min = min;
	_max = max;
	_inc = inc;

	if( size > 0 ){
		_valid = true;
	}else{
		_valid = false;
	}

	for( unsigned i=0; i<_cnt.size(); i++){
		if( _inc ){
			_cnt[i] = _min+i;
		}else{
			_cnt[i] = _min;
		}
	}
	if( _cnt.back() >= _max ){
		invalidate();
	}
}

vector<unsigned>::iterator counter::begin(){
	return _cnt.begin();
}
vector<unsigned>::iterator counter::end(){
	return _cnt.end();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<unsigned> counter::get_counter( ){
	return _cnt;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// manually invalidate the counter
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void counter::invalidate(){
	_valid = false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool counter::isvalid( ){
	return _valid;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// prefix increment
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

counter& counter::operator++(){
    // berechne den neuen Bruch
    // geht auch: *this = *this + 1;
	int x;
	unsigned inject = 0;


	if( inject >= _cnt.size() ){
		invalidate();
		return *this;
	}
	if( ! this->isvalid() ){
		return *this;
	}

//	cout << "++counter"<<endl;
	for( int i = (int)_cnt.size()-1; i>=(int)inject; i-- ){
		_cnt[i]++;
//		copy(pos.begin(), pos.end(), ostream_iterator<int>(cout," ")); cout <<endl;
		if(_inc){
			x = ((int)_cnt.size()-1-i );
		}else{
			x = 0;
		}
		if( _cnt[i] + x < _max ){
//			cout << pos[i]<<"+"<<x<<"<"<<max<<endl;
			for(unsigned j=i; j<_cnt.size()-1; j++){
				if(_inc){
					_cnt[j+1] = _cnt[j]+1;
				}else{
					_cnt[j+1] = _min;
				}
			}
//			cout << "radj ";copy(pos.begin(), pos.end(), ostream_iterator<int>(cout," ")); cout <<endl;
			break;
		}
	}
//	cout << "++end "<<pos.back()<<" >= "<<max<<endl;
	if( _cnt.back() >= _max ){
		invalidate();
	}
    return *this;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// postfix increment
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

counter counter::operator++(int){
    counter oldcounter =*this;	// save old
    ++(*this);					// add
    return oldcounter;			// return old
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned &counter::operator[](unsigned i){
  if(i<0 || i>= _cnt.size() ) {
    cout << "counter: range violation"<< endl;
    exit(1);
  }
  return _cnt[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const unsigned& counter::operator [] (unsigned i) const {
	  if(i<0 || i>= _cnt.size() ) {
	    cout << "counter: range violation"<< endl;
	    exit(1);
	  }
	  return _cnt[i];
}

ostream& operator<<(ostream& out, counter c){
	out << "counter(";
	copy( c._cnt.begin(), c._cnt.end(), ostream_iterator<unsigned>(out, " ") );
	out << ") ";
	if(c.isvalid())
		out << "valid";
	return out;
}
