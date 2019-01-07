/**
 * rearrangements.cpp
 * @todo simplify of ordered, unordered, and alternative: if a scenario contains a scenario of the same type -> flatten
 */

#include <iterator>
#include <algorithm>
#include <vector>

#include "rearrangements.hpp"
#include "costfoo.hpp"

//#define DEBUG_INTERSECTION
//#define DEBUG_ATOMINTERSECTION
//#define DEBUG_DEAPPLY
//#define DEBUG_APPLY

using namespace std;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// n length of the genomes
// x probability of a breakpoint to be fragile
// y probability of choosing a fragile breakpoint
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned random_breakpoint( unsigned n, float x, float y, vector<vector<int> > &bp){

	unsigned bpt; // breakpoint type (fragile, non-fragile)

	if( x<0 ){
		return ask_rng() % n;
	}

	if( bp.size( ) != n ){
		bp = vector<vector<int> >(2);
		for( unsigned i=0; i<n; i++ ){
			bp[ (ask_rng_f() <= x) ? 0 : 1 ].push_back(i);
		}
	}

	bpt = (ask_rng_f() <= y) ? 0 : 1;
	return bp[ bpt ][ ask_rng(0, bp[ bpt ].size()) ];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt::~rrrmt(){
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rrrmt::append( vector<vector<rrrmt*> > &v ) const{
	v.push_back( vector<rrrmt*>((unsigned)1, clone()) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rrrmt::apply( genom &g, set<genom> &gs )const{
	apply(g);
	gs.insert(g);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rrrmt *rrrmt::create() const{
	return new emptyscen;
}


bool rrrmt::equal( const emptyscen *r ) const{
	return false;
}
bool rrrmt::equal( const indel *r ) const{
	return false;
}
bool rrrmt::equal( const rev *r ) const{
	return false;
}
bool rrrmt::equal( const transp *r) const{
	return false;
}
bool rrrmt::equal( const revtransp *r) const{
	return false;
}
bool rrrmt::equal( const tdrl *r) const{
	return false;
}
bool rrrmt::equal( const unordered *r) const{
	return false;
}
bool rrrmt::equal( const ordered *r) const{
	return false;
}
bool rrrmt::equal( const alternative *r) const{
	return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void rrrmt::getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const{
	//do nothing
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool rrrmt::hascomposite() const{
	return false;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection functions of the basis class. all return empty. they have to be overwritten
// if another behaviour is wished
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt* rrrmt::intersect( const emptyscen *rr, bool prefix )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const indel *rr, bool prefix  )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const rev *rr, bool prefix  )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const transp *rr, bool prefix  )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const revtransp *rr, bool prefix  )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const tdrl *rr, bool prefix  )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const unordered *rr, bool prefix  )const{
	return new emptyscen();
}
rrrmt* rrrmt::intersect( const ordered *rr, bool prefix  )const{
	return new emptyscen();
}

rrrmt* rrrmt::intersect( const alternative *rr, bool prefix  )const{
	return new emptyscen();
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool rrrmt::iscomplete() const{
	return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool rrrmt::isempty() const{
	return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned rrrmt::length( int mode ) const{
	return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *rrrmt::mkealt(){
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//rrrmt *rrrmt::mkepar() const{
//	return clone();
//}

unsigned rrrmt::nralt() const{
	return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

float rrrmt::getmincost(const costfoo *cst) const{
	//return corresponding cost value of costfoo
	return (*cst)[ this ];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &operator<<(ostream &os, const rrrmt &r){
        return r.output(os);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rrrmt::output(unsigned l, unsigned d, string pquot){
	this->output( cout, l, d, pquot );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ostream & rrrmt::output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx ){
	out << typestrg( 1 )<< dynamic_ridx( this, ridx );
	return out;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *rrrmt::simplify(){
	return clone();
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//int rrrmt::type() const{
//	return RRRMT;
//}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************

void emptyscen::apply( genom &g )const{}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

emptyscen* emptyscen::clone() const{
	return new emptyscen(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void emptyscen::deapply(genom &g) const{}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void emptyscen::elements( vector<set<int> > &elms, bool ab ) const{}

bool emptyscen::equal( const rrrmt *r ) const{
	return r->equal(this);
}

bool emptyscen::equal( const emptyscen *r ) const{
	return true;
}


void emptyscen::getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode ){}
void emptyscen::getrrrmt( vector<rrrmt*> &out, int mode ){}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection with empty will always return emty .. so we dont need double dispatch here
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *emptyscen::intersect(const rrrmt *r, bool prefix ) const{
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool emptyscen::isempty() const{
	return true;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned emptyscen::length( int mode ) const{
	return 0;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool emptyscen::operator<(const rrrmt *rrrmtp) const{
	return type() < rrrmtp->type();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &emptyscen::output(ostream &out, unsigned l, unsigned d, string pquot ) const{
	if( d == 0 ){
		out<<"Empty"<<endl;
	}else if( d == 1 ){
		out<<"Empty";
	}else{
		out << " ";
	}
	return out;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int emptyscen::type() const{
	return EMPTYSCE;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string emptyscen::typestrg(unsigned d) const{
	if(d==0)
		return "empty";
	else
		return "";
}



// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************


indel::indel( int e, bool del, const genom &g ){
	_nmap = g.get_nmap();
	_d = vector<set<int> >( 3 );
	_del = del;

//	cerr << "new indel "<<e<<endl;

	unsigned x=0;
	for( unsigned i=0; i<g.size(); i++ ){
		if( abs(e) == abs( g[i] ) ){
			x++;
		}
		_d[x].insert(abs(g[i]));
		if( x==1 ){
			_inv = (g[i]<0) ? true : false;
		}

		if( abs(e) == abs( g[i] ) ){
			x++;
		}
//		if( k<i ){
//			_d[0].insert(abs(g[k]));
//		}else if( k>i ){
//			_d[2].insert(abs(g[k]));
//		}else{
//			_d[1].insert(abs(g[k]));
//		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

indel::~indel(){
	_d.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a indel
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void indel::apply(genom &g) const{
#ifdef DEBUG_APPLY
	cerr << "apply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_APPLY
	insertdelete( g, _d, _del, _inv );
#ifdef DEBUG_APPLY
	cerr <<"->"<< g<< endl;
#endif//DEBUG_APPLY
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

indel* indel::clone() const{
	return new indel(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a indel
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void indel::deapply(genom &g) const{
	insertdelete( g, _d, ! _del, _inv );
#ifdef DEBUG_DEAPPLY
	cout << "deapply "<< *this << endl;
	cout << g<< endl;
#endif//DEBUG_DEAPPLY
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void indel::elements( vector<set<int> > &elms, bool ab ) const{
	// add the deleted elements
	elms.push_back( _d[1] );
}

bool indel::equal( const rrrmt *r ) const{
	return r->equal(this);
}
bool indel::equal( const indel *r ) const{
	return (_d == r->_d) ;
}


void indel::getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode ){
	indel *i = this->clone();
	if( ! out.insert( i ).second ){
		delete i;
	}
}
void indel::getrrrmt( vector<rrrmt*> &out, int mode ){
	out.push_back( clone() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for indel
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *indel::intersect(const rrrmt *rrrmtp, bool prefix ) const{
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection indel - indel
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *indel::intersect(const indel *rr, bool prefix ) const{
	#ifdef DEBUG_ATOMINTERSECTION
	cout << "indel::intersect "<<*this<< " vs. "<<rr;
	#endif//DEBUG_ATOMINTERSECTION
	if(_d[1] == rr->_d[1] ){
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> copy"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return clone();
	}else{
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> empty"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return new emptyscen();
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection indel - unordered .. implemented in unordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *indel::intersect(const unordered *rr, bool prefix ) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection indel - ordered .. implemented in ordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *indel::intersect(const ordered *rr, bool prefix ) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool indel::operator<(const rrrmt *rrrmtp) const{
	if( type() != rrrmtp->type() ){
		return type() < rrrmtp->type();
	}
	const indel* revp = dynamic_cast<const indel*>(rrrmtp);
	return _d[1] < revp->_d[1];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// output foo
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &indel::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) out << "\t";
	}
	out << typestrg(d);

	if( d < 2 ){
		out <<pquot<<"(";
		for( unsigned i=0; i<_d.size(); i++ ){
			for( set<int>::iterator it = _d[i].begin(); it!=_d[i].end(); it++ ){
				if( _nmap == NULL ){
					out << *it<<" ";
				}else{
					if( *it < 0 ){
						cout << "-";
					}
					out << (*_nmap)[ abs(*it) ]<<" ";
				}
			}
			out << ",";
		}
		out <<pquot<< ")";
		if( d == 0 ){
			out << endl;
		}
	}

	return out;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void indel::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	vector<set<int> > nd(3);

	for( unsigned j=0; j<_d.size(); j++ ){
		for( set<int>::iterator it = _d[j].begin(); it!=_d[j].end(); it++ ){
			for(unsigned i=0; i<mapping[*it].size(); i++){
				nd[j].insert(abs(mapping[*it][i]));
			}
		}
	}
	_d = nd;
	_nmap = nmap;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int indel::type() const{
	return INDEL;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string indel::typestrg(unsigned d) const{
	if(d > 0){
		if(_del)
			return "iD";
		else
			return "Id";
	}else{
		if(_del)
			return "inDel";
		else
			return "Indel";
	}
}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rev::rev(int i, int j, const genom &g){
//	cout << "rev::rev("<<i<<","<<j<<")"<<endl;

	if( i < 0 || j < 0 ){
		cerr << "cannot construct inversion with negative indicees: "<<i<<","<<j<<endl;
		exit(1);
	}

	if(j<i){
		cerr << "cannot construct inversion with end index "<<j<<" < start index "<<i<<endl;
		exit(1);
	}

	if( j >= (int)g.size() ){
		cerr << "cannot construct inversion exeeding the genome boundaries ("<<g.size()<<"): "<<i<<","<<j<<endl;
		exit(1);
	}


	for(int x=i; x<=j; x++)
		r.insert( abs(g[x]) );

	_nmap = g.get_nmap();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rev::rev( const genom &g, unsigned ml ){
	vector<unsigned> iv;	// random interval

	if( ml == 0 ){
		cerr << "rev::rev ml=0"<<endl;
		exit(EXIT_FAILURE);
	}

	iv = rng_inc_seq( 2, g.size(), ml );

//	cout <<"rev "<< iv[1]-iv[0]<<" "<<iv[1]<<" "<<iv[0]<<" "<<ml<<endl;
	for(unsigned i=iv[0]; i<=iv[1]; i++)
		r.insert( abs(g[i]) );
	_nmap = g.get_nmap();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rev::rev(const set<int> &rset, vector<string> *nm){
	r = rset;
	_nmap = nm;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rev::~rev(){
	r.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a reversal
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rev::apply(genom &g) const{
#ifdef DEBUG_APPLY
	cerr << "apply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_APPLY
	reverse(g, r);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rev* rev::clone() const{
	return new rev(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a reversal
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rev::deapply(genom &g) const{
	reverse(g, r);
#ifdef DEBUG_DEAPPLY
	cerr << "deapply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_DEAPPLY
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rev::elements( vector<set<int> > &elms, bool ab ) const{
	elms.push_back( set<int>() );
//	vector<int> e;
//	e.insert(e.end(), r.begin(), r.end());
	for( set<int>::const_iterator it=r.begin(); it!=r.end(); it++  ){
		elms.back().insert( ab ? abs(*it) : *it  );
	}

//	elms.push_back(e);
}

bool rev::equal( const rrrmt *r ) const{
	return r->equal(this);
}

bool rev::equal( const rev *r ) const{
	return (this->r == r->r);
}

void rev::getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode ){
	rev *t = this->clone();
	if( ! out.insert( t ).second ){
		delete t;
	}
}
void rev::getrrrmt( vector<rrrmt*> &out, int mode ){
	out.push_back( this->clone() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for reversals
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *rev::intersect(const rrrmt *rrrmtp, bool prefix ) const{
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection reversal - reversal
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *rev::intersect(const rev *rr, bool prefix ) const{
	#ifdef DEBUG_ATOMINTERSECTION
	cout << "rev::intersect "<<*this<< " vs. "<<*rr;
	#endif//DEBUG_ATOMINTERSECTION
	if(r == rr->r){
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> copy"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return clone();
	}else{
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> empty"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return new emptyscen;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection reversal - unordered .. implemented in unordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *rev::intersect(const unordered *rr, bool prefix ) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection reversal - ordered .. implemented in ordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *rev::intersect(const ordered *rr, bool prefix ) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//bool rev::operator==(const rrrmt *rrrmtp) const{
//	if( const rev* revp = dynamic_cast<const rev*>(rrrmtp) ){
//		return r == revp->r ;
//	}
//	return false;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool rev::operator<(const rrrmt *rrrmtp) const{
	if( type() != rrrmtp->type() ){
		return type() < rrrmtp->type();
	}
	const rev* revp = dynamic_cast<const rev*>(rrrmtp);
	return r < revp->r ;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// output foo
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &rev::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) out << "\t";
	}
	out << typestrg(d);

	if( d < 2 ){
		out <<pquot<<"(";
		for( set<int>::iterator it = r.begin(); it!=r.end(); it++ ){
			if( _nmap == NULL ){
				out << *it<<" ";
			}else{
				if( *it < 0 ){
					cout << "-";
				}
				out << (*_nmap)[ abs(*it) ]<<" ";
			}
		}
		out <<pquot<< ")";
		if( d == 0 ){
			out << endl;
		}
	}

	return out;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void rev::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	set<int> nr;

	for( set<int>::iterator it = r.begin(); it!=r.end(); it++ ){
		for(unsigned i=0; i<mapping[*it].size(); i++){
			nr.insert(abs(mapping[*it][i]));
		}
//		nr.insert( mapping[*it].begin(), mapping[*it].end() );
	}
	r = nr;
	_nmap = nmap;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int rev::type() const{
	return REV;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string rev::typestrg(unsigned d) const{
	if(d > 0){
		return REVNM;
	}else{
		return "inversion";
	}
}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

transp::transp(int i, int j, int k, const genom &g){
	set<int> tmp;

	if ( i < 0 || j < 0 || k < 0 ){
		cerr << "cannot construct transposition with negative indices: "<<i<<","<<j<<","<<k<<endl;
		exit(1);
	}

	if( !(i<=j && j<=k) ){
		cerr << "cannot construct transposition with unsorted indices: "<<i<<","<<j<<","<<k<<endl;
		exit(1);
	}

	if( k > (int)g.size() ){
		cerr << "cannot construct transposition exceeding the genome boundaries: ("<<g.size()<<"): "<<i<<","<<j<<","<<k<<endl;
		exit(1);
	}

	for(int x=i; x<j; x++){
		tmp.insert( abs(g[x]) );
	}
	t.insert(tmp);
	tmp.clear();

	for(int x=j; x<k; x++){
		tmp.insert( abs(g[x]) );
	}
	t.insert(tmp);
	tmp.clear();
	_nmap = g.get_nmap();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

transp::transp(const genom &g, unsigned ml){
	vector<unsigned> iv;	// random increasing seq of length 3 with iv[2]-iv[0]<=d
	set<int> tmp;

	iv = rng_inc_seq( 3, g.size(), ml );
		// assign elements to the transp
	for(unsigned i=iv[0]; i<iv[1]; i++){
		tmp.insert( abs(g[i]) );
	}
	t.insert(tmp);
	tmp.clear();
	for(unsigned i=iv[1]; i<=iv[2]; i++){
		tmp.insert( abs(g[i]) );
	}
	t.insert(tmp);
	tmp.clear();
	_nmap = g.get_nmap();
	//	cout << "transp::transp():       "<<*this << endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

transp::transp( const vector<set<int> > &tra, vector<string> *nm ){
	if( tra.size() != 2 ){
		cerr << "internal error: transp::transp() called with "<<tra.size()<<" sets"<<endl;
		exit(1);
	}
	t.insert( tra[0] );
	t.insert( tra[1] );

	_nmap = nm;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

transp::~transp(){
	t.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void transp::apply(genom &g) const{
#ifdef DEBUG_APPLY
	cerr << "apply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_APPLY
	transpose(g, t);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

transp* transp::clone() const{
	return new transp(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void transp::deapply(genom &g) const{
	transpose(g, t);
#ifdef DEBUG_DEAPPLY
	cerr << "deapply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_DEAPPLY
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void transp::elements( vector<set<int> > &elms, bool ab )const{
//	vector<int> e;
	for(set<set<int> >::const_iterator it = t.begin(); it!= t.end(); it++){
		elms.push_back(set<int>());
		for( set<int>::const_iterator jt=it->begin(); jt!=it->end(); jt++ ){
			elms.back().insert( ab?abs(*jt):*jt );
		}
	}
//	elms.push_back(e);
}

bool transp::equal( const rrrmt *r ) const{
	return r->equal(this);
}
bool transp::equal( const transp *r ) const{
	return ( this->t == r->t );
}

void transp::getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode ){
	transp *t = this->clone();
	if( ! out.insert( t ).second ){
		delete t;
	}
}
void transp::getrrrmt( vector<rrrmt*> &out, int mode ){
	out.push_back( clone() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for transpositions
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *transp::intersect(const rrrmt *rrrmtp, bool prefix) const{
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection transpositions - transpositions
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *transp::intersect(const transp *rr, bool prefix) const{
	#ifdef DEBUG_ATOMINTERSECTION
	cout << "transp::intersect "<<*this<< " vs. "<<*rr;
	#endif//DEBUG_ATOMINTERSECTION
	if(t == rr->t){
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> copy"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return clone();
	}else{
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> empty"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return new emptyscen;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection transpositions - unordered .. implemented in unordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *transp::intersect(const unordered *rr, bool prefix) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection transpositions - ordered .. implemented in ordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *transp::intersect(const ordered *rr, bool prefix) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *transp::mkealt() {
	alternative *alt = new alternative();
	unordered *uor = NULL;
	vector< set<int> > tmp;

	for( set<set<int> >::iterator it = t.begin(); it != t.end(); it++ ){
		tmp.push_back(*it);
	}
	tmp.push_back( set<int>() );
	set_union(tmp[0].begin(), tmp[0].end(), tmp[1].begin(), tmp[1].end(), inserter(tmp[2], tmp[2].begin()));

		// get alternative reversals
	uor = new unordered();
	for(unsigned i=0; i<tmp.size(); i++)
		uor->insert( new rev(tmp[i], _nmap) );
	alt->insert(uor);

		// get alternative reverse transposition + reversal scenarios
	uor = new unordered();
	uor->insert( new revtransp( tmp[0], tmp[1], _nmap ) );
	uor->insert( new rev( tmp[0], _nmap ) );
	alt->insert(uor);

	uor = new unordered();
	uor->insert( new revtransp( tmp[1], tmp[0], _nmap ) );
	uor->insert( new rev( tmp[1], _nmap ) );
	alt->insert(uor);

		// and finally the transposition itself is also an 'alternative'
	alt->insert( clone() );

	return alt;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//bool transp::operator==(const rrrmt *rrrmtp) const{
//	if( const transp* transpp = dynamic_cast<const transp*>(rrrmtp) ){
//		return t == transpp->t ;
//	}
//	return false;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool transp::operator<(const rrrmt *rrrmtp) const{
	if( type() != rrrmtp->type() ){
		return type() < rrrmtp->type();
	}
	const transp* transpp = dynamic_cast<const transp*>(rrrmtp);
	return t < transpp->t;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// output foo
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &transp::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) out << "\t";
	}
	out << typestrg(d);
	if( d < 2 ){
		out <<pquot<<"(";
		for(set<set<int> >::iterator it = t.begin(); it != t.end(); it++){
			for( set<int>::iterator jt = it->begin(); jt!= it->end(); jt++ ){
				if( _nmap == NULL ){
					out << *jt<<" ";
				}else{
					if( *jt < 0){
						out << "-";
					}
//					out << abs(*jt)<<":";
					out << (*_nmap)[abs(*jt)]<<" ";
				}
			}
			out << ",";
		}
		out<<pquot << ")";
		if( d == 0 ){
			out <<endl;
		}
	}
	return out;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void transp::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	set<int> tmp;
	set<set<int> > nt;

	for(set<set<int> >::iterator it = t.begin(); it!=t.end(); it++){
		for( set<int>::iterator jt = it->begin(); jt!=it->end(); jt++ ){
			for(unsigned i=0; i<mapping[*jt].size(); i++){
				tmp.insert(abs(mapping[*jt][i]));
			}
//			tmp.insert( mapping[*jt].begin(), mapping[*jt].end());
		}
		nt.insert(tmp);
		tmp.clear();
	}
	t = nt;
	_nmap = nmap;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// get the type identifier
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int transp::type() const{
	return TRA;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string transp::typestrg(unsigned d) const{
	if(d > 0){
		return TRANM;
	}else{
		return "transposition";
	}
}






// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

revtransp::revtransp(int i, int j, int k, const genom &g){
	set<int> tmp;

	if ( i < 0 || j < 0 || k < 0 ){
		cerr << "cannot construct inverse transposition with negative indices: "<<i<<","<<j<<","<<k<<endl;
		exit(1);
	}

	if( i>j || !(k<i || k>j) ){
		cerr << "cannot construct inverse transposition with unsorted indices: "<<i<<","<<j<<endl;
		exit(1);
	}

	if( k > (int)g.size() || j >= (int)g.size() ){
		cerr << "cannot construct inverse transposition exeeding the genome boundaries: ("<<g.size()<<"): "<<i<<","<<j<<","<<k<<endl;
		exit(1);
	}

	for(int x=i; x<=j; x++){
		tmp.insert( abs(g[x]) );
	}
	rt.push_back(tmp);
	tmp.clear();

	for(int x=((k<i)?k:j+1); x<=((k<i)?i-1:k-1); x++){
		tmp.insert( abs(g[x]) );
	}
	rt.push_back(tmp);
	tmp.clear();

	_nmap = g.get_nmap();
}


revtransp::revtransp(int i, int j, int k, int l, const genom &g){
	set<int> tmp;

	if ( i < 0 || j < 0 || k < 0 || l < 0 ){
		cerr << "cannot construct inverse transposition with negative indices: "<<i<<","<<j<<","<<k<<","<<l<<endl;
		exit(1);
	}

	if( !( i<=j && k<=l ) ){
		cerr << "cannot construct inverse transposition with unsorted indices: "<<i<<","<<j<<","<<k<<","<<l<<endl;
		exit(1);
	}

	if( l >= (int)g.size() || j >= (int)g.size() ){
		cerr << "cannot construct inverse transposition exeeding the genome boundaries: ("<<g.size()<<"): "<<i<<","<<j<<","<<k<<","<<l<<endl;
		exit(1);
	}

	for(int x=i; x<=j; x++){
		tmp.insert( abs(g[x]) );
	}
	rt.push_back(tmp);
	tmp.clear();

	for(int x=k; x<=l; x++){
		tmp.insert( abs(g[x]) );
	}
	rt.push_back(tmp);
	tmp.clear();

	_nmap = g.get_nmap();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

revtransp::revtransp( const genom &g, unsigned ml ){

	vector<unsigned> iv;
	set<int> tmp;

	iv = rng_inc_seq(3, g.size(), ml);

	for(unsigned i=iv[0]; i<iv[1]; i++){
		tmp.insert( abs(g[i]) );
	}
	rt.push_back(tmp);
	tmp.clear();

	for(unsigned i=iv[1]; i<=iv[2]; i++){
		tmp.insert( abs(g[i]) );
	}
	rt.push_back(tmp);
	tmp.clear();

	if( ask_rng() % 2 == 1){
		swap( rt[0], rt[1] );
	}

	_nmap = g.get_nmap();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

revtransp::revtransp( const set<int> &rtra, const set<int> &tra, vector<string> *nm ){
	rt.push_back(rtra);
	rt.push_back(tra);

	_nmap = nm;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

revtransp::~revtransp(){
	rt.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a reverse transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void revtransp::apply(genom &g) const{
#ifdef DEBUG_APPLY
	cerr << "apply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_APPLY
	reverse(g, rt[0]);
	transpose(g, rt);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

revtransp* revtransp::clone() const{
	return new revtransp(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a reverse transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void revtransp::deapply(genom &g) const{
	reverse(g, rt[0]);
	transpose(g, rt);
#ifdef DEBUG_DEAPPLY
	cerr << "deapply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_DEAPPLY
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void revtransp::elements( vector<set<int> > &elms, bool ab ) const{
//	vector<int> e;
	for( unsigned i=0;i<rt.size(); i++ ){
		elms.push_back(set<int>());
		for( set<int>::const_iterator it = rt[i].begin(); it!=rt[i].end(); it++ ){
			elms.back().insert( ab?abs(*it):*it );
		}
	}
//	elms.push_back(e);
}

bool revtransp::equal( const rrrmt *r ) const{
	return r->equal(this);
}

bool revtransp::equal( const revtransp *r ) const{
	return ( this->rt == r->rt );
}

void revtransp::getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode ){
	revtransp *r = this->clone();
	if( ! out.insert( r ).second ){
		delete r;
	}

}
void revtransp::getrrrmt( vector<rrrmt*> &out, int mode ){
	out.push_back( clone() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for reverse transpositions
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *revtransp::intersect(const rrrmt *rrrmtp, bool prefix) const{
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection reverse transposition - reverse transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *revtransp::intersect(const revtransp *rr, bool prefix) const{
	#ifdef DEBUG_ATOMINTERSECTION
	cout << "revtransp::intersect "<<*this<< " vs. "<<*rr;
	#endif//DEBUG_ATOMINTERSECTION
	if(rt == rr->rt){
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> copy"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return clone();
	}else{
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> empty"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return new emptyscen;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection reverse transposition - unordered .. implemented in unordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *revtransp::intersect(const unordered *rr, bool prefix) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection reverse transposition - ordered .. implemented in ordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *revtransp::intersect(const ordered *rr, bool prefix) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *revtransp::mkealt(){
	alternative *alt = new alternative();
	unordered *uor = NULL;
	set<int> tmp;
		// add the transposition + reversal scenarios
	uor = new unordered();
	uor->insert( new transp(rt, _nmap) );
	uor->insert( new rev(rt[0], _nmap) );
	alt->insert(uor);

		// add the reversal + reversal scenarios
	uor = new unordered();
	set_union(rt[0].begin(), rt[0].end(), rt[1].begin(), rt[1].end(), inserter(tmp, tmp.begin()));
	uor->insert( new rev(tmp, _nmap) );
	uor->insert( new rev(rt[1], _nmap) );
	alt->insert(uor);

		// add the reverse transposition itsef
	alt->insert( clone() );
	return alt;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//bool revtransp::operator==(const rrrmt *rrrmtp) const{
//	if( const revtransp* revtranspp = dynamic_cast<const revtransp*>(rrrmtp) ){
//		return rt == revtranspp->rt ;
//	}
//	return false;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool revtransp::operator<(const rrrmt *rrrmtp) const{
	if( type() != rrrmtp->type() ){
		return type() < rrrmtp->type();
	}
	const revtransp* revtranspp = dynamic_cast<const revtransp*>(rrrmtp);
	return rt < revtranspp->rt ;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &revtransp::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) out << "\t";
	}
	out << typestrg(d);

	if( d < 2 ){
		out<<pquot <<"(";
		for(unsigned i=0; i<rt.size(); i++){
			for( set<int>::iterator it = rt[i].begin(); it != rt[i].end(); it++ ){
				if( _nmap != NULL ) {
					if( *it < 0){
						out << "-";
					}
					out << (*_nmap)[ abs(*it) ]<<" ";
				}else{
					out << *it<<" ";
				}
			}
			out << ",";
		}
		out <<pquot<< ")";
		if (d==0){
			out <<endl;
		}
	}
	return out;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void revtransp::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	vector<set<int> > nrt(2);

	for(unsigned i=0; i<rt.size(); i++){
		for( set<int>::iterator jt = rt[i].begin(); jt!=rt[i].end(); jt++ ){
			for(unsigned j=0; j<mapping[*jt].size(); j++){
				nrt[i].insert( abs(mapping[*jt][j]) );
			}
//			nrt[i].insert( mapping[*jt].begin(), mapping[*jt].end());
		}
	}

	rt = nrt;
	_nmap = nmap;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int revtransp::type() const{
	return RTRA;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string revtransp::typestrg(unsigned d) const{
	if(d > 0){
		return RTRANM;
	}else{
		return "inverse transposition";
	}
}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************

void tdrl::_init(const vector<bool> &tdl, const genom &g){
	int s = 0, e = 0;
	unsigned blockidx=0;

	_t = vector<set<int> >(2, set<int>());
	_orgord = vector<int>( g.size()+1, 0 );

//	cout <<"tdrl::tdrl (";copy(tdl.begin(), tdl.end(), ostream_iterator<bool>(cout,"") ); cout<<", "<<g<<")" <<endl;
//	cout << g <<endl;

	// get normal form, i.e. get start and end of the interval which has to be copied
	// and check if tdrl is a "real tdrl", i.e. not a transposition
	tdl_normalform( tdl, s, e );
	realtdl = is_tdl( tdl, s, e);
//	cout << "s "<<s<<" e "<<e<<" t "<<realtdl<<endl;

	// insert the elements of the genome in first / second copy
	// and determine the blocks, each max. interval which is in the same
	// copy gets the same index, indices are counted from left to right,
	// elements not in the copy get index 0
 	for(int i=s; i<e; i++){
		if( tdl[i] ){
			_t[0].insert( abs(g[i]) );
		}else{
			_t[1].insert( abs(g[i]) );
		}
		if( i==s || tdl[i] != tdl[i-1]){
			blockidx++;
		}
		_orgord[ abs(g[i]) ] = blockidx;
	}
 	_nmap = g.get_nmap();
// 	cout << "tdrl::tdrl -> "<<*this<<endl;
//	cout << "orgord "; copy(orgord.begin(), orgord.end(), ostream_iterator<int>(cout, " ")); cout << endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tdrl::tdrl(const vector<bool> &tdl, const genom &g){

	_init( tdl, g );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tdrl::tdrl(const vector<set<int> > &set, vector<string> *nm, bool rtdrl){
	_t=set;
	_nmap=nm;
	realtdl=rtdrl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// contruct a random tdrl applied to a given genome
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tdrl::tdrl(const genom &g, bool randrange, unsigned ml) {
	vector<unsigned> iv;	// interval
	vector<bool> rt;

	if( ml == 0 ){
		cerr << "tdrl::tdrl ml = 0"<<endl;
		exit(EXIT_FAILURE);
	}

	if( randrange ){
		iv = rng_inc_seq( 2, g.size(), ml );
		iv[1]++;
	}else{
		iv = vector<unsigned>( 2, 0 );
		iv[1] = g.size();
	}

	rt = vector<bool>(g.size());
	for(unsigned i=0; i<g.size(); i++){
		if( i<iv[0] ){			// keep all < s in the 1st copy
			rt[i] = true;
		}else if( i>=iv[1] ){	// keep all > e in the 2nd copy
			rt[i] = false;
		}else{				// otherwise choose randomly
			rt[i] = ask_rng(0,1);
		}
	}

//	cout << "s "<<iv[0]<<" e "<<iv[1]<<endl;
//	cout << " cand "; copy(rt.begin(), rt.end(), ostream_iterator<bool>(cout,"")); cout << endl;

	_init( rt, g );

	////	cout << "tdrl::tdrl random"<<endl;
//	realtdl = true;
//	do{
//		rt = random_tdl( g.size() );
////		cout << " cand "; copy(rt.begin(), rt.end(), ostream_iterator<bool>(cout,"")); cout << endl;
////		cout << is_tdl(rt)<<endl;;
//	}while( !is_tdl( rt ) );
//	cout << "  take this "<<endl;
//	cout <<"tdrl::tdrl (";copy(rt.begin(), rt.end(), ostream_iterator<bool>(cout,"") );  cout<<", "<<g<<")" <<endl;

//	unsigned s = std::numeric_limits< unsigned >::max(),
//		e = std::numeric_limits< unsigned >::max();
//	tdl_normalform( rt, s, e );
//	realtdl = is_tdl( rt, s, e );
//
////	realtdl = is_tdl( rt );
////	s=0;
////	e=g.size()-1;
//
//	for(unsigned i=s; i<=e; i++){
//		if( rt[i] ){
//			t[0].insert( abs(g[i]) );
//		}else{
//			t[1].insert( abs(g[i]) );
//		}
//		if( i==s || rt[i] != rt[i-1]){
//			blockidx++;
//		}
//		orgord[abs(g[i])] = blockidx;
//	}
//
//	nmap = g.get_nmap();
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tdrl::~tdrl(){
	_t.clear();
	_orgord.clear();
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a tdrl
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl::apply(genom &g) const{
	genom tmp;
	unsigned f=0, 	// index of the first
		s=1;		// and second copy
#ifdef DEBUG_APPLY
	cerr << "apply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_APPLY

	tmp.set_nmap( g.get_nmap() );


	// first check if the tdrl is applied to an interval
	unsigned bc = 0;
	for( unsigned i=0; i<g.size()-1; i++ ){
		if( (_t[f].find(abs(g[i])) == _t[f].end() && _t[s].find(abs(g[i])) == _t[s].end())
			!= (_t[f].find(abs(g[i+1])) == _t[f].end() && _t[s].find(abs(g[i+1])) == _t[s].end()) ){

			bc++;
		}
	}
	if( bc > 2 ){
		set<int> tmp;
		set_union(_t[0].begin(), _t[0].end(), _t[1].begin(), _t[1].end(), inserter(tmp, tmp.begin()));
		throw NoIntervalException( g, tmp);
//		cerr << "tdrl::apply(): tdrl is no interval "<<bc<<endl;
//		for(unsigned j=0; j<g.size(); j++){
//			if( _t[f].find(abs(g[j])) != _t[f].end() ){cout << "<";}
//			if( _t[s].find(abs(g[j])) != _t[s].end() ){cout << ">";}
//			cout << g[j]<<" ";
//		}
		cout << endl;
		exit(EXIT_FAILURE);
	}

	// take all elements which are not in the duplication (i.e. not in t[0] or t[1])
	// until an element is found which is in the  duplication
	unsigned i=0, j=0, m=0;
	for(j=0; j<g.size(); j++){
		if( _t[f].find(abs(g[j])) == _t[f].end() && _t[s].find(abs(g[j])) == _t[s].end() ){
			tmp.push_back( g[j] );
		}else{
			break;
		}
	}

	m = j;


	// if the stored original order does not match the first element of the
	// genome that is affected by the tdrl then it can be assumed that there
	// was a reversal that included the elements of the tdrl, thus the
	// roles of 1st and 2nd copy have to be changed
	if( j<g.size() &&  _orgord[abs(g[j])] != 1 ){
		//cout << "swap" << endl;
		//cout << "swap "<<abs(g[i])<<" : "<< orgord[abs(g[i])]<<endl;
		swap( f, s );
	}

		// take all elements in the first copy
	for( i=j; i<g.size(); i++){
		if( _t[f].find(abs(g[i])) != _t[f].end() ){
			tmp.push_back( g[i] );
			m = max(m,i);
		}
	}
		// take all elements in the second copy
	for( i=j; i<g.size(); i++){
		if( _t[s].find(abs(g[i])) != _t[s].end() ){
			tmp.push_back( g[i] );
			m = max(m,i);
		}
	}
		// take the remaining elements
	for(unsigned i=m+1; i<g.size(); i++){
		tmp.push_back( g[i] );
	}
	g = tmp;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

tdrl* tdrl::clone() const{
	return new tdrl(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a tdrl
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl::deapply(genom &g) const{
//	cout << "deapply "<< *this << endl;
	genom oldg,
		oldginv;
	int i=0,
		maxidx=0;
	unsigned d, invd;
	vector<vector<int> > ord,
		invord;

#ifdef DEBUG_DEAPPLY
	cerr << "deapply "<<*this;
	cerr << "orgord "; copy(_orgord.begin(), _orgord.end(), ostream_iterator<int>(cout, " ")); cout << "size: "<<_orgord.size()<< endl;
#endif//DEBUG_DEAPPLY

	// deapply a tdrl.
	// - make k decks of cards where k is the max value found in orgord
	// - all elements in the prefix of g that have orgord=0 are just appended
	//   to the new genome
	// - then the elements with orgord != 0 are pushed to the corresponding stack.
	// - the elements from the stacks are appended to the new genome
	// - the remaining elements are appended to the genome

	oldg.set_nmap(g.get_nmap());
	oldginv.set_nmap(g.get_nmap());

	// get the maximal block
	for(unsigned i=0; i<_orgord.size(); i++){
		maxidx = max( _orgord[i], maxidx );
	}

	// initialize stacks
	ord = vector<vector<int> >( maxidx );
	invord = vector<vector<int> >( maxidx );

	// append elements with orgord==0
	for( i=0; i<(int)g.size(); i++ ){
		if(_orgord[abs(g[i])] == 0){
			oldg.push_back(g[i]);
			oldginv.push_back(g[i]);
		}else{
			break;
		}
	}

	for( ; i<(int)g.size(); i++){
		if( _orgord[abs(g[i])]==0 ){
			break;
		}
		ord[ _orgord[ abs(g[i]) ]-1 ].push_back( g[i] );
		invord[ maxidx-_orgord[ abs(g[i]) ] ].push_back( g[i] );
	}

	for(unsigned j=0; j<ord.size(); j++){
		for( vector<int>::iterator it=ord[j].begin(); it!=ord[j].end(); it++){
			oldg.push_back( *it );
		}
		for( vector<int>::iterator it=invord[j].begin(); it!=invord[j].end(); it++){
			oldginv.push_back( *it );
		}
	}
//	cout << i <<" "<< g[i]<<endl;
	for( ; i<(int)g.size(); i++){
		if( _orgord[abs(g[i])] != 0 ){
			cerr << "internal error: tdrl::deapply "<<_orgord[abs(g[i])]<<" != 0 "<<endl;
			cerr << g<<endl;
			cerr << "tdrl::deapply "<<*this;
			cerr << "tdrl::deapply orgord ";copy(_orgord.begin(), _orgord.end(), ostream_iterator<int>(cout, " "));cout << endl;
			exit(EXIT_FAILURE);
		}
		oldg.push_back(g[i]);
		oldginv.push_back(g[i]);
	}


	d = tdrl_distance(oldg, g);
	invd = tdrl_distance(oldginv, g);

#ifdef DEBUG_DEAPPLY
	cerr << "oldg "<<oldg<<" d "<<d<<endl;
	cerr << "oldg b) "<<oldginv<<" d "<<invd<<endl;
#endif//DEBUG_DEAPPLY

	if( min(d, invd) != 1 ){
		cout << "internal error: tdrl::deapply: could not reconstruct tdrl with distance 1"<<endl;
		cout << *this;
		exit(1);
	}

	if( d == 1 && invd != 1 ){
		g = oldg;
	}else if( d!= 1 && invd == 1 ){
		g = oldginv;
	}else{
		cerr << "internal error: tdrl::deapply equal distance"<<endl;
		exit(1);
	}
#ifdef DEBUG_DEAPPLY
	cerr << "deapply "<< *this << endl;
	cerr << g<< endl;
#endif//DEBUG_DEAPPLY
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl::elements( vector<set<int> > &elms, bool ab) const{
//	set<int> e;
	for( unsigned i=0;i<_t.size(); i++ ){
		elms.push_back(set<int>());
		for( set<int>::const_iterator it = _t[i].begin(); it!=_t[i].end(); it++ ){
			elms.back().insert( ab?abs(*it):*it );
		}
	}
//	elms.push_back(e);
}

bool tdrl::equal( const rrrmt *r ) const{
	return r->equal(this);
}

bool tdrl::equal( const tdrl *r ) const{
	return (_t == r->_t);
}

void tdrl::getrrrmt( set<rrrmt*, HDereferenceLess> &out, int mode ){
	tdrl *t = this->clone();
	if( ! out.insert( t ).second ){
		delete t;
	}
}
void tdrl::getrrrmt( vector<rrrmt*> &out, int mode ){
	out.push_back( clone() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for tdrls
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *tdrl::intersect(const rrrmt *rrrmtp, bool prefix) const{
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection tdrl - tdrl
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *tdrl::intersect(const tdrl *rr, bool prefix) const{
	#ifdef DEBUG_ATOMINTERSECTION
	cout << "tdrl::intersect "<<*this<< " vs. "<<*rr;
	#endif//DEBUG_ATOMINTERSECTION
	if(_t == rr->_t){
		// @todo handling of orgord ??
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> copy"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return clone();
	}else{
		#ifdef DEBUG_ATOMINTERSECTION
		cout << " -> empty"<<endl;
		#endif//DEBUG_ATOMINTERSECTION
		return new emptyscen;
	}
	// @todo handling of orgord ??
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection tdrl - unordered .. implemented in unordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *tdrl::intersect(const unordered *rr, bool prefix) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection tdrl - ordered .. implemented in ordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *tdrl::intersect(const ordered *rr, bool prefix) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//bool tdrl::operator==(const rrrmt *rrrmtp) const{
//	if( const tdrl* tdrlp = dynamic_cast<const tdrl*>(rrrmtp) ){
//		return t == tdrlp->t ;
//	}
//	return false;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool tdrl::operator<(const rrrmt *rrrmtp) const{
	if( type() != rrrmtp->type() ){
		return type() < rrrmtp->type();
	}
	const tdrl* tdrlp = dynamic_cast<const tdrl*>(rrrmtp);
	return _t < tdrlp->_t ;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &tdrl::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) out << "\t";
	}
	out << typestrg(d);
	if( d < 2 ){
		out<<pquot<<"(";
		for(unsigned i=0; i<_t.size(); i++){
			for( set<int>::iterator it = _t[i].begin(); it!=_t[i].end(); it++ ){
				if( _nmap == NULL ){
					out << *it<<" ";
				}else{
					if( *it < 0 ){
						out << "-";
					}
					out <<(*_nmap)[abs(*it)]<<" ";
				}
			}
			out << ",";
		}
		out<<pquot<<")";
		if(d==0){
			out<<endl;
		}
	}
	return out;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	vector<set<int> > nt(2);
	vector<int> noo;

//	cerr << "tdrl::rename("<<*this<<")"<<endl;
//	cerr << "mapping "<<endl;
//	for(unsigned i=0; i<mapping.size(); i++){
//		cerr<<i<<" -> "; copy(mapping[i].begin(), mapping[i].end(), ostream_iterator<int>(cerr," ")); cerr<<endl;
//	}
//	copy( _orgord.begin(), _orgord.end(), ostream_iterator<int>(cerr, " ") ); cerr<<endl;

	for(unsigned i=0; i<_t.size(); i++){
		for( set<int>::iterator jt = _t[i].begin(); jt!=_t[i].end(); jt++ ){
			for(unsigned j=0; j<mapping[*jt].size(); j++){
				nt[i].insert( abs(mapping[*jt][j]) );
			}
		}
	}
	_t = nt;

//	cerr << "max "<<mx<<endl;
	noo = vector<int>(mx+1, 0);
	for( unsigned i=0; i<_orgord.size(); i++ ){
		if(_orgord[i] != 0){
			for( unsigned j=0; j<mapping[i].size(); j++ ){
				noo[abs(mapping[i][j])] = _orgord[i];
			}
		}
	}
	_orgord = noo;
	_nmap = nmap;
//	cerr << "tdrl::rename ->"<<*this<<endl;
//	copy( _orgord.begin(), _orgord.end(), ostream_iterator<int>(cerr, " ") ); cerr<<endl;
}

void tdrl::set_orgord(vector<int > &bla){
	_orgord=bla;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *tdrl::simplify(){
	if( _t[0].size() == 0 && _t[1].size() == 0 ){
		return new emptyscen();
	}else if( !realtdl ){
		return new transp(_t, _nmap);
	}else{
		return clone();
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int tdrl::type() const{
	return TDRL;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string tdrl::typestrg(unsigned d) const{
	if( realtdl ){
		//return TDRLNM;
		return "tdrl";
	}else{
		return "tdrl";
	}
}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unordered::unordered(){
	_complete = true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unordered::unordered( const unordered &c ){
	_complete = c._complete;
	for(set<rrrmt*, HDereferenceLess>::iterator it=c._sce.begin(); it != c._sce.end(); it++){
		insert( (*it)->clone() );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unordered::unordered( const vector<rrrmt*> &c ){
	_complete = true;
	for(vector<rrrmt*>::const_iterator it=c.begin(); it != c.end(); it++){
		insert( (*it)->clone() );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// construct a random scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unordered::~unordered(){
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		delete *it;
	}
	_sce.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unordered::append( vector<vector<rrrmt*> > &v ) const{
	for(set<rrrmt*, HDereferenceLess>::const_iterator it=_sce.begin(); it!=_sce.end(); it++){
		(*it)->append(v);
//		v.push_back( vector<rrrmt*>((unsigned)1, (*it)->clone()) );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a ordered scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unordered::apply(genom &g) const{
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		(*it)->apply(g);
	}
}

void unordered::apply( genom &g, set<genom> &gs )const{

	vector<rrrmt*> sceperm;
	sceperm.insert(sceperm.begin(), _sce.begin(), _sce.end());
	do{
		genom tmp = g;
		for( unsigned i=0; i<sceperm.size(); i++){
				sceperm[i]->apply( tmp, gs );
		}
	}while( next_permutation( sceperm.begin(), sceperm.end() ));

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unordered* unordered::clone() const{
	return new unordered(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unordered* unordered::create() const{
//	cout << "unordered::create()"<<endl;
	return new unordered();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a unordered scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unordered::deapply(genom &g) const{
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		if( (*it)->iscomplete() ){
			(*it)->deapply(g);
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool unordered::equal( const rrrmt *r ) const{
	return r->equal(this);
}

bool unordered::equal( const unordered *r ) const{
	if( _sce.size() != r->_sce.size() ){
		return false;
	}

	set<rrrmt*, DereferenceLess>::iterator it = _sce.begin(),
			jt = r->_sce.begin();

	while( it!=_sce.end() && jt!=r->_sce.end() ){
		if( ! (*it)->equal( *jt ) ){
			return false;
		}
	}

	return true;
}


void unordered::elements( vector<set<int> > &elms, bool ab ) const{
	for( set<rrrmt*, DereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		(*it)->elements(elms, ab);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unordered::getrrrmt( set<rrrmt*,HDereferenceLess> &rset, int mode ) {
	for( set<rrrmt*,HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		(*it)->getrrrmt(rset, mode);
	}
}
void unordered::getrrrmt( vector<rrrmt*> &rset, int mode ) {
	for( set<rrrmt*,HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		(*it)->getrrrmt(rset, mode);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void unordered::getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const{
	for( set<rrrmt*,HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		if((*it)->type()==1 || (*it)->type()==2 || (*it)->type()== 3 || (*it)->type()==7 || (*it)->type()==8){
				scen.push_back((*it)->clone());
		}else{
				(*it)->getminrrrmt(scen, cst);
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool unordered::hascomposite() const{
	for( set<rrrmt*,HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		if( (*it)->type() >= UORD ){
			return true;
		}
	}
	return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool unordered::insert(rrrmt *r){
	bool suc = false;

	// try to insert the element into the scenario
	if( !r->isempty() ){
		suc = (_sce.insert(r)).second;
	}
	// if its already inside -> free memory
	if( !suc ){
		delete r;
	}
	return suc;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for unordered scenarios
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *unordered::intersect( const rrrmt *rrrmtp, bool prefix) const{
//#ifdef DEBUG_INTERSECTION
//cerr << "start   unordered::intersect (***)"; output(cerr, 0, 1); cerr<< " vs. "; rrrmtp->output(cerr, 0, 1); cerr<<endl;
//#endif//DEBUG_INTERSECTION
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - indel
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt* unordered::intersect( const indel *rr, bool prefix ) const{
//#ifdef DEBUG_INTERSECTION
//cerr << "start   unordered::intersect (ind)"; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//#endif//DEBUG_INTERSECTION
	rrrmt *tmp = NULL;
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		tmp = (*it)->intersect( rr, prefix );
		if( !tmp->isempty() ){
//			#ifdef DEBUG_INTERSECTION
//			cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr <<endl;
//			cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
//			#endif//DEBUG_INTERSECTION
			return tmp;
		}else{
			delete tmp;
		}
	}
//	#ifdef DEBUG_INTERSECTION
//	cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr <<endl;
//	cerr << "-> empty"<<endl;
//	#endif//DEBUG_INTERSECTION
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - reversal
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt* unordered::intersect( const rev *rr, bool prefix ) const{
//#ifdef DEBUG_INTERSECTION
//cerr << "start   unordered::intersect (rev)"; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//#endif//DEBUG_INTERSECTION

	rrrmt *tmp = NULL;
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		tmp = (*it)->intersect( rr, prefix );
		if( !tmp->isempty() ){
//			#ifdef DEBUG_INTERSECTION
//			cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//			cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
//			#endif//DEBUG_INTERSECTION
			return tmp;
		}else{
			delete tmp;
		}
	}
//	#ifdef DEBUG_INTERSECTION
//	cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//	cerr << "-> empty"<<endl;
//	#endif//DEBUG_INTERSECTION
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt* unordered::intersect( const transp *rr, bool prefix ) const{
//#ifdef DEBUG_INTERSECTION
//cerr << "start   unordered::intersect (tra) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//#endif//DEBUG_INTERSECTION
	rrrmt *tmp = NULL;
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		tmp = (*it)->intersect( rr, prefix );
		if( !tmp->isempty() ){
//			#ifdef DEBUG_INTERSECTION
//			cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//			cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
//			#endif//DEBUG_INTERSECTION
			return tmp;
		}else{
			delete tmp;
		}
	}
//	#ifdef DEBUG_INTERSECTION
//	cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//	cerr << "-> empty"<<endl;
//	#endif//DEBUG_INTERSECTION
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - reverse transposition
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt* unordered::intersect( const revtransp *rr, bool prefix ) const{
//#ifdef DEBUG_INTERSECTION
//cerr << "start   unordered::intersect (rtr) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//#endif//DEBUG_INTERSECTION

	rrrmt *tmp = NULL;
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		tmp = (*it)->intersect( rr, prefix );
		if( !tmp->isempty() ){
//			#ifdef DEBUG_INTERSECTION
//			cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//			cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
//			#endif//DEBUG_INTERSECTION
			return tmp;
		}else{
			delete tmp;
		}
	}
//	#ifdef DEBUG_INTERSECTION
//	cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//	cerr << "-> empty"<<endl;
//	#endif//DEBUG_INTERSECTION
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - tdrl
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt* unordered::intersect( const tdrl *rr, bool prefix ) const{
//#ifdef DEBUG_INTERSECTION
//cerr << "start   unordered::intersect (tdl) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//#endif//DEBUG_INTERSECTION


	rrrmt *tmp = NULL;
	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		tmp = (*it)->intersect( rr, prefix );
		if( !tmp->isempty() ){
//			#ifdef DEBUG_INTERSECTION
//			cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//			cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
//			#endif//DEBUG_INTERSECTION
			return tmp;
		}else{
			delete tmp;
		}
	}
//	#ifdef DEBUG_INTERSECTION
//	cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
//	cerr << "-> empty"<<endl;
//	#endif//DEBUG_INTERSECTION
	return new emptyscen();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - ordered scenario .. implemented in ordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *unordered::intersect( const ordered *rr, bool prefix ) const{
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - alternative scenario .. implemented in alternative
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *unordered::intersect( const alternative *rr, bool prefix ) const{
#ifdef DEBUG_INTERSECTION
cerr << "start   unordered::intersect (alt) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
#endif//DEBUG_INTERSECTION
	return rr->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// intersection unordered - unordered
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *unordered::intersect( const unordered *rr, bool prefix ) const{
	rrrmt *tmp, *ret;
	set<rrrmt*,HDereferenceLess> its;// set of potential intersections
	vector<vector<rrrmt*> > incpy; 	// temporary copy of the unordered scenarios
	unordered *tempu;				// stores the result of all-vs-all intersection
	unsigned max = 0;

#ifdef DEBUG_INTERSECTION
cerr << "start   unordered::intersect (uno) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
#endif//DEBUG_INTERSECTION

	// copy input scenarios
	incpy = vector<vector<rrrmt*> >( 2 );
	incpy[0].insert(incpy[0].begin(), _sce.begin(), _sce.end());
	incpy[1].insert(incpy[1].begin(), rr->_sce.begin(), rr->_sce.end());

	// intersection of the complete scenario with every element of the other
	// if non-empty add to its
	for( set<rrrmt*,DereferenceLess>::const_iterator it=this->_sce.begin(); it!=this->_sce.end(); it++ ){
		tmp = (*it)->intersect( rr, prefix );
		if( ! its.insert( tmp ).second )
			delete tmp;
	}
	for( set<rrrmt*,DereferenceLess>::const_iterator it=rr->_sce.begin(); it!=rr->_sce.end(); it++ ){
		tmp = (*it)->intersect( this, prefix );
		if( ! its.insert( tmp ).second )
			delete tmp;
	}

	// intersection of every vs every
	tempu = create();
	for( set<rrrmt*,DereferenceLess>::const_iterator it=this->_sce.begin(); it!=this->_sce.end(); it++ ){
		for( set<rrrmt*,DereferenceLess>::const_iterator jt=rr->_sce.begin(); jt!=rr->_sce.end(); jt++ ){
			tmp = (*it)->intersect( *jt, prefix);
			tempu->insert( tmp );
		}
	}
	tmp = tempu->simplify();
	delete tempu;
	if( ! its.insert( tmp ).second )
		delete tmp;

	// determine the max of its[0]...its[2]
	max = getmax(its);
	#ifdef DEBUG_INTERSECTION
	cerr << "unordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION
	if( its.size() == 0 ||  max == 0 ){
		tmp = new emptyscen;
	}else if( its.size() == 1){
		tmp = (*its.begin())->clone();
	}else{
		alternative *atmp = new alternative();
		for( set<rrrmt*, DereferenceLess>::const_iterator it=_sce.begin(); it!=_sce.end(); it++){
			atmp->insert( (*it)->simplify() );
		}
		tmp = atmp;
	}
	#ifdef DEBUG_INTERSECTION
	cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION

	for( set<rrrmt*,HDereferenceLess>::iterator it=its.begin(); it!=its.end(); it++ ){
		delete *it;
	}
	its.clear();

	// simplify result if possible
	ret = tmp->simplify();
	delete tmp;

	if(ret->isempty()){
		delete ret;
		ret = new emptyscen();
	}
	return ret;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool unordered::iscomplete() const{
	if(!_complete){
		return false;
	}
	for(set<rrrmt*,HDereferenceLess>::const_iterator it=_sce.begin(); it!=_sce.end(); it++){
		if( !(*it)->iscomplete() ){
			return false;
		}
	}
	return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool unordered::isempty() const{
	return (_sce.size() == 0);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned unordered::length( int mode ) const{
	unsigned s = 0;
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		s += (*it)->length( mode );
	}
	return s;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *unordered::mkealt(){
	rrrmt *tmp;
	set<rrrmt*, HDereferenceLess> update;

	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end();){
		if( !(*it)->iscomplete() ){
			it++;
			continue;
		}
		tmp = (*it)->mkealt();
		if( !tmp->isempty() ){
			delete *it;
			_sce.erase(it++);
			update.insert( tmp );
		}else{
			delete tmp;
			it++;
		}
	}
	_sce.insert( update.begin(), update.end() );
	return new emptyscen;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//rrrmt *unordered::mkepar() const{
//	unordered *ret = new unordered();
//	rrrmt *t;
//
//	for( set<rrrmt*,DereferenceLess>::const_iterator it=sce.begin(); it!=sce.end(); it++ ){
//		t = (*it)->mkepar();
//		if( !t->isempty() ){
//			ret->insert( t );
//		}else{
//			delete t;
//		}
//	}
//	return ret;
//}

unsigned unordered::nralt( ) const{
	unsigned a = 1;
	for(set<rrrmt*,DereferenceLess>::const_iterator it = _sce.begin(); it!=_sce.end(); it++){
		a *= (*it)->nralt();
	}
	return a;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

float unordered::getmincost(const costfoo *cst) const{
	float sum=0.0;
	for(set<rrrmt*,DereferenceLess>::const_iterator it = _sce.begin(); it!=_sce.end(); it++){
		sum += (*it)->getmincost(cst);
	}
	return sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool unordered::operator<(const rrrmt *rrrmtp) const{
	if( type() != rrrmtp->type() ){
		return type() < rrrmtp->type();
	}

	const unordered* uop = dynamic_cast<const unordered*>(rrrmtp);
	set<rrrmt *, HDereferenceLess>::iterator ita, itb;
	for( ita=_sce.begin(), itb=uop->_sce.begin(); ita!=_sce.end() && itb != uop->_sce.end(); ita++, itb++ ){
		if( *(*ita) < *itb ){
			return true;
		}else if( *(*itb) < *ita ){
			return false;
		}
	}
	return _sce.size() < uop->_sce.size();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &unordered::output(ostream &out, unsigned l, unsigned d, string pquot ) const{
	string indent = "";
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) indent+="\t";
	}

	out << indent;
	out << typestrg(d);
	out <<pquot<<"{";
	if (d==0){
		out<<endl;
	}
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		(*it)->output(out, l+1, d, pquot);
		out << ",";
	}
	out<<indent;
	out<<pquot<<"}";

	if( d==0 ){
		out<<" complete="<<_complete<<endl;
	}
	return out;
}

ostream &unordered::output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx ) {
	out<<"{";
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		if(it!=_sce.begin())
			out << ",";
		(*it)->output(out, ridx);
	}
	out<<"}";
	return out;
}

//void unordered::shortdesc(string &d) const{
//	d+="\\(";
//	for(set<rrrmt*, DereferenceLess>::iterator it=sce.begin(); it != sce.end(); it++){
//		(*it)->shortdesc(d);
//		d+=",";
//	}
//	d+="\\)";
//}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unordered::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	set<rrrmt *, HDereferenceLess> tmp;	// temporary copy of the rearrangements
	tmp = _sce;
	_sce.clear();

	for(set<rrrmt *, HDereferenceLess>::iterator it=tmp.begin(); it != tmp.end(); it++){
		(*it)->rename(mapping, nmap, mx);
		insert(*it);
	}
	tmp.clear();
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// replace all elements of the ordered scenario by a simplified version (if it exist)
// if the ordered scenario has length one return a copy of this one rearrangement
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *unordered::simplify(){
	rrrmt *tmp;
	unordered *ret = create();

	// - get a simplified version of each child
	// - lift the childs of unordered children one level up
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		tmp = (*it)->simplify();

		unordered *tmpu = dynamic_cast<unordered*>(tmp);
		if( tmpu && ret->type() == tmpu->type() ){
			for(set<rrrmt*, HDereferenceLess>::iterator jt=tmpu->_sce.begin(); jt != tmpu->_sce.end(); jt++){
				ret->insert( *jt );
			}
			tmpu->_sce.clear();
			delete tmpu;
		}else{
			ret->insert( tmp );
		}
	}
//	cerr << "simplify";this->output(cerr, 0, 1); cerr<< endl;

	if( ret->_sce.size() == 0 ){
		delete ret;
		tmp = new emptyscen();
//		cerr << "->"; tmp->output(cerr, 0, 1);cerr<<endl;
		return tmp;
	}else if( ret->_sce.size() == 1 ){
		tmp = (*ret->_sce.begin())->clone();
		delete ret;
//		cerr << "->"; tmp->output(cerr, 0, 1);cerr<<endl;
		return tmp;
	}else{
//		cerr << "->"; ret->output(cerr, 0, 1);cerr<<endl;
		return ret;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int unordered::type() const{
	return UORD;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string unordered::typestrg(unsigned d) const{
	if(d == 0){
		return "unordered";
	}else{
		return "";
	}
}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered::ordered(){
	_complete = true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered::ordered( const ordered &c ){
	_complete = c._complete;
	for(unsigned i=0; i<c.sce.size(); i++){
		sce.push_back( c.sce[i]->clone() );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered::ordered( const vector<rrrmt*> &c, bool comp ){
	_complete = comp;
	for(unsigned i=0; i<c.size(); i++){
		sce.push_back(c[i]->clone());
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// construct a random scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered::ordered( int k, vector<float> prob, genom &g, vector<genom> &trace, bool randrange, unsigned ml){
	int evtype;
	rrrmt *tmprandev = NULL;
	_complete = true;

	trace.clear();
	trace.push_back( g );
	weighted_choice_init(prob);
	for(int i=0; i<k; i++){
		evtype = weighted_choice( prob );
		switch ( evtype ){
			case 0: {	// reversal
				tmprandev = new rev(g, ml);
				break;
			}
			case 1: {	// transposition
				tmprandev = new transp(g, ml);
				break;
			}
			case 2: {	// reverse transposition
				tmprandev = new revtransp(g, ml);
				break;
			}
			case 3: {	// tdrl
				tmprandev = new tdrl(g, randrange, ml);
				break;
			}
			default:{
				cerr << "internal error: unknown event type choosen"<<endl;
			}
		}
//		cout << *tmprandev<<endl;
		tmprandev->apply(g);
		trace.push_back(g);
		sce.push_back( tmprandev );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered::~ordered(){
	for(unsigned i=0; i<sce.size(); i++){
		delete sce[i];
	}
	sce.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ordered::append( vector<vector<rrrmt*> > &v ) const{
	v.push_back( vector<rrrmt*>() );
	for( unsigned i=0; i<sce.size(); i++ ){
//	for( int i=sce.size()-1; i>=0; i-- ){
//	for(vector<rrrmt*>::reverse_const_iterator it=sce.rbegin(); it!=sce.rend(); it++){
		(v.back()).push_back( sce[i]->clone());
	}
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a ordered scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ordered::apply(genom &g) const{
	for(unsigned i=0; i<sce.size(); i++){
		sce[i]->apply(g);
	}
}

void ordered::apply( genom &g, set<genom> &gs )const{
	genom tmp = g;
	for( unsigned i=0; i<sce.size(); i++){
		sce[i]->apply( tmp, gs );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered* ordered::clone() const{
	return new ordered(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ordered* ordered::create() const{
	return new ordered();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply a ordered scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ordered::deapply(genom &g) const{
	for(int i=sce.size()-1; i>=0; i--){
		if(sce[i]->iscomplete())
			sce[i]->deapply(g);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void ordered::elements( vector<set<int> > &elms, bool ab ) const{
	for(unsigned i=0; i<sce.size(); i++){
		sce[i]->elements(elms, ab);
	}
}

bool ordered::equal( const rrrmt *r ) const{
	return r->equal(this);
}

bool ordered::equal( const ordered *r ) const{
	if( sce.size() != r->sce.size() ){
		return false;
	}
	for( unsigned i=0; i<sce.size(); i++ ){
		if( ! sce[i]->equal( r->sce[i] ) ){
			return false;
		}
	}
	return true;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ordered::getrrrmt( set<rrrmt*,HDereferenceLess> &rset, int mode ) {
	for( vector<rrrmt*>::iterator it=sce.begin(); it!=sce.end(); it++ ){
		(*it)->getrrrmt(rset, mode);
	}
}
void ordered::getrrrmt( vector<rrrmt*> &rset, int mode ) {
	for( vector<rrrmt*>::iterator it=sce.begin(); it!=sce.end(); it++ ){
		(*it)->getrrrmt(rset, mode);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ordered::getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const{
	for(vector<rrrmt*>::const_iterator it = sce.begin(); it!=sce.end(); it++){
		if((*it)->type()==1 || (*it)->type()==2 || (*it)->type()==3 || (*it)->type()==7 || (*it)->type()==8){
			scen.push_back((*it)->clone());
		}else{
			(*it)->getminrrrmt(scen, cst);
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ordered::hascomposite() const{
	for( vector<rrrmt*>::const_iterator it=sce.begin(); it!=sce.end(); it++ ){
		if( (*it)->type() >= UORD ){
			return true;
		}
	}
	return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for ordered scenarios
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *ordered::intersect(const rrrmt *rrrmtp, bool prefix ) const{
	return rrrmtp->intersect( this, prefix);
}
rrrmt *ordered::intersect( const indel *rr, bool prefix ) const{
//	#ifdef DEBUG_INTERSECTION
//	cout << "ordered::intersect "<<*this<< " <-> "<<rr<<endl;
//	#endif//DEBUG_INTERSECTION
	if( sce.size() > 0 ){
		if(prefix){
			return sce.front()->intersect(rr, prefix);
		}else{
			return sce.back()->intersect( rr, prefix ) ;
		}
	}else{
		return new emptyscen();
	}
}
rrrmt *ordered::intersect( const rev *rr, bool prefix ) const{
//	#ifdef DEBUG_INTERSECTION
//	cout << "ordered::intersect "<<*this<< " <-> "<<rr<<endl;
//	#endif//DEBUG_INTERSECTION
	if( sce.size() > 0 ){
		if(prefix){
			return sce.front()->intersect( rr, prefix ) ;
		}else{
			return sce.back()->intersect( rr, prefix ) ;
		}
	}else{
		return new emptyscen();
	}
}
rrrmt *ordered::intersect( const transp *rr, bool prefix ) const{
//	#ifdef DEBUG_INTERSECTION
//	cout << "ordered::intersect "<<*this<< " <-> "<<rr<<endl;
//	#endif//DEBUG_INTERSECTION
	if( sce.size() > 0 ){
		if(prefix){
			return sce.front()->intersect( rr, prefix ) ;
		}else{
			return sce.back()->intersect( rr, prefix ) ;
		}
	}else{
		return new emptyscen();
	}
}
rrrmt *ordered::intersect( const revtransp *rr, bool prefix ) const{
//	#ifdef DEBUG_INTERSECTION
//	cout << "ordered::intersect "<<*this<< " <-> "<<rr<<endl;
//	#endif//DEBUG_INTERSECTION
	if( sce.size() > 0 ){
		if(prefix){
			return sce.front()->intersect( rr, prefix ) ;
		}else{
			return sce.back()->intersect( rr, prefix ) ;
		}
	}else{
		return new emptyscen();
	}
}
rrrmt *ordered::intersect( const tdrl *rr, bool prefix ) const{
//	#ifdef DEBUG_INTERSECTION
//	cout << "ordered::intersect "<<*this<< " <-> "<<rr<<endl;
//	#endif//DEBUG_INTERSECTION
	if( sce.size() > 0 ){
		if(prefix){
			return sce.front()->intersect( rr, prefix ) ;
		}else{
			return sce.back()->intersect( rr, prefix ) ;
		}
	}else{
		return new emptyscen();
	}
}
rrrmt *ordered::intersect(const ordered *rr, bool prefix) const{

	ordered *tempo = create();	// intersection of scenarios, i.e. biggest suffix
	rrrmt *tmp = NULL, *ret;
	set<rrrmt*,HDereferenceLess> its;
	vector<rrrmt*>::iterator it;
	unsigned max = 0,
			i, j;

	// compute the intersection of the last (resp. first) element of one
	// of the scenarios with the other scenario
	if(prefix && sce.size()>0){
		tmp = sce.front()->intersect(rr, prefix);
	}else if(!prefix && sce.size()>0){
		tmp = sce.back()->intersect(rr, prefix);
	}
	if( tmp !=NULL && ! its.insert(tmp).second )
		delete tmp;

	if(prefix && sce.size()>0){
		tmp = rr->sce.front()->intersect( this, prefix );
	}else if(!prefix && sce.size()>0){
		tmp = rr->sce.back()->intersect( this, prefix );
	}
	if( tmp != NULL && !its.insert(tmp).second )
		delete tmp;

	// intersection of scenarios, i.e. biggest suffix of complete intersections
	// i.e. sce[i] = rr.sce[i] plus one nonempty incomplete
	if( prefix ){
		i = 0;
		j = 0;
	}else{
		i = sce.size()-1;
		j = rr->sce.size()-1;
	}
	do{
		tmp = sce[ i ]->intersect( rr->sce[ j ], prefix );
		if( !tmp->isempty() ){
			if(prefix){
				tempo->sce.push_back( tmp );
			}else{
				tempo->sce.insert( tempo->sce.begin(), tmp );
			}
			if( tmp->length( ALTMAX ) < sce[ i ]->length(ALTMAX) || tmp->length( ALTMAX ) < rr->sce[ j ]->length(ALTMAX) ){
				break;
			}
		}else{
			delete tmp;
			break;
		}

		if(prefix){
			i++;
			j++;
		}else{
			i--;
			j--;
		}
	}while(i >= 0 && j >= 0 && i < sce.size() && j < rr->sce.size());
	tmp = tempo->simplify();
	delete tempo;
	if( !its.insert( tmp ).second)
		delete tmp;

	max = getmax(its);
#ifdef DEBUG_INTERSECTION
	cerr << "ordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	cerr << "   max "<<max<<" " << endl;
#endif//DEBUG_INTERSECTION
	if( max == 0 ){
		tmp = new emptyscen;
	}else{
		if(its.size() != 1){
			cerr << "warning: found "<<its.size()<<"!= 1 maximal intersections in ordered::intersect "<<*this<< " vs. "<<*rr<<endl;
			for( set<rrrmt*, DereferenceLess>::const_iterator it=its.begin(); it!=its.end(); it++){
				(*it)->output(cerr, 0, 1); cerr<<endl;
			}
		}


		tmp = (*its.begin())->clone();
	}
	#ifdef DEBUG_INTERSECTION
	cerr << " -> "; tmp->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION

	for(set<rrrmt*,HDereferenceLess>::iterator it = its.begin(); it!=its.end(); it++){
		delete *it;
	}
	its.clear();

	// simplify result if possible
	ret = tmp->simplify();
	delete tmp;

	if(ret->isempty()){
		delete ret;
		ret = new emptyscen();
	}

	return ret;
}
rrrmt *ordered::intersect(const unordered *rr, bool prefix) const{
	ordered *tempo;			// temporary ordered rearrangement
	rrrmt *tmp = NULL, *ret;				// the rearrangment to return
	set<rrrmt*,HDereferenceLess> its,
		itstmp;
	vector<rrrmt*> rrsce;				// temporary copy of rr.sce (for using it as a vector)
	unsigned max = 0;

	rrsce.insert(rrsce.begin(), rr->_sce.begin(), rr->_sce.end());

	// intersection of last (resp. first) element of ordered with the complete unordered
	if(prefix && sce.size()>0 ){
		tmp = sce.front()->intersect( rr, prefix );
	}else if( !prefix && sce.size()>0){
		tmp = sce.back()->intersect( rr, prefix );
	}
	if( tmp != NULL && ! its.insert(tmp).second )
		delete tmp;

	// intersection of the complete ordered with every element of the unordered
	for(unsigned i=0; i<rrsce.size(); i++){
		tmp = rrsce[i]->intersect( this, prefix );
		if( ! its.insert(tmp).second )
			delete tmp;
	}

	// intersection of every vs every
	tempo = create();
	for(unsigned i=0; i<sce.size(); i++){
		for( unsigned j=0; j<rrsce.size(); j++ ){
			if(prefix){
				tmp = sce[i]->intersect(rrsce[j], prefix);
			}else{
				tmp = sce[sce.size()-1-i]->intersect(rrsce[j], prefix);
			}
			if( ! itstmp.insert(tmp).second )
				delete tmp;
		}
		max = getmax( itstmp );
		if( max > 0 ){
			tempo->sce.insert( tempo->sce.begin(), (*itstmp.begin())->clone() );
			if(itstmp.size() != 1){
				cerr << "warning: found != 1 maximal intersections in ordered::intersect "<<*this<< " vs. "<<rr<<endl;
				cerr << "         during intersection of every vs every"<<endl;
			}
		}
		for( set<rrrmt*, HDereferenceLess>::iterator it=itstmp.begin(); it!=itstmp.end(); it++)
			delete *it;
		itstmp.clear();
			// check for incomplete intersection
		if( max < sce[ sce.size()-1-i]->length(ALTMAX) ){
			break;
		}
		max = 0;
	}
	tmp = tempo->simplify();
	delete tempo;
	if( ! its.insert( tmp ).second )
		delete tmp;
	tmp = NULL;

	// determine the max of its[0]...its[2]
	max = getmax(its);
	#ifdef DEBUG_INTERSECTION
	cerr << "ordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	cerr << "   max "<<max<< endl;
	#endif//DEBUG_INTERSECTION
	if( max == 0 ){
		tmp = new emptyscen;
	}else{
		if(its.size() != 1){
			cerr << "warning: found != 1 maximal intersections in ordered::intersect "<<*this<< " vs. "<<rr<<endl;
			for( set<rrrmt*, DereferenceLess>::const_iterator it=its.begin(); it!=its.end(); it++){
				(*it)->output(cerr, 0, 1); cerr<<endl;
			}
		}
		tmp = (*its.begin())->clone();
	}
	#ifdef DEBUG_INTERSECTION
	cerr << " -> "; tmp->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION

	for( set<rrrmt*,HDereferenceLess>::iterator it=its.begin(); it!=its.end(); it++){
		delete *it;
	}
	its.clear();
	rrsce.clear();

	// simplify result if possible
	ret = tmp->simplify();
	delete tmp;

	if(ret->isempty()){
		delete ret;
		ret = new emptyscen();
	}

	return ret;
}

rrrmt *ordered::intersect( const alternative *rr, bool prefix ) const{
	return rr->intersect( this, prefix);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ordered::iscomplete() const{
	if(!_complete){
		return false;
	}
	for(vector<rrrmt*>::const_iterator it=sce.begin(); it!=sce.end(); it++){
		if( !(*it)->iscomplete() ){
			return false;
		}
	}
	return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ordered::isempty() const{
	return (sce.size() == 0);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned ordered::length( int mode ) const{
	unsigned s = 0;
	for(unsigned i=0; i<sce.size(); i++){
		s += sce[i]->length(mode);
	}
	return s;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *ordered::mkealt(){
	rrrmt *tmp;
	for(unsigned i=0; i<sce.size(); i++){
		if( !sce[i]->iscomplete() )
			continue;
		tmp = sce[i]->mkealt();
		if( !tmp->isempty() ){
			delete sce[i];
			sce[i] = tmp;
		}
		else{
			delete tmp;
		}
	}
	return new emptyscen;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//rrrmt *ordered::mkepar() const{
//	ordered *ret = new ordered();
//	rrrmt *t;
//
//	for( vector<rrrmt*>::const_iterator it=sce.begin(); it!=sce.end(); it++ ){
//		t = (*it)->mkepar();
//		if( !t->isempty() ){
//			ret->push_back( t );
//		}else{
//			delete t;
//		}
//	}
//	return ret;
//}
unsigned ordered::nralt( ) const{
	unsigned a = 1;
	for(vector<rrrmt*>::const_iterator it = sce.begin(); it!=sce.end(); it++){
		a *= (*it)->nralt();
	}
	return a;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

float ordered::getmincost(const costfoo *cst) const{
	float sum=0.0;
	for(vector<rrrmt*>::const_iterator it = sce.begin(); it!=sce.end(); it++){
		sum += (*it)->getmincost(cst);
	}
	return sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool ordered::operator<(const rrrmt *rrrmtp) const{
	if(type() != rrrmtp->type()){
		return type() < rrrmtp->type();
	}
	const ordered* ordp = dynamic_cast<const ordered*>(rrrmtp);
	for(unsigned i=0; i<sce.size() && i<ordp->sce.size(); i++){
		if( *(sce[i]) < ordp->sce[i] ){
			return true;
		}else if( *(ordp->sce[i]) < sce[i] ){
			return false;
		}
	}
	return ( sce.size() < ordp->sce.size() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &ordered::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	string indent = "";
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) indent+="\t";
	}

	out << indent;
	out << typestrg(d);
	out <<pquot<<"[";
	if (d==0){
		out<<endl;
	}
	for(unsigned i=0; i<sce.size(); i++){
		sce[i]->output(out, l+1, d, pquot);
		if( i!=sce.size()-1 )
			out << ",";
	}
	out<<indent;
	out<<pquot<<"]";

	if( d==0 ){
		out<<" complete="<<_complete<<endl;
	}
	return out;
}

ostream &ordered::output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx ) {
	out<<"[";
	for(vector<rrrmt*>::iterator it=sce.begin(); it != sce.end(); it++){
		if(it!=sce.begin())
			out << ",";
		(*it)->output(out, ridx);

	}
	out<<"]";
	return out;
}

//void ordered::shortdesc(string &d) const{
//	d+="\\[";
//	for(vector<rrrmt*>::const_iterator it=sce.begin(); it != sce.end(); it++){
//		(*it)->shortdesc(d);
//		d+=",";
//	}
//	d+="\\]";
//}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ordered::push_back(rrrmt *r){
	sce.push_back(r);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ordered::rename(const vector<vector<int> > &mapping, vector<string> *nmap, unsigned mx){
	for(unsigned i=0; i<sce.size(); i++){
		sce[i]->rename(mapping, nmap, mx);
	}
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// replace all elements of the ordered scenario by a simplified version (if it exist)
// if the ordered scenario has length one return a copy of this one rearrangement
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rrrmt *ordered::simplify(){
	ordered *ret = create();
	rrrmt *tmp;

	// - get a simplified version of the child nodes
	// - lift all children of ordered children one level up
	for(unsigned i=0; i<sce.size(); i++){
		tmp = sce[i]->simplify();

		ordered *tmpo = dynamic_cast<ordered*>(tmp);
		if( tmpo && ret->type() == tmpo->type() ){
			for(vector<rrrmt*>::iterator jt=tmpo->sce.begin(); jt != tmpo->sce.end(); jt++){
				ret->push_back( *jt );
			}
			tmpo->sce.clear();
			delete tmpo;
		}else{
			ret->push_back( tmp );
		}

	}

	if( sce.size() == 0 ){
		delete ret;
		tmp = new emptyscen();
		return tmp;
	}else if( sce.size() == 1 ){
		tmp = sce[0]->clone();
		delete ret;
		return tmp;
	}else{
		return ret;
	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int ordered::type() const{
	return ORD;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string ordered::typestrg(unsigned d) const{
	if(d == 0){
		return "ordered";
	}else{
		return "";
	}
}

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void alternative::append( vector<vector<rrrmt*> > &v) const{
	unsigned min = length(ALTMIN);

	for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
		if( (*it)->iscomplete() && (*it)->length(ALTMIN) == min ){
			(*it)->append( v );
//			(v.back()).push_back( (*it)->clone() );
			break;
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply an alternative scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void alternative::apply(genom &g) const{
	if( _sce.size() == 0){
		cerr << "error: apply called for zero size alternative"<<endl;
		exit(0);
	}else{
		unsigned min = length( ALTMIN );
		for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
			if( (*it)->length(ALTMIN) == min ){
				(*it)->apply(g);
				break;
			}
		}
	}
}



void alternative::apply( genom &g, set<genom> &gs )const{

		for( set<rrrmt*,HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++){
			genom tmp = g;
			(*it)->apply( tmp, gs );
		}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

alternative* alternative::clone() const{
	return new alternative(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

alternative* alternative::create() const{
//	cout << "alternative::create()"<<endl;
	return new alternative();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

set<rrrmt*, HDereferenceLess> alternative::get_alternatives() const{
	set<rrrmt*, HDereferenceLess> s;

	for(set<rrrmt*, HDereferenceLess>::const_iterator it=_sce.begin(); it!=_sce.end(); it++){
		s.insert( (*it)->clone() );
	}
	return s;
}

rrrmt * alternative::get_random_alternative() const{
	if( _sce.size() == 0 ){
		return new emptyscen();
	}else{
		set<rrrmt*, HDereferenceLess>::const_iterator it(_sce.begin());
		advance(it,ask_rng( 0, _sce.size()-1 ));
		return (*it)->clone();
	}
}


void alternative::getrrrmt( set<rrrmt*,HDereferenceLess> &rset, int mode ) {
	int m;
	set<rrrmt*,HDereferenceLess> mset;

	switch(mode){
	case ALTALL:
		break;
	case ALTMIN:
		m = std::numeric_limits< int >::max();
		break;
	case ALTMAX:
		m = -1;
		break;
	default:
		cerr << "invalid mode "<< mode <<" given to alternative::getrrrmt()"<<endl;
		exit(1);
	}

	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		switch(mode){
		case ALTMIN:
			if((int)(*it)->length( mode ) < m ){
				m = (*it)->length( mode );
				for(set<rrrmt*,HDereferenceLess>::iterator jt=mset.begin(); jt!=mset.end(); jt++){delete *jt;} mset.clear();
				(*it)->getrrrmt(mset, mode);
			}
			break;
		case ALTMAX:
			if((int)(*it)->length( mode ) > m ){
				m = (*it)->length( mode );
				for(set<rrrmt*,HDereferenceLess>::iterator jt=mset.begin(); jt!=mset.end(); jt++){delete *jt;} mset.clear();
				(*it)->getrrrmt(mset, mode);
			}
			break;
		case ALTALL:
			(*it)->getrrrmt(mset, mode);
			break;
		}
	}

	rset.insert( mset.begin(), mset.end() );
	mset.clear();
}
void alternative::getrrrmt( vector<rrrmt*> &rset, int mode ) {
	int m;
	set<rrrmt*,HDereferenceLess> mset;

	switch(mode){
	case ALTMIN:
	case ALTALL:
		m = std::numeric_limits< unsigned >::max();
		break;
	case ALTMAX:
		m = -1;
		break;
	default:
		cerr << "invalid mode "<< mode <<" given to alternative::getrrrmt()"<<endl;
		exit(1);
	}

	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		switch(mode){
		case ALTMIN:
			if((int)(*it)->length( mode ) < m ){
				m = (*it)->length( mode );
				for(set<rrrmt*,HDereferenceLess>::iterator jt=mset.begin(); jt!=mset.end(); jt++){delete *jt;} mset.clear();
				(*it)->getrrrmt(mset, mode);
			}
			break;
		case ALTMAX:
			if((int)(*it)->length( mode ) > m ){
				m = (*it)->length( mode );
				for(set<rrrmt*,HDereferenceLess>::iterator jt=mset.begin(); jt!=mset.end(); jt++){delete *jt;} mset.clear();
				(*it)->getrrrmt(mset, mode);
			}
			break;

		case ALTALL:
			(*it)->getrrrmt(rset, mode);
			break;
		}
	}

	//rset.insert(rset.end(), mset.begin(), mset.end() );
	mset.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void alternative::getminrrrmt(vector<rrrmt*> &scen, const costfoo *cst) const{
	float tempcost=999999.9;
	rrrmt* tempscen;
	for(set<rrrmt*,DereferenceLess>::const_iterator it = _sce.begin(); it!=_sce.end(); it++){
		if(tempcost <= (*it)->getmincost(cst)){
			//do nothing
		}else{
			tempscen=(*it)->clone();
			tempcost=(*it)->getmincost(cst);
		}
	}
	if(tempscen->type()==1 || tempscen->type()==2 || tempscen->type()==3 || tempscen->type()==7 || tempscen->type()==8){
		// elementar rrrmt
		scen.push_back(tempscen->clone());
	}else{
		//ordered, unordered or alternative
		tempscen->getminrrrmt(scen, cst);
	}
	delete tempscen;
}


bool alternative::equal( const rrrmt *r ) const{
	return r->equal(this);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// double dispatch for ordered scenarios
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



rrrmt *alternative::intersect(const rrrmt *rrrmtp, bool prefix ) const{

#ifdef DEBUG_INTERSECTION
cerr << "start alternative::intersect (***) "; output(cerr, 0, 1); cerr<< " vs. "; rrrmtp->output(cerr, 0, 1); cerr<<endl;
#endif//DEBUG_INTERSECTION
	return rrrmtp->intersect( this, prefix);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


rrrmt *alternative::intersect( const unordered *rr, bool prefix ) const{

	alternative *atmp;
	rrrmt *tmp, *ret;						// temp variable for rearrangements
	set<rrrmt*, HDereferenceLess> its;				// candidate intersections
	vector<vector<rrrmt*> > incpy; 	// temporary copy of the input scenarios

	#ifdef DEBUG_INTERSECTION
	cerr << "start alternative::intersect (uno) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION


//	// make a copy of the input, for using it as a vector
//	incpy = vector<vector<rrrmt*> >( 2, vector<rrrmt*>() );
//	incpy[0].insert(incpy[0].begin(), _sce.begin(), _sce.end());
//	incpy[1].insert(incpy[1].begin(), rr->_sce.begin(), rr->_sce.end());

	// intersection of the complete unordered scenario with every element of the other
	for( set<rrrmt*, DereferenceLess>::const_iterator it=this->_sce.begin(); it != this->_sce.end(); it++ ){
		tmp = (*it)->intersect( rr, prefix );
		if( ! its.insert( tmp ).second )
			delete tmp;
	}
	for( set<rrrmt*, DereferenceLess>::const_iterator it=rr->_sce.begin(); it != rr->_sce.end(); it++ ){
		tmp = (*it)->intersect( this, prefix );
		if( ! its.insert( tmp ).second )
			delete tmp;
	}

	// determine the maxima of its[0],its[1]
	getmax(its);
	atmp = new alternative();
	for( set<rrrmt*,DereferenceLess>::const_iterator it=its.begin(); it!=its.end(); it++ ){
		atmp->insert( (*it)->simplify() );
	}
	#ifdef DEBUG_INTERSECTION
	cerr << "alternative::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	cerr << " -> "; atmp->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION

	for( set<rrrmt*,HDereferenceLess>::iterator it=its.begin(); it!=its.end(); it++){
		delete *it;
	}
	its.clear();
	incpy.clear();

	// simplify result if possible
	ret = atmp->simplify();
	delete atmp;

	if(ret->isempty()){
		delete ret;
		ret = new emptyscen();
	}

	return ret;
}

rrrmt *alternative::intersect( const ordered *rr, bool prefix ) const{
	alternative *atmp;
	rrrmt *tmp=NULL, *ret;
	set<rrrmt*, HDereferenceLess> its;
	vector<rrrmt*> itstmp,
		altsce;				// temporary copy of the alternative (for using it as a vector)

	#ifdef DEBUG_INTERSECTION
	cerr << "start alternative::intersect (ord) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION


	altsce.insert(altsce.begin(), _sce.begin(), _sce.end());

	// intersection of last element of ordered with the complete alternative
	if(prefix && rr->sce.size()>0){
		tmp = rr->sce.front()->intersect( this, prefix );
	}else if(!prefix && rr->sce.size()>0){
		tmp = rr->sce.back()->intersect( this, prefix );
	}
	if( tmp !=NULL && ! its.insert( tmp ).second)
		delete tmp;

	// intersection of the complete ordered with every element of the alternative
	for(unsigned i=0; i<altsce.size(); i++){
		tmp = altsce[i]->intersect( rr, prefix );
		if(! its.insert( tmp ).second )
			delete tmp;
	}

		// determine the max of its[0], its[1]
	getmax(its);
	#ifdef DEBUG_INTERSECTION
	cerr << "ordered::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION
	atmp = new alternative();
	for( set<rrrmt*, DereferenceLess>::const_iterator it = its.begin(); it!=its.end(); it++ ){
		atmp->insert( (*it)->simplify() );
	}

	#ifdef DEBUG_INTERSECTION
	cerr << "-> "; tmp->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION

	for(set<rrrmt*,HDereferenceLess>::iterator it=its.begin(); it != its.end(); it++){
		delete *it;
	}
	its.clear();
	altsce.clear();

	// simplify result if possible
	ret = atmp->simplify();
	delete atmp;
	if(ret->isempty()){
		delete ret;
		ret = new emptyscen();
	}

	return ret;
}

rrrmt *alternative::intersect( const alternative *rr, bool prefix ) const{
	alternative *atmp;
	rrrmt *tmp, *ret;
	set<rrrmt*,HDereferenceLess> its;
	vector<rrrmt*>::iterator it;
	vector<vector<rrrmt*> > itstmp,
		incpy;

	#ifdef DEBUG_INTERSECTION
	cerr << "start alternative::intersect (alt) "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION


	incpy = vector<vector<rrrmt*> >( 2, vector<rrrmt*>() );
	incpy[0].insert(incpy[0].begin(), this->_sce.begin(), this->_sce.end());
	incpy[1].insert(incpy[1].begin(), rr->_sce.begin(), rr->_sce.end());

	// intersection of the complete scenario with every element of the other
	for( unsigned i=0; i<incpy[0].size(); i++ ){
		tmp = incpy[0][i]->intersect( rr, prefix );
		if( ! its.insert( tmp ).second )
			delete tmp;
	}

	for( unsigned i=0; i<incpy[1].size(); i++ ){
		tmp = incpy[1][i]->intersect( this, prefix );
		if( ! its.insert( tmp ).second )
			delete tmp;
	}

//	// intersection of every vs every
////	itstmp = vector<vector<rrrmt*> >( incpy[0].size(), vector<rrrmt*>(incpy[1].size()) );
//	for( unsigned i=0; i<incpy[0].size(); i++ ){
//		for( unsigned j=0; j<incpy[1].size(); j++ ){
//			tmp = incpy[0][i]->intersect( incpy[1][j], prefix );
//			if( ! its.insert( tmp ).second )
//				delete tmp;
//		}
//	}


//
//	maxrc = vector<vector<unsigned> >(2 );
//	maxidxrc = vector<vector<vector<unsigned> > >(2);
//	// get maxima of every row
//	maxrc[0] = vector<unsigned>( incpy[0].size() );
//	maxidxrc[0] = vector<vector<unsigned> >( incpy[0].size() );
//	for( unsigned i=0; i<itstmp.size(); i++ ){
//		getmax( itstmp[i], maxrc[0][i], maxidxrc[0][i] );
//	}
//
//	// get maxima of every column
//	maxrc[1] = vector<unsigned>( incpy[1].size() );
//	maxidxrc[1] = vector<vector<unsigned> >( incpy[1].size() );
//	for( unsigned i=0; i<itstmp[0].size(); i++ ){
//		transpitstmp = vector<rrrmt*>( itstmp.size());
//		for( unsigned j=0; j<itstmp.size(); j++ ){
//			transpitstmp[j] = itstmp[j][i];
//		}
//		getmax( transpitstmp, maxrc[1][i], maxidxrc[1][i] );
//		transpitstmp.clear();
//	}
//
//	tempa = create();
//	for( unsigned i=0; i<itstmp.size(); i++ ){
//		if( maxrc[0][i] > 0 ){
//			if(maxidxrc[0][i].size() != 1){
//				cerr << "warning: found != 1 maximal intersections in alternative::intersect "; this->output(cerr, 0, 0); cerr<< " vs. "; rr->output(cerr, 0, 0); cerr<<endl;
//				cerr << "         during intersection of the complete scenario with every element of the other"<<endl;
//				exit(EXIT_FAILURE);
//			}
//
//			if( maxidxrc[1][ maxidxrc[0][i][0] ].size() != 1 ){
//				cerr << "warning: found != 1 maximal intersections in alternative::intersect "; this->output(cerr, 0, 0); cerr<< " vs. "; rr->output(cerr, 0, 0); cerr<<endl;
//				cerr << "         During intersection of the complete scenario with every element of the other"<<endl;
//				cerr << itstmp.size()<<"x"<<itstmp[0].size()<<endl;
//				for( unsigned j=0; j<itstmp.size(); j++ ){
//					for( unsigned k=0; k<itstmp[j].size(); k++){
//						cerr << itstmp[j][k]->length(ALTMAX);
//					}cerr << endl;
//				}
//				exit(EXIT_FAILURE);
//			}
//
//			if( maxidxrc[1][ maxidxrc[0][i][0] ][0] != i ){
//				cerr << "warning: max intersection mismatch in alternative::intersect "; this->output(cerr, 0, 0); cerr<< " vs. "; rr->output(cerr, 0, 0); cerr<<endl;
//				cerr << "         do not match "<<endl;
//				exit(EXIT_FAILURE);
//			}
//			tempa->insert( itstmp[i][maxidxrc[0][i][0]]->clone() );
//		}
//	}
//	tmp = tempa->simplify();
//	delete tempa;
//	ir = its.insert( tmp );
//	if( ! ir.second )
//		delete tmp;

	// determine the max of its[0]...its[2]
	getmax(its);
	#ifdef DEBUG_INTERSECTION
	cerr << "alternative::intersect "; output(cerr, 0, 1); cerr<< " vs. "; rr->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION
	atmp = new alternative();
	for( set<rrrmt*, DereferenceLess>::const_iterator it=_sce.begin(); it!=_sce.end(); it++){
		atmp->insert( (*it)->simplify() );
	}
	#ifdef DEBUG_INTERSECTION
	cerr << " -> "; atmp->output(cerr, 0, 1); cerr<<endl;
	#endif//DEBUG_INTERSECTION
//	for( unsigned i=0; i<itstmp.size(); i++ ){
//		for( unsigned j=0; j<itstmp[i].size(); j++ ){
//			delete itstmp[i][j];
//		}
//	}
//	itstmp.clear();
	for(set<rrrmt*,HDereferenceLess>::iterator it = its.begin(); it!=its.end(); it++){
		delete *it;
	}
	its.clear();
	incpy.clear();

	// simplify result if possible
	ret = atmp->simplify();
	delete atmp;
	if(ret->isempty()){
		delete ret;
		ret = new emptyscen();
	}

	return ret;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// de-apply an alternative scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void alternative::deapply(genom &g) const{
	if( _sce.size() == 0){
		cerr << "error: apply called for zero size alternative"<<endl;
		exit(0);
	}else{
		unsigned min = length( ALTMIN );

		for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
			if( (*it)->iscomplete() && (*it)->length(ALTMIN) == min ){
				(*it)->deapply(g);
				break;
			}
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool alternative::iscomplete() const{
	if(!_complete){
		return false;
	}
	for(set<rrrmt*,HDereferenceLess>::const_iterator it=_sce.begin(); it!=_sce.end(); it++){
		if( (*it)->iscomplete() ){
			return true;
		}
	}
	return false;
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool alternative::isempty() const{
	if (_sce.size() == 0){
		return true;
	}else{
		for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++){
			if( !(*it)->isempty() ){
				return false;
			}
		}
		return true;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned alternative::length( int mode ) const{
//	bool c = false;

	unsigned s;

	switch(mode){
	case ALTMIN:
		s = std::numeric_limits< unsigned >::max();
		break;
	case ALTMAX:
		s = 0;
		break;
	default:
		cerr << "invalid mode "<< mode <<" given to alternative::length()"<<endl;
		exit(1);
	}



	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
//		if( (*it)->iscomplete() ){
			switch(mode){
			case ALTMIN:
				s = min( s, (*it)->length( mode ) );
				break;
			case ALTMAX:
				s = max( s, (*it)->length( mode ) );
				break;
			}
//		}
	}
	return s;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//rrrmt *alternative::mkepar() const{
//	rrrmt *ret = new emptyscen(),
//		*t = NULL;
//
//	unsigned min = std::numeric_limits< unsigned >::max();
//	for( set<rrrmt*, DereferenceLess>::iterator it=sce.begin(); it!=sce.end(); it++ ){
//		if( (*it)->size() < min ){
//			delete ret;
//			t = (*it)->mkepar();
//			ret = t;
//			min = (*it)->size();
//		}
//	}
//
//	return ret;
//}

unsigned alternative::nralt( ) const{
	unsigned a = 0;
	for(set<rrrmt*,DereferenceLess>::const_iterator it = _sce.begin(); it!=_sce.end(); it++){
		a += (*it)->nralt();
	}
	return a;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


float alternative::getmincost(const costfoo *cst) const{
	float sum=999999.9;
	for(set<rrrmt*,DereferenceLess>::const_iterator it = _sce.begin(); it!=_sce.end(); it++){
		if(sum <= (*it)->getmincost(cst)){
			//do nothing
		}else{
			sum=(*it)->getmincost(cst);
		}
	}
	return sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ostream &alternative::output(ostream &out, unsigned l, unsigned d, string pquot  ) const{
	string indent = "";
//	unsigned min = length( ALTMIN );
	if( d == 0 ){
		for(unsigned i=0; i<l; i++) indent+="\t";
	}
	out << indent;
	out << typestrg(d);
	out <<pquot<<"(";
	if (d==0){
		out<<endl;
	}
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		(*it)->output(out, l+1, d);
		out << ",";
	}
	out<<indent;
	out<<pquot<<")";

	if( d==0 ){
		out<<" complete="<<_complete<<endl;
	}
	return out;
}
ostream &alternative::output( ostream &out, map<rrrmt*, unsigned, HDereferenceLess> &ridx ) {
	out<<"(";
	for(set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it != _sce.end(); it++){
		if(it!=_sce.begin())
			out << ",";
		(*it)->output(out, ridx);

	}
	out<<")";
	return out;
}

//void alternative::shortdesc(string &d) const{
//	unsigned min = size();
//	d+="\\{";
//	for( set<rrrmt*, DereferenceLess>::iterator it=sce.begin(); it!=sce.end(); it++ ){
//		if( (*it)->size() == min ){
//			(*it)->shortdesc(d);
//			if( sce.size() > 1 ){
//				d+="..";
//			}
//			break;
//		}
//	}
//	d+="\\}";
//}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int alternative::type() const{
	return ALT;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string alternative::typestrg(unsigned d) const{
	if(d == 0){
		return "alternative";
	}else{
		return "";
	}

}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned dynamic_ridx( rrrmt *r, map<rrrmt*, unsigned, HDereferenceLess> &ridx ){

	if( ridx.find( r ) ==
			ridx.end() ){
		unsigned m = 0;
		for(map<rrrmt*, unsigned, HDereferenceLess>::const_iterator it=ridx.begin();
				it!=ridx.end(); it++){

			if( r->type() != (it->first)->type() ){
				continue;
			}
			m = max( m, it->second );
		}

		ridx[r] = m+1;
	}

	return ridx[r];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

map<rrrmt*, unsigned, HDereferenceLess> getrrrmtidx( const set<rrrmt*, HDereferenceLess> &rset ){

	map<rrrmt*, unsigned, HDereferenceLess> rrrmtidx;
	map<string, int> typecnt;
	string type;

	for( set<rrrmt*, HDereferenceLess>::const_iterator it=rset.begin(); it!=rset.end(); it++ ){
		if( rrrmtidx.find( *it ) != rrrmtidx.end() ){
			continue;
		}
		type = (*it)->typestrg( 1 );
		if( typecnt.find(type) == typecnt.end() ){
			typecnt[ type ] = 0;
		}

		typecnt[ type ]++;

		rrrmtidx[ *it ] = typecnt[ type ];
	}
	return rrrmtidx;
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void getmax( const vector<rrrmt*> &its, unsigned &max, vector<unsigned> &maxidx){
//	unsigned s = 0;
//
//	max = 0;
//	maxidx.clear();
//	for(unsigned i=0; i<its.size(); i++){
//		s = its[i]->length( ALTMAX );
//		if( s > 0 && s >= max){
//			if( s > max ){
//				max = s;
//				maxidx.clear();
//			}
//			maxidx.push_back(i);
//		}
//	}
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned getmax( set<rrrmt*,HDereferenceLess> &its ){
	unsigned l = 0, a = 0,
		mxl = 0, mxa = 0;
	vector< set<rrrmt*,HDereferenceLess>::iterator > todel, tokeep;

	for(set<rrrmt*,HDereferenceLess>::iterator it=its.begin(); it!=its.end(); ){
		l = (*it)->length( ALTMAX );
		a = (*it)->nralt(  );

		if( l>0 && ( l > mxl || (l == mxl && a >= mxa ) ) ){
			if( l > mxl || (l == mxl && a > mxa ) ){
				mxl = l;
				mxa = a;
				its.erase( its.begin(), it ); // remove everything up to (not including it)
			}
			++it;
		}else{
			delete *it;
			its.erase(it++);
		}
	}
	for( vector< set<rrrmt*,HDereferenceLess>::iterator >::iterator it = todel.begin(); it != todel.end(); ){
//		delete **it;
		its.erase(*it);
	}
//	for(set<rrrmt*,HDereferenceLess>::iterator it = its.begin(); it != its.end(); ){
//		l = (*it)->length( ALTMAX );
//		if ( l == 0 || l < mxl){
//			delete *it;
//			its.erase(it++);
//		}else{
//			++it;
//		}
//	}
	return mxl;
}
