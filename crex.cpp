/**
 * author:  T. Hartmann and Matthias Bernt
 * (coauthors: D.Merkle, M. Middendorf, K. Ramsch)
 */
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <sstream>

#include "crex.hpp"
#include "sb4type.hpp"


#define MARK_RED "#ff6f00"
#define MARK_GREEN "#6fff00"
#define MARK_BLUE "006fff"

//#define DEBUG_BPSCENARIO
//#define DEBUG_CREX
//#define DEBUG_CREX2
//#define DEBUG_CREX2_PRIM
//#define DEBUG_CREX2_PRIM_APPROX
//#define DEBUG_CREXTRA
//#define DEBUG_CREXREV
//#define DEBUG_CREXRTRA
//#define DEBUG_CREXTDRL
//#define DEBUG_CREX_INDELS

using namespace std;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void bpscenario::bpscenario_add( genom g1, const genom &g2, vector<rrrmt*> &current,
		unsigned maxalt, unordered *check, bool complete ){

	bool success  = false,	// found some rearrangement(s)
		comp = false;		// final scenario is complete ?
	vector<int> pi,
		invpi;

	if( maxalt > 0 && _sce.size() >= maxalt )
		return;

		// apply the last operation on g1 (if there is one)
		// and get the identified and their inverse of the result
	if(current.size() > 0){
		(current.back())->apply(g1);
	}
	pi = g2.identify(g1);
	invpi = vector<int>(pi.size(), 0);
	for(unsigned i=0; i<pi.size(); i++){
		invpi[abs(pi[i])] = i;
	}
	pi.push_back( pi.size() );
	invpi.push_back( invpi.size() );

#ifdef DEBUG_BPSCENARIO
	cout << g1<<endl;
	cout << g2<<endl;
	copy(pi.begin(), pi.end(), ostream_iterator<int>(cout," ")); cout << endl;
	copy(invpi.begin(), invpi.end(), ostream_iterator<int>(cout," ")); cout << endl;
#endif//DEBUG_BPSCENARIO


	for(unsigned i=0; i<pi.size()-1; i++){
		#ifdef DEBUG_BPSCENARIO
		cout << "adjacency "<<pi[i]<<","<<pi[i+1]<<endl;
		#endif//DEBUG_BPSCENARIO
			// ignore conserved adjacencies
		if( pi[i+1] - pi[i] == 1 ){
			#ifdef DEBUG_BPSCENARIO
			cout << "   conserved "<<endl;
			#endif//DEBUG_BPSCENARIO
			continue;
		}

			//    a(g_i, g_{i+1})   a(g_j, g_{j+1})
			// -> b(g_i, -g_j)      b(-g_{i+1}, g_{j+1})

			// check for b(g_i, -g_j), i.e. start of a reversal
			// and   for b(-g_{j+1}, g_{i+1})), i.e. end of a reversed reversal
		if( (pi[i] >= 0 && pi[i+1] < 0) ){
			#ifdef DEBUG_BPSCENARIO
			cout << "   potential reversal start size"<<endl;
			cout << "   endpositions " << invpi[ pi[i]+1 ] <<","<<invpi[ abs(pi[i+1])+1 ]<<endl;
			#endif//DEBUG_BPSCENARIO

			// check if the potential reversal end elements are neighboured
			// and have the correct signs
			if( invpi[abs(pi[i+1])+1]-invpi[pi[i]+1]==1 &&
				pi[invpi[ pi[i]+1 ]] < 0 &&
				pi[invpi[ abs(pi[i+1])+1 ]] >= 0 ){

				int rs=i+1-1,
					re=invpi[pi[i]+1]-1;;
				if( rs > re ){
					swap(rs, re);
					rs++;
					re--;
				}
				#ifdef DEBUG_BPSCENARIO
				cout << "   reversal("<< rs <<","<< re <<")"<<endl;
				#endif//DEBUG_BPSCENARIO

				success = true;
				current.push_back( new rev( rs, re , g1 ) );
				bpscenario_add( g1, g2, current, maxalt, check, complete );
				delete current.back();
				current.pop_back();
			}
		}
		// (g_i,g_{i+1}) ... (g_j,g_{j+1}) ... (g_k, g_{k+1})
		// (g_i,g_{j+1}) ... (g_k,g_{i+1}) ... (g_j, g_{k+1})

		// check for a potential transposition start site
		// i.e. adjacency with same signs
		if( (pi[i] >= 0 && pi[i+1] >= 0) || (pi[i] < 0 && pi[i+1] < 0) ){
			#ifdef DEBUG_BPSCENARIO
			cout << "   potential transposition start site"<<endl;
			#endif//DEBUG_BPSCENARIO
			unsigned x = i, 	// current position
				nx;		// next position
			vector<int> tra;
			for(unsigned j=0; j<3; j++){
				nx = ( pi[x] < 0 ) ? invpi[(-1*pi[x])-1] : invpi[pi[x]+1];
				#ifdef DEBUG_BPSCENARIO
				cout << "   x "<<x<<" nx "<<nx<< " pi[x] "<<pi[x]<<" pi[nx] "<<pi[nx]<<endl;
				#endif//DEBUG_BPSCENARIO
				// check for:
				// - a sing change to the next adjacency
				// - that next pos. is bigger (or the cycle is completed)
				if( ((pi[x] < 0 && pi[nx] < 0) || (pi[x]>=0 && pi[nx]>=0)) &&
						(x < nx || j==2) ){
					x = nx-1;
					tra.push_back(x);
				}else{
					break;
				}
			}
			if ( tra.size()==3 && x == i ){

				#ifdef DEBUG_BPSCENARIO
				cout << "   transp: ";copy(tra.begin(), tra.end(), ostream_iterator<int>(cout," ")); cout << endl;
				#endif//DEBUG_BPSCENARIO

				success = true;
				current.push_back( new transp( tra[2], tra[0], tra[1], g1 ) );
				bpscenario_add( g1, g2, current, maxalt, check, complete );
				delete current.back();
				current.pop_back();
			}
			tra.clear();
		}
	}

	if( !success && current.size()>0 ){
		// insert current in the check structure .. if its not there
		// allready then insert it in the scenario
		if( check->insert( new unordered(current) ) ){
			comp = ( g1==g2 ) ? true : false;
			if( !complete || comp ){
				_sce.insert( new ordered(current, comp) );
			}
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// construct an emptyscen scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bpscenario::bpscenario() : alternative(){}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// construct an scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bpscenario::bpscenario(const genom &g1, const genom &g2, unsigned maxalt, bool complete): alternative(){
	vector<rrrmt*> current;
	unordered *check;
	check = new unordered();
	// start recursion

	bpscenario_add( g1, g2, current, maxalt, check, complete );
	delete check;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// copy construct an scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bpscenario::bpscenario(const bpscenario &c ) : alternative(c){}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// destruct
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bpscenario::~bpscenario(){}

string bpscenario::typestr() const{
	return "bpscenario";
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//bpscenario* bpscenario::clone() const{
//	return new bpscenario(*this);
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//bpscenario* bpscenario::create() const{
//	cout << "bpscenario::create()"<<endl;
//	return new bpscenario();
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//ostream &bpscenario::output(ostream &out, unsigned l) const{
//	out << typestrg()<<"[";
//	for(set<rrrmt*, DereferenceLess>::iterator it=sce.begin(); it != sce.end(); it++){
//		cout << *(*it) <<" ";
//	}
//	out<<"]";
//	return out;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// construct an emptyscen scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crex2::crex2() : alternative(){}

crex2::crex2( const genom &g1, const genom &g2, bool oriented, const costfoo *cst, bool &is_optimal,
		bool complete, bool all, const unsigned time_bound, const unsigned distance, const bool approx){

#ifdef DEBUG_CREX2
	cout << "CREX2 called with:" << endl;
	cout << "source (g1): ";
	for(int z=0; z<g1.size(); z++){
		cout << g1[z] << " ";
	}
	cout << endl;
	cout << "target (g2): ";
	for(int z=0; z<g2.size(); z++){
		cout << g2[z] << " ";
	}
	cout << endl;
#endif//DEBUG_CREX2
	int n;
	itnode *iroot;			// interval tree relative to g1
	list<itnode *> nodes; 	// nodes of the interval tree
	// DP matrix containing the scores and corresponding rearrangement scenarios
	vector<float> dps;
	float mn;
	vector<rrrmt *> dpr;
	n = g1.size();
	genom g; // identified genom
	// identify the permutations
	g = g2.identify_g(g1);
	g.set_nmap(g1.get_nmap());
	// get the SIT and its nodes
	interval_tree(g1, g2, n, &iroot);
	#ifdef DEBUG_CREX2
		cout << "CREX2 constructor: " << endl;
		cout << "resulting g1: ";
		for(int z=0; z<g1.size(); z++){
			cout << g1[z] << " ";
		}
		cout << endl;
		cout << "resulting g: ";
		for(int z=0; z<g.size(); z++){
			cout << g[z] << " ";
		}
		cout << endl;
		cout << "g mapped: " << g << endl;
		interval_tree_print(iroot, g, cout); cout << endl;
	#endif//DEBUG_CREX2

	dps = vector<float>(2, 0);
	dpr = vector<rrrmt*>(2, NULL);
	// note that g1 is the original genome and g is the identified genome such that g2=id
	_crex2_main(iroot, g1, g, dps, dpr, cst, complete, all, time_bound, is_optimal, distance, approx);

	if(oriented){
		if(distance!=0 && dpr[INC]->length(1)>distance){
			//do nothing
		}else{
			_sce.insert( dpr[ INC ]->simplify() );
		}
		delete dpr[ INC ];
		delete dpr[ DEC ];
	}else{
		// get the better of the two
		mn = min( dps[0], dps[1] );
		for( unsigned i=0; i<dps.size(); i++ ){
			if ( dps[i] == mn && ( all || _sce.size() == 0 ) ){
				if(distance!=0 && dpr[i]->length(1)>distance){
					//do nothing
				}else{
					_sce.insert( dpr[i]->simplify() );
				}
			}
			delete dpr[i];
		}
	}
	interval_tree_free(iroot);
}

crex2::crex2(const crex2 &c ) : alternative(c){}

crex2::~crex2(){}

void crex2::_crex2_leaf( itnode *n, const genom &g1, vector<float> &dps,
		vector<rrrmt*> &dpr, const costfoo *cst ){

	dps = vector<float>(2, 0);
	dpr = vector<rrrmt*>(2, NULL);

	// leaf has sign x
	// to get sign x -> score 0 & empty scenario
	dpr[n->sign] = new emptyscen(  );
	dps[n->sign] = 0;

	// leaf has sign x
	// to get sign -x -> score reversal & rev scenario
	dpr[1 - n->sign] = new rev( abs((n->i).first), abs((n->i).second), g1);
	dps[1 - n->sign] = (*cst)[ dpr[1 - n->sign] ];
}

void crex2::_crex2_lin_k0( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr){

	int sign = n->sign;		// target sign = sign of the node
	unordered *tmp;
	// k0: option 1 to keep the sign
	//     not doing anything at the node and make the childs get the same sign
	tmp = new unordered(  );

	ms = 0;
	for (unsigned i=0; i<cdps.size(); i++){
		ms += cdps[i][sign];
		tmp->insert( cdpr[i][sign]->clone() );
	}
	mr = tmp;
}

void crex2::_crex2_lin_k1( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst, bool all){

	int sign = 1-n->sign,		// target sign = - sign of the node
		nsign = n->sign;		// inverse target sign

	vector<float> s(4, 0);
	vector<unordered *> r(4, NULL);
	alternative *ar;
	rrrmt *tmp;

	// I
	r[0] = new unordered();
	tmp = new rev( abs((n->i).first), abs((n->i).second), g1 );
	s[0] = (*cst)[tmp];
	r[0]->insert( tmp );
	for( unsigned i=0; i<n->children.size(); i++ ){
		s[0] += cdps[i][ nsign ];
		r[0]->insert( cdpr[i][nsign]->clone() );
	}
//#ifdef DEBUG_CREX2
//	cout << "\tI "<<s[0]<<" "; r[0]->output(cout, 0, 1); cout <<endl;
//#endif//DEBUG_CREX2

	// T
	r[1] = new unordered();
	if( n->children.size() == 2 ){
		tmp = new transp(n->children[0]->i.first, n->children[0]->i.second+1, n->children[1]->i.second+1, g1);
		s[1] = (*cst)[tmp] + cdps[0][sign] + cdps[1][sign];
		r[1]->insert( tmp );
		r[1]->insert( cdpr[0][sign]->clone() );
		r[1]->insert( cdpr[1][sign]->clone() );
	}else{
		s[1] = std::numeric_limits< float >::max();
		// r[1] is left uninitialized since by using infty as weight it will not be chosen
	}
//#ifdef DEBUG_CREX2
//	cout << "\tT "<<s[1]<<" "; r[1]->output(cout, 0, 1); cout <<endl;
//#endif//DEBUG_CREX2
	// siT
	r[2] = new unordered();
	tmp = new revtransp(n->children[1]->i.first, n->children.back()->i.second,
				n->children[0]->i.first, n->children[0]->i.second, g1);
	r[2]->insert( tmp );
	r[2]->insert( cdpr[0][sign]->clone() );
	s[2] = (*cst)[tmp] + cdps[0][sign] ;
	for( unsigned i=1; i<n->children.size(); i+=1 ){
		s[2] += cdps[i][nsign];
		r[2]->insert( cdpr[i][nsign]->clone() );
	}
//#ifdef DEBUG_CREX2
//	cout << "\tsiT "<<s[2]<<" "; r[2]->output(cout, 0, 1); cout <<endl;
//#endif//DEBUG_CREX2
	// piT
	r[3] = new unordered();
	tmp = new revtransp(n->children[0]->i.first, n->children[n->children.size()-2]->i.second,
				n->children.back()->i.first, n->children.back()->i.second, g1);
	r[3]->insert( tmp );
	r[3]->insert( cdpr.back()[sign]->clone() );
	s[3] = (*cst)[tmp] + cdps.back()[sign] ;

	for( unsigned i=0; i<n->children.size()-1; i+=1 ){
		s[3] += cdps[i][nsign];
		r[3]->insert( cdpr[i][nsign]->clone() );
	}
//#ifdef DEBUG_CREX2
//	cout << "\tpiT "<<s[3]<<" "; r[3]->output(cout, 0, 1); cout <<endl;
//#endif//DEBUG_CREX2
	// get the best
	ar = new alternative();
	ms = *(min_element( s.begin(), s.end() ));


	//	random_shuffle( s.begin(), s.end() );
	for( unsigned i=0; i<s.size(); i++ ){
		if( s[i] == ms && ( all || ar->nralt() == 0 )){
			ar->insert( r[i] );
		}else{
			delete r[i];
		}
	}
	mr = ar->simplify();
	delete ar;
}

void crex2::_crex2_lin_k2( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst, bool all ){

	int sign = n->sign,			// target sign = sign of the node
			nsign = 1-sign;		// inverse target sign
	vector<float> s(4, std::numeric_limits< float >::max());
	vector<unordered *> r(4, NULL);
	alternative *ar;
	rrrmt *tmp;
	rrrmt *I, *T, *IT;
	I = new rev( 0, 1, g1 );
	T = new transp( 0, 1, 2, g1 );
	IT = new revtransp(0, 1, 2, g1);

	// I + T
	r[0] = new unordered();
	if( n->children.size() == 2 and (*cst)[T] <= (*cst)[I]  ){
		tmp = new rev(abs((n->i).first), abs((n->i).second), g1);
		s[0] = (*cst)[ tmp ];
		r[0]->insert( tmp );
		tmp = new transp(n->children[0]->i.first, n->children[0]->i.second+1, n->children[1]->i.second+1, g1);
		s[0] += (*cst)[ tmp ];
		r[0]->insert( tmp );
		s[0] += cdps[0][nsign] + cdps[1][nsign];
		r[0]->insert( cdpr[0][nsign]->clone() );
		r[0]->insert( cdpr[1][nsign]->clone() );
	}
//#ifdef DEBUG_CREX2
//	cout << "\tI+T "<<s[0]<<" ";
//	r[0]->output(cout, 0, 1);
//	cout <<endl;
//#endif//DEBUG_CREX2
	// T + piT
	// piT + T
	r[1] = new unordered();
	r[2] = new unordered();
	if( n->children.size() == 2  and (*cst)[IT] + (*cst)[T] <= (*cst)[I]  ){

		tmp = new transp(n->children[0]->i.first, n->children[0]->i.second+1, n->children[1]->i.second+1, g1);
		s[1] = (*cst)[ tmp ];
		r[1]->insert( tmp );
		tmp = new revtransp(n->children[1]->i.first, n->children[1]->i.second,
				n->children[0]->i.first, n->children[0]->i.second, g1);
		s[1] += (*cst)[tmp];
		r[1]->insert( tmp );
		s[1] += cdps[0][sign] + cdps[1][nsign];
		r[1]->insert(cdpr[0][sign]->clone());
		r[1]->insert(cdpr[1][nsign]->clone());

		tmp = new transp(n->children[0]->i.first, n->children[0]->i.second+1, n->children[1]->i.second+1, g1);
		s[2] = (*cst)[tmp];
		r[2]->insert( tmp );
		tmp = new revtransp(n->children[0]->i.first, n->children[0]->i.second,
				n->children[1]->i.first, n->children[1]->i.second, g1);
		s[2] += (*cst)[tmp];
		r[2]->insert( tmp );
		s[2] += cdps[0][nsign] + cdps[1][sign];
		r[2]->insert(cdpr[0][nsign]->clone());
		r[2]->insert(cdpr[1][sign]->clone());
	}
//#ifdef DEBUG_CREX2
//	cout << "\tT+piT "<<s[1]<<" ";
//	r[1]->output(cout, 0, 1);
//	cout <<endl;
//	cout << "\tpiT+T "<<s[2]<<" ";
//	r[2]->output(cout, 0, 1);
//	cout <<endl;
//#endif//DEBUG_CREX2
	// 2iT
	r[3] = new unordered();
	if( (*cst)[IT] <= (*cst)[I] ){ //TODO cst[RTRANM] <= cst[REVNM]

		tmp = new revtransp( n->children[0]->i.first, n->children[n->children.size()-2]->i.second,
				n->children.back()->i.first, n->children.back()->i.second, g1 );
		s[3] = (*cst)[tmp];
		r[3]->insert( tmp );

		tmp = new revtransp( n->children[1]->i.first, n->children.back()->i.second,
				n->children[0]->i.first, n->children[0]->i.second, g1 );
		s[3] += (*cst)[tmp];
		r[3]->insert( tmp );

		s[3] += cdps[0][nsign] + cdps.back()[nsign];
		r[3]->insert( cdpr[0][nsign]->clone() );
		r[3]->insert( cdpr.back()[nsign]->clone() );
		for( unsigned i=1; i<n->children.size()-1; i++ ){
			s[3] += cdps[i][sign];
			r[3]->insert( cdpr[i][sign]->clone() );
		}
	}
//#ifdef DEBUG_CREX2
//	cout << "\t2iT "<<s[3]<<" ";
//	r[3]->output(cout, 0, 1);
//	cout <<endl;
//#endif//DEBUG_CREX2
	// get the best
	ar = new alternative();
	ms = *(min_element( s.begin(), s.end() ));
//	random_shuffle( s.begin(), s.end() );
	for( unsigned i=0; i<s.size(); i++ ){
		if( s[i] == ms && ( all || ar->nralt() == 0 )){
			ar->insert( r[i] );
		}else{
			delete r[i];
		}
	}
	mr = ar->simplify();
	delete ar;

	delete I;
	delete T;
	delete IT;
}

void crex2::_crex2_lin_k3( itnode *n, const genom &g1, float &ms, rrrmt* &mr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst, bool all){

	int //sign = 1-n->sign,			// target sign = - sign of the node
			nsign = n->sign;		// inverse target sign
	vector<float> s(1, std::numeric_limits< float >::max());
	vector<unordered *> r(1, NULL);
	alternative *ar;
	rrrmt *tmp;
	rrrmt *I, *T, *IT;
		I = new rev( 0, 1, g1 );
		T = new transp( 0, 1, 2, g1 );
		IT = new revtransp(0, 1, 2, g1);
	// T + 2iT
	r[0] = new unordered();
	if(n->children.size() == 2 && 2*(*cst)[IT]+(*cst)[T]<=(*cst)[I]  ){ // TODO and 2*cst[RTRANM]+cst[TRANM]<=cst[REVNM]

		tmp = new transp(n->children[0]->i.first, n->children[0]->i.second+1, n->children[1]->i.second+1, g1);
		s[0] = (*cst)[tmp];
		r[0]->insert( tmp );

		tmp = new revtransp(n->children[0]->i.first, n->children[0]->i.second,
				n->children[1]->i.first, n->children[1]->i.second, g1);
		s[0] += (*cst)[tmp];
		r[0]->insert( tmp );

		tmp = new revtransp(n->children[1]->i.first, n->children[1]->i.second,
						n->children[0]->i.first, n->children[0]->i.second, g1);
		s[0] += (*cst)[tmp];
		r[0]->insert( tmp );

		s[0] += cdps[0][nsign] + cdps[1][nsign];
		r[0]->insert(cdpr[0][nsign]->clone());
		r[0]->insert(cdpr[1][nsign]->clone());
	}
//#ifdef DEBUG_CREX2
//	cout << "\tT+2iT "<<s[0]<<" "; r[0]->output(cout, 0, 1); cout <<endl;
//#endif//DEBUG_CREX2
	// get the best
	ar = new alternative();
	ms = *(min_element( s.begin(), s.end() ));
//	random_shuffle( s.begin(), s.end() );
	for( unsigned i=0; i<s.size(); i++ ){
		if( s[i] == ms && s[i] != std::numeric_limits< float >::max() && ( all || ar->nralt() == 0 )){
			ar->insert( r[i] );
		}else{
			delete r[i];
		}
	}
	mr = ar->simplify();
	delete ar;
	delete I;
	delete T;
	delete IT;
}

void crex2::_crex2_lin_all(itnode *n, const genom &g1, vector<float> &ms, vector<rrrmt*> &mr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst, bool all){

	int tsign;  // target sign
	genom id,   // quotient permutation
		tmpg;   // temp copy of the quotient permutation for applying rearrangements to it
	unordered *tu;  // tmp unordered scenario
	float ts;       // tmp score
	vector<rrrmt*> pr, prq;  // rearrangement wrt g and quotient permutation
	vector<alternative*> result(2, NULL);

	// construct the id
	id = genom( n->children.size(), 0 );

	// init the scores to maximum in order to be able to minimize
	ms[0] =std::numeric_limits< float >::max();
	ms[1] =std::numeric_limits< float >::max();

	// init the rearrangements as empty alternative (to be able to apply a delete)
	result[0] = new alternative();
	result[1] = new alternative();

	// get the potential rearrangements I, T, siT, piT
	// wrt the quotient permutation
	prq.push_back( new rev( 0, id.size()-1, id ) );
	if( id.size() == 2  ){
		prq.push_back( new transp(0, 1, 2, id) );
	}
	prq.push_back( new revtransp(1, id.size()-1, 0, 0, id) );
	prq.push_back( new revtransp(0, id.size()-2, id.size()-1,id.size()-1, id) );

	// get the potential rearrangements I, T, siT, piT
	// wrt the (real) permutation
	pr.push_back( new rev( abs((n->i).first), abs((n->i).second), g1 ) );
	if( n->children.size() == 2  ){
		pr.push_back( new transp(n->children[0]->i.first, n->children[0]->i.second+1, n->children[1]->i.second+1, g1) );
	}
	pr.push_back( new revtransp(n->children[1]->i.first, n->children.back()->i.second,
			n->children[0]->i.first, n->children[0]->i.second, g1) );
	pr.push_back( new revtransp(n->children[0]->i.first, n->children[n->children.size()-2]->i.second,
			n->children.back()->i.first, n->children.back()->i.second, g1) );

	// evaluate all subsets of the potential rearrangements
	for( int i=0; i<ppow( pr.size() ); i++ ){
		tsign = n->sign;
		tmpg = id;
		ts = 0;
		tu = new unordered();

		/*
		// for each subset of potential rearrangements
		// - apply them to the id permutation in order to determine the
		//   signs that need to be reached for the children
		// - get the score of the rearrangement
		// - get the target sign

		// Idea:
		// we have one of the 4 possible cases (each pic shows the given configuration (left)
		// and the target configuration on the right)
		//
		// left) odd number of rearrangement: target sign is != original node sign
		// right) even number of rearrangements: target sign is == original node sign
		// top) target sign is +
		// bottom) target sign is -
		//
		//    -  ===>  +         +  ===>  +
		//   / \      / \       / \      / \
		//  -...+    +...+     -...+    +...+
		//
		//    +  ===>  -         -  ===>  -
		//   / \      / \       / \      / \
		//  -...+    -...-     -...+    -...-
		//
		// assume we apply the set of rearrangements on the id that we would apply on the
		// quotient permutation. let pi be the result. then (for target sign +)
		// pi(i) > 0 indicates that the sign of the subtree at child i is
		//           not modified by the rearrangements
		// pi(i) < 0 indicates that the sign of the subtree at child i is
		//           switched by the rearrangements
		//
		// this is because if the quotient permutation is signed in this way the result
		// of applying the rearrangements will be "all +".
		// for target sign = -1 the orientation needs to be switched.
		 */


		for( int j=0; j<(int)pr.size(); j++ ){
			if( i & ppow(j) ){
				prq[j]->apply( tmpg );
				ts += (*cst)[ pr[j] ];
				tu->insert( pr[j]->clone() );
				tsign = 1-tsign;	// each rrrmt flips the target sign
			}
		}
//#ifdef DEBUG_CREX2
//		cout << "\ttesting "; tu->output(cout, 0, 1); cout << endl;
//#endif//DEBUG_CREX2
		// add score for the children
		for( unsigned j=0; j<tmpg.size(); j++ ){
			if( (tmpg[j]<0 && tsign == INC) || ( tmpg[j]>0 && tsign == DEC ) ){
				ts += cdps[ abs(tmpg[j])-1 ][DEC];
				tu->insert( cdpr[ abs(tmpg[j])-1 ][DEC]->simplify() );
			}else{ // (tmpg[j]>0 && tsign == INC) || ( tmpg[j]<0 && tsign == DEC )
				ts += cdps[ abs(tmpg[j])-1 ][INC];
				tu->insert( cdpr[ abs(tmpg[j])-1 ][INC]->simplify() );
			}
		}
//#ifdef DEBUG_CREX2
//		cout << "\t        "<<tsign<<" "<< ts<<" "; tu->output(cout, 0, 1); cout<< endl;
//#endif//DEBUG_CREX2
		if( ts > ms[tsign] ){
//#ifdef DEBUG_CREX2
//			cout << "\t        worse -> ignore"<<endl;
//#endif//DEBUG_CREX2
			delete tu;
			continue;
		}else if( ts < ms[tsign] ){
//#ifdef DEBUG_CREX2
//			cout << "\t        reset"<<endl;
//#endif//DEBUG_CREX2
			delete result[tsign];
			result[tsign] = new alternative();
			ms[tsign] = ts;
			result[tsign]->insert( tu );
		}else { // if( ts == ms[tsign] )
			result[tsign]->insert( tu );
		}
	}

	for(unsigned i=0; i<2; i++){
		if( all ){
			mr[i] = result[i];
		}else{
			mr[i] = result[i]->get_random_alternative();
			delete result[i];
		}
//#ifdef DEBUG_CREX2
//		cout << "==>"<< ms[i]<<" "; mr[i]->output(cout, 0, 1); cout<< endl;
//#endif//DEBUG_CREX2
	}

	// clean up memory for the rearrangements
	for(unsigned i=0; i<pr.size(); i++){
		delete pr[i];
		delete prq[i];
	}
}

void crex2::_crex2_lin( itnode *n, const genom &g1, vector<float> &dps, vector<rrrmt*> &dpr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst, bool all){

	alternative *ar;
	int sign = n->sign,			// sign of the node
		nsign = 1-sign;		// inverse of the sign of the node

	vector<float> sk01, sk23;
	vector<rrrmt*> rk01, rk23;

	sk01 = vector<float>(2, 0);
	rk01 = vector<rrrmt*>(2, NULL);
	sk23 = vector<float>(2, 0);
	rk23 = vector<rrrmt*>(2, NULL);

	_crex2_lin_k0(n, g1, sk01[sign], rk01[sign], cdps, cdpr);
//#ifdef DEBUG_CREX2
//	cerr << "kappa0 -> "<< sk01[sign]<<" ";
//	rk01[sign]->output(cerr, 0, 1);
//	cerr<< endl;
//#endif//DEBUG_CREX2
//	k1
	_crex2_lin_k1(n, g1, sk01[nsign], rk01[nsign], cdps, cdpr, cst, all);
//#ifdef DEBUG_CREX2
//	cerr << "kappa1 -> "<< sk01[nsign]<<" "; \
//	rk01[nsign]->output(cerr, 0, 1);
//	cerr<< endl;
//#endif//DEBUG_CREX2
//	k2
	_crex2_lin_k2(n, g1, sk23[sign], rk23[sign], cdps, cdpr, cst, all);
//#ifdef DEBUG_CREX2
//	cerr << "kappa2 -> "<< sk23[sign]<<" ";
//	rk23[sign]->output(cerr, 0, 1);
//	cerr<< endl;
//#endif//DEBUG_CREX2

//	k3
	_crex2_lin_k3(n, g1, sk23[nsign], rk23[nsign], cdps, cdpr, cst, all);
//#ifdef DEBUG_CREX2
//	cerr << "kappa3 -> "<< sk23[nsign]<<" "; rk23[nsign]->output(cerr, 0, 1); cerr<< endl;
//#endif//DEBUG_CREX2
	for(unsigned s=0; s<2; s++){
		dps[s] = min( sk01[s], sk23[s] );

		ar = new alternative();
		if( sk01[s] == dps[s] ){
			ar->insert( rk01[s] );
		}else{
			delete rk01[s];
		}
		if( sk23[s] == dps[s] && (all || sk01[s] > dps[s] ) ){
			ar->insert( rk23[s] );
		}else{
			delete rk23[s];
		}
		dpr[s]  = ar;
	}
//#ifdef DEBUG_CREX2
//	cerr << "lin 0" << dps[0]<<" "; dpr[0]->output(cerr, 0, 1); cerr<< endl; //<< *() << endl;
//	cerr << "lin 1" << dps[1]<<" "; dpr[1]->output(cerr, 0, 1); cerr<< endl; //<< *() << endl;
//#endif//DEBUG_CREX2

}


void crex2::_crex2_pri( itnode *n, const genom &g1, const genom &g, vector<float> &dps, vector<rrrmt*> &dpr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst, bool all, const unsigned time_bound, bool &is_optimal, const unsigned distance){

	// build the prime quotient permutation
	genom qperm;	// quotient permutation
	qperm = genom(n->children.size(), 0);

	for(unsigned k=0; k<n->children.size(); k++){
		qperm[k] = g[ abs(n->children[k]->i.first) ];
	}

	quotient_permutation(qperm, g.size());
	// for the SP$TYPE version, we need to set some signs: I choose the signs depending on the minimum cost
	for(unsigned k=0; k<cdps.size();k++){
		if(cdps[k][INC]<cdps[k][DEC]){
			//do nothing since qperm element is already positive
		}else if(cdps[k][INC]>cdps[k][DEC]){
			// make element negative
			qperm[k] *= -1;
		}else{
			//decide based on the SIT
			if(n->children[k]->sign == DEC){
				qperm[k] *= -1;
			}
		}
	}

	//
	// sort qperm to id and -id
	//
	genom id;					// identity permutation
	genom inv_id;				// inverse identity permutation
	rrrmt *sort_plu_id;			// sorting scenario qperm -> +id
	rrrmt *sort_min_id;			// sorting scenario qperm -> -id
	double plu_cost=0;			// cost for sorting qperm -> +id
	double min_cost=0;			// cost for sorting qperm -> -id

	//	build +id and -id permutation
	id = genom(qperm.size(), 0);
	inv_id = genom(qperm.size(), 0);
	reverse (inv_id.chromosom.begin(), inv_id.chromosom.end());

	for(int i=0; i<inv_id.size();i++){
		inv_id[i] *=-1;
	}

	// initialize the mapping from quotient permutation elements
	// to permutation elements
	vector< vector<int> > quotgmap;	// mapping from elements of the quotient permutation
	quotgmap = vector<vector<int> >( qperm.size()+1 );
	for( unsigned i=0; i<qperm.size(); i++ ){
		for(int j= abs(n->children[ i ]->i.first); j<=abs(n->children[ i ]->i.second); j++ ){
			quotgmap[ abs(qperm[i]) ].push_back( abs( g1[j] ) );
		}
	}

	// run SP4TYPE 
	vector<rrrmt*> cm(2);
	sort_by_4type(qperm, id, cm[0]);
	sort_by_4type(qperm, inv_id, cm[1]);
	float cnt_id=0.0;
	float cnt_inv_id=0.0;

	// we also have to add the cost of the child sortings
	float sorting_childs=0.0;
	for(int k=0; k<qperm.size();k++){
		if(qperm[k]<0){
			sorting_childs+=cdps[k][DEC];
		}else{
			sorting_childs+=cdps[k][INC];
		}
	}
	cnt_id=cm[0]->getmincost(cst)+sorting_childs;
	cnt_inv_id=cm[1]->getmincost(cst)+sorting_childs;


	#ifdef DEBUG_CREX2_PRIM
		cout << "min scen_cost id: " << cnt_id << endl;
		cout << "min scen_cost inv_id: " << cnt_inv_id << endl;
		cout << (*cm[0]) << endl;
		cout << (*cm[1]) << endl;
	#endif//DEBUG_CREX2_PRIM

	vector<rrrmt*> min_scen_plu;
	cm[0]->getminrrrmt(min_scen_plu,cst);
	ordered* scen_pl;
	scen_pl = new ordered(min_scen_plu, true);
	scen_pl->rename(quotgmap, g1.get_nmap(), g1.size());

	// add scenario
	dps[INC]=cnt_id;
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			dps[INC] += cdps[b][DEC];
		}else{
			dps[INC] += cdps[b][INC];
		}
	}
	unordered *final_sp4type_scenario;
	final_sp4type_scenario = new unordered();
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			final_sp4type_scenario->insert(cdpr[b][DEC]->clone());
		}else{
			final_sp4type_scenario->insert(cdpr[b][INC]->clone());
		}
	}
	final_sp4type_scenario->insert(scen_pl->clone());
	//correct mapping of scenario
	dpr[INC]=final_sp4type_scenario->clone();

	delete final_sp4type_scenario;
	delete scen_pl;

	// calculate parsimonious SP4TYPE scenario old
	vector<rrrmt*> min_scen_min;
	cm[1]->getminrrrmt(min_scen_min, cst);
	ordered* scen_min;
	scen_min = new ordered(min_scen_min, true);
	scen_min->rename(quotgmap, g1.get_nmap(), g1.size());

	// add scenario
	dps[DEC]=cnt_inv_id;
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			dps[DEC] += cdps[b][DEC];
		}else{
			dps[DEC] += cdps[b][INC];
		}
	}

	unordered *final_sp4type_minus_scen;
	final_sp4type_minus_scen = new unordered();
	for(int b=0; b < cdpr.size(); b++){
		if( qperm[b] < 0){
			final_sp4type_minus_scen->insert( cdpr[b][DEC]->clone() );
		}else{
			final_sp4type_minus_scen->insert( cdpr[b][INC]->clone() );
		}
	}
	final_sp4type_minus_scen->insert(scen_min->clone());

	dpr[DEC]=final_sp4type_minus_scen->clone();

	delete final_sp4type_minus_scen;
	delete scen_min;
	delete sort_min_id;
	delete sort_plu_id;
}

void crex2::_crex2_pri_approx(itnode *n, const genom &g1, const genom &g, vector<float> &ms, vector<rrrmt*> &mr,
		const vector<vector<float> >&cdps, const vector<vector<rrrmt*> > &cdpr,
		const costfoo *cst){
#ifdef DEBUG_CREX2_PRIM_APPROX
	cout << "CREX2_PRIM_APPROX CALLED FOR:" << endl;
	cout << "g1: " << g1 << endl;
	cout << "g: " << g << endl;
#endif//DEBUG_CREX2_PRIM_APPROX
	// build the prime quotient permutation
	genom qperm;	// quotient permutation
	qperm = genom(n->children.size(), 0);
	for(unsigned k=0; k<n->children.size(); k++){
		qperm[k] = g[ abs(n->children[k]->i.first) ];
	}
	quotient_permutation(qperm, g.size());

	// usually we have to test for all different sign combinations every value of qperm
	//	-> instead the weight minimum sign-combination is used
	for(unsigned k=0; k<cdps.size();k++){
		if(cdps[k][INC]<cdps[k][DEC]){
			//do nothing since qperm element is already positive
		}else if(cdps[k][INC]>cdps[k][DEC]){
			// make element negative
			qperm[k] *= -1;
		}else{
			//decide based on the SIT
			if(n->children[k]->sign == DEC){
				qperm[k] *= -1;
			}
		}
	}
#ifdef DEBUG_CREX2_PRIM_APPROX
	cout << "qperm: " << qperm << endl;
#endif//DEBUG_CREX2_PRIM_APPROX
	genom id;					// identity permutation
	genom inv_id;				// inverse identity permutation
	rrrmt *sort_plu_id;			// sorting scenario qperm -> +id
	rrrmt *sort_min_id;			// sorting scenario qperm -> -id
	double plu_cost=0;			// cost for sorting qperm -> +id
	double min_cost=0;			// cost for sorting qperm -> -id

	//	build +id and -id permutation
	id = genom(qperm.size(), 0);
	inv_id = genom(qperm.size(), 0);
	reverse (inv_id.chromosom.begin(), inv_id.chromosom.end());
	for(int i=0; i<inv_id.size();i++){
		inv_id[i] *=-1;
	}

	// initialize the mapping from quotient permutation elements
	// to permutation elements
	vector< vector<int> > quotgmap;	// mapping from elements of the quotient permutation
	quotgmap = vector<vector<int> >( qperm.size()+1 );
	for( unsigned i=0; i<qperm.size(); i++ ){
		for(int j= abs(n->children[ i ]->i.first); j<=abs(n->children[ i ]->i.second); j++ ){
			quotgmap[ abs(qperm[i]) ].push_back( abs( g1[j] ) );
		}
	}

//	cout << "quotgmap: " << endl;
//	for(int i=0; i<quotgmap.size();i++){
//		cout << "i: " << i << "-> ";
//		for(int j=0; j<quotgmap[i].size(); j++){
//			cout << quotgmap[i][j] << " ";
//		}
//		cout << endl;
//	}

	// calculate sorting scenarios
	sort_by_4type(qperm, id, sort_plu_id);
	sort_by_4type(qperm, inv_id, sort_min_id);

#ifdef DEBUG_CREX2_PRIM_APPROX
	cout << "SIT--> "; interval_tree_print(n, g1, cerr); cerr << endl;
	for( unsigned i=0; i<2; i++ ){
		for( unsigned j=0; j<cdps.size(); j++ ){
			cout << cdps[j][i]<< " ";
		}
		cout << endl;
	}

	cout << "scenario qperm ---> iota: "<< endl;
	cout << *sort_plu_id << endl;
	cout << endl;
	cout << "scenario qperm ---> -iota: "<< endl;
	cout << *sort_min_id << endl;
	cout << endl;
#endif//DEBUG_CREX2_PRIM_APPROX

	// rename elements of scenario according to quotient permutation mapping
	sort_plu_id->rename(quotgmap, g1.get_nmap(), g1.size());
	sort_min_id->rename(quotgmap, g1.get_nmap(), g1.size());

#ifdef DEBUG_CREX2_PRIM_APPROX
	cout << "scenario after renaming qperm ---> iota: "<< endl;
	cout << *sort_plu_id << endl;
	cout << endl;
	cout << "scenario after renaming qperm ---> -iota: "<< endl;
	cout << *sort_min_id << endl;
	cout << endl;
#endif//DEBUG_CREX2_PRIM_APPROX

	// calculate the weight of both scenarios using costfoo and store in ms[0]<->INC and ms[1]<->DEC
	//		-> weights for scenario to iota
	vector<rrrmt*> scen_plus;
	sort_plu_id->getrrrmt(scen_plus, ALTMIN);
	for(int z=0; z<scen_plus.size();z++){
		if(scen_plus[z]->type()==1 || scen_plus[z]->type()==4 || scen_plus[z]->type()==5 || scen_plus[z]->type()==6){
			cerr << "Internal error: unknown type of rrrmt: " << scen_plus[z]->type() << " ---> exiting!"<< endl;
			exit(EXIT_FAILURE);
		}
		ms[INC] += (*cst)[scen_plus[z]];
	}
	// add cost for child sorting scenarios depending on sign used for qperm
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			ms[INC] += cdps[b][DEC];
		}else{
			ms[INC] += cdps[b][INC];
		}
	}
	// 		-> weights for scenario to -iota
	vector<rrrmt*> scen_min;
	sort_min_id->getrrrmt(scen_min, ALTMIN);
	for(int z=0; z<scen_min.size(); z++){
		if(scen_min[z]->type()==1 || scen_min[z]->type()==4 || scen_min[z]->type()==5 || scen_min[z]->type()==6){
			cerr << "Internal error: unknown type of rrrmt: " << scen_min[z]->type() << " ---> exiting!"<< endl;
			exit(EXIT_FAILURE);
		}
		ms[DEC] += (*cst)[scen_min[z]];
	}
	// add cost for child sorting scenarios depending on sign used for qperm
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			ms[DEC] += cdps[b][DEC];
		}else{
			ms[DEC] += cdps[b][INC];
		}
	}

#ifdef DEBUG_CREX2_PRIM_APPROX
	cout << "weight for scenario tree rooted at qperm -> +sign: " << ms[INC] << endl;
	cout << "weight for scenario tree rooted at qperm -> -sign: " << ms[DEC] << endl;
	cout << endl;
#endif//DEBUG_CREX2_PRIM_APPROX

	// calculate scenarios
	unordered* final_plus_scenario = new unordered();
	unordered* final_min_scenario = new unordered();
	//			--> scenario for iota
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			final_plus_scenario->insert(cdpr[b][DEC]->clone());
		}else{
			final_plus_scenario->insert(cdpr[b][INC]->clone());
		}
	}
	final_plus_scenario->insert(sort_plu_id);
	mr[INC]=final_plus_scenario->clone();
	//			--> scenario for -iota
	for(int b=0; b < cdps.size(); b++){
		if(qperm[b]<0){
			final_min_scenario->insert(cdpr[b][DEC]->clone());
		}else{
			final_min_scenario->insert(cdpr[b][INC]->clone());
		}
	}
	final_min_scenario->insert(sort_min_id);
	mr[DEC]=final_min_scenario->clone();

#ifdef DEBUG_CREX2_PRIM_APPROX
	cout << "scenario for tree rooted at qperm -> +sign: " << endl;
	cout << *mr[INC] << endl;
	cout << "scenario tree rooted at qperm -> -sign: " << endl;
	cout << *mr[DEC] << endl;
	cout << endl;
#endif//DEBUG_CREX2_PRIM_APPROX

//	delete final_min_scenario;
//	delete final_plus_scenario;
}

void crex2::_crex2_main( itnode *n, const genom &g1, const genom &g, vector<float> &dps, vector<rrrmt*> &dpr,
		const costfoo *cst, bool complete, bool all, const unsigned time_bound, bool &is_optimal, const unsigned distance, const bool approx){

	// determine cost of a leaf node
	if( n->children.size() == 0 ){
		_crex2_leaf( n, g1, dps, dpr, cst );
		return;
	}

	// DP scores and rrrmts of the childs
	// cdps[i][j]: score of child i and sign j (j either 0=INC or 1=DEC)
	// cdpr[i][j]: corresponding rearrangements for this child and sign
	vector<vector<float> > cdps;
	vector<vector<rrrmt*> > cdpr;

	cdps = vector<vector<float> >( n->children.size(), vector<float>(2, 0) );
	cdpr = vector<vector<rrrmt*> >( n->children.size(), vector<rrrmt*>(2,NULL) );
	for( unsigned i=0; i<n->children.size(); i++ ){
		_crex2_main( n->children[i], g1, g, cdps[i], cdpr[i], cst, complete, all, time_bound, is_optimal, distance, approx);
	}

#ifdef DEBUG_CREX2
	cout << "--> "; interval_tree_print(n, g1, cerr); cerr << endl;
	cout << "scores:" <<endl;
	for( unsigned i=0; i<2; i++ ){
		if(i==0){
			cout << "+ ";
		}else{
			cout << "- ";
		}
		for( unsigned j=0; j<cdps.size(); j++ ){
			cout << cdps[j][i]<< " ";
		}
		cout << endl;
	}
	cout << "scenarios:" << endl;
	for( unsigned j=0; j<cdps.size(); j++ ){
		for( unsigned i=0; i<2; i++ ){
			cout << "c(" << j<<") s(";
			if(i==0){
				cout << "+) ";
			}else{
				cout << "-) ";
			}
			cdpr[j][i]->output(cout, 0, 1);
			cout << endl;
		}
		cout << endl;
	}
#endif//DEBUG_CREX2
	if( n->type == LIN ){
		if(complete){
			_crex2_lin_all( n, g1, dps, dpr, cdps, cdpr, cst, all );
		}else{
			_crex2_lin( n, g1, dps, dpr, cdps, cdpr, cst, all );
		}
	}else{
		if(approx){// chose approximation algorithm instead of exact one
			_crex2_pri_approx(n, g1, g, dps, dpr, cdps, cdpr, cst);
		}else{// calculate exact solutions or use timebound
			_crex2_pri( n, g1, g, dps, dpr, cdps, cdpr, cst, all, time_bound, is_optimal, distance);
		}
	}

	// child data not needed anymore -> free
	for( unsigned i=0; i<cdpr.size(); i++ ){
		delete cdpr[i][0];
		delete cdpr[i][1];
	}
}

string crex2::typestr() const{
	return "crex2";
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// construct an emptyscen scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crex::crex() : unordered(){
	_complete = false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// get the crex scenario for a pair of genomes
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crex::crex(const genom &src, const genom &tgt, bool bps, unsigned maxszen, bool mkalt) : unordered(){
	bool tmatch, rtmatch;// did the (reverse) transposition pattern match?
	genom srcp, tgtp;	// genomes with equalized gene content
	int n;			//
	itnode *iroot;		// interval tree relative to g1
	list<itnode *> nodes; 	// nodes of the interval tree
	rrrmt *tmp = NULL;
	pair<vector<rrrmt*>,vector<rrrmt*> > indels; 	// temporary storage for indels
	unsigned mxel = 0; 	// maximum element;
	vector<vector<int> > rnm;	// mapping from new to old name (for back renaming)
	vector<string> orgnmap;	// temporary nmap for that will apply to srcp and tgtp
	vector<string> *nmap;	//

#ifdef DEBUG_CREX
	cout << "crex::crex():"<<endl;
	cout << ">src "<<endl<<src<<endl;
	cout << ">tgt "<<endl<<tgt<<endl;
#endif//DEBUG_CREX

	if( src.size() < 2 || tgt.size() < 2 ){
		cerr << "internal error: crex received a genome of size < 2"<<endl;
		cerr << "  src: "<<src<<endl;
		cerr << "  tgt: "<<tgt<<endl;
		exit(EXIT_FAILURE);
	}

	// determine indels and get a genome with same gene set
	_crex_indels( src, tgt, srcp, tgtp, indels, rnm, orgnmap );

//	srcp = src;
//	tgtp = tgt;
	n = srcp.size();

//	cerr << "START " << endl << srcp << endl << tgtp << endl;

		// assign the genome relative to which the intervals are computed
//	genome = g1;
		// get the interval tree and a list of its nodes (with the leaves)
	interval_tree(srcp, tgtp, n, &iroot);
#ifdef DEBUG_CREX
	interval_tree_print(iroot, srcp, cout); cout << endl;
#endif//DEBUG_CREX

	std::insert_iterator<std::list<itnode_struct*> > x = inserter(nodes, nodes.begin());
	interval_tree_nodes(iroot, x, true);


	// construct the the scenario in a single iteration over the nodes
	for( list<itnode*>::iterator it=nodes.begin(); it!=nodes.end(); it++ ){
//		#ifdef DEBUG_CREX
//		cout << "   node ["<<(*it)->i.first<<","<<(*it)->i.second<<"]"<<endl;
//		#endif//DEBUG_CREX
		if( (*it)->type == LIN ){

			tmatch = rtmatch = false;

			// try to match the transposition pattern
			tmatch = _crex_transpositions(*it, srcp);
			// try to match the reverse transposition pattern
			// (not necessary if the transposition pattern already matched)
			if( !tmatch )
				rtmatch = _crex_reversetranspositions(*it, srcp);

			// if nothing helps the reversal pattern has to be applied
			// - to the children
			// - and also to root nodes
			if( !tmatch && !rtmatch ){
				for( unsigned i=0; i<(*it)->children.size(); i++ )
					_crex_reversals((*it)->children[i], srcp);
			}

			if( (*it)->parent==NULL )
				_crex_reversals(*it, srcp);
		}else{
			_crex_tdls(*it, srcp, tgtp, bps, maxszen);
		}
	}

	// make alternatives for some rearrangements
	if( mkalt == true ){
		tmp = mkealt();
		if( tmp->isempty() ){
			delete tmp;
		}else{
			cerr << "internal error: crex:crex: mkealt() did not return empty"<<endl;
			exit(EXIT_FAILURE);
		}
	}


		// determine the maximum element
	for( unsigned i=0; i<src.size(); i++ ){
		mxel = max(mxel, (unsigned)abs(src[i]));
	}
	for( unsigned i=0; i<tgt.size(); i++ ){
		mxel = max(mxel, (unsigned)abs(tgt[i]));
	}
	this->rename( rnm, src.get_nmap(), mxel );

	// if there are (single gene) insertions or deletion
	// we need to ensure that there come first resp. last
	// so we embed the crex scenario into an ordered scenario
	// framed by the deletion and insertions
	if( (indels.first.size() + indels.second.size()) > 0){

//		cerr << "RECoNSTRUCTED "<<endl << *this<<endl;

		ordered *tmpo = new ordered();
		for(unsigned i=0; i<indels.second.size(); i++){
			tmpo->push_back(  indels.second[i] );
		}

		tmp = this->simplify();
		if( tmp->isempty() ){
			delete tmp;
		}else{
			tmpo->push_back( tmp );
		}

		for(unsigned i=0; i<indels.first.size(); i++){
			tmpo->push_back( indels.first[i] );
		}
		for( set<rrrmt*, HDereferenceLess>::iterator it=_sce.begin(); it!=_sce.end(); it++ ){
			delete *it;
		}
		_sce.clear();

		tmp = tmpo->simplify();
		delete tmpo;
		insert( tmp );



	}
	// restore nmap
	nmap = src.get_nmap();
	for(unsigned i=0; i<orgnmap.size(); i++){
		(*nmap)[i] = orgnmap[i];
	}

#ifdef DEBUG_CREX
	cerr << "reconstruction "<<endl << *this<<endl;
#endif//DEBUG_CREX
	// check the reconstruction


	genom t = src;
//	cerr << "apply "<<t<<endl;
	apply(t);
	if (tgt != t){
		cerr << "crex: apply failed"<<endl;
		cerr << ">tgt"<<endl<<tgt<<endl<<">src"<<endl<<src<<endl;
		cerr << ">t"<<endl<<t<<endl;
		cerr << *this<<endl;
		exit(EXIT_FAILURE);
	}

		// check the reconstruction
	t = tgt;
//	cerr << "deapply"<<t<<endl;
	deapply(t);
	if (src != t){
		cerr << "crex: deapply failed"<<endl;
		cerr << ">tgt"<<endl<<tgt<<endl<<">src"<<endl<<src<<endl;
		cerr << ">t"<<endl<<t<<endl;
		cerr << *this<<endl;
		exit(EXIT_FAILURE);
	}


	interval_tree_free(iroot);

//#ifdef DEBUG_CREX
//	cout << *this <<endl;
//#endif//DEBUG_CREX

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// copy constructor of crex scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crex::crex(const crex &c ) : unordered(c){
	_complete = c._complete;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// destructor for a crex scenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crex::~crex(){
	_complete = false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

crex* crex::clone() const{
	return new crex(*this);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void crex::_crex_indels( const genom &src, const genom &tgt,
			genom &srcp, genom &tgtp, pair<vector<rrrmt*>, vector<rrrmt*> > &indels,
			vector<vector<int> > &rnm, vector<string> &orgnmap){

	int n, N;
	rrrmt *tmp;	// tmp for inserting
	vector<int> _srcinv, _tgtinv;	// inverse permutations
	vector<string> *nmap = src.get_nmap();

	if( nmap != NULL ){
		orgnmap =  (* src.get_nmap() );
	}

	// equalize gene content and add corresponding single
	// gene indels to the scenario. take care that the resulting
	// permutations are 1..x, i.e. no missing element
	_srcinv = src.inverse();
	_tgtinv = tgt.inverse();

//	cerr << "inverser permutations "<<endl;
//	copy(_srcinv.begin(), _srcinv.end(), ostream_iterator<int>(cerr," "));cerr<<endl;
//	copy(_tgtinv.begin(), _tgtinv.end(), ostream_iterator<int>(cerr," "));cerr<<endl;

	// 1. check for elements from src missing in tgt
	//    these elements are deleted on the way from src to tgt

	srcp = src;
	for( unsigned i=0; i < src.size(); i++ ){
		if( (int)_tgtinv.size() < abs(src[i])+1 || _tgtinv[abs(src[i])] == std::numeric_limits< int >::max() ){
			tmp = new indel( src[i], true, srcp );
			indels.second.push_back(tmp);
#ifdef DEBUG_CREX_INDELS
			cout <<  src[i]<<" missing in tgt -> adding  "<< *tmp <<endl;
#endif//DEBUG_CREX_INDELS
			tmp->apply( srcp );
		}

//		else{
//			srcp.chromosom.push_back( src[i] );
//		}
	}
	// 2. check for elements from tgt missing in src
	//    these elements are inserted on the way from src to tgt
	tgtp = tgt;
	for( unsigned i=0; i < tgt.size(); i++ ){
		if( (int)_srcinv.size() < abs(tgt[i])+1 || _srcinv[abs(tgt[i])] == std::numeric_limits< int >::max() ){
			tmp = new indel( tgt[i], false, tgtp );
			indels.first.push_back(tmp);
#ifdef DEBUG_CREX_INDELS
			cout <<  tgt[i]<<" missing in src -> adding  "<< *tmp <<endl;
#endif//DEBUG_CREX_INDELS
			tmp->deapply( tgtp );
		}
//		else{
//			tgtp.chromosom.push_back( tgt[i] );
//		}
	}
	std::reverse( indels.first.begin(), indels.first.end() );


	if( srcp.size() != tgtp.size() ){
		cerr << "internal error: gene content not equalized"<<endl;
		cerr << "  src: "<<src<<endl;
		cerr << "  tgt: "<<tgt<<endl;
		cerr << "  src': "<<srcp<<endl;
		cerr << "  tgt': "<<tgtp<<endl;
		exit(EXIT_FAILURE);
	}

//#ifdef DEBUG_CREX_INDELS
//	cerr << "equalized gene content"<<endl;
//	cerr << srcp<<endl<<tgtp<<endl;
//#endif//DEBUG_CREX_INDELS

	n = srcp.size();

	// rename the elements of src and tgt such that we have a
	// permutation of {1,...,n}
	N = max( _srcinv.size(), _tgtinv.size() );


	rnm = vector<vector<int> >(n+1, vector<int>());
//	cerr << "N "<<N<<endl;

	for( unsigned o=1, nn=1; o < (unsigned)N; o++ ){
		if( ( o < _srcinv.size() && _srcinv[o] != std::numeric_limits< int >::max()) &&
				( o < _tgtinv.size() && _tgtinv[o] != std::numeric_limits< int >::max() ) ){
//			cerr << "o "<<o<<" nn "<<nn<<endl;
			rnm[nn].push_back( o );
			if( nmap != NULL && orgnmap.size()>0){
				(*nmap)[nn] = (orgnmap)[o];
			}
			for( unsigned i=0; i<srcp.size(); i++ ){
				if( (unsigned)abs(srcp[i]) == o ){
					srcp.chromosom[i] = (srcp[i]<0) ? -nn : nn;
				}
				if( (unsigned)abs(tgtp[i]) == o ){
					tgtp.chromosom[i] = (tgtp[i]<0) ? -nn : nn;
				}
			}
			nn++;
		}
	}

#ifdef DEBUG_CREX_INDELS
	for(unsigned i=0; i<rnm.size(); i++){
		cerr << "rename("<<i<<") ";
		copy( rnm[i].begin(), rnm[i].end(), ostream_iterator<int>(cerr, " ") ); cerr<<endl;
	}
#endif//DEBUG_CREX_INDELS

	if( src.get_nmap() != NULL ){
		srcp.set_nmap( nmap );
		tgtp.set_nmap( nmap );
	}else{
		srcp.set_nmap( NULL );
		tgtp.set_nmap( NULL );
	}
#ifdef DEBUG_CREX_INDELS
	if( nmap != NULL ){
		cerr << "new name maps"<<endl;
		copy( nmap->begin(), nmap->end(), ostream_iterator<string>(cerr, ",") );cerr <<endl;
	}
	cerr << "equalized gene content"<<endl;
	cerr << srcp<<endl<<tgtp<<endl;
#endif//DEBUG_CREX_INDELS

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// find reversal patterns in the remaining nodes of the interval tree
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool crex::_crex_reversals(itnode* node, const genom &g1){
	rev *tmp;

#ifdef DEBUG_CREXREV
	cout << "crex_reversals"<<endl;
#endif//DEBUG_CREXREV

	if( node->type == PRI || (node->parent != NULL && node->parent->type == PRI)){
		return false;
	}

	// negative root node or nodes with sign different from the linear parent have to be reversed
	if( !( ( node->parent == NULL && node->sign == DEC ) || ( node->parent != NULL &&  node->parent->sign != node->sign )) ){
		return false;
	}
//	cout << "crex found reversal ["<<abs(((*it)->i).first)-1<<","<<abs(((*it)->i).second)-1<<"]"<<endl;
	// construct the reversal
	tmp = new rev( abs((node->i).first), abs((node->i).second), g1);
#ifdef DEBUG_CREXREV
	cout << " adding  "<< *tmp <<endl;
#endif//DEBUG_CREXREV
	// apply the reversal to g1 and to the tree
	insert(tmp);


	return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// find reverse transposition patterns in the remaining nodes of the interval tree
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool crex::_crex_reversetranspositions(itnode* node, const genom &g1){

	int ptype,	// type and ...
		psign;  // ... sign of the parent
	pair<int, int> r;
	unsigned nesign = 0,	// number of sign differences between the current node and its children
		nepos = 0;	// last position where a sign difference occured
//	vector< vector<itnode *> > newchildren(2);
	revtransp *tmp = NULL;

	#ifdef DEBUG_CREXRTRA
	cout << "crex_reversetranspositions"<<endl;
	#endif//DEBUG_CREXRTRA

//	for(list<itnode *>::const_iterator it = node.begin(); it != node.end(); ){
			// only look for linear nodes
	if( node->type == PRI){
		return false;
	}
		// get the parent sign and type
	if( node->parent == NULL || (node->parent->type == PRI && node->parent->sign == UNK ) ){
		ptype = LIN;
		psign = INC;
	}else{
		ptype = node->parent->type;
		psign = node->parent->sign;
	}

	// nodes must have a sign different from the parent
	// TODO WHY sing == INC
	if( ptype == UNK || node->sign == INC || psign == node->sign ){
		#ifdef DEBUG_CREXRTRA
		cout << " -> no rT "<<ptype<<" "<<psign<<endl;
		#endif//DEBUG_CREXRTRA
		return false;
	}
	#ifdef DEBUG_CREXRTRA
	cout << " psign  "<<psign<<" childs ";
	#endif//DEBUG_CREXRTRA

	// determine the number of child nodes that have a sign
	// different from the sign of the node; and store the position of the
	// last child with a different sign
	for( unsigned j=0; j<node->children.size(); j++ ){
		#ifdef DEBUG_CREXRTRA
		cout << node->children[j]->sign;
		#endif//DEBUG_CREXRTRA
		if( node->children[j]->type == PRI){
			continue;
		}
		if( node->sign != node->children[j]->sign ){
			nesign++;
			nepos = j;
		}
	}
	#ifdef DEBUG_CREXRTRA
		cout << endl;
	#endif//DEBUG_CREXRTRA

	// a reverse transposition must have exactly one child that has a different sign
	// this child must be the first or last child
	if( !(nesign == 1 && (nepos == 0 || nepos == node->children.size()-1 ))){
		#ifdef DEBUG_CREXRTRA
		cout << " -> no rT nesign "<<nesign<<" nepos "<<nepos<<endl;
		#endif//DEBUG_CREXRTRA
		return false;
	}
		// put the first child and all with the same sign into newchildren[0]
		// the remaining into newchildren[1]
//	for(unsigned j=0; j<node->children.size(); j++){
//
//
//		if( newchildren[0].size() == 0 || newchildren[0][0]->sign == node->children[j]->sign ){
//			newchildren[0].push_back( node->children[j] );
//		}else{
//			newchildren[1].push_back( node->children[j] );
//		}
//	}

		// save the reverse transposition
		// the child 0 is the transposed, all remaining children are reverse transposed part
	if ( nepos == 0 ){
		tmp = new revtransp(node->children[1]->i.first,
				node->children.back()->i.second,
				node->children[0]->i.first,
				node->children[0]->i.second, g1);
	}
		// the last child is the transposed and all remaining children are the reverse transposed part
	else{
		tmp = new revtransp(node->children[0]->i.first,
				node->children[node->children.size()-2]->i.second,
				node->children.back()->i.first,
				node->children.back()->i.second, g1);
	}
	#ifdef DEBUG_CREXRTRA
		cout << " adding  "<< *tmp <<endl;
	#endif//DEBUG_CREXRTRA
	insert( tmp );

////			interval_tree_print(*it, g1, cout); cout << endl;
//		// reverse the reverse-transposition nodes
//	for (unsigned j=0; j<newchildren.size(); j++){
//		if( newchildren[j][0]->sign != psign ){
//			for (unsigned k=0; k< newchildren[j].size(); k++){
//				interval_tree_reverse( newchildren[j][k] );
//			}
//			reverse( newchildren[j].begin(), newchildren[j].end() );
//		}
//	}
//	swap( newchildren[0], newchildren[1] );
//
//		// insert the nodes in the new (transposed) order
//	node->children.clear();
//	node->children.insert(node->children.end(), newchildren[0].begin(), newchildren[0].end());
//	node->children.insert(node->children.end(), newchildren[1].begin(), newchildren[1].end());
//
		// correct the sign of the node
	if( node->sign == INC ){
		node->sign = DEC;
	}else if( node->sign == DEC ){
		node->sign = INC;
	}

////			interval_tree_print(*it, g1, cout); cout << endl;
//	newchildren[0].clear();
//	newchildren[1].clear();
	return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// find transposition patterns in the remaining nodes of the interval tree
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool crex::_crex_transpositions( itnode* node, const genom &g1 ){

	pair<int,int> r1, r2;
	itnode *c0, *c1;
	transp *tmp = NULL;

		// only linear nodes with two children can be transpositions ..
	if( node->children.size() != 2 || node->type != LIN ){
		return false;
	}

	c0 = node->children[0];
	c1 = node->children[1];

#ifdef DEBUG_CREXTRA
	cout << "crex_transpositions"<<endl;
	cout << "     signs  n:"<<((node->sign == INC )?'+':'-') <<" children:"<< ((c0->sign == INC )?'+':'-')<<" "<<((c1->sign == INC )?'+':'-')<<endl;
#endif//DEBUG_CREXTRA

	if( ! (( node->sign == DEC && c0->sign == INC && c1->sign == INC ) ||
		( node->sign == INC && c0->sign == DEC && c1->sign == DEC ) ) ) {
#ifdef DEBUG_CREXTRA
		cout << "     no transposition"<<endl;
#endif//DEBUG_CREXTRA


		return false;
	}
	r1 = node->children[0]->i;
	r1.first = abs(r1.first);
	r1.second = abs(r1.second);
	r2 = node->children[1]->i;
	r2.first = abs(r2.first);
	r2.second = abs(r2.second);

	tmp = new transp(r1.first, r1.second+1, r2.second+1, g1);
#ifdef DEBUG_CREXTRA
		cout << "       add  "<< *tmp <<endl;
#endif//DEBUG_CREXTRA

	insert( tmp );

	if ( node->sign == INC ){
		node->sign = DEC;
	}else{
		node->sign = INC;
	}
//		if( c1->type == PRI ){
//			c1->sign = (*it)->sign;
//		}
//		if( c0->type == PRI ){
//			c0->sign = (*it)->sign;
//		}
//		swap((*it)->children[0], (*it)->children[1]);

	return true;

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// find tdrl patterns in the remaining nodes of the interval tree
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool crex::_crex_tdls(itnode* node, const genom &g1, const genom &g2, bool bps, unsigned maxszen ){
	alternative *asce = NULL;	// the resulting alternative scenarios
	genom g,			// identified genome
		quot,			// quotient permutation for a prime node
		tgt;
	int n = g1.size(),	// length of the genomes
		bestd,			// number of rearrangements in the shortest alternative
		tdld,
		psign;			// the sign of the parent
	map<genom, bool> used;			// for avoiding recomputations in the reversal process (temporary memory for reversal function)
	rrrmt *tmp;
	vector< pair<int,int> > blocks,	// negative/posiztive blocks of the qoutient permutation of a prime node
		blocks_inv,
		currz;							// temporary memory for the reversal funtion
	vector<vector<bool> > tdls;			// one tdl szenario
	vector< vector<int> > quotgmap;	// mapping from elements of the quotient permutation
										// to elements of the permutation (0 is empty!)
	vector<vector<pair<int, int> > > revscen;	// reversal szenarios leading to the positive / negative  permutations
	vector<vector<vector<bool> > > tdlscen; 	// tdl szenarios (each tdl is a bool vector, a szenario contains
												// possibly several tdls and we store all szenarios ->3x vector)

	// only look at pnodes
	if(node->type != PRI){
		return false;
	}

	// identify the permutations
	g = g2.identify_g(g1);

	#ifdef DEBUG_CREXTDRL
	cout << "crex_tdls"<<endl;
	cout <<"g: "<< g << endl;
	#endif//DEBUG_CREXTDRL

	// first alternative scenarios are constructed for the quotient
	// permutations, the scenarios for g1 are derived afterwards
	asce = new alternative();

	// initialise the:
	// - best distance to std::numeric_limits< int >::max()
	// - the quotient permutation
	// - the sign of the parent
	bestd = std::numeric_limits< int >::max();

	quot = genom(node->children.size(), 0);
	for(unsigned k=0; k<node->children.size(); k++){
		quot[k] = g[ abs(node->children[k]->i.first) ];
	}
	quotient_permutation(quot, n);
	for(unsigned k=0; k<node->children.size(); k++){
		if(node->children[k]->sign == DEC)
			quot[k] *= -1;
	}
		// initialise the mapping from quotient permutation elements
		// to permutation elements
	quotgmap = vector<vector<int> >( quot.size()+1 );
	for( unsigned i=0; i<quot.size(); i++ ){
		for(int j= abs(node->children[ i ]->i.first); j<=abs(node->children[ i ]->i.second); j++ ){
			quotgmap[ abs(quot[i]) ].push_back( abs( g1[j] ) );
		}
#ifdef DEBUG_CREXTDRL
		cout << "qm("<<abs(quot[i])<<") = "; copy(quotgmap[ abs(quot[i]) ].begin(), quotgmap[ abs(quot[i]) ].end(), ostream_iterator<int>(cout," ")); cout << endl;
#endif//DEBUG_CREXTDRL
	}

	if( node->parent != NULL  && node->parent->sign == DEC ){
		tgt = genom_nid(quot.size());
		psign = DEC;
	}else{
		tgt = genom(quot.size(), 0 );
		psign = INC;
	}

#ifdef DEBUG_CREXTDRL
	cout << "quot "<<quot<<endl;
	cout << "psig "<<psign<<endl;
#endif//DEBUG_CREXTDRL

	if( bps ){
		bpscenario *bps = NULL;
		set<rrrmt*, HDereferenceLess> bpa;
		bps = new bpscenario( quot, tgt, maxszen, true );
		if( !bps->isempty() ){
			bpa = bps->get_alternatives();
			for( set<rrrmt*, HDereferenceLess>::const_iterator it=bpa.begin(); it!=bpa.end(); it++ ){
				asce->insert( *it );
			}
			bpa.clear();
			delete bps;
//			asce->insert( bps );
		}else{
			delete bps;
		}
	}

	// get the positive blocks of the quotient permutation
	get_blocks(quot, tgt, blocks);
	get_blocks(tgt, quot, blocks_inv);
#ifdef DEBUG_CREXTDRL
	cout << "TO "<<tgt<<endl;
	cout << "blocks    : ";	for(unsigned j=0; j<blocks.size(); j++){cout << "("<<blocks[j].first<<","<<blocks[j].second<<") ";}cout << endl;
	cout << "blocks_inv: "; for(unsigned j=0; j<blocks_inv.size(); j++){cout << "("<<blocks_inv[j].first<<","<<blocks_inv[j].second<<") ";}cout << endl;
#endif//DEBUG_CREXTDRL
	if( blocks.size() == 0 ){
		tdl_sort( quot, tgt, tdls );
		add_quot_rrrmt(tdls, quot, tgt, asce);
	}else{
		tdld = std::numeric_limits< int >::max();
		reverse_blocks(quot, tgt, blocks, tdlscen, revscen, used, maxszen, currz, tdld);
		#ifdef DEBUG_CREXTDRL
		cout << "        -> ("<<blocks.size()<<" rev + "<<tdld<<" tdrls)"<<endl;
		#endif//DEBUG_CREXTDRL
		if(tdld + (int)blocks.size() <= bestd){
			add_quot_rrrmt(tdlscen, revscen, quot, tgt, true, asce);
			bestd = tdld + blocks.size();
		}
		currz.clear(); used.clear(); tdlscen.clear(); revscen.clear();

		/*tdld = std::numeric_limits< int >::max();
		reverse_blocks_inv(quot, tgt, blocks_inv, tdlscen, revscen, used, maxszen, currz, tdld);
		#ifdef DEBUG_CREXTDRL
		cout << "            -> ("<<blocks_inv.size()<<" rev + "<<tdld<<" tdrls)"<<endl;
		#endif//DEBUG_CREXTDRL
		if(tdld + (int)blocks_inv.size() <= bestd){
			add_quot_rrrmt(tdlscen, revscen, quot, tgt, false, asce);
			bestd = tdld + blocks.size();
		}
		currz.clear(); used.clear(); tdlscen.clear(); revscen.clear();*/
	}
	blocks.clear();
	blocks_inv.clear();

//	cerr <<"pnode scenario "<<*asce<<endl;

		// up to now alternative rearrangement sceanarios for the quotient
		// permutation have been constructed. by renaming the elements that are
		// affected by the rearrangements with the element sets  corresponding
		// to the child nodes of the pnode scenarios for the permutations
		// are gained

	asce->rename( quotgmap, g1.get_nmap(), g1.size() );

		// finally simplify and insert the scenario
	tmp = asce->simplify();
	delete asce;
	insert( tmp );

	// apply the consequences of the rearrangement
	// i.e. make the node linear and fix the sign
	node->type = LIN;
	node->sign = psign;

	quot.clear();
	return true;
}

/*
void crex::crex_tdls(list<itnode *> node, const genom &g1, const genom &g2, bool bps, unsigned maxszen ){

	alternative *asce = NULL;	// the resulting alternative scenarios
	bpscenario *bpsce = NULL;	// the scenario computed with the breakpoint method
	rrrmt *tmp;

	genom g,		// identified genome
		quot,		// quotient permutation for a prime node
		id,			// identity and ..
		nid;		// .. negative identity
	int n = g1.size(),	// length of the genomes
		psign,			// the sign of the parent
		tdld,
		bestd;
	map<genom, bool> used;			// for avoiding recomputations in the reversal process (temporary memory for reversal function)
	vector< pair<int,int> > blocks,	// negative/posiztive blocks of the qoutient permutation of a prime node
		blocks_inv,
		currz;						// temporary memory for the reversal funtion
	vector<vector<pair<int, int> > > revscen;	// reversal szenarios leading to the positive / negative  permutations
	vector<vector<bool> > tdls;					// one tdl szenario
	vector<vector<vector<bool> > > tdlscen; 	// tdl szenarios (each tdl is a bool vector, a szenario contains
												// possibly several tdls and we store all szenarios ->3x vector)
	vector< vector< int > > quotgmap;	// mapping from elements of the quotient permutation
										// to elements of the permutation (0 is empty!)

		// identify the permutations
	g = g2.identify_g(g1);
#ifdef DEBUG_CREXTDRL
	cout <<"g: "<< g << endl;
#endif//DEBUG_CREXTDRL
		// search the nodes
	for(list<itnode*>::iterator it = node.begin(); it!= node.end(); it++){
			// only look at pnodes
		if((*it)->type != PRI){
			continue;
		}

		asce = new alternative();
			// initialise the
			// - best distance to std::numeric_limits< int >::max()
			// - the quotient permutation
			// - the identity and negative identity (of len(quot) length)
			// - the positive and negative blocks
			// - the sign of the parent
		bestd = std::numeric_limits< int >::max();

		quot = genom((*it)->children.size(), 0);
		for(unsigned k=0; k<(*it)->children.size(); k++){
			quot[k] = g[ abs((*it)->children[k]->i.first)-1 ];
		}

		quotient_permutation(quot, n);
		for(unsigned k=0; k<(*it)->children.size(); k++){
			if((*it)->children[k]->sign == DEC)
				quot[k] *= -1;
		}
			// initialise the mapping from quotient permutation elements
			// to permutation elements
		quotgmap = vector<vector<int> >( quot.size()+1 );
		for( unsigned i=0; i<quot.size(); i++ ){
//			cout << "("<< ((*it)->children[ i ]->i.first-1) <<","<< (*it)->children[ i ]->i.second<<") "<<g1.size()<<endl;
			for(int j= abs((*it)->children[ i ]->i.first)-1; j<abs((*it)->children[ i ]->i.second); j++ ){
				quotgmap[ abs(quot[i]) ].push_back( abs( g1[j] ) );
			}
//			cout << "qm("<<abs(quot[i])<<") = "; copy(quotgmap[ abs(quot[i]) ].begin(), quotgmap[ abs(quot[i]) ].end(), ostream_iterator<int>(cout," ")); cout << endl;
		}

		id = genom(quot.size(), 0 );
		nid = genom_nid(quot.size());

		if( (*it)->parent != NULL ){
			psign = (*it)->parent->sign;
			cout << "take parent sign "<<endl;
		}else{
			psign = INC;
		}

#ifdef DEBUG_CREXTDRL
		cout << "quot "<<quot<<endl;
		cout << "psig "<<psign<<endl;
#endif//DEBUG_CREXTDRL

			// try all possibilities to reverse the blocks (with respect to the sign of the parent node)
			// and select outcomes with minimal tdl distance (stored in pbestd)
			// if the sign is unknown -> compute both
		if( psign == INC || psign == UNK ){

			if( bps ){
				bpsce = new bpscenario( quot, id, maxszen, true );
					if( !bpsce->isempty() ){
						asce->insert( bpsce );
					}else{
						delete bpsce;
				}
			}

			// get the positive blocks of the quotient permutation
			get_blocks(quot, id, blocks);
			get_blocks(id, quot, blocks_inv);
#ifdef DEBUG_CREXTDRL
			cout << "TO POSITIVE"<<endl;
			cout << "blocks    : ";	for(unsigned j=0; j<blocks.size(); j++){cout << "("<<blocks[j].first<<","<<blocks[j].second<<") ";}cout << endl;
			cout << "blocks_inv: "; for(unsigned j=0; j<blocks_inv.size(); j++){cout << "("<<blocks_inv[j].first<<","<<blocks_inv[j].second<<") ";}cout << endl;
#endif//DEBUG_CREXTDRL
			if( blocks.size() == 0 ){
				tdl_sort( quot, id, tdls );
				add_quot_rrrmt(tdls, quot, id, asce);
			}else{
				tdld = std::numeric_limits< int >::max();
				reverse_blocks(quot, id, blocks, tdlscen, revscen, used, maxszen, currz, tdld);
				#ifdef DEBUG_CREXTDRL
				cout << "positive : ("<<blocks.size()<<" rev + "<<tdld<<" tdrls)"<<endl;
				#endif//DEBUG_CREXTDRL
				if(tdld + (int)blocks.size() <= bestd){
					add_quot_rrrmt(tdlscen, revscen, quot, id, true, asce);
					bestd = tdld + blocks.size();
				}
				currz.clear(); used.clear(); tdlscen.clear(); revscen.clear();

				tdld = std::numeric_limits< int >::max();
				reverse_blocks_inv(quot, id, blocks_inv, tdlscen, revscen, used, maxszen, currz, tdld);
				#ifdef DEBUG_CREXTDRL
				cout << "inv positive : ("<<blocks_inv.size()<<" rev + "<<tdld<<" tdrls)"<<endl;
				#endif//DEBUG_CREXTDRL
				if(tdld + (int)blocks_inv.size() <= bestd){
					add_quot_rrrmt(tdlscen, revscen, quot, id, false, asce);
					bestd = tdld + blocks.size();
				}
				currz.clear(); used.clear(); tdlscen.clear(); revscen.clear();

			}
			blocks.clear();
			blocks_inv.clear();
		}

		if( psign == DEC ){
			if( bps ){
				bpsce = new bpscenario( quot, nid, maxszen, true );
					if( !bpsce->isempty() ){
						asce->insert( bpsce );
					}else{
						delete bpsce;
				}
			}

				// get the negative blocks of the quotient permutation
			get_blocks(quot, nid, blocks);
			get_blocks(nid, quot, blocks_inv);
#ifdef DEBUG_CREXTDRL
			cout << "TO NEGATIVE"<<endl;
			cout << "blocks    : "; for(unsigned j=0; j<blocks.size(); j++){cout << "("<<blocks[j].first<<","<<blocks[j].second<<") ";}cout << endl;
			cout << "blocks_inv: ";	for(unsigned j=0; j<blocks_inv.size(); j++){cout << "("<<blocks_inv[j].first<<","<<blocks_inv[j].second<<") ";}cout << endl;
#endif//DEBUG_CREXTDRL
			if( blocks.size() == 0 ){
				tdl_sort( quot, nid, tdls );
				add_quot_rrrmt(tdls, quot, nid, asce);
			}else{
				tdld = std::numeric_limits< int >::max();
				reverse_blocks(quot, nid, blocks, tdlscen, revscen, used, maxszen, currz, tdld);
				#ifdef DEBUG_CREXTDRL
				cout << "negative: ("<< blocks.size() << " rev + "<<tdld<<" tdrls )"<<endl;
				#endif//DEBUG_CREXTDRL
				if(tdld + (int)blocks.size() <= bestd){
					add_quot_rrrmt(tdlscen, revscen, quot, nid, true, asce);
					bestd = tdld + blocks.size();
				}
				currz.clear(); used.clear(); tdlscen.clear(); revscen.clear();

				tdld = std::numeric_limits< int >::max();
				reverse_blocks_inv(quot, nid, blocks_inv, tdlscen, revscen, used, maxszen, currz, tdld);
				#ifdef DEBUG_CREXTDRL
				cout << "inv negative: ("<< blocks_inv.size() << " rev + "<<tdld<<" tdrls )"<<endl;
				#endif//DEBUG_CREXTDRL
				if(tdld + (int)blocks_inv.size() <= bestd){
					add_quot_rrrmt(tdlscen, revscen, quot, nid, false, asce);
					bestd = tdld + blocks.size();
				}
				currz.clear(); used.clear(); tdlscen.clear(); revscen.clear();
			}
			blocks.clear();
			blocks_inv.clear();
		}

//		cout << *asce<<endl;

			// rename the elements of the rearrangements from quot. perm. to the real permutation
		asce->rename( quotgmap, g1 );
//		insert(asce);
			// simplify the scenario
		tmp = asce->simplify();
		if( tmp->isempty() ){	// impossible -> take as is
			delete tmp;
			insert( asce );
		}else{					// otherwise -> take the simplified version
			delete asce;
			insert( tmp );
		}

		quot.clear();
		id.clear();
		nid.clear();

			// if sorting to negative (resp. positive) identity was more parsimonious
			// store the reversal tdl szenario to the negative (resp. positive) id
//		if ( pbestd <= nbestd ){
//			for(unsigned j=0; j<pgenomes.size(); j++){
//				tdl_sort( pgenomes[j], id, tdls );
//				tdl_szenarios.push_back( tdls );
//				revszens.push_back( prevszen[j] );
//				tdls.clear();
//			}
//		}
//
//		if ( nbestd <= pbestd ){
//			if ( nbestd < pbestd ){
//				tdl_szenarios.clear();
//				revszens.clear();
//			}
//			for(unsigned j=0; j<ngenomes.size(); j++){
//				tdl_sort( ngenomes[j], nid, tdls );
//				tdl_szenarios.push_back( tdls );
//				revszens.push_back( nrevszen[j] );
//				tdls.clear();
//			}
//		}

			// sort (low numbers of real tdl first, i.e. high numbers of transpositions first)
//		for(unsigned j=0; j<tdl_szenarios.size(); j++){
//			ttc.push_back( make_pair(tdl_transposition_count(tdl_szenarios[j]), j) );
//		}
//		sort( ttc.begin(), ttc.end() );

			// store the sorted alternative scenarios
//		for(unsigned j=0; j<ttc.size(); j++){
//
//			for(unsigned k=0; k<tdl_szenarios[ttc[j].second].size(); k++){
//				tdrl *tmptdrl = new tdrl( abs((*it)->i.first)-1, abs((*it)->i.second)-1, g1 );
//				for(unsigned l=0; l<(*it)->children.size(); l++){
//					tmptdrl->add_elements(abs((*it)->children[l]->i.first)-1,
//						abs((*it)->children[l]->i.second)-1,
//						(tdl_szenarios[ttc[j].second][k][l])? 0 : 1 ,g1);
//				}
//				rearrangements.insert(tmptdrl);
//			}
//
//			for( unsigned k=0; k<revszens[ttc[j].second].size(); k++ ){
//				rearrangements.insert(new rev(abs(revszens[ttc[j].second][k].first)-1, abs(revszens[ttc[j].second][k].second)-1, g1 ));
//			}
//			// get_tdl_szenario( id, tdl_szenarios[ttc[j].second], revszens[ttc[j].second], (*it), g1, g2, nmap, szenarios );
//		}
//
//		r = (*it)->i;
//		r.first = abs(r.first);
//		r.second = abs(r.second);
//		q = translate_interval(g1, r, g2 );
//
//			// if crex could identify a unique szenario -> apply it
//			// if not then the scenario will be marked as incomplete
//			// because g1 is not transformed into g2 if the quotient
//			// permutation is not sorted
//		if( tdl_szenarios.size() == 1 ){
//			genom q;
//			vector<itnode *> newchildren;
//			q = genom(quot.size(), 0 );
//				// apply the reversals
//			for(unsigned i=0; i<revszens[0].size(); i++){
//				reverse(q, revszens[0][i]);
//			}
//			for(unsigned i=0; i<tdl_szenarios[0].size(); i++){
//				tdl( q, tdl_szenarios[0][i] );
//			}
//				// reverse the children which are to reverse
//				// and store the new child order
//			for( unsigned i=0; i<q.size(); i++){
//				if (q[i] < 0){
//					q[i] *= -1;
//					interval_tree_reverse((*it)->children[ q[i]-1 ]);
//				}
//				newchildren.push_back( (*it)->children[ q[i]-1 ] );
//			}
//			(*it)->children = newchildren;
//		}
//
//			// free memory for the next iteration
	}
}*/

string crex::typestr() const{
	string r = "crex(";
	if(_complete){
		r+="full)";
	}else{
		r+="part)";
	}
	return r;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//ostream &crex::output(ostream &out) const{
//	out << "crex(";
//	if( complete ) cout << "complete";
//	else cout << "incomplete";
//	cout <<") : [";
//	for(set<rrrmt*, DereferenceLess>::iterator it=sce.begin(); it != sce.end(); it++){
//		cout << *(*it) <<endl;
//	}
//	out<<"]";
//	return out;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper function computing the correct ranges for a tdrl which was computed for a quotient permutation
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void adjust_tdrl_ranges( vector<vector<vector<bool> > > &tdlscen, itnode *nd, int n){
//	for(unsigned i=0; i<tdlscen.size(); i++){
//		adjust_tdrl_ranges( tdlscen[i], nd, n);
//	}
//}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper function computing the correct ranges for a reversal which was computed for a quotient permutation
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//void adjust_tdrl_ranges( vector<bool> &tdl, itnode *nd, const genom &quot, int n){
//	vector<bool> tmp;
//
////	copy(tdl.begin(), tdl.end(), ostream_iterator<bool>(cout, ""));cout << endl;
//	tmp.insert( tmp.end(), nd->i.first-1, true );
//	for(unsigned i=0; i<tdl.size(); i++){
//			tmp.insert(tmp.end(),
//						abs(nd->children[abs(quot[i])-1]->i.second)-abs(nd->children[abs(quot[i])-1]->i.first)+1,
//						tdl[i]);
//	}
//	tmp.insert( tmp.end(), n-nd->i.second, false );
//	tdl = tmp;
////	copy(tdl.begin(), tdl.end(), ostream_iterator<bool>(cout, ""));cout << endl;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void adjust_rev_ranges( pair<int, int> &r, const genom &quot, itnode *nd){
////	cout << "q "<<quot<<endl;
////	cout << "rev ("<<r.first<<"-"<<r.second<<") -> ";
////	for(vector<itnode*>::iterator it = nd->children.begin(); it!=nd->children.end(); it++){
////		cout <<"("<< (*it)->i.first <<","<< (*it)->i.first<<"] ";
////	}cout << endl;
//
//	int s = nd->i.first - 1,
//		l = 0;
//	for( int i=0; i<r.first; i++ ){
//		s += abs(nd->children[ abs(quot[i])-1 ]->i.second) - abs(nd->children[ abs(quot[i])-1 ]->i.first)+1;
//	}
//
//	for(int i=r.first; i<=r.second; i++){
//		l += abs(nd->children[ abs(quot[i])-1 ]->i.second) - abs(nd->children[ abs(quot[i])-1 ]->i.first)+1;
//	}
//	r.first = s;
//	r.second = s + l - 1;
////	cout <<r.first<<"-"<<r.second<<endl;
//}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//void adjust_ranges( vector<vector<bool> > &tdlscen, itnode *nd, int n){
//	genom tmpg(nd->children.size(), 0);
////	cout << "adjusting to intervals ";
////	for(unsigned i=0; i<nd->children.size(); i++)
////		cout << "("<<nd->children[i]->i.first<<","<<nd->children[i]->i.second<<") ";
////	cout << endl;
//	for(unsigned i=0; i<tdlscen.size(); i++){
////		copy(tdlscen[i].begin(), tdlscen[i].end(), ostream_iterator<bool>(cout, "") );cout << endl;
////		cout << tmpg<<endl;
//		adjust_tdrl_ranges( tdlscen[i], nd, tmpg, n);
//		tdl( tmpg, tdlscen[i] );
////		copy(tdlscen[i].begin(), tdlscen[i].end(), ostream_iterator<bool>(cout, "") );cout << endl<<endl;;
//	}
////	cout << "done adjusting"<<endl;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void add_quot_rrrmt( vector<vector<bool> > &ts, genom q, const genom &tgt, alternative *alt ){
	ordered *ret = 0;	// the ordered scenario to return
	rrrmt *tmp;			// temp for the rrrmts

	ret = new ordered();
	for(unsigned i=0; i<ts.size(); i++){
		tmp = new tdrl(ts[i], q);
		tmp->apply(q);
		ret->push_back( tmp );
	}
	alt->insert( ret );

	if( q != tgt ){
		cerr << "add_quot_rrrmt(): target quotient permutation not reached"<<endl;
		exit(1);
	}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void add_quot_rrrmt( vector<vector<vector<bool> > > &ts, vector<vector<pair<int, int> > > &rs,
		const genom &q, const genom &tgt, bool revfirst, alternative *alt ){

	genom tmpq;
	ordered *tmpo = NULL;
	rrrmt *tmpr = NULL;

	if( rs.size() != ts.size() ){
		cerr << "internal error: add_quot_rrrmt unequal number of reversal tdrl scenarios"<<endl;
	}

//	cout << "--------------"<<endl;
	for(unsigned i=0; i<rs.size(); i++){
		tmpo = new ordered();
		tmpq = q;

		if(revfirst){
			for(unsigned j=0; j<rs[i].size(); j++){
//				cout << "tmp quot "<<tmpq<<endl;
//				cout << "+rev("<< rs[i][j].first<< ","<< rs[i][j].second <<")"<<endl;
				tmpr = new rev( rs[i][j].first, rs[i][j].second, tmpq );
//				cout << *tmpr<<endl;
				tmpr->apply(tmpq);
//				cout << "-> "<<tmpq<<endl;
				tmpo->push_back(tmpr);

			}
			for(unsigned j=0; j<ts[i].size(); j++){
//				cout << "tmp quot "<<tmpq<<endl;
//				cout << "+tdrl("; copy(ts[i][j].begin(), ts[i][j].end(), ostream_iterator<bool>(cout,""));cout << ")"<<endl;
				tmpr = new tdrl( ts[i][j], tmpq );
//				cout << *tmpr<<endl;
				tmpr->apply(tmpq);
//				cout << "-> "<<tmpq<<endl;
				tmpo->push_back(tmpr);
			}
		}else{
			for(unsigned j=0; j<ts[i].size(); j++){
//				cout << "tmp quot "<<tmpq<<endl;
//				cout << "+tdrl("; copy(ts[i][j].begin(), ts[i][j].end(), ostream_iterator<bool>(cout,""));cout << ")"<<endl;
				tmpr = new tdrl( ts[i][j], tmpq );
//				cout << *tmpr<<endl;
				tmpr->apply(tmpq);
//				cout << "-> "<<tmpq<<endl;
				tmpo->push_back(tmpr);
			}
			for(unsigned j=0; j<rs[i].size(); j++){
//				cout << "tmp quot "<<tmpq<<endl;
//				cout << "+rev("<< rs[i][j].first<< ","<< rs[i][j].second <<")"<<endl;
				tmpr = new rev( rs[i][j].first, rs[i][j].second, tmpq );
//				cout << *tmpr<<endl;
				tmpr->apply(tmpq);
//				cout << "-> "<<tmpq<<endl;
				tmpo->push_back(tmpr);
			}
		}
#ifdef DEBUG_CREXTDRL
		cout << "adding "<<*tmpo<<endl;
#endif//DEBUG_CREXTDRL
		alt->insert(tmpo);
		if( tmpq != tgt ){
			cerr << "add_quot_rrrmt(): target quotient permutation not reached"<<endl;
			cerr << q << " -> "<<tmpq<<" != "<<tgt<<endl;
			exit(1);
		}
	}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//ordered *adjust_ranges( vector<vector<bool> > &ts, itnode *nd, const genom &g, int n){
//
//	genom tmpq(nd->children.size(), 0),
//		tmpg = g, t;
//
//	ordered *ctmp = NULL;
//	vector<bool> tmps;
//
//	ctmp = new ordered();
//	for(unsigned i=0; i<ts.size(); i++){
//		tmps = ts[i];
////		copy(ts[i].begin(), ts[i].end(), ostream_iterator<bool>(cout, "") );cout << " -> ";
//		adjust_tdrl_ranges( ts[i], nd, tmpq, n);
////		copy(ts[i].begin(), ts[i].end(), ostream_iterator<bool>(cout, "") );cout << endl;
////		cout << " : "<<tmpq<<" . "<<tmpg<<endl;
//		ctmp->push_back( new tdrl( ts[i], tmpg ) );
//		tdl( tmpq, tmps );
//		tdl( tmpg, ts[i] );
//	}
//
//	return ctmp;
//
////	genom tmp = g;
////
////	for(unsigned i=0; i<ts.size(); i++){
////		ctmp->push_back( new tdrl( ts[i], tmp ) );
//////		tdl( tmp, ts[i] );
////	}
////	sce.insert(ctmp);
//}
//
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//void adjust_ranges( vector<vector<vector<bool> > > &ts,
//	vector<vector<pair<int, int> > > &rs, itnode *nd, bool revfirst, const genom &g, int n, alternative *alts){
//
//	genom tmpg,
//		tmpq;
//	ordered *ctmp = NULL;
//	pair<int,int> tmpr;
//	vector<bool> tmpt;
//
//	if( rs.size() != ts.size() ){
//		cerr << "add_tdlrev_scenario: unequal number of reversal tdrl scenarios"<<endl;
//	}
//
//	for(unsigned i=0; i<rs.size(); i++){
//		ctmp = new ordered();
//		tmpq = genom(nd->children.size(), 0);
//		tmpg = g;
//
//		if(revfirst){
//			for(unsigned j=0; j<rs[i].size(); j++){
//				tmpr = rs[i][j];
//				adjust_rev_ranges( rs[i][j], tmpq, nd );
//				ctmp->push_back( new rev(rs[i][j].first, rs[i][j].second, tmpg) );
//				reverse(tmpq, tmpr);
//				reverse(tmpg, rs[i][j]);
//			}
//			for(unsigned j=0; j<ts[i].size(); j++){
//				tmpt = ts[i][j];
//				adjust_tdrl_ranges( ts[i][j], nd, tmpq, n);
//				ctmp->push_back( new tdrl( ts[i][j], tmpg ) );
//				tdl( tmpq, tmpt );
//				tdl( tmpg, ts[i][j] );
//			}
//		}else{
//			for(unsigned j=0; j<ts[i].size(); j++){
//				tmpt = ts[i][j];
//				adjust_tdrl_ranges( ts[i][j], nd, tmpq, n);
//				ctmp->push_back( new tdrl( ts[i][j], tmpg ) );
//				tdl( tmpq, tmpt );
//				tdl( tmpg, ts[i][j] );
//			}
//			for(unsigned j=0; j<rs[i].size(); j++){
//				tmpr = rs[i][j];
//				adjust_rev_ranges( rs[i][j], tmpq, nd );
//				ctmp->push_back( new rev(rs[i][j].first, rs[i][j].second, tmpg) );
//				reverse(tmpq, tmpr);
//				reverse(tmpg, rs[i][j]);
//			}
//		}
//		alts->insert(ctmp);
//	}
//
////	if( revscen.size() != tdlscen.size() ){
////		cerr << "adjust_ranges: unequal number of reversal tdrl scenarios"<<endl;
////	}
////
////	for(unsigned i=0; i<revscen.size(); i++){
////		genom tmpg(nd->children.size(), 0);
////		if(revfirst){
////			for(unsigned j=0; j<revscen[i].size(); j++){
////				adjust_rev_ranges(revscen[i][j], tmpg, nd);
////				reverse(tmpg, revscen[i][j]);
////			}
////			for(unsigned j=0; j<tdlscen[i].size(); j++){
////				adjust_tdrl_ranges( tdlscen[i][j], nd, tmpg, n);
////				tdl( tmpg, tdlscen[i][j] );
////			}
////		}else{
////			for(unsigned j=0; j<tdlscen[i].size(); j++){
////				adjust_tdrl_ranges( tdlscen[i][j], nd, tmpg, n);
////				tdl( tmpg, tdlscen[i][j] );
////			}
////			for(unsigned j=0; j<revscen[i].size(); j++){
////				adjust_rev_ranges(revscen[i][j], tmpg, nd);
////				reverse(tmpg, revscen[i][j]);
////			}
////		}
////	}
//
//}







































/* @todo web crex funtions commented out -> move to separate file

















// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// wrapper for computing the number of common intervals from python
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int cnt_intervals( vector<int> p1, vector<int> p2, int circular ){
	genom  g1(p1, circular),
		g2(p2, circular);
	vector<pair<int, int> > comint;

	common_intervals(g1, g2, circular, 0, 0, comint);

	return comint.size();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// wrapper for computing the number of breakpoints from python
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int cnt_breakpoints( vector<int> p1, vector<int> p2, int circular ){
	genom  g1(p1, circular),
		g2(p2, circular);
	hdata hd;

	init_data(hd, p1.size(), circular);

	return breakpoints(g1, g2, hd);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// wrapper for computing the number of breakpoints from python
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int cnt_reversals( vector<int> p1, vector<int> p2, int circular ){
	genom  g1(p1, circular),
		g2(p2, circular);
	dstnc_inv rd(p1.size(), circular );
	return rd.calc( g1, g2 );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex ..
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool compare( genom g1, genom g2, vector<string> namemap, int circular, vector<string> &szenarios, int maxprev ){
	genom frnt;
	int n = g1.size();	// length of the genomes
	itnode *iroot;		// interval tree relative to g1
	list<itnode *> pnode,  	// prime nodes of the two trees
		node;
	stringstream out;

	szenarios.clear();

	interval_tree(g1, g2, n, &iroot);
	interval_tree_primenodes(iroot, pnode);
	interval_tree_nodes(iroot, node);

	find_transpositions(node, g1, g2, namemap, szenarios);
	find_reversetranspositions(node, g1, g2, namemap, szenarios);
	find_reversals(node, g1, g2, namemap, szenarios);
	find_tdls(node, g1, g2, namemap, szenarios, maxprev);

	front( iroot, g1, frnt );

	interval_tree_free(iroot);

//	cout << g1 << endl;
//	cout << g2 << endl;
//	cout << frnt << endl;
	if( g2 == frnt ){
		return true;
	}else{
		return false;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex wrapper for python
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<string> compare( vector<int> p1, vector<int> p2, vector<string> namemap, int circular, int norm_to, int maxprev){
	genom  g1(p1, circular, &namemap, norm_to),// the two genomes
		g2(p2, circular, &namemap, norm_to);

	vector<string> szenario;
	compare( g1, g2, namemap, circular, szenario, maxprev );
	return szenario;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// check if two gene orders (of same length) are equal (wrapper for python)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool equal(vector<int> p1, vector<int> p2, int circular){
	genom  g1(p1, circular ),// the two genomes
		g2(p2, circular );

	return (g1 == g2);
}




// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex function: reversal localization
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void find_reversals( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios ){

	pair<int,int> r, q;
	stringstream treestream;

	for(list<itnode *>::iterator it = node.begin(); it != node.end(); ){
			// do not consider prime nodes or nodes with prime parent
		if( (*it)->type == PRI || ((*it)->parent != NULL && (*it)->parent->type == PRI)){
			it++;
			continue;
		}

			// negative root node or nodes with sign different from the linear parent have to be reversed
		if( ( (*it)->parent == NULL && (*it)->sign == DEC ) || ( (*it)->parent != NULL &&  (*it)->parent->sign != (*it)->sign ) ){
			r = (*it)->i;
			r.first = abs(r.first);
			r.second = abs(r.second);
			q = translate_interval(g1, r, g2 );
			szenarios.push_back( "#next"+int2string((*it)->i.first)+"_"+int2string((*it)->i.second) );
			szenarios.push_back( "rev" );

			treestream << "|rev1|";
			interval_tree_print( (*it) , g1, treestream);
			szenarios.push_back( treestream.str() );
			treestream.str("");

			szenarios.push_back("r|rev1|"+int2string(r.first)+","+int2string(r.second) );
			szenarios.push_back("q|rev1|"+int2string(q.first)+","+int2string(q.second) );

			interval_tree_reverse( (*it) );

			treestream << "|rev1|";
			interval_tree_print( (*it), g1, treestream);
			szenarios.push_back( treestream.str() );
			treestream.str("");

			it = node.erase(it);
		}else{
			it++;
		}
	}
}




// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex function: reverse transposition localization
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void find_reversetranspositions( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios){

	genom f;
	unsigned nesign = 0,	// number of sign differences between the current node and its children
		nepos = 0;	// last position where a sign difference occured
	int ptype,	// type and ...
		psign;  // ... sign of the parent
	stringstream click,
		mout,
		out2,
		treestream;
	pair<int, int> r, q;
	vector< vector<itnode *> > newchildren(2);
	vector<string> par_tree;
	string op_str,
		op_com;

	for(list<itnode *>::iterator it = node.begin(); it != node.end(); ){

			// only look for linear nodes
		if( (*it)->type == PRI){
			it++;
			continue;
		}
			// get the parent sign and type
		if( (*it)->parent == NULL || ((*it)->parent->type == PRI && (*it)->parent->sign == UNK) ){
			ptype = LIN;
			psign = INC;
		}else{
			ptype = (*it)->parent->type;
			psign = (*it)->parent->sign;
		}
		cout <<"rT? "<< (*it)->i.first<<","<< (*it)->i.second<<endl;
			 // nodes have to have a sign different from the parent
		if( ptype == UNK || (*it)->sign==INC ||  psign == (*it)->sign ){
			cout << "no rT1"<<endl;
			it++;
			continue;
		}

		for(unsigned j=0; j<(*it)->children.size(); j++){

			if( (*it)->children[j]->type == PRI ){
				it++;
				continue;
			}

			if( (*it)->sign != (*it)->children[j]->sign ){
				cout << j<< " " << (*it)->sign <<"!="<< (*it)->children[j]->sign<< " nesign "<<nesign<<endl;
				nesign++;
				nepos = j;
			}
		}

		if( !(nesign == 1 && (nepos == 0 || nepos == (*it)->children.size()-1 ))){
			it++;
			cout  << "abort nesign" << nesign<<"  nepos "<<nepos << " size "<<(*it)->children.size()-1<<endl;
			continue;
		}

		cout << "parent "<<ptype<< " "<< psign<<endl;
		cout << "nesign "<< nesign <<" nepos "<< nepos<< endl;

//		out<<signchanges<<"</span></li>"<<endl;
		szenarios.push_back( "#next"+int2string((*it)->i.first)+"_"+int2string((*it)->i.second) );
		szenarios.push_back("revtra");

		interval_tree_print_prefix(*it, treestream);
		for(unsigned j=0; j<(*it)->children.size(); j++){
			r = (*it)->children[j]->i;
			r.first = abs(r.first);
			r.second = abs(r.second);
			q = translate_interval(g1, r, g2 );
			if( nepos == 0 ){	//(*it)->sign
				if( j>0 ){
					op_str = "|rev1|";
					op_com = ",";
					newchildren[0].push_back( (*it)->children[j] );
				}else{
					op_str = "|tra2|";
					op_com = ";";
					newchildren[1].push_back( (*it)->children[j] );
				}
			}else{
				if( j<(*it)->children.size()-1 ){
					op_str = "|rev1|";
					op_com = ",";
					newchildren[1].push_back( (*it)->children[j] );
				}else{
					op_str = "|tra2|";
					op_com = ";";
					newchildren[0].push_back( (*it)->children[j] );
				}
			}

			interval_tree_print((*it)->children[j], g1, treestream, op_str);if(j<(*it)->children.size()-1){treestream<<",";}
			szenarios.push_back("r"+op_str+int2string(r.first)+op_com+int2string(r.second));
			szenarios.push_back("q"+op_str+int2string(q.first)+op_com+int2string(q.second));
		}

		interval_tree_print_suffix(*it, treestream);
		szenarios.push_back( treestream.str() );
		treestream.str("");

			// reverse the reverse-transposition nodes
		for (unsigned j=0; j<newchildren.size(); j++){
			if( (nepos == 0 && j==0) || (nepos != 0 && j==newchildren.size()-1)){
				for (unsigned k=0; k< newchildren[j].size(); k++){
					interval_tree_reverse( newchildren[j][k] );
				}
				reverse( newchildren[j].begin(), newchildren[j].end() );
				op_str = "|rev1|";
			}else{
				op_str = "|tra2|";
			}
			for (unsigned k=0; k< newchildren[j].size(); k++){
				interval_tree_print( newchildren[j][k] , g1, treestream, op_str);if(k<newchildren[j].size()-1){treestream<<",";}
			}
			if(j<newchildren.size()-1)
				treestream<<",";
			par_tree.push_back( treestream.str() );
			treestream.str("");
		}
//		swap( newchildren[0], newchildren[1] );
//		swap( par_tree[0], par_tree[1] );

			// transpose the reverse transposition nodes
		if( (*it)->sign == INC ){
			(*it)->sign = DEC;
		}else if( (*it)->sign == DEC ){
			(*it)->sign = INC;
		}

		(*it)->children.clear();
		interval_tree_print_prefix(*it, treestream);
		for(unsigned j=0; j<par_tree.size(); j++){
			treestream << par_tree[j];
		}
		interval_tree_print_suffix(*it, treestream);
		szenarios.push_back( treestream.str() );
		treestream.str("");
		par_tree.clear();

		(*it)->children.insert((*it)->children.end(), newchildren[0].begin(), newchildren[0].end());
		(*it)->children.insert((*it)->children.end(), newchildren[1].begin(), newchildren[1].end());

		newchildren[0].clear();
		newchildren[1].clear();
		nesign = 0;
		nepos = 0;
		it = node.erase(it);

	}
	return;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex function: tdl localization
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void find_tdls( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios, int maxprev){

	genom g,		// identified genome
		quot,		// quotient permutation for a prime node
		id,			// identity and ..
		nid;		// .. negative identity
	int n = g1.size(),	// length of the genomes
		psign,			// the sign of the parent
		pbestd = std::numeric_limits< int >::max(),			// best tdl distance if blocks reversed to positive
		nbestd = std::numeric_limits< int >::max(); 		// ... or negative
	map<genom, bool> used;		// for avoiding recomputations in the reversal process (temporary memory for reversal function)
	pair<int, int> r, q;
	vector<genom> pgenomes,		// best results of the reversals of the negative blocks -> positive genomes
		ngenomes;				// best results of the reversals of the positive blocks -> negative genomes
	vector< pair<int,int> > neg_blocks,	// negative blocks of the qoutient permutation of a prime node
		pos_blocks,						// positive blocks of the qoutient permutation of a prime node
		currz;							// temporary memory for the reversal funtion
	vector< pair< pair<int,int>, int > > ttc;	// 1st: pair of (#tdls, #transp), 2nd: index of tdl-&reversal-szenario
	vector<vector<pair<int, int> > > prevszen,	// reversal szenarios leading to the positive / negative
		nrevszen,				// permutations
		revszens;				// all (combination of p&n) reversal szenarios for the tdl szenarios
	vector<vector<bool> > tdls;			// one tdl szenario
	vector<vector<vector<bool> > > tdl_szenarios; 	// the tdl szenarios (each tdl is a bool vector, a szenario contains
							// possibly several tdls and we store all szenarios ->3x vector)

		// identify the permutations
	g = g2.identify_g(g1);
//	out <<"g: "<< g << "<br>";
		// search the nodes
	for(list<itnode*>::iterator it = node.begin(); it!= node.end(); it++){
			// only look at pnodes
		if((*it)->type != PRI){
			continue;
		}
			// initialise the
			// - tdl distance to std::numeric_limits< int >::max()
			// - the quotient permutation
			// - the identity and negative identity (of len(quot) length)
			// - the positive and negative blocks
			// - the sign of the parent
		pbestd = std::numeric_limits< int >::max();
		nbestd = std::numeric_limits< int >::max();
		quot = genom((*it)->children.size(), 0);
		for(unsigned k=0; k<(*it)->children.size(); k++){
			quot[k] = g[ abs((*it)->children[k]->i.first)-1 ];
		}

		quotient_permutation(quot, n);
		for(unsigned k=0; k<(*it)->children.size(); k++){
			if((*it)->children[k]->sign == DEC)
				quot[k] *= -1;
		}

		id = genom(quot.size(), 0 );
		nid = genom_nid(quot.size());

		if( (*it)->parent != NULL ){
			psign = (*it)->parent->sign;
		}else{
			psign = INC;
		}

			// get the positive and negative blocks of the quotient permutation
		get_blocks(quot, id, neg_blocks);

//		cout << "blocks "<< neg_blocks.size()<<endl;

		if( neg_blocks.size() == 0 ){
			tdl_sort( quot, id, tdls );
			tdl_szenarios.push_back( tdls );
			tdls.clear();
			revszens.push_back( vector<pair<int, int> >() );
		}else{
//			cout << "maxprev "<<maxprev<<endl;
			reverse_blocks(quot, id,    neg_blocks, tdl_szenarios,  revszens, used, maxprev, currz, pbestd);
		}

//		cout << tdl_szenarios.size()<< "tdl_szenarios"<<endl;
//		cout << prevszen.size()<<"prevszen"<<endl;

//		get_blocks(quot, nid, neg_blocks);

			// try all possibilities to reverse the blocks (with respect to the sign of the parent node)
			// and select outcomes with minimal tdl distance (stored in pbestd)
			// if the sign is unknown -> compute both
//		if( psign == INC || psign == UNK ){
//			if( maxprev == 0 || maxprev >= (int) neg_blocks.size() ){
//				reverse_blocks(quot, id, 0, neg_blocks, pgenomes, prevszen, used, currz, pbestd);
//				currz.clear();
//				used.clear();
//			}
//		}
//		if( psign == DEC || psign == UNK ){
//			if( maxprev == 0 || maxprev >= (int) neg_blocks.size() ){
//				reverse_blocks(quot, nid, 1, pos_blocks, ngenomes, nrevszen, used, currz, nbestd);
//				currz.clear();
//				used.clear();
//			}
//		}

			// if sorting to negative (resp. positive) identity was more parsimonious
			// store the reversal tdl szenario to the negative (resp. positive) id
//		if ( pbestd <= nbestd ){
//			cout << "POS better"<<endl;
//			for(unsigned j=0; j<pgenomes.size(); j++){
//				tdl_sort( pgenomes[j], id, tdls );
//				tdl_szenarios.push_back( tdls );
//				revszens.push_back( prevszen[j] );
//				tdls.clear();
//			}
//		}
//
//		if ( nbestd <= pbestd ){
//			cout << "NEG better"<<endl;
//			if ( nbestd < pbestd ){
//				tdl_szenarios.clear();
//				revszens.clear();
//			}
//			for(unsigned j=0; j<ngenomes.size(); j++){
//				tdl_sort( ngenomes[j], nid, tdls );
//				tdl_szenarios.push_back( tdls );
//				revszens.push_back( nrevszen[j] );
//				tdls.clear();
//			}
//		}

			// sort (low numbers of real tdl first, i.e. high numbers of transpositions first)
		for(unsigned j=0; j<tdl_szenarios.size(); j++){
			ttc.push_back( make_pair(tdl_transposition_count(tdl_szenarios[j]), j) );
		}

//		cerr << "ttc "<<ttc.size() <<endl;
		sort( ttc.begin(), ttc.end() );

		r = (*it)->i;
		r.first = abs(r.first);
		r.second = abs(r.second);
		q = translate_interval(g1, r, g2 );
		for(unsigned j=0; j<ttc.size(); j++){
			get_tdl_szenario( id, tdl_szenarios[ttc[j].second], revszens[ttc[j].second], (*it), g1, g2, nmap, szenarios );
		}
//		cerr << szenarios.size()<<"szenarios" <<endl;
//		for(unsigned j=0; j<szenarios.size(); j++)
//			cerr <<"szen "<<j<<" " <<szenarios[j]<<endl;

			// if crex could identify a unique szenario -> apply it
		if( tdl_szenarios.size() == 1 ){
			genom q;
			vector<itnode *> newchildren;
			q = genom(quot.size(), 0 );
				// apply the reversals
			for(unsigned i=0; i<revszens[0].size(); i++){
				reverse(q, revszens[0][i]);
			}
			for(unsigned i=0; i<tdl_szenarios[0].size(); i++){
				tdl( q, tdl_szenarios[0][i] );
			}
				// reverse the children which are to reverse
				// and store the new child order
			for( unsigned i=0; i<q.size(); i++){
				if (q[i] < 0){
					q[i] *= -1;
					interval_tree_reverse((*it)->children[ q[i]-1 ]);
				}
				newchildren.push_back( (*it)->children[ q[i]-1 ] );
			}
			(*it)->children = newchildren;
		}

			// free memory for the next iteration
		quot.clear();
		id.clear();
		nid.clear();
		pgenomes.clear();
		ngenomes.clear();
		pos_blocks.clear();
		neg_blocks.clear();
		tdl_szenarios.clear();
		prevszen.clear();
		nrevszen.clear();
		revszens.clear();
		ttc.clear();

	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex function: transposition localization
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void find_transpositions( list<itnode *> node, const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios ){

	genom f;
//	int ptype, // type of the parent node
//		psign;
	pair<int,int> q1, q2, r1, r2;
	stringstream treestream;

	for(list<itnode *>::iterator it = node.begin(); it != node.end(); ){

			// only nodes with two children can be transpositions ..
			//
		if( (*it)->children.size() == 2 && (*it)->type == LIN &&
			(( (*it)->sign == DEC && (*it)->children[0]->sign == INC && (*it)->children[1]->sign == INC) ||
			( (*it)->sign == INC && (*it)->children[0]->sign == DEC && (*it)->children[1]->sign == DEC)) ){
			szenarios.push_back( "#next"+int2string((*it)->i.first)+"_"+int2string((*it)->i.second) );
			szenarios.push_back("tra");

			interval_tree_print_prefix(*it, treestream);
				interval_tree_print( (*it)->children[0] , g1, treestream, "|tra2|");
				treestream<<",";
				interval_tree_print( (*it)->children[1] , g1, treestream, "|tra1|");
			interval_tree_print_suffix(*it, treestream);
			szenarios.push_back( treestream.str() );
			treestream.str("");

			r1 = (*it)->children[0]->i;
			r1.first = abs(r1.first);
			r1.second = abs(r1.second);
			r2 = (*it)->children[1]->i;
			r2.first = abs(r2.first);
			r2.second = abs(r2.second);
			q1 = translate_interval(g1, r1, g2 );
			q2 = translate_interval(g1, r2, g2 );

			szenarios.push_back( "r|tra2|"+int2string(r1.first)+";"+int2string(r1.second) );
			szenarios.push_back( "r|tra1|"+int2string(r2.first)+","+int2string(r2.second) );
			szenarios.push_back( "q|tra2|"+int2string(q1.first)+";"+int2string(q1.second) );
			szenarios.push_back( "q|tra1|"+int2string(q2.first)+","+int2string(q2.second) );

			if ( (*it)->sign == INC )
				(*it)->sign = DEC;
			else
				(*it)->sign = INC;
			swap((*it)->children[0], (*it)->children[1]);

			interval_tree_print_prefix(*it, treestream);
				interval_tree_print( (*it)->children[0] , g1,  treestream, "|tra1|");
				treestream<<",";
				interval_tree_print( (*it)->children[1] , g1,  treestream, "|tra2|");
			interval_tree_print_suffix(*it, treestream);
			szenarios.push_back( treestream.str() );
			treestream.str("");

			it = node.erase(it);

		}else{
			it++;
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// crex function: determine tdl szenario
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void get_tdl_szenario( genom g, const vector<vector<bool> > &tdls, const vector<pair<int, int> > &revs, itnode *node,
	const genom &g1, const genom &g2, const vector<string> &nmap, vector<string> &szenarios){

	genom ng;		// the next genom
	stringstream treestream;
	int keepcnt = 0;
	pair<int, int> r, q;
	string op;

	szenarios.push_back( "#next"+int2string(node->i.first)+"_"+int2string(node->i.second) );

	if(revs.size() > 0 || tdls.size() > 1){
		r = node->i;
		r.first = abs(r.first);
		r.second = abs(r.second);
		q = translate_interval(g1, r, g2 );
		szenarios.push_back("r|pnode1|"+int2string(r.first)+","+int2string(r.second) );
		szenarios.push_back("q|pnode1|"+int2string(q.first)+","+int2string(q.second) );
	}

	for(unsigned i=0; i<revs.size(); i++){
		ng = g;
		reverse(ng, revs[i]);

		szenarios.push_back( "rev" );

		interval_tree_print_prefix(node, treestream);
		for(int j=0; j<(int)g.size(); j++){
			if (revs[i].first <= j && j <= revs[i].second )
				treestream << "|rev1|";
			if( g[j] < 0 ){
				interval_tree_reverse(node->children[ abs(g[j])-1 ]);
			}
			interval_tree_print( node->children[ abs(g[j])-1 ] , g1, treestream );if(j<(int)g.size()-1){treestream<<",";}
			if( g[j] < 0 ){
				interval_tree_reverse(node->children[ abs(g[j])-1 ]);
			}
		}
		interval_tree_print_suffix(node, treestream);
		szenarios.push_back( treestream.str() );
		treestream.str("");

		interval_tree_print_prefix(node, treestream);
		for(int j=0; j<(int)ng.size(); j++){
			if (revs[i].first <= j && j <= revs[i].second )
				treestream << "|rev1|";
			if( ng[j] < 0 ){
				interval_tree_reverse(node->children[ abs(ng[j])-1 ]);
			}
			interval_tree_print( node->children[ abs(ng[j])-1 ] , g1, treestream );if(j<(int)g.size()-1){treestream<<",";}
			if( ng[j] < 0 ){
				interval_tree_reverse(node->children[ abs(ng[j])-1 ]);
			}

		}
		interval_tree_print_suffix(node, treestream);
		szenarios.push_back( treestream.str() );
		treestream.str("");

		g=ng;
	}

	for(unsigned i=0; i<tdls.size(); i++){
		ng = g;
		tdl( ng, tdls[i] );
		keepcnt = 0;

		if( is_tdl(tdls[i]) ){
			szenarios.push_back( "tdl" );
		}else{
			szenarios.push_back( "tra" );
		}

		if(tdls.size() == 1 && revs.size() == 0){
			for(unsigned j=0; j<tdls[i].size(); j++){
				r = node->children[ j ]->i;
				r.first = abs(r.first);
				r.second = abs(r.second);
				q = translate_interval(g1, r, g2 );

				if(tdls[i][ j ] == false){
					szenarios.push_back("r|tra2|"+int2string(r.first)+";"+int2string(r.second) );
					szenarios.push_back("q|tra2|"+int2string(q.first)+";"+int2string(q.second) );
				}else{
					szenarios.push_back("r|tra1|"+int2string(r.first)+","+int2string(r.second) );
					szenarios.push_back("q|tra1|"+int2string(q.first)+","+int2string(q.second) );
				}
			}
		}

		interval_tree_print_prefix(node, treestream);
		for(unsigned j=0; j<tdls[i].size(); j++){
			if(tdls[i][j]){
				treestream<< "|tra1|"; // del
				keepcnt++;
			}else{
				treestream<< "|tra2|"; // kep
			}
			if( g[j] < 0 ){
				interval_tree_reverse(node->children[ abs(g[j])-1 ]);
			}
			interval_tree_print( node->children[ abs(g[j])-1 ] , g1, treestream );if(j<tdls[i].size()-1){treestream<<",";}
			if( g[j] < 0 ){
				interval_tree_reverse(node->children[ abs(g[j])-1 ]);
			}
		}
		interval_tree_print_suffix(node, treestream);
		szenarios.push_back( treestream.str() );
		treestream.str("");

		interval_tree_print_prefix(node, treestream);
		for(unsigned j=0; j<tdls[i].size(); j++){
			keepcnt--;
			if (keepcnt>=0){
				treestream<< "|tra1|"; // kep
			}else{
				treestream<< "|tra2|"; // del
			}

			if( ng[j] < 0 ){
				interval_tree_reverse(node->children[ abs(ng[j])-1 ]);
			}
			interval_tree_print( node->children[ abs(ng[j])-1 ] , g1, treestream );if(j<tdls[i].size()-1){treestream<<",";}
			if( ng[j] < 0 ){
				interval_tree_reverse(node->children[ abs(ng[j])-1 ]);
			}
		}
		interval_tree_print_suffix(node, treestream);
		szenarios.push_back( treestream.str() );
		treestream.str("");

		g = ng;
	}
}

//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// crex function:
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////void print_tdl_szenario_html( genom g, const vector<vector<bool> > &tdls, const vector<pair<int, int> > &revs, itnode *node,
////	const genom &g1, const vector<string> &nmap, stringstream &out){
////
////	genom ng;		// the next genom
////	string cls="rev";
////	int keepcnt = 0;
////
////	if (revs.size() == 1){
////		out <<"1 Reversal <br>";
////	}else if (revs.size() > 1){
////		out << revs.size()<<" Reversals<br>";
////	}
////	for(unsigned i=0; i<revs.size(); i++){
////		ng = g;
////		reverse(ng, revs[i]);
////		out << "<table border=0><tr><td>Rev.</td><td><table border=1 rules=\"none\" class=\"pri\"><tr>";
////		for(unsigned j=0; j<g.size(); j++){
////			out << "<td><div class=\""<<cls<<"\">"; interval_tree_print_table( node->children[ abs(g[j])-1 ] , g1, nmap, "", out, g[j]<0);	out << "</div></td>";
////		}
////		out << "</tr></table>"<<endl;
////		out << "</td></tr><tr><td> &rarr; </td><td>";
////		out << "<table border=1 rules=\"none\" class=\"pri\"><tr>"<<endl;
////		for(unsigned j=0; j<g.size(); j++){
////			out << "<td><div class=\""<<cls<<"\">"; interval_tree_print_table( node->children[ abs(ng[j])-1 ] , g1, nmap, "", out, ng[j]<0);	out << "</div></td>";
////		}
////		out << "</tr></table>"<<endl;
////		out << "</td></tr></table>";
////		g=ng;
////	}
////
////	for(unsigned i=0; i<tdls.size(); i++){
////		ng = g;
////		tdl( ng, tdls[i] );
////		keepcnt = 0;
////		out << "<table border=0><tr><td>";
////		if( is_tdl(tdls[i]) ){
////			out << "TDL";
////		}else{
////			out << "Transp.";
////		}
////
////		out << "</td><td><table border=1 rules=\"none\" class=\"pri\"><tr>"<<endl;
////		for(unsigned j=0; j<tdls[i].size(); j++){
////			if(tdls[i][j]){
////				cls = "tdlkep";
////				keepcnt++;
////			}else{
////				cls = "tdldel";
////			}
////			out << "<td><div class=\""<<cls<<"\">"; interval_tree_print_table( node->children[ abs(g[j])-1 ] , g1, nmap, "", out, g[j]<0);	out << "</div></td>";
////		}
////		out << "</tr></table>"<<endl;
////		out << "</td></tr><tr><td> &rarr; </td><td>";
////		out << "<table border=1 rules=\"none\" class=\"pri\"><tr>"<<endl;
////		for(unsigned j=0; j<tdls[i].size(); j++){
////			keepcnt--;
////			if (keepcnt>=0){
////				cls = "tdlkep";
////			}else{
////				cls = "tdldel";
////			}
////			out << "<td><div class=\""<<cls<<"\">"; interval_tree_print_table( node->children[ abs(ng[j])-1 ] , g1, nmap, "", out, ng[j]<0);	out << "</div></td>";
////		}
////		out << "</tr></table>"<<endl;
////		out << "</td></tr></table>";
////		g = ng;
////	}
////
//////	for(unsigned j=0; j<szenario.size() -1; j++){
//////		keepcnt = 0;
////////		for(unsigned k=0; k<tdls[j].size(); k++)
////////			out << tdls[j][k]<<" ";
////////		out << endl;
//////		out << "<table border=0><tr><td>";
//////		if( is_tdl(tdls[j]) ){
//////			out << "TDL";
//////		}else{
//////			out << "Transp.";
//////		}
//////
//////		out << "</td><td><table border=1 rules=\"none\" class=\"pri\"><tr>"<<endl;
//////		for(unsigned k=0; k<szenario[j].size(); k++){
//////			if(tdls[j][k]){
//////				cls = "tdlkep";
//////				keepcnt++;
//////			}else{
//////				cls = "tdldel";
//////			}
//////			out << "<td><div class=\""<<cls<<"\">"; interval_tree_print_table( node->children[ abs(szenario[j][k])-1 ] , g1, nmap, "", out);	out << "</div></td>";
//////		}
//////		out << "</tr></table>"<<endl;
//////		out << "</td></tr><tr><td> &rarr; </td><td>";
//////		out << "<table border=1 rules=\"none\" class=\"pri\"><tr>"<<endl;
//////		for(unsigned k=0; k<szenario[j+1].size(); k++){
//////			keepcnt--;
//////			if (keepcnt>=0){
//////				cls = "tdlkep";
//////			}else{
//////				cls = "tdldel";
//////			}
//////			out << "<td><div class=\""<<cls<<"\">"; interval_tree_print_table( node->children[ abs(szenario[j+1][k])-1 ] , g1, nmap, "", out);	out << "</div></td>";
//////		}
//////		out << "</tr></table>"<<endl;
//////		out << "</td></tr></table>";
//////	}
////}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pair<int, int> translate_interval(const genom &g1, const pair<int,int> &i1, const genom &g2 ){
	pair<int,int> i2 = make_pair( std::numeric_limits< int >::max(), 0 );
	vector<int> inv_g2(g1.size()+1, 0);

	for(unsigned i=0; i<g2.size(); i++){
		inv_g2[ abs(g2[i]) ] = i+1;
	}

	for(int i=abs(i1.first)-1; i<abs(i1.second); i++){
//		cout << "looking for "<<abs(g1[i]) << " -> pos "<< inv_g2[abs(g1[i])];
		if( inv_g2[abs(g1[i])] > i2.second )
			i2.second = inv_g2[abs(g1[i])];
		if( inv_g2[abs(g1[i])] < i2.first )
			i2.first = inv_g2[abs(g1[i])];
//		cout << " -> ["<<i2.first << ","<<i2.second<<"]"<<endl;
	}

	return i2;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string tree(vector<int> p1, vector<int> p2, vector<string> namemap, int circular, int norm_to){
	genom  g1(p1, circular, &namemap, norm_to),// the two genomes
		g2(p2, circular, &namemap, norm_to);
	int n = p1.size();		// length of the genomes
	stringstream out;
	itnode *iroot;			 	// interval tree relative to g1

		// compute the interval tree and get the prime nodes
	interval_tree(g1, g2, n, &iroot);

	interval_tree_print(iroot, g1, out);
	interval_tree_free(iroot);
	return out.str();
}
//
////// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////
////string tree_dot( vector<int> p1, vector<int> p2, vector<string> namemap, int circular, int norm_to ){
////
////	genom  g1(p1, circular, norm_to),// the two genomes
////		g2(p2, circular, norm_to);
////	int n = p1.size();		// length of the genomes
////	stringstream out;
////	itnode *iroot;			 	// interval tree relative to g1
////
////		// compute the interval tree and get the prime nodes
////	interval_tree(g1, g2, n, &iroot);
////
////	interval_tree_print_dot(iroot, g1, out, namemap);
////	interval_tree_free(iroot);
////	return out.str();
////}
////
////// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////
////string tree_table( vector<int> p1, vector<int> p2, vector<string> namemap, int circular, string prefix, int norm_to ){
////
////	genom  g1(p1, circular, norm_to),	// the two genomes
////		g2(p2, circular, norm_to);
////	int n = p1.size();					// length of the genomes
////	stringstream out;
////	itnode *iroot;			 			// interval tree relative to g1
////
////		// compute the interval tree print it and delete it
////	interval_tree(g1, g2, n, &iroot);
////	interval_tree_print_table(iroot, g1, namemap, "r"+prefix, out);
////	interval_tree_free(iroot);
////
////
////		// the same in the other direction
////	interval_tree(g2, g1, n, &iroot);
////	interval_tree_print_table(iroot, g2, namemap, "q"+prefix, out);
////	interval_tree_free(iroot);
////
////	return out.str();
////
////}

*/
