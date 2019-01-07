/**
 * @file tdl.cpp tdl functions
 * @author: Matthias Bernt, Kai Ramsch
 */
#include <algorithm>
#include <iostream>
#include <iterator>
#include <math.h>
#include <numeric>

#include "helpers.hpp"
#include "tdl.hpp"
#include "genom.hpp"

#define THETA 0
#define PHI 1


//#define DEBUG_TDLBFMEDIAN

using namespace std;



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// get the minimum of d(p1,p2), d(p2,p1)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int min_tdrl_distance(const genom &origin,const genom &goal, int &dir){

	int d1,
		d2;

	d1 = tdrl_distance(origin, goal);
	d2 = tdrl_distance(goal, origin);

	if( d1 < d2 ){
		dir = 1;
		return d1;
	}else if(d2 < d1){
		dir = -1;
		return d2;
	}else{
		dir = 0;
		return d1;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// print chains
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void print_chains( const vector<vector<int> > &chains ){
	for(unsigned i=0; i<chains.size(); i++){
		cout << "chain "<<i<<" : ";
		copy(chains[i].begin(), chains[i].end(), ostream_iterator<int>(cout," "));
		cout << endl;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper/degubbing function: print dp matrix
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void print_dp( const vector<vector<vector<int> > > &a ){
	for( unsigned i=0; i<a.size(); i++ ){
		for( unsigned j=0; j<a[i].size(); j++ ){
			cout << "("<<a[i][j][0]<<","<<a[i][j][1]<<")\t";
		}
		cout << endl;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper/degubbing function: print tdrl matrix
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void print_tdrl( const vector<int> &tdrl, vector<string> *nmap ){
	for( unsigned i=1; i<tdrl.size(); i++ ){
		if(nmap == NULL){
			cout << i;
		}else{
			cout << (*nmap)[ i ];
		}
		cout << ":"<<tdrl[i]<<" ";
	}
	cout << endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper/degubbing function: print tdrl matrix
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void print_tdrls( const vector<vector<int> > &tdrls ){
	for( unsigned i=0; i<tdrls.size(); i++ ){
		print_tdrl(tdrls[i]);
	}
}

void print_tdrl( const vector<int> &tdrl, const genom &g, vector<string> *nmap){
	for( unsigned i=0; i<g.size(); i++ ){
		print_element(g[i], cout, 1, "", nmap);
		if( tdrl[abs(g[i])] == 0 ) cout << "*";
		cout << " ";
	}
	cout << endl;
}

void print_tdrl_scenario( const vector<vector<int> > &scenario, const genom &g, vector<string> *nmap){
	genom tmpg = g;
	for( unsigned i=0; i<scenario.size(); i++ ){
		print_tdrl( scenario[i], tmpg, nmap );
		tdrl_apply(tmpg, scenario[i]);
		cout << tmpg << endl<<endl;
	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a tdrl; the tdrl vector specifies for each position in g if it is
// kept in the first copy (0), or in the second copy (1)
// tdrl is of size n
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void tdrl_apply( genom &g, const vector<int> &tdrl ){
	genom res;
	int j=0;

	if( tdrl.size() != g.size()+1 ){
		cerr << "error: |tdrl| != n"<<endl;
		exit(1);
	}

	res = g;

	for( int copy=0; copy<=1; copy++ ){
		for(unsigned i=0; i<g.size(); i++){
			if( tdrl[ abs(g[i]) ] == copy ){
				res[j] = g[i];
				j++;
			}
		}
	}
	g = res;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute number of chains in g wrt id
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_chaincnt( const genom &g ){
	genom id(g.size(), g.circular);
	return tdrl_chaincnt(g, id);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute number of chains in g wrt h
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_chaincnt( const genom &g, const genom &h ){
	unsigned cc=1;		// chain count
	vector<int> ig;		// inverse of g

	ig = g.inverse();	// get the inverse

	if( h[0] - g[ abs(ig[ abs(h[0]) ] )] != 0 ){
		return std::numeric_limits< unsigned >::max();
	}

	for(unsigned i=1; i<g.size(); i++){				// for all adjacencies in h
		if( h[i] - g[abs(ig[ abs(h[i] ) ])] != 0 ){
			return std::numeric_limits< unsigned >::max();
		}

		if( abs(ig[ abs(h[i] ) ]) < abs(ig[ abs(h[i-1]) ]) ){	// if \g^{-1}(h(i)) < \g^{-1}(h(i-1))
			cc++;									// .. i.e. g(h(i)) left of g(h(i-1))
		}											// .. start a new chain
	}
	return cc;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// the chainindex, i.e. for each element 1..n in which chain it is
// @param[in] chains the chains
// @return the chain index (0 is undefined)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<int> tdrl_chain_index( const vector<vector<int> > &chain, unsigned n ){
	vector<int> chainidx = vector<int>( n+1, std::numeric_limits< int >::max() );

	for( unsigned i=0; i<chain.size(); i++ ){
		for(unsigned j=0; j<chain[i].size(); j++){
			chainidx[ chain[i][j] ] = i;
		}
	}
	return chainidx;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute chains in g wrt id
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl_chains( const genom &g, vector<vector<int> > &chains ){
	genom id(g.size(), g.circular);
	tdrl_chains(g, id, chains);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute chains in g wrt h
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl_chains( const genom &g, const genom &h, vector<vector<int> > &chains ){
	vector<int> ig;		// inverse of g

	chains = vector<vector<int> >( 1, vector<int>(1, abs(h[0])) );	// start the first chain
	ig = g.inverse();	// get the inverse
	for(unsigned i=1; i<h.size(); i++){				// for all adjacencies in h
		if( abs(ig[ abs(h[i] ) ]) < abs(ig[ abs(h[i-1]) ]) ){	// if \g^{-1}(h(i)) < \g^{-1}(h(i-1))
			chains.push_back( vector<int>() );		// .. i.e. g(h(i)) left of g(h(i-1))
		}											// .. start a new chain
		chains.back().push_back(abs(h[i]));
	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute tdrl distance from g to id
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_distance( const genom &g ){
	return tdrl_distance( tdrl_chaincnt(g) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute tdrl distance from g to h
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_distance( const genom &g, const genom &h ){
	return tdrl_distance( tdrl_chaincnt(g, h) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute tdrl distance from chain count
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_distance( unsigned cc ){
	if( cc == std::numeric_limits< unsigned >::max() ){
		return std::numeric_limits< unsigned >::max();
	}else if(cc > 0){
		return (int)ceil(log2((float)cc));
	}else{
		return 0;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the number of sorting tdrls given number of chains (i.e. wo. breaking chains)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_sorting_cnt( unsigned c ){
	unsigned d,			// tdrl distance
		cnt = 0;

	if( c < 2 ){
		return 0;
	}

	d = tdrl_distance(c);

	for( int i=c-(ppow(d-1)); i<=(int)floor(c/2.0); i++){
		cnt += binom(c+1, 2*i+1);
	}
	return cnt;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the number of sorting tdrls for g wrt. id (wo. breaking chains)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_sorting_cnt( const genom &g ){
	return tdrl_sorting_cnt( tdrl_chaincnt(g) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the number of sorting tdrls for g wrt. h (wo. breaking chains)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_sorting_cnt( const genom &g, const genom &h ){
	return tdrl_sorting_cnt( tdrl_chaincnt(g, h) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute a tdrl given transition positions
// @param[in] length of the permutation
// @param[in] chain the chains
// @param[in] transitions ordered sequence of transitions
// @param[in] state start state 0/1
// @return the resulting tdrl (specifying for each element in which copy it is)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<int> tdrl_from_transitions( int n, const vector<vector<int> > &chain, const vector<int> &transition, char state){
	vector<int> tdrl(n+1, std::numeric_limits< int >::max());

	for( unsigned i=0, t=0; i<chain.size(); i++){
//		cout << (int)state;
		for( unsigned j=0; j<chain[i].size(); j++ ){
			tdrl[ chain[i][j] ] = state;
		}
		if( t < transition.size() && i == (unsigned)transition[t] ){
			state = (state == 0) ? 1 : 0;
			t++;
		}
	}
//	cout << endl;


//	cout << "stdrl ";
//	for(unsigned j=0; j<tdrl.size(); j++){
//		cout << (int)tdrl[j];
//	}
//	cout << endl;
	return tdrl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the sorting tdrls given number of chains (i.e. wo. breaking chains)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<vector<int> > tdrl_sorting( const vector<vector<int> > &chain, unsigned n ){
	unsigned c = chain.size(),
		d = tdrl_distance(c);
	vector<int> transition;
	vector<vector<int> > tdrl;	// the sorting tdrls


	for( unsigned i=c-(ppow(d-1)); i<=(unsigned)floor(c/2.0); i++){
//		cout <<"i " <<i << endl;

		// t_1 = t_n = 0/1 -> 2k transitions at n-1 positions
//		cout << "choose("<<2*i<<","<<c-1<<") = "<<binom(c-1, 2*i)<<endl;
		if( binom(c-1, 2*i) > 0){
			counter_init( transition, 2*i, 0, true );
			do{
//				cout << "A transitions ";copy(transition.begin(), transition.end(), ostream_iterator<int>(cout," ")); cout << endl;
				tdrl.push_back( tdrl_from_transitions( n, chain, transition, 0) );
				tdrl.push_back( tdrl_from_transitions( n, chain, transition, 1) );
				counteradd( transition, 0, c-1, 0, true );
			}while(counter_valid(transition));
		}
		// t_1 = 0 t_n = 1 -> 2k-1 transitions at n-1 positions
//		cout << "choose("<<2*i-1<<","<<c-1<<") = "<<binom(c-1, 2*i-1)<<endl;
		if( binom(c-1, 2*i-1) > 0 ){
			counter_init( transition, 2*i-1, 0, true );
			do{
//				cout << "B transitions ";copy(transition.begin(), transition.end(), ostream_iterator<int>(cout," ")); cout << endl;
				tdrl.push_back( tdrl_from_transitions( n, chain, transition, 0) );
				counteradd( transition, 0, c-1, 0, true );
			}while(counter_valid(transition));
		}
		// t_1 = 1 t_n = 0 -> 2k+1 transitions at n-1 positions
//		cout << "choose("<<2*i+1<<","<<c-1<<") = "<<binom(c-1, 2*i+1)<<endl;
		if( binom(c-1, 2*i+1) > 0 ){
			counter_init( transition, 2*i+1, 0, true );
			do{
//				cout << "C transitions ";copy(transition.begin(), transition.end(), ostream_iterator<int>(cout," ")); cout << endl;
				tdrl.push_back( tdrl_from_transitions( n, chain, transition, 1) );
				counteradd( transition, 0, c-1, 0, true );
			}while(counter_valid(transition));
		}
	}
//	cout << tdrl.size()<<" sorting"<<endl;
	return tdrl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the sorting tdrls for g wrt. id (wo. breaking chains)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<vector<int> > tdrl_sorting( const genom &g ){
	vector<vector<int> > chain;
	tdrl_chains(g, chain);	// get the chains
	return tdrl_sorting(chain, g.size());
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the sorting tdrls for g wrt. h (wo. breaking chains)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<vector<int> > tdrl_sorting( const genom &g, const genom &h ){
	vector<vector<int> > chain;

	tdrl_chains(g, h, chain);	// get the chains
	return tdrl_sorting(chain, g.size());

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the dynamic programming matrix for the (really) all sorting tdrls
// case
// @param[in] chain the chains
// @param[in] n the length of the genomes
// @param[out] a the dynamic programming matrix ..
//    - access [n][k][0/1]
//    - for k: 0 is at position maxred, i.e. to get k access [maxred+k]
//      e.g. maxred = 3
//      a[i]  0  1  2 3 4 5 6 7 ...
//        k  -3 -2 -1 0 1 2 3 4 ...
// @param[out] pos the position string
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl_sorting_all_dp( const vector<vector<int> > &chain, unsigned n, vector<vector<vector<int> > > &a, vector<int> &pos ){
	pair<int,int> krng;				// range of ks which have to be considered
	int maxred = 0,				// maximal possible reduction of k
		maxinc = 0,				// maximal possible increase of k
		c = chain.size();		// number of chains
//	vector<int> chainidx;		// the chainindex

	//chainidx = tdrl_chain_index(chain, n);	// get the chain index
	//pos = vector<int>(n-1, std::numeric_limits< int >::max());		// init the position string
	for(unsigned i=0; i<chain.size(); i++){
		for( unsigned j=0; j<chain[i].size()-1; j++ ){
			pos.push_back( THETA );
		}
		if( i != chain.size()-1 ){
			pos.push_back( PHI );
		}

		maxinc += (int)floor( chain[i].size() / 2.0 );
	}

	maxred = (unsigned)floor(c/2.0);	// maximal chain reduction
//	maxinc = c*2;						// maximal increase in the number of chains

//	cout << "chains   "<<endl;
//	print_chains(chain);
////	cout << "chainidx "; copy( chainidx.begin()+1, chainidx.end(), ostream_iterator<int>(cout, " ") ); cout << endl;
//	cout << "pos      ";
//	for(unsigned i=0; i<pos.size(); i++){
//		if( pos[i] == THETA )
//			cout << "T";
//		else if ( pos[i] == PHI )
//			cout << "P";
//		else
//			cout << "?";
//	}
//	cout << endl;
//	cout << "maxred "<<maxred<<endl;
//	cout << "maxinc "<<maxinc<<endl;

		// init dp matrix
	a = vector<vector<vector<int> > >( n, vector<vector<int> >( maxred+maxinc+1, vector<int>( 2, 0 ) ) );
	a[0][maxred][0] = 1;	// set a_{1,0}^0 = 1
	a[0][maxred][1] = 1;	// set a_{1,0}^1 = 1

	krng = make_pair(0,0);
	for( unsigned i=0; i<n-1; i++ ){
		for( int k = krng.first; k<= krng.second; k++ ){
//			cout << "calculate ("<<i<<","<<k<<")"<<endl;

			// perpendicular update
			a[i+1][maxred+k][0] += a[i][maxred+k][0];
			a[i+1][maxred+k][1] += a[i][maxred+k][1];
			// "diagonal with same k" update
			if( pos[i] == PHI ){
				a[i+1][maxred+k][0] += a[i][maxred+k][1];
			}else{
				a[i+1][maxred+k][1] += a[i][maxred+k][0];
			}

			// "diagonal with different" update
			if( pos[i] == PHI and a[i][maxred+k][0] != 0 ){
				a[i+1][maxred+k-1][1] += a[i][maxred+k][0];
				if( krng.first == k )
					krng.first--;
			}else if (pos[i] == THETA and a[i][maxred+k][1] != 0){
				a[i+1][maxred+k+1][0] += a[i][maxred+k][1];
				if( krng.second == k )
					krng.second++;
			}
		}

//		cout << "range "<<krng.first <<".."<<krng.second<<endl;
	}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helper function: construct sorting tdrls from the dp matrix
// @param i,k,c the indices where to start
// @param a the dp
// @param pos the position string
// @param[in] inv inverse of g
// @param[in] ctdrl temp memory fot a tdrl (must be allocated allredy)
// @param[out] the resulting tdrls
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl_from_dp( int n, int k, int c, const vector<vector<vector<int> > > &a, const vector<int> &pos,
	int maxred, const genom &h, vector<int> &ctdrl, vector<vector<int> > &tdrls  ){


//	cout << "("<<n<<","<<k<<"["<<maxred+k<<"]"<<","<<c<<")";cout.flush();
	if( (maxred+k < 0) or (maxred+k >= (int)a[0].size()) ){
//		cout << "out of range"<<endl;
		return;
	}

	if( a[n][maxred+k][c] <= 0 ){
//		cout <<" "<< c<<" <= 0 ("<< a[n][maxred+k][0] << a[n][maxred+k][1]<<")" <<endl;
		return;
	}
//	cout << endl;
	ctdrl[ abs(h[n]) ] = c;

	if( n == 0 ){
		tdrls.push_back( ctdrl );
	}else{

//		ctdrl.push_back( c );
		tdrl_from_dp( n-1, k, c, a, pos, maxred, h, ctdrl, tdrls );
//		ctdrl.pop_back();

		if( pos[n-1] == PHI and c == 0 ){
//			ctdrl.push_back( 1 );
			tdrl_from_dp( n-1, k, 1, a, pos, maxred, h, ctdrl, tdrls );
//			ctdrl.pop_back();
		}else if( pos[n-1] == PHI and c == 1 ){
//			ctdrl.push_back( 0 );
			tdrl_from_dp( n-1, k+1, 0, a, pos, maxred, h, ctdrl, tdrls );
//			ctdrl.pop_back();
		}else if(pos[n-1] == THETA and c == 0 ){
//			ctdrl.push_back( 1 );
			tdrl_from_dp( n-1, k-1, 1, a, pos, maxred, h, ctdrl, tdrls );
//			ctdrl.pop_back();
		}else if( pos[n-1] == THETA and c == 1 ){
//			ctdrl.push_back( 1 );
			tdrl_from_dp( n-1, k, 0, a, pos, maxred, h, ctdrl, tdrls );
//			ctdrl.pop_back();
		}

	}
//	ctdrl.pop_back();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the number of sorting tdrls for g wrt. h with the dp approach
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_sorting_cnt_all_dp( const genom &g, const genom &h ){
	int n = g.size(),		// length of the genome
		c = 0,					// number of chains
		d = 0,					// tdrl distance
		cnt = 0,				// number of sorting tdrls
		maxred = 0;				// maximal possible chain reduction
	vector<int> pos;
	vector<vector<int> > chain;	// the chains
	vector<vector<vector<int> > > a; // the dp matrix

	if( g==h ){
		return 0;
	}

	tdrl_chains(g, h, chain);				// get the chains
	c = chain.size();
	d = tdrl_distance(c);

	tdrl_sorting_all_dp(chain, n, a, pos);		// do the dp


	maxred = (int)floor(c/2.0);			// maximal chain reduction

	// fetch the results
	for( int k= -1*( c-(ppow(d-1)) ) ; k>= -1*( (int)floor(c/2.0) ); k--){
//		cout << "add "<<k<<endl;
		cnt += a[n-1][maxred+k][0] + a[n-1][maxred+k][1];
	}

	return cnt;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the number of sorting tdrls given n and c
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_sorting_cnt_all( unsigned n, unsigned c ){

	unsigned cnt = 0;
	for( unsigned i=0; i<=(unsigned)ppow((unsigned)ceil(log2(c)))-c; i++ ){
		cnt += binom(n,i);
	}
	return cnt;

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the number of sorting tdrls as sum of binomials
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned tdrl_sorting_cnt_all( const genom &g, const genom &h ){
	return tdrl_sorting_cnt_all(g.size(), tdrl_chaincnt(g));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// compute the sorting tdrls for g wrt. h
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<vector<int> > tdrl_sorting_all( const genom &g, const genom &h ){

	int n = g.size(),		// length of the genome
		c = 0,					// number of chains
		d = 0,					// tdrl distance
		maxred = 0;				// maximal possible chain reduction
	vector<int> tmpt;			// tmp tdrl
	vector<int> pos;			// the position string
	vector<vector<int> > tdrl;	// the sorting tdrls
	vector<vector<int> > chain;	// the chains
	vector<vector<vector<int> > > a; // the dp matrix

	if(g==h){
		return tdrl;
	}

	tdrl_chains(g, h, chain);				// get the chains
	c = chain.size();
	d = tdrl_distance(c);

	tdrl_sorting_all_dp(chain, n, a, pos);	// do the dp
	maxred = (int)floor(c/2.0);	// maximal chain reduction

//	cout << "maxred " <<maxred<<endl;
//	cout << endl;
//	print_dp(a);
//	cout << endl;

	tmpt = vector<int>( g.size()+1, std::numeric_limits< int >::max() );
	// fetch the results
	for( int k= -1*( c-(ppow(d-1)) ) ; k>= -1*( (int)floor(c/2.0) ); k--){
		tdrl_from_dp( n-1, k, 0, a, pos, maxred, h, tmpt, tdrl);
		tdrl_from_dp( n-1, k, 1, a, pos, maxred, h, tmpt, tdrl);
	}

//	cout << "found "<<tdrl.size() <<" sorting tdrls"<<endl;
//	for(unsigned i=0; i<tdrl.size(); i++){
//		for(unsigned j=0; j<tdrl[i].size(); j++)
//			cout << (int)tdrl[i][j];
//		cout << endl;
//	}
	return tdrl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// internal function: recursively get all sorting tdrl scenarios
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl_sort_recursion(const genom &g, const genom &h, vector<vector<vector<int> > > &scenarios,
		vector<vector<int> > &ctdrl, bool restricted, unsigned maxscen){
	genom tmp; 					// tmp genom for applying tdrls
	vector<vector<int> > stdrl;	// the sorting tdrls

	if( maxscen > 0 && scenarios.size() >= maxscen ){
		return;
	}

	if( g == h ){
		scenarios.push_back(ctdrl);
		return;
	}

	if ( restricted ){
		stdrl = tdrl_sorting(g, h);
	}else{
		stdrl = tdrl_sorting_all(g, h);
	}

	// shuffle the sorting tdrls .. just in case the number of sorting tdrls
	// is bigger then maxscen
	random_shuffle(stdrl.begin(), stdrl.end());
	for(unsigned i=0; i<stdrl.size(); i++){
		tmp = g;
		tdrl_apply(tmp, stdrl[i]);
		ctdrl.push_back(stdrl[i]);
		tdrl_sort_recursion(tmp, h, scenarios, ctdrl, restricted, maxscen);
		ctdrl.pop_back();
	}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// get all sorting tdrl scenarios
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdrl_sort( const genom &g, const genom &h, vector<vector<vector<int> > > &scenarios, bool restricted, unsigned maxscen){
	if(g==h){
		return;
	}

	vector<vector<int> > ctdrl;
	tdrl_sort_recursion(g, h, scenarios, ctdrl, restricted, maxscen);

//	for( unsigned i=0; i<scenarios.size(); i++ ){
//		cout << "scenario "<< i <<" size "<<scenarios[i].size()<<endl;
//	}
}

















/**
 * internal function transforming a median order into a set of permutations
 * which are consistent with the median order
 * @param[in] mo median orders mo[i] stores all nodes which come after node i
 * @param[out] cperm consistent permutations are appended(!) here
 * @param[in] curinv for internal usage should be initialized to length of mo
 * @param[in] npos for internal usage .. current possition
 */
void mo2permutation( vector<set<int> > &mo, set<vector<int> > &cperm, vector<int> &curinv, unsigned npos );

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//int cnt_blocks(const genom &quot){
//	int b = 0;
//
//	if( quot[0] < 0){
//		b++;
//	}
//	for(unsigned j=1; j<quot.size(); j++){
//		if(quot[j-1] > 0 &&  quot[j] < 0){
//			b++;
//		}
//	}
//	return b;
//
////	for(unsigned j=1; j<quot.size(); j++){
////		if(quot[j-1] * quot[j] < 0){
////			b++;
////		}
////	}
////
////	return b/2+b%2;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void get_blocks(const genom &quot, const genom &tgt, vector<pair<int,int> > &blocks){
	int start = 0;
	vector<bool> same( quot.size(), true );
	vector<int> tinv(tgt.size()+1, std::numeric_limits< int >::max());	// inverse permutation of tgt

//	cout << "get blocks : "<<quot<<" & "<<tgt << endl;

	// preprocessing:
	// - get inverse of tgt
	// - get a bool vector telling us if element i in quot has the same sign in tgt
	for( unsigned i=0; i<tgt.size(); i++ ){
		if( tgt[i] < 0 )
			tinv[ abs(tgt[i]) ] = -1;
		else
			tinv[ abs(tgt[i]) ] = 1;
	}
	for( unsigned i=0; i<quot.size(); i++ ){
		if( (quot[i] < 0 && tinv[ abs(quot[i]) ] > 0) || (quot[i] > 0 && tinv[ abs(quot[i]) ] < 0) ){
			same[i] = false;
		}
	}
//	cout << "tinv "; copy(tinv.begin(), tinv.end(), ostream_iterator<int>(cout," "));cout << endl;
//	cout << "same "; copy(same.begin(), same.end(), ostream_iterator<bool>(cout," "));cout << endl;


	blocks.clear();
	if( ! same[0] ){
		start = 0;
	}

	for(unsigned j=1; j<quot.size(); j++){
		if( !same[j-1] &&  same[j] ){
			blocks.push_back( make_pair(start, j-1) );
		}
		if( same[j-1] &&  !same[j] ){
			start = j;
		}
	}

	if( ! same.back() ){
		blocks.push_back( make_pair(start, quot.size()-1) );
	}
//	for(unsigned j=1; j<quot.size(); j++){
//		if(quot[j-1] * quot[j] < 0){
//			b++;
//		}
//	}
//
//	return b/2+b%2;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void get_blocks_inv(const genom &quot, int pos, vector<pair<int,int> > &blocks){
//	genom id(quot.size(), 0);
//
//	// mark the elements of the identity permutation
//	// negative which are negative in quot
//	for(unsigned i=0; i<quot.size(); i++){
//		if( quot[i] < 0 ){
//			id[ abs( quot[i] )-1 ] *= -1;
//		}
//	}
//	cout << "get_blocks_inv("<<quot<<")"<<endl;
//	cout << "    ID "<<id << endl;
//	// get the blocks in the (now signed) identity permutation
//	get_blocks( id, pos, blocks );
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool is_tdl( const vector<bool> &tdl, unsigned s, unsigned e ){
	int zo = 0;	// number of 0->1 changes

//	copy(tdl.begin(), tdl.end(), ostream_iterator<bool>(cout, "")); cout<<" " << s<<" "<<e<<endl;
// 	cout << "checking ";
	for(unsigned i=s; (i<tdl.size()-1) && (i<=e) ; i++){
//		cout << i<<tdl[i]<<" ";
		if( !tdl[i] && tdl[i+1] ){
			zo++;
			if(zo > 1){
//				cout << " zo "<<zo<<endl;
				return true;
			}
		}
	}
//	cout << " zo "<<zo<<endl;
	return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void mark_tdl(const genom &g, vector<int> &marks){

	int take = 1;
	marks[0] = 1;
	for(unsigned i=1; i<g.size(); i++){
		if(abs(g[i]) < abs(g[i-1])){
			take*=-1;
		}
		if(take > 0){
			marks[i] = 1;
		}else{
			marks[i] = 0;
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void mark_tdl(const genom &origin,const genom &goal, vector<int> &marks){

    genom id=origin.identify_g(goal);
    return mark_tdl(id, marks);

}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// internal function transforming a median order into a set of permutations
// which are consistent with the median order
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void mo2permutation( vector<set<int> > &mo, set<vector<int> > &cperm, vector<int> &curinv, unsigned npos ){
	vector<unsigned> iecnt;	// indegree for each edge

	if( npos == curinv.size()-1 ){
		vector<int> nperm( npos, 0 );
//		copy( curinv.begin(), curinv.end(), ostream_iterator<int>(cout, " ") ); cout << endl;
		for( unsigned i=1; i<curinv.size(); i++ ){
			nperm[ curinv[i] ] = i;
		}
//		copy( nperm.begin(), nperm.end(), ostream_iterator<int>(cout, " ") ); cout << endl;
		cperm.insert(nperm);
		return;
	}
//	cout << "npos "<<npos << " curinv "; copy(curinv.begin(), curinv.end(), ostream_iterator<int>(cout," ")); cout <<endl;
		// determine the indegree for each node
	iecnt=vector<unsigned>(mo.size(), 0);
	for(unsigned i=1; i<mo.size(); i++){
		if( curinv[i] != std::numeric_limits< int >::max() ){
			continue;
		}
		for( set<int>::iterator it=mo[i].begin(); it!=mo[i].end(); it++ ){
			iecnt[*it]++;
		}
	}

		// start recursion for the node with indegree = 0
	for( unsigned i=1; i<iecnt.size(); i++ ){
		if( iecnt[i] > 0 || curinv[i] != std::numeric_limits< int >::max() ){
			continue;
		}
//		cout << "add "<<i<< " to pos "<<npos<<endl;
		curinv[ i ] = npos;
		mo2permutation( mo, cperm, curinv, npos+1 );
		curinv[ i ] = std::numeric_limits< int >::max();
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdl_normalform( const vector<bool> &tdrl, int &start, int &end){

		// determine the 1st element that is kept in the 2nd copy
		// and the last element that is kept in the 1st copy
	end = -1;
	start = tdrl.size();
	for(int i=0; i<(int)tdrl.size(); i++){
		if( tdrl[i]==false ){
			start = min(start, i);
		}
		if( tdrl[i] == true ){
			end = max(end, i+1);
		}
	}
//	if(start == std::numeric_limits< unsigned >::max() || end == std::numeric_limits< unsigned >::max()){
//		cerr << "internal error: could not determine NF of tdrl: ";
//		copy(tdrl.begin(), tdrl.end(), ostream_iterator<bool>(cerr, "") ); cerr << endl;
//		exit(0);
//	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<bool> random_tdl( unsigned n){
	unsigned s, l;
	vector<bool> t;
	l = ask_rng( 2, n );
	s = ask_rng(0, n-l);
//	cout << "random_tdl start "<<s<<" length "<<l<<endl;
	t = vector<bool>( n, false );
	for(unsigned i=0; i<s; i++){
		t[i] = true;
	}
	for( unsigned i=s; i<s+l; i++){
		if( ask_rng(0,1)==0 ){
			t[i] = false;
		}else{
			t[i] = true;
		}
	}
	for(unsigned i=s+l; i<n; i++){
		t[i] = false;
	}
	if( t.size() != n ){
		cerr << "internal error: random tdrl has a size different from the size of the genome"<<endl;
		exit(1);
	}

	return t;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void reverse_blocks(genom g, const genom &tgt, vector<pair<int, int> > blocks,
	vector<vector<vector<bool> > > &tdlsce,	vector<vector<pair<int, int> > > &revszen,
	map<genom, bool> &used, unsigned maxszen, vector<pair<int, int> > &currz, int &bestd ){

//	cout << "reverse_blocks "<< blocks.size()<<endl;
	if( tdlsce.size() >= maxszen && maxszen > 0 ){
//		cout <<"return "<<tdlsce.size()<<" >= "<<maxszen<<endl;
		return;
	}

	map<genom, bool>::iterator used_it;	// tmp iterator for searching in the map
	pair<int, int> rev;					// temp reversal
	vector<pair<int,int> > new_blocks; 	// negative blocks after a reversal
	vector<vector<bool> > tdls;			// a tdrl scenario

	for(unsigned i=0; i<blocks.size(); i++){
		for(unsigned j=i; j< blocks.size(); j++){

			rev.first = abs(blocks[i].first);
			rev.second =  abs(blocks[j].second);
			reverse(g, rev);
//			cout << "reverse "<< rev.first << ","<<rev.second<<endl;
				// check if the genome was visited already, if then undo the
				// reversal and go to the next
			used_it = used.find(g);
			if(used_it != used.end()){
				reverse(g, rev);
				continue;
			}
				// if the genome was NOT visited before, mark it and process it
			used[g] = true;
			currz.push_back( rev );
			get_blocks(g, tgt, new_blocks);

			if(new_blocks.size() > 0){
				reverse_blocks( g, tgt, new_blocks, tdlsce, revszen, used, maxszen, currz, bestd );
			}else{
				int d =tdrl_distance(g, tgt);

//				cout << "D "<<d<<endl;

				if( d <= bestd ){
					if( d < bestd ){
						bestd = d;
						tdlsce.clear();
						revszen.clear();
					}

					if( maxszen == 0 || tdls.size() <= maxszen ){
						tdl_sort( g, tgt, tdls );
						tdlsce.push_back( tdls );
						tdls.clear();
						revszen.push_back(currz);
					}
				}
			}
			new_blocks.clear();
			reverse(g, rev);
			currz.pop_back();
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void reverse_blocks_inv(const genom &g, genom tgt, vector<pair<int, int> > blocks,
	vector<vector<vector<bool> > > &tdlsce,	vector<vector<pair<int, int> > > &revszen,
	map<genom, bool> &used, unsigned maxszen, vector<pair<int, int> > &currz, int &bestd ){

	if( tdlsce.size() >= maxszen && maxszen > 0 ){
		return;
	}

//	cout << "reverse_blocks_inv"<<endl;
//	cout << g<<endl;
//	cout << tgt<<endl;

	map<genom, bool>::iterator used_it;
	pair<int, int> rev;
	vector<pair<int,int> > new_blocks; 	// negative blocks after a reversal
	vector<vector<bool> > tdls;			// a tdrl scenario

	for(unsigned i=0; i<blocks.size(); i++){
		for(unsigned j=i; j< blocks.size(); j++){

			rev.first = abs(blocks[i].first);
			rev.second =  abs(blocks[j].second);

			reverse(tgt, rev);
//			cout << "rev "<<rev.first<<","<<rev.second<<" -> "<<tgt << endl;

				// check if the genome was visited already, if then undo the
				// reversal and go to the next
			used_it = used.find(tgt);
			if(used_it != used.end()){
				reverse(tgt, rev);
				continue;
			}
				// if the genome was NOT visited before, mark it and process it
			used[tgt] = true;

			currz.insert(currz.begin(), rev);
			get_blocks(tgt, g, new_blocks);

			if(new_blocks.size() > 0){
				reverse_blocks_inv( g, tgt, new_blocks, tdlsce, revszen, used, maxszen, currz, bestd );
			}else{
				int d =tdrl_distance(g, tgt);

				if( d <= bestd ){
					if( d < bestd ){
						bestd = d;
						tdlsce.clear();
						revszen.clear();
//						out << "-----------------------<br>";
					}
					if( maxszen == 0 || tdls.size() < maxszen ){
						tdl_sort( g, tgt, tdls );
						tdlsce.push_back( tdls );
//						cout << "append "<<tdls.size()<<"tdls and "<<endl;
						tdls.clear();
						revszen.push_back(currz);
//						for(unsigned k=0; k<currz.size(); k++)
//							cout << "rev("<< currz[k].first<<","<<currz[k].second<<")"<<endl;
					}
				}
			}

			new_blocks.clear();
			reverse(tgt, rev);
			currz.erase( currz.begin() );
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// apply a tdl step; the keep vector specifies which genes are kept in the
// first copy (1), and which are kept in the second copy (0)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void tdl( genom &g, const vector<bool> &keep ){
	genom res;
	int j=0;

	if( keep.size() < g.size() ){
		cerr << "tdl: keep vector smaller than genome"<<endl;
		exit(1);
	}

	res = genom(g.size(), 0);
	for(unsigned i=0; i<g.size(); i++){
		if( keep[i] == true ){
			res[j] = g[i];
			j++;
		}
	}
	for(unsigned i=0; i<g.size(); i++){
		if( keep[i] == false ){
			res[j] = g[i];
			j++;
		}
	}
	g = res;
}



/*void tdl_bfmedian( const genom &g, const genom &h, const counter &distances, int n, unsigned cyccnt,
	const vector<vector<int> > &gecm, const vector<vector<int> > &hecm, set<vector<int> > &cperm ){

	vector<int> curinv; 		// inverse for solution construction
	vector<counter> cpos(2); 	// the cut positions
	vector<bool> cyccover;			// for determining if all cycles are covered
	vector<set<int> > tmo;
	vector<vector<set<int> > > mo;	// median orders

	cyccover = vector<bool>(cyccnt, false);
	curinv = vector<int>( g.size()+1, std::numeric_limits< int >::max() );

	#ifdef DEBUG_TDLBFMEDIAN
	cout << "cut ranges "<< max(0, ppow(distances[0]-1)) <<" - " << min(n-1, ppow(distances[0])-1)<<endl;
	cout << "           "<< max(0, ppow(distances[1]-1)) <<" - " << min(n-1, ppow(distances[1])-1)<<endl;
	#endif//DEBUG_TDLBFMEDIAN
	for( int gcuts = max(0, ppow(distances[0]-1)); gcuts <= min(n-1, ppow(distances[0])-1); gcuts++ ){
		for( int hcuts = max(0, ppow(distances[1]-1)); hcuts <= min(n-1, ppow(distances[1])-1); hcuts++ ){
			#ifdef DEBUG_TDLBFMEDIAN
			cout << "gcuts "<< gcuts<<" hcuts "<<hcuts<<endl;
			#endif//DEBUG_TDLBFMEDIAN

			cpos[0] = counter(gcuts, n-1, 0, true);
//			counter_init(cpos[0], gcuts, 0, true);
			do{
				cpos[1] = counter(hcuts, n-1, 0, true);
//				counter_init(cpos[1], hcuts, 0, true);
				do{
//					#ifdef DEBUG_TDLBFMEDIAN
//					cout << "counters ";
//					copy( cpos[0].begin(), cpos[0].end(), ostream_iterator<int>(cout, " ") ); cout << " .. ";
//					copy( cpos[1].begin(), cpos[1].end(), ostream_iterator<int>(cout, " ") ); cout << endl;
//					#endif//DEBUG_TDLBFMEDIAN
					for(unsigned i=0; i<gcuts; i++){
						for(unsigned k=0; k<gecm[cpos[0][i]].size(); k++){
							cyccover[ gecm[cpos[0][i]][k] ] = true;
						}
					}
					for(unsigned i=0; i<hcuts; i++){
						for(unsigned k=0; k<hecm[cpos[1][i]].size(); k++){
							cyccover[ hecm[cpos[1][i]][k] ] = true;
						}
					}

					if( sum(cyccover) == cyccnt ){
						tmo = vector<set<int> >( n+1 );
						#ifdef DEBUG_TDLBFMEDIAN
						cout << "complete covered by ("; copy(cpos[0].begin(), cpos[0].end(), ostream_iterator<int>(cout,","));cout << ") ("; copy(cpos[1].begin(), cpos[1].end(), ostream_iterator<int>(cout,","));cout<<")" <<endl;

						cout << "c "<<cpos[0].size()<<" " <<cpos[1].size()<<endl;
						#endif//DEBUG_TDLBFMEDIAN
						for( int k=0, p=0; k<n-1; k++ ){
							if( p >= (int)gcuts || k != cpos[0][p] ){
								tmo[g[k]].insert( g[k+1] );
							}else{
								p++;
							}
						}
						for( int k=0, p=0; k<n-1; k++ ){
							if( p >= (int)hcuts || k != cpos[1][p] ){
								tmo[h[k]].insert( h[k+1] );
							}else{
								p++;
							}
						}
						mo.push_back(tmo);
					}
					cyccover.assign(cyccnt, false);

					cpos[1]++;
//					counter( cpos[1], 0, n-1, 0, true);
				}while( cpos[1].isvalid());

				cpos[0]++;
			}while( cpos[0].isvalid() );
		}
	}

	for( unsigned i=0; i<mo.size(); i++ ){
//		cout << "median order "<<i<<endl;
//		for( unsigned j=1; j<mo[i].size(); j++ ){
//			cout << j << " -> "; copy(mo[i][j].begin(), mo[i][j].end(), ostream_iterator<int>(cout, " ")); cout << endl;
//		}
		mo2permutation( mo[i], cperm, curinv, 0 );
	}
	mo.clear();

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdl_bfmedian(const genom &g, const genom &h, vector<genom> &median ){
	vector<int> ginv,	// inverse permutation of g
		hinv;			// .. and of h
	unsigned cyccnt = 0;		// number of cycles
	int dist = 0,				// smaller of the two distances
		dir;					// direction of the smaller distance
	vector<vector<int> > gecm,	// edge -> cycle (index) for g and h
		hecm; 					// .. means: edge i->i+1 (g[i]->g[i+1]) cuts all the cycles in gecm[i]
	set<vector<int> > cperm;

	#ifdef DEBUG_TDLBFMEDIAN
	cout << "tdl_bfmedian("<<g<<","<<h<<")"<<endl;
	#endif//DEBUG_TDLBFMEDIAN
	if( !g.ispos() || !h.ispos() ){
		cerr << "tdl_1median: permutations must be positive "<<endl;
		cerr << g << endl;
		cerr << h << endl;
		exit(1);
	}

	median.clear();
		// get the inverse permutations
	ginv = g.inverse();
	hinv = h.inverse();
		// init gecm and hecm
	gecm = vector<vector<int> >(g.size()-1);
	hecm = vector<vector<int> >(h.size()-1);
	#ifdef DEBUG_TDLBFMEDIAN
	cout<<"g "; for(unsigned i=0; i<g.size(); i++){cout << g[i]<<" ";} cout << " ginv ";copy(ginv.begin()+1, ginv.end(), ostream_iterator<int>(cout, " "));cout <<endl;
	cout<<"h "; for(unsigned i=0; i<h.size(); i++){cout << h[i]<<" ";} cout << " hinv ";copy(hinv.begin()+1, hinv.end(), ostream_iterator<int>(cout, " "));cout <<endl;
	#endif//DEBUG_TDLBFMEDIAN

		// determine cycles which are defined by decents in g,
		// i.e. adjacencies in g which are in the other order in h
		// if there is a conserved adjaceny in h -> dont save in the cycles
	for( unsigned i=0; i<g.size()-1; i++ ){
//		cout << "adjacency "<<g[i]<<","<<g[i+1]<<" h-> "<<hinv[ g[i] ] <<","<< hinv[g[i+1]]<<endl;
		if( hinv[ g[i] ] > hinv[g[i+1]] ){
//			cout << "cycle"<<endl;
			gecm[ i ].push_back( cyccnt );
			for( int j=hinv[ g[i+1] ]; j<hinv[g[i]]; j++ ){
//				cout << h[j]<<".."<<h[j+1]<<" ? "<< ginv[h[j]]<<".."<<ginv[h[j+1]] <<endl;
				if( ginv[h[j+1]]-ginv[h[j]] != 1 ){	// conserved adjacencies dont cover anything
					hecm[j].push_back( cyccnt );
				}
			}
			cyccnt++;
		}
	}

		// determine cycles which are defined by decents in h,
		// i.e. adjacencies in h which are in the other order in g
		// note if the elements are adjacent in g then they can be skipped
		// as they are allready have been defined by the corresp. decent in g
	for( unsigned i=0; i<h.size()-1; i++ ){
//		cout << "adjacency "<<h[i]<<","<<h[i+1]<<" g-> "<<ginv[ h[i] ] <<","<< ginv[h[i+1]]<<endl;
		if( ginv[ h[i] ] > ginv[h[i+1]]+1 ){
//			cout << "cycle"<<endl;
			hecm[i].push_back( cyccnt );
			for( int j=ginv[ h[i+1] ]; j<ginv[h[i]]; j++ ){
//				cout << g[j]<<".."<<g[j+1]<<" ? "<< hinv[g[j]]<<".."<<hinv[g[j+1]] <<endl;
				if( hinv[g[j+1]]-hinv[g[j]]!=1 ){	// conserved adjacencies dont cover anything
					gecm[j].push_back( cyccnt );
				}
			}
			cyccnt++;
		}
	}

	#ifdef DEBUG_TDLBFMEDIAN
	cout << cyccnt << " cycles "<<endl;
	cout << "g edge -> cycle mapping"<<endl;
	for( unsigned i=0; i<g.size()-1; i++ ){
		cout << g[i]<<"->"<<g[i+1]<<" : ";
		for( unsigned j=0; j<gecm[i].size(); j++ )
			cout << gecm[i][j]<<" ";
		cout << endl;
	}
	cout << "h edge -> cycle mapping"<<endl;
	for( unsigned i=0; i<h.size()-1; i++ ){
		cout << h[i]<<"->"<<h[i+1]<<" : ";
		for( unsigned j=0; j<hecm[i].size(); j++ )
			cout << hecm[i][j]<<" ";
		cout << endl;
	}
	#endif//DEBUG_TDLBFMEDIAN

	counter distances;				// number of cuts to make in the two permutations
	// try to find a median with target score tscore
	// where target score is iterated from to minimal score 1 to its max (= min(dist(g,h),dist(h,g)))
	// @todo <= -> < (if <= then sorting has the same cost)
	dist = min_tdrl_distance(g, h, dir);
	for( int tscore=1; tscore <= dist; tscore++ ){

//		counter_init( distances, 2, 0, false );
		distances = counter( 2, g.size()-1, 0, false );
		do{
			if( accumulate(distances.begin(), distances.end(), 0) != tscore ){
//				cout << "counter != tscore "<<endl;
//				counter( distances, 0, g.size()-1, 0, false );
				distances++;
				continue;
			}
			#ifdef DEBUG_TDLBFMEDIAN
			cout << "score "<< tscore<<" distances "; copy(distances.begin(), distances.end(), ostream_iterator<int>(cout, " ")); cout <<endl;
			#endif//DEBUG_TDLBFMEDIAN
			tdl_bfmedian(g, h, distances, g.size(), cyccnt, gecm, hecm, cperm);
			distances++;
//			counter( distances, 0, g.size()-1, 0, false );
		}while( distances.isvalid() );
		if( cperm.size() > 0 ){
			break;
		}
	}

	for( set<vector<int> >::iterator it = cperm.begin(); it!=cperm.end(); ){
		median.push_back( genom(*it, g.getCircular(), g.get_nmap()) );
		cperm.erase(it++);
	}
	cperm.clear();
//	cyccover = vector<bool>(cyccnt, false);
//	for( unsigned i=0; i<gecm.size(); i++ ){
//		for(unsigned j=0; j<hecm.size(); j++){
//			for(unsigned k=0; k<gecm[i].size(); k++){
//				cyccover[ gecm[i][k] ] = true;
//			}
//			for(unsigned k=0; k<hecm[j].size(); k++){
//				cyccover[ hecm[j][k] ] = true;
//			}
//
////			cout << "("<< g[i]<<","<<g[i+1]<<") ("<<h[j]<<","<<h[j+1]<<") -> "; copy(cyccover.begin(), cyccover.end(), ostream_iterator<bool>(cout, "")); cout<<" = "<<sum(cyccover)<<"/"<<cyccnt<<endl;
//
//			if( sum(cyccover) == cyccnt ){
//				tmo = vector<set<int> >( g.size()+1 );
//				cout << "complete covered"<<endl;
//				for( unsigned k=0; k<g.size()-1; k++ ){
//					if( k != i ){
//						tmo[g[k]].insert( g[k+1] );
//					}
//				}
//				for( unsigned k=0; k<h.size()-1; k++ ){
//					if( k != j ){
//						tmo[h[k]].insert( h[k+1] );
//					}
//				}
//				mo.push_back(tmo);
//			}
//			cyccover.assign(cyccnt, false);
//		}
//	}
//
//	for( unsigned i=0; i<mo.size(); i++ ){
////		cout << "median order "<<i<<endl;
////		for( unsigned j=1; j<mo[i].size(); j++ ){
////			cout << j << " -> "; copy(mo[i][j].begin(), mo[i][j].end(), ostream_iterator<int>(cout, " ")); cout << endl;
////		}
//		mo2permutation( mo[i], cperm, curinv, 0 );
//	}
//	for( set<vector<int> >::iterator it = cperm.begin(); it!=cperm.end(); ){
//		median.push_back( genom(*it, g.getCircular(), g.get_nmap()) );
//		cperm.erase(it++);
//	}
}*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//int tdl_distance(const genom &g){
//	int descents = 0;
////	cout << "TDL DIST --------------"<<endl;
////	cout << g << endl;
//
//	for(unsigned i=1; i<g.size(); i++){
//		if(abs(g[i]) < abs(g[i-1])){
//			descents++;
//		}
//	}
//
////	cout << descents+1<<endl;
////	cout << "--------------"<<endl;
//	if(descents > 0){
//		descents++;
//		return (int)ceil(log2((float)descents));
//	}else{
//		return 0;
//	}
//
//}

//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//int tdl_distance(const genom &origin,const genom &goal){
//    genom id=origin.identify_g(goal);
//    return tdl_distance(id);
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdl_sort( const genom &g, vector<vector<bool> > &tdls ){
	genom id;			// the start of the tdl ..
	int idx = 0,		// current increasing substring index
		d = 0,			// the tdl distance
		p = 1;			// a power of two
	vector<bool> keep;	// the vector specifying which genes to keep in the 1st resp. 2nd copy
	vector<int> isi;	// increasing substring index of the elements

//	cout <<"tdl_sort "<< g<<endl;

	tdls.clear();
	isi = vector<int>(g.size()+1, 0);
	keep = vector<bool>(g.size(), false);
	id = genom(g.size(), 0);

		// compute increasing substring indices
	isi[ abs(g[0]) ] = 0;
	for(unsigned i=1; i<g.size(); i++){
		if(abs(g[i]) < abs(g[i-1])){
			idx++;
		}
		isi[ abs(g[i]) ] = idx;
	}
		// get the tdrl distance
	d = (int)ceil(log2((float)(idx+1)));

//	cout << "tdl distance "<<d<<endl;
//	cout << "isi ";copy(isi.begin(), isi.end(), ostream_iterator<int>(cout, " "));cout << endl;

	for( int i=0; i<d; i++ ){
//		cout << i+1<<" th step"<<endl;
		keep.assign( g.size(), false );
		p = ppow(i);
//		cout << "pow "<<p<<endl;
		for(unsigned k=0; k<id.size(); k++){
//			cout << isi[id[k]]<<" & "<<p << " = "<<(isi[ id[k] ] & p)<<endl;;
			if( (isi[ id[k] ] & p) == 0 ){	// kth bit is 0
				keep[k] = true;
			}
		}
//		copy( keep.begin(), keep.end(), ostream_iterator<bool>(cout, " ") ); cout << endl;
		tdl(id, keep);
//		cout << "id "<<id<<endl;
		tdls.push_back(keep);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void tdl_sort(const genom &origin, const genom &goal, vector<vector<bool> > &tdls) {
	genom id;

	tdls.clear();

	id = origin.identify_g(goal);

//	cout <<"o "<< origin<<endl;
//	cout <<"g "<< goal<<endl;
//	cout << "i "<<id<<endl;

	tdl_sort( id, tdls );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pair<int, int> tdl_transposition_count( const vector<vector<bool> > &tdls ){
	pair<int, int> ttc(0,0);	// count of transpositions / tdls

	for( unsigned i=0; i<tdls.size(); i++){
		if( is_tdl(tdls[i]) ){
			ttc.first++;
		}else{
			ttc.second++;
		}
	}

	return ttc;
}

