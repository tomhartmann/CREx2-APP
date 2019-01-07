/**
 * @file genom.cpp
 * genom functions
 * @author M. Bernt
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>


//#include "conserved.hpp"
#include "counter.hpp"
#include "genom.hpp"
#include "helpers.hpp"

//	// includes for condensation (from grappa)
//#include "condense.h"

//#define DEBUG_CYCLES
//#define DEBUG_DISTANCE

using namespace std;






FILE *outfile;

///////////////////////////////////////////////////////////////////////////////
/*!todo false as preset ? */
genom::genom() {
	//~ name = "";
	chromosom.clear();
	circular = false;
	nmap = NULL;
	return ;
}

///////////////////////////////////////////////////////////////////////////////

genom::genom(int cnt, char c, bool inv) {
	chromosom.clear();
	chromosom = vector<int>(cnt, 0);
	for (int i = 1; i <= cnt; i++) {
		chromosom[i-1] = i;
	}
	nmap = NULL;
	circular = c;

	if( inv ){
		reverse (chromosom.begin(), chromosom.end());
		for(unsigned i=0; i<chromosom.size(); i++){
			chromosom[i] *= -1;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

genom::genom(int *genes, int cnt, char c, vector<string> *nm) {
	//~ name = n;
	chromosom.clear();
	for (int i = 0; i < cnt; i++) {
		chromosom.push_back(genes[i]);
	}
	nmap = nm;
	circular = c;

	if(circular)
		normalize();
}

///////////////////////////////////////////////////////////////////////////////

genom::genom(const vector<int> &p, int c, vector<string> *nm, int norm_to){
	chromosom = p;
	circular = c;
	if( c ){
		normalize(norm_to);
	}

	nmap = nm;
}

///////////////////////////////////////////////////////////////////////////////

//genom::genom(int cnt, int dist, int start_len, int end_len, bool steps, bool random, char c, bool unsig, hdata &hd) {
//	nmap = NULL;
//
//	if(random){
//		int idx;
//		vector<int> elements;
//
//			// init the genom at idendity
//		chromosom.clear();
//		circular = c;
//
//		for (int i = 1; i <= cnt; i++)
//			elements.push_back(i);
//
//		while (elements.size()){
//			idx = ask_rng() % elements.size();
//			if(!unsig){
//				if(ask_rng() % 2)
//					elements[idx] *= -1;
//			}
//
//			chromosom.push_back(elements[idx]);
//			elements.erase(elements.begin()+idx);
//		}
//
//	}else{
//		//~ vector<pair<int,int> > reversals;	// vector of possible reversals
//		genom id(cnt, c);			// identity permutation
//		int len,
//			start,
//			end,
//			rounds=0,
//			currDist = 0,		// current distance
//			newDist = 0;		// new distance
//
//		if(dist > cnt){
//			cout << "Could not create a genome with length "<<cnt<<" and distance "<<dist << " to ID"<<endl;
//			exit(1);
//		}
//
//			// init the genom at idendity
//		chromosom.clear();
//		circular = c;
//		chromosom = id.chromosom;
//
//		currDist = 0;			// init the distance
//		while (currDist < dist){							// iterate until we are at the goal distance
//			if(rounds > 10*dist){
//				currDist = 0;
//				newDist = 0;
//				chromosom = id.chromosom;
//				rounds = 0;
//			}
//
//			len = start_len + (ask_rng() % (end_len - start_len));
//
//			if(circular)
//				start = 1+ask_rng() % (size() - len-1);
//			else
//				start = ask_rng() % (size() - len);
//			end = start + len;
//			reverse(*this, start, end);
//
//			if(steps){
//				newDist++;
//			}else{
//				newDist = distance(id, hd);						// calculate the new distance
//			}
//
//			if (newDist <= currDist){						// if it's not better re-reverse
//				reverse(*this, start, end);
//			}
//			else{										// else remember the new permutation
//				currDist = newDist;
//			}
//			rounds++;
//		}
//	}
//}

///////////////////////////////////////////////////////////////////////////////

genom genom_nid(int n){
	genom g;
	g.chromosom = vector<int>(n);

	for (int i = 0; i < n; i++) {
		g.chromosom[i] = -1*(n-i);
	}
	return g;
}

///////////////////////////////////////////////////////////////////////////////

vector<int>::iterator genom::begin(){
	return chromosom.begin();
}

///////////////////////////////////////////////////////////////////////////////

//int breakpoints(const genom &g1, const genom &g2, hdata &hd){
//	int bp = 0;
//	vector<int> g = g1.identify(g2, hd);
//
////	copy(g.begin(), g.end(), ostream_iterator<int>(cout, " ")); cout << endl;
//
//	for(unsigned j=0; j<g.size()-1; j++){
//		if(g[j] >= 0){
//			if(g[j+1] - g[j] != 1){
//				bp++;
//			}
//		}else if(g[j] < 0){
//			if( g[j] - g[j-1] != 1){
//				bp++;
//			}
//		}
//	}
//	return bp;
//}
//
/////////////////////////////////////////////////////////////////////////////////
//
//int breakpoints(const vector<genom> &genomes, hdata &hd){
//
//	vector<bool> breakpoints(genomes[0].size()+1, 0);
//	vector<int> g;
//
//	for( unsigned i=1; i<genomes.size(); i++){
//		g = genomes[0].identify(genomes[i], hd);
//
//		for(unsigned j=0; j<g.size()-1; j++){
//			if( breakpoints[abs(g[j])] == 1)
//				continue;
//
//			if(g[j] >= 0){
//				if(g[j+1] - g[j] != 1){
//					breakpoints[ g[j] ] = 1;
//				}
//			}else if(g[j] < 0){
//				if( g[j] - g[j-1] != 1){
//					breakpoints[ -1*g[j] ] = 1;
//				}
//			}
//		}
//	}
//	//~ cout << "->"<<endl;
//	//~ for (unsigned i = 1; i < breakpoints.size(); i++) {
//		//~ cout.width(3);
//		//~ cout << i << " ";
//	//~ }
//	//~ cout << endl;
//
//	//~ for (unsigned i = 1; i < breakpoints.size(); i++) {
//		//~ cout.width(3);
//		//~ cout << breakpoints[i] << " ";
//	//~ }
//	//~ cout << endl;
//
//	return  sum(breakpoints);
//}

///////////////////////////////////////////////////////////////////////////////

void genom::clear() {
	//~ name = "";
	chromosom.clear();
	circular = false;
}

///////////////////////////////////////////////////////////////////////////////

void deidentify( const genom &map, vector<genom> &genomes ){

	for(unsigned i=0; i<genomes.size(); i++){
		for(unsigned j=0; j<genomes[i].size(); j++){
			if(genomes[i][j] > 0){
				genomes[i][j] = map[ genomes[i][j]-1 ];
			}else{
				genomes[i][j] = -1 * map[ (-1*genomes[i][j])-1 ];
			}
		}
	}

}

///////////////////////////////////////////////////////////////////////////////

void genom::decomposition(vector< set<int> > &deco){
	set<int> temp;
	for(unsigned int i=0; i<chromosom.size()-1; i++){
		// i is end of maximum increasing string
		if(chromosom[i]>chromosom[i+1]){
			temp.insert(chromosom[i]);
			deco.push_back(temp);
			temp.clear();
		}else{ // inside of max increasing substring
			temp.insert(chromosom[i]);
		}
	}
	// last element:
	temp.insert(chromosom[chromosom.size()-1]);
	deco.push_back(temp);
}

///////////////////////////////////////////////////////////////////////////////

void genom::decomposition_strict(vector< set<int> > &deco){
	set<int> temp;
	for(unsigned int i=0; i<chromosom.size()-1; i++){
		// i is end of maximum increasing string
		if(chromosom[i+1] - chromosom[i] != 1){
			temp.insert(chromosom[i]);
			deco.push_back(temp);
			temp.clear();
		}else{ // inside of max increasing substring
			temp.insert(chromosom[i]);
		}
	}
	// last element:
	temp.insert(chromosom[chromosom.size()-1]);
	deco.push_back(temp);
}

///////////////////////////////////////////////////////////////////////////////

unsigned int genom::max_inc_substrings(){
	int count=1;
	for(unsigned int i=0; i<chromosom.size()-1; i++){
		// i is end of maximum increasing string
		if(chromosom[i]>chromosom[i+1]){
			count++;
		}
	}
	return count;
}

////////////////////////////////////////////////////////////////////////////////

//int genom::distance(const genom &g, hdata &hd) const{
//	struct genome_struct g1, g2;
//	//~ static distmem_t * distmem;		// structure for grappa distance computations
//	//~ static bool initialised = false;	// true if the structure is already initialised
//
//		// if distmem is uninitialised
//		// -> initialise it (once per program run)
//	//~ if(!initialised){
//		//~ distmem = new_distmem ( size() );
//		//~ initialised = true;
//	//~ }
//
//	// Test for valid chromosom
//	/*!todo validity check beim einlesen machen*/
//	//~ if (!this->size() || !g.size()) {
//		//~ return 0;
//	//~ }
//
//		// it seem that the vector is direct usable as c-array
//		// it's ugly to cast to int*, but the genes should not
//		// get altered during the distance calculation
//	g1.genes = get_pointer();
//	g2.genes = g.get_pointer();
//
//		// (let) calculate the distance
//	if (circular) {
//		return invdist_circular_mb ( &g1, &g2, size(), hd.distmem );
//		//~ return invdist_circular_nomem(&g1, &g2, this->size());
//	} else{
//		return invdist_noncircular_mb ( &g1, &g2, 0, size(), hd.distmem );
//		//~ return invdist_noncircular_nomem(&g1, &g2, 0, this->size());
//	}
//
//}

///////////////////////////////////////////////////////////////////////////////

//int genom::distance(const vector<genom> &respect, hdata &hd) const{
//	int dist = 0;
//
//	for (unsigned int i = 0; i < respect.size(); i++)
//		dist += distance(respect[i], hd);
//
//	return dist;
//}

///////////////////////////////////////////////////////////////////////////////

vector<int>::iterator genom::end(){
	return chromosom.end();
}

///////////////////////////////////////////////////////////////////////////////

void genom::erase(int i){
	for(int j=chromosom.size()-1; j>=0; j--){
		if(i == abs(chromosom[j]) ){
//			cout << "erase "<< chromosom[j]<<" -> " << *this<<endl;
			chromosom.erase(chromosom.begin() + j);
		}
	}

	for(int j=chromosom.size()-1; j>=0; j--){
		if(abs(chromosom[j]) > i ){
			if(chromosom[j] > 0)
				chromosom[j]--;
			else
				chromosom[j]++;
		}
//		if(chromosom[j] > 0 && i < chromosom[j]){
//			chromosom[j]--;
//		}
//
//		if(chromosom[j] < 0 && i > chromosom[j]){
//			chromosom[j]++;
//		}
//		cout << " -> " << *this<<endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

void genom::evolve(int r, int t, int l){
	int idx0, idx1, idx2;

	if (l==0){
		for (int j = 0; j<r; j++){
			if(t>0 && ask_rng()%t == 0){
				idx0 = ask_rng()%size();
				do{
					idx1= ask_rng()%size();
				}while(idx0 == idx1);
				do{
					idx2= ask_rng()%size();
				}while(idx0 == idx2 || idx1 == idx2);
				transpose(*this, idx0,idx1,idx2);
			}else{

				idx0 = ask_rng()%size();
				idx1= ask_rng()%size();
				reverse(*this, idx0, idx1);
			}
		}
	}else if(l > 0){
		for (int j = 0; j<r; j++){
			if(t>0 && ask_rng()%t == 0){
				idx0 = ask_rng()%(size()-2*l);
				idx1 = idx0 + l;
				idx2 = idx1 + l;
				transpose(*this, idx0,idx1,idx2);
			}else{
				idx0 = ask_rng()%(size()-l);
				idx1 = idx0+l;
				reverse(*this, idx0, idx1);
			}
		}
	}else{
		cerr << "evolve: length < 0"<<endl;
		exit(0);
	}


}

///////////////////////////////////////////////////////////////////////////////

vector<unsigned> genom::genset() const{
	vector<unsigned> gset;
	for(unsigned i=0; i<chromosom.size(); i++){
		gset.push_back( abs(chromosom[i]) );
	}

	sort(gset.begin(), gset.end());
	return gset;
}

///////////////////////////////////////////////////////////////////////////////

vector<string> *genom::get_nmap() const{
	return nmap;
}

///////////////////////////////////////////////////////////////////////////////

int * genom::get_pointer() const{
	return (int*) &(chromosom[0]);
}

///////////////////////////////////////////////////////////////////////////////

vector<int> genom::getChromosom() const{
	return chromosom;
}

///////////////////////////////////////////////////////////////////////////////

char genom::getCircular() const{
	return circular;
}

///////////////////////////////////////////////////////////////////////////////

//vector<int> genom::getComponents(const genom &g, vector<int> &component_orientation, unsigned &component_cnt) const{
//	vector<int> component(size()+1,-1),
//		orientation;
//	vector< pair<int,int> > ici;
//	vector<int> cior;
//	int cur_component_idx = 0;
//
//
//	ici = getIrreducibleConservedIntervals(*this, g, cior);
//
//		// the irreducible conserved intervals (ici) are the spans of the components
//		// - each point belongs to exactly one component
//		// - each point belongs to the ici which contains the point
//		// -> for each ici mark all point in the range from its start to its end with the index of the ici
//		// 	if a point is already marked by another ici check if the current ici is contained in the other ici
//		// 	if so: then overwrite; else leave as is
//	//~ for(unsigned i=0; i<ici.size(); i++){
//	for(int i=ici.size()-1; i>=0; i--){
//			// if component has length 1 -> index 0
//		if(ici[i].first == ici[i].second-1){
//			component[ici[i].first] = 0;
//			continue;
//		}
//		cur_component_idx++;
//			// else label the point
//		for(int j = ici[i].first; j < ici[i].second; j++){
//			component[j] = cur_component_idx;
//		}
//	}
//
//	//~ for( unsigned i=0; i<component.size(); i++)
//		//~ cout << component[i]<<" ";
//	//~ cout << endl;
//
//		// get the orientation of the points of the permutation
//	orientation = getPointOrientations(g);
//	component_cnt = cur_component_idx+1;
//	//~ cout << "cc"<<component_cnt<<endl;
//		// init the orientations of the components (all unoriented)
//		// component 0 (all components of length 1) is always oriented
//	component_orientation = vector<int>(component_cnt, std::numeric_limits< int >::max());
//	component_orientation[0] = 0;
//
//		// get the orientation of each component
//	for(unsigned i=0; i<component.size(); i++){
//			// dont change orientation of components of length 1
//		if(component[i] == 0){
//			continue;
//		}
//
//			// init if first element of the component
//		if(component_orientation[ component[i] ] == std::numeric_limits< int >::max()){
//			component_orientation[ component[i] ] = orientation[i];
//			continue;
//		}
//			// if the orientation of the of the current point is different to the orientation of the points
//			// belonging to same component - mark it as oriented (0)
//		if(component_orientation[ component[i] ] != orientation[i]){
//			component_orientation[ component[i] ] = 0;
//		}
//	}
//		// for convenience - rename orientation entries -1 & 1 to 0 (unoriented) and 0 to 1 (oriented)
//	for(unsigned i=0; i< component_orientation.size(); i++){
//		if(component_orientation[i]==0){
//			component_orientation[i] = 1;
//		}else{
//			component_orientation[i] = 0;
//		}
//	}
//	//~ for (unsigned i=0; i<component_orientation.size(); i++)
//		//~ cout << "i "<<i<<" : "<<component_orientation[i]<<endl;
//
//
//	return component;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool genom::hasduplicates( ) const{
	map<unsigned, unsigned> cnt;

//	vector<unsigned> cnt( nmap->size(), 0 );
	for( unsigned i=0; i<chromosom.size(); i++ ){
		if( cnt.find( abs( chromosom[i] ) ) == cnt.end() ){
			cnt[abs( chromosom[i] )]=0;
		}
		cnt[abs( chromosom[i] )]++;
		if(cnt[abs( chromosom[i] )] > 1){
			return true;
		}
	}
	return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int cycles(const vector<int> &pi, int n, vector<int> &cycs){
	int c = 0, 	// the number of cycles
		nei = 0;	// next elementary interval
	static vector<bool> marked;		// status of the points (was there or not)
	static vector<pair<int,int> > ei;	// the elementary intervals
	static vector<vector<int> > mi;		// the intervals meeting at a point p

	if((int)ei.size() != n+1){
		ei = vector<pair<int, int> > (n+1);
		marked = vector<bool>(n+1, false);
		mi = vector<vector<int> >(n+1, vector<int>(2, -1));
	}else{
		marked.assign(n+1, false);
		for(unsigned i=0; i<mi.size(); i++)
			mi[i].assign(2, -1);
	}

	elementary_intervals(pi, n, ei);
#ifdef DEBUG_CYCLES
	cout << "ei"<<endl;
	for(unsigned i=0; i<ei.size(); i++){
		cout <<"I"<<i<<" "<< ei[i].first <<" "<<ei[i].second<<endl;
	}
#endif
		// get the intervals meeting at a point
	for(unsigned i=0; i<ei.size(); i++){
		if(mi[ ei[i].first ][0] == -1 ){
			mi[ ei[i].first ][0] = i;
		}else{
			mi[ ei[i].first ][1] = i;
		}

		if(mi[ ei[i].second ][0] == -1 ){
			mi[ ei[i].second ][0] = i;
		}else{
			mi[ ei[i].second ][1] = i;
		}
	}
#ifdef DEBUG_CYCLES
	cout << "mi"<<endl;
	for(unsigned i=0; i<mi.size(); i++){
		cout << i<< " -> "<<mi[i][0]<<" "<<mi[i][1]<<endl;
	}
#endif
	for(int p=0; p<n+1; p++){
		if( !marked[p] ){
			while( !marked[p] ){
				marked[p] = true;
				cycs[p] = c+1;
					// determine next edge
				if( ei[ mi[ p ][0] ] == ei[nei]){
					nei = mi[ p ][1];
				}else{
					nei = mi[ p ][0];
				}
					// determine the next point
				if( p != ei[nei].first){
					p = ei[nei].first;
				}else if( p != ei[nei].second){
					p = ei[nei].second;
				}
			}
			c++;
		}
	}

	return c;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void genom::de_identify_g( const genom &origin ){
	genom tmp( size(), 0 ); 	// temporary genom



	for( unsigned i=0; i<size(); i++ ){
		tmp[ i ] = origin[ abs(chromosom[i])-1 ];
		if( chromosom[i] < 0)
			tmp[ i ] *= -1;
	}

	chromosom = tmp.chromosom;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//int distance(const genom &g1, const genom &g2, int n, hdata &hd){
//	int d = 0,  // the reversal distance
//		c = 0, 	// number of cycles
//		t = 0, 	// the ominous parameter t
//		len = 0;	// variable for the determination of the branch lengths in the pqtree
//	pqnode *pqroot;	// the root of the pq-tree
//	vector<int> pi,
//		comp_or,
//		branch_lens;
//	vector<pair<int,int> > comp;
//	static bool initialised = false;
//	static vector<int>	cycs;	// the cylces of the permutation
//
//	if(!initialised){
//		cycs = vector<int>(n+1, 0);
//		initialised = true;
//	}
//
//	pi = g2.identify(g1);
//#ifdef DEBUG_DISTANCE
//	cout << "pi"<<endl;
//	for(unsigned i=0; i<pi.size(); i++)
//		cout << pi[i]<<" ";
//	cout << endl;
//#endif
//	c = cycles( pi, n, cycs );
//
//	comp = getIrreducibleConservedIntervals(g1, g2, comp_or);
//#ifdef DEBUG_DISTANCE
//	for(unsigned i=0; i<comp.size(); i++)
//		cout <<"C"<< comp[i].first<<","<<comp[i].second<<" "<<comp_or[i]<<endl;
//#endif
//
//	pqtree(comp, comp_or, n, &pqroot);
////	pqtree_print(pqroot);
//	pqtree_defoliate(&pqroot);
//
//	pqtree_branches(pqroot, len, branch_lens);
//
//
//	if(branch_lens.size() > 0){
//		t = branch_lens.size();
//
//		if( branch_lens.size() % 2 != 0 ){	// odd number of leaves
//			for (unsigned i=0; i<branch_lens.size(); i++){	// search one short branch
//				if( branch_lens[i] == 1 ){
//					branch_lens[0] = std::numeric_limits< int >::max();
//					break;
//				}
//			}
//			if(branch_lens[0] != std::numeric_limits< int >::max()){	// if there is no short branch -> fortress
//				t++;
//			}
//		}
//	}
//
//	pqtree_free(pqroot);
//#ifdef DEBUG_DISTANCE
//	cout <<"c : "<<c<<endl;
//	cout <<"t : "<<t<<endl;
//#endif
//	d = n + 1 - c + t;
//	return d;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void double_genome( genom &g, vector<string> *nmap, bool directed ){
	genom lg;

	lg = genom();

	for(unsigned j=0; j<g.size(); j++){
		if( g[j] > 0 ){
			lg.push_back((2 * g[j]) - 1);
			lg.push_back(2 * g[j]);
		}else{
			lg.push_back(-2 * g[j]);
			lg.push_back((-2 * g[j]) - 1);
		}
	}


	if( directed ){
		lg.chromosom.insert(lg.chromosom.begin(), 0);
		lg.push_back( lg.size() );
	}
	lg.set_nmap( nmap );

	g = lg;

//	for(unsigned i=0; i<genomes.size(); i++){
//		genomes[i].set_nmap( NULL );
////		genomes[i].set_nmap( &nmap );
//	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void double_genomes( vector<genom> &genomes, vector<string> *nmap, bool directed ){
	for(unsigned i=0; i<genomes.size(); i++){
		double_genome(genomes[i], nmap, directed);
	}
}

vector<string> double_nmap( const vector<string> &nmap, bool directed ){
	vector<string> lnmap;

	lnmap = vector<string>( 2*nmap.size()-1 );
	for( unsigned i=1; i<nmap.size(); i++ ){
		lnmap[ 2*i-1 ] = "-"+nmap[i];
		lnmap[ 2*i ] = "+"+nmap[i];
	}
	if( directed ){
		lnmap[0] = "STA";
		lnmap.push_back("END");
	}
	return lnmap;
}

///////////////////////////////////////////////////////////////////////////////

void elementary_intervals(const vector<int> &pi, int n, vector<pair<int,int> > &ei){

		// permutation: p =  p_0   p_1   p_2   p_3  ....   p_n
		// points          0     1     2     3     4     n
	static bool initialised = false;
	static vector<int> pi_inv;			// the inverse permutation (big enough to store the frame elements)

	if(!initialised){
		pi_inv = vector<int>(n+2);
		initialised = true;
	}

	for(unsigned i=0; i<pi.size(); i++)	// get the inverse permutation (here it's not necessary to treat the signs)
		pi_inv[abs(pi[i])] = i;				// -> now it's possible to get the index of elements in O(1)

		// for every element i of the permution there is an elementary reversal
		// which start right of i if i is positive and left otherwise
		// 	and ends left of the element (i+1) if it's positive and right of it else
		// see Anne Bergeron: "Reversal distance without hurdles and fortresses"
	for (unsigned i=0; i< pi.size()-1; i++){
		if(pi[pi_inv[i]]>=0)
			ei[i].first = pi_inv[i];
		else
			ei[i].first = pi_inv[i]-1;

		if(pi[pi_inv[i+1]]>=0)
			ei[i].second = pi_inv[i+1]-1;
		else
			ei[i].second = pi_inv[i+1];
	}
		// the result contains all elementary reversals
		// they are at the positions of the elements where they start

}


///////////////////////////////////////////////////////////////////////////////


void get_all_permutations(unsigned n, int circular, int sig, vector<genom> &all){
	genom g;
	int offset = 0;
	vector<vector<unsigned> > signs;
	vector<int> chromosom(n);

	for (unsigned i=0; i<n; i++){
		chromosom[i] = i+1;
	}

	if (circular != 0){
		offset = 1;
	}

	if(sig){
		counter s(n, 2, 0, false);
//		s = vector<int>(n, 0);

		do{
			signs.push_back( s.get_counter() );
			s++;
		}while( s.isvalid() );
//		for(int i=0; i<pow(2, (int) n); i++){
//			signs.push_back(s);
//			s[0]++;
//			for(unsigned j=0; j<s.size(); j++){
//				if( s[j] > 1 ){
//					s[j] = 0;
//					if( j+1 < s.size() ){
//						s[j+1]++;
//					}
//				}else{
//					break;
//				}
//			}
//		}
//		cout << "signs "<<signs.size()<<endl;
	}else{
		signs.push_back( vector<unsigned>(n, 0) );
	}
		// construct all possible permutations
	do{
		for(unsigned i=0; i<signs.size(); i++){
			g.setChromosom(chromosom);
			for(unsigned j=0; j<signs[i].size(); j++){
				if(signs[i][j])
					g[j] *= -1;
			}
			all.push_back(g);
		}
	}while (next_permutation(chromosom.begin()+offset, chromosom.end()));

	return;
}


/*string genom::getName() {
	return name;
}*/

vector<int> genom::getPointOrientations(const genom &g) const{
	vector<int> orientations(size()+1,0);
	vector<int> pi;
		// a point p.q is:
		// 	positiv (+1) if p and q are positive
		// 	negative (-1) if p and q are negative
		// 	0 else
	pi = g.identify(*this );		// make 2 permutions to 1 and id
	//~ pi.insert(pi.begin(),0);		// frame pi by 0 and size+1
	//~ pi.insert(pi.end(),size()+1);

	for(unsigned i=0; i<orientations.size(); i++){
		if(pi[i]>=0 && pi[i+1]>=0){
			orientations[i] = 1;
		}else if(pi[i]<0 && pi[i+1]<0){
			orientations[i] = -1;
		}
	}
	return orientations;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



vector< pair<int, int> > genom::getReversals() const{
	vector< pair<int, int> > reversals;

	// linear case :
	// 	* run through all possiblilities,
	// circular case
	// 	* except all including the first element

	unsigned start = 0;
	if (circular)
		start = 1;

	for (; start < size(); start++) {
		for (unsigned int end = start; end < size(); end++) {
			reversals.push_back(pair<int,int>(start, end));
		}
	}

	return reversals;
}

///////////////////////////////////////////////////////////////////////////////

vector< pair<int, int> > genom::getReversalsSameCycle(const genom &g) const{
	unsigned cycle_cnt;
	static vector<int> cycs;
	vector<int>	pi;
	vector< vector<int> > cycleIndices;
	vector<pair<int,int> > cycleReversals;

	if( cycs.size() != g.size()+1 ){
		cycs = vector<int>(g.size()+1, 0);
	}

	pi = g.identify(*this);			// make 2 permutions to 1 and id
	cycle_cnt = cycles(pi, g.size(), cycs);			// get the cycles

		// in the cycles vector for every point the index of the cycle is stored
		// --
		// check for each point: if there are other points of the same cycle already known
		//      - IF: make pairs of the current point and the already known points of the cycle
		// store the current point on the "cycleIndexArray"

	for(unsigned i=0; i<cycs.size(); i++){		// for every cycle
		if(cycleIndices.size() <= (unsigned) cycs[i]){
			cycleIndices.resize(cycs[i]+1);
		}
			// check if there was earlier points with the same index
		if(cycleIndices[cycs[i]].size()){
				// make pairs with all of them
			for(unsigned j=0; j<cycleIndices[cycs[i]].size(); j++){
				cycleReversals.push_back(pair<int,int>(cycleIndices[cycs[i]][j], i-1));
			}
		}
			// put the current point on the Array of its cycle
		cycleIndices[cycs[i]].push_back(i);
	}
	return cycleReversals;
}

///////////////////////////////////////////////////////////////////////////////


vector< pair<int, int> > genom::getReversalsSameUnorientedComponent(const genom &g) const{
//	unsigned component_cnt,			// the number of components
//		cycle_cnt;
	vector< pair<int, int> > rsuc;			//
//	vector<int> pi,
//		component,					// component index of each point
//		cycs,								// cycle index of each point
//		orientation;						// orientation of each point
//	vector<int> component_orientation;		// orientation of each component
//	vector< vector<int> > component_indices;
//
//	if( cycs.size() != g.size()+1 ){
//		cycs = vector<int>(g.size()+1, 0);
//	}
// @todo: reimplement
	cerr << "getReversalsSameUnorientedComponent: only stub"<<endl;
	exit(EXIT_FAILURE);

//		// get the component indices,  cycle indices and orientation for all points
//	component = getComponents(g, component_orientation, component_cnt);
//
//	pi = g.identify(*this);		// make 2 permutions to 1 and id
//	cycle_cnt = cycles(pi, g.size(),cycs);
//	//~ orientation = getPointOrientations(g);
//
//		//~ // resize the vectors to the number of component std::numeric_limits< int >::max() marks as unwritten
//	//~ component_orientation = vector<int>(number_of_components, std::numeric_limits< int >::max());
//	component_indices = vector< vector<int> >(component_cnt, vector<int>());
//
//		// get the reversals
//	for(unsigned i=0; i<cycs.size(); i++){
//			// if the current point (and thus the current cycle) belongs to an oriented component it isn't interresting
//		if(component_orientation[ component[i] ] == 1)
//			continue;
//			// check if there was earlier points with the component same index
//			// if the points belong to different cycles make pairs with all of them
//		if(component_indices[ component[i]-1 ].size()){
//			for(unsigned j=0; j<component_indices[component[i]-1].size(); j++){
//				if(cycs[i] != cycs[ component_indices[component[i]-1][j] ])
//					rsuc.push_back(pair<int,int>(component_indices[component[i]-1][j], i-1));
//			}
//		}
//			// mark the index as already visited (should be concerned while processing points of the same component)
//		component_indices[component[i]-1].push_back(i);
//	}
//
	return rsuc;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< pair<int, int> > genom::getReversalsSameUnorientedComponentAndSameCycle(const genom &g) const{
	vector< pair<int, int> > add,
		rsucsc;	//
// @todo: reimplement
	cerr << "getReversalsSameUnorientedComponentAndSameCycle: only stub"<<endl;
	exit(EXIT_FAILURE);
//	rsucsc = getReversalsSameCycle(g);
//	add = getReversalsSameUnorientedComponent(g);
//	rsucsc.insert(rsucsc.end(), add.begin(), add.end());
//
	return rsucsc;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


vector< pair<int,int> > genom::getReversals_nonAdjacencyBreaking(const genom &respect) const{
	vector<int> pi;
	vector< pair<int, int> > reversals;

	//~ cout << *this << endl;
	//~ cout << respect << endl;

	pi = respect.identify(*this);			// make 2 permutions to 1 and id
	//~ pi.insert(pi.begin(),0);		// frame pi by 0 and size+1 to get the points at front and
	//~ pi.insert(pi.end(),size()+1);	// 	back of the permutation

	//~ cout<<"PI "<<endl;
	//~ for(unsigned i=0; i<pi.size(); i++){
			//~ cout << pi[i]<<" ";
	//~ }
	//~ cout << endl;

	for(unsigned i=0; i<pi.size()-1; i++){	// go through all possible startpoints
		//~ cout << "i"<<i<<" : " << pi[i+1] <<"-"<<pi[i] <<endl;
		if(pi[i+1] - pi[i] == 1)		// Adjacency -> dont start reversal here
			continue;						// try the next
		//~ cout << "start "<<i<< endl;
		for(unsigned j=i+1; j<pi.size()-1;j++){ // go through all possible endpoints
			//~ cout << "j"<<i<<" : " << pi[j+1] <<"-"<<pi[j] <<endl;
			if(pi[j+1] - pi[j] == 1)	// Adjacency -> dont start reversal here
				continue;					// try the next
			//~ cout << "("<<i<<","<<j-1<<")"<<endl;
			reversals.push_back(pair<int,int>(i,j-1) );
		}
	}
	return reversals;
}


///////////////////////////////////////////////////////////////////////////////


vector<vector<unsigned> > genom::getTranspositions(){
	/*!@todo implement circular */

	vector<unsigned> t(3);
	vector<vector<unsigned> > transpositions;

	for(unsigned i=0; i<chromosom.size()-1; i++){
		t[0]=i;
		for(unsigned j=i+1; j<chromosom.size(); j++){
			t[1]=j;
			for(unsigned k=j+1; k<chromosom.size()+1; k++){
				t[2]=k;
				transpositions.push_back(t);
			}
		}
	}

	return transpositions;
}

///////////////////////////////////////////////////////////////////////////////

genom genom::identify_g(const genom &q) const{
	unsigned s = size();
	static genom g;
	static vector<int> pi_inv;

	if(g.size() != s){
		g.chromosom.resize(s);
		pi_inv.resize(s);
	}

		// get the inverse permutation of p [this] -> p^-1
		// this is the permutation where the i-th element
		// is the position of the i-th element in p
		// (+ possible sign swaps)
	for(unsigned i=0; i<s; i++){
		pi_inv[ abs(chromosom[i]) - 1 ] = i+1;
		if(chromosom[i]<0)
			pi_inv[abs(chromosom[i]) - 1] *= -1;
	}

		// compute (p^-1) * q
		// this is the permutation where the i-th element is the
		// element of p^-1 at the position of the i-th element of q
	for(unsigned i=0; i<s;i++){
		if(q[i]<0)
			g[i] = -1 * pi_inv[abs(q[i]) - 1];
		else
			g[i] = pi_inv[abs(q[i]) - 1];
	}

	return g;
}

///////////////////////////////////////////////////////////////////////////////

vector<int> genom::identify( const genom &q ) const{
	//!@todo ?? is it enough to begin with 1 for circular genomes because they always begin with 1

	unsigned s = size()+1;
	vector<int> pi_inv( q.size() ),
		identified_chr( q.size() + 1);

	// get the inverse permutation of p [this] -> p^-1
	// this is the permutation where the i-th element
	// is the position of the i-th element in p
	// (+ possible sign swaps)
	for(unsigned i=1; i<s; i++){
		pi_inv[abs(chromosom[i-1])-1] = i;
		if(chromosom[i-1]<0)
			pi_inv[abs(chromosom[i-1])-1] *= -1;
	}

	// compute (p^-1) * q
	// this is the permutation where the i-th element is the
	// element of p^-1 at the position of the i-th element of q
 	for(unsigned i=1; i<s;i++){
 		identified_chr[i] = pi_inv[abs(q.chromosom[i-1])-1];
 		if(q.chromosom[i-1]<0)
 			identified_chr[i] *= -1;
	}

	return identified_chr;

}

///////////////////////////////////////////////////////////////////////////////

void identify(vector<genom> &genomes){
	unsigned s = genomes[0].size();
	vector<int> pi_inv( s, 0);
		// get the inverse permutation of p [this] -> p^-1
		// this is the permutation where the i-th element
		// is the position of the i-th element in p
		// (+ possible sign swaps)
	for(unsigned i=0; i<s; i++){
		pi_inv[ abs(genomes[0][i]) - 1] = i+1;
		if(genomes[0][i]<0)
			pi_inv[abs(genomes[0][i]) - 1] *= -1;
	}

		// compute (p^-1) * q
		// this is the permutation where the i-th element is the
		// element of p^-1 at the position of the i-th element of q
	for(unsigned g=0; g<genomes.size(); g++){
		for(unsigned i=0; i<s;i++){
			//~ genomes[g][i] = d.pi_inv[abs(genomes[g][i])];
			if(genomes[g][i]<0)
				genomes[g][i] = -1 * pi_inv[abs(genomes[g][i]) - 1];
			else
				genomes[g][i] = pi_inv[abs(genomes[g][i]) - 1];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void insertdelete(genom &g, const vector<set<int> > &id, bool del, bool inv){


	if( id.size() != 3 ){
		cerr << "error: indel called with "<<id.size()<<" set(s)"<<endl;
		exit(EXIT_FAILURE);
	}

	vector<pair<int,int> > ints(3);
	int s;

	if( del ){
		ints[1] = g.interval( id[1] );
//		cout << "int1 "<<ints[1].first<<" "<<ints[1].second<<endl;
		g.chromosom.erase( g.chromosom.begin()+ints[1].first, g.chromosom.begin()+ints[1].second+1 );
	}else{
		set<int> ins;
		if( inv ){
			for( set<int>::const_iterator it=id[1].begin(); it!=id[1].end(); it++ ){
				ins.insert( -1*(*it) );
			}
		}else{
			ins = id[1];
		}


		if( id[0].size() != 0 && id[2].size() != 0 ){
			ints[0] = g.interval( id[0] );
			ints[2] = g.interval( id[2] );
//			cout << "int0 "<<ints[0].first<<" "<<ints[0].second<<endl;
//			cout << "int2 "<<ints[2].first<<" "<<ints[2].second<<endl;
			ints[0].first = max(ints[0].first, ints[0].second );
			ints[2].first = max(ints[2].first, ints[2].second );
			s = min( ints[0].first, ints[2].first );
//			cerr << "s "<<s<<endl;
			g.chromosom.insert( g.chromosom.begin()+s+1, ins.begin(), ins.end() );
		}else if( id[0].size() == 0 && id[2].size() != 0 ){
			g.chromosom.insert( g.chromosom.begin(), ins.begin(), ins.end() );
		}else if( id[2].size() == 0 && id[0].size() != 0 ){
			g.chromosom.insert( g.chromosom.end(), ins.begin(), ins.end() );
		}else {
			cerr<<"error: insertdelete with empty frame"<<endl;
			exit(EXIT_FAILURE);
		}
//		cerr << "iserted "<<g<<endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

pair<int,int> genom::interval(const set<int> &elements) const{
	pair<int,int> i(std::numeric_limits< int >::max(), std::numeric_limits< int >::min());
	vector<int> inv;		// the inverse permutation

//	cerr << "find interval "; copy(elements.begin(), elements.end(), ostream_iterator<int>(cerr, " ")); cerr << endl;
//	cerr <<"in "<<*this<<endl;

		// get the inverse permutation
	inv = inverse( );
//	copy(inv.begin(), inv.end(), ostream_iterator<int>(cerr, " ")); cerr << endl;

	for(set<int>::iterator it=elements.begin(); it!=elements.end(); it++){
//		cerr  << "element "<<*it<<" pos "<< abs( inv[*it] )<<endl;
		i.first  = min( i.first , abs( inv[*it] ) );
		i.second = max( i.second, abs( inv[*it] ) );
	}

	if( ((int)elements.size()-1) != (i.second - i.first) ){
//		cerr << "genom::interval error: invalid interval ("<< i.first<<","<<i.second <<") ";
//		for( int j=i.first; j<=i.second; j++){
//			cerr << chromosom[j]<<" ";
//		}cerr << endl;
//		cerr << "for set: ";
//		copy(elements.begin(), elements.end(), ostream_iterator<int>(cerr, " ")); cerr <<endl;
//		cerr << "                       genom: "<<*this<<endl;
		throw NoIntervalException( *this, elements );
//		exit(1);
	}

	return i;
}

///////////////////////////////////////////////////////////////////////////////

vector<int> genom::inverse( bool sign ) const{

//	cerr << "genom::inverse() "<< *this << endl;
//	cerr << "size() "<<size()<<endl;

	unsigned s=0;
	for(unsigned i=0; i<chromosom.size(); i++){
		s = max( s, (unsigned)abs(chromosom[i]) );
	}

	vector<int> inv(s+1, std::numeric_limits< int >::max());	// the inverse permutation

	for(unsigned i=0; i<chromosom.size(); i++){
		//cerr <<"i "<<i<<" " <<abs(chromosom[i]) <<endl;
		inv[ abs(chromosom[i]) ] = i;
		if( sign && chromosom[i] < 0){
			inv[ abs(chromosom[i]) ] *= -1;
		}
	}
	return inv;
}

///////////////////////////////////////////////////////////////////////////////

bool genom::isneg() const{
	for(unsigned i=0; i<chromosom.size(); i++){
		if( chromosom[i] > 0 ){
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

bool genom::ispos() const{
	for(unsigned i=0; i<chromosom.size(); i++){
		if( chromosom[i] < 0 ){
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

void genom::normalize(int to){
	//!@todo distinction between undirected & directed linear chromosomes. til now all linear are directed
	//!@todo seems to be called a little bit too often
	unsigned onePos = 0;

		// search element to
	for(unsigned i=0; i<chromosom.size(); i++){
		if(abs(chromosom[i]) == abs(to)){
			onePos = i;
			break;
		}
	}

		// if its in the wrong direction -> rotate
	if(chromosom[onePos] == -1 * to){
		std::reverse(chromosom.begin(), chromosom.end());
		for(unsigned i=0; i<chromosom.size(); i++)
			chromosom[i] *= -1;
	}

		// bring the element to the front
	rotate(chromosom.begin(), find(chromosom.begin(), chromosom.end(), to), chromosom.end());
}

///////////////////////////////////////////////////////////////////////////////

int& genom::operator[](unsigned i){
	return chromosom[i];
}

///////////////////////////////////////////////////////////////////////////////

int genom::operator[](unsigned i) const{
	return chromosom[i];
}

///////////////////////////////////////////////////////////////////////////////

bool genom::operator==(const genom& g) const {
	return equal(chromosom.begin(), chromosom.end(), g.chromosom.begin());
}

///////////////////////////////////////////////////////////////////////////////

bool genom::operator!=(const genom& g) const {
	if (*this == g)
		return false;
	else
		return true;
}

///////////////////////////////////////////////////////////////////////////////

bool genom::operator<(const genom& g) const {
//	cerr << *this << " < "<<g<<endl;

	unsigned s = min(this->size(), g.size());
//cerr << s<<endl;
	for (unsigned i = 0; i < s; i++) {
//		cerr << chromosom[i]<<" < "<<g.chromosom[i]<<endl;

		if ( chromosom[i] > g.chromosom[i]) {
			return false;
		} else if (chromosom[i] < g.chromosom[i]){
			return true;
		}
	}

	// -> up to s both are equal
	// -> decide by length

	if( chromosom.size() < g.size() )
		return true;
	else
		return false;
}


///////////////////////////////////////////////////////////////////////////////

ostream &operator<<(ostream &out, const genom &g) {
	for (unsigned i = 0; i < g.size(); i++) {
		print_element( g.chromosom[i], out, 1, "", g.nmap );
		if( i < g.size() -1  )
			out << " ";
	}
	return out;
}

///////////////////////////////////////////////////////////////////////////////

ostream & operator<<(ostream &out, const vector<genom> &genomes){
	for (unsigned i=0; i<genomes.size(); i++){
		out << genomes[i]<<endl;
	}

	return out;

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void print_element( int e, ostream &o, vector<string> namemap, string plus ){
//	if( namemap.size() == 0 ){
//		o.width(3);
//		o << e;
//	}else{
//		if( e < 0 ){
//			o << "-";
//		}else{
//			o<<plus;
//		}
//		o << namemap[ abs(e) ];
//	}
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void genom::print_element( unsigned idx, ostream &out, int inv, string plus ) const{
//	print_element( chromosom[idx], out, inv, plus, nmap );
//}

void print_element( int e, ostream &out, int inv, string plus, const vector<string> *nmap ){

	if( nmap == NULL ){
		out.width(3);
		out << e;
	}else{
//		cerr << "print "<<e<<" ";
	//	cerr << (*nmap)[ abs(e) ]<<endl;

//		if(nmap->size()<=(unsigned)abs(e)){
//			cerr << "error in print_element "<<e<<" not in namemap of size "<<nmap->size() <<endl;
//			for(unsigned i=0;i<nmap->size(); i++)
//				cerr << (*nmap)[i]<<" ";
//			cerr << endl;
//			exit(0);
//		}


		if( inv*e < 0 ){
			out << "-";
		}else{
			out<<plus;
		}
		out << (*nmap)[ abs(e) ];
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void print_genom(ostream &out, const genom &g, const vector<string> &nmap) {
//	for (unsigned i = 0; i < g.size(); i++) {
//		print_element( g[i], out, nmap);
//		if(i < g.size()-1)
//			out << " ";
//	}
//	out.flush();
//}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void genom::push_back(int gene) {
	chromosom.push_back(gene);
}


///////////////////////////////////////////////////////////////////////////////

void genom::randomise(){
	random_shuffle( chromosom.begin(), chromosom.end() );
	for( unsigned i=0; i<chromosom.size(); i++ ){
		if(ask_rng() % 2 != 0){
			chromosom[i] *= -1;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////

//pair<unsigned, unsigned> genom::rank(hdata &d){
//	pair<unsigned, unsigned> rank;
//	rank.second = 0;
//
//		// for every element that is greater than 0 push it
//		// on the unsigned vector ; for smaller values : push
//		// the abs value on the unsigned vector
//		// additional decrease by 1 to get a permutation from
//		// 0..n-1 (for my convenience)
//	for (unsigned i=0; i< chromosom.size(); i++){
//		if(chromosom[i]>0){
//			d.u_pi[i]= chromosom[i]-1;
//		}else{
//			d.u_pi[i] = (-1*chromosom[i])-1;
//			rank.second += pow(2,i);
//		}
//	}
//
//		// get the inverse permutation of the unsigned permutation
//	//~ d.inv_u_pi.resize(u_pi.size());
//	for(unsigned i=0; i<d.u_pi.size(); i++){
//		d.inv_u_pi[ d.u_pi[ i ] ] = i;
//	}
//
//		// compute the rank
//	rank.first = rank_perm(d.u_pi.size(), d.u_pi, d.inv_u_pi);
//
//	return rank;
//}

///////////////////////////////////////////////////////////////////////////////

void reverse(genom &g, int start, int end) {
		// test if start > end -> swap
	if (start > end) {
		swap(start, end);
	}

		// reverse the interval
	while ((end - start) > 0) {
		swap (g.chromosom[start], g.chromosom[end]);
		g.chromosom[start] *= -1;
		g.chromosom[end] *= -1;

		start++;
		end--;
	}
	if (start == end){ // odd size
		g.chromosom[start] *= -1;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void reverse(genom &g, pair<int, int> r) {
	reverse(g, r.first , r.second);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void reverse(genom &g, const set<int> &elements){
	reverse( g, g.interval(elements) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void reverse(vector<genom> &genomes, int s, int e){
	for(unsigned i=0; i<genomes.size(); i++){
		reverse(genomes[i], s, e );
	}
}

///////////////////////////////////////////////////////////////////////////////

void reverse(genom g, vector< pair<int, int> > &r, map<genom, bool> &used, int dist, unsigned &out){
	// vector<genom> result;
	vector<pair<int,int> > todo;

	/*cout << "d"<< dist<<endl;
	if (dist == 0) {
		return;
	}
	cout << dist << endl;*/
		// do all remaining reversals
	for (unsigned i=0; i<r.size(); i++){

		genom tmp = g;				// reverse the genome (a temporary copy)
		reverse(tmp, r[i]);
		if(used[tmp]){				// test if the genome was already generated
			//cout << "USED"<<endl;
			continue;				// if then continue with the next iteration
		}


		used[tmp] = 1;			// else mark the genome
		//result.push_back(g);	// an store it in the result vector

			// check if other reversals got disturbed
		for(unsigned j=0; j<r.size(); j++){
			if(i!=j){
				pair<int,int> temp(r[j]);

					// test if r[j] don't overlaps r[i] from right and test if r[j] don't overlaps r[i] from left
				if(!(r[i].first < temp.first && temp.first <= r[i].second && r[i].second < temp.second) &&
					!(r[i].second > temp.second && temp.second >= r[i].first && r[i].first > temp.first)){
						// if the intevals are nested the indices need a update
					if(r[i].first <= temp.first && temp.second <= r[i].second){
						temp.first = r[i].first+(r[i].second - temp.first);
						temp.second = r[i].first+(r[i].second - temp.second);
						if (temp.first > temp.second){	// let the indices be in the right order
							int t = temp.first;
							temp.first = temp.second;
							temp.second = t;
						}
					}
					if(temp.first < 0 || temp.second < 0){
						cout << "SMALLER "<<endl;

						cout << "ri "<< r[i].first << ","<<r[i].second << endl;
						cout << "rj "<< r[j].first << ","<<r[j].second << endl;

					}

					todo.push_back(temp);		// store the reversal on a todo vector
				}else{
					out++;
				}
			}
		}
			// if there is work left
		if(todo.size()){
			//vector<genom> gg =
			reverse(g, todo, used, dist-1, out);				// do it
			//result.insert(result.end(), gg.begin(), gg.end());
			todo.clear();
		}
	}
	//cout << result.size()<<endl;

	return;
}


///////////////////////////////////////////////////////////////////////////////

vector<genom> genom::reverseHamiltonian(vector<pair<int,int> > &r){
	vector<unsigned> horder;
	vector<genom> hgenomes;
	genom temp = *this;

	//~ cout << "start H"<<endl;
	//~ cout << *this<<endl;
	//~ cout << temp<<endl;

	horder.push_back(0);
	//~ cout << "H 1"<<endl;
	for(unsigned i=1; i<r.size() && i<8; i++){
		//~ horder.push_back(i);
		//~ cout << horder.size()<<" new size " <<horder.size()*2-1<<endl;
		//~ back_insert_iterator<vector<unsigned> > ii(horder);

		vector<unsigned>::iterator start = horder.begin(),
			end=horder.end();
		horder.resize(horder.size()*2);
		copy(start, end,  end);
		//~ horder.pop_back();

		cout << i << ": ";
		for(unsigned j=0; j<horder.size(); j++){
			cout << horder[j]<<" " ;
		}
		cout << endl;

		//~ for(int j=horder.size()-2; j>=0; j--){
			//~ horder[horder.size()-j] = horder[j];
		//~ }

		cout << horder.size()<<endl;;
	}
	cout << "H 2"<<endl;
	for(unsigned i=0; i<horder.size(); i++){
		reverse(temp, r[horder[i]]);
		hgenomes.push_back(temp);
	}

	cout << "end H"<<endl;
	return hgenomes;
}

///////////////////////////////////////////////////////////////////////////////

/*vector<genom> genom::reverse(vector< pair<int, int> > &r, map<genom, bool> &used, int dist){
	vector<genom> result;
	vector<pair<int,int> > todo;
	//cout << "d"<< dist<<endl;
	if (dist == 0)
		return result;

		// do all remaining reversals
	for (unsigned i=0; i<r.size(); i++){

		genom g = *this;		// reverse the genome (a temporary copy)
		g.reverse(r[i]);
		if(used[g]){				// test if the genome was already generated
			//cout << "USED"<<endl;
			continue;			// if then continue with the next iteration
		}
		used[g] = 1;			// else mark the genome
		result.push_back(g);	// an store it in the result vector

			// check if other reversals got disturbed
		for(unsigned j=0; j<r.size(); j++){
			if(i!=j){
				pair<int,int> temp(r[j]);

					// test if r[j] don't overlaps r[i] from right and test if r[j] don't overlaps r[i] from left
				if(!(r[i].first < temp.first && temp.first <= r[i].second && r[i].second < temp.second) &&
					!(r[i].second > temp.second && temp.second >= r[i].first && r[i].first > temp.first)){
						// if the intevals are nested the indices need a update
					if(r[i].first <= temp.first && temp.second <= r[i].second){
						temp.first = r[i].first+(r[i].second - temp.first);
						temp.second = r[i].first+(r[i].second - temp.second);
						if (temp.first > temp.second){	// let the indices be in the right order
							int t = temp.first;
							temp.first = temp.second;
							temp.second = t;
						}
					}
					if(temp.first < 0 || temp.second < 0){
						cout << "SMALLER "<<endl;

						cout << "ri "<< r[i].first << ","<<r[i].second << endl;
						cout << "rj "<< r[j].first << ","<<r[j].second << endl;

					}

					todo.push_back(temp);		// store the reversal on a todo vector
				}
			}
		}
			// if there is work left
		if(todo.size()){
			vector<genom> gg = g.reverse(todo, used, dist-1);				// do it
			if(dist == 1000){
				cout << "+"<<gg.size()<<endl;
			}

			result.insert(result.end(), gg.begin(), gg.end());
			todo.clear();
		}
	}
	//cout << result.size()<<endl;

	return result;
}*/

///////////////////////////////////////////////////////////////////////////////

void genom::setChromosom(vector<int> c) {
	chromosom = c;
}

///////////////////////////////////////////////////////////////////////////////

void genom::setCircular(char c) {
	circular = c;
	if(circular)
		normalize();
}

///////////////////////////////////////////////////////////////////////////////
void genom::set_nmap( vector<string> *nm ){
	nmap = nm;
}

///////////////////////////////////////////////////////////////////////////////

//~ void genom::setName(string n) {
	//~ name = n;
//~ }

///////////////////////////////////////////////////////////////////////////////

unsigned int genom::size() const{
	return chromosom.size();
}


///////////////////////////////////////////////////////////////////////////////

void transpose(genom &g, unsigned i, unsigned j, unsigned k){
	vector<unsigned> t;
	t.push_back(i);
	t.push_back(j);
	t.push_back(k);
	transpose(g, t);
}

///////////////////////////////////////////////////////////////////////////////

void transpose(genom &g, vector<unsigned> &t){
	sort(t.begin(), t.end());
	g.chromosom.insert(g.chromosom.begin()+t[2], g.chromosom.begin()+t[0], g.chromosom.begin()+t[1]);
	g.chromosom.erase(g.chromosom.begin()+t[0], g.chromosom.begin()+t[1]);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void transpose(genom &g, const set<set<int> > &t){
	vector<pair<int,int> > tp;
	if( t.size() != 2 ){
		cerr << "error: transpose called with "<<t.size()<<" set(s)"<<endl;
		exit(1);
	}

	for(set<set<int> >::iterator it=t.begin(); it!=t.end(); it++){
		tp.push_back( g.interval(*it) );
	}

	sort(tp.begin(), tp.end());
//	cout << "transpose "<< tp[0].first << ","<<tp[0].second<< " - "<<tp[1].first << ","<<tp[1].second<<endl;
//	cout << "       in "<<g<<endl;

	if( (tp[0].second + 1) != tp[1].first ){
		throw NonConsecutiveIntervalsException(g, t);
	}
	transpose(g, tp[0].first, tp[1].first, tp[1].second+1);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void transpose(genom &g, const vector<set<int> > &transp){
	vector<pair<int,int> > tp;
	if( transp.size() != 2 ){
		cerr << "error: transpose called with "<<transp.size()<<" set(s)"<<endl;
		exit(1);
	}

	for(unsigned i=0; i<transp.size(); i++){
		tp.push_back( g.interval(transp[i]) );
	}

	sort(tp.begin(), tp.end());
//	cout << "transpose "<< tp[0].first << ","<<tp[0].second<< " - "<<tp[1].first << ","<<tp[1].second<<endl;
//	cout << "       in "<<g<<endl;
	if( (tp[0].second + 1) != tp[1].first ){
		cerr << "error: transpose non tandem intervals"<<endl;
		for(unsigned i=0; i<transp.size(); i++){
			cout << "{"; copy(transp[i].begin(), transp[i].end(), ostream_iterator<int>(cout, " "));cout << "}"<<endl;
		}
		cout << "genom "<<g<<endl;
		exit(1);
	}
	transpose(g, tp[0].first, tp[1].first, tp[1].second+1);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unique(vector<genom> &genomes){
	vector<genom>::iterator it;

	sort(genomes.begin(), genomes.end());
	it = unique(genomes.begin(), genomes.end());
	genomes.erase(it, genomes.end());
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unique( vector<genom> &genomes, vector<string> &names,
		vector<vector<string> > &tax, vector<unsigned> &cnt, bool report ){

	vector<bool> remove;
	vector<genom> tmpg;
	vector<string> tmps;
	vector<vector<string > > tmpt;

	cnt.clear();

	for(unsigned i=0; i<genomes.size(); i++){
		tmpg.push_back( genomes[i] );
		reverse( tmpg[i], 0, genomes[i].size()-1);
	}

	remove = vector<bool> (genomes.size(), false);
	for(unsigned i=0; i<genomes.size(); i++){
		if(remove[i])
			continue;

		cnt.push_back( 1 );

		for(unsigned j=i+1; j<genomes.size(); j++){
			if(genomes[i] == genomes[j] || genomes[i] == tmpg[j]){
				if( report ){
					cerr << names[i] << "==" <<names[j]<<endl;
				}else{
					names[i]+=","+names[j];
				}
				if( tax.size() > 0 ){
//						copy(tax[i].begin(), tax[i].end(), ostream_iterator<string>(cerr, " ")); cerr << endl;
//						copy(tax[j].begin(), tax[j].end(), ostream_iterator<string>(cerr, " ")); cerr << endl;
					pair<vector<string>::iterator, vector<string>::iterator> mm;
					mm = mismatch(tax[i].begin(), tax[i].begin()+min( tax[i].size(), tax[j].size() ), tax[j].begin() );

					tax[i].erase(mm.first, tax[i].end());
					tax[j].erase(mm.second, tax[j].end());
//						copy(tax[i].begin(), tax[i].end(), ostream_iterator<string>(cerr, " ")); cerr << endl;
//						copy(tax[j].begin(), tax[j].end(), ostream_iterator<string>(cerr, " ")); cerr << endl;
//						cerr << "---------"<<endl;
				}

				cnt[ cnt.size()-1 ]++;
				remove[j] = true;
			}
		}
	}

	tmpg.clear();
	for(unsigned i=0; i<genomes.size(); i++){
		if(remove[i])
			continue;
		tmpg.push_back(genomes[i]);
		tmps.push_back(names[i]);
		if( tax.size() > 0 )
			tmpt.push_back( tax[i] );
	}
	genomes = tmpg;
	names = tmps;
	tax = tmpt;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

genom undouble_genom( const genom &g, vector<string> *snmap, bool directed ){
	genom sg;	// the undoubled genom to be constructed

	sg = genom();
	for( unsigned j = ( directed ) ? 1 : 0;\
		j< ((directed) ? g.size()-1 : g.size()); j+=2 ){
		if( abs(g[j] - g[j+1]) != 1 ){
			cerr << "error: undouble_genomes invalid pair "<<g[j] << " "<<g[j+1]<<endl;
		}

		if( g[j] < g[j+1] ){
			sg.push_back( g[j+1]/2 );
		}else{
			sg.push_back( -1 * g[j]/2 );
		}
	}

#ifdef DEBUG_NONAMES
	sg.set_nmap( NULL );
#else
	sg.set_nmap( snmap );
#endif


	return sg;
}

void undouble_genomes( vector<genom> &genomes, vector<string> *snmap,
		bool directed){

	for( unsigned i=0; i<genomes.size(); i++ ){
		genomes[i] = undouble_genom( genomes[i], snmap, directed );
	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void unsign_genomes( vector<genom> &genomes, vector<string> *nmap, bool directed  ){
	for( unsigned i=0; i<genomes.size(); i++ ){
		for( unsigned j=0; j<genomes[i].size(); j++ ){
			genomes[i][j] = abs( genomes[i][j] );
		}
	}
	if( directed ){
		for( unsigned i=0; i<genomes.size(); i++ ){
			genomes[i].chromosom.insert( genomes[i].chromosom.begin(), 0 );
			genomes[i].chromosom.push_back( genomes[i].size() );

			if( nmap != NULL ){
				(*nmap)[0]="STA";
				nmap->push_back("END");
			}
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned unsign_other( unsigned e ){
	if ( e==0 ){
		return std::numeric_limits< unsigned>::max();
	}else{
		return ( ( e % 2 == 0 ) ? e-1 : e+1 );
	}
}

//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//void condense (vector<genom> &genomes, char circ, int *condense_succ, int *condense_decode, int sign ){
//	struct genome_struct *genome_list;
//	int n;		// number of genes
//
//	if(genomes.size()<3){
//		cerr << "ERROR: condense .. need at least 3 genomes"<<endl;
//		exit(1);
//	}
//
//	n = genomes[0].size();
//
//	genome_list = ( struct genome_struct * ) malloc ( genomes.size() * sizeof ( struct genome_struct ) );
//    if ( genome_list == ( struct genome_struct * ) NULL )
//        cerr << "ERROR: genome_list NULL"<<endl;
//
//	for(unsigned i=0; i<genomes.size(); i++){
//		genome_list[i].genes = (int *) malloc(genomes[i].size() * sizeof(int) );
//		for(unsigned j=0; j<genomes[i].size(); j++){
//			genome_list[i].genes[j] = genomes[i][j];
//		}
//	}
//
//	condense_genes ( genome_list, genomes.size(), &n, circ,
//			condense_succ + n,
//			condense_decode + n);
//
//	for(unsigned i=0; i<genomes.size(); i++){
//		genomes[i] = genom(genome_list[i].genes, n, circ);
//		free(genome_list[i].genes);
//	}
//
//	free(genome_list);
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void uncondense ( vector<genom> &genomes, char circ, int *condense_succ, int *condense_decode, int orig_num_genes ){
//    int j,
//		*succ,
//		*decode,
//		curr;
//
//	vector<genom> uncond_g;
//    succ = condense_succ + orig_num_genes;
//    decode = condense_decode + orig_num_genes;
//
//	for(unsigned g=0; g<genomes.size(); g++){
//		uncond_g.push_back( genom(orig_num_genes, circ) );
//		curr = 0;
//		for (unsigned i = 0; i < genomes[g].size(); i++ ){
//
//			j = decode[ genomes[g][i] ];
//
//
//			while ( j != BREAK )
//			{
//				//~ cout << " "<<j;
//				uncond_g[g][curr] = j;
//				curr ++;
//				j = succ[j];
//			}
//
//		}
//		//~ cout << endl;
//	}
//
////~ cout <<"uncondensed "<<endl<< uncond_g<<endl;
//	genomes = uncond_g;
//	//~ uncond_g.clear();
//    return;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pair<unsigned, unsigned> matching_edge(int pii, int piip1, int len){
	pair<unsigned, unsigned> edge;
	//~ cerr <<"matching_edge("<< pii <<","<<piip1<<","<< len << ") -> ";

	if(abs(pii) < len){
		if( pii >= 0 ){
			edge.first = 2*pii;
		}else{
			edge.first = 2 * (-1*pii)-1;
		}
	}else{
		edge.first = 2*(len-1)+1;
	}

	if(abs(piip1) <= len){
		if( piip1 > 0 ){
			edge.second = 2*piip1 - 1;
		}else if(piip1 < 0){
			edge.second = 2* (-1*piip1);
		}else{
			edge.second = 0;
		}
	}else{
		edge.second = 2*len+1;
	}
	//~ cerr << edge.first<<","<<edge.second <<endl;
	return edge;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//bool is_initialised(hdata &hd){
//	return hd.initialised;
//}
//
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//void init_data(hdata &d, unsigned len, char circ){
//		// init memory for identify and identify_g
//	d.g = genom(len, circ);
//	d.pi_inv.resize(len, 0);
//	d.identified_chr.resize(len+2, 0);
//	d.identified_chr[0] = 0;
//	d.identified_chr[len+1] = len+1;
//		// init memory for rank
//	d.u_pi.resize(len, 0);
//	d.inv_u_pi.resize(len, 0);
//
//	/* @todo remove */
//	// init memory for distance
//	d.distmem = new_distmem ( len );
//		// init memory for wrapper to capraras median solver
//	d.genes = (int *) malloc(len * sizeof(int));
//	for(unsigned i=0; i<3; i++){
//		d.gs[i] = (genome_struct *) malloc(sizeof(genome_struct));
//	}
//
//	d.initialised = true;
//}
//
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//extern "C" void init_data_c(hdata &d, unsigned len, char circ){
//	init_data( d, len, circ);
//}
//
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//void free_data(hdata &d){
//		// clear memory for identify and identify_g
//	d.g.clear();
//	d.pi_inv.clear();
//	d.identified_chr.clear();
//		// free memory for rank
//	d.u_pi.clear();
//	d.inv_u_pi.clear();
///* @todo remove */
//	// free memory for distance
//	free_distmem ( d.distmem );
//		// free memory for wrapper to capraras median solver
//	free(d.genes);
//	for(unsigned i=0; i<3; i++){
//		free(d.gs[i]);
//	}
//
//	d.initialised = false;
//}
//
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//extern "C" void free_data_c(hdata &d){
//	free_data(d);
//}

