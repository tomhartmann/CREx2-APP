#include <algorithm>
#include <iostream>
#include <limits>
#include <list>
#include <stack>

#include "common.hpp"
#include "helpers.hpp"

//#define DEBUG_SORTINGREV

//#define DEBUG_COMMUTINGINT
//#define DEBUG_GENERATOR
//#define DEBUG_PERFECTD
//#define DEBUG_PRIMED
//#define DEBUG_SORTINTERVALS

using namespace std;

void canonical_generator(vector<int> &r, vector<int> &l){
	int n = r.size()-1;
	static vector<int> rsup, lsup,
		rc, lc;	// canonical generators;

	if((int)rsup.size() != n+1){
		rsup = vector<int>(n+1, 0);
		lsup = vector<int>(n+1, 0);
		rc = vector<int>(r.size(), 0);
		lc = vector<int>(l.size(), 0);
	}

	support(r, l, rsup, lsup);

	rc[1] = n;
	for(int k=2; k<=n; k++){
		rc[k] = k;
	}
	for(int k=n; k>=2; k--){
		if( l[rc[k]] <= rsup[k] && rsup[k] <= rc[k] && rc[k] <= r[rsup[k]]){
			rc[ rsup[k] ] = max(rc[k], rc[ rsup[k] ]);
		}
	}

	lc[n] = 1;
	for(int k=1; k<n; k++){
		lc[k] = k;
	}
	for(int k=1; k<n; k++){
		if( l[ lsup[k] ] <= lc[k] && lc[k] <= lsup[k] && lsup[k] <= r[lc[k]]){
			lc[ lsup[k] ] = min(lc[k], lc[ lsup[k] ]);
		}
	}

	r = rc;
	l = lc;

#ifdef DEBUG_GENERATOR
	cout << "canonical generator "<<endl;
	for(int i=1; i<=n; i++)
		cout << "("<< rc[i]<<","<<lc[i]<<")";
	cout <<endl;
#endif
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool cmp_interval_first( pair<int, int> a, pair<int, int> b ){
	return a.first < b.first;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool cmp_interval_second( pair<int, int> a, pair<int, int> b ){
	return a.second < b.second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void combine_generator(vector<int> &r, vector<int> &l, const vector<int> &radd, const vector<int> &ladd){
	for(unsigned j=0; j<r.size(); j++){
		r[j] = min(r[j], radd[j]);
		l[j] = max(l[j], ladd[j]);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void combined_generator(const vector<genom> &genomes, int n, vector<int> &r, vector<int> &l){
	static vector<int> 	radd,		// aditional generator
		ladd;

	if((int)radd.size() != n+1){
		radd = vector<int>( n+1, 0 );
		ladd = vector<int>( n+1, 0 );
	}

	generator(genomes[1], r, l);

	for(unsigned i=2; i<genomes.size(); i++){
		generator(genomes[i], radd, ladd);
		combine_generator(r, l, radd, ladd);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void common_intervals(const vector<int> &r, const vector<int> &l, int circular, int trivial, int sign, vector<pair<int, int> > &comint){
	int i=0,
		n = r.size()-1,
		x = 1;
	static vector<int> rsup,
		lsup;

	if( (int) rsup.size() != n+1 ){
		rsup = vector<int>(n+1, 0);
		lsup = vector<int>(n+1, 0);
	}
	comint.clear();

	if(circular)
		x = 2;

	support(r, l, rsup, lsup);

	for(int j=n; j>0; j--){
		i=j;
		while( i >= l[j] ){
//			cerr << "X"<<i-1<<" "<<j-1<<endl;
				// if the trivial intervals are wanted -> output them
				// otherwise exclude:
				// - intervals (i..j) with i==j
				// - the interval (1..n), i.e. length = n
				// - in the circular case: the intervals with length n-1, i.e. the complementary
				//	intervals of (i..j) with i==j
			if( trivial || ( i<j && j-i < n-x) ){
				comint.push_back( make_pair(i-1,j-1) );
//				cout << i-1 << ":"<<j-1<<endl;
				if(circular){
//					if( j%n > (n+i-2)%n ){
						comint.push_back( make_pair( j%n, (n+i-2)%n) );
//						cout << j%n << ";"<< (n+i-2)%n <<endl;
//					}
//					comint.push_back( make_pair( (j+1<=n) ? j+1 : 1 , (i-1 > 0) ? i-1 : n) );
				}
			}
//			else{
//				cerr <<"nope "<<trivial<<" "<<(i<j)<<" "<<j-i<<endl;
//			}
			i = rsup[i];
		}
	}


	if(circular){
		vector<pair<int, int> >::iterator it;
		sort( comint.begin(), comint.end() );
		it = unique(comint.begin(), comint.end() );
		comint.erase(it, comint.end());
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void common_intervals(const genom &g1, const genom &g2, int circular, int trivial, int sign, vector<pair<int, int> > &comint){
	genom g;			// the 'identified' genomes (g->id is the same as g1->g2)
	int n = g1.size();
	static vector<int> r,		// the generator
		l,
		s;						// the sign arrays

	if((int) r.size() != n+1){
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
		if(sign != 0)
			s = vector<int>( n, 0 );
	}

	g = g1.identify_g(g2);
//	cout << g<< endl;
		// compute the sign array .. every interval which has the same sign gets an unique number
	if(sign != 0){
		for(unsigned i=1; i<g.size(); i++){
			if( ( g[i] < 0 && g[i-1] < 0 ) || ( g[i] > 0 && g[i-1] > 0 ) ){
				s[i] = s[i-1];
			}else{
				s[i] = s[i-1] + 1;
			}
		}
	}

	generator(g, r, l);
	common_intervals(r, l, circular, trivial, sign, comint);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void common_intervals(const vector<genom> &genomes, int circular, int trivial, int sign, vector<pair<int, int> > &comint ){
	vector<genom> gi = genomes;
	int n = genomes[0].size();
	static vector<int> r,		// the generator
		l;

	if((int)r.size() != n+1){
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
	}

	if(sign != 0){
		cerr << "signed common intervals for more than 2 permutations are unimplemented"<<endl;
		exit(0);
	}

	identify(gi);
	combined_generator(gi, n, r, l);
	common_intervals(r, l, circular, trivial, sign, comint);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
void common_intervals_diff(const vector<genom> &genomes, const genom &gadd,
	int n, int circular, int sign, vector<int> &dest_elements, hdata &hd){

	genom g;
	static vector<int> r,		// the generator
		l,
		radd,		// aditional generator
		ladd;
	vector<pair<int, int> > comint_pre,
		comint_post,
		comint_dif;
	vector<genom> genomesid;
	vector<vector<int> > components;
	vector<int> inner,
		small_components;
	int last_found,
		cnt_found_elms,
		max;

	if((int) r.size() != n+1){
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
		radd = vector<int>( n+1, 0 );
		ladd = vector<int>( n+1, 0 );
	}

	if(sign != 0){
		cerr << "signed common intervals for more than 2 permutations are unimplemented"<<endl;
		exit(0);
	}

		// compute the common intervals of 'genomes'
	combined_generator(genomes, n, r, l);
	common_intervals(r, l, 0, 0, sign, comint_pre);

		// compute the common intervals of 'genomes' \cup 'gadd'
	g = genomes[0].identify_g(gadd);
	generator(g, radd, ladd);
	combine_generator(radd, ladd, r, l);
	common_intervals(radd, ladd, 0, 0, sign, comint_post);

		// compute the set difference, i.e. the destroyed intervals
	sort (comint_pre.begin(), comint_pre.end());
	sort (comint_post.begin(), comint_post.end());

	unique(comint_pre.begin(), comint_pre.end());
	unique(comint_post.begin(), comint_post.end());

	set_difference(comint_pre.begin(), comint_pre.end(),
					comint_post.begin(), comint_post.end(), back_inserter(comint_dif));

	//~ cout <<"comint_dif"<< endl;
	//~ for(unsigned i=0; i<comint_dif.size(); i++){
		//~ cout << "{"<<comint_dif[i].first<<","<<comint_dif[i].second<<"} ";
	//~ }
	//~ cout << endl<<endl;

		// for each destroyed element the elements which destroyed these are sought.
		// the idea is the in gadd some of the elements which were in the destroyed interval
		// are in gadd not in the common interval, i.e. they are outside. This can be seen
		// otherway round, too. i.e. some of the elements which were formerly not in the common
		// interval, are now in it (and destroying it therefore)
		// The problem is to decide which is the more probable. Here the desicion is made by the
		// sizes of the sets A = elements which entered the destroyed common interval
		// and B = elements which moved away

	for(unsigned i=0; i<comint_dif.size(); i++){
		//~ cout << "["<<comint_dif[i].first<<","<<comint_dif[i].second<<"]"<<endl;
		last_found = std::numeric_limits< int >::max();
		cnt_found_elms = 0;
		components.clear();
		inner.clear();
		small_components.clear();
		max = 0;
		for(int j=0; j< (int)g.size(); j++){
			//~ cout << abs(gadd[j])<<" ";
			if(abs(g[j]) >= comint_dif[i].first &&  abs(g[j]) <= comint_dif[i].second){
				//~ cout << "c ";
				if(j == 0 || last_found != j-1){
					//~ cout << "nc ";
					components.push_back( vector<int>() );
				}
				components.back().push_back(abs(gadd[j]));
				last_found = j;
				cnt_found_elms++;
			}else if(cnt_found_elms > 0 && cnt_found_elms < (1 + comint_dif[i].second - comint_dif[i].first)){
				//~ cout << "i";
				inner.push_back(abs(gadd[j]));
			}
		}
		//~ cout << endl;
			// decided which elements destroyed the interval
		for(unsigned j=0; j<components.size(); j++){
			if((int)components[j].size() > max){
				max = components[j].size();
			}
		}
		for(unsigned j=0; j<components.size(); j++){
			//~ if((int)components[j].size() < max){
				small_components.insert(small_components.end(), components[j].begin(), components[j].end());
			//~ }
		}
		//~ cout << inner.size()<< " inner ";
		//~ for(unsigned j=0; j<inner.size(); j++)
			//~ cout << inner[j] <<" ";
		//~ cout << endl;
		//~ cout << small_components.size()<<" small comp ";
		//~ for(unsigned j=0; j<small_components.size(); j++)
			//~ cout << small_components[j] <<" ";
		//~ cout << endl;
		//~ cout << "components"<<endl;
		//~ for(unsigned j=0; j<components.size(); j++){
			//~ for(unsigned k=0; k<components[j].size(); k++)
				//~ cout << components[j][k]<<" ";
			//~ cout <<endl;
		//~ }
		if(small_components.size() < inner.size()){
			for(unsigned j=0; j<small_components.size(); j++)
				dest_elements[ small_components[j] ]++;
		}
		if(inner.size() < small_components.size()){
			for(unsigned j=0; j<inner.size(); j++)
				dest_elements[ inner[j] ]++;
		}
	}

}
*/
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*
extern "C" int common_intervals_diff_c(
	struct genome_struct **genomes_ptr, int m,
	int *g_ptr, int n, int circular, int sign, int add){

	genom g;
	vector<genom> genomes;

	for(int i=0; i<m; i++){
		genomes.push_back( genom(genomes_ptr[i]->genes, n, 0) );
	}
	g = genom(g_ptr, n, 0);

	return common_intervals_diff(genomes, g, n, circular, sign, add);
}
*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int common_intervals_diff(const vector<genom> &genomes, const genom &gadd,
	int n, int circular, int sign, int add){

	genom g;
	static vector<int> r,		// the generator
		l,
		radd,		// aditional generator
		ladd;
	static vector<pair<int, int> > comint;
	vector<genom> gi; 	// identified genomes, i.e. gi[0] = id
	vector<pair<int, int> > comint2;

	if((int) r.size() != n+1){
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
		radd = vector<int>( n+1, 0 );
		ladd = vector<int>( n+1, 0 );
	}

	if(sign != 0){
		cerr << "signed common intervals for more than 2 permutations are unimplemented"<<endl;
		exit(0);
	}

	if(add == 0){
		gi = genomes;
		identify(gi);
		combined_generator(gi, n, r, l);
		common_intervals(r, l, circular, WITHOUTTRIV, sign, comint);
//		for(unsigned i=0; i<comint.size(); i++)
//			cout << "["<<comint[i].first<<","<<comint[i].second<<"] ";
//		cout<<endl;
	}

	if(add == 1){
//		cerr << "gadd "<< gadd<<endl;
//		cerr << "genomes "<< genomes<<endl;
		g = genomes[0].identify_g(gadd);
		generator(g, radd, ladd);
		combine_generator(radd, ladd, r, l);
		common_intervals(radd, ladd, circular, WITHOUTTRIV, sign, comint2);
//		for(unsigned i=0; i<comint2.size(); i++)
//			cout << "["<<comint2[i].first<<","<<comint2[i].second<<"] ";
//		cout<<endl;

	}

	return comint.size() - comint2.size();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void generator(const genom &g, vector<int> &sup, vector<int> &inf){
	int n = g.size();

	stack<int> s;
	static vector<int> lmax,	// left end right bounds of imax
		rmax,
		lmin,			// left end right bounds of imin
		rmin,
		m,				// helping vars for generator computations
		M,
		p,				// the augmented genome (0 ... 0) / (n+1 ... n+1)
		p_inv;			// inverse permutation

	if((int)p.size() != n+2) {
		p = vector<int>(n+2,0);
		p_inv = vector<int>(n+2, 0);
		lmax = vector<int>( n+1, 0 );
		rmax = vector<int>( n+1, 0 );
		lmin = vector<int>( n+1, 0 );
		rmin = vector<int>( n+1, 0 );
		m = vector<int>( n+1, 0 );
		M = vector<int>( n+1, 0 );
	}

		// augment g -> produce p
	for(int i=0; i<n; i++) {
		p[i+1] = abs(g[i]);
	}
	p[0] = 0;
	p[n+1] = 0;
		// compute the inverse
	for(int i=1; i<=n; i++){
		p_inv[p[i]] = i;
	}
	p_inv[n+1] = n+1;
#ifdef DEBUG_GENERATOR
	cout << "inverse"<<endl;
	for(int i=1; i<=n; i++)
		cout << p_inv[i] << " ";
	cout << endl;
#endif
		// compute the left bound lmax[p_i] of the Intervals imax[p_i] for all p_i
	s.push(0);
	for(int i=1; i<=n; i++){
		while( p[i] < p[ s.top() ] )
			s.pop();
		lmax[ p[i] ] = s.top() + 1;
		s.push(i);
	}
	while(s.size() > 0){
		s.pop();
	}
		// compute the right bound rmax[p_i] of the Intervals imax[p_i] for all p_i
	s.push(n+1);
	for(int i=n; i>0; i--){
		while( p[i] < p[ s.top() ] )
			s.pop();
		rmax[ p[i] ] = s.top() - 1;
		s.push(i);
	}
	while(s.size() > 0){
		s.pop();
	}

	p[0] = n+1;
	p[n+1] = n+1;
		// compute the left bound lmin[p_i] of the Intervals imin[p_i] for all p_i
	s.push(0);
	for(int i=1; i<=n; i++){
		while( p[i] > p[ s.top() ] )
			s.pop();
		lmin[ p[i] ] = s.top() + 1;
		s.push(i);
	}
	while(s.size() > 0){
		s.pop();
	}

		// compute the right bound rmin[p_i] of the Intervals imin[p_i] for all p_i
	s.push(n+1);
	for(int i=n; i>0; i--){
		while( p[i] > p[ s.top() ] )
			s.pop();
		rmin[ p[i] ] = s.top() - 1;
		s.push(i);
	}
	while(s.size() > 0){
		s.pop();
	}

#ifdef DEBUG_GENERATOR
	cout << "LMin ";
	for(int i=1; i<=n; i++){
		cout << "["<<lmin[i]<<","<<rmin[i]<<"] ";
	}
	cout << endl;
	cout << "LMax ";
	for(int i=1; i<=n; i++){
		cout << "["<<lmax[i]<<","<<rmax[i]<<"] ";
	}
	cout << endl;
#endif

		// compute the generator (sup, inf)
	inf[1] = 1;
	sup[n] = n;

	for(int k=1; k<=n; k++){
		m[k] = k;
		M[k] = k;
	}
	for(int k=2; k<=n; k++){
		while( lmin[k] <= p_inv[ m[k]-1 ] && p_inv[ m[k]-1 ] <= rmin[k] )
			m[k] = m[ m[k]-1 ];
		inf[k] = m[k];
	}

	for(int k=n-1; k>=1; k--){
		while((lmax[k] <= p_inv[ M[k]+1 ]) && (p_inv[ M[k]+1 ] <= rmax[k]))
			M[k] = M[ M[k]+1 ];
		sup[k] = M[k];
	}

#ifdef DEBUG_GENERATOR
	cout << "generator "<<endl;
	for(int i=1; i<=n; i++)
		cout << "("<< sup[i]<<","<<inf[i]<<")";
	cout <<endl;
#endif
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void front( itnode *n, const genom &g,  genom &f ){

	for(unsigned i=0; i<n->children.size(); i++){
		front( n->children[i], g, f );
	}

	if( n->children.size() == 0 ){
		if( n->i.first < 0 || n->sign == DEC){
			f.push_back( -1*g[ abs(n->i.first) ] );
		}else{
			f.push_back( g[ abs(n->i.first) ] );
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void init_itnode(pair<int, int> &i, itnode* parent, itnode** n){

	*n = new(itnode);
	(*n)->i = i;
	(*n)->parent = parent;

	(*n)->children = vector<itnode*>();
	(*n)->type = UNK;
	(*n)->sign = UNK;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree(vector<pair<int, int> > &strong_int, int n, const genom &g, itnode **itroot){
	unsigned k;
	itnode *itcur,		// current node of the interval tree
		*ittmp;			// temporary for new nodes

	sort_intervals(n, strong_int);

	init_itnode(strong_int[0], NULL, itroot);
	itcur = *itroot;

	k=1;
	while( k< strong_int.size() ){
		if(itcur->i.first <= strong_int[k].first && strong_int[k].second <= itcur->i.second){
			if(strong_int[k] != itcur->i){	// filter duplicate strong intervals
				init_itnode(strong_int[k], itcur, &ittmp);
				itcur->children.push_back(ittmp);
				itcur = ittmp;
			}
			k++;
		}else{
			itcur = itcur->parent;
		}
	}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree(const genom &g1, const genom &g2, int n, itnode **tree_root){
	genom g;			// the 'identified' genomes (g->id is the same as g1->g2)
	static vector<int> r,	// the generator
		l;
	static vector<pair<int,int> > strong_int;	// strong intervals

	if((int) strong_int.size() != 2*n){
		strong_int = vector<pair<int,int> >(2*n, pair<int,int>());
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
	}

		// get the strong intervals of the genomes
	strong_intervals(g1, g2, n, r, l, strong_int);

		// now identify in the correct order ( make the g1 the target and g2 the source
		// -- because the strong intervals .. and therefore the interval tree is computed
		// relative to g2 ) to get the interval tree
	g = g2.identify_g(g1);

	interval_tree( strong_int, n, g, tree_root );
	interval_tree_type( *tree_root, r, l, g );
	interval_tree_sign( *tree_root, g );

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree(vector<genom> genomes, int n, itnode **tree_root ){
	static vector<int> r,
		l;
	static vector<pair<int,int> > strong_int;	// strong intervals

	if((int) strong_int.size() != 2*n){
		strong_int = vector<pair<int,int> >(2*n, pair<int,int>());
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
	}

//	cout << genomes << endl;

	identify(genomes);

//	for(unsigned i=0; i<genomes.size(); i++)
//		cout << genomes[i] << "<br>"<<endl;

	strong_intervals(genomes, n, r, l, strong_int);

//	cout << strong_int.size()<<" strong intervals "<<endl;

	interval_tree( strong_int, n, genomes[0], tree_root );
	interval_tree_type( *tree_root, r, l, genomes[0] );
	interval_tree_sign( *tree_root, genomes[0] );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_copy( itnode *orig, itnode **copy ){
	itnode *tmp;

	init_itnode( orig->i , orig->parent, copy);
	(*copy)->sign = orig->sign;
	(*copy)->type = orig->type;

	for( unsigned i=0; i<orig->children.size(); i++ ){
		interval_tree_copy( orig->children[i], &tmp );
		(*copy)->children.push_back( tmp );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int interval_tree_depth( itnode *p ){
	int cd = 0;
	for( unsigned i=0; i< p->children.size(); i++ ){
		cd = max( cd, interval_tree_depth( p->children[i] ) );
	}
	if(p->children.size() == 0){
		return 0;
	}else{
		return cd + 1;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_free(itnode *p){
	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_free( p->children[i] );
	}
	p->children.clear();
	delete( p );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_get_siblings( itnode *p, vector<itnode *> &siblings ){

	siblings.clear();
	if( p->parent == NULL ){
		return;
	}
	for(unsigned i=0; i<p->parent->children.size(); i++){
		if(p->parent->children[i] == p)
			continue;
		siblings.push_back(p->parent->children[i]);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool interval_tree_is_prime( const genom &g1, const genom &g2 ){
	int n = g1.size();
	itnode *itroot = NULL;
	list<itnode *> pnode;

	interval_tree(g1, g2, n, &itroot);
	interval_tree_primenodes(itroot, pnode);
	interval_tree_free(itroot);
	// note that the content of pnode must not be accessed after this point
	// but we still check if the list of prime nodes is empty
	return (not pnode.empty());
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_leafs(itnode *p, vector<itnode *> &node){

	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_leafs( p->children[i], node);
	}

	if(p->children.size() == 0)
		node.push_back( &(*p) );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_linearnodes(itnode *p, list<itnode *> &lnode){
	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_linearnodes( p->children[i], lnode);
	}

	if( p->children.size()> 0 && p->parent != NULL && p->type == LIN ){
		lnode.push_back( &(*p) );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_primenodes(itnode *p, list<itnode *> &pnode){
	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_primenodes( p->children[i], pnode);
	}

	if( p->type == PRI ){
		pnode.push_back( p );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_pnode_components(itnode *p, vector<vector<itnode *> > &pcomp){
	vector<itnode *> ptop;

	interval_tree_pnode_top(p, ptop);
	pcomp = vector<vector<itnode *> >(ptop.size());

	for(unsigned i=0; i<ptop.size(); i++){
		interval_tree_pnode_fill( ptop[i], pcomp[i] );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_pnode_fill(itnode *p, vector<itnode *> &pcomp){

	if(p->type == PRI){
		pcomp.push_back( p );
	}

	for( unsigned i=0; i<p->children.size(); i++ ){
		if(p->children[i]->type == PRI)
			interval_tree_pnode_fill( p->children[i], pcomp);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_pnode_top(itnode *p, vector<itnode *> &ptop){
	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_pnode_top( p->children[i], ptop);
	}

	if(p->type == PRI && (p->parent == NULL || p->parent->type == LIN) ){
		ptop.push_back( p );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_print(itnode *p, const genom &g, ostream &o, string prefix){

	interval_tree_print_prefix(p, o, prefix);

	if( p->children.size() == 0 ){
		if ( p->i.first < 0 ){
			print_element(g[(-1*p->i.first)], o, -1, "", g.get_nmap());
		}else{
			print_element(g[(p->i.first)],  o,  1, "", g.get_nmap());
		}
//		cout << p->i.first ;
	}else{
		for(unsigned i=0; i< p->children.size(); i++){
			interval_tree_print( p->children[i], g, o);
			if( i < p->children.size()-1 )
				o << ",";
		}
	}

	interval_tree_print_suffix(p, o);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_print_prefix( itnode *p, ostream &o, string prefix){
	o << prefix;
	switch( p->sign ){
		case INC:{o << "+";break;}
		case DEC:{o << "-";break;}
		case UNK:{o << "?";break;}
	}
	o << "<"<<abs(p->i.first) <<"-"<<abs(p->i.second)<<">";
	switch (p->type){
		case LIN: {o << "["; break;}
		case PRI: {o << "("; break;}
		case UNK: {o << "{"; break;}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_print_suffix( itnode *p, ostream &o ){
	switch (p->type){
		case LIN: {o << "]"; break;}
		case PRI: {o << ")"; break;}
		case UNK: {o << "}"; break;}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void interval_tree_print_dot(itnode *p, const genom &g, string fname, string ext, const vector<string> &namemap){
//	ofstream o;
//	string dfname,
//		call;
//
//	dfname = fname + ".dot";
//	o.open(dfname.c_str());
//	interval_tree_print_dot(p, g, o, namemap);
//	o.close();
//
//	call = "dot -T"+ext+" "+dfname+" -o "+fname+"."+ext;
//	system(call.c_str());
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_print_dot(itnode *p, const genom &g, ostream &o, const vector<string> &namemap){
	if( p->parent == NULL ){

		o << "digraph g{"<<endl;
		o << "ordering=out;"<<endl;
		o << "nodesep=0.05;"<<endl;
		o << "fontname=Helvetica;"<<endl;
		o << "graph[fontsize=20];"<<endl;
		o << "size=\"32,8\";";
		o << "ranksep=0.2;"<<endl;

		o << "{ rank = same; ";

		for( unsigned i=0; i<g.size(); i++){
			 o << "\"f"<< i+1 <<"t"<<i+1<<"\"; ";
		}
		o << "}"<<endl;
	}

	o << "f"<< p->i.first<<"t"<<p->i.second<<" [";
	switch (p->type){
		case LIN: {
			o<<"shape=box, ";
			if(p->sign == INC)
				o<<"color=red";
			else
				o<<"color=green";
			break;
		}
		case PRI: {o<<"shape=ellipse, color=blue"; break;}
	}

	if(p->children.size() == 0){
		int len;
//		len = int2string(g[ p->i.first - 1]).size();
		len = namemap[abs(g[ p->i.first])].size();
		if(g[ p->i.first ] < 0)
			len++;

		o << ", fixedsize=true, width="<<len*0.15<<", height=0.25";
		o <<",label=\"";
		print_element(g[ p->i.first ], o, 1, "", &namemap);
	}else{
		o << ",height=0.25,label=\"";
		for(int i=p->i.first  ; i<=p->i.second; i++){
			print_element(g[i], o, 1, "", &namemap);
		}
	}
//	o <<", fixedsize=true, width=0.0001, height=0.0001, label=\"\"]"<< endl;
	o << "\"];";
	if( p->parent != NULL ){
		o << "f"<<p->parent->i.first<<"t"<<p->parent->i.second <<" -> "<< "f"<< p->i.first<<"t"<<p->i.second<< "[dir=none]"<<endl;
	}

	for(unsigned i=0; i<p->children.size(); i++){
		interval_tree_print_dot( p->children[i], g, o, namemap );
	}
	if( p->parent == NULL ){
		o << "}"<<endl;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_print_table(itnode *p, const genom &g, const vector<string> &namemap, const string prefix, ostream &o, bool inv ){

	string cls="";	// css class
	int start=0,
		end=p->children.size(),
		dir=1;

	if(inv){
		start = p->children.size()-1;
		end = -1;
		dir = -1;
	}

	if(p->type == PRI){
		cls="pri";
	}else{
		if((!inv && p->sign == INC) || (inv && p->sign == DEC)){
			cls = "inc";
		}else if((!inv && p->sign == DEC) || (inv && p->sign == INC) ){
			cls = "dec";
		}
	}

		// print a table with a class (for css formating): pri/inc/dec
		// and an id for accessing an interval in an permutation: prefix_start_end
	o << "<table border=1 rules=\"none\" class=\""<<cls<<"\"";
	o << " id=\""<<prefix<<"_"<< p->i.first<<"_"<<p->i.second <<"\">";

	o<<"<tr>"<<endl;
	for(int i=start; i!=end; i+=dir){
		o << "<td>"<<endl;
//		o << "<td style=\"padding:3px;\">"<<endl;
		interval_tree_print_table( p->children[i], g, namemap, prefix, o, inv );
		o<< "</td>"<<endl;
	}

	if(p->children.size() == 0){
		o << "<td>"<<endl;
//		o << "<td style=\"padding:3px;\">"<<endl;
		if (inv){
			print_element(g[p->i.first], o, -1, "&nbsp", &namemap);
		}else{
			print_element(g[p->i.first], o, 1, "&nbsp;", &namemap);
		}
		o<< "</td>"<<endl;
	}

	o<< "</tr>"<<endl;
	o<< "</table>"<<endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//void interval_tree_properties(itnode *p, const genom &g, unsigned n, unsigned &dist, dstnc_inv *rd){
//
//	for( unsigned i=0; i<p->children.size(); i++ ){
//		interval_tree_properties( p->children[i], g, n, dist, rd);
//	}
//
//		// prime node with linear parent (i.e. known sign)
//		// compute the distance depending on the sign of the parent
//	if(p->type == PRI){
//		dist += prime_distance(p, n, g, rd);
//	}
//		// linear node
//		// count the number of child nodes with different sign
//		// (so far the orientation of the root does not matter, i.e. linear undirected genome model)
//	else {
//		if( p->parent == NULL ){
//				// currently i want the linear directed genome model
//				// so check if the (linear) root is decreasing, if so add a reversal
//				// (comment out to get linear undirected)
//			if ( p->sign == DEC ){
//#ifdef DEBUG_PERFECTD
//				cout << "sign diff root"<<endl;
//#endif//DEBUG_PERFECTD
//				dist ++;
//			}
//		}else{
//			if( p->parent->type != PRI && ( p->sign != p->parent->sign ) ){
//#ifdef DEBUG_PERFECTD
//				cout << "sign diff "<<"["<<p->i.first <<","<<p->i.second<<"] - ["<<p->parent->i.first<<","<<p->parent->i.second<<"]"<<endl;
//#endif//DEBUG_PERFECTD
//				dist++;
//			}
//		}
//	}
////	cout << "["<<p->i.first <<","<<p->i.second<<"] :"<<dist<<endl;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_random( itnode **r, itnode *p, pair<int,int> i, unsigned dm, unsigned dM ){
	itnode *tmp;
	unsigned deg;
	vector<unsigned> split;

	dM = min( dM, (unsigned)(i.second - i.first + 1));
	dm = min( dm, (unsigned)(i.second - i.first + 1));

	deg = ask_rng( dm, dM );

	init_itnode( i, p, r );

	(*r)->type = LIN;
	(*r)->sign = ask_rng(0,1);
	if( i.second <= i.first ){
		return;
	}

	deg = min( deg, (unsigned)(i.second - i.first + 1) );
	split = rng_inc_seq( deg-1, i.second-i.first, std::numeric_limits< unsigned >::max());
	for(unsigned j=0; j<split.size(); j++){
		split[j] += i.first+1;
	}
	split.insert(split.begin(), i.first);
	split.push_back( i.second+1 );
	for(unsigned j=0; j<split.size()-1; j++){
		interval_tree_random( &tmp, (*r), make_pair( split[j], split[j+1]-1), dm, dM );
		(*r)->children.push_back( tmp );
	}

	// if all children have the same sign the node must have opposite sign
	bool esign = true;
	for(unsigned i=0; i<(*r)->children.size()-1; i++){
		if( (*r)->children[i]->sign != (*r)->children[i+1]->sign ){
			esign = false;
			break;
		}
	}
	if( esign ){
		(*r)->sign = 1-(*r)->children[0]->sign;
	}

//	if( (*r)->sign == DEC ){
//		reverse( (*r)->children.begin(), (*r)->children.end() );
//	}


}

void interval_tree_reverse(itnode *n){
	for( unsigned i=0; i<n->children.size(); i++){
		interval_tree_reverse(n->children[i]);
	}

	if(n->sign == INC){
		n->sign = DEC;
	}else if(n->sign == DEC){
		n->sign = INC;
	}

	if(n->children.size() > 0){
		reverse(n->children.begin(), n->children.end());
	}
	else{
		n->i.first *= -1;
		n->i.second *= -1;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_sign( itnode *n, const genom &g){

	if( n->type == PRI ){
		if(n->parent == NULL){
			if ( n->parent == NULL ){		// ... the root is positive to implement the
				n->sign = INC;				// linear directed model (comment out for linear undirected)
			}
		}else{
			if(n->parent->type == PRI){		//...  without linear parent have an
				n->sign = UNK;				// 	unknown sign
			}else if(n->parent != NULL){	// ... with linear parent have the sign
				n->sign = n->parent->sign;	// 	of the parent
			}
		}
	}

	for( unsigned i=0; i<n->children.size(); i++ ){
		interval_tree_sign( n->children[i], g);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*void interval_tree_sorting_reversals(itnode *p, const genom &g, int n,
		vector<pair<int,int> > &reversals, dstnc_inv *rd){

	for( unsigned i=0; i<p->children.size(); i++ ){
		interval_tree_sorting_reversals( p->children[i], g, n, reversals, rd);
	}

		// prime node with linear parent (i.e. known sign)
		// compute the distance depending on the sign of the parent
	if(p->type == PRI){
		prime_sorting_reversals(p, n, g, reversals, rd);
//		dist += prime_distance(p, n, g, hd);
	}
		// linear node
		// count the number of childnodes with different sign
		// (so far the orientation of the root does not matter, i.e. linear undirected genome model)
	else {
		if( p->parent == NULL ){
				// currently i want the linear directed genome model
				// so check if the (linear) root is decreasing, if so add a reversal
				// (comment out to get linear undirected)
			if ( p->sign == DEC ){
#ifdef DEBUG_SORTINGREV
				cout << "sign diff root"<<endl;
#endif//DEBUG_SORTINGREV
				reversals.push_back( p->i );
			}
		}else{
			if( p->parent->type != PRI && ( p->sign != p->parent->sign ) ){
#ifdef DEBUG_SORTINGREV
				cout << "sign diff "<<"["<<p->i.first <<","<<p->i.second<<"] - ["<<p->parent->i.first<<","<<p->parent->i.second<<"]"<<endl;
#endif//DEBUG_SORTINGREV
				reversals.push_back( p->i );
			}
		}
	}
//	cout << "["<<p->i.first <<","<<p->i.second<<"] :"<<dist<<endl;

}*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*void interval_tree_sorting_reversals(const genom g1, const genom g2,
		vector<pair<int,int> > &reversals){

	int n = (int)g1.size();
	itnode* tree_root;	// the tree (accesible by its root)
	genom g;			// the 'identified' genomes (g->id is the same as g1->g2)
	static vector<int> r,	// the generator
		l;
	static vector<pair<int,int> > strong_int;	// strong intervals
	dstnc_inv *rd = new dstnc_inv(n, false);

	if((int) strong_int.size() != 2*n){
		strong_int = vector<pair<int,int> >(2*n, pair<int,int>());
		r = vector<int>( n+1, 0 );
		l = vector<int>( n+1, 0 );
	}

		// clear the reversals vector, just in case that it is not empty
	reversals.clear();

		// get the strong intervals of the genomes
	strong_intervals(g1, g2, n, r, l, strong_int);

		// now identify in the correct order ( make the g1 the target and g2 the source
		// -- because the strong intervals .. and therefore the interval tree is computed
		// relative to g2 ) to get the interval tree
	g = g2.identify_g(g1);

	interval_tree( strong_int, n, g, &tree_root );
	interval_tree_type( tree_root, r, l, g );
	interval_tree_sign( tree_root, g );
#ifdef DEBUG_SORTINGREV
	cout << "g "<<g<<endl;
	cout << "strong intervals"<<endl;
	for(unsigned i=0; i<strong_int.size(); i++)
		cout << "["<< strong_int[i].first <<","<<strong_int[i].second<<"] ";
	cout << endl;
	cout << "interval tree"<<endl;
	interval_tree_print(tree_root);
#endif

	interval_tree_sorting_reversals(tree_root, g, n, reversals, rd);

		// free memory
	interval_tree_free(tree_root);

#ifdef DEBUG_SORTINGREV
	cout << "distance "<<dist<<endl;
#endif//DEBUG_PERFECTD
	delete rd;

}*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_switch_sign( itnode *p ){

	if( p->sign == INC ){
		p->sign = DEC;
	}else if( p->sign == DEC ){
		p->sign = INC;
	}

	for(unsigned i=0; i<p->children.size(); i++){
		interval_tree_switch_sign( p->children[i] );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void interval_tree_type( itnode *n,  const vector<int> &r, const vector<int> &l, const genom &g ){

		// init node as P node
	n->type = PRI;
		// determine the type of nodes
		// if the node has no children (leaf) or 2 children -> linear (just mark as INC)
		// else if the union of the first two child nodes is an interval -> linear
	if( n->children.size() > 2){
		int i = n->children[0]->i.first  +1,
			j = n->children[1]->i.second +1;

			// if (i..j) is an interval -> it is an Q node
		if(l[j] <= i && i <= j && j <= r[i]){
			n->type = LIN;
		}
	}
		// nodes with 2 or 0 childs are autimatically linear nodes (1 child is impossible)
	else if ( n->children.size() <= 2 ){
		n->type = LIN;
	}

		// decide orientation of linear nodes
	if(n->type == LIN){
		if( n->children.size() > 0 ){
			if( abs(g[ n->children[0]->i.first ]) < abs(g[ n->children[1]->i.first ]) ){
				n->sign = INC;
			}else{
				n->sign = DEC;
			}
		}else{							// leaf nodes have the sign of the
			if( g[n->i.first] > 0 )	// 	corresponding element
				n->sign = INC;
			else
				n->sign = DEC;
		}
	}

		// recursion
	for(unsigned i = 0; i< n->children.size(); i++){
		interval_tree_type( n->children[i], r, l, g );
	}

	/*	// if the node has children :
		// - check for the first 2 children if inc / dec
		// - check for subsequent neighbours if inc / dec
		// - if a change occurs -> set to pri
		// the check can be done in O(1) by checking two arbitrary elements
		// from each interval, because the id is in the genom set

	//~ cout << "n ["<<n->i.first<<","<<n->i.second<<"] "<< n->children.size() <<" children" <<endl;
	if(n->children.size()>1){
		//~ cout << " "<<abs(g[ n->children[0]->i.first-1 ]) <<"?"<< abs(g[ n->children[1]->i.first-1 ]);
		if( abs(g[ n->children[0]->i.first-1 ]) < abs(g[ n->children[1]->i.first-1 ]) ){
			n->type = INC;
		}else{
			n->type = DEC;
		}
		//~ cout << " "<<n->type;

		for(unsigned k=2; k<n->children.size(); k++){
			if( n->type == INC && abs(g[ n->children[k-1]->i.first-1 ]) > abs(g[ n->children[k]->i.first-1 ])
				|| n->type == DEC && abs(g[ n->children[k-1]->i.first-1 ]) < abs(g[ n->children[k]->i.first-1 ]) ) {
				n->type = PRI;
					break;
			}
		}
		//~ cout << endl;
	}else{
		if( g[ n->i.first-1 ] > 0){
			n->type = INC;
		}else{
			n->type = DEC;
		}
	}
		// recursion
	for( unsigned i=0; i<n->children.size(); i++ ){
		interval_tree_type( n->children[i], g );
	}*/
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void interval_tree_get_subpermutation( itnode *p, const genom start, genom &g, const unsigned insert_pos){
	if(p->children.size()==0){
		//case leaf in SIT
		if ( p->i.first < 0 ){
			g[insert_pos]=start[ (-1)*p->i.first ];
			//cout << "add element: " << start[ (-1)*p->i.first ] << " on position: " << (-1)*p->i.first << endl;
		}else{
			g[insert_pos]=start[ p->i.first ];
			//cout << "add element: " << start[ p->i.first ] << " on position: " << p->i.first << endl;
		}
	}else{
		//case inner node
		genom test;
		int size=0;
		for(unsigned k=0; k<p->children.size(); k++){
				size+=p->children[k]->i.second - p->children[k]->i.first +1;
		}
		test=genom(size, 0);
		unsigned cur_pos=0;
		for(unsigned k=0; k<p->children.size(); k++){
			interval_tree_get_subpermutation(p->children[k], start, test, cur_pos);

			for(unsigned copy=0; copy<p->children[k]->i.second - p->children[k]->i.first +1; copy++){
				g[cur_pos+copy+insert_pos]=test[cur_pos+copy];
			}
			cur_pos+=p->children[k]->i.second - p->children[k]->i.first +1;
		}
	}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void irreducible_intervals(itnode *nd, vector<pair<int,int> > &ii){

	for( unsigned i=0; i<nd->children.size(); i++ ){
		irreducible_intervals( nd->children[i], ii );
	}

	if( nd->type == LIN && nd->children.size()>0 ){
		for( unsigned j=0; j<nd->children.size()-1; j++ ){
			ii.push_back( pair<int,int>( nd->children[j]->i.first, nd->children[j+1]->i.second) );
		}
	}else if( nd->type == PRI ){
		ii.push_back( nd->i );
	}
//	(nd->type == LIN && (nd->i.second - nd->i.first == 1 )) ||
//	(nd->type == PRI && (nd->parent == NULL || nd->parent->type == PRI) )
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void irreducible_intervals(const genom &g1, const genom &g2, int n, vector<pair<int,int> > &ii){
	itnode *tree_root;
	interval_tree(g1, g2, n, &tree_root);
	ii.clear();
	irreducible_intervals(tree_root, ii);
	interval_tree_free(tree_root);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void irreducible_intervals(vector<genom> genomes, int n, vector<pair<int,int> > &ii){
	itnode *tree_root;
	interval_tree(genomes, n, &tree_root);
	ii.clear();
	irreducible_intervals(tree_root, ii);
	interval_tree_free(tree_root);
}


//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// compute the perfect distance between g1 and g2, i.e. the minimal number
//// of reversals, which do not destroy the common intervals of g1 and g1,
//// needed to transform g1 into g2
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//int perfect_distance(const genom g1, const genom g2){
//	unsigned n = g1.size(),
//		dist = 0;		// the computed distance
//	itnode* tree_root;	// the tree (accesible by its root)
//	genom g;			// the 'identified' genomes (g->id is the same as g1->g2)
//	dstnc_inv *rd;
//	static vector<int> r,	// the generator
//		l;
//	static vector<pair<int,int> > strong_int;	// strong intervals
//
//	rd = new dstnc_inv(n, false);
//	if((int) strong_int.size() != 2*n){
//		strong_int = vector<pair<int,int> >(2*n, pair<int,int>());
//		r = vector<int>( n+1, 0 );
//		l = vector<int>( n+1, 0 );
//	}
//
//		// get the strong intervals of the genomes
//	strong_intervals(g1, g2, n, r, l, strong_int);
//
//		// now identify in the correct order ( make the g1 the target and g2 the source
//		// -- because the strong intervals .. and therefore the interval tree is computed
//		// relative to g2 ) to get the interval tree
//	g = g2.identify_g(g1);
//	interval_tree( strong_int, n, g, &tree_root );
//	interval_tree_type( tree_root, r, l, g );
//	interval_tree_sign( tree_root, g );
//#ifdef DEBUG_PERFECTD
//	cout << "g "<<g<<endl;
//	cout << "strong intervals"<<endl;
//	for(unsigned i=0; i<strong_int.size(); i++)
//		cout << "["<< strong_int[i].first <<","<<strong_int[i].second<<"] ";
//	cout << endl;
//	cout << "interval tree"<<endl;
//	interval_tree_print(tree_root, g, cout, "");
//#endif
//	interval_tree_properties(tree_root, g, n, dist, rd);
//		// free memory
//	interval_tree_free(tree_root);
//	delete rd;
//#ifdef DEBUG_PERFECTD
//	cout << "distance "<<dist<<endl;
//#endif//DEBUG_PERFECTD
//	return dist;
//}



//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// compute the perfect distance between g1 and g2 where the common intervals
//// of genomes have to be preserved
//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//int perfect_distance_glob(const genom &g1, const genom &g2, const vector<genom> &genomes, int n){
//	genom g;									// the 'identified' genomes (g->id is the same as g1->g2)
//	unsigned dist = 0;								// the distance to compute
//	itnode* tree_root;	// the tree (accesible by its root)
//	static vector<int> r,
//		l;
//	static vector<pair<int,int> > strong_int;	// strong intervals
//	vector<genom> gi;
//
//	dstnc_inv *rd;
//
//	rd = new dstnc_inv(n, false);
//
//	if( (int)strong_int.size() != 2*n ){
//		strong_int = vector<pair<int,int> >(2*n, pair<int,int>());
//		r = vector<int>( n+1, 0 );
//		l = vector<int>( n+1, 0 );
//	}
//
//		// construct the genome vector consisting of g1, g2, genome_1...genome_m
//		// in order to compute the intervals relative to g1
//	gi.push_back(g1);
//	gi.push_back(g2);
//	gi.insert(gi.end(), genomes.begin(), genomes.end());
//
//		// get the strong intervals relative to gi[0], i.e. g1
//	strong_intervals(gi, n, r, l, strong_int );
//
//			// now identify in the correct order ( make the g1 the target and g2 the source
//			// -- because the strong intervals .. and therefore the interval tree is computed
//			// relative to g2 ) to get the interval tree
//	g = gi[1].identify_g(gi[0]);
//
//	interval_tree( strong_int, n, g, &tree_root );
//	interval_tree_type( tree_root, r, l, g );
//	interval_tree_sign( tree_root, g );
//#ifdef DEBUG_PERFECTD
//	cout << "g "<<g<<endl;
//	cout << "strong intervals"<<endl;
//	for(unsigned i=0; i<strong_int.size(); i++)
//		cout << "["<< strong_int[i].first <<","<<strong_int[i].second<<"] ";
//	cout << endl;
//	cout << "interval tree"<<endl;
//	interval_tree_print(tree_root, g, cout, "");
//#endif
//
//	interval_tree_properties(tree_root, g, n, dist, rd);
//			// free memory
//	interval_tree_free(tree_root);
//	delete rd;
//	return dist;
//
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// perfect distance computation for c programs

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// computes the sum of the perfect distances between g1 and the genomes in
// genomes, where the common intervals of \f ${g_1, genomes}$ \f are
// preserved, the normal use of this function is to compute the distance sum
// of some permutation to their (common interval preserving) median
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//int perfect_distance_glob_sum(const genom &g1, const vector<genom> &genomes, int n ){
//
//	genom g;			// the 'identified' genomes (g->id is the same as g1->g2)
//	unsigned dist = 0;		// the distance to compute
//	itnode* tree_root;	// the tree (accesible by its root)
//	dstnc_inv *rd;
//	static vector<int> r,
//		l;
//	static vector<pair<int,int> > strong_int;	// strong intervals
//	vector<genom> gi;
//
//	rd = new dstnc_inv(n, false);
//	if( (int)strong_int.size() != 2*n ){
//		strong_int = vector<pair<int,int> >(2*n, pair<int,int>());
//		r = vector<int>( n+1, 0 );
//		l = vector<int>( n+1, 0 );
//	}
//
//		// construct the genome vector consisting of g1, genome_1...genome_m
//		// in order to compute the intervals relative to g1
//	gi.push_back(g1);
//	gi.insert(gi.end(), genomes.begin(), genomes.end());
//
//		// get the strong intervals relative to gi[0], i.e. g1
//	strong_intervals(gi, n, r, l, strong_int);
//
//	for(unsigned i=1; i<gi.size(); i++){
//			// now identify in the correct order ( make the g1 the target and g2 the source
//			// -- because the strong intervals .. and therefore the interval tree is computed
//			// relative to g2 ) to get the interval tree
//		g = gi[i].identify_g(gi[0]);
//
//		interval_tree( strong_int, n, g, &tree_root );
//		interval_tree_type( tree_root, r, l, g );
//		interval_tree_sign( tree_root, g );
//#ifdef DEBUG_PERFECTD
//		cout << "g "<<g<<endl;
//		cout << "strong intervals"<<endl;
//		for(unsigned i=0; i<strong_int.size(); i++)
//			cout << "["<< strong_int[i].first <<","<<strong_int[i].second<<"] ";
//		cout << endl;
//		cout << "interval tree"<<endl;
//		interval_tree_print(tree_root, g, cout, "");
//#endif
//		interval_tree_properties(tree_root, g, n, dist, rd);
////		cout << "dist "<<dist<<endl;
//			// free memory
//		interval_tree_free(tree_root);
//
//	}
//	delete rd;
//	return dist;
//}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*void prime_sorting_reversals(itnode *p, int n, const genom &g, vector<pair<int,int> > &reversals, dstnc_inv *rdist){
	genom id,
		nid;
	unsigned pdist = std::numeric_limits< unsigned >::max(),		// best distance to the positive identity
		ndist = std::numeric_limits< unsigned >::max(),		// best distance to the negative identity
		d;
	int options,
		dir = 0;
	static vector<int> ambiguous_signs,	// current sign assignment to children with unknown sign
		ambiguous_pos;					// indices of children with unknown sign
	int ambiguous_cnt = 0;				// number of childnodes with unknown sign
	vector<pair<int,int> > rev;			// the best reversals
	vector<bool> pgoodrev,
		ngoodrev;

	if((int) ambiguous_signs.size() != n){
		ambiguous_signs = vector<int>(n, INC);
		ambiguous_pos = vector<int>(n, 0);
	}

#ifdef DEBUG_SORTINGREV
	cout << "sorting reversals for ["<<p->i.first<<","<<p->i.second<<"] typ "<< p->type<<" sig "<< p->sign <<endl;
#endif

		// get quotient permutation
		// - take an arbitrary element of each children (e.g. the first)
		// - compute the quotient permutation
		// - assign the signs of the children to the corresponding elements of the
		// 	quotient permutation (if unknown add the to the ambiguos positions)
	genom quot;
	quotient_permutation(p, g, quot, 0, 0);

	for(unsigned i=0; i<p->children.size(); i++){
		if(p->children[i]->sign == DEC)
			quot[i] *= -1;
		else if( p->children[i]->sign == UNK ){
			ambiguous_signs[ambiguous_cnt] = INC;
			ambiguous_pos[ambiguous_cnt] = i;
			ambiguous_cnt++;
		}
	}

	id = genom(quot.size(), 0);
	nid = genom_nid(quot.size());

	rev = id.getReversals();
	pgoodrev = vector<bool>(rev.size(), false);
	ngoodrev = vector<bool>(rev.size(), false);

	options = pow(2, ambiguous_cnt);
#ifdef DEBUG_PRIMED
	cout << ambiguous_cnt<<" ambig. positions -> #options = "<<options <<endl;
	cout << "ambiguous positions: ";
	for(int j=0; j<ambiguous_cnt; j++){
		cout << ambiguous_pos[j] <<" ";
	}cout << endl;
#endif

	for(int i=0; i<options; i++){
			// assign the current sign assignment to the ambiguous positions
		for(int j=0; j<ambiguous_cnt; j++){
			if(quot[ ambiguous_pos[j] ] > 0 && ambiguous_signs[j] == DEC)
				quot[ ambiguous_pos[j] ] *= -1;
			if(quot[ ambiguous_pos[j] ] < 0 && ambiguous_signs[j] == INC)
				quot[ ambiguous_pos[j] ] *= -1;
		}

		if(p->sign == INC || p->sign == UNK){
			d = rdist->calc(quot, id);
			if(d<=pdist){
				if(d<pdist){
					pdist = d;
					pgoodrev.assign( rev.size(), false );
				}
				for( unsigned r=0; r<rev.size(); r++ ){
					reverse(quot, rev[r]);
					if( rdist->calc( quot, id) < d )
						pgoodrev[r] = true;
					reverse(quot, rev[r]);
				}
			}
		}
			// if negative -> sort to nedative identity d( q, -id ) = d( -q, id )
		if(p->sign == DEC || p->sign == UNK){
			d = rdist->calc(quot, nid);
			if(d<=ndist){
				if(d<ndist){
					ndist = d;
					pgoodrev.assign( rev.size(), false );
				}
				for( unsigned r=0; r<rev.size(); r++ ){
					reverse(quot, rev[r]);
					if( rdist->calc(quot,nid) < d )
						ngoodrev[r] = true;
					reverse(quot, rev[r]);
				}
			}
		}

			// compute the next possible sign assignment
		ambiguous_signs[0]++;
		for(int j=0; j<ambiguous_cnt; j++){
			if(ambiguous_signs[j] == DEC || ambiguous_signs[j] == INC){
				break;
			}else{
				ambiguous_signs[j] = INC;
				if(j+1<ambiguous_cnt)
					ambiguous_signs[j+1]++;
			}
		}
	}
#ifdef DEBUG_SORTINGREV
	cout << "pdist  = "<<pdist<<" ndist = "<<ndist <<endl;
#endif
		// if the sign of the node is known (only possible if it is the child of a linear node)
		// than take the minimal distance to the positive/negative identity, depending on the
		// sign of the parent node
		// if the sign is unknown take the minimal distance and set the sign of the
		// node if the distances are unequal, else the sign stays unknown
	if(p->sign == INC){
		d = pdist;
		dir = 1;
	}else if(p->sign == DEC){
		d = ndist;
		dir = -1;
	}else{
		if(pdist < ndist){
			p->sign = INC;
			d = pdist;
			dir = 1;
		}else if(ndist < pdist){
			p->sign = DEC;
			d = ndist;
			dir = -1;
		}else{
			d = pdist;
			dir = 0;
		}
	}
#ifdef DEBUG_SORTINGREV
	cout << "prime distance = "<<d <<endl;
#endif

	for(unsigned i=0; i<rev.size(); i++){
		if( (dir >= 0 && pgoodrev[i] == true) || (dir <= 0 && ngoodrev[i] == true) ){
			reversals.push_back( make_pair(p->children[ rev[i].first ]->i.first, p->children[ rev[i].second ]->i.second) );
		}
	}
}*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void quotient_permutation(itnode *p, const genom &g, genom &qperm, int circular, int sign){
	vector<itnode *> sib;

	qperm.chromosom.clear();

	qperm.chromosom = vector<int>(p->children.size(), 0);
		// construct the quotient permutation of the children
	for(unsigned i=0; i<p->children.size(); i++){
		qperm[i] = abs(g[ p->children[i]->i.first ]);
	}

		// if circular -> 1 sibling has to be added
	if (circular != 0){
		interval_tree_get_siblings(p, sib);
//		cout <<sib.size()<< "siblings"<<endl;
		if (sib.size() > 0){	// happens if the node is the root .. is ok
			qperm.push_back(abs( g[ sib[0]->i.first ] ));
		}
	}

//	genom f;
//	front( p, g,  f );
	quotient_permutation(qperm, g.size());

		// make it signed
	if(sign == 0){
		return;
	}
	for (unsigned i=0; i<p->children.size(); i++){
		if( p->children[i]->type == LIN && p->children[i]->sign == DEC ){
			qperm[i] *= -1;
		}
	}
	if (circular != 0 && sib.size() > 0){
		if ( p->parent->type == LIN && p->parent->sign == DEC){
			qperm[qperm.size()-1] *= -1;
		}
	}
//	cout << " -> "<<qperm<<endl;
	if(circular != 0){
		qperm.normalize(1);
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void quotient_permutation(genom &quot, int n){
	int k = 0;
	static vector<int> bucket;

	if((int)bucket.size() != n+1){
		bucket = vector<int>(n+1, std::numeric_limits< int >::max());
	}

		// put the index of each element from pi into the bucket of the
		// corresponding bucket
	for(unsigned i=0; i<quot.size(); i++){
		bucket[ abs(quot[i]) ] = (i+1);
	}

		// go through the buckets:
		// - if a bucket is empty -> just continue
		// - if there is something -> overwrite the element in quot at the position
		// 	which is stored in the bucket, with the number of
		// 	nonempty buckets +1 which are found so far
	for(unsigned i=1; i<bucket.size(); i++){
		if(bucket[i] == std::numeric_limits< int >::max()){
			continue;
		}
		quot[ bucket[i]-1 ] = k+1;
		bucket[i] = std::numeric_limits< int >::max();
		k++;
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sort_intervals(int n, vector<pair<int, int> > &intervals){

		// decreasing bucket sort the intervals according to their right bound
	stable_sort(intervals.begin(), intervals.end(), cmp_interval_second);
	reverse( intervals.begin(), intervals.end() );
		// increasing bucket sort the intervals according to their left bound
	stable_sort(intervals.begin(), intervals.end(), cmp_interval_first);

	return;

	int j;
	static vector<vector<pair<int, int> > > bucket;
	static vector<int> bucket_size;

#ifdef DEBUG_SORTINTERVALS
	cout << "unsorted :"<<endl;
	for(unsigned i=0; i<intervals.size(); i++)
		cout << "("<<intervals[i].first<<","<<intervals[i].second<<")";
	cout << endl;
#endif

	if((int) bucket.size() != n+1){
//		bucket = vector<vector<pair<int,int> > >(n+1, vector<pair<int,int> >(((n+1)*n) / 2 , pair<int,int>()));
		bucket = vector<vector<pair<int,int> > >(n+1, vector<pair<int,int> >(n+1 , pair<int,int>()));
		bucket_size = vector<int>(n+1, 0);
	}

		// decreasing bucket sort the intervals according to their right bound
	for(unsigned i=0; i<intervals.size(); i++){
//		if (i>0 && intervals[i].first == intervals[i-1].first && intervals[i].second == intervals[i-1].second)
//			continue;

		bucket[ intervals[i].second ][ bucket_size[intervals[i].second] ] = intervals[i];
		bucket_size[intervals[i].second]++;
	}
	j=0;
	for(int i=n; i>0; i--){
		for(int k=bucket_size[i]-1; k>=0; k--){
			intervals[j] = bucket[i][ k ];
			j++;
		}
		bucket_size[i] = 0;
	}

#ifdef DEBUG_SORTINTERVALS
	cout << "decreasing sorted by right bound :"<<endl;
	for(unsigned i=0; i<intervals.size(); i++)
		cout << "("<<intervals[i].first<<","<<intervals[i].second<<")";
	cout << endl;
#endif
			// increasing bucket sort the intervals according to their left bound
	for(unsigned i=0; i<intervals.size(); i++){
		bucket[ intervals[i].first ][ bucket_size[intervals[i].first] ] = intervals[i];
		bucket_size[intervals[i].first]++;
	}

	j=0;
	for(int i=1; i<=n; i++){
		for(int k=0; k<bucket_size[i]; k++){
			intervals[j] = bucket[i][k];
			j++;
		}
		bucket_size[i] = 0;
	}
#ifdef DEBUG_SORTINTERVALS
	cout << "increasing sorted by left bound :"<<endl;
	for(unsigned i=0; i<intervals.size(); i++)
		cout << "("<<intervals[i].first<<","<<intervals[i].second<<")";
	cout << endl;
#endif
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void strong_intervals(const genom &g1, const genom &g2, int n, vector<pair<int,int> > &strong_int, bool uniq, bool trivial){
	genom g;				// the 'identified' genomes (g->id is the same as g1->g2)
	static vector<int> r,
		l;

	if((int) r.size() != n+1){
		r = vector<int>(n+1,0);
		l = vector<int>(n+1,0);
	}
		// get the strong intervals of the genomes (first generator and canonical generator)
		// identify the genomes (the strong intervals are computed with respect to the
		// second permutation (the identity) -> so identify in the reverse order)
	g = g1.identify_g(g2);
	//~ cout << "g "<<g<<endl;
	generator(g, r, l);
	canonical_generator(r,l);
	strong_intervals(r, l, strong_int, uniq, trivial);

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void strong_intervals(const genom &g1, const genom &g2, int n, vector<int> &r, vector<int> &l, vector<pair<int,int> > &strong_int, bool uniq, bool trivial ){
	genom g;				// the 'identified' genomes (g->id is the same as g1->g2)

		// get the strong intervals of the genomes (first generator and canonical generator)
		// identify the genomes (the strong intervals are computed with respect to the
		// second permutation (the identity) -> so identify in the reverse order)
	g = g1.identify_g(g2);
	//~ cout << "g "<<g<<endl;
	generator(g, r, l);
	canonical_generator(r,l);
	strong_intervals(r, l, strong_int, uniq, trivial);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void strong_intervals( const vector<genom> &genomes, int n, vector<pair<int,int> > &strong_int, bool uniq, bool trivial ){
	vector<genom> gi = genomes;
	static vector<int> r,
		l;

	if((int) r.size() != n+1){
		r = vector<int>(n+1,0);
		l = vector<int>(n+1,0);
	}

	identify( gi );
	combined_generator(gi, n, r, l);
	canonical_generator(r, l);
	strong_intervals(r, l, strong_int, uniq, trivial);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void strong_intervals( const vector<genom> &genomes, int n, vector<int> &r, vector<int> &l, vector<pair<int,int> > &strong_int, bool uniq, bool trivial ){
	vector<genom> gi = genomes;

	identify( gi );
	combined_generator(gi, n, r, l);
	canonical_generator(r, l);
	strong_intervals(r, l, strong_int, uniq, trivial);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void strong_intervals(const vector<int> &r, const vector<int> &l, vector<pair<int,int> > &strong_int, bool uniq, bool trivial){
	int n = r.size()-1;
	int si = 0;					// a counter
	stack<int> s;
	static vector<int> a,	// the sorted bounds of (i,R[i]) and (L[j],j)
		alb, 		// 1 if the bound in left
		tlb, 		// how often is a element left bound
		trb;		// how often is a element right bound

	if((int) trb.size() != n+1){
		a  = vector<int>(4*n + 1, 0);
		alb  = vector<int>(4*n + 1, 0);
		tlb = vector<int>(n+1, 1);
		trb = vector<int>(n+1, 1);
	}

		// count hof often i is a left border in (L[i],i)
		// +1 [from (i, R(i))] is the final number (done
		// in the initialisation of tlb)
 	for(unsigned i=1; i<l.size(); i++){
		tlb[ l[i] ]++;
		trb[ r[i] ]++;
	}

#ifdef DEBUG_COMMUTINGINT
	cout << "left  ";
	for(unsigned i=1; i<tlb.size(); i++)
		cout << tlb[i]<<" ";
	cout << endl;
	cout << "right ";
	for(unsigned i=1; i<trb.size(); i++)
		cout << trb[i]<<" ";
	cout << endl;
#endif

		// build vector a
	int k = 1;
	for(int i=1; i<=n; i++){
		for(int j=0; j<tlb[i]; j++){
			a[k] = i;
			alb[k] = 1;
			k++;
		}
		for(int j=0; j<trb[i]; j++){
			a[k] = i;
			k++;
		}
		tlb[i] = trb[i] = 1;
	}

#ifdef DEBUG_COMMUTINGINT
	for(unsigned i=1; i<a.size(); i++)
		cout << a[i];
	cout << endl;
	for(unsigned i=1; i<a.size(); i++)
		cout << alb[i];
	cout << endl;
#endif

	for(int i=1; i<=4*n; i++){
		if(alb[i] == 1){
			s.push(a[i]);
		}else{
			strong_int[ si ].first = s.top() - 1;
			strong_int[ si ].second = a[i] - 1;
			si ++;
			s.pop();
		}
		a[i] = alb[i] = 0;
	}
#ifdef DEBUG_COMMUTINGINT
	cout << "strong ints "<<endl;
	for(unsigned i=0; i<strong_int.size(); i++){
		cout << "("<<strong_int[i].first<<","<<strong_int[i].second<<")";
	}
	cout <<endl;
#endif

	if(uniq){
		vector<pair<int,int> >::iterator new_end;
		new_end = unique(strong_int.begin(), strong_int.end());
		strong_int.erase(new_end, strong_int.end());
	}

	if(trivial == false){
		vector<pair<int,int> >::iterator it;
		it=strong_int.begin();
		while(it != strong_int.end()){
			if(it->second == it->first || it->second - it->first == n-1 )
				it = strong_int.erase(it);
			else
				++it;
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void support(const vector<int> &r, const vector<int> &l, vector<int> &rsup, vector<int> &lsup){
	int n=r.size()-1;
	stack<int> s;

	//~ rsup = vector<int>(n+1, 0);
	//~ lsup = vector<int>(n+1, 0);

	s.push(1);
	for(int i=2; i<=n; i++){
		while(r[s.top()] < i)
			s.pop();
		rsup[i] = s.top();
		s.push(i);
	}
	while(s.size() > 0)
		s.pop();

	s.push(n);
	for(int i=n-1; i>0; i--){
		while(l[s.top()] > i)
			s.pop();
		lsup[i] = s.top();
		s.push(i);
	}
	while(s.size() > 0)
		s.pop();
#ifdef DEBUG_GENERATOR
	cout << "r support ";
	for(unsigned i=1; i< rsup.size(); i++)
		cout << rsup[i]<<" ";
	cout <<endl;
	cout << "l support ";
	for(unsigned i=1; i< lsup.size(); i++)
		cout << lsup[i]<<" ";
	cout <<endl;

#endif
}

//// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//void normalize_common(const genom &g1, const genom &g2,
//	genom &g1n, genom &g2n, int sign, hdata &hd){
//
//	vector<genom> genomes,
//		genomes_norm;
//
//	genomes.push_back(g1);
//	genomes.push_back(g2);
//
//	genomes_norm = normalize_common(genomes, sign, hd);
//
//	g1n = genomes_norm[0];
//	g2n = genomes_norm[1];
//
//	genomes.clear();
//	genomes_norm.clear();
//
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<genom> normalize_common(vector<genom> &genomes, int sign){
	int min = 0,
		max = 0;
	vector<genom> genomes_norm;
	vector<pair<int, int> > ici;
	vector<bool> border_range_elements;

#ifdef DEBUG_NORMALIZE
	for(unsigned i=0; i<genomes.size(); i++)
		cout << genomes[i]<<endl;
#endif// DEBUG_NORMALIZE

	genomes_norm = genomes;
	return genomes_norm;
	border_range_elements = vector<bool>(genomes_norm[0].size()+1, false);
		// get the irreducible common intervals for circular permutations
		// because normalize_common is only called for circular genomes -> circular
		// replacing should also be done
//	ici = find_common_intervals(genomes_norm, sign, 1, 1, 1, hd);
	irreducible_intervals(genomes_norm, genomes[0].size(), ici);
#ifdef DEBUG_NORMALIZE
	cout << "NORM: "<< ici.size() <<" irreducible common intervals "<<endl<<"NORM: ";
	for (unsigned i=0; i<ici.size(); i++)
		cout <<"["<< ici[i].first<<","<<ici[i].second<<"] ";
	cout << endl;
#endif

		// find the longest common intervall startting at the smallest (not 0) index
	for(unsigned i=0; i<ici.size(); i++){
		if( ici[i].first > 0 && min == 0 ){
			min = ici[i].first;
			max = ici[i].second + 1;
			continue;
		}

		if (min > 0 && (ici[i].second + 1) > max){
			max = ici[i].second + 1;
		}

		if(min > 0 && ici[i].first > min)
			break;
	}

#ifdef DEBUG_NORMALIZE
	cout << "NORM: border range elements "<<endl<<"NORM: ";
#endif
		// find the elements between the first two marked borders
	for(int i=min; i<max; i++){
		border_range_elements[ abs(genomes_norm[0][ i ]) ] = true;
#ifdef DEBUG_NORMALIZE
		cout << abs(genomes_norm[0][ i ])<<" ";
#endif
	}
#ifdef DEBUG_NORMALIZE
	cout << endl;
#endif

		// if there are no border elements, i.e. there is only one common interval starting
		// at index 0, then there is no need for rotating...
	if(border_range_elements.size() != 0){
			// renormalize each genome, by searching the first element from the set
			// border_range_elements and cutting there
		for(unsigned i=0; i<genomes_norm.size(); i++){
			vector<int> chromosom = genomes_norm[i].getChromosom();

			for(unsigned j=0; j<chromosom.size(); j++){
				if( border_range_elements[ abs(chromosom[j]) ] ){
					rotate(chromosom.begin(), chromosom.begin()+j , chromosom.end());
					break;
				}
			}
			genomes_norm[i].setChromosom(chromosom);
		}
	}

#ifdef DEBUG_NORMALIZE
	cout << "NORM: normalized to "<<endl<<"NORM: ";
	for(unsigned i=0; i<genomes_norm.size(); i++)
		cout << genomes_norm[i]<<endl;
#endif// DEBUG_NORMALIZE

	return genomes_norm;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int overlapping_ici(vector<pair<int,int> > &com_int,
	vector<vector<unsigned> > &overlapping_int){

	vector<int> overlap_idx; 	// stores for the common intervals to which overlap component they belong
	vector< unsigned > todo;
	int cur_idx = 0;
	unsigned cur_ci;

#ifdef DEBUG_OVERLAPPING
	cout << "irreducible common intervals "<<endl;
	for (unsigned i=0; i<com_int.size(); i++)
		cout <<i<<" "<< com_int[i].first<<","<<com_int[i].second<<" "<<endl;
#endif//DEBUG_OVERLAPPING

	overlap_idx = vector<int>(com_int.size(), -1);
	for(unsigned i=0; i<com_int.size(); i++){
		if (overlap_idx[i] != -1)
			continue;

		todo.push_back( i );
		overlap_idx[i] = cur_idx;

		while(todo.size() > 0){
			cur_ci = todo.back();
			todo.pop_back();
			for(unsigned j=0; j<com_int.size(); j++){
				if(cur_ci ==j || overlap_idx[j] == overlap_idx[cur_ci])
					continue;

				if((com_int[j].first<com_int[ cur_ci ].first && com_int[ cur_ci ].first<=com_int[j].second && com_int[j].second<com_int[ cur_ci ].second )||
					(com_int[ cur_ci ].first<com_int[j].first && com_int[j].first<=com_int[ cur_ci ].second && com_int[ cur_ci ].second<com_int[j].second )){
						todo.push_back( j );
						overlap_idx[j] = overlap_idx[ cur_ci ];
				}
			}
		}

		cur_idx++;
	}
	overlapping_int = vector<vector<unsigned> >( cur_idx ) ;
	for(unsigned i=0; i<overlap_idx.size(); i++){
		overlapping_int[ overlap_idx[i] ].push_back( i );
	}


#ifdef DEBUG_OVERLAPPING
	cout << "overlapping_int"<<endl;
	for(unsigned i=0; i<overlapping_int.size(); i++){
		cout << i << " -> " ;
		for(unsigned j=0; j<overlapping_int[i].size(); j++)
			cout << "["<< com_int[ overlapping_int[i][j] ].first<<","<<com_int[ overlapping_int[i][j] ].second <<"] ";
		cout << endl;
	}
#endif//DEBUG_OVERLAPPING


	return 1;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< vector<int> > valid_adjacencies_common(vector< genom > &genomes,
		int sig, int circular,
		vector<vector<int> > &element_range_map,
		vector<vector<int> > &range_size,
		vector<int> &range_max){

	vector<genom> genomes_norm;
	vector<pair<int,int> > com_int;	// the common intervals of the genomes
	vector<vector<unsigned> > overlapping_int;	// indices of overlapping intervals

	vector<bool> range_borders;	// 1range if one of the overlapping interval starts or ends
	int start_idx = 0, 	// start and end of ranges
		end_idx = 0;
	unsigned n = 0;					// length of the genomes
	vector< vector<int> > genesets;	// the sets of genes in the ranges of overlapping common intervals
	pair<unsigned, unsigned> e;		// a matching edge
	vector<int> connected_int_idx;	// store where connected (overlapping / including) common intervals are
									// -1: no connected common interval;
									// all other connected common intervals get a unique number
	int max_connected_int_idx = -1,		// maximal number given for a connected common intervals
		curr_connected_int_idx = -1;	// current number for a connected common intervals
	vector< vector<int> > va;			// the valid adjacencies matrix
	int range_nr = 0;

		// if the genomes are circular -> normalize them correctly
		// then the restrictions (va, rangemap) are computed in a way that caprara's
		// median solver produces one circular representative which preserves
		// the common intervals
	if(circular != 0){
		genomes_norm = normalize_common(genomes, sig);
	}else{
		genomes_norm = genomes;
	}

		// some datastructures must be initialised
	n = genomes_norm[0].size();
	va = vector< vector<int> >(2*n+2, vector<int>(2*n+2, 1));
	element_range_map.clear();
	range_size.clear();
	range_max.clear();
		// get the common intervals of the genomes
	common_intervals(genomes_norm, circular, 0, sig, com_int);
//	com_int = find_common_intervals(genomes_norm, sig, 1, circular, 1, hd);
#ifdef DEBUG
	for(unsigned i=0; i<genomes_norm.size(); i++)
		cout << genomes_norm[i]<<endl;
	cout << "irred common intervals"<<endl;
	for(unsigned i=0; i<com_int.size(); i++)
		cout << "["<<com_int[i].first<<","<<com_int[i].second<<"] ";
	cout << endl;
#endif

		// find all overlapping common intervals (which are sorted in ascending order)
	overlapping_ici(com_int, overlapping_int);

	range_borders = vector<bool>(n+1, false);

#ifdef DEBUG
	cout << "overlapping_int"<<endl;
	for(unsigned i=0; i<overlapping_int.size(); i++){
		cout << i << " -> " ;
		for(unsigned j=0; j<overlapping_int[i].size(); j++)
			cout << overlapping_int[i][j]<<" ";
		cout << endl;
	}
#endif

		// for each set of overlapping common intervals
		// - mark the start and end indices in the range_borders vector
		// 	and get the start and the end of the component
		// - get the sets of genes between the marked indices genesets[1..] and the
		//	set of genes outside the component genesets[0]
		// - then forbid all ajacencies +a+b, -a-b, +b+a, -b-a with a\in genesets[0]
		//		and b \in genesets[2..size-1]
		//	then forbid all adjacencies +a+b, -a-b, +b+a, -b-a with a\in genesets[j]
		//		and b \in genesets[k] and j+1<k
		//	(note if sig>0 (true) then -ab, a-b, -ba, b-a are forbidden in the initialisation
		//		of va; if sig==0 (false) then this has to be done here)
		//	(note that +a+b and -b-a imply the matching edges (x,y) and (y,x))
		// - setup the rangemaps for each overlap component (also the trivial overlaps
		//	consisting of only one intervals)
	for(unsigned i=0; i<overlapping_int.size(); i++){
		//~ if(overlapping_int[i].size() < 2)
			//~ continue;

		genesets.clear();
			// - mark the start and end indices in the range_borders vector
			// 	and get the start and the end of the component
		range_borders.assign(n+1, false);
		start_idx = com_int[ overlapping_int[i][0] ].first;
		end_idx = 0;

		for(unsigned j=0; j<overlapping_int[i].size(); j++){
		range_borders[ com_int[ overlapping_int[i][j] ].first ] = true;
			range_borders[ com_int[ overlapping_int[i][j] ].second+1 ] = true;
			if( com_int[ overlapping_int[i][j] ].second+1 > end_idx)
				end_idx = com_int[ overlapping_int[i][j] ].second+1;
		}
#ifdef DEBUG
		cout << "range_borders"<<endl;
		for (unsigned j=0; j<range_borders.size(); j++){
			cout << range_borders[j]<<" ";
		}
		cout <<endl;
#endif
			// - get the sets of genes between the marked indices genesets[1..] and the
			//	set of genes outside the component genesets[0]
		genesets.push_back(vector<int>(1,0));	//  element 0 -> outside elements
		for(unsigned j=0; j<range_borders.size() - 1; j++){
			if((int)j < start_idx || (int)j > end_idx - 1){
				genesets[0].push_back(genomes_norm[0][j]);
			}else{
				if(range_borders[j]){
					genesets.push_back(vector<int>());
				}
				genesets.back().push_back( genomes_norm[0][j] );
			}
		}
		genesets[0].push_back(n+1);	//  element n+1 -> outside elements
#ifdef DEBUG
		cout << "genesets"<<endl;
		for (unsigned j=0; j<genesets.size(); j++){
			for(unsigned k=0; k<genesets[j].size(); k++)
				cout << genesets[j][k] <<" ";
			cout <<endl;
		}
#endif
			// - then forbid all ajacencies ab, ba with a\in genesets[0] and b \in genesets[2..size-1]
			// 	(only needed if there is really an overlap)
		if(overlapping_int[i].size() > 1){
			for(unsigned j=0; j<genesets[0].size(); j++){
				for(unsigned k=2; k<genesets.size()-1; k++){
					for(unsigned l=0; l<genesets[k].size(); l++){
						e = matching_edge(genesets[0][j], genesets[k][l], n+1);	// +a+b & -b-a
						va[e.first][e.second] = va[e.second][e.first] = 0;

						e = matching_edge(genesets[k][l], genesets[0][j], n+1);	// +b+a && -a-b
						va[e.first][e.second] = va[e.second][e.first] = 0;

						//~ if(sig==0){
							//~ e = matching_edge(-1*genesets[0][j], genesets[k][l], n+1);	// -a+b & -ba
							//~ cout << "a exclude "<<-1*genesets[0][j]<<","<<genesets[k][l]<<" -> "<< e.first<<","<<e.second<<endl;
							//~ va[e.first][e.second] = va[e.second][e.first] = 0;

							//~ e = matching_edge(genesets[k][l], -1*genesets[0][j], n+1);	// +b-a && a-b
							//~ cout << "b exclude "<<genesets[k][l]<<","<<-1*genesets[0][j]<<" -> "<< e.first<<","<<e.second<<endl;
							//~ va[e.first][e.second] = va[e.second][e.first] = 0;
						//~ }
					}
				}
			}
				//	then forbid all adjacencies +a+b, -a-b, +b+a, -b-a with a\in genesets[j] and b \in genesets[k] and j+1<k
			for(unsigned j=1; j<genesets.size()-1; j++){
				for(unsigned k=j+2; k<genesets.size(); k++){
					for(unsigned l=0; l<genesets[j].size(); l++){
						for(unsigned m=0; m<genesets[k].size(); m++){
							e = matching_edge(genesets[j][l], genesets[k][m], n+1);	// +a+b & -b-a
							va[e.first][e.second] = va[e.second][e.first] = 0;

							e = matching_edge(genesets[k][m], genesets[j][l], n+1);	// +b+a && -a-b
							va[e.first][e.second] = va[e.second][e.first] = 0;

							if(sig==0){
								e = matching_edge(-1*genesets[j][l], genesets[k][m], n+1);	// -a+b & -b+a
								//~ cout << "c exclude "<<-1*(int)genesets[j][l]<<","<<genesets[k][m]<<" -> "<< e.first<<","<<e.second<<endl;
								va[e.first][e.second] = va[e.second][e.first] = 0;


								e = matching_edge(genesets[k][m], -1*genesets[j][l], n+1);	// +b-a && a-b
								//~ cout << "d exclude "<<genesets[k][m]<<","<<-1*(int)genesets[j][l]<<" -> "<< e.first<<","<<e.second<<endl;
								va[e.first][e.second] = va[e.second][e.first] = 0;
							}
						}
					}
				}
			}
		}

			// element 0 and n+1 -> -1
		element_range_map.push_back( vector<int>(n+2, -1) );

		range_size.push_back( vector<int>(1, 0) );
		range_nr = 0;
		for(unsigned j=1; j<genesets.size(); j++){

			for(unsigned k=0; k<genesets[j].size(); k++){
				(element_range_map.back())[ abs(genesets[j][k]) ] = range_nr;
				(range_size.back())[ range_nr ]++;
			}
			range_size.back().push_back(0);
			range_nr++;
		}
		range_size.back().pop_back();
		range_max.push_back( range_nr-1 );
	}

		// if signed variant -> all adjacencies of elements a and b, where a and b are in the same
		// connected component and a and b have different signs can be excluded
	if(sig>0){
			// get the regions of connected common intervals
			//
			// check for all common intervals (ascending sorted) if they start in an
			// area already of a connected common interval; if so -> the interval from
			// start to end of the common interval belongs to this connected component;
			// else a new connected component has to be started and the interval from
			// start to end of the common interval belongs to this new component

			// for circular genomes the irreducible genomes are needed with unreplaced
			// common intervals .. because the sign restriction only hold for the original ci
		if(circular > 0){
			com_int.clear();
//			com_int = find_common_intervals(genomes_norm, sig, 1, circular, 0, hd);

			common_intervals(genomes_norm, circular, 0, sig, com_int);
		}

		connected_int_idx = vector<int>(n, -1);
		for(unsigned i=0; i<com_int.size(); i++){
				// the check if the current interval is connected has to respect the circularity
				// for linear the 1st and 2nd test are enough
				// circular may need the other tests too
			if(connected_int_idx[ com_int[i].first ] == -1 && connected_int_idx[ com_int[i].second ] == -1){
				max_connected_int_idx++;
				curr_connected_int_idx = max_connected_int_idx;
			}else if(connected_int_idx[ com_int[i].first ] > -1 && connected_int_idx[ com_int[i].second ] == -1){
				curr_connected_int_idx = connected_int_idx[ com_int[i].first ];
			}else if(connected_int_idx[ com_int[i].first ] == -1 && connected_int_idx[ com_int[i].second ] > -1){
				curr_connected_int_idx = connected_int_idx[ com_int[i].second ];
			}else{
				curr_connected_int_idx = connected_int_idx[ com_int[i].first ];
				for(unsigned j=0; j<connected_int_idx.size(); j++){
					if(connected_int_idx[j] == connected_int_idx[ com_int[i].second ])
						connected_int_idx[j]  = curr_connected_int_idx;
				}
			}

			for(unsigned j=com_int[i].first; j<=(unsigned)com_int[i].second; j++){
				j %= n;	// for circular genomes with border crossing intervals
				connected_int_idx[ j ] = curr_connected_int_idx;
			}
		}

#ifdef DEBUG
		cout << "connected intervals"<<endl;
		for(unsigned i=0; i<connected_int_idx.size(); i++){
			cout << connected_int_idx[i] << " ";
		}
		cout << endl;
#endif

		start_idx = end_idx = 0;
		curr_connected_int_idx = -1;
		for(unsigned i=0; i <= connected_int_idx.size(); i++){
				// end of a connected component
			if((i==connected_int_idx.size() && curr_connected_int_idx != -1) ||
					curr_connected_int_idx != connected_int_idx[i]){

				end_idx = i;
				if(curr_connected_int_idx != -1 && start_idx < end_idx){
					for(unsigned j=start_idx; j<(unsigned)end_idx; j++){
						for(unsigned k=j+1; k<(unsigned)end_idx; k++){
							// a-b; b-a;
							e = matching_edge(genomes_norm[0][j], -1*genomes_norm[0][k], n+1);
							va[e.first][e.second] = va[e.second][e.first] = 0;

							//-a,b;  -ba
							e = matching_edge(-1*genomes_norm[0][j], genomes_norm[0][k], n+1);
							va[e.first][e.second] = va[e.second][e.first] = 0;

						}
					}
				}
				start_idx = i;
				if(i<connected_int_idx.size())
					curr_connected_int_idx = connected_int_idx[i];
			}
		}
	}


		// forbid nonsense adjacencied:
		// - adjacencies of the element with itself or its reverse are nonsense -> reset the diagonal to 0
		// - adjacencies (x,0)
		// - adjacencies (n+1, x)
	for(unsigned i=0; i<va.size()-1; i++){
		va[i][i] = 0;
		if(i>0 && i<va.size() && i%2){
			va[i][i+1] = va[i+1][i] = 0;
		}
		va[i][0] = 0;
	}
	va.back().assign(2*n+2, 0);

	if(circular != 0){
		for(unsigned i=0; i<2*n+2; i++){
			va[0][i] = 0;
		}
		va[0][1] = 1;
	}

#ifdef DEBUG
	cout << "valid adjacencies"<<endl;
	for(unsigned i=0; i<va.size(); i++){
		for(unsigned j=0; j<va[i].size(); j++){
			cout << va[i][j]<<" ";
		}
		cout << endl;
	}

	cout << "element_range_maps"<<endl;
	for (unsigned i=0; i<element_range_map.size(); i++){
		for(unsigned j=0; j<genomes[0].size(); j++){
		//~ for(unsigned j=0; j<element_range_map[i].size(); j++){
			cout.width(2);
			cout <<element_range_map[i][ abs(genomes[0][j]) ]<<" ";
		}
		cout << endl;
	}
	cout << "range_sizes"<<endl;
	for(unsigned i=0; i<range_size.size(); i++){
		cout << "........................"<<endl;
		for(unsigned j=0; j<range_size[i].size(); j++)
			cout << j << " : "<<range_size[i][j]<<endl;
	}
	cout << "range max"<<endl;

	for(unsigned i=0; i<range_max.size(); i++){
		cout << i<<" - "<<range_max[i]<<endl;
	}
#endif
	return va;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*extern "C" int ** valid_adjacencies_common(
		struct genome_struct **passgenome,
		int ngenomes,
		int ngenes,
		int sig, int circular,
		int ***element_range_map,
		int ***range_size,
		int **range_max,
		int *overlap_cmp_cnt){

	int **va;

	vector<genom> genomes;
	vector<vector<int> > va_vec,
		element_range_map_vec,
		range_size_vec;
	vector<int> range_max_vec;

	for(unsigned i=0; i<(unsigned)ngenomes; i++){
		genomes.push_back(genom(passgenome[i]->genes, ngenes, circular));
	}

	va_vec = valid_adjacencies_common(genomes, sig, circular,
		element_range_map_vec,
		range_size_vec,
		range_max_vec);

		// make va array from vector
	va = (int**)malloc(va_vec.size() * sizeof(int *));
	for(unsigned i=0; i<va_vec.size(); i++){
		va[i] = (int*) malloc(va_vec[i].size() * sizeof(int));
		for(unsigned j=0; j<va_vec[i].size(); j++){
			va[i][j] = va_vec[i][j];
		}
	}

		// make c-arary from element_range..
	*overlap_cmp_cnt = element_range_map_vec.size();
	*element_range_map = (int**) malloc(element_range_map_vec.size() * sizeof(int *));
	*range_size = (int**) malloc(element_range_map_vec.size() * sizeof(int *));
	for(unsigned i=0; i<element_range_map_vec.size(); i++){
		(*element_range_map)[i] = (int*) malloc(element_range_map_vec[i].size() * sizeof(int));
		(*range_size)[i] = (int*) malloc(range_size_vec[i].size() * sizeof(int));

		for(unsigned j=0; j<element_range_map_vec[i].size(); j++){
			(*element_range_map)[i][j] = element_range_map_vec[i][j];
		}
		for(unsigned j=0; j<range_size_vec[i].size(); j++){
			(*range_size)[i][j] = range_size_vec[i][j];
		}
	}

	*range_max = (int *) malloc(range_max_vec.size() * sizeof(int));
	for( unsigned i=0; i<range_max_vec.size(); i++ ){
		(*range_max)[i] = range_max_vec[i];
	}
	return va;
}*/
