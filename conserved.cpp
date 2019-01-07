#include <algorithm>
#include <iostream>
#include <limits>
#include <stack>

#include "conserved.hpp"

using namespace std;

//~ #define DEBUG_CONSERVED_INTERVALS
//~ #define DEBUG_VALID_ADJACENCIES



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void init_pqnode(pqnode* parent, pqnode** n, int type, pair<int,int> interval, int orientation){

	*n = new(pqnode);
	(*n)->parent = parent;

	if(parent != NULL){
		parent->children.push_back(*n);
	}


	(*n)->i = interval;
	(*n)->type = type;
	(*n)->orientation = orientation;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void pqtree(const vector<pair<int, int> > &comp, const vector<int> &comp_or, int n, pqnode **pqroot){
	pqnode *p, *q;
	vector<int> start,
		end;

	start = vector<int>(n+2, -1);
	end = vector<int>(n+2, -1);

		// which intervals start/end at idx i
	for(unsigned i=0; i<comp.size(); i++){
		start[ comp[i].first ] = i;
		end[ comp[i].second ] = i;
	}

	init_pqnode( NULL, pqroot, Q, make_pair(0,0), PQUNORIENTED);
	p = *pqroot;
	init_pqnode( *pqroot, &q, P, comp[ start[0] ], comp_or[ start[0] ] );

	for(int i=1; i<n; i++){
		if(start[i] >= 0){
			if(end[i] < 0){
				init_pqnode( q, &p, Q, make_pair(0,0), PQUNORIENTED);
			}
			init_pqnode( p, &q, P, comp[ start[i] ], comp_or[ start[i] ]);
		}else if(end[i] >= 0){
			q = p->parent;
			p = q->parent;
		}
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void pqtree_branches(pqnode *pqroot, int len, vector<int> &branch_lens){

	if(pqroot->type == P &&  pqroot->orientation == PQUNORIENTED){
		len++;
	}

		// reached a leaf
	if(pqroot->children.size() == 0 && len > 0){
		branch_lens.push_back(len);
	}

		// reset the length at nodes with degree >= 3 (i.e. >= 2 leaves)
	if(pqroot->children.size() >= 2){
		len = 0;
	}

	for(unsigned i=0; i<pqroot->children.size(); i++){
		pqtree_branches(pqroot->children[i], len, branch_lens);
	}


}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool pqtree_defoliate(pqnode **p){
	bool d;

	for(int i=0; i<(int)(*p)->children.size(); i++){
		d = pqtree_defoliate(&(*p)->children[i]);
		if( d ){
			delete( (*p)->children[i] );
			(*p)->children.erase( (*p)->children.begin() + i );
			i--;
		}
	}

			// handling of the root node
	if((*p)->parent == NULL ){
		pqtree_reroot(p);
	}

		// mark oriented leaves and square nodes as 'to delete'
	if( (*p)->children.size() == 0 && ((*p)->orientation == PQORIENTED || (*p)->type == Q)){
//		cout << "remove "<<(*p)->i.first <<".."<<(*p)->i.second<<endl;
		return true;
	}else
		return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void pqtree_free(pqnode *p){
	for( unsigned i=0; i<p->children.size(); i++ ){
		pqtree_free( p->children[i] );
		delete( p->children[i] );
	}
	p->children.clear();

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void pqtree_print(pqnode *p){

	if(p == NULL){
		cout << " is NULL"<<endl;
		return;
	}

	cout << "("<< p->i.first<<","<<p->i.second<<") ";
	switch (p->type){
		case Q: {cout << "Q "; break;}
		case P: {cout << "P "; break;}
	}
	switch( p->orientation ){
		case PQORIENTED:{cout << "  ORIENTED ";break;}
		case PQUNORIENTED:{cout << "UNORIENTED ";break;}
	}
	cout << "-> ";
	for( unsigned i=0; i<p->children.size(); i++ ){
		cout << "("<<p->children[i]->i.first<<","<<p->children[i]->i.second<<") ";
	}
	cout << endl;
	for( unsigned i=0; i<p->children.size(); i++ ){
		pqtree_print( p->children[i] );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void pqtree_reroot(pqnode **root){
	pqnode *newroot;

		// T consists of only 1 node -> return
	if((*root)->children.size() == 0){
		return;
	}

	newroot = *root;
	while( newroot->children.size() == 1 ){
		newroot = newroot->children[0];
	}

		// if the newroot pointer points to a leaf -> the tree is just a line
		// than take the next square node 'under' the old root as new root
		// (its less work, compared to the strategy which takes the last found square node)
		// if there is no square node under the root -> leave it as it is
	if(newroot->children.size() == 0){
		if((*root)->children[0]->children.size() > 0 )
			newroot = (*root)->children[0]->children[0];
		else
			return;
	}

	while( *root != newroot ){
		(*root)->parent = (*root)->children[0];
		(*root)->children[0]->children.push_back( *root );
		(*root)->children.clear();
		(*root) = (*root)->parent;
		(*root)->parent = NULL;

//		pqtree_print(*root);
//		cout << "--------"<<endl;
	}
}



int changes_conserved_intervals(
		const vector<genom> &genomes, const genom &g,
		vector<pair<unsigned, unsigned> > &ci_pre,
		vector<pair<unsigned, unsigned> > &ci_dif,
		int circular){

	int difference;
	vector<pair<unsigned, unsigned> > ci_post;
	vector<genom> genomes_post;

	ci_dif.clear();
		// get the conserved intervals without g
	if(ci_pre.size() == 0){
		ci_pre = conserved_intervals(genomes, circular);
	}
		// get the conserved intervals with g
	genomes_post = genomes;
	genomes_post.push_back(g);
	ci_post = conserved_intervals(genomes_post, circular);
	genomes_post.clear();

	difference = (int)ci_pre.size() - (int)ci_post.size();

	if(difference > 0)
		set_difference(ci_pre.begin(), ci_pre.end(),
					ci_post.begin(), ci_post.end(), back_inserter(ci_dif));

	return difference;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*extern "C" int changes_conserved_intervals(
	struct genome_struct **genomes_ptr, int m,
	int *g_ptr, int n,
	unsigned ***ci_pre_ptr, unsigned *ci_pre_size,
	unsigned ***ci_dif_ptr, unsigned *ci_dif_size,
	int circular){

	genom g;
	int difference;
	vector<genom> genomes;
	vector<pair<unsigned, unsigned> > ci_dif,
		ci_pre;

	g = genom(g_ptr, n, 0);
	for(int i=0; i<m; i++){
		genomes.push_back( genom(genomes_ptr[i]->genes, n, 0) );
	}

	if(*ci_pre_size != 0){
		for(unsigned i=0; i<*ci_pre_size; i++){
			ci_pre.push_back(pair<unsigned, unsigned>((*ci_pre_ptr)[0][i], (*ci_pre_ptr)[1][i]) );
		}
	}

	if(*ci_dif_size != 0){
		free((*ci_dif_ptr)[0]);
		free((*ci_dif_ptr)[1]);
		free((*ci_dif_ptr));
		*ci_dif_size = 0;
	}

	difference = changes_conserved_intervals(genomes, g, ci_pre, ci_dif, circular);

	if(*ci_pre_size == 0){
		*ci_pre_size = ci_pre.size();
		*ci_pre_ptr = (unsigned **) malloc( 2 * sizeof(unsigned *));
		(*ci_pre_ptr)[0] = (unsigned *) malloc( *ci_pre_size * sizeof(unsigned));
		(*ci_pre_ptr)[1] = (unsigned *) malloc( *ci_pre_size * sizeof(unsigned));


		for(unsigned i=0; i<*ci_pre_size; i++){
			(*ci_pre_ptr)[0][i] = ci_pre[i].first;
			(*ci_pre_ptr)[1][i] = ci_pre[i].second;
		}
	}

	if(ci_dif.size() > 0){
		*ci_dif_size = ci_dif.size();

		*ci_dif_ptr = (unsigned **) malloc( 2 * sizeof(unsigned *));
		(*ci_dif_ptr)[0] = (unsigned *) malloc( ci_dif.size() * sizeof(unsigned));
		(*ci_dif_ptr)[1] = (unsigned *) malloc( ci_dif.size() * sizeof(unsigned));


		for(unsigned i=0; i<*ci_dif_size; i++){
			(*ci_dif_ptr)[0][i] = ci_dif[i].first;
			(*ci_dif_ptr)[1][i] = ci_dif[i].second;
		}
	}

	return difference;
}*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<pair<unsigned, unsigned> > conserved_intervals(
	const vector< vector<int> > &mc, int circular, int n){

	unsigned last_index;
	vector<pair<unsigned, unsigned> > ci;

#ifdef DEBUG_CONSERVED_INTERVALS
	cout << "maximal chains"<<endl;
	for(unsigned i=0; i<mc.size(); i++){
		for(unsigned j=0; j<mc[i].size(); j++){
			cout << mc[i][j]<<" ";
		}
		cout << endl;
	}
#endif

		// enumerate all conserved intervalls by enumaration of all pairs of maximal
		// chain elements
		// for linear genomes this are simply all (directed) pairs (i.e. ab but not ba)
		// for circular genomes there is an exception :
		// - pairs with (n+1,*) or (*, n+1) are not wanted
	for(unsigned i=0; i<mc.size(); i++){
		last_index = mc[i].size();
		if(circular == 1 && i == 0){
			last_index--;
		}
		for(unsigned j=0; j<last_index; j++){
			for(unsigned k=j+1; k<last_index; k++){
				if(circular == 1 && mc[i][j] != n+1){
					ci.push_back(pair<unsigned, unsigned>( mc[i][k], mc[i][j]));
					//~ cout << ci.back().first << ","<<ci.back().second<<endl;
					ci.push_back(pair<unsigned, unsigned>( mc[i][j], mc[i][k]));
					//~ cout << ci.back().first << ","<<ci.back().second<<endl;
				}else if(circular == 0){
					ci.push_back(pair<unsigned, unsigned>( mc[i][k], mc[i][j]));
				}
				//~ cout << "----------"<<endl;
			}
		}
	}
	return ci;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<pair<unsigned, unsigned> > conserved_intervals(
		const genom &g1, const genom &g2, int circular){

	vector<pair<unsigned, unsigned> > ci;
	vector< vector<int> > mc;

	mc = getMaximalChains(g1, g2);

	return conserved_intervals(mc, circular, g1.size());
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<pair<unsigned, unsigned> > conserved_intervals(
		const vector<genom> &genomes, int circular){

	vector<pair<unsigned, unsigned> > ci;
	vector< vector<int> > mc;

	mc = getMaximalChains(genomes);

	return conserved_intervals(mc, circular, genomes[0].size());
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



vector< pair<int, int> > getIrreducibleConservedIntervals(const genom &source,
		const genom &target,vector<int> &cior){
	int n = source.size()+1,
		s,r;
	vector< pair<int,int> > ci;	// the irreducible conserved intervals
	vector<int> pi;			// temp memory for the permutation

		// stacks ... for positive conserved intervals
	stack<int> M1;				// bigger elements left from the element
	stack<int> S1;				// next possible start index
	vector<int> M( n+1);		// the next smaller element at the right

	stack<int> M2;				// smaller elements left from the element
	stack<int> S2;				// next possible start index
	vector<int> m( n+1);			// the next smaller element at the left
	static bool initialised = false;
	static vector<int> o;		// the orientation of the points

	if(!initialised){
		o = vector<int>(n+1,0);
		initialised = true;
	}

	pi = target.identify( source );	// use a little bit of extra memory for the chromosome

		// init the point orientations
	for(int i=0; i<n; i++){
		if(pi[i] >= 0 && pi[i+1] >= 0)
			o[i] = 1;
		else if(pi[i] < 0 && pi[i+1] < 0)
			o[i] = -1;
		else
			o[i] = 0;
	}

	S1.push(0);						// init the stacks ...
	M1.push(n);
	M[0] = n;

	S2.push(0);
	M2.push(0);
	m[0]=0;

	for (int i=1; i<=n; i++){
			// get the positive conserved intervals [a,b]
		while(M1.top() < abs(pi[i]))				// get the M_i
			M1.pop();
		M[i] = M1.top();
		M1.push(abs(pi[i]));

		s = S1.top();
		while(abs(pi[i]) < pi[s] || abs(pi[i]) > M[s])	{ // exclude impossible startpoints
			S1.pop();
			r = S1.top();
			if(o[s] != o[r])
				o[r] = 0;
			s = r;
		}
		if(pi[i]>0 && i-s == pi[i] - pi[s] && M[i] == M[s]){	// test if they are realy startpoints
			ci.push_back(pair<int,int>(s,i));
			if(s+1 != i){			// all had the same sign -> unoriented
				if(o[s] != 0)
					cior.push_back(0);
				else						// there were sign differences -> oriented
					cior.push_back(1);
			}else{
				cior.push_back(1);		// trivial components are oriented
			}
		}
		if(pi[i]>0)		// possible startpoints for the next round
			S1.push(i);

			// get the negative conserved intervals [-b,-a]
		while(M2.top() > abs(pi[i]))				// get the m_i
			M2.pop();
		m[i] = M2.top();
		M2.push(abs(pi[i]));

		s = S2.top();
		while((abs(pi[i]) > abs(pi[s]) || abs(pi[i]) < m[s]) && s > 0){
			S2.pop();
			r = S2.top();
			if(o[s] != o[r])
				o[r] = 0;
			s = r;
		}
		if(pi[i]<0 && i-s == pi[i]-pi[s] && m[i] == m[s]){
			ci.push_back(pair<int,int>(s,i));
			if(s+1 != i){
				if(o[s] != 0)		// all had the same sign -> unoriented
					cior.push_back(0);
				else						// there were sign differences -> oriented
					cior.push_back(1);
			}else{						// trivial components are oriented
				cior.push_back(1);
			}
		}
		if(pi[i]<0)		// possible startpoints for the next round
			S2.push(i);
	}

	return ci;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< vector<int> > getMaximalChains(const genom &source, const genom &target){
	vector<int> cior;
	vector< pair<int,int> > ci;

	ci = getIrreducibleConservedIntervals(source, target, cior);

	return getMaximalChains(source.size(), ci);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< vector<int> > getMaximalChains(const vector<genom> &targets) {

	vector< vector<int> > I1,
		I2;
	vector<bool> I1_starts(targets[0].size()+1, false),	// starts: if there is an irreducible interval that starts
		I2_starts(targets[0].size()+1, false);			// at at index i -> starts[i] = true
	vector<unsigned> I1_ends(targets[0].size()+2, std::numeric_limits< unsigned >::max()),		// ends: if there is an irreducible interval [x,i] ending at
		I2_ends(targets[0].size()+2,std::numeric_limits< unsigned >::max()),				// index i -> ends[i] gives x
		I1_chainnr(targets[0].size()+2,std::numeric_limits< unsigned >::max()),			// chainnr: gives the index of the maximal chain in which
		I2_chainnr(targets[0].size()+2,std::numeric_limits< unsigned >::max());			// index i is contained in
	vector< pair<int,int> > ci;

	vector<unsigned> S;							// stack S
	unsigned s;									// top of S

	for(unsigned t=2; t<targets.size(); t++){

			// init the maximal chains I1 and I2 ;
			// I1 are in the first iteration the maximal chains of 0 and 1,
			// in all furter iterations it is the result of the last iteration; I2 are the maximal cahins of 0 and i
		if(t==2){
			I1 = getMaximalChains(targets[0], targets[1]);
		}else{
			I1 = getMaximalChains(targets[0].size(), ci);
			ci.clear();
		}
		I2 = getMaximalChains(targets[0], targets[t]);

			// init the helping datastructures
		I1_starts.assign(targets[0].size()+1, false);
		I2_starts.assign(targets[0].size()+1, false);
		I1_ends.assign(targets[0].size()+2,std::numeric_limits< unsigned >::max());
		I2_ends.assign(targets[0].size()+2,std::numeric_limits< unsigned >::max());
		I1_chainnr.assign(targets[0].size()+2,std::numeric_limits< unsigned >::max());
		I2_chainnr.assign(targets[0].size()+2,std::numeric_limits< unsigned >::max());

		for(unsigned i=0; i<I1.size(); i++){
			for(unsigned j=0; j<I1[i].size(); j++){
				if(j>0)
					I1_starts[I1[i][j]] = true;

				if(j<I1[i].size()-1)
					I1_ends[I1[i][j]] = I1[i][j+1];
				I1_chainnr[I1[i][j]] = i;
			}
		}
		for(unsigned i=0; i<I2.size(); i++){
			for(unsigned j=0; j<I2[i].size(); j++){
				if(j>0)
					I2_starts[I2[i][j]] = true;
				if(j<I2[i].size()-1)
					I2_ends[I2[i][j]] = I2[i][j+1];
				I2_chainnr[I2[i][j]] = i;
			}
		}

		S.clear();
		S.push_back(0);

		for(unsigned i=1; i<targets[0].size()+2; i++){

			if(I1_ends[i] != std::numeric_limits< unsigned >::max()){
				while(S.back() > I1_ends[i]){
					S.pop_back();
				}
			}

			if(I2_ends[i] != std::numeric_limits< unsigned >::max()){
				while(S.back() > I2_ends[i]){
					S.pop_back();
				}
			}
			s = S.back();

			if(I1_chainnr[s] == I1_chainnr[i] && I2_chainnr[s] == I2_chainnr[i]){
				S.pop_back();
				ci.push_back(pair<int,int>(s,i));
			}

			if(I1_starts[i] && I2_starts[i]){
				S.push_back(i);
			}
		}

	}

		// return the maximal chains of the irreducible intervals determined in the last iteration
	return getMaximalChains(targets[0].size(), ci);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< vector<int> > getMaximalChains(unsigned n, const vector< pair<int,int> > &ici){
	vector< vector<int> > maxChains;		// the maximal chains to compute
	vector<int> maxChainsStart(n+2, -1);
	int maxChainsCnt=0;

		// get the maximal chains; store it in reverse direction (lexical order)
		// e.g.: 10 9 5 4 3 2 1 0 then 8 7 6
	for (int i=ici.size()-1; i>=0; i--){
			// end of the irreducible conserved interval was the start of an other -> append
		if(maxChainsStart[ici[i].second]>=0){
			maxChains[maxChainsStart[ici[i].second]].push_back(ici[i].first);
			maxChainsStart[ici[i].first] = maxChainsStart[ici[i].second];
		}else{	// begin new maximal chain
			maxChains.push_back(vector<int>());
			maxChains[maxChainsCnt].push_back(ici[i].second);
			maxChains[maxChainsCnt].push_back(ici[i].first);
			maxChainsStart[ici[i].first] = maxChainsCnt;
			maxChainsCnt++;
		}
	}

	return maxChains;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< pair<int,int> > getPreservingReversals(const vector<genom> &G){
	vector< vector<int> > mc;		// maximal chains

	mc = getMaximalChains(G);
	return getPreservingReversalsFromMaxChains(mc);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< pair<int,int> > getPreservingReversals(const genom &source, const genom &target){
	vector< vector<int> > mc;		// maximal chains

	mc = getMaximalChains(source, target);
	return getPreservingReversalsFromMaxChains(mc);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< pair<int,int> > getPreservingReversalsFromMaxChains(const vector< vector<int> > &mc){
	vector< pair<int,int> > pr;	// the preserving reversals

	vector< vector<pair<int,int> > > groups;
	vector<pair<int,int> > groupsStack;
	vector<pair<int,int> > groupsBorder;

	int groupCnt=-1;

		// group the maximal chains together, if they are nonoverlaping
		// neighbours in the same nesting level e.g.
		//          --
		//  --    ------
		// ---------------
		// 1 2 3 4 5 6 7 8
		// will result in 4 groups

	// each maximal chains in mc is sorted backwards ... so it could look like this:
	// 	chain1: 11,9,7,1,0,
	// 	chain2: 6,2,
	// 	chain3: 4,3,
	for (unsigned i = 0; i<mc.size(); i++){		// iterate through all maximal chains
		if(groupsStack.size()){							// if there is a stack
				// search the right hole (read comment below)
				// the endindex of this group need to be bigger than the startindex of the right hole
				// so pop all false holes and insert in the right one
			while(mc[i][ 0 ] <  groupsStack[groupsStack.size()-1].first){	// if out of the border this belongs to the next group
				groupsStack.pop_back();
			}
			groups[groupsStack[groupsStack.size()-1].second].push_back(pair<int,int>(mc[i][0], mc[i][ mc[i].size()-1] ));	// store it
		}
			// search the chain for 'holes' (nested 'sub-chains') (index difference greater than 1)
			// in the example: chain 1 has the 'holes' 11-9; 9-7; 7-1 ...
			// - each hole could include other chains so start a new group and store the borders
			// 	of this group (the borders of this 'hole') -> corresponding groups and groupsBorder has the same index
			// - because the groups included in this hole come in later max.chains (backward sorted) -> need to
			// 	store the startindex (for assigning to the right hole) and the index of the hole in the group array
		for(unsigned j = mc[i].size()-1; j>0; j--){
			if(mc[i][j-1] - mc[i][j] > 1){
				groups.push_back(vector<pair<int,int> >());					// append a new group
				groupsBorder.push_back(pair<int,int>(mc[i][j-1], mc[i][j]));	// store the border indices
				groupCnt++;
				groupsStack.push_back(pair<int,int>(mc[i][j], groupCnt));		// store the start and the groupNr on a Stack
			}
		}
	}

		// include the missing chains (the one-element chains) at the corect position
	for(unsigned i=0; i<groups.size(); i++){
		int curSubGroup=0;
			// we have stored the borders of each group
			// now iterate through all indices between them
		for (int j=groupsBorder[i].first-1; j>groupsBorder[i].second; j--){
				// if a chain is already included then jump over it (and test for empty groups)
			if(groups[i].size() && curSubGroup < (int)groups[i].size() && j == groups[i][curSubGroup].first){
				j = groups[i][curSubGroup].second;
				curSubGroup++;
				continue;
			}
				// the other elements should be included
			groups[i].insert(groups[i].begin()+curSubGroup, pair<int,int>(j,j));
			curSubGroup++;	// don't be tricked by the new element
		}
	}
		// finaly get the preserving reversals
		// in every group all nonoverlaping neighbours could be reversed together
		// ! subtract 1 to start with index 0
	for(unsigned i=0;i<groups.size();i++)
		for(unsigned s=0; s<groups[i].size();s++)
			for(unsigned l=s; l<groups[i].size();l++)
				pr.push_back(pair<int,int>(groups[i][l].second-1, groups[i][s].first-1));

	return pr;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int intervalDistance(const genom &source, const genom &target){
	vector<vector<int> > mc = getMaximalChains(source, target);
	int s = source.size() + 2,
		n = 0;

		// interval distance = N1 + N2 - 2N
		//							= (s * s-1)/2 + (s * s-1)/2 - 2*n
		//							= (s * s-1)/2 - 2*n
		// s is the size of the genome (+2 (frame elements))
		// n is the count of the common conserved intervals

		// each maximal chain of irreducible intervals
		// contribute k*(k-1)/2 to the count of conserved intervals
		// where k is the count of the elements
	for (unsigned i=0; i<mc.size(); i++)
		n += (mc[i].size()*(mc[i].size()-1))/2;

	return s*(s-1) - (2*n);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int intervalDistance(const genom &source, const vector<genom> &targets){
	int dist=0;

	for(unsigned i=0; i<targets.size(); i++){
		dist += intervalDistance(source, targets[i]);
	}
	return dist;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//int preserving_reversal_distance(const genom &src,
//		const genom &tgt, int circular ){
//	unsigned cycle_cnt,
//		component_cnt,
//		unoriented_component_cnt = 0;
//	vector<int> pi,
//		component,			// component index of each point
//		cycle,						// cycle index of each point
//		component_orientation,		// orientation of each component
//		component_length;
//
//		// get the cycles and components
//	pi = src.identify(tgt);		// make 2 permutions to 1 and id
//	cycle_cnt = cycles(pi, src.size(), cycle);
//	component = src.getComponents(tgt, component_orientation, component_cnt);
//	component_length = vector<int>(component_cnt, 0);
//
//		// count unoriented components
//	for(unsigned i=0; i<component_orientation.size(); i++){
//		if(component_orientation[i]==0){
//			unoriented_component_cnt++;
//		}
//	}
//
//	return src.size()+1-cycle_cnt+unoriented_component_cnt;
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//extern "C" int preserving_reversal_distance(
//	int *src_p, int *tgt_p,
//	int n, char circular){
//
//	genom src(src_p, n, circular),
//		tgt(tgt_p, n, circular);
//
//	return preserving_reversal_distance(src, tgt, circular);
//}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector< vector<int> > valid_adjacencies(const vector< genom > &G,
		vector< vector< int > > &dependencies){

	pair<unsigned, unsigned> e;
	unsigned n = G[0].size(),
		ci_cnt = 0,
		cur_mc_end,
		end_of_outer_mc;
	vector< bool > single_element(G[0].size(), true);;
	vector< int > g0chr;
	vector< vector<int> > mc;
	vector< vector< pair< unsigned , unsigned > > > ci_elements;
	vector< vector<int> > va(2*G[0].size()+2, vector<int>(2*n+2, 0));

	mc = getMaximalChains(G);
#ifdef DEBUG_VALID_ADJACENCIES
	cout << "maximal chains"<<endl;
	for(unsigned i=0; i<mc.size(); i++){
		for(unsigned j=0; j<mc[i].size(); j++)
			cout << mc[i][j] << " ";
		cout << endl;
	}
#endif
		// iterate over all maximal chains and for each through the irreducible intervals it contains
	for(unsigned i=0; i<mc.size(); i++){
			// store the end of the outermost maximal chain, i.e. irreducible intervals at an index (in ci_elements)
			// which is greater than the stored endindex does not belong to the outermost maximal chain (and can be reversed)
		if(i==1){
			end_of_outer_mc = ci_elements.size()-1;
		}
		for(unsigned j=0; j<mc[i].size()-1; j++){
				// reset the no_single_element vector
			single_element.assign(n, true);
				// push the new frame
			ci_elements.push_back( vector< pair< unsigned , unsigned > >() );
				// store the bounding elements as first element
			ci_elements[ci_cnt].push_back( pair< unsigned , unsigned >(mc[i][j+1] , mc[i][j] ) );
				// check if one of the following maximal chains is contained in the
				// current irreducible interval
			cur_mc_end = mc[i][j];
			for(unsigned k=i+1; k<mc.size(); k++){
					// if the currently checked maximal chain starts with an outside the boundings of
					// the current irreducible interval -> it (and the following) could not be included
					// in the current irreducible interval
				if(mc[k][0] < mc[i][j+1] || mc[k].back() > mc[i][j] ){
					continue;
				}

					// if there are more than one maximal chains in the current irreducible interval;
					// then they are maybe nested; therefore the endindex of the latest included maximal
					// chain is stored; if the start index of the current maximal chain is greater than
					// this index it is nested;
				if((unsigned) mc[k][0] > cur_mc_end){
					continue;
				}

					// if none of the two criteria above were true -> append the start- end-index
					// pair to the elements of the current ci. and store the endindex of the latest
					// included maximal chain (i.e. this chain)
				ci_elements[ci_cnt].push_back( pair< unsigned , unsigned >(mc[k][ mc[k].size()-1 ] , mc[k][0] ) );
				cur_mc_end = mc[k][ mc[k].size()-1 ];

					// mark all elements in the added maximal chain as 'this is not a single element'
				for(unsigned l=mc[k][ mc[k].size()-1 ]; l<= (unsigned) mc[k][0]; l++){
					single_element[l] = false;
				}
			}
				// go through the interior of the irreducible interval
			for(unsigned k=mc[i][j+1]+1; k<(unsigned) mc[i][j]; k++ ){
				if(single_element[k])
					ci_elements[ci_cnt].push_back( pair< unsigned , unsigned >( k, k ) );
			}
			ci_cnt++;
		}
	}

#ifdef DEBUG_VALID_ADJACENCIES
	for(unsigned i=0; i<ci_elements.size(); i++){
		cout << "ci_elements" <<i<<endl;
		for (unsigned j=0; j<ci_elements[i].size(); j++){
			cout <<"("<<ci_elements[i][j].first<<","<<ci_elements[i][j].second<<") ";
		}
		cout << endl;
	}
#endif
		// get G[0] and frame it with 0 and n+1 ; to get the elements behing the indices
	g0chr = G[0].getChromosom();
	g0chr.insert(g0chr.begin(), 0);
	g0chr.push_back(g0chr.size());

		// construct the valid adjacencies
	for(unsigned i=0; i<ci_elements.size(); i++){
			// if this one contains no elements -> then only the elements at start and end indices form a valid pair
			// and the reverse if its not the outermost maximal chain
		if(ci_elements[i].size() == 1){
			e = matching_edge( g0chr[ ci_elements[i][0].first ], g0chr[ ci_elements[i][0].second ], n + 1 );
			va[e.first][e.second] = 1;
			//~ if(i > end_of_outer_mc){
				va[e.second][e.first] = 1;
			//~ }
		}

			// if there are elements -> form all valid pairs
			// consider an irreducible interval with bounding elements b1, b2 and internal maximal chains bounded
			// by i1, i2 resp. j1, j2 (note that a single element a is modeled as a1 = a2 = a).
			// then the valid pairs are :
			// - bounding elements and all inner maximal chains: (b1,i1) ; (b1-i2) ; (i2,b2), (-i1,b2)
			// 	because all (but not the outermost maximal chain 0..n+1) can be reversed also the pairs
			// 	(-b2,i1) ; (-b2-i2) ; (i2,-b1), (-i1,-b1) are valid
			// - pairs of inner nodes: (i2,j1); (i2, -j2); (-i1,j1); (-i1,-j2); (j2,i1); (j2,-i2); (-j1,i1); (-j1,-i2)
			// from these pairs of elements the edges in the MBGraph are computed
			// note: if (x,y) and (-y,-x) are valid pairs. and the matching edge coresponding to (x,y) is (e,f) then
			// (f,e) is the matching edge coresponding to (-y,-x).
			// -> the eight cases for pairs of inner maximal chains reduce to four cases.
		else{
			for(unsigned j=1; j<ci_elements[i].size(); j++){
					// between the bounding elements and the inner elements
					// (b1,i1)
				e = matching_edge(g0chr[ ci_elements[i][0].first ], g0chr[ ci_elements[i][j].first ], n+1);
				va[e.first][e.second] = 1;
				//~ if(i > end_of_outer_mc){
					va[e.second][e.first] = 1;
				//~ }
					// (b1-i2)
				e = matching_edge(g0chr[ ci_elements[i][0].first ], -1* g0chr[ ci_elements[i][j].second ], n+1);
				va[e.first][e.second] = 1;
				//~ if(i> end_of_outer_mc){
					va[e.second][e.first] = 1;
				//~ }
					// (i2,b2)
				e = matching_edge(g0chr[ ci_elements[i][j].second ], g0chr[ ci_elements[i][0].second ], n+1);
				va[e.first][e.second] = 1;
				//~ if(i> end_of_outer_mc){
					va[e.second][e.first] = 1;
				//~ }
					// (-i1,b2)
				e = matching_edge(-1 * g0chr[ ci_elements[i][j].first ], g0chr[ ci_elements[i][0].second ], n+1);
				va[e.first][e.second] = 1;
				//~ if(i> end_of_outer_mc){
					va[e.second][e.first] = 1;
				//~ }

					// and between each air of inner elements
				for(unsigned k=j+1; k<ci_elements[i].size(); k++){
						// (i2,j1)
					e = matching_edge(g0chr[ ci_elements[i][j].second ], g0chr[ ci_elements[i][k].first ], n+1);
					va[e.first][e.second] = 1;
						// (-j1,-i2)
					va[e.second][e.first] = 1;

						// (i2, -j2)
					e = matching_edge(g0chr[ ci_elements[i][j].second ], -1 * g0chr[ ci_elements[i][k].second ], n+1);
					va[e.first][e.second] = 1;
						// (j2,-i2)
					va[e.second][e.first] = 1;

						// (-i1,j1)
					e = matching_edge(-1 * g0chr[ ci_elements[i][j].first ], g0chr[ ci_elements[i][k].first ], n+1);
					va[e.first][e.second] = 1;
						// (-j1,i1)
					va[e.second][e.first] = 1;

						// (-i1,-j2)
					e = matching_edge(-1 * g0chr[ ci_elements[i][j].first ], -1 * g0chr[ ci_elements[i][k].second ], n+1);
					va[e.first][e.second] = 1;
						// (j2,i1)
					va[e.second][e.first] = 1;
				}
			}
		}
	}

		// forbid the impossible adjacencies (x,0) and (n+1,x)
	for(unsigned i=0; i<va.size(); i++){
		va[i][0] = 0;
		va[i].back() = 1;
	}
	va.back().assign(2*n+2, 0);


		/*!@todo dont know why but this seems important */
	va[0].assign(2*n+2, 1);

		// finaly get the element dependencies. if the direction of a gene, which is a linking element of a maximal
		// chain, is fixed then the signs of the other linking elements of the chain are determined
		// so we store for each element which is part of a maximal chain the depending elements
		//
		// e.g. maximal chain linking elements: -3 5 1 ; n=5; then
		// dependencies[2] = {10, 6}	 3 -> 5,1
		// dependencies[8] = {0,4}	-3 -> -5,-1
		// dependencies[10] = {2,6}	 5 -> -3,1
		// dependencies[0] = {4,8}	-5 -> 3,-1
		// dependencies[6] = {2,10}	 1 -> -3,5
		// dependencies[4] = {0,8}	-1 -> 3,-5

		// first init or clear the memory
	if(dependencies.size() != 2*n+2){
		dependencies.resize(2*n+2);
	}
	for(unsigned i=0; i<dependencies.size(); i++){
		dependencies[i].clear();
	}

	int gene_j;
	int gene_k;

		// assign the dependencies
	for(unsigned i=0; i<mc.size(); i++){
		for(unsigned j=0; j<mc[i].size(); j++){
			if(mc[i][j] == (int)0 || mc[i][j] == (int)n+1)	// the bounding elements at 0 and n+1 aren't relevant
				continue;

			gene_j = G[0][ mc[i][j] - 1 ];

			for(unsigned k=0; k<mc[i].size(); k++){
				if(k == j)
					continue;

				if (mc[i][k] == (int)0 || mc[i][k] == (int)n+1)
					continue;

				gene_k = G[0][ mc[i][k] - 1 ];

					// same sign -> could not be different signs
				if((gene_j > 0 && gene_k > 0) || (gene_j < 0 && gene_k < 0) ) {
					dependencies[ 2*abs(gene_j) ].push_back( 2*abs(gene_k) - 1 );	// gene_j negative -> gene_k can not be positive
					dependencies[ 2*abs(gene_j) -1 ].push_back(2*abs(gene_k));		// gene_j positive   -> gene_k can not be negative
				}
					// different signs -> could not be same signs
				else if( (gene_j > 0 && gene_k < 0 ) || (gene_j < 0 && gene_k > 0) ){
					dependencies[ 2*abs(gene_j) - 1 ].push_back( 2*abs(gene_k) - 1);	// gene_j positive -> gene_k can not be positive
					dependencies[ 2*abs(gene_j) ].push_back( 2*abs(gene_k) );	// gene_j negative -> gene_k can not be negative
				}
			}
		}
	}
#ifdef DEBUG_VALID_ADJACENCIES
	for(unsigned i=0; i<va.size(); i++){
		for(unsigned j=0; j<va[i].size(); j++)
			cout << va[i][j]<<" ";
		cout << endl;
	}
#endif
	return va;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*extern "C" int ** valid_adjacencies(
		struct genome_struct **passgenome,
		int ngenomes,
		int ngenes,
		int circular,
		int ***dep,
		int **dep_len){

	vector<genom> G;
	vector<vector<int> > va_vec;
	int **va;
	vector< vector< int > > dep_vec;

	for(unsigned i=0; i<(unsigned)ngenomes; i++){
		G.push_back(genom(passgenome[i]->genes, ngenes, circular));
	}

	va_vec = valid_adjacencies(G, dep_vec);

		// make c arrays for valid adj
	va = (int**)malloc(va_vec.size() * sizeof(int *));
	for(unsigned i=0; i<va_vec.size(); i++){
		va[i] = (int*) malloc(va_vec[i].size() * sizeof(int));
		for(unsigned j=0; j<va_vec[i].size(); j++){
			va[i][j] = va_vec[i][j];
		}
	}

		// make c arrays for dependencies
	*dep = (int**) malloc(dep_vec.size() * sizeof(int *));
	*dep_len = (int*) malloc(dep_vec.size() * sizeof(int));
	for(unsigned i=0; i<dep_vec.size(); i++){
		(*dep_len)[i] = dep_vec[i].size();
		(*dep)[i] = (int*) malloc(dep_vec[i].size() * sizeof(int));
		for(unsigned j=0; j<dep_vec[i].size(); j++){
			(*dep)[i][j] = dep_vec[i][j];
		}
	}

	return va;
}*/
