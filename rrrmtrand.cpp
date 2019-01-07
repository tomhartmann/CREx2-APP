/*
 * rrrmtrand.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: maze
 */

#include "common.hpp"
#include "rrrmtrand.hpp"

namespace std {

rrrmt_rand::rrrmt_rand( int k, vector<float> prob, genom &g, vector<genom> &trace, bool randrange, unsigned ml){
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
		insert( tmprandev );
	}
}

rrrmt_rand::~rrrmt_rand() {
}

rrrmt_rand_sit::rrrmt_rand_sit(  unsigned k, vector<float> prob, genom &g, vector<genom> &trace, unsigned dm, unsigned dM  ) {
	genom gp = g;
	itnode *itroot, *nd;
	unsigned i = 0,
			evtype;
	vector<itnode*> nodes;
	rrrmt *tmp;

	trace.clear();
	trace.push_back( g );

	weighted_choice_init(prob);

	interval_tree_random( &itroot, NULL, make_pair(0,g.size()-1), dm, dM );

//	cerr << "simulated itree "<< endl;
//	interval_tree_print( itroot, g, cerr ); cerr<<endl;

	insert_iterator<vector<itnode*> > x = inserter(nodes, nodes.begin());
	interval_tree_nodes( itroot, x, true );

	do{
		// get a random node
		nd = nodes[ask_rng( 0, nodes.size()-1 )];
//		interval_tree_print(nd, gp, cerr); cerr<<endl;
		// get a random rearrangement
		evtype = weighted_choice( prob );
		switch ( evtype ){
			case 0: {	// reversal
				tmp = new rev( nd->i.first, nd->i.second, gp );
				break;
			}
			case 1: {	// transposition
				if( nd->children.size() !=2 ) continue;
				tmp = new transp(nd->children[0]->i.first, nd->children[0]->i.second+1, nd->children[1]->i.second+1, gp);
				break;
			}
			case 2: {	// reverse transposition
				if( nd->children.size() < 2 ) continue;

				if( ask_rng(0,1) == 0 ){
					tmp = new revtransp(nd->children[0]->i.first, nd->children[nd->children.size()-2]->i.second,
						nd->children.back()->i.first, nd->children.back()->i.second, gp);
				}else{
					tmp = new revtransp(nd->children[1]->i.first, nd->children.back()->i.second,
							nd->children[0]->i.first, nd->children[0]->i.second, gp);
				}

				break;
			}
			case 3: {	// tdrl
				//cerr << "skipping TDRLs"<<endl;
				continue;
				break;
			}
			default:{
				cerr << "internal error: unknown event type choosen"<<endl;
			}
		}
		// insert the rearrangement do not insert twice
		if( ! insert( tmp ) ){
			continue;
		}

//		cout << *tmp<<endl;
		tmp->apply(g);
		trace.push_back(g);

		i++;
	}while( i < k );

	interval_tree_free( itroot );
}

rrrmt_rand_sit::~rrrmt_rand_sit() {
}
} /* namespace std */
