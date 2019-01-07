#include <algorithm>
#include <limits>
#include <math.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "sb4type.hpp"

//#define DEBUG

using namespace std;

///////////////////////////////////////////////////////////////////////////////

genom multiply(genom pi, genom sigma){
	if(pi.size() != sigma.size()){
		cerr << "ERROR: pi and sigma do not have the same size! ---> exiting!" << endl;
		exit(EXIT_FAILURE);
	}
	genom product;
	product.setCircular(sigma.getCircular());
	product.set_nmap(sigma.get_nmap());

	//π°σ is the function that maps any element x of the set to π(σ(x))
	vector<int> tmp;
	for(unsigned int x=0; x< pi.size(); x++){
		int value;
		value=abs(pi[abs(sigma[x])-1]);
		// take care about the signs
		if((sigma[x]<0 && pi[abs(sigma[x])-1]>0) || (sigma[x]>0 && pi[abs(sigma[x])-1]<0)){
			value*=-1;
		}
		tmp.push_back( value );
	}
	product.setChromosom(tmp);
	return product;
}

///////////////////////////////////////////////////////////////////////////////

void sort_by_4type(genom pi, genom sigma, rrrmt* &scenario){
	if(pi.size()!=sigma.size()){
		cerr << "ERROR: genomes do not have the same length: " << pi << " <-> " << sigma << "! ---> exiting!" << endl;
		exit(EXIT_FAILURE);
	}
	if(pi!=sigma){
		// Recall that pi -> sigma <=> iota -> pi^{-1}*sigma
		// STEP 1: obtain pi^{-1}*sigma
		// STEP 2: Calculate scenario iota -> pi^{-1}*sigma
		// STEP 3: transform scenario to pi->sigma

		// STEP 1:
		// 			-> get pi^{-1]
		vector<int > cr_ipi(pi.size(), std::numeric_limits< int >::max());	// chromosome of inverse permutation
		for(unsigned i=0; i<pi.size(); i++){
			cr_ipi[ abs(pi[i])-1 ] = i+1;
			if(pi[i] < 0){
				cr_ipi[ abs(pi[i])-1 ] *= -1;
			}
		}
//		cout << "SIZE: " << cr_ipi.size() << "; SIZE sigma=" << sigma.size() << endl;
//		cout << "pi_chromo:\t ";
//		for(int i=0; i<pi.size(); i++){
//			cout << pi[i]<< " ";
//		}
//		cout << endl;
//		cout << "si_chromo: \t";
//		for(int i=0; i<sigma.size(); i++){
//			cout << sigma[i]<< " ";
//		}
//		cout << endl;
//		cout << "inv_chromo: \t";
//		for(int i=0; i<cr_ipi.size(); i++){
//			cout << cr_ipi[i]<< " ";
//		}
//		cout << endl;
		genom pii=pi;
		pii.setChromosom(cr_ipi);
		pii.setCircular(pi.getCircular());
		pii.set_nmap(pi.get_nmap());
		//			-> apply pi^{-1} to sigma
//		cout << "pi: " << pi << endl;
//		cout << "sigma: " << sigma << endl;
//		cout << "pii: " << pii << endl;
		genom temp=multiply(pii, sigma);
//		cout << "temp: " << temp << endl;

		// STEP 2:
		vector <rrrmt*> rhos;
		sort_by_4type(temp, rhos);
//		for(int i=0; i<rhos.size(); i++){
//			cout << *(rhos[i]) << endl;
//		}

		// build renaming mapping
		vector< vector<int> > temp_to_pi_map;	// mapping from temp to pi
		temp_to_pi_map = vector<vector<int> >( pi.size()+1 );
		for( unsigned i=0; i<temp.size(); i++ ){
				temp_to_pi_map[ abs(temp[i]) ].push_back( abs(pi[abs(temp[i])-1]) );
		}
//		cout << "pi_chromo:\t ";
//		for(int i=0; i<pi.size(); i++){
//			cout << pi[i]<< " ";
//		}
//		cout << "temp_chromo:\t ";
//		for(int i=0; i<temp.size(); i++){
//			cout << temp[i]<< " ";
//		}
//		cout << "mapping temp to pi:" << endl;
//		for(int i=0; i<temp_to_pi_map.size();i++){
//			cout << "i: " << i << "-> ";
//			for(int j=0; j<temp_to_pi_map[i].size(); j++){
//				cout << temp_to_pi_map[i][j] << " ";
//			}
//			cout << endl;
//		}

		// STEP 3:
		//			-> reverse the order of rhos and rename rrrmts
		vector <rrrmt*> new_rhos;
		for(int i=rhos.size()-1; i>=0; i--){
			rhos[i]->rename(temp_to_pi_map, pi.get_nmap(), pi.size());
			new_rhos.push_back(rhos[i]->clone());
		}


//		vector <rrrmt*> new_rhos;
//
//			vector< set<int> > elems;
//			rhos[i]->elements(elems, false);
////			cout << "elems_size=" << elems.size() << endl;
//			// if rhos[i]==I; elems[0] reversed elements;
//			// if rhos[i]==T; elems[0] interval 1; elems[1] interval 2
//			// if rhos[i]==iT: elems[0} are the reversed elements; elems[1] are the trnasposed elements
//			// if rhos[i]==TDRL: elems[0] first copy; elems[1] second copy
//
//			//#define TRA 2
//			//#define TDRL 3
//			//#define RTRA 8
//
//			if(rhos[i]->type()==7){// case inversion
//				set<int> new_elems;
//				set<int>::iterator iter;
//				for(iter=elems[0].begin();iter!=elems[0].end();++iter){
//					//cout << *iter << endl;
////					if( ( *iter>0 && pi[abs(*iter)-1]>0 ) || ( *iter<0 && pi[abs(*iter)-1]<0 ) ){
//						new_elems.insert(abs(pi[abs(*iter)-1]));
////					}else{
////						new_elems.insert((-1)*abs(pi[abs(*iter)-1]));
////					}
//				}
//				rrrmt* re= new rev(new_elems, sigma.get_nmap());
////				cout << "new rev:" <<  *re <<endl;
//				new_rhos.push_back(re);
//			}
//			if(rhos[i]->type()==2 || rhos[i]->type()==8 || rhos[i]->type()==3 ){//case transposition or inverse transposition or TDRL
// 				set<int> new_first;
//				set<int> new_second;
//				set<int>::iterator iter;
////				cout << "first set: ";
//				for(iter=elems[0].begin();iter!=elems[0].end();++iter){
////					cout << *iter << " ";
////					if( ( *iter>0 && pi[abs(*iter)-1]>0 ) || ( *iter<0 && pi[abs(*iter)-1]<0 ) ){
//						new_first.insert(abs(pi[abs(*iter)-1]));
////					}else{
////						new_first.insert((-1)*abs(pi[abs(*iter)-1]));
////					}
//				}
////				cout << endl;
////				cout << "size:" << elems[1].size() << endl;
////				for(iter=elems[1].begin(); iter!=elems[1].end(); ++iter){
////					cout << *iter << " " <<endl;
////				}
////				cout << "second set: ";
//				for(iter=elems[1].begin();iter!=elems[1].end();++iter){
////					cout << *iter << " ";
////					if( ( *iter>0 && pi[abs(*iter)-1]>0 ) || ( *iter<0 && pi[abs(*iter)-1]<0 ) ){
//						new_second.insert(abs(pi[abs(*iter)-1]));
////					}else{
////						new_second.insert((-1)*abs(pi[abs(*iter)]));
////					}
//				}
////				cout << endl;
//				if(rhos[i]->type()==2){//case transposition
//					vector <set<int> > bla;
//					bla.push_back(new_first);
//					bla.push_back(new_second);
//					rrrmt* t=new transp(bla, sigma.get_nmap());
//					new_rhos.push_back(t);
//				}else if(rhos[i]->type()==8){// case inverse transposition
//					rrrmt* rt= new revtransp( new_first, new_second, sigma.get_nmap() );
//					new_rhos.push_back(rt);
//				}else{// case TDRL
//					vector <set<int> > bla;
////					set<int> new_first2;
////					set<int> new_second2;
////					for(iter=new_first.begin();iter!=new_first.end();++iter){
////						new_first2.insert(abs(*iter));
////					}
////					for(iter=new_second.begin();iter!=new_second.end();++iter){
////						new_second2.insert(abs(*iter));
////					}
////					bla.push_back(new_first2);
////					bla.push_back(new_second2);
//					bla.push_back(new_first);
//					bla.push_back(new_second);
//					tdrl* tdl= new tdrl(bla, sigma.get_nmap(), true);
//					vector< int> vec;
//					vec.push_back(-1);
//					tdl->set_orgord(vec);
//					new_rhos.push_back(tdl);
//				}
//			}//END CASE TRA or REV or TDRL
//			if((rhos[i]->type() != 3) && (rhos[i]->type() != 7) && (rhos[i]->type() != 8) && (rhos[i]->type() != 2) ){
//				cerr << "Error: Unknown type of rrrrmt ---> exiting!"<< endl;
//				cerr << "type=" << rhos[i]->type() << endl;
//				exit(EXIT_FAILURE);
//			}
//		}//END forall rhos
		rrrmt* scen= new ordered(new_rhos, false);
		scenario=scen->clone();
	}// END pi != sigma
}

///////////////////////////////////////////////////////////////////////////////

void sort_by_4type(genom g, vector<rrrmt* > &rrrmts){
	// check g != iota
	genom iota=genom( g.size(), 0 );

	// just for debugging
//	cout << "*************************************************************************" <<endl;
//	vector<int> g2=g.getChromosom();
//	for(int i=0; i<g2.size();i++){
//		cout << g2[i] << " ";
//	}
//	cout << endl;

	if(g != iota){// scenario cannot be empty; calculate scenario
		genom temp=g;
//		int test=0;

		while(temp != iota){
//			test++;
#ifdef DEBUG
			cout << "temp: ";
			for (int z=0; z< temp.size(); z++){
				cout << temp[z] << " ";
			}
			cout << endl;
#endif //DEBUG
			rrrmt *rho; // new rrrmt of the sorting scen in reversed order (\rho_j)
			unsigned int mis=temp.max_inc_substrings();
			// calculate all strict maximal increasing substring
			vector< set<int> > decomp;
			temp.decomposition_strict(decomp);
			// find position of first/last negative element
			int max_neg_pos=0; // position of last negative element
			int min_neg_pos=temp.size()-1; // position of first negative element
			for(unsigned int i=0; i<temp.size();i++){
				if(temp[i] < 0){
					if(max_neg_pos<i){
						max_neg_pos=i;
					}
					if(min_neg_pos>i){
						min_neg_pos=i;
					}
				}
			}
			if(mis==1){// one one maximum increasing substring
				// check that all negative elements of temp form a substring
				for (int i=0; i<=max_neg_pos; i++){
					if(temp[i]>0){
						cerr << "ERROR: temp is not of form (---- ++++) ---> exiting!" << endl;
						exit(EXIT_FAILURE);
					}
				}

//				cout << abs(temp[min_neg_pos]) << " - " << abs(temp[max_neg_pos]) << "=?=" << max_neg_pos << " - " << min_neg_pos <<endl;
				bool iTcond1= (abs(temp[min_neg_pos])-abs(temp[max_neg_pos])) == (max_neg_pos-min_neg_pos);
				bool iTcond2= decomp.size()>1;
				bool iTcond3= decomp.size()<4;
				bool iTcond4= false;
				if(decomp.size()>1){
					if(temp[max_neg_pos+1]==1){
						iTcond4=true;
					}
				}
//				cout << iTcond1 << " " << iTcond2 << " " << iTcond3 << " " << iTcond4 << endl;
 				if( iTcond1 && iTcond2 && iTcond3 && iTcond4){
					// temp is  (-l...-m 1...m-1 l+1...n) -> case inverse transposition
					// find position of element l+1
					int k=0;
					for(int i=max_neg_pos+1; i<temp.size(); i++){
						if(temp[i]==abs(temp[max_neg_pos])-1){
							k=i;
						}
					}
					if(k==0){
						cerr << "End of second interval not found ---> exiting!" << endl;
						exit(EXIT_FAILURE);
					}

					// define inverse transposition
					rho= new revtransp(0,max_neg_pos,max_neg_pos+1, k, temp);
					// apply
					rho->apply(temp);
					// store
					rrrmts.push_back(rho);
					// just for debugging
//					cout << "here" << endl;
//					cout << g << endl;
//					cout << temp << endl;
				}else{
					// temp is (-1...-l l+1...n) or (---- +++++) -> case inversion
					// define reversal
					rho= new rev(0, max_neg_pos, temp);
					// apply
//					for(int z=0; z< temp.size();z++){
//						cout << temp[z] << " ";
//					}
//					cout << endl;
					rho->apply(temp);
//					for(int z=0; z< temp.size();z++){
//						cout << temp[z] << " ";
//					}
//					cout << endl;
					// store
					rrrmts.push_back(rho);
					// just for debugging
					//cout << g << endl;
				} // case inversion
			}else{ // more than one maximal increasing substring
				// just for debugging : show strict max inc substrings
//				for(int i=0; i<decomp.size();i++){
//					set<int>::iterator it;
//					for(it=decomp[i].begin(); it!=decomp[i].end(); ++it){
//						cout << *it << " ";
//					}
//					cout << endl;
//				}
//				cout << endl;
				// necessary conditions for the one rrrmt cases
				bool cond1=(mis==2);			// only two maximal increasing substrings
				bool cond2=(temp[0]>0);		// first MIS contains only positive elements -> temp is of form (+++++ -----(++++))
				bool cond3=(decomp.size()<=4);	// maximal 4 strict maximal increasing substrings
				// check if all negative elements form a maximal increasing substring
				bool cond4=true;
				for (int i=min_neg_pos; i<=max_neg_pos; i++){
					if(temp[i]>0){
						cond4=false;
					}
				}
				// check if for temp exists l,m such that \Upsilon(temp)=[l:m]
				bool cond5=true;
				set<int> Upsilon;
				for(int i=0; i<temp.size(); i++){
					if(temp[i]<0){
						Upsilon.insert(abs(temp[i]));
					}
				}
				if(Upsilon.size()>1){
					set<int>::iterator iter;
					set<int>::iterator iterend = Upsilon.end();
					--iterend;
					for(iter=Upsilon.begin(); iter!=iterend; ++iter){
						int b = *iter;
						iter++;
						int a = *iter;
						iter--;
						if(a - b != 1){
							cond5=false;
						}
					}
				}
				// check if Upsilon is empty
				bool cond6 = !(Upsilon.empty());

				// just for debugging
//				cout << cond1 << cond2 << cond3 << cond4 << cond5 << cond6 << endl;
				// maybe case one inversion or inverse transposition
				if(cond1 && cond2 && cond3 && cond4 && cond5 && cond6){
					// simply try a single inversion \rho; if \rho(temp)=iota; it was the right one;
					if(decomp.size()<=3){
//						cout << min_neg_pos << " " << max_neg_pos << endl;
						// try inversion
						rrrmt* test;
						genom testg = temp;
						test= new rev(min_neg_pos,max_neg_pos,testg);
						test->apply(testg);
						if(testg==iota){
							rho= new rev(min_neg_pos,max_neg_pos,temp);
							rho->apply(temp);
							rrrmts.push_back(rho);
							delete test;
							continue;
						}
						delete test;
					}
					// try every possible inverse transposition
					if(decomp.size()==2){
						rrrmt* test;
						genom testg = temp;
						test= new revtransp(min_neg_pos,max_neg_pos, 0, decomp[0].size()-1, testg);
						test->apply(testg);
						if(testg==iota){
							rho= new revtransp(min_neg_pos,max_neg_pos, 0, decomp[0].size()-1, temp);
							rho->apply(temp);
							rrrmts.push_back(rho);
							delete test;
							continue;
						}
						delete test;
					}
					if(decomp.size()==3){
						rrrmt* test;
						genom testg=temp;
						// there are two cases Upsilon is the second or the third strict MIS
						if(temp[temp.size()-1]<0){// case that upsilon is the third
							test= new revtransp(min_neg_pos, max_neg_pos, decomp[0].size(), decomp[0].size()+decomp[1].size()-1, testg);
							test->apply(testg);
							if(testg==iota){
								rho= new revtransp(min_neg_pos, max_neg_pos, decomp[0].size(), decomp[0].size()+decomp[1].size()-1, temp);
								rho->apply(temp);
								rrrmts.push_back(rho);
								delete test;
								continue;
							}
							test->deapply(testg);
							delete test;
						}else{// case that upsilon the second
							// there are two cases: that transpose upsilon to the left or right
							if(temp[0]>abs(temp[min_neg_pos])){ // case to the left
								test= new revtransp(min_neg_pos, max_neg_pos, 0, decomp[0].size()-1, testg);
								test->apply(testg);
								if(testg==iota){
									rho= new revtransp(min_neg_pos,max_neg_pos, 0, decomp[0].size()-1, temp);
									rho->apply(temp);
									rrrmts.push_back(rho);
									delete test;
									continue;
								}
								test->deapply(testg);
								delete test;
							}else{// case to the right
								test= new revtransp(min_neg_pos, max_neg_pos, decomp[0].size() + decomp[1].size(), testg.size()-1, testg);
								test->apply(testg);
								if(testg==iota){
									rho= new revtransp(min_neg_pos,max_neg_pos, decomp[0].size()+decomp[1].size(), testg.size()-1, temp);
									rho->apply(temp);
									rrrmts.push_back(rho);
									delete test;
									continue;
								}
								test->deapply(testg);
								delete test;
							}
						}
					}// end decomp.size()==3
					if(decomp.size()==4){
						rrrmt* test;
						genom testg=temp;
						// there are two cases 1 <= o < p < q <= n:
						// Case 1: temp = (1...o-1 p+1...q -p...-o q+1....n) -> iT(-p...-o; p+1...q)
						// Case 2: temp = (1...o-1 -q...-(p+1) o...p q+1...n) -> iT(-q...-(p+1); o...p)
						if(temp[decomp[0].size()]>0){ // CASE 1
							test = new revtransp(min_neg_pos, max_neg_pos, decomp[0].size(), min_neg_pos-1, testg);
							test->apply(testg);
							if(testg==iota){
								rho = new revtransp(min_neg_pos, max_neg_pos, decomp[0].size(), min_neg_pos-1, temp);
								rho->apply(temp);
								rrrmts.push_back(rho);
								delete test;
								continue;
							}
							test->deapply(testg);
							delete test;
						}else{ // CASE 2
							test = new revtransp(min_neg_pos, max_neg_pos, max_neg_pos+1,  max_neg_pos+decomp[2].size(), testg);
							test->apply(testg);
							if (testg==iota){
								rho = new revtransp(min_neg_pos, max_neg_pos, max_neg_pos+1,  max_neg_pos+decomp[2].size(), temp);
								rho->apply(temp);
								rrrmts.push_back(rho);
								delete test;
								continue;
							}
							test->deapply(testg);
							delete test;
						}
					} // end decomp.size()=4
				}// END maybe sort-able by application of a single inversion or inverse transposition

				//
				// apply transformation T and store corresponding TDRL:
				//
				// -> step1: build MIS decomposition of TEMP
				vector < set<int> > decomposition;
				temp.decomposition(decomposition);
				// just for debugging only
//				for(int i=0; i<decomposition.size();i++){
//					set<int>::iterator it;
//					for(it=decomposition[i].begin(); it!=decomposition[i].end(); ++it){
//						cout << *it << " ";
//					}
//					cout << endl;
//				}
//				cout << endl;
				vector<int> new_chromo; // chromosome that replaces the one of temp
				// apply transformation T to temp
				set<int> result_set;
				for(int i=0; i<floor((float)decomposition.size()/2); i++){
					oplus(result_set, decomposition[i], decomposition[ (int)ceil((float)decomposition.size()/2) +i ]);
					set<int>::iterator ita;
					for(ita=result_set.begin(); ita!=result_set.end(); ++ita){
						new_chromo.push_back(*ita);
					}
					result_set.clear();
				}
				if(decomposition.size()%2==1){
					int middle= ceil((float)decomposition.size()/2)-1;
					set<int>::iterator ita;
					for(ita=decomposition[middle].begin(); ita!=decomposition[middle].end(); ++ita){
						new_chromo.push_back(*ita);
					}
				}
				//just for debugging
//				for (int i=0; i<new_chromo.size(); i++){
//					cout << new_chromo[i] << " ";
//				}
//				cout << endl;
				temp.setChromosom(new_chromo);

				set <int> first;
				set <int> second;
				for(int i=0; i < decomposition.size(); i++){
					set<int>::iterator ita;
					for(ita=decomposition[i].begin(); ita!=decomposition[i].end(); ++ita){
						if(i < ceil((float)decomposition.size()/2)){
							first.insert(*ita);
						}else{
							second.insert(*ita);
						}
					}
				}
				// define TDRL; check whether it is a transposition or not; if yes -> store transposition
				vector< bool > SF;
				for(int i=0; i<temp.size(); i++){
					if( first.find(temp[i]) != first.end() ){ //found element in first copy
						SF.push_back(true);
					}else{ // element must be in second copy
						SF.push_back(false);
					}
				}
				rrrmt* test;
				test = new tdrl(SF, temp);
				//just for debugging
//				cout << *test << endl;
				// check if transposition
				rrrmt* rho=test->simplify();
				rrrmts.push_back(rho);
				continue;
			}// END more that one maximal increasing substring
		} // END WHILE LOOP (temp != iota)
	}// if g=iota, then scenario is empty -> do nothing
}

///////////////////////////////////////////////////////////////////////////////

void test_sorting_algorithm(genom g){
	genom test=g;

	// CASE inversion
	cout << "Test all inversions ...";
	for(int i=0; i<test.size(); i++){
		for( int j=i; j<test.size();j++){
			rrrmt* inv;
			cout << "i=" << i << "; j=" << j << endl;
			inv = new rev(i, j, test);
			inv->apply(test);
			vector<rrrmt* > rhos;
			sort_by_4type(test, rhos);
			// test 1: does sequence rhos have length 1
			if(rhos.size()!=1){
				cerr << "Error: Sequence is longer than one for: " << endl;
				cerr << g << endl;
				cerr << test << endl;
				cerr << "length is: " << rhos.size() << endl;
				exit(EXIT_FAILURE);
			}
			// test 2: does the sequence transforms test into g?
			for(int k=rhos.size()-1; k>=0; k--){
				rhos[k]->apply(test);
			}
			if(test != g){
				cerr << "Error: Sequence does not sort test to g: " << endl;
				cerr << test << endl;
				for(int k=rhos.size()-1; k>=0; k--){
					cerr << "\t" << *rhos[k];
				}
				cerr << g << endl;
				exit(EXIT_FAILURE);
			}
			delete inv;
		}
	}
	cout << " success!" << endl;

	// CASE transposition
	cout << "Test all transpositions ...";
	for(int i=0; i < test.size(); i++){
		for( int j=i+1; j < test.size(); j++){
			for( int k=j+1; k < test.size(); k++){
				rrrmt* tra;
				if( (i != j) || (j != k)){
//					cout << "i=" <<  i << "; j=" << j << "; k=" << k << endl;
					tra = new transp(i, j, k, test);
					tra->apply(test);
					vector<rrrmt* > rhos;
					sort_by_4type(test, rhos);
					// test 1: does sequence rhos have length 1
					if(rhos.size()!=1){
						cerr << "Error: Sequence is longer than one for: " << endl;
						cerr << g << endl;
						cerr << test << endl;
						cerr << "length is: " << rhos.size() << endl;
						exit(EXIT_FAILURE);
					}
					// test 2: does the sequence transforms test into g?
					for(int k=rhos.size()-1; k>=0; k--){
						rhos[k]->apply(test);
					}
					if(test != g){
						cerr << "Error: Sequence does not sort test to g: " << endl;
						cerr << test << endl;
						for(int k=rhos.size()-1; k>=0; k--){
							cerr << "\t" << *rhos[k];
						}
						cerr << g << endl;
						exit(EXIT_FAILURE);
					}
					delete tra;
				}
			}
		}
	}
	cout << " success!" << endl;

	// CASE inverse transposition
	cout << "Test all inverse transpositions to the left ...";
	for(int i=0; i < test.size(); i++){
		for( int j=i; j < test.size(); j++){
			for( int k=j+1; k < test.size(); k++){
//				cout << "i=" <<  i << "; j=" << j << "; k=" << k << endl;
				rrrmt* itra;
				itra = new revtransp(i, j, j+1, k, test);
//				cout << *itra << endl;
				itra->apply(test);
				vector<rrrmt* > rhos;
				sort_by_4type(test, rhos);
				// test 1: does sequence rhos have length 1
				if(rhos.size()!=1){
					cerr << "Error: Sequence is longer than one for: " << endl;
					cerr << g << endl;
					cerr << test << endl;
					cerr << "length is: " << rhos.size() << endl;
					exit(EXIT_FAILURE);
				}
				// test 2: does the sequence transforms test into g?
				for(int k=rhos.size()-1; k>=0; k--){
					rhos[k]->apply(test);
				}
				if(test != g){
					cerr << "Error: Sequence does not sort test to g: " << endl;
					cerr << test << endl;
					for(int k=rhos.size()-1; k>=0; k--){
						cerr << "\t" << *rhos[k];
					}
					cerr << g << endl;
					exit(EXIT_FAILURE);
				}
				delete itra;
			}
		}
	}
	cout << " success!" << endl;

	cout << "Test all inverse transpositions to the right ...";
	for(int i=0; i < test.size(); i++){
		for( int j=i; j < test.size(); j++){
			for( int k=0; k < test.size(); k++){
				if((k<i)){
//					cout << "i=" <<  i << "; j=" << j << "; k=" << k << endl;
					rrrmt* itra;
					itra = new revtransp(i, j, k, i-1, test);
//					cout << *itra << endl;
					itra->apply(test);
					vector<rrrmt* > rhos;
					sort_by_4type(test, rhos);
					// test 1: does sequence rhos have length 1
					if(rhos.size()!=1){
						cerr << "Error: Sequence is longer than one for: " << endl;
						cerr << g << endl;
						cerr << test << endl;
						cerr << "length is: " << rhos.size() << endl;
						exit(EXIT_FAILURE);
					}
					// test 2: does the sequence transforms test into g?
					for(int k=rhos.size()-1; k>=0; k--){
						rhos[k]->apply(test);
					}
					if(test != g){
						cerr << "Error: Sequence does not sort test to g: " << endl;
						cerr << test << endl;
						for(int k=rhos.size()-1; k>=0; k--){
							cerr << "\t" << *rhos[k];
						}
						cerr << g << endl;
						exit(EXIT_FAILURE);
					}
					delete itra;
				}
			}
		}
	}
	cout << " success!" << endl;


	cout << "Test 10.000 random TDRLs ...";
	for(int i=0; i< 10000; i++){
		rrrmt* tdl= new tdrl(g, false);
		genom test1=g;
		genom test2=g;
		tdl->apply(test1);
		if(test1 != test2){
			vector<rrrmt* > rhos;
			sort_by_4type(test1, rhos);
			// test 1: does sequence rhos have length 1
			if(rhos.size()!=1){
				cerr << "Error: Sequence is longer than one for: " << endl;
				cerr << g << endl;
				cerr << test1 << endl;
				cerr << "length is: " << rhos.size() << endl;
				exit(EXIT_FAILURE);
			}
			// test 2: does the sequence transforms test into g?
			for(int k=rhos.size()-1; k>=0; k--){
				rhos[k]->apply(test2);
			}
			if(test1 != test2){
				cerr << "Error: Sequence does not sort test to g: " << endl;
				cerr << test1 << endl;
				for(int k=rhos.size()-1; k>=0; k--){
					cerr << "\t" << *rhos[k];
				}
				cerr << test2 << endl;
				exit(EXIT_FAILURE);
			}
			delete tdl;
		}
	}
	cout << " success!" << endl;

}

///////////////////////////////////////////////////////////////////////////////

void test_sorting_algorithm_pairwise(genom g){

	// CASE inversion
	cout << "Test all pairwise inversions ...";
	for(int i=0; i<g.size(); i++){
		for( int j=i; j<g.size();j++){
			genom test=g;
			rrrmt* inv;
//			cout << "i=" << i << "; j=" << j << endl;
			inv = new rev(i, j, test);
			inv->apply(test);
			rrrmt* rhos;
//			cout << "start: " << g << endl;
//			cout << "target: " << test << endl;
			sort_by_4type(g, test, rhos);
			// test 1: does sequence rhos have length 1
			if(rhos->length(ALTMIN)!=1){
				cerr << "Error: Sequence is unequal to one!" << endl;
				cerr << g << endl;
				cerr << test << endl;
				cerr << "length is: " << rhos->length(ALTMIN) << endl;
				exit(EXIT_FAILURE);
			}
			// test 2: does the sequence transforms test into g?
			genom test2=g;
			rhos->apply(test2);
			if(test != test2){
				cerr << "Error: Sequence does not sort test to g: " << endl;
				cerr << test2 << endl;
				cerr << "\t" << *rhos;
				cerr << test << endl;
				exit(EXIT_FAILURE);
			}
			delete inv;
		}
	}
	cout << " success!" << endl;

	// CASE transposition
	cout << "Test all pairwise transpositions ...";
	for(int i=0; i < g.size(); i++){
		for( int j=i+1; j < g.size(); j++){
			for( int k=j+1; k < g.size(); k++){
				genom test=g;
				genom test2=g;
				rrrmt* tra;
				if( (i != j) || (j != k)){
//						cout << "i=" <<  i << "; j=" << j << "; k=" << k << endl;
					tra = new transp(i, j, k, test);
					tra->apply(test);
					rrrmt* rhos;
					sort_by_4type(test2, test, rhos);
					// test 1: does sequence rhos have length 1
					if(rhos->length(ALTMIN)!=1){
						cerr << "Error: Sequence does not have length one! " << endl;
						cerr << test2 << endl;
						cerr << test << endl;
						cerr << "length is: " << rhos->length(ALTMIN) << endl;
						exit(EXIT_FAILURE);
					}
					// test 2: does the sequence transforms test into g?
					rhos->apply(test2);
					if(test != test2){
						cerr << "Error: Sequence does not sort test to g: " << endl;
						cerr << test2 << endl;
						cerr << "\t" << *rhos;
						cerr << test << endl;
						exit(EXIT_FAILURE);
					}
					delete tra;
				}
			}
		}
	}
	cout << " success!" << endl;

	// CASE inverse transposition
	cout << "Test all pairwise inverse transpositions to the left ...";
	for(int i=0; i < g.size(); i++){
		for( int j=i; j < g.size(); j++){
			for( int k=j+1; k < g.size(); k++){
				genom test=g;
				genom test2=g;
//				cout << "i=" <<  i << "; j=" << j << "; k=" << k << endl;
				rrrmt* itra;
				itra = new revtransp(i, j, j+1, k, test);
//				cout << *itra << endl;
				itra->apply(test);
				rrrmt* rhos;
				sort_by_4type(test2, test, rhos);
				// test 1: does sequence rhos have length 1
				if(rhos->length(ALTMIN)!=1){
					cerr << "Error: Sequence does not have lenght one! " << endl;
					cerr << test2 << endl;
					cerr << test << endl;
					cerr << "length is: " << rhos->length(ALTMIN) << endl;
					exit(EXIT_FAILURE);
				}
				// test 2: does the sequence transforms test into g?
				rhos->apply(test2);
				if(test != test2){
					cerr << "Error: Sequence does not sort test to g: " << endl;
					cerr << test2 << endl;
					cerr << "\t" << *rhos;
					cerr << test << endl;
					exit(EXIT_FAILURE);
				}
				delete itra;
			}
		}
	}
	cout << " success!" << endl;

	cout << "Test all pairwise inverse transpositions to the right ...";
	for(int i=0; i < g.size(); i++){
		for( int j=i; j < g.size(); j++){
			for( int k=0; k < g.size(); k++){
				if((k<i)){
					genom test=g;
					genom test2=g;
//					cout << "i=" <<  i << "; j=" << j << "; k=" << k << endl;
					rrrmt* itra;
					itra = new revtransp(i, j, k, i-1, test);
//					cout << *itra << endl;
					itra->apply(test);
					rrrmt* rhos;
					sort_by_4type(test2, test, rhos);
					// test 1: does sequence rhos have length 1
					if(rhos->length(ALTMIN)!=1){
						cerr << "Error: Sequence does not have length one! " << endl;
						cerr << test2 << endl;
						cerr << test << endl;
						cerr << "length is: " << rhos->length(ALTMIN) << endl;
						exit(EXIT_FAILURE);
					}
					// test 2: does the sequence transforms test into g?
					rhos->apply(test2);
					if(test != test2){
						cerr << "Error: Sequence does not sort test to g: " << endl;
						cerr << test2 << endl;
						cerr << "\t" << *rhos;
						cerr << test << endl;
						exit(EXIT_FAILURE);
					}
					delete itra;
				}
			}
		}
	}
	cout << " success!" << endl;

	cout << "Test 10.000 pairwise random TDRLs ...";
	for(int i=0; i< 10000; i++){
		rrrmt* tdl= new tdrl(g, false);
		genom test=g;
		genom test2=g;
		tdl->apply(test);
		if(test != test2){
			rrrmt* rhos;
			sort_by_4type(test2, test, rhos);
//			cout << "TEST2: " << test2 << endl;
//			cout << "\t" << *rhos << endl;
//			cout << "TEST: " << test << endl;
			// test 1: does sequence rhos have length 1
			if(rhos->length(ALTMIN)!=1){
				cerr << "Error: Sequence does not have length one! " << endl;
				cerr << test2 << endl;
				cerr << test << endl;
				cerr << "length is: " << rhos->length(ALTMIN) << endl;
				exit(EXIT_FAILURE);
			}
			// test 2: does the sequence transforms test into g?
			rhos->apply(test2);
			if(test != test2){
				cerr << "Error: Sequence does not sort test to g: " << endl;
				cerr << test2 << endl;
				cerr << "\t" << *rhos;
				cerr << test << endl;
				exit(EXIT_FAILURE);
			}
			delete tdl;
		}
	}
	cout << " success!" << endl;

}

///////////////////////////////////////////////////////////////////////////////

void oplus(set<int> &result, set<int> first, set<int> second){
	set<int>::iterator it;
	for(it=first.begin(); it!=first.end(); ++it){
		result.insert(*it);
	}
	for(it=second.begin(); it!=second.end(); ++it){
		result.insert(*it);
	}
}

