/**
 * program to calculate a parsimonious (w.r.t. given weighting), preserving scenario rearranging one given gene order into another using
 * 			Inversions, Transpositions, inverse Transpositions, and Tandem duplication random loss.
 * 			Combination of DP + ILP
 * @author: T. Hartmann + M. Bernt
 */

#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <queue>
#include <set>
#include <string>
#include <cstddef>

#include "crex.hpp"

#include "genom.hpp"
#include "costfoo.hpp"
#include "helpers.hpp"
#include "io.hpp"
#include "rrrmtrand.hpp"
#include "sb4type.hpp"

//#define DEBUG_CREX
//#define DEBUG_CREX2

using namespace std;

/**
 * get program parameters
 * @param[in] argc argument count
 * @param[in] argv arguments
 * @param[in/out] fname filename
 * @param[in/out] circluar circularity
 * @param[in/out] mkalt create alternatives for transpositions and inverse transpositions
 * @param[in/out] bps compute breakpoint scenario for the given perm
 * @param[in/out] crexone compute crex1 scenario for the given perm
 * @param[in/out] crex add breakpoint scenario as alternative for prime nodes
 * @param[in/out] maxalt maximum number of alternatives constructed for prime nodes (inv+tdrl, bp)
 * @param[in/out] wI weight of inversions
 * @param[in/out] wT weight of transpositions
 * @param[in/out] wiT weight of inverse transpositions
 * @param[in/out] wTDRL weight of TDRLs
 * @param[in/out] distance: print only distance table
 */
void getoptions(int argc, char *argv[], string &fname, int &circular, bool &mkalt,
		bool &bps, bool &crexone, bool &crexbps, unsigned &maxalt, float &wI, float &wT, float &wiT, float &wTDRL,
		bool &distance);

/**
 * print usage info and exit
 */
void usage();

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char *argv[]){
	bool bps = false,		// compute rearrangements with breakpoint method from [ZhaoBourque07]
		crexone = false,	// use crex1 to compute rearrangements
		crexbps = true,		// use breakpoint method for prime nodes in crex
		mkalt = false;		// generate alternatives
	costfoo_by_type *cst; 		// cost function for crex2
	int circular = true;		// are the genomes circular?
	unsigned maxszen = 2,		// maximal number of reversals for combined reversal+tdrl scenarios
		time_bound = std::numeric_limits<unsigned>::max();
	string filename;		// input filename
	vector<genom> genomes; 		// the genomes
	vector<float> weights(4,1);	// weights for reconstruction
	vector<string> names, 		// names of the genomes
		nmap;			// names of the elements
	bool use_approx_alg=true;
	bool dist=false;
	unsigned distance=0;
	vector<vector<string> > tax;

	// get and check the parameters
	getoptions( argc, argv, filename, circular, mkalt, bps, crexone, crexbps, maxszen, weights[0], weights[1], weights[2], weights[3], dist);

	cst = new costfoo_by_type();
	cst->set( REVNM, weights[0] );
	cst->set( TRANM, weights[1] );
	cst->set( RTRANM, weights[2] );
	cst->set( "tdrl", weights[3] );


	if(filename == ""){
		cout << "CREx2 needs a gene order file!" << endl << endl;
		usage();
	}else{
		vector<vector<rrrmt *> > cm;	// pairwise crex scenarios

		read_taxonomy(filename, tax);
		read_genomes(filename, genomes, names, circular, nmap, true);	// read the genom file
		if(genomes.size() == 0){
			cout << "no genomes in the file"<<endl;
			usage();
		}

		cm = vector<vector<rrrmt*> >(genomes.size(), vector<rrrmt*>(genomes.size(), NULL));
		for( unsigned i=0; i<genomes.size(); i++ ){
			for(unsigned j=0; j<genomes.size(); j++ ){
				bool is_optimal=true;
				if(i==j){
					continue;
				}
				if( crexone ){

					cm[i][j] = new crex(genomes[i], genomes[j], crexbps, maxszen, mkalt);
			
					//print only first name of names[i] and names[j]
					string i_name;
					size_t found = names[i].find_first_of(",");
					if(found==std::string::npos){
						i_name=names[i];
					}else{
						i_name=names[i].substr(0,found);
					}
					string j_name;
					found = names[j].find_first_of(",");
					if(found==std::string::npos){
						j_name=names[j];
					}else{
						j_name=names[j].substr(0,found);
					}

					if(dist){ //solve distance problem
						cout << names[i] << "\t" << names[j] << "\t" << cm[i][j]->length(ALTMIN) << endl;
					}else{// solve sorting problem
						cout << names[i] << " -> " << names[j] << ":" << endl;
						cout << "scenario: " << (*cm[i][j]) << endl;
					}
				}else{

					// calulate crex2 solution
					cm[i][j] = new crex2(genomes[i], genomes[j], true, cst, is_optimal, false, mkalt, time_bound, distance);

					//only first name of names[i] and names[j]
					string i_name;
					size_t found = names[i].find_first_of(",");
					if(found==std::string::npos){
						i_name=names[i];
					}else{
						i_name=names[i].substr(0,found);
					}
					string j_name;
					found = names[j].find_first_of(",");
					if(found==std::string::npos){
						j_name=names[j];
					}else{
						j_name=names[j].substr(0,found);
					}

					if(cm[i][j]->length(ALTMIN)!=0){
						if(dist){//solve distance problem only
							cout << i_name << "\t" << j_name << "\t" << cm[i][j]->length(ALTMIN) << endl;
						}else{// solve sorting problem
							cout << i_name << " --> " << j_name << ":" << endl;
							cout << genomes[i] << endl;
							cout << "scenario: " << (*cm[i][j]);
							cout << genomes[j] << endl;
						}
						// TEST
						genom source=genomes[i];
						genom target=genomes[j];
						cm[i][j]->apply(source);
						if(source != target){
							cerr << "SCENARIO IS NO SORTING SCENARIO ---> exiting!" << endl;
							cerr << "source= " << genomes[i] << endl;
							cerr << "reached genome= " << source << endl;
							cerr << "target= " << target << endl;
							exit(EXIT_FAILURE);
						}
					}else{
						if(dist){//solve distance problem only
							cout << i_name << "\t" << j_name << "\t" << 0 << endl;
						}else{// solve sorting problem
							cout << i_name << " == " << j_name << endl;
						}
					}
				}
			}
		}

		for(unsigned i=0; i<cm.size(); i++){
			for(unsigned j=0; j<cm[i].size(); j++){
				delete cm[i][j];
			}
		}
		cm.clear();
	}
	delete cst;

	return 0;
}
void getoptions(int argc, char *argv[], string &fname, int &circular, bool &mkalt,
		bool &bps, bool &crexone, bool &crexbps, unsigned &maxalt, 
		float &wI, float &wT, float &wiT, float &wTDRL, bool &distance){

	int c;

	while (1) {
		static struct option long_options[] = {
			{"noalt",   no_argument, 0, 'a'},
			{"bp",      no_argument, 0, 'b'},
			{"crex1",   no_argument, 0, 'c'},
			{"distance", no_argument, 0, 'd'},
			{"file",    required_argument, 0, 'f'},
			{"linear",  no_argument, 0, 'l'},
			{"primxalt",required_argument, 0, 'm'},
			{"prinobp", no_argument, 0, 'o'},
			{"wI", required_argument, 0, 'I'},
			{"wT", required_argument, 0, 'T'},
			{"wiT", required_argument, 0, 'N'},
			{"wTDRL", required_argument, 0, 'D'},
			{"help", no_argument, 0, 'h'},
			{0, 0, 0, 0}
		};
        int option_index = 0;			// getopt_long stores the option index here.
        c = getopt_long (argc, argv, "abcdD:f:hI:lmN:oT:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
        	break;
        switch (c){
			case 0: // If this option set a flag, do nothing else now.
				if (long_options[option_index].flag != 0)
					break;
				cout << "option "<< long_options[option_index].name;
				if (optarg)
					cout << " with arg " << optarg;
				cout << endl;
				break;
			case 'a':
				mkalt = false;
				break;
			case 'b':
				bps = true;
				break;
			case 'c':
				crexone = true;
				break;
			case 'd':
				distance = true;
				break;
			case 'D':
				wTDRL = atof( optarg );
				break;
			case 'f':
				fname = optarg;
				break;
			case 'h':
				usage();
				exit(EXIT_FAILURE);
			case 'I':
				wI = atof( optarg );
				break;
			case 'l':
				circular = false;
				break;
			case 'm':
				maxalt = atoi( optarg );
				break;
			case 'N':
				wiT = atof( optarg );
				break;
			case 'o':
				crexbps = false;
				break;
			case 'T':
				wT = atof( optarg );
				break;
			case '?':
				exit(EXIT_FAILURE);
				break; /* getopt_long already printed an error message. */
			default:
				usage();
        }
	}

	if( crexone and bps ){
		cerr << "the options --crex1 and --bp must not be set simultaneously"<<endl;
		exit( EXIT_FAILURE );
	}
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void usage(){
cout << "CREx2 is an algorithm that aims to compute a weight-minimum sequence of rearrangements for arbitrary mitochondrial gene orders and the following types of weighted rearrangement operations: inversions, transpositions, inverse transpositions, and tandem duplication random loss. CREx2 considers only rearrangement operations that preserve common intervals, i.e., groups of genes that form an interval in both given gene orders. If the common intervals of a problem instance are organized in a linear structure, then CREx2 has a linear runtime. Otherwise, there are two modes for CREx2. The first mode, called CREx2-ILP, computes an exact weight-minimum rearrangement scenario within an exponential runtime (in worst case). The second mode, called CREx2-APP, computes approximated solutions efficiently." << endl << endl;

cout << "Caution: This implementation does not contain the first mode of CREx2! Therefore, only approximated solutions can be computed by either CREx or CREx2-APP." << endl<< endl;

	cout << "###### CREx2 usage ######"<<endl;
	cout << "call:" << endl;
	cout << "./crex2 -f [gene_order_file] [OPTIONS]" << endl;
	cout << endl;
	cout << "general options:"<<endl;
	cout << "--crex1:	-c: compute with CREx1 (default CREx2)"<<endl;
	cout << "--distance	-d: calculate genomic distances/weights only "<< endl;
	cout << "--file  -f: specify a filename (mandatory)"<<endl;
	cout << "--linear 	-l: handle genomes as linear (default: circular)"<<endl;
	cout << endl;
	cout << "CREx2 options:"<<endl;
	cout << "--wI 		-I: weight of an inversion"<<endl;
	cout << "--wT 		-T: weight of a transposition"<<endl;
	cout << "--wiT 		-N: weight of an inverse transposition"<<endl;
	cout << "--wTDRL		-D: weight of a TDRL"<<endl;
	cout << endl;
	cout << "CREx1 options:"<<endl;
	cout << "--bp:		-b: compute with breakpoint scenario [ZhaoBourque07]"<<endl;
	cout << "--noalt:	-a: don't compute alternatives for T+iT"<<endl;
	cout << "--prinobp:	-m: don't construct breakpoint scenario for prime nodes"<<endl;
	cout << "--primxalt:	-o: maximal number of alternatives for prime nodes (I+TDRL, bp) (default: 2)"<<endl;
	cout << endl;
	exit(1);
}
