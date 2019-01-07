#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ostream>
#include <regex.h>
#include <string>
#include <vector>

#include "helpers.hpp"
#include "io.hpp"

//#define DEBUG_NONAMES

using namespace std;

int parseFilename(const char *filename, unsigned &m, unsigned &n, unsigned &k, unsigned &r_start, unsigned &r_end){
	int status;
	regex_t re_fname;
	regmatch_t match[6];
	string fname(filename);

	status = regcomp(&re_fname,"^[-_//.[:alnum:]]+[[:digit:]]+_m([[:digit:]]+)_n([[:digit:]]+)_k([[:digit:]]+)$",REG_ICASE|REG_EXTENDED);
	if( status != 0 ) {
		cout << "read_taxa: Could not compile regex pattern."<<endl;
		exit(-1);
	}

	if(regexec (&re_fname, &filename[0], 4, &match[0], 0)==0){
		string s;
		s.assign(fname.begin()+match[1].rm_so, fname.begin()+match[1].rm_so+match[1].rm_eo);
		m = string2int(s);
		s.assign(fname.begin()+match[2].rm_so, fname.begin()+match[3].rm_so+match[3].rm_eo);
		n = string2int(s);
		s.assign(fname.begin()+match[3].rm_so, fname.begin()+match[3].rm_so+match[3].rm_eo);
		k = string2int(s);
		r_start = 0;
		r_end = n;
		return 1;
	}
	regfree(&re_fname);

	status = regcomp(&re_fname,"^[-_//.[:alnum:]]+[[:digit:]]+_m([[:digit:]]+)_n([[:digit:]]+)_k([[:digit:]]+)_w([[:digit:]]+)-([[:digit:]]+)$",REG_ICASE|REG_EXTENDED);
	if( status != 0 ) {
		cout << "read_taxa: Could not compile regex pattern."<<endl;
		exit(-1);
	}

	if(regexec (&re_fname, &filename[0], 6, &match[0], 0)==0){
		string s;
		s.assign(fname.begin()+match[1].rm_so, fname.begin()+match[1].rm_so+match[1].rm_eo);
		m = string2int(s);
		s.assign(fname.begin()+match[2].rm_so, fname.begin()+match[3].rm_so+match[3].rm_eo);
		n = string2int(s);
		s.assign(fname.begin()+match[3].rm_so, fname.begin()+match[3].rm_so+match[3].rm_eo);
		k = string2int(s);

		s.assign(fname.begin()+match[4].rm_so, fname.begin()+match[4].rm_so+match[4].rm_eo);
		r_start = string2int(s);
		s.assign(fname.begin()+match[5].rm_so, fname.begin()+match[5].rm_so+match[5].rm_eo);
		r_end = string2int(s);
		return 1;
	}
	regfree(&re_fname);
	return 0;

}

///////////////////////////////////////////////////////////////////////////////

void print_name(string name, ostream &out, bool appdots){
	string::iterator end;

	for( string::iterator it = name.begin(); it!=name.end(); it++ ){
		if( *it == ',' ){
			if( appdots )
//				out << "...";
//			else
				out << "\\("<<count(name.begin(), name.end(), ',')+1 << "\\)";
			break;
		}
		out << *it;
	}
}

//void read_gb(const char *filename, int &length, vector<pair<int, int> > &genes, vector<int> &orientation){
//	int status,			// status flag for regex compile
//		start,
//		end;
//	ifstream gbfile;
//	string line,
//		s;
//	regex_t re_genecds,	// regular expression for a gene or cds
//		re_source,
//		re_complement;
//	regmatch_t match[4];
//
//		// open file
//	gbfile.open(filename);
//	if(!gbfile.is_open()){
//		cerr << "error: read_gb could not read " << filename<< endl;
//		exit(1);
//	}
//
//	status = regcomp(&re_genecds,"([[:digit:]]+)\\.\\.([[:digit:]]+)",REG_ICASE|REG_EXTENDED);
//	if( status != 0 ) {
//		cout << "read_gb: Could not compile regex pattern re_genecds."<<endl;
//		exit(-1);
//	}
//
//	status = regcomp(&re_source,"source",REG_ICASE|REG_EXTENDED);
//	if( status != 0 ) {
//		cout << "read_gb: Could not compile regex pattern re_source."<<endl;
//		exit(-1);
//	}
//
//	status = regcomp(&re_complement,"complement",REG_ICASE|REG_EXTENDED);
//	if( status != 0 ) {
//		cout << "read_gb: Could not compile regex pattern re_complement."<<endl;
//		exit(-1);
//	}
//
//	genes.clear();
//	orientation.clear();
//
//		// read the file line for line
//	while(gbfile){
//		getline(gbfile, line);
//		if(!gbfile)
//			break;
//
//		if( regexec (&re_source, &line[0], 1, &match[0], 0)==0 ){
//			if( regexec (&re_genecds, &line[0], 3, &match[0], 0) == 0){;
//				s.assign(line.begin()+match[2].rm_so, line.begin()+match[2].rm_eo);
//				length = atoi( s.c_str() );
//			}
//			continue;
//		}
//		if( regexec (&re_genecds, &line[0], 3, &match[0], 0)==0 ){
//			s.assign(line.begin()+match[1].rm_so, line.begin()+match[1].rm_eo);
//			start = atoi( s.c_str() );
//			s.assign(line.begin()+match[2].rm_so, line.begin()+match[2].rm_eo);
//			end = atoi( s.c_str() );
//
//			if(genes.size()>0 && start == genes.back().first && end == genes.back().second)
//				continue;
//
//			if( regexec (&re_complement, &line[0], 3, &match[0], 0)==0 )
//				orientation.push_back(1);
//			else
//				orientation.push_back(0);
//			genes.push_back( make_pair(start, end) );
//
//		}
//	}
//
//	gbfile.close();
//}

///////////////////////////////////////////////////////////////////////////////

void read_taxa(const string filename, vector<genom> &genomes, vector<string> &names, char circular, int normalize_to){

	string line;
	vector<string> data;

		// open file
	ifstream taxafile(filename.c_str());
	if(!taxafile){
		cout << "Could not read: "<< filename<<endl;
		exit(1);
	}

		// read the file line for line
	while(taxafile){
		getline(taxafile, line);
		if (!taxafile)
			break;
		data.push_back(line);

	}

	read_taxa(data, genomes, names, circular, normalize_to);

	taxafile.close();

}

///////////////////////////////////////////////////////////////////////////////

void read_taxonomy( const string &fname, vector<vector<string> > &tax ){
	ifstream taxafile;						// input file handle
	string line;

	taxafile.open(fname.c_str(), ios::in);
	if(!taxafile.is_open()){
		cerr << "read_genomes: could not read: "<< fname.c_str()<<endl;
		exit(1);
	}
	while(taxafile){
		vector<string> t;
		getline(taxafile, line);
		if (!taxafile)
			break;

		strip_leading_ws(line);
		strip_trailing_ws(line);

		if(line.size() == 0 || line[0] != ']')
			continue;

		line = line.substr(1,line.length());

		split_string(line, t, " ");

		tax.push_back(t);
	}
}

void read_genomes(const string &fname, vector<genom> &genomes, vector<string> &names, int circular,
		vector<string> &nmap, bool allowdel, bool allowdup, int normto){

	bool err = false;						// error flag
	ifstream taxafile;						// input file handle
	int s;									// sign
	map<string, int> nm,					// mapping from genenames to integers
		cnt;								// number of occ. of genes in one genome
	string line,							// one line of the text
		g;
	unsigned n = 1;							// element number
	vector<int> chr;						// one chromosome
	vector<string> tg;						// one text genome
	vector<vector<string> > txtgenomes;		// all the textgenomes from the file
	vector<vector<int> > signs;

//	cerr << "adel "<<allowdel<<endl;
//	cerr << "adup "<<allowdup<<endl;
		// open file
	taxafile.open(fname.c_str(), ios::in);
	if(!taxafile.is_open()){
		cerr << "read_genomes: could not read: "<< fname.c_str()<<endl;
		exit(1);
	}

	// recycle existing name map
	for(unsigned i=1; i<nmap.size(); i++){
//		if( nm.find(nmap[i]) == nm.end() ){
			nm[ nmap[i] ] = i;
//		}
	}
	n = nmap.size()+1;

		// read the file line for line
		// and save the names of the genomes as well as the genomes as text
	while(taxafile){

		getline(taxafile, line);

//		cerr <<"LINE "<< line << endl;
		if (!taxafile)
			break;

//		cout << line<< endl;
		strip_leading_ws(line);
		strip_trailing_ws(line);
			// skip empty or comment lines as well as lines containing trees
		if(line.size() == 0 || line[0] == '#' || line[0] == '(' || line[0] == ']')
			continue;

			// get the species name from 'name lines'
			// and start a new genome
		if(line[0] == '>'){
			line = line.erase(line.find_last_not_of(" ") + 1 ) ;
			names.push_back( line.substr(1,line.length()) );
			txtgenomes.push_back( vector<string>() );
			signs.push_back( vector<int>() );
		}
			// and the genomes from the other lines
			// add them to the last added genome (so multiline fasta should be possible)
		else{
			if( signs.size() == 0 ){
				cerr << "error: missing header"<<endl;
				exit(EXIT_FAILURE);
			}

			split_string(line, tg);
			for(unsigned i=0; i<tg.size(); i++){
				if( tg[i][0] == '-' ){
					s = -1;
					g = tg[i].substr(1);
				}else if( tg[i][0] == '+' ){
					s = 1;
					g = tg[i].substr(1);
				}else{
					s = 1;
					g = tg[i];
				}
				signs.back().push_back( s );
				txtgenomes.back().push_back( g );

				if( nm.find(g) == nm.end() ){
					nm[ g ] = n;
					n++;
				}
			}
			tg.clear();
		}
	}
	taxafile.close();

	// check if the number of name lines is equal to the number of genome lines
	if(names.size() != txtgenomes.size()){
		cout << "unequal number of names and genomes in the file"<<endl;
		cout << "names "<<endl;
		for(unsigned i=0; i<names.size(); i++){
		 	cout << i<< " : "<<names[i]<<endl;
		}
		cout << "genomes "<<endl;
		for(unsigned i=0; i<txtgenomes.size(); i++){
			cout << i<< " : "; copy(txtgenomes[i].begin(), txtgenomes[i].end(), ostream_iterator<string>(cout," ")); cout << endl;
		}
		exit(1);
	}

		// init the gene occurence counter map (needed for checking for
		// duplications and deletions)
	for( map<string,int>::iterator it=nm.begin(); it!=nm.end(); it++ ){
		cnt[ it->first ] = 0;
	}

	nmap.resize(n);
	nmap[0] = "INVALID_GENE_0";
	for( map<string,int>::iterator it=nm.begin(); it!=nm.end(); it++ ){
		nmap[it->second] = it->first;
	}

	for(unsigned i=0; i<txtgenomes.size(); i++){
		for(unsigned j=0; j<txtgenomes[i].size(); j++){
			chr.push_back( signs[i][j] * nm[ txtgenomes[i][j] ] );
			cnt[ txtgenomes[i][j] ]++;
		}
#ifdef DEBUG_NONAMES
		genomes.push_back( genom(chr, circular, NULL, normto) );
#else
		genomes.push_back( genom(chr, circular, &nmap, normto) );
#endif//DEBUG_NONAMES
		for( map<string, int>::iterator it = cnt.begin(); it!=cnt.end(); it++ ){
			if( allowdel == false && it->second == 0 ){
				cerr << names[i]<<": missing gene "<< it->first <<endl;
				err = true;
			}else if( allowdup == false && it->second > 1 ){
				cerr <<i<< names[i]<<": "<<it->second<<" gene "<<it->first <<endl;
				err = true;
			}
			cnt[it->first] = 0;
		}
		chr.clear();
	}
	cnt.clear();

//	cout << "nmap "<<nmap.size()<<endl;
//	copy(nmap.begin(), nmap.end(), ostream_iterator<string>(cout, " ")); cout << endl;
//
//	chr.resize(n+1, 0);
//	for( unsigned i=0; i<genomes.size(); i++ ){
//		chr.assign(n+1, 0);
//		for(unsigned j=0; j<genomes[i].size(); j++){
//			chr[abs(genomes[i][j])]++;
//		}
//
//		for(unsigned j=0; j<chr.size(); j++){
//			if(chr[j] != 1){
//				cerr << names[i]<<" : "<<chr[j] <<"x" << nmap[j]<<endl;
//				err = true;
//			}
//		}
//	}
	if( err ){
		cout << "errors occured while reading "<< fname << " -> exiting"<<endl;
		exit(1);
	}
		// strip leading and trailing white spaces from the species names
	for(unsigned i=0; i<names.size(); i++){
		strip_leading_ws(names[i]);
		strip_trailing_ws(names[i]);
	}

	txtgenomes.clear();
}

///////////////////////////////////////////////////////////////////////////////

void read_taxa( const vector<string> &data, vector<genom> &genomes, vector<string> &names, int circular, int normalize_to){
	int status;
	string line, name;
	regex_t re_gene;
	regmatch_t match_gene;
	genom tempGenome;		// gene to read in

	status = regcomp(&re_gene,"^([[:space:]]*-{0,1}[[:digit:]]+)",REG_ICASE|REG_EXTENDED);
	if( status != 0 ) {
		cout << "read_taxa: Could not compile regex pattern."<<endl;
		exit(-1);
	}

		// if neccesary reinit genomes
	for (unsigned i=0; i<genomes.size(); i++)
		genomes[i].clear();
	genomes.clear();
		// if neccesary reinit names
	for (unsigned i=0; i<names.size(); i++)
		names[i].clear();
	names.clear();

		// read the file line for line
		// lines beginning with # are comments
		// lines beginning with > are name lines and mark the beginning of a new genome
		// the other lines are the genomes
	for(unsigned i=0; i<data.size(); i++){
		line = data[i];
		if(line.size() == 0 || line[0] == '#')
			continue;
		if(line[0] == '>'){
			if(tempGenome.getChromosom().size()){
				tempGenome.setCircular(circular);
				genomes.push_back(tempGenome);
				names.push_back(name);
				tempGenome.clear();
				name.clear();
			}
			line = line.erase(line.find_last_not_of(" ") + 1 ) ;
			name = line.substr(1,line.length());
		}else{
			while(regexec (&re_gene, &line[0], 1, &match_gene, 0)==0){
				string gene_str;
				gene_str.assign(line.begin()+match_gene.rm_so, line.begin()+match_gene.rm_so+match_gene.rm_eo);
				line = line.substr(match_gene.rm_eo - match_gene.rm_so ,line.length());
				tempGenome.push_back(string2int(gene_str));
			}
		}
	}
	if(tempGenome.getChromosom().size()){
		tempGenome.setCircular(circular);
		genomes.push_back(tempGenome);
		names.push_back(name);
		tempGenome.clear();
		name.clear();
	}

	if(circular == 1){
		for(unsigned i=0; i< genomes.size(); i++){
			genomes[i].normalize( normalize_to );
		}
	}

		// strip leading and trailing white spaces from the names
	for(unsigned i=0; i<names.size(); i++){
		strip_leading_ws(names[i]);
		strip_trailing_ws(names[i]);
	}

	regfree(&re_gene);
}

///////////////////////////////////////////////////////////////////////////////

void read_trees( const string &fname, vector<string> &trees ){
	ifstream file;
	string line;

	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "error: read_trees could not read" << fname.c_str() << endl;
		exit(1);
	}

	while(! file.eof() ){
		getline(file, line);
		if(!file)
			break;

			// skip everything what does not look like a newick tree description
		if(line.size() == 0 || line[0] != '(')
			continue;

		trees.push_back( line );
	}

}

///////////////////////////////////////////////////////////////////////////////

vector<string> read_namemapping(string fname){
	ifstream f;
	int status;
	unsigned idx;
	regex_t re_map;
	regmatch_t match[3];
	string line,
		s;
	vector<string> m;

	f.open(fname.c_str());
	if( !f.is_open() ){
		cerr << "read_namemapping: can not open file: "<<fname <<endl;
		exit(1);
	}

	status = regcomp(&re_map,"^(.+) ([[:digit:]]+)",REG_ICASE|REG_EXTENDED);
	if( status != 0 ) {
		cout << "read_namemapping: Could not compile regex pattern."<<endl;
		exit(1);
	}



	while(f){
		getline(f, line);
		if (!f)
			break;
		if(regexec (&re_map, &line[0], 3, &match[0], 0) == 0){
			s.assign(line.begin()+match[2].rm_so, line.begin()+match[2].rm_eo);
			idx = atoi( s.c_str() );
			s.assign(line.begin()+match[1].rm_so, line.begin()+match[1].rm_eo);
			if( m.size() <= idx ){
				m.resize( idx+1 );
			}
			m[idx] = s;
		}

	}

	f.close();

	return m;
}

///////////////////////////////////////////////////////////////////////////////




//int read_grappaClusters(const char *filename, vector< pair<string,int> > &newick_trees, int &min_score){
//	return 1;
//	//~ ifstream grappa_out;									// file-stream for grappa-output
//	//~ string line,											// string for line-wise file reading
//		//~ s;														// string for regexp matches
//	//~ boost::regex re("^(.*);$"),						// regexp for newick-trees
//		//~ re_distance("(?:NJ )*inversion score  =\\s*(\\d+)$");	// regexp for inversion scores
//	//~ boost::smatch match;									// regexp-match
//	//~ // vector< pair<string, int> > trees;
//	//~ pair<string, int> tree;
//
//		//~ // reset the given parameters
//	//~ min_score=0;
//	//~ newick_trees.clear();
//		//~ // open the input-file
//	//~ grappa_out.open(filename);
//	//~ if(!grappa_out){
//		//~ cout << "Kann Inputdatei nicht lesen : " << filename << endl;
//		//~ return 0;
//	//~ }
//		//~ // read the grappa output and search for trees in the newick-tree
//		//~ // and inversion scores
//	//~ while(grappa_out){
//		//~ getline(grappa_out, line);
//		//~ if (line[0]=='#')
//			//~ continue;
//		//~ if (boost::regex_match(line, match, re)){
//			//~ s.assign(match[1].first, match[1].second);
//			//~ // newick_trees.push_back(s);
//			//~ tree.first = s;
//		//~ }
//		//~ if (boost::regex_match(line, match, re_distance)){
//			//~ s.assign(match[1].first, match[1].second);
//			//~ if(!min_score || min_score > string2int(s))
//				//~ min_score = string2int(s);
//			//~ tree.second = string2int(s);
//		//~ }
//		//~ if(tree.second && tree.first.size()){
//			//~ newick_trees.push_back(tree);
//			//~ tree.first.clear();
//			//~ tree.second = 0;
//		//~ }
//	//~ }
//		//~ // close the file
//	//~ grappa_out.close();
//		//~ // check if results
//	//~ if(!min_score || !newick_trees.size())
//		//~ return 0;
//	//~ else
//		//~ return 1;
//}

