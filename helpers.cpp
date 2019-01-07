#include <algorithm>
#include <limits>

#include "helpers.hpp"

using namespace std;


//void error(char *e, const char *ee = ""){
//	cerr << e << ee << endl;
//	exit(EXIT_FAILURE);
//}
///////////////////////////////////////////////////////////////////////////////

string output(vector<pair<unsigned, unsigned> > &v){
	string o;
	for (unsigned i=0; i<v.size(); i++)
		o += "[" + int2string(v[i].first) +","+ int2string(v[i].second) + "] ";
	return o;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++ counting ++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void counter_init( vector<int> &pos, unsigned size, unsigned min, bool inc ){
	pos.resize(size, 0);
	for( unsigned i=0; i<pos.size(); i++){
		if( inc ){
			pos[i] = min+i;
		}else{
			pos[i] = min;
		}
	}
}

void counter_invalidate( vector<int> &pos ){
	pos.assign( pos.size(), std::numeric_limits< int >::max() );
}

bool counter_valid( const vector<int> &pos ){
	if( pos.size() > 0 && pos.back() != std::numeric_limits< int >::max() ){
		return true;
	}else{
		return false;
	}
}

void counteradd( vector<int> &pos, int min, int max, unsigned inject, bool inc ){

	//cerr << "use of deprecated counteradd function"<<endl;

	int x;
	if( inject >= pos.size() ){
		pos.assign(pos.size(), std::numeric_limits< int >::max());
		return;
	}
	if( !counter_valid(pos) ){
		return;
	}

//	cout << "++counter"<<endl;
	for( int i = (int)pos.size()-1; i>=(int)inject; i-- ){
		pos[i]++;
//		copy(pos.begin(), pos.end(), ostream_iterator<int>(cout," ")); cout <<endl;
		if(inc){
			x = ((int)pos.size()-1-i );
		}else{
			x = 0;
		}
		if( pos[i] + x < max ){
//			cout << pos[i]<<"+"<<x<<"<"<<max<<endl;
			for(unsigned j=i; j<pos.size()-1; j++){
				if(inc){
					pos[j+1] = pos[j]+1;
				}else{
					pos[j+1] = min;
				}
			}
//			cout << "radj ";copy(pos.begin(), pos.end(), ostream_iterator<int>(cout," ")); cout <<endl;
			break;
		}
	}
//	cout << "++end "<<pos.back()<<" >= "<<max<<endl;
	if( pos.back() >= max ){
		counter_invalidate(pos);
	}
//	cout << "-> ";copy(pos.begin(), pos.end(), ostream_iterator<int>(cout," ")); cout <<endl;
}

//void counter( vector<int> &pos, int min, int max, unsigned inject, bool inc ){
//
//	if( inject >= pos.size() ){
//		pos.assign(pos.size(), std::numeric_limits< int >::max());
//		return;
//	}
//
//	pos[inject]++;
//	for(unsigned i=inject; i<pos.size(); i++){
//		if(pos[i] < max){
//			break;
//		}else{
//			pos[i] = min;
//			if(i+1 < pos.size())
//				pos[i+1]++;
//			else
//				pos.assign(pos.size(), std::numeric_limits< int >::max());
//		}
//	}
//	return;
//}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+ generalised counting ++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void countt( unsigned inject, const vector<int> &lengths, vector<int> &pos){

	if( inject >= lengths.size() ){
		counter_invalidate(pos);
		return;
	}

	pos[inject]++;
	for(unsigned i=inject; i<pos.size(); i++){
		if(pos[i] < lengths[i]){
			break;
		}else{
			pos[i] = 0;
			if(i+1 < pos.size())
				pos[i+1]++;
			else
				counter_invalidate(pos);
		}
	}
	return;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*Permutation ranking  ++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

unsigned rank_perm(unsigned n, vector<unsigned> &pi, vector<unsigned> &pi_inv){
	int s = pi[n-1];

	if(n==1)
		return 0;

	swap( pi[n-1], pi[ pi_inv[n-1] ] );
	swap( pi_inv[s], pi_inv[n-1]);

	unsigned ret;
	ret = s + n * rank_perm(n-1, pi, pi_inv);

	return ret;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* STRING CONVERTER FUNCTIONS  +++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int string2int(string str){
	istringstream ins;
	ins.str(str);
	int cv=0;
	ins >> cv;
	return cv;
}

string int2string(int i, int m){
	std::stringstream ss ;
	ss << i ;
	return ss.str();
}

string float2string(float f, int m){
	std::stringstream ss;
	ss.precision(m);
	ss << f ;
	return ss.str();
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* STRING FUNCTIONS  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void split_string(const string& str,
                      vector<string>& tokens,
                      const string& delimiters){
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void strip_leading_ws(string &str){
//	const char *whitespace = " \b\f\n\r\t\v";
//	size_t p = str.find_first_not_of(whitespace);
//	if( p!= string::npos )
//		str.erase(str.begin(), str.begin()+p);//Strip trailing whitespace

    str.erase(str.begin(), std::find_if(str.begin(), str.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
//    return str;


//	unsigned i=0;
//	for( i=0;i<=str.size()&&str.find_first_of(whitespace, i)==i&&str.find_first_of(whitespace, i)!=string::npos;i++);
//	str.erase(0,i);//Strip leading whitespace
}

void strip_trailing_ws(string &str){
//	const char *whitespace = " \b\f\n\r\t\v";
//	size_t p = str.find_last_not_of(whitespace);
//
//	if( p != string::npos )
//		str.erase(str.begin()+p+1, str.end());//Strip trailing whitespace

	//	unsigned i=0;
	//	for(i=str.size()-1;(i>=0)&&(str.find_last_of(whitespace, i)==i)&&(str.find_last_of(whitespace, i)!=string::npos);i--);
	//	str.erase(i+1,str.size()-(i+1));//Strip trailing whitespace

	str.erase(std::find_if(str.rbegin(), str.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), str.end());

}


string lower( const string &s){
	string cp(s.size(), ' ');
	for(unsigned i=0; i<s.size(); i++){
		cp[i] = tolower(s[i]);
	}
	return cp;
}

string upper( const string &s){
	string cp(s.size(), ' ');
	for(unsigned i=0; i<s.size(); i++){
		cp[i] = toupper(s[i]);
	}
	return cp;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* MATH FUNCTIONS    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

unsigned binom(unsigned n, unsigned r){
	if ( n < r ){
		return 0;
	}
	if  (r == 0){
		return 1;
	}

	if (r > (n - r)){
		return binom(n, n-r);
	}


	unsigned long b = n-r+1;
	for( unsigned i=2; i<=r; i++ ){
		if ( std::numeric_limits< unsigned long >::max() / (n-r+i) < b ){
			throw overflow_error("binom");
		}
		b *= (n-r+i);
		b /= i;
	}

	if ( b > std::numeric_limits< unsigned>::max() ){
		throw overflow_error("binom");
	}

	return b;
}


//vector<unsigned> pascal_row(unsigned n){
////       1        1
////      1 1       2
////     1 2 1      3
////    1 3 3 1     4
////   1 4 6 4 1    5
////  1 5 10105 1   6
////  ...
//
//	vector<unsigned> row(n+1, 1);
//
//	for(unsigned i=1; i<=n; i++){
////		cerr << row[i-1] << " * ("<<n+1-i<<"/"<<i<<")"<<endl;
//		row[i] = (unsigned)(row[i-1] * (float)(n+1-i)/(float)(i));
//	}
//	return row;
//}

vector<double> pascal_row(unsigned n){
	vector<double> row(n+1, 1);

	for(unsigned i=1; i<=n; i++){
//		cerr << row[i-1] << " * ("<<n+1-i<<"/"<<i<<")"<<endl;
		row[i] = (row[i-1] * (double)(n+1-i)/(double)(i));
	}
	return row;
}


unsigned digits(const unsigned number){
	unsigned d=1, q=10;

	while(number/q!=0){
		d++;
		q *= 10;
	}

	return d;
}

unsigned fact(unsigned a){
	unsigned ar = 1;

	if (a == 0)
		return 1;

	for (unsigned i=2; i <= a; i++){
		if ( std::numeric_limits< unsigned >::max() / i < ar ){
			return 0;
		}
		ar *= i;

	}
	return ar;
}

int pow(int b, int p){
	int ret = b;

	if (!p)
		return 1;

	for (int i=1; i<p; i++){
		ret *= b;
	}

	return ret;
}

int ppow(int p){
	return 1 << p;
}

int sign(const int &a) {
	return (a == 0) ? 0 : (a<0 ? -1 : 1);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* BIT_VECTOR HANDLING FUNNCTIONS  +++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool conjunction(const vector<bool> &b, unsigned maxSize){
		// if not all elements are allocated jet
		// means that those are false
	if(b.size() != maxSize){
		return false;
	}
		// check if there is a false -> then the result is false also
	for(unsigned i=0; i<b.size(); i++){
		if(!b[i])
			return false;
	}
	return true;
}

bool disjunction(const vector<bool> &b){
	for(unsigned i=0; i<b.size(); i++){
		if(b[i])
			return true;
	}
	return false;
}

void mark(vector<bool> &b, unsigned pos, bool val){
	//~ cout << "mark "<<pos << " with "<<val;

	if(pos >= b.size())
		b.resize(pos+1);
	b[pos] = val;
}

bool isMarked(const vector<bool> &b, unsigned pos){
	if(pos < b.size()){
		if(b[pos])
			return true;
	}
	return false;
}

unsigned sum(const vector<bool> &b){
	unsigned s=0;
	for(unsigned i=0; i<b.size(); i++){
		if(b[i])
			s++;
	}
	return s;
}

void print(const vector<bool> &b){
	for(unsigned i=0; i<b.size(); i++){
		if(b[i])
			cout << "1";
		else
			cout << "0";
	}
	//~ cout << endl;
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* STACK HANDLING FUNNCTIONS  ++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void outputStack(stack<int> &s){
	stack<int> temp = s;

	while(temp.size()){
		cout << temp.top()<< " ";
		temp.pop();
	}
	cout << endl;
	return;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* RANDOM NUMBER GENERATOR FUNKTIONS   ++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void init_rng(unsigned init){
// man drand48:
// """
//	NOTES
//       These functions are declared obsolete by SVID 3, which states that rand(3) should be used instead.
// """
#ifdef _GLIBCPP_HAVE_DRAND48
	srand48(init);
//	cout << "SRAND48"<<endl;
#else
//	cout << "SRAND"<<endl;
	srand(init);
#endif

}

int ask_rng(){
#ifdef _GLIBCPP_HAVE_DRAND48
	return drand48();
#else
	return rand();
#endif
}

float ask_rng_f(){
#ifdef _GLIBCPP_HAVE_DRAND48
	return lrand48();
#else
	return ((float) rand()/(float)RAND_MAX);
#endif
}

int ask_rng( int s, int e ){
	return (ask_rng()%(e-s+1)) + s;
}

vector<unsigned> rng_inc_seq( unsigned l, unsigned m, unsigned d ){
	set<unsigned> seqset;	// a set containing the seq elements
	vector<unsigned> seq;	// the seq to be constructed
	set<unsigned>::const_iterator mx, mn; 		// two temporary variables for getting min and max
	unsigned max, min;

	if( m == 0 && l > 1 ){
		cerr << "rng_inc_seq() impossible increasing sequence m=0 and l="<<l<<endl;
		exit(EXIT_FAILURE);
	}
	if( d < l ){
		cerr << "rng_inc_seq() impossible increasing sequence d="<<d<<" and l="<<l<<endl;
		exit(EXIT_FAILURE);
	}
//	cout << "incseq "<<l<<" "<<m<<" "<<d<<endl;

	do{
		seqset.clear();
		while( seqset.size() < l ){
			seqset.insert( ask_rng() % m );
		}
		mx = seqset.end(); mx--;
		mn = seqset.begin();
		max = *mx;
		min = *mn;
	}while( (max-min) > d );
//	copy(seqset.begin(), seqset.end(), ostream_iterator<int>(cout, " ")); cout << endl;
//	cout<<"d " << *max_element(seqset.begin(), seqset.end())-*min_element(seqset.begin(), seqset.end())<<" = "<<*max_element(seqset.begin(), seqset.end())<<"-"<<*min_element(seqset.begin(), seqset.end())<<endl;
	seq.insert( seq.begin(), seqset.begin(), seqset.end() );
	return seq;
}

void weighted_choice_init( vector<float> &prob){
	for(unsigned i=0; i<prob.size(); i++){
		if(i>0)
			prob[i] += prob[i-1];
	}
}

/**
 * choose a random integer between 0 and prob.size()-1 given different 'probabilities'
 * of each integer
 * @param[in] prob probabilities of each integer (index)
 * @param[in] probsum sum of the probabilities (compute it with random_choice_init())
 */
int weighted_choice( const vector<float> &prob){
	float r;
	r = ask_rng_f();	// rand float in [0:1]

	//	cout << "r"<<r<<endl;
	for( unsigned i=0; i<prob.size(); i++){
		if( r <= prob[i] / prob.back() )
			return i;
	}
	return prob.size()-1;

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* TIME FUNNCTIONS   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double get_time(){
#ifdef USEMPI
	return MPI::Wtime();
#else
	struct timespec t;
	#if _POSIX_TIMERS > 0
	clock_gettime(CLOCK_REALTIME, &t);
	#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	t.tv_sec = tv.tv_sec;
	t.tv_nsec = tv.tv_usec*1000;
	#endif
	return (t.tv_sec) + (t.tv_nsec)/1000000000.0;
#endif//USEMPI
}

bool scoreTableCmp(const vector<unsigned> &a, const vector<unsigned> &b){
	return a[0] < b[0];
}

#ifdef LAMMPI

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* GENOME VECTOR HANDLING FUNNCTIONS +++++++++++++++++++++++++++++++++++++++ */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void genomVector2intVector(const vector<genom> &gv, vector<int> &iv, const unsigned &n){
	iv.clear();
	iv.resize(gv.size()*n);

	//~ vector<int>::iterator curIvPos=iv.begin();
	unsigned curOffset = 0;

	for(unsigned i=0; i<gv.size(); i++){
		vector<int> c = gv[i].getChromosom();
		//~ iv.insert(curIvPos, c.begin(), c.end());
		memcpy(&iv[curOffset], &c[0], n*sizeof(int));
		//~ curIvPos += n;
		curOffset += n;
	}
}

void intVector2genomVector(const vector<int> &iv, vector<genom> &gv, const unsigned &n, const char &circ){
	gv.clear();
	gv.reserve(iv.size()/n);

	genom temp(n, circ);
	vector<int> c(n);
	int curStart = 0;

	for(unsigned i=0; i<iv.size()/n; i++){
		memcpy(&c[0], &iv[curStart], n*sizeof(int));

		//~ c.insert(c.begin(), iv.begin()+curStart, iv.begin()+curStart+n);
		temp.setChromosom(c);
		gv.push_back(temp);
		curStart += n;
		//~ c.clear();
	}
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* MPI HELPER & WRAPPER FUNNCTIONS +++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void SendGenomeVector(const vector<genom> &g, unsigned from, unsigned to,
				vector<int> &i, unsigned &uBuf, const unsigned &n,
				int dest, MPI::Request &req){
	if(to > g.size())
		to = g.size();

	vector<genom> part;
	part.insert(part.begin(), g.begin()+from, g.begin()+to);

		// construct an integer vector from the genom vector
	genomVector2intVector(part, i, n);
		// transfer the size to the dest
	uBuf = i.size();
	MPI::COMM_WORLD.Isend( &uBuf, 1, MPI::UNSIGNED, dest, 0);
	if(uBuf>0){
			// transfer the integer vector to the dest
		req = MPI::COMM_WORLD.Isend( &i[0], i.size(), MPI::INT, dest, 2);
	}
}

//~ void SendGenomeMap(const map<genom, vector<bool> > &gm, vector<genom> &gv, vector<int> &i, const unsigned &n, int dest, MPI::Request &req){
	//~ gv.clear();
	//~ gv.reserve(gm.size());

	//~ for(map<genom, vector<bool> >::const_iterator it=gm.begin(); it!=gm.end(); it++){
		//~ gv.push_back(it->first);
	//~ }
	//~ SendGenomeVector(gv, 0, gv.size(), i, n, dest, req);
//~ }

void RecvGenomeVector(vector<int> &i, unsigned size, const unsigned &n, int source, MPI::Request &req){
		// receive the size of the comming integervector
	//~ int size=0;
	//~ MPI::COMM_WORLD.Irecv( &size, 1, MPI::INT, source, 0);
	i.resize(size);
	if(size>0){
			// make the receive buffer large enough and receive the int vector
		req = MPI::COMM_WORLD.Irecv( &i[0], size, MPI::INT, source, 2);
			// transform the int vector to a genom vector
		//~ intVector2genomVector(i, g, n, circ);
		//~ cout << "received "<<g.size();
	}
}

void BcastGenomeVector(vector<genom> &g, char circ, vector<int> &i, const unsigned &n){

		// the root node has to construct the integer vector
	int size=0;

	if(MPI::COMM_WORLD.Get_rank() == 0){
		genomVector2intVector(g, i, n);
		size = i.size();
		//~ for(unsigned i=0; i< g.size(); i++)
			//~ cout << "bi"<< i<<" c "<<(int)g[i].getCircular()<<endl;
	}

		// broadcast the size of the following message
	//~ MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Bcast( &size, 1, MPI::INT, 0);
		// the receiving nodes have to resize the receive buffer
	if(MPI::COMM_WORLD.Get_rank() != 0){
		i.resize(size);
	}
		// receive

	MPI::COMM_WORLD.Bcast( &i[0], i.size(), MPI_INT, 0);

	// get the genom vector
	if(MPI::COMM_WORLD.Get_rank() != 0){
		intVector2genomVector(i, g, n, circ);
		//~ for(unsigned i=0; i< g.size(); i++)
			//~ cout << "ai"<< i<<" c "<<(int)g[i].getCircular()<<endl;
	}
}

void GathervGenomeVector(vector<genom> &g, char circ, vector<int> &i, const unsigned &n){
	vector<int> sizes(MPI::COMM_WORLD.Get_size(),0);
	vector<int> displ(MPI::COMM_WORLD.Get_size(),0);
	int size=0;

	vector<int> ii;

		// construct the int-vector on every node and get its size
	genomVector2intVector(g, ii, n);
	size = ii.size();

		// gather the size of each intvector on node 0 in the sizes vector
	MPI::COMM_WORLD.Gather(&size, 1, MPI::INT, &sizes[0], 1, MPI_INT, 0);
		// let node 0 calculate the sum of the sizes, allocate enough memory
		// and calculate the displacements
	if(MPI::COMM_WORLD.Get_rank() == 0){
		unsigned rsize=0;
		for(unsigned k=0; k<sizes.size(); k++){
			//~ cout << "k "<<k<<"  size "<<sizes[k]<<endl;
			rsize += sizes[k];
			if(k<sizes.size()-1){
				displ[k+1] = rsize;
			}
		}
		i.resize(rsize);
	}

		// gather the int vector on node 0
	MPI::COMM_WORLD.Gatherv(&ii[0], size, MPI_INT, &i[0], &sizes[0], &displ[0], MPI_INT, 0);

		// get the genome vectors
	if(MPI::COMM_WORLD.Get_rank() == 0){
		intVector2genomVector(i, g, n, circ);
	}
}

void GathervVector(vector<unsigned> &u, const int &proc_size, const int &proc_rank){
	vector<int> sizes(proc_size,0);
	vector<int> displ(proc_size,0);
	vector<unsigned> uu;
	int size = u.size();

	MPI::COMM_WORLD.Gather(&size, 1, MPI::INT, &sizes[0], 1, MPI_INT, 0);

	if(proc_rank == 0){
		unsigned rsize=0;
		for(unsigned k=0; k<sizes.size(); k++){
			rsize += sizes[k];
			if(k<sizes.size()-1){
				displ[k+1] = rsize;
			}
		}
		uu.resize(rsize);
	}
	MPI::COMM_WORLD.Gatherv(&u[0], size, MPI::UNSIGNED, &uu[0], &sizes[0], &displ[0], MPI::UNSIGNED, 0);

	if(proc_rank==0){
		u = uu;
	}
}

void mpiSpeedTest(int MPI_rank, int MPI_size, unsigned n, unsigned size){
	double start_time = 0, end_time = 0;
	vector<int> testI;
	vector<genom> testG, testG2;

	if(MPI_rank == 0){
		for (unsigned i=0; i< size; i++){
			genom tmpGenom = genom(n, 5, 0, 0, true, true, 0);
			testG.push_back(tmpGenom);
		}
	}

	if(MPI_rank == 0){
		cout << "broadcast test "<< testG.size()<<" genomes "<<endl;
	}
	MPI::COMM_WORLD.Barrier();

	start_time = MPI::Wtime();
	for(unsigned i=0; i< 1000; i++){
		BcastGenomeVector(testG, 1, testI, n);
	}
	end_time = MPI::Wtime();
	cout << MPI_rank<<" btime "<<end_time - start_time<<endl;


	//~ start_time = MPI::Wtime();
	//~ XMPI_Buoy("G");
	//~ GathervGenomeVector(temp, testI, n);
	//~ end_time = MPI::Wtime();
	//}
	//~ cout << MPI_rank<<" gtime "<<end_time - start_time<<endl;

	//~ cout <<MPI_rank<<" : " << testG.size()<<endl;

}
#endif

