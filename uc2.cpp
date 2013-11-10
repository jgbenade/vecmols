// This file is part of BOINC.
// http://boinc.berkeley.edu
// Copyright (C) 2008 University of California
//
// BOINC is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// BOINC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with BOINC.  If not, see <http://www.gnu.org/licenses/>.

// This program serves as both
// - An example BOINC application, illustrating the use of the BOINC API
// - A program for testing various features of BOINC
//
// NOTE: this file exists as both
// boinc/apps/upper_case.cpp
// and
// boinc_samples/example_app/uc2.cpp
// If you update one, please update the other!

// The program converts a mixed-case file to upper case:
// read "in", convert to upper case, write to "out"
//
// command line options
// --run_slow: sleep 1 second after each character
// --cpu_time N: use about N CPU seconds after copying files
// --early_exit: exit(10) after 30 chars
// --early_crash: crash after 30 chars
// --trickle_up: sent a trickle-up message
// --trickle_down: receive a trickle-up message
//

#ifdef _WIN32
#include "boinc_win.h"
#else
#include "config.h"
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <time.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#endif

#include "str_util.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"
#include "graphics2.h"


#ifdef APP_GRAPHICS
#include "uc2.h"
UC_SHMEM* shmem;
#endif

using std::string;
using namespace std;

typedef vector<int> permutation ;
typedef vector<permutation> lsquare;

#define CHECKPOINT_FILE "vecmols_state"
#define INPUT_FILENAME "in"
#define OUTPUT_FILENAME "out"

bool run_slow = false;
bool early_exit = false;
bool early_crash = false;
bool early_sleep = false;
bool trickle_up = false;
bool trickle_down = false;
bool critical_section = false;    // run most of the time in a critical section
double cpu_time =2, comp_result;
double start_time, end_time;

//MFILE out;
//FILE* state, *infile;

class MOLS {
public:
	typedef vector<int> permutation ;
	typedef vector<permutation> lsquare;
	MOLS(int n, int k) ;
	MOLS(FILE* infile, MFILE out);
	MOLS(FILE* infile, char* output_p);
	~MOLS();
	int enumerateMOLS();
	void  printAllStatics();
	vector<vector<vector< int> > > sqSymCurrUni;
	int restart(FILE* state,bool t);
	void printMOLSPerms( vector<lsquare> mols  , FILE* outfile );

	struct My
	{
	    static bool permutationComparator (permutation p1, permutation    p2);
	};

private:
	int n;
	int k;
	vector<lsquare> partMOLS;
	vector<lsquare> tempMOLS;
	vector<lsquare> root;
	vector<vector<vector<int> > > currLS;
	vector<permutation> cycleStructureReps;
	char* output_path;
	//MFILE out;
	//vector<vector <lsquare> > completedMOLS;

	int currSquare;
	int count_MOLS;
	vector<long long int> branchCount_  ;
	bool printOut;
	vector<permutation> possibleShuffles;
	vector<int> currentCS;
	string filename;
	vector<vector<int> > RCSsquares;
	vector<vector<int> > RCSperms;
	vector<int> z, diffs, sizes;
	vector<int> identity;
	vector<int> identityIt;
	vector<vector<vector<permutation> > > sqSymPossPerms;
	//vector<vector<long long int> > detailedCount;
	vector<int> numIsSmallest;
	vector<int>  numIsSmallestTrue;
	clock_t dtime, prev_time;
	int send_checkpoint; //0 means no, 1 means next time, 2 means done up to limit, 3 means written out, done checkpointing
	long nof_calls_made;
	const static long nof_calls_limit = 5000000000; // 5.10e5 is so .7 sek, 5.10e8 is 300sek, 5.109 so 3100 sek (op 83.9)
	//const static long nof_calls_limit = 50000000;

	int printUniversals(lsquare &l);
	void printPerm(permutation p);
	void printLS(lsquare &l);
	void printMOLSPerms( vector<lsquare> mols   );

	void printMOLS( vector<lsquare> mols  );
	bool isOrthogonal(permutation perm);
	permutation composition(permutation  &p1, permutation &p2);
	permutation inverse(permutation &p);
	permutation copyPerm(permutation p);
	permutation rcs(permutation &p1, permutation &p2);
	void permuteSymbols(vector<lsquare> &partMOLS , vector<lsquare> &newMOLS );
	void permuteMOLS(vector<lsquare> &partMOLS , permutation &rowPerm, permutation &colPerm  );
	bool testPermuteMOLS(vector<lsquare> &partMOLS , permutation &rowPerm, permutation &colPerm  );
	//void permuteMOLSRCS(vector<lsquare> partMOLS , permutation rowPerm, permutation colPerm, vector<lsquare> newMOLS );
	//void testpermuteMOLS(vector<lsquare> partMOLS , permutation rowPerm, permutation colPerm, vector<lsquare> newMOLS , vector<lsquare> newpMOLS);
	void printOA(vector<vector<int> > OA);
	bool  comparePartialMOLS(vector<lsquare> &partMOLS1, vector<lsquare> &partMOLS2);
	void  MOLStoOA(vector<lsquare> partMOLS , vector<vector<int> > OA) ;
	void  OAtoMOLS(vector<vector<int> > OA , vector<lsquare> newMOLS  ) ;
	void  permuteOA(vector<vector<int> > OA , int p[],vector<vector<int> > OANew ) ;
	void  getSmallest(vector<lsquare> partMOLS, vector<lsquare> newMOLS ) ;
	vector<lsquare>  transposeMOLS(vector<lsquare> &partMOLS ) ;
	permutation  CSToPerm(vector<int> &curr) ;
	void  genCRS(vector<int> &curr , int ctr, int rem) ;
	void  getCycleStructure(permutation &p, vector<int>  &nofCycles, map<int, vector<vector<int> > > &cycleLength_cycles_map ) ;
	void  getCycleStructure(permutation &p, vector<int>  &nofCycles) ;
	vector<permutation>  genRelevantPermutations(permutation &p ) ;
	bool  isIdentity(permutation p) ;
//	bool  isSmallestPermSymbols(lsquare ls, int nofPerms) ;
	void  changeOrder(vector<lsquare> &partMOLS, vector<int> &order , vector<lsquare> &newMOLS) ;
	void  changeOrder(vector<lsquare> &partMOLS, int first, int second , vector<lsquare> &newMOLS) ;
	vector<lsquare>  changeOrder(vector<lsquare> &partMOLS, int first, int second ) ;
	void  rollUp(vector<lsquare> &partMOLS,  int pos,  permutation &origP);
	void  standardForm(vector<lsquare> &partMOLS, int rowUp, int pos ,permutation &origP ) ;
	vector<int>  rowMeets(permutation &p1, permutation &p2) ;
	int  compareCS(permutation &p , vector<int> &targetCS) ;
	int  compareCS(vector<int> &pCS , vector<int> &targetCS, bool t) ;
	list<permutation>  getShuffles(permutation &pOrig, permutation &pNow) ;
	bool  noSmallerRCS(permutation &smallestRCS ) ;
	bool  noSmallerRCS(permutation &smallestRCS, vector<lsquare> &pMOLS ) ;
	bool  isSmallest4() ;
	bool  isSmallest4(bool t) ;
	bool isSmallest4(permutation &smallestRCS, vector<lsquare> &pMOLS);
	bool  isSmallestConjugateMOLS(vector<lsquare>& pMOLS) ;
	bool  checkFit( permutation &p) ;
	bool  checkOrthogonal(permutation &P) ;
	bool  checkRCS(permutation &P) ;
	void  printDots(int cursquare, int size, int k) ;
	void  printDots(int num) ;
	void  buildPossibleSquare( vector<vector<int> > &square ) ;
	void  generatePossiblePermsToInsertRec(vector<vector<int> > &square, permutation &p, int row,  list<permutation> &possiblePermList ) ;
	void  generatePossiblePermsToInsert(  list<permutation> &possiblePermList ) ;
	void  findMOLS4( ) ;
	void buildCurrentLS();
	list<vector<int> >  getSmallRelCS( vector<lsquare> &pMOLS, list<permutation> &listP);
	void removeUniversal(permutation &p	);
	void addUniversal(permutation &p	);
	void removeFromLS(permutation p);
	void addToLS();
	void generatePossiblePermsLSRec(permutation &p, int row,  list<permutation> &possiblePermList, vector<bool> &left );
	void  generatePossiblePermsToInsert2(  list<permutation> &possiblePermList ) ;
	void generatePossiblePermsLSRec2(permutation &p, int row,  list<permutation> &possiblePermList, vector<vector<int> > &poss, vector<bool> &left );
	vector<permutation>  generatePossiblePermsToInsert3(int square   ) ;
	//void generatePossiblePermsToInsert3(int square   ) ;
	void generatePossiblePermsLSRec3(permutation &p, int row,  vector<permutation> &possiblePermList, vector<vector<int> > &poss, vector<bool> &left );
	std::vector<std::string> & split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string>   split(const std::string &s, char delim );
	permutation  strtoPerm(string s);
	void updateLS();
	permutation flatten(permutation &pNow);
	bool getSmallRelCS( vector<lsquare>& pMOLS );
	bool rcsOrthog(permutation &p1, permutation &p2);
	static bool permutationComp (permutation p1, permutation    p2);
	void updatePossiblePerms();
	bool  getSmallRelCS( vector<lsquare>& pMOLS, permutation& mapTo );
	double calculate_fraction();
} ;



void printMOLSPerms( vector<lsquare> mols, FILE* outfile, int n, int k  ){
	//int k = mols;
	int m =0;
	//cout << endl;
	//int n = mols[0].front().size();
	stringstream outp[k][n];
	for (m=0; m<k; m++){
			vector<vector<int> >  L(n, vector<int>(n,0));
			int uni_ctr = 0, row_ctr = 0;
			int i,j;
			for ( i =0; i<n; i++){
				for ( j=0; j<n; j++){
					L[i][j]=-1;
				}
			}

			for(unsigned int i =0; i< mols[m].size(); i++){
					for(unsigned int j= 0; j< mols[m][i].size(); j++){
					outp[m][i] << mols[m][i][j];
				}
			}

	}
	int i=0, j=0;
	for ( i =0; i<n; i++){
		string s;
		for ( j =0; j<k; j++){
			s =s+ outp[j][i].str()+ (j<k-1? " ":"") ;
		}
		if (s.size()>k)
			fprintf(outfile, "%s\n", s.c_str());
	}
}

// do about .5 seconds of computing
// (note: I needed to add an arg to this;
// otherwise the MS C++ compiler optimizes away
// all but the first call to it!)
//
static double do_some_computing(int foo) {
    double x = 3.14159*foo;
    int i;
    for (i=0; i<50000; i++) {
        x += 5.12313123;
        x *= 0.5398394834;
    }
    return x;
}
//do_checkpoint(out, sqSymCurrUni);
//int do_checkpoint(MFILE& mf, int k, int n, vector<vector<vector<int> > > sqSymCurrUni) {
int do_checkpoint(int k, int n, vector<vector<vector< int> > > sqSymCurrUni, vector<long long int>& branchCount, clock_t dtime, int countMOLS, vector<vector<vector< int> > > root, int& send_checkpoint, int nof_calls) {

    int retval;
    string resolved_name;
    char output_path[512];
    boinc_resolve_filename(OUTPUT_FILENAME, output_path, sizeof(output_path));
    string outsign = "B";
    //send_checkpoint = true; //!@ 
    FILE* f;	
    if (send_checkpoint==1){
    	cout <<  "sending checkpoint as result -  nofcallsmade"<<endl;
    	f = boinc_fopen(output_path, "ab");

    	send_checkpoint = 3;
    }
    else{
    	if (send_checkpoint==2){
    	    	cout <<  "sending checkpoint as result - limit"<<endl;
    	    	f = boinc_fopen(output_path, "ab");
    	    	send_checkpoint = 3;
    	    	outsign = "A";
    	    }
    	else
    		f = boinc_fopen("temp", "ab");
    }
    	
    if (!f) {cout<< "file open failed"<<endl;return 1;}
    fprintf(f, "@%s\n%d %d\n", outsign.c_str(), n, k);
    printMOLSPerms(root, f,n,k);
    fprintf(f, "#Positions\n");
    for (int j=0; j<n; j++){
    	for (int i=0; i<k; i++ ){
    		fprintf(f, "%d %d %d ", sqSymCurrUni[i][j][0], sqSymCurrUni[i][j][2], sqSymCurrUni[i][j][1]);
    	}
    	fprintf(f, "\n");
    }
    fprintf(f, "%d \n",nof_calls);
    fprintf(f, "#Branchcounts\n");
    for (int i=0; i<branchCount.size(); i++ ){
    	//cout<< branchCount[i] << " ";
        fprintf(f, "%ld ", branchCount[i]);
    }
    //	cout<<endl;
    fprintf(f, "\n");
    fprintf(f, "#Total time\n");
   // "%.3f seconds \n", (clock() - dtime+prev_time)/1000000.0
    fprintf(f, "%.3f seconds\n", dtime/1000000.0 );
    fprintf(f, "#MOLS found\n");
    fprintf(f, "%d\n", countMOLS);

    //fprintf(f, "%d", nchars);
    fclose(f);

    // retval = mf.flush();
    //if (retval) return retval;
    
    boinc_resolve_filename_s(CHECKPOINT_FILE, resolved_name);
    retval = boinc_rename("temp", resolved_name.c_str());
    //cout <<"on return send_valis="<< send_checkpoint<<endl;
    //this->out.flush();
    // retval = boinc_rename("temp", resolved_name);
    if (retval) return retval;

    return 0;
}

#ifdef APP_GRAPHICS
void update_shmem() {
    if (!shmem) return;

    // always do this; otherwise a graphics app will immediately
    // assume we're not alive
    shmem->update_time = dtime();

    // Check whether a graphics app is running,
    // and don't bother updating shmem if so.
    // This doesn't matter here,
    // but may be worth doing if updating shmem is expensive.
    //
    if (shmem->countdown > 0) {
        // the graphics app sets this to 5 every time it renders a frame
        shmem->countdown--;
    } else {
        return;
    }
    shmem->fraction_done = boinc_get_fraction_done();
    shmem->cpu_time = boinc_worker_thread_cpu_time();;
    boinc_get_status(&shmem->status);
}
#endif

int main(int argc, char **argv) {
	start_time=dtime();
    int i;
    int c, nchars = 0, retval, n;
    double fsize, fd;
    char line[64];
    char input_path[512], output_path[512], chkpt_path[512], buf[256];
    MFILE out;
    FILE* state, *infile;

    for (i=0; i<argc; i++) {
        if (strstr(argv[i], "early_exit")) early_exit = true;
        if (strstr(argv[i], "early_crash")) early_crash = true;
        if (strstr(argv[i], "early_sleep")) early_sleep = true;
        if (strstr(argv[i], "run_slow")) run_slow = true;
        if (strstr(argv[i], "critical_section")) critical_section = true;
        if (strstr(argv[i], "cpu_time")) {
            cpu_time = atof(argv[++i]);
        }
        if (strstr(argv[i], "trickle_up")) trickle_up = true;
        if (strstr(argv[i], "trickle_down")) trickle_down = true;
    }
    fprintf(stderr, "%s app started; CPU time %f, flags:%s%s%s%s%s%s%s\n",
        boinc_msg_prefix(buf, sizeof(buf)),
        cpu_time,
        early_exit?" early_exit":"",
        early_crash?" early_crash":"",
        early_sleep?" early_sleep":"",
        run_slow?" run_slow":"",
        critical_section?" critical_section":"",
        trickle_up?" trickle_up":"",
        trickle_down?" trickle_down":""
    );

    retval = boinc_init();
    if (retval) {
        fprintf(stderr, "%s boinc_init returned %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit(retval);
    }
    boinc_set_min_checkpoint_period(160);
    // open the input file (resolve logical name first)
    //
    boinc_resolve_filename(INPUT_FILENAME, input_path, sizeof(input_path));
    infile = boinc_fopen(input_path, "r");
    if (!infile) {
        fprintf(stderr,
            "%s Couldn't find input file, resolved name %s.\n",
            boinc_msg_prefix(buf, sizeof(buf)), input_path
        );
        exit(-1);
    }

    // get size of input file (used to compute fraction done)
    file_size(input_path, fsize);

    boinc_resolve_filename(OUTPUT_FILENAME, output_path, sizeof(output_path));
	cout<<"out"<<endl;
    // See if there's a valid checkpoint file.
    // If so seek input file and truncate output file
    //
    boinc_resolve_filename(CHECKPOINT_FILE, chkpt_path, sizeof(chkpt_path));
    state = boinc_fopen(chkpt_path, "r");
    /*if (state) {
        n = fscanf(state, "%d", &nchars);
        fclose(state);
    }
    if (state && n==1) {
        fseek(infile, nchars, SEEK_SET);
        boinc_truncate(output_path, nchars);
        retval = out.open(output_path, "ab");
        fprintf(stderr,"worker: restarting at %d \n", i_start);
    } else {
    	fprintf(stderr,"worker: Starting from scratch.  \n");
        retval = out.open(output_path, "wb");
    }*/
    if (state) {
    	//n = fscanf(state, "%d", &nchars);
    	//read checkpoint file and do something
    	//fclose(state);
    	retval = out.open(output_path, "ab");
    	out.close();
    	fprintf(stderr,"worker: restarting at %d \n", 1);
    	cout << "restarting from state"<<endl;
    	MOLS threemols(infile, output_path); //change def to make it a checkpoint-start
    	threemols.restart(state, true);

    	//threemols.enumerateMOLS(out);
    }
    else {
    	fprintf(stderr,"worker: Starting from scratch.  \n");

    	retval = out.open(output_path, "wb");
    	out.close(); //destroy current outfile

    	MOLS threemols(infile, output_path);
    	cout<<"Starting"<<endl;
    	threemols.enumerateMOLS( );

    }
    if (retval) {
        fprintf(stderr, "%s APP: upper_case output open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
        fprintf(stderr, "%s resolved name %s, retval %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), output_path, retval
        );
        perror("open");
        exit(1);
    }

#ifdef APP_GRAPHICS
    // create shared mem segment for graphics, and arrange to update it
    //
    shmem = (UC_SHMEM*)boinc_graphics_make_shmem("uppercase", sizeof(UC_SHMEM));
    if (!shmem) {
        fprintf(stderr, "%s failed to create shared mem segment\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
    }
    update_shmem();
    boinc_register_timer_callback(update_shmem);
#endif

    /*for (i=0; ; i++) {
	std::cout<<i;
	fgets(line,64,infile);
	if (feof(infile)) break;
	out.printf(line);
    	out.printf( "lll \n");
    }*/
	//MOLS threemols(infile, out);
    //MOLS threemols(7,3);
	//string filename = argv[1];

	//threemols.enumerateMOLS(out);


   /* for (i=0; ; i++) {
        c = fgetc(infile);

        if (c == EOF) break;
        c = toupper(c);
        out._putchar(c);
	out._putchar('G');
        nchars++;
        if (run_slow) {
            boinc_sleep(1.);
        }

        if (early_exit && i>30) {
            exit(-10);
        }

        if (early_crash && i>30) {
            boinc_crash();
        }
        if (early_sleep && i>30) {
            boinc_disable_timer_thread = true;
            while (1) boinc_sleep(1);
        }

        if (boinc_time_to_checkpoint()) {
            retval = do_checkpoint(out, nchars);
            if (retval) {
                fprintf(stderr, "%s APP: upper_case checkpoint failed %d\n",
                    boinc_msg_prefix(buf, sizeof(buf)), retval
                );
                exit(retval);
            }
            boinc_checkpoint_completed();
        }

        fd = nchars/fsize;
        if (cpu_time) fd /= 2;
        boinc_fraction_done(fd);
    }*/

    /*retval = out.flush();
    if (retval) {
        fprintf(stderr, "%s APP: upper_case flush failed %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit(1);
    }*/

    if (trickle_up) {
        boinc_send_trickle_up(
            const_cast<char*>("example_app GG"),
            const_cast<char*>("sample trickle message GG")
        );
    }

    if (trickle_down) {
        boinc_sleep(10);
        retval = boinc_receive_trickle_down(buf, sizeof(buf));
        if (!retval) {
            fprintf(stderr, "Got trickle-down message: %s\n", buf);
        }
    }

    // burn up some CPU time if needed
    //
    /*if (cpu_time) {
        double start = dtime();
        for (i=0; ; i++) {
            double e = dtime()-start;
            if (e > cpu_time) break;
            fd = .5 + .5*(e/cpu_time);
            boinc_fraction_done(fd);

            if (boinc_time_to_checkpoint()) {
                retval = do_checkpoint(out, nchars);
                if (retval) {
                    fprintf(stderr, "%s APP: upper_case checkpoint failed %d\n",
                        boinc_msg_prefix(buf, sizeof(buf)), retval
                    );
                    exit(1);
                }
                boinc_checkpoint_completed();
            }
            if (critical_section) {
                boinc_begin_critical_section();
            }
            comp_result = do_some_computing(i);
            if (critical_section) {
                boinc_end_critical_section();
            }
        }
    }*/
    boinc_fraction_done(1);
#ifdef APP_GRAPHICS
    update_shmem();
#endif
    boinc_finish(0);
}

#ifdef _WIN32
int WINAPI WinMain(
    HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode
) {
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line(command_line, argv);
    return main(argc, argv);
}
#endif







inline bool MOLS::My::permutationComparator    (permutation p1, permutation    p2) { return p1[0] < p2[0]; }

inline bool MOLS::permutationComp    (permutation p1, permutation    p2) { return p1[0] < p2[0]; }

MOLS::~MOLS(){}

MOLS::MOLS( int n, int k){
	//nof_calls_limit = 500000;
    nof_calls_made=0;
	send_checkpoint = 0;
	printOut = true;
	this->n = n;
	this->k = k;
	currSquare = 0;
 	count_MOLS = 0;
	branchCount_.resize(n*k, 0);
	currentCS.resize(n+1, 0);
	RCSsquares.clear();
	RCSperms.clear();
	z.resize(n);
	diffs.resize(n);
	sizes.resize(n);
	sqSymPossPerms.resize(k, vector<vector<permutation> >(n) );
	identityIt.resize(k, 0 );
	identityIt[0]=1;
	//detailedCount.resize(n*k+1);
	 numIsSmallest.resize(n, 0) ;
	 numIsSmallestTrue.resize(n, 0) ;
	 sqSymCurrUni.resize(k, vector<vector<int> >(n, vector<int>(3, -1)));
	 prev_time = (clock()-clock())/1000000.0;

    for (unsigned int i=0; i<n; i++){
		identity.push_back(i);
	}

	lsquare emptyLS;

	int i,j, ii ;
	for (i =0; i<k; i++){
 		partMOLS.push_back(emptyLS);
		tempMOLS.push_back(emptyLS);
		root.push_back(emptyLS);
	}
	//Get cycle structure representatives for U_0^1
	vector<int> cycles(n+1, 0);
	cycles[1] =1;
	genCRS(cycles, 2, n-1);
	vector<permutation>::iterator CSRit;

	if (printOut){
		cout << "u0^1 list has "<< cycleStructureReps.size()<< " universals - ";
		for (CSRit=cycleStructureReps.begin(); CSRit!= cycleStructureReps.end(); ++CSRit){
			printPerm(*CSRit);
		}
		cout<<endl;
	}

	//Initialize currLS
	currLS.resize(k, vector<vector<int> > (n , vector<int>(n, this->n)));
 	cout<< "constructor completed"<<endl;

	return;
}

std::vector<std::string> & MOLS::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> MOLS::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

permutation MOLS::strtoPerm(string s){
	//vector<string> vec = split(s, '-');
	int i;
	permutation p;
	for(i=0; i<n; i++){
		p.push_back(  s[i]-'0');
	}
	return p;
}
MOLS::MOLS(FILE* infile, MFILE out){
	//nof_calls_limit = 500000;
        nof_calls_made=0;
	printOut = true;
	send_checkpoint = 0;

	char sc[100];
	string s;
	int i;
	//Readall the existing mols untill we get to the start of the checkpoint
	fgets(sc,100,infile);
	s= (string) sc;
	while(sc[0] !='@'){		
		out.printf( sc);
		fgets(sc,100,infile);
	}
	out.flush(); out.close();
	
	for (i=0; ; i++) {
		if (i==0){
			fgets(sc,100,infile);
			s = (string) sc;
			vector<string> vec = split(s,' ');
			n =  vec[0][0]- '0';
			if (n==1) n=10;
			k = vec[1][0] - '0';

			int i;
			for (i =0; i<k; i++){
				lsquare l1, l2, l3;
				partMOLS.push_back(l1);
				tempMOLS.push_back(l2);
				root.push_back(l3);
			}
			sqSymCurrUni.resize(k, vector<vector<int> >(n, vector<int>(3, -1)));
		}
		else{
			fgets(sc,100,infile);
			if (feof(infile)) break;
			s = (string) sc;
			if (s[0]=='#') break;
			cout<< s;
			vector<string> vec = split(s,' ');
			int numUni = vec.size();

			for (int j=0; j< numUni; j++){
				partMOLS[j].push_back(strtoPerm(vec[j]));
				root[j].push_back(strtoPerm(vec[j]));
			}
			currSquare = numUni%k;
		}
	}

	prev_time = 0;

	/*int i;
	string s;
	ifstream infile;
	infile.open(filename.c_str(), ios::in );

	 if (infile.is_open())
	  {
		 if  ( infile.good() ){
			 getline (infile,s);
			 vector<string> vec = split(s,' ');
			 n =  vec[0][0]- '0';
			 if (n==1) n=10;
			 k = vec[1][0] - '0';
			 int i,j, ii ;
			 	for (i =0; i<k; i++){
			 		lsquare l1, l2;
			 		partMOLS.push_back(l1);
			 		tempMOLS.push_back(l2);
			 	}
		 }

		cout<< filename <<endl;
	    while ( infile.good() )
	    {
	      getline (infile,s);
	       if (s.size()>0){
	    	   vector<string> vec = split(s,' ');
	    	   int numUni = vec.size();

	    	  for (i=0; i< numUni; i++){
	    		  partMOLS[i].push_back(strtoPerm(vec[i]));
	    	  }
	    	  currSquare = numUni%k;
	      }

	    }
	    infile.close();
	  }
	 else
		 cout<< "File not found"<<endl;*/

	printMOLSPerms(partMOLS);
	//cout<< n<< k;
	//possibleShuffles = genRelevantPermutations(partMOLS[1][0]);
	count_MOLS = 0;
	branchCount_.resize(n*k, 0);
	currentCS.resize(n+1, 0);
	RCSsquares.clear();
	RCSperms.clear();
	z.resize(n);
	diffs.resize(n);
	sizes.resize(n);
	identityIt.resize(k, 0 );
	identityIt[0]=1;
	for (unsigned int i=0; i<n; i++){
			identity.push_back(i);
		}

	sqSymPossPerms.resize(k, vector<vector<permutation> >(n) );

	//detailedCount.resize(n*k+1);
	numIsSmallest.resize(n, 0) ;
	numIsSmallestTrue.resize(n, 0) ;

	//Get cycle structure representatives for U_0^1
	vector<int> cycles(n+1, 0);
	cycles[1] =1;
	genCRS(cycles, 2, n-1);
	vector<permutation>::iterator CSRit;
	if (printOut){
		cout << "u0^1 list has "<< cycleStructureReps.size()<< " universals - ";
		for (CSRit=cycleStructureReps.begin(); CSRit!= cycleStructureReps.end(); ++CSRit){
			printPerm(*CSRit);
		}
		cout<<endl;
	}

	//Initialize currLS
 	updateLS();
  	updatePossiblePerms();
  	printAllStatics();
  	rewind(infile);//move thepointerbackout of the positions
  	restart(infile,false);
 	//Initialise symbSqRC
//	symbSqRC.resize(n, vector<vector<vector<int> > >(k, vector<vector<int> >(n, vector<int>(n, 1))));
	//completedMOLS.clear();
 	//cout<<"done"<<endl;
	return;
}

MOLS::MOLS(FILE* infile, char* output_p){
	//nof_calls_limit = 500000;
        nof_calls_made=0;
	//this->filename = filename;
	printOut = true;
	send_checkpoint = 0;
	this->output_path = output_p;

	char sc[100];
	string s;
	int i;	
	MFILE out;
	//out.open(output_path, "wb");
	cout<<"a"<<endl;
	//skip all the first lines, copyall the found mols over to outfile
	fgets(sc,200,infile);
	//s= (string) sc;
	while(sc[0]!='@'){
		cout<<"XX"<<sc<<endl;
		//out.printf( sc);
		fgets(sc,200,infile);		
	}
	//cout<<"done while"<<endl;
	//out.flush(); out.close();
	//cout<<"a"<<endl;
	for (i=0; ; i++) {
		if (i==0){
			fgets(sc,200,infile);
			s = (string) sc;
			cout << "SR"<<s<<endl;
			vector<string> vec = split(s,' ');
			n =  vec[0][0]- '0';
			if (n==1) n=10;
			k = vec[1][0] - '0';
			cout<<"SR"<<n<<k<<endl;
			for (unsigned int j =0; j<k; j++){
				lsquare l1, l2, l3;
				partMOLS.push_back(l1);
				tempMOLS.push_back(l2);
				root.push_back(l3);
			}
			sqSymCurrUni.resize(k, vector<vector<int> >(n, vector<int>(3, -1)));
			cout<<"Done reading"<<endl;
		}
		else{
			fgets(sc,200,infile);
			if (feof(infile)) break;
			s = (string) sc;
			if (sc[0]=='#') break;
			cout<< "Comp " <<s<<endl;
			vector<string> vec = split(s,' ');
			int numUni = vec.size();

			for (int j=0; j< numUni; j++){
				partMOLS[j].push_back(strtoPerm(vec[j]));
				root[j].push_back(strtoPerm(vec[j]));
			}
			currSquare = numUni%k;
		}
	}
	cout<<"a"<<endl;
	prev_time = 0;

	printMOLSPerms(partMOLS);
	//cout<< n<< k;
	//possibleShuffles = genRelevantPermutations(partMOLS[1][0]);
	count_MOLS = 0;
	branchCount_.resize(n*k, 0);
	currentCS.resize(n+1, 0);
	RCSsquares.clear();
	RCSperms.clear();
	z.resize(n);
	diffs.resize(n);
	sizes.resize(n);
	identityIt.resize(k, 0 );
	identityIt[0]=1;
	for (unsigned int i=0; i<n; i++){
			identity.push_back(i);
		}

	sqSymPossPerms.resize(k, vector<vector<permutation> >(n) );

	//detailedCount.resize(n*k+1);
	numIsSmallest.resize(n, 0) ;
	numIsSmallestTrue.resize(n, 0) ;

	//Get cycle structure representatives for U_0^1
	vector<int> cycles(n+1, 0);
	cycles[1] =1;
	genCRS(cycles, 2, n-1);
	vector<permutation>::iterator CSRit;
	if (printOut){
		cout << "u0^1 list has "<< cycleStructureReps.size()<< " universals - ";
		for (CSRit=cycleStructureReps.begin(); CSRit!= cycleStructureReps.end(); ++CSRit){
			printPerm(*CSRit);
		}
		cout<<endl;
	}

	//Initialize currLS
 	updateLS();
  	updatePossiblePerms();
  	printAllStatics();
  	rewind(infile);//move thepointerbackout of the positions
  	restart(infile,false);
 	//Initialise symbSqRC
//	symbSqRC.resize(n, vector<vector<vector<int> > >(k, vector<vector<int> >(n, vector<int>(n, 1))));
	//completedMOLS.clear();
 	//cout<<"done"<<endl;
	return;
}

void MOLS::updateLS(){
 	currLS.resize(k, vector<vector<int> > (n , vector<int>(n, n)));
 	int   i;
	for (i=0; i<k; i++){
 		//		for(universalIt=partMOLS[i].begin(); universalIt != partMOLS[i].end(); ++universalIt){
		for(unsigned int j =0; j< partMOLS[i].size(); j++){
  			for(unsigned int jj= 0; jj< partMOLS[i][j].size(); jj++){//(partMOLS[i][j]).begin(); perm_it != (partMOLS[i][j]).end(); ++perm_it){
				//cout<<partMOLS[i][j][jj]<<"-";				cout.flush();
				currLS[i][jj][partMOLS[i][j][jj]] = j;
			}
		}
	}

}

int MOLS::enumerateMOLS(){
/*	string outfile = filename+ ".out.txt";
	std::ofstream out(outfile);
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf());*/
	//this->out = out;
	cout<<"starting enumeration "<< prev_time<< endl;
	dtime = clock();
	cout<<   k << " MOLS of order "<< n<<endl;
	int i;
	//count_MOLS =0;

	if (partMOLS[0].size()==0){
		permutation ss;
		for (i=0;i<n; i++){
			ss.push_back(i);
		}
		addUniversal((ss));
		currSquare++;
		findMOLS4();
	}
	else{
		if (partMOLS[1].size()>0){
			map<int, vector<vector<int> > > dummy_cycles_map;
			getCycleStructure(partMOLS[1].front(),  currentCS, dummy_cycles_map );
			possibleShuffles = genRelevantPermutations(partMOLS[1].front());
		}

 		//accout for funny behaviour on first branch
		int dec = -1;
		for (i=0; i<k; i++)
			dec = dec+root[i].size();
		branchCount_[dec] = -1;
		/////

		cout<<"about to call findMOLS - "<< dec<< " send_scheckpoint:"<< send_checkpoint<<", "<<nof_calls_made<<endl;
		findMOLS4();
	}
	cout<<"Backfrom findmols"<<endl;
	if (send_checkpoint<1){ //ifwe are sending the checkpoint, then the out file is created by do checkpoint
	MFILE out;
	out.open(output_path, "ab");
	out.printf("@A\n");
	out.printf("#Branchcounts\n");
	for (i=0; i<n*k;i++){
		cout<<branchCount_[i]<<" ";
		out.printf( "%d ", branchCount_[i] );
	} out.printf("\n");
	out.printf("#Total time\n");
	out.printf( "%.3f seconds \n", (clock() - dtime+prev_time)/1000000.0  ); 
	cout<< (clock() - dtime+prev_time)/1000000.0<<endl;
	out.printf("#MOLS found\n");
	out.printf("%d" , count_MOLS)	;
	out.flush(); out.close();
	}
	/*	
	out.printf( "# %d %d MOLS of order %d found \n",count_MOLS, k,n );//this->out.flush();
	cout<< "# "<<count_MOLS << " MOLS found"<<endl;
    cout<< "# ";
    out.printf( "# " );
	for (i=0; i<n*k;i++){
		cout<<branchCount_[i]<<" ";
		out.printf( "%d ", branchCount_[i] );
	}
	cout<<endl;out.printf( "\n" );
	for (i=0; i< n; i++){
		cout<<"# "<< branchCount_[k-1 +i*k] <<" takke op vlak "<<i<< " | "<< numIsSmallest[i]<<" "<<numIsSmallestTrue[i]<< " "<< (numIsSmallestTrue[i]*1.)/numIsSmallest[i]<<endl;
		out.printf( "# %d takke op vlak %d | %d %d %f \n", branchCount_[k-1 +i*k], i,numIsSmallest[i],numIsSmallestTrue[i], (numIsSmallestTrue[i]*1.)/numIsSmallest[i]  );
	}
	//dtime = clock() - dtime;
	cout<<"# "<<  (clock() - dtime+prev_time)/1000000.0 << " seconds"<< endl;
	out.printf( "# %.3f seconds \n", (clock() - dtime+prev_time)/1000000.0  ); out.flush();
	out.flush();out.close();
	cout<< count_MOLS;
	*/
	return count_MOLS;

}

int MOLS::printUniversals(lsquare& l){
	list<list<int> >::iterator universalIt; //iterates over universals
	list<int>::iterator perm_it; //iterates over permutation
	for (unsigned int i=0; i<l.size(); i++){
	//for(universalIt=l.begin(); universalIt != l.end(); ++universalIt){

		//list<int>& universalIt = *i;
		cout<<"< ";
		//for(perm_it=(*universalIt).begin(); perm_it != (*universalIt).end(); ++perm_it){
		for (unsigned int j=0; j<l[i].size(); j++){
			cout<<l[i][j];
		}
		cout<<" >";
	}
	 return 0;
}

void MOLS::printPerm(permutation p){
	list<int>::iterator perm_it; //iterates over permutation
	cout<<" ";
	for (unsigned int j=0; j<p.size(); j++){
		cout<<p[j];
	}
	cout<<" ";
}
/*string MOLS::printPerm(permutation p , boolean t){
	string s ="";
	for (unsigned int j=0; j<p.size(); j++){
		s=s+p[j];
	}
	return s;
}*/

void MOLS::printLS(lsquare & l){
	//int n = l.front().size();
	vector<vector<int> >  L(n, vector<int>(n,0));
	int uni_ctr = 0, row_ctr = 0;
	int i,j;
	for ( i =0; i<n; i++){
		for ( j=0; j<n; j++){
			L[i][j]=-1;
		}
	}

	list<list<int> >::iterator universalIt; //iterates over universals
	list<int>::iterator perm_it; //iterates over permutation
	uni_ctr=0;
	for(unsigned int j =0; j< l.size(); j++){
		for(unsigned int jj= 0; jj< l[j].size(); jj++){
//	for(universalIt=l.begin(); universalIt != l.end(); ++universalIt){
		//for(perm_it=(*universalIt).begin(); perm_it != (*universalIt).end(); ++perm_it){
			L[jj][l[j][jj]] = j;
			//row_ctr++;
		}
		//row_ctr=0;
		//uni_ctr++;
	}
	uni_ctr = 0;
	//int i,j;
	for ( i =0; i<n; i++){
		for ( j=0; j<n; j++){
			if (L[i][j]>= 0)
				//cout <<  L[i][j]   <<" "<< (j==n-1? "\\\\" : "& ");
				cout <<  L[i][j]   <<" ";
			else
				cout <<    "- ";
				//cout <<    "- "<< (j==n-1? "\\\\" : "& ");
		}
		cout<< endl;
	}
	cout<< endl;

}

void MOLS::printMOLSPerms( vector<lsquare> mols, FILE* outfile  ){
	//int k = mols;
	int m =0;
	//cout << endl;
	//int n = mols[0].front().size();
	stringstream outp[k][n];
	for (m=0; m<k; m++){
			vector<vector<int> >  L(n, vector<int>(n,0));
			int uni_ctr = 0, row_ctr = 0;
			int i,j;
			for ( i =0; i<n; i++){
				for ( j=0; j<n; j++){
					L[i][j]=-1;
				}
			}

			for(unsigned int i =0; i< mols[m].size(); i++){
					for(unsigned int j= 0; j< mols[m][i].size(); j++){
					outp[m][i] << mols[m][i][j];
				}
			}

	}
	int i=0, j=0;
	for ( i =0; i<n; i++){
		string s;
		for ( j =0; j<k; j++){
			s =s+ outp[j][i].str()+ " ";
		}
		if (s.size()>k)
			fprintf(outfile, "%s\n", s.c_str());
	}
}

void MOLS::printMOLSPerms( vector<lsquare> mols  ){
	//int k = mols;
	int m =0;
	//cout << endl;
	//int n = mols[0].front().size();
	stringstream outp[k][n];
	for (m=0; m<k; m++){
		//if (mols[m].size()>0){

			vector<vector<int> >  L(n, vector<int>(n,0));
			int uni_ctr = 0, row_ctr = 0;
			int i,j;
			for ( i =0; i<n; i++){
				for ( j=0; j<n; j++){
					L[i][j]=-1;
				}
			}

			//list<list<int> >::iterator universalIt; //iterates over universals
			//list<int>::iterator perm_it; //iterates over permutation
			//uni_ctr=0;
			for(unsigned int i =0; i< mols[m].size(); i++){
					for(unsigned int j= 0; j< mols[m][i].size(); j++){
			//for(universalIt=mols[m].begin(); universalIt != mols[m].end(); ++universalIt){
			//	for(perm_it=(*universalIt).begin(); perm_it != (*universalIt).end(); ++perm_it){
					//L[row_ctr][*perm_it] = uni_ctr;
					//row_ctr++;

					outp[m][i] << mols[m][i][j];
				}
				//row_ctr=0;
				//uni_ctr++;
			}

		//}
		//else cout<< "Empty"<< endl;
	}
	int i=0, j=0;
	for ( i =0; i<n; i++){
		string s;
		for ( j =0; j<k; j++){
			//s =s+ outp[j][i].str()+(outp[j][i].str().size()==0? new string(k+1, " "): " ");
			s =s+ outp[j][i].str()+ " ";
		}
		if (s.size()>k)
			cout<<s<<endl;
	}
}

void MOLS::printMOLS( vector<lsquare> mols  ){
	int m;
	stringstream outp[k][n];
	for (m=0; m<k; m++){
		//if (mols[m].size()>0){

		int i,j;
		for ( i =0; i<n; i++){
			for ( j=0; j<n; j++){
				if (currLS[m][i][j]>=  0 )
					//cout <<  L[i][j]   <<" "<< (j==n-1? "\\\\" : "& ");
					outp[m][i]<<currLS[m][i][j]<<" ";
				else
					outp[m][i]<<  "- ";
				//cout <<    "- "<< (j==n-1? "\\\\" : "& ");
			}

		}
		//}
		//else cout<< "Empty"<< endl;
	}

	MFILE out;
	out.open(output_path, "ab");
	int i=0, j=0;
	for ( i =0; i<n; i++){
		for ( j =0; j<k; j++){
			cout<<outp[j][i].str()<<"\t";
		    out.printf("%s \t",outp[j][i].str().c_str());
		}
		cout<<endl;
		out.printf("\n");
	}
	cout <<"*******************************"<<endl;
	out.printf("*****************************\n");
	out.flush(); out.close();

}

/*is Orthogonal accepts a relative cycle structure in the form of
 * a permutation and checks wheter the two original permutations
 * are orthonal, ie, whether the rcs has exactly one fixed point.
 * Input: perm - the permutation of the relative cycle structure to be checked
 * Output: true if orthogonal (only one fixed point)
 */
 bool MOLS::isOrthogonal(permutation perm){

	// cout<< "Enter"<<endl;
	int  fixed_ctr = 0;

	for(unsigned int i =0; i< perm.size(); i++){
		if (i == perm[i]) fixed_ctr++;

	}
	 //cout<< "exit"<<endl;
	if (fixed_ctr==1)
		return true;
	else
		return false;
}

permutation MOLS::composition(permutation &p1, permutation &p2){
	// cout<< "Enter composite"<<endl;

	permutation composite(p1.size()) ; //= {0,0,0,0};
	for(unsigned int i =0; i< n; i++)
        composite[i] = p1[p2[i]];

	return composite;
}

permutation MOLS::inverse(permutation &p){
	// cout<< "Enter inverse"<<endl;
 	permutation inverse(p.size()) ;

	for (unsigned int i=0; i<p.size(); i++){
		inverse[p[i]] = i;
	}
	return inverse;
}
permutation MOLS::copyPerm(permutation p){
	permutation pNew(p) ;
	return pNew;
}

permutation MOLS::rcs(permutation &p1, permutation &p2){
	// cout<< "Enter RCS"<<endl;
//TODO test whether permutations are of same length
	// actually just write the two permutations above each other and sort. has to be faster. although..nlogn
	permutation pNew(p1.size()) ;
	for (unsigned int i=0; i<p1.size(); i++)
		pNew[p1[i]] = p2[i];
	return pNew;

	//return composition(p2, inverse(p1));
}

//void MOLS::permuteSymbols(vector<lsquare> partMOLS, vector<lsquare> &newMOLS){
void MOLS::permuteSymbols(vector<lsquare> &partMOLS, vector<lsquare> &newMOLS){
	int i;
	for (i =0;i<k;i++){
		newMOLS[i].clear();
		if (partMOLS[i].size()>0){ //if the square is still empty it cant be permuted
			lsquare::iterator permIt ;
			for (permIt = partMOLS[i].begin(); permIt != partMOLS[i].end(); ++permIt){
				newMOLS[i].push_back(*permIt);
 				//newMOLS[i].push_back(copyPerm(*permIt));
			} //TODO
		/*	for (unsigned int j=0; j< partMOLS[i].size(); j++)
				newMOLS[i][partMOLS[i][j][0]] = partMOLS[i][j];*/
			//newMOLS[i].std::sort(newMOLS[i].begin(), newMOLS[i].end(), My::permutationComparator);
			std::sort(newMOLS[i].begin(), newMOLS[i].end(), My::permutationComparator);
		}
	}

}

//true if partMOLS (1) is smaller or equal (weak ordering)to newMOLS(pMOLS after permutation (2))
 // 1 < 2 - TRUE
// 1 = 2 - TRUE
// 1 > 2 - FALSE
bool MOLS::testPermuteMOLS(vector<lsquare> &pMOLS, permutation &rowPerm, permutation &colPerm ){
	vector<int> permIt1(k,0);//arrays of iterators
	permIt1[0]++;

	permutation invRowPerm = inverse(rowPerm);

	int indOfOne =-1; bool found = false;
	unsigned int  i=0;
	unsigned int cs = 0;  //the number of the square currently being compared.

	for (cs =0; cs<k; cs++){
		sizes[cs] = pMOLS[cs].size();
	}

	indOfOne = invRowPerm[0];

	unsigned int  ci=0;
	int diff;
	for (cs = 1%k; true; cs=(cs+1)%k){//this goes on forever, break on when iterato reaches end

		if (permIt1[cs] == partMOLS[cs].size() || permIt1[cs] == sizes[cs])
			return true;
		ci=0;
		while(ci<sizes[cs]){
			if (colPerm[pMOLS[cs][ci][indOfOne]]==permIt1[cs]){
				//cout<< "ci = "<< ci;
				break;}
			ci++;
		}

		if (ci< sizes[cs]){ //so weve found the one that gets a whatever in the first position
			//permutation z(n, 0);
			//z.resize(n);
			//diffs.resize(n);

//			for (unsigned int  jj=0; jj<n; jj++){
//				diffs[invRowPerm[rowPerm[jj]]] = partMOLS[cs][permIt1[cs]][ invRowPerm[rowPerm[jj]] ] - colPerm[pMOLS[cs][ci][rowPerm[jj]]];
//				if (diffs[invRowPerm[rowPerm[jj]]] !=0){
//                    //cout<< partMOLS[cs][permIt1[cs]][ rowPerm[jj]]<<"|"<<colPerm[pMOLS[cs][ci][jj]];
//                    //cout<<"R"<<(diffs[invRowPerm[jj]]<0)<< " ";
//                    return (diffs[invRowPerm[rowPerm[jj]]]<0);
//                }
//			}//cout<<"X"<< " ";

			for (unsigned int  jj=0; jj<n; jj++){

				//diffs[rowPerm[jj]] = partMOLS[cs][permIt1[cs]][ rowPerm[jj]] - colPerm[pMOLS[cs][ci][jj]];
				//diffs[rowPerm[invRowPerm[jj]]] = partMOLS[cs][permIt1[cs]][ rowPerm[invRowPerm[jj]]] - colPerm[pMOLS[cs][ci][invRowPerm[jj]]];

                if ((partMOLS[cs][permIt1[cs]][ rowPerm[invRowPerm[jj]]] - colPerm[pMOLS[cs][ci][invRowPerm[jj]]]) !=0){
					//cout<< partMOLS[cs][permIt1[cs]][ rowPerm[jj]]<<" | "<<colPerm[pMOLS[cs][ci][jj]]<<"Return "<<(diff<0);
					//cout<<"Return "<<(diffs[jj]<0);
					//return (diffs[rowPerm[invRowPerm[jj]]]<0);
					return ((partMOLS[cs][permIt1[cs]][ rowPerm[invRowPerm[jj]]] - colPerm[pMOLS[cs][ci][invRowPerm[jj]]])<0);
				}
				//z[ rowPerm[jj] ] = colPerm[pMOLS[cs][ci][jj]];
				//printPerm(z);


			}
			//printPerm(diffs);
//			for (unsigned int  jj=0; jj<n; jj++){
//				if (diffs[jj] !=0){
//					//cout<< partMOLS[cs][permIt1[cs]][ rowPerm[jj]]<<" | "<<colPerm[pMOLS[cs][ci][jj]]<<"Return "<<(diff<0);
//					//cout<<"Return "<<(diffs[jj]<0);
//					return (diffs[jj]<0);
//				}
//			}
			//cout<< "Equal "<<cs<<ci<<endl;
		}
		else return true;

		permIt1[cs]++;

	}

	return true;
}

//permutes i
void MOLS::permuteMOLS(vector<lsquare> &pMOLS, permutation &rowPerm, permutation &colPerm ){
	//permutation z(n);
	//vector<lsquare> newMOLS(k);

	for (unsigned int i=0; i<k; i++	){
 		if (pMOLS[i].size()>0){ //if the square is not empty (otherwise it cant be permuted)
 			unsigned int sz = pMOLS[i].size();
			for (unsigned int j=0; j<sz; j++){
				for (unsigned int  jj=0; jj<n; jj++)
					z[ rowPerm[jj] ] = colPerm[pMOLS[i][j][jj]];
				//newMOLS[i].push_back(z);
				pMOLS[i][j] = z;
 			}

 		}
		std::sort(pMOLS[i].begin(), pMOLS[i].end(), My::permutationComparator );
	}
	//pMOLS= newMOLS;


//	std::list<int> third (second.begin(),second.end());  // iterating through second
}




void MOLS::printOA(vector<vector<int> > OA){
	int i; int j; int w = n*n;
	cout << "ORTHOGONAL ARRAY"<< endl;
		for (i=0; i<k+2; i++){
			for (j=0; j<w; j++){
					cout<<	OA[i][j] << " ";
			}
			cout<< endl;
		}
}

//true if first is smaller or equal (weak ordering)
// this function has been tested and works with 1-mols and 2mols
//at the very least.  2 mols means by extension also with kmols
// 1 < 2 - TRUE
// 1 = 2 - TRUE
// 1 > 2 - FALSE
bool MOLS::comparePartialMOLS(vector<lsquare> &partMOLS1, vector<lsquare> &partMOLS2){
	vector<int> permIt1(k,0);//arrays of iterators
 	vector<int> doneArr(k, 0);

 	//printMOLSPerms(partMOLS1);
 	//printMOLSPerms(partMOLS2);
 	//cout<<"enter"<<endl;
	permIt1[0]++;
	permIt1[1]++;

	//doneArr keeps track of when the end of each array is reached to prevent the iterator from enumerating an array which is already on its end.

	int currSquare = 0;  //the number of the square currently being compared.
	for (currSquare =0; currSquare<k; currSquare++){
		//if the square is empty the iterator obv wont function properly
		//Assume for now that if one doesnt have anythiong in a square neither does the other, as they are permutations
		if (permIt1[currSquare] >= partMOLS1[currSquare].size() || permIt1[currSquare]>= partMOLS2[currSquare].size()  ){
			doneArr[currSquare] = 1;

		}
	}
	//printPerm(doneArr);printPerm(permIt1); cout.flush();
	int diff;
	unsigned int i;
	/*for (i=0;i<2; i++){

		permIt1[i]++;
		permIt2[i]++;

		//if (permIt1[currSquare]== partMOLS1[currSquare].end() ){
		if (   partMOLS1[i].size()==1 ||  partMOLS2[i].size()==1 ){
			// so if we can move along on this square lets move along
			// remember k cycles between the squares
			doneArr[i] = 1;

			int prod = 1;

			for (i=0; i<k; i++){
				prod=prod*doneArr[i];
			}
			if (prod==1)//then all permutation in all the squares have been enumerated, still equal
				return true;

		}
	}
*/
	//for (currSquare = (2%k); currSquare<k; currSquare=(currSquare+1)%k){//this goes on forever, break on when iterato reaches end
	for (currSquare = 2%k; doneArr[currSquare]==0; currSquare=(currSquare+1)%k){//this goes on forever, break on when iterato reaches end
		for (i=0; i<n; i++	){
		    diff = partMOLS1[currSquare][permIt1[currSquare]][i] - partMOLS2[currSquare][permIt1[currSquare]][i];
			if (diff !=0)
				return (diff<0);
		}
		//if we havent returned anything by now the elements compared are equivalent, move on to next.
		permIt1[currSquare]++;

		//if (permIt1[currSquare]== partMOLS1[currSquare].end() ){
		if ( (permIt1[currSquare]== partMOLS1[currSquare].size() || permIt1[currSquare]== partMOLS2[currSquare].size())
			|| doneArr[(currSquare-1+k)%k]==1){
			// so if we can move along on this square lets move along
			// remember k cycles between the squares
			doneArr[currSquare] = 1;

		}
	}
	return true;
 }

//true if first is smaller or equal (weak ordering)
// this function has been tested and works with 1-mols and 2mols
//at the very least.  2 mols means by extension also with kmols
// 1 < 2 - TRUE
// 1 = 2 - TRUE
// 1 > 2 - FALSE
/*bool MOLS::comparePartialMOLS(vector<lsquare> partMOLS1, vector<lsquare> partMOLS2){
	lsquare::iterator permIt1[k];//arrays of iterators
	lsquare::iterator permIt2[k];
	vector<int> doneArr(k, 0);
	int i=0;

	//doneArr keeps track of when the end of each array is reached to prevent the iterator from enumerating an array which is already on its end.



	int currSquare = 0;  //the number of the square currently being compared.
	for (currSquare =0; currSquare<k; currSquare++){
		//if the square is empty the iterator obv wont function properly
		//Assume for now that if one doesnt have anythiong in a square neither does the other, as they are permutations
		if (partMOLS1[currSquare].size()>0){
			permIt1[currSquare] = partMOLS1[currSquare].begin();
			permIt2[currSquare] = partMOLS2[currSquare].begin();
		}
		else {
			doneArr[currSquare] = 1;
		}

	}

	for (i=0;i<2; i++){

		permIt1[i]++;
		permIt2[i]++;

		//if (permIt1[currSquare]== partMOLS1[currSquare].end() ){
		if (   partMOLS1[i].size()==1 ||  partMOLS2[i].size()==1 ){
			// so if we can move along on this square lets move along
			// remember k cycles between the squares
			doneArr[i] = 1;

			int prod = 1;

			for (i=0; i<k; i++){
				prod=prod*doneArr[i];
			}
			if (prod==1)//then all permutation in all the squares have been enumerated, still equal
				return true;

		}
	}

	//for (currSquare = (2%k); currSquare<k; currSquare=(currSquare+1)%k){//this goes on forever, break on when iterato reaches end
	for (currSquare = 0; currSquare<k; currSquare=(currSquare+1)%k){//this goes on forever, break on when iterato reaches end
		 //printPerm(*permIt1[currSquare]);
		 //printPerm(*permIt2[currSquare]);
		//cout<< lexicographical_compare((*permIt1[currSquare]).begin(), (*permIt1[currSquare]).end(), (*permIt2[currSquare]).begin(), (*permIt2[currSquare]).end() )<<endl;

		if ( lexicographical_compare((*permIt1[currSquare]).begin(), (*permIt1[currSquare]).end(), (*permIt2[currSquare]).begin(), (*permIt2[currSquare]).end() ))
			return true;
		else if ( lexicographical_compare((*permIt2[currSquare]).begin(), (*permIt2[currSquare]).end(), (*permIt1[currSquare]).begin(), (*permIt1[currSquare]).end() ))
			return false;

		//if we havent returned anything by now the elements compared are equivalent, move on to next.

		permIt1[currSquare]++;
		permIt2[currSquare]++;

		//if (permIt1[currSquare]== partMOLS1[currSquare].end() ){
		if ( (permIt1[currSquare]== partMOLS1[currSquare].end() || permIt2[currSquare]== partMOLS2[currSquare].end())
			|| doneArr[(currSquare-1+k)%k]==1){
			// so if we can move along on this square lets move along
			// remember k cycles between the squares
			doneArr[currSquare] = 1;
			int prod = 1;
			for (i=0; i<k; i++){
				prod=prod*doneArr[i];
			}
			if (prod==1)//then all permutation in all the squares have been enumerated, still equal
				return true;
		}
	}

	return true;
	//lexicographical_compare(partMOLS[0].front(), partMOLS[0].back(), ))
}*/

void MOLS::MOLStoOA(vector<lsquare> partMOLS, vector<vector<int> > OA){
	//int n = partMOLS[0].front().size();
	int i=0;
	int j=0;
	int w = n*n;
	for (i=0; i<w; i++){
		OA[0][i] = i / n;
		OA[1][i] = (j++) %n;
	}

	int uni_ctr=0, row_ctr=0;

	/*for (i=0; i<k; i++){
			for    (j=0; j<w; j++){
				OA[i+2][j] = k-1;
			}
		}*/


	for (i=0; i<k; i++){
		uni_ctr = 0;
		lsquare::iterator universalIt; //iterates over universals
		permutation::iterator perm_it; //iterates over permutation

		for(universalIt=partMOLS[i].begin(); universalIt != partMOLS[i].end(); ++universalIt){
			row_ctr=0;
			for(perm_it=(*universalIt).begin(); perm_it != (*universalIt).end(); ++perm_it){
				OA[i+2][n*row_ctr + (*perm_it)] = uni_ctr;
				row_ctr++;
			}
			uni_ctr++;
		}
	}
	uni_ctr = 0; row_ctr = 0;
	//printOA<k,n>(OA);

return;

}

void MOLS::OAtoMOLS(vector<vector<int> > OA, vector<lsquare> newMOLS ){

	map<int, int> pos_to_colidx;
	int i=0;
	int w=n*n;
	permutation parr[k][n];
	lsquare tempMOLS[k];

	for (i=0;i<k; i++){
				newMOLS[i].clear();
		}

	for (i=0; i<w; i++){
		pos_to_colidx[OA[0][i]*n + OA[1][i]] = i ;
		//map.insert(std::make_pair(OA[0][i]*n + OA[1][i], i));
	}

	int j;
	for (i=0;i<w; i++){
		for (j=0; j<k; j++){
			int col = pos_to_colidx[i];
			//cout << col<< " , "<< OA[j+2][col] << " , "<< OA[1][col]<< " " ;cout.flush();
			parr[j][ OA[j+2][col] ].push_back( OA[1][col] );
		}

	}

	for (i=0;i<k; i++){
		for (j=0; j<n; j++){
			newMOLS[i].push_back(parr[i][j]);
		}
	}

	//newMOLS=tempMOLS;
	return;
}

void MOLS::permuteOA(vector<vector<int> > OA, int p[],vector<vector<int> > OANew){
	int i=0;
	int j=0;
	int w = n*n;
	for (i=0; i<k+2; i++){
		for (j=0; j<w; j++){
			OANew[p[i]][j]	= OA[i][j];
		}
	}
	return;
}

void MOLS::getSmallest(vector<lsquare> partMOLS, vector<lsquare> newMOLS){


	permutation colPerm;
	permutation rowPerm;

	//lsquare tempMOLS[2];
int x;
	int i=0;
	for (i=0; i<n; i++){
		rowPerm.push_back(i);
		colPerm.push_back(i);
	}

	while(next_permutation(rowPerm.begin(), rowPerm.end())){
		while(next_permutation(colPerm.begin(), colPerm.end())){
			lsquare tempMOLS[k];

			//lsquare l1; lsquare l2; tempMOLS[0] = l1; tempMOLS[1] = l2;
			//printPerm(rowPerm) ;printPerm(colPerm);

			//transposeMOLS<2>(partMOLS, transpMOLS);
			permuteMOLS(partMOLS, rowPerm, colPerm);
			//tempmols has the postpermutation, first make this even smaller by permutin symbols



			if (!comparePartialMOLS(partMOLS, newMOLS) ){// if  part > final
				return;
			}
			else{
				for (x=0; x<k; x++)
					newMOLS[x].clear();
			}


		}
	}

	for (x=0; x<k; x++)
		newMOLS[x]=partMOLS[x];

	return ;
}



vector<lsquare>  MOLS::transposeMOLS(vector<lsquare> &pMOLS){

	vector<lsquare> newMOLS(k );

	for (unsigned int i =0;i<k;i++){

		unsigned int sz = pMOLS[i].size();
		for (unsigned int j =0; j<sz; j++)
			newMOLS[i].push_back(inverse(pMOLS[i][j]));

	}
	return newMOLS;

	//return newMOLS;
	//TODO shuffle squares
}


permutation MOLS::CSToPerm(vector<int> &curr ){
	int i,j, currCycleLen, ind=0;
	int ctr =0;
	int cycleStart;
	permutation p(n) ;

	for (i=0; i<n+1; i++){
		int j=curr[i]; //j = number of cycles of length i

		while (j-- > 0){ //while there is another cycle remaining of lenght i
			currCycleLen = i;
			cycleStart = ctr;

			while (currCycleLen>0){ //while we are still in this cycle
				//p.push_back((currCycleLen==1? cycleStart: ctr+1	 ));
				 p[ind++] = ((currCycleLen==1? cycleStart: ctr+1	 ));
				currCycleLen--;
				ctr++;
			}
		}
	}

	return p;
}

void MOLS::genCRS(vector<int> &curr, int ctr, int rem ){
	if (ctr==0){
		genCRS(curr, ctr+1, rem );
		return;
	}
	if (rem==0){
		cycleStructureReps.push_back( CSToPerm(curr) );
		return;
	}

	if (ctr<= rem){
		//first we select it and see what we get
		curr[ctr]++;
		rem = rem- ctr;
		genCRS(curr, ctr, rem );
		//then we dont
		curr[ctr]--;
		rem = rem+ctr;
		genCRS(curr, ctr+1, rem );
	}
	return;
}

void MOLS::getCycleStructure(permutation &p, vector<int> &nofCycles ){
	vector<bool> visited (n+1, false);
	int i=0, j, front, count=0;

	for (i=0; i<n; i++){
		//cout<<i<<endl;
		if (i==p[i]){
			nofCycles[1]++;
 		}
		else
		{
			if (!visited[i]){
				j=i+0;
				front =j;
				count=0;
				do{
					count++;
					visited[j] = true;
					j = p[j];
				}while (front !=j);
				//visited[j] = true;
				nofCycles[count]++;
 			}
		}
	}
}

void MOLS::getCycleStructure(permutation &p, vector<int> &nofCycles, map<int, vector<vector<int> > > &cycleLength_cycles_map ){
	//bool visited[n+1] ={0};
	vector<bool> visited (n+1, false);
	permutation::const_iterator permIt;
	vector<int> permArr(n+1, 0);
	//int permArr[n]={0};
	int i=0;

	for (i=0; i<n; i++){
		//cout<<i<<endl;
		if (i==p[i]){
			nofCycles[1]++;
			vector<int> thisCycle;
			thisCycle.push_back(i);
			cycleLength_cycles_map[thisCycle.size()].push_back(thisCycle);
		}
		else
		{
			if (!visited[i]){
				vector<int> thisCycle;
				int j=i+0;
				do{
					thisCycle.push_back(j);
					visited[j] = true;
					j = p[j];
				}while (thisCycle.front() !=j);

				nofCycles[thisCycle.size()]++;
				cycleLength_cycles_map[thisCycle.size()].push_back(thisCycle);
			}
		}
	}

	//for (i=0; i<n+1; i++)
	//	cout<< nofCycles[i]<< " ";
}

/*list<permutation> MOLS::genRelevantPermutations(permutation p){
	list<permutation> result;
	int i;
	vector<int> nofCycles(n+1, 0);
	//int nofCycles[n+1] ={0};
	map<int, vector<vector<int> > > cycleLength_cycles_map;

	for (i=0; i<n+1; i++){
		list<list<int> > l;
		cycleLength_cycles_map[i] = l;
	}
	getCycleStructure(p, nofCycles, cycleLength_cycles_map);

	//for some reason this doesnt work if there are more than  one cycle length 1
	//TODO fis, for now just return empty list
	if (nofCycles[1]>1)
		return result ;

	int nofCycleLengths = 0;
	for (i=0; i<n+1; i++)
		if (nofCycles[i] >0)
			nofCycleLengths++	;

	list<permutation> arrPossPerLength[nofCycleLengths];

}*/

vector<permutation> MOLS::genRelevantPermutations(permutation &p){
	vector<permutation> result;
	int i;
	vector<int> nofCycles(n+1, 0);
	//int nofCycles[n+1] ={0};
	map<int, vector<vector<int> > > cycleLength_cycles_map;

	/*for (i=0; i<n+1; i++){
		list<list<int> > l;
		cycleLength_cycles_map[i] = l;
	}*/
	getCycleStructure(p, nofCycles, cycleLength_cycles_map);

	//for some reason this doesnt work if there are more than  one cycle length 1
	//TODO fis, for now just return empty list
	if (nofCycles[1]>1)
		return result ;

	int nofCycleLengths = 0;
	for (i=0; i<n+1; i++)
		if (nofCycles[i] >0)
			nofCycleLengths++	;

	vector<permutation> arrPossPerLength[nofCycleLengths];

/*	list<permutation> l;
	for (i=0; i<nofCycleLengths; i++){
		arrPossPerLength[i]=l;
	}*/

	int length=0;
	int ctr =0;
	for (length = 1; length<n+1; length++){
		vector<permutation> arrPossibilities[nofCycles[length]];

		if (nofCycles[length]>0){
			int ind=0;
			//list<permutation> arrPossibilities[nofCycles[length]];

			lsquare::iterator cycleIt;
			permutation::iterator permIt;
			//all cycles of this length
			//cout<< "Cycles of length "<< length <<": "<<endl;
			for (cycleIt = cycleLength_cycles_map[length].begin(); cycleIt !=  cycleLength_cycles_map[length].end(); ++cycleIt){
					int a=0;
 					for (permIt = (*cycleIt).begin(); a<length; ++permIt){
						std::rotate((*cycleIt).begin(), permIt, (*cycleIt).end()) ;
						// printPerm(*cycleIt);

						arrPossibilities[ind].push_back(copyPerm(*cycleIt)			);
						if (a++>0)
							--permIt;
					}
					ind++;
			}

			if (nofCycles[length] == 1){ //no need to worry about shuffling permutations, there is only one
				for (cycleIt = arrPossibilities[0].begin(); cycleIt != arrPossibilities[0].end(); ++cycleIt){
					arrPossPerLength[ctr].push_back(copyPerm(*cycleIt));
				}
				ctr++;

			}
			else{ //now we need to both shuffle and select

				vector<int> orderOfCycle;
				for (int i=0; i<nofCycles[length]; i++)
					orderOfCycle.push_back(i);

				map<int, list<permutation> > cperm_to_possibilities_map;
				int pno = 0;

				do{
					vector<permutation> templ;
					vector<permutation> templ2;

					for (ind=0; ind<nofCycles[length]; ind++){
						//blist<permutation> templ;
						if (ind==0){
							for (cycleIt = arrPossibilities[orderOfCycle[ind]].begin(); cycleIt != arrPossibilities[orderOfCycle[ind]].end(); ++cycleIt){
								permutation a = copyPerm(*cycleIt);
								templ.push_back(a);
								//printPerm(a);
							}
						}
						else{
							//cout<<endl<< "W"<< cperm_to_possibilities_map[pno-1].size();

							lsquare::const_iterator partialpermIt;
							//while nextpermutation
 							for (partialpermIt = templ.begin(); partialpermIt != templ.end(); ++partialpermIt){

								for (cycleIt = arrPossibilities[orderOfCycle[ind]].begin(); cycleIt != arrPossibilities[orderOfCycle[ind]].end(); ++cycleIt){
									permutation a(*partialpermIt);
								 	a.insert(a.end(), (*cycleIt).begin(), (*cycleIt).end());
									templ2.push_back(a);
 								}
							}

							templ.clear();
							templ.assign(templ2.begin(), templ2.end());
							templ2.clear();
						}//end else

					}
					cperm_to_possibilities_map[pno++].assign(templ.begin(), templ.end());
				}
				while (next_permutation(orderOfCycle.begin(), orderOfCycle.end()	));

				//merge all the different list , one for each permutation, into one list
				int y;
				for (y=0; y<pno; y++	)
					arrPossPerLength[ctr].insert(arrPossPerLength[ctr].end(),cperm_to_possibilities_map[y].begin(), cperm_to_possibilities_map[y].end());

				ctr++;
			}

			//now create options for this length and every other length
			// options[length]
			//then pick one of the options per lenght

		} //end if there are cycles of this length


	} //end for every length

	/*now we assume we have all the combinations of every length in an array, for example
	 * <0> 	| <1234> 	| <567>
	 * 		| <1243>	| <675>
	 * 		| <2134>	| <756>
	 * 		| <2143>
	 * 		| <3412>
	 * 		| <4312>
	 * 		| <3421>
	 * 		| <4321>
	 *
	 * 	Select one from each column of arrposperlength
	 */

/*	int ind = 0;
	while (ind<nofCycleLengths){
		cout<<ind<<endl;
		list<list<int> >::iterator cycleIt;

		for (cycleIt = arrPossPerLength[ind].begin(); cycleIt != arrPossPerLength[ind].end(); ++cycleIt){

			printPerm(*cycleIt);
		}

		ind++;
	}*/

	////endprint
	//cout<<"==========";

	list<permutation> tempres;
	permutation b;
	result.push_back(b);

	vector<vector<int> >::const_iterator partialpermIt;
	vector<vector<int> >::const_iterator cycleIt;
	vector<permutation> temp;
	int ind = 0;
	while (ind<nofCycleLengths){
		temp.clear();
		for (partialpermIt = result.begin(); partialpermIt != result.end(); ++partialpermIt){

			for (cycleIt = arrPossPerLength[ind].begin(); cycleIt != arrPossPerLength[ind].end(); ++cycleIt){
				permutation a(*partialpermIt);
				a.insert(a.end(), (*cycleIt).begin(), (*cycleIt).end());
				temp.push_back(a);
				//printPerm(a);
			}
		}
		ind++;
		result = temp;
	}
	//possibleShuffles = result;
	return result;


}

bool MOLS::isIdentity(permutation p){
	//int n = p.size();
	permutation::iterator permIt ;

	int checkValue = 0;
	for (permIt = p.begin(); permIt != p.end(); ++permIt){
		if (checkValue++ != *permIt)
			return false;
	}
	cout<< "ID TRUE";
	return true;

}


void MOLS::changeOrder(vector<lsquare> &partMOLS , vector<int> &order , vector<lsquare> &newMOLS){
	int i;
	for (i=0; i<k; i++){
		newMOLS[i]= partMOLS[order[i]];
	}
}

vector<lsquare> MOLS::changeOrder(vector<lsquare> &pMOLS , int first, int second  ){
	lsquare emptyLS;
	vector<lsquare> newMOLS(k, emptyLS);
	newMOLS[0]=pMOLS[first];
	newMOLS[1]= pMOLS[second];
 	int i=0;
	int j=2;
	for (i=0; i<k; i++){
		if (i!=first && i!=second)
			 newMOLS[j++]= pMOLS[i];
	}

	return newMOLS;
}

void MOLS::rollUp(vector<lsquare> &pMOLS ,  int pos,  permutation &origP ){

	int i;
	for (i=0; i<k; i++){
 		for (unsigned int j=0; j< pMOLS[i].size(); j++){
			for (unsigned int jj=0 ; jj<n ; jj++){
				pMOLS[i][j][jj] = (pMOLS[i][j][jj]-pos+n)%n;
			}

		}
	}
	return;
}

void MOLS::standardForm(vector<lsquare> &pMOLS , int rowUp, int pos ,permutation &origP){

	rollUp(pMOLS, pos, origP);
	/*cout<<"after rollup2";
	printMOLSPerms(pMOLS);*/

	permutation p(n);
	for (unsigned int i=0;i<n;i++){
 			p[i] = (origP[i]-pos+n)%n;
	}

	//TODO maybe unify this with testpermute
	permuteMOLS(pMOLS, p, identity);

}

vector<int> MOLS::rowMeets(permutation &p1, permutation &p2){

 	vector<int>	res(2);
 	for (unsigned int i=0; i<p1.size(); i++) {
		if (p1[i]==p2[i]){
			res[0] = i; //row
			res[1] = p1[i];
			return res;
		}

	}
	return res;
}

/*
 * accept  a permutation p, of which the cycle structure is to be found and comparded to
 * another cycle structure in the form of an integer arrayu where
 * the number stored at index i is the number of cycles of length i
 * Returns 0 if the same number of cycles of every length
 *
 * Returns <0 if p < target
 * Returns >0 if p> target
 */
int MOLS::compareCS(permutation &p , vector<int> &targetCS){
	vector<int> thisCS(n+1, 0);
	//map<int, vector<vector<int> > > dummy_cycles_map;
	getCycleStructure(p, thisCS);
 	/*cout<<endl;printPerm(p);
	for(j=0; j<n;j++){
		cout<<thisCS[j]<<" ";}
	cout<<endl;
	for(j=0; j<n;j++){
			cout<<targetCS[j]<<" ";}
	cout<<endl;*/

	for (unsigned int i=0; i<n; i++){
		//cout<<targetCS[i]<<"<>" <<thisCS[i];
		if (! (targetCS[i] == thisCS[i])	)
			return targetCS[i]-thisCS[i];
	}
	return 0; //the same
}

//boolean just says that p is alreasy a CS, not just a perm
int MOLS::compareCS(permutation &pCS , vector<int> &targetCS, bool t){
	//printPerm(pCS); printPerm(targetCS);
	for (unsigned int i=0; i<n; i++){
		//cout<<targetCS[i]<<"<>" <<thisCS[i];
		if (! (targetCS[i] == pCS[i])	)
			return targetCS[i]-pCS[i];
	}

	return 0; //the same

}

list<vector<int> > MOLS::getSmallRelCS( vector<lsquare> &pMOLS, list<permutation> &listP){
	permutation targetRCS(pMOLS[1].front());

	/*vector<int> targetCS(n+1, 0);
	map<int, list<list<int> > > dummy_cycles_map;
	getCycleStructure(targetRCS, targetCS, dummy_cycles_map );*/

	list<vector<int> > squares;
 	list<permutation>::iterator permIt1;
	list<permutation>::iterator permIt2;
	vector<int> data(4) ;
	vector<int> data1(4) ;
	vector<int> dataPos(2);
 	for(unsigned int i=0; i<k-1; i++){
		for(unsigned int j=i+1; j<k; j++){

		 	for (unsigned int ii=0; ii< pMOLS[i].size(); ii++){
				for (unsigned int jj=0; jj< pMOLS[j].size(); jj++){

					permutation rcsV = rcs(pMOLS[i][ii], pMOLS[j][jj]);

					//if this rcs is equal to smallest store in list
					//if (!lexicographical_compare(rcsV.begin(), rcsV.end(), targetRCS.begin(), targetRCS.end())
					//		and !lexicographical_compare(targetRCS.begin(), targetRCS.end(), rcsV.begin(), rcsV.end())){
					int comparison = compareCS(rcsV, currentCS);

					if (comparison==0){
					// if (true){

						//data.assign(4);
						data[0]=i;
						data[1] =j;
						dataPos = rowMeets(pMOLS[i][ii], pMOLS[j][jj]);
						data[2] = dataPos[0];
						data[3] = dataPos[1];
						//if (!(data[0]==0&&data[1]==1&&data[2]==0) ){
						listP.push_back(pMOLS[i][ii]);
						//cout<<  data[0]<<" "<<data[1]<<" "<<data[2]<<" "<<data[3]<< " perm "; printPerm(*permIt1);cout<<", "; printPerm(*permIt2); cout<<endl;
						squares.push_back(data);

						 //and in reverse
						data[0]=j;
						data[1] =i;
						dataPos = rowMeets(pMOLS[j][jj],pMOLS[i][ii]);
						data[2] = dataPos[0];
						data[3] = dataPos[1];
						// if (!(data1[0]==0&&data1[1]==1&&data1[2]==0) ){
						listP.push_back(pMOLS[j][jj]);
						//cout<<  data1[0]<<" "<<data1[1]<<" "<<data1[2]<<" "<<data1[3]<< " perm "; printPerm(*permIt2);cout<<", "; printPerm(*permIt1); cout<<endl;

						squares.push_back(data);

					}



				}
			}


		}
	}
 }
bool MOLS::getSmallRelCS( vector<lsquare>& pMOLS, permutation& mapTo ){
 	//permutation targetRCS(pMOLS[1].front());
	RCSperms.clear();
	RCSsquares.clear();
	int i,j;

 	vector<vector<vector<vector<int> > > >  tempPerms(k, vector<vector<vector<int> > >(k  ));
	vector<vector<vector<vector<int> > > >  tempSquares(k, vector<vector<vector<int> > >(k  ));
	vector<int> data(4) ;
	vector<int> dataPos(2);

 	for(i=0; i<k-1; i++){
		for(j=i+1; j<k; j++){

		 	for (unsigned int ii=0; ii< pMOLS[i].size(); ii++){
				for (unsigned int jj=0; jj< pMOLS[j].size(); jj++){


					permutation rcsV = rcs(pMOLS[i][ii], pMOLS[j][jj]);

					//if this rcs is equal to smallest store in list
					//if (!lexicographical_compare(rcsV.begin(), rcsV.end(), targetRCS.begin(), targetRCS.end())
					//		and !lexicographical_compare(targetRCS.begin(), targetRCS.end(), rcsV.begin(), rcsV.end())){
					int comparison = compareCS(rcsV, mapTo);
					if (comparison<0){
 						return false;
					}
					if (comparison==0){
					// if (true){

						data[0]=i;
						data[1] =j;
					    dataPos = rowMeets(pMOLS[i][ii], pMOLS[j][jj]);
						data[2] = dataPos[0];
						data[3] = dataPos[1];
						//if (!(data[0]==0&&data[1]==1&&data[2]==0) ){
						tempPerms[data[0]][data[1]].push_back(pMOLS[i][ii]);
						//RCSperms.push_back(pMOLS[i][ii]);
						//cout<<  data[0]<<" "<<data[1]<<" "<<data[2]<<" "<<data[3]<< " perm "; printPerm(*permIt1);cout<<", "; printPerm(*permIt2); cout<<endl;
						tempSquares[data[0]][data[1]].push_back(data);
						//RCSsquares.push_back(data);

						//and in reverse

						data[0]=j;
						data[1] =i;
						dataPos = rowMeets(pMOLS[j][jj],pMOLS[i][ii]);
						data[2] = dataPos[0];
						data[3] = dataPos[1];
						// if (!(data1[0]==0&&data1[1]==1&&data1[2]==0) ){
						//RCSperms.push_back(pMOLS[j][jj]);
						//cout<<  data1[0]<<" "<<data1[1]<<" "<<data1[2]<<" "<<data1[3]<< " perm "; printPerm(*permIt2);cout<<", "; printPerm(*permIt1); cout<<endl;

						//RCSsquares.push_back(data1);

						tempPerms[data[0]][data[1]].push_back(pMOLS[j][jj]);
						tempSquares[data[0]][data[1]].push_back(data);
					}

					/*rcsV = rcs(pMOLS[j][jj],pMOLS[i][ii]);

					if (compareCS(rcsV, currentCS)==0){
					// if(true){
						vector<int> data1 ;
						data1.assign(4, 0);
						data1[0]=j;
						data1[1] =i;
						vector<int> dataPos = rowMeets(pMOLS[j][jj],pMOLS[i][ii]);
						data1[2] = dataPos[0];
						data1[3] = dataPos[1];
						// if (!(data1[0]==0&&data1[1]==1&&data1[2]==0) ){
							listP.push_back(pMOLS[j][jj]);
							//cout<<  data1[0]<<" "<<data1[1]<<" "<<data1[2]<<" "<<data1[3]<< " perm "; printPerm(*permIt2);cout<<", "; printPerm(*permIt1); cout<<endl;

							squares.push_back(data1);
						// }
					}*/
				}
			}
		}
	}
 	//cout<<"begin fors"<<endl;
 	for(i=0; i<k ; i++){
 		for(j=0; j<k; j++){
 			if (tempPerms[i][j].size()>0){
 				RCSperms.insert( RCSperms.end(), tempPerms[i][j].begin(), tempPerms[i][j].end());
 				RCSsquares.insert(RCSsquares.end(),tempSquares[i][j].begin(), tempSquares[i][j].end() );
 			}
 		}
 	} 	//cout<<"end fors"<<RCSperms.size()<<" "<<RCSsquares.size()<<endl;
 /*	for(i=0; i<RCSsquares.size() ; i++){
 	 	 cout<<RCSsquares[i][0]<<" "<<RCSsquares[i][1]<<endl;
 	 	}*/

 	return true;
}

bool MOLS::getSmallRelCS( vector<lsquare>& pMOLS ){
 	//permutation targetRCS(pMOLS[1].front());
	RCSperms.clear();
	RCSsquares.clear();
	int i,j;

 	vector<vector<vector<vector<int> > > >  tempPerms(k, vector<vector<vector<int> > >(k  ));
	vector<vector<vector<vector<int> > > >  tempSquares(k, vector<vector<vector<int> > >(k  ));
	vector<int> data(4) ;
	vector<int> dataPos(2);

 	for(i=0; i<k-1; i++){
		for(j=i+1; j<k; j++){

		 	for (unsigned int ii=0; ii< pMOLS[i].size(); ii++){
				for (unsigned int jj=0; jj< pMOLS[j].size(); jj++){


					permutation rcsV = rcs(pMOLS[i][ii], pMOLS[j][jj]);

					//if this rcs is equal to smallest store in list
					//if (!lexicographical_compare(rcsV.begin(), rcsV.end(), targetRCS.begin(), targetRCS.end())
					//		and !lexicographical_compare(targetRCS.begin(), targetRCS.end(), rcsV.begin(), rcsV.end())){
					int comparison = compareCS(rcsV, currentCS);
					if (comparison<0){
 						return false;
					}
					if (comparison==0){
					// if (true){

						data[0]=i;
						data[1] =j;
					    dataPos = rowMeets(pMOLS[i][ii], pMOLS[j][jj]);
						data[2] = dataPos[0];
						data[3] = dataPos[1];
						//if (!(data[0]==0&&data[1]==1&&data[2]==0) ){
						tempPerms[data[0]][data[1]].push_back(pMOLS[i][ii]);
						//RCSperms.push_back(pMOLS[i][ii]);
						//cout<<  data[0]<<" "<<data[1]<<" "<<data[2]<<" "<<data[3]<< " perm "; printPerm(*permIt1);cout<<", "; printPerm(*permIt2); cout<<endl;
						tempSquares[data[0]][data[1]].push_back(data);
						//RCSsquares.push_back(data);

						//and in reverse

						data[0]=j;
						data[1] =i;
						dataPos = rowMeets(pMOLS[j][jj],pMOLS[i][ii]);
						data[2] = dataPos[0];
						data[3] = dataPos[1];
						// if (!(data1[0]==0&&data1[1]==1&&data1[2]==0) ){
						//RCSperms.push_back(pMOLS[j][jj]);
						//cout<<  data1[0]<<" "<<data1[1]<<" "<<data1[2]<<" "<<data1[3]<< " perm "; printPerm(*permIt2);cout<<", "; printPerm(*permIt1); cout<<endl;

						//RCSsquares.push_back(data1);

						tempPerms[data[0]][data[1]].push_back(pMOLS[j][jj]);
						tempSquares[data[0]][data[1]].push_back(data);
					}

					/*rcsV = rcs(pMOLS[j][jj],pMOLS[i][ii]);

					if (compareCS(rcsV, currentCS)==0){
					// if(true){
						vector<int> data1 ;
						data1.assign(4, 0);
						data1[0]=j;
						data1[1] =i;
						vector<int> dataPos = rowMeets(pMOLS[j][jj],pMOLS[i][ii]);
						data1[2] = dataPos[0];
						data1[3] = dataPos[1];
						// if (!(data1[0]==0&&data1[1]==1&&data1[2]==0) ){
							listP.push_back(pMOLS[j][jj]);
							//cout<<  data1[0]<<" "<<data1[1]<<" "<<data1[2]<<" "<<data1[3]<< " perm "; printPerm(*permIt2);cout<<", "; printPerm(*permIt1); cout<<endl;

							squares.push_back(data1);
						// }
					}*/
				}
			}
		}
	}
 	//cout<<"begin fors"<<endl;
 	for(i=0; i<k ; i++){
 		for(j=0; j<k; j++){
 			if (tempPerms[i][j].size()>0){
 				RCSperms.insert( RCSperms.end(), tempPerms[i][j].begin(), tempPerms[i][j].end());
 				RCSsquares.insert(RCSsquares.end(),tempSquares[i][j].begin(), tempSquares[i][j].end() );
 			}
 		}
 	} 	//cout<<"end fors"<<RCSperms.size()<<" "<<RCSsquares.size()<<endl;
 /*	for(i=0; i<RCSsquares.size() ; i++){
 	 	 cout<<RCSsquares[i][0]<<" "<<RCSsquares[i][1]<<endl;
 	 	}*/

 	return true;
}

permutation MOLS::flatten(permutation &pNow){

	//cout<<" Flatten in - "; printPerm(pNow);cout.flush();
	vector<bool> visited (n+1, false);
	permutation z, y;
	vector<permutation> size_cycles(n
                                 );

	int j;
	for (unsigned int i=0; i<n; i++){
		if (!visited[i]){
			//bcout<<i<<endl;
			vector<int > thiscycle;
			j=i;
			do{

				visited[j] = true;
				//cout<<"push "<<pNow[j]<<endl;
				thiscycle.push_back(pNow[j]);

				//z.push_back(pNow[j]);
				j = pNow[j];
			}while (!visited[j]);
			size_cycles[thiscycle.size()].insert(size_cycles[thiscycle.size()].begin(), thiscycle.begin(), thiscycle.end()	);

		}
	}
	for (unsigned int i=1; i<n; i++){
			 if (size_cycles[i].size()>0)
				 y.insert(y.end(), size_cycles[i].begin(), size_cycles[i].end());
		}



	//cout<<" Flatten outut - "; printPerm(z);cout<<endl;
	return y;



}

list<permutation> MOLS::getShuffles(permutation &pOrig, permutation &pNow){

	//vector<permutation> listShuffles ;
	// listShuffles=genRelevantPermutations(pNow);
	list<permutation> listtemp;

	vector<permutation>::const_iterator pi;
	//cout<< "listShuffles size - "<< listShuffles.size()<<endl;
	permutation f = flatten(pNow);
	//for (pi=listShuffles.begin(); pi!=listShuffles.end(); ++pi){
	//permutation z(n);
	for (pi=possibleShuffles.begin(); pi!=possibleShuffles.end(); ++pi){
		 for (unsigned int i=0; i<n; i++)
			 z[f[(*pi)[i]]] =i;

		 listtemp.push_back( z);
	}

	return listtemp;
}

bool MOLS::noSmallerRCS(permutation &smallestRCS, vector<lsquare> &pMOLS){
 	vector< vector<int> >::const_iterator testRel;
	vector< permutation >::iterator testPermsit;
	testPermsit= RCSperms.begin();

	list<permutation>::iterator permIt; int i;
	vector<int> currentOrder(2, 0);

	vector<lsquare> ordMOLS ;

	for (testRel=RCSsquares.begin(); testRel!=RCSsquares.end(); testRel++){
		// cout<<endl<<" "<<(*testRel)[0]<<" "<<(*testRel)[1]<<" "<<(*testRel)[2]<<" "<<(*testRel)[3] <<" ";

		if ((*testRel)[0]!= currentOrder[0] || currentOrder[1]!=(*testRel)[1]){
			currentOrder[0]=(*testRel)[0];
			currentOrder[1]=(*testRel)[1];
			ordMOLS = changeOrder(pMOLS, (*testRel)[0],(*testRel)[1]);
		}


		vector<lsquare> roMOLS(ordMOLS);

		standardForm(roMOLS, (*testRel)[2],(*testRel)[3], *testPermsit);
		/*cout<<"Before";
		printMOLSPerms(partMOLS );
		cout<<"Stadard"<<endl;
		printMOLSPerms(roMOLS );

		printMOLS(roMOLS);
		printMOLSPerms(roMOLS);

		cout<<"STD Form - "<<endl;printMOLSPerms(roMOLS);cout.flush();*/

		list<permutation> allShuffles = getShuffles(smallestRCS, roMOLS[1].front());
		for (permIt = allShuffles.begin(); permIt != allShuffles.end(); ++permIt) {
			//vector<lsquare> tMOLS(roMOLS);
			/*//permutation realShuffle = composition(inverse( *permIt),targetP );
			//printPerm( *permIt);
			//permuteMOLS(tMOLS, *permIt, *permIt);
 			b*/

			//if (!comparePartialMOLS(partMOLS, tMOLS) )// if  part > temp
			if (!testPermuteMOLS(roMOLS, *permIt, *permIt) )// if  part > temp
			{
				/*for (i=0; i<k; i++){
					printPerm(partMOLS[i].front());cout<<" ";}*/
				/*cout<<"***********************"<<endl;
				printMOLSPerms(partMOLS);
				cout<<"v"<<endl;
				printMOLSPerms(ttMOLS);*/
				//cout<< "rel CS smaller.."<< (*testRel)[0]<<" "<<(*testRel)[1]<<" "<<(*testRel)[2]<<" "<<(*testRel)[3]<< " perm "; printPerm(*permIt);cout<<endl;

				return false;
			}

		}
		testPermsit++;
	}


	//==================\
	//TRANSPOSES??
 	 vector<lsquare> transpMOLS = transposeMOLS(pMOLS);
 	//list<permutation> listPerms1;
	//list <vector<int> > relCS1 = getSmallRelCS(transpMOLS, listPerms1);
	vector< vector<int> >::const_iterator testRel1;
	vector< permutation >::iterator testPermsit1;
	testPermsit1= RCSperms.begin();
	//list<permutation>::iterator permIt;int i;
	vector<lsquare> ordMOLS1;
	currentOrder[0] =0;
	currentOrder[1] = 0;
	permutation inv;

	for (testRel1=RCSsquares.begin(); testRel1!=RCSsquares.end(); testRel1++){

		if ((*testRel1)[0]!= currentOrder[0] || currentOrder[1]!=(*testRel1)[1]){
					currentOrder[0]=(*testRel1)[0];
					currentOrder[1]=(*testRel1)[1];
					ordMOLS1 = changeOrder(transpMOLS, (*testRel1)[0],(*testRel1)[1]);

		}

		vector<lsquare> roMOLS(ordMOLS1);


		//vector<lsquare> roMOLS(transpMOLS);

		//changeOrder(roMOLS, (*testRel1)[0],(*testRel1)[1]);
		inv = inverse(*testPermsit1);
		standardForm(roMOLS, (*testRel1)[3],(*testRel1)[2], inv );

		permutation currPerm = roMOLS[1].front();
		list<permutation> allShuffles = getShuffles(smallestRCS,currPerm  );

		for (permIt = allShuffles.begin(); permIt != allShuffles.end(); ++permIt) {
		//	vector<lsquare> tMOLS(roMOLS) ;
			//permuteMOLS(tMOLS, *permIt, *permIt );


			//if (!comparePartialMOLS(partMOLS, tMOLS) )// if  part > temp
			if (!testPermuteMOLS(roMOLS, *permIt, *permIt ))
			{
				/*for (i=0; i<k; i++){
					printPerm(compMOLS[i].back());cout<<" ";}
*/
				//cout<< "Initial";printMOLSPerms(compMOLS, k);
				//cout<< "New after transpose";
				//printMOLSPerms(tMOLS, k);
				//cout<<"***********************"<<endl;

				//cout<< "rel CS smaller after transpose.."<< (*testRel1)[0]<<" "<<(*testRel1)[1]<<" "<<(*testRel1)[2]<<" "<<(*testRel1)[3]<< " perm "; printPerm(*permIt);cout<<endl;
								//printMOLS(tMOLS, k);
				return false;
			}
		}
		testPermsit1++;
	}

	return true;
}

bool MOLS::noSmallerRCS(permutation &smallestRCS){
	return noSmallerRCS(smallestRCS, partMOLS);
}


bool MOLS::isSmallest4( ){

if (currSquare == k-1	){
	numIsSmallest[partMOLS[k-1].size()]++;

	if (!noSmallerRCS( partMOLS[1].front() )){
		/*for (i=0; i<k; i++)
				printPerm(partMOLS[i].back());
		cout << "smaller cs"<<endl;*/
			return false;
		}
numIsSmallestTrue[partMOLS[k-1].size()]++;
}
	return true;

}

bool MOLS::isSmallest4(bool t ){
	return true;

	/*if (partMOLS[currSquare].size()>2)
		return true;*/
	if (currSquare == k-1	)	{
		numIsSmallest[partMOLS[k-1].size()]++;
		if (!noSmallerRCS( partMOLS[1].front() )){
		/*for (i=0; i<k; i++)
				printPerm(partMOLS[i].back());
		cout << "smaller cs"<<endl;*/
			return false;
		}

	numIsSmallestTrue[partMOLS[k-1].size()]++;
	}
	return true;

}

bool MOLS::isSmallest4(permutation &smallestRCS, vector<lsquare> &pMOLS){

	numIsSmallest[partMOLS[k-1].size()]++;
	if (!noSmallerRCS(smallestRCS, pMOLS ) ){
		return false;
	}
	numIsSmallestTrue[partMOLS[k-1].size()]++;

	return true;
}

bool MOLS::isSmallestConjugateMOLS(vector<lsquare> &pMOLS ){
	int maxCS=cycleStructureReps.size();
	currentCS.clear();currentCS.resize(n+1, 0);
	vector<lsquare> nMOLS(k, vector<vector<int> >(n, vector<int>(n,0)));
	vector<vector<int> > oa(k+2, vector<int>(n*n, 0));
	//map<int, vector<vector<int> > > dummy_cycles_map;
	//map<int, vector<vector<int> > > dummy_cycles_map2;
	vector<int> thisCS(n+1, 0);
	getCycleStructure(partMOLS[1].front(),  currentCS );
	//cout<<endl<<"Current cs ";printPerm(partMOLS[1].front());cout <<": "; printPerm(currentCS);

	//populate the orthoganal array from the current MOLS
	for (int i =0; i<n*n; i++){
		 oa[0][i] = i/n;
		 oa[1][i] = i%n;
		 for (int j =0; j<k; j++){
			 oa[2+j][i] = currLS[j][oa[0][i]][oa[1][i]];
		 }
	}
	/*for (int i =0; i<k+2; i++){
		for (int j =0; j<n*n; j++)
			cout<< oa[i][j]<<" ";
		cout<<endl;
	}*/
	permutation oaPerm;
	for (int i=0; i<k+2; i++)
		oaPerm.push_back(i);

	for (unsigned int i =0; i< cycleStructureReps.size(); i++){
		getCycleStructure(cycleStructureReps[i],  thisCS );
		if (compareCS(thisCS, currentCS, true)> 0){
			 maxCS=i; break;
		}
		thisCS.clear();thisCS.resize(n+1, 0);
	}

	/*vector<vector<int> >listPerms((k*k+3*k)/2, oaPerm );
	int temp;
	for (int i =0; i<k; i++){
		temp = listPerms[i][0];
		listPerms[i][0] = listPerms[i][i+2];
		listPerms[i][i+2]= temp;

		temp = listPerms[i+k][1];
		listPerms[i+k][1] = listPerms[i+k][i+2];
		listPerms[i+k][i+2]= temp;
	}
	int counter=2*k;
	for (int i =0; i<k-1; i++){
		for (int j =i+1; j<k; j++){
			temp = listPerms[counter][0];
			listPerms[counter][0] = listPerms[counter][i+2];
			listPerms[counter][i+2]= temp;

			temp = listPerms[counter][1];
			listPerms[counter][1] = listPerms[counter][j+2];
			listPerms[counter][j+2]= temp;
			counter++;
		}
	}
	for (int ii =0; ii<listPerms.size(); ii++){
		oaPerm = listPerms[ii];
		for (int j =0; j<k; j++){
					for (int i =0; i<n*n; i++){
						// the new value, number of universal is	oa[oaPerm[k]][i]
						// the new row (position in universal ) is oa[oaPerm[0]][i]
						// the new column, value in that position in oa[oaPerm[1]][i]
						nMOLS[j][oa[oaPerm[j+2]][i]][oa[oaPerm[0]][i]] = oa[oaPerm[1]][i] ;
					}
				}
				//printMOLSPerms(nMOLS);

				for (unsigned int i =0; i< maxCS; i++){
					//we want to try and map it to all of the smaller and equal csreps
					getCycleStructure(cycleStructureReps[i],  thisCS, dummy_cycles_map );
					//int cnt = compareCS(thisCS, currentCS, true); cout<<"Z"<<cnt;
					if (getSmallRelCS(nMOLS,  thisCS )){
						if (!isSmallest4(cycleStructureReps[i], nMOLS)){
							//currentCS.clear();currentCS.resize(n+1, 0);
							return false;
						}
					}
					thisCS.clear(); thisCS.resize(n+1, 0);

					//dummy_cycles_map.clear();
				}

	}*/


	while (next_permutation(oaPerm.begin(), oaPerm.end())){
		//cout<<"permutations "; printPerm(oaPerm);cout<<endl;
		for (int j =0; j<k; j++){
			for (int i =0; i<n*n; i++){
				// the new value, number of universal is	oa[oaPerm[k]][i]
				// the new row (position in universal ) is oa[oaPerm[0]][i]
				// the new column, value in that position in oa[oaPerm[1]][i]
				nMOLS[j][oa[oaPerm[j+2]][i]][oa[oaPerm[0]][i]] = oa[oaPerm[1]][i] ;
			}
		}
		//printMOLSPerms(nMOLS);

		for (unsigned int i =0; i< maxCS; i++){
			//we want to try and map it to all of the smaller and equal csreps
			getCycleStructure(cycleStructureReps[i],  thisCS );
			//int cnt = compareCS(thisCS, currentCS, true); cout<<"Z"<<cnt;
			if (getSmallRelCS(nMOLS,  thisCS )){
				if (!isSmallest4(cycleStructureReps[i], nMOLS)){
					//currentCS.clear();currentCS.resize(n+1, 0);
					return false;
				}
			}
			thisCS.clear(); thisCS.resize(n+1, 0);

			//dummy_cycles_map.clear();
		} //end for csreps
	}//end while
	//currentCS.clear();currentCS.resize(n+1, 0);
	return true;
}

void MOLS::buildCurrentLS(){
		//currLS[]
		int x, y;
		for(x = 0; x < n; x ++) {
		    for(y = 0; y < n; y ++)
				currLS[currSquare][x][y] = -1;
		}


		list<permutation >::iterator universalIt; //iterates over universals
		permutation::iterator perm_it; //iterates over permutation

		for (unsigned int i=0;i< partMOLS[currSquare].size(); i++){
			for (unsigned int j=0; j<partMOLS[currSquare][i].size(); j++){
				currLS[currSquare][j][partMOLS[currSquare][i][j]] = i;
			}
		}
}

bool MOLS::checkFit( permutation &p){
	if (n==0) return true;
	int i= 0; //row



	i=0;
//	printMOLS(partMOLS); cout<<"Try to fit perm "; printPerm(p);cout<<endl;
	for (unsigned int i=0;i< p.size(); i++){
	//for(perm_it=p.begin(); perm_it != p.end(); ++perm_it){
		if (currLS[currSquare][i][p[i]] < n && currLS[currSquare][i][p[i]] >=0)
			return false; //if the cell is already filled
	}
	return true;

}

bool MOLS::rcsOrthog(permutation &p1, permutation &p2){
 	int count=0;
	for (unsigned int i=0; i<p1.size(); i++)
		if (p1[i] == p2[i])
			count++;

	return (count==1);

}

bool MOLS::checkOrthogonal(permutation &P){
	if (nof_calls_made++ > nof_calls_limit&& send_checkpoint<1) send_checkpoint = 1;
	
	unsigned int sz, count;
	for (unsigned int i=0; i<k; i++){
		if (i!=currSquare){ //for every square except the current one
			sz = partMOLS[i].size();
			for (unsigned int j =0; j<sz; j++){
                 count=0;
                for (unsigned int l=0; l<n; l++)
                    if (partMOLS[i][j][l] == P[l])
                        count++;
                if (count != 1)
                    return false;

			}

//            for (unsigned int j =0; j<sz; j++)
//				if(!rcsOrthog(partMOLS[i][j], P))
//					return false;

		}
	}
	return true;

//    unsigned int sz, count;
//    vector<int> met(n,-1);
//    for (unsigned int i=0; i<k; i++){
//
//		if (i!=currSquare){ //for every square except the current one
//			sz = partMOLS[i].size();
//			count = 0;
//			//vector<bool> met(n, false);
//			for (unsigned int j =0; j<n; j++)
//                if (currLS[i][j][P[j]] != n){
//                   if  (met[currLS[i][j][P[j]]] == i)
//                        return false;
//                   else{
//                        met[currLS[i][j][P[j]]]= i;
//                        count++;
//                   }
//                }
//            if (count!=sz   )
//                return false;
//		}
//	}
//    return true;
}

bool MOLS::checkRCS( permutation &P){
	if (k<2) return true;
	/*if (partMOLS[1].size()<1) // ie if this is where we insert u_0^(0)
		return true;//return isCSR(P);*/

	//int targetCS[n+1]={0};
	vector<int> targetCS(n+1);
	//map<int, vector<vector<int> > > dummy_cycles_map;
	getCycleStructure(partMOLS[1].front(), targetCS  );

	int i=0;
	for (i=0; i<k; i++){
		lsquare::iterator permIt1;
		for (permIt1=partMOLS[i].begin(); permIt1!= partMOLS[i].end(); ++permIt1){
			permutation tempPerm = rcs(*permIt1, P);
			int test =compareCS(tempPerm,targetCS);
			//cout<<test<<endl;

			if (test<0 ){
				//cout << "Rel cycle structure of ";printPerm(*permIt1);cout<< " and ";printPerm(P);cout<< "is "; printPerm(tempPerm); cout <<"Smaller";
				return false;
			}
		}
	}
	return true;
}



void MOLS::printDots(int cursquare, int size, int k){
	int p = size*k+cursquare;
	int i=0;
	for (i=0;i<p; i++)
		cout<<".";
}
void MOLS::printDots(int num){
	int i=0;
	for (i=0;i<num; i++)
		cout<<".";
}

void MOLS::buildPossibleSquare( vector<vector<int> > &square ){
	int r,c;
	//square = {{1}};

		for (c=0;c<n;c++){
			square[0][c]=0;
		}


	square[0][partMOLS[currSquare].size()]=1;
	list<permutation>::iterator permIt;
	permutation::iterator pit;

	for (unsigned int i=0;i< partMOLS[currSquare].size(); i++){
			for (unsigned int j=0; j<partMOLS[currSquare][i].size(); j++){

				square[j][partMOLS[currSquare][i][j]]=0;
			}

		}



/*	int q,w;
	for (q=0;q<n;q++)	{
			for (w=0;w<n;w++){
				cout<< square[q][w]<<" ";
			}
			cout<<endl;
		}*/
}


void MOLS::generatePossiblePermsToInsertRec(vector<vector<int> > &square, permutation &p, int row,  list<permutation> &possiblePermList ){

	if (row==n){
		possiblePermList.push_back((p));
		//printPerm(p);
		return;
	}

	int i,j;
	for (i=0; i<n; i++){
		if (square[row][i] ==1){
			p.push_back(i);
			for(j=0; j<n; j++){
				if (square[j][i]==1)
					square[j][i]=-row;
			}

			generatePossiblePermsToInsertRec(square, p, row+1,   possiblePermList	);
			p.pop_back();
			for(j=0; j<n; j++){
				if (square[j][i]==-row)
					square[j][i]=1;
			}
		}
	}
	return;
}

void MOLS::generatePossiblePermsLSRec(permutation& p, int row,  list<permutation> &possiblePermList, vector<bool > &left ){
	if (row==n){
		possiblePermList.push_back(p);
		//printPerm(p);
		return;
	}

	//cout<< endl<< row<<endl; printPerm(p);cout<<endl;
	//int q,w;
		/*for (q=0;q<n;q++)	{
			for (w=0;w<n;w++){
				cout<< currLS[currSquare][q][w]<<" ";
			}
			cout<<endl;
		}*/

	//int i,j;
	for (unsigned int i=0; i<n; i++){
		if ((left[i]) && (currLS[currSquare][row][i]==n)){
			left[i] = false;
			p[row] = i;
			/*for(j=row+1; j<n; j++){
				if (currLS[currSquare][j][i]==n)
					currLS[currSquare][j][i]= (-row-1);
			}*/

			generatePossiblePermsLSRec( p, row+1,   possiblePermList, left	);
			left[i] = true;
			//bp.pop_back();
			/*for(j=row+1; j<n; j++){
				if (currLS[currSquare][j][i]== (-row-1) )
					currLS[currSquare][j][i]=n;
			}*/
		}
	}
	return;
}

void MOLS::generatePossiblePermsLSRec2(permutation& p, int row,  list<permutation> &possiblePermList, vector<vector<int> > &poss, vector<bool > &left ){
	if (row==n){
		possiblePermList.push_back(p);
		//printPerm(p);
		return;
	}
	//printPerm(p);

	/*if (p[row]>=0)
		generatePossiblePermsLSRec2( p, row+1,   possiblePermList, poss, left	);
	else{*/
		int j = poss[row].size();
		for (unsigned int i=0; i<j; i++){
			if (left[poss[row][i]]  ){
				p[row] = poss[row][i];

				left[poss[row][i]] = false;
				generatePossiblePermsLSRec2( p, row+1,   possiblePermList, poss, left	);
				left[poss[row][i]] = true;
			}
		}
	//}
	return;
}

void MOLS::generatePossiblePermsLSRec3(permutation& p, int row,  vector<permutation> &possiblePermList, vector<vector<int> > &poss, vector<bool > &left ){
 	if (row==n){
		possiblePermList.push_back(p);
		//printPerm(p);
		return;
	}
	//printPerm(p);

	/*if (p[row]>=0)
		generatePossiblePermsLSRec2( p, row+1,   possiblePermList, poss, left	);
	else{*/
		int j = poss[row].size();
		for (unsigned int i=0; i<j; i++){
			if (left[poss[row][i]]  ){
				p[row] = poss[row][i];

				left[poss[row][i]] = false;
				generatePossiblePermsLSRec3( p, row+1,   possiblePermList, poss, left	);
				left[poss[row][i]] = true;
			}
		}
	//}
	return;
}

vector<permutation> MOLS::generatePossiblePermsToInsert3(int square){
//void MOLS::generatePossiblePermsToInsert3(int square){
	//int cs = partMOLS[k-1].size();

	permutation p(n);
	int j=partMOLS[currSquare].size();
	vector<permutation> currList;
	p[0] =j; // because i in the first row is always in the ith position
	vector<int> evec;
	vector<vector<int> > poss(n);
	vector<bool> left(n, true);
	left[j] = false;

	for( unsigned int i=1; i<n; i++){
		for(unsigned int j=0; j<n; j++){
			if (currLS[square][i][j] == n)
				poss[i].push_back(j);
		}
		/*if (poss[i].size()==1){
			p[i] = poss[i].front();
			left[p[i]] = false;
		}*/
	}


	generatePossiblePermsLSRec3( p, 1,   currList, poss, left	) ;
	//sqSymPossPerms[square][cs] = currList;
  	return currList ;
}

void MOLS::generatePossiblePermsToInsert2( list<permutation> &possiblePermList ){
	permutation p(n);
	int j=partMOLS[currSquare].size();

	p[0] =j; // because i in the first row is always in the ith position
	vector<int> evec;
	vector<vector<int> > poss(n);
	vector<bool> left(n, true);
	left[j] = false;

	for( unsigned int i=1; i<n; i++){
		for(unsigned int j=0; j<n; j++){
			if (currLS[currSquare][i][j] == n)
				poss[i].push_back(j);
		}
		/*if (poss[i].size()==1){
			p[i] = poss[i].front();
			left[p[i]] = false;
		}*/
	}

	generatePossiblePermsLSRec2( p, 1,   possiblePermList, poss, left	) ;


	//cout<< "Size"<<possiblePermList.size()<<endl;

	return;
}


void MOLS::generatePossiblePermsToInsert( list<permutation> &possiblePermList ){
	int i; permutation p(n);
	int j=partMOLS[currSquare].size();
	p[0] = j; // because i in the first row is always in the ith position
	vector<bool> left(n, true);
	left[j] = false;
	/*for(i=1; i<n; i++){
		if (currLS[currSquare][i][j]==n)
			currLS[currSquare][i][j]=-1;
	}*/
	generatePossiblePermsLSRec( p, 1,   possiblePermList, left	) ;
	/*for(i=1; i<n; i++){
		if (currLS[currSquare][i][j]==-1)
			currLS[currSquare][i][j]=n;
	}*/
	//cout<< "Size"<<possiblePermList.size()<<endl;

	return;
}


void addToCompletedMOLS(){

}

void MOLS::addToLS(){
	//cout<<"enter add LS"; cout.flush();
	//list<permutation>::iterator permIt;

	//permutation::const_iterator pit;
	//cout<<partMOLS[currSquare].size(); cout.flush();
	//printPerm(partMOLS[currSquare].back()); cout.flush();

	//int r=0;

	//for (pit=partMOLS[currSquare].back().begin(); pit != partMOLS[currSquare].back().end(); ++pit){
	//for (unsigned int i=0; i<n; i++)
	//	currLS[currSquare][i][partMOLS[currSquare][i]] = partMOLS[currSquare].size()-1;

	//cout<<"exit add LS"; cout.flush();
}

void MOLS::removeFromLS(permutation p){

	for (unsigned int i=0; i<n; i++)
		currLS[currSquare][i][p[i]] = n;

}

void MOLS::addUniversal(permutation &p){
	//cout<<"enter add"; cout.flush();
	//printPerm(p);
	//cout<< currSquare; cout.flush();

	partMOLS[currSquare].push_back(p);

	for (unsigned int i=0; i<n; i++)
			currLS[currSquare][i][p[i]] = partMOLS[currSquare].size()-1;

	//cout<< partMOLS[currSquare].size(); cout.flush();
	//addToLS();
	//cout<<"universal added"; cout.flush();


}

void MOLS::removeUniversal(permutation &p){
	//removeFromLS(p);
	for (unsigned int i=0; i<n; i++)
	currLS[currSquare][i][p[i]] = n;
	partMOLS[currSquare].pop_back();
}

void MOLS::updatePossiblePerms(){
	int csymb = partMOLS[k-1].size();

	for (unsigned int i=0; i<k; i++){
		sqSymPossPerms[i][csymb] = (generatePossiblePermsToInsert3(i));
		//sqSymCurrUni[i][csymb][0] = 0;
		sqSymCurrUni[i][csymb][1] = sqSymPossPerms[i][csymb].size();

		if (sqSymCurrUni[i][csymb][0] < 0  ) //if i have no position in the current list then i also have no
			sqSymCurrUni[i][csymb][2] = sqSymPossPerms[i][csymb].size()+1; //limit one more, <
	}


	/*for (unsigned int i=0; i<k; i++)
		generatePossiblePermsToInsert3(i);*/

 }

double MOLS::calculate_fraction(){
	double fd = (1.*nof_calls_made) / nof_calls_limit;

	//double prev =1;
	/*for (int j=0; j<3; j++){
		for (int i=0; i<k; i++){
			if (sqSymCurrUni[i][j][0] >0){
				prev = 1./sqSymCurrUni[i][j][1]*prev;
				fd = fd + sqSymCurrUni[i][j][0]*prev;
			}
		}
	}*/
	return fd;
}

void MOLS::findMOLS4(){
	int counter = (partMOLS[currSquare].size())*k+currSquare-1;
	//cout<< sqSymCurrUni[currSquare][partMOLS[currSquare].size()][0]<<" "<< sqSymCurrUni[currSquare][partMOLS[currSquare].size()][2]<<" "<< sqSymCurrUni[currSquare][partMOLS[currSquare].size()][1]<<endl;
 	//branchCount_[(partMOLS[currSquare].size())*k+currSquare-1]++;
 	// cout<<".";
	int i=0;

	if (partMOLS[k-1].size()==n) //if the last square is filled in completely, contains n permutations
		{
			if (isSmallestConjugateMOLS(partMOLS)){
				//if (isSmallest2<k>(partMOLS)){
				branchCount_[(partMOLS[currSquare].size())*k+currSquare-1]++;

				printMOLS(partMOLS);
				count_MOLS++;
				//addToCompletedMOLS();
				return;

				//}
			}
			else	return;

		}
		else
			branchCount_[(partMOLS[currSquare].size())*k+currSquare-1]++;

	/*printDots(2*partMOLS[0].size());
	cout<<branchCount_[partMOLS[currSquare].size()*k+currSquare-1];
	printDots(4);

	for (i=0; i<k; i++){
		if (partMOLS[i].size()>0){
			printPerm(partMOLS[i].back()	);
			cout<<" ";}
	}
	cout<<endl;*/

	/*if (currSquare==0&& partMOLS[k-1].size()==2)
		return;*/
	//else
		//detailedCount[counter].push_back(0);
	if (currSquare==0){
		if (partMOLS[0].size()>0){
			/*printDots(2*partMOLS[0].size());
			cout<<branchCount_[partMOLS[currSquare].size()*k+currSquare-1];
			printDots(4);

			for (i=0; i<k; i++){
				printPerm(partMOLS[i].back()	);
				cout<<" ";
			}
			cout<<endl;*/
			updatePossiblePerms();
  		}
	}
/*
	if (currSquare==0&& partMOLS[k-1].size()==2)
		return;*/
	//if this is the first universal iln the second square we look at the class representatives

		//}
		if (partMOLS[currSquare].size()>0)
		{//in general
 			int currUni = partMOLS[currSquare].size();

			for (unsigned int possPermIt=max(0,sqSymCurrUni[currSquare][currUni][0]); possPermIt< sqSymPossPerms[currSquare][currUni].size() && send_checkpoint<3; ++possPermIt){
				//sqSymCurrUni[currSquare][currUni][0] = possPermIt; //keep track of current one for progress and checkpoint
				sqSymCurrUni[currSquare][currUni][0] = possPermIt; //keep track of current one for progress and checkpoint
 				//if (boinc_time_to_checkpoint()  || nof_calls_made > nof_calls_limit) {
				if ((possPermIt+1)%600==0 || nof_calls_made > nof_calls_limit	|| possPermIt == sqSymCurrUni[currSquare][currUni][2]) { //doing it before so that I have the 'prior' branchcounts
					//if (partMOLS[0].size()<4){
						//double fd = calculate_fraction();
						//cout<< fd << endl;
						//if (cpu_time) fd /= 2;
						boinc_fraction_done(calculate_fraction());
					//}
					//cout<<" going into checkpoint sendval is="<<send_checkpoint<<endl;
					if (possPermIt == sqSymCurrUni[currSquare][currUni][2]) {
						send_checkpoint = 2;
						cout<< "up against limit square: uni, curruni, size[1], limit[2], size() "<< currSquare<<" "<<currUni<<" "<<possPermIt<<" "<<sqSymCurrUni[currSquare][currUni][1]<<" "<<sqSymCurrUni[currSquare][currUni][2]<<", "<<sqSymPossPerms[currSquare][currUni].size()<<endl;
					}
					int retval = do_checkpoint(k,n, sqSymCurrUni, branchCount_, (clock()-dtime) + prev_time, count_MOLS,root, send_checkpoint, nof_calls_made);
					//cout<<"Coming backsendval is="<<send_checkpoint<<endl;
/*					if (retval) {
				 						fprintf(stderr, "%s APP: upper_case checkpoint failed %d\n",
				 								boinc_msg_prefix(buf, sizeof(buf)), retval
				 						);
				 						exit(retval);
				 					}*/
					if (send_checkpoint>0) send_checkpoint =3;
					boinc_checkpoint_completed();
				}


				if (checkOrthogonal(sqSymPossPerms[currSquare][currUni][possPermIt])){
					//printDots(currSquare, partMOLS[currSquare].size(), k);  printPerm((*possPermIt));
					//if (checkRCS( (*possPermIt))){
 						addUniversal(sqSymPossPerms[currSquare][currUni][possPermIt]);
 						if (getSmallRelCS(partMOLS)){
 							if (isSmallest4()){
 								currSquare = (currSquare+1)%k;
 								//detailedCount[counter].back()++;
								findMOLS4( );
								currSquare = (currSquare-1+k)%k;
							}
						}
						removeUniversal(sqSymPossPerms[currSquare][currUni][possPermIt]);
 					//}
				/*	else{

												for (i=0; i<currSquare; i++){
													printPerm(partMOLS[i].back());cout<<"\t";}
												printPerm((*possPermIt));
												cout<< "checkrcs "; cout<<endl;
											}*/
				}
				/*else{

					for (i=0; i<currSquare; i++){
						printPerm(partMOLS[i].back());cout<<"\t";}
					printPerm((*possPermIt));
					cout<< "orthog "; cout<<endl;
				}*/
			}
			for (int ii=currSquare; ii<k; ii++){
				sqSymCurrUni[ii][currUni][0] = -1;
				//cout<<"clear from "<<currSquare<<endl;

			}
			//sqSymCurrUni[currSquare][currUni][1] = -1;


		}
    else {
        if(currSquare==1&& partMOLS[currSquare].size()==0	){
        vector<permutation>::iterator CSRit;
		int j=1;
		for (CSRit=cycleStructureReps.begin(); CSRit!= cycleStructureReps.end(); ++CSRit){
			//++CSRit; ++CSRit;++CSRit;
			if (printOut){
				cout << j++ << ". u0^1 = ";
				printPerm(*CSRit);
				cout<<endl;

			}
			//No need for any of the checks, as this is only the first perm in the second square.
			//Automatically orthog, will fit (empty), forced to be RCS..

			//map<int, vector<vector<int> > > dummy_cycles_map;
			getCycleStructure(*CSRit,  currentCS);
			addUniversal(*CSRit);
			/*
			cout<<"Added cs"<<endl;
			printMOLS(partMOLS);*/
			possibleShuffles = genRelevantPermutations(partMOLS[currSquare].front());
			currSquare = (currSquare+1)%k;
			//detailedCount[counter].back()++;
			findMOLS4( );
			currSquare = (currSquare-1+k)%k;
			/*cout<<"Return from findmols structure"<<endl;
			printMOLS(partMOLS);*/
			removeUniversal(*CSRit);

			//CSRit = --cycleStructureReps.end();
			/*cout<<"cs removed"<<endl;
			printMOLS(partMOLS);*/
			currentCS.clear();currentCS.resize(n+1, 0);
			//partMOLS[currSquare].pop_back();
		}
	}
	else{ //not u_0^(1) but still  0universal
		permutation P(identity) ;
		if (partMOLS[currSquare].size()==0 ){
		//	if (!testPerms.size()==0){
				//cout << " sq "<< (currSquare+1)<< ", "<< testPerms.size()<<endl;
				do{
					if (P.front()<1){
						if (checkFit(P)){
							if (checkOrthogonal(P)){
									addUniversal(P);
									if (getSmallRelCS(partMOLS)){
										if (isSmallest4()){
											currSquare = (currSquare+1)%k;
											//detailedCount[counter].back()++;
											findMOLS4( );
											currSquare = (currSquare-1+k)%k;
										}
									}
									removeUniversal(P);

								/*else{

									for (i=0; i<currSquare; i++){
										printPerm(partMOLS[i].back());cout<<" ";}
									printPerm(P);
									cout<< "checkrcs "; cout<<endl;
								}*/
							}
							/*else{

								for (i=0; i<currSquare; i++){
									printPerm(partMOLS[i].back());cout<<" ";}
								printPerm(P);
								cout<< "orthog "; cout<<endl;
							}*/
						}
						/*else{

							for (i=0; i<currSquare; i++){
								printPerm(partMOLS[i].back());cout<<" ";}
							printPerm(P);
							cout<< "checkfit "; cout<<endl;
						}*/
					}
					else break;
				}while(next_permutation(P.begin(), P.end()));
			}
}
	}
	return ;
}

void MOLS::printAllStatics(){
	cout<<"***************Statics*****************"<<endl	;
	cout<< "n, k = "<< n<<", "<<k<<endl;
	cout<<"Currsquare = "<<currSquare<<endl;
	cout<<"CurrCS = ";
	for (unsigned int i = 0; i<n; i++)
		cout<< currentCS[i]<<" ";
	cout<<endl;
	cout<<"partMOLS = ";
	printMOLSPerms(partMOLS);
	cout<<"tempMOLS = ";
	printMOLSPerms(tempMOLS);
	cout<<"possibleShuffles "<<possibleShuffles.size()<<"= ";
	vector<permutation>::const_iterator it;
	for (it = possibleShuffles.begin(); it != possibleShuffles.end(); it++)
		printPerm(*it);
	cout<< "sqSymbCurrUni \n";
	for (int i=0; i<n; i++){
		for (int j=0; j<k; j++){
			cout << sqSymCurrUni[j][i][0]<<" "<<sqSymCurrUni[j][i][2] <<" "<<sqSymCurrUni[j][i][1] <<" ";
		}
		cout <<endl;
	}
	cout<<nof_calls_made<<endl;
	for( int i=0; i <n*k; i++)
		cout << branchCount_[i]<<" ";
	cout <<endl;
	cout<<endl<<"************DONE***************"<<endl;


}

//enumerate says whether it is called from a checkpoint restart, if it is we are interested in the branchcounts and nofcalls made
int MOLS::restart(FILE* state, bool enumerate){
	cout<<"into restart: "<<nof_calls_made<<endl;
	//this->out = out;
	char sc[100];
	string s;
	//skipthe first few lines, all found MOLS onlyinterested from the positions onwards
	fgets(sc,100,state);
	s= (string) sc;
	cout<<"YY"<<s<<endl;
	unsigned int i=0;
	while(s[0] !='#' )
	{
		if (feof(state))	{
			break;
		}
		fgets(sc,100,state);
		s= (string) sc;
	} //so the current line after the while is "#Positions"
		
	for (int i =0; i<1;i++){

			for (int j=0; j<n; j++){
				fgets(sc,100,state);
				if (feof(state)) break;
				s = (string) sc;
				cout<<s<<endl;
				vector<string> vec = split(s,' ');
				for (int jj =0; jj<k; jj++){
					sqSymCurrUni[jj][j][0] = atoi(vec[3*jj].c_str());
					sqSymCurrUni[jj][j][2] = atoi(vec[3*jj+1].c_str());
					if (sqSymCurrUni[jj][j][0] >= 0 ) {
						sqSymCurrUni[jj][j][1] = atoi(vec[3*jj+2].c_str());
						sqSymCurrUni[jj][j][2] = atoi(vec[3*jj+1].c_str());
					}
					else
						sqSymCurrUni[jj][j][2] = 10000000;
				}

			}// done reading positinos, assume we are not interested in branchcounts
			//cout<<"done reading positions"<<endl;
			for (int j=0; j< k*n; j++){
				if (sqSymCurrUni[j%k][j/k][0]>=0 &&sqSymCurrUni[(j+1)%k][(j+1)/k][0]>=0){
					//cout<<"decreasing branches"<<endl;
					branchCount_[j] =  -1; //count every piece seperately, sum from resultfiles.
					nof_calls_made--;
				}
				else
					branchCount_[j] = 0;
				//branchCount_[j] = atol(vec[j].c_str());
			}

			fgets(sc,100,state); //number of calls ie 454769976
			if (feof(state)) break;
			s = (string) sc;
			if (enumerate) nof_calls_made = atol(s.c_str()); //not interestedin thiswhen restart iscalled form MOLS() with infile, only when called wiuth checkpoint
			// when called with a checkpoint sent a s a starting position it breaks afer this, only up to callsmade is sent.
			fgets(sc,100,state); //#Branchcounts
			if (feof(state)) break;
			//s = (string) sc;
 			fgets(sc,100,state); // actual branches
			s = (string) sc;
			cout<< s<<endl;
			vector<string> vec = split(s,' ');
			cout <<vec.size()<<vec[vec.size()-1]<< endl;
			if (enumerate){
				cout<< "I am reading branches from vecmols_state"<<endl;
				for (int j=0; j< k*n; j++){
					if (sqSymCurrUni[j%k][j/k][0]>=0 &&sqSymCurrUni[(j+1)%k][(j+1)/k][0]>=0){
						branchCount_[j] = atol(vec[j].c_str())-1;
					}
					else
						branchCount_[j] = atol(vec[j].c_str());
				}
			}

			fgets(sc,100,state); // # time
			fgets(sc,100,state); // actual time
			if (feof(state)) break;
			s = (string) sc;
			vec = split(s,' ');
			cout<< s;
			if (enumerate) prev_time=(long) (1000000*atof(vec[0].c_str()));

			fgets(sc,100,state); // # number of mols
			fgets(sc,100,state); // actual number
			if (feof(state)) break;
			s = (string) sc;

			cout<< s;
			if (enumerate) count_MOLS=atoi(s.c_str());

			// out.printf(s.c_str());
			// out.printf( "restarting \n");

		}
		//out.flush();
	cout<< "done restartingX: "<<nof_calls_made<<endl ;
	printAllStatics();
	
	//return count_MOLS;
	//return 0;

	if (enumerate)
		return enumerateMOLS();
	else 
		return 0;		
}


/*int main(int argc,char *argv[]){

	MOLS threemols(7,3);
	string filename = argv[1];
   //MOLS threemols(filename);
/*
	string outfile = "out103_1.txt";//+filename;
	std::ofstream out(outfile.c_str());
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf());
*/




	//threemols.enumerateMOLS();
	// std::cout.rdbuf(coutbuf); //reset to standard output again

//}

const char *BOINC_RCSID_33ac47a071 = "$Id: upper_case.cpp 20315 2010-01-29 15:50:47Z davea $";

