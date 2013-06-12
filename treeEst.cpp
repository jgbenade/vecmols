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
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include "../bigintlib/BigIntegerLibrary.hh"

using namespace std;

typedef vector<int> permutation ;
typedef vector<permutation> lsquare;
int main(int argc,char *argv[]);

class MOLS {
public:
	typedef vector<int> permutation ;
	typedef vector<permutation> lsquare;
	MOLS(int n, int k) ;
	MOLS( std::string filename);
	~MOLS();
	int enumerateMOLS(void);
	void  printAllStatics();
	struct My
	{
	    static bool permutationComparator (permutation p1, permutation    p2);
	};

private:
	int n;
	int k;
	vector<lsquare> partMOLS;
	vector<lsquare> tempMOLS;
	vector<vector<vector<int> > > currLS;
	list<permutation> cycleStructureReps;
	//list<vector <lsquare> > completedMOLS;

	int currSquare;
	int count_MOLS;
	vector<int> branchCount_  ;
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
	vector<vector<int> > detailedCount;
	vector<BigInteger> numIsSmallest;
	vector<BigInteger>  numIsSmallestTrue;
	vector<BigInteger> estimates;
	vector<vector<BigInteger> > generalSummary;
	int positionsFound;

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
	list<permutation>  getShuffles(permutation &pOrig, permutation &pNow) ;
	bool  noSmallerRCS(permutation &smallestRCS ) ;
	bool  noSmallerRCS(permutation &smallestRCS, vector<lsquare> &pMOLS ) ;
	bool  isSmallest4() ;
	bool  isSmallest4(bool t) ;
	bool  isSmallestConjugateMOLS(vector<lsquare> partMOLS) ;
	bool  checkFit( permutation &p) ;
	bool  checkOrthogonal(permutation &P) ;
	bool  checkRCS(permutation &P) ;
	void  printDots(int cursquare, int size, int k) ;
	void  printDots(int num) ;
	void  buildPossibleSquare( vector<vector<int> > &square ) ;
	void  generatePossiblePermsToInsertRec(vector<vector<int> > &square, permutation &p, int row,  list<permutation> &possiblePermList ) ;
	void  generatePossiblePermsToInsert(  list<permutation> &possiblePermList ) ;
	int  findMOLS4(float scaleFactor, bool Visit) ;
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
} ;

inline bool MOLS::My::permutationComparator    (permutation p1, permutation    p2) { return p1[0] < p2[0]; }

inline bool MOLS::permutationComp    (permutation p1, permutation    p2) { return p1[0] < p2[0]; }

MOLS::~MOLS(){}

MOLS::MOLS( int n, int k){
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
	detailedCount.resize(n*k+1);
	numIsSmallest.resize(n ) ;
	numIsSmallestTrue.resize(n ) ;
	estimates.resize(n*k );
	generalSummary.resize(n*k, vector<BigInteger>(2 ));
	srand(time(0))	;
	positionsFound=0;

    for (unsigned int i=0; i<n; i++){
		identity.push_back(i);
	}

	lsquare emptyLS;

	int i,j, ii ;
	for (i =0; i<k; i++){
 		partMOLS.push_back(emptyLS);
		tempMOLS.push_back(emptyLS);
	}
	//Get cycle structure representatives for U_0^1
	vector<int> cycles(n+1, 0);
	cycles[1] =1;
	genCRS(cycles, 2, n-1);
	list<permutation>::const_iterator CSRit;

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
	int i;
	permutation p;
	for(i=0; i<s.length(); i++){
		p.push_back(  s[i]-'0');
	}
	return p;
}

MOLS::MOLS( string filename){

	//this->filename = filename;
	printOut = true;
	/*this->n = 5;
	this->k = 3;
	currSquare = 0;
	count_MOLS = 0;
	branchCount_[100] = {0} ;*/

	int i;
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
			 cout<<   k << " MOLS of order "<< n<<endl;
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
		 cout<< "File not found"<<endl;

	printMOLSPerms(partMOLS);
	cout<< n<< k;

	//possibleShuffles = genRelevantPermutations(partMOLS[1][0]);
	positionsFound=0;
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
	detailedCount.resize(n*k+1);
	numIsSmallest.resize(n) ;
	numIsSmallestTrue.resize(n ) ;
	estimates.resize(n*k );
	generalSummary.resize(n*k, vector<BigInteger>(2 ));
	srand(time(NULL))	;



	//Get cycle structure representatives for U_0^1
	vector<int> cycles(n+1, 0);
	cycles[1] =1;
	genCRS(cycles, 2, n-1);
	list<permutation>::const_iterator CSRit;
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

	//Initialise symbSqRC
//	symbSqRC.resize(n, vector<vector<vector<int> > >(k, vector<vector<int> >(n, vector<int>(n, 1))));
	//completedMOLS.clear();

	return;
}

void MOLS::updateLS(){
	currLS.resize(k, vector<vector<int> > (n , vector<int>(n, n)));
	int   i;
	for (i=0; i<k; i++){
 		//		for(universalIt=partMOLS[i].begin(); universalIt != partMOLS[i].end(); ++universalIt){
		for(unsigned int j =0; j< partMOLS[i].size(); j++){
			for(unsigned int jj= 0; jj< partMOLS[i][j].size(); jj++){//(partMOLS[i][j]).begin(); perm_it != (partMOLS[i][j]).end(); ++perm_it){
				currLS[i][jj][partMOLS[i][j][jj]] = j;
			}

		}
	}

}

int MOLS::enumerateMOLS(void	){
/*	string outfile = filename+ ".out.txt";
	std::ofstream out(outfile);
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf());*/

	clock_t t;
	t = clock();
	cout<<   k << " MOLS of order "<< n<<endl;

	int i;
	count_MOLS =0;

	if (partMOLS[0].size()==0){
		permutation ss;
		for (i=0;i<n; i++){
			ss.push_back(i);
		}
		addUniversal((ss));
		currSquare++;
		findMOLS4(1, true);
	}
	else{
		if (partMOLS[1].size()>0){
			map<int, vector<vector<int> > > dummy_cycles_map;
			getCycleStructure(partMOLS[1].front(),  currentCS, dummy_cycles_map );
			possibleShuffles = genRelevantPermutations(partMOLS[1].front());
		}
		//if (isSmallest4())
		printDots(2);
		cout<<1;
		printDots(4);

		for (i=0; i<k; i++){
			printPerm(partMOLS[i].back()	);
			cout<<" ";
		}
		cout<<endl;
		findMOLS4(1, true);
	}


	/*vector<int>::const_iterator vit;
	for (i=0; i< k*n; i++){
		cout <<"@ ";
		for (vit = detailedCount[i].begin(); vit!= detailedCount[i].end(); vit++)
			cout<< *vit<< " ";
		cout<< " opsies na "<<i/k<<"."<<i%k<<endl;
	}*/

	cout<< "# "<<count_MOLS << " MOLS found"<<endl;
    cout<< "# Visited ";
	for (i=0; i<n*k;i++)
		cout<<branchCount_[i]<<" ";
	cout<<endl;

	cout<< "# Estimated ";
	for (i=0; i<n*k;i++)
			cout<<estimates[i]<<" ";
		cout<<endl;

		cout<< "# Summ ";
			for (i=0; i<n*k;i++)
					cout<<bigIntegerToString(generalSummary[i][0])<<"/"<<bigIntegerToString(generalSummary[i][1])<<" ";
				cout<<endl;

	/*for (i=0; i< n; i++){
		cout<<"# "<< branchCount_[k-1 +i*k] <<" takke op vlak "<<i<< " | "<< numIsSmallest[i]<<" "<<numIsSmallestTrue[i]<< " "<< (numIsSmallestTrue[i]*1.)/numIsSmallest[i]<<endl;
	}*/


	t = clock() - t;
	cout<<"# "<<  t/1000000.0 << " seconds"<< endl;

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

void MOLS::printMOLSPerms( vector<lsquare> mols  ){
	//int k = mols;
	int m =0;
	cout << endl;
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
	int m;stringstream outp[k][n];
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
		int i=0, j=0;
		for ( i =0; i<n; i++){
			for ( j =0; j<k; j++){
				cout<<outp[j][i].str()<<"\t";
			}
			cout<<endl;
		}
		cout <<"*******************************"<<endl;

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
	//int n = p.size();

	//int thisCS[n+1] ={0};
	vector<int> thisCS(n+1, 0);
	map<int, vector<vector<int> > > dummy_cycles_map;
	getCycleStructure(p, thisCS ,dummy_cycles_map );

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
	numIsSmallest[partMOLS[k-1].size()]+=1;

	if (!noSmallerRCS( partMOLS[1].front() )){
		/*for (i=0; i<k; i++)
				printPerm(partMOLS[i].back());
		cout << "smaller cs"<<endl;*/
			return false;
		}
	numIsSmallestTrue[partMOLS[k-1].size()]=numIsSmallestTrue[partMOLS[k-1].size()]+1;
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


bool MOLS::isSmallestConjugateMOLS(vector<lsquare> partMOLS ){
	return true;

 	int i;
 	vector<lsquare> smallestMOLS ;
	for (i=0; i<k; i++){
		smallestMOLS[i] = partMOLS[i];
	}

	list<permutation> allCSR; //list of all smalle s u0^1 cycle structures
	//int cycles[n+1];
	vector<int> cycles(n+1, 0);
	/*for (i=0; i<=n+1; i++)
		cycles[i]=0;*/
	cycles[1] =1;
	genCRS(cycles, 2, n-1); //we may only have one overlap

	list<permutation>::iterator CSRit;

	permutation smallestRCS = partMOLS[1].front();
	vector<lsquare> t1MOLS , newMOLS;



    vector<vector<int> > oa (k+2, vector<int>(n*n)	);
    vector<vector<int> > oanew(k+2, vector<int>(n*n));

   	MOLStoOA(partMOLS, oa);
	//printOA<k,n>(oa);
	//int rowls[] = {0,1,2,3,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
	permutation p[3] ;//(rowls, rowls + n / sizeof(int) );

	for (i=0;i<k+2; i++){
			p[0].push_back(i);
			p[1].push_back(i);
			p[2].push_back(i);
	}


	p[0][0]=2; p[0][2] = 0; //(i,j,k) - (k,, j, i...)
	p[1][1]=2; p[1][2] = 1; //(i,j,k) - (k,, j, i...)
	p[2][0]=2; p[2][1] = 3;p[2][2]=0; p[2][3] = 1; //(i,j,k) - (k,l, i, j...)

	int parr[k+2];
	int perm_ctr = 0;

	while(perm_ctr<3 ){
		i=0;

		permutation::iterator permIt;
		for (permIt = p[perm_ctr].begin(); permIt !=  p[perm_ctr].end(); ++permIt){
			parr[i++] =  *permIt;
		}

		permuteOA(oa, parr, oanew);
		//cout<< "  -  oanew "; cout.flush();
		//printOA<k,n>(oanew);
		OAtoMOLS(oanew, t1MOLS);
		//cout<< "  -  t1MOLS "; cout.flush();
		//OAtoMOLS<k,n>(oanew, (flag1? t2MOLS : t1MOLS));
		//printMOLS(t1MOLS, k);

		//getSmallest<k,n>(t1MOLS, newMOLS);
		//newmols comes bac empty
		//list<permutation> testPs =
		//if( !noSmallerRCS<k,n>(smallestRCS, partMOLS, partMOLS){
		for (CSRit=allCSR.begin(); CSRit!= allCSR.end(); ++CSRit){
			//printPerm(*CSRit); cout.flush();
			if (!noSmallerRCS((*CSRit), t1MOLS )){
				/*for (i=0; i<k; i++)
		    				printPerm(partMOLS[i].back());
		    		cout << "smaller cs"<<endl;*/
				cout <<"Not smallest conjugate";
				return false;
			}

		}
		int x;//if this wasnt smalleer, wip the mols so that a new one may be tried
		for (x=0; x<k; x++)
			t1MOLS[x].clear();

		perm_ctr++;
	}

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
	map<int, vector<vector<int> > > dummy_cycles_map;
	getCycleStructure(partMOLS[1].front(), targetCS,dummy_cycles_map );

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
	for (unsigned int i=0; i<k; i++)
		sqSymPossPerms[i][csymb] = (generatePossiblePermsToInsert3(i));

	/*for (unsigned int i=0; i<k; i++)
		generatePossiblePermsToInsert3(i);*/

 }

 int MOLS::findMOLS4(float scaleFactor, bool Visit){
	int counter = (partMOLS[currSquare].size())*k+currSquare-1;
 	branchCount_[counter]++;
 	if (positionsFound>100)
 		return 0;

	int i=0;
	int nofFeasible;

	if (partMOLS[k-1].size()==n) //if the last square is filled in completely, contains n permutations
		{
			if (isSmallestConjugateMOLS(partMOLS)){
				//if (isSmallest2<k>(partMOLS)){
					printMOLS(partMOLS);
					count_MOLS++;
					//addToCompletedMOLS();
					return 0;

				//}
			}
			else return 1;

		}
	/*printDots(2*partMOLS[0].size());
	cout<<branchCount_[partMOLS[currSquare].size()*k+currSquare-1];
	printDots(4);

	for (i=0; i<k; i++){
		if (partMOLS[i].size()>0){
			printPerm(partMOLS[i].back()	);
			cout<<" ";}
	}
	cout<<endl;*/

	/*if (currSquare==1&& partMOLS[k-1].size()==2)
		return 0;
	else*/
		detailedCount[counter].push_back(0);
		detailedCount[counter].push_back(0);

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
		int totalBranches = sqSymPossPerms[currSquare][currUni].size();
		int nofFeasible = 0;
		//count numfeasible completions and store indices
		vector<int> feasiblePoss;
		feasiblePoss.reserve(totalBranches);
		for (unsigned int possPermIt=0;    possPermIt< sqSymPossPerms[currSquare][currUni].size(); ++possPermIt ){
			if (checkOrthogonal(sqSymPossPerms[currSquare][currUni][possPermIt])){
				//printDots(currSquare, partMOLS[currSquare].size(), k);  printPerm((*possPermIt));
				//if (checkRCS( (*possPermIt))){
				addUniversal(sqSymPossPerms[currSquare][currUni][possPermIt]);
				if (getSmallRelCS(partMOLS)){
					nofFeasible++;
					feasiblePoss.push_back(possPermIt);
				}
				removeUniversal(sqSymPossPerms[currSquare][currUni][possPermIt]);
			}
		}

		if (Visit && nofFeasible>0){

			int numVisited=0;
			int numPolled=0;
			if (nofFeasible>3){
				numVisited = min(int(.05*nofFeasible), 10);
				if (numVisited<3)
					numVisited = min(3, nofFeasible);
				if (numVisited < .05*nofFeasible){
					numPolled = min(5, nofFeasible); //repolling already visited one, especially when there are few
				}
			}
			else
				numPolled = nofFeasible;

			//this if stops the full estimation, use dfor example for getting startrting point on level i  etc
			if (currSquare==k-1 && partMOLS[k-1].size() ==2){
				numPolled = min(numVisited+numPolled, nofFeasible);
				numVisited=0;
				int numtries = 0;

				for (unsigned int possPermIt=0;  numtries <1 && possPermIt<nofFeasible ;
						possPermIt= rand()%nofFeasible){
					numtries++;
					addUniversal(sqSymPossPerms[currSquare][currUni][feasiblePoss[possPermIt]]);
					//if (getSmallRelCS(partMOLS)){
					//cout<<" inc numdone";
					if (isSmallest4()){
						positionsFound++;
						cout<<positionsFound;cout<<" ";
						for (unsigned int q =0; q<3; q++){

							printDots(2*q);

							//printDots(4);

							for (i=0; i<k; i++){
								printPerm(partMOLS[i][q]	);
								cout<<" ";
							}
							cout<<endl;
						}

					}
					//}
					removeUniversal(sqSymPossPerms[currSquare][currUni][feasiblePoss[possPermIt]]);
				}
			}
			else{

				//int jump = (numVisited>0 ? nofFeasible /numVisited: 0);

				int numtries = 0;
				//now visited some number of them
				BigInteger sum = 0;
				for (unsigned int possPermIt=0;  numtries <numVisited /*&& possPermIt<numdone possPermIt< sqSymPossPerms[currSquare][currUni].size()*/; /*++possPermIt*/
						possPermIt= rand()%nofFeasible){
					numtries++;
					addUniversal(sqSymPossPerms[currSquare][currUni][feasiblePoss[possPermIt]]);
					//if (getSmallRelCS(partMOLS)){
					//cout<<" inc numdone";
					if (isSmallest4(true)){
						currSquare = (currSquare+1)%k;
						//detailedCount[counter].back()++;
						sum = sum+findMOLS4(scaleFactor*(float(nofFeasible)/numVisited), true);
						currSquare = (currSquare-1+k)%k;
					}
					//}
					removeUniversal(sqSymPossPerms[currSquare][currUni][feasiblePoss[possPermIt]]);
				}
				for (unsigned int possPermIt=0;  numtries <numVisited+numPolled /*&& possPermIt<numdone possPermIt< sqSymPossPerms[currSquare][currUni].size()*/; /*++possPermIt*/
						possPermIt= rand()%nofFeasible){
					numtries++;
					addUniversal(sqSymPossPerms[currSquare][currUni][feasiblePoss[possPermIt]]);
					//if (getSmallRelCS(partMOLS)){
					//cout<<" inc numdone";
					if (isSmallest4(true)){
						currSquare = (currSquare+1)%k;
						//detailedCount[counter].back()++;
						sum = sum+findMOLS4(scaleFactor*(float(nofFeasible)/numVisited), false);
						currSquare = (currSquare-1+k)%k;
					}
					//}
					removeUniversal(sqSymPossPerms[currSquare][currUni][feasiblePoss[possPermIt]]);
				}
				BigInteger avgFeasible(sum);
				avgFeasible/=numtries;
				//BigInteger avgF(avgFeasible);
				long int prod = ( long int)(nofFeasible * scaleFactor);
				BigInteger est = avgFeasible*prod ;
				//if (currUni<2&& currSquare<2) cout<<currUni<<"."<< currSquare<<" visiting "<< numVisited<< " / "<<numdone<<" / "<<totalBranches<<" Avg: "<<avgFeasible<< " Est: "<< est <<" , " << scaleFactor<<endl;
				estimates[counter] = estimates[counter] + est;
				//return (avgFeasible *(float(nofFeasible)/numVisited));
				generalSummary[counter][0] = generalSummary[counter][0]+nofFeasible;
				generalSummary[counter][1] = generalSummary[counter][1]+numVisited;
			}
		}
		return nofFeasible;

	}//cout<< currUni<<"."<< currSquare<<"V "<< numdone<< "/"<<numVisited <<"/"<<totalBranches<<endl;


    else {
        if(currSquare==1&& partMOLS[currSquare].size()==0	){
		list<permutation>::iterator CSRit;
		int j=1;
		BigInteger sum ;
		for (CSRit=cycleStructureReps.begin(); CSRit!= cycleStructureReps.end(); ++CSRit){
			//++CSRit; ++CSRit;++CSRit;
			if (printOut){
				cout << j++ << ". u0^1 = ";
				printPerm(*CSRit);
				cout<<endl;
			}
			//No need for any of the checks, as this is only the first perm in the second square.
			//Automatically orthog, will fit (empty), forced to be RCS..

			map<int, vector<vector<int> > > dummy_cycles_map;
			getCycleStructure(*CSRit,  currentCS, dummy_cycles_map );
			addUniversal(*CSRit);
			/*
			cout<<"Added cs"<<endl;
			printMOLS(partMOLS);*/

			possibleShuffles = genRelevantPermutations(partMOLS[currSquare].front());
			currSquare = (currSquare+1)%k;
			detailedCount[counter].back()++;
			sum = sum+findMOLS4(1, true);
			currSquare = (currSquare-1+k)%k;
			/*cout<<"Return from findmols structure"<<endl;
			printMOLS(partMOLS);*/
			removeUniversal(*CSRit);
			/*cout<<"cs removed"<<endl;
			printMOLS(partMOLS);*/
			currentCS.clear();currentCS.resize(n+1, 0);
			//partMOLS[currSquare].pop_back();
		}

		estimates[counter] = estimates[counter]+ sum;
		return cycleStructureReps.size();

	}
	else{ //not u_0^(1)
		permutation P(identity) ;
		int currUni = partMOLS[currSquare].size();
		int totalBranches = sqSymPossPerms[currSquare][currUni].size();
		int numdone = 0;

		BigInteger sum;
		do{
			if (P.front()<1){
				if (checkFit(P)){
					if (checkOrthogonal(P)){
						addUniversal(P);
						if (getSmallRelCS(partMOLS)){
							if (isSmallest4( )){
								currSquare = (currSquare+1)%k;
								numdone++;
								sum = sum+findMOLS4(1, true);
								currSquare = (currSquare-1+k)%k;
							}
						}
						removeUniversal(P);
					}
				}
			}
			else break;
		}while(next_permutation(P.begin(), P.end()));


			//int avgFeasible = sum/float(numdone);
		estimates[counter] = estimates[counter]+ sum;
			//return (avgFeasible *(float(nofFeasible)/numVisited));
		return numdone;



		}


    }


        return nofFeasible ;

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

	cout<<endl<<"************DONE***************"<<endl;


}

int main(int argc,char *argv[]){
 	// MOLS threemols(8,3);
	string filename = argv[1];
    MOLS threemols(filename);
/*
	string outfile = "out103_1.txt";//+filename;
	std::ofstream out(outfile.c_str());
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf());
*/
	threemols.enumerateMOLS();
	// std::cout.rdbuf(coutbuf); //reset to standard output again
}
