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

// A sample validator that accepts all results

#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include "error_numbers.h"
#include "boinc_db.h"
#include "sched_util.h"
#include "validate_util2.h"
#include "validate_util.h"

using std::string;
using std::vector;

struct DATA {
    int num_MOLS;
    vector<long long int> branch_counts;
};

std::vector<std::string> & split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

int init_result(RESULT& result, void*& data) {
	FILE* f;
    OUTPUT_FILE_INFO fi;
    int  retval;

    char sc[10000], sc1[10000];
    string s;
    int num_MOLS, k, n;
    bool read = true;
    vector<long long int> branch_count;

    retval = get_output_file_path(result, fi.path);
    if (retval) return retval;
    retval = try_fopen(fi.path.c_str(), f, "r");
    if (retval) return retval;

    while (read){
    	 if (feof(f)) break;
    	 fgets(sc, 10000, f);
    	 if (sc[0]=='#'&& sc[1]=='B'){ //#Branchcounts
        	 //fgets(sc, 10000, f);
    		 //vector<string> data  = split((string)sc, ' ');
    	 	 fgets(sc1, 10000, f);
    		 vector<string> branches  = split((string) sc1, ' ');
    		 branch_count.resize(branches.size(), 0);
    		 for (int i =0; i<branches.size(); i++)
    			 branch_count[i] = atoll(branches[i].c_str()); //first element is a #

    		 //read = false;
    	 }
    	 if (sc[0]=='#'&& sc[1]=='M'){ //#Branchcounts
    	         	 fgets(sc, 10000, f);
    	     		 vector<string> data  = split((string)sc, ' ');
    	     		 num_MOLS = atol(data[0].c_str());
    	     		 read = false;
    	}

    }
    if (read) return ERR_XML_PARSE;
    //now read branches everywhere

    retval = fclose(f);
    if (retval) return retval;

    DATA* dp = new DATA;
    dp->num_MOLS = num_MOLS;
    dp->branch_counts = branch_count;
    data = (void*) dp;
    return 0;
}

int compare_results(RESULT& r1, void* _data1, RESULT const& r2, void* _data2, bool& match) {
    DATA* data1 = (DATA*)_data1;
    DATA* data2 = (DATA*)_data2;
    match = true;
    if (data1->num_MOLS != data2->num_MOLS) match = false;
    if (data1->branch_counts.size() != data2->branch_counts.size()) match = false;
    else
    	for (int i =0; i<data1->branch_counts.size(); i++)
    		if (data1->branch_counts[i] != data2->branch_counts[i]) match = false;
    //if (fabs(data1->x - data2->x) > 0.01) match = false;
    return 0;
}

int cleanup_result(RESULT const& r, void* data) {
    if (data) delete (DATA*) data;
    return 0;
}

/*int main(int argc, char** argv) {
	return 0;
}*/


const char *BOINC_RCSID_f3a7a34795 = "$Id$";
