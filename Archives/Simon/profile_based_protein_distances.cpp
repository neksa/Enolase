#include <bitset>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>

/**
   == Input format ==
   tab-separated data to STDIN
   col1: profile index
   col2: protein name
   col3: subgroup assignment
 */


using namespace std;
using namespace boost;

void usage(char* name) {
    cerr << "Input taken from STDIN\n";
    cerr << "Usage: " << name << " SFLD|CDD" << endl;
}

int main (int argc, char** argv) {
    if (argc != 2) {
	usage(argv[0]);
	exit(0);
    }
    string type_string(argv[1]);
    enum {CDD,SFLD} type;
    //db_type type;
    if (type_string == string("CDD")){
	type = CDD;
    } else if (type_string == string("SFLD")) {
	type = SFLD;
    } else {
	usage(argv[0]);
	exit(0);
    }
    
    map<string, string> group2abbrev;
    map<string, size_t> protein2index, profile2index;
    map<string, set<string> > protein2profiles, protein2group, protein2abbrev;
    vector<string> index2profile, index2protein, index2label;

    if (type == CDD) {
	group2abbrev["enolase"] = "ENO";
	group2abbrev["galactarate_dehydratase"] = "GALD";
	group2abbrev["glucarate_dehydratase"] = "GLUD";
	group2abbrev["mandelate_racemase"] = "MR";
	group2abbrev["mannonate_dehydratase"] = "MD";
	group2abbrev["methylaspartate_ammonia-lysase"] = "MA";
	group2abbrev["muconate_cycloisomerase"] = "MC";
    } else if (type == SFLD) {
	// don't do anything now
    }
    
    /** read input (from matches2subgroup.pl or compatible) **/
    do {
	string profile, protein, group;
	cin >> profile >> protein >> group;
	if (profile == "" || protein == "") {
	    continue;
	}
	protein2group[protein].insert(group);
	protein2profiles[protein].insert(profile);

	//assign an index to each unique profile
	if (profile2index.count(profile) == 0) {
	    size_t n = profile2index.size(); 
	    profile2index[profile] = n;
	    index2profile.push_back(profile);
	    assert(index2profile.size() == profile2index.size());
	}
    } while (cin.good());

    /** get rid of proteins that match only one profile **/
    vector<map<string, set<string> >::iterator> remove_these;
    for (map<string, set<string> >::iterator it = protein2profiles.begin();
	 it != protein2profiles.end(); ++it) {
	if (it->second.size() <= 1) {
	    remove_these.push_back(it);
	}
    }
    typedef map<string, set<string> >::iterator p2p_it;
    BOOST_FOREACH(const p2p_it& it, remove_these) {
	protein2profiles.erase(it);
    }
    
    /** assign an index to each unique protien **/
    size_t n_proteins = protein2profiles.size();
    size_t index = 0;
    index2protein.resize(n_proteins,"");
    typedef pair<string, set<string> > p2p_entry;
    BOOST_FOREACH(const p2p_entry& kv_pair, protein2profiles) {
	protein2index[kv_pair.first] = index;
	index2protein[index++] = kv_pair.first;
    }

    /** describe each protein as a binary string representing the
	profiles it matches **/
    size_t n_profiles = profile2index.size();
    index2label.resize(n_proteins);
    vector<dynamic_bitset<> > conf(n_proteins,dynamic_bitset<>(n_profiles));
    cout << "Analyzing " << n_proteins << " proteins and " 
	 << n_profiles << " profiles" << endl;
    for (size_t p = 0; p < n_proteins; ++p) {
	string protein = index2protein[p];
	set<string> group = protein2group[protein];
	stringstream label;
	
	//each protein needs a unique label too (strictly required by R
	//in subsequent analysis)
	for (set<string>::const_iterator it = protein2group[protein].begin();
	     it != protein2group[protein].end();) {
	    if (type == SFLD) {
		if (group2abbrev.count(*it)) {
		    label << group2abbrev[*it];
		} else {
		    label << "??";
		}
	    } else if (type == CDD) {
		label << *it;
	    }
	    if (++it != protein2group[protein].end()) {
		label << ",";
	    }
	}
	label << ":" << protein;


	for (size_t q = 0; q < n_profiles; ++q) {
	    string profile = index2profile[q];
	    if (protein2profiles[protein].count(profile)) {
		conf[p][q] = 1;
		label << "_" << profile;
	    } else {
		conf[p][q] = 0;
	    }
	}
	index2label[p] = label.str();
    }

    /** calculate and output distance matrix **/
    ofstream fout("dist_matrix");
    for (size_t p = 0; p < n_proteins; ++p) {
	fout << index2label[p] << "\t";
    }
    fout << endl;
    for (size_t p = 0; p < n_proteins; ++p) {
	fout << index2label[p] << "\t";
	for (size_t q = 0; q < n_proteins; ++q) {
	    /** hamming distance **/
	    //fout << (conf[p] ^ conf[q]).count() << "\t";
	    /** jaccard distance **/
	    fout << 1. - (conf[p] & conf[q]).count()/(double)(conf[p] | conf[q]).count() << "\t";
	}
	fout << endl;
    }
    fout.close();

    
    /** output protein to group assignments **/
    /*
    fout.open("group_map");
    for (size_t p = 0; p < n_proteins; ++p) {
	string protein = index2protein[p];
	fout << protein << "\t" << protein2group[protein] << endl;
    }    
    fout.close();
    */
    return 0;
}
