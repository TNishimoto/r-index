// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>

#include "internal/r_index.hpp"

#include "internal/utils.hpp"

using namespace ri;
using namespace std;

string check = string();//check occurrences on this text
bool hyb=false;
string ofile;

void help(){
	cout << "ri-locate: locate all occurrences of the input patterns." << endl << endl;

	cout << "Usage: ri-locate [options] <index> <patterns>" << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file (must be the same indexed)" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   -o <ofile>   write pattern occurrences to this file (ASCII)" << endl;
	cout << "   <index>      index file (with extension .ri)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-c")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		ofile = string(argv[ptr]);
		ptr++;

	}/*else if(s.compare("-h")==0){

		hyb=true;

	}*/else{

		cout << "Error: unknown option " << s << endl;
		help();

	}

}


template<class idx_t>
void locate(std::ifstream& in, string patterns){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    string text;
    bool c = false;

    ofstream out;

    if(ofile.compare(string())!=0){

    	out = ofstream(ofile);

    }

    if(check.compare(string()) != 0){

    	c = true;

		ifstream ifs1(check);
		stringstream ss;
		ss << ifs1.rdbuf();//read the file
		text = ss.str();

    }

    auto t1 = high_resolution_clock::now();

    idx_t idx;

	idx.load(in);


	cout << "searching patterns ... " << endl;
	ifstream ifs(patterns);

	//read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
	string header;
	std::getline(ifs, header);

	ulint n = get_number_of_patterns(header);
	ulint m = get_patterns_length(header);

	uint last_perc = 0;

	ulint occ_tot=0;


	std::vector<uint64_t> backward_search_time_vector;
	std::vector<uint64_t> computing_sa_time_vector;

	auto t2 = high_resolution_clock::now();

	//extract patterns from file and search them in the index
	for(ulint i=0;i<n;++i){

		uint perc = (100*i)/n;
		if(perc>last_perc){
			cout << perc << "% done ..." << endl;
			last_perc=perc;
		}

		string p = string();

		for(ulint j=0;j<m;++j){
			char c;
			ifs.get(c);
			p+=c;
		}

		//cout << "locating " << idx.occ(p) << " occurrences of "<< p << " ... " << flush;

		auto t3 = high_resolution_clock::now();
		auto res = idx.count_and_get_occ(p);
		auto t4 = high_resolution_clock::now();
		auto OCC = idx.locate_all(res);	//occurrences
		auto t5 = high_resolution_clock::now();


		uint64_t backward_search_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();
		uint64_t computing_sa_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count();
		backward_search_time_vector.push_back(backward_search_time);
		computing_sa_time_vector.push_back(computing_sa_time);

		if(ofile.compare(string())!=0){

			sort(OCC.begin(),OCC.end());

			for(auto x:OCC) out << (int)x << endl;

		}

		occ_tot += OCC.size();

		if(c){//check occurrences

			//remove duplicates, if any (there shouldn't be!)
			sort(OCC.begin(),OCC.end());
			auto it = unique(OCC.begin(),OCC.end());
			OCC.resize(std::distance(OCC.begin(),it));

			if(OCC.size() != idx.occ(p)){

				cout << "Error: wrong number of located occurrences: " << OCC.size() << "/" << idx.occ(p) << endl;
				exit(0);

			}

			for(auto o:OCC){

				//for(auto c : p) cout << int(c) << " ";

				if(text.substr(o,p.size()).compare(p) != 0){

					cout << "Error: wrong occurrence: " << o << " ("  << occ_tot << " occurrences"  << ") "<< endl;
					for(auto c : text.substr(o,p.size())) cout << int(c) << " ";
					cout << "  /  ";
					for(auto c : p) cout << int(c) << " ";
					cout << endl;

					break;

					//exit(0);

				}

			}

		}

	}

	double occ_avg = (double)occ_tot / n;

	cout << endl << occ_avg << " average occurrences per pattern" << endl;

	ifs.close();

	//auto t3 = high_resolution_clock::now();

	//printRSSstat();

	uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Load time : " << load << " milliseconds" << endl;

	for(uint64_t i = 0; i < backward_search_time_vector.size();i++){
		std::cout << "LOCATE, " << i << ", backward search (nanoseconds), " << backward_search_time_vector[i] << ", computing sa (nanoseconds), " << computing_sa_time_vector[i]  << std::endl;
	}


    uint64_t _total_backward_search_time = std::reduce(std::begin(backward_search_time_vector), std::end(backward_search_time_vector));
    uint64_t _average_backward_search_time = _total_backward_search_time / backward_search_time_vector.size();
    std::cout << "Total backward search time" << ": \t\t\t\t\t" << (_total_backward_search_time/1000) << " " << "microseconds" << std::endl;
    std::cout << "\t" << "Average: \t\t\t\t\t" << (_average_backward_search_time/1000) << " " << "microseconds" << std::endl;

    uint64_t _total_computing_sa_time = std::reduce(std::begin(computing_sa_time_vector), std::end(computing_sa_time_vector));
    uint64_t _average_computing_sa_time = _total_computing_sa_time / computing_sa_time_vector.size();
    std::cout << "Total computing sa time" << ": \t\t\t\t\t" << (_total_computing_sa_time/1000) << " " << "microseconds" << std::endl;
    std::cout << "\t" << "Average: \t\t\t\t\t" << (_average_computing_sa_time/1000) << " " << "microseconds" << std::endl;

	uint64_t  _total_search_time = _total_backward_search_time + _total_computing_sa_time;
    uint64_t _average_search_time = _total_backward_search_time / computing_sa_time_vector.size();
    std::cout << "Total search time" << ": \t\t\t\t\t" << (_total_search_time/1000) << " " << "microseconds" << std::endl;
    std::cout << "\t" << "Average: \t\t\t\t\t" << (_average_search_time/1000) << " " << "microseconds" << std::endl;






	cout << "number of patterns n = " << n << endl;
	cout << "pattern length m = " << m << endl;
	cout << "total number of occurrences  occ_t = " << occ_tot << endl;

	//uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

	//cout << "Total time : " << search << " milliseconds" << endl;
	//cout << "Search time : " << (double)search/n << " milliseconds/pattern (total: " << n << " patterns)" << endl;
	//cout << "Search time : " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;

}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	string idx_file(argv[ptr]);
	string patt_file(argv[ptr+1]);

	std::cout << "file: " << argv[ptr+1] << std::endl;

	std::ifstream in(idx_file);

	bool fast;

	//fast or small index?
	in.read((char*)&fast,sizeof(fast));

	cout << "Loading r-index" << endl;

	if(hyb){

		//locate<r_index<sparse_hyb_vector,rle_string_hyb> >(in, patt_file);

	}else{

		locate<r_index<> >(in, patt_file);

	}

	in.close();

}
