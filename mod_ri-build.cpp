
#include <iostream>

#include "internal/r_index.hpp"
#include "utils.hpp"
#include "internal/r_index.hpp"
#include "added_files/sampling_sa_builder.hpp"

using namespace ri;
using namespace std;

string out_basename=string();
string input_file=string();
int sa_rate = 512;
bool sais=true;
ulint T = 0;//Build fast index with SA rate = T
bool fast = false;//build fast index
bool hyb = false; //use hybrid bitvectors instead of sd_vectors?

void help(){
	cout << "ri-build: builds the r-index. Extension .ri is automatically added to output index file" << endl << endl;
	cout << "Usage: ri-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>        use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
	cout << "   -divsufsort          use divsufsort algorithm to build the BWT (fast, 7.5n Bytes of RAM). By default,"<<endl;
	cout << "                        SE-SAIS is used (about 4 time slower than divsufsort, 4n Bytes of RAM)."<<endl;
	cout << "   <input_file_name>    input text file." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_basename = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-divsufsort")==0){

		sais = false;

	}else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

	//parse options

    out_basename=string();
    input_file=string();
	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr]);

	if(out_basename.compare("")==0)
		out_basename = string(input_file);

	string idx_file = out_basename;
	idx_file.append(".ri");


	cout << "Building r-index of input file " << input_file << endl;
	cout << "Index will be saved to " << idx_file << endl;

	/*
	string input;

	{

		std::ifstream fs(input_file);
		std::stringstream buffer;
		buffer << fs.rdbuf();

		input = buffer.str();

	}
	*/

	string path = string(out_basename).append(".ri");
	std::ofstream out(path);

	//save flag storing whether index is fast or small
	out.write((char*)&fast,sizeof(fast));


	if(hyb){

		//auto idx = r_index<sparse_hyb_vector,rle_string_hyb>(input,sais);
		//idx.serialize(out);

	}else{
		auto tuple_bwt_and_samples = stool::r_index::SamplingSATBuilder::build(input_file); 
		auto idx = r_index<>(tuple_bwt_and_samples);
		idx.serialize(out);

	}


	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;


	out.close();

}
