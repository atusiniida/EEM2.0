//
// wgetopt.h : simple getopt for eemCPP
//
//
#ifndef Wgetopt_H
#define Wgetopt_H

#include <iostream>
#include <string>
using namespace std;

#define no_argument       0
#define required_argument 1
#define optional_argument 2

struct option{
	const char *name;
	int         has_arg;
	int        *flag;
	int         val;
};

option options[] = {
	{"absrad",     required_argument, 0, 'A'},
	{"relrad",     required_argument, 0, 'R'},
	//{"cor",        required_argument, 0, 'C'},
	{"pcut",       required_argument, 0, 'p'},
	{"itr",        required_argument, 0, 'i'},
	//{"absolute",   no_argument,       0, 'a'},
	{"outfile",    required_argument, 0, 'o'},
	{"mingeneset", required_argument, 0, 'm'},
	{"maxgeneset", required_argument, 0, 'M'},
	//{"expmodset",  required_argument, 0, 'e'},
	{"log",        required_argument, 0, 'l'},
	{"recnull",    no_argument,       0, 'r'},
	{0,            0,                 0, 0  }
};

char*  optarg;
int optind, optopt;// opterr;

static void usage(const string &cmdline, const struct option *options)
{
	cerr << cmdline << endl;
	cerr << "Options:" << endl;
	while (options->name) {
		cerr << "--" << options->name;
		if (options->val) {
			cerr << ", -" << char(options->val);
		}
		if (options->has_arg == required_argument) {
			cerr << " arg";
		}
		//
		cerr << endl;
		options++;
	}
}

int search_optstring(string& optstring, string& str)
{
	int flag;
	//
	// Ex. optstring == "A:R:C:p:i:ao:m:M:e:l:rc:"
	//
	string::size_type pos;
	pos = optstring.find(str);

	if(string::npos != pos){
		if(optstring.length() < pos+1){
			if(':' == optstring[pos+1]){
				flag=1;//':' exist 
			}else{
				flag=2;
			}
		}else{
			flag=2;
		}
	}else{
		optopt = str[0];
		flag=0;// not option 
	}

	return flag;
}

int getopt_long(int argc, char* argv[], string optstring, option opt[], int* option_index)
{
	int ret(-1);
	bool existOption(false);

	for(int i=0; i < argc; i++){

		if(argv[i][0]=='-'){
			existOption=true;

			string str = argv[i];
			str.erase(0, 1);
			
			if(str.length() > 0){
				ret = str[0];

				int flag = search_optstring(optstring, str);
				
				if(flag==1){// exist ? option argument
					++i;
					optarg = argv[i];
				}else if(flag==0){
					ret = '?';
					optarg = 0;
				}else if(flag==2){
					optarg = 0;
				}
			}else{
				ret = '?';
				optarg = 0;
			}
		}
		optind = i+1;// optind : next argv index

		if(existOption) break;
	};

	if(!existOption){
		optind = 1;
	}

	return ret;
}

#endif // Wgetopt_H
