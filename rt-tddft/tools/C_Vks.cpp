#include "C_Vks.h"

using namespace TD_Post;

C_Vks::C_Vks(std::string dump_dir){
	this->dump_dir = dump_dir;
	
	// check the number of processors
	num_proc = 0;
	while(1){
		// determine a filename
		/*std::string filename = "rho"+std::to_string(num_proc);
		if(dump_dir.back()=='/')*/
		std::ostringstream os;
		os << num_proc;
		std::string filename = "vks_"+os.str();
		if(dump_dir[dump_dir.length()-1]=='/')
			filename = dump_dir+filename;
		else
			filename = dump_dir+"/"+filename;

		std::ifstream f(filename.c_str());
		if(f.good())
			num_proc++;
		else
			break;
	}

	// report number of processors and dump file location
	std::cout << num_proc << " Kohn-Sham potential files are found under " << dump_dir << std::endl;
}

C_Vks::~C_Vks(){}
