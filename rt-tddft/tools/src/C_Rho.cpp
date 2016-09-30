#include "C_Rho.h"

using namespace TD_Post;

C_Rho::C_Rho(std::string dump_dir){
	this->dump_dir = dump_dir;
	
	// check the number of processors
	num_proc = 0;
	while(1){
		// determine a filename
		/*std::string filename = "rho"+std::to_string(num_proc);
		if(dump_dir.back()=='/')*/
		std::ostringstream os;
		os << num_proc;
		std::string filename = "rho_"+os.str();
		if(dump_dir[dump_dir.length()-1]=='/')
			filename = dump_dir+filename;
		else
			filename = dump_dir+"/"+filename;
		
		std::ifstream f(filename.c_str());
		if(f.good())
			num_proc++;
		else
			break;
		f.close();
	}
	
	// report number of processors and dump file location
	std::cout << num_proc << " charge density files are found under " << dump_dir << std::endl;
	
	// return if num_proc == 0;
	if(num_proc == 0) return;
	
	// open all files
	std::ifstream *infile = new std::ifstream[num_proc];
	for(int idx=0; idx<num_proc; idx++){
		std::ostringstream os;
		os << idx;
		std::string filename = "rho_"+os.str();
		if(dump_dir[dump_dir.length()-1]=='/')
			filename = dump_dir+filename;
		else
			filename = dump_dir+"/"+filename;
		
		infile[idx].open(filename.c_str());
	}
	
	// init the arrarys
	idz_start.resize(num_proc);
	idz_end.resize(num_proc);
	idzs_length.resize(num_proc);
	
	// get nspin, dt,  and c
	std::string tmpLine;
	std::stringstream spliter;
	std::getline(infile[0], tmpLine);
	spliter.clear();
	spliter.str(tmpLine);
	spliter >> tmpLine >> tmpLine \
	        >> nspin >> dt \
	        >> idz_start[0] >> idz_end[num_proc-1] \
	        >> c;
	
	// fill in the idz
	for(int idx=1; idx<num_proc; idx++){
		std::getline(infile[idx], tmpLine);
		spliter.clear();
		spliter.str(tmpLine);
		spliter >> tmpLine >> tmpLine \
		        >> tmpLine >> tmpLine \
		        >> idz_start[idx] >> idz_end[idx-1] \
		        >> tmpLine;
	}
	for(int idx=0; idx<num_proc; idx++){
		std::getline(infile[idx], tmpLine);
		spliter.clear();
		spliter.str(tmpLine);
		spliter >> tmpLine >> init_step;
		int tmpLen = 0;
		while(infile[idx].good()){
			std::getline(infile[idx], tmpLine);
			if(tmpLine[0] != 'T')
				tmpLen++;
			else
				break;
		}
		idzs_length[idx] = tmpLen;
	}
	
	// get num_step, init_step, num_mids
	init_step++;
	if(infile[0].good())
		num_step = 1;
	else
		num_step = 0;
	while(infile[0].good()){
		std::getline(infile[0], tmpLine);
		if(tmpLine[0] == 'T')
			num_step++;
	}
	num_mids = idz_end[num_proc-1];
	
	// report simulation summary
	std::cout << "Number of spin: " << nspin << std::endl;
	std::cout << "Delta t (atto): " << dt << std::endl;
	std::cout << "Lattice c (A):  " << c << std::endl;
	std::cout << "Initial step:   " << init_step << std::endl;
	std::cout << "Number of step: " << num_step << std::endl;
	std::cout << "Number of mids: " << num_mids << std::endl;
	
	// finalize
	for(int idx=0; idx<num_proc; idx++)
		infile[idx].close();
	delete[] infile;
}

C_Rho::~C_Rho(){
}

void C_Rho::load(double *vbias, double *rho){
	// passing values to local pointers
	this->vbias = vbias;
	this->rho = rho;
	
	// open all files
	std::ifstream *infile = new std::ifstream[num_proc];
	for(int idx=0; idx<num_proc; idx++){
		std::ostringstream os;
		os << idx;
		std::string filename = "rho_"+os.str();
		if(dump_dir[dump_dir.length()-1]=='/')
			filename = dump_dir+filename;
		else
			filename = dump_dir+"/"+filename;
		
		infile[idx].open(filename.c_str());
	}
	
	// read data
	std::string tmpLine;
	std::stringstream spliter;
	for(int idx=0; idx<num_proc; idx++){
		// remove info line
		std::getline(infile[idx], tmpLine);
		
		// get vbias and rho
		for(int tdx=0; tdx<=num_step; tdx++){
			// get vbias
			std::getline(infile[idx], tmpLine);
			spliter.clear();
			spliter.str(tmpLine);
			spliter >> tmpLine >> tmpLine >> vbias[tdx];
			
			// get rhos
			int idz_true = idz_end[idx] - idz_start[idx];
			for(int idzs=0; idzs<idzs_length[idx]; idzs++){
				int tmps, tmpz;
				double tmp;
				std::getline(infile[idx], tmpLine);
				spliter.clear();
				spliter.str(tmpLine);
				spliter >> tmps >> tmp >> tmpz;
				tmps--;
				
				// check tmpz range
				if(tmpz < idz_end[idx])
					rho[tmps*(num_step+1)*num_mids + tdx*num_mids + tmpz] = tmp;
			}
		}
	}
	
	// finalize
	for(int idx=0; idx<num_proc; idx++)
		infile[idx].close();
	delete[] infile;
}

void C_Rho::comp_cur(double* cur){
	// passing value to local pointer
	this->cur = cur;
	
	// start from init_step
	for(int is=0; is<nspin; is++)
		for(int tdx=0; tdx<num_step; tdx++){
			// left edge
			cur[is*num_step*num_mids + tdx*num_mids] = \
				      rho[is*(num_step+1)*num_mids + (tdx+1)*num_mids] \ 
				    - rho[is*(num_step+1)*num_mids + tdx*num_mids];
			// rest values
			for(int idz=1; idz<num_mids; idz++)
				cur[is*num_step*num_mids + tdx*num_mids + idz] = \
				      cur[is*num_step*num_mids + tdx*num_mids + idz-1] \
				    + rho[is*(num_step+1)*num_mids + (tdx+1)*num_mids + idz] \ 
				    - rho[is*(num_step+1)*num_mids + tdx*num_mids +idz];
		}
	
	// with coefficients
	const double ePerAtto2muAmp = 0.1602176462e+6;
	for(int is=0; is<nspin; is++)
		for(int tdx=0; tdx<num_step; tdx++)
			for(int idz=1; idz<num_mids; idz++)
				cur[is*num_step*num_mids + tdx*num_mids + idz] *= ePerAtto2muAmp;
}
