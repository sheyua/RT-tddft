#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#ifndef C_RHO_H
#define C_RHO_H
namespace TD_Post{
	class C_Rho{
		private:
			int num_proc;
			std::vector<int> idz_start;
			std::vector<int> idz_end;
			std::vector<int> idzs_length;
		
		public:
			std::string dump_dir;
			int nspin;
			double dt;
			double c;
			int init_step;
			int num_step;
			int num_mids;
			double *vbias;
			double *rho;
			double *cur;
			
			C_Rho(std::string);
			~C_Rho();
			
			void load(double*, double*);
			void comp_cur(double*);
	};
}
#endif
