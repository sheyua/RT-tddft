#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#ifndef C_VKS_H
#define C_VKS_H
namespace TD_Post{
	class C_Vks{
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
			double *vks;
			
			C_Vks(std::string);
			~C_Vks();
			
			void load(double*, double*);
	};
}
#endif
