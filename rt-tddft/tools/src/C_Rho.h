#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifndef C_RHO_H
#define C_RHO_H
namespace TD_Post{
	class C_Rho{
		private:
			int num_proc;
			std::vector<int> idz_start;
			std::vector<int> idz_end;
			std::vector<int> idz_length;
		
		public:
			std::string dump_dir;
			int nspin;
			double dt;
			double c;
			int num_step;
			int init_step;
			
			C_Rho(std::string);
			~C_Rho();
			
			void load();
	};
}
#endif
