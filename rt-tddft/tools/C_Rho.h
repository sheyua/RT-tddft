#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#ifndef C_RHO_H
#define C_RHO_H
namespace TD_Post{
	class C_Rho{
		public:
			int num_proc;
			std::string dump_dir;
	
			C_Rho(std::string);
			~C_Rho();

			void load();
	};
}
#endif
