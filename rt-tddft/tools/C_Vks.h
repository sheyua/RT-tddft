#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#ifndef C_VKS_H
#define C_VKS_H
namespace TD_Post{
	class C_Vks{
		public:
			int num_proc;
			std::string dump_dir;
	
			C_Vks(std::string);
			~C_Vks();
	};
}
#endif
