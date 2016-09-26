#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>
using namespace std;

int main(int argc, char *argv[]){
	// open files
	int numOfFile = atoi(argv[1]);
	double lfrac = atof(argv[2]);
	double rfrac = atof(argv[3]);
	const double ePerAtto2Amp = 0.1602176462;

	ifstream *infile = new ifstream[2*numOfFile];
	for(int idx=0; idx<numOfFile; idx++){
		char infilename[64];
		sprintf(infilename,"rho_%d",idx);
		infile[idx].open(infilename);
	}
	for(int idx=numOfFile; idx<2*numOfFile; idx++){
		char infilename[64];
		sprintf(infilename,"v_%d",idx-numOfFile);
		infile[idx].open(infilename);
	}


	// read the header
	int *nr3x = new int[numOfFile];
	string tmpline;
	stringstream spliter;
	for(int idx=0; idx<2*numOfFile; idx++){
		getline(infile[idx],tmpline);
		spliter.clear();
		spliter.str(tmpline);
		
		// #Procid=, me_pool, nspin
		spliter >> tmpline >> tmpline >> tmpline;
		
		// nnr, nr1x, nr2x
		if(idx<numOfFile){
			int nnr, nr1x, nr2x;
			spliter >> nnr;
			spliter >> nr1x;
			spliter >> nr2x;
			nr3x[idx] = nnr/(nr1x*nr2x);
		}
	}

	// form the 2D time-dependent charge density
	vector< vector<double> > trho;
	const double minrho = 1e-20;
	while(1){
		bool flag = false;
		vector<double> thisrho;
		vector<int> thisindex;
		double tmprho;
		int tmpindex;

		// read each rho file
		for(int idx=0; idx<numOfFile; idx++){
			// time tag
			if( infile[idx].eof() ){
				flag = true;
				continue;
			}
			getline(infile[idx],tmpline);

			// values
			for(int jdx=0; jdx<nr3x[idx]; jdx++){
				if( infile[idx].eof() ){
					flag = true;
					break;
				}
				getline(infile[idx],tmpline);
				spliter.clear();
				spliter.str(tmpline);
				spliter >> tmpindex >> tmpline >> tmprho;

				// judge whether this value is valid
				if( tmpindex >= thisindex.size() ){
					// if(idx==numOfFile-1 && jdx==nr3x[idx]-1 && tmprho<minrho){
					if(jdx==nr3x[idx]-1 && fabs(tmprho)<minrho){
					} else{
						thisrho.push_back(tmprho);
						thisindex.push_back(tmpindex);
					}
				}
			}
		}

		// decide whether or not to add
		if(flag)
			break;
		else
			trho.push_back(thisrho);
	}

	// report the trho matrix
	ofstream outfile;
	outfile.open("rho.dat");
	outfile.precision(15);
	for(int idx=0; idx<trho.size(); idx++){
		for(int jdx=0; jdx<trho[idx].size()-1; jdx++)
			outfile << trho[idx][jdx] << "\t";
		outfile << trho[idx][trho[idx].size()-1] << endl;
	}
	cout << trho.size() << " time steps for the charge density are successfully found!" << endl;

	// compute the current
	int numOfGrid = trho[0].size();
	int lInt = static_cast<int>( numOfGrid*lfrac );
	int rInt = static_cast<int>( numOfGrid*rfrac );
	
	vector< double > lrho;
	vector< double > rrho;
	vector< double > tmpcur;
	vector< double > cur;

	for(int idx=0; idx<trho.size(); idx++){
		double tmplrho=0.0;
		double tmprrho=0.0;
		for(int jdx=0; jdx<lInt; jdx++){
			tmplrho += trho[idx][jdx];
		}

		for(int jdx=rInt; jdx<numOfGrid; jdx++){
			tmprrho += trho[idx][jdx];
		}

		lrho.push_back(tmplrho);
		rrho.push_back(tmprrho);
		tmpcur.push_back(0.5*(tmplrho-tmprrho));
	}
	cur.push_back(0.0);
	for(int idx=1; idx<tmpcur.size(); idx++){
		cur.push_back(ePerAtto2Amp * (tmpcur[idx]-tmpcur[idx-1]) );
	}

	// report current
	ofstream outfile2;
	outfile2.open("cur.dat");
	outfile2.precision(15);
	for(int idx=0; idx<cur.size(); idx++){
		outfile2 << idx<< "\t" << cur[idx] << "\t" << lrho[idx] << "\t" << rrho[idx] << endl;
        }

	// form the 2D time-dependent KS potential
	vector< vector<double> > tv;
	const double minv = 1e-20;
	while(1){
		bool flag = false;
		vector<double> thisv;
		vector<int> thisindex;
		double tmpv;
		int tmpindex;

		// read each rho file
		for(int idx=numOfFile; idx<2*numOfFile; idx++){
			// time tag
			if( infile[idx].eof() ){
				flag = true;
				continue;
			}
			getline(infile[idx],tmpline);

			// values
			for(int jdx=0; jdx<nr3x[idx-numOfFile]; jdx++){
				if( infile[idx].eof() ){
					flag = true;
					break;
				}
				getline(infile[idx],tmpline);
				spliter.clear();
				spliter.str(tmpline);
				spliter >> tmpindex >> tmpline >> tmpv;

				// judge whether this value is valid
				if( tmpindex >= thisindex.size() ){
					// if(idx==2*numOfFile-1 && jdx==nr3x[idx-numOfFile]-1 && tmpv<minv){
					if(jdx==nr3x[idx-numOfFile]-1 && fabs(tmpv)<minv){
					} else{
						thisv.push_back(13.6057*tmpv);
						thisindex.push_back(tmpindex);
					}
				}
			}
		}

		// decide whether or not to add
		if(flag)
			break;
		else
			tv.push_back(thisv);
	}

	// report the tv matrix
	ofstream outfile3;
	outfile3.open("v.dat");
	outfile3.precision(15);
	for(int idx=0; idx<tv.size(); idx++){
		for(int jdx=0; jdx<tv[idx].size()-1; jdx++)
			outfile3 << tv[idx][jdx] << "\t";
		outfile3 << tv[idx][tv[idx].size()-1] << endl;
	}
	cout << tv.size() << " time steps for the KS potential are successfully found!" << endl;

	// finalize
	for(int idx=0; idx<2*numOfFile; idx++){
		infile[idx].close();
	}
	outfile.close();
	outfile2.close();
	outfile3.close();
	delete[] infile;
	delete[] nr3x;
}
