#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>  
#include <cmath>
#include <omp.h>

using namespace std;

vector<double> get_col (vector< vector<int> > &ab, int idx){
	vector<double> vec;
	for(int i = 0;i<ab.size();i++){
		vec.push_back(ab[i][idx]);
	}
	return vec;
}

void print_vec (vector<double> &vec){
	for(auto& elem : vec){
		cout<<elem<<',';	
	}
	cout<<'\n';
}

vector< vector<double> > sorenson(vector< vector<int> > &ab){
	int colsize = ab[0].size();
	vector<vector<double> > sorenson_mat(colsize, std::vector<double>(colsize));
	omp_set_dynamic(0);
	omp_set_num_threads(64);
	#pragma omp parallel for
	for(int i=0;i<colsize;i++){
		vector<double> samp1 = get_col(ab,i);
		for(int j=0+i;j<colsize;j++){
			vector<double> samp2 = get_col(ab,j);
			//print_vec(samp1);
			//print_vec(samp2);
			double a_pw = 0;
			double b = 0;
			double c = 0;
			for(int k=0;k<samp1.size();k++){
				double s1 = samp1[k];
				double s2 = samp2[k];
				int action = 0;
				if( s1 > 0 && s2 > 0 ){
					a_pw+=1;
					action = 1;
				}
				if( s1 > 0 && action == 0 ){
					b+=1;
					action = 1;
				}
				if( s2 > 0 && action == 0){
					c+=1;
				}
			}
			//cout << "a: "<<a_pw;
			double a = 2*a_pw;
			//cout << a << ',' << b << ',' << c << '\n' << '\n';
			double new_entry = 1-(a/(a+b+c));
			sorenson_mat[i][j] = new_entry;
			sorenson_mat[j][i] = new_entry;
		}
	}	
	return sorenson_mat;
}

vector< vector<int> > to_matrix (vector< vector<string> > &ia){
	vector<vector<int> > output_mat;
	for(int i = 1;i<ia.size();i++){
		vector<int> rowv;
		for(int j = 1;j<ia[0].size();j++){
			rowv.push_back(stoi(ia[i][j]));
		}
		output_mat.push_back(rowv);
	}
	return output_mat;
}

void to_csv(vector< vector<double> > &sm){
	ofstream outputfile;
	outputfile.open("sorenson.csv");
	int count;
	for (auto& elem: sm){
		count = 0;
		for (auto it = elem.begin(); it != elem.end();++it){
			if(count == 0){
				outputfile << *it;
				count =1;
			}
			else	
				outputfile << ',' << *it;
		}
		outputfile << '\n';
	}
	outputfile.close();	
}

bool String2Int(const std::string& str, int& result)
{
    std::string::const_iterator i = str.begin();
 
    if (i == str.end())
        return false;
 
    bool negative = false;
 
    if (*i == '-')
    {
        negative = true;
        ++i;
 
        if (i == str.end())
            return false;
    }
 
    result = 0;
 
    for (; i != str.end(); ++i)
    {
        if (*i < '0' || *i > '9')
            return false;
 
        result *= 10;
        result += *i - '0';
    }
 
    if (negative)
    {
        result = -result;
    }
 
    return true;
}

int main (int argc, char **argv) {

    std::ifstream f;

    if (argc > 1) {         
        f.open (argv[1]);   
        if (! f.is_open()) {    
            std::cerr << "error: file open failed '" << argv[1] << "'.\n";
            return 1;
        }
    }
    else {  
        std::cerr << "error: insufficient input. <filename> required.\n";
        return 1;
    }

    std::string line, val;                  
    std::vector<std::vector<int> > array;    
    int tester = 0;
    while (std::getline (f, line)) {        
        if(tester == 0){
		tester = 1;
		continue;
	}
	std::vector<int> v;                 
        std::stringstream s (line);
	int t2 = 0;        
	int holder; 
	while (getline (s, val, ',')){       
		if(t2 != 0){
			String2Int(val,holder);
			v.push_back (holder); 
		}
		else 
			t2=1;
	} 
        array.push_back (v);                
    }
    //vector< vector<int> > binabm = to_matrix(array);
    //array.clear();
    auto sormat = sorenson(array);
    to_csv(sormat);
}