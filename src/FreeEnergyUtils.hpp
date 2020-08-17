#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <unistd.h>

using namespace std;

inline bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

string reverseString(string a){
	string tmp="";
	for(int i=a.length()-1;i>=0;i--){
		tmp = tmp + a[i];
	}
	return tmp;
}

class FreeEnergyUtils{

	//parameters
	private: double Init_GC = 0.98;
	private: double Init_AT = 1.03;
	private: double Init = 2.01;
	private: int window_size = 2;
	private: int primer_size;
	private: int minDGscore;
	private: string forward_primer = "";
	private: string reverse_primer = "";
	private: vector<string> string_arr_list;

	public: FreeEnergyUtils(string last_five_primer, int minDG, int DGlen){
		
		primer_size = DGlen; 
		minDGscore = minDG; 
		
		forward_primer = last_five_primer;

		for (int i = 0; i < forward_primer.length(); i++)
		{
			reverse_primer += comp_gen(forward_primer[i]);
		}
		if (!forward_primer.empty() && !reverse_primer.empty()) {
			build_free_energy_parimeters();
		}
//		else
//			cout << "EMPTY CASE!!" << endl;
	}

	public: vector<string> get_string_arr_list () {
		return string_arr_list;
	}

	public: void build_free_energy_parimeters() {
		for (int i = primer_size-1; i > 0 ; i--) {
			string forward_temp = "";
			string reverse_temp = "";
			for (int j = forward_primer.length()-1-i; j <  forward_primer.length()-1 -i + window_size ; j++) {
				forward_temp += forward_primer[j];
				reverse_temp += reverse_primer[j];
			}
			string_arr_list.push_back(forward_temp + reverse_temp);
		}
	}

	public: char comp_gen(char a) {
		
		switch (a) {
		case 'T':
			return 'A';
		case 'A':
			return 'T';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		default: ;
		
		}
			return a;
	}

	public: bool comp_char(char a, char b) {

		switch (a) {
		case 'T':
			if (b == 'A') {
				return true;
			} else {
				return false;
			}
		case 'A':
			if (b == 'T') {
				return true;
			} else {
				return false;
			}
		case 'C':
			if (b == 'G') {
				return true;
			} else {
				return false;
			}
		case 'G':
			if (b == 'C') {
				return true;
			} else {
				return false;
			}
		default: ;
			break;
		}
		return false;

	}

	public: bool self_complementary() {

		if ((comp_char(forward_primer[1], forward_primer[2]))
				&& comp_char(forward_primer[0], forward_primer[3])) {
			return true;
		}
		if ((comp_char(forward_primer[1], forward_primer[4]))
				&& comp_char(forward_primer[2], forward_primer[3])) {
			return true;
		}
		if ((comp_char(forward_primer[0], forward_primer[4]))
				&& comp_char(forward_primer[1], forward_primer[3])) {
			return true;
		}

		return false;
	}

	public: int selfCmp() {
		
		string tmpP = forward_primer;
		string P = tmpP.substr(tmpP.length()-5, 5);
		int max = 0;

		string rvsP = reverseString(P);
		for (int k = 1; k <= P.length(); k++) {
			int tmpscore = 0;
			for (int i = 0, j = rvsP.length()-k; i < k; i++, j++) {

				if (P[i] == 'A') {
					if (rvsP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'T') {
					if (rvsP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'C') {
					if (rvsP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'G') {
					if (rvsP[j] == 'C')
						tmpscore ++;
					else 
						tmpscore --;
				}
				if (tmpscore < 0) 
					tmpscore = 0;
				if (tmpscore > max)
					max = tmpscore;
			}
		}

		for (int k = P.length()-1; k > 0; k--) {
			int tmpscore = 0;
			for (int i = P.length()-k, j = 0; j < k; i++, j++) {
				if (P[i] == 'A') {
					if (rvsP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'T') {
					if (rvsP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'C') {
					if (rvsP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'G') {
					if (rvsP[j] == 'C')
						tmpscore ++;
					else 
						tmpscore --;
				}
				if (tmpscore < 0) 
					tmpscore = 0;
				if (tmpscore > max)
					max = tmpscore;
			}
		}
		return max;
	}
	

	public: double get_score() {
		double score = 0.0;
		for (int i = 0; i < string_arr_list.size(); i++) {
			string param = string_arr_list[i];
			//strcpy(param,string_arr_list[i].c_str());
			if (param=="AATT")
				score += -1.00;
			else if (param=="ATTA")
				score += -0.88;
			else if (param=="TAAT")
				score += -0.58;
			else if (param=="CAGT")
				score += -1.45;
			else if (param=="GTCA")
				score += -1.44;
			else if (param=="CTGA")
				score += -1.28;
			else if (param=="GACT")
				score += -1.30;
			else if (param=="CGGC")
				score += -2.17;
			else if (param=="GCCG")
				score += -2.24;
			else if (param=="GGCC")
				score += -1.84;
			else if (param=="ACTG")
				score += -1.44;
			else if (param=="AGTC")
				score += -1.28;
			else if (param=="CCGG")
				score += -1.84;
			else if (param=="TTAA")
				score += -1.00;
			else if (param=="TCAG")
				score += -1.30;
			else if (param=="TGAC")
				score += -1.45;
			
		}
		score = score + Init; //Init_AT + Init_GC;
///		String endP = forward_primer.substring(forward_primer.length()-5, forward_primer.length());
//		Constraints Pconstraints = new Constraints(endP);
		
		if (selfCmp() != 0) {
			score += 0.43;
		}
		
//		if (ends_with(string_arr_list[string_arr_list.size()-1],"A") || ends_with(string_arr_list[string_arr_list.size()-1],"T"))
//		{
//			score += 0.05;
//		}
		return score;
	}
	public: bool get_free_energy() {

		if ( minDGscore < get_score() ) {
			return true;
		} else {
			return false;
		}		
	}
};
