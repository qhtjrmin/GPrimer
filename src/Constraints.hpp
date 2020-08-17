#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* reverseChar(char* a){
	int len = strlen(a);
	char* tmp = new char[len + 1];
	memset(tmp, 0, len + 1);
	int j = len - 1;
	for(int i = 0; i < len; i++){
		tmp[i] = a[j];
		j--;
	}
	return tmp;
}

template <typename T>
T absol(T a){
	if(a >= 0)
		return a;
	else
		return (-a);
}

class PrimerTm {

	private: char* Primer;
	private: int P_len;
	private: double dna_conc;	// Annealing oligo concentration  (DNA concentration) (nanomolar)
	private: double salt_conc;	// Salt concentration (concentration of monovalent cations) (millimolar)
	private: double divalent_conc; // Concentration of divalent cations (millimolar)
	private: double dntp_conc ; // Concentration of dNTPs (millimolar)


	// Primer3 default values are following
	// Annealing oligo concentration: 50.0
	// Salt concentration: 50.0
	// Concentration of divalent cations: 1.5
	// Concentration of dNTPs: 0.6
	public: PrimerTm (char* pri) {//, double dnaConc, double saltConc, double divalentConc, double dntpConc) {
		P_len = strlen(pri);
		Primer = new char[strlen(pri) + 1];
		strcpy(Primer, pri);
		dna_conc = 50.0;
		salt_conc = 50.0;
		divalent_conc = 1.5;
		dntp_conc = 0.6;
	}

	public: ~PrimerTm(){
		delete[] Primer;
		Primer = NULL;
	}

	//public oligotm(String pri) {
		// TODO Auto-generated constructor stub
	//}

	public: double getEnthalpy () {
		double h = 0.0;

		h += 0.2;

		char lStr = Primer[P_len - 1];

		if (lStr == 'A' || lStr == 'T') {
			h += 2.2;
		}

		char* subtmp = new char[3];
		memset(subtmp, 0, 3);
		for (int i = 0; i < P_len - 1; i++){
			memcpy(subtmp, Primer + i, 2);
//			cout << subtmp << endl;
			if (!strcmp(subtmp, "AA")) {
				h += -7.6;
			}
			else if (!strcmp(subtmp, "AC")) {
				h += -8.4;
			}
			else if (!strcmp(subtmp, "AG")) {
				h += -7.8;
			}
			else if (!strcmp(subtmp, "AT")) {
				h += -7.2;
			}
			else if (!strcmp(subtmp, "CA")) {
				h += -8.5;
			}
			else if (!strcmp(subtmp, "CC")) {
				h += -8.0;
			}
			else if (!strcmp(subtmp, "CG")) {
				h += -10.6;
			}
			else if (!strcmp(subtmp, "CT")) {
				h += -7.8;
			}
			else if (!strcmp(subtmp, "GA")) {
				h += -8.2;
			}
			else if (!strcmp(subtmp,"GC")) {
				h += -9.8;
			}
			else if (!strcmp(subtmp, "GG")) {
				h += -8.0;
			}
			else if (!strcmp(subtmp, "GT")) {
				h += -8.4;
			}
			else if (!strcmp(subtmp, "TA")) {
				h += -7.2;
			}
			else if (!strcmp(subtmp, "TC")) {
				h += -8.2;
			}
			else if (!strcmp(subtmp, "TG")) {
				h += -8.5;
			}
			else if (!strcmp(subtmp, "TT")) {
				h += -7.6;
			}
		}
		delete[] subtmp;
		subtmp = NULL;
		return h;
	}

	public: static double divalent_to_monovalent(double divalent, double dntp) {
		// TODO Auto-generated method stub
		if(divalent == 0) dntp = 0;
		if(divalent < dntp)
		/* According to theory, melting temperature does not depend on	divalent cations */
		   	divalent = dntp;
		return 120 * (sqrt(divalent - dntp));
	}

	public: double getEntropy () {
		double s = 0.0;

		s = s + (-5.7);

		salt_conc = salt_conc + divalent_to_monovalent(divalent_conc, dntp_conc);

		s += 0.368 * (double)(P_len - 1) * log(salt_conc / 1000.0);

		char lStr = Primer[P_len - 1];
		if (lStr == 'A' || lStr == 'T') {
			s += 6.9;
		}

		char* subtmp = new char[3];
		memset(subtmp, 0, 3);
		for (int i = 0; i < P_len-1; i++) {

			memcpy(subtmp, Primer + i, 2);

			if (!strcmp(subtmp, "AA")) {
				s += -21.3;
			}
			else if (!strcmp(subtmp, "AC")) {
				s += -22.4;
			}
			else if (!strcmp(subtmp, "AG")) {
				s += -21.0;
			}
			else if (!strcmp(subtmp, "AT")) {
				s += -20.4;
			}
			else if (!strcmp(subtmp, "CA")) {
				s += -22.7;
			}
			else if (!strcmp(subtmp, "CC")) {
				s += -19.9;
			}
			else if (!strcmp(subtmp, "CG")) {
				s += -27.2;
			}
			else if (!strcmp(subtmp, "CT")) {
				s += -21.0;
			}
			else if (!strcmp(subtmp, "GA")) {
				s += -22.2;
			}
			else if (!strcmp(subtmp, "GC")) {
				s += -24.4;
			}
			else if (!strcmp(subtmp, "GG")) {
				s += -19.9;
			}
			else if (!strcmp(subtmp, "GT")) {
				s += -22.4;
			}
			else if (!strcmp(subtmp, "TA")) {
				s += -21.3;
			}
			else if (!strcmp(subtmp, "TC")) {
				s += -22.2;
			}
			else if (!strcmp(subtmp, "TG")) {
				s += -22.7;
			}
			else if (!strcmp(subtmp, "TT")) {
				s += -21.3;
			}
		}
		delete[] subtmp;
		subtmp = NULL;
		return s;
	}


	public: double getTM () {
		double tm = 0.0;

		double h = getEnthalpy();
		double s = getEntropy();

		dna_conc = dna_conc * 0.000000001;
		tm = ((1000 * h) / (s + (1.987 * log((dna_conc) / 4)))) - 273.15;

		return tm;
	}

	private: char* Pos;

	public: char* GetPrimer() {
		// TODO Auto-generated method stub
		return Primer;
	}

	public: char* GetPos() {
		// TODO Auto-generated method stub
		return Pos;
	}
};

class Constraints {

	private: char* Primer;
	private: char* rPrimer;
	private: int PrimerPos;
	private: int rPrimerPos;
	private: int P_len;
	private: int rP_len;

	//for single 
	public: Constraints (char* P) {
		P_len = strlen(P);
		Primer = new char[P_len + 1];
		strcpy(Primer, P);
		rPrimer = new char[1];
        PrimerPos = 0;
        rPrimerPos = 0;
	}

	//for pair
	public: Constraints (char* fP, char* rP, int fPos, int rPos) {
		P_len = strlen(fP);
		rP_len = strlen(rP);
		Primer = new char[P_len + 1];
		rPrimer = new char[rP_len + 1];
		PrimerPos = fPos;
		rPrimerPos = rPos;
		strcpy(Primer, fP);
		strcpy(rPrimer, rP);
	}

	public: ~Constraints(){
		delete[] Primer;
		delete[] rPrimer;
		Primer = NULL;
		rPrimer = NULL;
	}	

	public: char* getPrimer () {
		return Primer;
	}
	
	public: char* getrPrimer() {
		return rPrimer;
	}
	
	public: int getfPos() {
		return PrimerPos;
	}

	public: int getrPos() {
		return rPrimerPos;
	}
	
	
	//single filtering constraints
	//melting temperature 
	public: double MeltTemp (char* P) {
		PrimerTm pritm = PrimerTm(P);		
		double tm = pritm.getTM();
		return tm ; 
	}

	//GC content 
	public: double GCcontent() {
		char* P = getPrimer();	
		int CountOfGC = 0;

		for (int cnt = 0; cnt < P_len; cnt++) {
			if ((P[cnt] == 'G') || (P[cnt] == 'C')) {
				CountOfGC++;
			}
		}
		double val = ((double) CountOfGC / (double) P_len) * 100.0;
		P = NULL;
		return val;
	}

	//self-complementary (not contiguous, local alignment)
	public: int selfCmp () {
		
		char* P = getPrimer();
		int max = 0;

		char* rvsP = reverseChar(P);
		for (int k = 1; k <= P_len; k++) {
			int tmpscore = 0;
			for (int i = 0, j = strlen(rvsP) - k; i < k; i++, j++) {
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

		for (int k = P_len - 1; k > 0; k--) {
			int tmpscore = 0;
			for (int i = P_len - k, j = 0; j < k; i++, j++) {
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
		delete[] rvsP;
		P = NULL;
		rvsP = NULL;
		return max;
	}

	//self-complementary (contiguous)
	public: int cont_selfCmp (int maxSC) {
		int max = 0;
		char* P = getPrimer();
		char* rptmp = new char[P_len + 1];
		memset(rptmp,0,P_len+1);
	
		for (int n = 0; n < P_len; n++) {
			if (P[n] == 'A') 
				rptmp[n] = 'T';
			else if (P[n] == 'T') 
				rptmp[n] = 'A';
			else if (P[n] == 'C') 
				rptmp[n] = 'G';
			else if (P[n] == 'G') 
				rptmp[n] = 'C';
		}
		
		char* rvsP = reverseChar(rptmp);
		
		for (int k = maxSC; k <= P_len - maxSC; k++) {
			int tmpscore = 0;
			for (int i = 0, j = strlen(rvsP) - k; i < k; i++, j++) {
				if (P[i] == 'A') {
					if (rvsP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore = 0;
				}
				else if (P[i] == 'T') {
					if (rvsP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore = 0;
				}
				else if (P[i] == 'C') {
					if (rvsP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore = 0;
				}
				else if (P[i] == 'G') {
					if (rvsP[j] == 'C')
						tmpscore ++;
					else 
						tmpscore = 0;
				}

				if (tmpscore > max)
					max = tmpscore; 
			}
		}

		if (max < maxSC) {			
			for (int k = maxSC; k <= P_len-maxSC; k++) {
				int tmpscore = 0;
				for (int i = P_len-k, j = 0; j < k; i++, j++) {
					if (P[i] == 'A') {
						if (rvsP[j] == 'T')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					else if (P[i] == 'T') {
						if (rvsP[j] == 'A')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					else if (P[i] == 'C') {
						if (rvsP[j] == 'G')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					else if (P[i] == 'G') {
						if (rvsP[j] == 'C')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					if (tmpscore > max)
						max = tmpscore; 
				}
			}
		}
		delete[] rvsP;
		delete[] rptmp;
		rvsP = NULL;
		rptmp = NULL;
		P = NULL;
		return max;
	}

	// end self-complementary (not contiguous)
	public: int endSelfCmp () {

		char* P = getPrimer();
	
		char* selfRvs = reverseChar(P);

		int max = 0;
		int len = P_len / 2;

		for (int k = P_len - 1; k > len; k--) {
			int tmpscore = 0;
			for (int i = P_len-k, j = 0; j < k; i++, j++) {
				if (P[i] == 'A') {
					if (selfRvs[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'T') {
					if (selfRvs[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'C') {
					if (selfRvs[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'G') {
					if (selfRvs[j] == 'C')
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
		delete[] selfRvs;
		selfRvs = NULL;
		P = NULL;
		return max;
	}

	//end self-complementary (3' 4bp) 
	public: int end4SelfCmp () {

		char* P = getPrimer();

		char* rptmp = new char[P_len + 1];
		memset(rptmp, 0, P_len+1);

		for (int n = 0; n < P_len; n++) {
			if (P[n] == 'A') 
				rptmp[n] = 'T';
			else if (P[n] == 'T') 
				rptmp[n] = 'A';
			else if (P[n] == 'C') 
				rptmp[n] = 'G';
			else if (P[n] == 'G') 
				rptmp[n] = 'C';
		}

		char* rvsP2 = reverseChar(rptmp);
		char* rvsP = new char[strlen(rvsP2)-3];
		memset(rvsP,0,strlen(rvsP2)-3);
		strcpy(rvsP,rvsP2+4);

		int max = 0; 
		
		char* end4 = new char[5];
		memset(end4,0,5);
		memcpy(end4,P+(P_len-4),4);
		
		for (int i = 0; i < strlen(rvsP) - 3; i++) {
			int num = 0;
			for (int j = 0; j < strlen(end4); j++) {
				if (rvsP[i + j] == 'A') {
					if (end4[j] == 'T')
						num ++;
				}
				else if (rvsP[i + j] == 'T') {
					if (end4[j] == 'A')
						num ++;
				}
				else if (rvsP[i + j] == 'C') {
					if (end4[j] == 'G')
						num ++;
				}
				else if (rvsP[i + j] == 'G') {
					if (end4[j] == 'C')
						num ++;
				}					
			}
			if (num > max)
				max = num; 
		}
		delete[] rptmp;
		delete[] rvsP;
		delete[] rvsP2;
		delete[] end4;
		P = NULL; rptmp = NULL; rvsP = NULL; rvsP2 = NULL; end4 = NULL;
		return max; 
	}

	// consiguous residue
	public: bool contiguous_resident(int contiguous) {

		char* P = getPrimer();
		
		int contiguous_num = contiguous;
		int count = 1, max = 1;
		if (P_len < contiguous_num) {
			return true;
		} else {
			for (int i = 1; i < P_len; i++) {				
				if (P[i] == P[i - 1]) {
					count ++;
					if (count > max)
						max = count; 
				}
				else 
					count = 1;
				if (max == contiguous_num)
					return false;
			}
		}
		P = NULL;
		return true;
	}

	//hairpin 
	public: int hairpin () {

		char* P = getPrimer();
		int max = 0;
		
		char* pre_s = new char[P_len / 2 + 1];
		memset(pre_s, 0, P_len / 2 + 1);
		char* temp_s = new char[P_len - (P_len / 2) + 1];
		memset(temp_s, 0, P_len - (P_len / 2) + 1);

		memcpy(pre_s, P, P_len / 2);
		strcpy(temp_s, P + P_len / 2);
		char* post_s = reverseChar(temp_s);

		for (int k = 1; k <= strlen(pre_s); k++) {
			int tmpscore = 0;
			for (int i = 0, j = strlen(post_s) - k; i < k; i++, j++) {

				if (pre_s[i] == 'A') {
					if (post_s[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (pre_s[i] == 'T') {
					if (post_s[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (pre_s[i] == 'C') {
					if (post_s[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (pre_s[i] == 'G') {
					if (post_s[j] == 'C')
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

		for (int k = strlen(pre_s) - 1; k > 0; k--) {
			int tmpscore = 0;
			for (int i = strlen(pre_s) - k, j = 0; j < k; i++, j++) {

				if (pre_s[i] == 'A') {
					if (post_s[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (pre_s[i] == 'T') {
					if (post_s[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (pre_s[i] == 'C') {
					if (post_s[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (pre_s[i] == 'G') {
					if (post_s[j] == 'C')
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
		delete[] post_s;
		delete[] temp_s;
		delete[] pre_s;
		P = NULL; post_s = NULL; temp_s = NULL; pre_s = NULL;
		return max;
	}

	// pair filtering 
	// length difference
	public: int LenDiff () {
		int num = absol(P_len - rP_len);
		return num;
	}

	// melting temperature difference 
	public: double TempDiff() {
		char* fP = getPrimer();
		char* rP = getrPrimer();
		double num = absol(double(MeltTemp(fP) - MeltTemp(rP)));
		fP = NULL; rP = NULL;
		return num;
	}

	// product size 
	public: int ProductSize() {		
		int fPos = getfPos();
		int rPos = getrPos();		
		int diff = rPos - fPos + 1; 
		return diff;
	}

	// pair complementary (not contiguous) 
	public: int pairCmp () {
		int max = 0;
		char* P = getPrimer();
		char* rP2 = getrPrimer();
		char* rP;

		rP = reverseChar(rP2);

		int minimum = min(P_len, rP_len);
		for (int k = 1; k <= minimum; k++) {
			int tmpscore = 0;
			for (int i = 0, j = rP_len-k; i < k ; i++, j++) {

				if (P[i] == 'A') {
					if (rP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'T') {
					if (rP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'C') {
					if (rP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'G') {
					if (rP[j] == 'C')
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

		for (int k = minimum -1; k > 0; k--) {
			int tmpscore = 0;
			for (int i = P_len-k, j = 0; j < k && i < k; i++, j++) {

				if (P[i] == 'A') {
					if (rP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'T') {
					if (rP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'C') {
					if (rP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'G') {
					if (rP[j] == 'C')
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
		delete[] rP;
		P = NULL; rP2 = NULL; rP = NULL;
		return max;
	}

	//pair complementary (contiguous) 
	public: int cont_pairCmp (int maxPC) {
		int max = 0;

		char* P = getPrimer();
		char* rP = getrPrimer();

		int minimum = min(P_len, rP_len);

		for (int k = maxPC; k <= minimum-maxPC; k++) {
			int tmpscore = 0;
			for (int i = 0, j = rP_len - k; i < k; i++, j++) {
				if (P[i] == 'A') {
					if (rP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore = 0;
				}
				else if (P[i] == 'T') {
					if (rP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore = 0;
				}
				else if (P[i] == 'C') {
					if (rP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore = 0;
				}
				else if (P[i] == 'G') {
					if (rP[j] == 'C')
						tmpscore ++;
					else 
						tmpscore = 0;
				}					

				if (tmpscore > max)
					max = tmpscore; 
			}
		}

		if (max < maxPC) {			
			for (int k = maxPC; k <= minimum - maxPC; k++) {
				int tmpscore = 0;
				for (int i = P_len - k, j = 0; j < k; i++, j++) {
					if (P[i] == 'A') {
						if (rP[j] == 'T')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					else if (P[i] == 'T') {
						if (rP[j] == 'A')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					else if (P[i] == 'C') {
						if (rP[j] == 'G')
							tmpscore ++;
						else 
							tmpscore = 0;
					}
					else if (P[i] == 'G') {
						if (rP[j] == 'C')
							tmpscore ++;
						else 
							tmpscore = 0;
					}					

					if (tmpscore > max)
						max = tmpscore; 	
				}
			}
		}
		P = NULL; rP = NULL;
		return max;
	}

	// 3' end pair complementary (not contiguous) 
	public: int endPairCmp () {


		int max = 0;
				
		char* P = getPrimer();		
		char* rP = getrPrimer();
		
		int minimum = min(P_len, rP_len);

		int len = P_len/2;

		for (int k = minimum-1; k > len; k--) {
			int tmpscore = 0;
			for (int i = P_len-k, j = 0; j < k && i < k; i++, j++) {
				if (P[i] == 'A') {
					if (rP[j] == 'T')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'T') {
					if (rP[j] == 'A')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'C') {
					if (rP[j] == 'G')
						tmpscore ++;
					else 
						tmpscore --;
				}
				else if (P[i] == 'G') {
					if (rP[j] == 'C')
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
		P = NULL; rP = NULL;
		return max;
	}

	// 3' end pair complementary (3' 4bp)
	public: int end4PairCmp () {
		
		char* P = getPrimer();
		char* rP = getrPrimer();

		char* fend4 = new char[5];
		char* rend4 = new char[5];
		memset(fend4, 0, 5);
		memset(rend4, 0, 5);
		memcpy(fend4, P + (P_len - 4), 4);
		memcpy(rend4, rP + (rP_len - 4), 4);
		
		int max = 0; 
		for (int i = 0; i < rP_len - 3; i++) {
			int num = 0;
			for (int j = 0; j < strlen(fend4); j++) {
				if (rP[i + j] == 'A') {
					if (fend4[j] == 'T')
						num ++;
				}
				else if (rP[i + j] == 'T') {
					if (fend4[j] == 'A')
						num ++;
				}
				else if (rP[i + j] == 'C') {
					if (fend4[j] == 'G')
						num ++;
				}
				else if (rP[i + j] == 'G') {
					if (fend4[j] == 'C')
						num ++;
				}							
			}
			if (num > max)
				max = num; 
		}

		for (int i = 0; i < P_len - 3; i++) {
			int num = 0;
			for (int j = 0; j < strlen(rend4); j++) {
				if (P[i + j] == 'A') {
					if (rend4[j] == 'T')
						num ++;
				}
				else if (P[i + j] == 'T') {
					if (rend4[j] == 'A')
						num ++;
				}
				else if (P[i + j] == 'C') {
					if (rend4[j] == 'G')
						num ++;
				}
				else if (P[i + j] == 'G') {
					if (rend4[j] == 'C')
						num ++;
				}						
			}
			if (num > max)
				max = num; 
		}
		delete[] fend4;
		delete[] rend4;
		P = NULL; rP = NULL; fend4 = NULL; rend4 = NULL;
		return max; 
	} 
};
