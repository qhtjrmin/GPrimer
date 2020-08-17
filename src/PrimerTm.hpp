/*
 * PrimerTm.hpp
 *
 *  Created on: Feb 10, 2020
 *      Author: jmbae
 */

#ifndef PRIMERTM_HPP
#define PRIMERTM_HPP

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



#endif /* PRIMERTM_HPP */
