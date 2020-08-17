/*
 * Step5Header.hpp
 *
 *  Created on: Nov 8, 2019
 *      Author: jmbae
 */

#include "Input.hpp"
#include "Global.hpp"
#include "DataStructures.hpp"

using namespace std;

typedef unordered_map<char*, set<char*, set_Cmp>, hash<string>, CmpChar> Hmap_set;

long write1(int my_threshold, float* result, long myRank, long passed,
		unsigned int* Ooffset, arraysForStep5 **myArraysStep5,
		Hmap_set *SIDSET, unsigned int *fidx, unsigned int *ridx,
		inputParameter* input,
		FILE *fout){

	int numOfGPUs = input->numOfGPUs;
	int numOfThreads = input->numOfThreads;
	int maxLen = input->maxLen;

	char* tmp_fP = new char[maxLen + 1];
	char* tmp_rP = new char[maxLen + 2];
	char* tmp_fpos = new char[6];
	char* tmp_rpos = new char[6];
	unsigned int index;
	int tmpIdx = myRank % numOfGPUs;
	int tmpRank = myRank / numOfGPUs;
	int tmpnumOfThreads = 0;

	unsigned int *fsid_cnt = myArraysStep5[tmpIdx]->FPoffset;
	unsigned int *rsid_cnt = myArraysStep5[tmpIdx]->RPoffset;
	unsigned char *fP = myArraysStep5[tmpIdx]->FP;
	unsigned char *rP = myArraysStep5[tmpIdx]->RP;
	unsigned int *fpos = myArraysStep5[tmpIdx]->Fpos;
	unsigned int *rpos = myArraysStep5[tmpIdx]->Rpos;

	if (numOfThreads % numOfGPUs == 0) {
		tmpnumOfThreads = numOfThreads / numOfGPUs;
	} else {
		for (int i = 0; i < numOfThreads; i++) {
			if (i % numOfGPUs == tmpIdx)
				tmpnumOfThreads++;
		}
	}

	for (int key = tmpRank; key < my_threshold; key += tmpnumOfThreads) {
		int tmp_rsid_cnt = rsid_cnt[key + 1]
				- rsid_cnt[key];
		for (long i = fsid_cnt[key];
				i < fsid_cnt[key + 1]; i++) {
			for (int j = 0; j < tmp_rsid_cnt; j++) {
				index = Ooffset[key]
						+ (i - fsid_cnt[key]) * tmp_rsid_cnt + j;
				if (result[index] > -1) {
					passed++;
					memset(tmp_fP, 0, maxLen + 1);
					memset(tmp_rP, 0, maxLen + 2);
					tmp_rP[0] = '*';
					for (int l = 0; l < maxLen; l++) {
						if (fP[maxLen * i + l] != '\0')
							tmp_fP[l] = fP[maxLen * i + l];
						else
							break;
					}
					long before = rsid_cnt[key] + j;
					for (int l = 0; l < maxLen; l++) {
						if (rP[before * maxLen + l] != '\0')
							tmp_rP[l + 1] = rP[before * maxLen
									+ l];
						else
							break;
					}
					set<char*, set_Cmp> set1, set2;
					if (SIDSET->find(tmp_fP) != SIDSET->end())
						set1 = SIDSET->find(tmp_fP)->second;
					else
						cout << "ERR_F : " << key << " " << tmp_fP << endl;
					if (SIDSET->find(tmp_rP) != SIDSET->end())
						set2 = SIDSET->find(tmp_rP)->second;
					else
						cout << "ERR_R : " << key << " " << tmp_rP << endl;
					vector<char*> commen_set;
					set_intersection(set1.begin(), set1.end(), set2.begin(),
							set2.end(), back_inserter(commen_set), set_Cmp());

					auto it = commen_set.begin();
					if (it != commen_set.end()) {
						fprintf(fout, "%s", (*it));
						it++;
						for (auto it2 = it; it2 != commen_set.end(); it2++)
							fprintf(fout, "-%s", (*it2));
						fprintf(fout, "\t%s+%s+", tmp_fP, tmp_rP + 1);
					} else {
						if (myRank == 0) {
							cout << "ERR_SID : " << key << " " << tmp_fP << " "
									<< tmp_rP << " ";
							for (auto it2 = set1.begin(); it2 != set1.end();
									it2++)
								cout << *it2 << " ";
							cout << "\t";
							for (auto it2 = set2.begin(); it2 != set2.end();
									it2++)
								cout << *it2 << " ";
							cout << endl;
						}
					}
					if (tmpIdx == 0) {
						int tmp_idx5 = tmpIdx + numOfGPUs * (key + 1);
						fprintf(fout, "%d+%d+%d\t",
								fidx[tmp_idx5] + 1,
								fpos[i],
								rpos[before]);
					} else {
						int tmp_idx5 = myRank + numOfGPUs * key;
						fprintf(fout, "%d+%d+%d\t",
								fidx[tmp_idx5] + 1,
								fpos[i],
								rpos[before]);
					}

					fprintf(fout, "%.15f\n", result[index]);
				}
			}
		}
	}
	return passed;

}

long write2(int i, int working_thread, float* result, long myRank, long passed,
		unsigned int* Ooffset, unsigned int max_off, arraysForStep5 **myArraysStep5,
		Hmap_set *SIDSET, unsigned int *fidx, unsigned int *ridx,
		inputParameter *input,
		FILE *fout, int threshold2_idx) {

	int numOfGPUs = input->numOfGPUs;
	int numOfThreads = input->numOfThreads;
	int maxLen = input->maxLen;

	char* tmp_fP = new char[maxLen + 1];
	char* tmp_rP = new char[maxLen + 2];
	char* tmp_fpos = new char[6];
	char* tmp_rpos = new char[6];
	unsigned int index;

	int tmpIdx = myRank % numOfGPUs;
	int tmpRank = myRank / numOfGPUs;
	int tmpnumOfThreads = 0;

	unsigned int *fsid_cnt = myArraysStep5[tmpIdx]->FPoffset;
	unsigned int *rsid_cnt = myArraysStep5[tmpIdx]->RPoffset;
	unsigned char *fP = myArraysStep5[tmpIdx]->FP;
	unsigned char *rP = myArraysStep5[tmpIdx]->RP;
	unsigned int *fpos = myArraysStep5[tmpIdx]->Fpos;
	unsigned int *rpos = myArraysStep5[tmpIdx]->Rpos;

	if (numOfThreads % numOfGPUs == 0) {
		tmpnumOfThreads = numOfThreads / numOfGPUs;
	} else {
		for (int i = 0; i < numOfThreads; i++) {
			if (i % numOfGPUs == tmpIdx)
				tmpnumOfThreads++;
		}
	}

	for (int Sid_idx = i + tmpRank;
			(Sid_idx < i + working_thread) && (Sid_idx < threshold2_idx);
			Sid_idx += tmpnumOfThreads) {
		int tmp_rsid_cnt = rsid_cnt[Sid_idx + 1]
				- rsid_cnt[Sid_idx];
		for (long j = fsid_cnt[Sid_idx];
				j < fsid_cnt[Sid_idx + 1]; j++) {
			//long result_idx = j - sorted_fsid_cnt[myRank][i];
			for (int k = 0; k < tmp_rsid_cnt; k++) {
				index = Ooffset[Sid_idx - i]
						+ (j - fsid_cnt[Sid_idx]) * tmp_rsid_cnt
						+ k;
				if (index > max_off)
					cout << "ERR!!" << myRank << " " << index << " is over max_off " << max_off
							<< endl;
				if (result[index] > -1) {
					passed++;
					memset(tmp_fP, 0, maxLen + 1);
					memset(tmp_rP, 0, maxLen + 2);
					tmp_rP[0] = '*';
					for (int l = 0; l < maxLen; l++) {
						if (fP[maxLen * j + l] != '\0')
							tmp_fP[l] = fP[maxLen * j + l];
						//                          fprintf(fout, "%c", sorted_fP[myRank][maxLen * j + l]);
						else
							break;
					}
					//                  fprintf(fout, "+");
					long before = rsid_cnt[Sid_idx] + k;
					for (int l = 0; l < maxLen; l++) {
						if (rP[before * maxLen + l] != '\0')
							tmp_rP[l + 1] = rP[before * maxLen
									+ l];
						//                          fprintf(fout, "%c", sorted_rP[myRank][before * maxLen + l]);
						else
							break;
					}
					set<char*, set_Cmp> set1, set2;
					if (SIDSET->find(tmp_fP) != SIDSET->end())
						set1 = SIDSET->find(tmp_fP)->second;
					else
						cout << "ERR_F : " << Sid_idx << " " << tmp_fP << endl;
					if (SIDSET->find(tmp_rP) != SIDSET->end())
						set2 = SIDSET->find(tmp_rP)->second;
					else
						cout << "ERR_R : " << Sid_idx << " " << tmp_rP << endl;
					vector<char*> commen_set;
					set_intersection(set1.begin(), set1.end(), set2.begin(),
							set2.end(), back_inserter(commen_set), set_Cmp());

					auto it = commen_set.begin();
					if (it != commen_set.end()) {
						fprintf(fout, "%s", (*it));
						it++;
						for (auto it2 = it; it2 != commen_set.end(); it2++)
							fprintf(fout, "-%s", (*it2));
						fprintf(fout, "\t%s+%s+", tmp_fP, tmp_rP + 1);
					} else {
						if (myRank == 0) {
							cout << "ERR_SID : " << Sid_idx << " " << tmp_fP
									<< " " << tmp_rP << " ";
							for (auto it2 = set1.begin(); it2 != set1.end();
									it2++)
								cout << *it2 << " ";
							cout << "\t";
							for (auto it2 = set2.begin(); it2 != set2.end();
									it2++)
								cout << *it2 << " ";
							cout << endl;
						}
					}
					if (tmpIdx == 0) {
						int tmp_idx5 = tmpIdx + numOfGPUs * (Sid_idx + 1);
						fprintf(fout, "%d+%d+%d\t",
								fidx[tmp_idx5] + 1,
								fpos[j],
								rpos[before]);
					} else {
						int tmp_idx5 = tmpIdx + numOfGPUs * Sid_idx;
						fprintf(fout, "%d+%d+%d\t",
								fidx[tmp_idx5] + 1,
								fpos[j],
								rpos[before]);
					}

					fprintf(fout, "%.15f\n", result[index]);
				}
			}
		}
	}
	return passed;
}

long write3(int i, int working_thread, float* result, long my_maxrP,
		long myRank, long passed, long fcnt, long rcnt, arraysForStep5 **myArraysStep5,
		Hmap_set *SIDSET, unsigned int *fidx, unsigned int *ridx,
		inputParameter *input, FILE *fout) {

	int numOfGPUs = input->numOfGPUs;
	int numOfThreads = input->numOfThreads;
	int maxLen = input->maxLen;

	char* tmp_fP = new char[maxLen + 1];
	char* tmp_rP = new char[maxLen + 2];
	char* tmp_fpos = new char[6];
	char* tmp_rpos = new char[6];

	int tmpIdx = myRank % numOfGPUs;
	int tmpRank = myRank / numOfGPUs;
	int tmpnumOfThreads = 0;

	unsigned int *fsid_cnt = myArraysStep5[tmpIdx]->FPoffset;
	unsigned int *rsid_cnt = myArraysStep5[tmpIdx]->RPoffset;
	unsigned char *fP = myArraysStep5[tmpIdx]->FP;
	unsigned char *rP = myArraysStep5[tmpIdx]->RP;
	unsigned int *fpos = myArraysStep5[tmpIdx]->Fpos;
	unsigned int *rpos = myArraysStep5[tmpIdx]->Rpos;

	if (numOfThreads % numOfGPUs == 0) {
		tmpnumOfThreads = numOfThreads / numOfGPUs;
	} else {
		for (int i = 0; i < numOfThreads; i++) {
			if (i % numOfGPUs == tmpIdx)
				tmpnumOfThreads++;
		}
	}

	for (long j = tmpRank; j < fcnt; j += tmpnumOfThreads) {
		for (long k = 0; k < rcnt; k++) {
			if (result[j * my_maxrP + k] > -1) {
				passed++;
				memset(tmp_fP, 0, maxLen + 1);
				memset(tmp_rP, 0, maxLen + 2);
				tmp_rP[0] = '*';
				long idx_fsid = fsid_cnt[i] + j;
				for (int l = 0; l < maxLen; l++) {
					if (fP[maxLen * idx_fsid + l] != '\0')
						tmp_fP[l] = fP[maxLen * idx_fsid + l];
//                          fprintf(fout, "%c", fP[maxLen * j + l]);
					else
						break;
				}
//                  fprintf(fout, "+");
				long before = rsid_cnt[i] + k;
				for (int l = 0; l < maxLen; l++) {
					if (rP[before * maxLen + l] != '\0')
						tmp_rP[l + 1] = rP[before * maxLen + l];
//                          fprintf(fout, "%c", rP[before * maxLen + l]);
					else
						break;
				}
				set<char*, set_Cmp> set1, set2;
				if (SIDSET->find(tmp_fP) != SIDSET->end())
					set1 = SIDSET->find(tmp_fP)->second;
				else
					cout << "ERR_F : " << i << " " << tmp_fP << endl;
				if (SIDSET->find(tmp_rP) != SIDSET->end())
					set2 = SIDSET->find(tmp_rP)->second;
				else
					cout << "ERR_R : " << i << " " << tmp_rP << endl;
				vector<char*> commen_set;
				set_intersection(set1.begin(), set1.end(), set2.begin(),
						set2.end(), back_inserter(commen_set), set_Cmp());

				auto it = commen_set.begin();
				if (it != commen_set.end()) {
					fprintf(fout, "%s", (*it));
					it++;
					for (auto it2 = it; it2 != commen_set.end(); it2++)
						fprintf(fout, "-%s", (*it2));
					fprintf(fout, "\t%s+%s+", tmp_fP, tmp_rP + 1);
				} else {
					if (myRank == 0) {
						cout << "ERR_SID : " << i << " " << tmp_fP << " "
								<< tmp_rP << " ";
						for (auto it2 = set1.begin(); it2 != set1.end(); it2++)
							cout << *it2 << " ";
						cout << "\t";
						for (auto it2 = set2.begin(); it2 != set2.end(); it2++)
							cout << *it2 << " ";
						cout << endl;
					}
				}
				if (tmpIdx == 0) {
					int tmp_idx5 = tmpIdx + numOfGPUs * (i + 1);
					fprintf(fout, "%d+%d+%d\t",
							fidx[tmp_idx5] + 1,
							fpos[idx_fsid],
							rpos[before]);
				} else {
					int tmp_idx5 = tmpIdx + numOfGPUs * i;
					fprintf(fout, "%d+%d+%d\t",
							fidx[tmp_idx5] + 1,
							fpos[idx_fsid],
							rpos[before]);
				}
				fprintf(fout, "%.15f\n", result[j * my_maxrP + k]);
			}
		}
	}
	return passed;
}

long fileReadingLine(char* fName) {
	ifstream fin;
	fileReadOpen(&fin, fName, -1);
	long totalLine = countLine(fin);
	fin.close();
	return totalLine;
}

int fileRadingSid(char* fName, int *Sid){
	ifstream fin;
	fileReadOpen(&fin, fName, -1);
	char *buf = new char[300];
	char *ptr, *ptr2;

	if(fin.is_open()){
		fin.seekg(-1, ios_base::end);
	    if(fin.peek() == '\n')
	    {
	      //Start searching for \n occurrences
	    	fin.seekg(-1, std::ios_base::cur);

	      for(int i = fin.tellg();i > 0; i--)
	      {
	        if(fin.peek() == '\n')
	        {
	          //Found
	        	fin.get();
	          break;
	        }
	        //Move one character back
	        fin.seekg(i, std::ios_base::beg);
	      }
	    }
        fin.getline(buf, 300);

        ptr = buf;
        if(strchr(ptr, '\t') != NULL)
        	ptr2 = strchr(ptr, '\t');
        *ptr2 = '\0';
        ptr = ++ptr2;
        if(strchr(ptr, '\t') != NULL)
        	ptr2 = strchr(ptr, '\t');
        *ptr2 = '\0';
        (*Sid) = stoi(ptr);
	}

	del(buf);
	return 0;
}
