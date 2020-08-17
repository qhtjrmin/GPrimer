/*
 * GPrimer.hpp
 *
 *  Created on: Feb 15, 2020
 *      Author: jmbae
 */

#ifndef GPRIMER_HPP_
#define GPRIMER_HPP_

#include "Constraints.hpp"
#include "FreeEnergyUtils.hpp"
#include "Input.hpp"
#include "RunPthread.hpp"
#include "Step5Header.hpp"
#include "KernelFunction.hpp"

#include "sys/types.h"
#include "sys/sysinfo.h"

using namespace std;

//#define DEBUG
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

unsigned long usedMemory;

//-------------parameters for filtering---------------
//single filtering
int minLen; int maxLen;
double minGC; double maxGC;
double minTM; double maxTM;
int maxSC;
int endMaxSC;
int endMinDG;
int DGlen;
int maxHP;
int contiguous;

//pair filtering
int lenDiff; int TMDiff;
int minPS; int maxPS;
int maxPC; int endMaxPC;
//-----------------------------------------------------

//type definition
using Hmap_vec = unordered_map<char*, vector<char*>, hash<string>, CmpChar>;
using Hmap_set = unordered_map<char*, set<char*, set_Cmp>, hash<string>, CmpChar>;

//global memory
const int MAX_BUF_SIZE = 300;
int numOfThreads, numOfGPUs, numOfStreams;
vector<char*>* vecData;
const int sidLen = 6;
inputParameter* myInput;
int memGPU; //total memory of GPU

//hash map
unordered_map<char*, bool, hash<string>, CmpChar> primerH;
Hmap_set sidsetH;
Hmap_vec suffixH;
Hmap_vec seedHf, seedHr;
unordered_map<char*, unsigned int, hash<string>, CmpChar> seedHfIdx;
unordered_map<char*, unsigned int, hash<string>, CmpChar> seedHrIdx;
vector<pair<unsigned int, unsigned int>> sortedseedHfIdx;
vector<pair<unsigned int, unsigned int>> sortedseedHrIdx;

//mutex
pthread_mutex_t primerH_lock;
pthread_mutex_t sidsetH_lock;
pthread_mutex_t suffixH_lock;
pthread_mutex_t seedHf_lock;
pthread_mutex_t seedHr_lock;
pthread_mutex_t *seedHf_lock2;
pthread_mutex_t *seedHr_lock2;

//barrier
pthread_barrier_t barrThreads;
pthread_barrier_t barrGPUs;
pthread_barrier_t barr2GPUs;

//global for Step 4
bool probingEnd;
bool rank0End, rank1End;
vector<char>* P1f, *P1r;
vector<unsigned int>* SID1f, *SID1r;
vector<unsigned long>* P1foffset, *P1roffset;
bool* finCheck;
long fPcnt, rPcnt;
long f_threshold, r_threshold;
unsigned int* myWorkloadThresholdF, *myWorkloadThresholdR;
unsigned int* myseedHfcnt, *myseedHrcnt;
unsigned int* myFPcnt, *myRPcnt;
unsigned long* Gf_off, *Gr_off;
char** Gf_key, **Gr_key;
int workloadC1;
arraysForStep4** myArraysC3;
arraysForStep4** myArraysC3_dev;
unsigned int *startPworkload;

//global for Step 5
double writeTime = 0;
float probingTime = 0;
long memSID;
unsigned int* sorted_fidx, *sorted_ridx;
unsigned int* sorted_rsid;
unsigned int* fsid, *rsid;
float **result;
unsigned int **Ooffset;
unsigned int *total_result3;
unsigned int *case3fcnt;
bool *finish;
int sid_workload;
int* fP_idx, *rP_idx;
int* sid_idx;
int* mySmallThreshold;
int* myLargeThreshold;
int smallThreshold = 10; //default decided by heuristic expeirments
int largeThreshold; //decided as distribution
pthread_mutex_t valid_lock;
long total_passed;
FILE** finalFile;
long total_input_workload;
int total_sid;
unsigned long total_input_line;
vector<unsigned int> sorted_fsid;
char *finName, *foutName;
arraysForStep5 *tmpArraysStep5;
arraysForStep5 **myArraysStep5;

int stage4Divide = 1;

void initialization(inputParameter* input); //initialization for Step 2~4
void initialization2(long input_workload); //initialization for Step 5

void* stage1PrimerGeneration(void* rank){
	long myRank = (long) rank;
	long bufSize = 1024 * 5000;
	char buf[bufSize];
	char *sid = new char[100];
	char *sequence = new char[bufSize];
	char *ptr, *ptr2;
	char *Primer = new char[maxLen + 1];
	char *rPrimer = new char[maxLen + 2];
	bool check = true;
	int line = 0;
	char *fName = new char[100];
	string tmpfName;

	ifstream fin;
	fileReadOpen(&fin, myInput->inputPath, -1);

	FILE *fout;
	tmpfName = string(myInput->dirPath) + "/tmpC1.txt";
	strcpy(fName, tmpfName.c_str());
	fileWriteOpen(&fout, fName, myRank);

	while(fin.good()){
		fin.getline(buf, bufSize);

		if(fin.gcount() > 0){

			if (line % numOfThreads == myRank) {

				ptr = buf;
				if(strchr(ptr, '\t') != NULL){
					ptr2 = strchr(ptr, '\t');
					*ptr2 = '\0';
					sid = ptr;
					sequence = ++ptr2;

					int sequenceLen = strlen(sequence);

					for (int i = 0; i < sequenceLen; i++) {
						//Primer, rPrimer initialization
						memset(Primer, 0, maxLen + 1);
						memset(rPrimer, 0, maxLen + 1);
						for (int k = minLen; k <= maxLen; k++) {
							if ((i + k) > sequenceLen) {
								break;
							}
							int t = 0;
							int r = k - 1;
							check = true;
							for (int j = 0; j < k; j++) {
								Primer[t] = sequence[i + j];
								if (sequence[i + r] == 'A')
									rPrimer[t] = 'T';
								else if (sequence[i + r] == 'T')
									rPrimer[t] = 'A';
								else if (sequence[i + r] == 'C')
									rPrimer[t] = 'G';
								else if (sequence[i + r] == 'G')
									rPrimer[t] = 'C';
								else
									check = false;
								t++;
								r--;
							}
							if (check) {
								fprintf(fout, "%s\t%s\t%i\n", Primer, sid, i + 1);
								fprintf(fout, "*%s\t%s\t%lu\n", rPrimer, sid,
										i + strlen(rPrimer));
							}
						}
					}
				}
			}
            line++;
		}
	}
	del(Primer); del(rPrimer);
	fin.close();
	fclose(fout);
	return NULL;
}

void stage1Sort(){
	string sortedFile, command;
	ifstream fin;

	sortedFile = string(myInput->dirPath) + "/sorted.txt";

	//sort by primer
	command = "sort -k1,1 -k2,2n -S " + to_string(myInput->bufferSize) + "% --parallel "
			+ to_string(myInput->numOfThreads) + " -T " + string(myInput->dirPath) + " "
			+ string(myInput->dirPath) + "/tmpC1.txt_*" + " -o " + sortedFile;
	sysCall(command);

	sysCall("rm " + string(myInput->dirPath) + "/tmpC1.txt_*");
}

void stage1FileDistribution(){

	char buf[MAX_BUF_SIZE];
	int count = 0;
	long line = 0;
	char *primer, *etc;
	char *ptr, *ptr2;
	char *fName = new char[100];
	string tmpfName;
	char *before = new char[maxLen + 2];
	memset(before, 0, maxLen + 2);

	tmpfName = string(myInput->dirPath) + "/sorted.txt";
	strcpy(fName, tmpfName.c_str());

	ifstream fin;
	fileReadOpen(&fin, fName, -1);

	ifstream fin2;
	fileReadOpen(&fin2, fName, -1);
    long total_input = countLine(fin2);

	FILE **fout = new FILE *[numOfThreads];
	for(int i = 0; i < numOfThreads; i++)
		fileWriteOpen(&fout[i], myInput->c1Path, i);

    while(fin.good()){

        fin.getline(buf,MAX_BUF_SIZE);

        if(fin.gcount() > 0){
            ptr = buf;
            if(strchr(ptr,'\t') != NULL){
				ptr2 = strchr(ptr, '\t');
				*ptr2 = '\0';
				primer = ptr;
				etc = ++ptr2;

				if (line >= (count + 1) * total_input / numOfThreads + 1) {
					if (strcmp(primer, before)) {
						fclose(fout[count]);
						count++;
					}
				}
				fprintf(fout[count], "%s\t%s\n", primer, etc);

				if (line >= ((count + 1) * total_input / numOfThreads))
					strcpy(before, primer);
				line++;
            }
        }
    }
    del(fout); del(fName); del(before);
}

void* stage2(void* rank) {

	long myRank = (long) rank;
	char buf[MAX_BUF_SIZE];
	char *suffix, *untagPrimer;
	char *primer, *sid;
	char *ptr, *ptr2;
	char **val;
	bool pass = false; //check wheter primer is filtered or not
	long line = 0, cnt = 0;
	int prefixLen = 4;

#ifdef DEBUG
	long count = 0;
#endif

	char* beforePrimer = new char[maxLen + 2];
	memset(beforePrimer, 0, maxLen + 2);

	char* beforeSid = new char[sidLen];
	memset(beforeSid, 0, sidLen);

	val = new char*[3];

	ifstream fin;
	fileReadOpen(&fin, myInput->c1Path, myRank);

	FILE *fout1, *fout2;
	fileWriteOpen(&fout1, myInput->c2Path, myRank); // C2.txt
	fileWriteOpen(&fout2, myInput->c1SidsetPath, myRank); //C1'.txt

	while (fin.good()) {

		fin.getline(buf, MAX_BUF_SIZE);

		if (fin.gcount() > 0) {

			//reading file
			ptr = buf;
			cnt = 0;
			while ((ptr2 = strchr(ptr, '\t')) != NULL) {
				*ptr2 = '\0';
				val[cnt] = ptr;
				ptr = ++ptr2;
				cnt++;
			}
			val[cnt] = ptr;

			if (val[0][0] == '*')
				untagPrimer = val[0] + 1;
			else
				untagPrimer = val[0];

#ifdef DEBUG
			if (myRank == 0) {
				if (line >= count * 10000000) {
					cout << beforePrimer << " " << val[0] << " " << untagPrimer
							<< " " << val[1] << " " << val[2] << " " << endl;
					count++;
				}
			}
#endif

			if (strcmp(val[0], beforePrimer)) { //if primer is changed (primer is sorted)
				FreeEnergyUtils FEU = FreeEnergyUtils(untagPrimer, endMinDG,
						DGlen);
				Constraints Pconstraints = Constraints(untagPrimer);

				//filtering
				if (((Pconstraints.MeltTemp(untagPrimer) >= minTM)
						&& (Pconstraints.MeltTemp(untagPrimer) <= maxTM))
						&& ((Pconstraints.GCcontent() >= minGC)
								&& (Pconstraints.GCcontent() <= maxGC))
						&& Pconstraints.selfCmp() < maxSC
						&& Pconstraints.endSelfCmp() < endMaxSC
						&& Pconstraints.hairpin() < maxHP
						&& FEU.get_free_energy()
						&& Pconstraints.contiguous_resident(contiguous)) { //the case of filtering pass

					primer = new char[maxLen + 2];
					sid = new char[sidLen];

					strcpy(primer, val[0]);
					strcpy(sid, val[1]);

					//primerH update
					pthread_mutex_lock(&primerH_lock);
					primerH.insert(pair<char*, bool>(primer, true));
					pthread_mutex_unlock(&primerH_lock);

					//sidsetH update
					set<char*, set_Cmp> tmp;
					tmp.insert(sid);
					pthread_mutex_lock(&sidsetH_lock);
					if (sidsetH.find(primer) == sidsetH.end()) {
						sidsetH.insert(
								pair<char*, set<char*, set_Cmp>>(primer, tmp));
					}
					pthread_mutex_unlock(&sidsetH_lock);
					freeContainer(tmp);

					//suffixH update
					suffix = new char[maxLen - 2];
					strcpy(suffix, untagPrimer + prefixLen);
					pthread_mutex_lock(&suffixH_lock);
					if (suffixH.find(suffix) == suffixH.end()) {
						vector<char*> tmp;
						tmp.push_back(primer);
						suffixH.insert(
								pair<char*, vector<char*>>(suffix, tmp));
						freeContainer(tmp);
					} else {
						(suffixH.find(suffix)->second).push_back(primer);
					}
					pthread_mutex_unlock(&suffixH_lock);

					fprintf(fout1, "%s\t%s\t%s\n", val[0], val[1], val[2]); //write to C2.txt
					pass = true; //check that the primer is passed
				}

				else
					// case of filtered out in single filtering
					pass = false;

				if (line > 0)
					fprintf(fout2, "\n%s\t%s", val[0], val[1]); //write to C1'
				else
					fprintf(fout2, "%s\t%s", val[0], val[1]); //write to C1'
			}

			else if (!strcmp(val[0], beforePrimer) && pass) {
				//beforePrimer is the same as current primer, and the before primer was pass

				fprintf(fout1, "%s\t%s\t%s\n", val[0], val[1], val[2]); //write to C2

				if (strcmp(beforeSid, val[1]))
					fprintf(fout2, "-%s", val[1]); //write to C1'

				primer = new char[maxLen + 2];
				sid = new char[6];

				strcpy(primer, val[0]);
				strcpy(sid, val[1]);

				//sidsetH update
				pthread_mutex_lock(&sidsetH_lock);
				if (sidsetH.find(primer) != sidsetH.end())
					(sidsetH.find(primer)->second).insert(sid);
				pthread_mutex_unlock(&sidsetH_lock);

			}

			else {
				//before primer is the same as current primer but before primer is filterd out
				if (strcmp(beforeSid, val[1]))
					fprintf(fout2, "-%s", val[1]); //write to C1'
			}
		}
		strcpy(beforePrimer, val[0]);
		strcpy(beforeSid, val[1]);
		line++;
	}
	fprintf(fout2, "\n"); //write to C1'

	del(beforePrimer);
	del(beforeSid);
	del(val);

	fin.close();
	fclose(fout1);
	fclose(fout2);

	return NULL;
}

void* stage3(void* rank) {
	long myRank = (long) rank;
	long line = 0;
	bool check;
	char *ptr, *ptr2;
	char *buf, *suffix, *primer;
	char *c1Primer, *c2Primer;
	char *c1Sidset, *oriSidset;
	int prefixLen = 4;

	set<char*, set_Cmp> c2Sidset;
	vector<char*> c2PrimerSet;
	vector<char*>::iterator it;

#ifdef DEBUG
	long count = 0;
#endif

	buf = new char[300000];
	c1Primer = new char[maxLen + 2];
	c1Sidset = new char[300000];
	oriSidset = new char[300000];

	ifstream fin;
	fileReadOpen(&fin, myInput->c1SidsetPath, myRank);

	//reading file C1'
	while (fin.good()) {
		fin.getline(buf, 300000);
		if (fin.gcount() > 0) {

			ptr = buf;
			ptr2 = strchr(ptr, '\t');
			*ptr2 = '\0';
			strcpy(c1Primer, ptr);
			strcpy(oriSidset, ++ptr2);

			if (c1Primer[0] == '*')
				primer = c1Primer + 1;
			else
				primer = c1Primer;

			suffix = primer + prefixLen;

#ifdef DEBUG
			if ((line >= count * 10000000) && (myRank == 1)) {
				if (myRank == 0) {
					cout << primer << " " << suffix << " " << oriSidset << endl;
					count++;
				}
			}
#endif

			//suffixH probing based on suffix (key)
			if (suffixH.find(suffix) != suffixH.end()) {
				c2PrimerSet = suffixH.find(suffix)->second;
				for (it = c2PrimerSet.begin(); it != c2PrimerSet.end(); it++) { // iteration to primer having the same suffix
					c2Primer = (*it);
					if (primerH.find(c2Primer) != primerH.end()) {
						if ((primerH.find(c2Primer)->second) == true) {

							check = true;
							c2Sidset = sidsetH.find(c2Primer)->second;

							strcpy(c1Sidset, oriSidset);
							ptr = c1Sidset;
							while ((ptr2 = strchr(ptr, '-')) != NULL) {
								*ptr2 = '\0';
								if (c2Sidset.find(ptr) == c2Sidset.end()) { //if c2Sidset does not include the sid (ptr) from c1_sidset
									primerH.find(c2Primer)->second = false; //primrH update
									check = false;
									break;
								}
								ptr = ++ptr2;
							}
							if (check) {
								if (c2Sidset.find(ptr) == c2Sidset.end()) {
									primerH.find(c2Primer)->second = false; //primerH update
								}
							}
							freeContainer(c2Sidset);
						}
					}
				}
			}
		}
		line++;
	}
	del(buf);
	del(c1Sidset);
	del(oriSidset);
	del(c1Primer);

	fin.close();

	return NULL;

}

void stage3Delete(){
	for(auto it = suffixH.begin(); it != suffixH.end(); it++){
		delete[] it->first;
		freeContainer(it->second);
	}
	freeContainer(suffixH);
}

int stage4Check(int k){
	long numOfTrues = 0;
	long memPredicted = 0;
	int myMemory;
	for(auto it = primerH.begin(); it != primerH.end(); it++){
		if((it->second) == true)
			numOfTrues++;
	}

	if(k == 1){
		memPredicted = ((numOfTrues / 1000000 + 1) * 90) / numOfGPUs;
		myMemory = 14;

	}
	else if(k == 2){
		memPredicted = ((numOfTrues / 1000000 + 1) * 130) / numOfGPUs;
		myMemory = 28;
	}
	else {
		cout << "ERR: wrong value of k. k should be one or two" << endl;
		return -1;
	}

	if((memGPU - memPredicted) - float(myMemory / numOfGPUs) > 0){
		stage4Divide = 1;
	}
	else{
		for(int i = 2; i < 100; i++){
			if((memGPU - (memPredicted / i) - float(myMemory/numOfGPUs)) > 0 ){
				stage4Divide = i;
				i = 100;
			}
		}
	}
	for(int i = 0; i < numOfThreads; i++)
		startPworkload[i] = 0;

	//cout << "divide: " << stage4Divide << " " << memPredicted << "/" << memGPU << endl;

	return 0;

}

void* stage4Building(void* param) {
	paramThread *writeInput = (paramThread *) param;
	long myRank = writeInput->rank;
	int k = writeInput->myParam;

	char *primer, *seed, *newKey, *untagPrimer;
	int seedLen, primerLen, keyLen;
	char *primerLenChar = new char[3];

	unsigned int primerHworkload = primerH.size() / numOfThreads + 1;
	unsigned int primerHworkload2;
	if(stage4Divide == 1)
		primerHworkload2 = primerHworkload;
	else
		primerHworkload2 = primerHworkload / stage4Divide + 1;

	Hmap_vec* seedH;
	pthread_mutex_t* seedH_lock;

#ifdef DEBUG
	long cnt = 0;
#endif

	if (k == 1)
		seedLen = minLen / 2;
	else if (k == 2)
		seedLen = minLen / 3;

	seed = new char[seedLen + 1];
	keyLen = seedLen + 6;

	unsigned int valid_cnt = 0;

	for (auto it_primerH = next(primerH.begin(), primerHworkload * myRank + startPworkload[myRank]);
			it_primerH != primerH.end(); it_primerH++) {
		if (valid_cnt < primerHworkload2) {
			if (it_primerH->second == true) {
				primer = new char[maxLen + 2];
				strcpy(primer, it_primerH->first);
				memset(seed, 0, seedLen + 1);

				if (primer[0] == '*')
					untagPrimer = primer + 1;
				else
					untagPrimer = primer;

				primerLen = strlen(untagPrimer);
				strcpy(primerLenChar, to_string(primerLen).c_str());

				for (int i = 0; i < primerLen - seedLen + 1; i++) {
					newKey = new char[keyLen];

					memcpy(seed, untagPrimer + i, seedLen);

					// make key(seed + i + len(primer))
					strcpy(newKey, seed);
					strcat(newKey, "+");
					strcat(newKey, to_string(i).c_str());
					strcat(newKey, "+");
					strcat(newKey, primerLenChar);


					#ifdef DEBUG
					if ((valid_cnt >= cnt * 100000) && (myRank == 0)) {
						cout << valid_cnt << " / " << primerHworkload << " : "
								<< primer << " " << newKey << endl;
						cnt++;
					}
					#endif

					if (primer[0] == '*') {
						seedH = &seedHr;
						seedH_lock = &seedHr_lock;
					} else {
						seedH = &seedHf;
						seedH_lock = &seedHf_lock;
					}

					//G update
					pthread_mutex_lock(&(*seedH_lock));
					if (seedH->find(newKey) == seedH->end()) {
						vector<char*> tmp;
						tmp.push_back(primer);
						seedH->insert(
								pair<char*, vector<char*>>(newKey, tmp));
						freeContainer(tmp);
					} else {
						((seedH->find(newKey))->second).push_back(primer);
					}
					pthread_mutex_unlock(&(*seedH_lock));
					i += (seedLen - 1);
				}
				vecData[myRank].push_back(primer);
			}
			valid_cnt++;
		} else
			break;
	}

	startPworkload[myRank] += primerHworkload2;
	del(seed);
	del(primerLenChar);
	return NULL;
}

void* stage4Prepare(void* param) { //prepare arraysC3 and copy it to GPU memory
	paramThread *writeInput = (paramThread *) param;
	long myRank = writeInput->rank;
	int k = writeInput->myParam;
	int seedLen, keyLen;

	unsigned int index = 0;
	unsigned long Pcnt = 0;
	unsigned int size;
	long myPFsize, myPRsize;
	long myseedHfsize, myseedHrsize;
	long mySIDFsize, mySIDRsize;
	double workloadThresholdF;
	double workloadThresholdR;

	long seedHfsize = seedHf.size();
	long seedHrsize = seedHr.size();

	int myMemory;

	if (k == 1) {
		seedLen = minLen / 2;
		workloadThresholdF = (double) ((double) 17 / 20 * (double) seedHfsize);
		workloadThresholdR = (double) ((double) 17 / 20 * (double) seedHrsize);
		myMemory = 14;
	} else if (k == 2) {
		seedLen = minLen / 3;
		workloadThresholdF = (double) ((double) 1 / 5 * (double) seedHfsize);
		workloadThresholdR = (double) ((double) 1 / 5 * (double) seedHrsize);
		myMemory = 28;
	}

	keyLen = seedLen + 6;

	char* key;
	char* key2 = new char[keyLen];
	char* ptr, *ptr2;

	pthread_barrier_wait(&barrThreads);

	if (myRank == 0) {

		Gf_off = new unsigned long[seedHfsize + 1];
		Gf_off[0] = 0;
		Gf_key = new char*[seedHfsize];
		for (Hmap_vec::iterator g_it = seedHf.begin(); g_it != seedHf.end();
				g_it++) {
			key = g_it->first;
			seedHfIdx.insert(pair<char*, unsigned long>(key, index));
			Gf_key[index] = key;
			size = (g_it->second).size(); // the number of primers in one key (group)
			sortedseedHfIdx.push_back(pair<int, unsigned long>(size, index));
			Pcnt += size;
			Gf_off[index + 1] = Pcnt;
			index++;
		}
		fPcnt = Pcnt; // the total number of primers in G_f

		sort(sortedseedHfIdx.begin(), sortedseedHfIdx.end()); // sort by the size

		cout << "minimum_seedHf: " << (*sortedseedHfIdx.begin()).first
				<< " maximum_seedHf: " << (*(--sortedseedHfIdx.end())).first << endl;

		myWorkloadThresholdF = new unsigned int[numOfGPUs];
		myseedHfcnt = new unsigned int[numOfGPUs];
		myFPcnt = new unsigned int[numOfGPUs];

		P1f = new vector<char> [seedHfsize];
		SID1f = new vector<unsigned int> [seedHfsize];
		P1foffset = new vector<unsigned long> [seedHfsize];
		seedHf_lock2 = new pthread_mutex_t[seedHfsize];
		for (long i = 0; i < seedHfsize; i++)
			pthread_mutex_init(&seedHf_lock2[i], NULL);

		probingEnd = false;
		rank0End = false;
		finCheck = new bool[numOfThreads];
	}
	else if(myRank == 1){
		index = 0;
		Pcnt = 0;
		Gr_off = new unsigned long[seedHrsize + 1];
		Gr_off[0] = 0;
		Gr_key = new char*[seedHrsize];

		for (Hmap_vec::iterator g_it = seedHr.begin(); g_it != seedHr.end();
				g_it++) {
			key = g_it->first;
			seedHrIdx.insert(pair<char*, unsigned long>(key, index));
			Gr_key[index] = key;
			size = (g_it->second).size();
			sortedseedHrIdx.push_back(pair<int, unsigned long>(size, index));
			Pcnt += size;
			Gr_off[index + 1] = Pcnt;
			index++;
		}
		rPcnt = Pcnt;

		sort(sortedseedHrIdx.begin(), sortedseedHrIdx.end());

		cout << "minimum_seedHr: " << (*sortedseedHrIdx.begin()).first
				<< " maximum_seedHr: " << (*(--sortedseedHrIdx.end())).first << endl;

		myWorkloadThresholdR = new unsigned int[numOfGPUs];
		myseedHrcnt = new unsigned int[numOfGPUs];
		myRPcnt = new unsigned int[numOfGPUs];

		P1r = new vector<char> [seedHrsize];
		SID1r = new vector<unsigned int> [seedHrsize];
		P1roffset = new vector<unsigned long> [seedHrsize];
		seedHr_lock2 = new pthread_mutex_t[seedHrsize];
		for (long i = 0; i < seedHrsize; i++)
			pthread_mutex_init(&seedHr_lock2[i], NULL);

	}

	pthread_barrier_wait(&barrThreads);

	if(myRank < 2 * numOfGPUs){
		int tmpRank = myRank / 2;
		if(myRank % 2 == 0){ //preparing arraysC3 for FPs
			myPFsize = fPcnt / numOfGPUs * 1.5;

			if (seedHfsize % numOfGPUs == 0)
				myseedHfsize = seedHfsize / numOfGPUs;
			else
				myseedHfsize = seedHfsize / numOfGPUs + 1;

			//malloc arraysC3 (forward primers)
			myArraysC3[tmpRank]->PFoffset = new unsigned long[myseedHfsize + 1];
			myArraysC3[tmpRank]->PF = new char[maxLen * myPFsize];
			myArraysC3[tmpRank]->SIDFoffset = new unsigned long[myPFsize + 1];
			myArraysC3[tmpRank]->SIDF = new unsigned int[3 * myPFsize];
			myArraysC3[tmpRank]->SEEDF = new unsigned int[myseedHfsize];
			myArraysC3[tmpRank]->SIDFoffset[0] = 0;
			myArraysC3[tmpRank]->PFoffset[0] = 0;
			myArraysC3[tmpRank]->outputF = new unsigned int[myPFsize];
			for (unsigned long i = 0; i < myPFsize; i++)
				myArraysC3[tmpRank]->outputF[i] = 0;

			index = 0;
			unsigned long sidIdx = 0;
			set<char*, set_Cmp> sidset;
			set<char*, set_Cmp>::iterator sid_it;
			Hmap_vec::iterator g_it;
			bool first_check = true;
			long left = 0;
			vector<char*> g_vec;
			vector<char*>::iterator it;
			char* primer;

			for (unsigned long i = tmpRank; i < seedHfsize; i += numOfGPUs) {
				int idx = sortedseedHfIdx[i].second;

				if (index > myPFsize || left > myseedHfsize || sidIdx > 2 * myPFsize) {
					cout << "ERRF!! " << index << "/" << myPFsize << " " << sidIdx
							<< "/" << 2 * myPFsize << " " << left << "/" << seedHfsize
							<< endl;
				}

				myArraysC3[tmpRank]->PFoffset[left + 1]
				                                        = myArraysC3[tmpRank]->PFoffset[left]
				                          + (Gf_off[idx + 1] - Gf_off[idx]);

				key = Gf_key[idx];
				g_it = seedHf.find(key);
				strcpy(key2, key);
				ptr = key2;
				ptr2 = strchr(ptr, '+');
				*ptr2 = '\0';
				ptr = ++ptr2;
				ptr2 = strchr(ptr, '+');
				*ptr2 = '\0';
				myArraysC3[tmpRank]->SEEDF[left] = stoi(ptr);

				if (i > workloadThresholdF && first_check) {
					f_threshold = left;
					myWorkloadThresholdF[tmpRank] = f_threshold;
					first_check = false;
				}

				g_vec = g_it->second;
				for (it = g_vec.begin(); it != g_vec.end(); it++) {
					primer = (*it);
					sidset = sidsetH.find(primer)->second;

					long old_idx_sid = sidIdx;
					for (sid_it = sidset.begin(); sid_it != sidset.end();
							sid_it++) {
						myArraysC3[tmpRank]->SIDF[sidIdx] = stoi(*sid_it);
						sidIdx++;
					}
					thrust::sort(myArraysC3[tmpRank]->SIDF + old_idx_sid, myArraysC3[tmpRank]->SIDF + sidIdx);

					myArraysC3[tmpRank]->SIDFoffset[index + 1]
					                     = myArraysC3[tmpRank]->SIDFoffset[index] + sidset.size();

					int plen = strlen(primer);
					for (int k = 0; k < maxLen; k++) {
						if (k < plen)
							myArraysC3[tmpRank]->PF[index * maxLen + k] = primer[k];
						else
							myArraysC3[tmpRank]->PF[index * maxLen + k] = '\0';
					}
					index++;
					freeContainer(sidset);
				}
				freeContainer(g_vec);
				left++;
			}
			myseedHfsize = left;
			mySIDFsize = sidIdx;
			myPFsize = index;
			myFPcnt[tmpRank] = myPFsize;
			myseedHfcnt[tmpRank] = myseedHfsize;

		}

		else{//preparing arraysC3 for RPs
			myPRsize = rPcnt / numOfGPUs * 1.5;

			if (seedHrsize % numOfGPUs == 0)
				myseedHrsize = seedHrsize / numOfGPUs;
			else
				myseedHrsize = seedHrsize / numOfGPUs + 1;

			//malloc arraysC3 (reverse primers)
			myArraysC3[tmpRank]->PRoffset = new unsigned long[myseedHrsize + 1];
			myArraysC3[tmpRank]->PR = new char[maxLen * myPRsize];
			myArraysC3[tmpRank]->SIDRoffset = new unsigned long[myPRsize + 1];
			myArraysC3[tmpRank]->SIDR = new unsigned int[3 * myPRsize];
			myArraysC3[tmpRank]->SEEDR = new unsigned int[myseedHrsize];
			myArraysC3[tmpRank]->SIDRoffset[0] = 0;
			myArraysC3[tmpRank]->PRoffset[0] = 0;
			myArraysC3[tmpRank]->outputR = new unsigned int[myPRsize];
			for (unsigned long i = 0; i < myPRsize; i++)
				myArraysC3[tmpRank]->outputR[i] = 0;

			index = 0;
			unsigned long sidIdx = 0;
			set<char*, set_Cmp> sidset;
			set<char*, set_Cmp>::iterator sid_it;
			Hmap_vec::iterator g_it;
			bool first_check = true;
			long left = 0;
			vector<char*>g_vec;
			vector<char*>::iterator it;
			char* primer;

			for (unsigned long i = tmpRank; i < seedHrsize; i += numOfGPUs) {
				int idx = sortedseedHrIdx[i].second;

				if (index > myPRsize || left > myseedHrsize || sidIdx > 2 * myPRsize) {
					cout << "ERRR!! " << index << "/" << myPRsize << " " << sidIdx
							<< "/" << 2 * myPRsize << " " << left << "/" << myseedHrsize
							<< endl;
				}

				myArraysC3[tmpRank]->PRoffset[left + 1] = myArraysC3[tmpRank]->PRoffset[left]
										  + (Gr_off[idx + 1] - Gr_off[idx]);

				key = Gr_key[idx];
				g_it = seedHr.find(key);
				strcpy(key2, key);
				ptr = key2;
				ptr2 = strchr(ptr, '+');
				*ptr2 = '\0';
				ptr = ++ptr2;
				ptr2 = strchr(ptr, '+');
				*ptr2 = '\0';
				myArraysC3[tmpRank]->SEEDR[left] = stoi(ptr);

				if (i > workloadThresholdR && first_check) {
					r_threshold = left;
					myWorkloadThresholdR[tmpRank] = r_threshold;
					first_check = false;
				}

				g_vec = g_it->second;
				for (it = g_vec.begin(); it != g_vec.end(); it++) {
					primer = (*it);
					sidset = sidsetH.find(primer)->second;

					long old_idx_sid = sidIdx;
					for (sid_it = sidset.begin(); sid_it != sidset.end();
							sid_it++) {
						myArraysC3[tmpRank]->SIDR[sidIdx] = stoi(*sid_it);
						sidIdx++;
					}
					thrust::sort(myArraysC3[tmpRank]->SIDR + old_idx_sid, myArraysC3[tmpRank]->SIDR + sidIdx);

					myArraysC3[tmpRank]->SIDRoffset[index + 1]
										 = myArraysC3[tmpRank]->SIDRoffset[index] + sidset.size();

					primer = primer + 1;
					int plen = strlen(primer);
					for (int k = 0; k < maxLen; k++) {
						if (k < plen)
							myArraysC3[tmpRank]->PR[index * maxLen + k] = primer[k];
						else
							myArraysC3[tmpRank]->PR[index * maxLen + k] = '\0';
					}
					index++;
					freeContainer(sidset);
				}
				freeContainer(g_vec);
				left++;
			}
			myseedHrsize = left;
			mySIDRsize = sidIdx;
			myPRsize = index;
			myRPcnt[tmpRank] = myPRsize;
			myseedHrcnt[tmpRank] = myseedHrsize;
		}
	}

	pthread_barrier_wait(&barrThreads);

	if (myRank < numOfGPUs * 2) {

		int tmpRank = myRank / 2;

		if (myRank % 2 == 0) {
			cudaSetDevice(tmpRank);

			//malloc arraysC3 for FPs
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->PFoffset,
							sizeof(unsigned long) * (myseedHfsize + 1))); //P3_offset
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->PF,
							sizeof(char) * maxLen * myPFsize)); //P3
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->SIDFoffset,
							sizeof(unsigned long) * (myPFsize + 1))); //SID3_offset
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->SIDF,
							sizeof(unsigned int) * mySIDFsize)); //SID3
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->SEEDF,
							sizeof(unsigned int) * myseedHfsize)); //SEED
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->outputF,
							sizeof(unsigned int) * myPFsize)); //output;
		}

		else {

			cudaSetDevice(tmpRank);
			//malloc arraysC3 for RPs
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->PRoffset,
							sizeof(unsigned long) * (myseedHrsize + 1))); //P3_offset
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->PR,
							sizeof(char) * maxLen * myPRsize)); //P3
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->SIDRoffset,
							sizeof(unsigned long) * (myPRsize + 1))); //SID3_offset
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->SIDR,
							sizeof(unsigned int) * mySIDRsize)); //SID3
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->SEEDR,
							sizeof(unsigned int) * myseedHrsize)); //SEED
			HANDLE_ERROR(
					cudaMalloc((void** ) &myArraysC3_dev[tmpRank]->outputR,
							sizeof(unsigned int) * myPRsize)); //output
		}

		if (myRank % 2 == 0) {
			cudaSetDevice(tmpRank);
			//memcpy arraysC3 for FPs
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->PFoffset, &myArraysC3[tmpRank]->PFoffset[0],
							sizeof(unsigned long) * (myseedHfsize + 1),
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->PF, &myArraysC3[tmpRank]->PF[0],
							sizeof(char) * maxLen * myPFsize,
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->SIDFoffset, &myArraysC3[tmpRank]->SIDFoffset[0],
							sizeof(unsigned long) * (myPFsize + 1),
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->SIDF, &myArraysC3[tmpRank]->SIDF[0],
							sizeof(unsigned int) * mySIDFsize,
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->SEEDF, &(myArraysC3[tmpRank]->SEEDF[0]),
							sizeof(unsigned int) * myseedHfsize,
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->outputF, &myArraysC3[tmpRank]->outputF[0],
							sizeof(unsigned int) * myPFsize,
							cudaMemcpyHostToDevice));
		} else {
			cudaSetDevice(tmpRank);
			//memcpy arraysC3 for RPs
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->PRoffset, &myArraysC3[tmpRank]->PRoffset[0],
							sizeof(unsigned long) * (myseedHrsize + 1),
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->PR, &myArraysC3[tmpRank]->PR[0],
							sizeof(char) * maxLen * myPRsize,
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->SIDRoffset, &myArraysC3[tmpRank]->SIDRoffset[0],
							sizeof(unsigned long) * (myPRsize + 1),
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->SIDR, &myArraysC3[tmpRank]->SIDR[0],
							sizeof(unsigned int) * mySIDRsize,
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->SEEDR, &myArraysC3[tmpRank]->SEEDR[0],
							sizeof(unsigned int) * myseedHrsize,
							cudaMemcpyHostToDevice));
			HANDLE_ERROR(
					cudaMemcpyAsync(myArraysC3_dev[tmpRank]->outputR, &myArraysC3[tmpRank]->outputR[0],
							sizeof(unsigned int) * myPRsize,
							cudaMemcpyHostToDevice));
		}
	}

	pthread_barrier_wait(&barrThreads);

	if(myRank == 0){
		usedMemory += sizeof(unsigned long) * (myseedHfsize + 1);
		usedMemory += sizeof(char) * maxLen * myPFsize;
		usedMemory += sizeof(unsigned long) * (myPFsize + 1);
		usedMemory += sizeof(unsigned int) * mySIDFsize;
		usedMemory += sizeof(unsigned int) * myseedHfsize;
		usedMemory += sizeof(unsigned int) * myPFsize;
	}

	pthread_barrier_wait(&barrThreads);

	if(myRank == 1){

		usedMemory += sizeof(unsigned long) * (myseedHrsize + 1);
		usedMemory += sizeof(char) * maxLen * myPRsize;
		usedMemory += sizeof(unsigned long) * (myPRsize + 1);
		usedMemory += sizeof(unsigned int) * mySIDRsize;
		usedMemory += sizeof(unsigned int) * myseedHrsize;
		usedMemory += sizeof(unsigned int) * myPRsize;

	}

	pthread_barrier_wait(&barrThreads);

	del(key2);

	if (myRank == 0) {
		del(Gf_off);
		del(Gf_key);
		for(int i = 0; i < numOfGPUs; i++){
			delete[] myArraysC3[i]->PFoffset;
			delete[] myArraysC3[i]->SIDFoffset;
			delete[] myArraysC3[i]->SIDF;
			delete[] myArraysC3[i]->SEEDF;
		}

		usedMemory = usedMemory / 1024;
		usedMemory = usedMemory / 1024;

		int remaining_mm;

		remaining_mm = memGPU - usedMemory;

		workloadC1 =  float(remaining_mm / float(myMemory / numOfGPUs)) * 10000;

		cout << usedMemory << " MB is utilized. "
				<< workloadC1 <<  " workload" << endl;

	} else if (myRank == 1) {
		del(Gr_off);
		del(Gr_key);
		for(int i = 0; i < numOfGPUs; i++){
			delete[] myArraysC3[i]->PRoffset;
			delete[] myArraysC3[i]->SIDRoffset;
			delete[] myArraysC3[i]->SIDR;
			delete[] myArraysC3[i]->SEEDR;
		}

	}



	return NULL;
}

void* stage4Probing(void* param) { //prepare arraysC1 and perform GPU operation
	paramThread *writeInput = (paramThread *) param;
	long myRank = writeInput->rank;
	int k = writeInput->myParam;

	bool f_check, g_check, mCheck;
	char *ptr, *ptr2;
	char buf[300000];
	char *primer, *c1_primer, *c1_sidset, *ori_sidset;
	long line = 0;
	char *seed;
	char *new_key;
	int seedLen;

	if (k == 1)
		seedLen = minLen / 2;
	else if (k == 2)
		seedLen = minLen / 3;
	else {
		cout << "ERR: wrong value of k. k should be one or two" << endl;
		return NULL;
	}

	vector<char*>::iterator it;

	c1_primer = new char[maxLen + 2];
	c1_sidset = new char[300000];
	ori_sidset = new char[300000];

	char* len_primer_char = new char[3];

	seed = new char[seedLen + 1];
	memset(seed, 0, seedLen + 1);
	new_key = new char[seedLen + 6];

	vector<char>* P;
	vector<unsigned int>* sid;
	vector<unsigned long>*sid_cnt;
	pthread_mutex_t* seedH_lock;
	unordered_map<char*, unsigned int, hash<string>, CmpChar>* G_idx;

	int tmp_idx;
	unsigned int index;
	unsigned long idx_sid;
	int loop = 0;
	finCheck[myRank] = false;

	unsigned long *fP_total = new unsigned long[numOfStreams];
	unsigned long *rP_total = new unsigned long[numOfStreams];
	long *fsid_total = new long[numOfStreams];
	long *rsid_total = new long[numOfStreams];
	int *Gf_size = new int[numOfStreams];
	int *Gr_size = new int[numOfStreams];

	char **P1f_dev = new char*[numOfStreams];
	char **P1r_dev = new char*[numOfStreams];
	unsigned int** SID1f_dev = new unsigned int*[numOfStreams];
	unsigned int **SID1r_dev = new unsigned int*[numOfStreams];
	unsigned long** P1foffset_dev = new unsigned long*[numOfStreams];
	unsigned long **P1roffset_dev = new unsigned long*[numOfStreams];
	unsigned long** c1_foff_dev = new unsigned long*[numOfStreams];
	unsigned long **c1_roff_dev = new unsigned long*[numOfStreams];

	cudaStream_t* streamF = new cudaStream_t[numOfStreams];
	cudaStream_t* streamR = new cudaStream_t[numOfStreams];

	char** myP1f = new char*[numOfStreams];
	unsigned int** mySID1f = new unsigned int*[numOfStreams];
	unsigned long** mySID1foffset = new unsigned long *[numOfStreams];
	unsigned long** myP1foffset = new unsigned long*[numOfStreams];

	char** myP1r = new char*[numOfStreams];
	unsigned int** mySID1r = new unsigned int*[numOfStreams];
	unsigned long** mySID1roffset = new unsigned long *[numOfStreams];
	unsigned long** myP1roffset = new unsigned long*[numOfStreams];

#ifdef DEBUG
	int cnt = 0;
#endif

	for (int i = 0; i < numOfStreams; i++) {
		myP1foffset[i] = new unsigned long[seedHf.size() + 1];
		myP1roffset[i] = new unsigned long[seedHr.size() + 1];
	}

	int myStream = 0;

	ifstream fin;
	fileReadOpen(&fin, myInput->c1SidsetPath, myRank);

	while (1) {

		if (line < workloadC1) {
			if (!finCheck[myRank])
				fin.getline(buf, 300000);

			if (fin.eof()) {
				finCheck[myRank] = true;
				line = workloadC1;
			}
		}

		if (line < workloadC1) {
			ptr = buf;
			ptr2 = strchr(ptr, '\t');
			*ptr2 = '\0';

			strcpy(c1_primer, ptr);
			strcpy(ori_sidset, ++ptr2);

			if (c1_primer[0] == '*') {
				primer = c1_primer + 1;
				f_check = false;
			} else {
				primer = c1_primer;
				f_check = true;
			}

#ifdef DEBUG
			if (myRank == 1 && line >= cnt * 1000000) {
				cout << loop << " " << line << endl;
				cnt++;
			}
#endif

			int len_primer = strlen(primer);
			strcpy(len_primer_char, to_string(len_primer).c_str());
			memset(seed, 0, seedLen + 1);
			for (int i = 0; i < len_primer - seedLen + 1; i++) {
				memcpy(seed, primer + i, seedLen);

				//make key(seed + i + len(primer)
				strcpy(new_key, seed);
				strcat(new_key, "+");
				strcat(new_key, to_string(i).c_str());
				strcat(new_key, "+");
				strcat(new_key, len_primer_char);

				if (!f_check) {
					G_idx = &seedHrIdx;
				} else {
					G_idx = &seedHfIdx;
				}

				if ((G_idx->find(new_key)) != G_idx->end()) {
					g_check = true;
					tmp_idx = (G_idx->find(new_key))->second;
				} else
					g_check = false;

				int total;
				if (g_check) {
					if (!f_check) {
						seedH_lock = &seedHr_lock2[tmp_idx];
						total = seedHr.size();
					} else {
						seedH_lock = &seedHf_lock2[tmp_idx];
						total = seedHf.size();
					}

					if (tmp_idx < total) {
						pthread_mutex_lock(seedH_lock);
						if (!f_check) {
							P = &P1r[tmp_idx];
							sid_cnt = &P1roffset[tmp_idx];
							sid = &SID1r[tmp_idx];
						} else {
							P = &P1f[tmp_idx];
							sid_cnt = &P1foffset[tmp_idx];
							sid = &SID1f[tmp_idx];
						}

						int plen = strlen(primer);
						for (int j = 0; j < maxLen; j++) {
							if (j < plen)
								P->push_back(primer[j]);
							else
								P->push_back('\0');
						}

						long sidset_cnt = 1;
						strcpy(c1_sidset, ori_sidset);
						ptr = c1_sidset;
						while ((ptr2 = strchr(ptr, '-')) != NULL) {
							*ptr2 = '\0';
							sid->push_back(stoi(ptr));
							sidset_cnt++;
							ptr = ++ptr2;
						}
						sid->push_back(stoi(ptr));
						sid_cnt->push_back(sidset_cnt);
						pthread_mutex_unlock(seedH_lock);
					} else
						cout << "ERR!!" << tmp_idx << " " << total << endl;
				}
				i += seedLen - 1;
			}
			line++;
		}

		else {

			mCheck = false;

			if (myRank == 0)
				rank0End = false;

			if (myRank == 1)
				rank1End = false;

			pthread_barrier_wait(&barrThreads);

			//change!
			if (myRank < numOfGPUs * 2) {

				int tmpRank = myRank / 2;

				if (myRank % 2 == 0) {

					memset(myP1foffset[myStream], 0, seedHf.size() + 1);
					fP_total[myStream] = 0;
					fsid_total[myStream] = 0;

					int off_idx = 0;
					myP1foffset[myStream][0] = 0;
					for (int i = tmpRank; i < sortedseedHfIdx.size(); i +=
							numOfGPUs) {
						int idx = sortedseedHfIdx[i].second;
						fP_total[myStream] += P1foffset[idx].size();
						fsid_total[myStream] += SID1f[idx].size();
						myP1foffset[myStream][off_idx + 1] = fP_total[myStream];
						off_idx++;
					}

					myP1f[myStream] = new char[maxLen * fP_total[myStream]];
					mySID1foffset[myStream] =
							new unsigned long[fP_total[myStream] + 1];
					mySID1f[myStream] = new unsigned int[fsid_total[myStream]];

					memset(myP1f[myStream], 0, maxLen * fP_total[myStream]);
					memset(mySID1foffset[myStream], 0, fP_total[myStream] + 1);
					memset(mySID1f[myStream], 0, fsid_total[myStream]);
					mySID1foffset[myStream][0] = 0;

					off_idx = 0;
					index = 0;
					idx_sid = 0;
					long old_sid_cnt = 0;
					for (int i = tmpRank; i < sortedseedHfIdx.size(); i +=
							numOfGPUs) {
						int idx = sortedseedHfIdx[i].second;
						int total_cnt = myP1foffset[myStream][off_idx + 1]
								- myP1foffset[myStream][off_idx];
						old_sid_cnt = 0;

						for (int j = 0; j < total_cnt; j++) {
							for (int k = 0; k < maxLen; k++) {
								myP1f[myStream][index * maxLen + k] =
										P1f[idx][j * maxLen + k];
							}
							for (long k = old_sid_cnt;
									k < old_sid_cnt + P1foffset[idx][j];
									k++) {
								mySID1f[myStream][idx_sid + k] =
										SID1f[idx][k];
							}
							old_sid_cnt += P1foffset[idx][j];
							mySID1foffset[myStream][index + 1] =
									mySID1foffset[myStream][index]
											+ P1foffset[idx][j];
							index++;
						}
						idx_sid += old_sid_cnt;
						off_idx++;
					}
					Gf_size[myStream] = off_idx;
				}

				else {
					memset(myP1roffset[myStream], 0, seedHr.size() + 1);
					rP_total[myStream] = 0;
					rsid_total[myStream] = 0;

					int off_idx = 0;
					myP1roffset[myStream][0] = 0;
					for (int i = tmpRank; i < sortedseedHrIdx.size(); i +=
							numOfGPUs) {
						int idx = sortedseedHrIdx[i].second;
						rP_total[myStream] += P1roffset[idx].size();
						rsid_total[myStream] += SID1r[idx].size();
						myP1roffset[myStream][off_idx + 1] = rP_total[myStream];
						off_idx++;
					}

					myP1r[myStream] = new char[maxLen * rP_total[myStream]];
					mySID1roffset[myStream] =
							new unsigned long[rP_total[myStream] + 1];
					mySID1r[myStream] = new unsigned int[rsid_total[myStream]];

					memset(myP1r[myStream], 0, maxLen * rP_total[myStream]);
					memset(mySID1roffset[myStream], 0, rP_total[myStream] + 1);
					memset(mySID1r[myStream], 0, rsid_total[myStream]);
					mySID1roffset[myStream][0] = 0;

					off_idx = 0;
					index = 0;
					idx_sid = 0;
					long old_sid_cnt = 0;
					for (int i = tmpRank; i < sortedseedHrIdx.size(); i +=
							numOfGPUs) {
						int idx = sortedseedHrIdx[i].second;
						int total_cnt = myP1roffset[myStream][off_idx + 1]
								- myP1roffset[myStream][off_idx];
						old_sid_cnt = 0;

						for (int j = 0; j < total_cnt; j++) {
							for (int k = 0; k < maxLen; k++) {
								myP1r[myStream][index * maxLen + k] =
										P1r[idx][j * maxLen + k];
							}
							for (long k = old_sid_cnt;
									k < old_sid_cnt + P1roffset[idx][j];
									k++) {
								mySID1r[myStream][idx_sid + k] =
										SID1r[idx][k];
							}
							old_sid_cnt += P1roffset[idx][j];
							mySID1roffset[myStream][index + 1] =
									mySID1roffset[myStream][index]
											+ P1roffset[idx][j];
							index++;
						}
						idx_sid += old_sid_cnt;
						off_idx++;
					}

					Gr_size[myStream] = off_idx;

				}

				if (myStream > 0) {
					int tmpStream = myStream - 1;

					if(myRank % 2 == 0){
						cudaSetDevice(tmpRank);
						cudaStreamSynchronize(streamF[tmpStream]);
						cudaStreamDestroy(streamF[tmpStream]);

						del(myP1f[tmpStream]);
						del(mySID1f[tmpStream]);
						del(mySID1foffset[tmpStream]);

						cudaSetDevice(tmpRank);

						HANDLE_ERROR(cudaFree(P1f_dev[tmpStream]));
						HANDLE_ERROR(cudaFree(SID1f_dev[tmpStream]));
						HANDLE_ERROR(cudaFree(P1foffset_dev[tmpStream]));
						HANDLE_ERROR(cudaFree(c1_foff_dev[tmpStream]));
					}

					else{
						cudaSetDevice(tmpRank);

						cudaStreamSynchronize(streamR[tmpStream]);
						cudaStreamDestroy(streamR[tmpStream]);

						del(myP1r[tmpStream]);
						del(mySID1r[tmpStream]);
						del(mySID1roffset[tmpStream]);

						cudaSetDevice(tmpRank);

						HANDLE_ERROR(cudaFree(P1r_dev[tmpStream]));
						HANDLE_ERROR(cudaFree(SID1r_dev[tmpStream]));
						HANDLE_ERROR(cudaFree(P1roffset_dev[tmpStream]));
						HANDLE_ERROR(cudaFree(c1_roff_dev[tmpStream]));
					}
				}

				int thread = 1024;
				int block_f, block_r;
				if (myWorkloadThresholdF[tmpRank] / thread > 1024)
					block_f = 1024;
				else
					block_f = myWorkloadThresholdF[tmpRank] / thread + 1;

				if (myWorkloadThresholdR[tmpRank] / thread > 1024)
					block_r = 1024;
				else
					block_r = myWorkloadThresholdR[tmpRank] / thread + 1;

				mCheck = true;

				if (myRank % 2 == 0) {

					cudaSetDevice(tmpRank);

					HANDLE_ERROR(cudaStreamCreate(&streamF[myStream]));

					HANDLE_ERROR(
							cudaMalloc((void** ) &P1f_dev[myStream],
									sizeof(char) * maxLen
											* fP_total[myStream]));
					HANDLE_ERROR(
							cudaMalloc((void** ) &SID1f_dev[myStream],
									sizeof(unsigned int)
											* fsid_total[myStream]));
					HANDLE_ERROR(
							cudaMalloc((void** ) &P1foffset_dev[myStream],
									sizeof(unsigned long)
											* (fP_total[myStream] + 1)));
					HANDLE_ERROR(
							cudaMalloc((void** ) &c1_foff_dev[myStream],
									sizeof(unsigned long)
											* (Gf_size[myStream] + 1)));

					HANDLE_ERROR(
							cudaMemcpyAsync(P1f_dev[myStream],
									&myP1f[myStream][0],
									sizeof(char) * maxLen * fP_total[myStream],
									cudaMemcpyHostToDevice, streamF[myStream]));
					HANDLE_ERROR(
							cudaMemcpyAsync(SID1f_dev[myStream],
									&mySID1f[myStream][0],
									sizeof(unsigned int) * fsid_total[myStream],
									cudaMemcpyHostToDevice, streamF[myStream]));
					HANDLE_ERROR(
							cudaMemcpyAsync(P1foffset_dev[myStream],
									&mySID1foffset[myStream][0],
									sizeof(unsigned long)
											* (fP_total[myStream] + 1),
									cudaMemcpyHostToDevice, streamF[myStream]));
					HANDLE_ERROR(
							cudaMemcpyAsync(c1_foff_dev[myStream],
									&myP1foffset[myStream][0],
									sizeof(unsigned long)
											* (Gf_size[myStream] + 1),
									cudaMemcpyHostToDevice, streamF[myStream]));

					stage4_probing<<<block_f, thread, 0, streamF[myStream]>>>(k,
							myWorkloadThresholdF[tmpRank], seedLen, myArraysC3_dev[tmpRank]->SEEDF,
							c1_foff_dev[myStream], P1f_dev[myStream],
							SID1f_dev[myStream], P1foffset_dev[myStream],
							myArraysC3_dev[tmpRank]->PFoffset, myArraysC3_dev[tmpRank]->PF,
							myArraysC3_dev[tmpRank]->SIDF, myArraysC3_dev[tmpRank]->SIDFoffset,
							myArraysC3_dev[tmpRank]->outputF, maxLen);

					stage4_probing2<<<1024, thread, 0, streamF[myStream]>>>(k,
							myWorkloadThresholdF[tmpRank], myseedHfcnt[tmpRank], seedLen,
							myArraysC3_dev[tmpRank]->SEEDF, c1_foff_dev[myStream],
							P1f_dev[myStream], SID1f_dev[myStream],
							P1foffset_dev[myStream], myArraysC3_dev[tmpRank]->PFoffset,
							myArraysC3_dev[tmpRank]->PF, myArraysC3_dev[tmpRank]->SIDF,
							myArraysC3_dev[tmpRank]->SIDFoffset, myArraysC3_dev[tmpRank]->outputF, maxLen);

				}

				else {

					cudaSetDevice(tmpRank);

					HANDLE_ERROR(cudaStreamCreate(&streamR[myStream]));

					HANDLE_ERROR(
							cudaMalloc((void** ) &P1r_dev[myStream],
									sizeof(char) * maxLen
											* rP_total[myStream]));
					HANDLE_ERROR(
							cudaMalloc((void** ) &SID1r_dev[myStream],
									sizeof(unsigned int)
											* rsid_total[myStream]));
					HANDLE_ERROR(
							cudaMalloc((void** ) &P1roffset_dev[myStream],
									sizeof(unsigned long)
											* (rP_total[myStream] + 1)));
					HANDLE_ERROR(
							cudaMalloc((void** ) &c1_roff_dev[myStream],
									sizeof(unsigned long)
											* (Gr_size[myStream] + 1)));

					HANDLE_ERROR(
							cudaMemcpyAsync(P1r_dev[myStream],
									&myP1r[myStream][0],
									sizeof(char) * maxLen * rP_total[myStream],
									cudaMemcpyHostToDevice, streamR[myStream]));
					HANDLE_ERROR(
							cudaMemcpyAsync(SID1r_dev[myStream],
									&mySID1r[myStream][0],
									sizeof(unsigned int) * rsid_total[myStream],
									cudaMemcpyHostToDevice, streamR[myStream]));
					HANDLE_ERROR(
							cudaMemcpyAsync(P1roffset_dev[myStream],
									&mySID1roffset[myStream][0],
									sizeof(unsigned long)
											* (rP_total[myStream] + 1),
									cudaMemcpyHostToDevice, streamR[myStream]));
					HANDLE_ERROR(
							cudaMemcpyAsync(c1_roff_dev[myStream],
									&myP1roffset[myStream][0],
									sizeof(unsigned long)
											* (Gr_size[myStream] + 1),
									cudaMemcpyHostToDevice, streamR[myStream]));

					stage4_probing<<<block_r, thread, 0, streamR[myStream]>>>(k,
							myWorkloadThresholdR[tmpRank], seedLen, myArraysC3_dev[tmpRank]->SEEDR,
							c1_roff_dev[myStream], P1r_dev[myStream],
							SID1r_dev[myStream], P1roffset_dev[myStream],
							myArraysC3_dev[tmpRank]->PRoffset, myArraysC3_dev[tmpRank]->PR,
							myArraysC3_dev[tmpRank]->SIDR, myArraysC3_dev[tmpRank]->SIDRoffset,
							myArraysC3_dev[tmpRank]->outputR, maxLen);

					stage4_probing2<<<1024, thread, 0, streamR[myStream]>>>(k,
							myWorkloadThresholdR[tmpRank], myseedHrcnt[tmpRank], seedLen,
							myArraysC3_dev[tmpRank]->SEEDR, c1_roff_dev[myStream],
							P1r_dev[myStream], SID1r_dev[myStream],
							P1roffset_dev[myStream], myArraysC3_dev[tmpRank]->PRoffset,
							myArraysC3_dev[tmpRank]->PR, myArraysC3_dev[tmpRank]->SIDR,
							myArraysC3_dev[tmpRank]->SIDRoffset, myArraysC3_dev[tmpRank]->outputR, maxLen);

				}

				pthread_barrier_wait(&barr2GPUs);

				if (myRank == 0) {
					for (int i = 0; i < seedHf.size(); i++) {
						freeContainer(P1f[i]);
						freeContainer(SID1f[i]);
						freeContainer(P1foffset[i]);
					}
				}

				if (myRank == 0) {
					for (int i = 0; i < seedHr.size(); i++) {
						freeContainer(P1r[i]);
						freeContainer(SID1r[i]);
						freeContainer(P1roffset[i]);
					}
					rank0End = true;
				}
			}//(if myRank < 2 * numOfGPUs)

			if (myRank == 1)
				rank1End = true;

			pthread_barrier_wait(&barrThreads);

			if (myRank == 0) {
				bool end = true;
				for (int i = 0; i < numOfThreads; i++) {
					if (finCheck[i] == false) {
						end = false;
						break;
					}
				}
				if (end) {
					probingEnd = true;
					rank0End = true;
				}
				rank0End = true;
			}

			pthread_barrier_wait(&barrThreads);

			if (probingEnd && rank0End && rank1End) {

				if (myRank < 2 * numOfGPUs) {

					int tmpRank = myRank / 2;

					cudaDeviceSynchronize();

					if(myRank % 2 == 0){

						cudaSetDevice(tmpRank);

						delete[] myP1f[myStream];
						delete[] mySID1f[myStream];
						delete[] mySID1foffset[myStream];

						HANDLE_ERROR(cudaStreamDestroy(streamF[myStream]));

						if(mCheck){
							cudaSetDevice(tmpRank);
							HANDLE_ERROR(cudaFree(P1f_dev[myStream]));
							HANDLE_ERROR(cudaFree(SID1f_dev[myStream]));
							HANDLE_ERROR(cudaFree(P1foffset_dev[myStream]));
							HANDLE_ERROR(cudaFree(c1_foff_dev[myStream]));
						}
					}

					else{

						cudaSetDevice(tmpRank);

						delete[] myP1r[myStream];
						delete[] mySID1r[myStream];
						delete[] mySID1roffset[myStream];

						HANDLE_ERROR(cudaStreamDestroy(streamR[myStream]));

						if(mCheck){
							cudaSetDevice(tmpRank);
							HANDLE_ERROR(cudaFree(P1r_dev[myStream]));
							HANDLE_ERROR(cudaFree(SID1r_dev[myStream]));
							HANDLE_ERROR(cudaFree(P1roffset_dev[myStream]));
							HANDLE_ERROR(cudaFree(c1_roff_dev[myStream]));
						}
					}
				}

				cudaDeviceSynchronize();

				for(int i = 0; i < numOfStreams; i++){
					del(myP1foffset[i]);
					del(myP1roffset[i]);
				}

				del(c1_primer);
				del(len_primer_char);
				del(new_key); del(seed);
				del(c1_sidset);
				del(ori_sidset);
				del(fP_total); del(rP_total);
				del(fsid_total); del(rsid_total);
				del(Gf_size); del(Gr_size);
				del(myP1foffset); del(myP1roffset);
				del(myP1f); del(myP1r);
				del(mySID1f); del(mySID1r);
				del(mySID1foffset); del(mySID1roffset);
				del(P1f_dev); del(P1r_dev);
				del(SID1f_dev); del(SID1r_dev);
				del(P1foffset_dev); del(P1roffset_dev);
				del(c1_foff_dev); del(c1_roff_dev);
				del(streamF); del(streamR);

				fin.close();
				return NULL;
			}

			line = 0;
			loop++;
			myStream++;

#ifdef DEBUG
			cnt = 0;
#endif

		}
	}
}

void* stage4Update(void* rank) { //update the validation of primerH
	long myRank = (long) rank;

	if (myRank < numOfGPUs * 2) {
		int tmpRank = myRank / 2;
		char* newP;
		if (myRank % 2 == 0) {

			newP = new char[maxLen + 1];
			newP[maxLen] = '\0';

			cudaSetDevice(tmpRank);

			cudaMemcpy(myArraysC3[tmpRank]->outputF, myArraysC3_dev[tmpRank]->outputF,
					sizeof(unsigned int) * myFPcnt[tmpRank], cudaMemcpyDeviceToHost);

			for (int j = 0; j < myFPcnt[tmpRank]; j++) {
				if (myArraysC3[tmpRank]->outputF[j] > 0) {
					memset(newP, 0, maxLen + 1);
					for (int k = 0; k < maxLen; k++) {
						newP[k] = myArraysC3[tmpRank]->PF[j * maxLen + k];
					}

					pthread_mutex_lock(&primerH_lock);
					if (primerH.find(newP)->second == true)
						primerH.find(newP)->second = false;
					pthread_mutex_unlock(&primerH_lock);
				}
			}
			cudaFree(myArraysC3_dev[tmpRank]->PF);
			cudaFree(myArraysC3_dev[tmpRank]->SIDF);
			cudaFree(myArraysC3_dev[tmpRank]->SIDFoffset);
			cudaFree(myArraysC3_dev[tmpRank]->PFoffset);
			cudaFree(myArraysC3_dev[tmpRank]->SEEDF);
			cudaFree(myArraysC3_dev[tmpRank]->outputF);
			delete[] myArraysC3[tmpRank]->outputF;
			delete[] myArraysC3[tmpRank]->PF;

		}

		else {
			newP = new char[maxLen + 2];
			newP[maxLen + 1] = '\0';
			newP[0] = '*';

			cudaSetDevice(tmpRank);

			cudaMemcpy(myArraysC3[tmpRank]->outputR, myArraysC3_dev[tmpRank]->outputR,
					sizeof(unsigned int) * myRPcnt[tmpRank], cudaMemcpyDeviceToHost);

			for (int j = 0; j < myRPcnt[tmpRank]; j++) {
				if (myArraysC3[tmpRank]->outputR[j] > 0) {
					memset(newP + 1, 0, maxLen + 1);
					for (int k = 0; k < maxLen; k++) {
						newP[k + 1] = myArraysC3[tmpRank]->PR[j * maxLen + k];
					}
					pthread_mutex_lock(&primerH_lock);
					if (primerH.find(newP)->second == true)
						primerH.find(newP)->second = false;
					pthread_mutex_unlock(&primerH_lock);
				}
			}
			cudaFree(myArraysC3_dev[tmpRank]->PR);
			cudaFree(myArraysC3_dev[tmpRank]->SIDR);
			cudaFree(myArraysC3_dev[tmpRank]->SIDRoffset);
			cudaFree(myArraysC3_dev[tmpRank]->PRoffset);
			cudaFree(myArraysC3_dev[tmpRank]->SEEDR);
			cudaFree(myArraysC3_dev[tmpRank]->outputR);
			delete[] myArraysC3[tmpRank]->outputR;
			delete[] myArraysC3[tmpRank]->PR;
		}

		del(newP);
	}

	pthread_barrier_wait(&barrThreads);

	if (myRank == 0) {

		del(P1f);
		del(SID1f);
		del(P1foffset);

		for (auto it = seedHf.begin(); it != seedHf.end(); it++) {
			delete[] it->first;
			freeContainer(it->second);
		}

		freeContainer(seedHf);
		freeContainer(seedHfIdx);
		freeContainer(sortedseedHfIdx);

		del(seedHf_lock2);
		del(finCheck);
		del(myFPcnt);
		del(myWorkloadThresholdF);
		del(myseedHfcnt);

	} else if (myRank == 1) {


		del(P1r);
		del(SID1r);
		del(P1roffset);

		for (auto it = seedHr.begin(); it != seedHr.end(); it++) {
			delete[] it->first;
			freeContainer(it->second);
		}

		freeContainer(seedHr);
		freeContainer(seedHrIdx);
		freeContainer(sortedseedHrIdx);

		del(seedHr_lock2);
		del(myRPcnt);
		del(myWorkloadThresholdR);
		del(myseedHrcnt);
	}

	for(auto it = vecData[myRank].begin(); it != vecData[myRank].end(); it++){
		delete[] (*it);
	}
	freeContainer(vecData[myRank]);

	return NULL;
}

void stage4Final(){// cat total result file for sorting
	del(vecData);
	del(myArraysC3);
	del(myArraysC3_dev);
}

void stage5Sort(){
	string sortedFile, command;
	ifstream fin;

	sortedFile = string(myInput->dirPath) + "/input_for_step5.txt";
	command = "sort -k2,2n -k1,1 -S " + to_string(myInput->bufferSize) + "% --parallel "
			+ to_string(myInput->numOfThreads) + " -T " + string(myInput->dirPath)
			+ " " + string(myInput->c4Path2) + "_*" + " -o " + sortedFile;
	sysCall(command);

	finName = new char[100];
	strcpy(finName, sortedFile.c_str());
	foutName = myInput->outputPath;

	fileRadingSid(finName, &total_sid);

	total_input_line = fileReadingLine(finName);
	long input_workload = total_input_line / numOfGPUs * 3;
	initialization2(input_workload);
	cout << "total sid: " << total_sid << " "<< "total line: " << total_input_line << endl;
	fin.close();
}

void *stage5Prepare(void *rank) {
	long myRank = (long) rank;

	if (myRank < 2) {
		ifstream fin;
		fileReadOpen(&fin, finName, -1);

		char buf[MAX_BUF_SIZE];
		long idx_rP = 0, idx_fP = 0;
		long idx_rpos = 0, idx_fpos = 0;
		char* ptr, *ptr2, *ptr3;
		char* tmp_primer;
		int tmp_sid, tmp_pos;
		int old_sid = 0;

		if (myRank == 0) {
			fsid[0] = 0;
			sorted_fsid.resize(total_sid + 1);
			for (int i = 0; i < total_sid + 1; i++) {
				sorted_fsid[i] = 0;
			}
		} else {
			rsid[0] = 0;
			for (int i = 0; i < total_sid + 1; i++) {
				sorted_rsid[i] = 0;
			}
		}

		bool myTurn = false;
		while (fin.good()) {
			fin.getline(buf, MAX_BUF_SIZE);
			myTurn = false;
			if (fin.gcount() > 0) {
				ptr = buf;
				ptr2 = strchr(ptr, '\t');
				*ptr2 = '\0';
				tmp_primer = ptr;

				if ((tmp_primer[0] != '*' && myRank == 0)
						|| (tmp_primer[0] == '*' && myRank == 1)) {
					myTurn = true;
				}

				if (myTurn) {
					//tmp_primer2 = ptr;
					ptr3 = strchr(++ptr2, '\t');
					*ptr3 = '\0';
					tmp_sid = stoi(ptr2);
					//tmp_sid2 = ptr2;
					tmp_pos = stoi(++ptr3);

					if (old_sid != tmp_sid) {
						if (tmp_sid != (old_sid + 1)) {
							if (myRank == 0) {
								for (int i = old_sid + 1; i < tmp_sid; i++)
									fsid[i] = fsid[old_sid];
							} else {
								for (int i = old_sid + 1; i < tmp_sid; i++)
									rsid[i] = rsid[old_sid];
							}
						}
						if (myRank == 0)
							fsid[tmp_sid] = fsid[old_sid];
						else
							rsid[tmp_sid] = rsid[old_sid];
						old_sid = tmp_sid;
					}

					if (myRank == 0) {
						fsid[tmp_sid]++; // sid_count++
						sorted_fsid[tmp_sid]++;

						PrimerTm temp = PrimerTm(tmp_primer);
						tmpArraysStep5->Ftemp[idx_fpos] = temp.getTM();
						FreeEnergyUtils deltag = FreeEnergyUtils(tmp_primer,
								endMinDG, DGlen);
						tmpArraysStep5->Fenergy[idx_fpos] = deltag.get_score();
						tmpArraysStep5->Fpos[idx_fpos] = tmp_pos;
						idx_fpos++;

						for (int i = 0; i < maxLen; i++) {
							if (tmp_primer[i] != '\0') {
								tmpArraysStep5->FP[idx_fP] = tmp_primer[i];
								idx_fP++;
							} else {
								for (int j = i; j < maxLen; j++) {
									tmpArraysStep5->FP[idx_fP] = '\0';
									idx_fP++;
								}
								break;
							}
						}
					} else {
						tmp_primer = ptr + 1;
						rsid[tmp_sid]++; // sid_count++ in reverse
						sorted_rsid[tmp_sid]++;

						PrimerTm temp = PrimerTm(tmp_primer);
						tmpArraysStep5->Rtemp[idx_rpos] = temp.getTM();
						FreeEnergyUtils deltag = FreeEnergyUtils(tmp_primer,
								endMinDG, DGlen);
						tmpArraysStep5->Renergy[idx_rpos] = deltag.get_score();
						tmpArraysStep5->Rpos[idx_rpos] = tmp_pos;
						idx_rpos++;

						for (int i = 0; i < maxLen; i++) {
							if (tmp_primer[i] != '\0') {
								tmpArraysStep5->RP[idx_rP] = tmp_primer[i];
								idx_rP++;
							} else {
								for (int j = i; j < maxLen; j++) {
									tmpArraysStep5->RP[idx_rP] = '\0';
									idx_rP++;
								}
								break;
							}
						}
					}
				}
			}
		}
		fin.close();
	}

	pthread_barrier_wait(&barrThreads);

	if (myRank == 0) {

		sort(sorted_fsid.begin(), sorted_fsid.end());

		int distribution[10000];

		for (int i = 0; i < 10000; i++)
			distribution[i] = 0;

		for (int i = 0; i < total_sid; i++) {
			if (sorted_fsid[i] < 10000)
				distribution[sorted_fsid[i]]++;
		}

		long dist_cnt = 0;
		for (int i = 0; i < 10; i++) {
			dist_cnt = 0;
			for (int j = 0; j < 1000; j++) {
				if (distribution[i * 1000 + j] > 0) {
					dist_cnt++;
				}
			}

			if (dist_cnt < 100) {
				largeThreshold = i * 1000;
				break;
			}
		}

		cout << "largeThreshold: " << largeThreshold << endl;

		bool* sid_check = new bool[total_sid];
		for (int i = 0; i < total_sid; i++)
			sid_check[i] = true;

		for (int i = 1; i < total_sid + 1; i++) {
			for (int j = 1; j < total_sid + 1; j++) {
				if ((sorted_fsid[i] == fsid[j] - fsid[j - 1])
						&& (sid_check[j - 1])) {
					sid_check[j - 1] = false;
					sorted_fidx[i] = j - 1;
					sorted_ridx[i] = sorted_fidx[i];
					break;
				}
			}
		}

		for (int i = 1; i < total_sid + 1; i++)
			sorted_rsid[i] = rsid[sorted_ridx[i] + 1] - rsid[sorted_ridx[i]];

		if (total_sid % numOfGPUs == 0)
			sid_workload = total_sid / numOfGPUs;
		else
			sid_workload = total_sid / numOfGPUs + 1;

		for (int i = 0; i < numOfGPUs; i++) {
			myArraysStep5[i]->FPoffset[0] = 0;
			myArraysStep5[i]->RPoffset[0] = 0;
			mySmallThreshold[i] = 0;
			myLargeThreshold[i] = 0;
		}
		int sid_cnt = 1;

		for (int i = 1; i < total_sid + 1; i++) {
			int idx = i % numOfGPUs;
			myArraysStep5[idx]->FPoffset[sid_cnt] =
					myArraysStep5[idx]->FPoffset[sid_cnt - 1]
					+ sorted_fsid[i];
			myArraysStep5[idx]->RPoffset[sid_cnt] =
					myArraysStep5[idx]->RPoffset[sid_cnt - 1]
					+ sorted_rsid[i];
			if ((sorted_fsid[i] > smallThreshold) && (mySmallThreshold[idx] == 0)) {
				mySmallThreshold[idx] = sid_cnt;
				cout << "threshold " << smallThreshold << " : " << i << " " << sorted_fsid[i]
						<< " " << sorted_rsid[i]  << endl;
			}
			if ((sorted_fsid[i] > largeThreshold) && (myLargeThreshold[idx] == 0)) {
				myLargeThreshold[idx] = sid_cnt;
				cout << "threshold " << largeThreshold << " : " << i << " " << sorted_fsid[i]
						<< " " << sorted_rsid[i] << endl;
			}
			sorted_fsid[i] += sorted_fsid[i - 1];
			sorted_rsid[i] += sorted_rsid[i - 1];
			sid_idx[idx] = sid_cnt;
			if (idx == 0)
				sid_cnt++;
		}
	}

	pthread_barrier_wait(&barrThreads);

	long idx_fpos = 0;
	long idx_rpos = 0;
	long idx_fP = 0;
	long idx_rP = 0;
	int cnt = 0;

	int tmpRank = myRank / 2;

	if (myRank < 2 * numOfGPUs) {
		if (myRank % 2 == 0) {
			for (int i = 1; i < total_sid + 1; i++) {
				if (i % numOfGPUs == tmpRank) {
					cnt = sorted_fsid[i] - sorted_fsid[i - 1];
					for (int j = 0; j < cnt; j++) {
						myArraysStep5[tmpRank]->Ftemp[idx_fpos] =
								tmpArraysStep5->Ftemp[fsid[sorted_fidx[i]] + j];
						myArraysStep5[tmpRank]->Fenergy[idx_fpos] =
								tmpArraysStep5->Fenergy[fsid[sorted_fidx[i]] + j];
						myArraysStep5[tmpRank]->Fpos[idx_fpos] =
								tmpArraysStep5->Fpos[fsid[sorted_fidx[i]] + j];
						idx_fpos++;
						for (int k = 0; k < maxLen; k++) {
							myArraysStep5[tmpRank]->FP[idx_fP] =
									tmpArraysStep5->FP[(fsid[sorted_fidx[i]] + j) * maxLen + k];
							idx_fP++;
						}
					}
					fP_idx[tmpRank] = idx_fpos;
				}
			}
		} else {
			for (int i = 1; i < total_sid + 1; i++) {
				if (i % numOfGPUs == tmpRank) {
					cnt = sorted_rsid[i] - sorted_rsid[i - 1];
					for (int j = 0; j < cnt; j++) {
						myArraysStep5[tmpRank]->Rpos[idx_rpos] =
								tmpArraysStep5->Rpos[rsid[sorted_ridx[i]] + j];
						myArraysStep5[tmpRank]->Rtemp[idx_rpos] =
								tmpArraysStep5->Rtemp[rsid[sorted_ridx[i]] + j]; //idx_rpos++;
						myArraysStep5[tmpRank]->Renergy[idx_rpos] =
								tmpArraysStep5->Renergy[rsid[sorted_ridx[i]] + j];
						idx_rpos++;
						for (int k = 0; k < maxLen; k++) {
							myArraysStep5[tmpRank]->RP[idx_rP] =
									tmpArraysStep5->RP[(rsid[sorted_ridx[i]] + j) * maxLen + k];
							idx_rP++;
						}
					}
					rP_idx[tmpRank] = idx_rpos;
				}
			}
		}
	}

	pthread_barrier_wait(&barrThreads);

	if (myRank == 0) {
		//delete
		del(tmpArraysStep5->FP);
		del(tmpArraysStep5->RP);
		del(tmpArraysStep5->Fpos);
		del(tmpArraysStep5->Rpos);
		del(tmpArraysStep5->Ftemp);
		del(tmpArraysStep5->Rtemp);
		del(tmpArraysStep5->Fenergy);
		del(tmpArraysStep5->Renergy);

		delete[] fsid;
		delete[] rsid;

	}
	return NULL;
}

void *stage5(void* rank) {
	long myRank = (long) rank;
	int tmpIdx = myRank % numOfGPUs;
	unsigned long my_threshold = mySmallThreshold[tmpIdx];
	int my_total_sid = sid_idx[tmpIdx];
	unsigned int total_result;
	unsigned char *fP_dev, *rP_dev;
	unsigned int *fpos_dev, *rpos_dev;
	unsigned int *fsid_dev, *rsid_dev;
	double *ftemp_dev, *rtemp_dev;
	double *fenergy_dev, *renergy_dev;
	float *result_dev;
	unsigned int *Ooffset_dev;

	int thread = 512;
	int block = my_threshold / thread + 1;

	int working_thread = thread * block;
	int old_fcnt = 0, old_rcnt = 0;
	int fcnt, rcnt;
	int tmp_fcnt, tmp_rcnt;

	double st1, et1, st2, et2;

	long passed = 0;

	inputParameter *myInput_dev;

#ifdef DEBUG
	unsigned long used_memory = 0;
#endif

	if (myRank < numOfGPUs) {

		cudaSetDevice(myRank);

		HANDLE_ERROR(cudaMalloc((void**) &fsid_dev,
				sizeof(unsigned int) * (my_total_sid + 1)));
		HANDLE_ERROR(cudaMalloc((void**) &rsid_dev,
				sizeof(unsigned int) * (my_total_sid + 1)));

		HANDLE_ERROR(cudaMemcpy(fsid_dev, &myArraysStep5[myRank]->FPoffset[0],
				sizeof(unsigned int) * (my_total_sid + 1),
				cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(rsid_dev, &myArraysStep5[myRank]->RPoffset[0],
				sizeof(unsigned int) * (my_total_sid + 1),
				cudaMemcpyHostToDevice));

		HANDLE_ERROR(cudaMalloc((void**) &myInput_dev, sizeof(inputParameter)));
		HANDLE_ERROR(cudaMemcpy(myInput_dev, &myInput[0],
				sizeof(inputParameter), cudaMemcpyHostToDevice));

#ifdef DEBUG
		used_memory = 0;
		used_memory += sizeof(unsigned int) * (my_total_sid + 1);
		used_memory += sizeof(unsigned int) * (my_total_sid + 1);
		memSID = used_memory;

		if (myRank == 0)
		cout << "SID memory: " << (used_memory / 1024) / 1024 << " MB"
		<< endl;
#endif

		st1 = getCurrentTime();
		st2 = getCurrentTime();

		//-----------------------------------CASE1---------------------------------------------

		Ooffset[myRank] = new unsigned int[my_threshold + 1];
		Ooffset[myRank][0] = 0;

		for (int i = 0; i < my_threshold; i++) {
			tmp_fcnt = myArraysStep5[myRank]->FPoffset[i + 1]
					- myArraysStep5[myRank]->FPoffset[i];
			tmp_rcnt = myArraysStep5[myRank]->RPoffset[i + 1]
					- myArraysStep5[myRank]->RPoffset[i];
			Ooffset[myRank][i + 1] = Ooffset[myRank][i] + tmp_fcnt * tmp_rcnt;
		}

		for (int i = 0; i < my_threshold;) {

			if ((i + working_thread) < my_threshold) {
				fcnt = myArraysStep5[myRank]->FPoffset[i + working_thread]
						- myArraysStep5[myRank]->FPoffset[i];
				rcnt = myArraysStep5[myRank]->RPoffset[i + working_thread]
						- myArraysStep5[myRank]->RPoffset[i];
			} else {
				//cout << "case2 : " << total_sid << endl;
				fcnt = myArraysStep5[myRank]->FPoffset[my_threshold]
						- myArraysStep5[myRank]->FPoffset[i];
				rcnt = myArraysStep5[myRank]->RPoffset[my_threshold]
						- myArraysStep5[myRank]->RPoffset[i];
			}

			total_result = Ooffset[myRank][my_threshold];
			result[myRank] = new float[total_result];

			memset(result[myRank], -1.0, sizeof(float) * total_result);

#ifdef DEBUF
			cout << i << " cnt: " << fcnt << " " << rcnt
					<< " " << total_result << endl;
#endif

			cudaSetDevice(myRank);
			HANDLE_ERROR(cudaMalloc((void**) &fP_dev, sizeof(unsigned char) * maxLen * fcnt));
			HANDLE_ERROR(cudaMalloc((void**) &rP_dev, sizeof(unsigned char) * maxLen * rcnt));

			HANDLE_ERROR(cudaMalloc((void**) &fpos_dev, sizeof(unsigned int) * fcnt));
			HANDLE_ERROR(cudaMalloc((void**) &rpos_dev, sizeof(unsigned int) * rcnt));

			HANDLE_ERROR(cudaMalloc((void**) &ftemp_dev, sizeof(double) * fcnt));
			HANDLE_ERROR(cudaMalloc((void**) &rtemp_dev, sizeof(double) * rcnt));

			HANDLE_ERROR(cudaMalloc((void**) &fenergy_dev, sizeof(double) * fcnt));
			HANDLE_ERROR(cudaMalloc((void**) &renergy_dev, sizeof(double) * rcnt));

			HANDLE_ERROR(cudaMalloc((void**) &Ooffset_dev,
					sizeof(unsigned int) * (my_threshold + 1)));
			HANDLE_ERROR(cudaMalloc((void**) &result_dev, sizeof(float) * total_result));

			HANDLE_ERROR(cudaMemcpy(fP_dev, &myArraysStep5[myRank]->FP[old_fcnt * maxLen],
					sizeof(char) * fcnt * maxLen, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(rP_dev, &myArraysStep5[myRank]->RP[old_rcnt * maxLen],
					sizeof(char) * rcnt * maxLen, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(fpos_dev, &myArraysStep5[myRank]->Fpos[old_fcnt],
					sizeof(unsigned int) * fcnt, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(rpos_dev, &myArraysStep5[myRank]->Rpos[old_rcnt],
					sizeof(unsigned int) * rcnt, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(ftemp_dev, &myArraysStep5[myRank]->Ftemp[old_fcnt],
					sizeof(double) * fcnt, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(rtemp_dev, &myArraysStep5[myRank]->Rtemp[old_rcnt],
					sizeof(double) * rcnt, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(fenergy_dev, &myArraysStep5[myRank]->Fenergy[old_fcnt],
					sizeof(double) * fcnt, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(renergy_dev, &myArraysStep5[myRank]->Renergy[old_rcnt],
					sizeof(double) * rcnt, cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(Ooffset_dev, &Ooffset[myRank][0],
					sizeof(unsigned int) * (my_threshold + 1),
					cudaMemcpyHostToDevice));
			HANDLE_ERROR(cudaMemcpy(result_dev, &result[myRank][0],
					sizeof(float) * total_result, cudaMemcpyHostToDevice));

#ifdef DEBUG
			used_memory += sizeof(unsigned char) * maxLen * fcnt;
			used_memory += sizeof(unsigned char) * maxLen * rcnt;
			used_memory += sizeof(unsigned int) * fcnt;
			used_memory += sizeof(unsigned int) * rcnt;
			used_memory += sizeof(double) * fcnt;
			used_memory += sizeof(double) * rcnt;
			used_memory += sizeof(double) * fcnt;
			used_memory += sizeof(double) * rcnt;
			used_memory += sizeof(unsigned int) * (my_threshold + 1);
			used_memory += sizeof(float) * total_result;

			if (myRank == 0)
			cout << "GPU memory: " << used_memory / (1024 * 1024) << " MB "
			<< "result: "
			<< (sizeof(float) * total_result) / (1024 * 1024)
			<< "MB" << endl;

#endif

			pair_filtering<<<block, thread>>>(i, my_threshold, fP_dev, rP_dev,
					fsid_dev, rsid_dev, fpos_dev, rpos_dev, ftemp_dev,
					rtemp_dev, fenergy_dev, renergy_dev, result_dev,
					Ooffset_dev, myInput_dev);

			i += working_thread;
			old_fcnt += fcnt;
			old_rcnt += rcnt;

			cudaFree(ftemp_dev);
			cudaFree(rtemp_dev);
			cudaFree(fenergy_dev);
			cudaFree(renergy_dev);
			cudaFree(fpos_dev);
			cudaFree(rpos_dev);
			cudaFree(fP_dev);
			cudaFree(rP_dev);

#ifdef DEBUG
			used_memory -= sizeof(unsigned char) * maxLen * fcnt;
			used_memory -= sizeof(unsigned char) * maxLen * rcnt;
			used_memory -= sizeof(unsigned int) * fcnt;
			used_memory -= sizeof(unsigned int) * rcnt;
			used_memory -= sizeof(double) * fcnt;
			used_memory -= sizeof(double) * rcnt;
			used_memory -= sizeof(double) * fcnt;
			used_memory -= sizeof(double) * rcnt;
#endif
		}

		cudaMemcpy(result[myRank], result_dev, sizeof(float) * total_result,
				cudaMemcpyDeviceToHost);

		cudaDeviceSynchronize();

	}

	pthread_barrier_wait(&barrThreads);

	et1 = getCurrentTime();

	if(myRank == 0)
		probingTime += (et1 - st1);

	st1 = getCurrentTime();

	passed = 0;
	passed = write1(mySmallThreshold[tmpIdx], result[tmpIdx], myRank, passed,
			Ooffset[tmpIdx], myArraysStep5, &sidsetH, sorted_fidx, sorted_ridx,
			myInput, finalFile[myRank]);

	pthread_barrier_wait(&barrThreads);

	et1 = getCurrentTime();
	et2 = getCurrentTime();

	if (myRank == 0)
		writeTime += (et1 - st1);

	pthread_mutex_lock(&valid_lock);
	total_passed += passed;
	pthread_mutex_unlock(&valid_lock);

	if (myRank < numOfGPUs) {
		cudaSetDevice(myRank);
		cudaFree(Ooffset_dev);
		cudaFree(result_dev);
		delete[] result[myRank];
		delete[] Ooffset[myRank];

#ifdef DEBUG
		used_memory -= sizeof(unsigned int) * (my_threshold + 1);
		used_memory -= sizeof(float) * total_result;
#endif
	}

	pthread_barrier_wait(&barrThreads);

	if (myRank == 0) {
		del(result);
		del(Ooffset);
		cout << "total_passed : " << total_passed << endl;
		cout << "each-thread mode : " << et2 - st2 << " sec"
				<< endl;
	}

	//-------------------------------CASE2------------------------------------

	passed = 0;
	int max_fcnt = 0;
	long used_memory2 = 0;
	numOfStreams = 2;
	unsigned int max_result = 0;

	long old_fcnt2[numOfStreams];
	long old_rcnt2[numOfStreams];
	long fcnt2[numOfStreams];
	long rcnt2[numOfStreams];

	unsigned char **fP_dev2, **rP_dev2;
	unsigned int **fpos_dev2, **rpos_dev2;
	double **ftemp_dev2, **rtemp_dev2;
	double **fenergy_dev2, **renergy_dev2;
	float **result_dev2;
	unsigned int **Ooffset_dev2;

	pthread_barrier_wait(&barrThreads);

#ifdef DEBUG
	if (myRank < numOfGPUs) {

		pthread_mutex_lock(&valid_lock);
		cout << "max : "
		<< myArraysStep5[myRank]->FPoffset[my_total_sid]
		- myArraysStep5[myRank]->FPoffset[my_total_sid - 1] << " "
		<< myArraysStep5[myRank]->FPoffset[my_total_sid + 1]
		- myArraysStep5[myRank]->FPoffset[my_total_sid] << endl;
		pthread_mutex_unlock(&valid_lock);
	}
#endif

	st1 = getCurrentTime();
	st2 = getCurrentTime();

	for (int k = 512; k >= 64;) {
		max_result = 0;
		working_thread = k;
		used_memory2 = 0;

		fcnt = myArraysStep5[tmpIdx]->FPoffset[myLargeThreshold[tmpIdx]]
				- myArraysStep5[tmpIdx]->FPoffset[myLargeThreshold[tmpIdx]
						- working_thread];
		rcnt = myArraysStep5[tmpIdx]->RPoffset[myLargeThreshold[tmpIdx]]
				- myArraysStep5[tmpIdx]->RPoffset[myLargeThreshold[tmpIdx]
						- working_thread];

		for (int i = (myLargeThreshold[tmpIdx] - working_thread);
				i < myLargeThreshold[tmpIdx]; i++) {
			int tmp_fcnt = myArraysStep5[tmpIdx]->FPoffset[i + 1]
					- myArraysStep5[tmpIdx]->FPoffset[i];
			int tmp_rcnt = myArraysStep5[tmpIdx]->RPoffset[i + 1]
					- myArraysStep5[tmpIdx]->RPoffset[i];
			max_result += tmp_fcnt * tmp_rcnt;
		}

		used_memory2 += memSID;
		used_memory2 += numOfStreams * sizeof(unsigned char) * maxLen * fcnt;
		used_memory2 += numOfStreams * sizeof(unsigned char) * maxLen * rcnt;
		used_memory2 += numOfStreams * sizeof(unsigned int) * fcnt;
		used_memory2 += numOfStreams * sizeof(unsigned int) * rcnt;
		used_memory2 += numOfStreams * sizeof(double) * fcnt;
		used_memory2 += numOfStreams * sizeof(double) * rcnt;
		used_memory2 += numOfStreams * sizeof(double) * fcnt;
		used_memory2 += numOfStreams * sizeof(double) * rcnt;
		used_memory2 += numOfStreams * sizeof(float) * max_result;
		used_memory2 += numOfStreams * sizeof(unsigned int) * working_thread;

		used_memory2 = used_memory2 / 1024;
		used_memory2 = used_memory2 / 1024;

		if (used_memory2 > memGPU)
			k = k / 2;
		else
			break;
	}

	pthread_barrier_wait(&barrThreads);

	if (myRank < numOfGPUs) {
		block = working_thread;

		cudaSetDevice(myRank);

		fP_dev2 = new unsigned char*[numOfStreams];
		rP_dev2 = new unsigned char*[numOfStreams];
		fpos_dev2 = new unsigned int*[numOfStreams];
		rpos_dev2 = new unsigned int*[numOfStreams];
		ftemp_dev2 = new double*[numOfStreams];
		rtemp_dev2 = new double*[numOfStreams];
		fenergy_dev2 = new double*[numOfStreams];
		renergy_dev2 = new double*[numOfStreams];
		result_dev2 = new float*[numOfStreams];
		Ooffset_dev2 = new unsigned int*[numOfStreams];

		for (int s = 0; s < numOfStreams; s++) {

			cudaSetDevice(myRank);

			HANDLE_ERROR(
					cudaMalloc((void** )&fP_dev2[s],
							sizeof(unsigned char) * maxLen * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&rP_dev2[s],
							sizeof(unsigned char) * maxLen * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&fpos_dev2[s],
							sizeof(unsigned int) * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&rpos_dev2[s],
							sizeof(unsigned int) * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&ftemp_dev2[s], sizeof(double) * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&rtemp_dev2[s], sizeof(double) * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&fenergy_dev2[s],
							sizeof(double) * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&renergy_dev2[s],
							sizeof(double) * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&Ooffset_dev2[s],
							sizeof(unsigned int) * (working_thread + 1)));
			HANDLE_ERROR(
					cudaMalloc((void** )&result_dev2[s],
							sizeof(float) * max_result));
		}
	}

	cudaDeviceSynchronize();
	pthread_barrier_wait(&barrThreads);

	if (myRank < numOfGPUs) {
#ifdef DEBUG
		used_memory += numOfStreams * sizeof(unsigned char) * maxLen * fcnt;
		used_memory += numOfStreams * sizeof(unsigned char) * maxLen * rcnt;
		used_memory += numOfStreams * sizeof(unsigned int) * fcnt;
		used_memory += numOfStreams * sizeof(unsigned int) * rcnt;
		used_memory += numOfStreams * sizeof(double) * fcnt;
		used_memory += numOfStreams * sizeof(double) * rcnt;
		used_memory += numOfStreams * sizeof(double) * fcnt;
		used_memory += numOfStreams * sizeof(double) * rcnt;
		used_memory += numOfStreams * sizeof(float) * max_result;
		used_memory += numOfStreams * sizeof(unsigned int)
		* (working_thread + 1);

		if (myRank == 0)
		cout << "GPU memory: " << used_memory / (1024 * 1024) << " MB "
		<< "result: "
		<< (sizeof(float) * max_result) / (1024 * 1024) << "MB"
		<< endl;
#endif

	}

	int tmpnumOfStreams = 1000;

	if (myRank == 0) {
		result = new float*[numOfGPUs * tmpnumOfStreams];
		Ooffset = new unsigned int*[numOfGPUs * tmpnumOfStreams];
		total_result3 = new unsigned int[numOfGPUs * tmpnumOfStreams];
		case3fcnt = new unsigned int[numOfGPUs * tmpnumOfStreams];
	}

	cudaStream_t stream[numOfStreams];

	pthread_barrier_wait(&barrThreads);

	int myStream = 0;
	int totalStream;

	if (myRank < numOfGPUs) {
		cudaSetDevice(myRank);
		for (int s = 0; s < numOfStreams; s++) {
			//cudaSetDevice(myRank);
			HANDLE_ERROR(cudaStreamCreate(&stream[s]));
		}

		for (int i = my_threshold; i < myLargeThreshold[tmpIdx];
				i += working_thread * numOfStreams) {

			for (int s = 0;
					(s < numOfStreams)
							&& ((i + s * working_thread)
									< myLargeThreshold[myRank]); s++) {

				int currentStream = myStream + s;
				int gIdx = tmpIdx * tmpnumOfStreams + currentStream;
				totalStream = currentStream;
				Ooffset[gIdx] = new unsigned int[working_thread + 1];
				Ooffset[gIdx][0] = 0;

				int tmpIdx = i + s * working_thread;
				if ((tmpIdx + working_thread) < myLargeThreshold[myRank]) {
					fcnt = myArraysStep5[myRank]->FPoffset[tmpIdx + working_thread]
							- myArraysStep5[myRank]->FPoffset[tmpIdx];
					rcnt = myArraysStep5[myRank]->RPoffset[tmpIdx + working_thread]
							- myArraysStep5[myRank]->RPoffset[tmpIdx];

					max_fcnt = myArraysStep5[myRank]->FPoffset[tmpIdx + working_thread
							+ 1]
							- myArraysStep5[myRank]->FPoffset[tmpIdx + working_thread];

					for (int j = 0; j < working_thread; j++) {
						tmp_fcnt = myArraysStep5[myRank]->FPoffset[tmpIdx + j + 1]
								- myArraysStep5[myRank]->FPoffset[tmpIdx + j];
						tmp_rcnt = myArraysStep5[myRank]->RPoffset[tmpIdx + j + 1]
								- myArraysStep5[myRank]->RPoffset[tmpIdx + j];
						Ooffset[gIdx][j + 1] = Ooffset[gIdx][j]
								+ tmp_fcnt * tmp_rcnt;
					}

					if (max_fcnt < 512)
						thread = max_fcnt;
					else
						thread = 512;

					total_result3[gIdx] = Ooffset[gIdx][working_thread];

				} else {
					fcnt = myArraysStep5[myRank]->FPoffset[myLargeThreshold[myRank]]
							- myArraysStep5[myRank]->FPoffset[tmpIdx];
					rcnt = myArraysStep5[myRank]->RPoffset[myLargeThreshold[myRank]]
							- myArraysStep5[myRank]->RPoffset[tmpIdx];

					for (int j = tmpIdx; j < myLargeThreshold[myRank]; j++) {
						tmp_fcnt = myArraysStep5[myRank]->FPoffset[j + 1]
								- myArraysStep5[myRank]->FPoffset[j];
						tmp_rcnt = myArraysStep5[myRank]->RPoffset[j + 1]
								- myArraysStep5[myRank]->RPoffset[j];
						Ooffset[gIdx][j - tmpIdx + 1] =
								Ooffset[gIdx][j - tmpIdx] + tmp_fcnt * tmp_rcnt;
					}
					thread = 512;
					total_result3[gIdx] = Ooffset[gIdx][myLargeThreshold[myRank]
							- tmpIdx];
				}

				cudaStreamSynchronize(stream[s]);

				result[gIdx] = new float[total_result3[gIdx]];
				memset(result[gIdx], -1.0, total_result3[gIdx] * sizeof(float));
				fcnt2[s] = fcnt;
				rcnt2[s] = rcnt;
				old_fcnt2[s] = old_fcnt;
				old_rcnt2[s] = old_rcnt;
				old_fcnt += fcnt2[s];
				old_rcnt += rcnt2[s];

				HANDLE_ERROR(
						cudaMemcpyAsync(fP_dev2[s],
								&myArraysStep5[myRank]->FP[old_fcnt2[s] * maxLen],
								sizeof(char) * fcnt2[s] * maxLen,
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(rP_dev2[s],
								&myArraysStep5[myRank]->RP[old_rcnt2[s] * maxLen],
								sizeof(char) * rcnt2[s] * maxLen,
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(fpos_dev2[s],
								&myArraysStep5[myRank]->Fpos[old_fcnt2[s]],
								sizeof(unsigned int) * fcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(rpos_dev2[s],
								&myArraysStep5[myRank]->Rpos[old_rcnt2[s]],
								sizeof(unsigned int) * rcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(ftemp_dev2[s],
								&myArraysStep5[myRank]->Ftemp[old_fcnt2[s]],
								sizeof(double) * fcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(rtemp_dev2[s],
								&myArraysStep5[myRank]->Rtemp[old_rcnt2[s]],
								sizeof(double) * rcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(fenergy_dev2[s],
								&myArraysStep5[myRank]->Fenergy[old_fcnt2[s]],
								sizeof(double) * fcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(renergy_dev2[s],
								&myArraysStep5[myRank]->Renergy[old_rcnt2[s]],
								sizeof(double) * rcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(Ooffset_dev2[s], &Ooffset[gIdx][0],
								sizeof(unsigned int) * (working_thread + 1),
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(result_dev2[s], &result[gIdx][0],
								sizeof(float) * total_result3[gIdx],
								cudaMemcpyHostToDevice, stream[s]));

			}

			for (int s = 0;
					(s < numOfStreams)
							&& ((i + s * working_thread)
									< myLargeThreshold[myRank]); s++) {

				//cudaSetDevice(myRank);
				pair_filtering2<<<block, thread, 0, stream[s]>>>(
						i + s * working_thread, myLargeThreshold[myRank],
						fP_dev2[s], rP_dev2[s], fsid_dev, rsid_dev,
						fpos_dev2[s], rpos_dev2[s], ftemp_dev2[s],
						rtemp_dev2[s], fenergy_dev2[s], renergy_dev2[s],
						result_dev2[s], Ooffset_dev2[s], myInput_dev);
			}

			for (int s = 0;
					(s < numOfStreams)
							&& ((i + s * working_thread)
									< myLargeThreshold[myRank]); s++) {

				int currentStream = myStream + s;
				int gIdx = tmpIdx * tmpnumOfStreams + currentStream;

				//cudaSetDevice(myRank);
				HANDLE_ERROR(
						cudaMemcpyAsync(result[gIdx], result_dev2[s],
								sizeof(float) * total_result3[gIdx],
								cudaMemcpyDeviceToHost, stream[s]));

			}

			myStream += numOfStreams;
		}

		cudaDeviceSynchronize();
		for (int s = 0; s < numOfStreams; s++) {
			//cudaSetDevice(myRank);
			HANDLE_ERROR(cudaStreamDestroy(stream[s]));
		}
	}

	pthread_barrier_wait(&barrThreads);

	et1 = getCurrentTime();

	if(myRank == 0)
		probingTime += (et1 - st1);

	st1 = getCurrentTime();

	int s = 0;
	for (int i = my_threshold; i < myLargeThreshold[tmpIdx]; i +=
			working_thread) {
		int gIdx = tmpIdx * tmpnumOfStreams + s;
		passed = write2(i, working_thread, result[gIdx], myRank, passed,
				Ooffset[gIdx], total_result3[gIdx], myArraysStep5, &sidsetH, sorted_fidx,
				sorted_ridx, myInput,
				finalFile[myRank], myLargeThreshold[tmpIdx]);
		s++;
	}

	pthread_barrier_wait(&barrThreads);

	et1 = getCurrentTime();
	if(myRank == 0)
		writeTime += (et1 - st1);

	if (myRank < numOfGPUs) {
		cudaSetDevice(myRank);
		for (int s = 0; s < numOfStreams; s++) {
			cudaFree(ftemp_dev2[s]);
			cudaFree(rtemp_dev2[s]);
			cudaFree(fenergy_dev2[s]);
			cudaFree(renergy_dev2[s]);
			cudaFree(fpos_dev2[s]);
			cudaFree(rpos_dev2[s]);
			cudaFree(fP_dev2[s]);
			cudaFree(rP_dev2[s]);
			cudaFree(result_dev2[s]);
			cudaFree(Ooffset_dev2[s]);
		}

		for(int s = 0; s < totalStream; s++){
			int gIdx = tmpIdx * tmpnumOfStreams + s;
			del(Ooffset[gIdx]);
			del(result[gIdx]);
		}

#ifdef DEBUG
		used_memory = memSID;
#endif

	}

	et2 = getCurrentTime();

	pthread_mutex_lock(&valid_lock);
	total_passed += passed;
	pthread_mutex_unlock(&valid_lock);

	pthread_barrier_wait(&barrThreads);

	if (myRank == 0) {
		cout << "passed : " << total_passed << endl;
		cout << "block-threads mode : " << et2 - st2 << " sec"
				<< endl;
	}

	//--------------------------CASE3----------------------------------------

	long my_maxrP = 0;
	block = 1024;
	thread = 512;
	working_thread = block;
	passed = 0;

	numOfStreams = 2;

	st1 = getCurrentTime();
	st2 = getCurrentTime();

	if (myRank < numOfGPUs) {

		cudaSetDevice(myRank);
		my_maxrP = 0;

		fcnt = myArraysStep5[myRank]->FPoffset[my_total_sid]
				- myArraysStep5[myRank]->FPoffset[my_total_sid - 1];
		rcnt = myArraysStep5[myRank]->RPoffset[my_total_sid]
				- myArraysStep5[myRank]->RPoffset[my_total_sid - 1];
		my_maxrP = rcnt;

		if (myRank == 0)
			cout << block << " " << thread << " " << fcnt
					<< " " << rcnt << endl;

		for (int s = 0; s < numOfStreams; s++) {


			//cudaSetDevice(myRank);
			HANDLE_ERROR(
					cudaMalloc((void** )&fP_dev2[s],
							sizeof(unsigned char) * maxLen * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&rP_dev2[s],
							sizeof(unsigned char) * maxLen * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&fpos_dev2[s],
							sizeof(unsigned int) * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&rpos_dev2[s],
							sizeof(unsigned int) * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&ftemp_dev2[s], sizeof(double) * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&rtemp_dev2[s], sizeof(double) * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&fenergy_dev2[s],
							sizeof(double) * fcnt));
			HANDLE_ERROR(
					cudaMalloc((void** )&renergy_dev2[s],
							sizeof(double) * rcnt));

			HANDLE_ERROR(
					cudaMalloc((void** )&result_dev2[s],
							sizeof(float) * my_maxrP * fcnt));
		}

#ifdef DEBUG
		used_memory += numOfStreams * sizeof(unsigned char) * maxLen * fcnt;
		used_memory += numOfStreams * sizeof(unsigned char) * maxLen * rcnt;
		used_memory += numOfStreams * sizeof(unsigned int) * fcnt;
		used_memory += numOfStreams * sizeof(unsigned int) * rcnt;
		used_memory += numOfStreams * sizeof(double) * fcnt;
		used_memory += numOfStreams * sizeof(double) * rcnt;
		used_memory += numOfStreams * sizeof(double) * fcnt;
		used_memory += numOfStreams * sizeof(double) * rcnt;
		used_memory += numOfStreams * sizeof(float) * my_maxrP * fcnt;


		if (myRank == 0)
		cout << "GPU memory: " << used_memory / (1024 * 1024) << " MB "
		<< "result: "
		<< (sizeof(float) * my_maxrP * fcnt) / (1024 * 1024)
		<< "MB" << endl;
#endif

		for (int s = 0; s < numOfStreams; s++) {
			//cudaSetDevice(myRank);
			HANDLE_ERROR(cudaStreamCreate(&stream[s]));
		}

		int myStream = 0;
		for (int i = myLargeThreshold[tmpIdx]; i < my_total_sid; i +=
				numOfStreams) {
			for (int s = 0; (s < numOfStreams) && ((i + s) < my_total_sid);
					s++) {

				int gIdx = tmpIdx * tmpnumOfStreams + myStream + s;
				fcnt = myArraysStep5[myRank]->FPoffset[(i + s) + 1]
						- myArraysStep5[myRank]->FPoffset[(i + s)];
				rcnt = myArraysStep5[myRank]->RPoffset[(i + s) + 1]
						- myArraysStep5[myRank]->RPoffset[(i + s)];

				case3fcnt[gIdx] = fcnt;
				total_result3[gIdx] = rcnt;
				result[gIdx] = new float[fcnt * total_result3[gIdx]];
				memset(result[gIdx], -1.0,
						sizeof(float) * fcnt * total_result3[gIdx]);


				//cudaSetDevice(myRank);
				cudaStreamSynchronize(stream[s]);

				fcnt2[s] = fcnt;
				rcnt2[s] = rcnt;
				old_fcnt2[s] = old_fcnt;
				old_rcnt2[s] = old_rcnt;
				old_fcnt += fcnt2[s];
				old_rcnt += rcnt2[s];

#ifdef DEBUG
				cout << myRank << "_largeThreshold : " << i + s << " " << block
				<< " " << thread << " " << fcnt2[s] << " " << rcnt2[s]
				<< " " << fcnt2[s] * total_result3[gIdx] << endl;
#endif

				//cudaSetDevice(myRank);
				HANDLE_ERROR(
						cudaMemcpyAsync(fP_dev2[s],
								&myArraysStep5[myRank]->FP[old_fcnt2[s] * maxLen],
								sizeof(char) * fcnt2[s] * maxLen,
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(rP_dev2[s],
								&myArraysStep5[myRank]->RP[old_rcnt2[s] * maxLen],
								sizeof(char) * rcnt2[s] * maxLen,
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(fpos_dev2[s],
								&myArraysStep5[myRank]->Fpos[old_fcnt2[s]],
								sizeof(unsigned int) * fcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(rpos_dev2[s],
								&myArraysStep5[myRank]->Rpos[old_rcnt2[s]],
								sizeof(unsigned int) * rcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(ftemp_dev2[s],
								&myArraysStep5[myRank]->Ftemp[old_fcnt2[s]],
								sizeof(double) * fcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(rtemp_dev2[s],
								&myArraysStep5[myRank]->Rtemp[old_rcnt2[s]],
								sizeof(double) * rcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(fenergy_dev2[s],
								&myArraysStep5[myRank]->Fenergy[old_fcnt2[s]],
								sizeof(double) * fcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(renergy_dev2[s],
								&myArraysStep5[myRank]->Renergy[old_rcnt2[s]],
								sizeof(double) * rcnt2[s],
								cudaMemcpyHostToDevice, stream[s]));
				HANDLE_ERROR(
						cudaMemcpyAsync(result_dev2[s], &result[gIdx][0],
								sizeof(float) * fcnt2[s] * total_result3[gIdx],
								cudaMemcpyHostToDevice, stream[s]));

			}
			for (int s = 0; (s < numOfStreams) && ((i + s) < my_total_sid);
					s++) {

				//cudaSetDevice(myRank);
				int gIdx = tmpIdx * tmpnumOfStreams + myStream + s;
				pair_filtering3<<<block, thread, 0, stream[s]>>>(i + s,
						fP_dev2[s], rP_dev2[s], fsid_dev, rsid_dev,
						fpos_dev2[s], rpos_dev2[s], ftemp_dev2[s],
						rtemp_dev2[s], fenergy_dev2[s], renergy_dev2[s],
						result_dev2[s], total_result3[gIdx], myInput_dev);
			}

			for (int s = 0; (s < numOfStreams) && ((i + s) < my_total_sid);
					s++) {

				//cudaSetDevice(myRank);
				int gIdx = tmpIdx * tmpnumOfStreams + myStream + s;
				cudaMemcpyAsync(result[gIdx], result_dev2[s],
						sizeof(float) * fcnt2[s] * total_result3[gIdx],
						cudaMemcpyDeviceToHost, stream[s]);
			}
			myStream += numOfStreams;
		}
		cudaDeviceSynchronize();
	}

	pthread_barrier_wait(&barrThreads);

	et1 = getCurrentTime();

	if(myRank == 0)
		probingTime += (et1 - st1);

	st1 = getCurrentTime();

	myStream = 0;
	for (int i = myLargeThreshold[tmpIdx]; i < my_total_sid; i++) {
		int gIdx = tmpIdx * tmpnumOfStreams + myStream;

		passed = write3(i, working_thread, result[gIdx], total_result3[gIdx],
				myRank, passed, case3fcnt[gIdx], total_result3[gIdx], myArraysStep5,
				&sidsetH, sorted_fidx,
				sorted_ridx, myInput,
				finalFile[myRank]);

		myStream++;
	}

	pthread_barrier_wait(&barrThreads);

	et1 = getCurrentTime();
	if(myRank == 0)
		writeTime += (et1 - st1);

	pthread_mutex_lock(&valid_lock);
	total_passed += passed;
	pthread_mutex_unlock(&valid_lock);

	if (myRank < numOfGPUs) {

		for (int s = 0; s < numOfStreams; s++) {
			cudaFree(ftemp_dev2[s]);
			cudaFree(rtemp_dev2[s]);
			cudaFree(fenergy_dev2[s]);
			cudaFree(renergy_dev2[s]);
			cudaFree(fpos_dev2[s]);
			cudaFree(rpos_dev2[s]);
			cudaFree(fP_dev2[s]);
			cudaFree(rP_dev2[s]);
			cudaFree(result_dev2[s]);
			cudaStreamDestroy(stream[s]);
		}

#ifdef DEBUG
		used_memory = memSID;
#endif


		pthread_barrier_wait(&barrGPUs);

		et2 = getCurrentTime();
		//#endif
		if (myRank == 0) {
			cout << "passed : " << total_passed << endl;
			cout << "all-threads mode : " << et2 - st2 << " sec"
					<< endl;
		}
	}
	fclose(finalFile[myRank]);
	return NULL;
}

void finalSort(){
	string command;
	ifstream fin;

	command = "sort -k1,1 -k3,3n -S " + to_string(myInput->bufferSize) + "% --parallel " + to_string(myInput->numOfThreads)
			+ " -T " + string(myInput->dirPath) + " " + 
			string(myInput->outputPath) + "_*" + " -o "
			+ string(myInput->outputPath);
	sysCall(command);
}

void primerHCheck(){
	long numOfTrues = 0;
	long numOfFalses = 0;

	for(auto it = primerH.begin(); it != primerH.end(); it++){
		if((it->second) == true)
			numOfTrues++;
		else
			numOfFalses++;
	}
	cout << "primerH true: " << numOfTrues << " false: " << numOfFalses << "\n\n";
}

void* writeOutput(void *param) {
	paramThread *writeInput = (paramThread *) param;
	long myRank = writeInput->rank;
	int stage = writeInput->myParam;

	char *primer, *val;
	char *ptr, *ptr2;
	long line = 0;

	char buf[MAX_BUF_SIZE];

#ifdef DEBUG
	long cnt = 0;
#endif

	ifstream fin;
	fileReadOpen(&fin, myInput->c2Path, myRank);

	FILE* fout;
	if (stage == 3)
		fileWriteOpen(&fout, myInput->c3Path, myRank);
	else if (stage == 4)
		fileWriteOpen(&fout, myInput->c4Path1, myRank);
	else
		fileWriteOpen(&fout, myInput->c4Path2, myRank);

	while (fin.good()) {
		fin.getline(buf, MAX_BUF_SIZE);
		if (fin.gcount() > 0) {
			ptr = buf;
			ptr2 = strchr(ptr, '\t');
			*ptr2 = '\0';
			primer = ptr;
			val = ++ptr2;

#ifdef DEBUG
			if ((line >= cnt * 1000000) && (myRank == 0)) {
				cout << primer << " " << val << endl;
				cnt++;
			}
#endif

			if (primerH.find(primer) != primerH.end()) {
				if ((primerH.find(primer)->second) == true)
					fprintf(fout, "%s\t%s\n", primer, val);
			}
		}
		line++;
	}

	fin.close();
	fclose(fout);
	return NULL;
}

void initialization(inputParameter* input) {

	cudaDeviceProp  prop;
	cudaGetDeviceProperties(&prop, 0);
    memGPU = prop.totalGlobalMem / (1024 * 1024);

	pthread_mutex_init(&primerH_lock, NULL);
	pthread_mutex_init(&sidsetH_lock, NULL);
	pthread_mutex_init(&seedHf_lock, NULL);
	pthread_mutex_init(&seedHr_lock, NULL);

	myInput = input;
	numOfThreads = input->numOfThreads;
	pthread_barrier_init(&barrThreads, 0, numOfThreads);
	numOfGPUs = input->numOfGPUs;
	pthread_barrier_init(&barrGPUs, 0, numOfGPUs);
	pthread_barrier_init(&barr2GPUs, 0, 2 * numOfGPUs);

	minLen = input->minLen;
	maxLen = input->maxLen;
	minGC = input->minGC;
	maxGC = input->maxGC;
	minTM = input->minTM;
	maxTM = input->maxTM;
	maxSC = input->maxSC;
	endMaxSC = input->endMaxSC;
	endMinDG = input->endMinDG;
	maxHP = input->maxHP;
	contiguous = input->contiguous;

	lenDiff = input->lenDiff;
	TMDiff = input->TMDiff;
	minPS = input->minPS;
	maxPS = input->maxPS;
	maxPC = input->maxPC;
	endMaxPC = input->endMaxPC;
	DGlen = 5;

	vecData = new vector<char*> [numOfThreads];
	myArraysC3 = new arraysForStep4* [numOfGPUs];
	myArraysC3_dev = new arraysForStep4* [numOfGPUs];
	for(int i = 0; i < numOfGPUs; i++){
		myArraysC3[i] = new arraysForStep4;
		myArraysC3_dev[i] = new arraysForStep4;
	}
	startPworkload = new unsigned int[numOfThreads];
	numOfStreams = 100;
}

void initialization2(long input_workload) {

	tmpArraysStep5 = new arraysForStep5;
	myArraysStep5 = new arraysForStep5 *[numOfGPUs];
	for(int i = 0; i < numOfGPUs; i++)
		myArraysStep5[i] = new arraysForStep5;

	sorted_fidx = new unsigned int[total_input_line];
	sorted_ridx = new unsigned int[total_input_line];

	fsid = new unsigned int[total_sid + 1];
	rsid = new unsigned int[total_sid + 1];
	sorted_rsid = new unsigned int[total_sid + 1];

	fP_idx = new int[numOfGPUs];
	rP_idx = new int[numOfGPUs];
	sid_idx = new int[numOfGPUs];
	mySmallThreshold = new int[numOfGPUs];
	myLargeThreshold = new int[numOfGPUs];

	result = new float*[numOfGPUs];
	Ooffset = new unsigned int*[numOfGPUs];

	tmpArraysStep5->FP = new unsigned char[maxLen * total_input_line];
	tmpArraysStep5->RP = new unsigned char[maxLen * total_input_line];
	tmpArraysStep5->Fpos = new unsigned int[total_input_line];
	tmpArraysStep5->Rpos = new unsigned int[total_input_line];
	tmpArraysStep5->Ftemp = new double[total_input_line];
	tmpArraysStep5->Fenergy = new double[total_input_line];
	tmpArraysStep5->Rtemp = new double[total_input_line];
	tmpArraysStep5->Renergy = new double[total_input_line];
	for (int i = 0; i < numOfGPUs; i++) {
		myArraysStep5[i]->FP = new unsigned char[input_workload * maxLen];
		myArraysStep5[i]->RP = new unsigned char[input_workload * maxLen];
		myArraysStep5[i]->Fpos = new unsigned int[input_workload];
		myArraysStep5[i]->Rpos = new unsigned int[input_workload];
		myArraysStep5[i]->FPoffset = new unsigned int[input_workload];
		myArraysStep5[i]->RPoffset = new unsigned int[input_workload];
		myArraysStep5[i]->Ftemp = new double[input_workload];
		myArraysStep5[i]->Rtemp = new double[input_workload];
		myArraysStep5[i]->Fenergy = new double[input_workload];
		myArraysStep5[i]->Renergy = new double[input_workload];
	}

	total_passed = 0;
}

#endif /* GPRIMER_HPP_ */
