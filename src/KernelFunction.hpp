/*
 * KernelFunction.h
 *
 *  Created on: Feb 11, 2020
 *      Author: jmbae
 */

#ifndef KERNELFUNCTION_H_
#define KERNELFUNCTION_H_


#include "Step5Kernel.hpp"
#include <thrust/sort.h>
#include <thrust/set_operations.h>

__global__ void stage4_probing(unsigned int mismatch, unsigned int threshold, unsigned int seedlen,
		unsigned int* seed_idx_dev, unsigned long* c1_off, char* c1_P_dev,
		unsigned int* c1_sid, unsigned long* c1_sid_cnt, unsigned long* c3_off,
		char* c3_P_dev, unsigned int* c3_sid, unsigned long* c3_sid_cnt,
		unsigned int* result_dev, int maxLen) {

	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int diff;
	int seed_idx;
	unsigned long c1_sp, c1_ep, c3_sp, c3_ep;
	unsigned long c1_sid_sp, c1_sid_ep, c3_sid_sp, c3_sid_ep;
	int Plen, count;
	int seedLen = seedlen;
	int mismatch_threshold = mismatch;
	char* c3_P, *c1_P;

	for (int idx = tid; idx < threshold; idx += gridDim.x * blockDim.x) {

		seed_idx = seed_idx_dev[idx];
		c1_sp = c1_off[idx];
		c1_ep = c1_off[idx + 1];
		c3_sp = c3_off[idx];
		c3_ep = c3_off[idx + 1];

		for (unsigned long i = c3_sp; i < c3_ep; i++) {
			if (result_dev[i] == 0) {
				c3_P = &c3_P_dev[i * maxLen];

#pragma unroll
				for (int k = maxLen - 1; k > 0; k--)
					if (c3_P[k] != '\0') {
						Plen = k + 1;
						break;
					}

//			if(Plen != maxLen && idx < 10)
//				printf("%d : %s\n", idx, c3_P);

				for (unsigned long j = c1_sp; j < c1_ep; j++) {
					c1_P = &c1_P_dev[j * maxLen];
					diff = -1;
					c3_sid_sp = c3_sid_cnt[i];
					c3_sid_ep = c3_sid_cnt[i + 1];
					c1_sid_sp = c1_sid_cnt[j];
					c1_sid_ep = c1_sid_cnt[j + 1];
#pragma unroll
					for (int k = c1_sid_sp; k < c1_sid_ep; k++) {
						if (diff == -1)
							thrust::set_difference(thrust::seq, c1_sid + k,
									c1_sid + (k + 1), c3_sid + c3_sid_sp,
									c3_sid + c3_sid_ep, &diff);
						else
							break;
					}

					if (diff != -1) {
						count = 0;
#pragma unroll
						for (int k = 0; k < seed_idx; k++)
							if (c3_P[k] != c1_P[k])
								count++;
						if (count <= mismatch_threshold) {
#pragma unroll
							for (int k = seed_idx + seedLen; k < Plen; k++)
								if (c3_P[k] != c1_P[k])
									count++;
						}
						if (count <= mismatch_threshold) {
//						if(Plen != maxLen){
//							printf("%s %d\n%s\n", c1_P, diff, c3_P);
//						}
							atomicAdd(&result_dev[i], 1);
							break;
						}
					}
				}
			}
		}
	}
}

__global__ void stage4_probing2(unsigned int mismatch, unsigned int threshold, unsigned int total_key,
		unsigned int seedlen, unsigned int* seed_idx_dev, unsigned long* c1_off,
		char* c1_P_dev, unsigned int* c1_sid, unsigned long* c1_sid_cnt,
		unsigned long* c3_off, char* c3_P_dev, unsigned int* c3_sid,
		unsigned long* c3_sid_cnt, unsigned int* result_dev, int maxLen) {

	int tid = threshold + blockIdx.x;
	int diff;
	int seed_idx;
	unsigned long c1_sp, c1_ep, c3_sp, c3_ep;
	unsigned long c1_sid_sp, c1_sid_ep, c3_sid_sp, c3_sid_ep;
	int Plen, count;
	int seedLen = seedlen;
	int mismatch_threshold = mismatch;
	char* c3_P, *c1_P;

	for (int idx = tid; idx < total_key; idx += gridDim.x) {

		seed_idx = seed_idx_dev[idx];
		c1_sp = c1_off[idx];
		c1_ep = c1_off[idx + 1];
		c3_sp = c3_off[idx] + threadIdx.x;
		c3_ep = c3_off[idx + 1];

		if (threadIdx.x < c3_off[idx + 1] - c3_off[idx])
			for (unsigned long i = c3_sp; i < c3_ep; i += blockDim.x) {
				if (result_dev[i] == 0) {
					c3_P = &c3_P_dev[i * maxLen];

#pragma unroll
					for (int k = maxLen - 1; k > 0; k--)
						if (c3_P[k] != '\0') {
							Plen = k + 1;
							break;
						}

					for (unsigned long j = c1_sp; j < c1_ep; j++) {
						c1_P = &c1_P_dev[j * maxLen];
						diff = -1;
						c3_sid_sp = c3_sid_cnt[i];
						c3_sid_ep = c3_sid_cnt[i + 1];
						c1_sid_sp = c1_sid_cnt[j];
						c1_sid_ep = c1_sid_cnt[j + 1];
#pragma unroll
						for (int k = c1_sid_sp; k < c1_sid_ep; k++) {
							if (diff == -1)
								thrust::set_difference(thrust::seq, c1_sid + k,
										c1_sid + (k + 1), c3_sid + c3_sid_sp,
										c3_sid + c3_sid_ep, &diff);
							else
								break;
						}

						if (diff != -1) {
							count = 0;
#pragma unroll
							for (int k = 0; k < seed_idx; k++)
								if (c3_P[k] != c1_P[k])
									count++;
							if (count <= mismatch_threshold) {
#pragma unroll
								for (int k = seed_idx + seedLen; k < Plen; k++)
									if (c3_P[k] != c1_P[k])
										count++;
							}
							if (count <= mismatch_threshold) {
								atomicAdd(&result_dev[i], 1);
								break;
							}
						}
					}
				}
			}
	}
}

__device__ double singlePenalty (unsigned char* P, int Len, double temp, double energy, inputParameter* input) {

	int minLen = input->minLen; int maxLen = input->maxLen;
	double minGC = input->minGC; double maxGC = input->maxGC;
	double minTM = input->minTM; double maxTM = input->maxTM;
	int maxSC = input->maxSC;
	int endMaxSC = input->endMaxSC;
	int endMinDG = input->endMinDG;
	int maxHP = input->maxHP;

	//default optimal values
	double penalty = 0.0;
	int optLen = (minLen + maxLen) / 2;
	double optTm = 60.0;
	double optGC = 50.0;

//	Constraints Pconstraints = Constraints (primer);

	int lenDiff = 0, len = 0;
	if (optLen - minLen > maxLen - optLen) {
		lenDiff = optLen - minLen;
	} else if (optLen - minLen < maxLen - optLen) {
		lenDiff = maxLen - optLen;
	}
	else
		lenDiff = maxLen - optLen;
	if (Len != optLen)
	{
		len = absol_dev(Len - optLen);
		penalty =  (double)((double)len / (double)lenDiff);
	}

	double TmDiff = 0.0;
	if (optTm - minTM > maxTM - optTm)
		TmDiff = optTm - minTM;
	else if (optTm - minTM < maxTM - optTm)
		TmDiff = maxTM - optTm;
	else
		TmDiff = maxTM - optTm;
	double tmd = absol_dev(temp - optTm);
	if (temp != optTm)
		penalty += tmd / TmDiff;

	double gcDiff = 0.0;
	if (optGC - minGC > maxGC - optGC)
		gcDiff = optGC - minGC;
	else if (optGC - minGC < maxGC - optGC)
		gcDiff = maxGC - optGC;
	else
		gcDiff = maxGC - optGC;
	double gcContent = absol_dev(GCcontent(P, Len) - optGC) ;
	if (gcContent != optGC)
		penalty += gcContent / gcDiff;

	int scDiff = maxSC - 0;
	int selfCmp_value = selfCmp(P, Len);
	double scNorm = 0.0;
	scNorm = ((double)selfCmp_value /(double) scDiff);
	penalty += scNorm;

	int escDiff = endMaxSC-0;
	int endSelfCmp_value = endSelfCmp(P, Len);
	double escNorm = 0.0;
	escNorm = ((double)endSelfCmp_value / (double)escDiff);
	penalty += escNorm;

	int hairpinDiff = maxHP - 0;
	int hairpin_value = hairpin(P, Len);
	double hpNorm = 0.0;
	hpNorm = ((double)hairpin_value / (double) hairpinDiff);
	penalty += hpNorm;


	int GCclampDiff = 5 - 2; // need to adjust
	int gcclamp = GCclamp(P, Len);
	double GCclampNorm = 0.0;
	if (gcclamp > 2)
		GCclampNorm = ((double) (gcclamp - 2) / (double) GCclampDiff);
	penalty += GCclampNorm;

	double stabilityDiff = endMinDG;
	double stability = absol_dev(energy);
	double stabilityNorm = absol_dev(stability / stabilityDiff);
	penalty += stabilityNorm;

	return penalty;
}

__device__ double pairPenalty(double fmScore, double rmScore,
		unsigned char* fPrimer, unsigned char* rPrimer,
		int fLen, int rLen, int fpos, int rpos, double ftemp, double rtemp,
		inputParameter* input){

	int lenDiff = input->lenDiff;
	int TMDiff = input->TMDiff;
	int minPS = input->minPS;
	int maxPS = input->maxPS;
	int maxPC = input->maxPC;
	int endMaxPC = input->endMaxPC;

	double penalty = 0.0;
	double optProductSize = (maxPS + minPS) / 2.0;

	penalty = fmScore + rmScore;

	double PSlenDiff=0, len=0;
	if (optProductSize - minPS > maxPS - optProductSize) {
		PSlenDiff = optProductSize - minPS;
	} else if (optProductSize - minPS < maxPS - optProductSize) {
		PSlenDiff = maxPS - optProductSize;
	}
	else
		PSlenDiff = maxPS - optProductSize;
	int productSize = rpos - fpos + 1;
	if (productSize != optProductSize)
	{
		len = absol_dev(productSize - optProductSize);
		penalty +=  (double)len / (double)PSlenDiff;
	}

	int frlenDiff = lenDiff;
	int frDiff = absol_dev(fLen - rLen);
	double lenDiffNorm = ((double) frDiff / (double) frlenDiff);
	penalty += lenDiffNorm;

	double tmdiffDiff = TMDiff;
	double tmdiff = absol_dev(ftemp - rtemp);
	double tmDiffNorm = (tmdiff / tmdiffDiff);
	penalty += tmDiffNorm;

	int pcDiff = maxPC;
	int pairCmp_value = pairCmp(fPrimer, rPrimer, fLen, rLen);
	double pcDiffNorm = ((double)pairCmp_value / (double)pcDiff);
	penalty += pcDiffNorm;

	int epcDiff = endMaxPC;
	int endPairCmp_value = endPairCmp(fPrimer, rPrimer, fLen, rLen);
	double epcDiffNorm = ((double)endPairCmp_value / (double)epcDiff);
	penalty += epcDiffNorm;

	return penalty;
}

__global__ void pair_filtering(unsigned int sid, unsigned int threshold, unsigned char* fP_dev, unsigned char* rP_dev, unsigned int* fsid_dev, unsigned int* rsid_dev,
	unsigned int* fpos_dev, unsigned int* rpos_dev, double* ftemp_dev, double* rtemp_dev, double* fenergy_dev, double* renergy_dev, float* result_dev, unsigned int* Ooffset_dev
	, inputParameter* input){

	int tid = sid + threadIdx.x + blockDim.x * blockIdx.x;

	if(tid < threshold){
	//	result_dev[tid] = 0;
	int fst = fsid_dev[tid] - fsid_dev[sid];
	int fcnt = fsid_dev[tid + 1] - fsid_dev[tid];
	int rst = rsid_dev[tid] - rsid_dev[sid];
	int rcnt = rsid_dev[tid + 1] - rsid_dev[tid];
	double fTemp, rTemp, fEnergy, rEnergy;
	double fScore, rScore, pairScore;
	int fPos, rPos; bool check;
	unsigned char* fP, *rP;
	int fLen, rLen;

	int lenDiff = input->lenDiff;
	int TMDiff = input->TMDiff;
	int minPS = input->minPS;
	int maxPS = input->maxPS;
	int maxPC = input->maxPC;
	int endMaxPC = input->endMaxPC;
	int maxLen = input->maxLen;

	#pragma unroll
	for(int i = fst; i < (fst + fcnt); i++){
		fPos = fpos_dev[i];
		fP = &fP_dev[i * maxLen];
		#pragma unroll
		for(int k = maxLen - 1; k > 0; k--) if(fP[k] != '\0'){fLen = k + 1; break;}
		#pragma unroll
		for(int j = rst; j < (rst + rcnt); j++){
			check = true;
			rPos = rpos_dev[j];
			rP = &rP_dev[j * maxLen];
			#pragma unroll
			for(int k = maxLen - 1; k > 0; k--) if(rP[k] != '\0'){rLen = k + 1; break;}
			if(fPos >= rPos) check = false;
			if(check) if(absol_dev(fLen - rLen) > lenDiff) check = false;
			if(check) if(absol_dev(ftemp_dev[i] - rtemp_dev[j]) > TMDiff) check = false;
			if(check) if(((rPos - fPos + 1) < minPS) || ((rPos - fPos + 1) > maxPS)) check = false;
			if(check) if((pairCmp(fP, rP, fLen, rLen)) >= maxPC) check = false;
			if(check) if((endPairCmp(fP, rP, fLen, rLen)) >= endMaxPC) check = false;
			if(check){
				fTemp = ftemp_dev[i]; rTemp = rtemp_dev[j];
				fEnergy = fenergy_dev[i]; rEnergy = renergy_dev[j];
				fScore =  singlePenalty(fP, fLen, fTemp, fEnergy, input);
				rScore = singlePenalty(rP, rLen, rTemp, rEnergy, input);
				pairScore = pairPenalty(fScore, rScore, fP, rP, fLen, rLen, fPos, rPos, fTemp, rTemp, input);
				result_dev[Ooffset_dev[tid] + (i - fst) * rcnt + (j - rst)] = pairScore;
				}
		}
	}
	}
}

__global__ void pair_filtering2(unsigned int sid, unsigned my_total_sid, unsigned char* fP_dev, unsigned char* rP_dev, unsigned int* fsid_dev, unsigned int* rsid_dev,
	unsigned int* fpos_dev, unsigned int* rpos_dev, double* ftemp_dev, double* rtemp_dev, double* fenergy_dev, double* renergy_dev, float* result_dev, unsigned int* Ooffset_dev
	, inputParameter* input){

	int tid = sid + blockIdx.x;

	int fst = fsid_dev[tid] - fsid_dev[sid];
	int fcnt = fsid_dev[tid + 1] - fsid_dev[tid];
	int rst = rsid_dev[tid] - rsid_dev[sid];
	int rcnt = rsid_dev[tid + 1] - rsid_dev[tid];
	int fLen, rLen;
	double fTemp, fEnergy, rTemp, rEnergy;
	double fScore, rScore, pairScore;
	int fPos, rPos;
	bool check;
	unsigned char* fP, *rP;

	int lenDiff = input->lenDiff;
	int TMDiff = input->TMDiff;
	int minPS = input->minPS;
	int maxPS = input->maxPS;
	int maxPC = input->maxPC;
	int endMaxPC = input->endMaxPC;
	int maxLen = input->maxLen;

	if((tid < sid + gridDim.x) && (tid < my_total_sid) && (tid >= sid)){
	#pragma unroll
	for(int i = fst + threadIdx.x; i < (fst + fcnt); i += blockDim.x){
		fPos = fpos_dev[i];
		fP = &fP_dev[i * maxLen];
		#pragma unroll
		for(int k = maxLen - 1; k > 0; k--) if(fP[k] != '\0'){fLen = k + 1; break;}
		#pragma unroll
		for(int j = rst; j < (rst + rcnt); j++){
			check = true;
			rPos = rpos_dev[j];
			rP = &rP_dev[j * maxLen];
			#pragma unroll
			for(int k = maxLen - 1; k > 0; k--) if(rP[k] != '\0'){rLen = k + 1; break;}
			if(fPos >= rPos) check = false;
			if(check) if(absol_dev(fLen - rLen) > lenDiff) check = false;
			if(check) if(absol_dev(ftemp_dev[i] - rtemp_dev[j]) > TMDiff) check = false;
			if(check) if(((rPos - fPos + 1) < minPS) || ((rPos - fPos + 1) > maxPS)) check = false;
			if(check) if((pairCmp(fP, rP, fLen, rLen)) >= maxPC) check = false;
			if(check) if((endPairCmp(fP, rP, fLen, rLen)) >= endMaxPC) check = false;
			if(check){
				fTemp = ftemp_dev[i]; rTemp = rtemp_dev[j];
                fEnergy = fenergy_dev[i]; rEnergy = renergy_dev[j];
                fScore =  singlePenalty(fP, fLen, fTemp, fEnergy, input);
                rScore = singlePenalty(rP, rLen, rTemp, rEnergy, input);
                pairScore = pairPenalty(fScore, rScore, fP, rP, fLen, rLen, fPos, rPos, fTemp, rTemp, input);
                result_dev[Ooffset_dev[tid - sid] + (i - fst) * rcnt + (j - rst)] = pairScore;
			}
		}
	}
	}
}

__global__ void pair_filtering3(unsigned int sid, unsigned char* fP_dev, unsigned char* rP_dev, unsigned int* fsid_dev, unsigned int* rsid_dev,
    unsigned int* fpos_dev, unsigned int* rpos_dev, double* ftemp_dev, double* rtemp_dev, double* fenergy_dev, double* renergy_dev,
    float* result_dev, unsigned long maxrP,inputParameter* input){

    int tid = sid;

    int fcnt = fsid_dev[tid + 1] - fsid_dev[tid];
    int rcnt = rsid_dev[tid + 1] - rsid_dev[tid];

    int fLen, rLen;
	double fTemp, rTemp, fEnergy, rEnergy;
	double fScore, rScore, pairScore;
	int fPos, rPos;
	bool check;
	unsigned char* fP, *rP;

	int lenDiff = input->lenDiff;
	int TMDiff = input->TMDiff;
	int minPS = input->minPS;
	int maxPS = input->maxPS;
	int maxPC = input->maxPC;
	int endMaxPC = input->endMaxPC;
	int maxLen = input->maxLen;

    #pragma unroll
    for(int i = blockIdx.x; i < fcnt; i += gridDim.x){
        fPos = fpos_dev[i];
        fP = &fP_dev[i * maxLen];
        #pragma unroll
        for(int k = maxLen - 1; k > 0; k--) if(fP[k] != '\0'){fLen = k + 1; break;}
        #pragma unroll
        for(int j = threadIdx.x; j < rcnt; j += blockDim.x){
            check = true;
			rPos = rpos_dev[j];
            rP = &rP_dev[j * maxLen];
            #pragma unroll
            for(int k = maxLen - 1; k > 0; k--) if(rP[k] != '\0'){rLen = k + 1; break;}
            if(fPos >= rPos) check = false;
            if(check) if(absol_dev(fLen - rLen) > lenDiff) check = false;
            if(check) if(absol_dev(ftemp_dev[i] - rtemp_dev[j]) > TMDiff) check = false;
            if(check) if(((rPos - fPos + 1) < minPS) || ((rPos - fPos + 1) > maxPS)) check = false;
            if(check) if((pairCmp(fP, rP, fLen, rLen)) >= maxPC) check = false;
            if(check) if((endPairCmp(fP, rP, fLen, rLen)) >= endMaxPC) check = false;
			if(check){//if(tid == 7846 && fcnt == 9446) printf("%d,%d ", i, j);
				fTemp = ftemp_dev[i]; rTemp = rtemp_dev[j];
                fEnergy = fenergy_dev[i]; rEnergy = renergy_dev[j];
                fScore =  singlePenalty(fP, fLen, fTemp, fEnergy, input);
                rScore = singlePenalty(rP, rLen, rTemp, rEnergy, input);
                pairScore = pairPenalty(fScore, rScore, fP, rP, fLen, rLen, fPos, rPos, fTemp, rTemp, input);
				result_dev[i * maxrP + j] = pairScore;
				if(j  > maxrP) printf("ERR: %d %d %ld\n", i, j, maxrP);}
        }
    }
}

#endif /* KERNELFUNCTION_H_ */
