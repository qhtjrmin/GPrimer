/*
 * Step5Kernel.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: jmbae
 */

#ifndef STEP5KERNEL_HPP_
#define STEP5KERNEL_HPP_

template <typename T>
__device__ T absol_dev(T a){
	if(a>=0) return a;
	else return (-a);
}

__device__ unsigned char* reverse(unsigned char* P, int PLen){
	unsigned char* tmp = (unsigned char*) malloc(sizeof(unsigned char) * (PLen));
	memset(tmp, 0, PLen);
	int j = PLen - 1;
	#pragma unroll
	for(int i = 0; i < PLen; i++){
		tmp[i] = P[j];
		j--;
	}
	return tmp;
}

__device__ int pairCmp(unsigned char* fP, unsigned char* rP, int fLen, int rLen){
	int max = 0;
	unsigned char* rP2 = reverse(rP, rLen);
	int minimum = min(fLen, rLen);

	#pragma unroll
	for (int k = 1; k <= minimum; k++) {
		int tmpscore = 0;
		#pragma unroll
		for (int i = 0, j = rLen - k; i < k ; i++, j++) {
			if (fP[i] == 'A') {
				if (rP2[j] == 'T') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'T') {
				if (rP2[j] == 'A') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'C') {
				if (rP2[j] == 'G') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'G') {
				if (rP2[j] == 'C') tmpscore ++;
				else tmpscore --;
			}
			if (tmpscore < 0) tmpscore = 0;
			if (tmpscore > max) max = tmpscore;
		}
	}

	#pragma unroll
	for (int k = minimum - 1; k > 0; k--) {
		int tmpscore = 0;
		#pragma unroll
		for (int i = fLen - k, j = 0; j < k && i < k; i++, j++) {
			if (fP[i] == 'A') {
				if (rP2[j] == 'T') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'T') {
				if (rP2[j] == 'A') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'C') {
				if (rP2[j] == 'G') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'G') {
				if (rP2[j] == 'C') tmpscore ++;
				else tmpscore --;
			}
			if (tmpscore < 0) tmpscore = 0;
			if (tmpscore > max) max = tmpscore;
		}
	}
	free(rP2);
	rP2 = NULL;
	return max;
}

__device__ int endPairCmp(unsigned char* fP, unsigned char* rP, int fLen, int rLen){
	int max = 0;
	int minimum = min(fLen, rLen);
	int len = fLen/2;

	#pragma unroll
	for (int k = minimum - 1; k > len; k--) {
		int tmpscore = 0;
		#pragma unroll
		for (int i = fLen - k, j = 0; j < k && i < k; i++, j++) {
			if (fP[i] == 'A') {
				if (rP[j] == 'T') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'T') {
				if (rP[j] == 'A') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'C') {
				if (rP[j] == 'G') tmpscore ++;
				else tmpscore --;
			}
			else if (fP[i] == 'G') {
				if (rP[j] == 'C') tmpscore ++;
				else tmpscore --;
			}
			if (tmpscore < 0) tmpscore = 0;
			if (tmpscore > max) max = tmpscore;
		}
	}
	return max;
}

//single
__device__ double GCcontent(unsigned char* P, int Len){
	int CountOfGC = 0;
	for(int cnt = 0; cnt < Len; cnt++){
		if((P[cnt] == 'G') || (P[cnt] == 'C')){
			CountOfGC++;
		}
	}
	double val = ((double)CountOfGC / (double) Len) * 100.0;
	return val;
}

__device__ int selfCmp(unsigned char* P, int P_len){

	int max = 0;

	unsigned char* rvsP = reverse(P, P_len);

	for (int k = 1; k <= P_len; k++) {
        int tmpscore = 0;
        for (int i = 0, j = P_len - k; i < k; i++, j++) {
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
	free(rvsP);
	rvsP = NULL;
	return max;
}

__device__ int endSelfCmp(unsigned char* P, int P_len){
	unsigned char* selfRvs = reverse(P, P_len);
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
	free(selfRvs);
	selfRvs = NULL;
	return max;
}

__device__ int hairpin(unsigned char* P, int P_len){

	int max = 0;
	unsigned char* pre_s = (unsigned char*) malloc(sizeof(unsigned char) * (P_len / 2 + 1));
	for(int i = 0; i < P_len/2; i++){
		pre_s[i] = P[i];
	}
	unsigned char* temp_s = (unsigned char*) malloc(sizeof(unsigned char) * (P_len - (P_len/2) + 1));
	for(int i = P_len/2; i < P_len; i++){
		temp_s[max] = P[i];
		max++;
	}
	max = 0;
	unsigned char* post_s = reverse(temp_s, P_len - (P_len/2));

	for (int k = 1; k <= P_len/2; k++) {
         int tmpscore = 0;
        for (int i = 0, j = P_len - (P_len/2) - k; i < k; i++, j++) {

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

	for (int k = (P_len / 2) - 1; k > 0; k--) {
        int tmpscore = 0;
        for (int i = P_len / 2 - k, j = 0; j < k; i++, j++) {

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
	free(post_s);
	free(temp_s);
	free(pre_s);
	post_s = NULL; temp_s = NULL; pre_s = NULL;
	return max;
}

__device__ int GCclamp(unsigned char* P, int P_len){
	int gcclamp = 0;
	for(int i = 0; i < 5; i++)
		if(P[P_len - i - 1] == 'G' || P[P_len - i - 1] == 'C')
			gcclamp++;
	return gcclamp;
}

#endif /* STEP5KERNEL_HPP_ */
