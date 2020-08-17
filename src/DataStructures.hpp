/*
 * DataStructures.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: jmbae
 */

#ifndef DATASTRUCTURES_HPP_
#define DATASTRUCTURES_HPP_

struct arraysForStep4 {
	unsigned long *PFoffset;
	unsigned long *PRoffset;
	char *PF;
	char *PR;
	unsigned long *SIDFoffset;
	unsigned long *SIDRoffset;
	unsigned int *SIDF;
	unsigned int *SIDR;
	unsigned int *SEEDF;
	unsigned int *SEEDR;
	unsigned int *outputF;
	unsigned int *outputR;
};

struct arraysForStep5 {
	unsigned int *FPoffset;
	unsigned int *RPoffset;
	unsigned char *FP;
	unsigned char *RP;
	unsigned int *Fpos;
	unsigned int *Rpos;
	double *Ftemp;
	double *Rtemp;
	double *Fenergy;
	double *Renergy;
};


#endif /* DATASTRUCTURES_HPP_ */
