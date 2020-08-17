/*
 * Input.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: jmbae
 */

#ifndef INPUT_HPP_
#define INPUT_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <unistd.h>

struct inputParameter{
	inputParameter(){
		dirPath = new char[300];
		c1Path = new char[300];
		c1SidsetPath = new char[300];
		c2Path = new char[300];
		c3Path = new char[300];
		c4Path1 = new char[300];
		c4Path2 = new char[300];
		numOfGPUs = 1;
		numOfThreads = 20;
		isWriteOutput = 0;

		minLen = 19;
		maxLen = 23;
		minGC = 40;
		maxGC = 60;
		minTM = 58.0;
		maxTM = 62.0;
		maxSC = 5;
		endMaxSC = 4;
		endMinDG = -9;
		maxHP = 3;
		contiguous = 6;

		bufferSize = 30;
		lenDiff = 5;
		TMDiff=3;
		minPS=100;
		maxPS=250;
		maxPC=5;
		endMaxPC=4;
	}

	~inputParameter(){
		delete[] dirPath;
		delete[] c1SidsetPath;
		delete[] c1Path;
		delete[] c2Path;
		delete[] c3Path;
		delete[] c4Path1;
		delete[] c4Path2;
	}

	void printInputParameter(){
		cout << "<Files>" << endl;
		cout << "input: " << inputPath << ", C1: " << c1Path <<
				", C1': " << c1SidsetPath << ", C2: " << c2Path << ", C3: " << c3Path
				<< ", C4(k=1): " << c4Path1 << ", C4(k=2): " << c4Path2
				<< ", output " << outputPath << "\n\n";
		cout << "<Settings>" << endl;
		cout << "numOfGPUs: " << numOfGPUs << ", numOfThreads: " << numOfThreads << "\n\n";
		cout << "<Parameters for single filtering>" << endl;
		cout << "minLen: " << minLen << ", maxLen: " << maxLen
				<< ", minGC: " << minGC << ", maxGC: " << maxGC
				<< ", minTM: " << minTM << ", maxTM: " << maxTM
				<< ", maxSC: " << maxSC << ", endMaxSC: " << endMaxSC
				<< ", endMinDG: " << endMinDG
				<< ", maxHP: " << maxHP << ", contiguous: " << contiguous << "\n\n";
		cout << "<Parameters for pair filtering>" << endl;
		cout << "lenDiff: " << lenDiff << ", TMDiff: " << TMDiff
				<< ", minPS: " << minPS << ", maxPS: " << maxPS
				<< ", maxPC: " << maxPC << ", endMaxPC: " << endMaxPC << "\n\n";
	}

	char* dirPath;
	char* inputPath;
	char* c1Path;
	char* c1SidsetPath;
	char* c2Path;
	char* c3Path;
	char* c4Path1;
	char* c4Path2;
	char* outputPath;
	int numOfGPUs;
	int numOfThreads;
	int isWriteOutput;

	int minLen;
	int maxLen;
	float minGC;
	float maxGC;
	float minTM;
	float maxTM;
	int maxSC;
	int endMaxSC;
	int endMinDG;
	int maxHP;
	int contiguous;

	int bufferSize;
	int lenDiff = 5;
	int TMDiff=3;
	int minPS=100;
	int maxPS=250;
	int maxPC=5;
	int endMaxPC=4;
};

bool readInputParameter(int argc, char*argv[], inputParameter* input) {
	if (argc == 1)
	{
 		cout << "-i <input sequence DB (necessary)>" <<endl;
		cout << "-o <final output path (necessary)>" << endl;
		cout << "-d <working and storage directory (necessary)" << endl;
		cout << "-t <number of CPU threads (in default 20)>" << endl;
		cout << "-g <number of GPUs (in default 1)>" << endl;
		cout << "-w <write intermediate output (0: no, 1: yes, in default 0)>"<<endl;
		cout << "-b <the percentage(%) of utilized buffer size for sort (in default 30)>" << endl;
		cout << "-p1 <change the parameter about single filtering (0: no, 1: yes, in default 0)>" << endl;
		cout << "-p2 <change the parameter about pair filtering (0: no, 1: yes, in default 0)>" << endl;
		return false;
	}
	int argnr = 0;
	bool inputCheck = false;
	bool dirCheck = false;
	bool outputCheck = false;
	while (++argnr < argc)
	{
		if (!strcmp(argv[argnr], "-i")) {
			input->inputPath = argv[++argnr];
			inputCheck = true;
		}
		else if(!strcmp(argv[argnr], "-d")){
			string dir = argv[++argnr];
			string tmpFname;

			strcpy(input->dirPath, dir.c_str());

			tmpFname = dir + "/C1.txt";
			strcpy(input->c1Path, tmpFname.c_str());

			tmpFname = dir + "/C1'.txt";
			strcpy(input->c1SidsetPath, tmpFname.c_str());

			tmpFname = dir + "/C2.txt";
			strcpy(input->c2Path, tmpFname.c_str());

			tmpFname = dir + "/C3.txt";
			strcpy(input->c3Path, tmpFname.c_str());

			tmpFname = dir + "/C4_1.txt";
			strcpy(input->c4Path1, tmpFname.c_str());

			tmpFname = dir + "/C4_2.txt";
			strcpy(input->c4Path2, tmpFname.c_str());

			dirCheck = true;
		}
		else if(!strcmp(argv[argnr], "-o")){
			input->outputPath = argv[++argnr];
			outputCheck = true;
		}
		else if (!strcmp(argv[argnr], "-t")){
			input->numOfThreads = stoi(argv[++argnr]);
		}
		else if (!strcmp(argv[argnr], "-g")){
			input->numOfGPUs = stoi(argv[++argnr]);
		}
		else if (!strcmp(argv[argnr], "-w")){
			input->isWriteOutput = stoi(argv[++argnr]);
		}
		else if (!strcmp(argv[argnr], "-b")){
			input->bufferSize = stoi(argv[++argnr]);
		}
		else if(!strcmp(argv[argnr], "-p1")){
			int tmpInt; float tmpFloat;
			if(stoi(argv[++argnr]) == 1){

				cout << "## Parameters for single filtering" << endl;
				cout << "The minimum length of primer (in default 19): ";
				cin >> tmpInt;
				input->minLen = tmpInt;

				cout << "The maximum length of primer (in default 23): ";
				cin >> tmpInt;
				input->maxLen = tmpInt;

				cout << "The minimum GC ratios (in default 40): ";
				cin >> tmpInt;
				input->minGC = tmpInt;

				cout << "The maximum GC ratios (in default 60): ";
				cin >> tmpInt;
				input->maxGC = tmpInt;

				cout << "The minimum primer melting temperatures (in default 58.0): ";
				cin >> tmpFloat;
				input->minTM = tmpFloat;

				cout << "The maximum primer melting temperatures (in default 62.0): ";
				cin >> tmpFloat;
				input->maxTM = tmpFloat;

				cout << "The maximum self-complementarity (in default 5): ";
				cin >> tmpInt;
				input->maxSC = tmpInt;

				cout << "The maximum 3' end self-complementarity (in default 4): ";
				cin >> tmpInt;
				input->endMaxSC = tmpInt;

				cout << "The maximum contiguous residues (in default 6): ";
				cin >> tmpInt;
				input->contiguous = tmpInt;

				cout << "The maximum end stability (in default -9): ";
				cin >> tmpInt;
				input->endMinDG = tmpInt;

				cout << "The maximum hairpin (in default 3): ";
				cin >> tmpInt;
				input->maxHP = tmpInt;

				cout << endl;
			}
		}
		else if(!strcmp(argv[argnr], "-p2")){
			int tmpInt;
			if(stoi(argv[++argnr]) == 1){

				cout << "## Parameters for pair filteirng " << endl;
				cout << "The maximum length difference (in default 5): ";
				cin >> tmpInt;
				input->lenDiff = tmpInt;

				cout << "The maximum temperature difference (in default 3): ";
				cin >> tmpInt;
				input->TMDiff = tmpInt;

				cout << "The minimum PCR amplicon size (in default 100): ";
				cin >> tmpInt;
				input->minPS = tmpInt;

				cout << "The maximum PCR amplicon size (indefault 250): ";
				cin >> tmpInt;
				input->maxPS = tmpInt;

				cout << "The maximum pair-complementarity (in default 5): ";
				cin >> tmpInt;
				input->maxPC = tmpInt;

				cout << "The maximum 3'end pair-complementarity (in default 4): ";
				cin >> tmpInt;
				input->endMaxPC = tmpInt;
			}
		}
	}

	if(!inputCheck){
		cout << "There is no input sequence file. It is essential." << endl;
		return false;
	}

	else if(!outputCheck){
		cout << "There is no output path. It is essential." << endl;
		return false;
	}

	else if(!dirCheck){
		cout << "There is no storage directory. It is essential." << endl;
		return false;
	}

	input->printInputParameter();
	return true;
}

#endif /* INPUT_HPP_ */
