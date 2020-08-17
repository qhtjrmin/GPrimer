
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include "GPrimer.hpp"

int main(int argc, char *argv[]) {
	inputParameter* input = new inputParameter();
	double elapsedTime, st, et;

	if (readInputParameter(argc, argv, input) == false) {
		return -1;
	}

	initialization(input);

	cout << "## Stage 1 : candidate generation" << endl;
	st = getCurrentTime();
	runPthread(stage1PrimerGeneration, input->numOfThreads, &elapsedTime);
	stage1Sort();
	stage1FileDistribution();
	et = getCurrentTime();
	cout << "Elapsed time(Stage 1) : " << (et - st) << " sec" << endl << endl;

	sysCall("rm " + string(myInput->dirPath) + "/sorted.txt");

	cout << "\n## Stage 2 : single filtering" << endl;
	runPthread(stage2, input->numOfThreads, &elapsedTime);
	cout << "Elapsed time (Stage 2) : " << elapsedTime << " sec"  << endl << endl;

	cout << "## Stage 3 : 5' cross hybridization probing" << endl;
	runPthread(stage3, input->numOfThreads, &elapsedTime);
	stage3Delete();
	cout << "Elapsed time (Stage 3) : " << elapsedTime << " sec"  << endl << endl;

	if (input->isWriteOutput) {
		cout << "## Stage 3 (wirte) : 5' cross hybridization writing" << endl;
		runPthreadParam(writeOutput, input->numOfThreads, 3, &elapsedTime);
		cout << "Elapsed time (Stage 3_write) : " << elapsedTime << " sec" << endl << endl;
	}

#ifdef DEBUG
	primerHCheck();
#endif

	st = getCurrentTime();

	for (int k = 1; k <= 2; k++) {

		stage4Check(k);

		for(int i = 0; i < stage4Divide; i++){

		cout << "## Stage 4_1 (k = " << k
				<< ") : general cross hybridization building" << endl;
		runPthreadParam(stage4Building, input->numOfThreads, k, &elapsedTime);
		cout << "Elapsed time (Stage 4_1, k = " << k << ") : " << elapsedTime
				<< " sec" << endl << endl;

		cout << "## Stage 4_2 (k = " << k
				<< ") : general cross hybridization preparing" << endl;
		runPthreadParam(stage4Prepare, input->numOfThreads, k, &elapsedTime);
		cout << "Elapsed time (Stage 4_2, k = " << k << ") : " << elapsedTime
				<< " sec"  << endl << endl;

		cout << "## Stage 4_3 (k = " << k
				<< ") : general cross hybridization probing" << endl;
		runPthreadParam(stage4Probing, input->numOfThreads, k, &elapsedTime);
		cout << "Elapsed time (Stage 4_3, k = " << k << ") : " << elapsedTime
				<< " sec"  << endl << endl;

		cout << "## Stage 4_4 (k = " << k
				<< ") : general cross hybridization updating" << endl;
		runPthread(stage4Update, input->numOfThreads, &elapsedTime);
		cout << "Elapsed time (Stage 4_4, k = " << k << ") : " << elapsedTime
				<< " sec" << endl << endl;

		}

		if(input->isWriteOutput || k == 2){
			cout << "## Stage 4 (wirte) : general cross hybridization writing" << endl;
				runPthreadParam(writeOutput, input->numOfThreads, 3 + k, &elapsedTime);
				cout << "Elapsed time (Stage 4_write) : " << elapsedTime << " sec"
						<< endl << endl;
		}
#ifdef DEBUG
		primerHCheck();
#endif
	}
	stage4Final();
	et = getCurrentTime();

	cout << "Elapsed time (Stage 4 total) : " << (et - st) << " sec" << endl << endl;

	st = getCurrentTime();

	stage5Sort();
	sysCall("rm " + string(myInput->c4Path2) + "_*");

	finalFile = new FILE*[numOfThreads];
	for (long i = 0; i < numOfThreads; i++) {
		fileWriteOpen(&finalFile[i], foutName, i);
	}

	cout << "\n## Stage 5_1 : pair filtering & scoring preparing" << endl;
	runPthread(stage5Prepare, input->numOfThreads, &elapsedTime);
	cout << "Stage 5_preparing : " << elapsedTime << " sec." << endl << endl;

	cout << "## Stage 5_2 : pair filtering & scoring" << endl;
	runPthread(stage5, input->numOfThreads, &elapsedTime);
	cout << "Stage5_pairfilteing&scoring : " << elapsedTime
			<< " sec." << endl << endl;

	finalSort();
	sysCall("rm " + string(myInput->outputPath) + "_*");
	et = getCurrentTime();
	cout << "Elapsed time (Stage 5 total) : " << (et - st) << " sec" << endl << endl;


}
