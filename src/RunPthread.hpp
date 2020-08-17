/*
 * RunPthread.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: jmbae
 */

#ifndef RUNPTHREAD_HPP_
#define RUNPTHREAD_HPP_

#include "Global.hpp"

struct paramThread
{
	long rank;
	int myParam;
};

void runPthread(void* step(void* rank), int numOfPthreads, double* elapsedTime){
	long thread;
	pthread_t* thread_handles;
	thread_handles = new pthread_t[numOfPthreads];

	double s1, e1;

	s1 = getCurrentTime();

	for (thread = 0; thread < numOfPthreads; thread++) {
		pthread_create(&thread_handles[thread], NULL, step, (void*) thread);
	}
	for (thread = 0; thread < numOfPthreads; thread++) {
		pthread_join(thread_handles[thread], NULL);
	}

	del(thread_handles);

	e1 = getCurrentTime();

	(*elapsedTime) =  (e1 - s1);
}

void runPthreadParam(void* step(void* rank), int numOfPthreads,  int stage, double* elapsedTime){
	long thread;
	pthread_t* thread_handles;
	thread_handles = new pthread_t[numOfPthreads];

	paramThread** param = new paramThread*[numOfPthreads];
	for(thread = 0; thread < numOfPthreads; thread++){
		param[thread] = new paramThread();
		param[thread]->myParam = stage;
	}

	double s1, e1;

	s1 = getCurrentTime();

	for (thread = 0; thread < numOfPthreads; thread++) {
		param[thread]->rank = thread;
		pthread_create(&thread_handles[thread], NULL, step, (void*) param[thread]);
	}
	for (thread = 0; thread < numOfPthreads; thread++) {
		pthread_join(thread_handles[thread], NULL);
		del(param[thread]);
	}

	del(thread_handles);
	del(param);

	e1 = getCurrentTime();

	(*elapsedTime) =  (e1 - s1);
}



#endif /* RUNPTHREAD_HPP_ */
