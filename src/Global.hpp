/*
 * Global.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: jmbae
 */

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <pthread.h>
#include <unordered_map>
#include <set>

#include <sys/time.h>
#include <fstream>

using namespace std;

//elasped time checking
double getCurrentTime() {
	struct timeval curr;
	struct timezone tz;
	gettimeofday(&curr, &tz);
	double tmp = static_cast<double>(curr.tv_sec) * static_cast<double>(1000000)
			+ static_cast<double>(curr.tv_usec);
	return tmp * 1e-6;
}

unsigned int GetTickCount() {
	struct timeval gettick;
	unsigned int tick;
	gettimeofday(&gettick, NULL);

	tick = gettick.tv_sec * 1000 + gettick.tv_usec / 1000;

	return tick;
}

//data structures free
template<typename T>
inline void freeContainer(T& p_container) {
	using std::swap;
	p_container.clear();
	T().swap(p_container);
}

template<typename T>
inline void del(T ptr) {
	delete[] ptr;
	ptr = NULL;
}

//unordered_map comparison function
class CmpChar {
public:
	bool operator()(const char* t1, const char* t2) const {
		return (!strcmp(t1, t2));
	}

};

//set comparison function
class set_Cmp {
public:
	bool operator()(const char* t1, const char* t2) const {
		return (strcmp(t1, t2) < 0);
	}
};

static void HandleError(cudaError_t err, const char * file, int line){
	if(err != cudaSuccess){
		printf( "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(EXIT_FAILURE);
	}
}

//count the total number of rows of file
size_t countLine(istream &is){
    if(is.bad()) return 0;
    std::istream::iostate state_backup = is.rdstate();
    is.clear();
    std::istream::streampos pos_backup = is.tellg();
    is.seekg(0);
    size_t line_cnt;
    size_t lf_cnt = std::count(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>(), '\n');
    line_cnt = lf_cnt;
    return line_cnt;
}

//open the file for reading
void fileReadOpen(ifstream* fin, char* fName, long my_rank){
	char* buf = new char[300];
	if(my_rank > -1){
		sprintf(buf, "%s_%ld", fName, my_rank);
		fin->open(buf);
	}
	else{
		fin->open(fName);
	}
	if(!fin){
		cout << "Can't open the " << buf << "file\n";
	}
	del(buf);
}

//open the file for writing
void fileWriteOpen(FILE** fout, char* fName, long my_rank){
	char* buf = new char[300];
	if(my_rank > -1){
		sprintf(buf, "%s_%ld", fName, my_rank);
		(*fout) = fopen(buf, "w");
	}
	else
		(*fout) = fopen(fName, "w");
	if(!(*fout)){
		cout << "Can't open the " << buf << "file\n";
	}
	del(buf);
}

int sysCall(string command){
	cout << "system call: " << command << endl;
	return system(command.c_str());
}
#endif /* GLOBAL_HPP_ */
