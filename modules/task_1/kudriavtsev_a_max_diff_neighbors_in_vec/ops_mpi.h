// Copyright 2020 Kudriavtsev Alexander
#ifndef MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_MPI_H_
#define MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_MPI_H_

#include <vector>
//#include <string>

struct MyPair{
	int diff;
	int indx;
};

std::vector<int> getRandomVector(int n);
int getParallelOperations(std::vector<int> global_vec, int count_size_vector);
int getSequentialOperations(std::vector<int> vec, int count_size_vector);

#endif  // MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_MPI_H_