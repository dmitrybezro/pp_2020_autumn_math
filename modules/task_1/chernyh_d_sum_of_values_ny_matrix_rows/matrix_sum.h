// Copyright 2018 Nesterov Alexander
#ifndef MODULES_TASK_1_NESTEROV_A_VECTOR_SUM_OPS_MPI_H_
#define MODULES_TASK_1_NESTEROV_A_VECTOR_SUM_OPS_MPI_H_

#include <vector>
#include <string>

std::vector<int> getRandomMat(int  count_row, int count_str);
std::vector<int> getParallelOperations(std::vector<int> global_mat,
	int count_row, int count_str);
std::vector<int> getSequentialOperations(std::vector<int> mat, int count_row, int count_str);

#endif  // MODULES_TASK_1_NESTEROV_A_VECTOR_SUM_OPS_MPI_H_
