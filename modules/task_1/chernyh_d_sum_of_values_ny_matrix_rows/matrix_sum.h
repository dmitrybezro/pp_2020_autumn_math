// Copyright 2018 Nesterov Alexander
#ifndef MODULES_TASK_1_CHERNYH_D_SUM_OF_VALUES_NY_MATRIX_ROWS_MATRIX_SUM_H_
#define MODULES_TASK_1_CHERNYH_D_SUM_OF_VALUES_NY_MATRIX_ROWS_MATRIX_SUM_H_

#include <vector>
#include <string>

std::vector<int> getRandomMat(int  count_row, int count_str);
std::vector<int> getParallelOperations(std::vector<int> global_mat,
  int count_row, int count_str);
std::vector<int> getSequentialOperations(std::vector<int> mat, int count_row, int count_str);

#endif  // MODULES_TASK_1_CHERNYH_D_SUM_OF_VALUES_NY_MATRIX_ROWS_MATRIX_SUM_H_
