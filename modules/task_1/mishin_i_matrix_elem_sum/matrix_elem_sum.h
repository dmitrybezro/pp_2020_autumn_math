// Copyright 2020 Mishin Ilya
#ifndef MODULES_TASK_1_MISHIN_I_MATRIX_ELEM_SUM_MATRIX_ELEM_SUM_H_
#define MODULES_TASK_1_MISHIN_I_MATRIX_ELEM_SUM_MATRIX_ELEM_SUM_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int sz);
int getParallelOperations(std::vector<int> vec,
    int count_size_rows, int count_size_cols, std::string ops);
int getSequentialOperations(std::vector<int> vec, std::string ops);

#endif  // MODULES_TASK_1_MISHIN_I_MATRIX_ELEM_SUM_MATRIX_ELEM_SUM_H_
