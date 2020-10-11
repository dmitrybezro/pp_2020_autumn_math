// Copyright 2020 Mishin Ilya
#ifndef MATRIX_ELEM_SUM_H_
#define MATRIX_ELEM_SUM_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int sz);
int getParallelOperations(std::vector<int> vec,
    int count_size_rows, int count_size_cols, std::string ops);
int getSequentialOperations(std::vector<int> vec, std::string ops);

#endif  // MATRIX_ELEM_SUM_H_