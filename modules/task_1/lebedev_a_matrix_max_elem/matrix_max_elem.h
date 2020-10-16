// Copyright 2020 Lebedev Andrew
#ifndef MODULES_TASK_1_LEBEDEV_A_MATRIX_MAX_ELEM_MATRIX_MAX_ELEM_H_
#define MODULES_TASK_1_LEBEDEV_A_MATRIX_MAX_ELEM_MATRIX_MAX_ELEM_H_

#include <vector>
#include <string>

std::vector<int> getRandomMatrix(int  mtx_size);
int getParallelOperations(std::vector<int> global_mtx, int mtx_size);
int getSequentialOperations(std::vector<int> mtx);

#endif  // MODULES_TASK_1_LEBEDEV_A_MATRIX_MAX_ELEM_MATRIX_MAX_ELEM_H_
