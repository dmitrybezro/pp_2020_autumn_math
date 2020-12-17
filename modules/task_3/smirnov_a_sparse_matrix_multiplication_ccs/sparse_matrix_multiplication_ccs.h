// Copyright 2018 Smirnov Alexander
#ifndef MODULES_TASK_3_SMIRNOV_A_SPARSE_MATRIX_MULTIPLICATION_CCS_SPARSE_MATRIX_MULTIPLICATION_CCS_H_
#define MODULES_TASK_3_SMIRNOV_A_SPARSE_MATRIX_MULTIPLICATION_CCS_SPARSE_MATRIX_MULTIPLICATION_CCS_H_

#include <vector>
#include <tuple>

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> getRandomMatrix(int sz);
std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
transpose(std::vector<double> values, std::vector<int> cols, std::vector<int> pointers);

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
getParallelOperations(int size, std::vector<double> Avalues, std::vector<int> Acols,std::vector<int> Apointers,
                                std::vector<double> Bvalues, std::vector<int> Bcols, std::vector<int> Bpointers);

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
getSequentialOperations(std::vector<double> Avalues, std::vector<int> Acols,std::vector<int> Apointers,
                        std::vector<double> Bvalues, std::vector<int> Bcols, std::vector<int> Bpointers);

#endif  // MODULES_TASK_3_SMIRNOV_A_SPARSE_MATRIX_MULTIPLICATION_CCS_SPARSE_MATRIX_MULTIPLICATION_CCS_H_
