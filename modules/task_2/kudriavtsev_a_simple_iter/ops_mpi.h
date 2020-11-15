// Copyright 2020 Kudriavtsev Alexander
#ifndef MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_OPS_MPI_H_
#define MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_OPS_MPI_H_

#include <vector>

std::vector<double> getRandomVector(int n);
std::vector<double> getRandomMatrix(int n);

std::vector<double> parallelIterMethod(std::vector<double> A, std::vector<double> b, int n, double eps, int NMax);
std::vector<double> sequentialIterMethod(std::vector<double> A, std::vector<double> b, int n, double eps, int NMax);

double discrepancyNorm(const std::vector<double>& x,
    const std::vector<double>& A, const std::vector<double>& b);

#endif  // MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_OPS_MPI_H_
