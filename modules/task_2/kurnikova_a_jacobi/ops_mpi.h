// Copyright 2020 Kurnikova Anastasia
#ifndef MODULES_TASK_2_KURNIKOVA_A_JACOBI_OPS_MPI_H_
#define MODULES_TASK_2_KURNIKOVA_A_JACOBI_OPS_MPI_H_

#include <math.h>
#include <vector>

int correctness(std::vector<double> a, int n);
int sequential(std::vector<double> a, std::vector<double> b,
                                   int n, const double e);
int parallel(std::vector<double> a, std::vector<double> b,
                                   int n, const double e);

#endif  // MODULES_TASK_2_KURNIKOVA_A_JACOBI_OPS_MPI_H_
