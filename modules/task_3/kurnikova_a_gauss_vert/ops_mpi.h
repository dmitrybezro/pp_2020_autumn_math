// Copyright 2020 Kurnikova Anastasia
#ifndef MODULES_TASK_3_KURNIKOVA_A_GAUSS_VERT_OPS_MPI_H_
#define MODULES_TASK_3_KURNIKOVA_A_GAUSS_VERT_OPS_MPI_H_

#include <vector>

const double gauss[9] = { 1, 2, 1,
                          2, 4, 2,
                          1, 2, 1 };

std::vector<double> sequential(std::vector<double> image, int w, int h);
std::vector<double> parallel(std::vector<double> image, int n);

#endif  // MODULES_TASK_3_KURNIKOVA_A_GAUSS_VERT_OPS_MPI_H_
