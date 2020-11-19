// Copyright 2020 Kudriavtsev Alexander
#ifndef MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_OPS_MPI_H_
#define MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_OPS_MPI_H_

#include <vector>

std::vector<int> getRandomVector(int n, int mod = 100'000);
int* radixSort(int vec[], int n);
int* merge(int* vec1, int n1, int* vec2, int n2);
int* parallelRadixSort(int* vec, int n);
int* radixSort(int vec[], int n);

#endif  // MODULES_TASK_1_KUDRIAVTSEV_A_MAX_DIFF_NEIGHBORS_IN_VEC_OPS_MPI_H_
