// Copyright 2020 Vodeneev Mikhail
#ifndef MODULES_TASK_3_VODENEEV_M_MERGESORT_OPS_MPI_H_
#define MODULES_TASK_3_VODENEEV_M_MERGESORT_OPS_MPI_H_
#include <vector>
#include <string>

std::vector<int> merge(std::vector<int> a, std::vector<int> b);
int compare(const void* x1, const void* x2);
void quicksort(std::vector<int>* Arr, int first, int last);
std::vector<int> mergesort(std::vector<int> vec);

#endif  // MODULES_TASK_3_VODENEEV_M_MERGESORT_OPS_MPI_H_
