// Copyright 2020 Kolesnikov Gleb
#ifndef MODULES_TASK_3_LEBEDEV_A_ODD_EVEN_QUICK_SORT_ODD_EVEN_QUICK_SORT_H_
#define MODULES_TASK_3_LEBEDEV_A_ODD_EVEN_QUICK_SORT_ODD_EVEN_QUICK_SORT_H_

#include <mpi.h>
#include <vector>

std::vector<int> getRandomVector(int size);
int compare(const void *el1, const void *el2);
void batchersort(int* a, int l, int r);
int*odd_even_merge(int* a, int* b, int size1, int size2);
int parallel_qsort(int* global_vec, int vec_size);

#endif  // MODULES_TASK_3_LEBEDEV_A_ODD_EVEN_QUICK_SORT_ODD_EVEN_QUICK_SORT_H_
