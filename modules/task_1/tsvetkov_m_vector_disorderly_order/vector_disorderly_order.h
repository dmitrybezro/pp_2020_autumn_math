  // Copyright 2020 Tsvetkov Maxim
#ifndef  MODULES_TASK_1_TSVETKOV_M_VECTOR_DISORDERLY_ORDER_VECTOR_DISORDERLY_ORDER_H_
#define MODULES_TASK_1_TSVETKOV_M_VECTOR_DISORDERLY_ORDER_VECTOR_DISORDERLY_ORDER_H_
#include <vector>
#include <string>

std::vector<int> getRandomVector(int sz);
int get_parallel_operations(std::vector<int> vec, int vec_size);
int get_sequential_operations(std::vector<int> vec);

#endif  // MODULES_TASK_1_TSVETKOV_M_VECTOR_DISORDERLY_ORDER_VECTOR_DISORDERLY_ORDER_H_
