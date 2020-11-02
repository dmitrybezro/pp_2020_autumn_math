// Copyright 2020 Smirnov Alexander
#ifndef MODULES_TASK_1_SMIRNOV_A_ALTERNATION_OS_SIGNS_ALTERNATION_OS_SIGNS_H_
#define MODULES_TASK_1_SMIRNOV_A_ALTERNATION_OS_SIGNS_ALTERNATION_OS_SIGNS_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int  sz);
int getParallelOperations(std::vector<int> global_vec,
    int count_size_vector);
int getSequentialOperations(std::vector<int> vec);

#endif  // MODULES_TASK_1_SMIRNOV_A_ALTERNATION_OS_SIGNS_ALTERNATION_OS_SIGNS_H_
