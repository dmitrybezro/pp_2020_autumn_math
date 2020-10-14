// Copyright 2020 Kolesnikov Gleb
#ifndef MODULES_TASK_1_KOLESNIKOV_G_APLHABETIC_COUNT_MPI_H_
#define MODULES_TASK_1_KOLESNIKOV_G_APLHABETIC_COUNT_MPI_H_

#include <vector>
#include <string>

std::vector<char> getRandomString(int  size);
int getParallelCount(std::vector<char> global_str,
	int vector_size);
int getSequentialCount(std::vector<char> str);

#endif  // MODULES_TASK_1_KOLESNIKOV_G_APLHABETIC_COUNT_MPI_H_