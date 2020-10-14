
#ifndef MODULES_TASK_1_KOLESNIKOV_G_APLHABETIC_COUNT_MPI_H_
#define MODULES_TASK_1_KOLESNIKOV_G_APLHABETIC_COUNT_MPI_H_

#include <vector>
#include <string>

std::vector<char> getRandomString(int sz);
int getSequentialCount(std::vector<char> str);
int getParallelCount(std::vector<char> str);


#endif  