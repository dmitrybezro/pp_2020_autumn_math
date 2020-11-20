// Copyright 2020 Vodeneev Mikhail
#ifndef MODULES_TASK_2_VODENEEV_M_BROADCAST_OPS_MPI_H_
#define MODULES_TASK_2_VODENEEV_M_BROADCAST_OPS_MPI_H_
#include <vector>
#include <string>

int Bcast(void* buffer, int count, MPI_Datatype datatype,
                          int root, MPI_Comm comm);

#endif  // MODULES_TASK_2_VODENEEV_M_BROADCAST_OPS_MPI_H_
