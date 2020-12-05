// Copyright 2020 Kolesnikov Gleb
#ifndef MODULES_TASK_2_KOLESNIKOV_G_REDUCE_REDUCE_H_
#define MODULES_TASK_2_KOLESNIKOV_G_REDUCE_REDUCE_H_

#include <mpi.h>

int MPI_Reduce_My_Own(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
    MPI_Op op, int root, MPI_Comm comm);

#endif  // MODULES_TASK_2_KOLESNIKOV_G_REDUCE_REDUCE_H_
