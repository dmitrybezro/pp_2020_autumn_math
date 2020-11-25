// Copyright 2020 Lebedev Andrew
#ifndef MODULES_TASK_2_LEBEDEV_A_GATHER_MY_GATHER_H_
#define MODULES_TASK_2_LEBEDEV_A_GATHER_MY_GATHER_H_

#include <string>

int my_gather(void *sbuf, int scount, MPI_Datatype stype, void *rbuf,
    int rcount, MPI_Datatype rtype, int dest, MPI_Comm comm);

#endif  // MODULES_TASK_2_LEBEDEV_A_GATHER_MY_GATHER_H_
