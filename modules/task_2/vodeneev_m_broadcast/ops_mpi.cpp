// Copyright 2020 Mikhail Vodeneev
#include <mpi.h>
#include <math.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include "../../../modules/test_tasks/test_mpi/ops_mpi.h"

int Bcast(void* buffer, int count,
              MPI_Datatype datatype, int root, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (root > size) {
        printf("Error");
        return 1;
    }


    if ((size == 1) && (root != 0)) {
        printf("Error");
        return 1;
    }


    if (root != 0) {
        if (rank == root)
           MPI_Send(buffer, count, datatype, 0, 0, MPI_COMM_WORLD);
        else if (rank == 0) {
            MPI_Status status;
            MPI_Recv(buffer, count, datatype, root, 0, MPI_COMM_WORLD, &status);
        }
    }

    if ((size == 1) && (root != 0)) {
        printf("Error");
        return 1;
    }

    if (size == 1) {
        return 0;
    }

    if (size == 2) {
        if (rank == 0) {
            MPI_Send(buffer, count, datatype, 1, 0, MPI_COMM_WORLD);
        } else {
            MPI_Status status;
            MPI_Recv(buffer, count, datatype, 0, 0, MPI_COMM_WORLD, &status);
        }
        return 0;
    }

    if (rank == 0) {
        MPI_Send(buffer, count, datatype, 1, 0, MPI_COMM_WORLD);
        MPI_Send(buffer, count, datatype, 2, 0, MPI_COMM_WORLD);
    } else {
        MPI_Status status;
        MPI_Recv(buffer, count, datatype,
                     (rank - 1) / 2, 0, MPI_COMM_WORLD, &status);
        if (2 * rank + 1 < size) {
            MPI_Send(buffer, count, datatype,
                      2 * rank + 1, 0, MPI_COMM_WORLD);
        }
        if (2 * rank + 2 < size) {
            MPI_Send(buffer, count, datatype, 2 * rank + 2, 0, MPI_COMM_WORLD);
        }
    }
    return 0;
}

